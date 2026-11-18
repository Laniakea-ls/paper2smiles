import asyncio
import base64
import hashlib
import hmac
import json
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from starlette.testclient import TestClient

from paper2smiles import server
from paper2smiles._pdf import fitz


class ServiceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.cfg = {
            "jobs_root": self.temp.name,
            "identity_secret": "test-only-" + ("x" * 48),
            "api_key": "test-key",
            "audience": "paper2smiles",
            "max_upload_bytes": 100000,
            "max_pages": 5,
        }
        self.app = server.create_app(lambda: self.cfg)
        self.repository = self.app.state.repository
        self.runner = self.app.state.runner
        self.client = TestClient(self.app)
        with fitz.open() as pdf:
            pdf.new_page().insert_text((30, 40), "Scientific paper")
            self.pdf = pdf.tobytes()

    def tearDown(self):
        self.client.close()
        self.temp.cleanup()

    def headers(self, path, method="GET", user="user-a"):
        def enc(x):
            return (
                base64.urlsafe_b64encode(json.dumps(x).encode())
                .rstrip(b"=")
                .decode()
            )

        h = enc({"alg": "HS256", "typ": "JWT"})
        p = enc(
            {
                "iss": "aidd-portal",
                "aud": "paper2smiles",
                "sub": user,
                "iat": int(time.time()),
                "exp": int(time.time()) + 30,
                "path": path,
                "method": method,
            }
        )
        sig = (
            base64.urlsafe_b64encode(
                hmac.new(
                    self.cfg["identity_secret"].encode(),
                    (h + "." + p).encode(),
                    hashlib.sha256,
                ).digest()
            )
            .rstrip(b"=")
            .decode()
        )
        return {
            "Authorization": "Bearer " + h + "." + p + "." + sig,
            "Content-Type": "application/pdf",
        }

    def create(self, user="user-a"):
        r = self.client.post(
            "/jobs?filename=test.pdf",
            headers=self.headers("/jobs", "POST", user),
            content=self.pdf,
        )
        self.assertEqual(r.status_code, 202, r.text)
        return r.json()["id"]

    def test_identity_is_method_and_path_bound(self):
        self.assertEqual(self.client.get("/jobs").status_code, 401)
        self.assertEqual(
            self.client.get(
                "/jobs", headers=self.headers("/health")
            ).status_code,
            401,
        )
        self.assertEqual(
            self.client.post(
                "/jobs", headers=self.headers("/jobs"), content=self.pdf
            ).status_code,
            401,
        )

    def test_upload_dedup_owner_and_private_artifacts(self):
        jid = self.create()
        again = self.client.post(
            "/jobs?filename=test.pdf",
            headers=self.headers("/jobs", "POST"),
            content=self.pdf,
        )
        self.assertEqual(again.json()["id"], jid)
        self.assertEqual(
            self.client.get(
                "/jobs", headers=self.headers("/jobs", user="user-b")
            ).json()["jobs"],
            [],
        )
        path = "/jobs/" + jid
        self.assertEqual(
            self.client.get(
                path, headers=self.headers(path, user="user-b")
            ).status_code,
            404,
        )
        folder = Path(self.temp.name) / jid / "report"
        (folder / "downloads").mkdir(parents=True)
        (folder / "report-data.json").write_text("{}")
        (folder.parent / "secret.json").write_text("private")
        with self.repository.connect() as con:
            con.execute(
                "UPDATE jobs SET status='completed' WHERE id=?", (jid,)
            )
        artifact = path + "/files/report-data.json"
        self.assertEqual(
            self.client.get(
                artifact, headers=self.headers(artifact)
            ).status_code,
            200,
        )
        self.assertEqual(
            self.client.get(
                artifact, headers=self.headers(artifact, user="user-b")
            ).status_code,
            404,
        )
        private = path + "/files/secret.json"
        self.assertEqual(
            self.client.get(
                private, headers=self.headers(private)
            ).status_code,
            404,
        )

    def test_invalid_pdf_and_size(self):
        self.assertEqual(
            self.client.post(
                "/jobs?filename=bad.pdf",
                headers=self.headers("/jobs", "POST"),
                content=b"not a PDF",
            ).status_code,
            400,
        )
        self.assertEqual(
            self.client.post(
                "/jobs?filename=x.txt",
                headers=self.headers("/jobs", "POST"),
                content=self.pdf,
            ).status_code,
            400,
        )
        self.cfg["max_upload_bytes"] = 5
        self.assertEqual(
            self.client.post(
                "/jobs?filename=test.pdf",
                headers=self.headers("/jobs", "POST"),
                content=self.pdf,
            ).status_code,
            413,
        )

    def test_cancel_and_retry(self):
        jid = self.create()
        base = "/jobs/" + jid
        p = base + "/cancel"
        self.assertEqual(
            self.client.post(p, headers=self.headers(p, "POST")).status_code,
            200,
        )
        self.assertEqual(
            self.client.get(base, headers=self.headers(base)).json()["status"],
            "cancelled",
        )
        p = base + "/retry"
        self.assertEqual(
            self.client.post(p, headers=self.headers(p, "POST")).status_code,
            200,
        )
        self.assertEqual(
            self.client.get(base, headers=self.headers(base)).json()["status"],
            "queued",
        )

    def test_retry_waits_for_stopping_worker(self):
        jid = self.create()
        base = "/jobs/" + jid
        marker = Path(self.temp.name) / jid / ".worker-running"
        marker.touch()
        p = base + "/cancel"
        self.client.post(p, headers=self.headers(p, "POST"))
        p = base + "/retry"
        self.assertEqual(
            self.client.post(p, headers=self.headers(p, "POST")).status_code,
            409,
        )
        marker.unlink()
        self.assertEqual(
            self.client.post(p, headers=self.headers(p, "POST")).status_code,
            200,
        )

    def test_history_owner_status_and_report_actions(self):
        jid = self.create()
        folder = Path(self.temp.name) / jid
        report = folder / "report"
        (report / "downloads").mkdir(parents=True)
        (report / "downloads/Article title.xlsx").write_bytes(b"PK")
        (report / "report-data.json").write_text(
            json.dumps({"excel_download": "downloads/Article title.xlsx"})
        )
        (folder / "summary.json").write_text(
            json.dumps({"title": "Article title"})
        )
        with self.repository.connect() as con:
            con.execute(
                "UPDATE jobs SET status='completed' WHERE id=?", (jid,)
            )
        r = self.client.get("/history", headers=self.headers("/history"))
        self.assertEqual(r.status_code, 200)
        record = r.json()["records"][0]
        self.assertEqual(record["job_id"], "paper_sar_" + jid)
        self.assertEqual(record["summary"], "Article title")
        self.assertTrue(
            record["download_path"].endswith("Article%20title.xlsx")
        )
        self.assertIn("/reports/literature-ui/", record["view_path"])
        self.assertTrue(
            all(s["status"] == "completed" for s in record["steps"])
        )
        self.assertFalse(record["can_stop"])
        self.assertEqual(
            self.client.get(
                "/history?user=user-a&scope=all",
                headers=self.headers("/history", user="user-b"),
            ).json()["records"],
            [],
        )
        with self.repository.connect() as con:
            con.execute(
                "UPDATE jobs SET status='cancelled' WHERE id=?", (jid,)
            )
        record = self.client.get(
            "/history", headers=self.headers("/history")
        ).json()["records"][0]
        self.assertEqual(record["status"], "stopped")
        self.assertFalse(record["view_path"])

    def test_history_includes_records_beyond_recent_jobs(self):
        with self.repository.connect() as con:
            for i in range(105):
                con.execute(
                    "INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?)",
                    (
                        f"{i:032x}",
                        "user-a",
                        f"{i}.pdf",
                        str(i),
                        "failed",
                        i,
                        i,
                        "",
                    ),
                )
        self.assertEqual(
            len(
                self.client.get(
                    "/history", headers=self.headers("/history")
                ).json()["records"]
            ),
            105,
        )
        self.assertEqual(
            len(
                self.client.get("/jobs", headers=self.headers("/jobs")).json()[
                    "jobs"
                ]
            ),
            100,
        )

    def test_shutdown_terminates_active_subprocess(self):
        jid = self.create()
        spawn = asyncio.create_subprocess_exec
        children = []

        async def scenario():
            async def fake_spawn(*args, **kwargs):
                child = await spawn(
                    sys.executable,
                    "-c",
                    "import time;time.sleep(60)",
                    start_new_session=True,
                )
                children.append(child)
                return child

            with patch.object(
                asyncio,
                "create_subprocess_exec",
                side_effect=fake_spawn,
            ):
                task = asyncio.create_task(self.runner.run())
                for _ in range(100):
                    if children:
                        break
                    await asyncio.sleep(0.02)
                self.assertTrue(children)
                task.cancel()
                with self.assertRaises(asyncio.CancelledError):
                    await task
                self.assertIsNotNone(children[0].returncode)
                self.assertFalse(
                    (Path(self.temp.name) / jid / ".worker-running").exists()
                )

        asyncio.run(scenario())


if __name__ == "__main__":
    unittest.main()
