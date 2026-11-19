"""Expose owner-scoped extraction endpoints through the Portal gateway."""

import hashlib
import re
import shutil
import time
import uuid
from pathlib import Path
from urllib.parse import unquote

from starlette.exceptions import HTTPException
from starlette.responses import FileResponse, JSONResponse

from paper2smiles.history import history_record


class ExtractionAPI:
    """Expose owner-scoped extraction endpoints through the Portal gateway.

    Parameters
    ----------
    repository : JobRepository
        Storage for job metadata and source files.
    authenticator : GatewayAuthenticator
        Request-bound gateway identity verifier.
    presenter : JobPresenter
        Projection of persisted jobs into client-visible state.
    """

    def __init__(self, repository, authenticator, presenter):
        """Bind storage, authentication, and presentation dependencies."""
        self.repository = repository
        self.authenticator = authenticator
        self.presenter = presenter

    def get_job(self, request):
        """Look up a task belonging to the authenticated request owner.

        Parameters
        ----------
        request : starlette.requests.Request
            Incoming authenticated Portal gateway request.

        Returns
        -------
        dict
            Persisted job row, including its storage identifier and status.

        Raises
        ------
        starlette.exceptions.HTTPException
            If the task does not exist for the authenticated owner.
        """
        user = self.authenticator.owner(request)
        jid = request.path_params["job_id"]
        if not re.fullmatch("[a-f0-9]{32}", jid):
            raise HTTPException(404, "任务不存在")
        with self.repository.connect() as con:
            row = con.execute(
                "SELECT * FROM jobs WHERE id=? AND owner=?", (jid, user)
            ).fetchone()
        if not row:
            raise HTTPException(404, "任务不存在")
        return dict(row)

    async def history(self, request):
        """List all literature jobs owned by the authenticated caller.

        Parameters
        ----------
        request : starlette.requests.Request
            Incoming authenticated Portal gateway request.

        Returns
        -------
        starlette.responses.JSONResponse
            Owner-scoped records in the Portal history schema.
        """
        user = self.authenticator.owner(request)
        with self.repository.connect() as con:
            rows = con.execute(
                "SELECT * FROM jobs WHERE owner=? ORDER BY created DESC",
                (user,),
            ).fetchall()
        return JSONResponse(
            {
                "records": [
                    history_record(r, self.presenter.present(r)) for r in rows
                ]
            },
            headers={"Cache-Control": "private, no-store"},
        )

    async def health(self, request):
        """Report service configuration status without disclosing credentials.

        Parameters
        ----------
        request : starlette.requests.Request
            Incoming authenticated Portal gateway request.

        Returns
        -------
        starlette.responses.JSONResponse
            Service version, upload limit, and model-configuration flag.
        """
        self.authenticator.owner(request)
        cfg = self.repository.settings()
        return JSONResponse(
            {
                "ok": True,
                "name": "文献活性提取",
                "configured": bool(cfg["api_key"]),
                "version": "1.0.0",
                "max_upload_bytes": cfg["max_upload_bytes"],
            }
        )

    async def jobs(self, request):
        """List recent jobs or validate and enqueue an uploaded PDF.

        Parameters
        ----------
        request : starlette.requests.Request
            Incoming authenticated Portal gateway request.

        Returns
        -------
        starlette.responses.JSONResponse
            Recent jobs, an existing matching job, or the new queued job.
        """
        user = self.authenticator.owner(request)
        if request.method == "GET":
            with self.repository.connect() as con:
                rows = con.execute(
                    (
                        "SELECT * FROM jobs WHERE owner=? ORDER BY created "
                        "DESC LIMIT 100"
                    ),
                    (user,),
                ).fetchall()
            return JSONResponse(
                {"jobs": [self.presenter.present(r) for r in rows]}
            )
        cfg = self.repository.settings()
        if not cfg["api_key"]:
            raise HTTPException(503, "模型服务尚未配置")
        filename = Path(
            unquote(request.query_params.get("filename", "paper.pdf"))
        ).name[:240]
        if not filename.lower().endswith(".pdf"):
            raise HTTPException(400, "请上传 PDF 文件")
        jid = uuid.uuid4().hex
        folder = self.repository.directory() / jid
        folder.mkdir(parents=True)
        pdf = folder / "upload.pdf"
        size = 0
        digest = hashlib.sha256()
        try:
            with pdf.open("wb") as stream:
                async for chunk in request.stream():
                    size += len(chunk)
                    if size > cfg["max_upload_bytes"]:
                        raise HTTPException(413, "PDF 不能超过 50 MB")
                    digest.update(chunk)
                    stream.write(chunk)
            with pdf.open("rb") as stream:
                signature = stream.read(5)
            if size < 5 or signature != b"%PDF-":
                raise HTTPException(400, "文件不是有效 PDF")
            from paper2smiles._pdf import fitz

            try:
                with fitz.open(pdf) as document:
                    if document.needs_pass:
                        raise HTTPException(400, "暂不支持加密 PDF")
                    if not 0 < len(document) <= cfg["max_pages"]:
                        raise HTTPException(400, "PDF 页数超出限制")
            except (RuntimeError, ValueError):
                raise HTTPException(400, "无法读取 PDF")
            with self.repository.connect() as con:
                con.execute("BEGIN IMMEDIATE")
                existing = con.execute(
                    (
                        "SELECT * FROM jobs WHERE owner=? AND pdf_sha=? AND "
                        "status IN ('queued','running','completed') ORDER BY "
                        "created DESC LIMIT 1"
                    ),
                    (user, digest.hexdigest()),
                ).fetchone()
                if existing:
                    pdf.unlink()
                    folder.rmdir()
                    return JSONResponse(self.presenter.present(existing))
                now = time.time()
                con.execute(
                    "INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?)",
                    (
                        jid,
                        user,
                        filename,
                        digest.hexdigest(),
                        "queued",
                        now,
                        now,
                        "",
                    ),
                )
            return JSONResponse(
                self.presenter.present(
                    {
                        "id": jid,
                        "filename": filename,
                        "status": "queued",
                        "created": now,
                        "updated": now,
                        "error": "",
                    }
                ),
                202,
            )
        except BaseException:
            if pdf.exists():
                pdf.unlink()
            if folder.exists() and (not any(folder.iterdir())):
                folder.rmdir()
            raise

    async def detail(self, request):
        """Return the status and available results of an owned job.

        Parameters
        ----------
        request : starlette.requests.Request
            Incoming authenticated Portal gateway request.

        Returns
        -------
        starlette.responses.JSONResponse
            Client-visible details for the requested job.
        """
        return JSONResponse(self.presenter.present(self.get_job(request)))

    async def delete(self, request):
        """Delete an owned terminal job and its private artifacts."""
        row = self.get_job(request)
        folder = self.repository.directory() / row["id"]
        with self.repository.connect() as con:
            con.execute("BEGIN IMMEDIATE")
            removed = con.execute(
                "DELETE FROM jobs WHERE id=? AND owner=? "
                "AND status IN ('completed','failed','cancelled')",
                (row["id"], row["owner"]),
            )
            if not removed.rowcount or (folder / ".worker-running").exists():
                raise HTTPException(409, "任务仍在运行或停止中，请稍后再试")
            if folder.exists():
                shutil.rmtree(folder)
        return JSONResponse({"deleted": True})

    async def cancel(self, request):
        """Mark an owned queued or running job for cancellation.

        Parameters
        ----------
        request : starlette.requests.Request
            Incoming authenticated Portal gateway request.

        Returns
        -------
        starlette.responses.JSONResponse
            Acknowledgement of the cancellation request.
        """
        row = self.get_job(request)
        with self.repository.connect() as con:
            con.execute(
                (
                    "UPDATE jobs SET status='cancelled',updated=? WHERE id=? "
                    "AND status IN ('queued','running')"
                ),
                (time.time(), row["id"]),
            )
        return JSONResponse({"ok": True})

    async def retry(self, request):
        """Requeue an owned terminal job after its worker has stopped.

        Parameters
        ----------
        request : starlette.requests.Request
            Incoming authenticated Portal gateway request.

        Returns
        -------
        starlette.responses.JSONResponse
            Acknowledgement that the job can resume its checkpoints.

        Raises
        ------
        starlette.exceptions.HTTPException
            If the job is active or its previous worker is still stopping.
        """
        row = self.get_job(request)
        if row["status"] not in ("failed", "cancelled"):
            raise HTTPException(409, "当前任务不能重试")
        with self.repository.connect() as con:
            con.execute("BEGIN IMMEDIATE")
            if (
                self.repository.directory() / row["id"] / ".worker-running"
            ).exists():
                raise HTTPException(409, "任务正在停止，请稍后重试")
            con.execute(
                (
                    "UPDATE jobs SET status='queued',error='',updated=? WHERE "
                    "id=? AND status IN ('failed','cancelled')"
                ),
                (time.time(), row["id"]),
            )
        return JSONResponse({"ok": True})

    async def artifact(self, request):
        """Serve a completed job artifact within the owner's report directory.

        Parameters
        ----------
        request : starlette.requests.Request
            Incoming authenticated Portal gateway request.

        Returns
        -------
        starlette.responses.FileResponse
            Private report data, evidence, or a downloadable workbook.
        """
        row = self.get_job(request)
        if row["status"] != "completed":
            raise HTTPException(409, "报告尚未完成")
        base = (self.repository.directory() / row["id"] / "report").resolve()
        relative = request.path_params["path"]
        path = (base / relative).resolve()
        allowed = relative == "report-data.json" or relative.startswith(
            ("assets/", "downloads/")
        )
        if (
            not allowed
            or not path.is_relative_to(base)
            or (not path.is_file())
            or (
                path.suffix.lower()
                not in {
                    ".json",
                    ".pdf",
                    ".xlsx",
                    ".csv",
                    ".tsv",
                    ".png",
                    ".svg",
                    ".jpg",
                }
            )
        ):
            raise HTTPException(404, "文件不存在")
        headers = {
            "Cache-Control": "private, no-store",
            "X-Content-Type-Options": "nosniff",
        }
        filename = (
            path.name
            if path.suffix.lower() in {".xlsx", ".csv", ".tsv"}
            else None
        )
        return FileResponse(path, filename=filename, headers=headers)

    async def error(self, request, exc):
        """Serialize an HTTP failure using the service error schema.

        Parameters
        ----------
        request : starlette.requests.Request
            Request associated with the HTTP failure.
        exc : starlette.exceptions.HTTPException
            HTTP failure to serialize for the client.

        Returns
        -------
        starlette.responses.JSONResponse
            Error detail with the original HTTP status code.
        """
        return JSONResponse({"error": exc.detail}, status_code=exc.status_code)
