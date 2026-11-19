"""Supervise queued extraction subprocesses and application shutdown."""

import asyncio
import json
import os
import signal
import sys
import time
from contextlib import asynccontextmanager

from paper2smiles.jobs.presentation import MODEL_BILLING_ERROR
from paper2smiles.settings import ROOT


class JobRunner:
    """Supervise queued extraction subprocesses and application shutdown.

    Parameters
    ----------
    repository : JobRepository
        Repository containing queued jobs and their checkpoint directories.
    """

    def __init__(self, repository):
        """Bind the repository used by this component."""
        self.repository = repository

    async def run(self):
        """Consume queued jobs and supervise one extraction subprocess."""
        for marker in self.repository.directory().glob("*/.worker-running"):
            marker.unlink()
        with self.repository.connect() as con:
            con.execute(
                "UPDATE jobs SET status='queued' WHERE status='running'"
            )
        process = None
        try:
            while True:
                with self.repository.connect() as con:
                    con.execute("BEGIN IMMEDIATE")
                    row = con.execute(
                        (
                            "SELECT * FROM jobs WHERE status='queued' ORDER "
                            "BY created LIMIT 1"
                        )
                    ).fetchone()
                    if row:
                        con.execute(
                            (
                                "UPDATE jobs SET status='running',updated=? "
                                "WHERE id=?"
                            ),
                            (time.time(), row["id"]),
                        )
                        (
                            self.repository.directory()
                            / row["id"]
                            / ".worker-running"
                        ).touch()
                if not row:
                    await asyncio.sleep(1)
                    continue
                jid = row["id"]
                folder = self.repository.directory() / jid
                try:
                    with (folder / "worker.log").open("ab") as log:
                        process = await asyncio.create_subprocess_exec(
                            sys.executable,
                            "-m",
                            "paper2smiles.worker",
                            jid,
                            cwd=ROOT,
                            stdout=log,
                            stderr=log,
                            start_new_session=True,
                        )
                        deadline = time.monotonic() + 12 * 3600
                        while process.returncode is None:
                            try:
                                await asyncio.wait_for(
                                    process.wait(), timeout=1
                                )
                            except asyncio.TimeoutError:
                                pass
                            with self.repository.connect() as con:
                                state = con.execute(
                                    "SELECT status FROM jobs WHERE id=?",
                                    (jid,),
                                ).fetchone()[0]
                            if (
                                state == "cancelled"
                                or time.monotonic() > deadline
                            ):
                                if process.returncode is None:
                                    os.killpg(process.pid, signal.SIGTERM)
                                try:
                                    await asyncio.wait_for(process.wait(), 10)
                                except asyncio.TimeoutError:
                                    os.killpg(process.pid, signal.SIGKILL)
                                    await process.wait()
                                break
                    with self.repository.connect() as con:
                        state = con.execute(
                            "SELECT status FROM jobs WHERE id=?", (jid,)
                        ).fetchone()[0]
                        if state != "cancelled":
                            completed = (
                                process.returncode == 0
                                and (
                                    folder / "report/report-data.json"
                                ).exists()
                            )
                            message = "提取未完成，可重试并继续检查点。"
                            try:
                                failure = json.loads(
                                    (folder / "status.json").read_text()
                                ).get("error", "")
                                if "ACCOUNT_BALANCE_ERROR" in failure:
                                    message = MODEL_BILLING_ERROR
                                elif "model_not_found" in failure:
                                    message = (
                                        "模型 API 暂无可用通道，请联系管理员。"
                                    )
                            except (OSError, ValueError):
                                pass
                            con.execute(
                                (
                                    "UPDATE jobs SET status=?,updated=?,error="
                                    "? WHERE id=?"
                                ),
                                (
                                    "completed" if completed else "failed",
                                    time.time(),
                                    "" if completed else message,
                                    jid,
                                ),
                            )
                except Exception:
                    with self.repository.connect() as con:
                        con.execute(
                            (
                                "UPDATE jobs SET status='failed',updated=?,err"
                                "or=? WHERE id=? AND status!='cancelled'"
                            ),
                            (time.time(), "任务执行失败，请重试。", jid),
                        )
                finally:
                    if process and process.returncode is None:
                        os.killpg(process.pid, signal.SIGTERM)
                        try:
                            await asyncio.wait_for(process.wait(), 10)
                        except asyncio.TimeoutError:
                            os.killpg(process.pid, signal.SIGKILL)
                            await process.wait()
                    process = None
                    (folder / ".worker-running").unlink(missing_ok=True)
        finally:
            if process and process.returncode is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    await asyncio.wait_for(process.wait(), 10)
                except asyncio.TimeoutError:
                    os.killpg(process.pid, signal.SIGKILL)
                    await process.wait()

    @asynccontextmanager
    async def lifespan(self, app):
        """Start the queue consumer and stop it when the application exits.

        Parameters
        ----------
        app : starlette.applications.Starlette
            Application whose background worker is supervised.

        Yields
        ------
        None
            Application lifespan context while the queue consumer is active.
        """
        task = asyncio.create_task(self.run())
        yield
        task.cancel()
        try:
            await task
        except asyncio.CancelledError:
            pass
