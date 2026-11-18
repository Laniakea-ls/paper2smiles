"""Store job metadata and resolve the private job directory."""

import sqlite3
from pathlib import Path


class JobRepository:
    """Store job metadata and resolve the private job directory.

    Parameters
    ----------
    config_provider : callable
        Return the current service configuration without exposing credentials.
    """

    def __init__(self, config_provider):
        """Bind the service configuration provider."""
        self.settings = config_provider

    def directory(self):
        """Resolve the directory containing persistent jobs and their database.

        Returns
        -------
        pathlib.Path
            Absolute job-storage directory.
        """
        return Path(self.settings()["jobs_root"]).resolve()

    def connect(self):
        """Open the SQLite job index and initialize its schema if needed.

        Returns
        -------
        sqlite3.Connection
            Connection with named rows and write-ahead logging enabled.
        """
        self.directory().mkdir(parents=True, exist_ok=True)
        con = sqlite3.connect(self.directory() / "jobs.sqlite3", timeout=20)
        con.row_factory = sqlite3.Row
        con.execute("PRAGMA journal_mode=WAL")
        con.execute(
            (
                "CREATE TABLE IF NOT EXISTS jobs (id TEXT PRIMARY KEY,owner "
                "TEXT NOT NULL,filename TEXT NOT NULL,pdf_sha TEXT NOT "
                "NULL,status TEXT NOT NULL,created REAL NOT NULL,updated REAL "
                "NOT NULL,error TEXT NOT NULL DEFAULT '')"
            )
        )
        return con
