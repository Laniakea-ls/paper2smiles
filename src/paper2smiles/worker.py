"""Execute one queued job using private model configuration."""

import re
import sys
from pathlib import Path

from paper2smiles.extraction.activity_workflow import ActivityWorkflow
from paper2smiles.extraction.model import Luna
from paper2smiles.settings import settings


class ExtractionWorker:
    """Execute a queued extraction with private model configuration.

    Parameters
    ----------
    job_id : str
        Validated hexadecimal queue identifier.
    config : dict
        Runtime paths and model configuration.
    """

    def __init__(self, job_id, config):
        """Validate the job and construct its extraction workflow."""
        if not re.fullmatch("[a-f0-9]{32}", job_id):
            raise ValueError("Invalid job id")
        folder = Path(config["jobs_root"]) / job_id
        self.api_key = config["api_key"]
        if not self.api_key:
            raise RuntimeError("Model API key not configured")
        self.workflow = ActivityWorkflow(
            folder / "upload.pdf",
            folder,
            Luna(self.api_key, config["api_base"], config["model"], folder),
        )

    def run(self):
        """Execute the workflow and record failures without exposing keys.

        Returns
        -------
        int
            Zero on success, or one when extraction fails.
        """
        try:
            self.workflow.run()
        except Exception as exc:
            self.workflow.status(
                "提取未完成",
                failed=True,
                error=str(exc).replace(self.api_key, "[REDACTED]")[:1200],
            )
            return 1
        return 0


def main():
    """Run the queued extraction identified by the command-line job ID.

    Returns
    -------
    int
        Worker exit status.
    """
    return ExtractionWorker(sys.argv[1], settings()).run()


if __name__ == "__main__":
    raise SystemExit(main())
