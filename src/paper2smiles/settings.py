"""Load private service configuration and resolve runtime storage."""

import json
import os
from pathlib import Path

ROOT = Path(os.environ.get("PAPER2SMILES_HOME", Path.cwd())).resolve()


def settings():
    """Load private service settings and apply environment overrides.

    Returns
    -------
    dict
        Resolved service, model, upload-limit, and storage settings.
    """
    path = Path(
        os.environ.get(
            "PAPER2SMILES_CONFIG", ROOT / "config/service.local.json"
        )
    )
    data = json.loads(path.read_text()) if path.exists() else {}
    data.setdefault("jobs_root", str(ROOT / "runs/jobs"))
    data.setdefault("model", "GPT-6-Luna")
    data.setdefault("api_base", "https://newapi.livzon.cn/v1")
    data["api_key"] = os.environ.get("PAPER_SI_API_KEY") or data.get(
        "api_key", ""
    )
    data.setdefault(
        "identity_secret", os.environ.get("PAPER2SMILES_IDENTITY_SECRET", "")
    )
    data.setdefault("audience", "paper2smiles")
    data.setdefault("max_upload_bytes", 50 * 1024 * 1024)
    data.setdefault("max_pages", 500)
    return data
