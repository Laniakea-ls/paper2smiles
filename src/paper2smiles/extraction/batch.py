"""Run an inventoried PDF batch with resumable per-paper evidence."""

import argparse
import getpass
import json
import time
import traceback
from pathlib import Path

from .activity_workflow import ActivityWorkflow
from .core import save
from .model import Luna
from .workflow import Workflow


def main():
    """Run the inventoried PDF batch with resumable per-paper checkpoints."""
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--ids", nargs="*")
    p.add_argument("--retry-failed", action="store_true")
    p.add_argument(
        "--v3",
        action="store_true",
        help=(
            "Compatibility flag: SI/name/scaffold workflow is now the default"
        ),
    )
    p.add_argument("--legacy-si-first", action="store_true")
    a = p.parse_args()
    root = a.root.resolve()
    key = getpass.getpass("API key (hidden, memory only): ")
    inventory = json.loads((root / "inventory.json").read_text())
    results = []
    for item in inventory:
        if (
            item["kind"] != "paper"
            or item.get("duplicate_of")
            or (a.ids and item["id"] not in a.ids)
        ):
            continue
        dest = root / item["id"]
        dest.mkdir(exist_ok=True)
        if (dest / "batch_result.json").exists() and not a.retry_failed:
            results.append(
                json.loads((dest / "batch_result.json").read_text())
            )
            continue
        print("START", item["id"], item["filename"], flush=True)
        started = time.time()
        try:
            llm = Luna(key, "https://newapi.livzon.cn/v1", "GPT-6-Luna", dest)
            (Workflow if a.legacy_si_first else ActivityWorkflow)(
                item["path"], dest, llm
            ).run()
            result = {
                "id": item["id"],
                "filename": item["filename"],
                "status": "completed",
                "summary": json.loads((dest / "summary.json").read_text()),
            }
        except Exception as e:
            result = {
                "id": item["id"],
                "filename": item["filename"],
                "status": "failed",
                "error": str(e).replace(key, "[REDACTED]"),
            }
            (dest / "failure.txt").write_text(
                traceback.format_exc().replace(key, "[REDACTED]")
            )
        result["seconds"] = round(time.time() - started, 1)
        save(dest / "batch_result.json", result)
        results.append(result)
        save(root / "batch_results.json", results)
        print("RESULT", json.dumps(result, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
