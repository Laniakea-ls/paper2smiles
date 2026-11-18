"""Adapt owned extraction jobs to the Portal history schema."""

import json
from datetime import datetime, timezone
from urllib.parse import quote

STEP_GROUPS = [
    ("plan", "读取文献 PDF", {"plan"}),
    ("activity_ocr", "识别活性编号与数值", {"activity_ocr"}),
    ("si_discovery", "查找作者 SI", {"si_discovery"}),
    ("names", "核对完整名称与结构", {"names", "opsin", "reference_names"}),
    ("markush", "通式连接与结构图识别", {"markush", "full_structure"}),
    ("report", "生成报告与活性数据总表", {"report", "completed"}),
]


def job_presentation(job, folder):
    """Build progress steps and validated report actions for one job.

    Parameters
    ----------
    job : dict
        Client-visible job fields and any available extraction summary.
    folder : pathlib.Path
        Directory for cached inputs, predictions, or generated files.

    Returns
    -------
    dict
        Step states and available same-origin report and workbook URLs.
    """
    status = {"cancelled": "stopped"}.get(job["status"], job["status"])
    stage = job.get("stage", "")
    index = next(
        (i for i, (_, _, stages) in enumerate(STEP_GROUPS) if stage in stages),
        0,
    )
    steps = []
    for i, (key, label, _) in enumerate(STEP_GROUPS):
        state = "pending"
        if status == "completed" or (status != "queued" and i < index):
            state = "completed"
        elif status != "queued" and i == index:
            state = status
        steps.append({"key": key, "label": label, "status": state})
    result = {"steps": steps, "report_url": "", "download_url": ""}
    if status != "completed":
        return result
    prefix = "/api/paper-sar/jobs/" + job["id"] + "/files/"
    result["report_url"] = (
        "/reports/literature-ui/literature-report.html?data="
        + quote(prefix + "report-data.json", safe="")
    )
    try:
        report = json.loads((folder / "report/report-data.json").read_text())
        relative = report.get("excel_download", "")
        base = (folder / "report").resolve()
        file = (base / relative).resolve()
        if (
            relative.startswith("downloads/")
            and file.is_relative_to(base)
            and file.is_file()
            and file.suffix == ".xlsx"
        ):
            result["download_url"] = prefix + quote(relative, safe="/")
    except (OSError, ValueError, TypeError):
        pass
    return result


def history_record(row, job):
    """Adapt an owned job to the Portal's unified history schema.

    Parameters
    ----------
    row : mapping
        Persisted SQLite job row, including its owner and status.
    job : dict
        Client-visible job fields and any available extraction summary.

    Returns
    -------
    dict
        Namespaced history row with owner, status, progress, and actions.
    """
    status = {"cancelled": "stopped"}.get(job["status"], job["status"])
    return {
        "job_id": "paper_sar_" + job["id"],
        "type": "paper_sar",
        "label": "文献活性提取",
        "icon": "file-text",
        "owner": row["owner"],
        "owner_name": row["owner"],
        "timestamp": datetime.fromtimestamp(
            job["created"], timezone.utc
        ).isoformat(),
        "summary": job.get("summary", {}).get("title") or job["filename"],
        "status": status,
        "status_detail": job["status"],
        "current_step": "排队中"
        if status == "queued"
        else job.get("message", ""),
        "steps": job["steps"],
        "error": job.get("error", ""),
        "can_stop": status in ("queued", "running"),
        "can_delete": False,
        "run_url": "/paper-sar?run=" + job["id"] + "&view=history",
        "report_url": job["report_url"],
        "view_path": job["report_url"] or None,
        "download_url": job["download_url"],
        "download_path": job["download_url"],
        "has_result_md": False,
        "result_md_path": None,
        "structures": None,
        "molecules": [],
    }
