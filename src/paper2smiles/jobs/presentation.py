"""Present persisted jobs without exposing private worker state."""

import json

from paper2smiles.history import job_presentation

STAGES = {
    "plan": "识别文献",
    "activity_ocr": "提取活性表",
    "si_discovery": "查找补充材料",
    "names": "提取完整名称",
    "opsin": "名称转结构",
    "reference_names": "核对名称",
    "markush": "通式与片段连接",
    "full_structure": "结构图识别",
    "report": "生成报告",
    "completed": "完成",
}


MODEL_BILLING_ERROR = "模型服务返回计费错误（ACCOUNT_BALANCE_ERROR）。"


class JobPresenter:
    """Present persisted jobs without exposing private worker state.

    Parameters
    ----------
    repository : JobRepository
        Repository owning the task directories and persisted metadata.
    """

    def __init__(self, repository):
        """Bind the repository used by this component."""
        self.repository = repository

    def present(self, row):
        """Project a persisted job into its client-visible status and actions.

        Parameters
        ----------
        row : mapping
            Persisted SQLite job row, including its owner and status.

        Returns
        -------
        dict
            Job status, available summary, progress steps, and artifact URLs.
        """
        row = dict(row)
        jid = row["id"]
        folder = self.repository.directory() / jid
        value = {
            k: row[k]
            for k in (
                "id",
                "filename",
                "status",
                "created",
                "updated",
                "error",
            )
        }
        if value["error"] == "模型 API 余额不足，请联系管理员充值后重试。":
            value["error"] = MODEL_BILLING_ERROR
        status = folder / "status.json"
        if status.exists():
            try:
                status = json.loads(status.read_text())
                value["stage"] = status.get("stage", "")
                value["message"] = STAGES.get(
                    status.get("stage"), str(status.get("message", ""))[:160]
                )
            except (OSError, ValueError):
                pass
        summary = folder / "summary.json"
        if summary.exists():
            try:
                s = json.loads(summary.read_text())
                value["summary"] = {
                    k: s[k]
                    for k in [
                        "title",
                        "doi",
                        "compounds",
                        "with_structure",
                        "measurements",
                    ]
                    if k in s
                }
            except (OSError, ValueError):
                pass
        if row["status"] == "completed":
            value["report_path"] = (
                f"/api/paper-sar/jobs/{jid}/files/report-data.json"
            )
        value.update(job_presentation(value, folder))
        return value
