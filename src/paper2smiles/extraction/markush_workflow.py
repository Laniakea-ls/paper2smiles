"""Source-image-only scaffold/fragment route, independent of SI and names."""

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

from rdkit import Chem

from paper2smiles.settings import ROOT

from .core import fitz, normalize_id, save, sha256
from .markush_chemistry import assemble

VERSION = "markush-v1"


def reviewed_examples(pdf_hash, targets):
    # Reuse an audited local trial only for the exact source PDF, with
    # immutable source provenance.
    """Load reviewed single-site assemblies only for an exact PDF match.

    Parameters
    ----------
    pdf_hash : str
        SHA-256 digest used to require an exact source PDF match.
    targets : set of str
        Active compound identifiers eligible for structure resolution.

    Returns
    -------
    list of dict
        Reviewed reconstruction records for the requested identifiers.
    """
    trial = ROOT / "runs/gpr88_table1_reconstruction"
    if not (trial / "summary.json").exists():
        return []
    summary = json.loads((trial / "summary.json").read_text())
    if summary["source_pdf_sha256"] != pdf_hash:
        return []
    out = []
    for r in json.loads((trial / "assembled_records.json").read_text()):
        if r["compound_id"] not in targets:
            continue
        out.append(
            {
                **r,
                "page": 2,
                "table_id": "page-002-t1",
                "review_status": "source_visually_reviewed",
                "review_note": "已完成的 Table 1 逐图复核；不是独立金标准。",
                "origin": "verified_local_trial",
                "source_pdf_sha256": pdf_hash,
            }
        )
    return out


def render_markush_crops(pdf, folder, page, layout, ids):
    """Render scaffold and fragment crops with the SAR cropping adapter.

    Parameters
    ----------
    pdf : pathlib.Path or str
        Original literature PDF used as source evidence.
    folder : pathlib.Path
        Directory for cached inputs, predictions, or generated files.
    page : int
        One-based PDF page number to render.
    layout : dict
        Model-located scaffold and fragment bounds on the source page.
    ids : set of str
        Allowed compound identifiers for the requested source table.

    Returns
    -------
    tuple
        Crop paths, accepted fragment records, and scaffold PDF bounds.
    """
    from .full_structure_workflow import SAR, render_with_edge_check, valid_box

    sys.path.insert(0, str(SAR))
    from pipeline.smiles_validation.ocsr_render import render_inputs

    box = layout.get("scaffold_bbox")
    if not valid_box(box):
        raise ValueError("Invalid scaffold crop")
    rows = [
        {
            "compound_id": "scaffold",
            "crop": "scaffold.png",
            "page": page,
            "bbox_pixels": list(box),
            "bbox_page_size": [1000, 1000],
        }
    ]
    accepted = []
    seen = set()
    for item in layout.get("rows", []):
        cid = normalize_id(item.get("compound_id"))
        box = item.get("fragment_bbox")
        if cid not in ids or cid in seen:
            continue
        if not valid_box(box):
            raise ValueError("Invalid fragment crop: " + cid)
        seen.add(cid)
        accepted.append({**item, "compound_id": cid})
        rows.append(
            {
                "compound_id": cid,
                "crop": "fragment-" + str(len(accepted) - 1) + ".png",
                "page": page,
                "bbox_pixels": list(box),
                "bbox_page_size": [1000, 1000],
            }
        )
    # Keep nearby long table rules outside the initial padding; they are
    # not molecular bonds.
    with fitz.open(pdf) as doc:
        source = doc[page - 1]
        height = source.rect.height
        width = source.rect.width
        for row in rows:
            x, y, u, v = row["bbox_pixels"]
            left, right = x * width / 1000, u * width / 1000
            for drawing in source.get_drawings():
                for item in drawing["items"]:
                    if item[0] == "re" and item[1].height <= 1:
                        a, b = item[1].tl, item[1].tr
                    elif item[0] == "l":
                        a, b = item[1:3]
                    else:
                        continue
                    if abs(a.y - b.y) > 0.1 or abs(a.x - b.x) < max(
                        50, (right - left) * 1.5
                    ):
                        continue
                    if min(a.x, b.x) > left or max(a.x, b.x) < right:
                        continue
                    if abs(a.y - v * height / 1000) < 3:
                        v = (a.y - 2) * 1000 / height
                    if abs(a.y - y * height / 1000) < 3:
                        y = (a.y + 2) * 1000 / height
            row["bbox_pixels"] = [x, y, u, v]
    mapping, audit = render_with_edge_check(pdf, rows, folder, render_inputs)
    crops = [Path(mapping[str((folder / r["crop"]).resolve())]) for r in rows]
    for item, entry in zip(accepted, audit[1:]):
        item["bbox_pdf_points"] = entry["pdf_bbox"]
    save(
        folder / "sar_render_manifest.json",
        {
            "renderer": "SAR ocsr_render.render_inputs",
            "dpi": 300,
            "crops": [
                {**r, **a, "crop_path": str(p), "sha256": sha256(p)}
                for r, a, p in zip(rows, audit, crops)
            ],
        },
    )
    return crops, accepted, audit[0]["pdf_bbox"]


def run(flow, plan, body):
    """Locate, reconstruct, and audit eligible scaffold-fragment tables.

    Parameters
    ----------
    flow : ActivityWorkflow
        Active extraction context with source pages and model access.
    plan : dict
        Document plan containing article metadata and relevant source pages.
    body : list of dict
        Activity measurements extracted from the main article.

    Returns
    -------
    list of dict
        Candidate and reviewed reconstructions with source evidence.
    """
    root = flow.root / "markush"
    root.mkdir(exist_ok=True)
    digest = sha256(flow.pdf)
    targets = set(flow.target_ids)
    from .multisite_reviewed import reviewed_multisite

    rows = reviewed_examples(digest, targets) + reviewed_multisite(
        digest, body, targets
    )
    reviewed_path = root / "source_reviewed_records.json"
    if reviewed_path.exists():
        reviewed = json.loads(reviewed_path.read_text())
        if reviewed.get("source_pdf_sha256") == digest:
            rows = list(
                {
                    (r["table_id"], r["compound_id"]): r
                    for r in rows
                    + [
                        r
                        for r in reviewed["records"]
                        if r["compound_id"] in targets
                        and r.get("review_status")
                        == "source_visually_reviewed"
                        and r.get("source_pdf_sha256") == digest
                    ]
                }.values()
            )
    reviewed_tables = {r["table_id"] for r in rows}
    source_tables = (
        json.loads((flow.root / "activity_ocr/tables.json").read_text())
        if (flow.root / "activity_ocr/tables.json").exists()
        else []
    )
    audit = []
    for t in source_tables:
        tid = t["table_id"]
        page = t["page"]
        ids = sorted(
            {r["compound_id"] for r in body if r.get("table_id") == tid}
            & targets
        )
        if not ids:
            continue
        if tid in reviewed_tables:
            audit.append({"table_id": tid, "status": "reused_reviewed_trial"})
            continue
        folder = root / tid
        folder.mkdir(exist_ok=True)
        event = {"table_id": tid, "page": page, "status": "pending"}
        try:
            prompt = (
                "Locate ONLY the source-image geometry for this "
                "bioactivity table. No SMILES, names, invented "
                "structures, or activity values. Find the generic "
                "scaffold belonging to the specified table, and each row "
                "fragment. The scaffold may be ABOVE the OCR table box. "
                "All bounding boxes use the FULL supplied PDF page "
                "normalized to [0,1000], [left,top,right,bottom]. Crop "
                "each complete molecular drawing with whitespace, "
                "including its attachment bond and R label but excluding "
                'other rows/text. Return {"supported":true,"table_label":"'
                'Table 1","variable_label":"R1","scaffold_bbox":[0,0,0,0],'
                '"rows":[{"compound_id":"4","fragment_bbox":[0,0,0,0],"ide'
                'ntifier_evidence":"literal printed row '
                'label","stereo_note":"literal stereochemistry evidence '
                'or unspecified"}],"reason":"source-based explanation"}. '
                "Support ONLY one variable attachment position and "
                "image-drawn substituents. Set supported=false for "
                "full-molecule tables, multiple variables, no visible "
                "shared scaffold, text-only substituents, or ambiguous "
                "table/row mapping. Never guess missing structures. Only "
                "return the requested IDs. The task is image "
                "localization, not molecular interpretation.\n"
            ) + json.dumps(
                {
                    "table_id": tid,
                    "page": page,
                    "caption": t["schema"].get("caption"),
                    "compound_ids": ids,
                }
            )
            layout = flow.llm.ask(
                VERSION + "_layout_" + tid,
                prompt,
                [flow.image(page)],
                max_tokens=4500,
            )
            save(folder / "layout.json", layout)
            if layout.get("supported") is not True:
                event.update(status="unsupported", reason=layout.get("reason"))
                graph_layout = flow.llm.ask(
                    VERSION + "_graph_columns_" + tid,
                    (
                        "Inspect this original table for MULTIPLE image-drawn "
                        "substituent columns (R, Ar1, Ar2, etc.) and a "
                        'visible ruled grid. Return {"supported":true,"columns'
                        '":{"1":"R","2":"Ar1","3":"Ar2"},"reason":"..."}. '
                        "Column indexes are zero-based INCLUDING the "
                        "compound-number column. Only identify actual "
                        "fragment column headers. No SMILES or guessed cell "
                        "inheritance. Return supported=false for full "
                        "molecules, unruled tables, or unclear headers. "
                        "Table: "
                    )
                    + tid,
                    [flow.image(page)],
                    max_tokens=1500,
                )
                save(folder / "graph_layout.json", graph_layout)
                if graph_layout.get("supported") is True:
                    from .fragment_graph_workflow import prepare

                    entries = json.loads(
                        (
                            flow.root / "activity_ocr/table_crops/inputs.json"
                        ).read_text()
                    )
                    entry = next(e for e in entries if e["id"] == tid)
                    ordered = list(
                        dict.fromkeys(
                            r["compound_id"]
                            for r in sorted(
                                body, key=lambda r: r.get("source_row", 0)
                            )
                            if r.get("table_id") == tid
                        )
                    )
                    columns = {
                        int(k): v for k, v in graph_layout["columns"].items()
                    }
                    if (
                        not columns
                        or any(k < 1 for k in columns)
                        or len(set(columns.values())) != len(columns)
                    ):
                        raise ValueError("Invalid fragment column schema")
                    event.update(
                        prepare(
                            flow.pdf,
                            folder / "fragment_graph",
                            entry,
                            columns,
                            ordered,
                        )
                    )
                audit.append(event)
                continue
            crops, accepted, core_box = render_markush_crops(
                flow.pdf, folder, page, layout, ids
            )
            if not accepted:
                raise ValueError("No allowed row mapping")
            signature = hashlib.sha256(
                "".join(sha256(p) for p in crops).encode()
            ).hexdigest()[:16]
            result = folder / ("molparser_" + signature + ".json")
            if not result.exists():
                env = os.environ.copy()
                env["HF_HOME"] = str(root / "hf_cache")
                cmd = [
                    "/data2/aidd/aidd_platform/tools/MolParser-Mobile/run.sh",
                    *[str(p) for p in crops],
                    "--device",
                    "cpu",
                    "--output",
                    str(result),
                ]
                proc = subprocess.run(
                    cmd, env=env, capture_output=True, text=True, timeout=240
                )
                (folder / "molparser.log").write_text(
                    proc.stdout + proc.stderr
                )
                if proc.returncode:
                    raise RuntimeError(
                        "Local MolParser failed; see per-table log"
                    )
            recognized = json.loads(result.read_text())
            if len(recognized) != len(crops):
                raise ValueError("Recognition row count mismatch")
            core_smi = recognized[0].get("smi", "")
            emitted = []
            for item, rec in zip(accepted, recognized[1:]):
                try:
                    m = assemble(core_smi, rec.get("smi", ""))
                    smi = Chem.MolToSmiles(m, isomericSmiles=True)
                    r = {
                        "compound_id": item["compound_id"],
                        "scaffold_smiles": core_smi,
                        "r1_smiles": rec["smi"],
                        "assembled_smiles": smi,
                        "page": page,
                        "table_id": tid,
                        "variable_label": layout.get("variable_label", "R"),
                        "source_pdf_sha256": digest,
                        "source_type": "main_pdf_scaffold_r1_reconstruction",
                        "review_status": "needs_visual_review",
                        ("review_note"): (
                            "自动定位与识图重建候选，连接点及立体化学待源图复"
                            "核。"
                        ),
                        "stereo_note": item.get("stereo_note", ""),
                        "unspecified_stereocenters": sum(
                            v == "?"
                            for _, v in Chem.FindMolChiralCenters(
                                m, includeUnassigned=True
                            )
                        ),
                        "scaffold_bbox_pdf_points": core_box,
                        "fragment_bbox_pdf_points": item["bbox_pdf_points"],
                        "recognition_file": str(result),
                        "crop_renderer": "SAR ocsr_render.render_inputs",
                        "crop_dpi": 300,
                        "crop_manifest": str(
                            folder / "sar_render_manifest.json"
                        ),
                    }
                    rows.append(r)
                    emitted.append(item["compound_id"])
                except Exception as exc:
                    event.setdefault("row_errors", []).append(
                        {"compound_id": item["compound_id"], "error": str(exc)}
                    )
            if event.get("row_errors"):
                flow.warnings.append(
                    "通式表 "
                    + tid
                    + " 有 "
                    + str(len(event["row_errors"]))
                    + " 行因连接点不明确未生成结构。"
                )
            event.update(status="completed", emitted_ids=emitted)
        except Exception as exc:
            event.update(status="failed", reason=str(exc))
            flow.warnings.append("通式重建 " + tid + " 未完成：" + str(exc))
        audit.append(event)
        save(root / "audit.json", audit)
    save(root / "audit.json", audit)
    save(flow.root / "markush_records.json", rows)
    return rows


def attach(cards, records):
    """Attach reconstruction evidence without promoting unreviewed candidates.

    Parameters
    ----------
    cards : list of dict
        Compound cards containing source-specific structure evidence.
    records : list of dict
        Source evidence records to audit, update, or attach.

    Returns
    -------
    list of dict
        Updated compound cards with reconstruction evidence and status.
    """
    byid = {}
    for r in records:
        byid.setdefault(r["compound_id"], []).append(r)
    for c in cards:
        c["markush"] = byid.get(c["compound_id"], [])
        c["structure_routes"] = {
            "si": bool(c.get("si_smiles")),
            "name": any(
                n["opsin"]["valid"]
                and n.get("eligible_for_structure_comparison", True)
                and (
                    n.get("name_text_supported")
                    or n.get("visual_name_review", {}).get("supported")
                )
                for n in c.get("names", [])
            ),
            "markush_candidate": bool(c["markush"]),
        }
        for r in c["markush"]:
            candidate = Chem.MolFromSmiles(r.get("assembled_smiles", ""))
            if candidate:
                smi = Chem.MolToSmiles(candidate, isomericSmiles=True)
                r["si_comparison"] = (
                    (
                        "一致"
                        if smi
                        == Chem.MolToSmiles(
                            Chem.MolFromSmiles(c["si_smiles"]),
                            isomericSmiles=True,
                        )
                        else "不同，需核对"
                    )
                    if c.get("si_smiles")
                    else "无 SI 结构可对照"
                )
        reviewed = any(
            r.get("review_status") == "source_visually_reviewed"
            for r in c["markush"]
        )
        c["has_structure"] = bool(c.get("has_structure") or reviewed)
        if reviewed and not (
            c["structure_routes"]["si"] or c["structure_routes"]["name"]
        ):
            c["status"] = "通式重建结构（已逐图复核）"
    return cards
