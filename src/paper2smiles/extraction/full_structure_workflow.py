"""Recognize complete molecular drawings with SAR crops and OCSR."""

import json
import sys
from collections import Counter
from pathlib import Path

from rdkit import Chem

from .core import normalize_id, save, sha256

SAR = Path("/data2/aidd/aidd_platform/pipeline/aidd_sar_extractor")
MODELS = ("MolParser-Mobile-V2", "MolSight-stereo", "OCSRGlyph")
VERSION = "full-structure-v1"


def complete_smiles(raw):
    """Canonicalize a complete molecular drawing prediction.

    Parameters
    ----------
    raw : str
        Unmodified text or SMILES returned by the source reader.

    Returns
    -------
    str
        Isomeric SMILES, or an empty string for invalid or partial
        structures.
    """
    m = Chem.MolFromSmiles(raw or "")
    if (
        m is None
        or not m.GetNumAtoms()
        or any(a.GetAtomicNum() == 0 for a in m.GetAtoms())
    ):
        return ""
    return Chem.MolToSmiles(m, isomericSmiles=True)


def consensus(predictions):
    """Accept complete structures when at least two model predictions agree.

    Parameters
    ----------
    predictions : dict
        Model names mapped to their SMILES predictions and errors.

    Returns
    -------
    dict
        Agreement status, selected SMILES, supporters, and review flags.
    """
    values = {
        k: complete_smiles(v.get("smiles", ""))
        for k, v in predictions.items()
        if not v.get("error")
    }
    count = Counter(v for v in values.values() if v)
    winner, n = count.most_common(1)[0] if count else ("", 0)
    supporters = [k for k, v in values.items() if v == winner and winner]
    valid = len([v for v in values.values() if v])
    status = (
        "three_agree"
        if n == 3
        else (
            "two_agree"
            if n == 2
            else ("needs_review" if valid == 3 else "incomplete")
        )
    )
    return {
        "status": status,
        "accepted": n >= 2,
        "requires_manual_review": status == "needs_review",
        "selected_smiles": winner if n >= 2 else "",
        "candidate_smiles": winner if n >= 2 else "",
        "agreeing_models": supporters if n >= 2 else [],
        "valid_models": valid,
        "basis": "exact_isomeric_smiles_no_salt_or_stereo_stripping",
        "policy": "majority_accept_v1",
        ("note"): (
            "至少两个有效模型的完整立体结构一致即采用；三个有效结果全不同时人"
            "工复核；结果缺失或失败单独重试。"
        ),
    }


def apply_policy(records):
    """Recompute three-model consensus for stored structure predictions.

    Parameters
    ----------
    records : list of dict
        Source evidence records to audit, update, or attach.

    Returns
    -------
    list of dict
        Prediction records with updated acceptance and review status.
    """
    for r in records:
        vote = consensus(r["predictions"])
        r["consensus"] = vote
        r["review_status"] = (
            "auto_accepted"
            if vote["accepted"]
            else (
                "manual_review"
                if vote["requires_manual_review"]
                else "recognition_incomplete"
            )
        )
        r["selected_smiles"] = vote["selected_smiles"]
        r["review_note"] = (
            "按多数一致规则自动采用，保留不同模型的原始结果。"
            if vote["accepted"]
            else (
                "三个有效模型结果各不相同，需人工复核。"
                if vote["requires_manual_review"]
                else "有效识别结果不足且无多数一致，等待补识别或重试。"
            )
        )
    return records


def valid_box(box):
    """Validate full-page bounding boxes normalized to zero through 1000.

    Parameters
    ----------
    box : list of float
        Full-page rectangle normalized independently to 0..1000.

    Returns
    -------
    bool
        Whether the coordinates describe a nonempty in-page rectangle.
    """
    return (
        isinstance(box, list)
        and len(box) == 4
        and all(isinstance(v, (int, float)) for v in box)
        and 0 <= box[0] < box[2] <= 1000
        and 0 <= box[1] < box[3] <= 1000
        and box[2] - box[0] >= 5
        and box[3] - box[1] >= 5
    )


def render_with_edge_check(pdf, rows, folder, renderer):
    # A cut terminal methyl or wedge can still produce a valid but wrong
    # SMILES.
    # Extend only edges with ink touching the crop boundary, then re-
    # render using SAR.
    """Expand clipped molecular crops before sending them to structure readers.

    Parameters
    ----------
    pdf : pathlib.Path or str
        Original literature PDF used as source evidence.
    rows : list
        Source rows carrying cell contents or extracted measurements.
    folder : pathlib.Path
        Directory for cached inputs, predictions, or generated files.
    renderer : callable
        SAR-compatible crop renderer returning paths and an audit.

    Returns
    -------
    tuple
        Rendered crop-path mapping and the SAR rendering audit.

    Raises
    ------
    ValueError
        If molecular ink still touches the crop after bounded expansion.
    """
    from PIL import Image

    from .core import fitz

    with fitz.open(pdf) as doc:
        for row in rows:
            row["initial_bbox_pixels"] = list(row["bbox_pixels"])
            page = doc[row["page"] - 1]
            x, y, u, v = row["bbox_pixels"]
            dx = 4 / page.rect.width * 1000
            dy = 1 / page.rect.height * 1000
            row["bbox_pixels"] = [
                max(0, x - dx),
                max(0, y - dy),
                min(1000, u + dx),
                min(1000, v + dy),
            ]
        for attempt in range(12):
            mapping, audit = renderer(
                pdf, rows, folder, folder / "final", dpi=300
            )
            dirty = []
            for row in rows:
                path = mapping[str((folder / row["crop"]).resolve())]
                with Image.open(path) as im:
                    g = im.convert("L")
                    w, h = g.size
                    bands = [
                        g.crop((0, 0, 1, h)),
                        g.crop((0, 0, w, 1)),
                        g.crop((w - 1, 0, w, h)),
                        g.crop((0, h - 1, w, h)),
                    ]
                    edges = [b.getextrema()[0] < 150 for b in bands]
                if not any(edges):
                    continue
                dirty.append(row["compound_id"])
                page = doc[row["page"] - 1]
                dx = 1 / page.rect.width * 1000
                dy = 1 / page.rect.height * 1000
                x, y, u, v = row["bbox_pixels"]
                row["bbox_pixels"] = [
                    max(0, x - dx) if edges[0] else x,
                    max(0, y - dy) if edges[1] else y,
                    min(1000, u + dx) if edges[2] else u,
                    min(1000, v + dy) if edges[3] else v,
                ]
            if not dirty:
                return mapping, audit
        raise ValueError(
            "Crop ink still touches border after expansion: " + ",".join(dirty)
        )


def run(flow, plan, body):
    """Locate full molecular drawings and apply the three-model readers.

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
        Source-reviewed crops, predictions, and consensus decisions.
    """
    sys.path.insert(0, str(SAR))
    from patent_sar.model_adapters import (
        glyph_iface,
        molparser_iface,
        ocsrglyph_iface,
    )
    from pipeline.smiles_validation.ocsr_render import render_inputs

    root = flow.root / "full_structure"
    root.mkdir(exist_ok=True)
    tables = (
        json.loads((flow.root / "activity_ocr/tables.json").read_text())
        if (flow.root / "activity_ocr/tables.json").exists()
        else []
    )
    crops = []
    audit = []
    for t in tables:
        tid = t["table_id"]
        page = t["page"]
        ids = sorted(
            {r["compound_id"] for r in body if r.get("table_id") == tid}
            & set(flow.target_ids)
        )
        if not ids:
            continue
        folder = root / tid
        folder.mkdir(exist_ok=True)
        event = {"table_id": tid, "page": page}
        try:
            prompt = (
                "Locate full molecular structures in the specified table "
                "on the attached original PDF page. This is geometric "
                "localization ONLY; do not output SMILES, names, "
                "activities, or reconstructed molecules. Return "
                '{"table_type":"full_structures|scaffold_fragments|text_on'
                'ly|uncertain","table_label":"actual '
                'label","rows":[{"compound_id":"29","structure_bbox":[0,0,'
                '0,0],"identifier_evidence":"printed identifier and '
                'location","is_complete_structure":true}],"reason":"eviden'
                'ce"}. Bounding boxes are [left,top,right,bottom] in FULL '
                "PAGE coordinates normalized independently to 0..1000. "
                "Include the ENTIRE single molecular drawing and all atom "
                "labels/wedges with small whitespace; exclude "
                "identifiers, unrelated text, other rows and adjacent "
                "molecules. Only include IDs from the supplied allowlist, "
                "using the printed row mapping. If a table mixes full "
                "molecules and fragments, include ONLY complete rows. "
                "R-group substituent fragments must NEVER be returned as "
                "full structures. Set rows=[] for uncertain or "
                "unsupported mappings.\n"
            ) + json.dumps(
                {
                    "table_id": tid,
                    "pdf_page": page,
                    "caption": t["schema"].get("caption"),
                    "allowed_compound_ids": ids,
                }
            )
            layout = flow.llm.ask(
                VERSION + "_locate_" + tid,
                prompt,
                [flow.image(page)],
                max_tokens=4500,
            )
            save(folder / "layout.json", layout)
            if layout.get("table_type") != "full_structures":
                event.update(
                    status="not_full_structure_table",
                    reason=layout.get("reason"),
                )
                audit.append(event)
                continue
            accepted = []
            seen = set()
            for row in layout.get("rows", []):
                cid = normalize_id(row.get("compound_id"))
                box = row.get("structure_bbox")
                if (
                    cid not in ids
                    or cid in seen
                    or row.get("is_complete_structure") is not True
                    or not valid_box(box)
                ):
                    continue
                seen.add(cid)
                accepted.append(
                    {
                        "compound_id": cid,
                        "crop": f"cell-{len(accepted)}.png",
                        "page": page,
                        "bbox_pixels": box,
                        "bbox_page_size": [1000, 1000],
                        "identifier_evidence": row.get(
                            "identifier_evidence", ""
                        ),
                    }
                )
            mapping, render_audit = render_with_edge_check(
                flow.pdf, accepted, folder, render_inputs
            )
            save(folder / "render_manifest.json", render_audit)
            for row in accepted:
                row["crop_path"] = mapping[
                    str((folder / row["crop"]).resolve())
                ]
            if accepted:
                review_prompt = (
                    "Check these proposed full-molecule crops against the "
                    "original PDF page (first image); subsequent images "
                    "are crops in the supplied order. Do not output "
                    "SMILES. For each crop check 1) correct printed "
                    "compound ID/table row, 2) entire molecule including "
                    "terminals/wedges/labels inside crop, 3) exactly one "
                    "complete molecular drawing, not an R substituent or "
                    'generic scaffold. Return {"checks":[{"compound_id":"2'
                    '9","supported":true,"reason":"specific source '
                    'evidence"}]}. Use supported=false if any point is '
                    "uncertain.\n"
                ) + json.dumps(
                    {
                        "table_id": tid,
                        "rows": [
                            {
                                "compound_id": r["compound_id"],
                                "bbox": r["bbox_pixels"],
                            }
                            for r in accepted
                        ],
                    }
                )
                review = flow.llm.ask(
                    VERSION + "_crop_review_" + tid,
                    review_prompt,
                    [flow.image(page)]
                    + [Path(r["crop_path"]) for r in accepted],
                    max_tokens=3500,
                )
                save(folder / "crop_review.json", review)
                checks = {
                    normalize_id(c.get("compound_id")): c
                    for c in review.get("checks", [])
                }
                for row in accepted:
                    check = checks.get(row["compound_id"], {})
                    if check.get("supported") is not True:
                        event.setdefault("rejected_crops", []).append(
                            {
                                "compound_id": row["compound_id"],
                                "reason": check.get(
                                    "reason", "No verification"
                                ),
                            }
                        )
                        continue
                    crops.append(
                        {
                            **row,
                            "table_id": tid,
                            "crop_review": check,
                            "crop_sha256": sha256(row["crop_path"]),
                            "dpi": 300,
                            "source_pdf_sha256": sha256(flow.pdf),
                        }
                    )
            event["status"] = "localized"
        except Exception as e:
            event.update(status="failed", error=str(e))
            flow.warnings.append(
                "完整结构截图定位 " + tid + " 失败：" + str(e)
            )
        audit.append(event)
        save(root / "audit.json", audit)
    save(root / "audit.json", audit)
    save(root / "crops.json", crops)
    paths = list(dict.fromkeys(r["crop_path"] for r in crops))
    predictions = {}
    # Use the exact SAR adapters and checkpoint-aware caches. Each model
    # is a separate subprocess.
    for name, reader in [
        (MODELS[0], molparser_iface.read_batch),
        (MODELS[1], glyph_iface.read_batch),
        (MODELS[2], ocsrglyph_iface.read_batch),
    ]:
        if not paths:
            predictions[name] = {}
            continue
        flow.status("完整结构截图识别：" + name, crops=len(paths))
        try:
            predictions[name] = reader(paths, timeout=600)
        except Exception as e:
            predictions[name] = {
                p: {"error": str(e), "smiles": ""} for p in paths
            }
        save(root / (name + ".json"), predictions[name])
    records = []
    for crop in crops:
        bymodel = {
            name: {
                **predictions[name].get(
                    crop["crop_path"], {"error": "No prediction", "smiles": ""}
                ),
                "display_model": name,
            }
            for name in MODELS
        }
        for pred in bymodel.values():
            pred["canonical_smiles"] = complete_smiles(pred.get("smiles", ""))
        vote = consensus(bymodel)
        records.append(
            {
                **crop,
                "source_type": "main_pdf_full_structure_three_models",
                "predictions": bymodel,
                "consensus": vote,
                "review_status": "model_consensus_needs_source_review"
                if vote["candidate_smiles"]
                else "model_disagreement_or_failure",
                ("review_note"): (
                    "完整结构截图；三模型逐一保留。候选不自动计作已人工确认结"
                    "构。"
                ),
            }
        )
    apply_policy(records)
    save(flow.root / "full_structure_records.json", records)
    return records


def attach(cards, records):
    """Attach accepted full-structure predictions and retain dissenting models.

    Parameters
    ----------
    cards : list of dict
        Compound cards containing source-specific structure evidence.
    records : list of dict
        Source evidence records to audit, update, or attach.

    Returns
    -------
    list of dict
        Updated compound cards with complete model evidence.
    """
    apply_policy(records)
    for c in cards:
        c["full_structure"] = [
            r for r in records if r["compound_id"] == c["compound_id"]
        ]
        c.setdefault("structure_routes", {})["full_structure_candidate"] = (
            bool(c["full_structure"])
        )
        accepted = sorted(
            [r for r in c["full_structure"] if r["consensus"]["accepted"]],
            key=lambda r: -len(r["consensus"]["agreeing_models"]),
        )
        c["full_structure_smiles"] = (
            accepted[0]["selected_smiles"] if accepted else ""
        )
        c["structure_routes"]["full_structure_accepted"] = bool(accepted)
        if accepted:
            c["has_structure"] = True
            c["status"] = (
                "完整结构识别已采用（"
                + str(len(accepted[0]["consensus"]["agreeing_models"]))
                + "/3 模型一致）"
            )
        elif c["full_structure"] and not c.get("has_structure"):
            c["status"] = (
                "三个模型结果不同，待人工复核"
                if any(
                    r["consensus"]["requires_manual_review"]
                    for r in c["full_structure"]
                )
                else "完整结构识别未完成，待重试"
            )
    return cards


def refresh_saved_policy(root):
    """Reapply the consensus policy to an existing extraction directory.

    Parameters
    ----------
    root : pathlib.Path or str
        Extraction directory containing checkpoints and source evidence.
    """
    root = Path(root)
    path = root / "full_structure_records.json"
    if not path.exists():
        return
    records = apply_policy(json.loads(path.read_text()))
    cards = attach(
        json.loads((root / "compound_records.json").read_text()), records
    )
    save(path, records)
    save(root / "compound_records.json", cards)
    summary = json.loads((root / "summary.json").read_text())
    summary.update(
        with_structure=sum(bool(c["has_structure"]) for c in cards),
        with_full_structure_accepted=sum(
            bool(c["full_structure_smiles"]) for c in cards
        ),
        full_structure_policy="majority_accept_v1",
    )
    save(root / "summary.json", summary)
    result = root / "batch_result.json"
    if result.exists():
        value = json.loads(result.read_text())
        value["summary"] = summary
        save(result, value)
