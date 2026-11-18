"""Build compound evidence with original molecule and activity crops."""

import argparse
import json
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from PIL import Image

from .activity_targets import select_targets
from .core import fitz, normalize_id, save
from .ocr_table_data import clean_text


def read(root, name):
    """Read a JSON evidence file from an extraction directory.

    Parameters
    ----------
    root : pathlib.Path or str
        Extraction directory containing checkpoints and source evidence.
    name : str
        Checkpoint, output, or asset name within the current directory.

    Returns
    -------
    object
        Decoded JSON evidence.
    """
    return json.loads((root / name).read_text())


def natural(value):
    """Build a sorting key that orders numeric identifier parts numerically.

    Parameters
    ----------
    value : object
        Source value to normalize or inspect.

    Returns
    -------
    list of tuple
        Comparable numeric and text components of the identifier.
    """
    return [
        (1, int(x)) if x.isdigit() else (0, x.lower())
        for x in re.split(r"(\d+)", str(value))
    ]


def supplemental_figure_activity(output, known_ids):
    """Recover known compound activities from supplemental figure evidence.

    Parameters
    ----------
    output : pathlib.Path
        Directory receiving the generated extraction evidence.
    known_ids : set of str
        Compound identifiers already supported by source evidence.

    Returns
    -------
    tuple
        Activity rows and their associated source crops.

    Notes
    -----
    Add numeric identifiers absent from tables, from a reviewed figure OCR
    crop.
    """
    source = output / "figure_ocr/paddle/figure-2.json"
    if not source.exists():
        return [], {}
    raw = read(source.parent, source.name)
    rows = []
    crops = {}
    image_path = output / "assets/figure-2.png"
    with Image.open(image_path) as image:
        for block in raw["parsing_res_list"]:
            cleaned = clean_text(clean_text(block["block_content"])).replace(
                "~", " "
            )
            match = re.search(
                r"(?<!\d)(\d+),\s*R=.*?cAMP\s*EC50\s*=\s*([\d.]+)\s*nM",
                cleaned,
            )
            if not match or match[1] in known_ids:
                continue
            cid, value = match.groups()
            record_id = "F2:" + cid
            row = {
                "record_id": record_id,
                "source_label": "Figure 2",
                "table_number": None,
                "page": 2,
                "source_row": None,
                "compound_id": cid,
                "raw_identifier": cid,
                "compound_markers": [],
                "compound_context": [
                    "正文 Figure 2 中的明确活性数值，补充表格未收录的编号。"
                ],
                "target": "GPR88",
                "assay": "cAMP",
                "measurements": [
                    {
                        "endpoint": "EC50",
                        "value_original": value,
                        "unit": "nM",
                        "status": "reported",
                        "footnotes": [],
                        "relation": "=",
                        "mean": float(value),
                        "sem": None,
                        "source_column": None,
                    }
                ],
                "raw_cells": [block["block_content"]],
                "source_type": "main_pdf_paddleocr_figure",
                ("value_origin"): (
                    "unmodified OCR figure annotation; deterministic parsing"
                ),
                "error_type": "unspecified",
            }
            rows.append(row)
            x0, y0, x1, y1 = block["block_bbox"]
            box = [
                max(0, x0 - 4),
                max(0, y0 - 2),
                min(image.width, x1 + 4),
                y1,
            ]
            path = output / "assets" / f"F2-{cid}.png"
            image.crop(box).save(path)
            crops[record_id] = {
                "path": str(path.resolve()),
                "header_path": None,
                "source_label": "Figure 2",
                "table_number": None,
                "page": 2,
                "compound_id": cid,
                "source_row": None,
                "bbox_figure_pixels": box,
                "source_figure_image": str(image_path.resolve()),
                ("method"): (
                    "PaddleOCR annotation block; original pixels; source "
                    "image visually reviewed"
                ),
            }
    save(output / "supplemental_activity.json", rows)
    return rows, crops


def activity_crops(activity_root, output):
    """Prepare contextual activity-table crops for compound evidence.

    Parameters
    ----------
    activity_root : pathlib.Path
        Directory containing parsed activity evidence.
    output : pathlib.Path
        Directory receiving the generated extraction evidence.

    Returns
    -------
    dict
        Activity crop references indexed by their source identifiers.
    """
    assets = output / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    result = {}
    tables = read(activity_root, "tables.json")
    pdf_path = Path(read(activity_root, "input.json")["pdf"])
    with fitz.open(pdf_path) as pdf:
        for table in tables:
            number = table["table_number"]
            page = table["page"]
            rows = table["rows"]
            raw = read(activity_root, f"paddle/table-{number}.json")
            block = next(
                b
                for b in raw["parsing_res_list"]
                if b["block_label"] == "table"
            )
            x0, y0, x1, y1 = block["block_bbox"]
            with Image.open(table["image"]) as original:
                image = np.asarray(original.convert("L"))
                if table["schema"]["header_rows"] == 1:
                    # Main-table rows are raster content.
                    # Locate their printed baselines
                    # in the identifier column, using
                    # horizontal rules to exclude the
                    # scaffold.
                    density = (image[:, x0 + 2 : x1 - 2] < 150).mean(axis=1)
                    rules = []
                    for y in np.where(density > 0.65)[0]:
                        if not rules or y - rules[-1][-1] > 3:
                            rules.append([int(y)])
                        else:
                            rules[-1].append(int(y))
                    rules = [g for g in rules if g[0] > 60]
                    if len(rules) != 3:
                        raise ValueError(
                            (
                                "Table "
                                f"{number}"
                                ": expected header/header-end/bottom rules"
                            )
                        )
                    header_start, header_end, bottom = (
                        rules[0][0],
                        rules[1][-1],
                        rules[2][0],
                    )
                    left = x0 + 10
                    right = round(x0 + (x1 - x0) * 0.14)
                    ink = (image[:, left:right] < 130).sum(axis=1)
                    on = (ink >= 2) & (ink < (right - left) * 0.8)
                    groups = []
                    for y in np.where(on)[0]:
                        if not header_end + 3 < y < bottom - 3:
                            continue
                        if not groups or y - groups[-1][-1] > 6:
                            groups.append([int(y)])
                        else:
                            groups[-1].append(int(y))
                    groups = [
                        g for g in groups if len(g) >= 5 and sum(ink[g]) >= 40
                    ]
                    if len(groups) != len(rows):
                        raise ValueError(
                            (
                                "Table "
                                f"{number}"
                                ": OCR rows and printed ID baselines differ: "
                                f"{len(rows)}"
                                " / "
                                f"{len(groups)}"
                            )
                        )
                    centers = [(g[0] + g[-1]) / 2 for g in groups]
                    boundaries = (
                        [header_end + 3]
                        + [(a + b) / 2 for a, b in zip(centers, centers[1:])]
                        + [bottom + 2]
                    )
                    boxes = [
                        [
                            max(0, x0 - 12),
                            int(boundaries[i]),
                            min(original.width, x1 + 12),
                            int(boundaries[i + 1]),
                        ]
                        for i in range(len(rows))
                    ]
                    header_box = [
                        max(0, x0 - 12),
                        header_start - 4,
                        min(original.width, x1 + 12),
                        header_end + 3,
                    ]
                    method = (
                        "raster identifier baselines; row count/order "
                        "validated against original OCR table"
                    )
                else:
                    # Table 5 is native PDF text; use
                    # only identifier locations to crop.
                    words = pdf[page - 1].get_text("words")
                    anchors = []
                    for row in rows:
                        marks = "".join(row["compound_markers"])
                        cid = row["compound_id"]
                        matches = [
                            w
                            for w in words
                            if 58 <= w[0] < 115
                            and 132 < w[1] < 180
                            and (
                                w[4] == cid
                                or (
                                    marks
                                    and re.fullmatch(
                                        re.escape(cid)
                                        + "["
                                        + re.escape(marks)
                                        + "]",
                                        w[4],
                                    )
                                )
                            )
                        ]
                        if len(matches) != 1:
                            raise ValueError(
                                "Ambiguous native row identifier " + cid
                            )
                        anchors.append(matches[0])
                    centers = [(w[1] + w[3]) / 2 for w in anchors]
                    boundaries = (
                        [min(anchors[0][1] - 3, centers[0] - 5)]
                        + [(a + b) / 2 for a, b in zip(centers, centers[1:])]
                        + [anchors[-1][3] + 2]
                    )
                    scale = table["dpi"] / 72
                    clip = table["bbox_pdf_points"]
                    # The remaining Table 5 columns are
                    # PK/properties. Show the ID and
                    # EC50 columns at readable size, with
                    # the complete table still linked.
                    ec_header = next(
                        w for w in words if w[4] == "EC50" and 100 < w[1] < 130
                    )
                    next_header = next(
                        w for w in words if w[4] == "clog" and 100 < w[1] < 130
                    )
                    right = int(
                        ((ec_header[2] + next_header[0]) / 2 - clip[0]) * scale
                    )
                    boxes = [
                        [
                            0,
                            int((boundaries[i] - clip[1]) * scale),
                            right,
                            int((boundaries[i + 1] - clip[1]) * scale),
                        ]
                        for i in range(len(rows))
                    ]
                    header_box = [
                        0,
                        int((ec_header[1] - 3 - clip[1]) * scale),
                        right,
                        boxes[0][1],
                    ]
                    method = (
                        "native PDF identifier coordinates; identifier/EC50 "
                        "columns only; numeric values remain original "
                        "PaddleOCR output"
                    )
                header_path = assets / f"table-{number}-header.png"
                original.crop(header_box).save(header_path)
                for row, box in zip(rows, boxes):
                    name = (
                        "T"
                        f"{number}"
                        "-R"
                        f"{row['source_row']}"
                        "-"
                        f"{row['compound_id']}"
                        ".png"
                    )
                    path = assets / name
                    original.crop(box).save(path)
                    result[row["record_id"]] = {
                        "path": str(path.resolve()),
                        "header_path": str(header_path.resolve()),
                        "table_number": number,
                        "page": page,
                        "compound_id": row["compound_id"],
                        "source_row": row["source_row"],
                        "bbox_table_pixels": box,
                        "source_table_image": table["image"],
                        "method": method,
                    }
    return result


def prepare(activity_root, structure_root, output):
    """Assemble compound evidence cards from activity and structure outputs.

    Parameters
    ----------
    activity_root : pathlib.Path
        Directory containing parsed activity evidence.
    structure_root : pathlib.Path
        Directory containing structure recognition evidence.
    output : pathlib.Path
        Directory receiving the generated extraction evidence.

    Returns
    -------
    dict
        Counts and summary of the prepared compound records.
    """
    activity_root = Path(activity_root).resolve()
    structure_root = Path(structure_root).resolve()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    assets = output / "assets"
    assets.mkdir(exist_ok=True)
    table_rows = read(activity_root, "activity_rows.json")
    known, _ = select_targets(table_rows)
    supplemental, figure_crops = supplemental_figure_activity(
        output, set(known)
    )
    targets, excluded = select_targets(table_rows + supplemental)
    save(
        output / "activity_targets.json",
        {"targets": list(targets.values()), "excluded": excluded},
    )
    source = read(activity_root, "input.json")
    metadata = read(structure_root, "metadata.json")
    if source["pdf_sha256"] != metadata["pdf_sha256"]:
        raise ValueError("Activity and structure sources are different PDFs")
    activity_images = activity_crops(activity_root, output)
    activity_images.update(figure_crops)
    save(output / "activity_crops.json", activity_images)
    candidates = {r["page"]: r for r in read(output, "si_candidates.json")}
    fallback = {r["page"]: r for r in read(output, "si_label_review.json")}
    by_id = defaultdict(list)
    for row in read(output, "si_ocr_labels.json"):
        review = fallback.get(row["page"])
        ids = (
            [review["compound_id"]]
            if review and review.get("compound_id")
            else row["ids"]
        )
        if len(ids) != 1:
            continue
        cid = normalize_id(ids[0])
        if cid not in targets:
            continue
        candidate = candidates[row["page"]]
        name = f"compound-{cid}-si-{row['page']}.png"
        shutil.copy2(candidate["image"], assets / name)
        by_id[cid].append(
            {
                "path": str((assets / name).resolve()),
                "source": "author_si_pdf",
                "page": row["page"],
                "bbox_pdf_points": candidate["candidate_bbox"],
                "pdf_rotation": candidate["rotation"],
                "shared": False,
                "printed_identifier": review["printed_label"]
                if review
                else "Compound " + cid,
                "form": review["form"] if review else "unspecified",
                "identifier_method": "Luna source-image review"
                if review
                else "local PaddleOCR label",
                "ocr_label_file": str(
                    (
                        output / "si_ocr/paddle" / f"si-{row['page']:03}.md"
                    ).resolve()
                ),
            }
        )
    names = read(structure_root, "names_opsin.json")
    names_by_id = defaultdict(list)
    for name in names:
        if normalize_id(name["compound_id"]) in targets:
            names_by_id[normalize_id(name["compound_id"])].append(name)
    si_by_id = {
        normalize_id(r["compound_id"]): r
        for r in read(structure_root, "merged_records.json")
        if normalize_id(r["compound_id"]) in targets
    }
    with fitz.open(source["pdf"]) as pdf:
        # Controls/early compounds share explicitly labelled figures.
        figures = [
            (["2-PCCA", "RTI-13951-33"], [59, 65, 300, 186], "Figure 1"),
            (["1", "2", "3", "4"], [59, 309, 301, 424], "Figure 2"),
        ]
        for ids, box, caption in figures:
            path = assets / (caption.replace(" ", "-").lower() + ".png")
            pdf[1].get_pixmap(dpi=250, clip=fitz.Rect(box)).save(path)
            for cid in ids:
                if cid not in targets:
                    continue
                by_id[cid].append(
                    {
                        "path": str(path.resolve()),
                        "source": "main_pdf",
                        "page": 2,
                        "bbox_pdf_points": box,
                        "shared": True,
                        "shared_compounds": ids,
                        "caption": caption,
                        "printed_identifier": cid,
                        ("identifier_method"): (
                            "visually checked original figure labels"
                        ),
                        ("note"): (
                            "共用骨架与 R / X 定义，按原文编号对应；活性值按活"
                            "性表分别保留。"
                        ),
                    }
                )
        for cid, items in names_by_id.items():
            for index, name in enumerate(items):
                crops = []
                for j, span in enumerate(name["source_ranges"]):
                    page = pdf[span["page"] - 1]
                    box = (
                        fitz.Rect(span["bbox"]) + (-5, -5, 5, 5)
                    ) & page.rect
                    path = assets / f"name-{cid}-{index}-{j}.png"
                    page.get_pixmap(dpi=200, clip=box).save(path)
                    crops.append(
                        {
                            "path": str(path.resolve()),
                            "page": span["page"],
                            "bbox_pdf_points": list(box),
                        }
                    )
                name["evidence_crops"] = crops
    cards = []
    for cid in sorted(targets, key=natural):
        target = targets[cid]
        si = si_by_id.get(cid)
        name_items = names_by_id.get(cid, [])
        if si and any(
            r["source_type"] == "main_pdf_paddleocr_figure"
            for r in target["activity_records"]
        ):
            # Preserve the older review as history, but use the
            # newly acquired figure
            # evidence when comparing an explicitly printed EC50
            # with the SI value.
            from .core import parse_number

            printed = [
                v
                for r in target["activity_records"]
                for v in r["measurements"]
                if v["endpoint"] == "EC50"
            ]
            si_number = parse_number(si["value"])
            if si_number and printed:
                si["prior_activity_comparison"] = si.get("activity_comparison")
                si["activity_comparison"] = (
                    "supported"
                    if any(
                        v["mean"] == si_number["value"]
                        and v["relation"] == si_number["relation"]
                        for v in printed
                    )
                    else "conflict"
                )
                si["activity_comparison_source"] = (
                    "current main-PDF figure OCR"
                )
        valid = [
            n
            for n in name_items
            if n["opsin"]["valid"]
            and (
                n["name_text_supported"]
                or n.get("visual_name_review", {}).get("supported")
            )
        ]
        if si:
            chosen = si["selected_smiles"]
            chosen_source = si["selected_structure_source"]
            comparison = si["structure_comparison"]
        elif len({n["opsin"]["canonical"] for n in valid}) == 1:
            chosen = valid[0]["opsin"]["smiles"]
            chosen_source = "main_name_opsin"
            comparison = "name_only"
        else:
            chosen = ""
            chosen_source = "unresolved"
            comparison = "structure_unresolved"
        molecule_images = sorted(
            by_id.get(cid, []),
            key=lambda x: (x.get("form") != "HCl", x["page"]),
        )
        cards.append(
            {
                **target,
                "names": name_items,
                "structure_record": si,
                "selected_smiles": chosen,
                "selected_structure_source": chosen_source,
                "structure_comparison": comparison,
                "molecule_crops": molecule_images,
                "activity_crops": [
                    activity_images[r["record_id"]]
                    for r in target["activity_records"]
                ],
            }
        )
    save(output / "compounds.json", cards)
    manifest = read(structure_root, "si_manifest.json")
    summary = {
        "activity_compounds": len(cards),
        "exact_numeric_compounds": sum(not c["only_limits"] for c in cards),
        "limit_only_compounds": sum(c["only_limits"] for c in cards),
        "activity_rows": sum(len(c["activity_records"]) for c in cards),
        "table_activity_rows": len(table_rows),
        "supplemental_figure_rows": len(supplemental),
        "molecule_screenshots": sum(bool(c["molecule_crops"]) for c in cards),
        "activity_row_screenshots": len(activity_images),
        "individual_si_structure_compounds": sum(
            any(not x["shared"] for x in c["molecule_crops"]) for c in cards
        ),
        "shared_main_figure_compounds": sum(
            any(x["shared"] for x in c["molecule_crops"]) for c in cards
        ),
        "valid_selected_structures": sum(
            bool(c["selected_smiles"]) for c in cards
        ),
        "structure_source_agreement": sum(
            c["structure_comparison"] in ["exact", "parent_only"]
            for c in cards
        ),
        "activity_tables": len(read(activity_root, "tables.json")),
        "selected_name_records": sum(len(c["names"]) for c in cards),
        "excluded_unselected_name_records": sum(
            normalize_id(n["compound_id"]) not in targets for n in names
        ),
        "structure_comparison": dict(
            Counter(c["structure_comparison"] for c in cards)
        ),
        ("selection_policy"): (
            "Select numeric bioactivity records first, including qualified "
            "bounds. Missing-only records and non-bioactivity properties do "
            "not create targets."
        ),
        "activity_root": str(activity_root),
        "structure_root": str(structure_root),
        "input_pdf": source["pdf"],
        "pdf_sha256": source["pdf_sha256"],
        "si_pdf": next(
            f["local_path"] for f in manifest if f["filename"].endswith(".pdf")
        ),
    }
    save(output / "summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False), flush=True)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activity-run", required=True)
    parser.add_argument("--structure-run", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    prepare(args.activity_run, args.structure_run, args.output)
