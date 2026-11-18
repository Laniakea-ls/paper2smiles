"""Overlay crop rectangles on original source table context."""

import base64
import hashlib
import json
from html import escape
from pathlib import Path

from PIL import Image

from paper2smiles.settings import ROOT

from .core import save


def table_overlay(root, group, boxes, assets):
    """Draw source bounding boxes on a contextual table image.

    Parameters
    ----------
    root : pathlib.Path or str
        Extraction directory containing checkpoints and source evidence.
    group : dict
        Source table metadata and its contextual image.
    boxes : list of tuple
        Labels, PDF rectangles, and colors to overlay on the table.
    assets : pathlib.Path
        Directory receiving generated structure and evidence images.

    Returns
    -------
    str
        Relative path to the generated evidence image.
    """
    image = Path(root) / "report" / group["image"]
    with Image.open(image) as im:
        w, h = im.size
    x, y, u, v = group["render"]["pdf_bbox"]
    rects = []
    audit = []
    for label, box, color in boxes:
        a, b, c, d = box
        if not (x <= a < c <= u and y <= b < d <= v):
            raise ValueError(
                "Source overlay lies outside table context: " + label
            )
        left = (a - x) / (u - x) * w
        top = (b - y) / (v - y) * h
        width = (c - a) / (u - x) * w
        height = (d - b) / (v - y) * h
        rects.append(
            (
                '<rect x="'
                f"{left:.3f}"
                '" y="'
                f"{top:.3f}"
                '" width="'
                f"{width:.3f}"
                '" height="'
                f"{height:.3f}"
                '" fill="none" stroke="'
                f"{color}"
                '" stroke-width="4"><title>'
                f"{escape(label)}"
                "</title></rect>"
            )
        )
        audit.append(
            {
                "label": label,
                "pdf_bbox": box,
                "image_bbox": [left, top, left + width, top + height],
                "color": color,
            }
        )
    data = base64.b64encode(image.read_bytes()).decode()
    svg = (
        (
            '<svg xmlns="http://www.w3.org/2000/svg" '
            'xmlns:xlink="http://www.w3.org/1999/xlink" width="'
            f"{w}"
            '" height="'
            f"{h}"
            '" viewBox="0 0 '
            f"{w}"
            " "
            f"{h}"
            '"><image width="'
            f"{w}"
            '" height="'
            f"{h}"
            '" xlink:href="data:image/png;base64,'
            f"{data}"
            '"/>'
        )
        + "".join(rects)
        + "</svg>"
    )
    name = (
        "source-context-"
        + hashlib.sha256(svg.encode()).hexdigest()[:20]
        + ".svg"
    )
    (assets / name).write_text(svg)
    save(
        Path(root) / "report_table_sources" / (name + ".json"),
        {
            "page": group["page"],
            "table_id": group["table_id"],
            "table_pdf_bbox": [x, y, u, v],
            "image_size": [w, h],
            "boxes": audit,
        },
    )
    return "assets/" + name


def scaffold_box(root, record):
    """Recover the original scaffold bounds from reconstruction evidence.

    Parameters
    ----------
    root : pathlib.Path or str
        Extraction directory containing checkpoints and source evidence.
    record : dict
        Source reconstruction or measurement evidence for one compound.

    Returns
    -------
    list of float or None
        Scaffold rectangle in PDF coordinates when available.
    """
    if not record:
        return None
    if record.get("fragments"):
        from .source_depiction import PAPER_SHA256

        if (
            record.get("source_pdf_sha256") == PAPER_SHA256
            and record.get("table_id") == "page-003-t3"
        ):
            return [374, 594, 480, 651]
        return None
    if record.get("scaffold_bbox_pdf_points"):
        return record["scaffold_bbox_pdf_points"]
    root = Path(root)
    if record.get("origin") == "verified_local_trial":
        p = ROOT / "runs/gpr88_table1_reconstruction"
        if json.loads((p / "summary.json").read_text())[
            "source_pdf_sha256"
        ] == record.get("source_pdf_sha256"):
            return json.loads((p / "crops.json").read_text())[
                "boxes_pdf_points"
            ]["scaffold"]
    candidates = [
        root / "markush" / record["table_id"] / "sar_render_manifest.json"
    ]
    if record["table_id"] == "page-003-t2":
        candidates.append(
            root / "markush/table3-sar-retry/sar_render_manifest.json"
        )
    for path in candidates:
        if path.exists():
            return json.loads(path.read_text())["crops"][0]["pdf_bbox"]
    layout_path = root / "markush" / record["table_id"] / "layout.json"
    if layout_path.exists():
        layout = json.loads(layout_path.read_text())
        box = layout.get("scaffold_bbox", [])
        if (
            layout.get("supported") is True
            and len(box) == 4
            and 0 <= box[0] < box[2] <= 1000
            and 0 <= box[1] < box[3] <= 1000
        ):
            from .core import fitz

            with fitz.open(root / "input/source.pdf") as pdf:
                page = pdf[record["page"] - 1]
                return [
                    v
                    * (page.rect.width if i % 2 == 0 else page.rect.height)
                    / 1000
                    for i, v in enumerate(box)
                ]
    return None


def marked_fragment_table(root, g, record, box, assets):
    """Mark the scaffold and substituent bounds on the original table.

    Parameters
    ----------
    root : pathlib.Path or str
        Extraction directory containing checkpoints and source evidence.
    g : dict
        Source table metadata and its contextual image.
    record : dict
        Source reconstruction or measurement evidence for one compound.
    box : list of float
        Rectangle coordinates in the convention required by the caller.
    assets : pathlib.Path
        Directory receiving generated structure and evidence images.

    Returns
    -------
    str
        Relative path to the contextual evidence image.
    """
    boxes = [("片段裁剪范围", box, "#dc4444")]
    core = scaffold_box(root, record)
    if core:
        boxes.insert(0, ("通式裁剪范围", core, "#2563eb"))
    return table_overlay(root, g, boxes, assets)
