"""Prepare merged fragment cells for explicit graph-port review.

Recognition is evidence, not an accepted structure. Empty cells never inherit
the preceding row implicitly; port maps require a reviewed graph plan.
"""

import json
import subprocess
import sys
from pathlib import Path

from PIL import Image

from .core import fitz, save
from .ruled_fragment_cells import ruled_cells


def ink_bounds(image, padding=2):
    """Find meaningful fragment ink while rejecting blank and dash-only cells.

    Parameters
    ----------
    image : PIL.Image.Image
        Original crop whose meaningful ink must be located.
    padding : int
        Pixel margin to include around detected fragment ink.

    Returns
    -------
    tuple
        Padded pixel bounds, or a rejection reason when no fragment exists.
    """
    gray = image.convert("L")
    bounds = gray.point(lambda v: 255 if v < 180 else 0).getbbox()
    if bounds is None:
        return None, "blank_cell"
    a, b, c, d = bounds
    if d - b <= 6 and c - a <= 35:
        return None, "dash_or_small_mark"
    return [
        max(0, a - padding),
        max(0, b - padding),
        min(image.width, c + padding),
        min(image.height, d + padding),
    ], None


def prepare(pdf_path, folder, entry, columns, ordered_ids):
    """Crop and recognize table fragments for source-reviewed graph assembly.

    Parameters
    ----------
    pdf_path : pathlib.Path or str
        Path to the original literature PDF.
    folder : pathlib.Path
        Directory for cached inputs, predictions, or generated files.
    entry : dict
        Table layout metadata used to locate source fragment cells.
    columns : dict
        Fragment column indices mapped to their substituent labels.
    ordered_ids : list of str
        Compound identifiers in their original table order.

    Returns
    -------
    dict
        Preparation counts and the remaining visual-review requirements.
    """
    from .full_structure_workflow import SAR

    sys.path.insert(0, str(SAR))
    from pipeline.smiles_validation.ocsr_render import render_inputs

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    with fitz.open(pdf_path) as pdf:
        page = pdf[entry["page"] - 1]
        size = [page.rect.width, page.rect.height]
        rows = ruled_cells(
            page, entry["bbox_pdf_points"], columns, ordered_ids
        )
        if {r["compound_id"] for r in rows} != set(ordered_ids):
            raise ValueError(
                "OCR rows and ruled table rows do not map completely"
            )
        cells = {}
        for row in rows:
            for label, part in row["parts"].items():
                cell = cells.setdefault(
                    part["cell_id"],
                    dict(
                        part,
                        node_label=label,
                        page=entry["page"],
                        compound_ids=[],
                    ),
                )
                cell["compound_ids"].append(row["compound_id"])
        specs = [
            dict(
                crop=k + ".png",
                page=entry["page"],
                bbox_pixels=c["bbox_pdf"],
                bbox_page_size=size,
            )
            for k, c in cells.items()
        ]
        mapping, audit = render_inputs(
            pdf_path, specs, folder, folder / "cells", dpi=300
        )
        active = []
        tight = []
        for cell, spec in zip(cells.values(), specs):
            path = mapping[str((folder / spec["crop"]).resolve())]
            cell["original_cell_crop"] = path
            with Image.open(path) as im:
                bounds, status = ink_bounds(im)
                w, h = im.size
            if status:
                cell["recognition_status"] = status
                continue
            a, b, c, d = bounds
            x, y, u, v = cell["bbox_pdf"]
            box = [
                x + a / w * (u - x),
                y + b / h * (v - y),
                x + c / w * (u - x),
                y + d / h * (v - y),
            ]
            cell["tight_bbox_pdf"] = box
            active.append(cell)
            tight.append(dict(spec, bbox_pixels=box))
        if tight:
            mapping, render_audit = render_inputs(
                pdf_path, tight, folder, folder / "tight", dpi=300
            )
            for cell, spec in zip(active, tight):
                cell["crop_path"] = mapping[
                    str((folder / spec["crop"]).resolve())
                ]
            save(folder / "tight_render_manifest.json", render_audit)
            result = folder / "recognition.json"
            subprocess.run(
                [
                    "/data2/aidd/aidd_platform/tools/MolParser-Mobile/run.sh",
                    *[c["crop_path"] for c in active],
                    "--device",
                    "cpu",
                    "--output",
                    str(result),
                ],
                check=True,
                capture_output=True,
                timeout=900,
            )
            recognized = json.loads(result.read_text())
            if len(recognized) != len(active):
                raise ValueError("Cell recognition count mismatch")
            for cell, prediction in zip(active, recognized):
                cell["recognition"] = prediction
    save(folder / "rows.json", rows)
    save(folder / "cells.json", list(cells.values()))
    save(folder / "cell_render_manifest.json", audit)
    return {
        "status": "fragment_graph_prepared",
        "rows": len(rows),
        "unique_cells": len(cells),
        "recognized_cells": len(active),
        ("review_required"): (
            "Explicit graph ports and fixed core require source review; no "
            "structure adopted automatically."
        ),
    }
