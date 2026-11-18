"""Interpret table schemas and parse numeric OCR cells."""

import importlib.util
import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .core import fitz, normalize_id, save
from .ocr_table_data import (
    compound_identifier,
    measurement,
    parse_tables,
    potency_cell,
    split_marks,
)

ENDPOINTS = {
    "EC50",
    "pEC50",
    "IC50",
    "pIC50",
    "Ki",
    "pKi",
    "Kd",
    "pKd",
    "Emax",
}


def activity_identifier(raw, known):
    """Separate a printed compound identifier from aliases and footnotes.

    Parameters
    ----------
    raw : str
        Unmodified text or SMILES returned by the source reader.
    known : set of str
        Known printed identifiers used to disambiguate OCR text.

    Returns
    -------
    tuple
        Normalized compound ID, footnotes, and any printed alias.
    """
    value, marks = split_marks(raw)
    reverse = re.fullmatch(
        r"(\d+[A-Za-z]?)\s*\(([A-Za-z][A-Za-z0-9._-]+)\)", value
    )
    if reverse:
        return reverse[1], marks, reverse[2]
    match = re.fullmatch(r"(.+?)\s*\((\d+[A-Za-z]?)\)", value)
    if match:
        cid, checked = compound_identifier(
            match[2] + "".join("^" + x for x in marks), known
        )
        return cid, checked, match[1].strip()
    cid, marks = compound_identifier(raw, known)
    return cid, marks, ""


def cell_values(cell, mapping, schema):
    """Parse an activity cell using its endpoint and encoding metadata.

    Parameters
    ----------
    cell : str
        Original OCR cell text containing the activity value.
    mapping : dict
        Column-to-endpoint mapping selected for the source table.
    schema : dict
        Reviewed interpretation of the source table columns.

    Returns
    -------
    list of dict
        Measurements retaining units, limits, and source uncertainty.
    """
    text, marks = split_marks(cell)
    if re.fullmatch(r"N\.?A\.?", text, re.I):
        notes = [
            n["text"]
            for n in schema.get("footnotes", [])
            if n["marker"] in marks
        ]
        bounds = []
        for note in notes:
            hit = re.search(
                (
                    "(EC50|IC50|Ki|Kd)\\s*(>=|<=|>|<|≥|≤)\\s*([\\d.]+)\\s*(nM|"
                    "[µμu]M|mM)"
                ),
                note,
                re.I,
            )
            if hit:
                endpoint = {
                    "ec50": "EC50",
                    "ic50": "IC50",
                    "ki": "Ki",
                    "kd": "Kd",
                }[hit[1].lower()]
                value = measurement(hit[2] + hit[3], endpoint, hit[4])
                value["value_evidence"] = note
                value["value_origin"] = (
                    "explicit numeric bound in linked footnote"
                )
                bounds.append(value)
        if len(bounds) != 1:
            raise ValueError(
                "NA has no unique explicit numeric footnote bound"
            )
        return bounds
    if mapping.get("encoding") == "EC50_with_Emax":
        pair = re.fullmatch(r"(.*?)\s*\(\s*([^()]+)\s*\)", text)
        if not pair:
            raise ValueError("Expected explicit EC50 (Emax) pair")
        values = [
            measurement(pair[1], "EC50", mapping.get("unit", "")),
            measurement(pair[2], "Emax", "%"),
        ]
        for v in values:
            v["footnotes"] = marks
        return values
    if mapping.get("encoding", "plain") == "plain" and mapping.get(
        "error_type"
    ) in ("SD", "SEM"):
        pair = re.fullmatch(r"([\d.,]+)\s*\(\s*([\d.,]+)\s*\)", text)
        if pair:
            value = measurement(
                pair[1] + " ± " + pair[2],
                mapping["endpoint"],
                mapping.get("unit", ""),
            )
            value["value_original"] = text
            value["value_evidence"] = (
                "Parenthesized error interpreted using explicit "
                + mapping["error_type"]
                + " column schema"
            )
            value["footnotes"] = marks
            return [value]
    return (
        potency_cell(cell, "pEC50_with_EC50")
        if mapping.get("encoding") == "pEC50_with_EC50"
        else [measurement(cell, mapping["endpoint"], mapping.get("unit", ""))]
    )


def extract_tables(flow, plan):
    """Extract and interpret activity tables using the local OCR services.

    Parameters
    ----------
    flow : ActivityWorkflow
        Active extraction context with source pages and model access.
    plan : dict
        Document plan containing article metadata and relevant source pages.

    Returns
    -------
    tuple
        Parsed activity records and per-table audit records.
    """
    root = flow.root / "activity_ocr"
    root.mkdir(exist_ok=True)
    inputs = [
        {"id": f"page-{p:03}", "page": p, "image": str(flow.image(p))}
        for p in plan["main_activity_pages"]
    ]
    save(root / "inputs.json", inputs)
    if not inputs:
        return [], []
    path = Path(
        (
            "/data2/aidd/aidd_platform/pipeline/aidd_sar_extractor/pipeline/ac"
            "tivity_table_extraction/paddle_service_ocr.py"
        )
    )
    spec = importlib.util.spec_from_file_location("batch_paddle", path)
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    adapter.CACHE_ROOT = root / "cache"
    adapter.run(
        root,
        workers=2,
        resume=True,
        urls=["http://127.0.0.1:18116", "http://127.0.0.1:18117"],
    )
    tables = []
    for inp in inputs:
        md = (root / "paddle" / (inp["id"] + ".md")).read_text()
        for index, grid in enumerate(parse_tables(md)):
            if len(grid) < 2:
                continue
            tables.append(
                {
                    **inp,
                    "table_id": inp["id"] + f"-t{index + 1}",
                    "grid": grid,
                    "markdown": md,
                }
            )
    # Second pass: small tables are rendered from the PDF at 300 dpi, using
    # first-pass layout coordinates. Keep the full-page footnotes as
    # schema context.
    zoom_root = root / "table_crops"
    zoom_root.mkdir(exist_ok=True)
    zoom_inputs = []
    with fitz.open(flow.pdf) as pdf:
        for inp in inputs:
            raw = json.loads(
                (root / "paddle" / (inp["id"] + ".json")).read_text()
            )
            blocks = [
                b
                for b in raw.get("parsing_res_list", [])
                if b.get("block_label") == "table"
            ]
            page_tables = [t for t in tables if t["page"] == inp["page"]]
            if len(blocks) != len(page_tables):
                continue
            page = pdf[inp["page"] - 1]
            for t, block in zip(page_tables, blocks):
                x0, y0, x1, y1 = block["block_bbox"]
                if x1 - x0 >= 1000:
                    continue
                sx = page.rect.width / raw["width"]
                sy = page.rect.height / raw["height"]
                box = (
                    fitz.Rect(x0 * sx, y0 * sy, x1 * sx, y1 * sy)
                    + (-4, -4, 4, 4)
                ) & page.rect
                image = zoom_root / (t["table_id"] + ".png")
                if not image.exists():
                    page.get_pixmap(dpi=300, clip=box).save(image)
                zoom_inputs.append(
                    {
                        "id": t["table_id"],
                        "image": str(image),
                        "page": inp["page"],
                        "bbox_pdf_points": list(box),
                        "dpi": 300,
                    }
                )
    if zoom_inputs:
        save(zoom_root / "inputs.json", zoom_inputs)
        adapter.CACHE_ROOT = zoom_root / "cache"
        adapter.run(
            zoom_root,
            workers=2,
            resume=True,
            urls=["http://127.0.0.1:18116", "http://127.0.0.1:18117"],
        )
        byid = {x["id"]: x for x in zoom_inputs}
        for t in tables:
            if t["table_id"] not in byid:
                continue
            md = (zoom_root / "paddle" / (t["table_id"] + ".md")).read_text()
            grids = parse_tables(md)
            if len(grids) == 1:
                t["grid"] = grids[0]
                t["image"] = byid[t["table_id"]]["image"]
                t["markdown"] = (
                    md + "\nFull-page context and footnotes:\n" + t["markdown"]
                )
                t["crop_source"] = byid[t["table_id"]]

    def schema(t):
        prompt = (
            "Identify ONLY numeric bioactivity columns in this source OCR "
            "table. Do not return body-row values. Ignore PK/exposure, animal "
            "doses, solubility, clogP, synthesis yields and statistical P "
            "values. Separate assays, species, and target columns. Return "
            '{"is_activity_table":true,"orientation":"compound_rows","caption"'
            ':"source caption","header_rows":1,"compound_id_column":0,"known_f'
            'ootnotes":["a"],"footnotes":[{"marker":"a","text":"literal '
            'source"}],"measurements":[{"column":2,"endpoint":"EC50","unit":"n'
            'M","assay":"specific assay","target":"GPR88","encoding":"plain","'
            'error_type":"SEM|SD|unspecified","reference":"Emax reference if '
            'applicable"}]}. Use encoding pEC50_with_EC50 ONLY when pEC50 and '
            "explicit EC50 appear in the same cell; its unit must describe "
            "the parenthesized EC50. Set orientation=compound_rows when each "
            "body ROW is a compound. If compounds are across COLUMNS and "
            "properties down rows, set orientation=compound_columns and "
            "return compound_header_row, compound_columns (one unique column "
            "per compound), and measurement_rows instead of measurements; "
            "each measurement_rows item uses row instead of column. Select "
            "only actual potency/efficacy rows, never PK/dose/property rows. "
            "For a cell EC50 (Emax), use encoding EC50_with_Emax with "
            "endpoint EC50 and concentration unit; emit one schema mapping "
            "for the pair, not a second Emax mapping. Distinguish this from "
            "uncertainty in parentheses according to headers/footnotes. "
            "Endpoint must be exactly one of EC50,pEC50,IC50,pIC50,Ki,pKi,Kd,p"
            "Kd,Emax (never combine names). Emax reference is the default "
            "comparator; optionally reference_by_marker maps printed cell "
            "footnotes to comparator overrides. All indices zero-based, "
            "header_rows counts all expanded header rows. Use "
            "is_activity_table=false when unsupported. No invented endpoint, "
            "units, or mappings.\n"
        )
        result = flow.llm.ask(
            "generic_schema_" + t["table_id"],
            prompt
            + json.dumps(
                {
                    "grid": [[c["text"] for c in row] for row in t["grid"]],
                    "ocr": t["markdown"],
                },
                ensure_ascii=False,
            ),
            [Path(t["image"])],
            max_tokens=6000,
        )
        save(root / "schemas" / (t["table_id"] + ".json"), result)
        return t, result

    rows = []
    audit = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(schema, t) for t in tables]):
            t, s = future.result()
            grid = t["grid"]
            coordinate_map = {}
            entry = {
                "table_id": t["table_id"],
                "page": t["page"],
                "schema": s,
                "rows": [],
                "rejected": [],
            }
            if not s.get("is_activity_table"):
                audit.append(entry)
                continue
            if s.get("orientation") == "compound_columns":
                original = grid
                cols = [
                    c.get("column") if isinstance(c, dict) else c
                    for c in s.get("compound_columns", [])
                ]
                mr = s.get("measurement_rows", [])
                hr = s.get("compound_header_row")
                if (
                    not isinstance(hr, int)
                    or not 0 <= hr < len(original)
                    or not cols
                    or not mr
                ):
                    entry["error"] = "Invalid transposed table mapping"
                    audit.append(entry)
                    continue
                if any(
                    not isinstance(c, int) or not 0 <= c < len(original[hr])
                    for c in cols
                ) or any(
                    not isinstance(m.get("row"), int)
                    or not 0 <= m["row"] < len(original)
                    for m in mr
                ):
                    entry["error"] = "Invalid transposed table coordinates"
                    audit.append(entry)
                    continue
                grid = [
                    [{"text": "compound"}]
                    + [{"text": m["endpoint"]} for m in mr]
                ]
                for ri, ci in enumerate(cols, 1):
                    grid.append(
                        [original[hr][ci]]
                        + [original[m["row"]][ci] for m in mr]
                    )
                    for mi, m in enumerate(mr, 1):
                        coordinate_map[(ri, mi)] = (m["row"] + 1, ci + 1)
                s = {
                    **s,
                    "header_rows": 1,
                    "compound_id_column": 0,
                    "measurements": [
                        {**m, "column": i + 1} for i, m in enumerate(mr)
                    ],
                }
            h = s.get("header_rows")
            col = s.get("compound_id_column")
            width = len(grid[0])
            if (
                not isinstance(h, int)
                or not 1 <= h < len(grid)
                or not isinstance(col, int)
                or not 0 <= col < width
            ):
                entry["error"] = "Invalid header/identifier mapping"
                audit.append(entry)
                continue
            for rid, cells in enumerate(grid[h:], h):
                raw = [c["text"] for c in cells]
                try:
                    cid, marks, alias = activity_identifier(
                        raw[col], set(s.get("known_footnotes", []))
                    )
                except ValueError as e:
                    entry["rejected"].append(
                        {"row": rid, "raw": raw, "reason": str(e)}
                    )
                    continue
                for m in s.get("measurements", []):
                    try:
                        if (
                            (
                                m["endpoint"] not in ENDPOINTS
                                and m.get("encoding")
                                not in ["pEC50_with_EC50", "EC50_with_Emax"]
                            )
                            or not isinstance(m["column"], int)
                            or not 0 <= m["column"] < len(raw)
                        ):
                            raise ValueError("Invalid endpoint/column")
                        cell = raw[m["column"]]
                        if not cell or cell in ["—", "–", "-"]:
                            continue
                        values = cell_values(cell, m, s)
                        for v in values:
                            if (
                                v["endpoint"] == "EC50"
                                and m.get("encoding") == "pEC50_with_EC50"
                                and not v.get("value_origin")
                            ):
                                v["unit"] = m.get("unit", "")
                            row = {
                                "compound_id": normalize_id(cid),
                                "raw_identifier": raw[col],
                                "alias": alias,
                                "footnote_marker": marks,
                                "target": m.get("target")
                                or (
                                    plan.get("targets", [None])[0]
                                    if len(plan.get("targets", [])) == 1
                                    else ""
                                ),
                                "assay": m.get("assay") or "",
                                "caption": s.get("caption", t["table_id"]),
                                "page": t["page"],
                                "endpoint": v["endpoint"],
                                "value": v["value_original"],
                                "unit": v["unit"],
                                "source_type": "main_pdf_paddleocr_table",
                                "table_id": t["table_id"],
                                "source_row": coordinate_map.get(
                                    (rid, m["column"]),
                                    (rid + 1, m["column"] + 1),
                                )[0],
                                "source_column": coordinate_map.get(
                                    (rid, m["column"]),
                                    (rid + 1, m["column"] + 1),
                                )[1],
                                "raw_cell": cell,
                                "value_origin": v.get(
                                    "value_origin", "original OCR cell"
                                ),
                                "value_evidence": v.get(
                                    "value_evidence", cell
                                ),
                                "raw_row": raw,
                                "relation": v["relation"],
                                "mean": v["mean"],
                                "sem": v["sem"]
                                if m.get("error_type") == "SEM"
                                else None,
                                "sd": v["sem"]
                                if m.get("error_type") == "SD"
                                else None,
                                "error_value": v["sem"],
                                "status": v["status"],
                                "error_type": m.get(
                                    "error_type", "unspecified"
                                ),
                                "reference": next(
                                    (
                                        m.get("reference_by_marker", {}).get(f)
                                        for f in v.get("footnotes", [])
                                        if f
                                        in m.get("reference_by_marker", {})
                                    ),
                                    m.get("reference", ""),
                                ),
                                "footnotes": s.get("footnotes", []),
                            }
                            rows.append(row)
                            entry["rows"].append(row)
                    except (ValueError, KeyError, IndexError) as e:
                        entry["rejected"].append(
                            {
                                "row": rid,
                                "raw": raw,
                                "measurement": m,
                                "reason": str(e),
                            }
                        )
            audit.append(entry)
    save(root / "tables.json", audit)
    save(root / "activity.json", rows)
    return rows, audit
