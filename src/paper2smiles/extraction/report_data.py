"""Export structured source evidence for the separate report frontend."""

import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

from .core import save, sha256


def display_activities(rows):
    """Remove duplicate display values without discarding raw evidence.

    Parameters
    ----------
    rows : list
        Source rows carrying cell contents or extracted measurements.

    Returns
    -------
    list of dict
        Activity records suitable for the compact compound view.
    """
    result = {}
    for row in rows:
        assay = re.sub(
            r"\s+assay\s*$", "", row.get("assay", ""), flags=re.I
        ).strip()
        key = (
            row["compound_id"],
            row.get("target", "").lower(),
            assay.lower(),
            row.get("endpoint"),
            row.get("value"),
            row.get("unit"),
            row.get("relation"),
            row.get("error_type"),
            row.get("reference"),
        )
        if key not in result:
            result[key] = {**row, "assay": assay, "display_sources": []}
        result[key]["display_sources"].append(
            {
                "page": row.get("page"),
                "caption": row.get("caption"),
                "table_id": row.get("table_id"),
            }
        )
    return list(result.values())


def assembly_source(root, record, images, contexts):
    """Locate source scaffold and fragment images for an assembled structure.

    Parameters
    ----------
    root : pathlib.Path or str
        Extraction directory containing checkpoints and source evidence.
    record : dict
        Source reconstruction or measurement evidence for one compound.
    images : sequence
        Source images used for visual review or report evidence.
    contexts : sequence
        Full source table images surrounding the selected crops.

    Returns
    -------
    dict
        Available scaffold and fragment image paths.
    """
    from .core import fitz
    from .source_overlay import scaffold_box

    box = scaffold_box(root, record)
    output = {
        "scaffold": "",
        "fragments": [i for i in images if i not in contexts],
    }
    if box:
        from .full_structure_workflow import SAR

        sys.path.insert(0, str(SAR))
        from pipeline.smiles_validation.ocsr_render import render_inputs

        key = hashlib.sha256(
            json.dumps([record["page"], box]).encode()
        ).hexdigest()[:16]
        folder = root / "report_table_sources" / ("assembly-core-" + key)
        folder.mkdir(parents=True, exist_ok=True)
        with fitz.open(root / "input/source.pdf") as pdf:
            page = pdf[record["page"] - 1]
            size = [page.rect.width, page.rect.height]
        mapping, audit = render_inputs(
            root / "input/source.pdf",
            [
                {
                    "crop": "scaffold.png",
                    "page": record["page"],
                    "bbox_pixels": box,
                    "bbox_page_size": size,
                }
            ],
            folder,
            folder / "final",
            dpi=300,
        )
        src = Path(next(iter(mapping.values())))
        name = "assembly-core-" + sha256(src)[:20] + ".png"
        shutil.copy2(src, root / "report/assets" / name)
        output["scaffold"] = "assets/" + name
        save(folder / "manifest.json", audit)
    return output


def export_report(
    root,
    cards,
    body,
    plan,
    selections,
    si_files,
    structure_image,
    source_renderer=None,
):
    """Export report JSON, structure depictions, and the titled Excel download.

    Parameters
    ----------
    root : pathlib.Path or str
        Extraction directory containing checkpoints and source evidence.
    cards : list of dict
        Compound cards containing source-specific structure evidence.
    body : list of dict
        Activity measurements extracted from the main article.
    plan : dict
        Document plan containing article metadata and relevant source pages.
    selections : dict
        Selected structure evidence indexed by compound identifier.
    si_files : dict
        Supplementary-file metadata and relative download references.
    structure_image : callable
        Render a selected structure and return its relative image path.
    source_renderer : callable
        Return structured source-image evidence for a compound.
    """
    root = Path(root)
    dest = root / "report"
    compounds = []
    for card in cards:
        cid = card["compound_id"]
        choice = selections[cid]
        evidence = choice["evidence"]
        source = {
            "filename": "",
            "href": "",
            "rows": [],
            "names": [],
            "images": [],
            "pages": [],
        }
        if choice["route"] == "si":
            file = si_files.get(card.get("si_source", ""))
            if file:
                source.update(filename=file["filename"], href=file["href"])
            seen = set()
            for r in evidence:
                key = (r.get("source_row"), r.get("si_smiles"))
                if key[0] and key[1] and key not in seen:
                    seen.add(key)
                    source["rows"].append({"number": key[0], "smiles": key[1]})
        elif choice["route"] != "missing":
            source.update(filename="正文 PDF", href="downloads/source.pdf")
            if choice["route"] == "name":
                source["names"] = list(
                    dict.fromkeys(n["name"] for n in evidence)
                )
                source["pages"] = sorted(
                    {
                        r["page"]
                        for n in evidence
                        for r in n.get("source_ranges", [])
                    }
                )
            else:
                source["pages"] = sorted({r["page"] for r in evidence})
                crop_records = []
                for r in evidence:
                    crop_records.extend(r.get("graph_evidence", []))
                    if r.get("crop_path"):
                        crop_records.append(r)
                for r in crop_records:
                    path = Path(r["crop_path"])
                    if not path.is_file():
                        continue
                    name = "evidence-" + sha256(path)[:20] + path.suffix
                    shutil.copy2(path, dest / "assets" / name)
                    src = "assets/" + name
                    if src not in source["images"]:
                        source["images"].append(src)
                if source_renderer:
                    rendered = source_renderer(card, choice)
                    # Reuse SAR-rendered original-table
                    # contexts with their actual crop
                    # overlays.
                    if rendered.get("images"):
                        source["images"] = rendered["images"]
                    if rendered.get("context_images"):
                        source["context_images"] = rendered["context_images"]
                if choice["route"] == "markush":
                    source["assembly"] = assembly_source(
                        root,
                        evidence[0],
                        source["images"],
                        source.get("context_images", []),
                    )
        activity = [r for r in body if r["compound_id"] == cid]
        models = []
        consensus = None
        if choice["route"] == "full_structure":
            from rdkit import Chem

            record = evidence[0]
            consensus = record.get("consensus", {})
            for name, prediction in record.get("predictions", {}).items():
                smiles = (
                    prediction.get("canonical_smiles")
                    or prediction.get("smiles")
                    or ""
                )
                valid = bool(smiles and Chem.MolFromSmiles(smiles))
                models.append(
                    {
                        "name": name,
                        "smiles": smiles,
                        "image": structure_image(
                            {"route": "model", "smiles": smiles}
                        )
                        if valid
                        else "",
                        "agrees": name in consensus.get("agreeing_models", []),
                        "error": prediction.get("error", ""),
                    }
                )
        compounds.append(
            {
                "id": cid,
                "route": choice["route"],
                "label": choice["label"],
                "smiles": choice["smiles"],
                "image": structure_image(choice) if choice["smiles"] else "",
                "source": source,
                "activities": activity,
                "display_activities": display_activities(activity),
                "model_predictions": models,
                "model_consensus": consensus,
            }
        )
    title = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", plan["title"]).strip(" .")
    while len(title.encode("utf-8")) > 220:
        title = title[:-1]
    filename = title + ".xlsx"
    shutil.copy2(
        dest / "downloads/results.xlsx", dest / "downloads" / filename
    )
    save(
        dest / "report-data.json",
        {
            "schema_version": 1,
            "title": plan["title"],
            "doi": plan.get("doi", ""),
            "excel_download": "downloads/" + filename,
            "excel_filename": filename,
            "compounds": compounds,
            "activities": body,
        },
    )
