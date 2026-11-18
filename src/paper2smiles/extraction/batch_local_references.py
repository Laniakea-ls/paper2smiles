"""Recover names from explicitly cited local papers and SI identity."""

import copy
import hashlib
import json
from pathlib import Path

from .core import mol_info, name_fold, save


def enrich(root):
    """Recover names from explicitly cited papers in a local batch.

    Parameters
    ----------
    root : pathlib.Path or str
        Extraction directory containing checkpoints and source evidence.

    Returns
    -------
    list of dict
        Audit records for source-supported name enrichment.
    """
    root = Path(root)
    papers = []
    for p in sorted(root.glob("paper-*")):
        if not (p / "compound_records.json").exists():
            continue
        plan = json.loads((p / "plan.json").read_text())
        doc = json.loads((p / "pdf_text.json").read_text())
        papers.append(
            (p, plan, doc, json.loads((p / "names_opsin.json").read_text()))
        )
    audit = []
    for p, plan, doc, names in papers:
        cards = json.loads((p / "compound_records.json").read_text())
        changed = False
        for c in cards:
            if not c["si_smiles"]:
                continue
            target = mol_info(c["si_smiles"])["canonical"]
            if any(
                n["opsin"].get("canonical") == target
                and n.get("eligible_for_structure_comparison", True)
                and (
                    n.get("name_text_supported")
                    or n.get("visual_name_review", {}).get("supported")
                )
                for n in c["names"]
            ):
                continue
            candidates = []
            for source, sourceplan, sourcedoc, sourcenames in papers:
                if source == p:
                    continue
                title = sourceplan["title"]
                pages = [
                    pg
                    for pg in doc["pages"]
                    if name_fold(title) in name_fold(pg["text"])
                ]
                if not pages:
                    continue
                for n in sourcenames:
                    if n.get("source_kind") == "cited_reference" or not n.get(
                        "source_ranges"
                    ):
                        continue
                    if not n.get(
                        "eligible_for_structure_comparison", True
                    ) or not (
                        n.get("name_text_supported")
                        or n.get("visual_name_review", {}).get("supported")
                    ):
                        continue
                    if n["opsin"].get("canonical") != target:
                        continue
                    # Structure must identify one current
                    # ID, rather than an ambiguous cross-
                    # paper mapping.
                    if (
                        sum(
                            mol_info(x["si_smiles"]).get("canonical") == target
                            for x in cards
                            if x["si_smiles"]
                        )
                        != 1
                    ):
                        continue
                    r = copy.deepcopy(n)
                    r.update(
                        compound_id=c["compound_id"],
                        original_compound_id=n["compound_id"],
                        source_kind="cited_reference",
                        source_reference={
                            "title": title,
                            "doi": sourceplan["doi"],
                            "local_pdf": str(source / "input/source.pdf"),
                            "pdf_sha256": hashlib.sha256(
                                (source / "input/source.pdf").read_bytes()
                            ).hexdigest(),
                            "citation_pages": [x["page"] for x in pages],
                        },
                        mapping_basis=(
                            "Exact isomeric canonical structure identity with "
                            "author SI; reference title explicitly present in "
                            "supplied paper"
                        ),
                        completeness_status="exact",
                        eligible_for_structure_comparison=True,
                    )
                    candidates.append(r)
            if not candidates:
                continue
            r = candidates[0]
            c["names"].append(r)
            names.append(r)
            c["has_structure"] = True
            c["status"] = "SI 与完整名称结构一致"
            changed = True
            audit.append(
                {
                    "paper": p.name,
                    "compound_id": c["compound_id"],
                    "reference_doi": r["source_reference"]["doi"],
                    "original_id": r["original_compound_id"],
                    "name": r["name"],
                }
            )
        if changed:
            save(p / "names_opsin.json", names)
            save(p / "compound_records.json", cards)
            s = json.loads((p / "summary.json").read_text())
            s["with_names"] = sum(bool(c["names"]) for c in cards)
            s["with_structure"] = sum(c["has_structure"] for c in cards)
            save(p / "summary.json", s)
            result = json.loads((p / "batch_result.json").read_text())
            result["summary"] = s
            save(p / "batch_result.json", result)
    if audit:
        save(root / "local_reference_recovery.json", audit)
    return audit
