"""Select one confirmed structure in the requested source priority order."""

from rdkit import Chem


def canonical(smiles):
    """Canonicalize a complete structure without unresolved attachment atoms.

    Parameters
    ----------
    smiles : str
        SMILES string to validate or normalize.

    Returns
    -------
    str or None
        Isomeric canonical SMILES, or None for an incomplete structure.
    """
    mol = Chem.MolFromSmiles(smiles or "")
    return (
        Chem.MolToSmiles(mol, isomericSmiles=True)
        if mol
        and mol.GetNumAtoms()
        and not any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())
        else None
    )


def select_structure(card):
    """Select one confirmed structure in the configured evidence priority.

    Parameters
    ----------
    card : dict
        Compound card with candidate structures and their evidence.

    Returns
    -------
    dict
        Selected route, label, SMILES, and the supporting source evidence.
    """
    if canonical(card.get("si_smiles")):
        return {
            "route": "si",
            "label": "作者 SI",
            "smiles": card["si_smiles"],
            "evidence": card.get("si_evidence", []),
        }
    names = [
        n
        for n in card.get("names", [])
        if n.get("eligible_for_structure_comparison", True)
        and (
            n.get("name_text_supported")
            or n.get("visual_name_review", {}).get("supported")
        )
        and n.get("opsin", {}).get("valid")
        and canonical(n["opsin"].get("smiles"))
    ]
    if len({canonical(n["opsin"]["smiles"]) for n in names}) == 1:
        return {
            "route": "name",
            "label": "正文名称",
            "smiles": names[0]["opsin"]["smiles"],
            "evidence": names,
        }
    graphs = [
        r
        for r in card.get("markush", [])
        if r.get("review_status") == "source_visually_reviewed"
        and canonical(r.get("assembled_smiles"))
    ]
    if len({canonical(r["assembled_smiles"]) for r in graphs}) == 1:
        return {
            "route": "markush",
            "label": "通式连接",
            "smiles": graphs[0]["assembled_smiles"],
            "evidence": graphs,
        }
    full = [
        r
        for r in card.get("full_structure", [])
        if r.get("consensus", {}).get("status") in ("two_agree", "three_agree")
        and canonical(r["consensus"].get("candidate_smiles"))
    ]
    if len({canonical(r["consensus"]["candidate_smiles"]) for r in full}) == 1:
        return {
            "route": "full_structure",
            "label": "原文结构图",
            "smiles": full[0]["consensus"]["candidate_smiles"],
            "evidence": full,
        }
    return {
        "route": "missing",
        "label": "待补充",
        "smiles": "",
        "evidence": [],
    }
