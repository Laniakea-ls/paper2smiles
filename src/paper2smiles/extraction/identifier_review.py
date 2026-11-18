"""Review ambiguous OCR identifiers against the original PDF."""

import json
from collections import defaultdict

from .core import fitz, identifier_present, normalize_id, save


def review_identifiers(flow, rows):
    """Review questionable compound IDs against their printed source evidence.

    Parameters
    ----------
    flow : ActivityWorkflow
        Active extraction context with source pages and model access.
    rows : list
        Source rows carrying cell contents or extracted measurements.

    Returns
    -------
    list of dict
        Activity rows with supported identifier corrections.
    """
    document = "\n".join(p["text"] for p in flow.pages.values())
    suspects = defaultdict(dict)
    for r in rows:
        cid = r["compound_id"]
        if identifier_present(cid, document):
            continue
        suggestions = {
            cid.replace("S", "5"),
            cid.replace("I", "1"),
            cid.replace("L", "1"),
            cid.replace("O", "0"),
        } - {cid}
        suggestions = {
            c for c in suggestions if identifier_present(c, document)
        }
        if suggestions:
            suspects[r["page"]][cid] = sorted(suggestions)
    audit = []
    for page, ids in suspects.items():
        image = flow.root / "identifier_review" / f"page-{page}.png"
        image.parent.mkdir(exist_ok=True)
        if not image.exists():
            with fitz.open(flow.pdf) as doc:
                doc[page - 1].get_pixmap(dpi=250).save(image)
        prompt = (
            "Verify compound IDENTIFIERS ONLY in the attached original "
            "PDF page. OCR may confuse 5/S, 1/I/l or 0/O. The suggestions "
            "merely locate possible errors and are not evidence. Read the "
            "printed row identifier directly; do not infer from sequence, "
            "chemical names or activity values. Return "
            '{"reviews":[{"ocr_identifier":"...","printed_identifier":"...'
            '","supported":true,"evidence":"specific table and printed '
            'label"}]}. Leave supported=false if the image is ambiguous. '
            "Do not output structures or activity numbers.\n"
        ) + json.dumps({"pdf_page": page, "suspects": ids})
        result = flow.llm.ask(
            "identifier_review_v1_" + str(page),
            prompt,
            [image],
            max_tokens=3000,
        )
        for review in result.get("reviews", []):
            old = normalize_id(review.get("ocr_identifier"))
            new = normalize_id(review.get("printed_identifier"))
            if (
                review.get("supported") is not True
                or old not in ids
                or new not in ids[old]
            ):
                continue
            for r in rows:
                if r["page"] == page and r["compound_id"] == old:
                    r["ocr_compound_id"] = old
                    r["compound_id"] = new
                    r["identifier_review"] = review
            audit.append({"page": page, **review})
    save(flow.root / "identifier_review/audit.json", audit)
    return rows
