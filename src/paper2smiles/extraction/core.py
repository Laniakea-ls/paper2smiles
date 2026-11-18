"""Source-preserving helpers for a PDF -> names/OPSIN -> SI workflow."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
import subprocess
import unicodedata
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# Import XML consumers before PyMuPDF (the installed builds otherwise clash).
from paper2smiles._pdf import fitz, openpyxl

DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", re.I)


def save(path: Path, value):
    """Atomically save JSON evidence without exposing a partially written file.

    Parameters
    ----------
    path : pathlib.Path or str
        Source or destination filesystem path.
    value : object
        Source value to normalize or inspect.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, default=str)
    )
    tmp.replace(path)


def sha256(path):
    """Compute the SHA-256 digest of a source file.

    Parameters
    ----------
    path : pathlib.Path or str
        Source or destination filesystem path.

    Returns
    -------
    str
        Hexadecimal content digest.
    """
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def normalize_id(value):
    """Normalize a compound identifier without merging distinct numeric tokens.

    Parameters
    ----------
    value : object
        Source value to normalize or inspect.

    Returns
    -------
    str
        Uppercase identifier with common compound prefixes removed.
    """
    value = unicodedata.normalize("NFKC", str(value or "")).strip()
    return (
        re.sub(r"^(?:compound|cmpd\.?|化合物)\s*", "", value, flags=re.I)
        .strip()
        .upper()
    )


def name_fold(value):
    """Normalize case, spacing, and hyphens for source-name comparison.

    Parameters
    ----------
    value : object
        Source value to normalize or inspect.

    Returns
    -------
    str
        Folded comparison key.
    """
    value = unicodedata.normalize("NFKC", str(value)).lower()
    return re.sub(r"[\s\-‐‑‒–—−\u00ad]", "", value)


def format_name(value):
    """Normalize printed chemical-name typography while retaining its text.

    Parameters
    ----------
    value : object
        Source value to normalize or inspect.

    Returns
    -------
    str
        Whitespace-normalized name with standard hyphens.
    """
    value = unicodedata.normalize("NFKC", str(value))
    value = value.translate(str.maketrans({x: "-" for x in "‐‑‒–—−"}))
    return " ".join(value.replace("\u00ad", "").split())


def extract_pdf(path):
    """Extract PDF text lines, bounding boxes, links, and candidate DOIs.

    Parameters
    ----------
    path : pathlib.Path or str
        Source or destination filesystem path.

    Returns
    -------
    dict
        Numbered pages, external links, and DOI occurrence evidence.
    """
    pages = []
    links = []
    with fitz.open(path) as doc:
        for page in doc:
            lines = []
            # Original text-block order preserves the two
            # separate text columns.
            for block in page.get_text("dict")["blocks"]:
                if block.get("type") != 0:
                    continue
                for line in block["lines"]:
                    text = "".join(
                        span["text"] for span in line["spans"]
                    ).strip()
                    if text:
                        lines.append(
                            {
                                "line": len(lines) + 1,
                                "text": text,
                                "bbox": list(line["bbox"]),
                            }
                        )
            pages.append(
                {
                    "page": page.number + 1,
                    "width": page.rect.width,
                    "height": page.rect.height,
                    "lines": lines,
                    "text": "\n".join(line["text"] for line in lines),
                }
            )
            links += [
                {"page": page.number + 1, "url": link["uri"]}
                for link in page.get_links()
                if link.get("uri")
            ]
    candidates = {}
    for page in pages:
        for match in DOI_RE.finditer(page["text"]):
            doi = match.group().rstrip(".,;:)").lower()
            candidates.setdefault(doi, [])
            if page["page"] not in candidates[doi]:
                candidates[doi].append(page["page"])
    return {
        "pages": pages,
        "links": links,
        "doi_candidates": [
            {"doi": k, "pages": v}
            for k, v in sorted(candidates.items(), key=lambda kv: -len(kv[1]))
        ],
    }


def numbered_page(page):
    """Format a page with explicit page and line references.

    Parameters
    ----------
    page : dict
        Extracted page record with text lines and page number.

    Returns
    -------
    str
        Text suitable for requesting source-grounded model citations.
    """
    return "\n".join(
        f"P{page['page']}:L{line['line']} {line['text']}"
        for line in page["lines"]
    )


def source_span(ranges, pages):
    """Validate cited line ranges and recover their text and bounds.

    Parameters
    ----------
    ranges : list of dict
        Printed page and line intervals supporting extracted evidence.
    pages : dict
        PDF page records indexed by one-based page number.

    Returns
    -------
    tuple
        Combined source text and corresponding PDF locations.

    Raises
    ------
    ValueError
        If the supplied page or line ranges are missing or outside the
        source.
    """
    snippets = []
    locations = []
    if not isinstance(ranges, list) or not ranges:
        raise ValueError("Missing source ranges")
    for r in ranges:
        n = int(r["page"])
        start = int(r["start_line"])
        end = int(r["end_line"])
        if (
            n not in pages
            or not 1 <= start <= end <= len(pages[n]["lines"])
            or end - start > 24
        ):
            raise ValueError("Invalid source line range")
        selected = pages[n]["lines"][start - 1 : end]
        snippets.append("\n".join(x["text"] for x in selected))
        locations.append(
            {
                "page": n,
                "start_line": start,
                "end_line": end,
                "bbox": [
                    min(x["bbox"][0] for x in selected),
                    min(x["bbox"][1] for x in selected),
                    max(x["bbox"][2] for x in selected),
                    max(x["bbox"][3] for x in selected),
                ],
            }
        )
    return "\n".join(snippets), locations


def identifier_present(identifier, text):
    # Match full tokens, never a substring of another compound number.
    """Check for a complete printed compound identifier in source text.

    Parameters
    ----------
    identifier : str
        Printed compound identifier to match as a complete token.
    text : str
        Source text to parse or compare.

    Returns
    -------
    bool
        Whether the identifier occurs without matching another ID's
        substring.
    """
    return bool(
        re.search(
            r"(?<![A-Za-z0-9])"
            + re.escape(str(identifier).strip())
            + r"(?![A-Za-z0-9])",
            text,
            re.I,
        )
    )


def mol_info(smiles):
    """Validate SMILES and derive canonical, connectivity, and parent forms.

    Parameters
    ----------
    smiles : str
        SMILES string to validate or normalize.

    Returns
    -------
    dict
        Validity flag and normalized structures, or empty invalid fields.
    """
    try:
        mol = Chem.MolFromSmiles(str(smiles)) if smiles else None
        if mol is None:
            return {
                "valid": False,
                "canonical": "",
                "connectivity": "",
                "parent": "",
            }
        fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        parent = max(
            fragments,
            key=lambda m: (
                sum(a.GetAtomicNum() == 6 for a in m.GetAtoms()),
                m.GetNumHeavyAtoms(),
            ),
        )
        parent = rdMolStandardize.Uncharger().uncharge(parent)
        return {
            "valid": True,
            "canonical": Chem.MolToSmiles(mol, isomericSmiles=True),
            "connectivity": Chem.MolToSmiles(mol, isomericSmiles=False),
            "parent": Chem.MolToSmiles(parent, isomericSmiles=True),
        }
    except Exception:
        return {
            "valid": False,
            "canonical": "",
            "connectivity": "",
            "parent": "",
        }


def opsin_batch(
    names,
    folder,
    java="/usr/bin/java",
    jar=(
        "/data2/aidd/aidd_platform/tools/opsin/opsin-cli-2.9.0-jar-with-depend"
        "encies.jar"
    ),
):
    """Resolve complete names with OPSIN while preserving input-row alignment.

    Parameters
    ----------
    names : list of str
        Complete chemical names in the required output-row order.
    folder : pathlib.Path
        Directory for cached inputs, predictions, or generated files.
    java : str
        Path to the Java executable running the OPSIN command line.
    jar : str
        Path to the locally installed OPSIN executable JAR.

    Returns
    -------
    list of dict
        Per-name parse status, SMILES, and normalized structure metadata.
    """
    folder.mkdir(parents=True, exist_ok=True)
    if not names:
        return []
    formatted = [format_name(name) for name in names]
    source = folder / "names.txt"
    output = folder / "smiles.txt"
    source.write_text("\n".join(formatted) + "\n")
    result = subprocess.run(
        [java, "-jar", jar, "-osmi", str(source), str(output)],
        capture_output=True,
        text=True,
        timeout=180,
    )
    (folder / "stderr.txt").write_text(result.stderr)
    if result.returncode:
        raise RuntimeError("OPSIN process failed; see stderr.txt")
    values = output.read_text().splitlines()
    if len(values) != len(names):
        raise RuntimeError("OPSIN output row count mismatch")
    return [
        {
            "input_name": name,
            "smiles": smi.strip(),
            "status": "parsed" if smi.strip() else "unparsed",
            **mol_info(smi.strip()),
        }
        for name, smi in zip(formatted, values)
    ]


def read_tables(path):
    """Read supported supplementary spreadsheets and delimited text files.

    Parameters
    ----------
    path : pathlib.Path or str
        Source or destination filesystem path.

    Returns
    -------
    list of dict
        Sheet names and textual row values with decoding metadata.
    """
    path = Path(path)
    if path.suffix.lower() in (".csv", ".tsv"):
        raw = path.read_bytes()
        encoding = (
            "utf-16"
            if raw.startswith((b"\xff\xfe", b"\xfe\xff"))
            else "utf-8-sig"
        )
        try:
            text = raw.decode(encoding)
        except UnicodeDecodeError:
            encoding = "cp1252"
            text = raw.decode(encoding)
        delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
        if path.suffix.lower() == ".csv":
            try:
                delimiter = (
                    csv.Sniffer()
                    .sniff(text[:65536], delimiters=",;\t")
                    .delimiter
                )
            except csv.Error:
                pass
        rows = list(
            csv.reader(io.StringIO(text, newline=""), delimiter=delimiter)
        )
        return [
            {
                "sheet": "",
                "rows": rows,
                "encoding": encoding,
                "delimiter": delimiter,
            }
        ]
    if path.suffix.lower() == ".xlsx":
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        try:
            return [
                {
                    "sheet": sheet.title,
                    "rows": [
                        ["" if x is None else str(x) for x in row]
                        for row in sheet.iter_rows(values_only=True)
                    ],
                }
                for sheet in wb
            ]
        finally:
            wb.close()
    return []


def parse_number(raw):
    """Parse a finite numeric value and its printed measurement relation.

    Parameters
    ----------
    raw : str
        Unmodified text or SMILES returned by the source reader.

    Returns
    -------
    dict or None
        Numeric value and relation, or None when parsing is unsupported.
    """
    text = str(raw).strip().replace(",", "").replace("−", "-")
    text = re.sub(r"\(\s*(±\s*[\d.eE+\-]+)\s*\)", r"\1", text)
    m = re.fullmatch(
        (
            "(<=|>=|<|>|≤|≥|=)?\\s*([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+"
            "-]?\\d+)?)(?:\\s*±\\s*[\\d.eE+-]+)?"
        ),
        text,
    )
    if not m:
        return None
    value = float(m.group(2))
    if not math.isfinite(value):
        return None
    return {
        "relation": {"≤": "<=", "≥": ">="}.get(m.group(1), m.group(1) or "="),
        "value": value,
    }


def potency_nm(value, endpoint, unit):
    """Convert potency to nanomolar units while preserving limits.

    Parameters
    ----------
    value : object
        Source value to normalize or inspect.
    endpoint : str
        Reported activity endpoint, such as EC50 or pEC50.
    unit : str
        Unit printed alongside the activity measurement.

    Returns
    -------
    dict or None
        Converted value, unit, relation, and derivation flag when valid.
    """
    number = parse_number(value)
    if not number:
        return None
    endpoint = re.sub(r"[^a-z0-9]", "", str(endpoint).lower())
    unit = str(unit).replace("μ", "u").replace("µ", "u").lower().strip()
    if endpoint == "pec50":
        if unit not in ("", "dimensionless", "unitless", "pec50"):
            return None
        relation = {"<": ">", ">": "<", "<=": ">=", ">=": "<=", "=": "="}[
            number["relation"]
        ]
        try:
            converted = 10 ** (9 - number["value"])
        except OverflowError:
            return None
    elif endpoint == "ec50":
        multiplier = {
            "nm": 1,
            "um": 1000,
            "mm": 1e6,
            "m": 1e9,
            "pm": 0.001,
        }.get(unit)
        if multiplier is None:
            return None
        relation = number["relation"]
        converted = number["value"] * multiplier
    else:
        return None
    return {
        "value": converted,
        "relation": relation,
        "unit": "nM",
        "derived": endpoint == "pec50",
    }


def assay_family(text):
    """Normalize common assay names for source-aware comparisons.

    Parameters
    ----------
    text : str
        Source text to parse or compare.

    Returns
    -------
    str
        Recognized assay family or a compact fallback comparison key.
    """
    text = str(text).lower()
    if "camp" in text:
        return "camp"
    if "gtp" in text:
        return "gtpgs"
    return re.sub(r"\s+", "", text)


def safe_filename(name):
    """Sanitize an attachment filename for local storage.

    Parameters
    ----------
    name : str
        Checkpoint, output, or asset name within the current directory.

    Returns
    -------
    str
        Bounded filename with unsupported filesystem characters removed.
    """
    name = Path(str(name)).name
    if name in ("", ".", "..") or "\\" in name:
        raise ValueError("Unsafe attachment filename")
    return name
