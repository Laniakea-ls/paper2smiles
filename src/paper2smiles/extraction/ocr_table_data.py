"""Parse OCR tables while retaining spans, limits, and footnotes."""

import html
import re
from html.parser import HTMLParser


def clean_text(value):
    """Normalize OCR cell whitespace and HTML text entities.

    Parameters
    ----------
    value : object
        Source value to normalize or inspect.

    Returns
    -------
    str
        Cleaned textual cell content.
    """
    text = html.unescape(str(value)).replace("$", "")
    for old, new in [
        (r"\pm", "±"),
        (r"\mu", "μ"),
        (r"\leq", "≤"),
        (r"\geq", "≥"),
        (r"\lt", "<"),
        (r"\gt", ">"),
        (r"\times", "×"),
    ]:
        text = text.replace(old, new)
    for _ in range(4):
        text = re.sub(
            r"\\(?:text|mathrm|mathbf|mathit|operatorname)\s*\{([^{}]*)\}",
            r"\1",
            text,
        )
    text = re.sub(r"_\s*\{([^{}]+)\}", r"\1", text)
    text = re.sub(r"\^\s*\{([^{}]+)\}", r"^\1", text)
    text = re.sub(r"\\(?:left|right|quad|qquad)\b", " ", text)
    text = text.replace("{", "").replace("}", "")
    text = re.sub(r"\\[:,;!]", " ", text)
    return " ".join(text.split()).strip()


class TableParser(HTMLParser):
    """Parse OCR HTML tables while retaining merged-cell provenance."""

    def __init__(self):
        """Initialize the HTML table parser and its span-tracking state."""
        super().__init__(convert_charrefs=True)
        self.tables = []
        self.table = None
        self.row = None
        self.cell = None

    def handle_starttag(self, tag, attrs):
        """Record HTML table boundaries, cells, and embedded image metadata.

        Parameters
        ----------
        tag : str
            HTML element name emitted by the parser.
        attrs : list of tuple
            HTML attribute names and values for the opening tag.
        """
        attrs = dict(attrs)
        if tag == "table":
            self.table = []
        elif tag == "tr" and self.table is not None:
            self.row = []
        elif tag in ["td", "th"] and self.row is not None:
            self.cell = {
                "raw_text": "",
                "images": [],
                "rowspan": int(attrs.get("rowspan", 1)),
                "colspan": int(attrs.get("colspan", 1)),
            }
        elif self.cell is not None:
            if tag == "img":
                self.cell["images"].append(attrs.get("src", ""))
            elif tag == "sup":
                self.cell["raw_text"] += "^"
            elif tag == "br":
                self.cell["raw_text"] += " "

    def handle_endtag(self, tag):
        """Finalize the table row or cell closed by an HTML tag.

        Parameters
        ----------
        tag : str
            HTML element name emitted by the parser.
        """
        if tag in ["td", "th"] and self.cell is not None:
            self.cell["text"] = clean_text(self.cell["raw_text"])
            self.row.append(self.cell)
            self.cell = None
        elif tag == "tr" and self.row is not None:
            self.table.append(self.row)
            self.row = None
        elif tag == "table" and self.table is not None:
            self.tables.append(expand_grid(self.table))
            self.table = None

    def handle_data(self, text):
        """Append source text to the current HTML table cell.

        Parameters
        ----------
        text : str
            Source text to parse or compare.
        """
        if self.cell is not None:
            self.cell["raw_text"] += text


def expand_grid(rows):
    """Expand row and column spans into a rectangular cell grid.

    Parameters
    ----------
    rows : list
        Source rows carrying cell contents or extracted measurements.

    Returns
    -------
    list of list
        Grid cells retaining their originating row and column.
    """
    slots = {}
    width = 0
    for i, row in enumerate(rows):
        j = 0
        for cell in row:
            while (i, j) in slots:
                j += 1
            for y in range(i, i + cell["rowspan"]):
                for x in range(j, j + cell["colspan"]):
                    if (y, x) in slots:
                        raise ValueError("Overlapping HTML table spans")
                    slots[(y, x)] = {**cell, "origin": [i, j]}
                    width = max(width, x + 1)
            j += cell["colspan"]
    return [
        [
            slots.get(
                (i, j),
                {"text": "", "raw_text": "", "images": [], "origin": [i, j]},
            )
            for j in range(width)
        ]
        for i in range(len(rows))
    ]


def parse_tables(markdown):
    """Parse HTML tables embedded in OCR Markdown output.

    Parameters
    ----------
    markdown : str
        OCR Markdown containing HTML table elements.

    Returns
    -------
    list of dict
        Structured tables with cell text, spans, and image references.
    """
    p = TableParser()
    p.feed(markdown)
    return p.tables


def split_marks(text):
    """Separate printed footnote markers from cell text.

    Parameters
    ----------
    text : str
        Source text to parse or compare.

    Returns
    -------
    tuple
        Clean cell text and deduplicated footnote labels.
    """
    text = clean_text(text)
    marks = re.findall(r"\^\s*([a-z]+)\b", text)
    return re.sub(r"\^\s*[a-z]+\b", "", text).strip(), list(
        dict.fromkeys("".join(marks))
    )


def compound_identifier(text, known_marks):
    """Normalize a printed compound identifier and its footnote markers.

    Parameters
    ----------
    text : str
        Source text to parse or compare.
    known_marks : sequence of str
        Footnote labels already identified in the table.

    Returns
    -------
    tuple
        Compound identifier and associated footnotes.
    """
    value, marks = split_marks(text)
    if any(mark not in known_marks for mark in marks):
        raise ValueError("Unknown compound footnote: " + str(marks))
    value = re.sub(r"^\(\s*±\s*\)\s*[-−]\s*", "(±)-", value)
    if not re.fullmatch(r"(?:\(±\)-)?[A-Za-z0-9][A-Za-z0-9._-]*", value):
        raise ValueError("Ambiguous compound identifier: " + value)
    return value, marks


NUMBER = r"(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"


def measurement(text, endpoint, unit):
    """Parse a source activity value without losing limits or uncertainty.

    Parameters
    ----------
    text : str
        Source text to parse or compare.
    endpoint : str
        Reported activity endpoint, such as EC50 or pEC50.
    unit : str
        Unit printed alongside the activity measurement.

    Returns
    -------
    dict
        Measurement fields, original text, units, and parsing status.
    """
    text, marks = split_marks(text)
    if re.fullmatch(r"N\.?\s*D\.?", text, re.I):
        return {
            "endpoint": endpoint,
            "value_original": text,
            "unit": unit,
            "status": "not_determined",
            "footnotes": marks,
            "relation": None,
            "mean": None,
            "sem": None,
        }
    numeric_text = re.sub(
        r"\(\s*(±\s*[\d.eE+\-,]+)\s*\)", r"\1", text
    ).replace(",", "")
    match = re.fullmatch(
        rf"(<=|>=|<|>|≤|≥|=)?\s*({NUMBER})(?:\s*±\s*({NUMBER}))?", numeric_text
    )
    if not match:
        raise ValueError("Unparsed measurement: " + repr(text))
    return {
        "endpoint": endpoint,
        "value_original": text,
        "unit": unit,
        "status": "reported",
        "footnotes": marks,
        "relation": match.group(1) or "=",
        "mean": float(match.group(2)),
        "sem": float(match.group(3)) if match.group(3) else None,
    }


def potency_cell(text, mode):
    """Decode potency cells containing EC50 or paired pEC50 and EC50 values.

    Parameters
    ----------
    text : str
        Source text to parse or compare.
    mode : str
        Encoding used by the potency cell, including paired values.

    Returns
    -------
    list of dict
        Source-preserving potency measurements.
    """
    if mode == "EC50":
        return [measurement(text, "EC50", "nM")]
    if mode != "pEC50_with_EC50":
        raise ValueError("Unsupported potency cell format: " + str(mode))
    clean, marks = split_marks(text)
    match = re.fullmatch(r"(.*?)\s*\(\s*([^()]*)\s*\)", clean)
    if match:
        values = (
            [measurement(match.group(1), "pEC50", "")]
            if match.group(1).strip()
            else []
        ) + [measurement(match.group(2), "EC50", "nM")]
    else:
        values = [measurement(clean, "pEC50", "")]
    for value in values:
        value["footnotes"] = marks
    # EC50 is emitted only if actually printed in parentheses, never
    # backfilled.
    return values
