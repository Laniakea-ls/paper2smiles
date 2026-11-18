"""Load PDF and workbook libraries in a compatible initialization order.

Notes
-----
The deployed XML libraries require openpyxl to be loaded before PyMuPDF.
Import PDF access through this module to preserve that initialization order.
"""

from importlib import import_module

openpyxl = import_module("openpyxl")
fitz = import_module("fitz")

__all__ = ["fitz", "openpyxl"]
