"""Extract source-grounded compound structures and literature activities."""

from importlib import import_module

# Initialize XML/PDF libraries before chemistry extensions on deployed builds.
import_module("paper2smiles._pdf")
