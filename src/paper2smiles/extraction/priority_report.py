"""One selected structure and its provenance per quantified compound."""

import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from .core import save, sha256
from .structure_selection import canonical, select_structure


class ReportBuilder:
    """Build a source-auditable report from one extraction directory.

    Parameters
    ----------
    root : pathlib.Path or str
        Directory containing extracted compound records and source evidence.
    """

    def __init__(self, root):
        """Bind the extraction directory for this report instance."""
        self.root = Path(root)

    def build(self):
        """Build structures, evidence, JSON, and Excel outputs."""
        self.dest = self.root / "report"
        self.assets = self.dest / "assets"
        self.downloads = self.dest / "downloads"
        for d in [self.dest, self.assets, self.downloads]:
            d.mkdir(parents=True, exist_ok=True)
        self.cards = self._load("compound_records.json", [])
        self.body = self._load("body_activity.json", [])
        self.plan = self._load("plan.json", {})
        self.summary = self._load("summary.json", {})
        self.si = self._load("si_records.json", [])
        self.selections = {
            c["compound_id"]: select_structure(c) for c in self.cards
        }
        manifest = self._load(
            "si_manifest.json", self._load("si_manifest.partial.json", [])
        )
        self.si_files = {}
        for item in manifest:
            p = Path(item["local_path"])
            if p.is_file() and p.suffix.lower() in (".csv", ".tsv", ".xlsx"):
                out = self.downloads / "si" / p.name
                out.parent.mkdir(exist_ok=True)
                shutil.copy2(p, out)
                self.si_files[item["filename"]] = {
                    **item,
                    "href": "downloads/si/" + p.name,
                }
        for name in [
            "summary.json",
            "plan.json",
            "compound_records.json",
            "body_activity.json",
            "names_opsin.json",
            "markush_records.json",
            "si_records.json",
        ]:
            if (self.root / name).exists():
                shutil.copy2(self.root / name, self.downloads / name)
        shutil.copy2(
            self.root / "input/source.pdf", self.downloads / "source.pdf"
        )
        selected_images = []
        for c in self.cards:
            choice = self.selections[c["compound_id"]]
            if choice["route"] in ("markush", "full_structure"):
                selected_images.append(
                    {
                        **c,
                        "markush": choice["evidence"]
                        if choice["route"] == "markush"
                        else [],
                        "full_structure": choice["evidence"]
                        if choice["route"] == "full_structure"
                        else [],
                    }
                )
        self.groups = []
        if selected_images:
            from .table_group_report import prepare

            self.groups, _ = prepare(
                self.root, self.body, selected_images, self.assets
            )
        from .excel_report import write_workbook

        write_workbook(
            self.root,
            self.body,
            self.selections,
            self.si,
            self.downloads / "results.xlsx",
        )
        from .report_data import export_report

        export_report(
            self.root,
            self.cards,
            self.body,
            self.plan,
            self.selections,
            self.si_files,
            self._structure_image,
            self._source,
        )
        self.summary["selected_structure_sources"] = dict(
            Counter((c["route"] for c in self.selections.values()))
        )
        self.summary["report_mode"] = "priority_one_structure"
        save(self.root / "summary.json", self.summary)
        save(self.downloads / "summary.json", self.summary)
        save(
            self.root / "selected_structures.json",
            [
                {
                    "compound_id": cid,
                    "selected_smiles": c["smiles"],
                    "selected_source": c["route"],
                }
                for cid, c in self.selections.items()
            ],
        )

    def _load(self, name, default):
        return (
            json.loads((self.root / name).read_text())
            if (self.root / name).exists()
            else default
        )

    def _structure_image(self, choice):
        mol = Chem.MolFromSmiles(choice["smiles"])
        if choice["route"] == "markush":
            from .markush_chemistry import assemble_record, write_colored_svg

            assembled = assemble_record(choice["evidence"][0])
            if Chem.MolToSmiles(assembled) != canonical(choice["smiles"]):
                raise ValueError("Reconstruction mismatch")
            return write_colored_svg(
                assembled,
                self.assets,
                "selected-"
                + hashlib.sha256(choice["smiles"].encode()).hexdigest()[:16],
                source_pdf_hash=sha256(self.root / "input/source.pdf"),
            )
        from .source_depiction import align_to_source

        aligned = align_to_source(mol, sha256(self.root / "input/source.pdf"))
        drawer = rdMolDraw2D.MolDraw2DSVG(650, 360)
        if aligned:
            drawer.DrawMolecule(
                rdMolDraw2D.PrepareMolForDrawing(
                    mol, addChiralHs=False, wedgeBonds=False
                )
            )
        else:
            rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        name = (
            "selected-"
            + hashlib.sha256(svg.encode()).hexdigest()[:20]
            + ".svg"
        )
        (self.assets / name).write_text(svg)
        return "assets/" + name

    def _source(self, card, choice):
        if choice["route"] == "markush":
            from .table_group_report import compound_source

            return compound_source(
                self.root,
                card,
                choice["evidence"][0],
                self.groups,
                self.body,
                self.assets,
            )
        if choice["route"] == "full_structure":
            from .core import fitz
            from .source_overlay import table_overlay

            r = choice["evidence"][0]
            src = Path(r["crop_path"])
            name = "full-" + sha256(src)[:20] + ".png"
            shutil.copy2(src, self.assets / name)
            g = next(
                (g for g in self.groups if g["table_id"] == r["table_id"]),
                None,
            )
            context = []
            if g:
                with fitz.open(self.root / "input/source.pdf") as doc:
                    pg = doc[r["page"] - 1]
                    size = r["bbox_page_size"]
                    box = [
                        v
                        * (
                            pg.rect.width / size[0]
                            if i % 2 == 0
                            else pg.rect.height / size[1]
                        )
                        for i, v in enumerate(r["bbox_pixels"])
                    ]
                context = [
                    table_overlay(
                        self.root,
                        g,
                        [("实际结构截图", box, "#dc4444")],
                        self.assets,
                    )
                ]
            return {
                "images": context + ["assets/" + name],
                "context_images": context,
            }
        return {"images": [], "context_images": []}
