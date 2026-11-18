import tempfile
import unittest
from pathlib import Path

from rdkit import Chem

from paper2smiles.extraction.markush_chemistry import (
    assemble,
    write_colored_svg,
)
from paper2smiles.extraction.markush_workflow import (
    attach,
    render_markush_crops,
    reviewed_examples,
)


class MarkushTests(unittest.TestCase):
    def test_chirality_and_origin_survive_join(self):
        m = assemble(
            "*c1ccc([C@H](CO)NC(=O)[C@@H](C)c2ccccc2)cc1", "*OC[C@@H](C)CC"
        )
        self.assertEqual(
            sorted(v for _, v in Chem.FindMolChiralCenters(m)), ["R", "S", "S"]
        )
        self.assertTrue(all(a.HasProp("source_part") for a in m.GetAtoms()))
        self.assertEqual(
            {a.GetProp("source_part") for a in m.GetAtoms()},
            {"scaffold", "r1"},
        )
        self.assertFalse(any(a.GetAtomicNum() == 0 for a in m.GetAtoms()))
        with tempfile.TemporaryDirectory() as d:
            path = write_colored_svg(m, Path(d), "test")
            svg = (Path(d) / Path(path).name).read_text()
            self.assertIn("#2563EB", svg)
            self.assertIn("#D97706", svg)

    def test_multiple_or_absent_attachments_rejected(self):
        for fragment in ["CO", "*O*"]:
            with self.assertRaises(ValueError):
                assemble("*c1ccccc1", fragment)

    def test_unknown_stereo_not_invented(self):
        m = assemble("*c1ccccc1", "*OCC(C)CCC")
        self.assertEqual(
            [
                v
                for _, v in Chem.FindMolChiralCenters(
                    m, includeUnassigned=True
                )
            ],
            ["?"],
        )

    def test_candidate_is_not_counted_as_confirmed(self):
        cards = [
            {
                "compound_id": "1",
                "si_smiles": "",
                "names": [],
                "has_structure": False,
                "status": "结构待补充",
            }
        ]
        attach(
            cards,
            [{"compound_id": "1", "review_status": "needs_visual_review"}],
        )
        self.assertFalse(cards[0]["has_structure"])
        self.assertTrue(cards[0]["structure_routes"]["markush_candidate"])

    def test_sar_crop_recovers_terminal_and_excludes_table_rule(self):
        import json

        from paper2smiles.extraction.core import fitz

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            pdf = root / "source.pdf"
            doc = fitz.open()
            page = doc.new_page(width=300, height=400)
            page.draw_line((50, 80), (50, 125), width=1)
            page.draw_line((60, 150), (60, 165), width=1)
            page.draw_rect(
                fitz.Rect(20, 170, 240, 170.4), color=None, fill=(0, 0, 0)
            )
            doc.save(pdf)
            doc.close()
            paths, rows, box = render_markush_crops(
                pdf,
                root,
                1,
                {
                    "scaffold_bbox": [100, 100, 300, 300],
                    "rows": [
                        {
                            "compound_id": "18",
                            "fragment_bbox": [100, 350, 300, 425],
                        }
                    ],
                },
                {"18"},
            )
            self.assertGreater(box[3], 125)
            self.assertLess(rows[0]["bbox_pdf_points"][3], 170)
            self.assertTrue(
                all(
                    p.name.endswith("300dpi.png") and p.exists() for p in paths
                )
            )
            self.assertEqual(
                json.loads((root / "sar_render_manifest.json").read_text())[
                    "dpi"
                ],
                300,
            )

    def test_stereo_wedge_at_attachment_can_be_colored(self):
        m = assemble(
            "*[C@H](NC(=O)[C@@H](C)c1ccccc1)c1ccc(OCC(C)CCC)cc1", "*CN"
        )
        before = Chem.MolToSmiles(m, isomericSmiles=True)
        with tempfile.TemporaryDirectory() as d:
            out = write_colored_svg(m, Path(d), "wedge")
            svg = (Path(d) / Path(out).name).read_text()
            self.assertIn("#2563EB", svg)
            self.assertIn("#D97706", svg)
        self.assertEqual(Chem.MolToSmiles(m, isomericSmiles=True), before)

    def test_trial_reuse_requires_exact_pdf(self):
        self.assertEqual(reviewed_examples("different-pdf", {"4"}), [])


if __name__ == "__main__":
    unittest.main()
