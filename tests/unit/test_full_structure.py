import tempfile
import unittest
from pathlib import Path

from paper2smiles.extraction.full_structure_workflow import (
    attach,
    complete_smiles,
    consensus,
    render_with_edge_check,
    valid_box,
)


class FullStructureTests(unittest.TestCase):
    def test_r_groups_are_not_full_structures(self):
        self.assertEqual(complete_smiles("*c1ccccc1"), "")
        self.assertEqual(complete_smiles(""), "")

    def test_two_exact_votes_keep_dissent(self):
        p = {
            "A": {"smiles": "CCO"},
            "B": {"smiles": "OCC"},
            "C": {"smiles": "COC"},
        }
        r = consensus(p)
        self.assertEqual(r["status"], "two_agree")
        self.assertEqual(r["agreeing_models"], ["A", "B"])
        self.assertEqual(p["C"]["smiles"], "COC")

    def test_stereo_difference_is_not_agreement(self):
        r = consensus(
            {
                "A": {"smiles": "N[C@H](C)O"},
                "B": {"smiles": "N[C@@H](C)O"},
                "C": {"smiles": ""},
            }
        )
        self.assertEqual(r["status"], "incomplete")
        self.assertFalse(r["requires_manual_review"])
        self.assertFalse(r["candidate_smiles"])

    def test_failures_and_unassigned_stereo_are_not_votes(self):
        r = consensus(
            {
                "A": {"smiles": "CCO", "error": "failed"},
                "B": {"smiles": "CCO"},
                "C": {"smiles": "*CCO"},
            }
        )
        self.assertEqual(r["valid_models"], 1)
        self.assertFalse(r["candidate_smiles"])

    def test_majority_policy(self):
        for values, status, accepted, manual in [
            (["CCO", "OCC", "COC"], "two_agree", True, False),
            (["CCO", "OCC", "CCO"], "three_agree", True, False),
            (["CCO", "OCC", ""], "two_agree", True, False),
            (["CCO", "COC", "CCC"], "needs_review", False, True),
            (["CCO", "", ""], "incomplete", False, False),
        ]:
            with self.subTest(values=values):
                vote = consensus(
                    {k: {"smiles": v} for k, v in zip("ABC", values)}
                )
                self.assertEqual(
                    (
                        vote["status"],
                        vote["accepted"],
                        vote["requires_manual_review"],
                    ),
                    (status, accepted, manual),
                )
                self.assertEqual(
                    vote["selected_smiles"], "CCO" if accepted else ""
                )

    def test_attach_adopts_majority_and_preserves_dissent(self):
        cards = [{"compound_id": "31", "has_structure": False}]
        records = [
            {
                "compound_id": "31",
                "predictions": {
                    "A": {"smiles": "CCO"},
                    "B": {"smiles": "OCC"},
                    "C": {"smiles": "COC"},
                },
            }
        ]
        attach(cards, records)
        self.assertTrue(cards[0]["has_structure"])
        self.assertEqual(cards[0]["full_structure_smiles"], "CCO")
        self.assertEqual(records[0]["review_status"], "auto_accepted")
        self.assertEqual(records[0]["predictions"]["C"]["smiles"], "COC")

    def test_crop_coordinates_reject_invalid_ranges(self):
        self.assertTrue(valid_box([10, 20, 100, 200]))
        self.assertFalse(valid_box([10, 20, 1001, 200]))
        self.assertFalse(valid_box([100, 20, 10, 200]))

    def test_crop_with_cut_bond_is_expanded_before_readers(self):
        from PIL import Image, ImageDraw

        from paper2smiles.extraction.core import fitz

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pdf = root / "source.pdf"
            doc = fitz.open()
            doc.new_page(width=600, height=800)
            doc.save(pdf)
            doc.close()
            calls = []

            def renderer(pdf, rows, folder, destination, dpi):
                destination.mkdir(exist_ok=True)
                image = Image.new("RGB", (100, 100), "white")
                if not calls:
                    ImageDraw.Draw(image).line(
                        (50, 50, 50, 99), fill="black", width=2
                    )
                path = destination / "crop.png"
                image.save(path)
                calls.append(list(rows[0]["bbox_pixels"]))
                return {
                    str((folder / rows[0]["crop"]).resolve()): str(path)
                }, []

            rows = [
                {
                    "page": 1,
                    "compound_id": "30",
                    "crop": "x.png",
                    "bbox_pixels": [100, 100, 200, 200],
                }
            ]
            render_with_edge_check(pdf, rows, root, renderer)
            self.assertEqual(len(calls), 2)
            self.assertGreater(calls[1][3], calls[0][3])
            self.assertEqual(calls[1][:3], calls[0][:3])


if __name__ == "__main__":
    unittest.main()
