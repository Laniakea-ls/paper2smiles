import tempfile
import unittest
from pathlib import Path

from paper2smiles.extraction.core import read_tables
from paper2smiles.extraction.structure_selection import select_structure
from paper2smiles.extraction.workflow import Workflow


class PriorityTests(unittest.TestCase):
    def test_first_confirmed_route_only(self):
        c = {
            "si_smiles": "CC",
            "names": [
                {
                    "opsin": {"valid": True, "smiles": "CO"},
                    "name_text_supported": True,
                }
            ],
            "markush": [
                {
                    "review_status": "source_visually_reviewed",
                    "assembled_smiles": "CN",
                }
            ],
            "full_structure": [
                {
                    "consensus": {
                        "status": "two_agree",
                        "candidate_smiles": "CF",
                    }
                }
            ],
        }
        self.assertEqual(select_structure(c)["route"], "si")
        c["si_smiles"] = ""
        self.assertEqual(select_structure(c)["route"], "name")
        c["names"][0]["eligible_for_structure_comparison"] = False
        self.assertEqual(select_structure(c)["route"], "markush")
        c["markush"][0]["review_status"] = "needs_visual_review"
        self.assertEqual(select_structure(c)["route"], "full_structure")
        c["full_structure"][0]["consensus"]["status"] = "needs_review"
        self.assertEqual(select_structure(c)["route"], "missing")

    def test_conflicting_names_not_arbitrarily_selected(self):
        c = {
            "names": [
                {
                    "opsin": {"valid": True, "smiles": s},
                    "name_text_supported": True,
                }
                for s in ["CO", "CN"]
            ]
        }
        self.assertEqual(select_structure(c)["route"], "missing")

    def test_semicolon_si_structure_without_activity(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "source.csv"
            p.write_text(
                'ID;SMILES;Value\n23;CCO;"1,2"\n', encoding="utf-8-sig"
            )
            self.assertEqual(
                read_tables(p)[0]["rows"][1], ["23", "CCO", "1,2"]
            )
            flow = object.__new__(Workflow)
            flow.warnings = []
            rows = flow.si_records(
                {
                    "tables": [
                        {
                            "file_id": "1",
                            "sheet": "",
                            "header_row": 0,
                            "data_start_row": 1,
                            "compound_id_column": 0,
                            "smiles_column": 1,
                            "measurements": [],
                        }
                    ]
                },
                [
                    {
                        "file_id": "1",
                        "local_path": str(p),
                        "filename": "source.csv",
                        "source_url": "https://example.org/file",
                    }
                ],
            )
            self.assertEqual(rows[0]["compound_id"], "23")
            self.assertTrue(rows[0]["si_structure"]["valid"])
            self.assertEqual(rows[0]["value"], "")


if __name__ == "__main__":
    unittest.main()
