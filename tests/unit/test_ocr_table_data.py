import unittest

from paper2smiles.extraction.ocr_table_data import (
    compound_identifier,
    measurement,
    parse_tables,
    potency_cell,
)


class OCRTableTests(unittest.TestCase):
    def test_nested_header_preserves_plasma_brain_alignment(self):
        grid = parse_tables(
            (
                '<table><tr><td rowspan="2">Compound</td><td '
                'colspan="2">PK</td></tr><tr><td>Plasma</td><td>Brain</td></tr'
                "><tr><td>90</td><td>270</td><td>39</td></tr></table>"
            )
        )[0]
        self.assertEqual(
            [[c["text"] for c in row] for row in grid],
            [
                ["Compound", "PK", "PK"],
                ["Compound", "Plasma", "Brain"],
                ["90", "270", "39"],
            ],
        )

    def test_bound_does_not_invent_ec50(self):
        values = potency_cell("<5.00", "pEC50_with_EC50")
        self.assertEqual(len(values), 1)
        self.assertEqual(values[0]["endpoint"], "pEC50")
        self.assertEqual(values[0]["relation"], "<")

    def test_sem_and_printed_ec50_remain_distinct(self):
        pec, ec = potency_cell("$ 6.05 \\pm 0.07 $ (891)", "pEC50_with_EC50")
        self.assertEqual(
            (pec["mean"], pec["sem"], ec["mean"]), (6.05, 0.07, 891)
        )
        self.assertIsNone(ec["sem"])

    def test_identifier_footnote_requires_definition(self):
        self.assertEqual(compound_identifier("$84^{e}$", {"e"}), ("84", ["e"]))
        with self.assertRaises(ValueError):
            compound_identifier("84'", {"e"})
        with self.assertRaises(ValueError):
            compound_identifier("84^e", set())

    def test_emax_marker_and_not_determined_are_not_values(self):
        data = measurement("$97 \\pm 7^{c}$", "Emax", "%")
        self.assertEqual(
            (data["mean"], data["sem"], data["footnotes"]), (97, 7, ["c"])
        )
        nd = measurement("N.D.", "Emax", "%")
        self.assertEqual(nd["status"], "not_determined")
        self.assertIsNone(nd["mean"])


if __name__ == "__main__":
    unittest.main()
