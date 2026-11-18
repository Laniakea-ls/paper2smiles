import unittest

from paper2smiles.extraction.activity_targets import select_targets


class ActivityTargetsTests(unittest.TestCase):
    def test_bioactivity_is_required_and_limits_are_preserved(self):
        targets, excluded = select_targets(
            [
                {"compound_id": "1", "endpoint": "EC50", "value": "59"},
                {"compound_id": "2", "endpoint": "pEC50", "value": "<5.00"},
                {"compound_id": "3", "endpoint": "Emax", "value": "N.D."},
                {"compound_id": "4", "endpoint": "clogP", "value": "3.4"},
                {"compound_id": "5", "measurements": []},
            ]
        )
        self.assertEqual(set(targets), {"1", "2"})
        self.assertFalse(targets["1"]["only_limits"])
        self.assertTrue(targets["2"]["only_limits"])
        self.assertEqual({r["compound_id"] for r in excluded}, {"3", "4", "5"})

    def test_repeated_assays_form_one_target_without_losing_rows(self):
        rows = [
            {
                "compound_id": "84",
                "measurements": [
                    {"endpoint": "EC50", "mean": 59, "relation": "="}
                ],
            },
            {
                "compound_id": "84",
                "measurements": [
                    {"endpoint": "EC50", "mean": 942, "relation": "="}
                ],
            },
        ]
        targets, _ = select_targets(rows)
        self.assertEqual(len(targets), 1)
        self.assertEqual(len(targets["84"]["activity_records"]), 2)


if __name__ == "__main__":
    unittest.main()
