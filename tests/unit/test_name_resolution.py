import unittest

from paper2smiles.extraction.core import mol_info
from paper2smiles.extraction.name_resolution import (
    audit_names,
    name_status,
    xml_names,
)


class NameResolutionTests(unittest.TestCase):
    def test_parseable_partial_name_is_not_complete(self):
        self.assertEqual(
            name_status(mol_info("NCCO"), mol_info("CC(=O)NCCO")),
            "partial_structure_name",
        )

    def test_salt_not_partial(self):
        self.assertEqual(
            name_status(mol_info("CCN.Cl"), mol_info("CCN")),
            "salt_form_difference",
        )

    def test_stereo_not_partial(self):
        self.assertEqual(
            name_status(mol_info("CC(O)C(=O)O"), mol_info("C[C@H](O)C(=O)O")),
            "stereo_difference",
        )

    def test_reference_names_preserve_local_id_and_full_heading(self):
        result = xml_names(
            (
                "<article><sec><title>5.1.9. (2S)-N-{(1R)-2-Hydroxy-1-[4-(2-me"
                "thylpentyloxy)phenyl]ethyl}-2-phenylpropanamide "
                "(3)</title></sec></article>"
            )
        )
        self.assertEqual(result[0]["original_compound_id"], "3")
        self.assertTrue(result[0]["name"].startswith("(2S)-N-"))
        self.assertTrue(result[0]["name"].endswith("phenylpropanamide"))

    def test_partial_quarantined_missing_queued(self):
        names = [
            {"compound_id": "2", "name": "partial", "opsin": mol_info("NCCO")}
        ]
        issues = audit_names(
            [
                {"compound_id": "2", "si_structure": mol_info("CC(=O)NCCO")},
                {"compound_id": "3", "si_structure": mol_info("CCO")},
            ],
            names,
        )
        self.assertFalse(names[0]["eligible_for_structure_comparison"])
        self.assertEqual(
            {i["reason"] for i in issues},
            {"partial_structure_name", "name_missing"},
        )


if __name__ == "__main__":
    unittest.main()
