import unittest

from rdkit import Chem

from paper2smiles.extraction.source_depiction import (
    PAPER_SHA256,
    align_to_source,
)


class SourceDepictionTests(unittest.TestCase):
    def test_wedges_round_trip_and_unknown_stereo_stays_unknown(self):
        for smi in [
            "CC(C)COc1ccc([C@H](CO)NC(=O)[C@@H](C)c2ccccc2)cc1",
            "CCCC(C)COc1ccc([C@H](C(C)(C)O)NC(=O)[C@@H](C)c2ccccc2)cc1",
            "CCCC(C)COc1ccc([C@@H](CN)NC(=O)[C@@H](C)c2ccccc2)cc1",
        ]:
            with self.subTest(smi=smi):
                m = Chem.MolFromSmiles(smi)
                expected = Chem.MolToSmiles(m)
                self.assertTrue(align_to_source(m, PAPER_SHA256))
                restored = Chem.MolFromMolBlock(Chem.MolToMolBlock(m))
                self.assertEqual(Chem.MolToSmiles(restored), expected)
                match = m.GetSubstructMatch(
                    Chem.MolFromSmiles("Oc1ccc(CNC=O)cc1")
                )
                conf = m.GetConformer()
                self.assertGreater(
                    conf.GetAtomPosition(match[6]).y,
                    conf.GetAtomPosition(match[5]).y,
                )
                bond = m.GetBondBetweenAtoms(match[5], match[6])
                self.assertEqual(bond.GetBeginAtomIdx(), match[5])
                self.assertIn(
                    bond.GetBondDir(),
                    [Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH],
                )

    def test_other_sources_not_forced_into_this_layout(self):
        m = Chem.MolFromSmiles("Oc1ccc(CNC=O)cc1")
        self.assertFalse(align_to_source(m, "other-pdf"))
        self.assertEqual(m.GetNumConformers(), 0)


if __name__ == "__main__":
    unittest.main()
