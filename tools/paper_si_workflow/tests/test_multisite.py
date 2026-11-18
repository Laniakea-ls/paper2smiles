import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from rdkit import Chem
from markush_chemistry import assemble_multiple,write_colored_svg
from multisite_reviewed import CORE
from source_depiction import align_to_source,PAPER_SHA256
class MultipleSitesTests(unittest.TestCase):
 def test_two_substituents_are_independent(self):
  a=assemble_multiple(CORE,{'4':'*C','5':'*F'},{'13':'N'})
  b=assemble_multiple(CORE,{'5':'*F','4':'*C'},{'13':'N'})
  self.assertEqual(Chem.MolToSmiles(a),Chem.MolToSmiles(b))
  swapped=assemble_multiple(CORE,{'4':'*F','5':'*C'},{'13':'N'})
  self.assertNotEqual(Chem.MolToSmiles(a),Chem.MolToSmiles(swapped))
 def test_ring_n_and_stereo_roundtrip(self):
  m=assemble_multiple(CORE,{'4':'*C','5':'*[H]'},{'13':'N'})
  self.assertEqual(sum(a.GetAtomicNum()==7 and a.GetIsAromatic() for a in m.GetAtoms()),1)
  self.assertFalse(any(a.GetAtomicNum()==0 or a.GetAtomMapNum() for a in m.GetAtoms()))
  expected=Chem.MolToSmiles(m);self.assertTrue(align_to_source(m,PAPER_SHA256))
  self.assertEqual(Chem.MolToSmiles(Chem.MolFromMolBlock(Chem.MolToMolBlock(m))),expected)
 def test_missing_site_is_rejected(self):
  with self.assertRaises(ValueError):assemble_multiple(CORE,{'4':'*C'},{'13':'N'})
if __name__=='__main__':unittest.main()
