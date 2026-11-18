import sys
from pathlib import Path
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from core import normalize_id,identifier_present,source_span,opsin_batch,potency_nm,mol_info

class EvidenceTests(unittest.TestCase):
    def test_identifiers_do_not_alias_footnotes_or_substrings(self):
        self.assertEqual(normalize_id('Compound 84'),'84')
        self.assertNotEqual(normalize_id("84'"),'84')
        self.assertTrue(identifier_present('84','amide (84).'))
        self.assertFalse(identifier_present('84','amide (184).'))
    def test_evidence_ranges_are_bounded(self):
        pages={1:{'lines':[{'text':'ethanol (1)','bbox':[1,2,3,4]}]}}
        self.assertEqual(source_span([{'page':1,'start_line':1,'end_line':1}],pages)[0],'ethanol (1)')
        with self.assertRaises(ValueError):source_span([{'page':1,'start_line':1,'end_line':2}],pages)
    def test_potency_bounds_and_units(self):
        self.assertEqual(potency_nm('<5.00','pEC50','')['relation'],'>')
        self.assertEqual(potency_nm('<5.00','pEC50','')['value'],10000)
        self.assertIsNone(potency_nm('6.2','pEC50','nM'))
        self.assertEqual(potency_nm('0.059','EC50','uM')['value'],59)
    def test_opsin_row_alignment_and_stereochemistry(self):
        with tempfile.TemporaryDirectory() as folder:
            rows=opsin_batch(['ethanol','not-a-valid-chemical-name-zzq','(R)-butan-2-ol','(S)-butan-2-ol'],Path(folder))
        self.assertEqual(len(rows),4)
        self.assertEqual(rows[0]['canonical'],mol_info('CCO')['canonical'])
        self.assertFalse(rows[1]['valid'])
        self.assertTrue(rows[2]['valid'] and rows[3]['valid'])
        self.assertNotEqual(rows[2]['canonical'],rows[3]['canonical'])

if __name__=='__main__':unittest.main()
