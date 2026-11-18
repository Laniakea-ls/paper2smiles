import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from generic_ocr import cell_values,activity_identifier
class GenericOCRTests(unittest.TestCase):
 def test_footnote_limit_does_not_become_exact_value(self):
  result=cell_values('NA^b',{'endpoint':'pEC50','encoding':'pEC50_with_EC50','unit':'nM'},{'footnotes':[{'marker':'b','text':'EC50 > 10 μM, two experiments.'}]})
  self.assertEqual([(r['endpoint'],r['relation'],r['mean'],r['unit']) for r in result],[('EC50','>',10,'μM')])
 def test_missing_footnote_does_not_guess_bound(self):
  with self.assertRaises(ValueError):cell_values('NA',{'endpoint':'EC50'},{'footnotes':[]})
 def test_alias_and_letter_suffix(self):
  self.assertEqual(activity_identifier('2-PCCA (1)',set()),('1',[],'2-PCCA'))
  self.assertEqual(activity_identifier('13e',{'e'}),('13e',[],''))
 def test_parenthesized_sd_and_thousands_separator(self):
  a=cell_values('500 ( ± 90)',{'endpoint':'EC50','encoding':'plain','unit':'nM'},{})[0]
  self.assertEqual(a['value_original'],'500 ( ± 90)');self.assertEqual(a['mean'],500)
  b=cell_values('>10,000',{'endpoint':'EC50','encoding':'plain','unit':'nM'},{})[0]
  self.assertEqual((b['relation'],b['mean']),('>',10000))
 def test_paired_potency_and_efficacy(self):
  rows=cell_values('301 (93)',{'endpoint':'EC50','encoding':'EC50_with_Emax','unit':'nM'},{})
  self.assertEqual([(v['endpoint'],v['mean'],v['unit']) for v in rows],[('EC50',301,'nM'),('Emax',93,'%')])
 def test_ocr_latex_grouping_preserves_numbers(self):
  rows=cell_values(r'{6.27} ± {0.01}\: ( {{537}} )',{'endpoint':'pEC50','encoding':'pEC50_with_EC50','unit':'nM'},{})
  self.assertEqual([(v['endpoint'],v['mean']) for v in rows],[('pEC50',6.27),('EC50',537)])
 def test_racemic_identifier_not_discarded(self):
  self.assertEqual(activity_identifier('( ± )-28',set()),('(±)-28',[],''))
 def test_parenthesized_bound_without_log_value(self):
  rows=cell_values('(>10,000)',{'endpoint':'pEC50','encoding':'pEC50_with_EC50','unit':'nM'},{})
  self.assertEqual([(v['endpoint'],v['mean'],v['relation']) for v in rows],[('EC50',10000,'>')])
 def test_plain_parentheses_require_error_schema(self):
  row=cell_values('270 (30)',{'endpoint':'EC50','encoding':'plain','unit':'nM','error_type':'SD'},{})[0]
  self.assertEqual((row['mean'],row['sem'],row['value_original']),(270,30,'270 (30)'))
  with self.assertRaises(ValueError):cell_values('270 (30)',{'endpoint':'EC50','encoding':'plain','unit':'nM'},{})
 def test_cell_unit_is_kept(self):
  result=cell_values('0.2 ± 0.1',{'endpoint':'EC50','encoding':'plain','unit':'μM'},{})
  self.assertEqual(result[0]['unit'],'μM');self.assertEqual(result[0]['mean'],.2)
if __name__=='__main__':unittest.main()
