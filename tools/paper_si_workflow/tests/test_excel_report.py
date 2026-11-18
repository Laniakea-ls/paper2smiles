import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from excel_report import merged_activity,number_or_text
class ExcelTests(unittest.TestCase):
 def test_pivot_preserves_units_assays_and_repeats(self):
  base={'compound_id':'8','target':'hGPR88','assay':'BRET','page':4,'caption':'Table 1','source_row':2}
  rows=[dict(base,endpoint='EC50',unit='nM',value='518'),dict(base,endpoint='Emax',unit='%',value='99'),dict(base,endpoint='EC50',unit='μM',value='0.518'),dict(base,endpoint='EC50',unit='nM',value='520'),dict(base,assay='other',endpoint='EC50',unit='nM',value='600')]
  cols,groups=merged_activity(rows)
  self.assertEqual(len(cols),3);self.assertEqual(len(groups),2)
  self.assertEqual(number_or_text(groups[0]['values'][('EC50','nM')]),'518\n520')
  self.assertEqual(sum(len(v) for g in groups for v in g['values'].values()),5)
 def test_limits_and_missing(self):
  _,groups=merged_activity([{'compound_id':'7','endpoint':'EC50','unit':'nM','value':'99000','relation':'>'}])
  self.assertEqual(number_or_text(groups[0]['values'][('EC50','nM')]),'>99000')
  self.assertEqual(number_or_text(None),'—');self.assertEqual(number_or_text(['518']),518)
if __name__=='__main__':unittest.main()
