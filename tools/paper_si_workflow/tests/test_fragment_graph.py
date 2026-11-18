import sys,unittest,tempfile,json,hashlib
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from core import fitz
from rdkit import Chem
from PIL import Image,ImageDraw
from markush_chemistry import assemble_graph,assemble_record
from ruled_fragment_cells import ruled_cells
from fragment_graph_workflow import ink_bounds
from workflow import Workflow

class FragmentGraphTests(unittest.TestCase):
 def test_bridge_and_stereo_preserved(self):
  nodes={'core':'C[C@H](N)[*:1]','bridge':'[*:1]c1ccc([*:2])cc1','tail':'[*:2]O'}
  m=assemble_record({'graph_nodes':nodes})
  self.assertEqual(Chem.MolToSmiles(m),Chem.MolToSmiles(Chem.MolFromSmiles('C[C@H](N)c1ccc(O)cc1')))
  self.assertEqual({a.GetProp('source_node') for a in m.GetAtoms()},set(nodes))
 def test_invalid_edges_rejected(self):
  for nodes in [{'core':'C[*:1]'}, {'core':'C[*:1]','a':'N[*:1]','b':'O[*:1]'}, {'core':'[*:1]CC[*:1]'}, {'core':'CC','b':'O'}]:
   with self.subTest(nodes=nodes),self.assertRaises(ValueError):assemble_graph(nodes)
 def test_merged_cells_and_ocr_count_guard(self):
  with fitz.open() as d:
   p=d.new_page(width=150,height=150)
   for x in [10,40,80,120]:p.draw_line((x,10),(x,100))
   for y in [10,40,100]:p.draw_line((10,y),(120,y))
   p.draw_line((10,70),(40,70));p.draw_line((80,70),(120,70))
   rows=ruled_cells(p,[10,10,120,100],{1:'Ar1',2:'Ar2'},['18','19'])
   self.assertEqual(len(rows),2)
   self.assertEqual(rows[0]['parts']['Ar1']['cell_id'],rows[1]['parts']['Ar1']['cell_id'])
   self.assertNotEqual(rows[0]['parts']['Ar2']['cell_id'],rows[1]['parts']['Ar2']['cell_id'])
   self.assertEqual(ruled_cells(p,[10,10,120,100],{1:'Ar1'},['18']),[])
 def test_blank_and_dash_never_recognized(self):
  im=Image.new('RGB',(100,100),'white');self.assertEqual(ink_bounds(im)[1],'blank_cell')
  ImageDraw.Draw(im).line((40,50,50,50),fill='black',width=2)
  self.assertEqual(ink_bounds(im)[1],'dash_or_small_mark')

class SIDiscoveryTests(unittest.TestCase):
 def test_global_csv_priority_and_download_failure_isolation(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);(root/'si').mkdir();pdf=root/'source.pdf';pdf.write_bytes(b'')
   content=b'id,SMILES\n1,CC\n';(root/'si/20_good.csv').write_bytes(content)
   doi='10.1/example';events=[]
   metas={1:{'resource_doi':doi,'files':[{'id':10,'name':'large.pdf','size':100,'download_url':'https://invalid/pdf'}]},2:{'resource_doi':doi,'files':[{'id':20,'name':'good.csv','size':len(content),'computed_md5':hashlib.md5(content).hexdigest(),'download_url':'https://invalid/csv'}]},3:{'resource_doi':doi,'files':[{'id':30,'name':'bad.csv','size':100,'download_url':'https://invalid/bad'}]}}
   class Response:
    def __init__(self,value):self.value=value
    def json(self):return self.value
   def request(method,url,**kwargs):
    if method=='POST':return Response([{'id':k,'resource_doi':doi} for k in metas])
    return Response(metas[int(url.rsplit('/',1)[1])])
   flow=object.__new__(Workflow);flow.root=root;flow.pdf=pdf;flow.warnings=[];flow.status=lambda *a,**k:events.append(k)
   with patch('workflow.public_request',side_effect=request),patch('workflow.requests.Session.get',side_effect=TimeoutError('test failure')):
    attachments=flow.discover_si({'doi':doi})
   self.assertEqual([a['filename'] for a in attachments],['good.csv'])
   self.assertTrue(any('bad.csv' in s for s in flow.warnings))
   self.assertTrue((root/'si_optional_downloads.json').exists())

if __name__=='__main__':unittest.main()
