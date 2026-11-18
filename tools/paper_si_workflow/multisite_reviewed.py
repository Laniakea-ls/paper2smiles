"""Source-reviewed Table 4 map, reused only for this exact PDF and printed rows."""
from rdkit import Chem
from markush_chemistry import assemble_multiple
from source_depiction import PAPER_SHA256
CORE='CC[C@H](C)COc1cc([*:5])c([C@H](CO)NC(=O)[C@@H](C)c2ccccc2)c([*:4])[cH:13]1'
FRAGMENTS={'H':'*[H]','Me':'*C','F':'*F','OMe':'*OC'}
EXPECTED={'5':['H','H','CH'],'24':['Me','H','CH'],'25':['F','H','CH'],'26':['OMe','H','CH'],'27':['Me','Me','CH'],'28':['Me','H','N']}
def reviewed_multisite(pdf_hash,body,targets):
 if pdf_hash!=PAPER_SHA256:return []
 out=[]
 for row in body:
  cid=row['compound_id'];raw=row.get('raw_row',[])
  if row.get('table_id')!='page-003-t3' or cid not in targets or cid not in EXPECTED:continue
  if raw[1:4]!=EXPECTED[cid]:raise ValueError('Table 4 printed row differs from reviewed site map: '+cid)
  fragments={'4':FRAGMENTS[raw[1]],'5':FRAGMENTS[raw[2]]};variants={'13':raw[3]}
  m=assemble_multiple(CORE,fragments,variants)
  out.append({'compound_id':cid,'table_id':'page-003-t3','page':3,'source_type':'main_pdf_multisite_reconstruction','source_pdf_sha256':pdf_hash,'scaffold_smiles':CORE,'fragments':fragments,'atom_variants':variants,'substituents':{'R4':raw[1],'R5':raw[2],'X':raw[3]},'assembled_smiles':Chem.MolToSmiles(m,isomericSmiles=True),'review_status':'source_visually_reviewed','review_note':'按原图 R4、R5 和环位 X 分别映射，按原表单元格连接；化合物 5 与 Table 1 校验一致。','stereo_note':'保留原图明确的三个立体中心。','source_row':row.get('source_row'),'variable_label':'R4 / R5 / X'})
 return out
