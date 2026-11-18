"""One selected structure and its provenance per quantified compound."""
import json,shutil,hashlib,re
from pathlib import Path
from collections import Counter
from core import save,sha256
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from structure_selection import select_structure,canonical

def build(root):
 root=Path(root);dest=root/'report';assets=dest/'assets';downloads=dest/'downloads'
 for d in [dest,assets,downloads]:d.mkdir(parents=True,exist_ok=True)
 def load(name,default):return json.loads((root/name).read_text()) if (root/name).exists() else default
 cards=load('compound_records.json',[]);body=load('body_activity.json',[]);plan=load('plan.json',{});summary=load('summary.json',{});si=load('si_records.json',[])
 selections={c['compound_id']:select_structure(c) for c in cards};manifest=load('si_manifest.json',load('si_manifest.partial.json',[]));si_files={}
 for item in manifest:
  p=Path(item['local_path'])
  if p.is_file() and p.suffix.lower() in ('.csv','.tsv','.xlsx'):
   out=downloads/'si'/p.name;out.parent.mkdir(exist_ok=True);shutil.copy2(p,out);si_files[item['filename']]={**item,'href':'downloads/si/'+p.name}
 for name in ['summary.json','plan.json','compound_records.json','body_activity.json','names_opsin.json','markush_records.json','si_records.json']:
  if (root/name).exists():shutil.copy2(root/name,downloads/name)
 shutil.copy2(root/'input/source.pdf',downloads/'source.pdf')
 selected_images=[]
 for c in cards:
  choice=selections[c['compound_id']]
  if choice['route'] in ('markush','full_structure'):
   selected_images.append({**c,'markush':choice['evidence'] if choice['route']=='markush' else [],'full_structure':choice['evidence'] if choice['route']=='full_structure' else []})
 groups=[]
 if selected_images:
  from table_group_report import prepare
  groups,_=prepare(root,body,selected_images,assets)
 def structure_image(choice):
  mol=Chem.MolFromSmiles(choice['smiles'])
  if choice['route']=='markush':
   from markush_chemistry import assemble_record,write_colored_svg
   assembled=assemble_record(choice['evidence'][0])
   if Chem.MolToSmiles(assembled)!=canonical(choice['smiles']):raise ValueError('Reconstruction mismatch')
   return write_colored_svg(assembled,assets,'selected-'+hashlib.sha256(choice['smiles'].encode()).hexdigest()[:16],source_pdf_hash=sha256(root/'input/source.pdf'))
  from source_depiction import align_to_source
  aligned=align_to_source(mol,sha256(root/'input/source.pdf'))
  drawer=rdMolDraw2D.MolDraw2DSVG(650,360)
  if aligned:drawer.DrawMolecule(rdMolDraw2D.PrepareMolForDrawing(mol,addChiralHs=False,wedgeBonds=False))
  else:rdMolDraw2D.PrepareAndDrawMolecule(drawer,mol)
  drawer.FinishDrawing();svg=drawer.GetDrawingText();name='selected-'+hashlib.sha256(svg.encode()).hexdigest()[:20]+'.svg';(assets/name).write_text(svg);return 'assets/'+name
 def source(card,choice):
  if choice['route']=='markush':
   from table_group_report import compound_source
   return compound_source(root,card,choice['evidence'][0],groups,body,assets)
  if choice['route']=='full_structure':
   from source_overlay import table_overlay
   from core import fitz
   r=choice['evidence'][0];src=Path(r['crop_path']);name='full-'+sha256(src)[:20]+'.png';shutil.copy2(src,assets/name)
   g=next((g for g in groups if g['table_id']==r['table_id']),None);context=[]
   if g:
    with fitz.open(root/'input/source.pdf') as doc:
     pg=doc[r['page']-1];size=r['bbox_page_size'];box=[v*(pg.rect.width/size[0] if i%2==0 else pg.rect.height/size[1]) for i,v in enumerate(r['bbox_pixels'])]
    context=[table_overlay(root,g,[('实际结构截图',box,'#dc4444')],assets)]
   return {'images':context+['assets/'+name],'context_images':context}
  return {'images':[],'context_images':[]}
 from excel_report import write_workbook
 write_workbook(root,body,selections,si,downloads/'results.xlsx')
 from report_data import export_report
 export_report(root,cards,body,plan,selections,si_files,structure_image,source)
 summary['selected_structure_sources']=dict(Counter(c['route'] for c in selections.values()))
 summary['report_mode']='priority_one_structure';save(root/'summary.json',summary);save(downloads/'summary.json',summary)
 save(root/'selected_structures.json',[{'compound_id':cid,'selected_smiles':c['smiles'],'selected_source':c['route']} for cid,c in selections.items()])
