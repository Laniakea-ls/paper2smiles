"""Finalize the source-reviewed BI-9508 experiment; not an automatic recognizer.
Requires the frozen pre-reference graph trial and DOI-verified author CSV.
"""
import sys,json,csv,hashlib,copy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from core import save,sha256,fitz,read_tables
from markush_chemistry import assemble_graph
from rdkit import Chem
from markush_workflow import attach
from full_structure_workflow import SAR
sys.path.insert(0,str(SAR))
from pipeline.smiles_validation.ocsr_render import render_inputs
p=Path('runs/gpr88_paper_batch/paper-13').resolve();out=p/'complex_trial/optimized'
digest=sha256(p/'input/source.pdf');frozen=json.loads((out/'source_graph_map_frozen.json').read_text());assert frozen['source_pdf_sha256']==digest
backup=p/'complex_trial/pre_update';backup.mkdir(exist_ok=True)
for name in ['compound_records.json','summary.json','names_opsin.json','markush_records.json']:
 if (p/name).exists() and not (backup/name).exists():(backup/name).write_bytes((p/name).read_bytes())
manifest=[]
for article,num in [(26140207,2),(26140210,3)]:
 meta=json.loads(Path(f'/tmp/bi9508-si{num}.json').read_text());assert meta['resource_doi']=='10.1021/acs.jmedchem.4c00665'
 save(p/'si'/f'article_{article}.json',meta)
 for f in meta['files']:
  path=p/'si'/f"{f['id']}_{f['name']}";md5=f.get('computed_md5') or f.get('supplied_md5')
  assert path.stat().st_size==f['size'] and hashlib.md5(path.read_bytes()).hexdigest()==md5
  tables=read_tables(path)
  manifest.append({'file_id':str(f['id']),'article_id':article,'filename':f['name'],'local_path':str(path),'article_doi':meta['resource_doi'],'supplement_doi':meta['doi'],'source_url':f['download_url'],'landing_url':meta['url_public_html'],'metadata_url':f'https://api.figshare.com/v2/articles/{article}','sha256':sha256(path),'md5':md5,'size':f['size'],'download_verified':True,'table_previews':[{'sheet':t['sheet'],'row_count':len(t['rows']),'preview_rows':t['rows'][:8]} for t in tables]})
save(p/'si_manifest.json',manifest)
file=next(x for x in manifest if x['filename'].endswith('_003.csv'))
with open(file['local_path'],encoding='utf-8-sig',newline='') as stream:csv_rows=list(csv.reader(stream,delimiter=';'))
header=csv_rows[0];refs={};si=[]
for rowno,row in enumerate(csv_rows[1:],2):
 if len(row)<7:continue
 cid=row[0];smi=Chem.MolToSmiles(Chem.MolFromSmiles(row[1]),isomericSmiles=True)
 refs.setdefault(cid,set()).add(smi)
 for col in range(6,len(header),2):
  if col>=len(row) or not row[col]:continue
  si.append({'compound_id':cid,'si_smiles':row[1],'canonical_smiles':smi,'source_type':'author_si_csv','source_file':file['filename'],'source_row':rowno,'source_column':col,'endpoint_header':header[col],'modifier':row[col-1],'raw_value':row[col],'value':row[col].replace(',','.'),'source_url':file['source_url'],'source_sha256':file['sha256']})
save(p/'si_records.json',si)
cards=json.loads((backup/'compound_records.json').read_text());names=json.loads((backup/'names_opsin.json').read_text())
for c in cards:
 values=refs.get(c['compound_id'],set());assert len(values)==1,(c['compound_id'],values)
 c.update(si_smiles=next(iter(values)),si_source=file['filename'],has_structure=True,status='作者 SI 结构')
 c['si_evidence']=[x for x in si if x['compound_id']==c['compound_id']]
 if c['compound_id']=='30':
  for n in c['names']:n.update(eligible_for_structure_comparison=False,exclusion_reason='PDF 27 的标题重复胺前体名称，缺失酰胺骨架；作者 SI、合成步骤及 HRMS 支持完整产物。')
  c['structure_note']='原文名称缺少酰胺骨架；采用作者 SI。'
for n in names:
 if n.get('compound_id')=='30':n.update(eligible_for_structure_comparison=False,exclusion_reason='Printed heading describes the amine precursor; author SI and experimental HRMS support the complete amide.')
save(p/'names_opsin.json',names)
cells=json.loads((out/'tight_cells.json').read_text());r14=next(c for c in cells if c['key']=='4-R-301.76-329.65')
with fitz.open(p/'input/source.pdf') as pdf:
 specs=[{'crop':'common-core.png','page':3,'bbox_pixels':[325,66,450,145],'bbox_page_size':[pdf[2].rect.width,pdf[2].rect.height]}]
 mapping,audit=render_inputs(p/'input/source.pdf',specs,out,out/'core',dpi=300)
core_evidence={'node_label':'通式','page':3,'tight_bbox_pdf':audit[0]['pdf_bbox'],'crop_path':next(iter(mapping.values()))}
save(out/'core_render_manifest.json',audit)
records=json.loads((out/'graph_pre_reference.json').read_text());comparisons=[]
for r in records:
 r['source_candidate_smiles']=r['assembled_smiles'];r['source_graph_nodes']=copy.deepcopy(r['graph_nodes'])
 if r['page']==5:
  r['graph_nodes']['R']='[*:1]CCOC(F)F'
  r['cross_source_resolution']={'node':'R','original':'[*:1]CCCOC(F)F','resolved':'[*:1]CCOC(F)F','source':'Table 1 compound 14, complete experimental names, independently checked against author SI','source_table':'page-004-t1','source_cell_key':r14['key']}
  r['graph_evidence'].append({**r14,'node_label':'R · Table 1 / 14'})
  r['display_note']='侧链按 Table 1 与 SI 校核；保留原图差异。'
 r['graph_evidence'].insert(0,core_evidence)
 r['assembled_smiles']=Chem.MolToSmiles(assemble_graph(r['graph_nodes']),isomericSmiles=True)
 exact=r['assembled_smiles'] in refs.get(r['compound_id'],set());assert exact,r['compound_id']
 r.update(review_status='source_visually_reviewed',review_note='Source-based local recognition plus manual attachment/ring review. Resolved Table 2 fixed-chain discrepancy using Table 1 and names; subsequently validated against author SI. Not fully automatic.',si_validation='exact_isomeric_match',table_context_bbox_pdf=[51,58,555,357] if r['page']==4 else [51,58,555,738])
 comparisons.append({'compound_id':r['compound_id'],'table_id':r['table_id'],'original_source_candidate_si_match':r['source_candidate_smiles'] in refs[r['compound_id']],'final_si_match':exact,'cross_source_resolution':r.get('cross_source_resolution')})
attach(cards,records);save(p/'compound_records.json',cards);save(p/'markush_records.json',records)
save(p/'markush/source_reviewed_records.json',{'source_pdf_sha256':digest,'records':records})
save(out/'independent_si_validation.json',comparisons)
summary=json.loads((backup/'summary.json').read_text());summary.update(with_structure=27,with_si_structure=27,with_names=23,with_markush_structure=len({r['compound_id'] for r in records}),si_measurements=len(si),warnings=['Table 2 公共侧链图与 Table 1、完整名称和作者 SI 不一致。','化合物 30 原文名称疑为胺前体，未作为完整结构采用。','SI 中化合物 6 有两条同结构、不同活性记录，原样保留；不覆盖正文。'])
summary['complex_trial']={'baseline_graph_structures':0,'assembled_unique_compounds':21,'assembled_table_rows':22,'exact_author_si_matches_after_source_resolution':21,'manual_source_review_required':True,'blocked_table_rows':json.loads((out/'graph_failures.json').read_text()),'not_fully_automatic':True}
save(p/'summary.json',summary)
from batch_report import build
build(p)
print(json.dumps(summary,ensure_ascii=False,indent=2))
