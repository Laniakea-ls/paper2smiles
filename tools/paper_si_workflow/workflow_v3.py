"""SI → source names → scaffold/fragment reconstruction; numeric activity gates cards."""
import json,re,shutil,time,hashlib
from collections import defaultdict,Counter
from pathlib import Path
from core import save,sha256,extract_pdf,parse_number,mol_info
from workflow import Workflow,public_request
from generic_ocr import extract_tables
from name_resolution import recover_cited_names,audit_names

class ActivityWorkflow(Workflow):
 def select_current_tables(self,plan,attachments):
  signature=hashlib.sha256(json.dumps([(a['file_id'],a['sha256']) for a in attachments],sort_keys=True).encode()).hexdigest()[:12]
  value=self.cached('si_selection_'+signature+'.json',lambda:self.select_tables(plan,attachments))
  save(self.root/'si_selection.json',value)
  return value
 def run(self):
  digest=sha256(self.pdf)
  if (self.root/'metadata.json').exists() and json.loads((self.root/'metadata.json').read_text())['pdf_sha256']!=digest:raise ValueError('Source mismatch')
  save(self.root/'metadata.json',{'workflow_version':'activity-first-priority-v7','pdf_sha256':digest,'input_pdf':str(self.pdf),'model':self.llm.model})
  (self.root/'input').mkdir(exist_ok=True);shutil.copy2(self.pdf,self.root/'input/source.pdf')
  document=self.cached('pdf_text.json',lambda:extract_pdf(self.pdf));self.pages={p['page']:p for p in document['pages']}
  self.stage='plan';self.status('识别正文与活性页')
  plan=self.cached('plan.json',lambda:self.plan_document(document))
  additional=[]
  for page in self.pages.values():
   if re.search(r'(?:EC\s*50|IC\s*50|pEC\s*50).{0,60}\d.{0,20}(?:nM|[µμu]M)',page['text'],re.I|re.S) and not re.search(r'^.*REFERENCES',page['text'][:100]):
    if page['page'] not in plan['main_activity_pages']:additional.append(page['page'])
  if additional:
   plan={**plan,'main_activity_pages':sorted(set(plan['main_activity_pages']+additional))}
   save(self.root/'plan_activity_refinement.json',{'additional_pages_with_printed_potency':additional})
  self.stage='activity_ocr';self.status('核对正文活性，确定需要展示的化合物编号')
  if (self.root/'body_activity_ocr_v8.json').exists():body=json.loads((self.root/'body_activity_ocr_v8.json').read_text())
  else:
   body,tables=extract_tables(self,plan)
   covered={t['page'] for t in tables if t['rows']}
   remaining=[p for p in plan['main_activity_pages'] if p not in covered]
   if remaining:
    extra=self.body_activity({**plan,'main_activity_pages':remaining})
    for row in extra:row['source_type']='main_pdf_vision_no_table';row['needs_visual_review']=True
    body+=extra
   save(self.root/'body_activity_ocr_v8.json',body)
  # Normalize explicitly parenthesized IDs also in vision-derived records.
  for r in body:
   match=re.fullmatch(r'(.+?)\s*\((\d+[A-Za-z]?)\)',r['compound_id'])
   if match:
    r['raw_identifier']=r['compound_id'];r['alias']=match[1].strip();r['compound_id']=match[2].upper()
  # A range with a pooled value is not an individual compound measurement.
  grouped=[r for r in body if re.fullmatch(r'\d+\s*[–—-]\s*\d+',r['compound_id'])]
  if grouped:
   save(self.root/'grouped_activity_review.json',grouped)
   self.warnings.append('合并编号范围的活性保留为待核对证据，不当作单个化合物。')
   body=[r for r in body if r not in grouped]
  from identifier_review import review_identifiers
  body=review_identifiers(self,body)
  alias_candidates=defaultdict(set)
  for r in body:
   if r.get('alias'):alias_candidates[r['alias'].upper()].add(r['compound_id'])
  for r in body:
   ids=alias_candidates.get(r['compound_id'].upper(),set())
   if len(ids)==1:r['compound_id']=next(iter(ids))
  targets={r['compound_id'] for r in body if (parse_number(r['value']) is not None or isinstance(r.get('mean'),(int,float))) and r['compound_id']}
  body=[r for r in body if r['compound_id'] in targets]
  save(self.root/'body_activity.json',body)
  self.status('正文活性初步完成',compounds=len(targets),measurements=len(body))
  if (self.root/'activity_ocr/tables.json').exists():
   table_audit=json.loads((self.root/'activity_ocr/tables.json').read_text())
   rejected=sum(len(t.get('rejected',[])) for t in table_audit)
   if rejected:self.warnings.append(f'OCR 表格有 {rejected} 个单元格或行未自动解析，详见原始 OCR 审计。')
  attachments=[];si=[]
  self.stage='si_discovery';self.status('第一步：检索并读取作者 SI')
  try:
   attachments=self.cached('si_manifest.json',lambda:self.discover_si(plan))
   if attachments:
    selection=self.select_current_tables(plan,attachments)
    si=self.si_records(selection,attachments);save(self.root/'si_records.json',si)
  except Exception as e:self.warnings.append('SI 获取或表格读取未完成：'+str(e))
  if not attachments and (self.root/'si_manifest.partial.json').exists():
   attachments=json.loads((self.root/'si_manifest.partial.json').read_text())
   try:
    selection=self.select_current_tables(plan,attachments);si=self.si_records(selection,attachments)
   except Exception as e:self.warnings.append('部分 SI 附件解析失败：'+str(e))
  if not attachments:self.warnings.append('当前公开附件接口未找到 SI；继续使用正文证据。')
  # The paper's quantified identifiers define scope; SI supplies structures for them.
  si_resolved=set()
  for cid in targets:
   candidates={r['si_structure']['canonical'] for r in si if r['compound_id']==cid and r['si_structure']['valid']}
   if len(candidates)==1:si_resolved.add(cid)
  self.target_ids=targets-si_resolved
  names=[]
  if self.target_ids:
   self.stage='names';self.status('第二步：为 SI 未解决的编号提取完整名称',targets=len(self.target_ids))
   signature=hashlib.sha256(json.dumps(sorted(self.target_ids)).encode()).hexdigest()[:12]
   raw=self.cached('names_extracted_v3_'+signature+'.json',lambda:self.names(plan))
   self.stage='opsin';names=self.cached('names_opsin_base_v3_'+signature+'.json',lambda:self.resolve_names(raw))
   from name_retry import retry_missing,normalize_explicit_n_dimethyl
   names=self.cached('names_retry_v3_'+signature+'.json',lambda:retry_missing(self,plan,names))
   names=normalize_explicit_n_dimethyl(names,self.root)
   if any(r['si_structure']['valid'] for r in si):
    self.stage='reference_names';names=self.cached('names_resolved_v3_'+signature+'.json',lambda:recover_cited_names(self,plan,si,names,public_request))
   audit_names(si,names)
  save(self.root/'names_opsin.json',names)
  byid=defaultdict(list)
  for n in names:byid[n['compound_id']].append(n)
  cards=[]
  for cid in sorted(targets,key=lambda x:([(0,int(x))] if x.isdigit() else [(1,x)])):
   own=[r for r in si if r['compound_id']==cid and r['si_structure']['valid']]
   smis={r['si_structure']['canonical'] for r in own};smi=own[0]['si_smiles'] if len(smis)==1 else ''
   ns=byid[cid];valid=[n for n in ns if n['opsin']['valid'] and n.get('eligible_for_structure_comparison',True) and (n['name_text_supported'] or n.get('visual_name_review',{}).get('supported'))]
   status='SI 与完整名称结构一致' if smi and any(n['opsin']['canonical']==own[0]['si_structure']['canonical'] for n in valid) else ('SI 与名称需核对' if smi and valid else ('仅 SI 结构' if smi else ('仅名称结构，缺少 SI 对照' if valid else '结构待补充')))
   cards.append({'compound_id':cid,'si_smiles':smi,'si_source':own[0]['source_filename'] if own else '', 'si_evidence':own,'names':ns,'status':status,'has_structure':bool(smi or valid)})
  from structure_selection import select_structure
  self.target_ids={c['compound_id'] for c in cards if not select_structure(c)['smiles']}
  self.stage='markush';self.status('第三步：通式与取代基识别、连接及结构候选生成')
  from markush_workflow import run as reconstruct,attach
  reconstruction=reconstruct(self,plan,body) if self.target_ids else []
  cards=attach(cards,reconstruction)
  self.target_ids={c['compound_id'] for c in cards if not select_structure(c)['smiles']}
  self.stage='full_structure';self.status('第四步：完整结构截图与 SAR 三模型识别')
  from full_structure_workflow import run as recognize_full,attach as attach_full
  full_records=recognize_full(self,plan,body) if self.target_ids else []
  cards=attach_full(cards,full_records)
  self.target_ids=targets
  # SI measurements remain in si_records.json and a separate Excel sheet.
  save(self.root/'body_activity.json',body);save(self.root/'compound_records.json',cards)
  summary={'doi':plan['doi'],'title':plan['title'],'compounds':len(cards),'with_structure':sum(c['has_structure'] for c in cards),'with_si_structure':sum(bool(c['si_smiles']) for c in cards),'with_full_structure_accepted':sum(bool(c.get('full_structure_smiles')) for c in cards),'full_structure_policy':'majority_accept_v1','with_full_structure_candidates':sum(bool(c.get('full_structure')) for c in cards),'with_full_structure_consensus':sum(any(r['consensus'].get('candidate_smiles') for r in c.get('full_structure',[])) for c in cards),'with_markush_candidates':sum(bool(c.get('markush')) for c in cards),'with_markush_reviewed':sum(any(r.get('review_status')=='source_visually_reviewed' for r in c.get('markush',[])) for c in cards),'structure_route_order':['si','name','markush','full_structure'],'with_names':sum(bool(c['names']) for c in cards),'measurements':len(body),'sources':dict(Counter(r['source_type'] for r in body)),'warnings':self.warnings,'status':'completed' if cards else 'no_quantified_activity','scope_note':'只保留有数值或测定限值的活性编号。SI 与正文数据分源保存；表格数值直接解析 PaddleOCR 单元格。无表格页的视觉提取须复核；结构来源一致率不等于独立准确率。'}
  save(self.root/'summary.json',summary)
  from batch_report import build
  build(self.root);self.stage='completed';self.status('报告生成完成',**summary)
