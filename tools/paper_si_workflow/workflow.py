"""PDF-only research input; Luna plans evidence extraction, OPSIN resolves names,
and verified publisher SI supplies auditable structure/activity records.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import getpass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import time
import traceback

from core import (save,sha256,extract_pdf,numbered_page,source_span,name_fold,format_name,
                  identifier_present,normalize_id,mol_info,opsin_batch,read_tables,parse_number,
                  potency_nm,assay_family,safe_filename,fitz)
from model import Luna
import requests

VERSION='pdf-opsin-si-v2'

def public_request(method,url,**kwargs):
    if not url.startswith('https://'):raise ValueError('Public attachment URLs must use HTTPS')
    for attempt in range(3):
        try:
            with requests.Session() as session:
                session.trust_env=False
                response=session.request(method,url,timeout=(15,180),**kwargs)
            response.raise_for_status();return response
        except requests.RequestException:
            if attempt==2:raise
            time.sleep(attempt+1)

class Workflow:
    def __init__(self,pdf,root,llm):
        self.pdf=Path(pdf).resolve();self.root=Path(root).resolve();self.llm=llm;self.started=time.time()
        self.root.mkdir(parents=True,exist_ok=True)
        self.pages={};self.stage='input';self.warnings=[]
    def status(self,message,**extra):
        record={'stage':self.stage,'message':message,'elapsed_seconds':round(time.time()-self.started,1),**extra}
        save(self.root/'status.json',record);print(json.dumps(record,ensure_ascii=False),flush=True)
    def cached(self,name,fn):
        path=self.root/name
        if path.exists():return json.loads(path.read_text())
        data=fn();save(path,data);return data
    def image(self,page):
        path=self.root/'pages'/f'page-{page:03}.png'
        if not path.exists():
            path.parent.mkdir(parents=True,exist_ok=True)
            with fitz.open(self.pdf) as doc:doc[page-1].get_pixmap(matrix=fitz.Matrix(1.8,1.8)).save(path)
        return path
    def plan_document(self,document):
        previews=[]
        for page in document['pages']:
            matches=[x['text'] for x in page['lines'] if re.search(r'experimental|supporting|associated content|references|table \d|\(\d{1,3}[a-z]?\)',x['text'],re.I)]
            previews.append({'page':page['page'],'first_text':page['text'][:550],'section_or_heading_lines':matches[:18]})
        prompt='''Plan a literature chemical-data extraction workflow from this PDF only. Select the paper's own DOI, not a cited reference DOI. Identify pages containing explicit systematic compound names linked to compound IDs, including experimental chemistry, and pages containing bioactivity tables/figures. Exclude bibliography-only pages. Include names on earlier introduction pages if explicit. All page numbers are 1-based PDF pages.
Return {"title":"...","doi":"one supplied DOI candidate","targets":["..."],"compound_name_pages":[1],"main_activity_pages":[1],"si_search_doi":"same DOI","evidence":"short source-based rationale"}. Do not provide any structures or activity values.
'''+json.dumps({'doi_candidates':document['doi_candidates'],'first_two_pages':[numbered_page(x) for x in document['pages'][:2]],'page_previews':previews},ensure_ascii=False)
        plan=self.llm.ask('document_plan',prompt)
        candidates={x['doi'].lower() for x in document['doi_candidates']}
        if str(plan.get('doi','')).lower() not in candidates:raise ValueError('Model-selected DOI is not present in PDF evidence')
        plan['doi']=plan['doi'].lower();plan['si_search_doi']=plan['doi']
        for key in ['compound_name_pages','main_activity_pages']:
            values=sorted(set(int(x) for x in plan.get(key,[])))
            if any(x not in self.pages for x in values):raise ValueError('Model selected an invalid PDF page')
            plan[key]=values
        return plan
    def names(self,plan):
        selected=set(plan['compound_name_pages'])
        in_experimental=False
        for page in self.pages.values():
            if not plan.get('targeted_name_search') and re.search(r'EXPERIMENTAL SECTION',page['text'],re.I):in_experimental=True
            if re.search(r'REFERENCES',page['text']):in_experimental=False
            if in_experimental:selected.add(page['page'])
        selected=sorted(selected);chunks=[selected[i:i+2] for i in range(0,len(selected),2)]
        def extract(index,batch):
            context=set(batch)
            for n in batch:
                if n>1:context.add(n-1)
                if n+1 in self.pages:context.add(n+1)
            texts=[]
            for n in sorted(context):
                page=self.pages[n]
                if n in batch:texts.append(numbered_page(page))
                elif n<min(batch):texts.append('\n'.join(numbered_page(page).splitlines()[-18:]))
                else:texts.append('\n'.join(numbered_page(page).splitlines()[:18]))
            prompt='''Extract ALL explicit systematic chemical names with their compound identifiers from the designated core pages. Only extract the allowed target IDs supplied below. Prioritize Experimental Section full systematic headings over introduction/discussion shorthand. A name must include every printed continuation across lines, columns, and pages; do not stop at a line break. Label roles. Do not extract unnamed reagents, synthesis solvents, bibliography titles, shorthand nicknames without a full chemical name, or fabricate names from structures. Preserve every locant, bracket, stereodescriptor and salt name. Repair only line wrapping/typographical spacing, not chemical content. Full names may span two columns or consecutive pages. Adjacent page context is supplied; include a record only if part of its name/ID evidence lies on a core page.
Return {"compounds":[{"compound_id":"84","name":"full systematic name, without the trailing compound ID","role":"target|intermediate|reference|unknown","source_ranges":[{"page":1,"start_line":1,"end_line":4}],"note":"any source ambiguity"}]}.
Each source range MUST use the supplied Pn:Lm line numbers and include the actual name and printed ID, normally 1–10 lines. Do not cite entire pages. If the same ID has multiple genuinely different printed names, preserve each. Return an empty list if none.
'''+json.dumps({'allowed_target_ids':sorted(getattr(self,'target_ids',[])),'core_pages':batch,'document_target':plan.get('targets',[]),'source':'\n\n'.join(texts)},ensure_ascii=False)
            output=self.llm.ask(f'names_{index:02}',prompt)
            results=[];rejected=[]
            for record in output.get('compounds',[]):
                try:
                    raw,locations=source_span(record.get('source_ranges'),self.pages)
                    if getattr(self,'target_ids',None) and normalize_id(record['compound_id']) not in self.target_ids:continue
                    if not any(p['page'] in batch for p in locations):raise ValueError('No core-page evidence')
                    if not identifier_present(record['compound_id'],raw):raise ValueError('Compound ID absent from source span')
                    name=format_name(record['name'])
                    if not name or len(name)>1600:raise ValueError('Invalid extracted name')
                    results.append({**record,'compound_id':normalize_id(record['compound_id']),'name':name,
                                    'source_ranges':locations,'source_text':raw,'name_text_supported':name_fold(name) in name_fold(raw)})
                except Exception as exc:rejected.append({'record':record,'reason':str(exc)})
            save(self.root/'raw'/f'names_batch_{index:02}.json',{'records':results,'rejected':rejected})
            return results
        results=[]
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures={pool.submit(extract,i,batch):batch for i,batch in enumerate(chunks)}
            for completed,future in enumerate(as_completed(futures),1):
                results+=future.result();self.status('正文名称批次完成',completed_batches=completed,total_batches=len(chunks),name_candidates=len(results))
        unique={}
        for record in results:
            key=(record['compound_id'],name_fold(record['name']))
            if key not in unique:unique[key]=record
            else:
                for span in record['source_ranges']:
                    if span not in unique[key]['source_ranges']:unique[key]['source_ranges'].append(span)
        return sorted(unique.values(),key=lambda r:(r['source_ranges'][0]['page'],r['source_ranges'][0]['start_line'],r['compound_id']))
    def resolve_names(self,names):
        parsed=opsin_batch([r['name'] for r in names],self.root/'opsin/initial')
        for record,result in zip(names,parsed):record['opsin']=result
        failed=[(i,r) for i,r in enumerate(names) if not r['opsin']['valid'] or not r['name_text_supported']]
        self.status('OPSIN 首次解析完成',names=len(names),valid=sum(r['opsin']['valid'] for r in names),source_or_parse_review=len(failed))
        # Correct source transcription only; SI has not been retrieved at this stage.
        for start in range(0,len(failed),5):
            group=failed[start:start+5]
            images=[self.image(p) for p in sorted({s['page'] for _,r in group for s in r['source_ranges']})]
            prompt='''Verify these chemical-name transcriptions against the original PDF page images and quoted text. OPSIN rejected the name or text-span matching failed. You may fix line-break hyphenation, Unicode ligatures, or transcription mistakes ONLY where the printed source supports it. Preserve stereochemistry and salt form. Do not construct names from chemical knowledge and do not omit unsupported stereochemistry to make OPSIN accept a name.
Return {"reviews":[{"index":0,"supported":true,"name":"complete corrected or unchanged source name","reason":"specific source evidence"}]}. If uncertain, supported=false and leave the name unchanged. No SMILES.
'''+json.dumps({'records':[{'index':i,'compound_id':r['compound_id'],'name':r['name'],'source_text':r['source_text'],'pages':[s['page'] for s in r['source_ranges']],'opsin_status':r['opsin']['status']} for i,r in group]},ensure_ascii=False)
            review=self.llm.ask(f'name_review_{start:03}',prompt,images)
            allowed={i for i,_ in group}
            for item in review.get('reviews',[]):
                index=item.get('index')
                if index not in allowed:continue
                record=names[index];record['visual_name_review']=item
                if item.get('supported') is True and item.get('name'):
                    revised=format_name(item['name']);record['original_name']=record['name'];record['name']=revised
                    record['opsin_initial']=record['opsin']
                    record['opsin']=opsin_batch([revised],self.root/'opsin'/f'review_{index:03}')[0]
            self.status('正文名称源图复核完成',reviewed=min(start+5,len(failed)),total_reviews=len(failed))
        return names
    def discover_si(self,plan):
        doi=plan['doi'];searches=[];found={}
        for query in [{'resource_doi':doi,'page_size':100},{'search_for':doi,'page_size':100}]:
            try:
                rows=public_request('POST','https://api.figshare.com/v2/articles/search',json=query).json()
                searches.append({'query':query,'returned_ids':[r.get('id') for r in rows]})
                for row in rows:
                    if str(row.get('resource_doi','')).lower()==doi:found[row['id']]=row
            except Exception as exc:searches.append({'query':query,'error':str(exc)})
            if found:break
        save(self.root/'si/search.json',{'doi':doi,'searches':searches,'matched_articles':list(found.values())})
        attachments=[];downloads=[]
        # Collect all metadata first: a slow PDF must not block later CSV articles.
        for article_id in sorted(found):
            metadata_url=f'https://api.figshare.com/v2/articles/{article_id}'
            try:
                meta=public_request('GET',metadata_url).json()
                save(self.root/'si'/f'article_{article_id}.json',meta)
                if str(meta.get('resource_doi','')).lower()!=doi:continue
                if meta.get('download_disabled') or meta.get('is_embargoed'):
                    self.warnings.append(f'SI article {article_id} is not openly downloadable');continue
                downloads.extend((article_id,metadata_url,meta,f) for f in meta.get('files',[]))
            except Exception as exc:
                self.warnings.append(f'SI metadata {article_id}: {exc}')
        downloads.sort(key=lambda item:(Path(item[3]['name']).suffix.lower() not in ('.csv','.tsv','.xlsx'),int(item[3].get('size',0))))
        for article_id,metadata_url,meta,file in downloads:
            try:
                name=safe_filename(file['name']);path=self.root/'si'/f"{file['id']}_{name}"
                expected=file.get('computed_md5') or file.get('supplied_md5')
                def verified():
                    return path.exists() and path.stat().st_size==file.get('size') and (not expected or hashlib.md5(path.read_bytes()).hexdigest()==expected)
                reused=verified()
                cache_source=None
                # A sibling SI archive is a cache, never a manually specified research input.
                # File identity must match fresh DOI-bound publisher metadata.
                cached=self.pdf.parent/'si'/doi.split('/',1)[1]/name
                shared=self.pdf.parent.parent/'si'/doi.split('/',1)[1]/name
                if not cached.exists() and shared.exists():cached=shared
                if not reused and expected and cached.is_file() and cached.stat().st_size==file.get('size') and hashlib.md5(cached.read_bytes()).hexdigest()==expected:
                    shutil.copy2(cached,path);reused=True;cache_source=str(cached)
                if not reused and path.suffix.lower()=='.pdf' and any(
                    any(str(cell).strip().lower() in ('smile','smiles','canonical smiles') for preview in a.get('table_previews',[]) for row in preview.get('preview_rows',[])[:3] for cell in row)
                    for a in attachments):
                    skipped=self.root/'si_optional_downloads.json'
                    pending=json.loads(skipped.read_text()) if skipped.exists() else []
                    item={'filename':name,'url':file['download_url'],'size':file.get('size'),'reason':'Verified structured SI with SMILES is available; large PDF deferred until needed.'}
                    if item not in pending:pending.append(item)
                    save(skipped,pending)
                    self.warnings.append('已获得作者 SMILES 表；附加 SI PDF 按需下载：'+name)
                    continue
                if not reused:
                    if int(file.get('size',0))>150*1024*1024:
                        self.warnings.append(f'SI attachment exceeds prototype size limit: {name}');continue
                    temp=path.with_suffix(path.suffix+'.part')
                    with requests.Session() as session:
                        session.trust_env=False
                        with session.get(file['download_url'],stream=True,timeout=(15,45)) as response:
                            response.raise_for_status();began=time.time();last_update=began;received=0
                            with temp.open('wb') as stream:
                                for chunk in response.iter_content(65536):
                                    stream.write(chunk);received+=len(chunk)
                                    if time.time()-began>240:raise TimeoutError('SI download exceeded four minutes; partial file retained')
                                    if time.time()-last_update>15:
                                        self.status('正在下载官方 SI 附件',filename=name,downloaded_bytes=received,total_bytes=file.get('size'))
                                        last_update=time.time()
                    temp.replace(path)
                if not verified():raise ValueError('Attachment checksum/size mismatch: '+name)
                item={'file_id':str(file['id']),'article_id':article_id,'filename':name,'local_path':str(path),
                      'article_doi':doi,'supplement_doi':meta.get('doi'),'source_url':file['download_url'],
                      'landing_url':meta.get('url_public_html'),'metadata_url':metadata_url,'sha256':sha256(path),
                      'md5':expected,'size':path.stat().st_size,'download_verified':True,'cache_reused':reused,'cache_source':cache_source}
                tables=read_tables(path)
                item['text_encoding']=tables[0].get('encoding') if tables else None
                item['table_previews']=[{'sheet':t['sheet'],'row_count':len(t['rows']),'preview_rows':t['rows'][:8]} for t in tables]
                if path.suffix.lower()=='.pdf':
                    with fitz.open(path) as doc:
                        item['pdf_pages']=len(doc);item['pdf_preview']=[{'page':n+1,'text':doc[n].get_text()[:2200]} for n in range(min(2,len(doc)))]
                attachments.append(item);save(self.root/'si_manifest.partial.json',attachments);self.status('已定位并校验 SI 附件',filename=name,files=len(attachments))
            except Exception as exc:
                self.warnings.append(f"SI attachment {file.get('name')}: {exc}")
                continue
        if not attachments:self.warnings.append('No downloadable DOI-matched SI found through the current Figshare provider. Other publisher adapters are not yet implemented.')
        return attachments
    def select_tables(self,plan,attachments):
        return self.llm.ask('si_table_selection','''Locate compound-to-structure-to-bioactivity tables in the discovered official SI attachments. Prefer directly machine-readable tables over spectra PDFs. Select actual table schemas, not values. Header row and data start row are ZERO-BASED; source CSV row numbers will later be retained. Column indexes are ZERO-BASED. Multiple header lines are possible. Only choose fields supported by supplied headers and article context.
Return {"tables":[{"file_id":"...","sheet":"","header_row":1,"data_start_row":2,"compound_id_column":0,"smiles_column":1,"name_column":null,"measurements":[{"column":2,"endpoint":"EC50","unit":"nM","target":"...","assay":"...","header_evidence":"exact heading text"}],"reason":"..."}],"pdf_tables":[{"file_id":"...","pages":[1],"reason":"..."}],"ignored_files":[{"file_id":"...","reason":"..."}]}.
Do not treat NMR/HPLC characterization spectra as potency assays. Use an empty list when unsupported. For PDF tables choose pages only when the provided preview gives actual table evidence; do not guess page numbers.
'''+json.dumps({'paper':plan,'attachments':attachments},ensure_ascii=False))
    def si_records(self,selection,attachments):
        files={a['file_id']:a for a in attachments};records=[]
        for schema in selection.get('tables',[]):
            file=files.get(str(schema.get('file_id')))
            if not file:raise ValueError('Model selected an unknown SI file')
            table=next((t for t in read_tables(file['local_path']) if t['sheet']==schema.get('sheet','')),None)
            if table is None:raise ValueError('Unknown SI worksheet')
            rows=table['rows'];header=int(schema['header_row']);start=int(schema['data_start_row'])
            if not 0<=header<start<=len(rows):raise ValueError('Invalid SI header/data row indexes')
            columns=[schema.get('compound_id_column'),schema.get('smiles_column'),schema.get('name_column')]+[m.get('column') for m in schema.get('measurements',[])]
            if any(not isinstance(c,int) or c<0 or c>=len(rows[header]) for c in columns if c is not None):raise ValueError('Invalid SI column mapping')
            def cell(row,col):return str(row[col]).strip() if col is not None and col<len(row) else ''
            for number,row in enumerate(rows[start:],start+1):
                compound=cell(row,schema['compound_id_column'])
                if not compound:continue
                smiles=cell(row,schema.get('smiles_column'));name=cell(row,schema.get('name_column'))
                count_before=len(records)
                for measurement in schema.get('measurements',[]):
                    value=cell(row,measurement['column'])
                    if not value:continue
                    records.append({'record_id':f"{file['file_id']}:{schema.get('sheet','')}:{number}:{measurement['column']}",
                        'compound_id':normalize_id(compound),'raw_compound_id':compound,'si_smiles':smiles,'si_name':name,
                        'si_structure':mol_info(smiles),'target':str(measurement.get('target','')),'assay':str(measurement.get('assay','')),
                        'endpoint':str(measurement.get('endpoint','')),'value':value,'unit':str(measurement.get('unit','')),
                        'source_type':'author_si_table_direct','file_id':file['file_id'],'source_filename':file['filename'],
                        'source_sheet':schema.get('sheet',''),'source_row':number,'source_column':measurement['column']+1,
                        'source_header':rows[header],'source_row_values':row,'source_url':file['source_url']})
                if len(records)==count_before and smiles:
                    records.append({'record_id':f"{file['file_id']}:{schema.get('sheet','')}:{number}:structure",
                        'compound_id':normalize_id(compound),'raw_compound_id':compound,'si_smiles':smiles,'si_name':name,
                        'si_structure':mol_info(smiles),'target':'','assay':'','endpoint':'','value':'','unit':'',
                        'source_type':'author_si_structure_table','file_id':file['file_id'],'source_filename':file['filename'],
                        'source_sheet':schema.get('sheet',''),'source_row':number,'source_header':rows[header],
                        'source_row_values':row,'source_url':file['source_url']})
        if selection.get('pdf_tables'):self.warnings.append('PDF-only SI activity tables require the PDF table adapter; this trial uses selected structured attachments.')
        return records
    def body_activity(self,plan):
        selected=plan['main_activity_pages']
        def extract(page):
            prompt='''Read the actual bioactivity measurements in this PDF page, especially compound activity tables and assay figures. Do not use external information or infer unprinted values. Keep distinct assays separate (e.g. cAMP versus tissue GTPgammaS). Extract EC50, pEC50, IC50, pIC50, Ki, Kd and Emax if present; do not treat PK exposure, cLogP or solubility as potency. Preserve inequalities, uncertainty, units and N.D. entries. If EC50 is not printed but pEC50 is, record only pEC50; do not calculate an EC50 here.
Compound identifiers must exclude a superscript footnote marker ONLY when the page footnote establishes it as an annotation. Preserve genuine letter suffixes denoting different compounds. Record the raw visible identifier and marker separately.
Return {"tables":[{"caption":"Table/Figure ...","source_role":"current_study|cited_prior_study|uncertain","citation_label":"printed reference number if any","target":"...","assay":"specific experimental assay","rows":[{"compound_id":"...","raw_identifier":"...","footnote_marker":"...","measurements":[{"endpoint":"EC50","value":"...","unit":"nM"}]}]}],"notes":"any uncertainty"}. No structure data. Only this page's results, not those of prior cited papers unless explicitly labelled as such.
'''+json.dumps({'page':page,'paper_targets':plan.get('targets',[]),'text':self.pages[page]['text']},ensure_ascii=False)
            data=self.llm.ask(f'body_activity_{page:03}',prompt,[self.image(page)])
            rows=[]
            for table in data.get('tables',[]):
                for row in table.get('rows',[]):
                    for m in row.get('measurements',[]):
                        rows.append({'compound_id':normalize_id(row.get('compound_id')),'raw_identifier':row.get('raw_identifier'),
                                     'footnote_marker':row.get('footnote_marker'),'target':table.get('target',''),
                                     'assay':table.get('assay',''),'caption':table.get('caption',''),'page':page,
                                     'endpoint':m.get('endpoint',''),'value':str(m.get('value','')),'unit':str(m.get('unit','')),
                                     'source_type':'main_pdf_luna_vision','source_role':table.get('source_role','uncertain'),'citation_label':table.get('citation_label','')})
            return rows
        result=[]
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures={pool.submit(extract,page):page for page in selected}
            for completed,future in enumerate(as_completed(futures),1):
                result+=future.result();self.status('正文活性核对页完成',completed_pages=completed,total_pages=len(selected))
        return sorted(result,key=lambda r:(r['page'],r['compound_id'],r['endpoint']))
    def merge(self,records,names,body):
        by_id={}
        for i,name in enumerate(names):
            if name.get('eligible_for_structure_comparison',True):by_id.setdefault(name['compound_id'],[]).append(i)
        output=[]
        for source in records:
            r=dict(source);indices=by_id.get(r['compound_id'],[]);r['name_indices']=indices
            supported=[names[i] for i in indices if names[i]['opsin']['valid'] and (names[i]['name_text_supported'] or names[i].get('visual_name_review',{}).get('supported') is True)]
            valid=[x['opsin'] for x in supported];si=r['si_structure']
            candidates={x['canonical'] for x in valid}
            if si['valid'] and any(x['canonical']==si['canonical'] for x in valid):match='exact'
            elif si['valid'] and any(x['parent']==si['parent'] for x in valid):match='parent_only'
            elif si['valid'] and any(x['connectivity']==si['connectivity'] for x in valid):match='stereo_difference'
            elif valid and si['valid']:match='different'
            elif indices:match='opsin_unresolved_or_evidence_review'
            else:match='name_not_found'
            r['structure_comparison']=match;r['multiple_name_structures']=len(candidates)>1
            if si['valid']:r['selected_smiles']=r['si_smiles'];r['selected_structure_source']='author_si'
            elif len(candidates)==1:r['selected_smiles']=next(iter(candidates));r['selected_structure_source']='main_name_opsin'
            else:r['selected_smiles']='';r['selected_structure_source']='unresolved'
            si_potency=potency_nm(r['value'],r['endpoint'],r['unit']);checks=[]
            for other in body:
                if other['compound_id']!=r['compound_id'] or assay_family(other['assay'])!=assay_family(r['assay']):continue
                if re.sub(r'\W','',other['target']).lower()!=re.sub(r'\W','',r['target']).lower():continue
                value=potency_nm(other['value'],other['endpoint'],other['unit'])
                if not value or not si_potency:continue
                tolerance=.025 if value['derived'] else 1e-9
                same=value['relation']==si_potency['relation'] and math.isclose(value['value'],si_potency['value'],rel_tol=tolerance,abs_tol=1e-8)
                checks.append({'page':other['page'],'caption':other['caption'],'endpoint':other['endpoint'],'value':other['value'],
                    'unit':other['unit'],'normalized':value,'agreement':same,'comparison':'rounded_pEC50_conversion' if value['derived'] else 'explicit_EC50'})
            r['body_activity_checks']=checks
            direct=[c for c in checks if c['comparison']=='explicit_EC50']
            r['activity_comparison']='conflict' if any(not c['agreement'] for c in direct) else ('supported' if any(c['agreement'] for c in checks) else ('conflict' if checks else 'not_found'))
            r['potency_normalized']=si_potency
            output.append(r)
        return output
    def review_outcome(self,plan,merged,names):
        conflicts=[]
        for r in merged:
            if r['structure_comparison'] in ['exact','parent_only'] and r['activity_comparison']!='conflict' and not r['multiple_name_structures']:
                continue
            conflicts.append({'record_id':r['record_id'],'compound_id':r['compound_id'],
                'structure_comparison':r['structure_comparison'],'activity_comparison':r['activity_comparison'],
                'si_smiles_original':r['si_smiles'],'si_value':r['value'],'si_unit':r['unit'],
                'body_checks':r['body_activity_checks'],
                'body_names':[{'name':names[i]['name'],'source_text':names[i]['source_text'],
                    'source_ranges':names[i]['source_ranges'],'opsin':names[i]['opsin'],
                    'visual_name_review':names[i].get('visual_name_review')} for i in r['name_indices']]})
        return self.llm.ask('workflow_review_v2','''Review source-specific differences in this PDF -> OPSIN -> official SI extraction. Return {"summary":"short factual Chinese summary","issues":[{"record_id":"an input record id","explanation":"Chinese explanation tied to the supplied evidence","recommendation":"review action, never an invented correction"}]}.
Do not change structures, IDs, or numerical values. Source names, OPSIN SMILES and original SI SMILES are supplied for interpretation, not for repairing the independently extracted names. A source-supported shorthand name can parse successfully yet omit parts of the actual compound. Distinguish this from an OPSIN parsing failure. Strict canonical-SMILES mismatch may reflect tautomer representations; do not call it a different scaffold without evidence. A missing stereodescriptor is different from an explicitly opposite stereodescriptor. Preserve salt forms and missing-data uncertainty. Mention absent name evidence and malformed printed names where relevant. Author SI is an extraction input, so agreement is not independent extraction accuracy.
'''+json.dumps({'paper':plan,'si_records':len(merged),'name_records':len(names),'conflicts':conflicts},ensure_ascii=False))
    def run(self):
        digest=sha256(self.pdf);metadata_path=self.root/'metadata.json'
        if metadata_path.exists() and json.loads(metadata_path.read_text())['pdf_sha256']!=digest:raise ValueError('Output folder belongs to a different PDF')
        save(metadata_path,{'workflow_version':VERSION,'input_pdf':str(self.pdf),'pdf_sha256':digest,
            'research_inputs':['PDF only'],'model':self.llm.model,'api_base':self.llm.base,
            'started_at_utc':datetime.now(timezone.utc).isoformat(),'no_previous_reference_outputs_used':True,
            'provider_scope':'DOI-matched public Figshare/ACS SI; structured CSV/TSV/XLSX tables supported in this prototype'})
        (self.root/'input').mkdir(exist_ok=True);shutil.copy2(self.pdf,self.root/'input/source.pdf')
        document=self.cached('pdf_text.json',lambda:extract_pdf(self.pdf))
        self.pages={p['page']:p for p in document['pages']}
        if not document['doi_candidates']:raise ValueError('No DOI found in PDF text; scanned PDFs/title search need an additional adapter')
        self.stage='plan';self.status('Luna 根据 PDF 识别 DOI 和处理范围')
        plan=self.cached('plan.json',lambda:self.plan_document(document));self.status('文献计划已生成',doi=plan['doi'],name_pages=plan['compound_name_pages'],activity_pages=plan['main_activity_pages'])
        self.stage='si_discovery';self.status('按 PDF DOI 自动检索官方 SI')
        attachments=self.cached('si_manifest.json',lambda:self.discover_si(plan))
        if not attachments:raise RuntimeError('No supported SI attachments found; partial name results retained')
        self.stage='si_tables';selection=self.cached('si_selection.json',lambda:self.select_tables(plan,attachments))
        records=self.cached('si_records.json',lambda:self.si_records(selection,attachments))
        if not records:raise RuntimeError('No structured SI activity rows selected; partial outputs retained')
        self.status('已从 SI 原始表格读取记录',records=len(records),valid_structures=sum(r['si_structure']['valid'] for r in records))
        self.target_ids={r['compound_id'] for r in records if parse_number(r['value']) is not None}
        self.stage='names';names=self.cached('names_extracted_v2.json',lambda:self.names(plan))
        self.stage='opsin';names=self.cached('names_opsin_base_v2.json',lambda:self.resolve_names(names))
        from name_resolution import recover_cited_names
        self.stage='name_completeness';self.status('检查完整名称，追溯有原文依据的参考文献')
        names=self.cached('names_resolved_v2.json',lambda:recover_cited_names(self,plan,records,names,public_request))
        save(self.root/'names_opsin.json',names)
        self.stage='body_activity';self.status('独立读取正文活性表，核对 SI 数值与实验')
        body=self.cached('body_activity.json',lambda:self.body_activity(plan))
        self.stage='merge';merged=self.merge(records,names,body);save(self.root/'merged_records.json',merged)
        review=self.review_outcome(plan,merged,names);save(self.root/'workflow_review.json',review)
        summary={'pdf_pages':len(self.pages),'doi':plan['doi'],'names':len(names),'opsin_valid':sum(r['opsin']['valid'] for r in names),
            'si_files':len(attachments),'si_records':len(records),'si_compounds':len({r['compound_id'] for r in records}),
            'valid_si_structures':sum(r['si_structure']['valid'] for r in records),'structure_comparison':dict(Counter(r['structure_comparison'] for r in merged)),
            'activity_comparison':dict(Counter(r['activity_comparison'] for r in merged)),
            'complete_structure_activity_records':sum(bool(r['selected_smiles']) and r['potency_normalized'] is not None for r in merged),
            'warnings':self.warnings,'scope_note':'SI is an extraction input, not an independent reference benchmark. Names are transcribed from source text without SI-driven rewriting. SI is used to detect partial names and verify cited-paper ID mappings; raw names and sources are retained.'}
        save(self.root/'summary.json',summary)
        self.stage='report';self.status('生成 SAR 风格报告',**summary)
        from priority_report import build
        build(self.root)
        self.stage='completed';self.status('PDF → OPSIN → SI → 报告流程完成',**summary)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--legacy-si-first',action='store_true',help='Use the previous SI-first workflow for comparison')
    parser.add_argument('--pdf',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--api-base',default='https://newapi.livzon.cn/v1')
    parser.add_argument('--model',default='GPT-6-Luna')
    args=parser.parse_args();args.output=args.output.resolve();args.output.mkdir(parents=True,exist_ok=True)
    key=os.environ.get('PAPER_SI_API_KEY') or getpass.getpass('API key (hidden; memory only): ')
    log=(args.output/'run.log').open('a',buffering=1)
    class Redacted:
        def write(self,text):return log.write(str(text).replace(key,'[REDACTED]'))
        def flush(self):log.flush()
        def isatty(self):return False
    sys.stdout=Redacted();sys.stderr=sys.stdout
    from workflow_v3 import ActivityWorkflow
    flow=(Workflow if args.legacy_si_first else ActivityWorkflow)(args.pdf,args.output,Luna(key,args.api_base,args.model,args.output))
    try:flow.run()
    except BaseException as exc:
        flow.status('流程未完成，已保留检查点',failed=True,error=str(exc).replace(key,'[REDACTED]'))
        traceback.print_exc();return 1
    return 0

if __name__=='__main__':raise SystemExit(main())
