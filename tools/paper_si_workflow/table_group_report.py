"""Source table groups rendered by the SAR PDF crop tool."""
import json,re,sys,hashlib,shutil,os
from pathlib import Path
from core import fitz,save

def prepare(root,body,cards,assets):
 root=Path(root);index={c['compound_id']:c for c in cards};groups=[];used=set()
 crops_path=root/'activity_ocr/table_crops/inputs.json'
 if not crops_path.exists():return [],[(None,c) for c in cards]
 crop_index={r['id']:r for r in json.loads(crops_path.read_text())}
 tids={r['table_id'] for c in cards for r in c.get('markush',[])+c.get('full_structure',[])}
 for path in (root/'full_structure').glob('*/layout.json'):
  if json.loads(path.read_text()).get('table_type')=='scaffold_fragments':tids.add(path.parent.name)
 with fitz.open(root/'input/source.pdf') as pdf:
  for tid in sorted(tids):
   if tid not in crop_index:continue
   evidence=[r for r in body if r.get('table_id')==tid]
   if not evidence:continue
   entry=crop_index[tid];page=pdf[entry['page']-1];box=list(entry['bbox_pdf_points'])
   caption=evidence[0].get('caption',tid);label=re.match(r'Table\s+\d+',caption,re.I)
   if not label:
    layout_path=root/'full_structure'/tid/'layout.json'
    if layout_path.exists():label=re.match(r'Table\s+\d+',json.loads(layout_path.read_text()).get('table_label',''),re.I)
   hits=page.search_for(label[0]) if label else []
   hits=[h for h in hits if box[0]-15<=h.x0<=box[2] and h.y0<=box[1] and box[1]-h.y0<230]
   reviewed_context=next((r.get('table_context_bbox_pdf') for c in cards for r in c.get('markush',[]) if r.get('table_id')==tid and r.get('review_status')=='source_visually_reviewed' and r.get('table_context_bbox_pdf')),None)
   if reviewed_context:box=list(reviewed_context)
   elif hits:
    heading=max(hits,key=lambda h:h.y0);box[1]=max(0,heading.y0-4)
   else:continue # Never show a cropped data grid as the complete original table.
   # Include adjacent footnotes, but avoid the next scheme or paragraph.
   if not reviewed_context:box[3]=min(page.rect.height,box[3]+34)
   ids=list(dict.fromkeys(r['compound_id'] for r in sorted(evidence,key=lambda r:r.get('source_row',0))))
   groups.append({'table_id':tid,'label':label[0] if label else tid,'page':entry['page'],'caption':caption,'ids':[cid for cid in ids if cid in index],'crop':tid+'.png','bbox_pixels':box,'bbox_page_size':[page.rect.width,page.rect.height]})
 if groups:
  from full_structure_workflow import SAR
  sys.path.insert(0,str(SAR))
  from pipeline.smiles_validation.ocsr_render import render_inputs
  folder=root/'report_table_sources';folder.mkdir(exist_ok=True)
  mapping,audit=render_inputs(root/'input/source.pdf',groups,folder,folder/'final',dpi=300)
  for g,a in zip(groups,audit):
   source=Path(mapping[str((folder/g['crop']).resolve())]);name='source-table-'+hashlib.sha256(source.read_bytes()).hexdigest()[:20]+'.png';shutil.copy2(source,assets/name);g['image']='assets/'+name;g['render']=a
  save(folder/'manifest.json',groups)
 contexts=[]
 for g in groups:
  for cid in g['ids']:contexts.append((g['table_id'],index[cid]));used.add(cid)
 contexts.extend((None,c) for c in cards if c['compound_id'] not in used)
 return groups,contexts

def compound_source(root,card,record,groups,body,assets):
 """Use the same table/row as the displayed reconstruction, never a different source."""
 from source_overlay import marked_fragment_table
 root=Path(root);lookup={g['table_id']:g for g in groups};cid=card['compound_id']
 tid=record.get('table_id') if record else next((r.get('table_id') for r in body if r['compound_id']==cid and r.get('table_id') in lookup),None)
 if tid not in lookup:return {'images':[],'context_images':[]}
 g=lookup[tid];box=record.get('fragment_bbox_pdf_points') if record else None
 if record and record.get('graph_nodes'):
  from graph_source_report import graph_source
  return graph_source(root,record,g,assets)
 if record and record.get('fragments'):
  from full_structure_workflow import SAR
  sys.path.insert(0,str(SAR))
  from pipeline.smiles_validation.ocsr_render import render_inputs
  with fitz.open(root/'input/source.pdf') as doc:
   page=doc[g['page']-1];left,top,right,bottom=g['bbox_pixels'];words=page.get_text('words')
   labels=[w for w in words if w[4]==cid and left<=w[0]<=left+35 and top<w[1]<bottom]
   if len(labels)!=1:raise ValueError('Ambiguous text substituent row: '+cid)
   label=labels[0];mid=(label[1]+label[3])/2
   cells=sorted([w for w in words if w[0]>label[2]+10 and w[2]<right and w[1]<=mid<=w[3]],key=lambda w:w[0])[:3]
   if [w[4] for w in cells]!=list(record['substituents'].values()):raise ValueError('Printed substituent row mismatch: '+cid)
   box=[min(w[0] for w in cells)-3,min(w[1] for w in cells)-0.5,max(w[2] for w in cells)+3,max(w[3] for w in cells)+0.5];size=[page.rect.width,page.rect.height]
  folder=root/'report_table_sources'/('text-sites-'+cid);folder.mkdir(exist_ok=True)
  mapping,audit=render_inputs(root/'input/source.pdf',[{'crop':'row.png','page':g['page'],'bbox_pixels':box,'bbox_page_size':size}],folder,folder/'final',dpi=300)
  src=Path(next(iter(mapping.values())));name='source-fragment-'+hashlib.sha256(src.read_bytes()).hexdigest()[:20]+'.png';shutil.copy2(src,assets/name);save(folder/'manifest.json',audit)
  g={**g,'image':marked_fragment_table(root,g,record,box,assets)}
  return {'images':['assets/'+name,g['image']],'context_images':[g['image']]}

 if not box and record and record.get('origin')=='verified_local_trial':
  trial=Path(__file__).resolve().parents[2]/'runs/gpr88_table1_reconstruction'
  if json.loads((trial/'summary.json').read_text())['source_pdf_sha256']==record.get('source_pdf_sha256'):
   box=json.loads((trial/'crops.json').read_text())['boxes_pdf_points'].get('r1_'+cid)
 if not box:
  layout_path=root/'markush'/tid/'layout.json'
  if layout_path.exists():
   for row in json.loads(layout_path.read_text()).get('rows',[]):
    if row.get('compound_id')==cid:
     with fitz.open(root/'input/source.pdf') as doc:
      page=doc[g['page']-1];box=[v*(page.rect.width if i%2==0 else page.rect.height)/1000 for i,v in enumerate(row['fragment_bbox'])]
 fragment=''
 if box:
  from full_structure_workflow import SAR
  sys.path.insert(0,str(SAR))
  from pipeline.smiles_validation.ocsr_render import render_inputs
  folder=root/'report_table_sources'/('fragment-'+tid+'-'+re.sub(r'[^A-Za-z0-9_-]','_',cid));folder.mkdir(exist_ok=True)
  with fitz.open(root/'input/source.pdf') as doc:size=list(doc[g['page']-1].rect)[2:]
  from pipeline.activity_table_extraction.table_subject_crops import unique_detector_box
  from pipeline.utils.crop_geometry import padded_fraction
  from PIL import Image
  detection_file=root/'report_table_sources/sar_detections.json'
  page_image=str((root/'pages'/f"page-{g['page']:03}.png").resolve())
  detections=json.loads(detection_file.read_text()) if detection_file.exists() else {}
  if page_image not in detections:
   from patent_sar.model_adapters.molparser_iface import detect_page
   previous=os.environ.get('PATENT_SAR_VL_CACHE');os.environ['PATENT_SAR_VL_CACHE']=str(root/'sar_detector_cache')
   try:detections.update(detect_page([page_image],timeout=600))
   finally:
    if previous is None:os.environ.pop('PATENT_SAR_VL_CACHE',None)
    else:os.environ['PATENT_SAR_VL_CACHE']=previous
   save(detection_file,detections)
  with Image.open(page_image) as im:pixel_size=im.size
  sx,sy=pixel_size[0]/size[0],pixel_size[1]/size[1]
  source_box=[v*(sx if i%2==0 else sy) for i,v in enumerate(box)]
  selected=unique_detector_box(source_box,detections.get(page_image,[]))
  if selected is None:
   table_path=str((root/'report'/g['image']).resolve());table_file=root/'report_table_sources/sar_table_detections.json'
   table_detections=json.loads(table_file.read_text()) if table_file.exists() else {}
   if table_path not in table_detections:
    from patent_sar.model_adapters.molparser_iface import detect_page
    previous=os.environ.get('PATENT_SAR_VL_CACHE');os.environ['PATENT_SAR_VL_CACHE']=str(root/'sar_detector_cache')
    try:table_detections.update(detect_page([table_path],timeout=600))
    finally:
     if previous is None:os.environ.pop('PATENT_SAR_VL_CACHE',None)
     else:os.environ['PATENT_SAR_VL_CACHE']=previous
    save(table_file,table_detections)
   with Image.open(table_path) as im:tw,th=im.size
   tx,ty,tu,tv=g['render']['pdf_bbox'];converted=[]
   for item in table_detections.get(table_path,[]):
    a,b,c,d=item['bbox'];converted.append({**item,'bbox':[(tx+a/tw*(tu-tx))*sx,(ty+b/th*(tv-ty))*sy,(tx+c/tw*(tu-tx))*sx,(ty+d/th*(tv-ty))*sy]})
   selected=unique_detector_box(source_box,converted)
  if selected is None:
   save(folder/'selection.json',{'status':'no_unique_detector_box','source_bbox_pdf':box})
   box=None
  else:
   final=padded_fraction(selected['bbox'],pixel_size,x_fraction=0,y_fraction=0,min_pad=2)
   save(folder/'selection.json',{'status':'selected','method':'SAR MolDet + unique_detector_box + padded_fraction','source_bbox_pdf':box,'detector_bbox_pixels':selected['bbox'],'final_bbox_pixels':final,'page_image_sha256':hashlib.sha256(Path(page_image).read_bytes()).hexdigest()})
   box=[v/(sx if i%2==0 else sy) for i,v in enumerate(final)]
  if box is None:
   return {'images':[g['image']],'context_images':[g['image']]}
  rows=[{'crop':'fragment.png','page':g['page'],'bbox_pixels':box,'bbox_page_size':size}]
  mapping,audit=render_inputs(root/'input/source.pdf',rows,folder,folder/'final',dpi=300)
  src=Path(next(iter(mapping.values())))
  # Remove only tiny peripheral row remnants separated from the main drawing by whitespace.
  with Image.open(src) as im:
   gray=im.convert('L');w,h=im.size
   ink=[y for y in range(h) if gray.crop((0,y,w,y+1)).getextrema()[0]<180]
  bands=[]
  for y in ink:
   if not bands or y-bands[-1][-1]>6:bands.append([y])
   else:bands[-1].append(y)
  if len(bands)>1:
   main=max(bands,key=lambda b:b[-1]-b[0]);span=main[-1]-main[0]+1
   if all(b is main or b[-1]-b[0]+1<=span*.12 for b in bands):
    top=max(0,main[0]-3);bottom=min(h,main[-1]+4)
    a,b,c,d=audit[0]['pdf_bbox'];refined=[a,b+top/h*(d-b),c,b+bottom/h*(d-b)]
    rows[0]['bbox_pixels']=refined
    original=audit
    mapping,audit=render_inputs(root/'input/source.pdf',rows,folder,folder/'final',dpi=300)
    audit[0]['initial_render']=original[0];audit[0]['refinement']='Removed tiny isolated neighboring-row remnants at whitespace boundaries'
    src=Path(next(iter(mapping.values())))
  name='source-fragment-'+hashlib.sha256(src.read_bytes()).hexdigest()[:20]+'.png';shutil.copy2(src,assets/name);save(folder/'manifest.json',audit)
  g={**g,'image':marked_fragment_table(root,g,record,audit[0]['pdf_bbox'],assets)}
  fragment='assets/'+name
 return {'images':([fragment] if fragment else [])+[g['image']],'context_images':[g['image']]}
