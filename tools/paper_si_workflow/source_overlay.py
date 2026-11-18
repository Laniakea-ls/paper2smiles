"""Non-destructive crop rectangles on a full source table, in PDF coordinates."""
import base64,hashlib,json
from pathlib import Path
from html import escape
from PIL import Image
from core import save

def table_overlay(root,group,boxes,assets):
 image=Path(root)/'report'/group['image']
 with Image.open(image) as im:w,h=im.size
 x,y,u,v=group['render']['pdf_bbox'];rects=[];audit=[]
 for label,box,color in boxes:
  a,b,c,d=box
  if not(x<=a<c<=u and y<=b<d<=v):raise ValueError('Source overlay lies outside table context: '+label)
  left=(a-x)/(u-x)*w;top=(b-y)/(v-y)*h;width=(c-a)/(u-x)*w;height=(d-b)/(v-y)*h
  rects.append(f'<rect x="{left:.3f}" y="{top:.3f}" width="{width:.3f}" height="{height:.3f}" fill="none" stroke="{color}" stroke-width="4"><title>{escape(label)}</title></rect>')
  audit.append({'label':label,'pdf_bbox':box,'image_bbox':[left,top,left+width,top+height],'color':color})
 data=base64.b64encode(image.read_bytes()).decode()
 svg=f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><image width="{w}" height="{h}" xlink:href="data:image/png;base64,{data}"/>'+''.join(rects)+'</svg>'
 name='source-context-'+hashlib.sha256(svg.encode()).hexdigest()[:20]+'.svg';(assets/name).write_text(svg)
 save(Path(root)/'report_table_sources'/(name+'.json'),{'page':group['page'],'table_id':group['table_id'],'table_pdf_bbox':[x,y,u,v],'image_size':[w,h],'boxes':audit})
 return 'assets/'+name

def scaffold_box(root,record):
 if not record:return None
 if record.get('fragments'):
  from source_depiction import PAPER_SHA256
  if record.get('source_pdf_sha256')==PAPER_SHA256 and record.get('table_id')=='page-003-t3':return [374,594,480,651]
  return None
 if record.get('scaffold_bbox_pdf_points'):return record['scaffold_bbox_pdf_points']
 root=Path(root)
 if record.get('origin')=='verified_local_trial':
  p=Path(__file__).resolve().parents[2]/'runs/gpr88_table1_reconstruction'
  if json.loads((p/'summary.json').read_text())['source_pdf_sha256']==record.get('source_pdf_sha256'):return json.loads((p/'crops.json').read_text())['boxes_pdf_points']['scaffold']
 candidates=[root/'markush'/record['table_id']/'sar_render_manifest.json']
 if record['table_id']=='page-003-t2':candidates.append(root/'markush/table3-sar-retry/sar_render_manifest.json')
 for path in candidates:
  if path.exists():return json.loads(path.read_text())['crops'][0]['pdf_bbox']
 layout_path=root/'markush'/record['table_id']/'layout.json'
 if layout_path.exists():
  layout=json.loads(layout_path.read_text());box=layout.get('scaffold_bbox',[])
  if layout.get('supported') is True and len(box)==4 and 0<=box[0]<box[2]<=1000 and 0<=box[1]<box[3]<=1000:
   from core import fitz
   with fitz.open(root/'input/source.pdf') as pdf:
    page=pdf[record['page']-1]
    return [v*(page.rect.width if i%2==0 else page.rect.height)/1000 for i,v in enumerate(box)]
 return None

def marked_fragment_table(root,g,record,box,assets):
 boxes=[('片段裁剪范围',box,'#dc4444')];core=scaffold_box(root,record)
 if core:boxes.insert(0,('通式裁剪范围',core,'#2563eb'))
 return table_overlay(root,g,boxes,assets)
