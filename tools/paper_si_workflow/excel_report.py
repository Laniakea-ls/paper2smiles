"""Human-facing activity workbook, with in-sheet SMILES depictions like SAR reports."""
import io,re
from collections import OrderedDict
from pathlib import Path
from openpyxl import Workbook
from openpyxl.styles import Alignment,Font,PatternFill,Border,Side
from openpyxl.drawing.image import Image
from openpyxl.drawing.spreadsheet_drawing import OneCellAnchor,AnchorMarker
from openpyxl.drawing.xdr import XDRPositiveSize2D
from openpyxl.utils import get_column_letter
from openpyxl.utils.units import pixels_to_EMU
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from core import sha256

def merged_activity(rows):
 columns=list(dict.fromkeys((r.get('endpoint',''),r.get('unit','')) for r in rows));groups=OrderedDict()
 for r in rows:
  key=tuple(r.get(k) for k in ('compound_id','target','assay','page','caption','table_id','source_row','reference','source_type'))
  group=groups.setdefault(key,{'source':r,'values':{}})
  column=(r.get('endpoint',''),r.get('unit',''));value=str(r.get('value',''));relation=r.get('relation','')
  if relation in ('<','>','<=','>=','≤','≥') and not re.match(r'^[<>=≤≥]',value.strip()):value=relation+value
  group['values'].setdefault(column,[]).append(value)
 return columns,list(groups.values())

def number_or_text(values):
 if not values:return '—'
 if len(values)>1:return '\n'.join(values)
 value=values[0]
 if re.fullmatch(r'[+-]?\d+',value):return int(value)
 if re.fullmatch(r'[+-]?(?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?',value):return float(value)
 return value

def write_workbook(root,body,selections,si,destination):
 root=Path(root);columns,groups=merged_activity(body);book=Workbook();sheet=book.active;sheet.title='化合物活性报告'
 headers=['编号','SMILES','2D 结构','实验条件']+[e+(' ('+u+')' if u else '') for e,u in columns]+['出处']
 sheet.append(headers);sheet.row_dimensions[1].height=30
 headerfill=PatternFill('solid',fgColor='E5EDF5');line=Side(style='thin',color='E3E9EF')
 for c in sheet[1]:c.font=Font(name='Calibri',size=12,bold=True,color='334A60');c.fill=headerfill;c.alignment=Alignment(horizontal='center',vertical='center',wrap_text=True)
 widths=[10,55,40,32]+[16]*len(columns)+[62]
 for i,width in enumerate(widths,1):sheet.column_dimensions[get_column_letter(i)].width=width
 image_cache={};buffers=[];digest=sha256(root/'input/source.pdf')
 for rownum,group in enumerate(groups,2):
  r=group['source'];cid=r['compound_id'];choice=selections[cid];smiles=choice['smiles'];target=r.get('target') or '';assay=r.get('assay') or ''
  condition=' · '.join(x for x in [target if target.lower() not in assay.lower() else '',assay] if x)
  source=(r.get('caption') or '')+(' · PDF '+str(r['page']) if r.get('page') else '')
  values=[cid,smiles,'',condition]+[number_or_text(group['values'].get(column)) for column in columns]+[source]
  sheet.append(values);sheet.row_dimensions[rownum].height=150
  for c in sheet[rownum]:
   if isinstance(c.value,str):c.data_type='s'
   c.font=Font(name='Calibri',size=11,color='334A60');c.alignment=Alignment(vertical='center',horizontal='left' if c.column in (2,4,len(headers)) else 'center',wrap_text=True);c.border=Border(bottom=line)
   if rownum%2==0:c.fill=PatternFill('solid',fgColor='FAFCFE')
  if smiles:
   if smiles not in image_cache:
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:raise ValueError('Invalid selected structure: '+cid)
    from source_depiction import align_to_source
    aligned=align_to_source(mol,digest);draw=rdMolDraw2D.MolDraw2DCairo(810,540)
    if aligned:draw.DrawMolecule(rdMolDraw2D.PrepareMolForDrawing(mol,addChiralHs=False,wedgeBonds=False))
    else:rdMolDraw2D.PrepareAndDrawMolecule(draw,mol)
    draw.FinishDrawing();image_cache[smiles]=draw.GetDrawingText()
   buffer=io.BytesIO(image_cache[smiles]);buffers.append(buffer);image=Image(buffer)
   image.anchor=OneCellAnchor(_from=AnchorMarker(col=2,row=rownum-1,colOff=pixels_to_EMU(5),rowOff=pixels_to_EMU(8)),ext=XDRPositiveSize2D(pixels_to_EMU(270),pixels_to_EMU(180)))
   sheet.add_image(image)
 sheet.freeze_panes='A2';sheet.auto_filter.ref=sheet.dimensions;sheet.sheet_view.showGridLines=False
 sheet.print_title_rows='1:1';sheet.sheet_properties.pageSetUpPr.fitToPage=True;sheet.page_setup.orientation='landscape';sheet.page_setup.paperSize=sheet.PAPERSIZE_A3;sheet.page_setup.fitToWidth=1;sheet.page_setup.fitToHeight=0
 # Keep raw measurements separately, without replacing article values with SI values.
 for title,rows in [('正文原始测定',body),('SI 原始测定',si)]:
  ws=book.create_sheet(title);fields=list(dict.fromkeys(k for r in rows for k in r));ws.append(fields or ['无数据'])
  import json
  for r in rows:
   ws.append([json.dumps(r.get(k),ensure_ascii=False) if isinstance(r.get(k),(list,dict)) else r.get(k) for k in fields])
   for c in ws[ws.max_row]:
    if isinstance(c.value,str):c.data_type='s'
  ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
 destination=Path(destination);temp=destination.with_suffix('.tmp.xlsx');book.save(temp);temp.replace(destination)
 return {'rows':len(groups),'structures':sum(bool(selections[g['source']['compound_id']]['smiles']) for g in groups),'measurements':sum(len(v) for g in groups for v in g['values'].values())}
