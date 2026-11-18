"""Extract row-specific molecular cells, preserving merged-cell ownership."""
import re
from core import fitz

def merge_positions(values,tolerance=.8):
 groups=[]
 for x in sorted(values):
  if not groups or x-groups[-1][-1]>tolerance:groups.append([x])
  else:groups[-1].append(x)
 return [sum(g)/len(g) for g in groups]

def ruled_cells(page,bounds,columns,allowed_ids):
 x0,y0,x1,y1=bounds;horizontal=[];vertical=[]
 for drawing in page.get_drawings():
  for item in drawing['items']:
   if item[0]=='l':a,b=item[1:3];left,right=sorted([a.x,b.x]);top,bottom=sorted([a.y,b.y])
   elif item[0]=='re':left,top,right,bottom=item[1]
   else:continue
   if left<x0-2 or right>x1+2 or top<y0-2 or bottom>y1+2:continue
   if bottom-top<1.5 and right-left>8:horizontal.append((left,(top+bottom)/2,right))
   if right-left<1.5 and bottom-top>8:vertical.append(((left+right)/2,top,bottom))
 candidates=merge_positions([r[0] for r in vertical]);xs=[]
 for x in candidates:
  intervals=sorted((a,b) for vx,a,b in vertical if abs(vx-x)<.8);merged=[]
  for a,b in intervals:
   if merged and a<=merged[-1][1]+1:merged[-1][1]=max(merged[-1][1],b)
   else:merged.append([a,b])
  if sum(b-a for a,b in merged)>(y1-y0)*.7:xs.append(x)
 if len(xs)<max(columns)+2:raise ValueError('Incomplete vertical table grid')
 horizontal=[(a,y,b) for a,y,b in horizontal if min(abs(a-x) for x in xs)<1.2 and min(abs(b-x) for x in xs)<1.2]
 def spans(col):
  center=(xs[col]+xs[col+1])/2
  ys=merge_positions([y for left,y,right in horizontal if left-.8<=center<=right+.8])
  return list(zip(ys,ys[1:]))
 idcells=spans(0);results=[];seen=set()
 ordered_ids=list(allowed_ids) if not isinstance(allowed_ids,set) else None
 for row_index,(top,bottom) in enumerate(idcells):
  text=page.get_textbox(fitz.Rect(xs[0]+1,top+1,xs[1]-1,bottom-1)).strip()
  cid=re.sub(r'\s+','',text)
  if not cid and ordered_ids is not None and len(idcells)==len(ordered_ids)+1 and row_index>0:cid=ordered_ids[row_index-1]
  if cid not in allowed_ids:continue
  if cid in seen:raise ValueError('Repeated row identifier: '+cid)
  seen.add(cid);middle=(top+bottom)/2;parts={}
  for col,label in columns.items():
   matches=[(a,b) for a,b in spans(col) if a-.1<=middle<=b+.1]
   if len(matches)!=1:raise ValueError('Ambiguous merged-cell ownership: '+cid+' '+label)
   a,b=matches[0];box=[xs[col]+1.5,a+1.5,xs[col+1]-1.5,b-1.5]
   cell_id=f'{label}-{a:.2f}-{b:.2f}'
   parts[label]={'cell_id':cell_id,'bbox_pdf':box,'text':page.get_textbox(fitz.Rect(box)).strip(),'row_span_pdf':[a,b]}
  results.append({'compound_id':cid,'row_bbox_pdf':[xs[0],top,xs[-1],bottom],'parts':parts})
 return results
