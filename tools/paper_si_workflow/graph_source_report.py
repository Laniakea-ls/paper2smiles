"""Source crop metadata for graph evidence; layout is owned by the frontend."""
import hashlib,shutil
from pathlib import Path
from source_overlay import table_overlay

def graph_source(root,record,group,assets):
    images=[];boxes=[]
    for node in record['graph_evidence']:
        src=Path(node['crop_path']);name='graph-source-'+hashlib.sha256(src.read_bytes()).hexdigest()[:20]+'.png'
        shutil.copy2(src,assets/name);images.append('assets/'+name)
        label=node.get('node_label',node.get('cell_id','片段').split('-')[0])
        if node['page']==group['page']:boxes.append((label,node['tight_bbox_pdf'],'#2563eb' if label=='通式' else '#dc4444'))
    context=table_overlay(root,group,boxes,assets)
    return {'images':images+[context],'context_images':[context]}
