"""Select compounds by actual bioactivity evidence before structure lookup."""
import math
from collections import defaultdict
from core import normalize_id,parse_number

BIOACTIVITY_ENDPOINTS={'EC50','PEC50','IC50','PIC50','KI','PKI','KD','PKD','EMAX','INHIBITION','ACTIVATION'}

def numeric_activity(value):
    if str(value.get('endpoint','')).replace('_','').upper() not in BIOACTIVITY_ENDPOINTS:return False
    if value.get('status') in ['not_determined','missing']:return False
    mean=value.get('mean')
    if isinstance(mean,(int,float)) and not isinstance(mean,bool):return math.isfinite(mean)
    return parse_number(value.get('value_original',value.get('value',''))) is not None

def select_targets(rows):
    grouped=defaultdict(list)
    for row in rows:
        cid=normalize_id(row.get('compound_id'))
        if not cid:continue
        grouped[cid].append(row)
    targets={};excluded=[]
    for cid,records in grouped.items():
        values=[v for r in records for v in r.get('measurements',[r]) if numeric_activity(v)]
        if not values:
            excluded.append({'compound_id':cid,'reason':'no numeric bioactivity endpoint'});continue
        def relation(v):
            parsed=parse_number(v.get('value_original',v.get('value','')))
            return v.get('relation') or (parsed or {}).get('relation') or '='
        targets[cid]={'compound_id':cid,'activity_records':records,
                      'only_limits':all(relation(v)!='=' for v in values),
                      'selection_basis':'at least one reported numeric bioactivity value or bound'}
    return targets,excluded
