"""Luna reads column/footnote schemas; values are parsed from PaddleOCR cells."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import getpass
import json
from pathlib import Path
import re

from core import save,sha256
from model import Luna
from ocr_table_data import parse_tables,compound_identifier,measurement,potency_cell,clean_text

def load_tables(root):
    tables=[]
    for source in json.loads((root/'inputs.json').read_text()):
        path=root/'paddle'/(source['id']+'.md');markdown=path.read_text()
        grids=parse_tables(markdown)
        if len(grids)!=1:raise ValueError('Expected exactly one OCR table in '+source['id'])
        tables.append({**source,'grid':grids[0],'markdown':markdown,'ocr_file':str(path),'ocr_sha256':sha256(path)})
    return tables

def schema_review(model,table):
    prompt='''Read the header and footnotes of this source table image and its local PaddleOCR transcription. This paper studies GPR88. Return a schema for deterministic extraction of every OCR row. Do NOT return, edit or fill any body-row values, identifiers or structures. Numbers will be taken directly from the OCR cells. No SI or prior extraction is supplied.
Return {"caption":"actual table caption","header_rows":1,"column_titles":["plain headers in original column order, retaining nested heading context and units"],"compound_id_column":0,"potency_column":3,"potency_format":"pEC50_with_EC50 or EC50","emax_column":4,"target":"GPR88","assay":"cAMP","default_emax_reference":"reference compound name or null","emax_reference_by_marker":{"d":"reference compound name"},"default_compound_markers":[],"compound_context_by_marker":{"e":"source meaning, e.g. HCl salt","c":"only if the compound marker cites earlier data"},"footnotes":[{"marker":"a","text":"literal English source footnote"}],"measurement_error":"standard error / standard deviation / unspecified, as printed","concerns":[{"row":0,"column":0,"issue":"possible OCR problem to check against the image"}]}.
All row and column indexes are ZERO-BASED in the supplied expanded grid. header_rows is the number of original header rows (including nested headers). Null is required for absent Emax. If a potency cell prints pEC50 with EC50 in parentheses use pEC50_with_EC50; if only EC50 is printed use EC50. pEC50 has no concentration unit; the parenthesized EC50 is nM here. Preserve blank cells; PK/solubility are distinct from potency. Emax reference comes from its column/header footnote, with row superscript overrides. A superscript on a compound identifier is separate from an Emax override. A compound-column header footnote applies to all rows. Concerns identify ambiguities only; never propose invented values or molecule names. Image-only R substituents are not OCR SMILES.
'''+json.dumps({'source_table':table['table_number'],'pdf_page':table['page'],
        'grid':[[c['text'] or ('[structure image]' if c['images'] else '') for c in row] for row in table['grid']],
        'ocr_markdown':table['markdown']},ensure_ascii=False)
    return model.ask('ocr_table_schema_'+str(table['table_number']),prompt,[Path(table['image'])],max_tokens=7000)

def compile_table(table,schema):
    grid=table['grid'];width=len(grid[0]);header=schema['header_rows']
    if not isinstance(header,int) or not 1<=header<len(grid):raise ValueError('Invalid header row count')
    if len(schema['column_titles'])!=width:raise ValueError('Column title count mismatch')
    for key in ['compound_id_column','potency_column','emax_column']:
        index=schema.get(key)
        if index is None and key=='emax_column':continue
        if not isinstance(index,int) or not 0<=index<width:raise ValueError('Invalid column '+key)
    known={f['marker'] for f in schema['footnotes']};rows=[]
    for index,cells in enumerate(grid[header:],header):
        raw=[c['text'] for c in cells]
        cid,marks=compound_identifier(raw[schema['compound_id_column']],known)
        annotations=list(dict.fromkeys(schema.get('default_compound_markers',[])+marks))
        if any(x not in known for x in annotations):raise ValueError('Unrecognized compound marker')
        values=potency_cell(raw[schema['potency_column']],schema['potency_format'])
        for value in values:value['source_column']=schema['potency_column']+1
        emax_reference=None
        if schema['emax_column'] is not None:
            emax=measurement(raw[schema['emax_column']],'Emax','%')
            emax_reference=schema['default_emax_reference']
            for mark in emax['footnotes']:
                if mark not in schema['emax_reference_by_marker']:raise ValueError('Unresolved Emax reference footnote '+mark)
                emax_reference=schema['emax_reference_by_marker'][mark]
            emax['source_column']=schema['emax_column']+1;emax['reference_compound']=emax_reference
            values.append(emax)
        rows.append({'record_id':f"T{table['table_number']}:R{index+1}",'table_number':table['table_number'],'page':table['page'],
            'source_row':index+1,'compound_id':cid,'raw_identifier':raw[schema['compound_id_column']],
            'compound_markers':annotations,'compound_context':[schema.get('compound_context_by_marker',{}).get(m,'') for m in annotations],
            'target':schema['target'],'assay':schema['assay'],'measurements':values,'raw_cells':raw,
            'source_type':'main_pdf_paddleocr_table','value_origin':'unmodified OCR cell; deterministic parsing',
            'error_type':schema['measurement_error']})
    return {k:v for k,v in table.items() if k!='markdown'}|{'schema':schema,'rows':rows}

def compile_all(root,tables):
    output=[]
    for table in tables:
        schema=json.loads((root/'schemas'/f"table-{table['table_number']}.json").read_text())
        output.append(compile_table(table,schema))
    save(root/'tables.json',output)
    rows=[r for table in output for r in table['rows']]
    values=[{k:v for k,v in row.items() if k!='measurements'}|measurement for row in rows for measurement in row['measurements']]
    save(root/'activity_rows.json',rows);save(root/'activity_values.json',values)
    summary={'tables':len(output),'table_rows':len(rows),'unique_compounds':len({r['compound_id'] for r in rows}),
        'rows_by_table':{str(t['table_number']):len(t['rows']) for t in output},
        'endpoints':dict(Counter(v['endpoint'] for v in values)),
        'reported_numeric_endpoints':sum(v['status']=='reported' for v in values),
        'not_determined_endpoints':sum(v['status']=='not_determined' for v in values),
        'qualified_endpoints':sum(v['relation'] not in ['=',None] for v in values),
        'ocr_engine':'PaddleOCR-VL-1.6-0.9B + PP-DocLayoutV3 (existing local services)',
        'schema_model':'GPT-6-Luna','numeric_value_source':'PaddleOCR cells; Luna only interprets headers/footnotes',
        'no_si_or_prior_extraction_values_used':True,
        'scope':'Main-PDF Tables 1–5. Table 5 contains four additional row occurrences, including repeats and previously cited comparator data.',
        'finished_utc':datetime.now(timezone.utc).isoformat()}
    save(root/'summary.json',summary);print(json.dumps(summary,ensure_ascii=False),flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-dir',required=True,type=Path)
    parser.add_argument('--compile-only',action='store_true');args=parser.parse_args();root=args.run_dir.resolve()
    tables=load_tables(root);save(root/'ocr_grid.json',[{k:v for k,v in t.items() if k!='markdown'} for t in tables])
    key=''
    try:
        if not args.compile_only:
            key=getpass.getpass('API key (hidden; memory only): ')
            model=Luna(key,'https://newapi.livzon.cn/v1','GPT-6-Luna',root)
            with ThreadPoolExecutor(max_workers=2) as pool:
                jobs={pool.submit(schema_review,model,t):t for t in tables}
                for future in as_completed(jobs):
                    table=jobs[future];schema=future.result();save(root/'schemas'/f"table-{table['table_number']}.json",schema)
                    print('Schema reviewed for table',table['table_number'],flush=True)
        compile_all(root,tables)
    except Exception as exc:
        print(str(exc).replace(key,'[REDACTED]') if key else str(exc),flush=True)
        raise SystemExit(1)
    finally:key=''

if __name__=='__main__':main()
