"""Targeted second search for missing names, using printed IDs to locate pages."""
import re
from pathlib import Path
from core import normalize_id,opsin_batch,save

def retry_missing(flow,plan,names):
    supported={n['compound_id'] for n in names if n['opsin']['valid'] and (n.get('name_text_supported') or n.get('visual_name_review',{}).get('supported'))}
    missing=flow.target_ids-supported
    pages=set()
    for page in flow.pages.values():
        for cid in missing:
            # Parenthesized synthetic compound IDs, not loose matches to numeric measurements.
            if re.search(r'\('+re.escape(cid)+r'\)',page['text'],re.I):
                pages.add(page['page'])
    if not pages or not missing:return names
    from workflow import Workflow
    helper=Workflow(flow.pdf,flow.root/'name_retry',flow.llm);helper.pages=flow.pages;helper.target_ids=missing
    helper.status('对缺失完整名称做编号定点复查',targets=sorted(missing),pages=sorted(pages))
    # Keep the search bounded to actual ID locations, rather than all experimental pages again.
    selected={p:flow.pages[p] for p in pages}
    helper.pages=flow.pages
    recovered=helper.names({**plan,'compound_name_pages':sorted(pages),'targeted_name_search':True})
    recovered=helper.resolve_names(recovered)
    existing={(n['compound_id'],n['name']) for n in names}
    for n in recovered:
        if (n['compound_id'],n['name']) not in existing:names.append(n)
    save(flow.root/'name_retry/results.json',recovered)
    return names

def normalize_explicit_n_dimethyl(names,root):
    for i,n in enumerate(names):
        if n['opsin']['valid']:continue
        if not (n.get('name_text_supported') or n.get('visual_name_review',{}).get('supported')):continue
        revised=re.sub(r'\bN-[Dd]imethyl\b','N,N-dimethyl',n['name'])
        if revised==n['name']:continue
        parsed=opsin_batch([revised],Path(root)/'opsin/explicit_n_dimethyl'/str(i))[0]
        if parsed['valid']:
            n['opsin_original']=n['opsin'];n['opsin']=parsed
            n['parser_normalization']={'source_name':n['name'],'parser_input':revised,'rule':'N-dimethyl explicitly denotes two methyl substituents on N; expand repeated N locant for parser, retain printed source name'}
    return names
