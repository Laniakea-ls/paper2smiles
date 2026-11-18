"""Evidence-gated name completeness and cited-paper name recovery."""
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from rdkit import Chem
from core import name_fold, normalize_id, opsin_batch, save

VERSION = 'name-resolution-v2'

def name_status(parsed, si):
    if not parsed.get('valid'): return 'unparsed'
    if not si.get('valid'): return 'no_structure_comparator'
    if parsed['canonical'] == si['canonical']: return 'exact'
    if parsed['parent'] == si['parent']: return 'salt_form_difference'
    if parsed['connectivity'] == si['connectivity']: return 'stereo_difference'
    a, b = Chem.MolFromSmiles(parsed['parent']), Chem.MolFromSmiles(si['parent'])
    if a and b and a.GetNumHeavyAtoms() < b.GetNumHeavyAtoms() and b.HasSubstructMatch(a):
        return 'partial_structure_name'
    return 'structure_difference_review'

def audit_names(records, names):
    by_id = {}
    for r in records: by_id.setdefault(r['compound_id'], []).append(r)
    issues = []
    for n in names:
        matches = by_id.get(n['compound_id'], [])
        if not matches: continue
        statuses = [name_status(n['opsin'], r['si_structure']) for r in matches]
        status = next((s for s in ['exact', 'salt_form_difference', 'stereo_difference'] if s in statuses), statuses[0])
        n['completeness_status'] = status
        n['eligible_for_structure_comparison'] = status != 'partial_structure_name'
        if status in ['partial_structure_name', 'unparsed', 'structure_difference_review']:
            issues.append({'compound_id': n['compound_id'], 'reason': status, 'name': n['name']})
    for cid in by_id:
        if not any(n['compound_id'] == cid for n in names): issues.append({'compound_id': cid, 'reason': 'name_missing'})
    return issues

def xml_names(xml):
    """Only names printed as section headings with explicit local IDs are candidates."""
    tree = ET.fromstring(xml)
    result = []
    for title in tree.findall('.//sec/title'):
        raw = ''.join(title.itertext()).strip()
        match = re.search(r'\s*\(([A-Za-z0-9-]+)\)\s*$', raw)
        if not match: continue
        name = re.sub(r'^\d+(?:\.\d+)*\.?\s+', '', raw[:match.start()]).strip()
        if len(name) < 20: continue
        result.append({'name': name, 'original_compound_id': match[1], 'source_text': raw})
    return result

def recover_cited_names(flow, plan, records, names, request):
    """Retrieve up to three source-quoted references; never ask a model to invent a name."""
    folder = flow.root / 'name_resolution'; folder.mkdir(exist_ok=True)
    issues = audit_names(records, names)
    save(folder / 'initial_audit.json', issues)
    targets = {x['compound_id'] for x in issues if x['reason'] in ['partial_structure_name', 'name_missing', 'unparsed']}
    if not targets: return names
    text = '\n'.join(p['text'] for p in flow.pages.values())
    prompt = '''Identify up to THREE prior publications cited by this paper that describe the synthesis or full names of the unresolved compounds. Return {"references":[{"title":"exact cited paper title","citation_quote":"verbatim bibliography passage","reason":"source-based relevance"}]}. Only quoted references present in the supplied paper are allowed. No chemical names, structures, invented DOIs or URLs. Return an empty list if unsupported.\n'''
    refs = flow.llm.ask('name_reference_selection_v2', prompt + str(sorted(targets)) + '\n' + text).get('references', [])[:3]
    audit = []
    for index, ref in enumerate(refs):
        event = {'reference': ref, 'accepted': []}
        try:
            quote = ref.get('citation_quote', '')
            if len(quote) < 30 or name_fold(quote) not in name_fold(text) or name_fold(ref['title']) not in name_fold(quote):
                raise ValueError('Reference not grounded in original bibliography')
            hits = request('GET', 'https://www.ebi.ac.uk/europepmc/webservices/rest/search', params={'query': 'TITLE:"' + ref['title'].replace('"', '') + '"', 'format': 'json', 'pageSize': 10}).json()['resultList']['result']
            hits = [h for h in hits if name_fold(h.get('title', '').rstrip('.')) == name_fold(ref['title'].rstrip('.')) and h.get('pmcid')]
            if len(hits) != 1: raise ValueError('No unique accessible exact-title reference')
            hit = hits[0]; pmcid = hit['pmcid']
            if not re.fullmatch(r'PMC\d+', pmcid): raise ValueError('Invalid PMC identifier')
            xml = request('GET', f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML').text
            (folder / (pmcid + '.xml')).write_text(xml)
            candidates = xml_names(xml)
            parsed = opsin_batch([c['name'] for c in candidates], folder / pmcid)
            for c, parsed_name in zip(candidates, parsed):
                if not parsed_name['valid']: continue
                matching = {r['compound_id'] for r in records if r['compound_id'] in targets and r['si_structure']['valid'] and parsed_name['canonical'] == r['si_structure']['canonical']}
                # A cross-paper ID mapping requires exact stereo-aware structure identity
                # and a unique current ID; similar names or equal activity are insufficient.
                if len(matching) != 1: continue
                cid = next(iter(matching))
                row = {**c, 'compound_id': cid, 'role': 'target', 'source_ranges': [], 'name_text_supported': True,
                       'opsin': parsed_name, 'source_kind': 'cited_reference', 'source_reference': {**hit, 'url': f'https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/', 'citation_quote': quote},
                       'mapping_basis': 'exact isomeric canonical structure identity with author SI', 'completeness_status': 'exact', 'eligible_for_structure_comparison': True}
                names.append(row); event['accepted'].append({'current_id': cid, 'original_id': c['original_compound_id'], 'name': c['name']})
            event['status'] = 'searched'
        except Exception as exc:
            event['status'] = 'unresolved'; event['error'] = str(exc)
        audit.append(event)
        save(folder / 'reference_audit.json', audit)
    audit_names(records, names)
    return names
