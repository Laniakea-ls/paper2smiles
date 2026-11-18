"""Parse OCR HTML tables without dropping merged headers, qualifiers or footnotes."""
from html.parser import HTMLParser
import html
import re

def clean_text(value):
    text=html.unescape(str(value)).replace('$','')
    for old,new in [(r'\pm','±'),(r'\mu','μ'),(r'\leq','≤'),(r'\geq','≥'),(r'\lt','<'),(r'\gt','>'),(r'\times','×')]:
        text=text.replace(old,new)
    for _ in range(4):
        text=re.sub(r'\\(?:text|mathrm|mathbf|mathit|operatorname)\s*\{([^{}]*)\}',r'\1',text)
    text=re.sub(r'_\s*\{([^{}]+)\}',r'\1',text)
    text=re.sub(r'\^\s*\{([^{}]+)\}',r'^\1',text)
    text=re.sub(r'\\(?:left|right|quad|qquad)\b',' ',text)
    text=text.replace('{','').replace('}','')
    text=re.sub(r'\\[:,;!]', ' ', text)
    return ' '.join(text.split()).strip()

class TableParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tables=[];self.table=None;self.row=None;self.cell=None
    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs)
        if tag=='table':self.table=[]
        elif tag=='tr' and self.table is not None:self.row=[]
        elif tag in ['td','th'] and self.row is not None:
            self.cell={'raw_text':'','images':[],'rowspan':int(attrs.get('rowspan',1)),'colspan':int(attrs.get('colspan',1))}
        elif self.cell is not None:
            if tag=='img':self.cell['images'].append(attrs.get('src',''))
            elif tag=='sup':self.cell['raw_text']+='^'
            elif tag=='br':self.cell['raw_text']+=' '
    def handle_endtag(self,tag):
        if tag in ['td','th'] and self.cell is not None:
            self.cell['text']=clean_text(self.cell['raw_text']);self.row.append(self.cell);self.cell=None
        elif tag=='tr' and self.row is not None:
            self.table.append(self.row);self.row=None
        elif tag=='table' and self.table is not None:
            self.tables.append(expand_grid(self.table));self.table=None
    def handle_data(self,text):
        if self.cell is not None:self.cell['raw_text']+=text

def expand_grid(rows):
    slots={};width=0
    for i,row in enumerate(rows):
        j=0
        for cell in row:
            while (i,j) in slots:j+=1
            for y in range(i,i+cell['rowspan']):
                for x in range(j,j+cell['colspan']):
                    if (y,x) in slots:raise ValueError('Overlapping HTML table spans')
                    slots[(y,x)]={**cell,'origin':[i,j]}
                    width=max(width,x+1)
            j+=cell['colspan']
    return [[slots.get((i,j),{'text':'','raw_text':'','images':[],'origin':[i,j]}) for j in range(width)] for i in range(len(rows))]

def parse_tables(markdown):
    p=TableParser();p.feed(markdown);return p.tables

def split_marks(text):
    text=clean_text(text)
    marks=re.findall(r'\^\s*([a-z]+)\b',text)
    return re.sub(r'\^\s*[a-z]+\b','',text).strip(),list(dict.fromkeys(''.join(marks)))

def compound_identifier(text,known_marks):
    value,marks=split_marks(text)
    if any(mark not in known_marks for mark in marks):raise ValueError('Unknown compound footnote: '+str(marks))
    value=re.sub(r'^\(\s*±\s*\)\s*[-−]\s*', '(±)-', value)
    if not re.fullmatch(r'(?:\(±\)-)?[A-Za-z0-9][A-Za-z0-9._-]*',value):raise ValueError('Ambiguous compound identifier: '+value)
    return value,marks

NUMBER=r'(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
def measurement(text,endpoint,unit):
    text,marks=split_marks(text)
    if re.fullmatch(r'N\.?\s*D\.?',text,re.I):
        return {'endpoint':endpoint,'value_original':text,'unit':unit,'status':'not_determined','footnotes':marks,'relation':None,'mean':None,'sem':None}
    numeric_text=re.sub(r'\(\s*(±\s*[\d.eE+\-,]+)\s*\)',r'\1',text).replace(',','')
    match=re.fullmatch(rf'(<=|>=|<|>|≤|≥|=)?\s*({NUMBER})(?:\s*±\s*({NUMBER}))?',numeric_text)
    if not match:raise ValueError('Unparsed measurement: '+repr(text))
    return {'endpoint':endpoint,'value_original':text,'unit':unit,'status':'reported','footnotes':marks,
            'relation':match.group(1) or '=','mean':float(match.group(2)),
            'sem':float(match.group(3)) if match.group(3) else None}

def potency_cell(text,mode):
    if mode=='EC50':return [measurement(text,'EC50','nM')]
    if mode!='pEC50_with_EC50':raise ValueError('Unsupported potency cell format: '+str(mode))
    clean,marks=split_marks(text)
    match=re.fullmatch(r'(.*?)\s*\(\s*([^()]*)\s*\)',clean)
    if match:
        values=([measurement(match.group(1),'pEC50','')] if match.group(1).strip() else [])+[measurement(match.group(2),'EC50','nM')]
    else:values=[measurement(clean,'pEC50','')]
    for value in values:value['footnotes']=marks
    # EC50 is emitted only if actually printed in parentheses, never backfilled.
    return values
