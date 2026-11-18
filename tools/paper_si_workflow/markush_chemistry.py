"""Source-tagged single-attachment assembly and two-color rendering."""
import re,copy,hashlib
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
def mapped(smi):
 m=Chem.MolFromSmiles(smi)
 if m is None:raise ValueError('Invalid recognized SMILES')
 ds=[a for a in m.GetAtoms() if a.GetAtomicNum()==0]
 if len(ds)!=1 or ds[0].GetDegree()!=1:raise ValueError('Need exactly one terminal attachment point')
 ds[0].SetAtomMapNum(1)
 return m
def assemble(core,fragment):
 c,f=mapped(core),mapped(fragment)
 for part,mol in [("scaffold",c),("r1",f)]:
  for atom in mol.GetAtoms():atom.SetProp("source_part",part)
 m=Chem.molzip(Chem.CombineMols(c,f));Chem.SanitizeMol(m);Chem.AssignStereochemistry(m,cleanIt=True,force=True)
 if any(a.GetAtomicNum()==0 for a in m.GetAtoms()) or len(Chem.GetMolFrags(m))!=1:raise ValueError('Unresolved attachment or disconnected product')
 assert m.GetNumHeavyAtoms()==c.GetNumHeavyAtoms()+f.GetNumHeavyAtoms()
 return m

def write_colored_svg(m,assets,name,source_pdf_hash=None):
 from source_depiction import align_to_source
 aligned=align_to_source(m,source_pdf_hash)
 colors={'scaffold':'#2563EB','r1':'#D97706'}
 d=rdMolDraw2D.MolDraw2DSVG(650,350);d.drawOptions().useBWAtomPalette();d.DrawMolecule(rdMolDraw2D.PrepareMolForDrawing(m,addChiralHs=False,wedgeBonds=False)) if aligned else rdMolDraw2D.PrepareAndDrawMolecule(d,m);d.FinishDrawing()
 tree=ET.fromstring(d.GetDrawingText())
 def color(atom):
  a=m.GetAtomWithIdx(atom);return colors[a.GetProp('source_part') if a.HasProp('source_part') else 'scaffold']
 def paint(el,c):
  style=el.get('style','')
  style=re.sub(r'(stroke|fill):(?!(?:none)\b)[^;]+',lambda x:x[1]+':'+c,style)
  el.set('style',style)
  if el.get('fill') and el.get('fill')!='none':el.set('fill',c)
  if el.get('stroke') and el.get('stroke')!='none':el.set('stroke',c)
 for parent in list(tree.iter()):
  for el in list(parent):
   ids=[int(x) for x in re.findall(r'atom-(\d+)',el.get('class',''))]
   if not ids:continue
   if len(ids)>1 and color(ids[0])!=color(ids[1]):
    match=re.fullmatch(r'M\s*([\d.-]+),([\d.-]+)\s*L\s*([\d.-]+),([\d.-]+)',el.get('d','').strip())
    if not match:
     triangle=re.fullmatch(r'M\s*([\d.-]+),([\d.-]+)\s*L\s*([\d.-]+),([\d.-]+)\s*L\s*([\d.-]+),([\d.-]+)\s*Z',el.get('d','').strip())
     if not triangle:raise ValueError('Unexpected attachment bond drawing')
     x,y,u,v,w,z=map(float,triangle.groups());a,b=(x+u)/2,(y+v)/2;c,e=(x+w)/2,(y+z)/2
     pos=d.GetDrawCoords(ids[0]);near=(pos.x-x)**2+(pos.y-y)**2 < (pos.x-(u+w)/2)**2+(pos.y-(v+z)/2)**2
     first,second=(color(ids[0]),color(ids[1])) if near else (color(ids[1]),color(ids[0]))
     other=copy.deepcopy(el);el.set('d',f'M {x},{y} L {a},{b} L {c},{e} Z');other.set('d',f'M {a},{b} L {u},{v} L {w},{z} L {c},{e} Z')
     paint(el,first);paint(other,second);parent.insert(list(parent).index(el)+1,other)
     continue
    x,y,u,v=map(float,match.groups());mx,my=(x+u)/2,(y+v)/2
    other=copy.deepcopy(el);el.set('d',f'M {x},{y} L {mx},{my}');other.set('d',f'M {mx},{my} L {u},{v}')
    paint(el,color(ids[0]));paint(other,color(ids[1]));parent.insert(list(parent).index(el)+1,other)
   else:paint(el,color(ids[0]))
 ET.register_namespace('','http://www.w3.org/2000/svg')
 content=ET.tostring(tree,encoding='unicode');version=hashlib.sha256(content.encode()).hexdigest()[:12]
 filename=name+'-'+version+'.svg';(assets/filename).write_text(content);return 'assets/'+filename

def assemble_multiple(core,fragments,atom_variants=None):
 """Join explicitly mapped independent sites and apply explicitly mapped ring variants."""
 m=Chem.MolFromSmiles(core)
 if m is None:raise ValueError('Invalid mapped scaffold')
 for a in m.GetAtoms():a.SetProp('source_part','scaffold')
 sites={a.GetAtomMapNum() for a in m.GetAtoms() if a.GetAtomicNum()==0}
 if 0 in sites or sites!={int(k) for k in fragments}:raise ValueError('Every attachment site must have exactly one supplied fragment')
 for key,symbol in (atom_variants or {}).items():
  hits=[a for a in m.GetAtoms() if a.GetAtomMapNum()==int(key) and a.GetAtomicNum()!=0]
  if len(hits)!=1 or symbol not in ('CH','N'):raise ValueError('Unsupported or ambiguous ring variant')
  a=hits[0]
  if not a.GetIsAromatic():raise ValueError('Ring variant must be aromatic')
  a.SetAtomicNum(7 if symbol=='N' else 6);a.SetNumExplicitHs(0 if symbol=='N' else 1);a.SetNoImplicit(True)
  if symbol=='N':a.SetProp('source_part','r1')
 for key,smi in fragments.items():
  params=Chem.SmilesParserParams();params.removeHs=False
  f=Chem.MolFromSmiles(smi,params)
  if f is None:raise ValueError('Invalid fragment')
  ds=[a for a in f.GetAtoms() if a.GetAtomicNum()==0]
  if len(ds)!=1 or ds[0].GetDegree()!=1:raise ValueError('Fragment must have one terminal attachment')
  ds[0].SetAtomMapNum(int(key))
  for a in f.GetAtoms():a.SetProp('source_part','r1')
  m=Chem.CombineMols(m,f)
 m=Chem.molzip(m);Chem.SanitizeMol(m);m=Chem.RemoveHs(m)
 for a in m.GetAtoms():a.SetAtomMapNum(0)
 Chem.AssignStereochemistry(m,cleanIt=True,force=True)
 if any(a.GetAtomicNum()==0 for a in m.GetAtoms()) or len(Chem.GetMolFrags(m))!=1:raise ValueError('Unresolved multiple-site assembly')
 return m

def assemble_record(record):
 if 'graph_nodes' in record:return assemble_graph(record['graph_nodes'])
 if 'fragments' in record:return assemble_multiple(record['scaffold_smiles'],record['fragments'],record.get('atom_variants'))
 return assemble(record['scaffold_smiles'],record['r1_smiles'])

def assemble_graph(nodes):
 """Connect a graph of fragments. Each labelled port occurs exactly twice, across nodes."""
 parts=[];owners={}
 for name,smi in nodes.items():
  m=Chem.MolFromSmiles(smi)
  if m is None:raise ValueError('Invalid graph node: '+name)
  for a in m.GetAtoms():
   a.SetProp('source_part','scaffold' if name=='core' else 'r1');a.SetProp('source_node',name)
   if a.GetAtomicNum()==0:
    key=a.GetAtomMapNum()
    if key<=0 or a.GetDegree()!=1:raise ValueError('Every graph attachment needs an explicit label and one bond')
    owners.setdefault(key,[]).append(name)
  parts.append(m)
 if not parts or any(len(v)!=2 or len(set(v))!=2 for v in owners.values()):raise ValueError('Each graph edge needs exactly two ports from different fragments')
 m=parts[0]
 for part in parts[1:]:m=Chem.CombineMols(m,part)
 m=Chem.molzip(m);Chem.SanitizeMol(m)
 if len(Chem.GetMolFrags(m))!=1 or any(a.GetAtomicNum()==0 for a in m.GetAtoms()):raise ValueError('Disconnected or unresolved fragment graph')
 for a in m.GetAtoms():a.SetAtomMapNum(0)
 Chem.AssignStereochemistry(m,cleanIt=True,force=True)
 return m
