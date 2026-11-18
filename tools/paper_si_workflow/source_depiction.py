"""Source-specific depiction anchors. Changes coordinates, never chemical identity."""
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Geometry import Point3D
PAPER_SHA256='72f5178fb70c813e8313ada31a5072f3281c55b6c092fb0ec84795b676c16ee7'

def align_to_source(m,pdf_hash):
 if pdf_hash!=PAPER_SHA256:return False
 # Shared para-alkoxy phenylglycinamide geometry, visually read from the source PDF.
 # O-aryl at lower left; benzylic C upper right; N above it; carbonyl up/right.
 template=Chem.MolFromSmiles('Oc1ccc(CNC=O)cc1')
 coords=[(-1.30,-.75),(0,0),(0,1.50),(1.30,2.25),(2.60,1.50),(3.90,2.25),(3.90,3.75),(5.20,4.50),(5.20,6.00),(2.60,0),(1.30,-.75)]
 match=m.GetSubstructMatch(Chem.MolFromSmarts('[#8]-c1ccc(CNC=O)c[c,n]1'))
 if not match:return False
 before=Chem.MolToSmiles(m,isomericSmiles=True)
 conf=Chem.Conformer(template.GetNumAtoms());conf.Set3D(False)
 for i,(x,y) in enumerate(coords):conf.SetAtomPosition(i,Point3D(x,y,0))
 template.AddConformer(conf)
 anchors={idx:coords[i] for i,idx in enumerate(match)}
 # Preserve the source's leftward acyclic alkoxy chain and downward methyl branch.
 oxygen=match[0];aryl=match[1]
 starts=[a.GetIdx() for a in m.GetAtomWithIdx(oxygen).GetNeighbors() if a.GetIdx()!=aryl]
 def branch_size(idx,parent,seen=None):
  seen=set() if seen is None else seen
  if idx in seen:return 0
  seen=seen|{idx}
  return 1+sum(branch_size(a.GetIdx(),idx,seen) for a in m.GetAtomWithIdx(idx).GetNeighbors() if a.GetIdx()!=parent and a.GetIdx() not in seen)
 if len(starts)==1:
  first=starts[0];second=[a.GetIdx() for a in m.GetAtomWithIdx(first).GetNeighbors() if a.GetIdx()!=oxygen]
  if not m.GetAtomWithIdx(first).IsInRing() and len(second)==1:
   center_left=second[0]
   if not m.GetAtomWithIdx(center_left).IsInRing():
    branches=[a.GetIdx() for a in m.GetAtomWithIdx(center_left).GetNeighbors() if a.GetIdx()!=first]
    if branches and all(m.GetAtomWithIdx(i).GetAtomicNum()==6 and not m.GetAtomWithIdx(i).IsInRing() for i in branches):
     anchors[first]=(-2.6,0);anchors[center_left]=(-3.9,-.75)
     branches.sort(key=lambda i:branch_size(i,center_left),reverse=True)
     if len(branches)==2 and branch_size(branches[1],center_left)==1:anchors[branches[1]]=(-3.9,-2.25)
     previous=center_left;current=branches[0];step=0
     while current not in anchors and not m.GetAtomWithIdx(current).IsInRing():
      anchors[current]=(-5.2-1.3*step,0 if step%2==0 else -.75)
      following=[a.GetIdx() for a in m.GetAtomWithIdx(current).GetNeighbors() if a.GetIdx()!=previous]
      if len(following)!=1:break
      previous,current=current,following[0];step+=1
 # Orient the variable polar side chain as drawn under R3.
 substituents=[a.GetIdx() for a in m.GetAtomWithIdx(match[5]).GetNeighbors() if a.GetIdx() not in (match[4],match[6])]
 if len(substituents)==1:
  c=substituents[0];neighbors=[a.GetIdx() for a in m.GetAtomWithIdx(c).GetNeighbors() if a.GetIdx()!=match[5]]
  if m.GetAtomWithIdx(c).GetAtomicNum()==6 and not m.GetAtomWithIdx(c).IsInRing():
   if len(neighbors)==1 and m.GetAtomWithIdx(neighbors[0]).GetAtomicNum() in (7,8):anchors[c]=(5.2,1.5);anchors[neighbors[0]]=(5.2,0)
   if len(neighbors)==3:
    oxy=[i for i in neighbors if m.GetAtomWithIdx(i).GetAtomicNum()==8];methyl=[i for i in neighbors if m.GetAtomWithIdx(i).GetAtomicNum()==6 and m.GetAtomWithIdx(i).GetDegree()==1]
    if len(oxy)==1 and len(methyl)==2:
     anchors[c]=(5.2,1.5);anchors[oxy[0]]=(5.2,0);anchors[methyl[0]]=(5.95,2.8);anchors[methyl[1]]=(6.7,1.5)
 reference=Chem.Mol(m);reference.RemoveAllConformers();reference_conf=Chem.Conformer(m.GetNumAtoms());reference_conf.Set3D(False)
 for idx,(x,y) in anchors.items():reference_conf.SetAtomPosition(idx,Point3D(x,y,0))
 reference.AddConformer(reference_conf)
 rdDepictor.GenerateDepictionMatching2DStructure(m,reference,[(i,i) for i in anchors])
 for b in m.GetBonds():b.SetBondDir(Chem.BondDir.NONE)
 Chem.WedgeMolBonds(m,m.GetConformer())
 center=match[5];bond=m.GetBondBetweenAtoms(center,match[6])
 if m.GetAtomWithIdx(center).GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED:
  for b in m.GetAtomWithIdx(center).GetBonds():
   if b.GetBondDir() in (Chem.BondDir.BEGINWEDGE,Chem.BondDir.BEGINDASH):b.SetBondDir(Chem.BondDir.NONE)
  if bond.GetBeginAtomIdx()!=center:
   other=bond.GetOtherAtomIdx(center);rw=Chem.RWMol(m)
   neighbors={i:[a.GetIdx() for a in m.GetAtomWithIdx(i).GetNeighbors()] for i in (center,other)}
   rw.RemoveBond(center,other);rw.AddBond(center,other,Chem.BondType.SINGLE)
   for idx,old in neighbors.items():
    atom=rw.GetAtomWithIdx(idx);order=[old.index(a.GetIdx()) for a in atom.GetNeighbors()]
    if sum(order[i]>order[j] for i in range(len(order)) for j in range(i+1,len(order)))%2:atom.InvertChirality()
   m.__init__(rw.GetMol());bond=m.GetBondBetweenAtoms(center,other)
  Chem.WedgeBond(bond,center,m.GetConformer())
 if Chem.MolToSmiles(m,isomericSmiles=True)!=before:raise ValueError('Source depiction changed molecular stereochemistry')
 return True
