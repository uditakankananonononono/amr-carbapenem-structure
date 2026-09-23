import re, json, math
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

def pdb_chain_seq(path, chain=None):
    seq={}
    resn3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M"}
    atoms=[]  # (resno, chain, atomname, x,y,z, resname)
    het=[]
    for line in open(path):
        if line.startswith("ATOM"):
            ch=line[21]
            if chain and ch!=chain: continue
            rn=int(line[22:26]); an=line[12:16].strip(); res=line[17:20].strip()
            x,y,z=float(line[30:38]),float(line[38:46]),float(line[46:54])
            atoms.append((rn,ch,an,x,y,z,res))
            if an=="CA": seq[rn]=resn3.get(res,"X")
        elif line.startswith("HETATM"):
            res=line[17:20].strip(); an=line[12:16].strip()
            if res in ("ZN",) or (res not in ("HOH","ZN","NA","CL","MG","CA","SO4","PO4","GOL","EDO","ACT","FMT") and an):
                het.append((res,an,float(line[30:38]),float(line[38:46]),float(line[46:54])))
    sseq="".join(seq[k] for k in sorted(seq))
    resnos=sorted(seq)
    return sseq,resnos,atoms,het

def aname(h):
    m=re.search(r"ARO:\d+\|([^ |]+)", h); return m.group(1)
refprot={}
for r in SeqIO.parse("data/card/all_ref.fasta","fasta"):
    refprot[aname(r.description)]=str(Seq(str(r.seq)).translate())

al=PairwiseAligner(); al.mode='global'; al.match_score=1; al.mismatch_score=0; al.open_gap_score=-2; al.extend_gap_score=-0.5

def build_map(ref_prot, pdb_seq, pdb_resnos):
    a=al.align(ref_prot, pdb_seq)[0]
    ra,pa=str(a[0]),str(a[1])
    m={}; ri=0; pi=0
    for x,y in zip(ra,pa):
        if x!='-': ri+=1
        if y!='-': pi+=1
        if x!='-' and y!='-':
            m[ri]=pdb_resnos[pi-1]
    return m

def dist_to_features(resno_pdb, atoms, het, want_het=("ZN",)):
    d_zn=999; d_lig=999
    ra=[a for a in atoms if a[0]==resno_pdb]
    if not ra: return d_zn,d_lig
    for res,an,x,y,z in het:
        for _,_,_,ax,ay,az,_ in ra:
            d=math.dist((x,y,z),(ax,ay,az))
            if res in want_het: d_zn=min(d_zn,d)
            elif res not in ("HOH",): d_lig=min(d_lig,d)
    return d_zn,d_lig

def zone(dz,dl,catdist):
    if dz<=4.0: return "zinc_first_shell"
    if dz<=10.0: return "active_site_near"
    if catdist<=8.0: return "catalytic_near"
    return "surface_other"

results={}

# ---- NDM ----
sseq,resnos,atoms,het=build=None,None,None,None
sseq,resnos,atoms,het = pdb_chain_seq("data/pdb/4EYL.pdb", chain="A")
ndm_map=build_map(refprot["NDM-1"], sseq, resnos)
print("NDM-1 map coverage:", len(ndm_map), "of", len(refprot["NDM-1"]))
# catalytic reference atoms: Zn ions from HETATM ZN
ndm_variants=defaultdict(set)
for r in SeqIO.parse("data/card/all_ref.fasta","fasta"):
    nm=aname(r.description)
    if not nm.startswith("NDM-"): continue
    p=str(Seq(str(r.seq)).translate())
    rp=refprot["NDM-1"]
    a=al.align(rp,p)[0]; ra,pa=str(a[0]),str(a[1]); ri=0
    for x,y in zip(ra,pa):
        if x!='-': ri+=1
        if x!=y and x!='-' and y!='-': ndm_variants[ri].add((nm,f"{x}{ri}{y}"))
# novel
nov=json.load(open("results/novel_variants.json"))
nov_ndm={}
for k,v in nov.items():
    if v["closest_ref"].startswith("NDM") and v["subs"]:
        for s in v["subs"]:
            m2=re.match(r"([A-Z])(\d+)([A-Z])",s)
            if m2: nov_ndm.setdefault(int(m2.group(2)),[]).append((k.split("|")[0],v["closest_ref"],s))
ndm_table=[]
for pos in sorted(set(ndm_variants)|set(nov_ndm)):
    pdbres=ndm_map.get(pos)
    if pdbres is None:
        ndm_table.append({"pos":pos,"pdb":None,"zone":"NOT_IN_STRUCTURE(signal/uncovered)","n_alleles":len(ndm_variants.get(pos,[])),"novel":nov_ndm.get(pos,[])})
        continue
    dz,dl=dist_to_features(pdbres,atoms,het)
    ndm_table.append({"pos":pos,"pdb":pdbres,"d_zn":round(dz,1),"d_lig":round(dl,1),"zone":zone(dz,dl,dl),
                      "alleles":sorted(ndm_variants.get(pos,[])),"novel":nov_ndm.get(pos,[])})
json.dump(ndm_table, open("results/ndm_structure_table.json","w"), indent=1)
zc=defaultdict(int)
for t in ndm_table: zc[t["zone"]]+=1
print("NDM variant positions:", len(ndm_table), dict(zc))

# ---- OXA-48 ----
sseq,resnos,atoms,het = pdb_chain_seq("data/pdb/6P97.pdb", chain="A")
oxa_map=build_map(refprot["OXA-48"], sseq, resnos)
print("OXA-48 map coverage:", len(oxa_map), "of", len(refprot["OXA-48"]))
oxa_variants=defaultdict(set)
for r in SeqIO.parse("data/card/all_ref.fasta","fasta"):
    nm=aname(r.description)
    if not nm.startswith("OXA-"): continue
    p=str(Seq(str(r.seq)).translate())
    rp=refprot["OXA-48"]
    a=al.align(rp,p)[0]; ra,pa=str(a[0]),str(a[1]); ri=0
    for x,y in zip(ra,pa):
        if x!='-': ri+=1
        if x!=y and x!='-' and y!='-': oxa_variants[ri].add((nm,f"{x}{ri}{y}"))
oxa_table=[]
for pos in sorted(oxa_variants):
    pdbres=oxa_map.get(pos)
    if pdbres is None:
        oxa_table.append({"pos":pos,"pdb":None,"zone":"NOT_IN_STRUCTURE","n_alleles":len(oxa_variants[pos])}); continue
    # distances: to imipenem ligand (non-Zn het), to catalytic S70/K73 (mapped)
    s70=oxa_map.get(70); k73=oxa_map.get(73)
    def dres(a,b):
        A=[x for x in atoms if x[0]==a]; B=[x for x in atoms if x[0]==b]
        return min((math.dist((x[3],x[4],x[5]),(y[3],y[4],y[5])) for x in A for y in B), default=999)
    _,dl=dist_to_features(pdbres,atoms,het,want_het=())
    ds70=dres(pdbres,s70); dk73=dres(pdbres,k73)
    cat=min(ds70,dk73)
    oxa_table.append({"pos":pos,"pdb":pdbres,"d_S70":round(ds70,1),"d_K73":round(dk73,1),"d_lig":round(dl,1),
                      "zone":("catalytic_near" if cat<=8 else ("ligand_near" if dl<=8 else "surface_other")),
                      "n_alleles":len(oxa_variants[pos]),"alleles":sorted(oxa_variants[pos])})
json.dump(oxa_table, open("results/oxa48_structure_table.json","w"), indent=1)
zc=defaultdict(int)
for t in oxa_table: zc[t["zone"]]+=1
print("OXA-48 variant positions:", len(oxa_table), dict(zc))
