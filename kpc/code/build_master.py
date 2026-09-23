#!/usr/bin/env python3
"""Master table: 296-allele union reference, Ambler annotation, structure, outbreak census."""
import sys, re, math, csv
sys.path.insert(0,"/home/sandbox/amr/code")
from collections import defaultdict
from Bio import SeqIO
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner
from kpc_caller import signature, REF

DATA="/home/sandbox/amr/data"; RES="/home/sandbox/amr/results"

AA = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(f"{DATA}/kpc_ref_union_aa.fasta","fasta")}

# Ambler mapping via 3DW0
p = MMCIFParser(QUIET=True)
st = p.get_structure("k", f"{DATA}/structures/3DW0.cif")
ch = list(next(iter(st)).get_chains())[0]
pdb_res = {r.id[1]: r for r in ch if r.id[0]==" "}
ambler_pos = list(pdb_res.keys())
pdb_seq = "".join(seq1(r.get_resname()) for r in pdb_res.values())
al = PairwiseAligner(); al.mode="global"; al.match_score=2; al.mismatch_score=-1
al.open_gap_score=-10; al.extend_gap_score=-0.5
aln = al.align(REF, pdb_seq)[0]
prot2ambler = {}
for (rs,re_),(ps,pe) in zip(*aln.aligned):
    for i in range(re_-rs): prot2ambler[rs+i+1]=ambler_pos[ps+i]
s70_ca = pdb_res[70]["CA"].coord
def d70(amb):
    r = pdb_res.get(amb)
    if not r or "CA" not in r: return None
    v = r["CA"].coord - s70_ca
    return round(float((v*v).sum()**0.5),1)
REGIONS=[("omega_loop",164,179),("loop237_243",237,243),("loop266_275",266,275)]
def regs_of(a): return [n for n,x,y in REGIONS if a and x<=a<=y]

# literature resistance (review Table 1)
resistant={}
for line in open(f"{DATA}/cazavi_table1_raw.md"):
    line=line.strip()
    if line.startswith("KPC-"):
        c=[x.strip() for x in line.split("|")]
        resistant[c[0].replace("*","").replace("\\","")]=" ".join(c[1:])

# outbreak census
census=defaultdict(lambda: dict(n=0,years=[],geos=defaultdict(int)))
for line in open(f"{DATA}/kpc_isolate_accessions.tsv"):
    f=line.rstrip("\n").split("\t")
    if len(f)<7: continue
    a=f[5]
    if not a.startswith("blaKPC-"): continue
    c=census[a]; c["n"]+=1
    y=f[3][:4]
    if y.isdigit() and 1990<=int(y)<=2026: c["years"].append(int(y))
    g=f[4] if f[4] and f[4]!="NULL" else "unknown"
    c["geos"][g.split(":")[0]]+=1

rows=[]
for name in sorted(AA, key=lambda x:int(x.split("-")[1])):
    sig=signature(AA[name])
    ops=[]
    for op in sig:
        m=re.match(r'^([A-Z])(\d+)([A-Z])$',op)
        if m:
            wt,pos,mut=m.group(1),int(m.group(2)),m.group(3)
            amb=prot2ambler.get(pos)
            ops.append((amb,f"{wt}{amb}{mut}" if amb else op,regs_of(amb),d70(amb) if amb else None))
        else:
            m=re.match(r'^(ins|del)(\d+)(?:-(\d+))?_(.*)$',op)
            if m:
                kind,a,b,seg=m.group(1),int(m.group(2)),m.group(3),m.group(4)
                b=int(b) if b else a
                ambs=[prot2ambler.get(pp) for pp in range(a,b+1) if prot2ambler.get(pp)]
                amb=ambs[0] if ambs else None
                rr=sorted({r for x in ambs for r in regs_of(x)})
                dd=min([d70(x) for x in ambs if d70(x)],default=None)
                ops.append((amb, f"{kind}{amb}_{seg}" if amb else op, rr, dd))
    sig_amb="; ".join(o[1] for o in ops) if ops else "identical to KPC-2"
    regs=sorted({r for o in ops for r in o[2]})
    dmin=min([o[3] for o in ops if o[3] is not None],default=None)
    known=name in resistant
    key="bla"+name
    c=census.get(key, dict(n=0,years=[],geos={}))
    top=max(c["geos"].items(),key=lambda kv:kv[1]) if c["geos"] else ("",0)
    if known: grade="documented CAZ-AVI resistance (literature)"
    elif regs or (dmin is not None and dmin<=10.0): grade="hotspot/active-site proximal - consistent with resistance (hypothesis)"
    else: grade="no resistance evidence"
    rows.append(dict(
        allele=name, signature_ambler=sig_amb, n_ops=len(ops),
        hotspot_regions=";".join(regs), min_dist_S70_A=dmin if dmin is not None else "",
        cazavi_resistant_lit=("yes" if known else "no"), lit_mutations=resistant.get(name,""),
        grade=grade, n_isolates=c["n"],
        year_first=min(c["years"]) if c["years"] else "",
        year_last=max(c["years"]) if c["years"] else "",
        n_2020plus=sum(1 for y in c["years"] if y>=2020),
        top_country=top[0], n_top_country=top[1], n_countries=len(c["geos"])))

with open(f"{RES}/master_allele_table.tsv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(rows)

obs=[r for r in rows if r["n_isolates"]>0]
res_obs=[r for r in obs if r["cazavi_resistant_lit"]=="yes"]
print("alleles in union:", len(rows))
print("observed in outbreaks:", len(obs))
print("observed + documented CAZ-AVI resistant:", len(res_obs), "isolates:", sum(r['n_isolates'] for r in res_obs))
print("observed + hotspot-proximal hypothesis:", sum(1 for r in obs if r['grade'].startswith('hotspot')))
print("\ntop resistant alleles by isolates:")
for r in sorted(res_obs,key=lambda x:-x['n_isolates'])[:12]:
    print(f"  {r['allele']:9s} n={r['n_isolates']:5d} {r['signature_ambler'][:60]:60s} {r['top_country']} first={r['year_first']} last={r['year_last']}")
