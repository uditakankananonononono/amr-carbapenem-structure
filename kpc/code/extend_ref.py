#!/usr/bin/env python3
"""Extended Gate G1 over current AMRFinderPlus db (2026-08-07.1) + combined set."""
import sys
sys.path.insert(0,"/home/sandbox/amr/code")
from Bio import SeqIO
from Bio.Seq import Seq
from kpc_caller import signature, AA as CARD_AA, NT as CARD_NT, REF

DATA="/home/sandbox/amr/data"
def load(path):
    d={}
    for r in SeqIO.parse(path,"fasta"):
        import re
        m = re.search(r'(KPC-\d+)', r.description)
        if m: d.setdefault(m.group(1), str(r.seq))
    return d
amf_aa = load(f"{DATA}/AMRProt.fa")
amf_nt = load(f"{DATA}/AMR_CDS.fa")
# remove trailing stops
amf_aa = {k: v.rstrip("*") for k,v in amf_aa.items()}
amf_nt = {k: v for k,v in amf_nt.items() if k in amf_aa}
print("AMRFinder 2026 KPC alleles with nt+aa:", len(amf_aa))

n_trans = n_sig = n_lookup = 0; fails=[]
for name, prot in amf_aa.items():
    nt = amf_nt[name]
    t = str(Seq(nt).translate(table=11)).rstrip("*")
    if t != prot:
        fails.append((name, f"trans {len(t)} vs {len(prot)}")); continue
    n_trans += 1
    if signature(prot) == signature(t): n_sig += 1
    else: fails.append((name,"sig"))
    hits = [n2 for n2,s2 in amf_nt.items() if s2.upper()==nt.upper()]
    if name in hits: n_lookup += 1
    else: fails.append((name, f"lookup {hits[:3]}"))
tot=len(amf_aa)
print(f"G1-ext translation: {n_trans}/{tot}  signature: {n_sig}/{tot}  nt-lookup: {n_lookup}/{tot}")
for f_ in fails[:10]: print("  fail:", f_)
print("GATE G1-extended:", "PASS (100%)" if (n_trans==tot and n_sig==tot and n_lookup==tot) else "FAIL")

# combined union
union_aa = dict(amf_aa); 
for k,v in CARD_AA.items(): union_aa.setdefault(k,v)
union_nt = dict(amf_nt)
for k,v in CARD_NT.items(): union_nt.setdefault(k,v)
print("union reference set:", len(union_aa), "alleles")
# cross-check: alleles in both - same protein?
both = set(amf_aa) & set(CARD_AA)
diff = [k for k in both if amf_aa[k].rstrip("*") != CARD_AA[k].rstrip("*")]
print("in both:", len(both), "protein disagreements:", len(diff), diff[:5])
only_card = sorted(set(CARD_AA)-set(amf_aa), key=lambda x:int(x.split('-')[1]))
print("CARD-only:", len(only_card), only_card[:10])
with open(f"{DATA}/kpc_ref_union_aa.fasta","w") as f:
    for k in sorted(union_aa, key=lambda x:int(x.split('-')[1])):
        f.write(f">{k}\n{union_aa[k].rstrip('*')}\n")
with open(f"{DATA}/kpc_ref_union_nt.fasta","w") as f:
    for k in sorted(union_nt, key=lambda x:int(x.split('-')[1])):
        f.write(f">{k}\n{union_nt[k]}\n")
