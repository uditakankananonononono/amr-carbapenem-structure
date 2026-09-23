import re, json, hashlib
from Bio import SeqIO
from Bio.Seq import Seq

CARD_NUC="data/card/nucleotide_fasta_protein_homolog_model.fasta"
CARD_PROT="data/card/protein_fasta_protein_homolog_model.fasta"

def name_of(h):
    m=re.search(r"\|ARO:\d+\|([^|]+)\[", h+"[")
    m2=re.search(r"ARO:\d+\|(.+?) \[", h)
    return m2.group(1) if m2 else h

# --- NDM set ---
ndm_nuc={}
for r in SeqIO.parse(CARD_NUC,"fasta"):
    nm=name_of(r.description)
    if re.fullmatch(r"NDM-\d+", nm):
        ndm_nuc[nm]=r
print("NDM nuc alleles:", len(ndm_nuc), sorted(ndm_nuc, key=lambda x:int(x.split('-')[1]))[:5], "...")

# --- OXA proteins: cluster vs OXA-48 ---
oxa_prot={}
for r in SeqIO.parse(CARD_PROT,"fasta"):
    nm=name_of(r.description)
    if re.fullmatch(r"OXA-\d+[a-z]?", nm):
        oxa_prot[nm]=str(r.seq)
print("OXA proteins total:", len(oxa_prot))
ref48=oxa_prot.get("OXA-48")
print("OXA-48 found:", ref48 is not None, "len:", len(ref48) if ref48 else None)

# quick identity via simple global alignment-free proxy won't work; use pairwise
from Bio.Align import PairwiseAligner
al=PairwiseAligner()
al.mode='global'; al.match_score=1; al.mismatch_score=0
al.open_gap_score=-0.5; al.extend_gap_score=-0.1
fam=[]
for nm,seq in oxa_prot.items():
    if not seq or not ref48: continue
    a=al.align(ref48,seq)[0]
    ident=a.score/max(len(ref48),len(seq))
    if ident>=0.90:
        fam.append((nm,round(ident,4)))
fam.sort(key=lambda x:(int(re.search(r"\d+",x[0]).group()),x[0]))
print("OXA-48-like family (>=90% aa identity):", len(fam))
for nm,i in fam: print(" ", nm, i)
with open("data/card/oxa48_family.json","w") as f:
    json.dump({nm:i for nm,i in fam}, f, indent=1)

# --- write subset fastas ---
with open("data/card/ndm_ref.fasta","w") as f:
    for nm in sorted(ndm_nuc, key=lambda x:int(x.split('-')[1])):
        SeqIO.write(ndm_nuc[nm], f, "fasta")
oxa48_fam_nuc=[]
for r in SeqIO.parse(CARD_NUC,"fasta"):
    nm=name_of(r.description)
    if nm in {x[0] for x in fam}:
        oxa48_fam_nuc.append(r)
with open("data/card/oxa48like_ref.fasta","w") as f:
    SeqIO.write(oxa48_fam_nuc, f, "fasta")
print("OXA-48-like nuc alleles written:", len(oxa48_fam_nuc))
for p in ["data/card/ndm_ref.fasta","data/card/oxa48like_ref.fasta","data/card/oxa48_family.json"]:
    print(p, hashlib.sha256(open(p,'rb').read()).hexdigest()[:16])
