import csv, re, json, glob, os
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def aname(h):
    m=re.search(r"ARO:\d+\|([^ |]+)", h); return m.group(1)
reflen={aname(r.description): len(r.seq) for r in SeqIO.parse("data/card/all_ref.fasta","fasta")}

# index genome fastas: acc -> path
idx={os.path.basename(p)[:-4]: p for p in glob.glob("data/asm/fna/*.fna")}

# parse hits
hits=defaultdict(list)
for line in open("results/blast/all_hits.tsv"):
    f=line.rstrip("\n").split("\t")
    q,s=f[0],f[1]
    acc=s.split("|",1)[0]
    qname = aname(q) if "ARO:" in q else q.split()[0].replace(">","")
    hits[acc].append({"q":qname,"s":s,"pident":float(f[2]),"length":int(f[3]),"mismatch":int(f[4]),
        "gapopen":int(f[5]),"qstart":int(f[6]),"qend":int(f[7]),"sstart":int(f[8]),"send":int(f[9]),
        "evalue":float(f[10]),"bitscore":float(f[11]),"qcovs":int(f[12])})

NDM=set(reflen)  # all CARD alleles in ref db
porins={"kp_ompK35","kp_ompK36","ec_ompC","ec_ompF"}

def call_alleles(hl):
    """best call per genome for carbapenemase genes"""
    bla=[h for h in hl if h["q"] in NDM]
    by={}
    for h in bla:
        fam="NDM" if h["q"].startswith("NDM") else "OXA48"
        cur=by.get(fam)
        if cur is None or (h["bitscore"],h["pident"],h["qcovs"])>(cur["bitscore"],cur["pident"],cur["qcovs"]):
            by[fam]=h
    calls={}
    for fam,h in by.items():
        L=reflen[h["q"]]
        exact = h["pident"]==100.0 and h["qstart"]==1 and h["qend"]==L
        calls[fam]={"allele":h["q"] if exact else h["q"]+"*NOVEL", "exact":exact, "hit":h}
    return calls

# NCBI calls for validation set
ncbi={}
for p in ["data/ncbi/ndm_calls.tsv","data/ncbi/oxa48_family_calls.tsv"]:
    for r in csv.DictReader(open(p), delimiter='\t'):
        ncbi.setdefault(r['asm_acc'],[]).append(r)

agree=disagree=only_ncbi=only_mine=0
details=[]
my_calls={}
for acc in sorted(idx):
    hl=hits.get(acc,[])
    c=call_alleles(hl)
    my_calls[acc]=c
    nb=ncbi.get(acc,[])
    for r in nb:
        sym=r['element_symbol']
        fam="NDM" if sym.startswith("blaNDM") else "OXA48"
        if fam not in c: only_ncbi+=1; details.append((acc,sym,"NCBI_ONLY")); continue
        my=c[fam]
        exp=sym.replace("bla","")
        if my["exact"] and my["allele"]==exp: agree+=1
        elif not my["exact"] and (r['pct_ref_identity']!='100.0' or r['pct_ref_coverage']!='100.0'): agree+=1  # both non-exact
        else:
            disagree+=1; details.append((acc,sym,f"mine={my['allele']} exact={my['exact']} ncbi={r['amr_method']} {r['pct_ref_identity']}/{r['pct_ref_coverage']}"))
n_with=sum(1 for a,c in my_calls.items() if c)
tot=agree+disagree+only_ncbi
print(f"G2: genomes with calls {n_with}/{len(idx)}; agree {agree}, disagree {disagree}, NCBI-only {only_ncbi} of {tot}")
if tot: print(f"G2 agreement: {100*agree/tot:.2f}%")
for d in details[:12]: print("  ", d)
json.dump({a:{f:{k:(v if k!='hit' else None) for k,v in c2.items()} for f,c2 in c.items()} for a,c in my_calls.items()},
          open("results/my_calls.json","w"))
