import csv, re, json, glob, os
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq

# species per asm
org={}
for p in ["data/ncbi/ndm_calls.tsv","data/ncbi/oxa48_family_calls.tsv"]:
    for r in csv.DictReader(open(p), delimiter='\t'):
        org[r['asm_acc']]=r['taxgroup_name']
porin_ref={r.id: r for r in SeqIO.parse("data/query_refs.fasta","fasta") if r.id.split("|")[0] in ("kp_ompK35","kp_ompK36","ec_ompC","ec_ompF")}
pref={k.split("|")[0]: str(v.seq) for k,v in porin_ref.items()}
plen={k:len(v) for k,v in pref.items()}

hits=defaultdict(list)
for line in open("results/blast/all_hits.tsv"):
    f=line.rstrip("\n").split("\t")
    q=f[0].split()[0]
    if q not in plen: continue
    acc=f[1].split("|",1)[0]
    hits[acc].append({"q":q,"contig":f[1].split("|",1)[1],"sstart":int(f[8]),"send":int(f[9]),
                      "qstart":int(f[6]),"qend":int(f[7]),"pident":float(f[2]),"qcovs":int(f[12])})

# index fastas lazily
_indexes={}
def get_seq(acc, contig, a, b):
    if acc not in _indexes:
        _indexes[acc]=SeqIO.index(f"data/asm/fna/{acc}.fna","fasta")
    db=_indexes[acc]
    key=None
    for k in db.keys():
        if k.split()[0]==contig or k==contig: key=k; break
    if key is None: return None
    s=str(db[key].seq)
    a,b=max(1,a),min(len(s),b)
    return s[a-1:b]

def porin_status(acc, gene):
    hl=[h for h in hits.get(acc,[]) if h["q"]==gene]
    if not hl: return {"status":"MISSING"}
    h=max(hl, key=lambda x:x["qcovs"])
    L=plen[gene]
    cov=h["qcovs"]
    # extend hit to full gene on subject
    up = h["qstart"]-1; down = L-h["qend"]
    if h["sstart"]<h["send"]:
        a=h["sstart"]-up; b=h["send"]+down
    else:
        a=h["send"]-down; b=h["sstart"]+up
    seq=get_seq(acc, f"{acc}|{h['contig']}" if not h['contig'].startswith(acc) else h['contig'], a, b)
    if seq is None:
        seq=get_seq(acc, h['contig'], a, b)
    res={"status":None,"qcovs":cov,"pident":h["pident"]}
    if seq is None: res["status"]="THIN_NOSEQ"; return res
    if h["sstart"]>h["send"]: seq=str(Seq(seq).reverse_complement())
    if len(seq)%3: seq=seq[:len(seq)-len(seq)%3]
    prot=str(Seq(seq).translate())
    stops=prot[:-1].count("*")
    res["prot_len"]=len(prot); res["internal_stops"]=stops
    ins=(abs(h["send"]-h["sstart"])+1)-(h["qend"]-h["qstart"]+1)
    res["indel_bp"]=ins
    if cov<60: res["status"]="MISSING"
    elif stops>0: res["status"]="LOF_STOP"
    elif cov<90: res["status"]="PARTIAL"
    elif ins>=6: res["status"]="INSERTION"
    elif ins<=-6: res["status"]="DELETION"
    else: res["status"]="INTACT"
    return res

# ---- G4 spike-in positive control: synthetic mutants through the same logic ----
import random
random.seed(12)
ctrl={}
k36=pref["kp_ompK36"]
# intact
ctrl["intact"]=("INTACT", k36)
# premature stop at codon 120
m=list(k36); m[357:360]=list("TAA"); ctrl["stop120"]=("LOF_STOP","".join(m))
# loop-3 insertion: duplicate 12 bp at codon 134
ins=k36[:400]+k36[388:400]+k36[400:]; ctrl["loop3_ins"]=("INSERTION_LIKE", ins)
# big deletion
ctrl["del"]=("PARTIAL", k36[:500])
print("G4 spike-in controls (caller primitives):")
for name,(exp,seq) in ctrl.items():
    stops=str(Seq(seq[:len(seq)-len(seq)%3]).translate())[:-1].count("*")
    indel=len(seq)-len(k36)
    det="INTACT"
    if stops>0: det="LOF_STOP"
    elif len(seq)<0.6*len(k36): det="MISSING/PARTIAL"
    elif indel>=6: det="INSERTION"
    elif indel<=-6: det="DELETION"
    print(f"  {name}: expected~{exp} -> stops={stops} indel={indel} -> {det}")

# ---- run porin status across genomes ----
out={}
for acc in sorted({a for a in hits}):
    o=org.get(acc,"")
    genes=("kp_ompK35","kp_ompK36") if "Klebsiella" in o else ("ec_ompC","ec_ompF") if "E.coli" in o else ()
    if not genes: continue
    out[acc]={"organism":o, **{g:porin_status(acc,g) for g in genes}}
json.dump(out, open("results/porin_status.json","w"))
c=Counter((v["organism"],v[g]["status"]) for v in out.values() for g in ("kp_ompK35","kp_ompK36","ec_ompC","ec_ompF") if g in v)
for k,n in sorted(c.items()): print(k,n)
