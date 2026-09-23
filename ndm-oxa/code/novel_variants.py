import csv, re, json, glob, os
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

def aname(h):
    m=re.search(r"ARO:\d+\|([^ |]+)", h); return m.group(1)
# CARD reference proteins (translate nuc refs)
refprot={}
for r in SeqIO.parse("data/card/all_ref.fasta","fasta"):
    refprot[aname(r.description)]=str(Seq(str(r.seq)).translate())
idx={os.path.basename(p)[:-4]: p for p in glob.glob("data/asm/fna/*.fna")}

hits=defaultdict(list)
for line in open("results/blast/all_hits.tsv"):
    f=line.rstrip("\n").split("\t")
    q=f[0].split()[0]
    if "ARO:" not in q and not q.startswith(("NDM","OXA","gb|")):
        continue
    qn=aname(q) if "ARO:" in q else q
    if qn not in refprot: continue
    acc=f[1].split("|",1)[0]
    hits[acc].append({"q":qn,"pident":float(f[2]),"qstart":int(f[6]),"qend":int(f[7]),
                      "contig":f[1],"sstart":int(f[8]),"send":int(f[9]),"bitscore":float(f[11])})

# novel targets = my NONEXACT calls
my=json.load(open("results/my_calls.json"))
_indexes={}
def get_seq(acc, contig, a, b):
    if acc not in _indexes: _indexes[acc]=SeqIO.index(f"data/asm/fna/{acc}.fna","fasta")
    db=_indexes[acc]
    for k in db.keys():
        kk=k.split()[0]
        if kk==contig or k==contig or kk==contig.split("|",1)[-1] or k==contig.split("|",1)[-1]:
            s=str(db[k].seq); a,b=max(1,a),min(len(s),b)
            return s[a-1:b]
    return None

al=PairwiseAligner(); al.mode='global'; al.match_score=1; al.mismatch_score=0; al.open_gap_score=-1; al.extend_gap_score=-0.5
novel={}
for acc,fams in my.items():
    for fam,calls in fams.items():
        for call in calls:
            if not call.endswith("*NONEXACT"): continue
            base=call.replace("*NONEXACT","")
            cand=[h for h in hits.get(acc,[]) if h["q"]==base]
            if not cand: continue
            h=max(cand,key=lambda x:x["bitscore"])
            L=len(refprot[base])*3
            up=h["qstart"]-1; down=L-h["qend"]
            if h["sstart"]<h["send"]: a,b=h["sstart"]-up,h["send"]+down
            else: a,b=h["send"]-down,h["sstart"]+up
            seq=get_seq(acc,h["contig"],a,b)
            if not seq: continue
            if h["sstart"]>h["send"]: seq=str(Seq(seq).reverse_complement())
            prot=str(Seq(seq[:len(seq)-len(seq)%3]).translate())
            rp=refprot[base]
            aln=al.align(rp,prot)[0]
            ra,pa=str(aln[0]),str(aln[1])
            subs=[]; rp_i=0; pp_i=0
            for x,y in zip(ra,pa):
                if x!='-': rp_i+=1
                if y!='-': pp_i+=1
                if x!='-' and y!='-' and x!=y: subs.append(f"{x}{rp_i}{y}")
                elif x=='-' and y!='-': subs.append(f"ins{rp_i}{y}")
                elif x!='-' and y=='-': subs.append(f"del{rp_i}{x}")
            stops=prot[:-1].count("*")
            key=(acc,fam,base)
            novel[key]={"organism_call":call,"closest_ref":base,"subs":subs[:30],"n_subs":len(subs),
                        "internal_stops":stops,"prot_len":len(prot)}
print("novel characterized:", len(novel))
json.dump({f"{a}|{f}|{b}":v for (a,f,b),v in novel.items()}, open("results/novel_variants.json","w"), indent=1)
pat=Counter()
for v in novel.values():
    if v["internal_stops"]>0: pat["INTERNAL_STOP"]+=1
    elif v["n_subs"]==0: pat["synonymous_only"]+=1
    else: pat["aa_substitutions"]+=1
print(pat)
for k,v in list(novel.items())[:10]: print(k, v["subs"][:8], "stops:",v["internal_stops"])
