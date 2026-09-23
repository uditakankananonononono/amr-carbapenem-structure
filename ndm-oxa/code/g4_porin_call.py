import subprocess, sys, json
from Bio import SeqIO
from Bio.Seq import Seq

asm=sys.argv[1]
# blast refs vs contigs
subprocess.run(["./blast/bin/makeblastdb","-in",asm,"-dbtype","nucl","-out","asmdb","-logfile","/dev/null"],check=True)
cat_fa="refs.fasta"
with open(cat_fa,"w") as o:
    for r in SeqIO.parse("ompK35.fasta","fasta"): SeqIO.write(r,o,"fasta")
    for r in SeqIO.parse("ompK36.fasta","fasta"): SeqIO.write(r,o,"fasta")
subprocess.run(["./blast/bin/blastn","-query",cat_fa,"-db","asmdb","-dust","no","-max_target_seqs","20000",
  "-outfmt","6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovs","-evalue","1e-50","-out","hits.tsv"],check=True)
plen={r.id.split()[0]: len(r.seq) for r in SeqIO.parse(cat_fa,"fasta")}
db=SeqIO.index(asm,"fasta")
out={}
for gene,L in plen.items():
    best=None
    for line in open("hits.tsv"):
        f=line.split("\t")
        if f[0].split()[0]!=gene: continue
        h={"contig":f[1],"pident":float(f[2]),"qstart":int(f[4]),"qend":int(f[5]),
           "sstart":int(f[6]),"send":int(f[7]),"qcovs":int(f[10]),"bitscore":float(f[9])}
        if best is None or h["qcovs"]>best["qcovs"] or (h["qcovs"]==best["qcovs"] and h["bitscore"]>best["bitscore"]):
            best=h
    if best is None:
        out[gene]={"status":"MISSING"}; continue
    up=best["qstart"]-1; down=L-best["qend"]
    if best["sstart"]<best["send"]: a,b=best["sstart"]-up,best["send"]+down
    else: a,b=best["send"]-down,best["sstart"]+up
    key=best["contig"]
    k=[k for k in db.keys() if k.split()[0]==key or k==key]
    if not k: out[gene]={"status":"THIN_NOSEQ"}; continue
    s=str(db[k[0]].seq)
    a,b=max(1,a),min(len(s),b)
    seq=s[a-1:b]
    if best["sstart"]>best["send"]: seq=str(Seq(seq).reverse_complement())
    if len(seq)%3: seq=seq[:len(seq)-len(seq)%3]
    prot=str(Seq(seq).translate())
    stops=prot[:-1].count("*")
    ins=(abs(best["send"]-best["sstart"])+1)-(best["qend"]-best["qstart"]+1)
    cov=best["qcovs"]
    if cov<60: st="MISSING"
    elif stops>0: st="LOF_STOP"
    elif cov<90: st="PARTIAL"
    elif ins>=6: st="INSERTION"
    elif ins<=-6: st="DELETION"
    else: st="INTACT"
    out[gene]={"status":st,"qcovs":cov,"pident":best["pident"],"prot_len":len(prot)-1,"internal_stops":stops,"indel_bp":ins,"stop_positions":[i for i,c in enumerate(prot) if c=='*' and i<len(prot)-1][:5]}
print(json.dumps(out, indent=1))
