#!/bin/bash
# chunk_blast.sh <chunk_idx_start> <chunk_idx_end> ; chunks of 100 genomes
export PATH=~/amr/blast/bin:$PATH
cd ~/amr
python3 - "$1" "$2" <<'PYEOF'
import os, sys, subprocess, glob
start,end=int(sys.argv[1]),int(sys.argv[2])
fnas=sorted(glob.glob("data/asm/fna/*.fna"))
chunks=[fnas[i:i+100] for i in range(0,len(fnas),100)]
for ci in range(start,min(end,len(chunks))):
    combined=f"data/asm/chunks/chunk{ci:02d}.fna"
    db=f"data/asm/chunks/chunk{ci:02d}"
    out=f"results/blast/chunk{ci:02d}.tsv"
    if os.path.exists(out): print(ci,"done"); continue
    with open(combined,"w") as o:
        for fp in chunks[ci]:
            acc=os.path.basename(fp)[:-4]
            for line in open(fp):
                if line.startswith(">"): o.write(">"+acc+"|"+line[1:])
                else: o.write(line)
    subprocess.run(["makeblastdb","-in",combined,"-dbtype","nucl","-out",db,"-logfile","/dev/null"],check=True)
    subprocess.run(["blastn","-query","data/query_refs.fasta","-db",db,"-dust","no",
        "-outfmt","6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs",
        "-max_target_seqs","20000","-evalue","1e-50","-out",out],check=True)
    os.remove(combined)
    for ext in (".ndb",".nhr",".nin",".nog",".nos",".not",".nsq",".ntf",".nto"):
        p=db+ext
        if os.path.exists(p): os.remove(p)
    print(ci,"ok",os.path.getsize(out))
PYEOF
