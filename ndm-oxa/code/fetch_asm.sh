#!/bin/bash
# usage: fetch_asm.sh <start_batch> <end_batch_exclusive>  (batches of 50)
cd ~/amr/data/asm
python3 - "$1" "$2" <<'PYEOF'
import json, sys, os, urllib.request, time, zipfile, shutil
start, end = int(sys.argv[1]), int(sys.argv[2])
targets = json.load(open("../ncbi/download_targets.json"))
B = 50
batches = [targets[i:i+B] for i in range(0, len(targets), B)]
os.makedirs("zips", exist_ok=True); os.makedirs("fna", exist_ok=True)
for bi in range(start, min(end, len(batches))):
    accs = batches[bi]
    zp = f"zips/batch{bi:03d}.zip"
    if os.path.exists(zp) and zipfile.is_zipfile(zp):
        print(f"batch {bi}: exists, skip"); continue
    body = json.dumps({"accessions": accs, "include_annotation_type": ["GENOME_FASTA"]}).encode()
    ok = False
    for t in range(3):
        try:
            req = urllib.request.Request("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/download",
                                         data=body, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=110) as r, open(zp, "wb") as f:
                shutil.copyfileobj(r, f)
            ok = zipfile.is_zipfile(zp); break
        except Exception as e:
            print(f"batch {bi} try {t}: {e}"); time.sleep(3+3*t)
    print(f"batch {bi}: {'ok' if ok else 'FAILED'} {os.path.getsize(zp) if os.path.exists(zp) else 0}")
    # extract fna members
    if ok:
        with zipfile.ZipFile(zp) as z:
            for n in z.namelist():
                if n.endswith("_genomic.fna"):
                    acc = n.split("/")[-2]
                    with z.open(n) as src, open(f"fna/{acc}.fna", "wb") as dst:
                        shutil.copyfileobj(src, dst)
    time.sleep(1)
PYEOF
