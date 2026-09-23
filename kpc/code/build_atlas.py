#!/usr/bin/env python3
"""KPC outbreak mutation atlas + structural effect table (milestone 4 core)."""
import sys, re, math, csv
sys.path.insert(0, "/home/sandbox/amr/code")
from collections import defaultdict
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner
from kpc_caller import signature, AA, NT, REF

DATA = "/home/sandbox/amr/data"

# ---------- Ambler mapping: align CARD KPC-2 protein to 3DW0 (Ambler-numbered) ----------
p = MMCIFParser(QUIET=True)
st = p.get_structure("k", f"{DATA}/structures/3DW0.cif")
model = next(iter(st))
chains = list(model.get_chains())
ch = chains[0]
pdb_res = {}   # ambler pos -> residue
for r in ch:
    if r.id[0] == " ":
        pdb_res[r.id[1]] = r
pdb_seq = "".join(seq1(r.get_resname()) for r in pdb_res.values())
ambler_pos = list(pdb_res.keys())

al = PairwiseAligner(); al.mode="global"; al.match_score=2; al.mismatch_score=-1
al.open_gap_score=-10; al.extend_gap_score=-0.5
aln = al.align(REF, pdb_seq)[0]
(rb, pb) = aln.aligned
prot2ambler = {}
for (rs, re_), (ps, pe) in zip(rb, pb):
    for i in range(re_-rs):
        prot2ambler[rs+i+1] = ambler_pos[ps+i]
# sanity
assert REF[prot2ambler and 0] is not None
s70_prot = [k for k,v in prot2ambler.items() if v==70][0]
print("KPC-2 protein pos of Ambler Ser70:", s70_prot, "residue:", REF[s70_prot-1])
assert REF[s70_prot-1] == "S"

# CA coords & distance to Ser70
s70_ca = pdb_res[70]["CA"].coord
def dist_to_s70(amb):
    r = pdb_res.get(amb)
    if r is None or "CA" not in r: return None
    d = r["CA"].coord - s70_ca
    return round(math.sqrt((d*d).sum()), 1)

REGIONS = [("omega_loop", 164, 179), ("loop237_243", 237, 243), ("loop266_275", 266, 275)]
def regions_of(amb):
    out = [n for n,a,b in REGIONS if a <= amb <= b]
    return ";".join(out) if out else ""

# ---------- literature ground truth: CAZ-AVI-resistant variants (review Table 1) ----------
resistant = {}
for line in open(f"{DATA}/cazavi_table1_raw.md"):
    line = line.strip()
    if not line.startswith("KPC-"): continue
    cells = [c.strip() for c in line.split("|")]
    name = cells[0].replace("*","").replace("\\","")
    muts = " ".join(cells[1:])
    resistant[name] = muts
print("literature CAZ-AVI-resistant variants parsed:", len(resistant))

# ---------- per-allele annotation ----------
def ops_annotated(sig):
    out = []
    for op in sig:
        m = re.match(r'^([A-Z])(\d+)([A-Z])$', op)
        if m:  # substitution
            wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
            amb = prot2ambler.get(pos)
            out.append(dict(op=op, kind="sub", prot_pos=pos, ambler=amb,
                            ambler_op=f"{wt}{amb}{mut}" if amb else op,
                            region=regions_of(amb) if amb else "", dist=dist_to_s70(amb) if amb else None))
        else:
            m = re.match(r'^(ins|del)(\d+)(?:-(\d+))?_(.*)$', op)
            if m:
                kind, a, b, seg = m.group(1), int(m.group(2)), m.group(3), m.group(4)
                b = int(b) if b else a
                ambs = [prot2ambler.get(pp) for pp in range(a, b+1)]
                ambs = [x for x in ambs if x]
                amb = ambs[0] if ambs else None
                amb_end = ambs[-1] if ambs else None
                regs = ";".join(sorted({r for x in ambs for r in regions_of(x).split(";") if r}))
                d = min([dist_to_s70(x) for x in ambs if dist_to_s70(x)], default=None)
                out.append(dict(op=op, kind=kind, prot_pos=a, ambler=amb, ambler_end=amb_end,
                                ambler_op=f"{kind}{amb}-{amb_end}_{seg}" if amb_end and amb_end!=amb else (f"{kind}{amb}_{seg}" if amb else op),
                                region=regs, dist=d))
    return out

rows = []
for name in sorted(AA, key=lambda x: int(x.split("-")[1])):
    sig = signature(AA[name])
    ops = ops_annotated(sig)
    sig_amb = "; ".join(o["ambler_op"] for o in ops) if ops else "WT (identical to KPC-2)"
    regs = ";".join(sorted({o["region"] for o in ops if o["region"]}))
    dmin = min([o["dist"] for o in ops if o["dist"] is not None], default=None)
    known = name in resistant
    # grading per G3
    if known:
        grade = "resistance documented (lit)"
    elif any(o["region"] for o in ops) or any((o["dist"] is not None and o["dist"] <= 10.0) for o in ops):
        grade = "consistent with resistance (hotspot/active-site proximity) - hypothesis"
    else:
        grade = "uncertain - no resistance evidence"
    rows.append(dict(allele=name, n_ops=len(ops), signature_ambler=sig_amb, regions=regs,
                     min_dist_to_s70_A=dmin, cazavi_resistant_lit=("yes" if known else "no"),
                     lit_muts=resistant.get(name,""), grade=grade))

with open(f"{DATA}/kpc_structural_effect_table.tsv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
print("structural effect table:", len(rows), "alleles")

# ---------- join with outbreak census ----------
census = defaultdict(lambda: dict(n=0, years=[], geos=defaultdict(int)))
for line in open(f"{DATA}/kpc_isolate_accessions.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 7: continue
    allele = f[5]
    if not allele.startswith("blaKPC-"): continue
    c = census[allele]
    c["n"] += 1
    y = f[3][:4]
    if y and y != "NULL" and y.isdigit(): c["years"].append(int(y))
    g = f[4] if f[4] and f[4] != "NULL" else "unknown"
    c["geos"][g.split(":")[0]] += 1

observed = set(census.keys())
card_set = {a.replace("KPC-","blaKPC-") for a in AA}
novel_in_outbreak = sorted(observed - card_set)
print("distinct exact alleles observed in outbreak data:", len(observed))
print("observed alleles NOT in CARD reference set (novel candidates):", len(novel_in_outbreak), novel_in_outbreak)

with open(f"{DATA}/kpc_outbreak_atlas.tsv","w",newline="") as f:
    fields = ["allele","n_isolates","year_min","year_max","n_2020plus","top_geo","n_top_geo","n_geos"]
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t"); w.writeheader()
    for a in sorted(census, key=lambda x:-census[x]["n"]):
        c = census[a]
        top = max(c["geos"].items(), key=lambda kv: kv[1]) if c["geos"] else ("unknown",0)
        w.writerow(dict(allele=a, n_isolates=c["n"],
                        year_min=min(c["years"]) if c["years"] else "",
                        year_max=max(c["years"]) if c["years"] else "",
                        n_2020plus=sum(1 for y in c["years"] if y>=2020),
                        top_geo=top[0], n_top_geo=top[1], n_geos=len(c["geos"])))
print("atlas written")
