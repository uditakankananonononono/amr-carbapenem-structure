#!/usr/bin/env python3
import sys, csv
sys.path.insert(0,"code")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

DATA="data"; RES="results"
rows=list(csv.DictReader(open(f"{RES}/master_allele_table.tsv"),delimiter="\t"))
obs=[r for r in rows if int(r["n_isolates"])>0]
res_lit={r["allele"] for r in rows if r["cazavi_resistant_lit"]=="yes"}

# ---- census timeseries from accessions ----
by_year=defaultdict(lambda: defaultdict(int))
by_geo=defaultdict(int); res_geo=defaultdict(int)
for line in open(f"{DATA}/kpc_isolate_accessions.tsv"):
    f=line.rstrip("\n").split("\t")
    if len(f)<7: continue
    a=f[5]; y=f[3][:4]; g=(f[4] if f[4] and f[4]!="NULL" else "unknown").split(":")[0]
    if not a.startswith("blaKPC-"): continue
    name=a.replace("bla","")
    if y.isdigit() and 2000<=int(y)<=2026:
        cls = "KPC-2" if name=="KPC-2" else "KPC-3" if name=="KPC-3" else ("CAZ-AVI-resistant allele" if name in res_lit else "other allele")
        by_year[int(y)][cls]+=1
    by_geo[g]+=1
    if name in res_lit: res_geo[g]+=1

plt.rcParams.update({"font.size":9,"figure.dpi":150})

# F1: top 25 alleles
top=sorted(obs,key=lambda r:-int(r["n_isolates"]))[:25]
fig,ax=plt.subplots(figsize=(7.5,4))
colors=["#c0392b" if r["allele"] in res_lit else ("#7f8c8d" if r["allele"] in ("KPC-2","KPC-3") else "#2980b9") for r in top]
ax.bar([r["allele"].replace("KPC-","") for r in top],[int(r["n_isolates"]) for r in top],color=colors)
ax.set_yscale("log"); ax.set_ylabel("isolates (log)"); ax.set_xlabel("KPC allele number")
ax.set_title("KPC allele census, NCBI Pathogen Detection Klebsiella snapshot PDG000000012.2532 (n=49,774 isolates)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#7f8c8d",label="KPC-2/3 (founder)"),Patch(color="#2980b9",label="other allele"),Patch(color="#c0392b",label="documented CAZ-AVI resistance")],frameon=False)
plt.tight_layout(); plt.savefig(f"{RES}/figures/F1_allele_census.png"); plt.close()

# F2: per-year stacked
years=sorted(by_year)
fig,ax=plt.subplots(figsize=(7.5,4))
bottom=np.zeros(len(years))
for cls,color in [("KPC-2","#95a5a6"),("KPC-3","#7f8c8d"),("other allele","#2980b9"),("CAZ-AVI-resistant allele","#c0392b")]:
    vals=np.array([by_year[y][cls] for y in years])
    ax.bar(years,vals,bottom=bottom,label=cls,color=color)
    bottom+=vals
ax.set_xlabel("collection year"); ax.set_ylabel("KPC+ isolates")
ax.set_title("KPC+ Klebsiella isolates per year by allele class")
ax.legend(frameon=False)
plt.tight_layout(); plt.savefig(f"{RES}/figures/F2_yearly_class.png"); plt.close()

# F3: variation map across observed alleles
import re
pos_hits=defaultdict(set)
for r in obs:
    for tok in r["signature_ambler"].split("; "):
        m=re.search(r'(\d+)',tok)
        if m and tok!="identical to KPC-2": pos_hits[int(m.group(1))].add(r["allele"])
fig,ax=plt.subplots(figsize=(7.5,4))
xs=sorted(pos_hits); ys=[len(pos_hits[x]) for x in xs]
ax.bar(xs,ys,color="#34495e",width=1.0)
for span,lab in [((164,179),"omega loop"),((237,243),"loop 237-243"),((266,275),"loop 266-275")]:
    ax.axvspan(span[0],span[1],color="#e74c3c",alpha=0.15)
    ax.text(sum(span)/2,max(ys)*1.02,lab,ha="center",fontsize=8,color="#c0392b")
for p,lab in [(179,"D179"),(243,"T243")]:
    ax.annotate(lab,(p,len(pos_hits.get(p,[]))),textcoords="offset points",xytext=(0,6),ha="center",fontsize=8)
ax.set_xlabel("Ambler position"); ax.set_ylabel("# observed alleles with variant at position")
ax.set_title(f"Convergent variation map across {len(obs)} outbreak-observed KPC alleles")
plt.tight_layout(); plt.savefig(f"{RES}/figures/F3_variation_map.png"); plt.close()

# F4: distance to S70
dres=[float(r["min_dist_S70_A"]) for r in obs if r["allele"] in res_lit and r["min_dist_S70_A"]]
doth=[float(r["min_dist_S70_A"]) for r in obs if r["allele"] not in res_lit and r["min_dist_S70_A"]]
fig,ax=plt.subplots(figsize=(6,4))
b=np.arange(0,40,2)
ax.hist(doth,bins=b,alpha=0.7,label=f"other observed (n={len(doth)})",color="#2980b9")
ax.hist(dres,bins=b,alpha=0.7,label=f"documented resistant (n={len(dres)})",color="#c0392b")
ax.set_xlabel("min CA distance of variant to catalytic Ser70 (A)"); ax.set_ylabel("# alleles")
ax.set_title("Variant proximity to catalytic Ser70 (3DW0)")
ax.legend(frameon=False)
plt.tight_layout(); plt.savefig(f"{RES}/figures/F4_dist_s70.png"); plt.close()

# F5: geography of resistant alleles
topg=sorted(res_geo.items(),key=lambda kv:-kv[1])[:15]
fig,ax=plt.subplots(figsize=(7,4))
ax.barh([g for g,_ in topg][::-1],[n for _,n in topg][::-1],color="#c0392b")
ax.set_xlabel("isolates with documented CAZ-AVI-resistant KPC allele")
ax.set_title("Geography of CAZ-AVI-resistant KPC alleles (top 15 countries)")
plt.tight_layout(); plt.savefig(f"{RES}/figures/F5_geo_resistant.png"); plt.close()

print("figures done")
print("resistant median dist:", np.median(dres), "other median:", np.median(doth))
from scipy import stats as st
print("Mann-Whitney:", st.mannwhitneyu(dres,doth,alternative="less"))
