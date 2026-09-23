import json, csv, re, hashlib
from collections import Counter

fam=json.load(open("data/card/oxa48_family.json"))
fam_syms={"bla"+n for n in fam}
# also allow letter-suffixed (e.g. blaOXA-48b) if base numeric is in family
def in_family(sym):
    if sym in fam_syms: return True
    m=re.fullmatch(r"blaOXA-(\d+)[a-z]", sym or "")
    return bool(m) and ("blaOXA-"+m.group(1)) in fam_syms

rows=json.load(open("data/ncbi/oxa48_calls_raw.json"))['ngout']['data']['content']
kept=[r for r in rows if in_family(r.get('element_symbol'))]
print(f"oxa48 raw {len(rows)} -> family-filtered {len(kept)}")
alleles=Counter(r['element_symbol'] for r in kept)
print("distinct family alleles observed:", len(alleles))
print("top 15:", alleles.most_common(15))
dropped=Counter(r['element_symbol'] for r in rows if not in_family(r.get('element_symbol')))
print("dropped (non-family):", dropped.most_common(6))

def toyear(cd):
    m=re.match(r'(\d{4})', cd or "")
    return int(m.group(1)) if m else None
recent=[r for r in kept if (toyear(r.get('collection_date')) or 0)>=2023]
print("2023+ family records:", len(recent), "unique asm:", len({r['asm_acc'] for r in recent}))
orgs=Counter(r.get('taxgroup_name') for r in recent)
print("2023+ top organisms:", orgs.most_common(8))

# unified field union writer
def write_tsv(recs, path):
    fields=sorted({k for r in recs for k in r})
    with open(path,"w",newline='') as f:
        w=csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for r in recs:
            w.writerow({k:(v if not isinstance(v,list) else ';'.join(map(str,v))) for k,v in r.items() if k in fields})
    return hashlib.sha256(open(path,'rb').read()).hexdigest()[:16]

print("oxa48_family_calls.tsv", write_tsv(kept,"data/ncbi/oxa48_family_calls.tsv"))
# rewrite ndm with union writer too
ndm=json.load(open("data/ncbi/ndm_calls_raw.json"))['ngout']['data']['content']
print("ndm_calls.tsv", write_tsv(ndm,"data/ncbi/ndm_calls.tsv"))

# novel/non-exact NDM calls (for targeted assembly download)
novel=[r for r in ndm if not (r.get('pct_ref_identity')==100.0 and r.get('pct_ref_coverage')==100.0)]
print("NDM non-exact calls:", len(novel), "methods:", Counter(r['amr_method'] for r in novel))
novel_recent=[r for r in novel if (toyear(r.get('collection_date')) or 0)>=2023]
print("NDM non-exact 2023+:", len(novel_recent), "unique asm:", len({r['asm_acc'] for r in novel_recent}))
bare=[r for r in ndm if r.get('element_symbol')=='blaNDM']
print("bare blaNDM calls:", len(bare), "recent:", len([r for r in bare if (toyear(r.get('collection_date')) or 0)>=2023]))
