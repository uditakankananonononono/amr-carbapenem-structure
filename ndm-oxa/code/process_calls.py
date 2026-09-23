import json, csv, re, hashlib
from collections import Counter

def load(p):
    d=json.load(open(p))
    return d['ngout']['data']['content']

def toyear(cd):
    if not cd: return None
    m=re.match(r'(\d{4})', cd)
    return int(m.group(1)) if m else None

for tag,path in [("ndm","data/ncbi/ndm_calls_raw.json"),("oxa48","data/ncbi/oxa48_calls_raw.json")]:
    rows=load(path)
    # dedupe by asm_acc x element_symbol x contig position; keep best identity then core scope
    best={}
    for r in rows:
        k=(r.get('asm_acc'),r.get('element_symbol'),r.get('contig_acc'),r.get('start_on_contig'))
        score=(r.get('pct_ref_identity') or 0, 1 if r.get('scope')=='core' else 0)
        if k not in best or score>best[k][0]:
            best[k]=(score,r)
    recs=[v[1] for v in best.values()]
    print(f"== {tag}: raw {len(rows)} -> dedup {len(recs)}")
    alleles=Counter(r['element_symbol'] for r in recs)
    print("distinct alleles:", len(alleles))
    print("top 12:", alleles.most_common(12))
    methods=Counter(r.get('amr_method') for r in recs)
    print("methods:", dict(methods))
    years=Counter(toyear(r.get('collection_date')) for r in recs)
    recent=[r for r in recs if (toyear(r.get('collection_date')) or 0)>=2023]
    print("2023+ records:", len(recent), "unique asm:", len({r['asm_acc'] for r in recent}))
    orgs=Counter(r.get('taxgroup_name') for r in recent)
    print("2023+ top organisms:", orgs.most_common(8))
    geo=Counter((r.get('geo_loc_name') or 'MISSING').split(':')[0] for r in recent)
    print("2023+ top geo:", geo.most_common(8))
    # exact calls (100/100) vs partial
    exact=[r for r in recs if (r.get('pct_ref_identity')==100.0 and r.get('pct_ref_coverage')==100.0)]
    print(f"exact 100/100 calls: {len(exact)} ({100*len(exact)/len(recs):.1f}%)")
    # write tidy csv
    with open(f"data/ncbi/{tag}_calls.tsv","w",newline='') as f:
        w=csv.DictWriter(f, fieldnames=sorted(recs[0].keys()), delimiter='\t')
        w.writeheader()
        for r in recs: w.writerow({k:(v if not isinstance(v,list) else ';'.join(map(str,v))) for k,v in r.items()})
    print("wrote", f"data/ncbi/{tag}_calls.tsv", hashlib.sha256(open(f"data/ncbi/{tag}_calls.tsv",'rb').read()).hexdigest()[:16])
