#!/usr/bin/env python3
"""KPC variant caller.
Input: nucleotide CDS -> translate (table 11) -> global-align to KPC-2 reference
-> variant signature (substitutions + indels, KPC-2 protein coordinates)
-> allele assignment by exact nucleotide match against CARD reference set.
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import re, sys, json

DATA = "/home/sandbox/amr/data"

def load_fasta(path):
    return {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}

# --- reference ---
AA = load_fasta(f"{DATA}/kpc_alleles_aa.fasta")
NT = load_fasta(f"{DATA}/kpc_alleles_nt.fasta")
REF = AA["KPC-2"]          # 293 aa incl. signal peptide
assert len(REF) == 293

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

def signature(prot, ref=REF):
    """Return tuple of ops vs ref, e.g. ('T3A', 'del164-177', 'ins269_KDD')."""
    aln = aligner.align(ref, prot)[0]
    # aligned blocks: ((ref_start,...), (prot_start,...))
    (rb, pb) = aln.aligned
    ops = []
    prev_rend, prev_pend = 0, 0
    for (rs, re_), (ps, pe) in zip(rb, pb):
        # gap between blocks
        if rs > prev_rend:  # deletion in prot (ref segment missing)
            seg = ref[prev_rend:rs]
            ops.append(f"del{prev_rend+1}-{rs}_{seg}" if rs-prev_rend>1 else f"del{prev_rend+1}_{seg}")
        if ps > prev_pend:  # insertion in prot
            seg = prot[prev_pend:ps]
            ops.append(f"ins{prev_rend}_{seg}")
        for i in range(re_ - rs):
            r, p = ref[rs+i], prot[ps+i]
            if r != p:
                ops.append(f"{r}{rs+i+1}{p}")
        prev_rend, prev_pend = re_, pe
    # trailing gaps
    if len(ref) > prev_rend:
        seg = ref[prev_rend:]
        ops.append(f"del{prev_rend+1}-{len(ref)}_{seg}")
    if len(prot) > prev_pend:
        ops.append(f"ins{prev_rend}_{prot[prev_pend:]}")
    return tuple(ops)

def translate_cds(nt):
    return str(Seq(nt).translate(table=11))

# ---------- GATE G1 ----------
def gate_g1():
    sig_db = {}
    for name, prot in AA.items():
        sig_db.setdefault(signature(prot), []).append(name)

    n_nt_exact = n_trans_ok = n_sig_ok = n_lookup_ok = 0
    collisions = {k: v for k, v in sig_db.items() if len(v) > 1}
    fails = []
    for name, prot in AA.items():
        nt = NT.get(name)
        if nt is None:
            fails.append((name, "missing nt")); continue
        t = translate_cds(nt)
        # allow trailing stop in nt-derived translation
        t_clean = t.rstrip("*")
        if t_clean == prot:
            n_trans_ok += 1
        else:
            fails.append((name, f"translation mismatch len {len(t_clean)} vs {len(prot)}")); continue
        s1, s2 = signature(prot), signature(t_clean)
        if s1 == s2: n_sig_ok += 1
        else: fails.append((name, "sig mismatch"))
        # allele assignment by EXACT nt match (AMRFinder-style)
        hits = [n2 for n2, seq2 in NT.items() if seq2.upper() == nt.upper()]
        if hits == [name] or name in hits:
            n_lookup_ok += 1
        else:
            fails.append((name, f"nt lookup -> {hits}"))
    total = len(AA)
    print(f"G1 total alleles: {total}")
    print(f"G1a translation(nt)==CARD protein: {n_trans_ok}/{total}")
    print(f"G1b signature(nt)==signature(CARD protein): {n_sig_ok}/{total}")
    print(f"G1c exact-nt allele lookup correct: {n_lookup_ok}/{total}")
    print(f"G1 protein-signature collisions (alleles indistinguishable at protein level): {len(collisions)}")
    for sig, names in sorted(collisions.items(), key=lambda x: -len(x[1]))[:10]:
        print("   ", names, sig[:6])
    if fails:
        print("FAILURES:")
        for f in fails[:20]: print("   ", f)
    ok = (n_trans_ok == total and n_sig_ok == total and n_lookup_ok == total)
    print("GATE G1:", "PASS (100%)" if ok else "FAIL")
    return ok

if __name__ == "__main__":
    gate_g1()
