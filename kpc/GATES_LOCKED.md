# KPC slice - locked success gates (builder 11)

Locked in writing BEFORE any outcome data was touched: milestone-1 report to orchestrator
2026-09-23 13:29 IST (agent_message), prior to all downloads and analysis.

NOTE (fleet QC rule): after an initial bundle, git history was rewritten (orchestrator QC,
2026-09-23) so this gates doc is the root commit, strictly before every results commit. The
gates themselves were locked and reported to the orchestrator before any download or analysis
(milestone-1 message, 2026-09-23 13:29 IST); the byte-lock manifest (DATA_MANIFEST.md) was
locked the same day before outcome analysis.

## G0 - data byte-lock
Every source with accession lists + sha256 checksums reported before analysis (milestone 2).
Honored: kpc/data/DATA_MANIFEST.md.

## G1 - positive control / validation (binding)
The variant caller must reproduce every blaKPC reference allele in NCBI AMRFinderPlus
ReferenceGeneCatalog + CARD exactly - nucleotide identity AND amino-acid substitution
signature vs KPC-2, 100% exact match, zero mismatches allowed. Any miss fixes the pipeline
before novel claims.
Outcome (reported milestone 3): CARD 231/231, extended AMRFinderPlus 2026-08-07.1 295/295,
zero protein-signature collisions, zero cross-database disagreements on 230 shared alleles.

## G2 - novelty threshold
An amino-acid variant is "novel" only if absent from the full CARD/AMRFinder allele set at
data-lock time; verified/thin/missing counts reported.
Outcome: 0 novel proteins among named outbreak calls; 18 outbreak-observed alleles absent
from CARD but present in AMRFinderPlus 2026-08-07.1 (reported as post-CARD, not novel).

## G3 - structural claim threshold
A binding-effect claim on any variant requires mapping onto >=1 experimental KPC structure
(KPC-2 apo 2OV5/3DW0 + avibactam complex 4ZBE) with residue context: distance to catalytic
Ser70, omega-loop (164-179), 240-loop (237-243), 266-275 hotspot membership, and cross-check
vs the literature-curated CAZ-AVI resistance position set. Claims graded: "documented" only
for literature-confirmed alleles; otherwise "consistent with resistance - hypothesis" at
documented positions/hotspots; else "uncertain/no evidence".

## G4 - honesty
Negative/uninformative results preserved and reported as such; no re-fishing.
Honored: Ser70-proximity null result (p=0.122) reported in paper section 3.5.
