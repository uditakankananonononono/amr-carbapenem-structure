# LOCKED SUCCESS GATES — amr-carbapenem-structure, Builder 12 slice
## NDM + OXA-48-like variants + porin-loss combination effects
Locked in writing 2026-09-23 13:30 IST, BEFORE any outcome data is inspected.

## Scope
- Enzyme classes: blaNDM (class B metallo-beta-lactamase) and blaOXA-48-like (class D carbapenemase) alleles.
- Porins: OmpK35/OmpK36 (Klebsiella pneumoniae) and OmpC/OmpF (E. coli) loss-of-function in carbapenemase-positive isolates.
- Data: NCBI Pathogen Detection isolate metadata + selected assemblies; CARD/ResFinder reference alleles (positive control); NDM/OXA-48 PDB structures.
- No wet lab. Open public data only. Sequence-level first; structure mapping onto existing PDB structures (homology-based positional mapping, not de novo folding of hundreds of proteins).

## GATES
G1 (positive control, PASS/FAIL): Variant-calling pipeline (blastn vs CARD nucleotide reference) must reproduce CARD reference allele calls at 100% identity/coverage on CARD's own curated NDM + OXA-48-like reference sequences. PASS requires 100% exact allele recovery (every reference sequence calls its own allele, no mismatches, no novel flags on references). If <100%, pipeline is fixed and re-gated before any novel claim.

G2 (independent positive control): On a set of public genomes with independently published allele calls (ResFinder/PubMLST-annotated), pipeline calls must agree with the published calls on >=95% of calls; every disagreement is manually resolved and reported by count. Disagreement count reported VERIFIED/THIN/MISSING style.

G3 (structural mapping): EVERY novel or outbreak-enriched amino-acid variant claim must be mapped onto the matching enzyme class structure (NDM -> NDM PDB structure; OXA-48-like -> OXA-48 PDB structure) with residue-number correspondence verified by alignment to the PDB sequence (100% coverage of the variant position, mismatches reported). Claims about effect on drug binding must be annotated by structural context (zinc-ligating, active-site loop, omega loop, dimer interface, surface/other) — NOT by sequence-only speculation.

G4 (porin-loss combination): Porin loss-of-function calls (premature stop, frameshift, IS insertion, large deletion, OmpK36 loop-3 insertion counted separately) must be validated on at least one published porin-loss genome (positive control) recovering the known porin disruption. Combination effect claims are limited to co-occurrence quantification in recent isolates (carbapenemase + porin state), with honest VERIFIED/THIN/MISSING counts. No MIC prediction claims without phenotype data; when NCBI Pathogen Detection phenotype fields are absent, state MISSING.

G5 (honest reporting): Every variant/effect table carries verified/thin/missing counts. Negative results (e.g., a variant class showing no structural clustering, or porin-loss rates not differing between enzyme classes) are preserved and reported, never re-fished. No stubs, simulations, or placeholders.

## Prior-art verdict: CROWDED (closest works named in milestone-1 report). Distinct angle: integrated, CARD-validated atlas of NDM + OXA-48-like allele variation in RECENT (2023+) outbreak-context genomes with per-variant structural mapping AND porin-loss combination states in one effect table.
