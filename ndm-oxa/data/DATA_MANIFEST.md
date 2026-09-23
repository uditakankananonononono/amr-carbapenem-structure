# NDM/OXA-48 slice - data provenance and re-derivation manifest

All payloads committed under `ndm-oxa/data/` are checksummed in `MANIFEST.sha256`
(regenerated 2026-09-23 ~14:25 IST to cover exactly the committed set).

## Sources

| Payload | Source | Notes |
|---|---|---|
| data/card/card-data.tar.bz2 | https://card.mcmaster.ca/latest/data (CARD "data" tarball) | sha256 94b89501...; downloaded 2026-08-11 (file dates), byte-identical re-download 2026-09-23 |
| data/card/ndm_ref.fasta, oxa48like_ref.fasta, oxa48_family.json | derived from card.json inside the tarball by code/build_refs.py | byte-identical regeneration (shas e37ebc63..., b1da8795..., 348a1acf...) |
| data/ncbi/ndm_calls_raw.json, oxa48_calls_raw.json | NCBI Pathogen Detection pathogens-srv API, collection=amr, filters element_symbol:blaNDM* and element_symbol:(blaOXA-48* OR ... 67 family wildcards) | analysis snapshot 2026-09-23 12:34 IST. A 14:24 IST re-pull returned 2/6 extra backend-appended rows (NCBI live update); the appended rows (ids listed below) were trimmed and the trimmed bytes reproduce the original snapshot shas exactly (b66b9609..., 77812230...) |
| data/ncbi/ndm_calls.tsv, oxa48_family_calls.tsv | derived from the raw JSONs by code/process2.py | byte-identical to analysis snapshot (9be8d5d9..., 6923a014...) |
| data/ncbi/download_targets.json, download_targets_stats.json, novel_targets.json | derived from the TSVs (recipe in commit history / process2 workflow) | byte-identical (1b4a8ecd..., 60ed4e71...) |
| data/ncbi/erd_map_recent.json | NCBI pathogens-srv collection=isolates, fq asm_acc==[...] 100-batches | DRIFT: NCBI re-clusters erd SNP groups continuously. Committed copy re-derived 2026-09-23 14:22 IST, sha 04b6c348...; the analysis used the 13:05 IST snapshot sha a08622a5... (not byte-restorable). results/cluster_spread.json was computed from the original snapshot and is byte-locked. |
| data/porin/*.fasta | NCBI eutils efetch: kp_ompK35/ompK36 from NZ_OZ547068 (fasta_cds_na), ec_ompC/ec_ompF from NC_000913.3 coords 2311646-2312749 / 985894-986982 (strand 2) | byte-identical regeneration |
| data/pdb/*.pdb | RCSB PDB: 3SPU (NDM-1), 4EYL, 4S2P (OXA-48), 6P97 | byte-identical re-download |
| data/g4/*, code/g4_porin_call.py | G4 porin-validation workbench; PLOS Pathogens 2022 S1 Table source: https://journals.plos.org/plospathogens/article/file?type=supplementary&id=10.1371/journal.ppat.1010334.s003 (xlsx committed at data/g4/ppat.1010334.s003.xlsx, sha256 4810d146ebc0cfa20e5af84967eb8de5e606bfa2c7f4959810de8dc6b0977d8b) | run DRR065574 (SRA), reads fetched via SRA toolkit |

## Not committed (re-derivable)

- data/asm/zips/batch000-014.zip (original shas in git history, manifest v1): NCBI assembly
  packages for the 722 accessions in data/ncbi/download_targets.json, fetched by
  code/fetch_asm.sh from the NCBI Datasets API. Re-derivable; large.
- data/query_refs.fasta (original sha 8f06868421e7...): concatenated BLAST query refs
  derived from the committed CARD + porin fastas by the chunk_blast.sh workflow.
- results/blast/all_hits.tsv (original sha 5e53bb05...): BLASTN output of
  query_refs vs the 722-genome database; regenerable with code/chunk_blast.sh
  (use -max_target_seqs >= 20000).

## Backend-drift record (2026-09-23)

NDM re-pull appended rows 022_PDT003397934.1, 041_PDT003397934.1.
OXA-48-family re-pull appended rows 017_/025_PDT003380424.1, 023_/042_PDT003397930.1,
039_/051_PDT003397934.1. All reported results use the trimmed snapshot.
