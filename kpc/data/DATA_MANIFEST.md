# KPC slice - data byte-lock manifest (locked 2026-09-23, builder 11)

All downloads performed 2026-09-23. sha256 recorded at lock time.

## Sources (bulk; not committed - re-fetch to reproduce)
| source | URL | bytes | sha256 |
|---|---|---|---|
| CARD data tarball | https://card.mcmaster.ca/latest/data | 4,366,637 | 94b89501fc26daafbe65e4151c17a8b007dc1272c4f8201737e44e44cdd4c065 |
| NCBI Pathogen Detection Klebsiella snapshot PDG000000012.2532 AMR metadata | https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/PDG000000012.2532/AMR/PDG000000012.2532.amr.metadata.tsv | 190,586,905 | 45d1f9d476333b61fc9e7b77ded0f38a6253e81d5082d64f6059ea8896a58e4a |
| AMRFinderPlus db 2026-08-07.1 proteins | https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt.fa | - | 41d5ebf4f807c9590f27de7dd23c9f037afb4d8efd613790c54847fb1fc2896b |
| AMRFinderPlus db 2026-08-07.1 CDS | .../AMR_CDS.fa | - | bcc57417c5a170fb4a664de6898455990fcc296dbfe214b6874d761dad63962c |
| CAZ-AVI resistance review (full text via web fetch) | https://pmc.ncbi.nlm.nih.gov/articles/PMC9487638/ | - | see cazavi_table1_raw.md sha 80771d3b203dfa5c8e7cfb7ca14207312cfe2182af6ae811e42327d22ec61a6c |

## PDB structures (mmCIF)
| id | description | sha256 |
|---|---|---|
| 2OV5 | WT KPC-2 apo 1.85A | 099700715f8d734a453b03182a3b97cbd85ee8049f043e79d1265755e3305b79 |
| 3DW0 | WT KPC-2 apo 1.6A | 15a340a9081d82512150cf59a5b6ff946014728cb0c690ce94558f668df62ebb |
| 4ZBE | WT KPC-2 + avibactam | dae1bc81c8311d3ad006a83900e5cf341d8938a418f88f34f303fb555d630a75 |
| 8AKK | KPC-2 E166Q + imipenem acyl-enzyme | 97d81c547d93cbd5da25111a58ff6508c3f0e4ebf0517e3829e2ed5a4967853d |
| 8AKL | KPC-2 E166Q + meropenem acyl-enzyme | 3dede707c60128dfe38edb450a7bedd6ddee397f60d877fd9f5d0992c81cd89f |
| 7TB7 | D179N KPC-2 apo 0.99A | 2b7255302c8454a4ca0b730e4e4261c995bdfe6e19f465dea3a62f3d7ff257ba |
| 8G2R | D179N KPC-2 + avibactam 1.28A | 93d5981931a6b78678e3d267b3e1b432acc5165c8a779e74d97db938f20d09ea |
| 8TMR | KPC-44 + avibactam 1.37A | 257328b4c8c31fa8222a1f9332a2cc342b2832e74a756b468edf9a6626b752e0 |

## Committed derived reference files
| file | sha256 |
|---|---|
| kpc_ref_union_aa.fasta (296-allele union: CARD 231 + AMRFinder 295, 230 overlap, 0 disagreements) | 527d1a3983679864e86bad8ca1207e254691fcb4b0e0643d2bcfefca798ac42f |
| kpc_ref_union_nt.fasta | eeb30475d7569346914158bb619c7fc240fa51b3fe7dfd7dd5597a7b96bf4633 |

## Derived census (from PD snapshot, quote/whitespace-normalized blaKPC calls)
- 49,774 blaKPC-bearing isolates; 50,238 isolate-allele pairs; 169 distinct exact alleles
- kpc_isolate_accessions.tsv sha256 50661ee8752d50f1d62061fd97a024d8f15091c20d283fb02187848496694762 (not committed; regenerable via code/build_master.py against the locked source)
- kpc_allele_census.txt (committed) sha256 094f3b4a09e699ef6c6a7ca873bd1c4858976921e37fb8b700f7a4b602f79063
