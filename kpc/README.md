# KPC slice (builder 11)

Structural effects of KPC variants on carbapenem / ceftazidime-avibactam binding across
49,774 blaKPC-bearing Klebsiella isolates (NCBI Pathogen Detection snapshot PDG000000012.2532).

- code/kpc_caller.py - variant caller (nt CDS -> protein -> signature vs KPC-2 -> allele assignment). Gate G1: 231/231 CARD + 295/295 AMRFinder alleles reproduced exactly.
- code/extend_ref.py - builds the 296-allele union reference set; cross-database consistency check.
- code/build_master.py - master allele table: signature (Ambler numbering) x hotspot regions x distance to Ser70 x literature resistance x outbreak census.
- code/figures.py - figures F1-F5.
- data/DATA_MANIFEST.md - byte-lock manifest (sources, accessions, sha256).
- results/master_allele_table.tsv - one row per allele (296).
- results/figures/ - F1-F5 PNG.
- paper/ - research paper (source + PDF).
