# ndm-oxa slice seal - 2026-10-10

Sealed commit: 75878a284c3e3aec72e0bdd8efcb8bd8f95eff0c
Root manifest: ndm-oxa/data/MANIFEST.sha256 (sha256 of file: 03b6b803a8f549b2ae1c5abecf9a791003c783fb48f51a34fa16a611f8c9dab1; 60 entries, sha256sum -c clean at the sealed commit)

Seal criteria (all met at the sealed commit):
1. Five G4-pending references replaced with completed status citing data/g4/g4_evidence.json and G4_MANIFEST_ADDENDUM.sha256 (commit 5d90dc99).
2. Gates path corrected to ndm-oxa/GATES_LOCKED.md (commit 5d90dc99).
3. Committed-set vs manifest-set distinction documented: 52-entry committed set vs 44-entry data/MANIFEST.json; manifest now 60 entries (commits 5d90dc99, bd076c4e, 6a788aad).
4. Precedence claim restated: current DAG places results/paper before standalone gate commit 2ca48cc; precedence is supported by that commit's attestation of the 13:29-13:30 lock plus matching content hash e74ab5a...070de, following the reported force-push. NOT DAG proof (commit 5d90dc99).
5. paper.pdf rebuilt from amended sources; result bytes untouched; sidecar refreshed (commits 3f14c8cc + 75878a284c3e3aec72e0bdd8efcb8bd8f95eff0c).
6. MANIFEST.sha256 covers 60 entries including GATES_LOCKED.md, code/g4_porin_call.py, data/DATA_MANIFEST.md, data/MANIFEST.json, data/download_targets_stats.json, data/oxa48_family.json, paper/paper.pdf, paper.pdf.sha256; manifest self-hash limitation and both duplicate pairs documented in DATA_MANIFEST.md (commits 6a788aad, bd076c4e).
7. kpc chronology note: equal author timestamps are not ordering evidence; DAG order is (commit 32b1f209).

Attestation artifacts: this SEAL.md and the tag seal-ndm-oxa-2026-10-10 are attestation records created after the sealed commit and are intentionally outside MANIFEST.sha256, documented as such here.
