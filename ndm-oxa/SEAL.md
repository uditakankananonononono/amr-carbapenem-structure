# ndm-oxa slice seal - 2026-10-10

Sealed commit: 75878a284c3e3aec72e0bdd8efcb8bd8f95eff0c
Root manifest: ndm-oxa/data/MANIFEST.sha256 (sha256 of file: 03b6b803a8f549b2ae1c5abecf9a791003c783fb48f51a34fa16a611f8c9dab1; 60 entries, sha256sum -c clean at the sealed commit)

Seal criteria (status at the sealed commit - qualified by the 2026-10-10 erratum below: two sub-manifests fail and the tag target differs from the sealed commit):
1. Five G4-pending references replaced with completed status citing data/g4/g4_evidence.json and G4_MANIFEST_ADDENDUM.sha256 (commit 5d90dc99).
2. Gates path corrected to ndm-oxa/GATES_LOCKED.md (commit 5d90dc99).
3. Committed-set vs manifest-set distinction documented: 52-entry committed set vs 44-entry data/MANIFEST.json; manifest now 60 entries (commits 5d90dc99, bd076c4e, 6a788aad).
4. Precedence claim restated: current DAG places results/paper before standalone gate commit 2ca48cc; precedence is supported by that commit's attestation of the 13:29-13:30 lock plus matching content hash e74ab5a...070de, following the reported force-push. NOT DAG proof (commit 5d90dc99).
5. paper.pdf rebuilt from amended sources; result bytes untouched; sidecar refreshed (commits 3f14c8cc + 75878a284c3e3aec72e0bdd8efcb8bd8f95eff0c).
6. MANIFEST.sha256 covers 60 entries including GATES_LOCKED.md, code/g4_porin_call.py, data/DATA_MANIFEST.md, data/MANIFEST.json, data/download_targets_stats.json, data/oxa48_family.json, paper/paper.pdf, paper.pdf.sha256; manifest self-hash limitation and both duplicate pairs documented in DATA_MANIFEST.md (commits 6a788aad, bd076c4e).
7. kpc chronology note: equal author timestamps are not ordering evidence; DAG order is (commit 32b1f209).

Attestation artifacts: this SEAL.md and the tag seal-ndm-oxa-2026-10-10 are attestation records created after the sealed commit and are intentionally outside MANIFEST.sha256, documented as such here.

## Erratum 2026-10-10 (record repair, post-seal disclosure)

Registry-level re-verification at main HEAD e36a8cccd403e168bf7e0d9c9603d9ffd9fabd28 found two stale sub-manifest entries and a tag-target mismatch. The root manifest result stands; the package must not be represented as cleanly sealed without this qualification.

1. Root manifest preserved: ndm-oxa/data/MANIFEST.sha256 verifies 60/60 clean at e36a8cc (sha256sum -c run from ndm-oxa/). Unaffected by the items below.
2. data/g4/G4_MANIFEST_ADDENDUM.sha256 FAILS on ndm-oxa/paper/paper.pdf: recorded d6f00845e0eac6dd0ec581e49e499eda3e12a1d79c2c69ebdad2d3ac0550b447 vs actual 5bbb6816b27e4a0f11f3f460473bc0efbad5464e0fed2ddf8cabed85c5cf5503. The addendum predates the paper rebuild (last touched in 4d963e8, before the sealed commit) and was never refreshed after it; the actual paper.pdf hash matches the committed sidecar paper.pdf.sha256 and the root manifest entry. This failure already existed at the sealed commit 75878a2; the "sha256sum -c clean" claim in this record covers the root manifest only.
3. kpc/MANIFEST.sha256 FAILS on ./README.md: recorded 004743230f0eec6fd935e303a40b6d405b258a409129b463163e9a0d835bd8a4 vs actual c961b2ead40c9325eb529e01628d08c352b23fec1b203786142dc6747949817a. Stale since 32b1f20 (kpc chronology note, review item A7) amended kpc/README.md after the kpc manifest was locked. Also already failing at the sealed commit.
4. Tag-target mismatch: annotated tag seal-ndm-oxa-2026-10-10 (tag object dbc93369cf47ed9ce0f49ca54c44318064a21b6a, unsigned) points to commit c6e9f6b5ffcf25a1f613d7c316eb82b63be6cee9 - the SEAL.md commit - NOT the sealed commit 75878a2 named above, and one commit behind main HEAD e36a8cc. The tag must not be described as sealing 75878a2. Tag retargeting is a maintainer decision; it has been flagged and no tag operation was performed.
5. Repair choice: the stale sub-manifest entries are disclosed, not refreshed, because G4_MANIFEST_ADDENDUM.sha256 is itself covered by the root manifest (entry 10), so rewriting it would cascade into the root manifest this seal preserves. Git history retains all prior bytes; nothing was rewritten.
