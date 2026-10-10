# amr-carbapenem-structure

Structural effects of carbapenemase variants on carbapenem / ceftazidime-avibactam binding across recent outbreak genomes.

## Slices
- `kpc/` - KPC variants (builder 11): outbreak mutation atlas + structural effect table + paper. See kpc/README.md.
- `ndm-oxa/` - Builder 12 NDM/OXA-48-like slice is built: locked gate document, code, source manifests, results and paper are present. Readiness review is pending; this is not a new analysis or a success claim.

Each slice is self-contained: `code/`, `data/` (manifests + checksums, not bulk data), `results/`, `paper/`.

## Documentary readiness audit (2026-10-10)

Addendum (2026-10-10, later same day): the amendment set landed - paper reference and gates-path corrections, manifest-count corrections (now 60 entries), and the rebuilt paper PDF; see ndm-oxa/DATA_MANIFEST.md and ndm-oxa/SEAL.md.

Verification only, with no new compute, predictions or nominations. All 52
entries of `ndm-oxa/data/MANIFEST.sha256` matched the current public bytes.
The G4 addendum records a published-genome control, while the existing paper
still says that control is pending. The paper also cites a nonexistent
`ndm-oxa/docs/GATES_LOCKED.md` path (actual path: `ndm-oxa/GATES_LOCKED.md`)
and an older 44-entry manifest count. No repository tags were present at audit.

The missing unit is documentary reconciliation and independent review of
gate evidence, provenance and the current portfolio success bar. No gate pass
or project success is inferred from file presence or checksums. Existing
results and paper are retained unchanged. 
