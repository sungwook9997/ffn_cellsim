# Outer-library built artifacts — where they live, and why they are not in git

**This directory holds no payload. Only this file.** That is the correction of 2026-08-09, made the
same day as the mistake it corrects, and it is the whole point of the page.

## The relocation, and why it was undone

The migration from `Project_Aleph`'s `codex/external-training-corpus` at `803ae3ec` moved 280 MB of
built `.npz` / `.jsonl` here, to `aleph/outputs/outer/`, on the correct principle that a rebuildable
artifact is an output. **The principle was right and the move was wrong**, because of a fact about
the code nobody checked first: every module in `aleph/outer/` resolves its own `results/` through
`Path(__file__).resolve().parent`. Moving the payload out of the package did not teach the pipeline
to read from two roots. It gave it one root with holes.

What that produced:

- **50 further artifacts (406 MB) were dropped entirely**, in neither location. They match the five
  patterns in `experiment_factory/outer_library/results/.gitignore`, inherited unchanged from the
  source repository — so a file-level copy honouring `.gitignore` read *"not versioned"* as *"not
  needed"*. They are neither.
- **14 of 132 tests in `test_outer_library.py` failed.** The worst was not a crash: the coverage
  audit read the holed directory and reported **1,024 observation rows where there are 329,303**,
  as an ordinary low-coverage finding. **A number that small, arrived at that quietly, is how a
  corpus comes to be described as a tenth of a percent of itself.**
- `verification/verify_corpus.py` — tracked source, 17.8 KB — was dropped too, and two more index
  files with it.

**Restored 2026-08-09: 15 files into `aleph/outer/`, and all 234 result artifacts back beside the
code that writes them.** The suite goes 14 failures → 4, and the 4 that remain are unrelated to this
(below). Artifacts now sit where they are produced and are ignored **by exact path** in the root
`.gitignore`, never by extension — `model/combined_content/weights_title_only.npz` is tracked
evidence while `model/oa_content/weights_title_only.npz` is a build artifact, so any glob wide
enough to catch the second silently untracks the first.

## Why gitignored rather than tracked

This tree's rule for a rebuildable artifact is the one the knowledge base already follows: **the
source of truth is the thing that makes it, never the artifact** — *"never hand-edit the vault or the
`.duckdb` — fix Notion and refresh."* A built tensor is the same category. That rule is unchanged;
only the location was wrong.

The rule was also tested today and failed once: a `git add -A` on a moved directory swept 42 MB of
regenerable `kb.duckdb` plus the bio-report PDFs into a commit, and it was reset before the blobs
reached history. Naming these here is what stops the next `add -A` from doing it again.

## The 4 remaining test failures — not from the migration

| test | what differs | class |
|---|---|---|
| `hpa_if_grouped_holdout…` | conformal radius `0.08959868616554265` vs `…437` — **1.3e-15 relative** | float replay across platforms; `REPORT.md` predicted it |
| `bbbc054_microglia_head…` | optimizer selection metrics | same class; `REPORT.md` predicted it |
| `observation_operator_registry_blocks_semantic_shortcuts` | `all_declared_source_evidence_resolves: false` | ⚠ real, below |
| `sweep_gate_resolves_preregistered_emt_markers…` | consumes the report above | cascade of the same |

⚠ **The registry cites a tree that is not this one.** Two entries pin simulator source evidence at
`aleph/vertical/connectors_transfer.py` and `aleph/vertical/cortex_filaments.py`. `aleph/vertical/`
is `Project_Aleph`'s layout and does not exist here.

**Not repointed, deliberately.** Each entry pins a `source_path`, a `source_symbol` and a digest, so
aiming it at this tree's nearest module asserts that module implements the same observation
operator — which is exactly the claim `forbidden_shortcut` exists to block (*"A force-bearing
callable is not a unit- and geometry-matched TFM operator"*). Merge record §1a is also live here: on
the axis these two share, this tree's `virtual_cell/observation_operator.py` was found the stronger.
**Repointing is a PI decision, and it is open.**

## What was NOT copied at all

**About 12.7 GB of untracked raw acquisition data** under
`corpus/external_training/experiment_factory/` in the source worktree — the downloaded microscopy
tiles, TIFFs and archives the pipeline reads. It was never in git there either. It stays at
`~/Project_Aleph/.claude/worktrees/external-training-corpus/`, and re-acquisition is what
`aleph/outer/acquisition/` and `ACQUISITION_POLICY.md` describe.

⚠ **That path is inside a tree decision A1 makes an archive.** Retiring `Project_Aleph` must not
delete that worktree until the raw data has somewhere else to live, or the tensors above become
unrebuildable and this file becomes a lie. **That is a PI decision and it is not scheduled.**

### What the 12.7 GB actually is — measured 2026-08-09

`aleph/outer/acquisition/manifest_raw_acquisition.py` hashed all of it. **89 files, 13.00 GB**, and
it is two blocks, not one:

| block | GB | files | digest pinned at acquisition |
|---|---:|---:|---|
| `data/allencell_label_free_confirmation/` | 11.00 | 64 | ✅ `results/allencell_label_free_acquisition.json` — per file: url, sha256, etag, size |
| `data/…/downloads/` (pbmc_two_lab, pbmc_cross_lab, gse158055) | 1.99 | 25 | ❌ **none, anywhere** |

**All 64 AllenCell payloads verify against their recorded digests. Nothing is contradicted.** That
block was never an unprovenanced pile — the receipts existed, they were simply not in
`acquisition_receipts.jsonl` (which has 53 entries and **zero** for AllenCell), which is why they
read as missing from outside.

The **25 GEO/SCP files have no digest recorded by anyone**. They carry accessions, so they can be
fetched again, but nothing says the bytes that were fetched are the bytes that trained anything.
`raw_acquisition_manifest.jsonl` now records their digests — ⚠ **as of today, not as of
acquisition.** It pins what is on disk; it cannot certify what the source served in August.

**A retracted intermediate result.** The first version of that module checked the AllenCell filename
prefix — 64 hex characters — as if it were a digest of the file, and reported **all 64 payloads
corrupt**. It is AllenCell's S3 logical key. No data was wrong; the check was. Recorded because a
verifier that fails on 100% of its inputs is a finding about the verifier, and because the number
64/64 would otherwise survive in a transcript with no retraction attached to it.

**What this does and does not unblock.** Deleting the source worktree is now *checkable* —
`--check` recomputes and exits non-zero on drift — but a manifest is not a copy. The 11 GB AllenCell
block is re-acquirable from a live S3 package with digests to verify against; the 1.99 GB is
re-acquirable by accession with nothing to verify against. **Whether that is enough to delete is
still the PI's call, and it is still not scheduled.**

## History

The 676-commit branch is fetched into this repository and reachable forever:

```
git log aleph-src/external-training-corpus
git show 803ae3ec
```

The commits carry evidence the files do not — `Seal IDR0168 source integrity refusal`,
`Validate IDR0168 uint16 values independent of byte order`, `Clarify IDR0168 tar stream diagnostic`.
That is a refusal being derived, tested and sealed across four commits, which is exactly the kind of
record this project keeps. Read the branch, not just the tree.
