# Independent external-corpus verification — 2026-08-05

Status: **AGENT-PROPOSED**.  This audit verifies corpus structure, byte and
digest integrity, retrieval authority boundaries, and current Git tracking.  It
does not establish scientific correctness, study quality, licence
interpretation, or physical validity.

## Independently recomputed inventory

The verifier reads the six original candidate manifests, the independent
second-pass candidate manifest, both acquisition receipt files, every referenced
OA object, and both ignored derived stores.  Counts are recomputed rather than
copied from producer summaries.

| Measure | First pass | Second pass | Combined |
|---|---:|---:|---:|
| Candidate rows | 7,218 | 3,000 | 10,218 raw |
| Unique candidate families | 7,213 | 3,000 | 10,208 |
| Acquisition receipts | 2,713 | 3,000 | 5,713 |
| Successful OA objects | 2,701 | 2,974 | 5,675 |
| Failed receipts | 12 | 26 | 38 |
| Verified XML bytes | 445,899,949 | 404,754,837 | 850,654,786 |
| Retrieval chunks | 139,573 | 145,047 | 284,620 |
| Sections | 79,071 | 82,918 | 161,989 |
| Figures located | 17,173 | 21,143 | 38,316 |
| Tables located | 2,213 | 2,237 | 4,450 |

There are five normalized candidate-family overlaps between the 7,213-family
initial metadata set and the 3,000-family second-pass metadata set.  There are
zero overlaps between the two successful OA object/derived family sets.

All 5,675 receipt-referenced object files were present, had the receipt byte
length, and reproduced the receipt SHA-256.  The object directory contained
exactly those 5,675 files: zero missing and zero unreferenced.  The verifier
also recomputed every chunk-file, text, and chunk-record hash over all 284,620
records; no digest, count, orphan-file, or `proposed`-authority error occurred.

The recomputed article-manifest set hashes are:

- first pass: `cf618563667281631b8df4f5453a83c65402e95ae619ba8788c72297607ac069`
- second pass: `9627db89509f0a2f363ec93e65a2da4b68d3ab1b0415195788145b54a30bcf0f`

## Git tracking hygiene

At the final audit, `git ls-files` contained no path beneath
`data/external_training/`, no `*.chunks.jsonl`, and no tracked PDF/XML/NXML/JATS
file.  A reachable-history path-name scan with `git rev-list --objects --all`
also returned none of those path classes.  Candidate metadata and digest-only
summaries remain tracked by design.

## Combined retrieval audit and boundary smoke

The combined index was independently decoded and checked against both derived
stores.  It contains 5,675 unique families and 284,620 chunks: 2,701 families
from the first pass and 2,974 from the second.  All 5,675 referenced chunk files
exist and reproduce their index digests.  The index contains no text, snippet,
abstract, body, or full-text field.

Exactly 5,669 index families are neural-routable.  Their family set equals the
combined model's split-manifest family set; the remaining six families are
retained in the structural index but have no usable model input.  The following
bindings all reproduced independently:

- index SHA-256 `ede8d2aff2a8786efa8e01bba6fbb398f46bdfdcc34eda0ef2662bdefe7a5887`;
- combined model weights, splits, and taxonomy digests from commit `5336672`;
- combined RAG validation receipt and ledger/TAG/projection hashes from commit
  `736153b`;
- combined retrieval index provenance at commit `3dbb3aa`.

Four allowed queries and six adversarial refusals were then exercised through
`CombinedExternalCorpusRetriever`.  Every allowed result remained `proposed`
and carried digest-verified exact locators.  The cortical-tension, TFM, and PIV
queries each returned results from **both** passes; the multimodal-state query
returned second-pass results.  Thus both passes were observed through the
actual query surface, not inferred from index membership.

The six refused cases covered three parameter estimation/calibration variants,
one training request, one authority/physics-proof request, and one unknown-token
query.  All returned their expected typed codes.  In particular, the mixed
instruction “search literature and then fit membrane tension” was refused.

The earlier first-pass-only bridge was also retained as a regression boundary.
Twelve historical adversarial cases were exercised through
`ExternalCorpusRetriever`:

- allowed, one proposed source returned: measurement-method and literature
  queries for cortical tension, Young's modulus, viscosity, and TFM — 5/5;
- refused: four parameter estimation/calibration variants, one training
  request, one authority/physics-proof request, and one unknown-token query —
  7/7 with the expected typed codes.

No returned snippet is stored in this report or any verification JSON.

## Commands and test results

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python \
  corpus/external_training/verification/verify_corpus.py --repo-root .

/Users/sw1/miniconda3/envs/aleph/bin/python \
  corpus/external_training/verification/smoke_retrieval.py

/Users/sw1/miniconda3/envs/aleph/bin/python \
  corpus/external_training/verification/smoke_combined_retrieval.py

/Users/sw1/miniconda3/envs/aleph/bin/python -m pytest \
  corpus/external_training -q

/Users/sw1/miniconda3/envs/aleph/bin/python -m pytest \
  tests/harness/test_cag.py tests/harness/test_ledger.py \
  tests/harness/test_objects.py tests/harness/test_snapshot.py \
  tests/harness/test_linkgrammar.py tests/harness/test_paths.py -q
```

The first default corpus collection encountered a foreign duplicate-module-name
collision between the two review lanes' `test_generated_results.py` files.  No
producer was edited by this audit.  The review owner independently landed
`dbc413d`, adding package isolation; the subsequent default full-corpus run was
green.  Exact final test counts are recorded in `test_results.json`.

## Limitations

- Structural and digest integrity do not validate extracted scientific claims.
- Access/licence fields are receipt provenance, not a new legal determination.
- Git hygiene checks current tracked paths and reachable object path names; they
  do not semantically classify arbitrary prose in all historical blobs.
- Six structurally indexed families are not neural-routable because the model
  input audit found no usable target content.  They remain visible in the
  5,675-family index rather than being silently discarded, but cannot be
  selected by the neural query candidate pool.
- Cross-pass query observations prove the current tested queries reach both
  stores; they do not guarantee that every possible query returns a balanced
  mixture of passes or that routing rank is scientifically correct.
