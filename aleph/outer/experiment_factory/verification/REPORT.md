# Independent ExperimentRecord factory verification report

Status: **PASS — AGENT-PROPOSED verification, no authority promotion**  
Date: 2026-08-05 KST

## Result

The complete external ExperimentRecord factory passed the independent structural, provenance,
privacy, leakage, visual-object and model-replay audit. The verifier writes only this directory and
does not import Aleph physics, learning authority, or producer mutation APIs.

| Surface | Independently observed |
|---|---:|
| frozen annotation candidates | 300 unique sources |
| local ExperimentRecords schema-validated | 52,725 |
| candidate / ExperimentRecord schema errors | 0 / 0 |
| extraction article / figure-join rows | 5,675 / 14,302 |
| source-accession bindings / dataset manifests | 1,982 / 1,801 |
| learning projection sources / leakage groups | 5,675 / 5,560 |
| cross-split group / accession / gold-component violations | 0 / 0 / 0 |
| visual figures / receipts / descriptors | 38,316 / 3,315 / 3,315 |
| visual duplicate IDs / receipt-descriptor mismatch | 0 / 0 |
| referenced visual objects fully hashed | 480 |
| missing / corrupt visual objects | 0 / 0 |
| tracked pilot/raw payloads | 0 |
| forbidden raw-text/pixel keys / `decided_by` | 0 / 0 |
| observed authority values | 99,145 `proposed`, no other value |

All 22 files named by producer SHA manifests reproduce their committed digests. Extraction and
learning digest bindings were separately recomputed rather than accepted from their summaries.

## Frozen 300 and leakage

The candidate set is exactly the frozen SHA-bound 300-source cohort. Its projection marker is
`frozen_300_gold_candidate_not_ground_truth`; the set is **not** treated as verified labels. Every
row keeps `gold_candidate_is_ground_truth=false`, and the model computes no performance metric on
the cohort or uses it for training, validation selection, or test scoring. Its ordinary split field
is retained only for future group-safe annotation workflows.

For the actual weak-route baseline, all 5,560 unioned source/accession components remain inside one
split. Independently checking the 1,982 source-accession rows and the candidate leakage components
found no cross-split component.

## Model audit and corrected reproducibility defect

The final model is a CPU NumPy linear baseline for the existing weak twelve-class discovery route,
not a biological outcome predictor, parameter estimator, or Aleph runtime model. It trains on the
5,375 non-cohort sources only, with train/validation/test counts 3,785 / 797 / 787 and seeds
11, 29, 47. All four views use the same split and test support; test and candidate-cohort results do
not select the model.

Mean test macro-F1 across the three seeds:

| View | Macro-F1 mean | std |
|---|---:|---:|
| metadata | 0.2846 | 0.0026 |
| metadata + numeric | 0.3086 | 0.0023 |
| metadata + numeric + visual | 0.3166 | 0.0036 |
| metadata + numeric + visual + dataset | 0.3148 | 0.0007 |

The first independent full replay found a real reproducibility defect: multithreaded BLAS changed
192–754 elements in each weight array by at most `3.33e-16`. Epochs and macro-F1 were unchanged,
but weights, metrics and cohort-audit bytes were not identical, so the initial byte-determinism
claim failed. The producer then forced OpenBLAS, OMP, MKL, Accelerate and NumExpr to one thread
before importing NumPy and regenerated the artifacts. Two fresh controlled launches and the final
independent replay now match all six artifact hashes byte-for-byte. The final audit records zero
model validation errors and an empty replay mismatch list. Unrestricted multithread BLAS remains
outside that byte-level guarantee.

## Data and authority boundaries

The Git tree contains no dataset pilot payload, PDF, XML, JATS, image, video, NumPy raw array, HDF5
or publisher payload under the factory. Structured projections were recursively scanned for raw
evidence/caption/full-text/pixel keys. None were present. Model weights are permitted model
artifacts and are not confused with raw datasets.

Every inspected `authority_status` is `proposed`; no `decided_by` field exists. Ambiguous
observations and candidate status are never promoted to labels. Visual metadata, acquired-asset
receipts and actual decoded descriptors remain distinct masks.

## Controls and limits

The full factory-focused test run includes schema adversaries, deterministic IDs, gold selection,
numeric false-positive controls, dataset registry and payload guards, visual host/panel/hash
adversaries, learning leakage/mask controls, model split/replay controls, and this verifier's own
recursive privacy/identity controls.

Final focused run: **48 passed**. Two complete verifier executions emitted byte-identical audit JSON
with SHA-256 `9277fb02adf2866b320d469075de5b7579d190f561d8035468c0a648ab84a1a6`.

This PASS establishes artifact integrity and declared-contract compliance. It does not validate
weak discovery labels as biological truth, certify image panel boxes, validate extracted numeric
observations by human annotation, or promote the baseline to Aleph model authority.
