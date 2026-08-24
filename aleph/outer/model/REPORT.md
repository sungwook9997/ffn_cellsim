# CPU weak-supervision baseline report

## Boundary and verdict

This run trained an external metadata organiser, not an Aleph physics model. It uses only title,
journal, publication year and corpus discovery labels. It neither imports nor changes
`aleph/learn`, `aleph/represent`, physics, guards, or the RAG–TAG–CAG lane. It does **not** unseal
MTG-PN, satisfy R5, validate an article, or train a physical parameter.

Within that narrow boundary, the learned model improves on both simple tag and retrieval baselines.
The result supports using a neural embedding to organise review queues; it is not evidence that the
embedding has learned cell mechanics.

## Dataset and leakage control

| Item | Measured value |
|---|---:|
| Snapshot rows | 7,200 |
| Seed rows | 18 |
| Seed/snapshot overlaps replaced by seed metadata | 5 |
| Unique source families after namespace/id case normalisation and deduplication | **7,213** |
| Weak labels | 29 |
| Discovery-query labels | 5,396 |
| Historical-title-rule labels | 1,800 |
| Seed-topic-rule labels | 18 |
| Train / validation / test | **5,025 / 1,042 / 1,146** |

The split unit is the namespace- and identifier-lowercased `source_family_id`. A seeded SHA-256
mapping assigns a family exactly once, before feature construction or training. The committed split
receipt contains 7,213 unique IDs and no duplicate family. This specifically collapses seed
`doi:10.7554/eLife.72381` and snapshot `doi:10.7554/elife.72381`. Temporal search lanes such as
`2010s_epithelial_states` are normalised to their
semantic domain; the 2010–2021 backfill is labelled by explicit title rules rather than publication
year. These choices prevent the model from receiving a synthetic year-as-class target, but they do
not turn weak labels into reviewed ground truth.

## Model and measured performance

The neural model is a one-hidden-layer ReLU MLP over a 4,096-bucket stable hashed bag of title and
journal tokens plus a five-year bucket. It has 64 hidden units, uses Adam and early stopping, and was
trained CPU-only with NumPy. It stopped after 11 of 35 allowed epochs.

| Test metric | Simple baseline | Neural model |
|---|---:|---:|
| Tag macro-F1 | majority: 0.0053 | **0.6055** |
| Tag macro-F1 | nearest centroid: 0.5053 | **0.6055** |
| Tag accuracy | — | 0.6396 |
| Retrieval hit-rate@5 | fixed hashed features: 0.6702 | **0.7618** |
| Retrieval precision@5 | fixed hashed features: 0.2894 | **0.6075** |

Retrieval was measured for all 1,146 test records against the 5,025-record training corpus. A hit
means at least one of the five neighbours shares the same weak discovery-domain label.

Per-domain F1 ranges from 0.3636 (`tissue_organoid_mechanics`, n=23 test) to 0.8793
(`neural_glial_states`, n=61 test). The full 29-domain table and supports are in `metrics.json`.
Weak domains with low F1 should remain review queues, not be auto-promoted.

## Reproducibility receipts

Two independent consecutive CPU runs produced the identical compressed weight digest:

`e0175c8ebcbbf8ada8d64c309bf1bdc9bd67024a48f58a1f0208c16f725b837b`

The deterministic split digest is:

`41e0598783a11d2fba7d2df65fcafb47e6508d93ad42d37cb4d8cf2b7efa8b94`

The two measured runs took 9.06 and 17.83 seconds on this shared CPU host. `SHA256SUMS.json` records the committed
configuration, source, weights, splits and metrics hashes. Runtime duration itself naturally changes
between runs; the weights and split are the reproducibility invariants.

## Appropriate next use

The embeddings may rank and diversify article-review batches or provide a non-authoritative
retrieval candidate list to the isolated corpus adapter. They must not select physical laws, invent
missing observations, promote evidence authority, or enter a released model bundle without a new
reviewed interface and the relevant Aleph gate decisions.

## Frozen inference surface

`infer.py` makes the measured baseline usable without retraining. It verifies the committed weights
and split digests, reconstructs the 7,213-family canonical catalogue, and exposes
`ExternalMetadataIndex.query(text, top_k)`. Results contain weak-domain probabilities and
cosine-ranked source-family/title metadata with their original evaluation split visible.

The CLI and API refuse four unsupported cases: empty queries, queries with no token in the frozen
catalogue, physical-parameter estimation requests, and authority/ground-truth claims. These are
typed `QueryRefusal` results rather than low-confidence fabricated answers. Twelve focused controls
cover determinism, digest validation through construction, identifier uniqueness, example-result
parity, refusals, and the absence of sealed-package imports.

The five-query machine-readable suite covers cortical tension, TFM, nuclear mechanics,
mechano-osmotic volume, and mixed cell-state routing. The last query routes too strongly toward
mechanotransduction; `example_results.json` records that failure visibly as a weak-label limitation.
