# External ExperimentRecord dataset factory — integrated report

Date: 2026-08-05 KST

Branch: `codex/external-training-corpus`

Authority: `proposed` only

## Outcome

The earlier paper-level retrieval corpus is now a reproducible experiment-data factory. It emits a
versioned ExperimentRecord contract, a stratified annotation queue, exact-locator text/table/numeric
candidates, official OA visual assets and descriptors, public-dataset manifests, a text-free learning
projection, and an external CPU routing baseline. It does not modify `aleph/**`, infer physical
parameters, or promote machine candidates to scientific truth.

| Surface | Measured result |
|---|---:|
| OA JATS sources processed | 5,675 / 5,675 |
| ExperimentRecord candidates | 52,725 |
| numeric observations | 238,489 |
| schema errors / observation-ID collisions | 0 / 0 |
| table / figure-linked records | 1,664 / 14,302 |
| frozen annotation candidates | 300 papers, 600 blinded A/B slots |
| public source-accession bindings | 1,982 |
| unique dataset manifests / providers | 1,801 / 18 |
| official metadata pilot | 81 resolved of 100 bounded requests |
| licensed raw-data pilot | 9 files, 3,487,812 bytes, 0 Git payloads |
| indexed figures | 38,316 |
| selected visual graphic identities | 3,315 |
| acquired graphic refs / unique objects | 627 / 480 |
| acquired visual bytes / sources | 174,796,895 / 116 |
| visual object digest failures | 0 |
| text-free learning sources / leakage groups | 5,675 / 5,560 |
| train / validation / test | 4,004 / 840 / 831 |
| split group/accession/gold-component violations | 0 / 0 / 0 |
| numeric candidate targets / masked observations | 39,344 / 199,145 |
| raw-text, caption or pixel fields copied to projection | 0 |

## Data states and human-review boundary

Automatic outputs remain `machine_candidate`. The frozen 300-paper set is `gold_candidate`, which
means selected for annotation rather than verified. Two independent annotations plus adjudication
are required to produce an `adjudicated` record; that human step has not been performed here. The
factory therefore supplies training infrastructure and weak supervision, not a biologically stable
external model or an Aleph parameter prior.

The 300-paper candidate set spans all 12 discovery domains, all 12 observed modality tags, all years
2010–2026, and 11 cell-state tags. Five multi-paper public-accession components containing 26 papers
were retained intact. The sparse modalities called out for deliberate annotation include magnetic
tweezers, optical tweezers, micropipette aspiration, PIV, single-cell sequencing, and TFM.

## Extraction controls

The final case-sensitive unit rules distinguish `nM` concentration, `nm` length, uppercase force
units, and sample-size notation. Regression controls removed `18S`, `W748S`, and bare `n` false
positives. Drug-context micrometre-like tokens that cannot be safely disambiguated are not guessed:
184 are retained as `length_or_concentration` ambiguity. Relative to the first rule set, the corrected
extractor removed 10,991 candidate records, 30,162 observations, and 9,734 ambiguous records.

All 52,725 records pass the canonical schema. Numeric observations without adequate context or with
ambiguity remain masked; zero ambiguous observations are promoted to model labels.

## Visual and public-dataset evidence

Only official PMC-page-listed NCBI CDN assets were attempted. The final visual ledger has one canonical
receipt per graphic identity, matching descriptor IDs exactly. Of 3,315 identities, 627 were described;
the others retain typed unresolved or deferred states. Descriptors cover technical image properties and
panel candidates only. The pipeline explicitly extracts zero scientific values from pixels.

The public dataset registry covers GEO, SRA, BioProject, ArrayExpress, PRIDE, Zenodo, Figshare, Dryad,
OSF, BioImage Archive and other providers. Accession co-membership is treated as a leakage boundary,
not proof that samples or processing are identical. The raw-data pilot downloaded only explicitly open,
size-bounded official files: TFM MAT, nuclear-mechanics XLSX, IF CZI, and PIV ZIP. The content-addressed
payloads remain ignored local data; none are committed.

## Weak-route baseline

The model target is the existing 12-class discovery route, not a cell state, experiment outcome, force,
physical parameter, evidence-quality score, or truth label. The frozen 300 candidates are excluded from
training, validation selection, and test metrics. Six `unmapped` sources are also excluded, leaving a
supervised split of 3,785 / 797 / 787. Values below are mean ± population standard deviation across seeds
11, 29, and 47 on the same 787-source test set.

| Evidence view | Macro-F1 | Micro-F1 | Balanced accuracy | ECE |
|---|---:|---:|---:|---:|
| metadata | 0.2846 ± 0.0026 | 0.3427 ± 0.0012 | 0.3231 ± 0.0055 | 0.0437 ± 0.0068 |
| + numeric | 0.3086 ± 0.0023 | 0.3647 ± 0.0018 | 0.3472 ± 0.0019 | 0.0309 ± 0.0028 |
| + visual | 0.3166 ± 0.0036 | 0.3744 ± 0.0032 | 0.3525 ± 0.0033 | 0.0388 ± 0.0087 |
| + dataset manifest | 0.3148 ± 0.0007 | 0.3761 ± 0.0010 | 0.3530 ± 0.0007 | 0.0492 ± 0.0031 |

Numeric summaries give the largest incremental improvement. Sparse visual evidence gives a smaller
improvement. Dataset features are mixed: micro-F1 rises by 0.0017, while macro-F1 falls by 0.0018 and
ECE worsens by 0.0104 relative to the visual view. No view is promoted as a production model.

Assay/cell tags and the discovery-route target are partly products of the same retrieval/query design.
Source/accession grouping cannot remove that semantic circularity. These scores demonstrate a stable
projection and a reproducible within-corpus routing benchmark, not real-world biological generalization.

## Reproducibility defect found and corrected

Independent replay initially found multithreaded-BLAS weight differences of `1.25e-16`–`3.33e-16`.
Epochs and metrics were unchanged, but the byte-determinism claim was false. The CLI now fixes five BLAS
thread controls to one before NumPy import. Two fresh launches and the independent verifier reproduce all
six model artifacts byte-for-byte. Unrestricted multithreaded or already-initialized BLAS remains outside
that guarantee.

## Independent acceptance

The independent verifier recomputed all 22 producer-manifest hashes, validated 300 annotation candidates
and 52,725 ExperimentRecords, hashed all 480 referenced visual objects, checked 1,982 accession bindings,
and replayed the model. It observed 99,145 `authority_status` values, all `proposed`; forbidden raw/pixel
keys and `decided_by` occurrences were both zero. Final audit status is PASS, its JSON is byte-identical
across two runs, and the factory-focused suite is **48 passed**.

This PASS establishes artifact integrity and contract compliance. It does not human-validate extracted
measurements, certify inferred panels, create an adjudicated gold set, or authorize scientific use.

## Commit chain

- `495980c` — ExperimentRecord schema and 300-paper annotation candidates
- `929cafc` — full OA JATS experiment/numeric extraction
- `834c94e` — public dataset registry and licensed pilot
- `31ef2c3` — visual index, OA assets, descriptors, and joins
- `d0e66f5` — text-free leakage-safe learning projection
- `d827962` — four-view weak-route CPU baseline
- `ff666a3` — independent end-to-end verification

Detailed reports live in each producer directory and in `verification/REPORT.md`.
