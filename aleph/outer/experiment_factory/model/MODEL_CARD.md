# Weak discovery-route projection baseline — model card

## Intended use

This proposed baseline checks whether the external-training projection can be converted into stable,
split-safe numeric matrices and whether successive evidence views carry information about the existing
twelve broad discovery routes. It may support corpus routing, retrieval diagnostics and future annotation
planning. It is not an Aleph physics model and is not approved for scientific inference.

The target is the source's weak discovery-route label. It is not a biological state, experimental
outcome, force, parameter, evidence-quality score or ground truth. Predictions must not be presented as
any of those quantities.

## Data and split contract

- Input: 5,675 records from `learning/results/source_projection.jsonl`, SHA-256
  `929adbbbd5b5ee7c59cabd5612c3961bd8e4a93edcdaafeb3e3eb654a7f42421`.
- The upstream source/accession-grouped assignments are preserved without reassignment: 4,004 train,
  840 validation and 831 test across all cohorts.
- The frozen 300 `gold_candidate` records are excluded from fitting, early stopping and test metrics.
- Six corpus records with an `unmapped` weak route are excluded: five train and one validation.
- Supervised counts are therefore 3,785 train, 797 validation and 787 test.
- Group leakage audit: 5,560 groups, zero split violations.

All vocabularies, feature columns, means and scales are learned from the corpus-only training split.
The same held-out test records are used for every feature view and seed. Validation macro-F1 alone
selects the stopping epoch; the test set and frozen cohort select nothing.

## Model and inputs

The model is a full-batch NumPy linear softmax classifier with inverse-frequency class weights, Adam,
L2 regularization and three fixed seeds (11, 29, 47). The four nested views are metadata; metadata plus
numeric summaries; those plus visual availability/descriptors; and those plus dataset-manifest signals.
Explicit missingness features accompany sparse evidence views. The target field, source identifiers and
artifact digests are excluded from model inputs.

## Performance

Values are mean ± population standard deviation across the three seeds on the same 787-record test set.

| View | Features incl. intercept | Macro-F1 | Micro-F1 | Balanced accuracy | ECE (10 bins) |
|---|---:|---:|---:|---:|---:|
| metadata | 27 | 0.2846 ± 0.0026 | 0.3427 ± 0.0012 | 0.3231 ± 0.0055 | 0.0437 ± 0.0068 |
| + numeric | 45 | 0.3086 ± 0.0023 | 0.3647 ± 0.0018 | 0.3472 ± 0.0019 | 0.0309 ± 0.0028 |
| + visual | 60 | 0.3166 ± 0.0036 | 0.3744 ± 0.0032 | 0.3525 ± 0.0033 | 0.0388 ± 0.0087 |
| + dataset manifest | 81 | 0.3148 ± 0.0007 | 0.3761 ± 0.0010 | 0.3530 ± 0.0007 | 0.0492 ± 0.0031 |

Adding dataset-manifest features does not improve macro-F1 and worsens ECE relative to the visual view;
it only yields a small micro-F1 increase. This is reported as a negative/mixed ablation, not an optimization
success. Per-seed metrics, per-class support/precision/recall/F1 and full confusion matrices are retained
in `results/metrics.json`.

## Frozen-cohort audit

The 300 frozen candidates receive coverage and prediction-distribution diagnostics only. Coverage is
298 numeric, 300 visual metadata, 116 visual asset receipts/descriptors and 98 dataset manifests. No
accuracy, F1 or calibration value is computed because `gold_candidate` identifies a sampling cohort,
not trusted labels.

## Limitations and prohibited interpretations

Zero source/accession split leakage does not imply semantic independence. The discovery-route target,
assay tags and cell tags all arise partly from the same retrieval/query design. A classifier can exploit
that construction and appear predictive through semantic leakage or circularity. Consequently these
numbers must not be interpreted as generalization to biology, experiments, cell states or physical laws.

The labels are weak, class support is imbalanced, descriptors exist for only 116 of 5,675 sources, and
dataset manifests for 885. A linear model cannot establish causal or mechanistic relations. External,
independently annotated evaluation is required before scientific use. The model has `proposed` authority
and must not update Aleph parameter priors or physics oracles.

## Reproducibility

The CLI fixes BLAS thread counts to one before NumPy initializes. Multi-threaded BLAS was independently
observed to change individual learned weights by up to `3.33e-16`, without changing macro-F1 or stopping
epochs, so unrestricted BLAS is not claimed to be byte deterministic. Under the enforced CLI condition,
two clean runs produced the same `SHA256SUMS.json`; deterministic archive timestamps make `weights.npz`
byte-identical. Importing this module after another library has already initialized BLAS is outside that
byte-reproducibility guarantee. The artifact manifest and exact configuration are in `results/`.
