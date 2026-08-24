# External experiment projection baseline report

## Result

The final 5,675-source projection was trained end to end with four incremental views and three fixed
seeds. The original grouped split was preserved exactly and the leakage audit found zero violations.
The strongest macro-F1 is the metadata + numeric + visual view (0.3166 ± 0.0036); the full dataset view
has the strongest micro-F1 (0.3761 ± 0.0010) but slightly lower macro-F1 and materially worse ECE.

This is a weak-route routing benchmark only. It does not demonstrate biological or physical prediction.

## Acceptance evidence

| Check | Result |
|---|---|
| Projection binding | SHA-256 `929adbbbd5b5ee7c59cabd5612c3961bd8e4a93edcdaafeb3e3eb654a7f42421` |
| Upstream split preserved | train 4,004 / validation 840 / test 831 |
| Supervised corpus split | train 3,785 / validation 797 / test 787 |
| Frozen cohort | 300, excluded from train/selection/test metrics |
| Unmapped weak labels | 6 excluded (train 5 / validation 1 / test 0) |
| Accession/source groups | 5,560 |
| Leakage violations | 0 |
| Seeds | 11, 29, 47 |
| Deterministic CLI rerun | artifact manifest diff exit 0 with enforced single-thread BLAS |
| Focused tests | 5 passed |

## Ablation

| Available evidence | Macro-F1 | Micro-F1 | Balanced accuracy | ECE |
|---|---:|---:|---:|---:|
| metadata | 0.2846 ± 0.0026 | 0.3427 ± 0.0012 | 0.3231 ± 0.0055 | 0.0437 ± 0.0068 |
| metadata + numeric | 0.3086 ± 0.0023 | 0.3647 ± 0.0018 | 0.3472 ± 0.0019 | 0.0309 ± 0.0028 |
| metadata + numeric + visual | 0.3166 ± 0.0036 | 0.3744 ± 0.0032 | 0.3525 ± 0.0033 | 0.0388 ± 0.0087 |
| metadata + numeric + visual + dataset | 0.3148 ± 0.0007 | 0.3761 ± 0.0010 | 0.3530 ± 0.0007 | 0.0492 ± 0.0031 |

Numeric summaries provide the largest incremental gain. Sparse visual signals provide a smaller gain.
Dataset signals are mixed: micro-F1 rises by 0.0017, macro-F1 falls by 0.0018 and ECE rises by 0.0104
relative to the visual view. No view is promoted as a production model.

## Interpretation boundary

The domain target is the corpus's twelve-class discovery route. Assay and cell tags are correlated with
how retrieval queries and routes were designed. This creates a semantic-leakage/circularity risk that a
source/accession grouping audit cannot detect. Therefore the test scores measure reproducibility within
this corpus construction, not real-world biological generalization.

The frozen 300 candidates are also not gold truth. They are used only to audit representation coverage
and prediction distribution. Model selection from that cohort is explicitly prohibited and no performance
metric is computed for it.

## Artifacts

`results/SHA256SUMS.json` binds configuration, feature schemas, weights, metrics, split audit and frozen
cohort audit. The model is deterministic on the recorded environment, while `results/metrics.json`
retains all per-class supports and confusion matrices needed to expose class imbalance and failure modes.
Independent unrestricted multi-thread BLAS replay exposed weight deltas up to `3.33e-16`; the CLI now
forces one BLAS thread before importing NumPy. Byte identity is claimed only for that controlled CLI path.
