# External projection weak-route CPU baseline

This directory trains a deterministic NumPy baseline for one target only: the existing twelve-class
weak discovery route attached to each source projection. It does not predict biology, experimental
outcomes, physical parameters, evidence quality, or truth.

Four incremental views are compared on the exact same source/accession-grouped split:

1. metadata excluding the route label itself;
2. metadata + numeric availability/count summaries and ambiguity masks;
3. the above + visual metadata, asset availability and source-level numeric descriptors;
4. the above + dataset-manifest availability/count/modality hints.

All category columns and standardization statistics are fitted on the corpus-only training split.
The 300 `gold_candidate` sources are excluded from training, validation model selection and test
metrics. They receive a coverage and prediction-distribution audit only.

Run with the Aleph interpreter:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python \
  corpus/external_training/experiment_factory/model/train.py \
  --projection corpus/external_training/experiment_factory/learning/results/source_projection.jsonl \
  --projection-summary corpus/external_training/experiment_factory/learning/results/summary.json \
  --output-dir corpus/external_training/experiment_factory/model/results
```

The model is a class-weighted linear softmax classifier trained with full-batch Adam. Seeds
11/29/47 are all reported. Early stopping observes validation macro-F1 only; test data never selects
an epoch, seed, feature view or hyperparameter. The weight archive has fixed ZIP metadata so a full
rerun is byte reproducible.

This benchmark is a pipeline and representation smoke test, not a biological generalization claim.
The weak target comes from discovery routes, while assay and cell tags are also influenced by the
retrieval/query design. Their correlation can therefore create semantic leakage or circularity even
when source/accession split leakage is exactly zero. See `MODEL_CARD.md` and `REPORT.md`.
