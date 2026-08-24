# Independent ExperimentRecord factory verification

`verify_factory.py` recomputes producer digests and cross-lane invariants without writing any
producer path. It verifies canonical ExperimentRecord schemas against ignored local records, exact
identity uniqueness, source/accession leakage groups, frozen-cohort semantics, raw payload and
authority boundaries, visual receipt/descriptor equality, every referenced local visual object,
and optional final model artifacts plus a clean replay directory.

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python -m \
  corpus.external_training.experiment_factory.verification.verify_factory \
  --root . \
  --model-replay-results data/external_training/experiment_factory/verification/model_replay \
  --output corpus/external_training/experiment_factory/verification/audit.json
```

The replay directory and source payloads remain ignored local data. A clean checkout can run unit
controls and producer digest/privacy checks, but full local ExperimentRecord/object verification
requires the corresponding ignored stores.

