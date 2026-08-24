# External metadata model experiment

This directory is an explicitly **non-authoritative** weak-supervision experiment over the external
literature catalogue. It is outside `aleph/learn` and `aleph/represent`, does not import either
package, does not read simulation trajectories, and cannot fit or validate a physical parameter.
It therefore does not unseal MTG-PN or satisfy gate R5.

The experiment asks two deliberately limited questions:

1. can title, journal and publication-year metadata recover the discovery-domain label attached by
   the corpus search process; and
2. does the learned hidden representation retrieve other records carrying the same weak label?

The label is search provenance, not biological truth. Metrics measure catalogue organisation, not
physics, article quality, or mechanistic validity.

Run on CPU with the Aleph interpreter:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python corpus/external_training/model/train.py \
  --config corpus/external_training/model/config.json
```

The trainer namespace- and case-normalises source-family identifiers, then deduplicates all five
snapshots plus the seed manifest before splitting.
It uses a deterministic SHA-256 group assignment, so no source family can cross train, validation,
or test. Outputs include the exact split manifest, NumPy weights, metrics, and SHA-256 receipts.

Query the frozen external index:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python corpus/external_training/model/infer.py \
  --query "actin cortical tension membrane pressure" --top-k 5
```

The API is `infer.ExternalMetadataIndex.query`. It verifies the frozen weights and split receipts,
reconstructs the case-normalised catalogue, and returns weak domain probabilities plus cosine-ranked
source families. Empty, catalogue-unknown, physical-parameter-seeking, and authority-seeking queries
raise the typed `QueryRefusal`. Example inputs and frozen machine-readable outputs are in
`example_queries.json` and `example_results.json`.
