# Proposed full-text and table ExperimentRecord extraction

This lane converts the two validated local Europe PMC JATS stores into conservative,
machine-generated `ExperimentRecord` candidates.  It is not a gold dataset and does
not promote literature statements to Aleph physical authority.

Only Methods/Results paragraphs, tables, figures, and explicitly experimental
abstract statements are eligible.  A record requires an explicit number with a
supported unit and an explicit biological-system or assay signal.  Missing systems,
protocols, comparison arms, uncertainty semantics, and replicate types remain
missing or ambiguous.  Background, Introduction, and Discussion prose is excluded.

Full records and bounded evidence text live only under ignored
`data/external_training/experiment_factory/extraction/`.  The committed `results/`
directory contains text-free article and figure join indices, hashes, counts, and an
audit receipt.

Run from the repository root with the Aleph interpreter:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python \
  corpus/external_training/experiment_factory/extraction/extract_experiments.py \
  --receipt first corpus/external_training/acquisition/europe_pmc_oa_2026-08-05.jsonl \
  --receipt second corpus/external_training/acquisition/second_pass/receipts_2026-08-05.jsonl \
  --object-root . \
  --output-root data/external_training/experiment_factory/extraction/full \
  --domain-map corpus/external_training/model/combined_content/splits.jsonl
```

`audit_and_publish.py` validates every local record against the canonical sibling
schema, recomputes record/evidence hashes, checks global ID uniqueness, and writes
the compact committed artifacts.
