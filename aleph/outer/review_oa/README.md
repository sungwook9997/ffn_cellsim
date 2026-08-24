# OA recoverability screen

**Authority:** `proposed`. This lane ranks manual review; it does not establish article quality,
Tier A status, physical truth, or training authority.

`screen_oa.py` deterministically screens the stable local JATS extraction. It records whether
article-level evidence appears *recoverable*, using exact section/figure/table locators and chunk
digests. A regex or vocabulary hit is only a retrieval aid. It never proves that a method was
adequate or that a result replicated.

Generated full records stay under ignored `data/external_training/review_oa/`. Git contains only:

- a compact per-source index with hashes and counts;
- a ranked manual-review queue;
- duplicate-dataset leakage groups identified from public accession identifiers;
- aggregate coverage, missingness and deterministic input/output hashes.

Decisions are limited to `hold`, `manual_priority`, and `exception_candidate`. Every record remains
`proposed`; there is deliberately no `tier_a` value.

Run from the repository root:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python \
  corpus/external_training/review_oa/screen_oa.py \
  --derived-root data/external_training/derived \
  --metadata-dir corpus/external_training/snapshots \
  --prior-review corpus/external_training/review/article_quality_screen.jsonl \
  --local-root data/external_training/review_oa \
  --output-root corpus/external_training/review_oa/results
```

The script refuses an extraction whose `final_summary.json` reports failures or whose manifest
count/hash does not reproduce. Re-running the same input must reproduce the same committed artifact
hashes.
