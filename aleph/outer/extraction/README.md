# OA JATS extraction lane

This directory defines a deterministic, resumable extractor for the successful
Europe PMC OA acquisition receipts.  It does **not** contain article text,
figures, tables, publisher payloads, or authoritative Aleph claims.

The extractor verifies every object against the SHA-256 in its receipt, parses
JATS XML with the Python standard library, and writes retrieval records beneath
the ignored local root `data/external_training/derived/`.  Local chunk records
contain text; committed result manifests contain digests and aggregate counts
only.  Every inferred tag has `authority_status: proposed`.

Run from the repository root:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python \
  corpus/external_training/extraction/extract_jats.py \
  --receipts corpus/external_training/acquisition/europe_pmc_oa_2026-08-05.jsonl \
  --object-root . \
  --derived-root data/external_training/derived
```

Re-running skips an article only when its local per-article completion record
matches the receipt payload hash and extractor version.  `--rebuild` forces a
fresh deterministic extraction.  The final summary hash is over canonical
digest-only article manifests, so two runs over the same successful receipt set
must agree even when file order differs.

This is an import proposal under CAG invariant `INV-IMPORT-1`.  Tag matches are
lexical retrieval aids, not findings, measurements, validation, or ground truth.

`extraction_record.schema.json` describes the ignored local chunk records;
`extraction_summary.schema.json` describes the digest-only committed result.
