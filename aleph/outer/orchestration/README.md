# Overnight orchestration

This directory contains only integration audits and the final measured report for the autonomous external-corpus run. It owns no source truth and promotes nothing beyond `proposed`.

The run has independent, append-only producers:

1. official Europe PMC OA acquisition into an ignored content-addressed local store;
2. an isolated RAG ledger, TAG snapshot, and CAG lint adapter using Aleph harness public APIs;
3. resumable JATS extraction and exact-locator retrieval indexes;
4. conservative recoverability screening with no automatic evidence promotion;
5. CPU-only external weak-supervision tagger/retriever baselines.

`audit_run.py` folds all acquisition passes and the committed producer summaries into one count
without reading publisher credentials or committing full text. It normalises source-family IDs,
merges duplicate metadata, distinguishes candidate files from acquisition receipts, and reports
missing lane outputs explicitly. `test_audit_run.py` covers case-normalised deduplication and
multi-pass receipt aggregation.
