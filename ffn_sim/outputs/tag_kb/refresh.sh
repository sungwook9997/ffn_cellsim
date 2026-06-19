#!/usr/bin/env bash
# Rebuild the TAG backend from source: Notion Contract-Graph (SoT) -> kb.duckdb
# table layer, then references/ PDFs -> content layer (paper_refs + paper_chunks
# + BM25 FTS). Run after Notion changes or new PDFs land. Mirrors the on-demand
# refresh of outputs/obsidian_rag_full/refresh.sh (Notion = source of truth).
set -euo pipefail
cd "$(dirname "$0")"

echo "[1/3] materializing Notion 8-DB Contract-Graph -> kb.duckdb ..."
python notion_to_duckdb.py

echo "[2/3] ingesting references/ PDFs -> paper_refs + paper_chunks + FTS ..."
python references_ingest.py

echo "[3/3] supersession / authoritative-chain layer (docs front-matter) ..."
python supersession.py

echo
echo "[sanity] TAG backend summary"
python references_ingest.py --check

echo
echo "[sanity] supersession / authoritative-chain (KU-3.5 test case)"
python supersession.py --check

echo
echo "[sanity] citation integrity (non-destructive; full re-verify = verify_sources.py)"
python verify_sources.py --check

echo
echo "[sanity] results integrity (disk-grounded; rebuilds run_audit table from results_manifest.yaml)"
python verify_runs.py

echo
echo "[sanity] ops-linkage drift (disk runs/modules not yet in the graph)"
python harvest_ops.py --check

echo
echo "TAG backend ready. Ask a question:"
echo "  python tag_query.py \"Which Parameters feed a failing ValidationGate?\""
