#!/usr/bin/env bash
# Rebuild the TAG backend from source: Notion Contract-Graph (SoT) -> kb.duckdb
# table layer, then references/ PDFs -> content layer (paper_refs + paper_chunks
# + BM25 FTS). Run after Notion changes or new PDFs land. Mirrors the on-demand
# refresh of outputs/obsidian_rag_full/refresh.sh (Notion = source of truth).
set -euo pipefail
cd "$(dirname "$0")"

echo "[1/2] materializing Notion 8-DB Contract-Graph -> kb.duckdb ..."
python notion_to_duckdb.py

echo "[2/2] ingesting references/ PDFs -> paper_refs + paper_chunks + FTS ..."
python references_ingest.py

echo
echo "TAG backend ready. Ask a question:"
echo "  python tag_query.py \"Which Parameters feed a failing ValidationGate?\""
