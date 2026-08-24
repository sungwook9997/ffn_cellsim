#!/usr/bin/env bash
# Rebuild the TAG backend from source: Notion Contract-Graph (SoT) -> kb.duckdb
# table layer, then references/ PDFs -> content layer (paper_refs + paper_chunks
# + BM25 FTS). Run after Notion changes or new PDFs land. Mirrors the on-demand
# refresh of outputs/obsidian_rag_full/refresh.sh (Notion = source of truth).
set -euo pipefail
cd "$(dirname "$0")"

# `python` is NOT on PATH in this project's environments — the conda env exposes its own bin dir and
# there is no bare `python` shim.  Every line below used it anyway, so this script died on its FIRST
# command with "python: command not found" — while `CLAUDE.md` tells every session to run it to check
# KB drift.  Same defect, same fix, as the Makefile and the pre-commit hook (2026-07-29).
#
# The predicate is "can this interpreter RUN the pipeline", not "does a python exist": resolving
# /usr/bin/python3 would trade a PATH error for a ModuleNotFoundError, which is the same failure in a
# different hat.  duckdb is the discriminating import here — it is what the first stage needs.
# Override with `PY=/path/to/python bash refresh.sh`.
if [ -z "${PY:-}" ]; then
  for c in "${CONDA_PREFIX:-}/bin/python" "$HOME/miniconda3/envs/ffn_sim/bin/python" \
           "$(command -v python 2>/dev/null || true)" "$(command -v python3 2>/dev/null || true)"; do
    [ -n "$c" ] && [ -x "$c" ] && "$c" -c 'import duckdb' >/dev/null 2>&1 && { PY="$c"; break; }
  done
fi
if [ -z "${PY:-}" ]; then
  echo "no interpreter on this machine can import duckdb — activate the env (conda activate ffn_sim)" >&2
  echo "or pass one explicitly: PY=/path/to/python bash refresh.sh" >&2
  exit 3
fi
echo "[env] interpreter: $PY"

echo "[1/3] materializing Notion 9-DB Contract-Graph -> kb.duckdb ..."
"$PY" notion_to_duckdb.py

echo "[1b] committing a git-durable snapshot of the SoT -> snapshots/notion_snapshot.json ..."
"$PY" dump_notion_snapshot.py

echo "[2/3] ingesting references/ PDFs -> paper_refs + paper_chunks + FTS ..."
"$PY" references_ingest.py

echo "[3/3] supersession / authoritative-chain layer (docs front-matter) ..."
"$PY" supersession.py

echo
echo "[sanity] TAG backend summary"
"$PY" references_ingest.py --check

echo
echo "[sanity] supersession / authoritative-chain (KU-3.5 test case)"
"$PY" supersession.py --check

echo
echo "[sanity] citation integrity (non-destructive; full re-verify = verify_sources.py)"
"$PY" verify_sources.py --check

echo
echo "[sanity] results integrity (disk-grounded; rebuilds run_audit table from results_manifest.yaml)"
"$PY" verify_runs.py

echo
echo "[sanity] parameter provenance (disk + citation-grounded; rebuilds param_audit from params_manifest.yaml)"
"$PY" verify_params.py

echo
echo "[sanity] ops-linkage drift (disk runs/modules not yet in the graph)"
"$PY" harvest_ops.py --check

echo
echo "TAG backend ready. Ask a question:"
echo "  python tag_query.py \"Which Parameters feed a failing ValidationGate?\""
