# One-off migrations (historical)

Apply-once / date-stamped scripts moved here 2026-06-08 (KB ops-linkage rec #2) to
keep the live TAG pipeline (notion_to_duckdb, references_ingest, supersession,
tag_query, harvest_ops) uncluttered. These ran once against Notion/duckdb and are
kept for provenance, NOT part of refresh.sh. Do not re-run without checking they
still match the current schema.
