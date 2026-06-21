#!/usr/bin/env python3
"""Dump the Notion Contract-Graph (8 DBs) to a committed JSON snapshot.

The durability fix (KB-audit 2026-06-18, finding #2): the knowledge layer was
reproducible only from a live Notion workspace behind a secret token — gitignored
duckdb, no offline copy, no diff/blame on the knowledge itself. This writes a
flattened, text-diffable snapshot of every row in all 8 databases into git, so the
SoT is recoverable from the repo alone (fresh clone / token rotated / Notion down)
and every change to the knowledge shows up in git history.

One file, `snapshots/notion_snapshot.json` (git history is the dated archive). Uses
the same Notion plumbing as notion_to_duckdb (token from
../obsidian_rag_full/.notion_token or $NOTION_TOKEN). Read-only against Notion.

RUN:
  conda activate ffn_sim
  python dump_notion_snapshot.py
"""
from __future__ import annotations

import datetime
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).parent
OBS = HERE.parent / "obsidian_rag_full"
sys.path.insert(0, str(OBS))
from notion_to_obsidian import (  # noqa: E402
    DATA_SOURCES, get_token, plain, query_all,
)

SNAP_DIR = HERE / "snapshots"
OUT = SNAP_DIR / "notion_snapshot.json"


def main() -> None:
    tok = get_token()
    SNAP_DIR.mkdir(exist_ok=True)
    out = {
        "source": "Notion Contract-Graph (ffn_cellsim KB) — single source of truth",
        "notion_version": "2022-06-28",
        "built_at": datetime.date.today().isoformat(),
        "note": "Flattened snapshot for durability + diff. Regenerate with "
                "dump_notion_snapshot.py; git history is the dated archive.",
        "row_counts": {},
        "databases": {},
    }
    for ds, dsid in DATA_SOURCES.items():
        rows = query_all(dsid, tok)
        recs = []
        for pg in rows:
            rec = {"_id": pg["id"].replace("-", "")}
            for pname, prop in pg["properties"].items():
                rec[pname] = plain(prop)
            recs.append(rec)
        # stable order so git diffs are meaningful
        recs.sort(key=lambda r: str(r.get("_id", "")))
        out["databases"][ds] = recs
        out["row_counts"][ds] = len(recs)
        print(f"  {ds:18s} {len(recs)} rows")
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=1))
    total = sum(out["row_counts"].values())
    print(f"\nwrote {OUT.relative_to(HERE.parents[2])} "
          f"({OUT.stat().st_size // 1024} KB, {total} rows across "
          f"{len(out['databases'])} databases)")


if __name__ == "__main__":
    main()
