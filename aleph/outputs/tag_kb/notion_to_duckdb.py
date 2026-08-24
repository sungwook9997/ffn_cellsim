#!/usr/bin/env python3
"""Materialize the ffn_cellsim Notion Contract-Graph (9 DBs) into a local DuckDB.

This is the **table layer** for the TAG (Table-Augmented Generation) engine
(`tag_query.py`). Notion = source of truth; `kb.duckdb` is a regenerable local
mirror that `syn`/`exec` run SQL against (no live-Notion querying, no rate limits).

Produces `kb.duckdb` with:
  • 9 node tables  — source_evidence, knowledge_claim, model_contract, parameter,
      validation_gate, code_mapping, run_result, decision_ledger, cell_state
      (one row per Notion row; columns = flattened scalar properties; relation
       columns kept as a JSON list of target node-ids for convenience)
  • edges          — long-form relation triples (src_id, src_type, rel, dst_id,
      dst_type), de-duplicated across Notion's dual two-way relations
  • _meta          — (table, n_rows, built_at) provenance

The references + paper_chunks content tables are added by `references_ingest.py`.

Reuses the verified Notion API plumbing from the Obsidian mirror
(`../obsidian_rag_full/notion_to_obsidian.py`): same token, same DATA_SOURCES
(database ids), same `query_all`/`plain` flatteners — single source of truth so
the two mirrors never drift.

RUN:
  conda activate ffn_sim
  python notion_to_duckdb.py            # -> kb.duckdb
  python notion_to_duckdb.py --stats    # print row counts + a sample join
"""
from __future__ import annotations

import json
import pathlib
import re
import sys

import duckdb

# --- reuse the Obsidian mirror's verified Notion plumbing (DRY) --------------
HERE = pathlib.Path(__file__).parent
OBS = HERE.parent / "obsidian_rag_full"
sys.path.insert(0, str(OBS))
from notion_to_obsidian import (  # noqa: E402  (path-dependent import)
    DATA_SOURCES,   # ds_name -> Notion database id
    get_token,      # reads ../obsidian_rag_full/.notion_token
    plain,          # flatten one Notion property to str | list[id]
    query_all,      # paginated /databases/{id}/query
)

DB_PATH = HERE / "kb.duckdb"

# Notion DB name -> SQL table name (snake_case)
TABLE = {
    "SourceEvidence": "source_evidence",
    "KnowledgeClaim": "knowledge_claim",
    "ModelContract":  "model_contract",
    "Parameter":      "parameter",
    "ValidationGate": "validation_gate",
    "CodeMapping":    "code_mapping",
    "RunResult":      "run_result",
    "DecisionLedger": "decision_ledger",
    # The 9th database (PI 2026-07-28). The (cell type x cell state) axis the graph had no way to
    # represent, which matters because stage 5 DEFINES project progress as
    # `cell types x states x parameters recovered` — a denominator that could not previously be
    # written down, let alone counted.
    "CellState":      "cell_state",
}


def col(name: str) -> str:
    """Notion property name -> snake_case SQL identifier."""
    out = []
    for ch in name.strip():
        out.append(ch.lower() if (ch.isalnum() or ch == "_") else "_")
    s = "".join(out)
    while "__" in s:
        s = s.replace("__", "_")
    s = s.strip("_") or "col"
    if s[0].isdigit():
        s = "c_" + s
    return s


def build():
    tok = get_token()

    # ---- pass 1: pull every row; build id -> node type (for edge dst typing) -
    raw: dict[str, list[dict]] = {}
    id2type: dict[str, str] = {}
    for ds, dsid in DATA_SOURCES.items():
        print(f"reading {ds} ...", flush=True)
        rows = query_all(dsid, tok)
        raw[ds] = rows
        for pg in rows:
            id2type[pg["id"].replace("-", "")] = TABLE[ds]
        print(f"  {len(rows)} rows", flush=True)

    # non-destructive: only the 9 node tables + edges + _meta are rebuilt
    # (CREATE OR REPLACE below); content/audit tables (paper_refs, paper_chunks,
    # source_audit) added by sibling scripts survive a refresh.
    con = duckdb.connect(str(DB_PATH))

    edges: list[tuple] = []
    seen_edges: set[tuple] = set()
    meta: list[tuple] = []

    # ---- pass 2: one table per DB + collect relation edges -------------------
    for ds, rows in raw.items():
        tname = TABLE[ds]
        records: list[dict] = []
        colset: list[str] = ["id"]          # ordered, stable
        for pg in rows:
            pid = pg["id"].replace("-", "")
            rec: dict[str, object] = {"id": pid}
            for pname, prop in pg["properties"].items():
                val = plain(prop)
                cname = col(pname)
                if prop["type"] == "relation":
                    # node table keeps a JSON list; normalized form -> edges
                    targets = [t for t in val if t in id2type]
                    rec[cname] = json.dumps(targets) if targets else None
                    for t in targets:
                        e = tuple(sorted((pid, t))) + (cname,)  # undirected dedupe key
                        if e in seen_edges:
                            continue
                        seen_edges.add(e)
                        edges.append((pid, tname, cname, t, id2type[t]))
                elif prop["type"] == "title":
                    rec["citation_key" if ds == "SourceEvidence" else "title"] = val
                    cname = "citation_key" if ds == "SourceEvidence" else "title"
                else:
                    rec[cname] = val if val != "" else None
                if cname not in colset:
                    colset.append(cname)
            records.append(rec)

        # CREATE TABLE (all TEXT — Notion props flatten to strings; CAST in SQL
        # when a query needs a number; keeps the schema robust to prop drift)
        cols_ddl = ", ".join(f'"{c}" TEXT' for c in colset)
        con.execute(f'CREATE OR REPLACE TABLE "{tname}" ({cols_ddl})')
        placeholders = ", ".join("?" for _ in colset)
        con.executemany(
            f'INSERT INTO "{tname}" VALUES ({placeholders})',
            [[r.get(c) for c in colset] for r in records],
        )
        meta.append((tname, len(records)))
        print(f"  wrote table {tname}: {len(records)} rows, {len(colset)} cols", flush=True)

    # ---- derived citation edges ----------------------------------------------
    # A KnowledgeClaim that names a paper's DOI (or a unique Surname+Year) in its
    # prose links to that SourceEvidence even when the Notion relation was never
    # wired. Faithful — the citing node itself states the reference. Skips pairs
    # already connected by a real relation so counts are not inflated.
    _norm_doi = lambda d: re.sub(r'^https?://(dx\.)?doi\.org/', '',
                                 (d or '').strip().lower()).rstrip('.')
    _DOI = re.compile(r'10\.\d{4,9}/[^\s;)\]]+')
    existing_pairs = {frozenset((s, d)) for s, st, rel, d, dt in edges}
    doi2se: dict[str, str] = {}
    ay_all: dict[tuple, set] = {}
    for pg in raw.get("SourceEvidence", []):
        sid = pg["id"].replace("-", "")
        ck = None
        for pr in pg["properties"].values():
            if pr["type"] == "title":
                ck = plain(pr)
            if pr["type"] == "url" and pr["url"]:
                doi2se.setdefault(_norm_doi(pr["url"]), sid)
        if ck:
            m = re.match(r'^([A-Za-zÀ-ſ]+?)(\d{4})', ck)
            if m and len(m.group(1)) >= 4:
                ay_all.setdefault((m.group(1).lower(), m.group(2)), set()).add(sid)
    ay2se = {k: next(iter(v)) for k, v in ay_all.items() if len(v) == 1}
    n_cite = 0
    for pg in raw.get("KnowledgeClaim", []):
        cid = pg["id"].replace("-", "")
        blob = " ".join(plain(p) for p in pg["properties"].values()
                        if p["type"] == "rich_text")
        if not blob.strip():
            continue
        hits = {doi2se[_norm_doi(m)] for m in _DOI.findall(blob)
                if _norm_doi(m) in doi2se}
        for (surn, yr), sid in ay2se.items():
            if re.search(r'\b' + re.escape(surn) + r'\b.{0,25}\b' + yr + r'\b',
                         blob, re.I):
                hits.add(sid)
        for sid in hits:
            if sid == cid or frozenset((cid, sid)) in existing_pairs:
                continue
            existing_pairs.add(frozenset((cid, sid)))
            edges.append((cid, "knowledge_claim", "cites", sid, "source_evidence"))
            n_cite += 1
    print(f"  derived {n_cite} citation (cites) edges from claim prose", flush=True)

    # ---- edges + meta --------------------------------------------------------
    con.execute(
        "CREATE OR REPLACE TABLE edges (src_id TEXT, src_type TEXT, rel TEXT, "
        "dst_id TEXT, dst_type TEXT)"
    )
    con.executemany("INSERT INTO edges VALUES (?, ?, ?, ?, ?)", edges)
    con.execute("CREATE OR REPLACE TABLE _meta (table_name TEXT, n_rows INTEGER)")
    con.executemany("INSERT INTO _meta VALUES (?, ?)", meta)

    n_nodes = sum(n for _, n in meta)
    print(f"\nkb.duckdb: {n_nodes} node rows across {len(meta)} tables, "
          f"{len(edges)} relation edges -> {DB_PATH}")
    con.close()


def stats():
    con = duckdb.connect(str(DB_PATH), read_only=True)
    print("=== row counts ===")
    for t, n in con.execute("SELECT table_name, n_rows FROM _meta ORDER BY 1").fetchall():
        print(f"  {t:18s} {n}")
    print(f"  {'edges':18s} {con.execute('SELECT count(*) FROM edges').fetchone()[0]}")
    print("\n=== sample 2-hop join: SourceEvidence -> KnowledgeClaim ===")
    rows = con.execute("""
        SELECT se.citation_key, e.rel, kc.title
        FROM source_evidence se
        JOIN edges e ON e.src_id = se.id AND e.dst_type = 'knowledge_claim'
        JOIN knowledge_claim kc ON kc.id = e.dst_id
        LIMIT 5
    """).fetchall()
    for ck, rel, title in rows:
        print(f"  {ck}  --{rel}-->  {(title or '')[:60]}")
    con.close()


if __name__ == "__main__":
    if "--stats" in sys.argv:
        stats()
    else:
        build()
        stats()
