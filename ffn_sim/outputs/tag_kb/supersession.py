#!/usr/bin/env python3
"""Supersession / authoritative-chain layer for the ffn_cellsim TAG+RAG KB.

WHY THIS EXISTS
---------------
Some claims in this project (KU-3.5 cortical tension is the canonical example)
have had their conclusion overturned multiple times. Older, since-refuted
hypotheses live on in docs / dev-logs / commit messages and — because they are
text-heavy and keyword-rich — can out-rank the *latest authoritative* conclusion
in a naive RAG/TAG retrieval. This module makes "X supersedes Y" / "Y is
superseded_by X" / "authoritative_as_of <date>" **first-class, queryable
metadata** so the newest authoritative record always wins.

SOURCE-OF-TRUTH DISCIPLINE (hard)
---------------------------------
Supersession is recorded ONLY where a human EXPLICITLY declared it — never
inferred locally from dates or text similarity:

  1. Docs front-matter   — a `kb_record:` YAML block at the top of a
     `ffn_sim/docs/**/*.md` doc. This is the live source today (e.g.
     `CORTICAL_TENSION_RECORD_2026-06-04.md`, which the gate-structure doc
     explicitly asked TAG/RAG to tag as authoritative over 5/31-6/3).
  2. DecisionLedger      — an explicit `@supersession ...` marker in a Notion
     DecisionLedger row's text (read from the already-materialized
     `decision_ledger` table; Notion stays the source of truth). No-op until
     a PI adds one.

A record referenced only as a target of someone's `supersedes:`/`superseded_by:`
(e.g. a fileless commit-message claim) becomes a *stub* record — its status is
taken from the explicit edge, not guessed.

Human-readable companion: `SUPERSEDED_FINDINGS_2026-06-04.md` is the PI-reviewable
candidate list ("WHY superseded"); this module is the engine that enacts those
calls in the read layers once they are recorded as `kb_record:` front-matter.

OUTPUTS (additive — never rebuilds the 8 node tables or paper_chunks)
---------------------------------------------------------------------
Into `kb.duckdb` (built first by notion_to_duckdb.py + references_ingest.py):
  • kb_record           — one row per record (topic, status, authoritative_as_of,
      conclusion, claim/gate anchors, aliases, source).
  • supersession        — directional edges (newer_id supersedes older_id) with
      the declaring source, de-duplicated.
  • authoritative_record (VIEW) — the single current authoritative record per
      topic (status='authoritative' preferred, then is_primary, then newest date).

RUN (after notion_to_duckdb.py + references_ingest.py):
  conda activate ffn_sim
  python supersession.py            # build the 3 objects into kb.duckdb
  python supersession.py --stats    # print every chain + the authoritative view
  python supersession.py --check    # KU-3.5 authoritative-status sanity (PASS/FAIL)
"""
from __future__ import annotations

import json
import pathlib
import re
import sys

import duckdb
import yaml

HERE = pathlib.Path(__file__).parent
DB_PATH = HERE / "kb.duckdb"
DOCS_DIR = HERE.parent.parent / "docs"          # .../ffn_sim/docs

VALID_STATUS = {"authoritative", "superseded", "open"}

# explicit DecisionLedger marker, e.g.
#   @supersession topic=KU-3.5-cortical-tension status=authoritative
#   as_of=2026-06-04 supersedes=foo,bar superseded_by=baz
LEDGER_MARKER = re.compile(r"@supersession\b(?P<body>[^\n]*)", re.I)
_KV = re.compile(r"(\w+)\s*=\s*([^\s]+)")


# --------------------------------------------------------------------------- #
# front-matter parsing
# --------------------------------------------------------------------------- #
def parse_front_matter(text: str) -> dict | None:
    """Return the leading YAML front-matter as a dict, or None if absent.

    A doc opts in by leading with a `---` ... `---` block; we only act on docs
    whose block carries a `kb_record:` key.
    """
    if not text.startswith("---"):
        return None
    # split on the first two `---` fence lines
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        return None
    end = None
    for i in range(1, len(lines)):
        if lines[i].strip() == "---":
            end = i
            break
    if end is None:
        return None
    try:
        data = yaml.safe_load("\n".join(lines[1:end]))
    except yaml.YAMLError:
        return None
    return data if isinstance(data, dict) else None


def _as_list(v) -> list[str]:
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return [str(x).strip() for x in v if str(x).strip()]
    return [s.strip() for s in str(v).split(",") if s.strip()]


def _norm_record(rec: dict, *, record_id: str, source_path: str | None,
                 source_kind: str) -> dict:
    """Coerce a raw kb_record mapping into the canonical record shape."""
    status = str(rec.get("status", "open")).strip().lower()
    if status not in VALID_STATUS:
        status = "open"
    return {
        "topic": str(rec.get("topic", "")).strip(),
        "record_id": str(rec.get("record_id") or record_id).strip(),
        "claim": (str(rec["claim"]).strip() if rec.get("claim") else None),
        "gate": (str(rec["gate"]).strip() if rec.get("gate") else None),
        "status": status,
        "is_primary": bool(rec.get("primary", False)),
        "authoritative_as_of": (str(rec["authoritative_as_of"]).strip()
                                if rec.get("authoritative_as_of") else None),
        "conclusion": (" ".join(str(rec["conclusion"]).split())
                       if rec.get("conclusion") else None),
        "aliases": _as_list(rec.get("aliases")),
        "supersedes": _as_list(rec.get("supersedes")),
        "superseded_by": _as_list(rec.get("superseded_by")),
        "source_path": source_path,
        "source_kind": source_kind,
    }


# --------------------------------------------------------------------------- #
# collectors (explicit sources only)
# --------------------------------------------------------------------------- #
def collect_doc_records(docs_dir: pathlib.Path = DOCS_DIR) -> list[dict]:
    """Scan docs/**/*.md for `kb_record:` front-matter. Offline; no Notion/duckdb."""
    out: list[dict] = []
    if not docs_dir.exists():
        return out
    for md in sorted(docs_dir.rglob("*.md")):
        try:
            fm = parse_front_matter(md.read_text(encoding="utf-8", errors="ignore"))
        except Exception:
            continue
        if not fm or "kb_record" not in fm:
            continue
        block = fm["kb_record"]
        if not isinstance(block, dict):
            continue
        rel = md.relative_to(docs_dir.parent.parent).as_posix()  # ffn_sim/docs/...
        out.append(_norm_record(block, record_id=md.stem,
                                source_path=rel, source_kind="doc"))
    return out


def collect_ledger_records(con) -> list[dict]:
    """Read explicit `@supersession` markers from the decision_ledger table.

    Notion is the source of truth; this only reads the already-materialized
    table. Returns [] when the table is absent or carries no marker (today).
    """
    have = con.execute(
        "SELECT count(*) FROM information_schema.tables "
        "WHERE table_name='decision_ledger'").fetchone()[0]
    if not have:
        return []
    cols = {r[0] for r in con.execute(
        "SELECT column_name FROM information_schema.columns "
        "WHERE table_name='decision_ledger'").fetchall()}
    text_cols = [c for c in ("decision_text", "rationale") if c in cols]
    id_col = "dec_id" if "dec_id" in cols else "title"
    sel = ", ".join([f'"{id_col}"', '"title"', *(f'"{c}"' for c in text_cols)])
    out: list[dict] = []
    for row in con.execute(f"SELECT {sel} FROM decision_ledger").fetchall():
        rid, title = row[0], row[1]
        blob = "\n".join(str(x) for x in row[2:] if x)
        m = LEDGER_MARKER.search(blob)
        if not m:
            continue
        kv = dict(_KV.findall(m.group("body")))
        if not kv.get("topic"):
            continue
        out.append(_norm_record(
            {"topic": kv["topic"],
             "status": kv.get("status", "authoritative"),
             "authoritative_as_of": kv.get("as_of"),
             "supersedes": kv.get("supersedes"),
             "superseded_by": kv.get("superseded_by"),
             "conclusion": (title or "").strip()},
            record_id=str(rid or title or "DEC"),
            source_path=None, source_kind="decision_ledger"))
    return out


def _stub_records(records: list[dict]) -> list[dict]:
    """Synthesize records for ids referenced only as supersession targets.

    Status comes from the explicit edge direction (an id in someone's
    `supersedes:` is superseded; in `superseded_by:` is the superseder) — never
    guessed. Topic inherits from the referencing record.
    """
    known = {r["record_id"] for r in records}
    stubs: dict[str, dict] = {}
    for r in records:
        for older in r["supersedes"]:
            if older not in known and older not in stubs:
                stubs[older] = {"status": "superseded", "topic": r["topic"]}
        for newer in r["superseded_by"]:
            if newer not in known and newer not in stubs:
                stubs[newer] = {"status": "authoritative", "topic": r["topic"]}
    return [_norm_record(
        {"topic": s["topic"], "status": s["status"],
         "conclusion": "(referenced supersession target; no standalone record)"},
        record_id=rid, source_path=None, source_kind="stub")
        for rid, s in stubs.items()]


def collect_records(con=None) -> tuple[list[dict], list[dict]]:
    """All explicit records + de-duplicated directional supersession edges.

    Edge = (topic, newer_id supersedes older_id, declaring source). Built from
    both `supersedes:` and `superseded_by:` so either declaration direction works.
    """
    records = collect_doc_records()
    if con is not None:
        records += collect_ledger_records(con)
    records += _stub_records(records)

    edges: dict[tuple, dict] = {}

    def add_edge(newer: str, older: str, topic: str, src: str, kind: str):
        if newer and older and newer != older:
            edges.setdefault((newer, older), {
                "topic": topic, "newer_id": newer, "older_id": older,
                "source_record": src, "source_kind": kind})

    for r in records:
        for older in r["supersedes"]:
            add_edge(r["record_id"], older, r["topic"], r["record_id"], r["source_kind"])
        for newer in r["superseded_by"]:
            add_edge(newer, r["record_id"], r["topic"], r["record_id"], r["source_kind"])
    return records, list(edges.values())


# --------------------------------------------------------------------------- #
# authoritative resolution (the same ordering the SQL view uses)
# --------------------------------------------------------------------------- #
def _auth_sort_key(r: dict):
    # higher tuple sorts first under reverse=True
    return (
        r["status"] == "authoritative",
        r["is_primary"],
        r["authoritative_as_of"] or "",
        r["record_id"],
    )


def resolve_authoritative(records: list[dict], topic: str) -> dict | None:
    cands = [r for r in records if r["topic"] == topic]
    if not cands:
        return None
    return sorted(cands, key=_auth_sort_key, reverse=True)[0]


def match_topic(records: list[dict], question: str) -> str | None:
    """Best topic for a free-text question via topic/claim/alias keyword hits.

    Used by tag_query.py to decide when to inject the authoritative chain. Pure
    substring/keyword matching — no inference about *which* record wins (that is
    governed entirely by the explicit metadata)."""
    q = question.lower()
    best, best_score = None, 0
    by_topic: dict[str, list[dict]] = {}
    for r in records:
        by_topic.setdefault(r["topic"], []).append(r)
    for topic, rs in by_topic.items():
        keys = {topic.lower(), topic.lower().replace("-", " ")}
        for r in rs:
            if r["claim"]:
                keys.add(r["claim"].lower())
            for a in r["aliases"]:
                keys.add(a.lower())
                keys.add(a.lower().replace("-", ""))
        score = sum(1 for k in keys if k and k in q)
        if score > best_score:
            best, best_score = topic, score
    return best if best_score else None


# --------------------------------------------------------------------------- #
# build into kb.duckdb (CREATE OR REPLACE only our own 3 objects)
# --------------------------------------------------------------------------- #
def build(db_path: pathlib.Path = DB_PATH) -> None:
    if not db_path.exists():
        sys.exit(f"{db_path} not found — run notion_to_duckdb.py first.")
    con = duckdb.connect(str(db_path))
    try:
        records, edges = collect_records(con)

        con.execute("""
            CREATE OR REPLACE TABLE kb_record (
                topic TEXT, record_id TEXT, claim TEXT, gate TEXT, status TEXT,
                is_primary TEXT, authoritative_as_of TEXT, conclusion TEXT,
                aliases TEXT, source_path TEXT, source_kind TEXT)
        """)
        con.executemany(
            "INSERT INTO kb_record VALUES (?,?,?,?,?,?,?,?,?,?,?)",
            [[r["topic"], r["record_id"], r["claim"], r["gate"], r["status"],
              "true" if r["is_primary"] else "false", r["authoritative_as_of"],
              r["conclusion"], json.dumps(r["aliases"]), r["source_path"],
              r["source_kind"]] for r in records])

        con.execute("""
            CREATE OR REPLACE TABLE supersession (
                topic TEXT, newer_id TEXT, older_id TEXT,
                source_record TEXT, source_kind TEXT)
        """)
        con.executemany(
            "INSERT INTO supersession VALUES (?,?,?,?,?)",
            [[e["topic"], e["newer_id"], e["older_id"],
              e["source_record"], e["source_kind"]] for e in edges])

        # one current authoritative record per topic (deterministic ordering)
        con.execute("""
            CREATE OR REPLACE VIEW authoritative_record AS
            WITH ranked AS (
                SELECT *, row_number() OVER (
                    PARTITION BY topic
                    ORDER BY (status='authoritative') DESC,
                             (is_primary='true') DESC,
                             authoritative_as_of DESC,
                             record_id DESC) AS rn
                FROM kb_record)
            SELECT topic, record_id, claim, gate, status, authoritative_as_of,
                   conclusion, source_path
            FROM ranked WHERE rn = 1
        """)
    finally:
        con.close()

    n_auth = sum(1 for r in records if r["status"] == "authoritative")
    n_sup = sum(1 for r in records if r["status"] == "superseded")
    topics = sorted({r["topic"] for r in records})
    print(f"supersession: {len(records)} records "
          f"({n_auth} authoritative, {n_sup} superseded) across "
          f"{len(topics)} topic(s), {len(edges)} supersession edge(s) -> {db_path}")
    for t in topics:
        auth = resolve_authoritative(records, t)
        print(f"  [{t}] authoritative = {auth['record_id'] if auth else '(none)'}"
              f" (as_of {auth['authoritative_as_of'] if auth else '?'})")


# --------------------------------------------------------------------------- #
# reporting + sanity
# --------------------------------------------------------------------------- #
def stats() -> None:
    records, edges = collect_records()
    by_topic: dict[str, list[dict]] = {}
    for r in records:
        by_topic.setdefault(r["topic"], []).append(r)
    for topic, rs in sorted(by_topic.items()):
        print(f"\n=== topic: {topic} ===")
        auth = resolve_authoritative(records, topic)
        print(f"  AUTHORITATIVE -> {auth['record_id']} (as_of {auth['authoritative_as_of']})")
        if auth.get("conclusion"):
            print(f"    {auth['conclusion'][:200]}")
        for r in sorted(rs, key=_auth_sort_key, reverse=True):
            tag = "AUTH" if r is auth else r["status"].upper()
            print(f"  - [{tag:12s}] {r['record_id']}  ({r['source_kind']})")
    print(f"\nsupersession edges ({len(edges)}):")
    for e in edges:
        print(f"  {e['newer_id']}  --supersedes-->  {e['older_id']}"
              f"   [{e['topic']}, src={e['source_record']}]")


# canonical KU-3.5 test case (request item #5)
KU35_TOPIC = "KU-3.5-cortical-tension"
KU35_AUTHORITATIVE = "CORTICAL_TENSION_RECORD_2026-06-04"
# the since-refuted claim that must NEVER read as the current conclusion
KU35_DEAD_MARKERS = ("v0-accel", "v0accel")


def sanity_ku35(records: list[dict] | None = None) -> tuple[bool, list[str]]:
    """KU-3.5 authoritative-status sanity (request item #5).

    PASS  : the resolved authoritative record is the 2026-06-04 force-aggregation
            / active-passive-split record, and its conclusion reflects that.
    FAIL  : the dead `--v0-accel resolved` claim resolves as authoritative.
    """
    if records is None:
        records, _ = collect_records()
    msgs: list[str] = []
    ok = True

    auth = resolve_authoritative(records, KU35_TOPIC)
    if auth is None:
        return False, [f"FAIL: no records for topic {KU35_TOPIC}"]

    if auth["record_id"] == KU35_AUTHORITATIVE:
        msgs.append(f"PASS: authoritative record = {auth['record_id']} "
                    f"(as_of {auth['authoritative_as_of']})")
    else:
        ok = False
        msgs.append(f"FAIL: authoritative record = {auth['record_id']}, "
                    f"expected {KU35_AUTHORITATIVE}")

    concl = (auth.get("conclusion") or "").lower()
    if "aggregation" in concl and ("active" in concl and "passive" in concl):
        msgs.append("PASS: conclusion reflects force-aggregation + active/passive split")
    else:
        ok = False
        msgs.append("FAIL: authoritative conclusion does not reflect the 6/4 "
                    "force-aggregation / active-passive split")

    # the dead --v0-accel "resolved" claim must be present AND not authoritative
    dead = [r for r in records
            if any(m in r["record_id"].lower() for m in KU35_DEAD_MARKERS)]
    if not dead:
        ok = False
        msgs.append("FAIL: the --v0-accel 'resolved' claim is not tracked at all "
                    "(cannot guarantee it stays demoted)")
    for r in dead:
        if r["status"] == "authoritative":
            ok = False
            msgs.append(f"FAIL: dead claim {r['record_id']} is marked authoritative")
        else:
            msgs.append(f"PASS: dead claim {r['record_id']} is '{r['status']}' (demoted)")

    return ok, msgs


def check() -> int:
    ok, msgs = sanity_ku35()
    print("=== KU-3.5 authoritative-status sanity ===")
    for m in msgs:
        print(f"  {m}")
    print(f"\n{'PASS' if ok else 'FAIL'}: KU-3.5 authoritative chain")
    return 0 if ok else 1


if __name__ == "__main__":
    if "--check" in sys.argv:
        sys.exit(check())
    elif "--stats" in sys.argv:
        stats()
    else:
        build()
