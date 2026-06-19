#!/usr/bin/env python3
"""Parameter-provenance integrity gate — the PARAMETERS-side twin of
verify_sources.py (citations) and verify_runs.py (results).

`verify_sources.py` audits the LITERATURE side (do cited papers exist / match).
`verify_runs.py` audits the RESULTS side (does each headline result hold on disk).
This audits the PARAMETERS side: does each headline SIMULATION CONSTANT the model
actually runs on (1) match its declared value in the real config YAML on disk,
(2) carry a KU anchor + a SourceEvidence citation_key, and (3) is that citation's
audit verdict OK — i.e. is the constant traceable to a *verified* source?

WHY THIS EXISTS (2026-06-20). The project's defining HARD rule is "no empirical
magic numbers — every tuning constant traces to a literature source." The KB
rigorously audits the *citation* (source_audit) and now the *results*
(run_audit), but had NO machine link between the two and the CONSTANTS the
simulator runs on. A constant could silently drift off its sourced value, or stay
anchored to a CHECK/DEAD/MISMATCH citation, and nothing failed. This closes that
gap: it binds constant -> KU -> citation_key -> verdict and makes the binding a
CI-checkable fact, not a YAML-comment honor-system. (See the 2026-06-18 KB audit,
finding C-1.) It runs on pure stdlib + PyYAML + the committed
source_audit_report.md, so the verdict chain is checkable WITHOUT kb.duckdb or a
Notion token (duckdb is used, fresher, only when present).

Verdicts (suspicion-ranked, worst first — mirrors verify_sources / verify_runs):
  VALUE_DRIFT       on-disk config value != declared (or constant absent)   <- worst
  UNSOURCED         no KU/citation_key — constant cannot be traced to a paper
  SOURCE_SUSPECT    citation verdict DOI_DEAD / DOI_MISMATCH / NO_DOI_NOMATCH (fabrication-risk)
  SOURCE_UNVERIFIED citation verdict CHECK / NO_DOI_FOUND (real paper, not confirmed-OK)
  VERIFIED          value matches disk + KU present + citation verdict OK    <- best

Writes the `param_audit` table into kb.duckdb (when duckdb is importable) +
`param_audit_report.md`.

RUN:
  conda activate ffn_sim
  python verify_params.py            # full audit -> param_audit table + report
  python verify_params.py --check    # read-only verdict distribution (refresh.sh)
  python verify_params.py --gate     # CI gate: exit 1 if any constant drifted
                                     # (disk verdict worse than declared)

The gate needs only PyYAML + the committed config YAMLs + source_audit_report.md
(duckdb optional), so it runs in a minimal CI env.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml

try:
    import duckdb
except Exception:  # duckdb is optional: the CI gate must run without it
    duckdb = None

HERE = Path(__file__).parent
DB_PATH = HERE / "kb.duckdb"
MANIFEST = HERE / "params_manifest.yaml"
AUDIT_MD = HERE / "source_audit_report.md"   # committed citation verdicts
REPORT = HERE / "param_audit_report.md"
REPO_ROOT = HERE.parents[2]  # tag_kb -> outputs -> ffn_sim -> repo root

# suspicion-ranked: LOWER = worse (mirror of verify_sources / verify_runs)
SUSPICION = {"VALUE_DRIFT": 0, "UNSOURCED": 1, "SOURCE_SUSPECT": 2,
             "SOURCE_UNVERIFIED": 3, "VERIFIED": 4}

MEANING = {
    "VALUE_DRIFT": "on-disk config value != declared (or constant absent on disk)",
    "UNSOURCED": "no KU / citation_key — constant not traceable to a paper",
    "SOURCE_SUSPECT": "citation verdict DOI_DEAD/MISMATCH/NO_DOI_NOMATCH (fabrication-risk)",
    "SOURCE_UNVERIFIED": "citation verdict CHECK/NO_DOI_FOUND (real paper, not confirmed-OK)",
    "VERIFIED": "value matches disk + KU present + citation verdict OK",
}

DECLARED_ALIASES = {
    "verified": "VERIFIED",
    "source-unverified": "SOURCE_UNVERIFIED", "unverified": "SOURCE_UNVERIFIED",
    "source-suspect": "SOURCE_SUSPECT", "suspect": "SOURCE_SUSPECT",
    "unsourced": "UNSOURCED",
    "value-drift": "VALUE_DRIFT", "drift": "VALUE_DRIFT",
}

# source_audit verdict (verify_sources.SUSPICION vocab) -> param verdict class
CITATION_VERDICT_CLASS = {
    "OK": "VERIFIED",
    "CHECK": "SOURCE_UNVERIFIED", "NO_DOI_FOUND": "SOURCE_UNVERIFIED",
    "DOI_DEAD": "SOURCE_SUSPECT", "DOI_MISMATCH": "SOURCE_SUSPECT",
    "NO_DOI_NOMATCH": "SOURCE_SUSPECT",
}
_KNOWN_VERDICTS = set(CITATION_VERDICT_CLASS)


def _load_config_value(config_rel: str, dotted_key: str):
    """Walk a dotted path into a config YAML on disk. Returns a float, or None if
    the file/key is absent or the value is non-numeric."""
    path = REPO_ROOT / config_rel
    if not path.exists():
        return None
    try:
        cur = yaml.safe_load(path.read_text())
    except Exception:
        return None
    for k in dotted_key.split("."):
        if isinstance(cur, dict) and k in cur:
            cur = cur[k]
        else:
            return None
    if isinstance(cur, bool):
        return None
    if isinstance(cur, (int, float)):
        return float(cur)
    try:
        return float(cur)
    except (TypeError, ValueError):
        return None


def _load_verdicts() -> dict:
    """citation_key -> verdict. Prefer the live kb.duckdb source_audit table
    (fresher); fall back to parsing the committed source_audit_report.md so the
    gate is fully CI-checkable without duckdb or a Notion token."""
    # 1. live duckdb table, if present
    if duckdb is not None and DB_PATH.exists():
        try:
            con = duckdb.connect(str(DB_PATH), read_only=True)
            have = con.execute("SELECT count(*) FROM information_schema.tables "
                               "WHERE table_name='source_audit'").fetchone()[0]
            if have:
                rows = con.execute(
                    "SELECT citation_key, verdict FROM source_audit").fetchall()
                con.close()
                if rows:
                    return {ck: v for ck, v in rows if ck}
            con.close()
        except Exception:
            pass
    # 2. committed markdown table (canonical, version-controlled fallback)
    if not AUDIT_MD.exists():
        return {}
    out: dict = {}
    for line in AUDIT_MD.read_text().splitlines():
        if not line.lstrip().startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 2:
            continue
        verdict, ck = cells[0], cells[1]
        # detail rows only: col0 is a known verdict, col1 is a citation_key
        # (NOT a pure count, which is how the summary table looks).
        if verdict in _KNOWN_VERDICTS and ck and not re.fullmatch(r"\d+", ck):
            out.setdefault(ck, verdict)
    return out


def audit_claim(c: dict, verdicts: dict) -> tuple[str, str]:
    """Compute the DISK verdict for one manifest constant. Returns (verdict, note).
    Reports only what disk + the citation audit support — never the optimistic
    declared label."""
    key = c.get("key", "?")
    declared_val = c.get("value")
    found = _load_config_value(c.get("config", ""), key)

    # 1. value must exist on disk and match the declared value
    if found is None:
        return ("VALUE_DRIFT", f"constant ABSENT on disk at {c.get('config','?')}:{key}")
    try:
        target = float(declared_val)
    except (TypeError, ValueError):
        return ("VALUE_DRIFT", f"manifest value not numeric: {declared_val!r}")
    tol = c.get("tol")
    tol = abs(target) * 1e-9 if tol is None else float(tol)
    if abs(found - target) > tol:
        return ("VALUE_DRIFT",
                f"value drift: disk {key}={found:g}, declared {target:g}")

    # 2. must carry a KU anchor + a citation_key to be traceable at all
    ku = (c.get("ku") or "").strip()
    ck = (c.get("citation_key") or "").strip() if c.get("citation_key") else ""
    if not ku and not ck:
        return ("UNSOURCED", f"value OK ({found:g}); no KU tag and no citation_key")
    if not ck:
        return ("UNSOURCED",
                f"value OK ({found:g}); KU {ku} present but KU->source link not committed")

    # 3. the citation's audit verdict
    cv = verdicts.get(ck)
    if cv is None:
        return ("SOURCE_UNVERIFIED",
                f"value OK; {ck} not found in citation audit (run verify_sources.py)")
    cls = CITATION_VERDICT_CLASS.get(cv, "SOURCE_UNVERIFIED")
    note = f"value OK ({found:g}); {ck} citation verdict={cv}"
    return (cls, note)


def load_claims() -> list[dict]:
    if not MANIFEST.exists():
        return []
    return (yaml.safe_load(MANIFEST.read_text()) or {}).get("claims", [])


def run_audit() -> list[tuple]:
    """Return rows: (id, config:key, ku, citation_key, verdict, declared, drift, note)."""
    verdicts = _load_verdicts()
    rows = []
    for c in load_claims():
        verdict, note = audit_claim(c, verdicts)
        declared = DECLARED_ALIASES.get((c.get("declared") or "").strip().lower(), "")
        # DRIFT = disk strictly worse than the declared (optimistic) level.
        drift = "DRIFT" if (declared and
                            SUSPICION.get(verdict, 9) < SUSPICION.get(declared, 9)) else ""
        loc = f"{c.get('config','?').split('/')[-1]}:{c.get('key','?')}"
        rows.append((c.get("id", "?"), loc, c.get("ku", ""),
                     c.get("citation_key") or "", verdict, declared, drift, note))
    return rows


def _counts(rows):
    from collections import Counter
    return Counter(r[4] for r in rows)


def _write_db(rows) -> bool:
    if duckdb is None:
        return False
    con = duckdb.connect(str(DB_PATH))
    con.execute("DROP TABLE IF EXISTS param_audit")
    con.execute("CREATE TABLE param_audit (claim_id TEXT, location TEXT, ku TEXT, "
                "citation_key TEXT, verdict TEXT, declared TEXT, drift TEXT, note TEXT)")
    con.executemany("INSERT INTO param_audit VALUES (?,?,?,?,?,?,?,?)", rows)
    con.close()
    return True


def _write_report(rows):
    counts = _counts(rows)
    ordered = sorted(rows, key=lambda r: (SUSPICION.get(r[4], 9), r[0]))
    drifts = [r for r in ordered if r[6]]
    lines = ["# Parameter-provenance audit (disk + citation-grounded)\n",
             f"Audited **{len(rows)}** headline simulation constants from "
             f"`params_manifest.yaml` (constant -> KU -> citation -> verdict).\n",
             "## Summary\n", "| verdict | n | meaning |", "|---|---|---|"]
    for v in sorted(counts, key=lambda v: SUSPICION.get(v, 9)):
        lines.append(f"| {v} | {counts[v]} | {MEANING.get(v, '')} |")
    lines += [f"\n**{len(drifts)} DRIFT rows** (disk worse than declared — these FAIL the CI gate).\n",
              "## All constants (suspicion-ranked)\n",
              "| verdict | drift | id | location | KU | citation | note |",
              "|---|---|---|---|---|---|---|"]
    for cid, loc, ku, ck, verdict, declared, drift, note in ordered:
        lines.append(f"| {verdict} | {drift} | {cid} | {loc} | {ku} | {ck} | "
                     f"{note.replace('|', '/')} |")
    REPORT.write_text("\n".join(lines))


def check() -> int:
    """Non-destructive verdict-distribution read for refresh.sh. Always exits 0."""
    rows = run_audit()
    if not rows:
        print("  [params] params_manifest.yaml empty/absent — no constants declared.")
        return 0
    print(f"  [params] param_audit: {len(rows)} headline simulation constants")
    counts = _counts(rows)
    for v in sorted(counts, key=lambda v: SUSPICION.get(v, 9)):
        flag = "" if v == "VERIFIED" else "  <-- not VERIFIED"
        print(f"           {v:17s} {counts[v]}{flag}")
    for cid, loc, ku, ck, verdict, declared, drift, note in rows:
        if drift:
            print(f"           ⚠️  DRIFT {cid}: disk={verdict} < declared={declared} — {note}")
    return 0


def gate() -> int:
    """CI gate: exit 1 if any constant's DISK verdict is worse than its declared
    level (over-claim). Honestly-declared pending constants (declared == disk:
    unsourced / source-unverified) stay GREEN — acknowledged-pending, not drift."""
    rows = run_audit()
    _write_report(rows)
    drifts = [r for r in rows if r[6]]
    if drifts:
        print(f"[params] GATE FAILED — {len(drifts)} constant(s) drifted "
              f"(disk worse than declared):")
        for cid, loc, ku, ck, verdict, declared, drift, note in drifts:
            print(f"   ✗ {cid} ({loc}): declared {declared}, disk {verdict} — {note}")
        print("Fix: restore the config value, re-anchor to a verdict=OK source, or "
              "honestly down-declare the constant in params_manifest.yaml.")
        return 1
    print(f"[params] gate OK — {len(rows)} constants, none drifted "
          f"(disk + citation audit support every declaration).")
    return 0


def main():
    if "--check" in sys.argv:
        sys.exit(check())
    if "--gate" in sys.argv:
        sys.exit(gate())
    rows = run_audit()
    wrote = _write_db(rows)
    _write_report(rows)
    print("=== param_audit verdicts ===")
    for v in sorted(_counts(rows), key=lambda v: SUSPICION.get(v, 9)):
        print(f"  {v:17s} {_counts(rows)[v]}")
    print(f"\ndrift (disk < declared): {sum(1 for r in rows if r[6])}  ->  report: {REPORT}")
    print(f"duckdb param_audit table: {'written' if wrote else 'SKIPPED (duckdb not importable)'}")


if __name__ == "__main__":
    main()
