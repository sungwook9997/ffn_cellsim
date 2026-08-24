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

WHAT THIS GATE COULD NOT DO UNTIL 2026-07-25 — and now does. Everything below
audits the rows the manifest LISTS; nothing could ever be wrong about what the
manifest OMITS. Measured that day: all 42 declared constants pointed at
`aleph/validation/oracles/configs/phase1_unit*.yaml` (retired v1 oracle configs
no runtime reads), while the live runtime tree `aleph/configs/**` carried 71
KU-tagged constants of which the manifest covered ZERO — and this gate printed
"gate OK — 42 constants, none drifted". `kb_coverage.py` now adds the OMISSION
side as a ratchet: the current undeclared set is written down item-by-item in the
checked-in `coverage_baseline.yaml`, its size is printed on every run, and it may
only SHRINK. A new KU-tagged constant that is neither declared nor in the
baseline, a silent value change to a baselined constant, or a stale baseline
entry FAILS this gate. A missing/unimportable kb_coverage.py also FAILS it.

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
                                     # (disk verdict worse than declared) OR the
                                     # coverage ratchet widened (see kb_coverage.py)

The gate needs only PyYAML + the committed config YAMLs + source_audit_report.md
(duckdb optional), so it runs in a minimal CI env.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml
from verify_manifest_common import (
    COVERAGE_IMPORT_ERROR,
    audit_manifest,
    check_coverage,
    duckdb,
    kb_coverage,
    run_coverage_gate,
    verdict_counts,
    write_audit_table,
)
from verify_manifest_common import (
    load_claims as load_manifest_claims,
)

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
    return load_manifest_claims(MANIFEST)


def run_audit() -> list[tuple]:
    """Return rows: (id, config:key, ku, citation_key, verdict, declared, drift, note)."""
    verdicts = _load_verdicts()
    return audit_manifest(
        MANIFEST,
        audit_claim=lambda claim: audit_claim(claim, verdicts),
        declared_aliases=DECLARED_ALIASES,
        suspicion=SUSPICION,
        row_builder=lambda claim, verdict, declared, drift, note: (
            claim.get("id", "?"),
            f"{claim.get('config', '?').split('/')[-1]}:{claim.get('key', '?')}",
            claim.get("ku", ""),
            claim.get("citation_key") or "",
            verdict,
            declared,
            drift,
            note,
        ),
    )


def _counts(rows):
    return verdict_counts(rows, 4)


def _write_db(rows) -> bool:
    return write_audit_table(
        DB_PATH,
        "param_audit",
        (
            "claim_id",
            "location",
            "ku",
            "citation_key",
            "verdict",
            "declared",
            "drift",
            "note",
        ),
        rows,
    )


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


def _coverage_gate() -> int:
    """Run the OMISSION-side ratchet. Missing module => FAIL, never skip."""
    return run_coverage_gate(
        kb_coverage,
        COVERAGE_IMPORT_ERROR,
        kind="params",
        label="params",
        manifest_name=MANIFEST.name,
    )


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
    check_coverage(
        kb_coverage, COVERAGE_IMPORT_ERROR, kind="params", label="params"
    )
    return 0


def gate() -> int:
    """CI gate: exit 1 if (a) any constant's DISK verdict is worse than its declared
    level (over-claim), or (b) the coverage ratchet widened — a KU-tagged runtime
    constant that no manifest row declares and no baseline entry admits, a silent
    value change under a baseline entry, or a stale baseline entry.
    Honestly-declared pending constants (declared == disk: unsourced /
    source-unverified) stay GREEN — acknowledged-pending, not drift."""
    rows = run_audit()
    _write_report(rows)
    drifts = [r for r in rows if r[6]]
    rc = 0
    if drifts:
        print(f"[params] GATE FAILED — {len(drifts)} constant(s) drifted "
              f"(disk worse than declared):")
        for cid, loc, ku, ck, verdict, declared, drift, note in drifts:
            print(f"   ✗ {cid} ({loc}): declared {declared}, disk {verdict} — {note}")
        print("Fix: restore the config value, re-anchor to a verdict=OK source, or "
              "honestly down-declare the constant in params_manifest.yaml.")
        rc = 1
    else:
        print(f"[params] declaration gate OK — {len(rows)} declared constants, none "
              f"drifted (disk + citation audit support every declaration).")
    rc |= _coverage_gate()
    return rc


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
    print("\n=== coverage ratchet (what the manifest OMITS) ===")
    if kb_coverage is None:
        print(f"  UNAVAILABLE: {COVERAGE_IMPORT_ERROR!r}  (this FAILS --gate)")
    else:
        kb_coverage.check("params", "params")


if __name__ == "__main__":
    main()
