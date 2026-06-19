#!/usr/bin/env python3
"""Results-integrity audit — the RESULTS-side twin of verify_sources.py.

`verify_sources.py` audits the LITERATURE side (do cited papers exist / match).
This audits the RESULTS side: does each headline "done / breakthrough / validated"
result-claim actually hold on DISK — artifact present? value matches? metric
sanctioned (top-down silhouette, NOT the PI-forbidden basal-contact / footprint)?
reproducible, or GPU-only with no committed build/CI trace?

WHY THIS EXISTS (2026-06-19). The RAG/TAG Contract-Graph rigorously audits
citations (`source_audit`) but had NO equivalent for results — and `harvest_ops.py`
launders unverified REPORT numbers straight into the graph as if true. The
2026-06-19 grounding pass found the worst drift exactly there: DCM REPORTs
headlined a PI-FORBIDDEN footprint A/A0 (~20x "spreading") while the only valid
top-down metric COMPACTS (A/A0 -> 0.597); "breakthrough" numbers (cleanball 1.87,
r^2=0.96) lived only in PNG titles with no committed data. This closes that hole.

The contract is `results_manifest.yaml`: every headline result-claim is declared
there, and verify_runs checks it against disk. DISK OVERRIDES THE CLAIM — this
can only DOWNGRADE an over-stated claim, never upgrade it. A claim whose disk
verdict is WORSE than its declared level is a DRIFT and FAILS the CI gate.

Verdicts (suspicion-ranked, worst first — mirrors verify_sources.SUSPICION):
  RETRACT          disk contradicts the claim: forbidden metric, or value drift  <- worst
  NEEDS_REGEN      claimed artifact ABSENT on disk (only prose/PNG) — regenerate
  GPU_UNREPRODUCED artifact present but GPU-only, no committed build/CI trace
  VERIFIED         artifact present, metric sanctioned, value matches            <- best

Writes the `run_audit` table into kb.duckdb (when duckdb is importable) +
`run_audit_report.md`.

RUN:
  conda activate ffn_sim
  python verify_runs.py            # full audit -> run_audit table + report
  python verify_runs.py --check    # read-only verdict distribution (refresh.sh)
  python verify_runs.py --gate     # CI gate: exit 1 if any claim drifted
                                   # (disk verdict worse than declared)

The gate needs only PyYAML + the committed artifacts (duckdb optional), so it
runs in a minimal CI env.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import yaml

try:
    import duckdb
except Exception:  # duckdb is optional: the CI gate must run without it
    duckdb = None

HERE = Path(__file__).parent
DB_PATH = HERE / "kb.duckdb"
MANIFEST = HERE / "results_manifest.yaml"
REPORT = HERE / "run_audit_report.md"
REPO_ROOT = HERE.parents[2]  # tag_kb -> outputs -> ffn_sim -> repo root

# suspicion-ranked: LOWER = worse (mirror of verify_sources.SUSPICION)
SUSPICION = {"RETRACT": 0, "NEEDS_REGEN": 1, "GPU_UNREPRODUCED": 2, "VERIFIED": 3}

# A/A0 / spreading MUST be the PI-mandated top-down xy silhouette (rule 2026-06-12).
# Any basal-contact / convex-hull / footprint metric inflates ~10-20x and inverts
# the sign -> automatic RETRACT.
FORBIDDEN_METRICS = {"basal-contact-hull", "basal-contact", "basal-footprint",
                     "footprint", "convex-hull", "basal"}

DECLARED_ALIASES = {"verified": "VERIFIED", "gpu": "GPU_UNREPRODUCED",
                    "gpu-unreproduced": "GPU_UNREPRODUCED",
                    "needs-regen": "NEEDS_REGEN", "regen": "NEEDS_REGEN",
                    "retract": "RETRACT", "retracted": "RETRACT"}

MEANING = {"RETRACT": "disk contradicts claim (forbidden metric / value drift)",
           "NEEDS_REGEN": "claimed artifact ABSENT — regenerate + commit",
           "GPU_UNREPRODUCED": "present but GPU-only, no committed build/CI trace",
           "VERIFIED": "artifact present, metric sanctioned, value matches"}


def _load_number(path: Path, field: str):
    """Best-effort pull of a numeric field from a .json artifact via a dotted
    path (e.g. 'a.b.0'). .pkl/.npz/.png are opaque to a pure-stdlib gate ->
    return None (existence-only check)."""
    if path.suffix != ".json":
        return None
    try:
        cur = json.loads(path.read_text())
    except Exception:
        return None
    for key in str(field).split("."):
        if isinstance(cur, list):
            try:
                cur = cur[int(key)]
            except Exception:
                return None
        elif isinstance(cur, dict) and key in cur:
            cur = cur[key]
        else:
            return None
    return cur if isinstance(cur, (int, float)) and not isinstance(cur, bool) else None


def audit_claim(c: dict) -> tuple[str, str]:
    """Compute the DISK verdict for one manifest claim. Returns (verdict, note).
    Reports only what disk supports — never the optimistic declared label."""
    metric = (c.get("metric") or "").strip().lower()
    art = c.get("artifact")
    path = (REPO_ROOT / art) if art else None

    if metric in FORBIDDEN_METRICS:
        return ("RETRACT",
                f"PI-FORBIDDEN metric '{metric}' — A/A0 must be top-down silhouette")
    if not art:
        return ("NEEDS_REGEN", "no artifact path declared")
    if not path.exists():
        return ("NEEDS_REGEN", f"claimed artifact ABSENT on disk: {art}")

    exp = c.get("expect")
    note = "artifact present"
    if isinstance(exp, dict) and "field" in exp and "value" in exp:
        found = _load_number(path, exp["field"])
        if found is None:
            note = f"present; '{exp['field']}' not machine-checkable ({path.suffix or 'dir'})"
        elif abs(found - float(exp["value"])) > float(exp.get("tol", 0)):
            return ("RETRACT",
                    f"value drift: {exp['field']}={found}, claimed {exp['value']}")
        else:
            note = f"value OK: {exp['field']}={found} ~= {exp['value']}"

    if (c.get("repro") or "").strip().lower() in ("gpu", "gpu-only"):
        return ("GPU_UNREPRODUCED", note + "; GPU-only, no committed build/CI trace")
    return ("VERIFIED", note)


def load_claims() -> list[dict]:
    if not MANIFEST.exists():
        return []
    return (yaml.safe_load(MANIFEST.read_text()) or {}).get("claims", [])


def run_audit() -> list[tuple]:
    """Return rows: (id, artifact, metric, verdict, declared, drift, note)."""
    rows = []
    for c in load_claims():
        verdict, note = audit_claim(c)
        declared = DECLARED_ALIASES.get((c.get("declared") or "").strip().lower(), "")
        # DRIFT = disk strictly worse than the declared (optimistic) level.
        drift = "DRIFT" if (declared and
                            SUSPICION.get(verdict, 9) < SUSPICION.get(declared, 9)) else ""
        rows.append((c.get("id", "?"), c.get("artifact", ""), c.get("metric", ""),
                     verdict, declared, drift, note))
    return rows


def _counts(rows):
    from collections import Counter
    return Counter(r[3] for r in rows)


def _write_db(rows) -> bool:
    if duckdb is None:
        return False
    con = duckdb.connect(str(DB_PATH))
    con.execute("DROP TABLE IF EXISTS run_audit")
    con.execute("CREATE TABLE run_audit (claim_id TEXT, artifact_path TEXT, "
                "metric TEXT, verdict TEXT, declared TEXT, drift TEXT, note TEXT)")
    con.executemany("INSERT INTO run_audit VALUES (?,?,?,?,?,?,?)", rows)
    con.close()
    return True


def _write_report(rows):
    counts = _counts(rows)
    ordered = sorted(rows, key=lambda r: (SUSPICION.get(r[3], 9), r[0]))
    drifts = [r for r in ordered if r[5]]
    lines = ["# Results-integrity audit (disk-grounded)\n",
             f"Audited **{len(rows)}** headline result-claims from `results_manifest.yaml`.\n",
             "## Summary\n", "| verdict | n | meaning |", "|---|---|---|"]
    for v in sorted(counts, key=lambda v: SUSPICION.get(v, 9)):
        lines.append(f"| {v} | {counts[v]} | {MEANING.get(v, '')} |")
    lines += [f"\n**{len(drifts)} DRIFT rows** (disk worse than declared — these FAIL the CI gate).\n",
              "## All claims (suspicion-ranked)\n",
              "| verdict | drift | claim | metric | note |", "|---|---|---|---|---|"]
    for cid, art, metric, verdict, declared, drift, note in ordered:
        lines.append(f"| {verdict} | {drift} | {cid} | {metric} | {note.replace('|', '/')} |")
    REPORT.write_text("\n".join(lines))


def check() -> int:
    """Non-destructive verdict-distribution read for refresh.sh. Always exits 0."""
    rows = run_audit()
    if not rows:
        print("  [runs] results_manifest.yaml empty/absent — no result-claims declared.")
        return 0
    print(f"  [runs] run_audit: {len(rows)} headline result-claims")
    for v in sorted(_counts(rows), key=lambda v: SUSPICION.get(v, 9)):
        flag = "" if v == "VERIFIED" else "  <-- not VERIFIED"
        print(f"           {v:16s} {_counts(rows)[v]}{flag}")
    for cid, art, metric, verdict, declared, drift, note in rows:
        if drift:
            print(f"           ⚠️  DRIFT {cid}: disk={verdict} < declared={declared} — {note}")
    return 0


def gate() -> int:
    """CI gate: exit 1 if any claim's DISK verdict is worse than its declared
    level (overclaim). Honestly-declared NEEDS_REGEN/RETRACT claims (declared ==
    disk) stay GREEN — they are acknowledged-pending, not drift."""
    rows = run_audit()
    _write_report(rows)
    drifts = [r for r in rows if r[5]]
    if drifts:
        print(f"[runs] GATE FAILED — {len(drifts)} result-claim(s) drifted (disk worse than declared):")
        for cid, art, metric, verdict, declared, drift, note in drifts:
            print(f"   ✗ {cid}: declared {declared}, disk {verdict} — {note}")
        print("Fix: regenerate+commit the artifact, correct the metric to top-down "
              "silhouette, or honestly down-declare the claim in results_manifest.yaml.")
        return 1
    print(f"[runs] gate OK — {len(rows)} claims, none drifted (disk supports every declaration).")
    return 0


def main():
    if "--check" in sys.argv:
        sys.exit(check())
    if "--gate" in sys.argv:
        sys.exit(gate())
    rows = run_audit()
    wrote = _write_db(rows)
    _write_report(rows)
    print("=== run_audit verdicts ===")
    for v in sorted(_counts(rows), key=lambda v: SUSPICION.get(v, 9)):
        print(f"  {v:16s} {_counts(rows)[v]}")
    print(f"\ndrift (disk < declared): {sum(1 for r in rows if r[5])}  ->  report: {REPORT}")
    print(f"duckdb run_audit table: {'written' if wrote else 'SKIPPED (duckdb not importable)'}")


if __name__ == "__main__":
    main()
