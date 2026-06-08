#!/usr/bin/env python3
"""harvest_ops.py — operational state (B) -> Contract-Graph (A) candidate harvester.

Populates the empty *contract* side of the graph (run_result=1, code_mapping=2)
from disk operational reality: `outputs/**/production/*.{json,md}` runs, per-unit
`REPORT.md` closeouts, and (P2) module docstrings.

Design: `docs/v2_audit/KB_OPS_LINKAGE_PLAN_2026-06-08.md`.

P0 (this file): deterministic + regex extraction only. NO Notion writes, NO LLM.
  - Deterministic fields  (path, commit, date, run_id, outcome heuristic, JSON
    metric snapshot) never touch an LLM.
  - Gate / contract links are regex-extracted from artifact text and matched
    ONLY against the existing ValidationGate / ModelContract id vocabulary read
    read-only from kb.duckdb. A reference with no matching row is FLAGGED in the
    manifest (e.g. H.7 Gate-A/Gate-B are not registered as ValidationGate rows),
    never auto-created.
  - Output is a dry-run manifest `OPS_HARVEST_CANDIDATES_<date>.md`. PI reviews;
    P1 adds the LLM snapshot + the Notion upsert behind `--apply`.

Usage:
    python harvest_ops.py                 # full scan -> manifest
    python harvest_ops.py --scope h7,layer2
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import duckdb

HERE = Path(__file__).resolve().parent
# tag_kb -> outputs -> ffn_sim -> repo root
REPO_ROOT = HERE.parent.parent.parent
OUTPUTS = HERE.parent
DB_PATH = HERE / "kb.duckdb"

# --------------------------------------------------------------------------- #
# id-vocabulary load (read-only mirror; matching target = Notion page id)
# --------------------------------------------------------------------------- #
KU_RE = re.compile(r"KU-?(\d+)\.(\d+)")          # KU-3.5 / KU3.5 (KnowledgeClaim id)
VG_RE = re.compile(r"VG-[A-Za-z0-9][\w-]+")      # explicit ValidationGate id
# H.7 Gate-A/B — '_' is a regex word-char so \b fails on gate_a_native; use lookarounds
GATE_AB_RE = re.compile(r"(?<![a-z])gate[-_ ]?([ab])(?![a-z])", re.I)
MC_RE = re.compile(r"MC-[A-Za-z0-9][\w-]+")


KU_COMPOUND_RE = re.compile(r"KU-?(\d+)\.(\d+)((?:/\d+(?:\.\d+)?)+)?")


def _ku_tokens(text: str) -> set[str]:
    """KU references -> {'ku<major><minor>'}, incl. compound 'KU-3.5/3.1' notation.

    A trailing '/3.1' carries its own major.minor; a bare '/1' inherits the major
    of the leading KU (so 'KU-3.5/1' -> {ku35, ku31}).
    """
    out: set[str] = set()
    for major, minor, rest in KU_COMPOUND_RE.findall(text):
        out.add(f"ku{major}{minor}")
        for part in re.findall(r"/(\d+(?:\.\d+)?)", rest or ""):
            if "." in part:
                a, b = part.split(".")
                out.add(f"ku{a}{b}")
            else:
                out.add(f"ku{major}{part}")  # bare continuation inherits major
    return out


def load_vocab(con) -> dict:
    """Build matchers from the materialised graph (ids ARE Notion page ids)."""
    gates = []
    for vg_id, title, pid in con.execute(
        "SELECT vg_id, title, id FROM validation_gate"
    ).fetchall():
        blob = f"{vg_id or ''} {title or ''}"
        gates.append(
            {
                "vg_id": vg_id,
                "title": title,
                "page_id": pid,
                "vg_lower": (vg_id or "").lower(),
                "ku": _ku_tokens(blob),
            }
        )
    contracts = []
    for mc_id, title, pid in con.execute(
        "SELECT mc_id, title, id FROM model_contract"
    ).fetchall():
        blob = f"{mc_id or ''} {title or ''}"
        contracts.append(
            {
                "mc_id": mc_id,
                "title": title,
                "page_id": pid,
                "mc_lower": (mc_id or "").lower(),
                "ku": _ku_tokens(blob),
            }
        )
    return {"gates": gates, "contracts": contracts}


def match_gates(text: str, vocab: dict) -> tuple[list[dict], list[str]]:
    """Return (matched gate rows, unmatched *gate-shaped* tokens) — regex only.

    Only explicit gate references (`VG-...`, `Gate-A/B`) that resolve to no row
    are reported as unmatched — those are genuinely-missing ValidationGate rows.
    Bare `KU-x.y` tokens are KnowledgeClaim ids, NOT gates; they are used only to
    POSITIVELY match a gate whose title carries that KU, and a non-matching KU is
    expected (it is just a knowledge claim) so it is never flagged as a missing gate.
    """
    matched: dict[str, dict] = {}
    unmatched: set[str] = set()

    # explicit VG-ids — a no-row hit is a genuine missing gate
    for vg in VG_RE.findall(text):
        hit = next((g for g in vocab["gates"] if g["vg_lower"] == vg.lower()), None)
        if hit:
            matched[hit["page_id"]] = hit
        else:
            unmatched.add(vg)

    # KU references -> gate whose id/title carries that KU token (positive only)
    for ku in _ku_tokens(text):
        for h in (g for g in vocab["gates"] if ku in g["ku"]):
            matched[h["page_id"]] = h

    # H.7 Gate-A/Gate-B — resolve against VG-H7-gate-{a,b} if registered, else flag
    for letter in GATE_AB_RE.findall(text):
        vgid = f"vg-h7-gate-{letter.lower()}"
        hit = next((g for g in vocab["gates"] if g["vg_lower"] == vgid), None)
        if hit:
            matched[hit["page_id"]] = hit
        else:
            unmatched.add(f"Gate-{letter.upper()}")

    return list(matched.values()), sorted(unmatched)


def match_contracts(text: str, vocab: dict) -> list[dict]:
    matched: dict[str, dict] = {}
    for mc in MC_RE.findall(text):
        hit = next((c for c in vocab["contracts"] if c["mc_lower"] == mc.lower()), None)
        if hit:
            matched[hit["page_id"]] = hit
    for ku in _ku_tokens(text):
        for c in vocab["contracts"]:
            if ku in c["ku"]:
                matched[c["page_id"]] = c
    return list(matched.values())


# --------------------------------------------------------------------------- #
# deterministic field extractors
# --------------------------------------------------------------------------- #
def git_last_commit(path: Path) -> tuple[str, str]:
    """(short_hash, author_date YYYY-MM-DD) of the last commit touching `path`."""
    try:
        out = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "log", "-1",
             "--format=%h\t%ad", "--date=short", "--", str(path)],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        if out:
            h, d = out.split("\t", 1)
            return h, d
    except Exception:
        pass
    # fall back to file mtime
    ts = datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc)
    return "", ts.strftime("%Y-%m-%d")


def run_id_for(path: Path) -> str:
    rel = path.relative_to(OUTPUTS).with_suffix("")
    parts = re.sub(r"[/_]+", "-", str(rel)).strip("-")
    # collapse "h7-production-..." -> "h7-..." (production is structural noise)
    parts = parts.replace("-production-", "-")
    return "RUN-" + parts.lower()


OUTCOME_RULES = [
    ("diagnosis", ("diagnos", "_diag", "triage")),
    ("feasibility", ("feas", "probe", "pilot", "smoke", "scaling")),
    ("refute", ("refute", "_refute")),
    ("confirm", ("confirm", "pass")),
    ("production", ("prod", "production", "tier", "result", "sweep", "go")),
]


def classify_outcome(path: Path, text: str) -> str:
    name = path.name.lower()
    body = text[:2000].lower()
    if "verdict" in body and "refute" in body:
        return "refute"
    if "verdict" in body and "confirm" in body:
        return "confirm"
    for outcome, keys in OUTCOME_RULES:
        if any(k in name for k in keys):
            return outcome
    return "production"


SALIENT_KEYS = re.compile(
    r"(verdict|outcome|gamma|gamma_total|gamma_soft|g_soft|band|speedup|"
    r"steps_per_s|steps/s|nonconv|s_grip|tension|status|result)", re.I)


def json_snapshot(obj, limit: int = 6) -> str:
    """Pull salient scalar metrics from a JSON artifact (deterministic)."""
    flat: list[str] = []

    def walk(o, prefix=""):
        if len(flat) >= limit:
            return
        if isinstance(o, dict):
            for k, v in o.items():
                if isinstance(v, (dict, list)):
                    walk(v, f"{prefix}{k}.")
                elif SALIENT_KEYS.search(str(k)) and isinstance(v, (int, float, str, bool)):
                    sv = v if not isinstance(v, float) else round(v, 6)
                    flat.append(f"{prefix}{k}={sv}")
        elif isinstance(o, list):
            for it in o[:3]:
                walk(it, prefix)

    walk(obj)
    return "; ".join(flat[:limit]) if flat else ""


def md_snapshot(text: str) -> str:
    for line in text.splitlines():
        s = line.strip().lstrip("#").strip()
        if s and not s.startswith(("---", "```", "|")):
            return s[:200]
    return ""


# --------------------------------------------------------------------------- #
# harvesters
# --------------------------------------------------------------------------- #
def iter_artifacts(scope: set[str] | None):
    """Production JSON/MD across units (+ layer2 + h1 root json)."""
    seen = set()
    for path in sorted(OUTPUTS.glob("**/production/*.json")) + \
            sorted(OUTPUTS.glob("**/production/*.md")) + \
            sorted(OUTPUTS.glob("layer2/**/*.json")) + \
            sorted(OUTPUTS.glob("h1/*_production.json")):
        if path in seen or "tag_kb" in path.parts or "obsidian_rag_full" in path.parts:
            continue
        seen.add(path)
        unit = path.relative_to(OUTPUTS).parts[0]
        if scope and unit not in scope:
            continue
        yield path, unit


def harvest_runs(vocab: dict, scope):
    rows = []
    for path, unit in iter_artifacts(scope):
        try:
            raw = path.read_text(errors="replace")
        except Exception:
            continue
        commit, date = git_last_commit(path)
        if path.suffix == ".json":
            try:
                obj = json.loads(raw)
                snap = json_snapshot(obj)
            except Exception:
                snap = ""
        else:
            snap = md_snapshot(raw)
        # scan content AND path — gate refs (gate_a) often live only in the filename
        gates, unmatched = match_gates(raw + "\n" + str(path.relative_to(OUTPUTS)), vocab)
        rows.append(
            {
                "run_id": run_id_for(path),
                "unit": unit,
                "title": path.stem.replace("_", " "),
                "outcome": classify_outcome(path, raw),
                "artifact_path": str(path.relative_to(REPO_ROOT)),
                "commit": commit,
                "date": date,
                "snapshot": snap,
                "gates": [g["vg_id"] for g in gates],
                "gate_pids": [g["page_id"] for g in gates],
                "unmatched": unmatched,
            }
        )
    return rows


def harvest_reports(vocab: dict, scope):
    rows = []
    for path in sorted(OUTPUTS.glob("**/REPORT.md")):
        unit = path.relative_to(OUTPUTS).parts[0]
        if scope and unit not in scope:
            continue
        raw = path.read_text(errors="replace")
        commit, date = git_last_commit(path)
        gates, unmatched = match_gates(raw, vocab)
        contracts = match_contracts(raw, vocab)
        rows.append(
            {
                "run_id": f"RUN-{unit}-closeout",
                "unit": unit,
                "title": f"{unit} closeout (REPORT.md)",
                "artifact_path": str(path.relative_to(REPO_ROOT)),
                "commit": commit,
                "date": date,
                "snapshot": md_snapshot(raw),
                "outcome": "",
                "gates": [g["vg_id"] for g in gates],
                "gate_pids": [g["page_id"] for g in gates],
                "contracts": [c["mc_id"] for c in contracts],
                "unmatched": unmatched,
            }
        )
    return rows


# --------------------------------------------------------------------------- #
# CodeMapping harvest — mechanism modules -> code_mapping (P2)
# --------------------------------------------------------------------------- #
FFN = REPO_ROOT / "ffn_sim"
PKG_DIRS = ["cortex", "cell", "ecm", "bridge", "junction",
            "integrator", "native", "common"]


def iter_modules():
    for d in PKG_DIRS:
        base = FFN / d
        if not base.exists():
            continue
        for p in sorted(base.rglob("*.py")):
            if (p.name == "__init__.py" or "__pycache__" in p.parts
                    or "test" in p.name):
                continue
            yield p


def _tests_blob() -> str:
    """Filenames + import lines of the test suite (status heuristic source)."""
    blob = []
    tdir = FFN / "tests"
    if tdir.exists():
        for t in sorted(tdir.glob("test_*.py")):
            blob.append(t.name)
            try:
                blob.append("\n".join(t.read_text(errors="replace").splitlines()[:40]))
            except Exception:
                pass
    return "\n".join(blob)


def harvest_code(vocab):
    """One CodeMapping candidate per mechanism module (deterministic + regex)."""
    tests = _tests_blob()
    rows = []
    for p in iter_modules():
        try:
            head = "\n".join(p.read_text(errors="replace").splitlines()[:45])
        except Exception:
            continue
        stem = p.stem
        contracts = match_contracts(head, vocab)
        # implemented if the suite references this module (import or filename)
        implemented = (f"import {stem}" in tests or f"/{stem}" in tests
                       or f"test_{stem}" in tests or f" {stem}" in tests)
        rel = str(p.relative_to(REPO_ROOT))
        rows.append(
            {
                "cm_id": "CM-" + str(p.relative_to(FFN)).replace("/", "-")[:-3],
                "path": rel,
                "title": str(p.relative_to(FFN)),
                "status": "implemented" if implemented else "draft",
                "contracts": [c["mc_id"] for c in contracts],
                "contract_pids": [c["page_id"] for c in contracts],
            }
        )
    return rows


# --------------------------------------------------------------------------- #
# manifest
# --------------------------------------------------------------------------- #
def write_manifest(runs, reports, today: str, code=None) -> Path:
    code = code or []
    out = HERE / f"OPS_HARVEST_CANDIDATES_{today}.md"
    n_gate = sum(1 for r in runs if r["gates"])
    all_unmatched: dict[str, int] = {}
    for r in runs + reports:
        for u in r["unmatched"]:
            all_unmatched[u] = all_unmatched.get(u, 0) + 1

    L = []
    L.append(f"# Ops-Harvest candidates (DRY RUN) — {today}\n")
    L.append("Generated by `harvest_ops.py` (P0: deterministic + regex, no LLM, "
             "no Notion writes). Review, then P1 adds LLM snapshot + `--apply`.\n")
    L.append(f"- RunResult candidates (production artifacts): **{len(runs)}** "
             f"(current graph: 1)")
    L.append(f"- RunResult candidates (REPORT.md closeouts): **{len(reports)}**")
    L.append(f"- ...with a matched ValidationGate link: **{n_gate}/{len(runs)}**\n")

    if all_unmatched:
        L.append("## ⚠️ Unmatched gate references (FLAGGED for PI — not auto-created)\n")
        L.append("Referenced in artifacts but no matching `ValidationGate` row "
                 "exists. Most important finding: H.7 Gate-A/Gate-B are unregistered.\n")
        L.append("| reference | # artifacts |")
        L.append("|---|---|")
        for tok, n in sorted(all_unmatched.items(), key=lambda x: -x[1]):
            L.append(f"| `{tok}` | {n} |")
        L.append("")

    L.append("## RunResult candidates — production artifacts\n")
    L.append("| run_id | outcome | commit | date | gate | snapshot |")
    L.append("|---|---|---|---|---|---|")
    for r in sorted(runs, key=lambda x: (x["unit"], x["run_id"])):
        gate = ", ".join(r["gates"]) or "—"
        snap = (r["snapshot"] or "").replace("|", "\\|")[:80]
        L.append(f"| `{r['run_id']}` | {r['outcome']} | {r['commit'] or '—'} "
                 f"| {r['date']} | {gate} | {snap} |")
    L.append("")

    L.append("## RunResult candidates — REPORT.md closeouts\n")
    L.append("| run_id | commit | date | gates | contracts |")
    L.append("|---|---|---|---|---|")
    for r in sorted(reports, key=lambda x: x["unit"]):
        L.append(f"| `{r['run_id']}` | {r['commit'] or '—'} | {r['date']} "
                 f"| {', '.join(r['gates']) or '—'} | {', '.join(r['contracts']) or '—'} |")
    L.append("")

    n_cmlinked = sum(1 for c in code if c["contracts"])
    L.append("## CodeMapping candidates — mechanism modules\n")
    L.append(f"{len(code)} modules (current graph: 2). {n_cmlinked} link to a "
             "ModelContract via docstring KU/MC tokens. Existing curated rows "
             "(matched by Path) are left untouched on --apply.\n")
    L.append("| cm_id | status | implements_contract |")
    L.append("|---|---|---|")
    for c in sorted(code, key=lambda x: x["path"]):
        L.append(f"| `{c['cm_id']}` | {c['status']} "
                 f"| {', '.join(c['contracts']) or '—'} |")
    L.append("")

    L.append("## Full snapshots (audit)\n")
    for r in sorted(runs, key=lambda x: (x["unit"], x["run_id"])):
        L.append(f"### `{r['run_id']}`")
        L.append(f"- artifact: `{r['artifact_path']}`")
        L.append(f"- commit `{r['commit'] or '?'}` · {r['date']} · outcome **{r['outcome']}**")
        if r["gates"]:
            L.append(f"- gate -> {', '.join(r['gates'])}")
        if r["unmatched"]:
            L.append(f"- ⚠️ unmatched: {', '.join(r['unmatched'])}")
        L.append(f"- snapshot: {r['snapshot'] or '(none)'}")
        L.append("")

    out.write_text("\n".join(L))
    return out


# --------------------------------------------------------------------------- #
# P1 — LLM result_snapshot (opt-in) + idempotent Notion upsert (--apply)
# --------------------------------------------------------------------------- #
# my deterministic outcome -> the RunResult.Outcome select vocabulary.
# Only confident verdicts are mapped; ambiguous ones stay None (Outcome omitted)
# so we never fabricate a PASS/FAIL the artifact did not state.
OUTCOME_TO_NOTION = {
    "refute": "FAIL",
    "confirm": "PASS",
    "feasibility": "smoke",
    "smoke": "smoke",
    "diagnosis": None,
    "production": None,
}


def llm_snapshot(text: str, run_id: str):
    """1-2 sentence result snapshot via the shared dual LLM backend (opt-in)."""
    try:
        from tag_query import llm  # same-dir; anthropic SDK or `claude -p` fallback
    except Exception:
        return None
    prompt = (
        "You are summarising a simulation run artifact for a database field. "
        "In ONE or TWO sentences, state the concrete result/verdict and the key "
        "number(s) only. No preamble, no markdown. Artifact "
        f"`{run_id}`:\n\n{text[:6000]}"
    )
    try:
        out = llm(prompt).strip()
        return out[:1800] or None
    except Exception:
        return None


def _notion():
    """(token, api_base, db_ids) reusing the existing read pipeline's auth."""
    import sys as _sys
    obs = OUTPUTS / "obsidian_rag_full"
    _sys.path.insert(0, str(obs))
    from notion_to_obsidian import (  # noqa: E402
        API, DATA_SOURCES, get_token, headers, query_all,
    )
    return get_token(), API, DATA_SOURCES, headers, query_all


def _rt(value: str):
    return {"rich_text": [{"text": {"content": (value or "")[:1900]}}]}


def upsert_runs(runs, reports, apply: bool, use_llm: bool):
    """Idempotent create/update of RunResult rows keyed on RUN ID."""
    import requests
    tok, API, DS, headers, query_all = _notion()
    db_id = DS["RunResult"]

    existing = {}
    for row in query_all(db_id, tok):
        props = row.get("properties", {})
        rid = props.get("RUN ID", {})
        key = "".join(x["plain_text"] for x in rid.get("rich_text", [])) if rid else ""
        if key:
            existing[key] = row["id"]

    n_create = n_update = 0
    allrows = runs + reports
    for r in allrows:
        snap = r.get("snapshot") or ""
        if use_llm:
            try:
                raw = (REPO_ROOT / r["artifact_path"]).read_text(errors="replace")
                snap = llm_snapshot(raw, r["run_id"]) or snap
            except Exception:
                pass
        props = {
            "Run": {"title": [{"text": {"content": r["title"][:1900]}}]},
            "RUN ID": _rt(r["run_id"]),
            "Result Snapshot": _rt(snap),
            "Notes": _rt(f"machine-harvested {datetime.now(tz=timezone.utc).date()} "
                         f"by harvest_ops.py (P1); artifact is SSOT"),
            "Commit": _rt(r.get("commit", "")),
            "Artifact Path": _rt(r["artifact_path"]),
        }
        if r.get("date"):
            props["Date"] = {"date": {"start": r["date"]}}
        outc = OUTCOME_TO_NOTION.get(r.get("outcome", ""))
        if outc:
            props["Outcome"] = {"select": {"name": outc}}
        if r.get("gate_pids"):
            props["Gate"] = {"relation": [{"id": p} for p in r["gate_pids"]]}

        if not apply:
            continue
        if r["run_id"] in existing:
            requests.patch(f"{API}/pages/{existing[r['run_id']]}",
                           headers=headers(tok), json={"properties": props},
                           timeout=30).raise_for_status()
            n_update += 1
        else:
            requests.post(f"{API}/pages", headers=headers(tok),
                          json={"parent": {"database_id": db_id}, "properties": props},
                          timeout=30).raise_for_status()
            n_create += 1
    return n_create, n_update


def upsert_code(code, apply: bool):
    """CREATE CodeMapping rows for modules not yet mapped (keyed on Path).

    Non-destructive: existing rows (matched by Path) are SKIPPED so manually
    curated CodeMapping entries are never clobbered.
    """
    import requests
    tok, API, DS, headers, query_all = _notion()
    db_id = DS["CodeMapping"]

    existing_paths = set()
    for row in query_all(db_id, tok):
        p = row.get("properties", {}).get("Path", {})
        val = "".join(x["plain_text"] for x in p.get("rich_text", [])) if p else ""
        if val:
            existing_paths.add(val.strip())

    n_create = n_skip = 0
    for c in code:
        # skip if this exact path (or a curated combo row containing it) exists
        if any(c["path"] == ep or c["path"] in ep for ep in existing_paths):
            n_skip += 1
            continue
        if not apply:
            n_create += 1
            continue
        props = {
            "Code Ref": {"title": [{"text": {"content": c["title"][:1900]}}]},
            "CM ID": _rt(c["cm_id"]),
            "Path": _rt(c["path"]),
            "Status": {"select": {"name": c["status"]}},
            "Notes": _rt(f"machine-harvested {datetime.now(tz=timezone.utc).date()} "
                         "by harvest_ops.py (P2); git is SSOT"),
        }
        if c.get("contract_pids"):
            props["Implements Contract"] = {
                "relation": [{"id": p} for p in c["contract_pids"]]}
        requests.post(f"{API}/pages", headers=headers(tok),
                      json={"parent": {"database_id": db_id}, "properties": props},
                      timeout=30).raise_for_status()
        n_create += 1
    return n_create, n_skip


def _drift_check() -> int:
    """Read-only: count disk run/code artifacts not yet in the Notion graph.

    Rides refresh.sh so each TAG refresh reports harvest drift without writing.
    Always exits 0 (informational).
    """
    con = duckdb.connect(str(DB_PATH), read_only=True)
    vocab = load_vocab(con)
    runs = harvest_runs(vocab, None)
    reports = harvest_reports(vocab, None)
    code = harvest_code(vocab)
    con.close()
    run_ids = {r["run_id"] for r in runs + reports}
    code_paths = {c["path"] for c in code}
    try:
        tok, API, DS, headers, query_all = _notion()
        have_runs = set()
        for row in query_all(DS["RunResult"], tok):
            rid = row.get("properties", {}).get("RUN ID", {})
            k = "".join(x["plain_text"] for x in rid.get("rich_text", [])) if rid else ""
            if k:
                have_runs.add(k)
        have_paths = set()
        for row in query_all(DS["CodeMapping"], tok):
            p = row.get("properties", {}).get("Path", {})
            v = "".join(x["plain_text"] for x in p.get("rich_text", [])) if p else ""
            if v:
                have_paths.add(v.strip())
    except Exception as e:
        print(f"  [ops-drift] Notion unreachable ({str(e)[:60]}); skipped")
        return 0
    new_runs = run_ids - have_runs
    new_code = {p for p in code_paths
                if not any(p == hp or p in hp for hp in have_paths)}
    print(f"  [ops-drift] RunResult: {len(run_ids)} disk / {len(have_runs)} in graph "
          f"-> {len(new_runs)} un-harvested")
    print(f"  [ops-drift] CodeMapping: {len(code_paths)} modules / {len(have_paths)} "
          f"mapped -> {len(new_code)} un-harvested")
    if new_runs or new_code:
        print("  [ops-drift] ⚠️ run `python harvest_ops.py --apply` to sync "
              "(PI-gated; dry-run first without --apply).")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scope", default="", help="comma unit filter e.g. h7,layer2")
    ap.add_argument("--date", default="", help="manifest date stamp (YYYY-MM-DD)")
    ap.add_argument("--llm", action="store_true", help="LLM-enrich result_snapshot")
    ap.add_argument("--apply", action="store_true",
                    help="WRITE to Notion RunResult (idempotent upsert on RUN ID)")
    ap.add_argument("--check", action="store_true",
                    help="drift report: how many disk artifacts are NOT yet "
                         "harvested into Notion (read-only, no writes)")
    args = ap.parse_args()

    if args.check:
        sys.exit(_drift_check())
    scope = {s.strip() for s in args.scope.split(",") if s.strip()} or None
    today = args.date or datetime.now(tz=timezone.utc).strftime("%Y-%m-%d")

    con = duckdb.connect(str(DB_PATH), read_only=True)
    vocab = load_vocab(con)
    runs = harvest_runs(vocab, scope)
    reports = harvest_reports(vocab, scope)
    code = harvest_code(vocab) if not scope else []  # code harvest is repo-wide
    con.close()
    # closeout report rows reuse the run upsert path; ensure gate_pids present
    for rep in reports:
        rep.setdefault("gate_pids", [])
        rep.setdefault("outcome", "")

    out = write_manifest(runs, reports, today, code)
    n_gate = sum(1 for r in runs if r["gates"])
    n_cm = sum(1 for c in code if c["contracts"])
    print(f"[harvest_ops] {len(runs)} run candidates, {len(reports)} report "
          f"closeouts, {n_gate} gate-linked; {len(code)} code modules "
          f"({n_cm} contract-linked) -> {out}")

    if args.apply:
        nc, nu = upsert_runs(runs, reports, apply=True, use_llm=args.llm)
        print(f"[harvest_ops] Notion RunResult upsert: {nc} created, {nu} updated"
              f"{' (LLM snapshots)' if args.llm else ''}")
        if code:
            cc, cs = upsert_code(code, apply=True)
            print(f"[harvest_ops] Notion CodeMapping upsert: {cc} created, "
                  f"{cs} skipped (already mapped)")
    else:
        print("[harvest_ops] dry-run (no --apply); manifest only")


if __name__ == "__main__":
    main()
