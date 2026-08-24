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
# non-result artifacts that the globs over-capture — intermediate state, not a run.
NON_RESULT = (".checkpoint", "_ckpt", ".ckpt", "_manifest")


def iter_artifacts(scope: set[str] | None):
    """Production JSON/MD across units (+ layer2 + h1 root json)."""
    seen = set()
    for path in sorted(OUTPUTS.glob("**/production/*.json")) + \
            sorted(OUTPUTS.glob("**/production/*.md")) + \
            sorted(OUTPUTS.glob("layer2/**/*.json")) + \
            sorted(OUTPUTS.glob("h1/*_production.json")):
        if path in seen or "tag_kb" in path.parts or "obsidian_rag_full" in path.parts:
            continue
        # skip intermediate checkpoints/manifests — they are not run RESULTS
        if any(k in path.name for k in NON_RESULT):
            continue
        seen.add(path)
        unit = path.relative_to(OUTPUTS).parts[0]
        if scope and unit not in scope:
            continue
        yield path, unit


# Unit -> representative ValidationGate, for runs that exercise a unit's model but
# carry no explicit gate reference in their text/filename (faithful "belongs-to-unit"
# link, not a fabricated specific-test claim). Multi-gate units use the headline gate.
UNIT_GATE = {
    "h1": "VG-U1-G0-linear",            # H.1 = ECM unit (KU-1.x)
    "h2": "VG-U1-wlc-hamiltonian",      # H.2 = polymer/WLC slab validation
    "h3": "VG-H3-KU35-cortex-tension",  # H.3 = cortex tension
    "h4": "VG-U2-biphasic-traction",    # H.4 = FA motor-clutch
    "h5": "VG-H3-KU35-cortex-tension",  # H.5 = cortex tension (myosin var-N)
    "layer2": "VG-U5-radial-expansion",  # collective spheroid A/A0 law
}
# run_id substrings that mark pure infra/perf/feasibility — NO science gate (honest leaf)
INFRA_KEYS = ("profile", "hotloop", "gpu-pilot", "gpu-probe", "native38k",
              "degcap", "fanin", "compartment-gpu", "native-fullcell-go",
              "force-stack")


def _rep_gate(unit: str, vocab: dict):
    """Resolve the unit's representative gate page-id, or None."""
    vgid = UNIT_GATE.get(unit)
    if not vgid:
        return None
    hit = next((g for g in vocab["gates"] if g["vg_lower"] == vgid.lower()), None)
    return hit["page_id"] if hit else None


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
        rid = run_id_for(path)
        gate_pids = [g["page_id"] for g in gates]
        gate_ids = [g["vg_id"] for g in gates]
        # fallback: science run with no explicit gate -> unit's representative gate.
        # infra/perf/feasibility runs test no science gate -> left as honest leaves.
        if not gate_pids and not any(k in rid for k in INFRA_KEYS):
            rep = _rep_gate(unit, vocab)
            if rep:
                gate_pids = [rep]
                gate_ids = [f"{UNIT_GATE[unit]} (unit-fallback)"]
        rows.append(
            {
                "run_id": rid,
                "unit": unit,
                "title": path.stem.replace("_", " "),
                "outcome": classify_outcome(path, raw),
                "artifact_path": str(path.relative_to(REPO_ROOT)),
                "commit": commit,
                "date": date,
                "snapshot": snap,
                "gates": gate_ids,
                "gate_pids": gate_pids,
                "unmatched": unmatched,
            }
        )
    return rows


def harvest_reports(vocab: dict, scope):
    rows = []
    for path in sorted(OUTPUTS.glob("**/REPORT.md")):
        # The glob is recursive, so keep the FULL directory chain: a nested
        # closeout (outputs/ac/gate_b_sf_motor/REPORT.md) must not collapse onto
        # its parent unit's id — every ac/* report used to land on RUN-ac-closeout
        # and only the last one survived. `unit` stays parts[0] because gate/
        # contract fallback is defined per unit, not per sub-directory.
        rel_dirs = path.relative_to(OUTPUTS).parts[:-1]
        unit = rel_dirs[0]
        if scope and unit not in scope:
            continue
        raw = path.read_text(errors="replace")
        commit, date = git_last_commit(path)
        gates, unmatched = match_gates(raw, vocab)
        contracts = match_contracts(raw, vocab)
        gate_ids = [g["vg_id"] for g in gates]
        gate_pids = [g["page_id"] for g in gates]
        # unit-fallback (same as harvest_runs): a unit closeout belongs to its
        # unit's gate. Foundation units (h_0_2, h1_baoab_freeze) aren't mapped -> leaf.
        if not gate_pids:
            rep = _rep_gate(unit, vocab)
            if rep:
                gate_pids = [rep]
                gate_ids = [f"{UNIT_GATE[unit]} (unit-fallback)"]
        rows.append(
            {
                # depth-1 keeps its historical id (RUN-<unit>-closeout) so the
                # idempotent RUN-ID upsert still matches existing Notion rows.
                "run_id": "RUN-" + "-".join(rel_dirs) + "-closeout",
                "unit": unit,
                "title": f"{'/'.join(rel_dirs)} closeout (REPORT.md)",
                "artifact_path": str(path.relative_to(REPO_ROOT)),
                "commit": commit,
                "date": date,
                "snapshot": md_snapshot(raw),
                "outcome": "",
                "gates": gate_ids,
                "gate_pids": gate_pids,
                "contracts": [c["mc_id"] for c in contracts],
                "unmatched": unmatched,
            }
        )
    return rows


# --------------------------------------------------------------------------- #
# CodeMapping harvest — mechanism modules -> code_mapping (P2)
# --------------------------------------------------------------------------- #
FFN = REPO_ROOT / "aleph"

#: Package roots scanned for mechanism modules.
#:
#: REPOINTED 2026-07-28 (PI: "제대로 수리").  The previous list was
#: ``["cortex", "cell", "ecm", "bridge", "junction", "integrator", "native", "common"]`` — seven of
#: those eight directories were retired on 2026-06-29 and only ``common/`` still existed, so the
#: harvester had been scanning a tree that was mostly gone for a month while reporting healthy drift.
#:
#: That is the root cause of the 84%-dead ``CodeMapping`` table, and it is why re-running ``--apply``
#: would NOT have repaired it: the 42 dead rows are this scanner's own past output, and re-running a
#: scanner pointed at the old layout re-affirms the old layout.  The audit's ordered proposal called
#: step 1 mechanical and sufficient; it was neither.
#:
#: ``iter_modules`` rglobs each entry, so ``ac`` covers ``ac/engine``, ``ac/cell``, ``ac/motor`` and
#: the rest.  ``dcm`` is deliberately absent: STATE.md (a) has it parked with no commits since
#: 2026-07-14, and harvesting a parked tree into the contract graph would present it as live.
#:
#: REPOINTED AGAIN 2026-08-20 — the SAME defect recurred, which is the point.  ``ac``, ``ff`` and
#: ``common`` were all gone by 2026-08-09 (``common/`` dissolved that day, per STATE.md (a)), so
#: ``iter_modules`` yielded ZERO modules and ``--check`` read "0 modules / 50 mapped -> 50 DEAD":
#: 100% of the contract graph's first hop pointing at absent files.  A repoint is therefore not a
#: one-off repair — a hard-coded directory list re-breaks on every tree reorganisation, silently,
#: because an empty scan looks exactly like a clean scan.  Running ``--apply`` against it would have
#: ARCHIVED all 50 rows and created nothing, deleting the table rather than repairing it.
#:
#: The three roots below are the mechanism stack CLAUDE.md §Stack names — "``aleph/engine/`` binds
#: ``aleph/laws/``; parts live in ``aleph/components/``" — so the list is now derived from the charter
#: rather than from whatever the tree happened to be called.  Deliberately ABSENT, each for a stated
#: reason, not by oversight: ``dcm`` (parked, above); ``validation/oracles`` (acceptance oracles, never
#: imported by a runtime path — mapping them would present an oracle as a mechanism); ``scripts``
#: (drivers and gate scripts, versioned as neither); ``inner``/``outer`` (the two LEARNED layers, not
#: mechanistic modules); ``world``/``observe``/``units``/``virtual_cell``/``infer`` — ⚠ these five are a
#: PI GOVERNANCE QUESTION surfaced, not decided here, exactly as the 43 unclaimed ``scripts/`` are in
#: ``ownership.yaml``.
PKG_DIRS = ["engine", "laws", "components"]


def iter_modules():
    # An empty scan is indistinguishable from a clean scan downstream: `--check` prints
    # "0 un-harvested" and `--apply` archives every machine row. That is how this scanner ran
    # against a retired tree for a month in 2026-07, and again for eleven days in 2026-08. The
    # list is hard-coded, so the only place the mistake is CHEAP to catch is here, before a
    # single row moves — hence a raise, not a warning.
    missing = [d for d in PKG_DIRS if not (FFN / d).is_dir()]
    if len(missing) == len(PKG_DIRS):
        raise SystemExit(
            f"[harvest_ops] FATAL: no PKG_DIRS root exists under {FFN} — {missing}.\n"
            "  The tree was reorganised and this scanner was not repointed. Fix PKG_DIRS;\n"
            "  do NOT run --apply first, it would archive every machine CodeMapping row."
        )
    if missing:
        print(f"[harvest_ops] ⚠️ PKG_DIRS entries absent, scanning without them: {missing}")
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


# package-path -> ModelContract fallback for modules whose docstring carries no KU that maps to a
# contract. Keyed on the longest matching path prefix, so "ac/cell" can differ from "ac/engine".
#
# DELIBERATELY SHORT (2026-07-28).  Only ``ac/cell`` is here, because it is the same package the old
# ``cell`` key pointed at, merely moved — a rename, not a judgement.  Every other current directory
# was left unmapped ON PURPOSE: choosing that, say, ``ac/motor`` implements ``MC-U2-fa-motor-clutch``
# (an ADHESION clutch contract, not an NMII one) would be authoring a contract relation, and
# CLAUDE.md reserves that for the PI — "new gates/contracts are PI-authored, never auto-created".
#
# A wrong edge is worse than a missing one here: a missing edge shows up as an honest leaf, while a
# wrong edge makes the traversal `run -> gate -> contract -> parameter -> source` return a confident
# answer about the wrong contract. The untagged modules are enumerated in the manifest instead, so
# the gap is visible and PI-resolvable rather than silently filled.
#: ⚠ EMPTIED 2026-08-20. The single key was ``ac/cell`` -> ``MC-U3-cell-integration``; ``aleph/ac/cell``
#: no longer exists and no current directory is its unambiguous successor. Choosing one would author a
#: contract relation, which CLAUDE.md reserves for the PI — and by this dict's own rule above, a wrong
#: edge is worse than a missing one. So the fallback is now empty and every module links by its explicit
#: ``Implements:`` line or not at all; re-establishing the prefix edge is a PI decision, surfaced here.
DIR_CONTRACT: dict[str, str] = {}


def _contract_by_mcid(mc_id: str, vocab: dict):
    return next((c for c in vocab["contracts"] if c["mc_lower"] == mc_id.lower()), None)


#: An explicit implements-declaration in a module docstring: ``Implements: KU-3.5, KU-3.1``.
#:
#: WHY AN EXPLICIT MARKER REPLACED "any KU in the docstring head" (2026-07-28).  The old rule linked a
#: module to every ModelContract carrying any KU its first 45 lines mentioned, which cannot tell an
#: implements-relation from a citation.  Measured on the repaired scanner, 2 of the 6 KU-derived edges
#: were plainly wrong: ``common/checkpoint.py`` (checkpoint I/O) linked to a cortical-tension contract
#: because its docstring says the checkpoints are used *during* the KU-3.5 sweep, and
#: ``ac/engine/sf_mechanics.py`` linked to the ECM-network contract because it cites KU-1.1 as the
#: SOURCE of a persistence length.  Both are prose about a KU, not an implementation of one.
#:
#: A missing edge is a visible leaf; a wrong edge makes the traversal answer confidently about the
#: wrong contract.  So mentions now generate a REVIEW SUGGESTION in the manifest and never an edge.
IMPLEMENTS_RE = re.compile(r"^\s*Implements:\s*(.+)$", re.MULTILINE)


def implements_contracts(head: str, vocab: dict) -> tuple[list[dict], list[str]]:
    """Return ``(contracts, mentioned_only)`` for one module docstring head.

    Args:
        head: The module's leading lines.
        vocab: The matcher vocabulary from :func:`load_vocab`.

    Returns:
        ``contracts`` are the ModelContracts the module DECLARES it implements, via an explicit
        ``Implements:`` line.  ``mentioned_only`` are KU tokens that appear elsewhere in the head and
        resolve to some contract but were not declared — surfaced for review, never written as edges.
    """
    declared = "\n".join(IMPLEMENTS_RE.findall(head))
    contracts = match_contracts(declared, vocab) if declared.strip() else []
    declared_pids = {c["page_id"] for c in contracts}
    mentioned = [
        c["mc_id"] for c in match_contracts(head, vocab) if c["page_id"] not in declared_pids
    ]
    return contracts, sorted(set(mentioned))


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
        contracts, mentioned_only = implements_contracts(head, vocab)
        has_ku = bool(IMPLEMENTS_RE.search(head))
        # fallback: no KU-matched contract -> link by longest matching package-path prefix
        if not contracts:
            parts = p.relative_to(FFN).parts
            for depth in (2, 1):
                mc = DIR_CONTRACT.get("/".join(parts[:depth]))
                hit = _contract_by_mcid(mc, vocab) if mc else None
                if hit:
                    contracts = [hit]
                    break
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
                "has_ku": has_ku,
                "mentioned_only": mentioned_only,
            }
        )
    return rows


def stale_code_rows(existing: dict) -> list[dict]:
    """Return graph ``CodeMapping`` rows whose ``Path`` no longer exists on disk.

    THE DIRECTION NOBODY WAS CHECKING.  ``_drift_check`` only ever computed ``disk - graph`` (modules
    not yet harvested), so a row pointing at a deleted file was invisible to every check that rode
    ``refresh.sh``.  That is how 42 of 50 rows could be dead for a month while the drift line read
    "1 un-harvested" and looked healthy.

    Args:
        existing: ``path -> (page_id, is_machine)``, as built by :func:`upsert_code`.

    Returns:
        One dict per stale row with ``path`` / ``page_id`` / ``is_machine``.  Combination paths
        (a curated row naming several files with ``+``) are resolved component-wise and count as
        live if ANY component still exists, because archiving a curated multi-file row on the
        strength of one moved file would destroy hand-written work.
    """
    stale = []
    for path, (page_id, is_machine) in existing.items():
        components = [c.strip() for c in path.split("+")] if "+" in path else [path]
        if any((REPO_ROOT / c).exists() for c in components if c):
            continue
        stale.append({"path": path, "page_id": page_id, "is_machine": is_machine})
    return stale


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
    untagged = [c for c in code if not c["has_ku"] and not c["contracts"]]
    L.append("## CodeMapping candidates — mechanism modules\n")
    L.append(f"{len(code)} modules. {n_cmlinked} link to a ModelContract via docstring KU/MC "
             "tokens or the package-path fallback. Existing curated rows (matched by Path) are left "
             "untouched on --apply; machine rows whose file no longer exists are ARCHIVED.\n")
    L.append("| cm_id | status | implements_contract |")
    L.append("|---|---|---|")
    for c in sorted(code, key=lambda x: x["path"]):
        L.append(f"| `{c['cm_id']}` | {c['status']} "
                 f"| {', '.join(c['contracts']) or '—'} |")
    L.append("")

    if untagged:
        L.append(f"### {len(untagged)} module(s) link to NOTHING — the real traversal gap\n")
        L.append("These carry no `KU-x.y` in their docstring head and match no package-path "
                 "fallback, so they harvest as honest leaves. Repointing the scanner (2026-07-28) "
                 "fixed WHERE it looks; it cannot invent what a module implements. Adding a KU to a "
                 "module's docstring head is what links it, and doing that is an authorship "
                 "decision per module, not a batch operation — a wrong KU manufactures a confident "
                 "edge to the wrong contract.\n")
        for c in sorted(untagged, key=lambda x: x["path"]):
            hint = f"  ← mentions {', '.join(c['mentioned_only'])}" if c["mentioned_only"] else ""
            L.append(f"- `{c['path']}`{hint}")
        L.append("")
        suggest = [c for c in untagged if c["mentioned_only"]]
        if suggest:
            L.append(f"{len(suggest)} of those mention a KU that resolves to a contract. Those are "
                     "REVIEW SUGGESTIONS, not edges — under the pre-2026-07-28 rule they would have "
                     "been written as links, which is how `common/checkpoint.py` came to point at a "
                     "cortical-tension contract for saying its checkpoints are used during that "
                     "sweep. Add `Implements: KU-x.y` to the docstring to make one real.\n")

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

    # prune spurious rows previously harvested from non-result files (checkpoints)
    n_prune = 0
    for key, pid in existing.items():
        if any(k in key for k in NON_RESULT):
            if apply:
                requests.patch(f"{API}/pages/{pid}", headers=headers(tok),
                               json={"archived": True}, timeout=30).raise_for_status()
            n_prune += 1

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
    return n_create, n_update, n_prune


def upsert_code(code, apply: bool):
    """Upsert CodeMapping rows keyed on Path, and archive machine rows whose file is gone.

    Non-destructive on CURATED rows: an existing row whose Notes do NOT carry the
    'machine-harvested ... harvest_ops.py' marker is left untouched (so manually
    curated entries are never clobbered). Machine-harvested rows are PATCHed so
    newly-resolved Implements-Contract links (e.g. after new ModelContracts land)
    propagate. Genuinely new paths are created.

    ARCHIVE PASS (added 2026-07-28).  Until now this function could only add, so a module that moved
    or was deleted left its row pointing into nothing forever — 42 of 50 rows by the time anyone
    counted, which broke the advertised ``run -> gate -> contract -> parameter -> source`` traversal
    on its first hop.  Machine-harvested rows whose path is absent from disk are now archived
    (Notion's reversible ``archived: true``, not a delete, because git is the SSOT for what the code
    is and the graph should be recoverable if this judgement is wrong).

    Curated stale rows are counted and returned but NEVER archived: a human wrote them, and a path
    that a person typed may be describing an intent rather than a file.  They surface in the manifest
    for PI instead.

    Returns:
        ``(n_create, n_update, n_skip, n_archived, curated_stale)`` where ``curated_stale`` is the
        list of hand-written rows that point at absent files.
    """
    import requests
    marker = "harvest_ops.py"
    tok, API, DS, headers, query_all = _notion()
    db_id = DS["CodeMapping"]

    # path -> (page_id, is_machine) for exact rows; curated_substr = combo paths
    existing = {}
    curated_substr = []
    for row in query_all(db_id, tok):
        props = row.get("properties", {})
        p = props.get("Path", {})
        val = "".join(x["plain_text"] for x in p.get("rich_text", [])) if p else ""
        nt = props.get("Notes", {})
        notes = "".join(x["plain_text"] for x in nt.get("rich_text", [])) if nt else ""
        if not val:
            continue
        val = val.strip()
        is_machine = marker in notes
        existing[val] = (row["id"], is_machine)
        if not is_machine:
            curated_substr.append(val)  # curated combo paths (e.g. "a + b + c")

    def _props(c, full: bool):
        d = {"Status": {"select": {"name": c["status"]}}}
        if c.get("contract_pids"):
            d["Implements Contract"] = {"relation": [{"id": p} for p in c["contract_pids"]]}
        if full:
            d["Code Ref"] = {"title": [{"text": {"content": c["title"][:1900]}}]}
            d["CM ID"] = _rt(c["cm_id"])
            d["Path"] = _rt(c["path"])
            d["Notes"] = _rt(f"machine-harvested {datetime.now(tz=timezone.utc).date()} "
                             "by harvest_ops.py (P2); git is SSOT")
        return d

    stale = stale_code_rows(existing)
    curated_stale = [s for s in stale if not s["is_machine"]]
    machine_stale = [s for s in stale if s["is_machine"]]
    n_archived = 0
    for s in machine_stale:
        if not apply:
            n_archived += 1
            continue
        requests.patch(f"{API}/pages/{s['page_id']}", headers=headers(tok),
                       json={"archived": True}, timeout=30).raise_for_status()
        n_archived += 1

    n_create = n_update = n_skip = 0
    for c in code:
        path = c["path"]
        # curated combo row that merely *contains* this path -> never touch
        if any(path != cs and path in cs for cs in curated_substr):
            n_skip += 1
            continue
        hit = existing.get(path)
        if hit and not hit[1]:          # exact curated row -> protect
            n_skip += 1
            continue
        if not apply:
            n_create += 1 if not hit else 0
            n_update += 1 if hit else 0
            continue
        if hit:                          # machine row -> PATCH (refresh contract links)
            requests.patch(f"{API}/pages/{hit[0]}", headers=headers(tok),
                           json={"properties": _props(c, full=False)},
                           timeout=30).raise_for_status()
            n_update += 1
        else:
            requests.post(f"{API}/pages", headers=headers(tok),
                          json={"parent": {"database_id": db_id},
                                "properties": _props(c, full=True)},
                          timeout=30).raise_for_status()
            n_create += 1
    return n_create, n_update, n_skip, n_archived, curated_stale


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
    # BOTH directions. Reporting only `disk - graph` is what let 42 dead rows sit unnoticed for a
    # month behind a drift line that read "1 un-harvested": a row pointing at a deleted file is not
    # un-harvested, it is WRONG, and the check had no way to say so.
    dead_code = [p for p in have_paths
                 if not any((REPO_ROOT / c.strip()).exists() for c in p.split("+") if c.strip())]
    print(f"  [ops-drift] RunResult: {len(run_ids)} disk / {len(have_runs)} in graph "
          f"-> {len(new_runs)} un-harvested")
    print(f"  [ops-drift] CodeMapping: {len(code_paths)} modules / {len(have_paths)} "
          f"mapped -> {len(new_code)} un-harvested, {len(dead_code)} DEAD (path absent from disk)")
    if dead_code:
        pct = 100.0 * len(dead_code) / max(1, len(have_paths))
        print(f"  [ops-drift] ⚠️ {pct:.0f}% of CodeMapping paths do not exist — the "
              "run->gate->contract->parameter->source traversal breaks on its first hop")
    untagged = [c["path"] for c in code if not c["has_ku"] and not c["contracts"]]
    if untagged:
        print(f"  [ops-drift] {len(untagged)} of {len(code)} modules carry no KU and match no "
              "contract -> they harvest as honest leaves, linked to nothing")
    if new_runs or new_code or dead_code:
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
        nc, nu, npr = upsert_runs(runs, reports, apply=True, use_llm=args.llm)
        print(f"[harvest_ops] Notion RunResult upsert: {nc} created, {nu} updated, "
              f"{npr} pruned (checkpoints){' (LLM snapshots)' if args.llm else ''}")
        if code:
            cc, cu, cs, ca, curated_stale = upsert_code(code, apply=True)
            print(f"[harvest_ops] Notion CodeMapping archived {ca} machine row(s) whose file is gone")
            if curated_stale:
                print(f"[harvest_ops] ⚠️ {len(curated_stale)} CURATED row(s) also point at absent "
                      "files and were left untouched (a person wrote them) — PI decides:")
                for s in curated_stale:
                    print(f"    - {s['path']}")
            print(f"[harvest_ops] Notion CodeMapping upsert: {cc} created, "
                  f"{cu} updated, {cs} skipped (curated)")
    else:
        print("[harvest_ops] dry-run (no --apply); manifest only")


if __name__ == "__main__":
    main()
