#!/usr/bin/env python3
"""One-time authoring helper: generate params_manifest.yaml for the FULL set of
KU-tagged physical constants in the oracle configs.

Walks the 6 config YAMLs, extracts every numeric constant whose inline comment
carries a `KU-x.y` tag, resolves KU -> SourceEvidence citation_key from the Notion
Contract-Graph (resolved 2026-06-20 via the Notion MCP — map below), joins the
committed source_audit_report.md verdict, and writes the manifest with each
constant's HONEST current status (declared == disk status, so the gate starts at
0 drift and guards future drift: a changed config value or a downgraded citation).

Excluded (not literature-anchored physical constants): discretization
(beads_per_fiber), validation/acceptance thresholds + test fixtures
(acceptance/rounding_test/pair_test), sim-control (n_steps), documentation labels
(gamma_drag_pa_s_um).

Run from tag_kb/_migrations:  python gen_params_manifest.py   (writes ../params_manifest.yaml)
"""
from __future__ import annotations

import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
TKB = HERE.parent
REPO = TKB.parents[2]
CFGDIR = REPO / "ffn_sim/validation/oracles/configs"
AUDIT = TKB / "source_audit_report.md"
OUT = TKB / "params_manifest.yaml"

CONFIGS = ["phase1_unit1.yaml", "phase1_unit2_1.yaml", "phase1_unit2_2.yaml",
           "phase1_unit3.yaml", "phase1_unit4_1.yaml"]

# KU -> primary SourceEvidence citation_key, read from each KnowledgeClaim's
# Citations field in Notion (2026-06-20). None = KnowledgeClaim has no single
# citation (empty Citations / a Phase-1 modelling assumption) or KU not resolved.
KU_SOURCE = {
    "KU-1.1": "Lindstrom2010_Biomaterials",     # collagen persistence length (measurement)
    "KU-1.2": "Jansen2018_BiophysJ",            # single-fibril Young's modulus
    "KU-1.3": "BroederszMacKintosh2014_RMP",    # network connectivity <z>
    "KU-1.5": None,                             # substrate-stiffness range claim, empty Citations
    "KU-1.21": "Discher2005_Science",           # PAA gel as model substrate
    "KU-1.22": None,                            # fiber contour length — KU not cleanly resolved
    "KU-1.26": None,                            # water viscosity / kBT — physical constant (NIST), no paper
    "KU-1.28": None,                            # crosslinker stiffness — KU not cleanly resolved
    "KU-2.6": "DelRio2009_Science",             # talin Bell-Evans (rod stretching)
    "KU-2.18": "Bangasser2013_BiophysJ",        # Odde motor-clutch base parameter set
    "KU-3.5": "Chugh2017_NatCellBiol",          # cortex surface tension band
    "KU-3.17": None,                            # cortex geometry — KU not cleanly resolved
    "KU-4.17": "Buckley2014_Science",           # cadherin catch-bond (KU-4.17 anchors to Buckley 2014)
}

# Per-constant overrides where the config comment names a MORE specific source.
KEY_SOURCE = {
    "junction.dx_star_phase1": "Bell1978_Science",          # "Bell 1978 lower bound"
    "bridge.catch_bond.k_off_slip": "Pereverzev2005_BiophysJ",  # "Pereverzev, KU-2.5/2.18"
    "bridge.catch_bond.F_s": "Pereverzev2005_BiophysJ",
    "bridge.catch_bond.k_off_catch": "Pereverzev2005_BiophysJ",
    "bridge.catch_bond.F_c": "Pereverzev2005_BiophysJ",
}


def excluded(key: str) -> bool:
    return ("beads_per_fiber" in key or "acceptance" in key
            or "rounding_test" in key or "pair_test" in key
            or key.endswith(".n_steps") or key == "cell.gamma_drag_pa_s_um")


# committed-md verdict -> declared status (mirror of verify_params.CITATION_VERDICT_CLASS)
def declared_for(citation_key, verdict):
    if not citation_key:
        return "unsourced"
    if verdict == "OK":
        return "verified"
    if verdict in ("DOI_DEAD", "DOI_MISMATCH", "NO_DOI_NOMATCH"):
        return "source-suspect"
    return "source-unverified"      # CHECK / NO_DOI_FOUND / not-in-audit


def load_verdicts() -> dict:
    out, known = {}, {"OK", "CHECK", "NO_DOI_FOUND", "NO_DOI_NOMATCH",
                      "DOI_MISMATCH", "DOI_DEAD"}
    for line in AUDIT.read_text().splitlines():
        if not line.lstrip().startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) >= 2 and cells[0] in known and cells[1] and not re.fullmatch(r"\d+", cells[1]):
            out.setdefault(cells[1], cells[0])
    return out


KU_RE = re.compile(r"KU-(\d+\.\d+|\d+)")
LINE = re.compile(r"^(\s*)([A-Za-z_]\w*):\s*([-+0-9.][-+0-9.eE]*)\s*#(.*)$")
KEYL = re.compile(r"^(\s*)([A-Za-z_]\w*):\s*(.*)$")


def extract(fname: str):
    stack, rows = [], []
    for ln in (CFGDIR / fname).read_text().splitlines():
        if not ln.strip() or ln.strip().startswith("#"):
            continue
        m = LINE.match(ln)
        if m:
            indent, key, val, cm = len(m.group(1)), m.group(2), m.group(3), m.group(4)
            ku = KU_RE.search(cm)
            while stack and stack[-1][0] >= indent:
                stack.pop()
            path = ".".join(k for _, k in stack) + ("." if stack else "") + key
            if ku:
                rows.append((path, val, "KU-" + ku.group(1), cm.strip()))
        else:
            mk = KEYL.match(ln)
            if mk and mk.group(3) in ("", "{}"):
                indent, key = len(mk.group(1)), mk.group(2)
                while stack and stack[-1][0] >= indent:
                    stack.pop()
                stack.append((indent, key))
    return rows


def clean(comment: str) -> str:
    c = re.sub(r"\s+", " ", comment).strip().strip("#").strip()
    return c.replace('"', "'")[:90]


HEADER = '''\
# params_manifest.yaml — the declarative contract of the simulation CONSTANTS.
#
# The PARAMETERS-side twin of the citation list (verify_sources) and the headline
# results list (verify_runs). Operationalizes the HARD rule "no empirical magic
# numbers": every KU-tagged physical constant the simulator runs on is declared
# here with WHERE it lives on disk, its value, its KU anchor, and the
# SourceEvidence citation_key that anchors it. verify_params.py checks each against
# the config value on disk + the committed citation audit; DISK OVERRIDES THE
# CLAIM (downgrade-only); a claim whose disk verdict is WORSE than declared is a
# DRIFT and FAILS the CI gate. Honestly-declared pending constants stay GREEN.
#
# GENERATED by _migrations/gen_params_manifest.py from the oracle configs +
# KU->source links resolved from the Notion Contract-Graph (2026-06-20). To widen
# coverage or re-anchor a KU, edit the map in that script and re-run, or hand-edit
# a single claim here. Excluded: discretization (beads_per_fiber), validation/
# acceptance thresholds + test fixtures, sim-control, documentation labels.
#
# Fields: id | desc | config | key (dotted) | value | unit | ku | citation_key |
#         declared (verified|source-unverified|source-suspect|unsourced)

claims:
'''


def main():
    verdicts = load_verdicts()
    blocks, tally = [], {}
    for fname in CONFIGS:
        for path, val, ku, comment in extract(fname):
            if excluded(path):
                continue
            ck = KEY_SOURCE.get(path) or KU_SOURCE.get(ku)
            verdict = verdicts.get(ck) if ck else None
            decl = declared_for(ck, verdict)
            tally[decl] = tally.get(decl, 0) + 1
            cid = path.replace(".", "-").replace("_", "-")
            ckline = f'    citation_key: {ck}\n' if ck else '    citation_key: null\n'
            blocks.append(
                f'  - id: {cid}\n'
                f'    desc: "{clean(comment)}"\n'
                f'    config: ffn_sim/validation/oracles/configs/{fname}\n'
                f'    key: {path}\n'
                f'    value: {val}\n'
                f'    ku: {ku}\n'
                f'{ckline}'
                f'    declared: {decl}\n'
            )
    OUT.write_text(HEADER + "\n".join(blocks) + "\n")
    print(f"wrote {OUT}: {len(blocks)} constants")
    for k in ("verified", "source-unverified", "source-suspect", "unsourced"):
        print(f"  {k:18} {tally.get(k, 0)}")


if __name__ == "__main__":
    main()
