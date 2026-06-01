#!/usr/bin/env python3
"""Build a MINI Obsidian vault from the 3 seeded Contract-Graph chains.

Demo / feasibility test only (PI 2026-06-01): does an Obsidian graph view of the
RAG v2 contract-graph add value over Notion's relation tables? This emits ~30
notes (the KU-3.5 cortex-tension chain, the KU-5.1 alias split, and the ECM
substrate chain) as markdown + YAML frontmatter + [[wikilinks]].

Data is hardcoded from the live Notion DBs (no API needed for the demo). The
full export (~240 nodes) would read Notion via MCP/API instead — this script is
the seed of that exporter.

Link convention: link text carries the relation TYPE so the graph stays legible,
e.g. body line "supports:: [[KB-3.5]]". Node type lives in frontmatter `type`.

Run:  python build_demo_vault.py
Then: open ffn_sim/outputs/obsidian_rag_demo/vault as an Obsidian vault, open
      Graph view, color by `type` (Settings > group: tag/path), and judge.
"""
from __future__ import annotations
import pathlib

VAULT = pathlib.Path(__file__).parent / "vault"
VAULT.mkdir(parents=True, exist_ok=True)

# ---- node tables: (filename, type, frontmatter dict, body lines) -------------
# body lines may contain typed wikilinks like "supports:: [[KB-3.5]]"

NOTES: list[tuple[str, dict, list[str]]] = []

def note(name, ntype, fm, links):
    fm = {"type": ntype, **fm}
    NOTES.append((name, fm, links))

# === SourceEvidence (papers) — link UP to the claims they support =============
SE = {
    "Salbreux2012_TCB":   ("Salbreux 2012 Trends Cell Biol 22(10):536", ["KB-3.5", "KB-3.1"]),
    "ChughPaluch2018_JCS":("Chugh & Paluch 2018 J Cell Sci 131:jcs186254", ["KB-3.5", "KB-3.1"]),
    "Tinevez2009_PNAS":   ("Tinevez 2009 PNAS 106(44):18581", ["KB-3.5"]),
    "Murrell2015_NRMCB":  ("Murrell 2015 Nat Rev Mol Cell Biol 16:486", ["KB-3.5", "KB-3.1"]),
    "Bi2016_NatPhys":     ("Bi 2016 Nat Phys 12:1074 (jamming p0*=3.81)", ["KB-5.1"]),
    "Park2015_NatMater":  ("Park 2015 Nat Mater 14:1040 (unjamming)", ["KB-5.1"]),
    "Bieling2016_Cell":   ("Bieling 2016 Cell 164:115 (Arp2/3 branching)", ["H5-LP-1"]),
    "Discher2005_Science":("Discher 2005 Science 310:1139 (tissue stiffness)", ["KB-1.5"]),
    "Janmey2019_PhysiolRev":("Janmey 2019 Physiol Rev (tissue mechanics)", ["KB-1.5"]),
    "Palchesko2012_PLoSOne":("Palchesko 2012 PLoS ONE 7:e51499 (PAA gel)", ["KB-1.5"]),
}
for key, (short, claims) in SE.items():
    links = [f"supports:: [[{c}]]" for c in claims]
    note(key, "SourceEvidence", {"short_source": short}, links)

# === KnowledgeClaim (literature atoms) — link to evidence + adopted-by =========
note("KB-3.5", "KnowledgeClaim",
     {"kb_id": "KB-3.5", "alias": "KU-3.5", "unit": "Unit 3 Cell",
      "confidence": "High", "status": "verified",
      "value": "gamma_cortex 0.1-1 mN/m; active+passive"},
     ["Cortical tension = active (myosin) + passive (elastic).",
      "evidence:: [[Salbreux2012_TCB]] [[ChughPaluch2018_JCS]] [[Tinevez2009_PNAS]] [[Murrell2015_NRMCB]]",
      "adopted-by:: [[MC-H3-cortex-tension-phase1]]"])
note("KB-3.1", "KnowledgeClaim",
     {"kb_id": "KB-3.1", "alias": "KU-3.1", "unit": "Unit 3 Cell",
      "confidence": "High", "status": "verified",
      "value": "cortex ~200nm, gamma 0.1-1 mN/m"},
     ["Cell cortex mechanical parameters.",
      "evidence:: [[Salbreux2012_TCB]] [[ChughPaluch2018_JCS]] [[Murrell2015_NRMCB]]"])
note("KB-5.1", "KnowledgeClaim",
     {"kb_id": "KB-5.1", "alias": "KU-5.1 (KB/collective)", "unit": "Unit 5 Collective",
      "confidence": "High", "status": "verified",
      "value": "p0=P/sqrt(A); threshold 3.81"},
     ["Jamming transition shape-index criterion (collective).",
      "NOTE: same legacy alias KU-5.1 as [[H5-LP-1]] but a DIFFERENT object.",
      "evidence:: [[Bi2016_NatPhys]] [[Park2015_NatMater]]"])
note("H5-LP-1", "KnowledgeClaim",
     {"kb_id": "H5-LP-1", "alias": "KU-5.1 (code/lamellipodium)", "unit": "Unit 3 Cell",
      "confidence": "High", "status": "verified",
      "value": "dendritic actin density (branched)"},
     ["Lamellipodium dendritic actin density.",
      "NOTE: same legacy alias KU-5.1 as [[KB-5.1]] but a DIFFERENT object (collision fix).",
      "evidence:: [[Bieling2016_Cell]]"])
note("KB-1.5", "KnowledgeClaim",
     {"kb_id": "KB-1.5", "alias": "KU-1.5", "unit": "Unit 1 ECM",
      "confidence": "High", "status": "verified",
      "value": "PAA 0.1-50 kPa; Phase-1 default E=5 kPa nu=0.45"},
     ["Substrate stiffness biological range + Phase-1 5 kPa default.",
      "evidence:: [[Discher2005_Science]] [[Janmey2019_PhysiolRev]] [[Palchesko2012_PLoSOne]]",
      "adopted-by:: [[MC-substrate-phase1]]"])

# === ModelContract — adopts claim, parameterized-by param, tested-by gate ======
note("MC-H3-cortex-tension-phase1", "ModelContract",
     {"mc_id": "MC-H3-cortex-tension-phase1", "phase": "Phase 1", "status": "implemented"},
     ["Phase-1 cortex-only tension gate (band 0.35-0.65 mN/m).",
      "Caveat: composite re-derive after H.8/H.9 -> [[MC-H3-composite-tension]]",
      "adopts:: [[KB-3.5]]",
      "parameterizes:: [[gamma_cortex]]",
      "tested-by:: [[VG-H3-KU35-cortex-tension]]",
      "implemented-in:: [[CM-myosin-gripwalk]] [[CM-cortex-tension]]",
      "governed-by:: [[DEC-grip-walk]]"])
note("MC-substrate-phase1", "ModelContract",
     {"mc_id": "MC-substrate-phase1", "phase": "Phase 1", "status": "implemented"},
     ["Phase-1 substrate = PAA-like linear elastic E=5 kPa, nu=0.45.",
      "adopts:: [[KB-1.5]]",
      "parameterizes:: [[E_substrate]]"])
note("MC-H3-composite-tension", "ModelContract",
     {"mc_id": "MC-H3-composite-tension", "phase": "Phase 1", "status": "draft"},
     ["KU-3.5/3.1 composite tension re-derivation (membrane+cortex+nucleus).",
      "Blocked on H.8/H.9 code.",
      "tested-by:: [[VG-H3-composite-tension]]"])

# === Parameter ================================================================
note("gamma_cortex", "Parameter",
     {"param_id": "PARAM-gamma_cortex", "default": "0.5 mN/m", "status": "implemented",
      "config": "ffn_sim/configs/phase1_h3.yaml"},
     ["Cortex tension parameter.",
      "source-claim:: [[KB-3.5]]", "contract:: [[MC-H3-cortex-tension-phase1]]"])
note("E_substrate", "Parameter",
     {"param_id": "PARAM-E_substrate", "default": "5 kPa", "status": "draft"},
     ["Substrate Young's modulus (sweep).",
      "source-claim:: [[KB-1.5]]", "contract:: [[MC-substrate-phase1]]"])

# === ValidationGate ===========================================================
note("VG-H3-KU35-cortex-tension", "ValidationGate",
     {"vg_id": "VG-H3-KU35-cortex-tension", "type": "Validation", "status": "failing",
      "band": "0.35-0.65 mN/m"},
     ["gamma_total ~3e-4 mN/m (~1600x under band).",
      "tests-contract:: [[MC-H3-cortex-tension-phase1]]",
      "observed-in:: [[RUN-ku35-v4-smoke]]"])
note("VG-H3-composite-tension", "ValidationGate",
     {"vg_id": "VG-H3-composite-tension", "type": "Validation", "status": "blocked"},
     ["Composite tension gate; blocked on H.8/H.9.",
      "tests-contract:: [[MC-H3-composite-tension]]"])

# === CodeMapping (pointers to git) ============================================
note("CM-myosin-gripwalk", "CodeMapping",
     {"cm_id": "CM-myosin-gripwalk", "path": "ffn_sim/cortex/myosin.py", "status": "draft"},
     ["myosin grip-walk (replaces binned_r0 proxy).",
      "implements:: [[MC-H3-cortex-tension-phase1]]",
      "gate:: [[VG-H3-KU35-cortex-tension]]"])
note("CM-cortex-tension", "CodeMapping",
     {"cm_id": "CM-cortex-tension", "path": "configs/phase1_h3.yaml + cortex/cortex.py", "status": "implemented"},
     ["cortex tension config + measurement.",
      "implements:: [[MC-H3-cortex-tension-phase1]]",
      "gate:: [[VG-H3-KU35-cortex-tension]]"])

# === RunResult (pointer to disk) ==============================================
note("RUN-ku35-v4-smoke", "RunResult",
     {"run_id": "RUN-ku35-v4-smoke", "outcome": "smoke", "commit": "6e1c6ac",
      "artifact": "outputs/h3/production/ku35_v4_smoke_diagnosis.md"},
     ["gamma_total ~3e-4 mN/m (~1600x under). FA clutch never loads + myosin step_advances=0.",
      "gate:: [[VG-H3-KU35-cortex-tension]]"])

# === DecisionLedger ===========================================================
note("DEC-grip-walk", "DecisionLedger",
     {"dec_id": "DEC-2026-05-31-grip-walk", "decided_by": "PI", "status": "PI-ratified"},
     ["Replace binned_r0 proxy with AFINES grip-walk.",
      "affects-contract:: [[MC-H3-cortex-tension-phase1]]",
      "affects-gate:: [[VG-H3-KU35-cortex-tension]]"])
note("DEC-ku5-optionA", "DecisionLedger",
     {"dec_id": "DEC-2026-05-31-ku5-optionA", "decided_by": "PI", "status": "PI-ratified"},
     ["KU-5.x collision -> Option A (alias split, no renumber).",
      "relates:: [[KB-5.1]] [[H5-LP-1]]"])

# ---- emit --------------------------------------------------------------------
def fm_block(fm: dict) -> str:
    lines = ["---"]
    for k, v in fm.items():
        lines.append(f"{k}: {v}")
    lines.append("---")
    return "\n".join(lines)

for name, fm, body in NOTES:
    text = fm_block(fm) + "\n\n# " + name + "\n\n" + "\n\n".join(body) + "\n"
    (VAULT / f"{name}.md").write_text(text, encoding="utf-8")

print(f"wrote {len(NOTES)} notes to {VAULT}")
# tally by type
from collections import Counter
c = Counter(fm["type"] for _, fm, _ in NOTES)
for t, n in sorted(c.items()):
    print(f"  {t}: {n}")
