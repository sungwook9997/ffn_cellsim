#!/usr/bin/env python3
"""Link the 30 pre-existing ORPHAN SourceEvidence rows (no Claims) to a KU.

These were created by earlier-session references_to_se runs (spheroid / cell-
spreading / FA corpus) and never claim-linked, so they show as orphan nodes in
Obsidian. PI-authorised (2026-06-04) content-based linkage, single best-fit KU
per paper. Reviews/landscape papers are linked low-confidence to the spheroid
hub KB-5.15 (they support the choice of spheroid as validation system) — marked
`# weak` below; PI can refine. The off-topic Balog2007 is NOT linked (flagged
for deletion in Notion UI instead — MCP cannot archive).

Sets SourceEvidence `Claims` (auto-mirrors to KC `evidence`). --commit to write.
"""
from __future__ import annotations
import sys, time, pathlib
import requests

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "obsidian_rag_full"))
from notion_to_obsidian import get_token, headers, API  # noqa: E402

KU = {
    "KB-5.15": "372120daec5d815ebfb3d7562eae7cb1",  # Spheroid radial expansion (primary observable)
    "KB-5.13": "372120daec5d816f94acc5e3e1813ed9",  # Spheroid surface tension
    "KB-5.3":  "372120daec5d81a6a498c03ea52dc0ef",  # Spheroid wetting (liquid-drop)
    "KB-5.14": "372120daec5d81e5bd1ff676f6296687",  # Compressive stress unjamming
    "KB-2.2":  "372120daec5d8176b628e26fce81d6f9",  # FA maturation stages
    "KB-2.4":  "372120daec5d81bead96e9d82fd9796d",  # Chan-Odde motor-clutch
    "KB-2.6":  "372120daec5d8130b752f235f64421ef",  # Talin rod unfolding
    "KB-2.12": "372120daec5d817b892ac25bcbf67c80",  # Traction force per FA
    "KB-2.14": "372120daec5d81db8e84f2d4b149e06d",  # Contact guidance via anisotropic FA
    "KB-2.15": "372120daec5d812da5c8eea09f76d056",  # Unified motility force
    "KB-3.11": "372120daec5d81c1b9a5d1dfd6dd8525",  # Migration modes by ECM/confinement
    "KB-3.22": "372120daec5d819b821ef539dcc5c6dd",  # Full-fidelity validation methodology
}

# SE page id -> (citation_key, KU key)
MAP = {
    "373120daec5d8185b4d9e92338351886": ("Arruda2026_NpjSystBiolAppl",          "KB-3.11"),  # sim-inference of migration
    "373120daec5d813b83f4f7de0b00cc84": ("Audoin2022_SciRep",                   "KB-5.14"),  # spheroids accelerate invasion (unjamming)
    "373120daec5d8119abc3ef94e4662da8": ("Berens2015_JoVE",                     "KB-5.15"),  # spheroid invasion assay  # weak/method
    "373120daec5d81d09700d9c44ccf5454": ("Betorz2023_JournalOfTheMechanicsAnd", "KB-2.15"),  # cell spreading/migration/taxis
    "373120daec5d81358ee3ec88e0a9748e": ("Boot2021_AdvancesInPhysicsX",         "KB-5.13"),  # spheroid mechanics / surface tension review
    "373120daec5d8186af85f71b8d6a0004": ("Botticelli2025_FrontBioengBiotechnol","KB-5.15"),  # spheroid growth model
    "373120daec5d813eaf43dad1320b8802": ("Brningk2019_SciRep",                  "KB-5.15"),  # CA spheroid model  # weak
    "373120daec5d8127b52ecc5bc8c63d1b": ("Bustamante2021_Biofabrication",       "KB-5.13"),  # spheroid fusion = surface tension
    "373120daec5d811f955ac83573194957": ("Chen2017_MBE",                        "KB-5.15"),  # multiscale tumor spheroid
    "373120daec5d81899178f41fafd56ef3": ("De2018_CommunBiol",                   "KB-2.14"),  # FA orientation / cyclic strain
    "373120daec5d8174822dfcfe7ca533ac": ("Dhandapani2023_AdvHealthcareMaterials","KB-5.15"), # 3D spheroid model  # weak/review
    "373120daec5d81c2b18be716745709c5": ("Fang2016_PhysRevE",                   "KB-2.12"),  # cell-spreading traction forces
    "373120daec5d81ba8277cc5a4c3a1ec7": ("Feng2025_BiochemicalEngineeringJo",   "KB-5.15"),  # spheroid maturity / invasion
    "373120daec5d81c2a3c7ff06a5eab2d4": ("Geiger2022_PLoSONE",                  "KB-5.15"),  # directed invasion of spheroids
    "373120daec5d819194efea3dfcde1f2c": ("Herold2023_PLoSComputBiol",           "KB-3.22"),  # sim-vs-exp scoring (validation method)
    "373120daec5d81c1b089e4ca75517b8f": ("Hou2018_SciRep",                      "KB-5.15"),  # spheroid quantification  # weak/method
    "373120daec5d814b8175faf54831be26": ("Kim2012_IntegrBiol",                  "KB-2.4"),   # FA + cytoskeleton + actin motor
    "373120daec5d810597a5ff9ded84bfb9": ("Mangani2025_Cancers",                 "KB-5.15"),  # spheroid breast matrix  # weak/review
    "374120daec5d8116a058d69f2b320b52": ("Mykuliak2020_BiophysicalJournal",     "KB-2.6"),   # mechanical unfolding of proteins (talin)
    "373120daec5d8186a497f2b6d0f3d382": ("Nayak2023_Cancers",                   "KB-5.15"),  # 3D spheroid review  # weak/review
    "373120daec5d81ff8793da83452a6663": ("Nieto2026_EuropeanJournalOfPharmac",  "KB-5.15"),  # nutrient-limited spheroid growth
    "373120daec5d81a8b096eccf8a860aaf": ("Odenthal2013_PLoSComputBiol",         "KB-2.2"),   # initial cell spreading mechanics
    "373120daec5d81bd99ced194a120d0d3": ("Rodenhizer2018_AdvHealthcareMaterials","KB-5.15"), # 3D tumor model landscape  # weak/review
    "373120daec5d81b8859ed4b79fe62ad4": ("Rosenbauer2023_JPhysChemB",           "KB-5.15"),  # multiscale spheroid nutrient
    "373120daec5d816e954afe05aa77516c": ("Shah2025_Cells",                      "KB-5.15"),  # TME spheroid review  # weak/review
    "373120daec5d81748061dfe69c4ba2fb": ("Wang2016_PLoSONE",                    "KB-5.3"),   # spheroid formation in microwells
    "373120daec5d815d87f4d2ee352a9675": ("Xiong2007_NatPrec",                   "KB-2.2"),   # 3D stochastic cell spreading
    "373120daec5d81dcbf34ddddb5d82f7f": ("Zhu2022_Organoids",                   "KB-5.15"),  # spheroid/organoid review  # weak/review
    "373120daec5d815da4eceb7755b87057": ("Zhuo2018_LightSciAppl",               "KB-2.2"),   # FA dynamics measurement
    "375120daec5d8132a0e1f7ed693858ca": ("Fu2019_ColloidsAndSurfacesBBioi",     "KB-5.3"),   # spontaneous spheroid formation
    # OFF-TOPIC — do NOT link (delete in Notion UI; yeast phosphoglycerate kinase):
    #   "373120daec5d81c2b583d0006ab61053": "ZZ-OFFTOPIC Balog2007_BiophysicalJournal"
}


def patch_claims(tok, se_id, ku_id):
    body = {"properties": {"Claims": {"relation": [{"id": ku_id}]}}}
    r = requests.patch(f"{API}/pages/{se_id}", headers=headers(tok), json=body, timeout=30)
    r.raise_for_status()


def main():
    commit = "--commit" in sys.argv
    tok = get_token()
    print(f"{'COMMIT' if commit else 'DRY-RUN'}: {len(MAP)} backlog SE rows "
          f"(Balog2007 off-topic excluded -> delete in Notion)\n")
    for se_id, (ck, ku) in MAP.items():
        if ku not in KU:
            print(f"  !! {ck}: unknown KU {ku}"); continue
        print(f"  {ck:42} -> {ku}")
        if commit:
            patch_claims(tok, se_id, KU[ku]); time.sleep(0.3)
    print("\n[dry-run]" if not commit else "\ndone.",
          "re-run with --commit." if not commit else "refresh kb.duckdb + Obsidian next.")


if __name__ == "__main__":
    main()
