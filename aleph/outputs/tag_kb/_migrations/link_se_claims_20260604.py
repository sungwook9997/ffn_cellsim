#!/usr/bin/env python3
"""One-shot: link the 14 SE rows created 2026-06-04 to their KnowledgeClaim(s).

PI-authorised (2026-06-04) claim linkage. Each SE -> KU mapping is content-based
(paper topic vs KU statement), high-confidence only; weak matches are left
UNLINKED (TAG-searchable but no graph edge) rather than fabricated. Sets the
SourceEvidence `Claims` relation (Notion auto-mirrors to the KC `evidence` side).

Dry-run by default; --commit to PATCH Notion.
"""
from __future__ import annotations
import sys, time, pathlib
import requests

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "obsidian_rag_full"))
from notion_to_obsidian import get_token, headers, API  # noqa: E402

KU = {
    "KB-3.5":   "371120daec5d815ab3e1ef281fc54cd4",  # Cortical tension = active+passive
    "KB-4.1":   "372120daec5d81af8e53dbfa4c8265a6",  # E-cadherin junction structure
    "KB-4.4":   "372120daec5d8123933ee62bd5424f8a",  # Junction tension vs cortex tension (Young)
    "KB-6.1.3": "372120daec5d8145a3d7fce3cafc9039",  # Cortical tension + contractility (breast)
    "KB-5.13":  "372120daec5d816f94acc5e3e1813ed9",  # Spheroid surface tension
}

# SE page id -> (citation_key, [KU keys]).  Rationale per row in the comment.
MAP = {
    "375120daec5d81abb10cfeef53c0965c": ("Tao2020_BiophysicalJournal",        ["KB-3.5"]),           # cell tension drives motility
    "375120daec5d815da433e98f3361743e": ("Tao2019_Jour",                      ["KB-3.5"]),           # preprint of Tao2020
    "375120daec5d81d0b10de3192e31bc40": ("Arslan2024_CurrentBiology",         ["KB-4.1", "KB-4.4"]), # E-cad contacts + cortical flow
    "375120daec5d8178a82bc20f3a7650e2": ("Samarage2015_DevelopmentalCell",    ["KB-3.5", "KB-4.4"]), # cortical tension allocates cells
    "375120daec5d81d6b2e7c34e5a48ae37": ("Miroshnikova2018_NatCellBiol",      ["KB-3.5", "KB-4.4"]), # adhesion+cortical tension couple
    "375120daec5d818ab266da5a01f44b49": ("Cornell2025_NatCellBiol",           ["KB-3.5"]),           # target-cell cortical tension
    "375120daec5d8159ad7bd7f845467373": ("Dmitrieff2017_ProcNatlAcadSciUSA",  ["KB-3.5"]),           # cortical tension vs MT stiffness
    "375120daec5d813aaeb7c9fad576c102": ("Warmt2021_NewJPhys",                ["KB-3.5", "KB-6.1.3"]),# breast-cell cortical contractility
    "375120daec5d81809f53c6a6bfbb20f6": ("Moazzeni2021_PhysRevE",             ["KB-3.5"]),           # single-cell tension measurement
    "375120daec5d816eac63c1dab967e524": ("Wang2021_IEEESensorsJ",             ["KB-3.5"]),           # cortical tension via constriction
    "375120daec5d8105ac7ce99c0c2e99b9": ("Winklbauer2015_JournalOfCellScience",["KB-3.5", "KB-4.4"]),# adhesion strength from cortical tension
    "375120daec5d8199a8b9ee68fd15fbe0": ("Bohec2025_Jour",                    ["KB-3.5"]),           # RhoGTPase -> cortical tension
    "375120daec5d8169a54dedec9e9409ca": ("Thiticharoentam2026_Jour",          ["KB-4.4", "KB-5.13"]),# cell-cell interfacial tension / sorting
    # Fu2019_ColloidsAndSurfacesBBioi (375120daec5d8132a0e1f7ed693858ca): spheroid-formation
    #   methods / CSC enrichment -> no high-confidence KU; left UNLINKED (TAG-only).
}


def patch_claims(tok, se_id, ku_ids):
    body = {"properties": {"Claims": {"relation": [{"id": k} for k in ku_ids]}}}
    r = requests.patch(f"{API}/pages/{se_id}", headers=headers(tok), json=body, timeout=30)
    r.raise_for_status()


def main():
    commit = "--commit" in sys.argv
    tok = get_token()
    print(f"{'COMMIT' if commit else 'DRY-RUN'}: linking {len(MAP)} SE rows "
          f"(Fu2019 intentionally left unlinked)\n")
    for se_id, (ck, kus) in MAP.items():
        missing = [k for k in kus if k not in KU]
        if missing:
            print(f"  !! {ck}: unknown KU {missing}"); continue
        print(f"  {ck:38} -> {', '.join(kus)}")
        if commit:
            patch_claims(tok, se_id, [KU[k] for k in kus])
            time.sleep(0.3)
    if not commit:
        print("\n[dry-run] re-run with --commit to write the Claims relations.")
    else:
        print("\ndone. next: refresh kb.duckdb + Obsidian to render the new edges.")


if __name__ == "__main__":
    main()
