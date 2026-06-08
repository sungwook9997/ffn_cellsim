"""Link the 2026-06-04 architecture-batch SE rows to KnowledgeClaims (single-cell cortex + Layer-2).

Mirrors link_se_claims_20260604.py: patches each SE page's "Claims" relation to the KU pages.
Dry-run by default; --commit to write. PI 2026-06-04: "layer2도 확인하고 업데이트" → architecture
papers get cortex KUs AND the Layer-2-relevant ones also get Layer-2 KUs.
"""
import sys
from notion_to_obsidian import get_token, headers, API  # noqa: E402
import requests

# KU kb_id -> Notion page id
KU = {
    "KB-3.23": "372120daec5d81749e14c8e04c451f3e",   # cortex fiber-network model
    "KB-3.5":  "371120daec5d815ab3e1ef281fc54cd4",   # cortical tension active+passive
    "KB-3.19": "372120daec5d8145af3edc1a5deb6946",   # crosslinker stochastic dynamics
    "KB-3.B1.4":"372120daec5d81fbae01c3b2e6a64f74",  # membrane tether / cortex-membrane adhesion
    "KB-5.13": "372120daec5d816f94acc5e3e1813ed9",   # spheroid surface tension (Layer-2)
    "KB-5.16": "372120daec5d81528ae2d936695c6799",   # Unit-5 collective spheroid (Layer-2)
    "KB-PIV-9":"372120daec5d81498a2bf3819ba7a382",   # Layer-2 collective framework
    "KB-4.18": "372120daec5d81929477ce4478d2d3ad",   # cortex-cortex collective coupling (Layer-2)
    "KB-6.1.3":"372120daec5d8145a3d7fce3cafc9039",   # cortical tension breast
}
# SE page id -> (label, [kb_ids])
MAP = {
    "375120da-ec5d-81a4-9ade-d683d4d34bcb": ("Flormann2024",   ["KB-3.23","KB-3.5","KB-6.1.3","KB-5.13"]),
    "375120da-ec5d-81c5-88ff-dabc08b19c55": ("Fritzsche2016",  ["KB-3.23","KB-3.5"]),
    "375120da-ec5d-814d-9a8e-f6cf6a629f60": ("Fritzsche2017",  ["KB-3.23"]),
    "375120da-ec5d-8189-9216-fe668126f214": ("Sakamoto2024",   ["KB-3.23"]),
    "375120da-ec5d-819b-87e2-cdb5c9733312": ("Li2022",         ["KB-3.23","KB-3.5"]),
    "375120da-ec5d-814a-8c0c-cb48c6f7b6e9": ("Bacher2021",     ["KB-3.23","KB-3.5"]),
    "375120da-ec5d-81f7-9d52-e5dcedae8d7a": ("Garlick2022",    ["KB-3.23"]),
    "375120da-ec5d-810a-ad24-fa9f4c347972": ("Chen2024NatPhys",["KB-3.23"]),
    "375120da-ec5d-8196-8f85-ea1e68f2d55c": ("Banerjee2021",   ["KB-3.23","KB-3.5"]),
    "375120da-ec5d-8135-af65-c4b92a501aa1": ("Chandrasekaran2024",["KB-3.19"]),
    "375120da-ec5d-81b0-8d7a-e6f2fab9478d": ("Serwas2021",     ["KB-3.B1.4"]),
    "375120da-ec5d-81fa-bda5-c31ca118c17c": ("Kadzik2026",     ["KB-3.23","KB-4.18","KB-5.13"]),
    "375120da-ec5d-8130-b971-d64d3c568d1d": ("MerinoCasallo2022",["KB-5.16","KB-PIV-9"]),
    "375120da-ec5d-8130-9d52-ec06dc547d61": ("Ray2024",        ["KB-3.5"]),
}


def patch(tok, se_id, ku_ids):
    body = {"properties": {"Claims": {"relation": [{"id": KU[k]} for k in ku_ids]}}}
    r = requests.patch(f"{API}/pages/{se_id}", headers=headers(tok), json=body, timeout=30)
    r.raise_for_status()


def main():
    commit = "--commit" in sys.argv
    tok = get_token()
    print(f"{'COMMIT' if commit else 'DRY-RUN'}: link {len(MAP)} architecture SE rows -> KUs\n")
    for se_id, (lab, kus) in MAP.items():
        miss = [k for k in kus if k not in KU]
        if miss:
            print(f"  !! {lab}: unknown KU {miss}"); continue
        l2 = [k for k in kus if k.startswith(("KB-5","KB-PIV","KB-4.18"))]
        print(f"  {lab:22} -> {kus}" + (f"   [Layer-2: {l2}]" if l2 else ""))
        if commit:
            patch(tok, se_id, kus)
    print("\n" + ("done: KC edges written (single-cell cortex + Layer-2)." if commit
                  else "[dry-run] re-run with --commit."))


if __name__ == "__main__":
    main()
