#!/usr/bin/env python3
"""Replace the 3 confirmed-hallucinated SourceEvidence rows with their real
substitutes, in place (the Claims relation to KB-4.1 / KB-3.18,3.19 / KB-6.2.3
is preserved — the claim now points at the correct paper).

Substitutes were web/PubMed-verified (see AUDIT_FINDINGS.md). The Yao->Ferrer
row additionally FLAGS the k_off value discrepancy: KB-3.19's "0.4 s⁻¹" is a
bulk no-load value, not Ferrer's single-molecule 0.066 s⁻¹ — left for PI to rule
on (the claimed literature value is NOT changed here, only the source is fixed).

Updates: Citation Key (title), DOI, Source Type, Short Source, Anchor Status,
Notes. Dry-run by default; `--apply` writes.
"""
from __future__ import annotations

import sys

import duckdb
import requests

HERE = __import__("pathlib").Path(__file__).parent
DB_PATH = HERE / "kb.duckdb"
sys.path.insert(0, str(HERE.parent / "obsidian_rag_full"))
from notion_to_obsidian import get_token, headers  # noqa: E402

API = "https://api.notion.com/v1"

# old citation_key -> replacement fields
REPL = {
    "Yao2011_NatCommun": {
        "key": "Ferrer2008_PNAS",
        "doi": "10.1073/pnas.0706124105",
        "type": "Direct measurement",
        "short": "Ferrer, Lee, Chen, Pelz, Nakamura, Kamm, Lang 2008, PNAS 105:9221 "
                 "(single-molecule alpha-actinin/actin rupture forces)",
        "anchor": "alpha-actinin/actin single-molecule: k_off0=0.066+/-0.028 s^-1 "
                  "(Bell-Evans), rupture 40-80 pN. WARNING prior '0.4 s^-1' is a BULK "
                  "no-load value, NOT single-molecule - KB-3.19 value pending PI.",
        "notes": "[2026-06-02 hallucination audit] Replaced fabricated 'Yao2011_NatCommun' "
                 "(no such Nat Commun alpha-actinin paper exists). Single-molecule anchor = "
                 "Ferrer 2008 PNAS; lifetime companion = Miyata 1996 BBA "
                 "10.1016/0304-4165(96)00003-7. Supports KB-3.18/KB-3.19.",
    },
    "YapKovacs_JCS": {
        "key": "Yap2015_DevCell",
        "doi": "10.1016/j.devcel.2015.09.012",
        "type": "Review consensus",
        "short": "Yap, Gomez, Parton 2015, Dev Cell 35:12 (adherens-junction nanoscale "
                 "organization; review)",
        "anchor": "adherens junction / cadherin nanoassembly organization",
        "notes": "[2026-06-02 hallucination audit] Replaced fabricated 'YapKovacs_JCS' "
                 "(no Yap & Kovacs JCS adherens-junction review exists). Supports KB-4.1.",
    },
    "NanoConvergence2021_Glioma": {
        "key": "Masud2025_SciRep",
        "doi": "10.1038/s41598-025-04841-4",
        "type": "Direct measurement",
        "short": "Masud et al 2025, Sci Rep 15:19302 (glioma whole-cell AFM: T98G 33.5, "
                 "U87 MG 24.6 kPa)",
        "anchor": "glioma AFM whole-cell Young's modulus: T98G 33.47+/-1.67 / "
                  "U87 MG 24.62+/-2.05 kPa",
        "notes": "[2026-06-02 hallucination audit] Replaced mis-cited "
                 "'NanoConvergence2021_Glioma' (PMC8253861 was Ketebo filamin-A/U87 pillar "
                 "paper, not AFM and no T98G). Supports KB-6.2.3.",
    },
}


def to_uuid(n):
    n = n.replace("-", "")
    return f"{n[:8]}-{n[8:12]}-{n[12:16]}-{n[16:20]}-{n[20:]}"


def rt(text):
    return [{"type": "text", "text": {"content": text[:1900]}}]


def main():
    apply = "--apply" in sys.argv
    con = duckdb.connect(str(DB_PATH), read_only=True)
    ids = {ck: nid for ck, nid in
           con.execute("SELECT citation_key, id FROM source_evidence").fetchall()}
    con.close()
    tok = get_token()

    for old_ck, r in REPL.items():
        if old_ck not in ids:
            print(f"  !! {old_ck} not found in kb.duckdb (already replaced?) — skip")
            continue
        pid = to_uuid(ids[old_ck])
        props = {
            "Citation Key": {"title": rt(r["key"])},
            "DOI": {"url": f"https://doi.org/{r['doi']}"},
            "Source Type": {"select": {"name": r["type"]}},
            "Short Source": {"rich_text": rt(r["short"])},
            "Anchor Status": {"rich_text": rt(r["anchor"])},
            "Notes": {"rich_text": rt(r["notes"])},
        }
        print(f"\n{'APPLY' if apply else 'DRY'}: {old_ck}  ->  {r['key']}  ({r['doi']})")
        print(f"    {r['anchor'][:100]}")
        if apply:
            resp = requests.patch(f"{API}/pages/{pid}", headers=headers(tok),
                                  json={"properties": props}, timeout=30)
            print(f"    HTTP {resp.status_code} "
                  + ("OK" if resp.status_code < 300 else resp.text[:200]))
    if not apply:
        print("\n(dry-run) re-run with --apply to write to Notion.")


if __name__ == "__main__":
    main()
