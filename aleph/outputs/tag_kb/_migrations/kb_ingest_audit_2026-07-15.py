#!/usr/bin/env python
"""Push the 2026-07-15 missing-inventory audit payload into the Notion KB (SoT).

Reads `kb_ingest_payload_2026-07-15.json` (from build_audit_payload_2026-07-15.py) and
creates SourceEvidence + KnowledgeClaim rows in the live Contract-Graph.

SAFETY (citation-integrity + SoT-hygiene, all HARD):
  * Every DOI is CrossRef-verified before upload. A DOI that does not resolve is SKIPPED
    (this auto-rejects the ~8 hallucinated/wrong DOIs the manifest flagged for DISCARD —
    Cheng2020 Science, Robert2015 FASEB, etc. — without needing a hand-maintained blocklist).
  * SE de-dup by normalized DOI AND citation key against the live SE DB (no dupes of the 376).
  * KU created as Status='draft' with a PROVISIONAL 'KB-DRAFT-*' id — PI renumbers at ingest.
    De-dup by normalized title against the live KU DB.
  * Every row carries the marker 'AUDIT-2026-07-15' so the whole batch is trivially findable
    and reversible.
  * --dry-run is the DEFAULT. --commit performs the writes. Rate-limited.

Usage:
    python kb_ingest_audit_2026-07-15.py            # dry-run: print the exact plan
    python kb_ingest_audit_2026-07-15.py --commit   # create the rows
    python kb_ingest_audit_2026-07-15.py --commit --se-only   # papers only
"""
from __future__ import annotations
import json, re, sys, time
from pathlib import Path
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "obsidian_rag_full"))
from notion_to_obsidian import get_token, headers, API, DATA_SOURCES  # noqa: E402

SE_DB = DATA_SOURCES["SourceEvidence"]
KC_DB = DATA_SOURCES["KnowledgeClaim"]
PAYLOAD = Path(__file__).with_name("kb_ingest_payload_2026-07-15.json")
MARKER = "AUDIT-2026-07-15"
TODAY = "2026-07-15"


def norm_doi(d: str | None) -> str:
    if not d:
        return ""
    return re.sub(r"^https?://(dx\.)?doi\.org/", "", d.strip(), flags=re.I).lower().rstrip(".")


def norm_title(t: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (t or "").lower())[:60]


def crossref_ok(doi: str) -> tuple[bool, str, str]:
    """Return (resolves, title, year). A non-resolving DOI is a DISCARD."""
    try:
        r = requests.get(f"https://api.crossref.org/works/{doi}",
                         headers={"User-Agent": "ffn_cellsim-kb/1.0 (mailto:sungwook999@gmail.com)"},
                         timeout=20)
        if r.status_code != 200:
            return False, "", ""
        m = r.json()["message"]
        title = (m.get("title") or [""])[0]
        year = ""
        for k in ("published-print", "published-online", "issued", "created"):
            dp = m.get(k, {}).get("date-parts", [[None]])
            if dp and dp[0] and dp[0][0]:
                year = str(dp[0][0]); break
        return True, title, year
    except Exception:
        return False, "", ""


def load_existing_se(tok):
    dois, keys, cursor = {}, set(), None
    while True:
        body = {"page_size": 100}
        if cursor:
            body["start_cursor"] = cursor
        r = requests.post(f"{API}/databases/{SE_DB}/query", headers=headers(tok), json=body, timeout=30)
        r.raise_for_status()
        data = r.json()
        for pg in data["results"]:
            p = pg["properties"]
            ck = "".join(x["plain_text"] for x in p.get("Citation Key", {}).get("title", []))
            if ck:
                keys.add(ck)
            d = norm_doi(p.get("DOI", {}).get("url"))
            if d:
                dois[d] = ck
        if not data.get("has_more"):
            break
        cursor = data["next_cursor"]; time.sleep(0.15)
    return dois, keys


def load_existing_ku_titles(tok):
    titles, cursor = set(), None
    while True:
        body = {"page_size": 100}
        if cursor:
            body["start_cursor"] = cursor
        r = requests.post(f"{API}/databases/{KC_DB}/query", headers=headers(tok), json=body, timeout=30)
        r.raise_for_status()
        data = r.json()
        for pg in data["results"]:
            t = "".join(x["plain_text"] for x in pg["properties"].get("Claim", {}).get("title", []))
            titles.add(norm_title(t))
        if not data.get("has_more"):
            break
        cursor = data["next_cursor"]; time.sleep(0.15)
    return titles


def rt(s):  # rich_text helper (Notion caps 2000 chars/text)
    return {"rich_text": [{"text": {"content": (s or "")[:1900]}}]}


def create_se(tok, key, doi, note):
    props = {
        "Citation Key": {"title": [{"text": {"content": key}}]},
        "Notes": rt(note),
        "Anchor Status": rt(f"auto-ingest {TODAY} ({MARKER} missing-inventory audit); "
                            f"Source Type + Claims unclassified — PI to classify"),
        "DOI": {"url": f"https://doi.org/{doi}"},
    }
    r = requests.post(f"{API}/pages", headers=headers(tok),
                      json={"parent": {"database_id": SE_DB}, "properties": props}, timeout=30)
    r.raise_for_status(); return r.json()["id"]


def create_ku(tok, k):
    props = {
        "Claim": {"title": [{"text": {"content": k["title"]}}]},
        "KB ID": rt(k["prov_kb_id"]),
        "Status": {"select": {"name": "draft"}},
        "Confidence": {"select": {"name": k.get("confidence", "Medium")}},
        "Unit": {"select": {"name": k["unit_partition"]}},
        "Value Range SI": rt(k.get("value_range_si", "")),
        "Subtopic": rt(k.get("compartment", "")),
        "Assumptions": rt(f"{MARKER} DRAFT ({k.get('kb_status','')}). {k.get('collision_flag','')}"),
        "Aliases": rt(f"{MARKER}; gap: {k.get('gap','')}"),
        "Citations (text)": rt(k.get("citations_text", "")),
    }
    r = requests.post(f"{API}/pages", headers=headers(tok),
                      json={"parent": {"database_id": KC_DB}, "properties": props}, timeout=30)
    r.raise_for_status(); return r.json()["id"]


def main():
    commit = "--commit" in sys.argv
    se_only = "--se-only" in sys.argv
    ku_only = "--ku-only" in sys.argv
    pay = json.loads(PAYLOAD.read_text())
    tok = get_token()

    # ---------- SourceEvidence ----------
    se_plan, se_skip, se_discard = [], [], []
    if not ku_only:
        print("querying live SourceEvidence for de-dup ...", flush=True)
        se_dois, se_keys = load_existing_se(tok)
        print(f"  {len(se_dois)} SE rows have DOIs ({len(se_keys)} total)")
        used_keys = set(se_keys)
        for s in pay["se_new"]:
            d = norm_doi(s["doi"])
            if not d:
                continue
            if d in se_dois:
                se_skip.append((s["citation_key"], d, "DOI already in SE")); continue
            ok, cr_title, cr_year = crossref_ok(d)
            time.sleep(0.2)
            if not ok:
                se_discard.append((s["citation_key"], d, "CrossRef 404 — hallucinated/wrong DOI (DISCARD)")); continue
            key = s["citation_key"] or "Ref"
            base, i = key, 2
            while key in used_keys:
                key = f"{base}-{i}"; i += 1
            used_keys.add(key)
            se_plan.append({**s, "citation_key": key, "cr_title": cr_title, "cr_year": cr_year})

    # ---------- KnowledgeClaim ----------
    ku_plan, ku_skip = [], []
    if not se_only:
        print("querying live KnowledgeClaim for de-dup ...", flush=True)
        ku_titles = load_existing_ku_titles(tok)
        print(f"  {len(ku_titles)} existing KU titles")
        seen_local = set()
        for k in pay["ku_draft"]:
            nt = norm_title(k["title"])
            if nt in ku_titles or nt in seen_local:
                ku_skip.append((k["prov_kb_id"], k["title"][:50], "title already in KU")); continue
            seen_local.add(nt)
            ku_plan.append(k)

    # ---------- report ----------
    print(f"\n=== PLAN ===")
    print(f"SourceEvidence: {len(se_plan)} create | {len(se_skip)} skip(dup) | {len(se_discard)} DISCARD(bad DOI)")
    for p in se_plan[:200]:
        print(f"  + SE {p['citation_key']:26} [{p['doi_confidence']:18}] {p['doi']}")
    if se_discard:
        print("  --- DISCARDED (DOI did not resolve) ---")
        for k, d, why in se_discard:
            print(f"  x {k:26} {d}  {why}")
    print(f"\nKnowledgeClaim: {len(ku_plan)} create(draft) | {len(ku_skip)} skip(dup)")
    for p in ku_plan[:200]:
        print(f"  + KU {p['prov_kb_id']:20} [{p['unit_partition']:20}] {p['title'][:60]}")

    if not commit:
        print(f"\n[dry-run] re-run with --commit to create "
              f"{len(se_plan)} SE + {len(ku_plan)} KU rows.")
        return

    if not se_only:
        print(f"\ncreating {len(se_plan)} SourceEvidence rows ...", flush=True)
        for p in se_plan:
            note = (f"Paper: {p.get('cr_title') or p['title']} ({p.get('cr_year','?')}). "
                    f"{MARKER} missing-inventory audit — unblocks: {p['gap'][:120]} "
                    f"[{p['compartment'][:40]}]. DOI CrossRef-verified.")
            create_se(tok, p["citation_key"], p["doi"], note)
            print(f"  + SE {p['citation_key']}"); time.sleep(0.34)
    if not ku_only and not se_only or ku_only:
        pass
    if not se_only:
        print(f"\ncreating {len(ku_plan)} KnowledgeClaim draft rows ...", flush=True)
        for k in ku_plan:
            create_ku(tok, k)
            print(f"  + KU {k['prov_kb_id']}  {k['title'][:50]}"); time.sleep(0.34)
    print(f"\nDONE. Marker='{MARKER}'. Next: bash outputs/tag_kb/refresh.sh (TAG picks up new rows).")


if __name__ == "__main__":
    main()
