#!/usr/bin/env python3
"""Create Notion SourceEvidence rows for the references/ TAG corpus papers.

This is the missing other half of references_ingest.py: that script LINKS a PDF
to an existing SourceEvidence row by DOI, but never CREATES the row when none
exists. As a result every corpus paper sat with linked_to_source_evidence=false,
so the Obsidian mirror (which reads only Notion) could not see them.

This tool: for the on-topic `source=="new"` corpus papers (default), enriches
each DOI via Crossref into a clean `Author+Year_Journal` Citation Key + real
title (matching the existing SE convention), de-duplicates against the live
SourceEvidence DB (by DOI), and CREATES a minimal SE row for the genuinely-new
ones. After this, `notion_to_obsidian.py` renders them as paper nodes.

Honest-by-construction:
  • Source Type is left BLANK (role classification needs a human) + a Note marks
    the provenance, so nothing is fabricated.
  • DOI de-dup against the live DB; never creates a row whose DOI already exists.
  • Dry-run by default; writes to Notion only with --commit.

RUN:
  conda activate ffn_sim
  python references_to_se.py                 # dry-run: print the plan
  python references_to_se.py --commit        # create the rows in Notion
  python references_to_se.py --commit --bundle   # also include cellpress bundle
"""
from __future__ import annotations
import json, re, sys, time, pathlib
import requests

# reuse the Obsidian exporter's token + API plumbing
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "obsidian_rag_full"))
from notion_to_obsidian import get_token, headers, API, DATA_SOURCES  # noqa: E402

REF_DIR = pathlib.Path("/Users/sw1/ffn_cellsim/ffn_sim/references")
CORPUS_JSON = REF_DIR / "tag_corpus.json"
SE_DB = DATA_SOURCES["SourceEvidence"]
MAILTO = "sungwook999@gmail.com"   # Crossref polite-pool identifier


def norm_doi(s):
    if not s:
        return ""
    s = s.strip().lower()
    s = re.sub(r"^https?://(dx\.)?doi\.org/", "", s)
    return s.rstrip(".)").strip()


def crossref(doi):
    """Return {title, family, year, journal, type} from Crossref, or {} on miss."""
    try:
        r = requests.get(f"https://api.crossref.org/works/{doi}",
                         params={"mailto": MAILTO}, timeout=20)
        if r.status_code != 200:
            return {}
        m = r.json()["message"]
    except Exception:
        return {}
    title = (m.get("title") or [""])[0]
    auth = m.get("author") or []
    family = (auth[0].get("family", "") if auth else "") or ""
    year = ""
    for k in ("published-print", "published-online", "published", "issued"):
        dp = (m.get(k) or {}).get("date-parts")
        if dp and dp[0] and dp[0][0]:
            year = str(dp[0][0]); break
    journal = (m.get("short-container-title") or m.get("container-title") or [""])[0]
    return {"title": title, "family": family, "year": year,
            "journal": journal, "type": m.get("type", "")}


def jabbr(journal):
    """'Nat. Phys.' -> 'NatPhys'; 'PLoS Comput. Biol.' -> 'PLoSComputBiol'."""
    j = re.sub(r"[^A-Za-z0-9 ]", "", journal)
    return "".join(w[:1].upper() + w[1:] for w in j.split())[:24] or "Jour"


def make_key(meta, fallback):
    fam = re.sub(r"[^A-Za-z]", "", meta.get("family", "")) or ""
    yr = meta.get("year", "")
    ab = jabbr(meta.get("journal", ""))
    if fam and yr:
        return f"{fam}{yr}_{ab}"
    return fallback  # corpus citation_key (slug-sha) as last resort


def load_existing_se(tok):
    """{norm_doi -> citation_key} and {citation_key} over the live SE DB."""
    dois, keys = {}, set()
    cursor = None
    while True:
        body = {"page_size": 100}
        if cursor:
            body["start_cursor"] = cursor
        r = requests.post(f"{API}/databases/{SE_DB}/query",
                          headers=headers(tok), json=body, timeout=30)
        r.raise_for_status()
        data = r.json()
        for pg in data["results"]:
            props = pg["properties"]
            ck = "".join(x["plain_text"] for x in props["Citation Key"]["title"])
            keys.add(ck)
            d = norm_doi(props.get("DOI", {}).get("url"))
            if d:
                dois[d] = ck
        if not data.get("has_more"):
            break
        cursor = data["next_cursor"]
        time.sleep(0.15)
    return dois, keys


def create_se_row(tok, key, doi, note):
    props = {
        "Citation Key": {"title": [{"text": {"content": key}}]},
        "Notes": {"rich_text": [{"text": {"content": note}}]},
        "Anchor Status": {"rich_text": [{"text": {"content": "auto-ingest 2026-06-02 (references TAG corpus); Source Type + Claims unclassified"}}]},
    }
    if doi:
        props["DOI"] = {"url": f"https://doi.org/{doi}"}
    r = requests.post(f"{API}/pages", headers=headers(tok),
                      json={"parent": {"database_id": SE_DB}, "properties": props},
                      timeout=30)
    r.raise_for_status()
    return r.json()["id"]


def main():
    commit = "--commit" in sys.argv
    include_bundle = "--bundle" in sys.argv
    tok = get_token()

    corpus = json.loads(CORPUS_JSON.read_text())
    want = [c for c in corpus
            if c.get("source") == "new" or (include_bundle and c.get("source") == "bundle")]

    # de-dup the corpus itself by DOI (the manifest has a few repeats)
    by_doi, no_doi = {}, []
    for c in want:
        d = norm_doi(c.get("doi"))
        if d:
            by_doi.setdefault(d, c)   # first wins
        else:
            no_doi.append(c)
    print(f"{len(want)} candidate corpus papers -> {len(by_doi)} unique DOIs "
          f"({len(want)-len(by_doi)} in-corpus dup), {len(no_doi)} without DOI")

    print("querying live SourceEvidence for de-dup ...", flush=True)
    se_dois, se_keys = load_existing_se(tok)
    print(f"  {len(se_dois)} SE rows have DOIs ({len(se_keys)} total rows)")

    plan_new, plan_skip = [], []
    used_keys = set(se_keys)
    for d, c in sorted(by_doi.items()):
        if d in se_dois:
            plan_skip.append((d, se_dois[d], "already in SE"))
            continue
        meta = crossref(d)
        time.sleep(0.2)
        key = make_key(meta, c["citation_key"])
        # ensure uniqueness
        base, i = key, 2
        while key in used_keys:
            key = f"{base}-{i}"; i += 1
        used_keys.add(key)
        title = meta.get("title") or c.get("title") or ""
        plan_new.append({"doi": d, "key": key, "title": title[:90],
                         "year": meta.get("year", "?"), "src": c["path"]})

    print(f"\n=== PLAN: {len(plan_new)} NEW SE rows, "
          f"{len(plan_skip)} skipped (already in SE) ===\n")
    for p in plan_new:
        print(f"  + {p['key']:28} {p['title']}")
    if plan_skip:
        print("\n  skipped (DOI already in SE):")
        for d, ck, _ in plan_skip:
            print(f"    = {ck}  ({d})")
    if no_doi:
        print(f"\n  {len(no_doi)} corpus papers without DOI (not created):")
        for c in no_doi:
            print(f"    ? {c['title'][:70]}")

    if not commit:
        print(f"\n[dry-run] re-run with --commit to create the {len(plan_new)} rows.")
        return

    print(f"\ncreating {len(plan_new)} SourceEvidence rows ...", flush=True)
    created = {}
    for p in plan_new:
        note = (f"Paper: {p['title']} ({p['year']}). Ingested from references/ "
                f"TAG corpus 2026-06-02 (spheroid / PI-exp validation track), "
                f"file {p['src']}.")
        pid = create_se_row(tok, p["key"], p["doi"], note)
        created[p["doi"]] = (p["key"], pid)
        print(f"  + {p['key']}")
        time.sleep(0.3)

    # reflect linkage back into tag_corpus.json
    linked = set(created) | {d for d, _, _ in plan_skip}
    for c in corpus:
        d = norm_doi(c.get("doi"))
        if d in created:
            c["linked_to_source_evidence"] = True
            c["se_citation_key"] = created[d][0]
            c["se_uid"] = created[d][1]
        elif d in {x[0] for x in plan_skip}:
            c["linked_to_source_evidence"] = True
            c["se_citation_key"] = se_dois[d]
    CORPUS_JSON.write_text(json.dumps(corpus, indent=2, ensure_ascii=False))
    print(f"\ndone: created {len(created)} SE rows; tag_corpus.json linkage updated.")
    print("next: re-run notion_to_obsidian.py (or refresh.sh) to render them as nodes.")


if __name__ == "__main__":
    main()
