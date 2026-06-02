#!/usr/bin/env python3
"""Apply the audit's DOI fixes to the Notion SourceEvidence DB (Notion = SoT).

Writes the `DOI` (url) property only — NO citation-key renames (keys are
contract-level; IDs never change meaning), NO relation edits. Two classes:

  • CORRECT   — 14 web-verified DOIs (curated; written directly).
  • BACKFILL  — 102 rows audited NO_DOI_FOUND with a CrossRef-suggested DOI.
      These auto-suggestions have a real false-positive rate (a top-hit can match
      author+year by coincidence — e.g. a French-humanities DOI for a "Cancers"
      paper), so each is GATED: the suggested DOI's CrossRef title must share a
      topical content word with the row's own description (anchor_status/notes).
      Failing rows are skipped and reported for manual DOI entry, never written.

Dry-run by default. `--test-one` PATCHes one web-verified row to confirm the
token has write access. `--apply` runs the gated bulk write.

After applying, re-run notion_to_duckdb.py so kb.duckdb reflects the fixes.
"""
from __future__ import annotations

import re
import sys
import time

import duckdb
import requests

HERE = __import__("pathlib").Path(__file__).parent
DB_PATH = HERE / "kb.duckdb"
sys.path.insert(0, str(HERE.parent / "obsidian_rag_full"))
from notion_to_obsidian import get_token, headers  # noqa: E402
from verify_sources import CR, msg_year, norm_doi   # noqa: E402

API = "https://api.notion.com/v1"
SESS = requests.Session()
SESS.headers["User-Agent"] = "ffn_cellsim-doi-fix/1.0 (mailto:sungwook999@gmail.com)"

CORRECT = {  # web-verified DOIs (citation_key -> DOI)
    "Cox_NatCommun": "10.1038/ncomms10366",
    "HuertaLopez2024_SciAdv": "10.1126/sciadv.adf9758",
    "Panzetta2019_PNAS": "10.1073/pnas.1904660116",
    "Andreu2021_NatCommun": "10.1038/s41467-021-24383-3",
    "K562_BBRC2019": "10.1016/j.bbrc.2019.06.054",
    "Nam2016_PNAS": "10.1073/pnas.1523906113",
    "Ray2017_NatCommun": "10.1038/ncomms14923",
    "Yang_NRMCB": "10.1038/s41580-020-0237-9",
    "Bi2015_PRX": "10.1103/PhysRevX.6.021011",
    "Mui2016_PNAS": "10.1242/jcs.183699",
    "MegeIshiyama_JBCReview": "10.1101/cshperspect.a028738",
    "Hosseini2021_BiophysJ": "10.1016/j.bpj.2021.05.006",
    "Liew2024_CMBE": "10.1007/s12195-024-00811-4",
    "Smelser2015_BMMB": "10.1007/s10237-015-0677-x",
}

STOP = set("about above actin actins after along among analysis based between cell "
           "cells cellular during dynamics effect effects forces human imaging living "
           "材料 mechanical mechanics model models models. role roles single studies "
           "study their these through tissue tissues toward using versus via within".split())


def to_uuid(n):
    n = n.replace("-", "")
    return f"{n[:8]}-{n[8:12]}-{n[12:16]}-{n[16:20]}-{n[20:]}"


def doi_url(doi):
    return f"https://doi.org/{norm_doi(doi)}"


def toks(s):
    return {w for w in re.findall(r"[a-z]{5,}", (s or "").lower()) if w not in STOP}


def verify_backfill(doi, ctx):
    """The suggested DOI must resolve AND its title must share a topical word
    with the row's own description — guards against coincidental author+year hits."""
    try:
        r = SESS.get(f"{CR}/{norm_doi(doi)}", timeout=20)
        if r.status_code != 200:
            return False, "DOI did not resolve"
        m = r.json()["message"]
    except Exception:
        return False, "fetch error"
    title = (m.get("title") or [""])[0]
    overlap = toks(title) & toks(ctx)
    if not overlap:
        return False, f"topic mismatch — CrossRef title: {title[:55]!r}"
    return True, "topic ok: " + ",".join(sorted(overlap)[:3])


def patch_doi(tok, page_id, doi):
    r = requests.patch(f"{API}/pages/{page_id}", headers=headers(tok),
                       json={"properties": {"DOI": {"url": doi_url(doi)}}}, timeout=30)
    return r.status_code, (r.text[:200] if r.status_code >= 300 else "")


def build_plan(con):
    plan = []  # (page_id, ck, old_doi, new_doi, kind, ctx)
    rows = con.execute(
        "SELECT se.id, se.citation_key, se.doi, se.anchor_status, se.notes, "
        "a.verdict, a.suggested_doi "
        "FROM source_evidence se LEFT JOIN source_audit a "
        "ON a.citation_key = se.citation_key").fetchall()
    for nid, ck, doi, anc, notes, verdict, sugg in rows:
        ctx = f"{ck} {anc or ''} {notes or ''}"
        if ck in CORRECT:
            new, kind = CORRECT[ck], "CORRECT"
        elif verdict == "NO_DOI_FOUND" and sugg and not doi:
            new, kind = sugg, "BACKFILL"
        else:
            continue
        if norm_doi(new) != norm_doi(doi or ""):
            plan.append((to_uuid(nid), ck, doi or "", norm_doi(new), kind, ctx))
    return plan


def main():
    con = duckdb.connect(str(DB_PATH), read_only=True)
    plan = build_plan(con)
    con.close()
    tok = get_token()

    if "--test-one" in sys.argv:
        pid, ck, old, new, kind, ctx = next(p for p in plan if p[4] == "CORRECT")
        print(f"test write (web-verified row): {ck}  DOI -> {doi_url(new)}")
        sc, err = patch_doi(tok, pid, new)
        print(f"  HTTP {sc} " + ("OK — token has write access" if sc < 300 else "FAIL " + err))
        return

    apply = "--apply" in sys.argv
    written = skipped = failed = 0
    skips = []
    print(f"plan: {len(plan)} candidate DOI writes "
          f"({sum(1 for p in plan if p[4]=='CORRECT')} correct, "
          f"{sum(1 for p in plan if p[4]=='BACKFILL')} backfill)\n")
    for pid, ck, old, new, kind, ctx in plan:
        if kind == "BACKFILL":
            ok, why = verify_backfill(new, ctx)
            time.sleep(0.1)
            if not ok:
                skipped += 1
                skips.append((ck, new, why))
                continue
        tag = "write" if apply else "WOULD-write"
        if apply:
            sc, err = patch_doi(tok, pid, new)
            if sc < 300:
                written += 1
            else:
                failed += 1
                print(f"  FAIL {ck}: HTTP {sc} {err}")
                continue
            time.sleep(0.2)
        else:
            written += 1
        if written <= 14 or kind == "CORRECT":
            print(f"  [{tag}/{kind}] {ck:32s} -> {new}")

    print(f"\n{'APPLIED' if apply else 'DRY-RUN'}: {written} "
          f"{'written' if apply else 'would write'}, {skipped} skipped (topic mismatch), "
          f"{failed} failed")
    if skips:
        print(f"\nSkipped backfills (need manual DOI) — {len(skips)}:")
        for ck, d, why in skips[:30]:
            print(f"  {ck:34s} {why}")
    if not apply:
        print("\n(dry-run) re-run with --apply to write verified DOIs to Notion.")


if __name__ == "__main__":
    main()
