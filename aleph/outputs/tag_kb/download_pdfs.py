#!/usr/bin/env python3
"""Download open-access PDFs for SourceEvidence DOIs via Unpaywall.

Only fetches legally-free OA copies (Unpaywall's best_oa_location). Closed/
paywalled DOIs are logged as "needs gbook" (KAIST institutional full-text, a
separate manual/SSH track) — never scraped. Downloaded PDFs land in
references/downloaded/<citation_key>.pdf so references_ingest.py picks them up
into paper_chunks on the next run.

  python download_pdfs.py --coverage      # just measure OA availability (no DL)
  python download_pdfs.py --download 40    # download up to 40 OA PDFs
  python download_pdfs.py --download 0     # download all available OA PDFs
"""
from __future__ import annotations

import json
import re
import sys
import time

import duckdb
import requests

HERE = __import__("pathlib").Path(__file__).parent
DB = HERE / "kb.duckdb"
OUT = __import__("pathlib").Path("/Users/sw1/ffn_cellsim/aleph/references/downloaded")
LOG = HERE / "download_log.json"
EMAIL = "sungwook999@gmail.com"
UPW = "https://api.unpaywall.org/v2/"
S = requests.Session()
S.headers["User-Agent"] = f"ffn_cellsim-oa-fetch/1.0 (mailto:{EMAIL})"


def norm_doi(s):
    return re.sub(r"^https?://(dx\.)?doi\.org/", "", (s or "").strip().lower()).rstrip(".)")


def oa_pdf_url(doi):
    try:
        r = S.get(f"{UPW}{doi}", params={"email": EMAIL}, timeout=25)
        if r.status_code != 200:
            return None, f"unpaywall {r.status_code}"
        j = r.json()
        if not j.get("is_oa"):
            return None, "closed (needs gbook)"
        loc = j.get("best_oa_location") or {}
        url = loc.get("url_for_pdf") or loc.get("url")
        return (url, "oa") if url else (None, "oa-no-pdf-url")
    except Exception as e:
        return None, f"err {type(e).__name__}"


def main():
    download = "--download" in sys.argv
    limit = int(sys.argv[sys.argv.index("--download") + 1]) if download else 0
    con = duckdb.connect(str(DB), read_only=True)
    rows = con.execute("SELECT citation_key, doi FROM source_evidence "
                       "WHERE doi IS NOT NULL ORDER BY citation_key").fetchall()
    con.close()

    have = {p.stem for p in OUT.glob("*.pdf")} if OUT.exists() else set()
    log = json.loads(LOG.read_text()) if LOG.exists() else {}
    n_oa = n_closed = n_dl = n_fail = 0
    for ck, doi in rows:
        d = norm_doi(doi)
        if ck in have or log.get(ck, {}).get("status") == "downloaded":
            continue
        url, why = oa_pdf_url(d)
        time.sleep(0.1)
        if not url:
            log[ck] = {"doi": d, "status": "closed" if "closed" in why else why}
            n_closed += "closed" in why
            continue
        n_oa += 1
        if not download:
            log[ck] = {"doi": d, "status": "oa-available", "url": url}
            continue
        if limit and n_dl >= limit:
            log[ck] = {"doi": d, "status": "oa-available", "url": url}
            continue
        try:
            OUT.mkdir(parents=True, exist_ok=True)
            pr = S.get(url, timeout=60)
            if pr.status_code == 200 and pr.content[:4] == b"%PDF":
                (OUT / f"{ck}.pdf").write_bytes(pr.content)
                log[ck] = {"doi": d, "status": "downloaded", "bytes": len(pr.content)}
                n_dl += 1
                print(f"  + {ck} ({len(pr.content)//1024} KB)", flush=True)
            else:
                log[ck] = {"doi": d, "status": "dl-fail", "http": pr.status_code}
                n_fail += 1
        except Exception as e:
            log[ck] = {"doi": d, "status": "dl-fail", "err": type(e).__name__}
            n_fail += 1
        time.sleep(0.3)

    LOG.write_text(json.dumps(log, indent=1))
    tot = len(rows)
    oa_total = sum(1 for v in log.values() if v.get("status") in ("oa-available", "downloaded"))
    closed = sum(1 for v in log.values() if v.get("status") == "closed")
    dled = sum(1 for v in log.values() if v.get("status") == "downloaded")
    print(f"\n{tot} DOIs | OA available: {oa_total} | closed(needs gbook): {closed} | "
          f"downloaded: {dled} | this run dl {n_dl} fail {n_fail}")
    print(f"log -> {LOG}")
    if download:
        print("next: python references_ingest.py  (ingest new PDFs into paper_chunks)")


if __name__ == "__main__":
    main()
