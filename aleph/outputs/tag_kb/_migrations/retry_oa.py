#!/usr/bin/env python3
"""Retry PDF downloads for dl-fail entries in download_log.json.

Tries ALL Unpaywall oa_locations (not just best_oa_location), parses HTML
landing pages for embedded PDF links, and applies publisher-specific URL
patterns (PMC, MDPI, Frontiers, bioRxiv). Only saves if bytes start with %PDF.

  conda activate ffn_sim
  python retry_oa.py           # dry-run: print what would be tried
  python retry_oa.py --dl      # actually download
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path
from urllib.parse import urljoin, urlparse

import html.parser

import requests

HERE = Path(__file__).parent
LOG = HERE / "download_log.json"
OUT = Path("/Users/sw1/ffn_cellsim/aleph/references/downloaded")
EMAIL = "sungwook999@gmail.com"
UPW = "https://api.unpaywall.org/v2/"

BROWSER_UA = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/124.0.0.0 Safari/537.36"
)

S = requests.Session()
S.headers["User-Agent"] = BROWSER_UA


def norm_doi(s: str) -> str:
    return re.sub(r"^https?://(dx\.)?doi\.org/", "", (s or "").strip().lower()).rstrip(".)")


def unpaywall_all_locations(doi: str) -> list[dict]:
    """Return list of oa_location dicts from Unpaywall (all, not just best)."""
    try:
        r = S.get(f"{UPW}{doi}", params={"email": EMAIL}, timeout=25)
        if r.status_code != 200:
            return []
        j = r.json()
        if not j.get("is_oa"):
            return []
        locs = j.get("oa_locations") or []
        # also include best_oa_location if somehow missing from list
        best = j.get("best_oa_location")
        if best and best not in locs:
            locs = [best] + locs
        return locs
    except Exception:
        return []


def candidate_urls(loc: dict) -> list[str]:
    """Return candidate PDF URLs from a single oa_location dict."""
    urls = []
    if loc.get("url_for_pdf"):
        urls.append(loc["url_for_pdf"])
    if loc.get("url") and loc["url"] not in urls:
        urls.append(loc["url"])
    return urls


def pmc_pdf_url(url: str) -> str | None:
    """Convert PMC article page to direct PDF URL."""
    m = re.search(r"PMC(\d+)", url)
    if m:
        return f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{m.group(1)}/pdf/"
    return None


def publisher_pdf_candidates(url: str) -> list[str]:
    """Generate publisher-specific PDF URL variants from a landing-page URL."""
    extras: list[str] = []
    p = urlparse(url)
    host = p.netloc.lower()
    path = p.path

    # PMC
    if "ncbi.nlm.nih.gov/pmc" in url:
        pmc = pmc_pdf_url(url)
        if pmc:
            extras.append(pmc)

    # bioRxiv / medRxiv
    if any(x in host for x in ("biorxiv.org", "medrxiv.org")):
        base = url.split("v")[0] if "/v" in url else url
        for suffix in (".full.pdf", ".pdf"):
            candidate = re.sub(r"(v\d+)?$", "", url.rstrip("/")) + suffix
            if candidate not in extras:
                extras.append(candidate)

    # MDPI
    if "mdpi.com" in host and not path.endswith("/pdf"):
        extras.append(url.rstrip("/") + "/pdf")

    # Frontiers
    if "frontiersin.org" in host and not path.endswith("/pdf"):
        extras.append(url.rstrip("/") + "/pdf")

    # Europe PMC
    if "europepmc.org" in host:
        m = re.search(r"PMC(\d+)", url)
        if m:
            extras.append(
                f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{m.group(1)}/pdf/"
            )

    return extras


class _PDFLinkParser(html.parser.HTMLParser):
    """Minimal HTML parser that collects PDF-related link/meta/a hrefs."""

    def __init__(self) -> None:
        super().__init__()
        self.citation_pdf_url: str | None = None
        self.link_pdf_url: str | None = None
        self.first_a_pdf: str | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        d = dict(attrs)
        if tag == "meta":
            name = (d.get("name") or "").lower()
            if name == "citation_pdf_url" and d.get("content"):
                self.citation_pdf_url = self.citation_pdf_url or d["content"]
        elif tag == "link":
            typ = (d.get("type") or "").lower()
            rel = (d.get("rel") or "").lower()
            if ("application/pdf" in typ or "pdf" in rel) and d.get("href"):
                self.link_pdf_url = self.link_pdf_url or d["href"]
        elif tag == "a":
            href = d.get("href") or ""
            if href.lower().endswith(".pdf") and not self.first_a_pdf:
                self.first_a_pdf = href


def parse_html_for_pdf(raw: bytes, base_url: str) -> str | None:
    """Try to extract a direct PDF URL from an HTML landing page."""
    try:
        text = raw.decode("utf-8", errors="replace")
    except Exception:
        return None
    parser = _PDFLinkParser()
    try:
        parser.feed(text)
    except Exception:
        pass
    for candidate in (parser.citation_pdf_url, parser.link_pdf_url, parser.first_a_pdf):
        if candidate:
            return urljoin(base_url, candidate)
    return None


def try_fetch_pdf(url: str) -> bytes | None:
    """Fetch URL; if PDF bytes, return them. If HTML, try to find PDF link inside."""
    try:
        r = S.get(url, timeout=60, allow_redirects=True)
        if r.status_code != 200:
            return None
        ct = r.headers.get("content-type", "")
        if r.content[:4] == b"%PDF":
            return r.content
        # HTML landing page — parse for PDF link
        if "html" in ct or r.content[:5] == b"<!DOC" or r.content[:9] == b"<!doctype":
            pdf_url = parse_html_for_pdf(r.content, r.url)
            if pdf_url and pdf_url != url:
                r2 = S.get(pdf_url, timeout=60, allow_redirects=True)
                if r2.status_code == 200 and r2.content[:4] == b"%PDF":
                    return r2.content
        return None
    except Exception:
        return None


def main() -> None:
    do_dl = "--dl" in sys.argv

    log = json.loads(LOG.read_text()) if LOG.exists() else {}
    fails = {k: v for k, v in log.items() if v.get("status") == "dl-fail"}
    print(f"dl-fail entries: {len(fails)}")

    if not do_dl:
        print("(dry-run mode — pass --dl to download)\n")

    OUT.mkdir(parents=True, exist_ok=True)
    n_recovered = 0

    for ck, entry in sorted(fails.items()):
        doi = norm_doi(entry.get("doi", ""))
        if not doi:
            print(f"  SKIP {ck}: no doi")
            continue

        print(f"  {ck} ({doi}) ...", end=" ", flush=True)

        locs = unpaywall_all_locations(doi)
        time.sleep(0.3)

        if not locs:
            print("no OA locs")
            continue

        # Collect all candidate URLs (Unpaywall + publisher extras)
        all_urls: list[str] = []
        for loc in locs:
            for u in candidate_urls(loc):
                if u not in all_urls:
                    all_urls.append(u)
            for u in publisher_pdf_candidates(candidate_urls(loc)[0] if candidate_urls(loc) else ""):
                if u and u not in all_urls:
                    all_urls.append(u)

        if not all_urls:
            print("no candidate URLs")
            continue

        if not do_dl:
            print(f"{len(all_urls)} candidate URLs")
            continue

        # Try each candidate
        pdf_bytes = None
        used_url = None
        for u in all_urls:
            pdf_bytes = try_fetch_pdf(u)
            time.sleep(0.3)
            if pdf_bytes:
                used_url = u
                break

        if pdf_bytes:
            dest = OUT / f"{ck}.pdf"
            dest.write_bytes(pdf_bytes)
            log[ck] = {
                "doi": doi,
                "status": "downloaded",
                "bytes": len(pdf_bytes),
                "url": used_url,
            }
            n_recovered += 1
            print(f"OK ({len(pdf_bytes)//1024} KB)")
        else:
            # Keep as dl-fail but record attempt
            log[ck]["tried_urls"] = len(all_urls)
            print("fail")

    if do_dl:
        LOG.write_text(json.dumps(log, indent=1))
        total_dl = sum(1 for v in log.values() if v.get("status") == "downloaded")
        still_fail = sum(1 for v in log.values() if v.get("status") == "dl-fail")
        print(f"\nRecovered: {n_recovered} | total downloaded: {total_dl} | still dl-fail: {still_fail}")
        print(f"log -> {LOG}")
        print("\nnext steps:")
        print("  cd /Users/sw1/ffn_cellsim/aleph/outputs/tag_kb")
        print("  python references_ingest.py")
        print("  python notion_to_duckdb.py")


if __name__ == "__main__":
    main()
