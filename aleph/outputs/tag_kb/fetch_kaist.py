#!/usr/bin/env python3
"""Fetch paywalled PDFs via KAIST institutional access (gbook SSH).

Strategy:
  1. PDFs already present locally → copy to downloaded/.
  2. Elsevier/ScienceDirect → flag "manual-browser" (bot-blocks institutional curl).
  3. Remaining publishers → batch-fetch on gbook via known PDF URL patterns,
     SCP results back, validate %PDF header.

Usage:
    cd aleph/outputs/tag_kb
    python fetch_kaist.py [--dry-run]
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
ROOT = Path("/Users/sw1/ffn_cellsim/ffn_sim")
OUT = ROOT / "references" / "downloaded"
LOG = HERE / "download_log.json"
REF_ROOT = ROOT / "references"
DRY_RUN = "--dry-run" in sys.argv

UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)

# ── 1. PDFs already on disk locally (not yet in downloaded/) ─────────────────

LOCAL_MAP: dict[str, Path] = {
    "Fang2016_PhysRevE": REF_ROOT / "PhysRevE.93.042404.pdf",
    "Kim2012_IntegrBiol": REF_ROOT / "c2ib20159c.pdf",
    "Dhandapani2023_AdvHealthcareMaterials": REF_ROOT / (
        "Adv Healthcare Materials - 2023 - Dhandapani - In Vitro 3D "
        "Spheroid Model Preserves Tumor Microenvironment of Hot and Cold.pdf"
    ),
    "Bustamante2021_Biofabrication": REF_ROOT / "Bustamante_2021_Biofabrication_13_035010.pdf",
}

# ── 2. Elsevier/ScienceDirect — flag manual-browser ──────────────────────────

ELSEVIER_KEYS: set[str] = {
    "Bastatas2012_BBA",          # 10.1016/j.bbagen.2012.02.006
    "Darling2008_JBiomech",      # 10.1016/j.jbiomech.2007.06.019
    "Dimova2014_ACIS",           # 10.1016/j.cis.2014.03.003
    "DizMunoz2013_TCB",          # 10.1016/j.tcb.2012.09.006
    "Feng2025_BiochemicalEngineeringJo",  # 10.1016/j.bej.2024.109567
    "Goldmann2002",              # 10.1006/cbir.2002.0900
    "Guilak2000_BBRC",           # 10.1006/bbrc.2000.2360
    "K562_BBRC2019",             # 10.1016/j.bbrc.2019.06.054
    "Li2008_BBRC",               # 10.1016/j.bbrc.2008.07.078
    "Omidvar2014_JBiomech",      # 10.1016/j.jbiomech.2014.08.002
    "PontesGauthier2017_SCDB",   # 10.1016/j.semcdb.2017.08.030
    "Salbreux2012_TCB",          # 10.1016/j.tcb.2012.07.001
}

# ── 3. Publisher-pattern PDF URLs (to fetch on gbook) ────────────────────────

PDF_URLS: dict[str, str] = {
    # APS
    "Marchetti2013_RMP": (
        "https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.86.995"
    ),
    # Science / AAAS
    "Bell1978_Science": (
        "https://www.science.org/doi/pdf/10.1126/science.347575"
    ),
    "Bustamante1994_Science": (
        "https://www.science.org/doi/pdf/10.1126/science.8079175"
    ),
    "ChanOdde2008_Science": (
        "https://www.science.org/doi/pdf/10.1126/science.1163595"
    ),
    "DenaisRaab2016_Science": (
        "https://www.science.org/doi/pdf/10.1126/science.aad7297"
    ),
    "Discher2005_Science": (
        "https://www.science.org/doi/pdf/10.1126/science.1116995"
    ),
    "Swift2013_Science": (
        "https://www.science.org/doi/pdf/10.1126/science.1240104"
    ),
    # Nature Publishing Group — article_id = doi suffix
    "Bertet2004_Nature": "https://www.nature.com/articles/nature02590.pdf",
    "Cavey2008_Cell":    "https://www.nature.com/articles/nature06953.pdf",
    "Kanchanawong2010_Nature": "https://www.nature.com/articles/nature09621.pdf",
    "Plodinec2012_NatNanotechnol": "https://www.nature.com/articles/nnano.2012.167.pdf",
    "Prost2015_NatPhys": "https://www.nature.com/articles/nphys3224.pdf",
    "Reffay2014_NCB":    "https://www.nature.com/articles/ncb2917.pdf",
    "SerraPicamal2012_NatPhys": "https://www.nature.com/articles/nphys2355.pdf",
    "Stewart2011_Nature": "https://www.nature.com/articles/nature09642.pdf",
    "Trepat2009_NatPhys": "https://www.nature.com/articles/nphys1269.pdf",
    "Yamada2019_NRMCB":  "https://www.nature.com/articles/s41580-019-0172-9.pdf",
    "Yonemura2010_NCB":  "https://www.nature.com/articles/ncb2055.pdf",
    # PNAS
    "Douezan2011_PNAS":  "https://www.pnas.org/doi/pdf/10.1073/pnas.1018057108",
    "Pajerowski2007_PNAS": "https://www.pnas.org/doi/pdf/10.1073/pnas.0702576104",
    # Wiley
    "MotteKaufman2013_Biopolymers": (
        "https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/bip.22133"
    ),
    "Rodenhizer2018_AdvHealthcareMaterials": (
        "https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/adhm.201701174"
    ),
    # ACS Publications
    "MarkoSiggia1995_Macromolecules": (
        "https://pubs.acs.org/doi/pdf/10.1021/ma00130a008"
    ),
    # RSC
    "Sharma2014_IntegrBiol": (
        "https://pubs.rsc.org/en/content/articlepdf/2014/ib/c3ib40246k"
    ),
    # MBoC (ASCB)
    "Stephens2017_MBoC": (
        "https://www.molbiolcell.org/doi/pdf/10.1091/mbc.e16-09-0653"
    ),
    # Springer Link
    "GalWeihs2012_CBB": (
        "https://link.springer.com/content/pdf/10.1007/s12013-012-9356-z.pdf"
    ),
    # SPIE (resolve via DOI redirect + PDF link extraction)
    "Bredfeldt2014_JBiomedOpt": "__resolve__:10.1117/1.jbo.19.1.016007",
    # ASME
    "Brodland2002": "__resolve__:10.1115/1.1449491",
    # Company of Biologists (Development)
    "Foty1996_PNAS": "__resolve__:10.1242/dev.122.5.1611",
    # Project MUSE
    "Evans2007_RBP": "__resolve__:10.1353/scu.2007.0015",
}


def _gbook(cmd: str, timeout: int = 60) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["ssh", "gbook", cmd],
        capture_output=True,
        timeout=timeout + 10,
    )


def resolve_doi_to_pdf(doi: str) -> str | None:
    """On gbook: follow DOI redirect, look for PDF URL in HTML."""
    script = f"""
python3 - <<'PYEOF'
import urllib.request, re, sys

UA = '{UA}'
doi = '{doi}'

req = urllib.request.Request(
    f'https://doi.org/{{doi}}',
    headers={{'User-Agent': UA, 'Accept': 'text/html,application/pdf'}},
)
try:
    with urllib.request.urlopen(req, timeout=20) as resp:
        url = resp.url
        ct = resp.headers.get('Content-Type', '')
        if 'pdf' in ct.lower():
            print(url)
            sys.exit(0)
        html = resp.read(65536).decode('utf-8', errors='replace')
except Exception as e:
    print(f'ERR {{e}}', file=sys.stderr)
    sys.exit(1)

# Look for PDF link patterns in HTML
patterns = [
    r'<meta[^>]+name=["\']citation_pdf_url["\'][^>]+content=["\']([^"\']+)["\']',
    r'content=["\']([^"\']+\.pdf(?:\?[^"\']*)?)["\'][^>]+name=["\']citation_pdf_url["\']',
    r'href=["\']([^"\']*\.pdf(?:\?[^"\']*)?)["\']',
    r'href=["\']([^"\']+/pdf(?:/|\?|$)[^"\']*)["\']',
    r'window\.location\s*=\s*["\']([^"\']+\.pdf[^"\']*)["\']',
]
for pat in patterns:
    m = re.search(pat, html, re.IGNORECASE)
    if m:
        u = m.group(1)
        if not u.startswith('http'):
            from urllib.parse import urljoin
            u = urljoin(url, u)
        print(u)
        sys.exit(0)

print(f'NOURL {{url}}', file=sys.stderr)
sys.exit(1)
PYEOF
"""
    r = _gbook(script, timeout=40)
    if r.returncode == 0 and r.stdout.strip():
        return r.stdout.decode(errors="replace").strip()
    return None


def fetch_pdf_on_gbook(key: str, url: str, remote_dir: str) -> bool:
    """Fetch URL on gbook, save to remote_dir/key.pdf. Return True on success."""
    cmd = (
        f"curl -sL -A '{UA}' --max-time 55 --retry 1 "
        f"-H 'Accept: application/pdf,*/*' "
        f"-o '{remote_dir}/{key}.pdf' '{url}' && "
        f"head -c 4 '{remote_dir}/{key}.pdf' 2>/dev/null | grep -q '%PDF' && "
        f"echo OK || (rm -f '{remote_dir}/{key}.pdf'; echo FAIL)"
    )
    r = _gbook(cmd, timeout=80)
    out = (r.stdout + r.stderr).decode(errors="replace").strip()
    return "OK" in out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log: dict = json.loads(LOG.read_text())

    closed = {
        ck: info
        for ck, info in log.items()
        if info.get("status") == "closed"
    }
    print(f"{len(closed)} closed entries to process")

    n_local = n_elsevier = n_fetched = n_fail = n_resolve_fail = 0
    remote_dir = "/tmp/kaist_pdfs"

    if not DRY_RUN:
        _gbook(f"mkdir -p {remote_dir}")

    # ── Step 1: handle local copies ──────────────────────────────────────────
    for ck, src in LOCAL_MAP.items():
        if ck not in closed:
            continue
        dest = OUT / f"{ck}.pdf"
        if dest.exists():
            print(f"  [skip local-exists] {ck}")
            continue
        if not src.exists():
            print(f"  [WARN local-missing] {ck}: {src}")
            continue
        if DRY_RUN:
            print(f"  [dry local] {ck}")
            continue
        shutil.copy2(src, dest)
        sz = dest.stat().st_size
        log[ck] = {"doi": log[ck]["doi"], "status": "downloaded", "bytes": sz}
        n_local += 1
        print(f"  + [local] {ck} ({sz // 1024} KB)")

    # ── Step 2: flag Elsevier ─────────────────────────────────────────────────
    for ck in ELSEVIER_KEYS:
        if ck not in closed:
            continue
        if DRY_RUN:
            print(f"  [dry elsevier] {ck}")
            continue
        log[ck]["status"] = "manual-browser"
        n_elsevier += 1
        print(f"  [manual-browser] {ck}")

    # ── Step 3: gbook batch fetch ─────────────────────────────────────────────
    to_fetch = {
        ck: url
        for ck, url in PDF_URLS.items()
        if ck in closed and ck not in LOCAL_MAP
    }

    for ck, url in to_fetch.items():
        dest = OUT / f"{ck}.pdf"
        if dest.exists():
            print(f"  [skip exists] {ck}")
            log[ck] = {"doi": log[ck]["doi"], "status": "downloaded",
                       "bytes": dest.stat().st_size}
            continue

        # Resolve DOIs for publishers with non-obvious PDF URLs
        if url.startswith("__resolve__:"):
            doi = url[len("__resolve__:"):]
            print(f"  [resolving] {ck} ({doi})", flush=True)
            if not DRY_RUN:
                url = resolve_doi_to_pdf(doi)
                if not url:
                    print(f"    -> resolve failed")
                    log[ck]["status"] = "gbook-resolve-fail"
                    n_resolve_fail += 1
                    continue
                print(f"    -> {url[:80]}")
            else:
                print(f"  [dry resolve] {ck}")
                continue

        if DRY_RUN:
            print(f"  [dry fetch] {ck}  {url[:80]}")
            continue

        print(f"  [fetch] {ck}", flush=True)
        ok = fetch_pdf_on_gbook(ck, url, remote_dir)
        if not ok:
            print(f"    -> FAIL")
            log[ck]["status"] = "gbook-fail"
            n_fail += 1
            continue

        # SCP back
        scp = subprocess.run(
            ["scp", f"gbook:{remote_dir}/{ck}.pdf", str(dest)],
            capture_output=True, timeout=120,
        )
        if scp.returncode != 0 or not dest.exists():
            print(f"    -> SCP fail")
            log[ck]["status"] = "gbook-scp-fail"
            n_fail += 1
            continue

        # Final %PDF validation
        if dest.read_bytes()[:4] != b"%PDF":
            print(f"    -> not PDF, removing")
            dest.unlink(missing_ok=True)
            log[ck]["status"] = "gbook-not-pdf"
            n_fail += 1
            continue

        sz = dest.stat().st_size
        log[ck] = {"doi": log[ck]["doi"], "status": "downloaded", "bytes": sz}
        n_fetched += 1
        print(f"    -> OK ({sz // 1024} KB)")

    # ── Write log ─────────────────────────────────────────────────────────────
    if not DRY_RUN:
        LOG.write_text(json.dumps(log, indent=1))
        print(f"\nlog -> {LOG}")

    total_dl = sum(1 for v in log.values() if v.get("status") == "downloaded")
    still_closed = sum(1 for v in log.values() if v.get("status") == "closed")
    print(
        f"\nThis run — local: {n_local}  elsevier→manual: {n_elsevier}  "
        f"gbook-fetched: {n_fetched}  fail: {n_fail}  resolve-fail: {n_resolve_fail}"
    )
    print(f"Totals — downloaded: {total_dl}  still-closed: {still_closed}")
    if not DRY_RUN:
        print("\nnext: python references_ingest.py && python notion_to_duckdb.py")


if __name__ == "__main__":
    main()
