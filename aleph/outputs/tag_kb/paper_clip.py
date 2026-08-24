#!/usr/bin/env python3
"""Discover candidate papers across open bibliographic sources and fetch their OA full text.

This is the DISCOVERY half of the corpus pipeline and it is deliberately a different concept
from its two neighbours:

* ``download_pdfs.py`` — BACKFILL.  Takes DOIs that are *already* SourceEvidence rows and fetches
  their OA copy.  It can never bring in a paper the KB has not already heard of.
* ``paper_clip.py`` (this file) — DISCOVERY.  Searches the open bibliographic world for papers the
  corpus does NOT have, in queries anchored to the engine's own component/connector vocabulary.
* ``references_ingest.py`` — INGEST.  Turns whatever PDFs are on disk into ``paper_refs`` +
  ``paper_chunks`` (BM25).  It reads ``references/<subdir>/*.pdf`` and takes the subdirectory name
  as the ``source`` column, which is why this script downloads into ``references/<tag>/``.

WHY THE QUERY SET IS NOT FREE TEXT.  Every query in ``clip_queries.yaml`` names the component or
connector it is trying to constrain, the parameter axis it is after, and — where one exists — the
PI-GAP card that recorded the previous search failing.  That does two things a keyword list cannot.
It makes each arriving PDF carry its *reason for arriving* (the seed of the paper -> connector edge,
recorded at clip time rather than retrofitted), and it lets the run target the four measured KINDS of
absence separately: a value that does not exist in the literature at all needs a different query than
one that exists only for another cell line, one that is paywalled, or one where a number exists but
is the WRONG PHYSICAL QUANTITY (PI-GAP card N1: the literature's myosin number is a step SIZE, not a
stall FORCE — the failure mode a relevance-ranked search cannot see).

OA ONLY, NEVER SCRAPED.  Only locations that Unpaywall / Europe PMC / arXiv / OpenAlex declare
openly available are fetched.  A closed DOI is recorded with ``status: no_oa`` and routed to the
existing KAIST institutional track (``fetch_kaist.py``); no paywall is ever circumvented.

DRY RUN IS THE DEFAULT, matching ``harvest_ops.py``: the first pass writes a candidate manifest for
review and downloads nothing.  ``--apply`` is the second, deliberate step.

Integrity Gate (the non-physics analogue of the Sanity Gate; this module computes no physics):
    * provenance: every hit records which provider returned it and which query id asked for it.
      A record whose provenance cannot be reconstructed is not written.
    * no invented metadata: title/DOI/year come from the provider response verbatim.  Nothing is
      inferred, completed, or guessed — the 2026-06-02 audit found three hallucinated sources in the
      KB and the only structural defence is that this layer never authors a field.
    * dedup: normalized DOI, then arXiv id (version-stripped), then normalized title, checked
      against ``paper_refs`` + ``source_evidence`` in ``kb.duckdb`` AND against previous clip runs.
      Byte-level duplicates are caught again downstream by ``references_ingest``'s sha1 pass.
    * content check: a fetched body must begin with the ``%PDF`` magic bytes.  Publisher "PDF" URLs
      routinely return an HTML landing page, and an HTML file ingested as a paper becomes silent
      garbage in the BM25 index.
    * rate limits: taken from each provider's published guidance and recorded in the manifest with
      their basis, never tuned to go faster.
    * reversibility: this script only ADDS files under ``references/<tag>/`` and appends to the
      analysis manifest.  It never deletes or rewrites an existing PDF.

Usage:
    python paper_clip.py                                  # dry run over clip_queries.yaml
    python paper_clip.py --id nmii-stall-force            # dry run, one query
    python paper_clip.py --query "cortical tension MCF7"  # dry run, ad-hoc
    python paper_clip.py --apply                          # fetch the OA subset
    python paper_clip.py --check                          # clipped but not yet ingested
"""

from __future__ import annotations

import argparse
import json
import logging
import pathlib
import re
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass, field
from datetime import date
from html import unescape
from typing import Any, Callable, Iterable, Sequence

import requests
import yaml

LOG = logging.getLogger("paper_clip")

HERE = pathlib.Path(__file__).resolve().parent
#: ``.../ffn_cellsim/ffn_sim`` — this file lives at ``aleph/outputs/tag_kb/``.
PKG_ROOT = HERE.parents[1]
REPO_ROOT = HERE.parents[2]
#: The corpus ``references_ingest.py`` reads.  It MUST be the package-level directory, not the
#: repo root: the ingest globs ``aleph/references/*/*.pdf`` and takes the subdirectory name as the
#: ``source`` column, so a corpus written one level up is invisible to it.
REF_DIR = PKG_ROOT / "references"
ANALYSIS_MANIFEST = REF_DIR / "analysis" / "_manifest.json"
DB_PATH = HERE / "kb.duckdb"
DEFAULT_QUERY_FILE = HERE / "clip_queries.yaml"

#: Contact address sent to every provider's polite pool, matching ``download_pdfs.py``.
EMAIL = "sungwook999@gmail.com"
USER_AGENT = f"ffn_cellsim-paper-clip/1.0 (mailto:{EMAIL})"

#: Bytes that a real PDF starts with.  Not a tunable — it is the format's magic number.
PDF_MAGIC = b"%PDF"

#: Largest body we will pull for a single paper.  A courtesy ceiling, not a physical constant;
#: recorded here so it is declared rather than hidden inside the fetch loop.
MAX_PDF_BYTES = 80 * 1024 * 1024

#: Statuses worth retrying: a throttle, and the transient server errors behind a load balancer.
RETRYABLE_STATUS = frozenset({429, 500, 502, 503, 504})
RETRY_ATTEMPTS = 4
#: First backoff step; doubles per attempt, and a server-sent ``Retry-After`` overrides it.
RETRY_BASE_DELAY_S = 2.0
#: Longest ``Retry-After`` we will ever honour.  A server may legitimately ask for far longer —
#: OpenAlex answered ``Retry-After: 35478`` (9.9 h) once its daily credit ran out — and obeying that
#: literally parks the process for the rest of the day with no output.  Past this ceiling the
#: provider is marked EXHAUSTED for the run instead, which is a state the caller can reason about.
RETRY_MAX_DELAY_S = 60.0
#: Requests to leave unspent when a provider reports a remaining quota, so a run does not consume
#: the last of a shared daily budget that another lane may need.
QUOTA_RESERVE = 25

DOI_PREFIX_RE = re.compile(r"^https?://(dx\.)?doi\.org/", re.IGNORECASE)
ARXIV_ID_RE = re.compile(r"(\d{4}\.\d{4,5}|[a-z-]+(\.[A-Z]{2})?/\d{7})(v\d+)?", re.IGNORECASE)
NON_ALNUM_RE = re.compile(r"[^a-z0-9]+")


@dataclass(frozen=True, slots=True)
class Provider:
    """One open bibliographic source, with the published basis for its rate limit.

    Attributes:
        name: Short identifier used on the command line and in the manifest.
        min_interval_s: Minimum seconds between requests to this provider.
        rate_basis: Where ``min_interval_s`` comes from.  Recorded in the manifest so a reviewer
            can see whether a delay is a provider requirement or our own courtesy choice.
        covers: One line on what this provider adds that the others do not.
    """

    name: str
    min_interval_s: float
    rate_basis: str
    covers: str


PROVIDERS: dict[str, Provider] = {
    "arxiv": Provider(
        name="arxiv",
        min_interval_s=3.0,
        rate_basis="arXiv API User Manual: 'no more than one request every three seconds'",
        covers="physics/q-bio preprints, incl. soft-matter and biophysics modelling not in PubMed",
    ),
    "europepmc": Provider(
        name="europepmc",
        min_interval_s=0.34,
        rate_basis="no published hard limit; courtesy ~3 req/s",
        covers="PubMed + PMC + bioRxiv/medRxiv preprints (SRC:PPR), with OA full-text links",
    ),
    "openalex": Provider(
        name="openalex",
        min_interval_s=0.1,
        rate_basis="OpenAlex docs: 10 requests/second, 100k/day in the polite pool (mailto)",
        covers="broadest metadata graph + best_oa_location; catches venues the others miss",
    ),
    "crossref": Provider(
        name="crossref",
        min_interval_s=0.2,
        rate_basis="Crossref polite pool (mailto); courtesy 5 req/s",
        covers="publisher-of-record DOI metadata; the authority when OpenAlex and EPMC disagree",
    ),
}

#: Resolver used to turn a DOI into an OA PDF location when the search provider gave none.
UNPAYWALL_INTERVAL_S = 0.1
UNPAYWALL_RATE_BASIS = "Unpaywall docs: 100,000 calls/day; same 0.1 s spacing as download_pdfs.py"


@dataclass
class Hit:
    """One candidate paper as returned by one provider, before dedup or fetching.

    Every field is copied verbatim from the provider response.  Nothing here is inferred; an
    unavailable field stays empty rather than being completed from another source.
    """

    provider: str
    query_id: str
    title: str = ""
    doi: str = ""
    arxiv_id: str = ""
    year: str = ""
    venue: str = ""
    authors: str = ""
    abstract: str = ""
    pdf_url: str = ""
    landing_url: str = ""
    is_oa: bool | None = None
    provider_id: str = ""
    status: str = "candidate"
    reason: str = ""
    local_path: str = ""
    work_type: str = ""
    cited_by_count: int | None = None
    fwci: float | None = None
    citation_percentile: float | None = None
    journal_2yr_mean_citedness: float | None = None
    screen: str = ""
    relevance: str = ""
    relevance_reason: str = ""
    relevance_matched: str = ""


    @property
    def dedup_keys(self) -> tuple[str, ...]:
        """The identity keys this hit may collide on, in descending authority."""
        keys: list[str] = []
        if self.doi:
            keys.append(f"doi:{norm_doi(self.doi)}")
        if self.arxiv_id:
            keys.append(f"arxiv:{norm_arxiv(self.arxiv_id)}")
        if self.title:
            keys.append(f"title:{norm_title(self.title)}")
        return tuple(keys)


@dataclass(frozen=True, slots=True)
class Query:
    """One clip query, anchored to the engine vocabulary rather than being free text.

    Attributes:
        qid: Stable identifier; also the provenance stamp carried by every hit it produces.
        all_terms: Phrases that must ALL appear.  Joined with AND where the provider supports it.
            These alone are the STRICT query, and they are what every provider runs first.
        any_terms: Alternatives to the LAST ``all_terms`` phrase, used only as a fallback when the
            strict query returns nothing.  They are never added as an extra required group: doing so
            narrows the query rather than widening it, which is how a first version of this file
            turned an 8-hit search into a 1-hit one.  Boolean providers (arXiv, Europe PMC) widen to
            ``AND(all[:-1]) AND OR(all[-1], *any)``; relevance providers (OpenAlex, Crossref) have
            no boolean syntax to widen with and fall back to the first phrase alone.
        component: Engine component this query is trying to constrain, if any.
        connector: Engine connector this query is trying to constrain, if any.
        connector_family: The ``ConnectorFamily`` member, when a connector is named.
        parameter_axis: The parameter this query is hunting.
        pi_gap_card: The PI-GAP evidence card recording the previous search, if there was one.
        absence: Which measured KIND of absence the previous pass hit — ``none`` (no value exists),
            ``proxy`` (other cell line only), ``paywalled``, ``wrong_quantity``, or empty.
        cell_state: The CellState row this is for, when the query is cell-line specific.
        note: Free text for the reviewer; never sent to a provider.
    """

    qid: str
    all_terms: tuple[str, ...] = ()
    any_terms: tuple[str, ...] = ()
    anchors: tuple[str, ...] = ()
    off_topic: tuple[str, ...] = ()
    component: str = ""
    connector: str = ""
    connector_family: str = ""
    parameter_axis: str = ""
    pi_gap_card: str = ""
    absence: str = ""
    cell_state: str = ""
    note: str = ""

    def provenance(self) -> dict[str, str]:
        """The subset of this query that should travel with every hit it produced."""
        return {
            k: v
            for k, v in {
                "query_id": self.qid,
                "component": self.component,
                "connector": self.connector,
                "connector_family": self.connector_family,
                "parameter_axis": self.parameter_axis,
                "pi_gap_card": self.pi_gap_card,
                "absence": self.absence,
                "cell_state": self.cell_state,
            }.items()
            if v
        }


# --------------------------------------------------------------------------------------
# normalisation
# --------------------------------------------------------------------------------------


def norm_doi(value: str | None) -> str:
    """Strip a DOI down to its bare lowercase form, matching ``download_pdfs.norm_doi``."""
    return DOI_PREFIX_RE.sub("", (value or "").strip().lower()).rstrip(".)")


def norm_arxiv(value: str | None) -> str:
    """Reduce an arXiv reference to its version-less id, so ``v1`` and ``v3`` collide."""
    if not value:
        return ""
    match = ARXIV_ID_RE.search(value)
    if not match:
        return ""
    return match.group(1).lower()


def norm_title(value: str | None) -> str:
    """Collapse a title to lowercase alphanumerics for duplicate detection."""
    return NON_ALNUM_RE.sub(" ", (value or "").lower()).strip()


def slugify(value: str, limit: int = 60) -> str:
    """Filename-safe slug, same shape as ``references_ingest.slug``."""
    out = NON_ALNUM_RE.sub("-", (value or "").lower()).strip("-")
    return out[:limit] or "untitled"


# --------------------------------------------------------------------------------------
# HTTP
# --------------------------------------------------------------------------------------


class Fetcher:
    """A requests session that honours each provider's published spacing."""

    def __init__(self, *, timeout: float = 30.0) -> None:
        self.session = requests.Session()
        self.session.headers["User-Agent"] = USER_AGENT
        self.timeout = timeout
        self._last_call: dict[str, float] = {}
        #: Providers that told us to come back later than :data:`RETRY_MAX_DELAY_S`.  Recorded rather
        #: than retried, so the caller can degrade honestly instead of the run stalling or, worse,
        #: treating an absent metric as a low one.
        self.exhausted: dict[str, str] = {}
        #: Last reported remaining quota per provider, when the provider publishes one.
        self.quota_remaining: dict[str, int] = {}

    def _note_quota(self, key: str, response: requests.Response) -> None:
        raw = response.headers.get("x-ratelimit-remaining", "")
        if raw.strip().lstrip("-").isdigit():
            remaining = int(raw)
            self.quota_remaining[key] = remaining
            if remaining <= QUOTA_RESERVE and key not in self.exhausted:
                reset = response.headers.get("x-ratelimit-reset", "?")
                self.exhausted[key] = f"quota reserve reached ({remaining} left, resets in {reset}s)"
                LOG.warning("%s: %s — no further calls this run", key, self.exhausted[key])

    def _wait(self, key: str, min_interval_s: float) -> None:
        last = self._last_call.get(key)
        if last is not None:
            gap = time.monotonic() - last
            if gap < min_interval_s:
                time.sleep(min_interval_s - gap)
        self._last_call[key] = time.monotonic()

    def get(self, key: str, min_interval_s: float, url: str, **kwargs: Any) -> requests.Response | None:
        """GET with per-provider spacing and backoff on a throttle or a transient server error.

        Retrying 429/5xx is not politeness theatre.  Without it a bulk run drops whole batches
        silently: a throttled OpenAlex enrichment request would leave 50 papers with no metrics, and
        they would then be screened out for "no fwci" — a rate limit reappearing as a quality
        verdict.  ``Retry-After`` is honoured when the server sends one.

        Returns:
            The response, or ``None`` once the retries are exhausted or the failure is permanent.
        """
        if key in self.exhausted:
            return None
        for attempt in range(RETRY_ATTEMPTS):
            self._wait(key, min_interval_s)
            try:
                response = self.session.get(url, timeout=self.timeout, **kwargs)
            except requests.RequestException as exc:
                LOG.warning("%s: request failed (%s)", key, type(exc).__name__)
                return None
            self._note_quota(key, response)
            if response.status_code == 200:
                return response
            if response.status_code in RETRYABLE_STATUS and attempt < RETRY_ATTEMPTS - 1:
                header = response.headers.get("Retry-After", "").strip()
                asked = float(header) if header.replace(".", "", 1).isdigit() else (
                    RETRY_BASE_DELAY_S * (2**attempt)
                )
                if asked > RETRY_MAX_DELAY_S:
                    self.exhausted[key] = (
                        f"HTTP {response.status_code}, Retry-After {asked:.0f}s "
                        f"exceeds the {RETRY_MAX_DELAY_S:.0f}s ceiling"
                    )
                    LOG.warning("%s: %s — treating as EXHAUSTED for this run",
                                key, self.exhausted[key])
                    return None
                LOG.info("%s: HTTP %s — backing off %.1fs (attempt %d/%d)",
                         key, response.status_code, asked, attempt + 1, RETRY_ATTEMPTS)
                time.sleep(asked)
                continue
            LOG.warning("%s: HTTP %s for %s", key, response.status_code, url)
            return None
        return None


# --------------------------------------------------------------------------------------
# providers
# --------------------------------------------------------------------------------------

_ATOM_NS = {"atom": "http://www.w3.org/2005/Atom", "arxiv": "http://arxiv.org/schemas/atom"}


def boolean_forms(query: Query, field: str = "") -> list[str]:
    """Build the strict query, then a genuinely WIDER fallback, for a boolean provider.

    Args:
        query: The query to render.
        field: Provider field prefix applied to each phrase (``all:`` for arXiv, empty for EPMC).

    Returns:
        One or two query strings, strict first.  The second is only present when ``any_terms``
        exist, and it relaxes the LAST required phrase into an alternation — that is the widening.
        Appending the alternation as an additional required group would narrow instead.
    """
    def phrase(term: str) -> str:
        return f'{field}"{term}"'

    if not query.all_terms:
        return [" OR ".join(phrase(t) for t in query.any_terms)] if query.any_terms else []

    forms = [" AND ".join(phrase(t) for t in query.all_terms)]
    if query.any_terms:
        head = [phrase(t) for t in query.all_terms[:-1]]
        alternation = "(" + " OR ".join(
            phrase(t) for t in (query.all_terms[-1], *query.any_terms)
        ) + ")"
        forms.append(" AND ".join([*head, alternation]))
    return forms


def relevance_forms(query: Query) -> list[str]:
    """Build the strict then fallback query for a relevance-ranked provider.

    OpenAlex and Crossref have no boolean operators here, so more words means a narrower result set.
    Widening therefore means using FEWER phrases, never adding ``any_terms``.
    """
    terms = list(query.all_terms) or list(query.any_terms)
    if not terms:
        return []
    forms = [" ".join(terms)]
    if len(terms) > 1:
        forms.append(terms[0])
    return forms


def search_arxiv(fetcher: Fetcher, query: Query, limit: int) -> list[Hit]:
    """Search arXiv's Atom API.  Preprints only; many carry a published DOI once accepted."""
    entries: list[ET.Element] = []
    for form in boolean_forms(query, field="all:"):
        response = fetcher.get(
            "arxiv",
            PROVIDERS["arxiv"].min_interval_s,
            "https://export.arxiv.org/api/query",
            params={"search_query": form, "max_results": limit, "sortBy": "relevance"},
        )
        if response is None:
            continue
        try:
            root = ET.fromstring(response.text)
        except ET.ParseError:
            LOG.warning("arxiv: unparseable Atom response for %s", query.qid)
            continue
        entries = root.findall("atom:entry", _ATOM_NS)
        if entries:
            break

    hits: list[Hit] = []
    for entry in entries:
        pdf_url = ""
        for link in entry.findall("atom:link", _ATOM_NS):
            if link.get("type") == "application/pdf":
                pdf_url = link.get("href", "")
        hits.append(
            Hit(
                provider="arxiv",
                query_id=query.qid,
                title=_text(entry, "atom:title"),
                doi=_text(entry, "arxiv:doi"),
                arxiv_id=norm_arxiv(_text(entry, "atom:id")),
                year=_text(entry, "atom:published")[:4],
                venue="arXiv",
                authors="; ".join(
                    _text(a, "atom:name") for a in entry.findall("atom:author", _ATOM_NS)
                ),
                abstract=_text(entry, "atom:summary"),
                pdf_url=pdf_url,
                landing_url=_text(entry, "atom:id"),
                is_oa=True,
                provider_id=norm_arxiv(_text(entry, "atom:id")),
            )
        )
    return hits


def _text(node: ET.Element, path: str) -> str:
    found = node.find(path, _ATOM_NS)
    return " ".join((found.text or "").split()) if found is not None else ""


def search_europepmc(fetcher: Fetcher, query: Query, limit: int) -> list[Hit]:
    """Search Europe PMC: PubMed + PMC + bioRxiv/medRxiv preprints in one index."""
    results: list[dict[str, Any]] = []
    for form in boolean_forms(query):
        response = fetcher.get(
            "europepmc",
            PROVIDERS["europepmc"].min_interval_s,
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
            params={"query": form, "format": "json", "pageSize": limit, "resultType": "core"},
        )
        if response is None:
            continue
        try:
            payload = response.json()
        except ValueError:
            continue
        results = payload.get("resultList", {}).get("result", [])
        if results:
            break

    hits: list[Hit] = []
    for record in results:
        pdf_url = ""
        for link in (record.get("fullTextUrlList") or {}).get("fullTextUrl", []):
            if link.get("availability") in {"Open access", "Free"} and link.get(
                "documentStyle"
            ) in {"pdf", "PDF"}:
                pdf_url = link.get("url", "")
                break
        hits.append(
            Hit(
                provider="europepmc",
                query_id=query.qid,
                title=(record.get("title") or "").strip().rstrip("."),
                doi=record.get("doi") or "",
                year=str(record.get("pubYear") or ""),
                venue=record.get("journalTitle") or record.get("source") or "",
                authors=record.get("authorString") or "",
                abstract=record.get("abstractText") or "",
                pdf_url=pdf_url,
                landing_url=f"https://europepmc.org/article/{record.get('source')}/{record.get('id')}",
                is_oa=(record.get("isOpenAccess") == "Y"),
                provider_id=f"{record.get('source')}:{record.get('id')}",
            )
        )
    return hits


def search_openalex(fetcher: Fetcher, query: Query, limit: int) -> list[Hit]:
    """Search OpenAlex.  Widest coverage, and carries ``best_oa_location`` directly."""
    results: list[dict[str, Any]] = []
    for form in relevance_forms(query):
        response = fetcher.get(
            "openalex",
            PROVIDERS["openalex"].min_interval_s,
            "https://api.openalex.org/works",
            params={"search": form, "per-page": limit, "mailto": EMAIL},
        )
        if response is None:
            continue
        try:
            payload = response.json()
        except ValueError:
            continue
        results = payload.get("results", [])
        if results:
            break

    hits: list[Hit] = []
    for record in results:
        oa_location = record.get("best_oa_location") or {}
        primary = record.get("primary_location") or {}
        source = (primary.get("source") or {}) if isinstance(primary, dict) else {}
        hits.append(
            Hit(
                provider="openalex",
                query_id=query.qid,
                title=(record.get("display_name") or "").strip(),
                doi=norm_doi(record.get("doi")),
                year=str(record.get("publication_year") or ""),
                venue=source.get("display_name") or "",
                authors="; ".join(
                    (a.get("author") or {}).get("display_name", "")
                    for a in (record.get("authorships") or [])[:12]
                ),
                pdf_url=oa_location.get("pdf_url") or "",
                landing_url=record.get("id") or "",
                is_oa=bool((record.get("open_access") or {}).get("is_oa")),
                provider_id=record.get("id") or "",
            )
        )
    return hits


def search_crossref(fetcher: Fetcher, query: Query, limit: int) -> list[Hit]:
    """Search Crossref bibliographic metadata — the publisher-of-record authority."""
    items: list[dict[str, Any]] = []
    for form in relevance_forms(query):
        response = fetcher.get(
            "crossref",
            PROVIDERS["crossref"].min_interval_s,
            "https://api.crossref.org/works",
            params={"query.bibliographic": form, "rows": limit, "mailto": EMAIL,
                    "select": "DOI,title,abstract,issued,container-title,author,URL"},
        )
        if response is None:
            continue
        try:
            payload = response.json()
        except ValueError:
            continue
        items = payload.get("message", {}).get("items", [])
        if items:
            break

    hits: list[Hit] = []
    for record in items:
        titles = record.get("title") or []
        issued = ((record.get("issued") or {}).get("date-parts") or [[""]])[0]
        containers = record.get("container-title") or []
        hits.append(
            Hit(
                provider="crossref",
                query_id=query.qid,
                title=(titles[0] if titles else "").strip(),
                doi=norm_doi(record.get("DOI")),
                year=str(issued[0] if issued else ""),
                venue=containers[0] if containers else "",
                authors="; ".join(
                    " ".join(filter(None, [a.get("given"), a.get("family")]))
                    for a in (record.get("author") or [])[:12]
                ),
                abstract=strip_jats(record.get("abstract")),
                landing_url=record.get("URL") or "",
                provider_id=norm_doi(record.get("DOI")),
            )
        )
    return hits


SEARCHERS: dict[str, Callable[[Fetcher, Query, int], list[Hit]]] = {
    "arxiv": search_arxiv,
    "europepmc": search_europepmc,
    "openalex": search_openalex,
    "crossref": search_crossref,
}


def resolve_oa_pdf(fetcher: Fetcher, doi: str) -> tuple[str, str]:
    """Ask Unpaywall for a legally-free PDF location.

    Returns:
        ``(url, reason)``.  ``url`` is empty when no OA copy exists, and ``reason`` then records
        why — which is what routes the DOI to the KAIST institutional track instead.
    """
    bare = norm_doi(doi)
    if not bare:
        return "", "no doi"
    response = fetcher.get(
        "unpaywall",
        UNPAYWALL_INTERVAL_S,
        f"https://api.unpaywall.org/v2/{bare}",
        params={"email": EMAIL},
    )
    if response is None:
        return "", "unpaywall unreachable"
    try:
        payload = response.json()
    except ValueError:
        return "", "unpaywall unparseable"
    if not payload.get("is_oa"):
        return "", "closed (route to KAIST full-text track)"
    location = payload.get("best_oa_location") or {}
    url = location.get("url_for_pdf") or location.get("url") or ""
    return (url, "oa") if url else ("", "oa but no pdf url")


# --------------------------------------------------------------------------------------
# relevance triage
# --------------------------------------------------------------------------------------

_JATS_RE = re.compile(r"<[^>]+>")


def strip_jats(text: str | None) -> str:
    """Flatten a Crossref JATS abstract to plain text and unescape provider HTML entities.

    Europe PMC titles arrive carrying ``&lt;i&gt;`` around species names and Crossref abstracts
    arrive as JATS XML.  Both would otherwise put markup between a search term and the word it is
    wrapped around, so an anchor or a disqualifier sitting inside italics would never match.
    """
    if not text:
        return ""
    return " ".join(unescape(_JATS_RE.sub(" ", unescape(text))).split())


def term_pattern(term: str) -> re.Pattern[str]:
    """Compile one anchor/disqualifier phrase with word boundaries and an optional plural.

    Boundaries matter more than they look.  ``plant`` must not fire on ``implant``, and ``wave``
    must fire on ``shock waves`` — the collision that put a ferromagnetic-soliton paper into the
    WAVE-complex query.  Internal whitespace spans any run of non-word characters, so
    ``type III secretion`` still matches ``type-III secretion``.

    THE OPTIONAL TRAILING ``s`` IS NOT COSMETIC.  Plain ``\\bplant\\b`` does not match ``plants``,
    so a disqualifier list written in the singular leaks every plural: ``axons``, ``dendrites`` and
    ``in plants`` all survived the first version, and the only reason the plant papers were caught
    at all is that they happened to also say ``Arabidopsis``.  A term already ending in ``s`` is
    left alone so ``mitosis`` does not become ``mitosiss``.
    """
    parts = [re.escape(p) for p in term.lower().split()]
    if parts and not parts[-1].endswith("s"):
        parts[-1] += "s?"
    return re.compile(r"\b" + r"\W+".join(parts) + r"\b", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class Triage:
    """The outcome of judging one hit against its query's subject matter.

    Attributes:
        relevance: ``candidate`` / ``needs_review`` / ``off_topic``.
        reason: Which term decided it, or why no decision was possible.
        matched: The specific phrase that fired, empty when nothing did.
    """

    relevance: str
    reason: str
    matched: str = ""


def triage(hit: Hit, query: Query, global_off_topic: Sequence[str]) -> Triage:
    """Decide whether a hit is about the object the query names, or merely shares its words.

    WHY A GATE THAT FWCI CANNOT REPLACE (PI 2026-07-28).  A citation metric ranks papers; it does
    not know what they are about.  An irrelevant but heavily cited paper passes any FWCI floor, and
    the measured noise in this corpus is exactly that shape: ``WAVE complex density`` returned
    many-body-localization and ferromagnetic-soliton physics plus a run of EEG and epileptology
    papers, ``filopodium tip complex`` returned bacterial type-III/VI secretion systems, and
    ``microtubule number cell`` returned mitotic spindle, plant and neuronal-axon work.  Every one
    of those shares a word with the query and none of them is about our object.

    THREE STATES, AND ``needs_review`` IS NOT A SOFT REJECT.  Auto-dropping what cannot be judged
    costs recall silently, so the only automatic exclusion is an EXPLICIT disqualifier.  Anything
    the available metadata cannot decide is preserved for a human, and only ``candidate`` is ever
    downloaded automatically.

    Order is deliberate: a disqualifier outranks an anchor.  A paper can be genuinely about
    microtubules and still be about the mitotic spindle, which is the wrong object for an
    interphase count — so seeing the wrong object is decisive even when the right word is present.

    Returns:
        A :class:`Triage`.  When a query declares no anchors the hit is ``needs_review``, never
        ``candidate``: a query with no stated subject cannot confirm one.
    """
    title = strip_jats(hit.title)
    abstract = strip_jats(hit.abstract)
    has_abstract = bool(abstract)

    for term in (*query.off_topic, *global_off_topic):
        if term_pattern(term).search(f"{title} {abstract}"):
            return Triage("off_topic", f"disqualifier {term!r} present", term)

    if not query.anchors:
        return Triage("needs_review", "query declares no anchors, so nothing can confirm it")

    # A PASSING MENTION IS NOT A SUBJECT.  Matching an anchor anywhere in a long abstract confirmed
    # the wrong things on the first run: an anaesthesia paper became an Arp2/3 candidate on one
    # clause of background, and — the case that settles it — a cell-motility paper matched on the
    # sentence "does NOT require Arp2/3", i.e. the anchor fired inside a denial.  A paper that is
    # ABOUT something says so in its title or returns to it; one recurrence is the least evidence
    # that distinguishes a subject from an aside, so that is the bar, and anything short of it is
    # preserved for review rather than fetched.
    single_mention = ""
    for term in query.anchors:
        pattern = term_pattern(term)
        if pattern.search(title):
            return Triage("candidate", f"anchor {term!r} in title", term)
        occurrences = len(pattern.findall(abstract))
        if occurrences >= 2:
            return Triage("candidate", f"anchor {term!r} recurs in abstract ({occurrences}x)", term)
        if occurrences == 1 and not single_mention:
            single_mention = term

    if single_mention:
        return Triage(
            "needs_review",
            f"anchor {single_mention!r} appears once in the abstract and not in the title — "
            "a passing mention, not a confirmed subject",
            single_mention,
        )
    if not has_abstract:
        return Triage("needs_review", "no abstract; the title alone did not confirm an anchor")
    return Triage("needs_review", "abstract present but no anchor matched")


def backfill_abstracts(fetcher: Fetcher, hits: Sequence[Hit]) -> int:
    """Fetch missing abstracts from Europe PMC by DOI, one request each.

    Crossref deposits an abstract for only a minority of records — one in six on a sample of this
    corpus — so without this step every Crossref-only hit lands in ``needs_review`` on metadata
    grounds rather than on its subject matter, and the reviewer inherits work the pipeline could
    have done.  Europe PMC is not credit-metered, unlike OpenAlex.

    Returns:
        How many abstracts were recovered.
    """
    targets = [h for h in hits if h.doi and not strip_jats(h.abstract)]
    recovered = 0
    for hit in targets:
        response = fetcher.get(
            "europepmc",
            PROVIDERS["europepmc"].min_interval_s,
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
            params={
                "query": f'DOI:"{norm_doi(hit.doi)}"',
                "format": "json",
                "pageSize": 1,
                "resultType": "core",
            },
        )
        if response is None:
            continue
        try:
            results = response.json().get("resultList", {}).get("result", [])
        except ValueError:
            continue
        if results and results[0].get("abstractText"):
            hit.abstract = results[0]["abstractText"]
            recovered += 1
    LOG.info("abstract backfill: %d recovered of %d missing", recovered, len(targets))
    return recovered


# --------------------------------------------------------------------------------------
# screening
# --------------------------------------------------------------------------------------

#: How many DOIs one OpenAlex ``filter=doi:a|b|c`` request may carry.  Provider-documented ceiling.
OPENALEX_BATCH = 50


def enrich_from_openalex(fetcher: Fetcher, hits: Sequence[Hit]) -> int:
    """Stamp year, work type and citation metrics onto hits that carry a DOI.

    One batched request per 50 DOIs.  Fields are copied verbatim; a work OpenAlex does not hold is
    left alone rather than being given defaults, because a zero written in place of "unknown" would
    be indistinguishable from a genuinely uncited paper at screening time.

    Returns:
        How many hits were enriched.
    """
    by_doi = {norm_doi(h.doi): h for h in hits if h.doi}
    dois = [d for d in by_doi if d]
    enriched = 0
    for start in range(0, len(dois), OPENALEX_BATCH):
        batch = dois[start : start + OPENALEX_BATCH]
        response = fetcher.get(
            "openalex",
            PROVIDERS["openalex"].min_interval_s,
            "https://api.openalex.org/works",
            params={
                "filter": "doi:" + "|".join(batch),
                "per-page": OPENALEX_BATCH,
                "select": "doi,publication_year,type,fwci,cited_by_count,"
                "citation_normalized_percentile,primary_location",
                "mailto": EMAIL,
            },
        )
        if response is None:
            continue
        try:
            payload = response.json()
        except ValueError:
            continue
        for work in payload.get("results", []):
            hit = by_doi.get(norm_doi(work.get("doi")))
            if hit is None:
                continue
            if work.get("publication_year"):
                hit.year = str(work["publication_year"])
            hit.work_type = work.get("type") or ""
            hit.cited_by_count = work.get("cited_by_count")
            hit.fwci = work.get("fwci")
            hit.citation_percentile = (work.get("citation_normalized_percentile") or {}).get("value")
            source = ((work.get("primary_location") or {}).get("source")) or {}
            if source.get("display_name") and not hit.venue:
                hit.venue = source["display_name"]
            enriched += 1
    return enriched


def screen(hit: Hit, *, min_year: int | None) -> str:
    """Apply the only HARD filter: publication year.  Returns "" to keep the hit.

    YEAR IS THE ONLY GATE, BY PI DECISION 2026-07-28 (option C).  Two things were tried as a quality
    floor and neither belongs here.

    Impact Factor cannot be applied at all: JIF is proprietary with no free API, and OpenAlex's
    journal-level ``2yr_mean_citedness`` — same formula — was measured unusable at a threshold of 5
    on this corpus.  It reads 0.47 for Biophysical Journal against a published ~3.2 and 1.70 for
    eLife against ~6.4, while resolving for only 34 of 254 candidate venues; applying it would have
    deleted 41 papers sitting at twice their field average, including a 278-citation Biophysical
    Journal paper and an 8,291-citation one in Scientific Reports.

    FWCI is real and is still collected, but it is NOT a gate either, because it answers the wrong
    question.  A citation metric ranks papers; it does not know what they are about, so an
    irrelevant but heavily cited paper clears any floor.  That is the defect the PI found in the
    first year-only dry run — ``WAVE complex density`` returning many-body-localization physics and
    EEG papers — and no FWCI threshold removes it.  Subject-matter fitness is :func:`triage`'s job;
    ``fwci`` survives only as the ORDER in which a human reads the results.
    """
    if min_year is None:
        return ""
    year = int(hit.year) if hit.year.isdigit() else None
    if year is None:
        return "no publication year"
    if year < min_year:
        return f"year {year} < {min_year}"
    return ""


# --------------------------------------------------------------------------------------
# corpus state
# --------------------------------------------------------------------------------------


def known_corpus_keys() -> set[str]:
    """Dedup keys for everything the corpus already has, from duckdb and from prior clip runs.

    A missing ``kb.duckdb`` is not an error: the clipper still works, it simply cannot dedup against
    the ingested corpus and says so.
    """
    keys: set[str] = set()
    if DB_PATH.exists():
        try:
            import duckdb

            con = duckdb.connect(str(DB_PATH), read_only=True)
            for table, has_title in (("paper_refs", True), ("source_evidence", False)):
                try:
                    columns = "doi, title" if has_title else "doi, NULL"
                    for doi, title in con.execute(f"SELECT {columns} FROM {table}").fetchall():
                        if doi:
                            keys.add(f"doi:{norm_doi(doi)}")
                        if title:
                            keys.add(f"title:{norm_title(title)}")
                except Exception as exc:  # table absent in an older db
                    LOG.debug("dedup: %s unavailable (%s)", table, exc)
            con.close()
        except Exception as exc:
            LOG.warning("dedup: kb.duckdb unreadable (%s) — deduping against disk only", exc)
    else:
        LOG.warning("dedup: %s absent — deduping against prior clip manifests only", DB_PATH)

    for manifest in REF_DIR.glob("*/_clip_manifest.json"):
        try:
            for record in json.loads(manifest.read_text()):
                if record.get("doi"):
                    keys.add(f"doi:{norm_doi(record['doi'])}")
                if record.get("arxiv_id"):
                    keys.add(f"arxiv:{norm_arxiv(record['arxiv_id'])}")
                if record.get("title"):
                    keys.add(f"title:{norm_title(record['title'])}")
        except (OSError, ValueError) as exc:
            LOG.warning("dedup: %s unreadable (%s)", manifest, exc)
    return keys


def merge_hits(hits: Iterable[Hit]) -> list[Hit]:
    """Collapse hits that name the same paper, keeping the record with the most fields filled.

    Providers are complementary, not redundant: Crossref has the authoritative DOI but no PDF,
    OpenAlex has the OA location, arXiv has the preprint.  Merging keeps whichever field is
    non-empty rather than letting the first provider seen win.
    """
    merged: dict[str, Hit] = {}
    for hit in hits:
        keys = hit.dedup_keys
        if not keys:
            continue
        existing_key = next((k for k in keys if k in merged), None)
        if existing_key is None:
            for key in keys:
                merged[key] = hit
            continue
        target = merged[existing_key]
        for name in ("title", "doi", "arxiv_id", "year", "venue", "authors", "abstract",
                     "pdf_url", "landing_url"):
            if not getattr(target, name) and getattr(hit, name):
                setattr(target, name, getattr(hit, name))
        if target.is_oa is None:
            target.is_oa = hit.is_oa
        if hit.provider not in target.provider.split("+"):
            target.provider = f"{target.provider}+{hit.provider}"
        for key in keys:
            merged.setdefault(key, target)

    seen: list[Hit] = []
    for hit in merged.values():
        if not any(hit is other for other in seen):
            seen.append(hit)
    return seen


# --------------------------------------------------------------------------------------
# fetching
# --------------------------------------------------------------------------------------


def fetch_pdf(fetcher: Fetcher, url: str, destination: pathlib.Path) -> tuple[bool, str]:
    """Download one PDF, rejecting anything that is not actually a PDF.

    Returns:
        ``(ok, reason)``.  On failure nothing is written.
    """
    response = fetcher.get("pdf", 0.5, url, stream=True, allow_redirects=True)
    if response is None:
        return False, "unreachable"
    body = bytearray()
    for block in response.iter_content(chunk_size=65536):
        body.extend(block)
        if len(body) > MAX_PDF_BYTES:
            return False, f"exceeds {MAX_PDF_BYTES} byte ceiling"
    if not bytes(body[:4]) == PDF_MAGIC:
        return False, "not a PDF (landing page or paywall interstitial)"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(bytes(body))
    return True, "ok"


def append_analysis_titles(records: Sequence[dict[str, Any]], tag: str) -> int:
    """Give ``references_ingest`` the real titles for the files this run added.

    ``references_ingest`` looks a PDF's title up in ``references/analysis/_manifest.json`` by
    FILENAME and otherwise falls back to the first long line on page 1 — which on a typeset paper is
    often the journal banner.  We know the real title from the provider, so we record it.  Existing
    entries are never rewritten.
    """
    if not records:
        return 0
    try:
        existing = json.loads(ANALYSIS_MANIFEST.read_text()) if ANALYSIS_MANIFEST.exists() else []
    except (OSError, ValueError) as exc:
        LOG.warning("analysis manifest unreadable (%s) — not appending titles", exc)
        return 0
    have = {entry.get("path") for entry in existing}
    next_n = max((entry.get("n", 0) for entry in existing), default=0) + 1
    added = 0
    for record in records:
        name = pathlib.Path(record["local_path"]).name
        if name in have:
            continue
        existing.append(
            {
                "n": next_n,
                "path": name,
                "pages": None,
                "title": record.get("title", ""),
                "slug": slugify(record.get("title", ""), 60),
                "clip_source": tag,
                "doi": record.get("doi", ""),
            }
        )
        next_n += 1
        added += 1
    if added:
        ANALYSIS_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        ANALYSIS_MANIFEST.write_text(json.dumps(existing, ensure_ascii=False, indent=1))
    return added


# --------------------------------------------------------------------------------------
# queries
# --------------------------------------------------------------------------------------


def load_queries(path: pathlib.Path) -> tuple[list[Query], tuple[str, ...]]:
    """Read ``clip_queries.yaml`` into ``Query`` objects, rejecting an entry with no terms."""
    payload = yaml.safe_load(path.read_text())
    if not isinstance(payload, dict) or "queries" not in payload:
        raise ValueError(f"{path} has no 'queries:' list")
    global_off_topic = tuple((payload.get("defaults") or {}).get("off_topic") or ())
    queries: list[Query] = []
    for raw in payload["queries"]:
        terms = raw.get("terms") or {}
        all_terms = tuple(terms.get("all") or ())
        any_terms = tuple(terms.get("any") or ())
        if not all_terms and not any_terms:
            raise ValueError(f"query {raw.get('id')!r} has no terms — it would match everything")
        queries.append(
            Query(
                qid=str(raw["id"]),
                all_terms=all_terms,
                any_terms=any_terms,
                component=raw.get("component", "") or "",
                connector=raw.get("connector", "") or "",
                connector_family=raw.get("connector_family", "") or "",
                parameter_axis=raw.get("parameter_axis", "") or "",
                pi_gap_card=raw.get("pi_gap_card", "") or "",
                absence=raw.get("absence", "") or "",
                cell_state=raw.get("cell_state", "") or "",
                anchors=tuple(raw.get("anchors") or ()),
                off_topic=tuple(raw.get("off_topic") or ()),
                note=raw.get("note", "") or "",
            )
        )
    return queries, global_off_topic


# --------------------------------------------------------------------------------------
# run
# --------------------------------------------------------------------------------------


def run(
    *,
    queries: Sequence[Query],
    providers: Sequence[str],
    limit: int,
    tag: str,
    apply: bool,
    min_year: int | None = None,
    global_off_topic: Sequence[str] = (),
    backfill: bool = True,
) -> dict[str, Any]:
    """Execute the clip: search, dedup, triage, screen, optionally fetch the confirmed subset."""
    fetcher = Fetcher()
    known = known_corpus_keys()
    LOG.info("dedup baseline: %d keys from corpus + prior clips", len(known))

    records: list[dict[str, Any]] = []
    counts: dict[str, int] = {}

    # Phase 1 — search and dedup every query BEFORE screening any of it.  Screening per query as we
    # went would let a provider that runs out of quota midway screen the later queries under a
    # different, harsher rule than the earlier ones, and the manifest would not show that the rule
    # had changed.  Doing all the searching first also lets enrichment batch across queries.
    paired: list[tuple[Hit, Query]] = []
    fresh: list[Hit] = []
    for query in queries:
        raw: list[Hit] = []
        for provider_name in providers:
            searcher = SEARCHERS[provider_name]
            found = searcher(fetcher, query, limit)
            LOG.info("[%s] %s -> %d hits", query.qid, provider_name, len(found))
            raw.extend(found)

        for hit in merge_hits(raw):
            if any(key in known for key in hit.dedup_keys):
                hit.status = "duplicate"
                hit.reason = "already in corpus or a prior clip run"
            else:
                fresh.append(hit)
            for key in hit.dedup_keys:
                known.add(key)
            paired.append((hit, query))

    # Phase 2 — enrich only what survived dedup.  Paying for metrics on a paper the corpus already
    # holds is a request nobody reads, and OpenAlex's quota is a shared daily budget.  ``fwci`` is
    # collected as READING ORDER, never as a gate — see screen() for why a citation metric cannot
    # do a relevance gate's job.  The abstract backfill exists so that a hit lands in needs_review
    # for what it is ABOUT rather than for which provider happened to return it.
    if fresh:
        n_enriched = enrich_from_openalex(fetcher, fresh)
        LOG.info("enriched %d/%d fresh hits from OpenAlex", n_enriched, len(fresh))
        if backfill:
            backfill_abstracts(fetcher, fresh)

    # Phase 3 — subject-matter triage over the WHOLE run, for the same reason phase 1 batches:
    # judging early queries on richer evidence than late ones is a rule change the manifest would
    # not show.
    query_of = {id(hit): query for hit, query in paired}
    for hit in fresh:
        verdict = triage(hit, query_of[id(hit)], global_off_topic)
        hit.relevance = verdict.relevance
        hit.relevance_reason = verdict.reason
        hit.relevance_matched = verdict.matched
    LOG.info("triage: %s", {
        state: sum(1 for h in fresh if h.relevance == state)
        for state in ("candidate", "needs_review", "off_topic")
    })

    out_dir = REF_DIR / tag
    for hit, query in paired:
        screened_out = "" if hit.status == "duplicate" else screen(hit, min_year=min_year)
        if hit.status == "duplicate":
            pass
        elif screened_out:
            hit.status = "screened"
            hit.reason = screened_out
        elif hit.relevance != "candidate":
            # needs_review is PRESERVED, never auto-dropped and never auto-fetched.  Only a
            # confirmed subject match is downloaded without a human in the loop (PI 2026-07-28):
            # auto-dropping what cannot be judged loses recall silently, and auto-fetching it
            # spends bandwidth on the charge-density-wave papers this gate exists to catch.
            hit.status = hit.relevance
            hit.reason = hit.relevance_reason
        else:
            if apply:
                url = hit.pdf_url
                reason = "provider-supplied OA pdf"
                if not url and hit.doi:
                    url, reason = resolve_oa_pdf(fetcher, hit.doi)
                if not url:
                    hit.status = "no_oa"
                    hit.reason = reason or "no open-access location"
                else:
                    stem = f"{slugify(hit.title, 50)}-{(norm_doi(hit.doi) or hit.arxiv_id or 'na')[-8:]}"
                    destination = out_dir / f"{slugify(stem, 64)}.pdf"
                    ok, why = fetch_pdf(fetcher, url, destination)
                    hit.status = "downloaded" if ok else "fetch_failed"
                    hit.reason = why
                    if ok:
                        hit.local_path = str(destination.relative_to(REPO_ROOT))
            else:
                hit.status = "candidate"
                hit.reason = "dry run — nothing fetched"

        counts[hit.status] = counts.get(hit.status, 0) + 1
        record = asdict(hit)
        record.update(query.provenance())
        records.append(record)

    payload = {
        "schema": "ffn-ac-paper-clip-v1",
        "tag": tag,
        "applied": apply,
        "n_queries": len(queries),
        "providers": {
            name: {"rate_basis": PROVIDERS[name].rate_basis, "covers": PROVIDERS[name].covers}
            for name in providers
        },
        "unpaywall_rate_basis": UNPAYWALL_RATE_BASIS,
        "screen": {
            "policy": "PI 2026-07-28 option C — year is the ONLY hard filter; fwci is reading "
                      "order; subject-matter fitness is a separate three-state relevance gate, "
                      "and only 'candidate' is fetched without a human in the loop",
            "min_year": min_year,
            "abstract_backfill": backfill,
            "n_global_off_topic_terms": len(global_off_topic),
            "provider_exhausted": dict(fetcher.exhausted),
            "provider_quota_remaining": dict(fetcher.quota_remaining),
            "quality_metric": "OpenAlex fwci — RECORDED FOR ORDERING ONLY, never a gate",
            "why_not_impact_factor": (
                "JIF is proprietary with no free API. OpenAlex's journal 2yr_mean_citedness uses the "
                "same formula but was measured unusable at a threshold of 5 on this corpus "
                "(eLife 1.70 vs published ~6.4; Biophysical Journal 0.47 vs ~3.2; journal name "
                "resolved for only 34 of 254 candidate venues)."
            ),
        },
        "counts": counts,
        "records": records,
    }

    if apply:
        downloaded = [r for r in records if r["status"] == "downloaded"]
        if downloaded:
            out_dir.mkdir(parents=True, exist_ok=True)
            n_titles = append_analysis_titles(downloaded, tag)
            LOG.info("appended %d titles to %s", n_titles, ANALYSIS_MANIFEST.name)
        manifest_path = out_dir / "_clip_manifest.json"
        if manifest_path.exists():
            try:
                previous = json.loads(manifest_path.read_text())
            except ValueError:
                previous = []
            records = previous + records
            payload["records"] = records
        out_dir.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(records, ensure_ascii=False, indent=1))
        LOG.info("manifest -> %s (%d records)", manifest_path, len(records))

    return payload


def _cell(value: Any, limit: int = 0) -> str:
    """Render one markdown table cell.

    Every cell goes through this.  A raw ``|`` anywhere in a row silently shifts every column after
    it, and the reviewer sees a table whose columns no longer mean what the header says — which is
    exactly what ``parameter_axis: "observable: total traction |F| vs substrate stiffness"`` did on
    the first full run.  Escaping only the title, as the first version did, is not enough: any field
    can carry a pipe.
    """
    text = " ".join(str(value or "").split()).replace("|", "\\|")
    return text[:limit] if limit else text


def write_report(payload: dict[str, Any], path: pathlib.Path) -> None:
    """Write the reviewable candidate manifest, in the shape ``harvest_ops`` uses.

    Also writes a sibling ``.json`` carrying the full records, so a dry run's result can be
    re-analysed without paying for the whole search again.
    """
    path.with_suffix(".json").write_text(json.dumps(payload, ensure_ascii=False, indent=1))
    lines = [
        f"# Paper clip candidates — tag `{payload['tag']}`",
        "",
        f"Applied: **{payload['applied']}** · queries: {payload['n_queries']} · "
        f"providers: {', '.join(payload['providers'])}",
        "",
        "| status | n |",
        "|---|---|",
    ]
    for status, n in sorted(payload["counts"].items(), key=lambda kv: -kv[1]):
        lines.append(f"| {status} | {n} |")
    screen_cfg = payload.get("screen", {})
    lines += [
        "",
        f"Policy: {screen_cfg.get('policy', 'n/a')}",
        "",
        f"min_year={screen_cfg.get('min_year')} · abstract_backfill="
        f"{screen_cfg.get('abstract_backfill')} · global off-topic terms="
        f"{screen_cfg.get('n_global_off_topic_terms')}",
        "",
        f"> **fwci is ordering only, not a filter**, and Impact Factor is not used at all. "
        f"{screen_cfg.get('why_not_impact_factor','')}",
        "",
        "## Candidates — grouped by relevance, ordered by fwci within each group",
        "",
        "`candidate` is auto-fetched. `needs_review` is PRESERVED for a human and never fetched "
        "automatically. `off_topic` matched an explicit disqualifier.",
        "",
        "| relevance | query | year | fwci | cites | venue | title | doi | status | why |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]

    # Ordering IS the product here: fwci lost its gate role, so its only remaining job is to put the
    # most-cited confirmed papers in front of the reviewer first.  Unknown fwci sorts last rather
    # than as zero, so "not scored" never reads as "poorly cited".
    # Group by the FINAL STATUS, not by the triage verdict.  Grouping by relevance put 51 rows under
    # a heading that says "auto-fetched" while the year filter had already removed them, so the
    # table's own caption was false for one row in eight.  `relevance` stays as a column because it
    # still explains WHY a surviving row is where it is.
    rank = {"candidate": 0, "needs_review": 1, "off_topic": 2, "screened": 3, "no_oa": 4}
    def sort_key(record: dict[str, Any]) -> tuple[Any, ...]:
        fwci = record.get("fwci")
        return (
            rank.get(record.get("status", ""), 5),
            record.get("query_id", ""),
            -fwci if isinstance(fwci, (int, float)) else float("inf"),
        )

    for record in sorted(
        (r for r in payload["records"] if r["status"] != "duplicate"), key=sort_key
    ):
        fwci = record.get("fwci")
        lines.append(
            "| " + " | ".join([
                _cell(record.get("relevance")),
                _cell(record.get("query_id")),
                _cell(record.get("year")),
                f"{fwci:.2f}" if isinstance(fwci, (int, float)) else "",
                _cell(record.get("cited_by_count")),
                _cell(record.get("venue"), 26),
                _cell(record.get("title"), 76),
                _cell(record.get("doi")),
                _cell(record.get("status")),
                _cell(record.get("relevance_reason") or record.get("reason"), 46),
            ]) + " |"
        )
    path.write_text("\n".join(lines) + "\n")
    LOG.info("report -> %s", path)


def check() -> int:
    """Report clipped PDFs that ``references_ingest`` has not yet turned into chunks."""
    pdfs = [
        pdf
        for directory in sorted(REF_DIR.iterdir())
        if directory.is_dir() and (directory / "_clip_manifest.json").exists()
        for pdf in sorted(directory.glob("*.pdf"))
    ]
    if not pdfs:
        print("no clipped PDFs on disk yet")
        return 0
    ingested: set[str] = set()
    if DB_PATH.exists():
        try:
            import duckdb

            con = duckdb.connect(str(DB_PATH), read_only=True)
            ingested = {row[0] for row in con.execute("SELECT path FROM paper_refs").fetchall()}
            con.close()
        except Exception as exc:
            print(f"kb.duckdb unreadable ({exc}) — cannot compare")
            return 1
    pending = [p for p in pdfs if f"{p.parent.name}/{p.name}" not in ingested]
    print(f"clipped PDFs: {len(pdfs)}   ingested: {len(pdfs) - len(pending)}   PENDING: {len(pending)}")
    if pending:
        print("run:  python references_ingest.py")
        for p in pending[:20]:
            print(f"   - {p.parent.name}/{p.name}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--queries", type=pathlib.Path, default=DEFAULT_QUERY_FILE,
                        help="query set anchored to the engine vocabulary")
    parser.add_argument("--query", help="ad-hoc free-text query (bypasses the query file)")
    parser.add_argument("--id", action="append", default=[], help="run only these query ids")
    parser.add_argument("--providers", default=",".join(PROVIDERS),
                        help=f"comma list from: {', '.join(PROVIDERS)}")
    parser.add_argument("--limit", type=int, default=25, help="max hits per provider per query")
    parser.add_argument("--tag", default="clipped", help="subdir under references/ (= the ingest 'source')")
    parser.add_argument("--min-year", type=int, default=None,
                        help="drop papers published before this year (exact, no proxy involved)")
    parser.add_argument("--no-backfill", action="store_true",
                        help="skip the Europe PMC abstract backfill; more hits then land in "
                             "needs_review for missing metadata rather than on their subject")
    parser.add_argument("--apply", action="store_true", help="actually fetch the OA subset")
    parser.add_argument("--check", action="store_true", help="report clipped-but-not-ingested")
    parser.add_argument("--report", type=pathlib.Path, default=HERE / "CLIP_CANDIDATES.md")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    if args.check:
        return check()

    providers = [p.strip() for p in args.providers.split(",") if p.strip()]
    unknown = [p for p in providers if p not in PROVIDERS]
    if unknown:
        parser.error(f"unknown provider(s): {unknown}; choose from {list(PROVIDERS)}")

    if args.query:
        queries = [Query(qid="adhoc", all_terms=(args.query,))]
        global_off_topic = ()
    else:
        if not args.queries.exists():
            parser.error(f"{args.queries} not found")
        queries, global_off_topic = load_queries(args.queries)
        if args.id:
            wanted = set(args.id)
            queries = [q for q in queries if q.qid in wanted]
            missing = wanted - {q.qid for q in queries}
            if missing:
                parser.error(f"no such query id(s): {sorted(missing)}")

    payload = run(queries=queries, providers=providers, limit=args.limit,
                  tag=args.tag, apply=args.apply, min_year=args.min_year,
                  global_off_topic=global_off_topic, backfill=not args.no_backfill)
    write_report(payload, args.report)

    print(f"\n{'APPLIED' if args.apply else 'DRY RUN'} — {payload['counts']}")
    if not args.apply:
        print(f"review {args.report}, then re-run with --apply")
    return 0


if __name__ == "__main__":
    sys.exit(main())
