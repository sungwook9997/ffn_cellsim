# Handoff prompts — TAG corpus PDF fetching (for Sonnet sessions)

Two independent, mechanical tasks (no heavy reasoning needed). They operate on
disjoint subsets of `tag_kb/download_log.json` — run them in parallel Sonnet
sessions if you like.

- **(b)** retries the ~90 `status: dl-fail` rows (Unpaywall said OA but the URL
  wasn't a direct PDF).
- **(a)** fetches the ~46 `status: closed` rows (paywalled → KAIST full-text via gbook).

Both end with the same ingest + commit steps. PDFs go to
`aleph/references/downloaded/<citation_key>.pdf`; they're gitignored — commit
only the scripts + `download_log.json` + `references/tag_corpus.json`.

---

## Prompt (b) — recover the OA download failures

```
You're picking up a mechanical PDF-fetching task in the ffn_cellsim repo
(/Users/sw1/ffn_cellsim). Read CLAUDE.md first.

Context: a TAG reference corpus lives in aleph/outputs/tag_kb/. An Unpaywall
downloader (download_pdfs.py) already saved open-access PDFs to
aleph/references/downloaded/<citation_key>.pdf. But ~90 DOIs that Unpaywall
reported as OA failed because best_oa_location returned an HTML landing page,
redirect, or 403 instead of a direct PDF — logged in tag_kb/download_log.json
with status "dl-fail".

Task: recover as many of those ~90 as possible. ONLY fetch the legally-free OA
copies Unpaywall points to — never scrape paywalled content.

Setup: conda activate ffn_sim (Python 3.13; requests, pypdf, fitz available).

Write tag_kb/retry_oa.py that, for each download_log.json entry with
status=="dl-fail":
  1. Re-query Unpaywall (https://api.unpaywall.org/v2/{doi}?email=sungwook999@gmail.com)
     and try ALL entries in oa_locations (not just best_oa_location), preferring
     url_for_pdf, then url.
  2. Fetch with a browser User-Agent and follow redirects. If the response is
     HTML (not %PDF), parse it for the PDF link: <meta name="citation_pdf_url">,
     <link type="application/pdf">, or <a href> ending in .pdf. Handle common
     publishers: PMC -> https://www.ncbi.nlm.nih.gov/pmc/articles/PMCxxxxxxx/pdf/,
     MDPI -> append /pdf, Frontiers -> /pdf, bioRxiv -> .full.pdf.
  3. Only save if the bytes start with b"%PDF". Write to
     references/downloaded/<citation_key>.pdf and set the log status to
     "downloaded"; else leave "dl-fail". Rate-limit ~0.3s.

Then: `python references_ingest.py` (ingests references/ incl downloaded/ into
paper_chunks + FTS) and `python notion_to_duckdb.py` (non-destructive refresh of
kb.duckdb node tables). Report recovered count + new paper_chunks total.

Commit on branch phase1/h3-cortex (do NOT push to ffn/foundation): retry_oa.py +
download_log.json + aleph/references/tag_corpus.json. PDFs are gitignored. Use
`git commit -F <msgfile>` — never heredoc a commit message (it corrupts stdin).
Other sessions commit to this branch concurrently; stage ONLY your own files.
```

---

## Prompt (a) — fetch the paywalled set via gbook (KAIST full-text)

```
You're picking up a mechanical PDF-fetching task in the ffn_cellsim repo
(/Users/sw1/ffn_cellsim). Read CLAUDE.md first.

Context: tag_kb/download_log.json lists ~46 SourceEvidence DOIs with status
"closed" — paywalled, no open-access copy. They need KAIST institutional
full-text, reachable from the gbook workstation: gbook's wifi is on the KAIST
network, so its IP unlocks publisher full-text. IMPORTANT (project memory): a
bare curl 403s; you must `ssh gbook` and curl with a browser User-Agent.

Task: for each "closed" DOI, fetch the institutional full-text PDF via gbook and
save it to aleph/references/downloaded/<citation_key>.pdf (same convention the
OA downloader uses). This is inherently partial — get what you can, log the rest.

Setup: conda activate ffn_sim. The list (citation_key, doi) is in
download_log.json (status=="closed"); also queryable via
`SELECT citation_key, doi FROM source_evidence` against tag_kb/kb.duckdb.

Per DOI: resolve https://doi.org/{doi} to the publisher article page, then on
gbook (`ssh gbook '...'`) curl the PDF with a browser UA and follow redirects
(-A 'Mozilla/5.0 (X11; Linux x86_64) ... Chrome/...' -L). PDF URL patterns vary
by publisher (Elsevier/ScienceDirect, Springer, Wiley, Nature, etc.) — derive
per publisher. Some (ScienceDirect) bot-block even institutional curl: flag those
"manual-browser" in the log rather than forcing it. Validate bytes start with
%PDF before saving. Update download_log.json status to "downloaded" or
"manual-browser".

Then: `python references_ingest.py` + `python notion_to_duckdb.py` to ingest +
refresh. Report fetched count.

Commit on branch phase1/h3-cortex (NOT ffn/foundation): your fetch script +
download_log.json + aleph/references/tag_corpus.json (PDFs gitignored). Use
`git commit -F <msgfile>` (never heredoc). Other sessions commit to this branch
concurrently; stage ONLY your own files.
```

---

Verify either result with the TAG engine:
`cd aleph/outputs/tag_kb && python tag_query.py "How many paper_refs are linked to SourceEvidence now?"`
