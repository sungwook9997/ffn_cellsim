# External training corpus

This directory is the versioned control plane for Aleph's external evidence corpus. It does **not** contain publisher PDFs or a claim that a paper's conclusions are already Aleph truth.

The corpus is organised by independent source family (normally a DOI), not by downloaded file. Git stores provenance, access and licence state, Aleph roles, acquisition priority, and extraction status. Full text and supplements live in the ignored local object store `data/external_training/objects/<sha256>`; institution-only material must never be committed or redistributed.

## Scope

The initial target is 1,200 independent source families:

- 300–500 Tier A sources with protocol- and observable-level extraction;
- 500–1,000 Tier B sources for priors, boundary conditions, failure modes, and domain coverage;
- explicit negative and contradictory evidence, rather than only canonical positive results.

Coverage is intentionally broader than cortical tension: cortex and membrane mechanics, adhesion and traction, nucleus, osmotic/poroelastic regulation, protrusion and migration, ECM and invasion, plus observation and inference methods for IF, WB, PCR, PIV, TFM, AFM, microscopy, and simulation-based inference.

## Files

- `manifest.schema.json`: contract for every source record.
- `sources/seed_primary_sources.jsonl`: verified or explicitly qualified starter records.
- `discovery_queries.json`: dated, reproducible Europe PMC discovery queries and quotas.
- `discovery_queries_expansion.json`: twelve complementary domains used to add 1,800 non-overlapping candidates.
- `discovery_queries_historical.json`: a publication-date-bounded 2010–2021 backfill that corrects recent-year bias.
- `discovery_queries_cell_states.json`: eight cell-type/state families for the Outer Library rather than a cortex-only corpus.
- `discovery_queries_cell_states_2010s.json`: the same cell-state families bounded to 2010–2019.
- `tools/discover_europe_pmc.py`: metadata-only retrieval and source-family deduplication.
- `coverage/targets.csv`: source-family quotas, not fabricated sample counts.
- `queues/acquisition.csv`: ordered acquisition work.
- `queues/journal_metrics_2026-08-05.csv`: deduplicated licensed-JCR lookup queue; blank values are intentionally unknown.
- `ACQUISITION_POLICY.md`: copyright, provenance, extraction, and leakage rules.
- `QUALITY_POLICY.md`: post-2010, JIF-screening, article-quality, and exception rules.

Records begin as `proposed`. Numeric values may enter training only after source-level extraction, unit normalisation, protocol capture, and an independent review gate elsewhere in Aleph.
