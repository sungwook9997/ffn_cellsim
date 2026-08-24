# Europe PMC OA acquisition receipts

This directory contains versioned receipts and summaries for payloads acquired
from Europe PMC's official `fullTextXML` endpoint. It contains no article XML,
PDF, publisher credential, cookie, or signed URL.

Eligibility is fail-closed: a discovery row must contain both
`is_open_access_provider_flag: true` and a valid PMCID. Successful payloads are
parsed as JATS-like article XML and stored outside Git at
`data/external_training/objects/<sha256>`. The receipt records the payload hash,
size, source URL, access classification, provider licence metadata, licence
elements found in the XML, and retrieval timestamp.

Run the collector from the repository root with explicit immutable inputs:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python \
  corpus/external_training/tools/acquire_europe_pmc_oa.py \
  --input corpus/external_training/snapshots/europe_pmc_2026-08-05.jsonl \
  --input corpus/external_training/snapshots/europe_pmc_expansion_2026-08-05.jsonl \
  --input corpus/external_training/snapshots/europe_pmc_historical_2010_2021.jsonl \
  --input corpus/external_training/snapshots/europe_pmc_cell_states_2026-08-05.jsonl \
  --input corpus/external_training/snapshots/europe_pmc_cell_states_2010s.jsonl \
  --object-root data/external_training/objects \
  --receipts corpus/external_training/acquisition/europe_pmc_oa_2026-08-05.jsonl \
  --summary corpus/external_training/acquisition/europe_pmc_oa_2026-08-05.summary.json \
  --workers 4
```

The receipt and summary are atomically replaced every checkpoint. Re-running
skips every PMCID with an existing receipt. Pass `--retry-failures` explicitly
to retry prior failures; successful records are never fetched again.

The default is one request lane. The measured run used four bounded lanes, each
with its own post-request delay; this remains well below a high-throughput bulk
harvester and avoids unbounded task submission.
