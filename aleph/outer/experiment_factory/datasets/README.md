# Public dataset/accession registry

This lane connects the two-pass OA JATS corpus to public repository identifiers. It is a proposed
discovery and leakage-control layer, not an evidence promotion and not proof that a paper generated
or reused a particular processed dataset.

Artifacts:

- `source_accession_registry.jsonl`: exact source-to-accession bindings with XML payload, locator and
  locator-text SHA-256 values;
- `dataset_manifests.jsonl`: one normalized proposed manifest per public identifier;
- `official_metadata_receipts.jsonl`: bounded official-API resolution receipts, without response
  bodies;
- `acquisition_queue.csv`: gold-candidate and modality-prioritized queue for deliberate raw-data
  review;
- `pilot_download_receipts.jsonl`: digest/licence/URL receipts for the conservative local-only pilot;
- `summary.json` and `SHA256SUMS.json`: measured counts and reproducibility bindings.

Rebuild the registry from local OA XML and existing official receipts:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python \
  corpus/external_training/experiment_factory/datasets/build_registry.py
```

Official metadata resolution is explicitly bounded:

```bash
.../python corpus/external_training/experiment_factory/datasets/build_registry.py \
  --resolve-official 100
```

The resolver uses only official NCBI, EMBL-EBI, Zenodo, Figshare, Dryad and OSF APIs implemented in
the script. It has no publisher, login, CAPTCHA or subscription fallback. `not_attempted`,
`not_found`, `failed`, and fields not exposed by an API remain typed states.

The pilot downloader has a hard 5,000,000-byte total cap and a fixed four-dataset plan. It accepts
only API-declared allowlisted open licences, official download URLs and files whose declared sizes
fit the cap. Payloads live under ignored `data/external_training/experiment_factory/datasets/**`;
Git receives receipts only.
