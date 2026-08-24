# Public dataset/accession registry — measured report

## Outcome

Both local OA JATS stores were scanned: 5,675/5,675 XML payloads, 850,654,786 input bytes. Exact
public repository identifiers were found in 885 sources (15.59%). The frozen 300-paper
`gold_candidate` set has 98 sources with at least one identifier (32.67%).

| Measure | Result |
|---|---:|
| Source↔accession bindings | 1,982 |
| Unique normalized dataset manifests | 1,801 |
| Gold-candidate bindings / unique manifests | 319 / 307 |
| Dataset manifests referenced by multiple papers | 123 |
| Maximum source count for one accession | 7 |
| Missing XML or digest locator | 0 |

Every binding retains the source-family ID, PMCID, exact XML SHA-256, structural XML locator and
locator-text SHA-256. No source text or XML is copied into the registry.

## Providers

| Provider | Unique manifests |
|---|---:|
| GEO | 1,110 |
| Zenodo | 233 |
| PRIDE | 103 |
| BioProject | 98 |
| Figshare | 62 |
| SRA | 39 |
| ArrayExpress | 30 |
| Dryad | 29 |
| BioSample | 25 |
| Mendeley Data | 24 |
| BioImage Archive | 11 |
| OSF | 11 |
| MassIVE | 10 |
| dbGaP | 5 |
| Dataverse / EMPIAR | 4 / 4 |
| MetaboLights / IDR | 2 / 1 |

The registry intentionally does not equate accession co-membership with identical samples,
processing, or measurements. The deterministic `ds:<hash>` is a split/leakage component, not a
scientific equivalence assertion.

## Official metadata resolution

A bounded 100-manifest pilot used official public APIs only. Candidates were prioritized before
non-candidates and non-GEO providers before GEO.

| API outcome | Count |
|---|---:|
| Resolved | 81 |
| Explicitly not found | 15 |
| No implemented official adapter | 4 |

Across all manifests, 71 titles, 59 file counts and 30 total-byte counts were exposed. Thirty
Figshare/Zenodo manifests supplied explicit open licence metadata; 1,771 manifests retain
`not_exposed`. Unknown fields were not inferred. The remaining 1,705 manifests are typed
`public_metadata_unresolved`; this is an acquisition backlog rather than evidence of inaccessibility.

## Modality-oriented acquisition queue

All 1,801 manifests appear in `acquisition_queue.csv`. Priority is deterministic: a manifest linked
to a gold candidate receives the largest increment, followed by target-modality coverage and source
count. Relevant manifest hints include:

- raw imaging 16;
- PIV 74;
- traction-force microscopy 65;
- AFM 250;
- immunofluorescence 1,001;
- western blot 1,041;
- PCR 1,359;
- sequencing 137 and broader omics 1,302;
- proteomics 113 and flow cytometry 715.

These are source/provider hints for review routing, not claims that each accession contains a
calibrated measurement of that modality.

## Conservative raw/source-data pilot

Four official, explicitly open datasets supplied nine selected files under the fixed 5,000,000-byte
cap. Total downloaded: **3,487,812 bytes**. Local SHA-256 audit: 9 files, 0 mismatches.

| Dataset | Files | Bytes | Licence | Purpose hint |
|---|---:|---:|---|---|
| Figshare `10.6084/m9.figshare.16826740` | 6 MAT | 31,953 | CC BY 4.0 | traction-force/actomyosin source data |
| Figshare `10.6084/m9.figshare.27241950` | 1 XLSX | 303,325 | CC BY 4.0 | nuclear-mechanics raw table |
| Figshare `10.6084/m9.figshare.28309337` | 1 CZI | 3,145,946 | CC BY 4.0 | raw fluorescence image |
| Zenodo `10.5281/zenodo.18779816` | 1 ZIP | 6,588 | GPL-3.0-or-later | PIV implementation artifact |

Payloads are content-addressed under ignored
`data/external_training/experiment_factory/datasets/pilot/**`. Only URL/licence/size/digest receipts
are committed. The remainder of each multi-gigabyte dataset was not downloaded.

## Reproducibility and limitations

- Focused tests: 7 passed.
- Registry generator replay: five generated artifacts byte-identical.
- Registry SHA-256: `c326a9c2c2cbfa2da0f14d20c3ee901c7571dc1221ca77c03793f32891bfe62c`.
- Dataset-manifest SHA-256 after metadata and pilot receipts:
  `d982dcb2f369533f04a5e19fa231e57935fd8ae7896fa687d1653461de818025`.
- Raw/source payload files tracked in Git: 0.

An exact identifier in a paper can describe generated data, reused data, a comparison dataset, or a
citation. Experiment-level annotation must adjudicate that relation before training. API resolution
only establishes public metadata reachability; it does not validate licence completeness,
measurement quality, biological truth, or fitness for Aleph parameter learning.
