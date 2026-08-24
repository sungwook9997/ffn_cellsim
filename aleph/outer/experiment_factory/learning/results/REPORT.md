# Learning projection report

Status: **PROPOSED projection, not a trained model, not adjudicated gold, and not physical authority**

## Coverage

| View | Sources present |
|---|---:|
| Text-free metadata | 5,675 |
| Conservative numeric candidates | 4,876 |
| Visual caption/index metadata | 5,635 |
| Acquired visual-asset receipts | 116 |
| Actual CPU visual descriptors | 116 |
| Public dataset manifests | 885 |

The visual distinction is intentional. There are 38,316 figure metadata rows and
627 described receipt rows over 480 unique assets, but actual descriptor coverage is
only 116 source families. Caption/index presence never sets the asset or descriptor
mask. Described sources carry deterministic source-level mean/std/min/max aggregates
for image dimensions, intensity, colorfulness, entropy, edge density, quantiles and
panel count, plus stable descriptor/asset/observation IDs and descriptor digests.

Numeric projection contains 238,489 structured observations. Only 39,344 have an
unmasked proposed-candidate target; 199,145 are masked because the observation or
its local biological/protocol context is ambiguous. **Zero ambiguous observations
are promoted to labels.** Local numeric rows contain structured values, units,
stable IDs and evidence digests but no evidence text.

Dataset coverage comprises 1,982 source-accession bindings across 885 source
families. Public-accession components and identical payload SHA-256 groups form the
leakage graph before splitting.

## Split and audit

| Split | Sources |
|---|---:|
| Train | 4,004 |
| Validation | 840 |
| Test | 831 |

- Sources: 5,675
- Leakage groups: 5,560
- Multi-source groups: 30; largest group: 69 sources
- Multi-source public-accession components: 123
- Frozen 300 evaluation cohort: exactly 300 sources, marked not-ground-truth
- Cross-split leakage-group violations: 0
- Duplicate projection/source IDs: 0
- Forbidden raw text/caption/pixel/publisher fields: 0
- Authority promotions: 0

The committed projection contains source-level text-free aggregates and stable join
references. Full per-asset descriptor rows remain in the visual producer; this
projection binds them by descriptor ID/SHA-256 and carries only source-level numeric
aggregates rather than copying pixels. Raw evidence,
captions, graphic hrefs, publisher payloads and image bytes are absent.

Final projection SHA-256 is
`929adbbbd5b5ee7c59cabd5612c3961bd8e4a93edcdaafeb3e3eb654a7f42421`.
Exact input hashes are recorded in `summary.json`; independent structural results
are in `audit.json`.
