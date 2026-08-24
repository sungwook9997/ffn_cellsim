# Targeted Europe PMC OA second pass — 2026-08-05

Status: measured discovery and acquisition receipt; no claim promotion

## Result

The second pass used Europe PMC's official search and `fullTextXML` APIs. The
search was constrained to `SRC:MED`, provider open access,
`PUB_TYPE:"research-article"`, PMCID-bearing records, publication dates from
2010 through 2026-08-05, and twelve mechanobiology, cell-state, or measurement
query families. It excluded all 2,713 first-pass PMCIDs and source families.

| Measure | Second pass | Combined first + second pass |
|---|---:|---:|
| Unique candidate receipts | 3,000 | 5,713 |
| Validated XML successes | 2,974 | 5,675 |
| Explicit failures | 26 | 38 |
| Unique SHA-256 payloads | 2,974 | 5,675 |
| Validated XML bytes | 404,754,837 | 850,654,786 |
| Payload duplicates | 0 | 0 |
| Local content store | — | 822 MiB on disk |

The second-pass acquisition span, including the explicit failure retry, was
1,872 seconds. The final measured success rate was **1.588675 payloads/s**.

## Relevance strata

The global cap was filled without broad fallback terms or padding. Sparse
strata remained sparse; balanced round-robin selection prevented early domains
from consuming the cap.

| Domain | Candidates |
|---|---:|
| adhesion / traction / ECM | 268 |
| cell-type and state mechanics | 268 |
| cortex / membrane / pressure | 235 |
| cytoskeleton / motor / rheology | 268 |
| division / morphogenesis / tissue | 262 |
| flow / shear / PIV | 238 |
| mechanical measurement | 268 |
| mechanotransduction signalling | 268 |
| migration / protrusion / invasion | 265 |
| multimodal observation / inference | 264 |
| nucleus / chromatin mechanics | 132 |
| osmotic volume / poroelasticity | 264 |

| Publication stratum | Candidates |
|---|---:|
| 2010–2013 | 670 |
| 2014–2017 | 756 |
| 2018–2021 | 774 |
| 2022–2026 | 800 |

The smallest underlying stratum was nucleus/chromatin mechanics in 2010–2013:
only four new eligible papers were available under the declared query. It was
not filled with unrelated records.

## Licence observations

Provider licence metadata for 2,974 successful XML payloads was: CC BY 2,104;
CC BY-NC-ND 379; CC BY-NC 173; CC BY-NC-SA 153; CC0 17; unknown 148. Embedded
licence text was present in 2,945 XML payloads. Unknown metadata never implies
permission; all downstream redistribution remains subject to the source
licence recorded in the per-paper receipt.

## Failures

All 26 failures returned HTTP 404 from the official full-text endpoint on the
initial request and explicit retry. No publisher, subscription, or non-OA
fallback was attempted. The exact failure records are:

- `PMC12295664` — DOI `10.3390/ijms26146705`
- `PMC2919753` — DOI `10.1038/nmat2732`
- `PMC3150339` — DOI `10.1038/nsmb.2084`
- `PMC3297676` — DOI `10.1038/nature10801`
- `PMC3310973` — DOI `10.1038/onc.2011.593`
- `PMC3536889` — DOI `10.1038/nature11693`
- `PMC3615085` — DOI `10.1038/pr.2013.3`
- `PMC3755030` — DOI `10.1038/jid.2013.184`
- `PMC3756671` — DOI `10.1038/ncomms3240`
- `PMC3777337` — DOI `10.1038/ncb2614`
- `PMC4662888` — DOI `10.1038/ncb3268`
- `PMC4666809` — DOI `10.1038/nmeth.3616`
- `PMC4720436` — DOI `10.1038/nature14215`
- `PMC4996707` — DOI `10.1038/nmat4654`
- `PMC5334365` — DOI `10.1038/nature21407`
- `PMC6344062` — DOI `10.1038/nphys4219`
- `PMC6345402` — DOI `10.1016/j.celrep.2018.10.101`
- `PMC7021530` — DOI `10.1016/j.celrep.2019.12.040`
- `PMC7210009` — DOI `10.1038/s41586-020-1998-1`
- `PMC7275893` — DOI `10.1038/s41592-020-0818-8`
- `PMC7680637` — DOI `10.3791/61433`
- `PMC7792532` — DOI `10.1016/j.celrep.2020.108409`
- `PMC7997775` — DOI `10.1016/j.celrep.2021.108816`
- `PMC8168789` — DOI `10.1016/j.celrep.2020.108117`
- `PMC8594876` — DOI `10.20517/jca.2021.25`
- `PMC9128075` — DOI `10.1101/sqb.2019.84.040360`

## Verification

- All 3,000 second-pass candidates are distinct by both PMCID and source
  family, and both intersections with the 2,713 first-pass records are empty.
- Every candidate is provider-OA, in Europe PMC, dated 2010–2026, and carries
  the `research-article` publication type.
- All 5,675 combined successful objects were re-read, SHA-256 hashed, size
  checked, and parsed as article XML.
- The 5,675 object filenames exactly equal the combined successful receipt
  digest set: no missing objects, orphan objects, or duplicate payloads.
- The combined byte total is exactly 850,654,786.
- All 38 combined failures have zero bytes, a null digest, and HTTP 404.
- A second-pass resume dry run reported 3,000 existing receipts and zero
  pending requests.

These checks prove retrieval integrity, access routing, and corpus separation.
They do not prove article quality, biological truth, or training eligibility.
