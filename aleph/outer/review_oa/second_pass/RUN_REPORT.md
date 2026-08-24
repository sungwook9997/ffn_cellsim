# Second-pass and combined OA recoverability screen — 2026-08-05

**Status:** AGENT-PROPOSED. This report ranks manual review. It is not Tier A promotion, a quality
verdict, a physical claim, or training authorization.

The independently selected second pass produced **2,974 / 2,974** stable JATS article manifests,
145,047 chunks, 21,143 figure captions and 2,237 table captions. The extraction replay agreed after
excluding only its expected resume counter. The manifest-set SHA-256 is
`9627db89509f0a2f363ec93e65a2da4b68d3ab1b0415195788145b54a30bcf0f`.

## Second-pass routing

| Proposed decision | Articles |
|---|---:|
| `manual_priority` | 507 |
| `exception_candidate` | 0 |
| `hold` | 2,467 |
| automatic Tier A | 0 |

All 2,974 candidate metadata records explicitly identify research articles, so no method or rare-
modality exception route was needed. A `manual_priority` value still means only that the paper has a
dense set of recoverable signals and exact locators.

## Recoverability

| Signal with locator | Articles | Coverage |
|---|---:|---:|
| units | 2,865 | 96.3% |
| uncertainty/error reporting | 2,599 | 87.4% |
| sample size | 2,357 | 79.3% |
| source/data availability | 1,080 | 36.3% |
| exclusions/outliers | 908 | 30.5% |
| calibration | 624 | 21.0% |
| biological replication | 588 | 19.8% |
| relevant figure/table observation locator | 2,803 | 94.3% |

Correction/retraction relation metadata was unavailable for 2,966 articles. Four had no-link status
from the earlier dated corpus check and four were flagged by title metadata for manual adjudication.
Independent-replication checks and affiliation-resolved laboratory groups remain missing for every
article.

## Combined first and second passes

- **5,675** distinct source families; cross-pass source-family overlap: **0**.
- Routing: 1,043 `manual_priority`, 74 exception candidates, 4,558 holds, zero automatic Tier A.
- Dataset-accession signals: 214 second-pass articles.
- Duplicate accession groups: 18 within the second pass, **120** across the combined corpus.
- **70** combined groups contain members from both passes. Group categories overlap: a cross-pass
  group may also have multiple members within one pass.

Accession co-membership is enforced as a leakage warning and split constraint. It does not prove
identical preprocessing or experimental independence.

Full second-pass records remain in ignored local storage. The committed compact-index SHA-256 is
`d9f494d552114f832130f4e8b90cfc714d7fe3641ab872a01926a854f96cbaa0`; the combined-leakage
artifact SHA-256 is
`16442a2b92b17fc0560fe7ed774f72d1a79ecde13cb3f50f85d759ad0152868c`.
