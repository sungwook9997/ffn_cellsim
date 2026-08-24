# Full OA recoverability screen — 2026-08-05

**Status:** AGENT-PROPOSED. This is a retrieval-priority report, not a quality verdict, Tier A
promotion, physical claim, or training authorization.

The stable extraction contained **2,701 / 2,701** successfully parsed OA JATS articles, **139,573**
chunks, 17,173 figure captions and 2,213 table captions. Its input manifest-set SHA-256 was
`cf618563667281631b8df4f5453a83c65402e95ae619ba8788c72297607ac069`.

## Decisions

| Proposed manual-routing decision | Articles |
|---|---:|
| `manual_priority` | 536 |
| `exception_candidate` | 74 |
| `hold` | 2,091 |

No article was called Tier A. A `manual_priority` article only has a dense combination of explicit
research-article metadata, sample-size, biological-replication and uncertainty signals, plus a
figure/table locator. A human must still decide applicability, adequacy, independence, corrections,
replication and biological meaning.

## Recoverability coverage

| Signal recovered with an exact locator | Articles | Coverage |
|---|---:|---:|
| units | 2,405 | 89.0% |
| uncertainty/error reporting | 1,958 | 72.5% |
| sample size | 1,918 | 71.0% |
| source/data availability | 1,324 | 49.0% |
| exclusions/outliers | 854 | 31.6% |
| biological replication | 610 | 22.6% |
| calibration | 594 | 22.0% |
| relevant figure/table observation locator | 2,260 | 83.7% |

The screen also found public dataset-accession signals in 319 articles and 61 accession groups used
by more than one source family. These are split constraints, not proof that two papers used the same
processed dataset.

## Critical uncertainty

- Affiliation-resolved laboratory groups and independent replication checks remain missing for all
  2,701 articles.
- Local correction/retraction relation metadata was unavailable for 2,578 articles; three more prior
  records were `unknown`. Five title/relation records were flagged for manual adjudication.
- Keyword and controlled-vocabulary signals can be false positives and do not establish adequacy.
- Full evidence records are ignored local data and contain locators plus hashes, not copied excerpts.

The compact committed index SHA-256 is
`0bd55f32829196ac65a13268988854aad6436dccdc0d79c8c8500f7b8a90e598`; the deterministic local
record-set SHA-256 is
`6b1ebd5895c29ec17dd321a3446144efdaa2ad6481bcb59d35fca7505bd1017f`.
