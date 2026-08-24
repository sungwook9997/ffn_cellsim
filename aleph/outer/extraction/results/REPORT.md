# OA JATS extraction result — 2026-08-05

Status: **AGENT-PROPOSED**.  This is structural extraction coverage, not
quantitative validity, a literature finding, or an Aleph claim.

The final frozen acquisition input contained 2,713 receipts: 2,701 successful
OA JATS XML objects and 12 acquisition failures.  Every successful object was
SHA-256 verified and parsed; extraction failures were zero.

| Measure | Count |
|---|---:|
| XML articles parsed | 2,701 |
| Validated payload bytes | 445,899,949 |
| Retrieval chunks | 139,573 |
| Sections including abstracts | 79,071 |
| Figures located | 17,173 |
| Tables located | 2,213 |
| Extraction failures | 0 |

Tag coverage counts articles with at least one lexical match in that proposed
dimension.  It does not say that the tag is a confirmed property of the study.

| Proposed tag dimension | Articles | Coverage |
|---|---:|---:|
| Cell type/state | 2,496 | 92.41% |
| Measurement modality | 2,118 | 78.42% |
| Mechanics observable | 1,564 | 57.90% |
| Perturbation | 2,425 | 89.78% |

The canonical digest-only article-manifest set SHA-256 was
`cf618563667281631b8df4f5453a83c65402e95ae619ba8788c72297607ac069`
on both final passes.  The first final pass parsed 1,631 new objects and resumed
1,070; the repeat resumed all 2,701.  A complete local audit independently
recomputed all 139,573 text hashes and chunk hashes: zero mismatches and zero
records above `proposed` authority.

Text-bearing records remain only under ignored
`data/external_training/derived/articles/` (2,701 manifest files and 2,701 chunk
files; 239 MB total derived footprint at closeout).  No article text, long
excerpt, figure, table, or publisher payload is committed here.
