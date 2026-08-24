# Second-pass OA JATS extraction — 2026-08-05

Status: **AGENT-PROPOSED**.  These are structural extraction and lexical-tag
coverage measurements, not quantitative validation or Aleph claims.

The independently discovered second pass contained 3,000 candidates.  Final
acquisition produced 2,974 successful OA XML objects and 26 HTTP 404 receipts.
The successful set contains 2,974 unique source families and 2,974 unique
payloads, with zero source-family overlap against the 2,701 successful first
pass objects.

| Measure | Count |
|---|---:|
| Successful XML families parsed | 2,974 |
| Validated payload bytes | 404,754,837 |
| Retrieval chunks | 145,047 |
| Sections including abstracts | 82,918 |
| Figures located | 21,143 |
| Tables located | 2,237 |
| Extraction failures | 0 |

| Proposed lexical tag dimension | Articles | Coverage |
|---|---:|---:|
| Cell type/state | 2,528 | 85.00% |
| Measurement modality | 2,348 | 78.95% |
| Mechanics observable | 1,793 | 60.29% |
| Perturbation | 2,574 | 86.55% |

The initial extraction and all-resume replay produced the identical canonical
article-manifest set SHA-256:
`9627db89509f0a2f363ec93e65a2da4b68d3ab1b0415195788145b54a30bcf0f`.
The independent local audit checked 2,974 manifest files, 2,974 chunk files and
all 145,047 chunk records: zero file-hash, text-hash, chunk-hash, count,
orphan-file, or authority errors.

The text-bearing derived store remains ignored at
`data/external_training/derived_second_pass/` and occupied 244 MB at audit time.
No full text, long excerpt, figure, table, or publisher payload is committed.
