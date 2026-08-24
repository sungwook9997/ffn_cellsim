# Full-text/table extraction report

Status: **PROPOSED automatic candidates; not human-verified gold and not physical authority**

## Final measured result

| Measure | Result |
|---|---:|
| Validated OA JATS articles processed | 5,675 / 5,675 |
| XML/payload/hash failures during rebuild | 0 |
| Articles with at least one conservative candidate | 4,876 (85.92%) |
| ExperimentRecord candidates | 52,725 |
| Explicit numeric observations | 238,489 |
| Exact or deterministic unit normalizations | 238,305 (99.923%) |
| Ambiguous records | 44,385 (84.18%) |
| Observations with reported uncertainty | 6,112 |
| Observations with sample-size/replicate evidence | 44,905 |
| Table / figure-caption records | 1,664 / 14,302 |
| Canonical-schema validation | 52,725 valid, 0 errors |
| Unique experiment / observation IDs | 52,725 / 238,489 |
| Full output hash-audit errors | 0 |

The high ambiguity rate is deliberate: a paragraph or caption frequently reports a
number without repeating its cell type, exact protocol, error semantics, or
comparison arms.  The extractor flags those gaps instead of filling them from a
distant paragraph.

Of the observations, 113,039 are protocol parameters and 125,450 are reported
measurement candidates.  These categories remain candidates until annotation.
No unresolved comparison arm was emitted as a structured comparison; all 52,725
records currently require downstream arm-resolution if a contrast is desired.

## False-positive correction audit

The first broad run produced 63,716 records, 268,651 observations, and 54,119
ambiguous records.  Restricting eligible evidence classes and hardening unit parsing
reduced that to 52,725 / 238,489 / 44,385 respectively: **−10,991 records,
−30,162 observations, and −9,734 ambiguous records**.

- `18 S rRNA` and residue names such as `W748S`: uppercase `S` is never seconds.
  Final reported-unit count for `S` is 0.
- sample notation: lowercase `n` is never newtons; final reported-unit count for
  `n` is 0. Bare uppercase `N` is accepted only with explicit force/load/newton/
  tension/thrust/traction context and is rejected for `N/m` and sample prose.
- `nM` is case-sensitively concentration (3,535 observations); `nm` is length
  (7,356 observations in the pre-final unit audit; the final rule is unchanged).
- lowercase `µm`/`μm` followed by a compound-like drug name is not silently treated
  as length: 184 observations are `length_or_concentration` with explicit ambiguity.
- ordinary Background/Introduction/Discussion mortality and literature-summary
  numbers are ineligible for automatic ExperimentRecords.

Nine focused controls cover exact conversion, no-unit refusal, background mortality,
`18 S`, `nM/nm/N/n`, lowercase-micrometre drug ambiguity, table locators, resume
determinism, and repeated-text ID uniqueness.

## Replay and integrity

A second execution resumed all 5,675 article manifests and reproduced manifest SHA
`05219037aa0fa3674ddffd44080b8bd10b706029997e488af6c6e2cdfdfb6e04`.
The final audit checked every record, observation, local evidence row, file digest,
and evidence binding.  Article-index SHA is
`bbc4bc81ee18a7b8d8fecdeed1f11dfb683131ac5bd6f6e6fb7c93a90ae9dc02`.

Full text and local candidates are ignored data.  No article text, JATS payload, or
scientific value promotion is committed.  The figure join index uses
`(source_family_id, figure_locator_id, document_sha256, caption_text_sha256)` for the
visual lane without embedding captions or assets.
