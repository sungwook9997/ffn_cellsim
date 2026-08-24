# External training corpus quality review — first batch

**Status:** AGENT-PROPOSED. No record in this review is an Aleph claim or training authorization.

## Counts

- Candidate-rich journals attempted: **50**
- JIF profiles actually verified in licensed Clarivate JCR: **15**
- Verified JIF >= 5: **8**
- JIF access failures left `UNKNOWN`: **35**
- Tier A candidates screened at article level: **120**
- Article records reviewed: **137**
- Tier A passed: **0**
- Holds: **135**
- Rejects: **0**
- Explicit exceptions: **2**
- Article full-text access failures or metadata-only reviews: **18**

## Method and limits

The JIF batch records the 2025 metric year, value, licensed JCR profile URL, and verification date. It does not substitute CiteScore, SJR, or inferred values. Chrome institutional access disconnected after 15 successful profile checks; the remaining 35 values are `UNKNOWN` rather than guessed.

Article screening is a conservative first pass over at least 100 research-article candidates from the candidate-rich top-50 journals, plus explicit seed exceptions. Recoverability signals come from Europe PMC full-text XML; correction/retraction links come from Europe PMC CORE `commentCorrectionList`. A keyword hit records a short source excerpt and section locator but does not itself prove methodological adequacy. Missing sample size, biological independence, calibration, units, uncertainty, exclusion, or source-data evidence causes a hold. Citation-level independent replication and exact figure/table observation locators remain `UNKNOWN` in batch 1, so JIF alone never produces Tier A.

## Leakage controls

Every paper, supplement, repository deposit, and replot remains in one `source_family_leakage_group` and is not counted as independent evidence. `lab_leakage_group` is only a last-author-surname proxy and therefore cannot authorize a split until affiliations are manually resolved. Extracted accession identifiers form dataset leakage groups; absent identifiers remain `UNKNOWN`.

## Next safe action

Restore Chrome extension connectivity and verify the remaining 35 top-journal JCR profiles. Then manually adjudicate exact figure/table locators, field applicability, affiliation-resolved laboratory groups, dataset identities, and citation-level replication/contradiction before any hold can become Tier A. All records remain `proposed`.
