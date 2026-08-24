# Europe PMC OA XML acquisition — 2026-08-05

Status: measured acquisition receipt; no claim promotion

## Scope and result

Five immutable discovery snapshots contained **2,713 unique PMCIDs** whose
records explicitly had `is_open_access_provider_flag: true`. Only those PMCIDs
were sent to Europe PMC's official `fullTextXML` REST endpoint.

| Measure | Result |
|---|---:|
| Eligible unique PMCIDs | 2,713 |
| Validated XML successes | 2,701 |
| Explicit failures | 12 |
| Success rate | 99.56% |
| Unique SHA-256 payloads | 2,701 |
| Duplicate payloads | 0 |
| Validated XML bytes | 445,899,949 |
| Acquisition span | 1,857 s (30 min 57 s) |
| Successful payload rate | 1.454496/s |
| XML with embedded licence text | 2,689 |

The payloads are ignored local objects at
`data/external_training/objects/<sha256>`. Git contains only this report, the
collector, the summary, and one receipt per candidate.

## Licence observations

Provider metadata among the 2,701 successful payloads was: CC BY 1,910; CC
BY-NC-ND 499; CC BY-NC 185; CC BY-NC-SA 47; CC0 6; provider licence unknown
54. The receipt also records licence elements extracted from each XML. Absence
of a recognised provider string does not grant redistribution; downstream use
remains subject to the recorded source licence.

## Failures

All 12 failures were HTTP 404 responses from the official full-text endpoint,
including a second explicit retry. This is recorded as a provider OA-flag / XML
availability mismatch, not as successful acquisition:

- `PMC12860937`, DOI `10.1016/j.jbc.2025.111126`
- `PMC3492504`, DOI `10.1038/nnano.2012.163`
- `PMC3981899`, DOI `10.1038/nm.3497`
- `PMC4189826`, DOI `10.1038/nphoton.2014.165`
- `PMC4452027`, DOI `10.1038/ncb3157`
- `PMC5857237`, DOI `10.1038/s41593-018-0083-7`
- `PMC6386196`, DOI `10.1016/j.celrep.2019.01.023`
- `PMC6486410`, DOI `10.1038/s41593-019-0369-4`
- `PMC6695368`, DOI `10.1007/s11912-019-0839-6`
- `PMC6788758`, DOI `10.1038/s41593-018-0316-9`
- `PMC6858575`, DOI `10.1038/s41588-019-0514-8`
- `PMC6899165`, DOI `10.1038/s41551-019-0420-5`

No publisher fallback, institutional download, CAPTCHA interaction, or non-OA
request was attempted.

## Verification

- All 2,713 receipts have unique PMCIDs and exactly match the eligible input
  set.
- Every successful object's full bytes were re-hashed and matched both its
  filename and receipt SHA-256.
- All 2,701 successful objects were parsed again as XML and had an `article` or
  `pmc-articleset` root.
- Object-store filenames exactly equalled the 2,701 successful receipt hashes;
  there were no missing or orphan objects after closeout.
- A resumability dry run reported 2,713 existing receipts and **0 pending
  attempts**.
- `bootstrap/scripts/bootstrap_check.py`: all checks passed.
- `python -m unittest discover -s bootstrap/tests -v`: 8 tests passed.

These checks establish retrieval integrity and provenance only. They do not
establish article quality, biological truth, or eligibility for model training.
