# Combined first- and second-pass RAG–TAG–CAG validation

**Authority:** every source remains `proposed`. Importing a paper or obtaining
its XML does not promote an Aleph claim.

The original six manifests contain 7,213 normalized source families. The 3,000
second-pass candidates overlap those families five times, giving **2,995 new
families** and **10,208 combined families**. The combined seven manifests hold
10,218 rows, of which ten are duplicate rows after case-normalized source-family
deduplication.

## Acquisition and object audit

| Item | Count |
|---|---:|
| First-pass receipts | 2,713 |
| Second-pass receipts | 3,000 |
| Successful XML payloads | 5,675 |
| Explicit failed receipts | 38 |
| Validated XML bytes | 850,654,786 |
| Content-addressed objects | 21,596 |
| Missing objects | 0 |
| Corrupt objects | 0 |

Every receipt is itself an immutable object. Each successful receipt references
payload bytes that were copied into the fresh external DataRoot and re-hashed.
Failed receipts remain referenced objects with their failure status and have no
fabricated payload.

## RAG–TAG–CAG result

| Check | Result |
|---|---:|
| First-run RAG events | 10,208 |
| TAG entities | 10,208 |
| Explicitly retracted entities | 21 |
| CAG lint violations | 0 |
| Second-run new events | 0 |
| Ledger chain | verified |

The idempotent second run reproduced ledger head
`662cea6bcf553d441b622de5a9c3d34745453f3586b8b187b7878728ea9e146d`,
TAG manifest
`f2de469010ce8de4bd234767d08fcbfb7c3c9676c1cbabfda28f2a65c9e9f203`,
and semantic projection
`15ff31f6719ed3df70a71556eb00a695eedb94d654f84dc326a178162aca963a`.
The ledger, snapshot and 850.7 MB of source payload bytes remain outside Git.
