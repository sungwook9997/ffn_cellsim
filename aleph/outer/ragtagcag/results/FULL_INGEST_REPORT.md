# Full external-corpus RAG–TAG–CAG validation

**Status:** `proposed`; this run imports literature evidence and promotes no
Aleph claim.

The six current manifests contain 7,218 rows. Normalising and deduplicating
their `source_family_id` values produces **7,213** source families and five
duplicate records. This is one fewer than the earlier conversational count of
7,214; the authoritative files have five seed/snapshot DOI overlaps, not four.
No placeholder source was invented to make the earlier count hold.

After the acquisition worker finished all 2,713 eligible PMCIDs, the final run
used an immutable 2,713-receipt copy with SHA-256
`2fdd7cd76c139d2d5619a42aab7b263cea2e89aae48ed147d9cee31cf6496d8c`.
It contained 2,701 successful XML acquisitions (445,899,949 validated bytes)
and 12 explicit failures. All 2,701 successful payload digests were found,
copied into the isolated Aleph ObjectStore, and re-hashed successfully. The
failed receipts remain provenance records and are not represented as payloads.

## Result

| Check | Result |
|---|---:|
| First-run RAG events | 7,213 |
| TAG entities | 7,213 |
| Acquisition receipts | 2,713 |
| Acquired XML payloads present | 2,701 |
| Content-addressed objects | 12,627 |
| Explicitly retracted entities | 8 |
| Missing snapshot objects | 0 |
| Broken object digests | 0 |
| CAG lint violations | 0 |
| Second-run new events | 0 |
| Ledger chain | verified |

The second run reproduced the same ledger head, TAG manifest hash, and semantic
projection hash. The full SQLite ledger, object bytes, and snapshot directory
remain outside Git; only this compact receipt is versioned.
