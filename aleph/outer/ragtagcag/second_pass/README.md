# Combined second-pass RAG–TAG–CAG validation

This directory contains only the deterministic wrapper, controls and compact
receipt for combining the original six manifests with the second-pass candidate
manifest and both OA acquisition receipt sets.

The actual ledger, immutable objects and TAG snapshot must use an Aleph
`DataRoot` outside the repository. Full text is never committed. Every source
remains `proposed`; failed acquisition receipts and explicit retractions remain
visible rather than being dropped.
