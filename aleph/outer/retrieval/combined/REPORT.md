# Combined external-corpus retrieval report

**Status:** implemented and locally validated. Every source, domain, chunk, and snippet remains
`proposed`; this bridge grants no training or scientific authority.

## Bound coverage

| Layer | Combined bridge | v1 bridge | Scope comparison only |
|---|---:|---:|---:|
| RAG–TAG–CAG families | 10,208 | 7,213 | +2,995 |
| local OA families indexed | 5,675 | 2,695 | +2,980 |
| neural-routable OA families | 5,669 | 2,695 | +2,974 |
| exact-locator chunks indexed | 284,620 | 139,557 | +145,063 |

The 145,063 chunk increase is 145,047 second-pass chunks plus 16 first-pass chunks that v1 omitted
because their six families lacked bounded v2 target-section content. The combined text-free index
records all 284,620 chunks. Neural results use only the 5,669-family v3 model/local intersection.

- index SHA-256: `ede8d2aff2a8786efa8e01bba6fbb398f46bdfdcc34eda0ef2662bdefe7a5887`
- v3 weights SHA-256: `ff2f0fb927ce484095bd881449b83ebf36b35e2414c7d7578f3718a173062bdb`
- taxonomy SHA-256: `2b0d69827c88eee5dcd8f4cc8956a05eb33dd93ed26362e7ddf4ccc98a2ae403`
- combined RAG ledger head: `662cea6bcf553d441b622de5a9c3d34745453f3586b8b187b7878728ea9e146d`
- combined RAG projection: `15ff31f6719ed3df70a71556eb00a695eedb94d654f84dc326a178162aca963a`
- committed article text: **0 bytes**

## CPU latency

Fresh-process import plus initialization and seven warm calls per query were measured separately for
combined and v1. These measurements compare operational cost only; v3's 12 coarse routes and v1's
29 weak labels are incompatible quality tasks.

| Query | Combined median / max | v1 median / max |
|---|---:|---:|
| cortical tension | 392 / 464 ms | 311 / 331 ms |
| TFM | 331 / 388 ms | 292 / 313 ms |
| IF/WB/PCR | 363 / 434 ms | 310 / 343 ms |
| PIV/cell state | 280 / 339 ms | 317 / 342 ms |

Cold import plus initialization was **5,113 ms combined** and **2,564 ms v1**. Exact deterministic
combined result hashes are recorded in `benchmark.json`. Timing is deliberately absent from query
JSON, preserving byte-stable results.

## Safety and exactness

Before emitting a maximum-280-character snippet, the bridge verifies the indexed file SHA-256, text
SHA-256, and canonical chunk-record SHA-256. Results carry source family, pass, PMCID, proposed
coarse domain, original weak label, split, exact section locator and chunk ordinal/hash.

The hardened v1 semantic boundary is reused unchanged. Empty/unknown queries, training requests,
authority claims, and estimation/inference/fitting/calibration/prediction of physical values return
typed refusals. Measurement-literature searches remain permitted. Focused controls cover exact
counts, deterministic index regeneration and query output, both-pass retrieval, bounded snippets,
the requested four query families, and adversarial refusal/allow cases.
