# Combined OA-content router v3 — measured report

## Input and taxonomy

The technically usable corpus contains **5,669 unique source families**: 2,695 from pass one and all
2,974 from pass two. The six omitted first-pass files comprise one empty chunk file and five families
without any bounded abstract/methods/results/conclusion text. Both final quality indices cover their
respective extracted inputs; no record is promoted beyond proposed authority.

The 29 first-pass routes and 12 second-pass routes are not directly compatible. `taxonomy.json`
therefore maps every observed weak label exactly once into 12 explicit coarse mechanobiology or
measurement routes. The broadest route, `cell_state_mechanics`, has 1,304 families; the smallest,
`organelle_general_mechanics`, has 222. This is a retrieval taxonomy, not biological ground truth.

## Leakage-safe split

The final combined quality artifact contains 120 proposed shared-accession groups. Transitive union
produces 27 connected components spanning 140 usable families; the largest contains 71 families and
14 usable components contain more than one coarse route. Every component is assigned as a single
unit. There are **0 group split violations**.

| Split | Families |
|---|---:|
| Train | 3,980 |
| Validation | 855 |
| Test | 834 |

All 12 routes have nonzero support in every split. The smallest test support is 20
(`nuclear_mechanics`).

## Fair same-split result

Both CPU NumPy MLPs use the same records, group split, 8,192 hash buckets, 64-unit ReLU layer,
optimizer and early-stopping rule.

| Test metric | Title/journal/year | Combined OA text + metadata | Delta |
|---|---:|---:|---:|
| Macro-F1 | 0.5677 | **0.6518** | +0.0842 |
| Accuracy | 0.6091 | **0.6990** | +0.0899 |
| Retrieval hit-rate@5 | 0.7218 | **0.7638** | +0.0420 |
| Retrieval precision@5 | 0.5700 | **0.6904** | +0.1204 |

The final receipt-bearing run trained title-only for 10 epochs in 6.96 seconds and combined content
for 13 epochs in 9.37 seconds; the whole pipeline took 21.06 seconds on CPU. Full per-route F1 and
support are in `metrics.json`. The weakest content route is `organelle_general_mechanics` at F1
0.3600 (30 test records), so coarse aggregation did not eliminate weak routing.

Two complete runs produced identical title weights
`a0e8606465ab8102a5cd5309e3c912638552aad38bd5e547e93939611605ae52`, content weights
`ff2f0fb927ce484095bd881449b83ebf36b35e2414c7d7578f3718a173062bdb`, and split receipt
`fcc7afd0bec12ebf199bef39039d8f0b694e494465a2e1a1ac59ff63166b0a5d`.

## Honest comparison to v2

V3 doubles technically usable family coverage from 2,695 to 5,669 and raises minimum test support
from 3 to 20 by using 12 coarse routes rather than v2's 29 narrow labels. Consequently v2 and v3
macro-F1 values answer different classification questions and **must not be compared as an accuracy
improvement**. The valid v3 improvement claim is only the same-split title-only versus content
comparison above.

## Inference and refusal boundary

The frozen example routes the second-pass OA article `doi:10.1038/ncb3525` to
`cortical_membrane_pressure` with probability 0.9968 and retrieves distinct first- and second-pass
families. This remains weak literature routing. It cannot validate a paper, estimate a physical
parameter, train Aleph physics, promote evidence, unseal R5/MTG-PN, or make an authority claim.
