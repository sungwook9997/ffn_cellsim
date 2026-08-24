# OA-content v2 measured report

## Input audit

| Item | Measured value |
|---|---:|
| `*.chunks.jsonl` files | 2,701 |
| Empty files | 1 |
| Non-empty unique source families | 2,700 |
| Families with bounded abstract/methods/results/conclusion text | **2,695** |
| Abstract / methods / results / conclusion coverage | 2,689 / 2,231 / 2,056 / 1,306 |
| Weak domains | 29 |
| Train / validation / test | **1,878 / 418 / 399** |

All 29 weak domains have nonzero support in every split. Splitting is deterministic by canonical
source family before feature construction. Input chunk files have aggregate SHA-256
`d7fcb6ebbc2156aa85d49e087bbbf1f6a54ddc370c8b32c5f5535403ae17c9d2`.

## Fair same-split comparison

Both models use the same 2,695 records, split, 8,192 hash buckets, 64-unit ReLU hidden layer,
optimizer, stopping rule and label vocabulary. The only comparison difference is bounded OA section
tokens.

| Test metric | Title/journal/year | OA content + metadata | Delta |
|---|---:|---:|---:|
| Macro-F1 | 0.5421 | **0.6349** | +0.0928 |
| Accuracy | 0.5965 | **0.6817** | +0.0852 |
| Retrieval hit-rate@5 | 0.7043 | **0.7920** | +0.0877 |
| Retrieval precision@5 | 0.5784 | **0.6807** | +0.1023 |

In the final receipt-bearing run, title-only trained for 18 epochs in 16.05 seconds; OA content
trained for 19 epochs in 19.75 seconds
on CPU. The full per-domain F1 and support table is in `metrics.json`. Low-support results remain
visible: for example `cancer_invasion_mechanics` has only 3 test records, and no per-domain result is
an authority claim.

Two consecutive complete runs produced identical artifacts: title-only weights
`7c5e8c30e127af8d27feb7cb39d3ab9dab2288cce00d698ae56b7726b0d935ba`, OA-content weights
`f13932828a0b4d69099247f4e6db58b20388699864b07fc7ed359e963ca4c262`, and split receipt
`c620b571c8e638cd74268c3938872a2562ef1c266a80ca3282337f0325501876`.

## Boundary

The improvement shows that bounded JATS text helps reproduce weak search-domain routing. It does not
show that the model learned mechanics, that any paper is reliable, or that a physical parameter can
be trained. All records remain proposed external evidence; MTG-PN and R5 remain sealed.

The frozen inference example uses the OA article `doi:10.1038/s41556-025-01807-6`, predicts its weak
`cortex_membrane_pressure` route at 0.9694, and retrieves five distinct non-self source families.
`example_inference.json` is machine-readable and checked against the live frozen index by a control.
