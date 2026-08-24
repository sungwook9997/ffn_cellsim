# External-corpus retrieval bridge report

**Status:** implemented and tested; every output remains `proposed` and non-authoritative.

## Indexed corpus

- RAG–TAG–CAG source families bound by receipt: **7,213**
- OA-content model/local source families: **2,695**
- digest-verified, locator-addressable chunks: **139,557**
- committed article text: **0 bytes**
- compact index SHA-256: `ada80630836a929eee13027d17e510fc9eabb0a9618743ff665fb279ee7c029e`
- OA-content weights SHA-256: `f13932828a0b4d69099247f4e6db58b20388699864b07fc7ed359e963ca4c262`
- RAG ledger head: `a93c1912ed4973a1fba090eb0a37da9564e93bfd7210ea0e3fb9912e9da978ce`

The extraction run contained 139,573 total chunks. Six families did not have the bounded target
sections required by the frozen OA model; the bridge intentionally indexes the model/local
intersection, leaving **16** chunks outside neural routing rather than silently changing the model
catalogue.

## CPU measurement

One index/model initialization took **4,736 ms**. After initialization, seven in-process runs per
query produced:

| Query family | median | observed maximum | top weak domain |
|---|---:|---:|---|
| cortical tension | 502 ms | 551 ms | cortex_membrane_pressure |
| TFM/focal adhesion | 536 ms | 632 ms | adhesion_traction |
| IF/WB/PCR/cell state | 471 ms | 522 ms | mechanotransduction_signalling |
| PIV/cell migration state | 490 ms | 609 ms | protrusion_migration |

These are routing measurements, not model-quality or physics-validation claims. Exact deterministic
result hashes and top source families are recorded in `benchmark.json`; timing is kept outside query
results so repeated query JSON is byte-stable.

## Controls and boundary

The focused suite verifies deterministic rebuilding and querying, source-family uniqueness, maximum
280-character snippets, exact PMCID/section/chunk locators, three levels of digest validation, and
separate typed refusals for empty queries, physical-parameter estimation, training, authority claims,
and unknown vocabulary. It imports neither sealed Aleph learning/representation packages nor any
physics/scenario module.

This is a usable literature retrieval layer. It is not permission to train the sealed Aleph model,
estimate cortical tension, validate a physical law, or promote a source to evidence authority.

## Adversarial refusal follow-up

An initial action/target pattern missed `estimate cortical tension parameter 0.5 nN per um`. The
boundary now requires the semantic conjunction of an estimation action (`estimate`, `infer`, `fit`,
`calibrate`, `predict`, `calculate`, `derive`, or `determine`, including inflections) and a physical
value target. Explicit parameter, tension, modulus, viscosity, stiffness, pressure, force,
coefficient, and constant estimation requests return
`physical_parameter_estimation_refused`. Measurement-literature queries such as `methods measuring
cortical tension` remain routable because they do not request estimation.
