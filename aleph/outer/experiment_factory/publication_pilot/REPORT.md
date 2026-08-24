# Publication-validation pilot report

Status: **implemented and frozen; human annotation not yet performed**

## Measured construction

| Item | Result |
|---|---:|
| parent candidate cohort | 300 papers |
| frozen pilot | 50 papers |
| A/B assignment slots | 100 |
| primary domains represented | 12 / 12 |
| observed modality tags represented | 12 / 12 |
| acquisition passes | first 27 / second 23 |
| year strata 2010–13 / 2014–17 / 2018–21 / 2022–26 | 7 / 11 / 12 / 20 |
| machine ExperimentRecord candidates exposed separately | 996 |
| split leakage components | 0 |
| completed human A / B / adjudication submissions | 0 / 0 / 0 |

Every domain has four papers except cortical/membrane-pressure and measurement/inference, which have
five. All twelve modality tags are present; the sparsest is magnetic tweezers with one paper. One
two-paper public-accession component is included whole. The 17-paper component is deferred rather than
partially sampled because including it would consume 34% of this pilot.

## Implemented gates

- Machine pre-annotations and blinded forms are separate artifacts.
- A/B forms contain no prediction, annotator identity, claim or adjudication result.
- The evaluator rejects the same A/B annotator and mismatched source sets.
- Adjudication preparation reports agreement and disagreement but chooses no claim.
- Machine evaluation accepts only explicit human-adjudicated records and claims.
- Ambiguous machine claims are masked from exact field metrics.
- Full-text, caption and pixel fields are absent from committed packets.
- The complete generator replays byte-for-byte, including deterministic gzip.
- The integrated factory and publication-pilot suite passes 58 focused tests.
- The committed pilot audit reports PASS with eight manifest files, zero digest mismatch, zero forbidden
  raw fields and zero split leakage components.

## Honest completion boundary

The software and frozen sampling protocol are complete. The scientific benchmark is not: no human has
filled either annotation slot or adjudicated a source. Consequently no inter-annotator agreement,
precision, recall, F1, failure prevalence or gold-label count is reported. `benchmark_readiness.json`
records this absence as data rather than allowing an empty form to look like a negative annotation.

This pilot is the quality-control gate for the broader externally supervised cell-state and observable
models. It is not itself that model and does not resolve the need for public raw datasets, author labels,
independent laboratory holdouts or prospective biological validation.
