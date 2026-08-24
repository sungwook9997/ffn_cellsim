# Experiment annotation candidate set — measured report

## Outcome

The input consists of exactly 1,043 `manual_priority` OA sources: 536 from the first acquisition
pass and 507 from the second. The deterministic sampler selected exactly 300 unique papers. Every
record is `authority_status: proposed` and `review_state: gold_candidate`; zero records are described
as human verified or as Aleph physics evidence.

Candidate JSONL SHA-256:
`f4be334ea3fce6659a7395e77a2107a3b6a3dfdb5547f8b163e535946ff83226`.

## Coverage

| Dimension | Measured selection |
|---|---:|
| First / second acquisition pass | 160 / 140 |
| Publication years | every year 2010–2026 |
| Year strata 2010–13 / 2014–17 / 2018–21 / 2022–26 | 44 / 65 / 73 / 118 |
| Primary domains represented | 12 / 12 |
| Observed modality tags represented | 12 / 12 |
| Observed cell-state tags represented | 11 / 11 |
| Papers with exact payload SHA-256 | 300 / 300 |

Ten domains contain 26 candidates each. `organelle_general_mechanics` and
`osmotic_volume_poroelasticity` contain 20 each, because the complete eligible pools contain only 21
and 28 sources respectively. This predeclared allocation sums to 300 while retaining broad coverage;
it does not pad a scarce domain with unrelated papers.

The most sparsely represented modalities remain magnetic tweezers (10), optical tweezers (18),
micropipette aspiration (19), PIV (24), single-cell sequencing (27), and TFM (28). They are present
but require deliberate attention during experiment-level annotation.

## Leakage components

The combined public-accession graph was unioned before sampling, including connections mediated by
out-of-pool sources. Among the 1,043 eligible papers this produces 1,022 components. Five components
contain multiple eligible papers, totalling 26 papers. All five components and all 26 papers were
selected as whole units. Selected components split: **0**.

Twenty-five selected papers also have an accession-connected relation to at least one source outside
the `manual_priority` pool. That count is retained on each candidate. Accession co-membership is a
leakage precaution, not a claim that measurements are identical.

## Explicitly unmet cross-strata

Five cross-strata are absent and are recorded in `summary.json` rather than hidden:

- cortical/membrane-pressure × 2010–2013;
- organelle-general-mechanics × second pass;
- organelle-general-mechanics × 2010–2013, 2014–2017, and 2018–2021.

The organelle class is a first-pass discovery route concentrated in 2022–2026. The candidate set
does not relabel unrelated papers to manufacture balance.

## Annotation workflow

`annotation_assignments_600.csv` creates two unassigned, blinded slots per source: A and B. The file
has 600 assignment rows. A source can reach `adjudicated` only after two distinct annotators and an
adjudicator. The schema separately represents exact source and locator hashes, units and conversion
evidence, uncertainty type, biological and technical replicates, visual/table/curve/vector-field
observations, comparison ambiguity and the full review state.

The checked-in guidelines cover IF, confocal/live imaging, WB, PCR/qPCR, PIV, TFM, AFM,
micropipette/tweezers, tables and curves. No automatic extraction is accepted as a gold label.

## Reproducibility and controls

- Candidate semantic validation: 300 records, 0 errors.
- Focused tests: 10 passed.
- Positive controls: 1 valid proposed ExperimentRecord.
- Planted negative controls: invalid evidence hashes/offsets and ungrounded ambiguity; both killed.
- Generator replay: byte-identical candidate set, assignment table, summary and hash manifest.
- Source queue SHA-256 values and the combined leakage graph SHA-256 are pinned in `summary.json`.

This lane does not commit XML, figures, captions or publisher payloads. It writes only schemas,
metadata, sampling records and annotation instructions.
