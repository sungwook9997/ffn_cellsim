# ExperimentRecord publication-validation pilot

This directory turns the frozen 300-paper annotation-candidate cohort into a deterministic 50-paper
human-validation pilot. It is a quality gate for later biological and multimodal learning, not the
final Outer Library and not a substitute for human review.

## Frozen pilot

The 50 papers cover all 12 primary domains and all 12 observed modality tags. The four publication
strata contain 7 / 11 / 12 / 20 papers and the two acquisition passes contain 27 / 23. Domain counts
are four each except cortical/membrane-pressure and measurement/inference, which contain five each.

Public-accession components are indivisible. One two-paper component is included whole. A 17-paper
component is explicitly deferred because it would consume 34% of the pilot; no component is split.

## Three separated artifacts

1. `forms_A.jsonl` and `forms_B.jsonl` are blank, blinded forms. They contain no machine prediction,
   other-annotator identity, or adjudication result.
2. `machine_preannotations.jsonl.gz` contains proposed structured candidates and exact evidence keys.
   It is hidden during independent annotation and is never called gold.
3. `adjudication_forms.jsonl` begins empty and pending. The comparison CLI can prepare disagreement
   packets, but it cannot choose a winning claim or fill an adjudicator identity.

`results/benchmark_readiness.json` deliberately reports that no accuracy or agreement metric exists:
zero human A submissions, zero human B submissions and zero adjudications have been completed.

## Rebuild

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python -m \
  corpus.external_training.experiment_factory.publication_pilot.build_pilot \
  --candidates corpus/external_training/experiment_factory/goldset/candidates_300.jsonl \
  --records-root data/external_training/experiment_factory/extraction/full/records \
  --output corpus/external_training/experiment_factory/publication_pilot/results
```

The committed machine packet is gzip-compressed with a zero timestamp. A clean rebuild must match
every hash in `results/SHA256SUMS.json` byte-for-byte.

## Human submission and evaluation

Human submissions follow `annotation_submission.schema.json`. A and B must have distinct non-empty
annotator identifiers and exact evidence keys. Prepare a blinded disagreement packet with:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python -m \
  corpus.external_training.experiment_factory.publication_pilot.evaluate_annotations \
  prepare-adjudication --a submissions_A.jsonl --b submissions_B.jsonl \
  --output adjudication_packets.jsonl
```

The packet contains agreed and side-specific claims, but `automatic_adjudication_performed` remains
false. After a human adjudicator submits records whose state and every accepted claim are explicitly
`adjudicated`, machine precision/recall/F1 can be computed:

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python -m \
  corpus.external_training.experiment_factory.publication_pilot.evaluate_annotations \
  evaluate --machine results/machine_preannotations.jsonl.gz \
  --adjudicated adjudicated_submissions.jsonl --output benchmark.json
```

The evaluator refuses single-annotator input, missing identities, non-adjudicated claims, duplicate
sources, mismatched A/B source sets, or the same person occupying both annotation slots.
