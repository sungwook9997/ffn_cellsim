# ExperimentRecord dataset-factory orchestration

This directory integrates the external-corpus experiment dataset factory. It owns no scientific
truth and promotes no record beyond `proposed`.

The measured end-to-end outcome, limitations, ablations, reproducibility correction and producer
commit chain are recorded in `FINAL_REPORT.md`.

The factory's completion boundary is broader than document routing:

1. a versioned, modality-neutral experiment schema and controlled vocabularies;
2. a stratified annotation-candidate set with leakage groups kept intact;
3. conservative full-corpus extraction of biological systems, conditions, protocols, comparisons,
   numerical observations, units, uncertainty and replicates;
4. local-only OA visual assets and deterministic image/panel descriptors with exact provenance;
5. a text-free, versioned learning projection and an external CPU baseline evaluated with grouped
   splits.

Human review remains a separate state transition. `gold_candidate` never means gold, and an
automatically extracted number never becomes a physical parameter authority.

## Promotion boundary

The factory has three deliberately separate states:

1. `machine_candidate`: deterministic extraction with exact source and locator provenance;
2. `gold_candidate`: a paper selected for blinded annotation, not a verified label;
3. `adjudicated`: two independent annotations plus adjudication, which is not produced by this
   automated pipeline.

Only structured, text-free projections may enter model training. Ambiguous observations retain an
explicit mask and cannot be promoted to physical targets. Raw JATS, figure pixels and publisher
payloads remain local OA evidence and are never copied into the committed learning artifacts.

## Required leakage boundary

Train, validation and test partitions are grouped by source payload and by connected components of
public dataset accessions. A paper pair connected through GEO, SRA, BioProject, ArrayExpress, PRIDE,
Zenodo, Figshare, Dryad, OSF or another registered public accession must remain in one partition.
The frozen 300-paper annotation cohort is reported separately and is not treated as a gold test set.

## Completion gates

The integrated report must publish measured counts for every modality, missing and ambiguous rates,
schema failures, identifier collisions, unresolved official metadata, unavailable OA assets and
split-component leakage. Model results must compare identical grouped splits for metadata-only,
numeric and visual/dataset-manifest views. These are routing and representation baselines only; they
do not establish biological validity or modify Aleph's physics modules.
