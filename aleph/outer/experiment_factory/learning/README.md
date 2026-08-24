# Text-free multimodal learning projection

This directory provides a deterministic **proposed** projection over all 5,675 OA
source families. It does not train or modify Aleph's sealed learning packages.

The projection joins five independently masked views by stable source, experiment,
observation, visual-descriptor, asset, and dataset-manifest IDs:

- text metadata: domain and structured assay/cell tags only, never article text;
- numeric: structured values and unit records in an ignored local observation view;
- visual metadata: caption-derived tags and receipt counts, without captions;
- visual descriptors: IDs/digests and source-level numeric aggregates for genuinely
  `described` assets; caption/index availability is not counted as descriptor coverage;
- dataset manifests: public accession/component IDs and modality hints.

Identical source-payload SHA-256 values and shared public-accession components are
unioned before a group-level deterministic train/validation/test split. The frozen
300-source annotation candidate set is marked as a separate evaluation cohort but
is explicitly not treated as ground truth.

Ambiguous numeric observations and records with missing cell/protocol context have
`target_mask=false` in the ignored local numeric view. Missing modalities remain
explicit; the projection never imputes them. Full numeric rows stay beneath ignored
`data/external_training/experiment_factory/learning/`.

The builder accepts every optional input independently. If the visual or dataset
producer is absent, it emits a stable missing-view projection and records the absent
input in the summary. Re-running with a completed producer updates only bound view
digests and masks.
