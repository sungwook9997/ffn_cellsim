# Experiment annotation guidelines v1

## Status and scope

The 300 papers are `gold_candidate` sources. They are not a human-verified gold set. A record becomes
`adjudicated` only after two independent annotations and a third conflict-resolution pass. Even then,
the record remains `authority_status: proposed`; dataset adjudication does not validate a biological
claim or authorize Aleph parameter inference.

The target is broad cell and tissue mechanobiology. Annotate cortical tension, but apply the same
rules to adhesion/ECM, cell state, cytoskeletal rheology, division/morphogenesis, flow/shear,
measurement methodology, mechanotransduction, migration/invasion, nuclear mechanics, organelle
mechanics, and osmotic/poroelastic measurements.

## What constitutes one experiment

Start a new `ExperimentRecord` when the paper changes the biological system, independent replicate
cohort, perturbation design, measurement protocol in a way that changes interpretation, or the
control/treatment comparison. Multiple outcomes from the same cells and condition design belong in
one record as multiple observations. A time course is one experiment when the same sampling design
and comparison remain intact; destructive measurements on different cohorts must retain their
separate replicate semantics.

Do not merge records merely because panels share a figure. Do not split a record merely because one
experiment appears in the text, a table, and a figure.

## Annotation order

1. Pin `source_family_id`, PMCID, acquisition pass, XML payload digest, and exact source locators.
2. Identify organism, tissue, cell type, cell line, state, genotype, donor/batch attributes.
3. Define control, treatment, baseline, and timepoint conditions before extracting outcomes.
4. Record modality, instrument, calibration, and analysis pipeline as the protocol.
5. Extract observations with reported value, reported unit, uncertainty, sample size, and replicate
   unit kept separate.
6. Link comparisons to condition and observation IDs; do not infer direction from a discussion-only
   statement when the underlying result is absent.
7. Mark every unresolved field `missing` or `ambiguous`, including candidates and a reason for the
   latter.

## Biological systems and states

Keep cell line, cell type, tissue, organoid, donor and organism distinct. A cell line name is not a
cell ontology identity. Record serum starvation, confluence, differentiation, disease, passage,
cell-cycle state and activation state only when reported. If a mixed population is assayed, retain
the mixture rather than selecting the most prominent cell type.

Normalize ontology terms only when the mapping is exact or defensible. Otherwise retain the source
label with `candidate`, `unmapped`, or `ambiguous` status.

## Conditions and comparisons

Encode dose, duration, substrate, ECM composition, stiffness, geometry, media, temperature, osmotic
condition and flow condition independently. A vehicle control and an untreated control are distinct
arms. Paired measurements must not be represented as independent arms. If the same panel compares
several doses or timepoints, preserve all arms and use `dose_response` or `time_course` rather than
creating unsupported pairwise claims.

## Values, units and uncertainty

The reported value and unit are immutable source facts. A normalized value is an additional field
and requires a formula plus evidence. Never overwrite `nN/µm` with `pN/µm`, percent with fraction,
or apparent modulus with Young's modulus. Dimensionless, missing-unit and not-applicable states are
different.

Record SD, SE, confidence interval, credible interval, IQR and range as different uncertainty types.
Do not convert error bars without a stated definition. Preserve whether `n` means cells, fields,
wells, cultures, animals, donors, organoids, sections, or experiment days. Biological and technical
replicates are separate; pooled cells are not independent biological replicates.

Derived values require the formula, input observation IDs and unit conversion evidence. Digitized
values must be marked extracted, not reported, and retain the panel calibration status.

## Modality-specific rules

- **IF/confocal/live imaging:** retain channel identity, panel, scale/calibration state, segmentation
  unit, normalization reference and whether values are per-cell, per-area or population summaries.
- **Western blot:** retain target band, loading control, lane/condition correspondence and whether
  the outcome is raw intensity, normalized ratio or fold change. Never infer absolute abundance.
- **PCR/qPCR:** distinguish endpoint PCR from qPCR; retain target, reference gene, Ct/ΔCt/ΔΔCt
  transformation and efficiency assumptions. Fold change without its reference is incomplete.
- **PIV:** retain spatial and temporal calibration, interrogation window, frame interval, coordinate
  system and vector-field mask. A velocity field is not a scalar mean.
- **TFM:** retain substrate modulus/model, Poisson ratio, bead imaging, displacement method,
  regularization, dimensionality and force reconstruction convention. A traction map is a vector
  field; integrated force and strain energy are separate observations.
- **AFM/micropipette/tweezers:** retain probe geometry, loading rate, indentation/deformation range,
  constitutive fit, calibration and contact model. Apparent stiffness is not interchangeable with a
  constitutive modulus or cortical tension.
- **Tables/curves:** pin row/column or axis/series identity. Preserve every transformation needed to
  obtain a scalar. A screenshot is not a numeric table.

## Evidence locators and digests

Every protocol, observation and comparison must point to an exact paragraph, figure/panel, table
cell, supplement or asset. Store the whole XML digest and the text/asset digest. Character offsets
may be null for figures and tables, but their locator ID and digest may not be invented. Do not
commit publisher payloads or long excerpts into this metadata set.

## Double annotation and adjudication

Annotators A and B work independently. Agreement is evaluated at the experiment boundary, entity
identity, condition graph, value/unit, uncertainty, replicate unit and evidence locator—not merely at
the paper level. The adjudicator sees both records, records the chosen resolution, and must not
silently average numeric disagreements. `adjudicated` requires two distinct annotator IDs and one
adjudicator ID. Rejected candidates retain a reason and provenance.

Shared-accession components are assigned as a unit when constructing train/validation/test splits.
Accession co-membership is a leakage precaution, not proof that two papers used identical data.

## Forbidden shortcuts

- Do not label machine output as human verified.
- Do not infer sample size from plotted dots without explicit confirmation.
- Do not treat citations or discussion claims as measurements from the current source.
- Do not use journal impact factor as record quality.
- Do not collapse missing, zero and below-detection-limit.
- Do not use this dataset to fit physical parameters before a separate authority and validation gate.
