# A provenance-preserving multimodal experiment-record factory for cell-mechanics literature

Manuscript state: internal methods/data draft; not submission-ready until the frozen human benchmark is
completed. Citation placeholders must be replaced from the verified bibliography before submission.

## Abstract

Cell-mechanics evidence is distributed across prose, tables, microscopy, force maps and public data
repositories, which makes conventional document-level retrieval insufficient for model training and
simulation calibration. We developed a provenance-preserving factory that converts open-access
literature into proposed ExperimentRecords while retaining exact document, locator and asset digests.
The system processed 5,675 OA JATS articles and emitted 52,725 schema-valid experiment candidates and
238,489 numerical observations. It indexed 38,316 figures, linked 14,302 figures to canonical experiment
records, registered 1,801 public datasets from 1,982 source-accession bindings, and constructed a
text-free 5,675-source multimodal learning projection. Source payloads and accession-connected
components were grouped before splitting, producing 5,560 leakage groups and zero observed cross-split
violations. In a weak twelve-route diagnostic benchmark, adding numeric summaries increased test
macro-F1 from 0.2846 to 0.3086 and adding sparse visual descriptors increased it to 0.3166; dataset
manifests produced mixed effects. These scores measure within-corpus routing rather than biology. We
freeze a 50-paper, 12-domain, 12-modality blinded dual-annotation pilot and an adjudication-aware
field-level evaluator. Human benchmark results are pending and no automatically extracted record is
called gold or used as physical authority.

## Introduction

Mechanobiology experiments combine biological identity, perturbation, geometry, instrumentation,
calibration, quantitative observations and visual evidence. Article-level labels erase this structure.
They also create leakage when multiple articles reuse the same public dataset or when captions,
supplements and repository records are independently sampled across train and test partitions.

We address the infrastructure problem rather than claiming an autonomous biological reader. The system
represents one interpretable experiment or comparison as a graph of source, biological system,
condition, protocol, observation, comparison and provenance. Missing and ambiguous states remain
explicit. Machine extraction is separated from human review, and all committed learning views omit
full text and pixels.

The present contribution is fourfold: (i) a versioned modality-neutral ExperimentRecord contract;
(ii) deterministic text, table, numeric, visual and public-dataset candidate generation with exact
provenance; (iii) source/accession-connected leakage control and multimodal missingness masks; and
(iv) a frozen human-validation protocol that mechanically prevents single-annotator or machine output
from being reported as adjudicated gold.

## Methods

### Corpus and authority boundary

The input comprised 5,675 digest-verified OA JATS payloads acquired through official Europe PMC paths.
Publisher-only and non-OA material was excluded from automated acquisition. Every derived record has
`authority_status: proposed`. OA payloads and public-dataset pilot files remain ignored local objects;
the repository contains only structured metadata, hashes, models and reports.

### ExperimentRecord representation

The schema models biological systems, condition arms, protocols, observations, comparisons and exact
evidence locators. Unit records distinguish reported and normalized forms and retain conversion evidence.
Replicate records separate biological and technical sample size. Visual observations distinguish asset,
panel, channel, axis and calibration state. Semantic validation checks referential integrity, evidence
digests, units, ambiguity and review-state requirements.

### Candidate extraction and negative controls

The extractor emits explicit numeric tokens from eligible evidence classes and records the exact locator
and character span. Case-sensitive rules distinguish concentration `nM`, length `nm`, force notation and
sample-size `n`. Regression controls reject mutation and biological-name confounders including `W748S`
and `18S`. Drug-context micrometre-like tokens that cannot be resolved between length and concentration
remain ambiguous.

### Visual evidence

The visual index records figure and graphic identities without copying captions into the learning view.
Only official PMC-page-listed NCBI CDN assets were attempted. Technical image descriptors and panel
candidates are computed for acquired objects. The pipeline does not extract scientific values from
pixels. Unresolved and deferred graphic identities remain typed records rather than disappearing from
the denominator.

### Public dataset registry

Exact accessions were detected for GEO, SRA, BioProject, ArrayExpress, PRIDE, Zenodo, Figshare, Dryad,
OSF, BioImage Archive and related repositories. Official APIs were queried for a bounded subset.
Accession co-membership defines a conservative leakage component but is not treated as proof of identical
samples or processing. A small acquisition pilot required an official URL, explicit open licence,
declared size and local digest.

### Learning projection and grouped split

The text-free projection joins structured metadata, numeric summaries, visual availability/descriptors
and dataset manifests by stable identifiers. Raw text, captions, pixels and publisher payloads are not
copied. Ambiguous or context-deficient numeric observations have a false target mask. Identical payload
hashes and accession-connected sources are unioned before deterministic group assignment.

### Diagnostic weak-route baseline

A class-weighted NumPy linear softmax classifier predicts only the corpus's twelve weak discovery routes.
Four nested views add numeric, visual and dataset features to metadata. The same grouped split and held-out
test support are used for every view and seeds 11, 29 and 47. The frozen annotation cohort is excluded
from fitting, early stopping and test metrics. Because discovery routes and assay/cell tags partly arise
from the retrieval design, results are interpreted as projection diagnostics rather than external
biological generalization.

### Frozen human-validation pilot

Fifty papers were selected from the frozen 300-paper cohort. The sample covers all twelve domains and
twelve modality tags, both acquisition passes and four publication strata. Public-accession components
are indivisible. Each source has independent blinded A and B forms. Machine pre-annotations are hidden
during those passes. A comparison tool prepares agreed and side-specific claims but cannot adjudicate.
Only a distinct human adjudicator may emit an adjudicated record. Field-level exact precision, recall and
F1 are unavailable until all human submissions are frozen.

## Results

### Factory yield and integrity

The factory processed 5,675/5,675 articles and emitted 52,725 ExperimentRecords and 238,489 observations.
All records passed semantic schema validation; experiment and observation identifier collisions were
zero. Correcting case and context false positives removed 10,991 records and 30,162 observations relative
to the first rule set. The final learning projection contains 5,675 sources in 5,560 leakage groups with
train/validation/test counts 4,004/840/831 and zero observed component violations.

### Visual and dataset coverage

The visual layer contains 38,316 figure records and 14,302 exact experiment joins. The frozen visual
attempt ledger contains 3,315 graphic identities; 627 references were described and resolved to 480
unique objects across 116 sources. All referenced object digests validated. The dataset layer contains
1,982 source-accession bindings and 1,801 unique manifests across 18 providers. A bounded 100-record
official metadata pilot resolved 81 records. Nine openly licensed files totalling 3,487,812 bytes were
downloaded locally and none were committed.

### Weak-route ablation

Across three seeds on the same 787-source test set, metadata achieved macro-F1 0.2846 ± 0.0026. Numeric
summaries increased macro-F1 to 0.3086 ± 0.0023, and visual features increased it to 0.3166 ± 0.0036.
Adding dataset-manifest features reduced macro-F1 to 0.3148 ± 0.0007 while slightly increasing micro-F1
and worsening calibration. The mixed dataset result was retained rather than selected away.

### Validation protocol readiness

The 50-paper pilot has 100 blank assignment slots and 996 separately stored machine ExperimentRecord
candidates. Domain counts are four each except cortical/membrane-pressure and measurement/inference,
which contain five. Year-stratum counts are 7/11/12/20 and acquisition-pass counts are 27/23. No public
accession component is split. Human A, B and adjudicated submissions are currently 0/0/0; therefore no
human agreement, extraction accuracy, failure prevalence or gold-label count is reported.

## Discussion

The measured results show that document-scale OA evidence can be converted into a deterministic,
traceable and leakage-aware experiment representation. Numeric and sparse visual features add route
information beyond metadata, but the modest scores and mixed dataset ablation underscore that this is
not yet a biological state model. The primary scientific limitation is absent human adjudication.

The system is designed to support the broader next stage: author-labelled public raw datasets,
modality-specific encoders, independently held-out laboratories and prediction of experimentally defined
cell states and observables. Those future models must pass the frozen annotation gate and out-of-dataset
evaluation before they can propose Aleph parameter constraints. No current output updates physics priors.

## Reproducibility and independent verification

An independent verifier checked 300 candidate and 52,725 ExperimentRecord schemas, 22 producer-manifest
files, all 480 referenced visual objects, accession-component split integrity and model replay. It first
identified multithreaded-BLAS weight differences up to `3.33e-16`; the controlled CLI was corrected to
initialize five thread controls before NumPy. Subsequent clean launches reproduced all six model artifacts
byte-for-byte. The final factory plus publication-pilot suite passed 58 focused tests.

## Data and code availability

Structured schemas, compact derived metadata, deterministic generators, model artifacts and audit reports
are maintained on the isolated external-training branch. OA article payloads, images and public raw-data
pilot files remain local content-addressed objects and are not redistributed in the repository. A release
must include per-provider licence receipts and exclude material not authorized for redistribution.

## Claims deliberately not made

- Machine candidates are not human gold.
- The weak-route baseline is not a cell-state or physical-parameter predictor.
- Technical visual descriptors are not scientific values extracted from images.
- Accession co-membership does not prove identical samples.
- Passing structural audits does not validate biological truth.

## References

`[VERIFIED REFERENCES TO BE INSERTED FROM THE PROJECT BIBLIOGRAPHY; DO NOT INVENT OR SUBMIT PLACEHOLDERS]`
