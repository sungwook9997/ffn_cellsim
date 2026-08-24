# Outer Library: evidence-bound multimodal training data

This package turns public author-labelled experimental source data into a
common, cell/sample-grouped observation schema. It is intentionally separate
from Aleph physics and has no permission to write physics parameters or decision
authority.

The near-term objective was deliberately changed because Aleph parameter sweeps
are not yet executable at the required scale. The current deliverable is an
independently useful external biophysical inference network: it must predict
observable cell state, uncertainty, mechanism compatibility, and the next useful
measurement without requiring an Aleph run. This is more than a scheduling
change and the network is not judged only as an Aleph helper.

The original sweep objective is nevertheless preserved as a future integration
path. Validated external outputs may later compress or condition a sweep only
through `latent_promotion.py`: independent dataset and lab holdouts, calibrated
uncertainty, resolved identifiability, an exact unit/geometry-matched Aleph
observation operator, and forward-sensitivity validation are all mandatory. The
promotion audit never chooses a parameter or emits a numeric range by itself.

## What is implemented

- A protocol-conditioned external mechanics network spanning 19 real cohort
  records from nine primary papers and 196 author-released individual magnetic
  wire measurements. It separates method, mechanical scope, adhesion state,
  probe scale, and unit; reconstructs effective viscosity on a sealed
  experimental-date holdout, rejects cell-line transfer after that head fails,
  and preregisters the `tau(a) ~ a^2` versus `a^0` discriminator without
  equating wire length to AFM contact radius.
- A protocol-conditioned cortical-tension reference head for MCF-7
  suspended/rounded interphase cells. It resolves the structured median/IQR
  summary against `KB-6.1.7`, `SE181`, and `source_audit=OK`, records the still-
  draft validation-gate status, and provides a value-redacting adapter for
  later full-native engine comparisons. MDA-MB-231 and adherent MCF-7 remain
  explicit coverage gaps rather than receiving transferred values.
- Exact MAT-index, XLSX-cell, CSV-cell, or archive-member provenance.
- Author labels only; missing units, replicate IDs, and cell types remain
  explicitly unknown.
- TFM-, cantilever-force-, IF pixel/table-, nuclear-shape-, RT-qPCR-, PIV-,
  FRET-, and calcium time-series adapters.
- Cell-grouped train/validation/test splitting.
- Calibration, semantic OOD, independent-lab OOD, and a full-lab task transfer
  evaluation.
- Lab-invariant within-study effect contrasts for mechanism-level observable
  evidence when absolute single-cell classification does not transfer.
- Multimodal mechanical-state coherence and same-lab cell-type holdout reports
  that cannot be mislabeled as independent-lab validation.
- A disjoint training/validation/calibration conformal audit that measures
  external-lab coverage and dataset shift without tuning on external labels.
- A three-lab paired migration-state path that retains person, independent
  experiment, and technical-replicate hierarchy: lab 1 learns the ROCK-
  inhibition direction, zero is the fixed paired no-effect threshold, lab 2
  validates/calibrates uncertainty, and lab 3 is an external retrospective
  test. A new lab is required for pristine confirmation of this revision.
- A CC0 H1299-PL multiscale migration corpus with 22,647 tracked-cell time
  observations, 60 image-derived features, explicit experimental-date/cell
  groups, and author Continuous/Discontinuous labels. A compact NPZ retains
  every row while the canonical manifest exposes cell-level robust summaries.
- A leakage-safe organizational-state model that excludes the five behavioral
  features used to define migration mode, selects on fibronectin dates,
  calibrates on disjoint dates, and retrospectively tests talin perturbation.
  This strengthens same-lab representation learning but is never counted as a
  new independent MULTIMOT lab.
- An independent eLife NIH3T3 3D cell-derived-matrix source with per-cell
  persistent speed and oscillation period under control, Y-27632, ML-7,
  CK666, and C8-BPA. It exposes a context-dependent ROCK effect rather than
  being forced into a universal speed constraint; experiment IDs are not
  reported, so cells are not promoted to biological replicates.
- An independent eLife HeLa mechano-osmotic source with 315 cell-spreading
  samples nested in seven author experiments and 314 AFM tether-force cell
  points nested in the published experiment groups. It preserves pN force,
  square-micrometre/minute spreading rate, and cubic-micrometre/minute volume
  flux, and evaluates treatment effects only over experiment means.
- A time-resolved operator guard: the strong early Y-27632 tether-force
  contrast is evidence for an apparent-membrane-tension proxy, not permission
  to equate AFM force with cortical tension or select a Helfrich parameter.
- A CC0 human dermal fibroblast collective-motion source with 7,645 raw
  single-cell velocity vectors, 3,154 canonical longitudinal/field
  observations, and 832 technical samples. All vectors and time points carry
  one mandatory lab/dataset split group because independent culture IDs are
  not reported. It is admitted for representation learning only, never as an
  external holdout or FAK-effect calibration source.
- An independent CC0 Dryad C2C12 supracellular-contractility source with 12
  author-labelled Fig. 2 LCN substrates and four Fig. 4 LCN substrates. It
  contributes 1,435 live-actin/PIV observations and 720 fibronectin IF angular
  frequencies. Spreadsheet fill is treated as data: all 112 red author-
  excluded cells are omitted, while nine green maximum-density alignment cells
  are retained and audited. The two control and two blebbistatin substrates
  reproduce the late order, nematic-correlation, velocity-correlation, and ECM-
  alignment direction, but n=2 per condition and missing aligned-lab/operator
  validation keep this as a directional pre-sweep constraint only.
- A 2026 CC0 Science/Dataverse NIH3T3 cellular-nematic source with 63 PIV
  frames (239,150 documented measured/interpolated vectors) and the complete
  12-field Figure S9 paired TFM/MSM set. The published before/after
  distributions are reproduced from `t0 -> t1`: mean traction falls from
  136.89 to 30.18 Pa and maximum principal monolayer tension from 25.40 to
  5.18 mN/m, with a decrease in 12/12 fields. Fields are nested technical
  views of one dated culture/batch, not twelve biological replicates. The
  ambiguous `t2` table and incomplete f13--f30 folders are withheld.
- A second independent CC0 eLife/Dryad opto-MDCK contour source with 68
  globally photoactivated doublets. RhoA activation yields positive author
  RSI in 47/68 x-axis and 46/68 y-axis objects. The released
  `stats_fullstim.csv` summarizes a different 35/15-object group than the
  22/30/16 groups in `fullstim_data.csv`; it is therefore excluded, as are all
  nine untrusted pickle members. Missing author sample/date identities force
  all 68 objects into one indivisible split group.
- A 2026 CC BY 4.0 Zenodo/KU Leuven MEF source that joins cell-level TFM,
  FLIM-FRET vinculin tension sensing, and focal-adhesion-level paired
  traction/FRET measurements. It contributes 847 canonical observations from
  358 cells or cell-grouped samples. Mean traction rises from a median 163.13
  Pa at 4.5 kPa to 220.22 Pa at 13 kPa, reproducing the stiffness direction of
  the independent NIH3T3 TFM source in the same unit. The released Figure 4
  tables also reproduce direct, inverse, and uncorrelated traction--FRET cell
  examples, preventing a false one-to-one load-transfer mapping. The paper
  reports at least three experiments, but released cells are not mapped to
  replicate IDs, so this remains a directional/operator constraint rather
  than an inferential holdout.
- A hard sweep gate that refuses proposals when validation is insufficient.
  It may emit a blocked calibration contract, but never an executable numeric
  axis before identifiability and external validation requirements pass.
- A pre-sweep latent-promotion contract connecting today's independent external
  inference layer to the eventual Aleph sweep. It explicitly distinguishes
  training-only, externally validated evidence, operator-ready constraints,
  and bounded-proposal eligibility.
- An applicability-gated theory registry with executable membrane-tether,
  single-mode viscoelastic-relaxation, and 2D persistent-walk experts. The
  tether expert returns the exactly identifiable product or apparent tension
  conditioned on bending rigidity while preserving the unresolved bilayer
  tension versus membrane-cortex adhesion family.
- A twenty-eight-head, task-scoped multimodal model contract spanning calcium,
  migration, IF/RNA fibrosis, mechano-osmotic observables, PIV, and collective
  velocity. It validates train/calibration/test lab separation and grants
  evidence authority per head; it does not claim that a shared neural encoder
  has already been trained.
- A preregistered PBMC pre-sweep context head. GSE96583 contributes 24,413
  author-labelled singlet cells split by whole donor; 508 fit-donor-selected
  genes align exactly to 22,594 retained author-annotated GSE132044 cells from
  a distinct laboratory. The sealed donor macro-F1 is 0.924 and the pristine
  external macro-F1 is 0.720, with ECE 0.046 and semantic-OOD AUROC 0.857.
  External 90% conformal coverage is only 0.502, so the head is explicitly
  blocked from deleting candidate axes or emitting parameter ranges. A frozen
  GSE164378 NYGC recalibration improved same-lab coverage to 0.900 without
  changing classifier weights. The subsequent single-use GSE150728 Stanford
  healthy-donor confirmation retained macro-F1 0.746, ECE 0.112, and
  semantic-OOD AUROC 0.725, but external coverage was only 0.694; B/DC/monocyte
  class-conditional coverage was 0.611/0.566/0.265. The head therefore remains
  uncertainty-blocked after a genuine second-lab calibration and third-lab
  confirmation. A separate
  GSE72502 three-donor paired test confirms the frozen IFN signature direction
  in all donors; its absolute scale is not transferred.
- A paper-linked M-TRACK A549 VIM-RFP temporal pilot with ten paired 1952x1952
  fluorescence/mask frames from one TGF-beta 4 ng/mL `xy01` trajectory. All
  frames are forced into one split group. Mask area rises from an early-three
  median of 34,768 to a late-three median of 63,167 pixels, while masked
  VIM-RFP mean falls by 15.123 gray units (Spearman time correlation -0.915).
  This conflicts with the ExIF fixed-IF TGF-beta/control Vimentin direction,
  so the path is an intensity-only sweep-conditioning blocker and temporal
  representation diagnostic—not an independent holdout or a numeric axis.
- A controlled CC BY 4.0 LINCS MCF10A IF phenotype source with 159 author
  well-level rows, 1,263 canonical observations, and seven released
  `collection + replicate` biological groups. Technical wells are median-folded
  before inference. A C1-trained TGFB+EGF-versus-EGF phenotype classifier scores
  8/8 at 24 h and 5/6 at 48 h on the later C2 collection, both with AUROC 1.0.
  KRT5 mean decreases, normalized neighbor distance increases, and mean cluster
  size decreases in every paired replicate in both collections at 48 h. This is
  a same-OHSU collection holdout, not an independent-lab holdout; semantic OOD
  remains unvalidated and the 10-versus-20 ng/mL companion-EGF discrepancy is
  retained as unresolved, so no numeric Aleph range is emitted.
- A pinned CC BY 4.0 Stanford HCC827 EMT CyTOF source with 90,066 raw
  single-cell events, 29,066 transformed cells in eight author CCAST
  EMT/MET states, six state markers, and seven treatment/withdrawal timepoints.
  The builder streams OOXML rather than materializing the 79 MB workbook as
  office objects. From 0 d to 10 d TGF-beta, raw median E-cadherin falls from
  41.996 to 1.549 and Vimentin rises from 2.136 to 2,877.259; both partially
  reverse after 10 d withdrawal. These directions reproduce the independent
  UNSW A549 fixed-IF observations across dataset, lab, cell line, and modality.
  The release contains no biological-replicate ID, however, and CCAST labels
  were derived from the same six markers. All 90,066 cells therefore remain
  one split group, marker-based state classification is explicitly circular,
  and the result is a pre-sweep direction filter rather than an inferential
  holdout, uncertainty calibration, OOD authority, or numeric parameter range.
- Five independent GEO EMT RNA studies add 72 genuinely biological samples
  from the University of Michigan (GSE17708), University of Virginia
  (GSE42373), Stony Brook University (GSE125369), Shanghai Advanced Research
  Institute (GSE69667), and Genentech (GSE49644). The ten declared
  study-by-gene endpoint contrasts all reproduce CDH1-down/VIM-up and every
  treated value is separated from every within-study control value. RNA-seq
  and microarray absolute scales are never pooled, and RNA is not equated to
  IF/CyTOF protein abundance. Five retrospective independent signatures give a
  nominal exact one-sided p=0.03125. Because the added studies were selected
  from literature reporting successful EMT, this is not an outcome-blind
  confirmatory holdout and receives no inferential authority. The head cannot
  emit a parameter axis or numeric range; a preregistered external study,
  numeric uncertainty, observation operators, and forward sensitivity remain
  absent. GSE42373's discordant 2D CDH1 response is correctly refused by the
  sealed direction diagnostic, but one OOD context does not validate general
  semantic-OOD refusal.
- A genuinely outcome-blind GSE325309 confirmation was protocol-frozen in
  commit `7984e5b` before its workbook was downloaded. All six biological
  samples, exact CDH1/VIM symbols, no exclusions, log2(FPKM+1), and the exact
  20-assignment statistic were fixed in advance. The observed directions are
  CDH1-down/VIM-up and the frozen composite test reaches its finite-sample
  minimum `p=0.05`. This promotes only the direction filter to a pristine
  external holdout; numeric scale, semantic OOD, and Aleph sweep-axis authority
  remain blocked.
- A three-lab, three-cell-type actomyosin pre-sweep constraint combines C2C12
  PIV/IF, NIH3T3 TFM/MSM, and opto-MDCK contour directions. This raises
  mechanism-family ranking confidence but is explicitly not an aligned
  same-observable holdout and emits neither an Aleph parameter nor a numeric
  range. Its job is to reduce a future forward sweep after matched observation
  operators, biological replication, uncertainty/OOD, identifiability, and
  monotonic sensitivity all pass.
- A verified CC BY 4.0 Zenodo 7432971 TFM/tether source. The 4.94 GB archive
  central directory and four selected XLSX member CRCs are verified by byte
  ranges, avoiding a full archive download. It contributes 60 NIH3T3 TFM
  cells, 132 Xenopus growth-cone TFM cells, and optical-tweezer tether
  measurements, with frames collapsed before analysis. Cells are not promoted
  to biological replicates. The head contributes to the independent external
  mechanical-state posterior; only an optional later adapter may use it to
  narrow an Aleph sweep after measurement-operator and forward-sensitivity
  validation.
- A focal-adhesion load-transfer sweep contract that combines the prior
  NIH3T3 stiffness--TFM data with the independent MEF TFM--FRET source. It
  requests cell-mask traction, FA traction, relative vinculin tension, and FA
  geometry operators while keeping `numeric_range: null` until replicate-
  complete external validation, calibrated uncertainty/OOD, and Aleph forward
  sensitivity all pass.
- Two verified CC BY 3.0 BBBC paired-channel IF plates: 96 U2OS
  FKHR-EGFP/DRAQ fields (BBBC013) and 96 MCF7/A549 NF-kB FITC/DAPI fields
  (BBBC014). A dependency-light BMP decoder pairs channels at the well before
  extracting 12 field observables. The sealed cross-dataset check retains all
  author dose/cell-type labels and detects dataset OOD, but fails the declared
  MCF7 rank threshold and conformal-width criterion; it is training-only and
  is not mislabeled as an independent-lab holdout because both plates share
  Ilya Ravkin provider lineage.
- The official HPA v25.1 CC BY 4.0 subcellular image-embedding export: 81,007
  IF images, 13,366 genes, 15,700 antibodies, 39 cell lines, 35 author
  localization labels, and 1,024 provider image features. A compact fixed
  128-D view is split by connected gene--antibody components, with zero gene
  or antibody overlap. Its grouped test remains training-only because all rows
  share one provider and annotation pipeline.
- A machine-generated adversarial coverage audit that separates observation
  rows from independent samples, quantifies unknown units and replicate IDs,
  and maintains an explicit missing-state/modality priority queue.
- A deterministic 512-image HPA raw-IF pilot spanning all 35 localization
  labels and four original channels while preserving the gene--antibody split.
- A frozen cross-provider raw-IF transfer diagnostic on 2,593 evaluable IDR0072
  confocal fields. The HPA checkpoint fails every preregistered performance and
  uncertainty endpoint (macro AP 0.317, macro F1 0.306, ECE 0.237, positive
  conformal coverage 0.806, mean set size 5.72/9). IDR0072 is consumed and
  cannot be reused to tune or validate a replacement model; the resulting head
  remains `blocked_external_domain_shift` with no evidence authority.
- A preregistered Light My Cells visual-observation branch using the official
  S-BIAD1047 metadata: 56,984 OME-TIFF files, 2,574 acquisition sets, 30
  independent studies, and 109.29 GB. Studies are frozen as 16 fit / 4
  selection / 4 calibration / 6 single-use confirmation. Its 102 nuclear and
  91 mitochondrial pairs are sufficient for the two mandatory cross-study
  targets. This branch predicts spatial
  fluorescence from BF/PC/DIC and is distinct from HPA protein-localization
  classification.
  A frozen OME-header resolver recovers 4,375/4,376 non-confirmation pairs:
  3,253 by physical `PositionZ`, 57 singleton inputs, and 1,065 explicitly
  marked lower-median fallbacks. One malformed target is refused. The selected
  7,933 payloads total 12.68 GB instead of downloading the full 109.29 GB.
  Provider integrity checking then refuses one metadata-only BigTIFF whose
  official FileList incorrectly claims 6.1 MB. The final training asset is
  4,374 pairs backed by 7,932 locally SHA-256-pinned files (12.673 GB) and a
  finite `[4374,64,64]` float16 tensor; no confirmation file is present.
  The first 26,641-parameter conditional U-Net is a frozen development
  diagnostic, not a completed model: cross-study selection gives Pearson 0.198,
  global SSIM 0.070, and 8.59% MAE improvement over the fit-only mean map. All
  miss the frozen readiness targets, so the six confirmation studies remain
  unopened while a stronger spatial architecture is selected.
  The frozen v2 two-level residual U-Net corrects that failure on the same
  study-disjoint development splits: Pearson 0.558, global SSIM 0.475, and
  30.09% macro-study MAE improvement. Nuclear and mitochondrial calibration
  coverage is 94.38% and 94.10%. All selection gates pass, while production
  authority remains false. The eligible confirmation was then consumed exactly
  once: 226 pairs across six studies and five experimenter groups, including
  two groups absent from training. SSIM 0.275, relative MAE improvement 11.92%,
  and interval coverage 92.31% pass, but Pearson 0.329 misses 0.35 and
  study-level embedding-OOD AUROC 0.385 misses 0.70. The resulting 30th head is
  `blocked_external_spatial_transfer`; the confirmation set cannot tune a
  replacement, and no evidence, production, uncertainty, or Aleph authority is
  granted.
- A verified BBBC054 LPS microglia time course with 58,186 exact-cell author
  morphophenotypes. Whole fields, never cells from one field, define the split;
  the weak handcrafted classifier is explicitly rejected and only the
  descriptive immune-response trajectory is retained for training context.
- The official BBBC048 Jurkat imaging-flow cell-cycle set with 32,266 cells and
  96,798 labelled channel images (brightfield, MPM2, and PI/DNA). All three
  channels remain in one cell split group. The seven author phases are highly
  imbalanced, and the single provider/batch is training-only rather than an
  external holdout.
- A deterministic four-way cell split and compact three-channel raw-image MLP.
  Selection, calibration, and test cells are disjoint. On the same-provider
  test it reaches 0.614 macro F1 and 0.651 balanced accuracy; 90% split-
  conformal marginal coverage is 0.895 with mean set size 1.28/7. Rare-phase
  counts are far too small for conditional coverage or production authority.
- Verified sci-Plex3 metadata for 799,317 single cells across A549, MCF7, and
  K562, two replicates, 24/72-hour contexts, 188 inhibitors plus vehicle, and
  five dose levels. Of these, 762,795 cells have complete context metadata;
  36,522 unassigned cells are excluded from supervised state training. The
  4,896 replicate-treatment-dose-time-cell-type groups are mandatory split
  units. The official 3,096,224,606-byte expression matrix passes gzip-stream,
  1,007,419,688 sparse-line, and SHA-256 verification and is training-ready.
- A source-derived 2,228-gene state feature basis from nine official human
  Reactome pathways. Every pathway API response is hashed; no hand-added marker
  genes are present. A streaming sparse-matrix builder aggregates only complete
  replicate-treatment-dose-time-cell-type contexts, keeping the full 110,983 x
  799,317 expression matrix out of dense memory.
- A treatment-disjoint sci-Plex baseline using 125/21/21/21 compound identities
  for fit/selection/calibration/sealed test, plus vehicle OOD. Fit-only nuisance
  residualization improves over raw expression but still fails general pathway
  transfer (sealed macro F1 0.198); nominal 90% sets average 11.52/16 classes
  and are uninformative. The tensor is retained while the model stays rejected.
- A treatment-disjoint sci-Plex baseline using 125/21/21/21 compound identities
  for fit/selection/calibration/sealed test, plus vehicle OOD. It fails general
  pathway transfer (sealed macro F1 0.139); its nominal 90% sets average
  13.33/16 classes and are therefore uninformative. The verified tensor is
  retained, while the model remains a rejected same-provider baseline.
- The official BBBC053 mitochondrial-stress IF archive: 29 DMSO and 29 FCCP
  murine CAD-cell TOM20 super-resolution fields, all CRC-verified and decoded
  from 2048x2048 16-bit TIFF. Compact raw-image and 13-feature morphology
  tensors are available for training, but missing biological-replicate IDs and
  the single Vitriol-lab source prohibit external or intervention authority.
- The CC0 S-BIAD2515 apoptosis resource: the official BioImage Archive indexes
  6,240 raw TIFFs (45.1 GB), while the pinned author-analysis commit supplies
  30 well profiles with 2,336 live morphology features and terminal Annexin-V.
  Fit, calibration, and sealed test each contain one well at all ten doses.
  Terminal features never enter the input. The same-provider test reaches R2
  0.878 and 90% well-level interval coverage, but has no dataset/lab authority.
- A verified independent programmed-cell-death design from Figshare 28202864:
  21,271 metadata sites, 15,430 QC sites, eight MCF7 plates, six author MoA
  classes, and three feature-extractor-specific well splits. Its Cell Painting
  mechanism labels are kept distinct from the S-BIAD2515 Annexin-V endpoint;
  the large feature archive must pass MD5 before training. Because 45--47
  compounds cross roles in the authors' well splits, the local evaluation uses
  a stricter zero-overlap 40/4/5/6 compound train/selection/calibration/test
  partition and refuses rare-class conditional calibration claims. The full
  2.452 GB feature archive is MD5-verified; only aggregate profiles are
  extracted. The resulting 954-well, 672-feature tensor reaches 0.500 balanced
  accuracy and 0.389 macro-F1 on six sealed compounds, with 0.667 marginal
  conformal coverage. It remains rejected same-provider representation
  training, not external evidence.
- An official NIH/NIA senescence design from GSE250041: 8,664 untreated
  proliferating and 4,949 10-Gy IR/day-10 senescent WI-38 barcodes share an
  identical 36,601-RNA plus eight-ADT feature axis. The state label is direct
  author metadata and the provider is independent of sci-Plex, but one pooled
  capture per condition means 13,613 cells are not 13,613 biological
  replicates. Both matrices are now gzip/SHA/shape verified and yield a compact
  13,613 x 43 RNA/ADT tensor. Capture is fully confounded with state, so it is
  representation training only and supplies no valid classification holdout.
- The CC0 SenSCOUT archive with 5,000 sampled cells and the authors' exact 87
  curated morphology features. Target biomarker intensities, author KMEANS,
  UMAP, labels, and identifiers are excluded from input. Bio1 alone supplies
  image-group-disjoint fit/selection/calibration; Bio2 is a sealed
  same-provider repeat and Bio3 is a frozen shift diagnostic. Direct BGAL
  balanced accuracy is 0.659 on Bio2, but image-group 90% interval coverage is
  only 0.444, so it is explicitly rejected for evidence or sweep authority.
  P16 is a promising candidate observable (Bio2/Bio3 Pearson r 0.770/0.808)
  pending independent dataset/lab validation and an Aleph observation operator.
- A frozen 35-marker transcriptomic senescence head fitted only on GSE63577
  cell lines, selected on IMR90, and calibrated on WI-38. GSE297406 reproduces
  all four paired passage directions but has sign-test `p=0.0625`, conformal
  coverage 0.125, and 100% OOD rejection. A genuinely post-freeze second lab,
  GSE301164, reproduces all three P33-vs-P3 and all three DOX-vs-DMSO directions
  (`6/6`, combined `p=0.015625`). Because GSE301164 remains 100% OOD, only the
  paired direction is admitted as external evidence; absolute state
  probabilities, uncertainty, and Aleph parameter selection remain refused.
- A dependency-light trainable 1D calcium CNN with exact gradients, deterministic
  L-BFGS optimization, validation-only architecture selection, and a hashed NPZ
  checkpoint. Its strong HeLa validation result fails EA.hy926 macro-F1, so the
  architecture is retained as a rejected diagnostic rather than promoted.
- A forward-sensitivity validator that refuses unit-mismatched total force,
  non-converged runs, under-replicated grids, non-monotonic responses, and
  interpolation outside the simulated response hull.

Generated observation manifests and compact row tensors are locally
reproducible but ignored by Git because they are generated data artifacts.
Their content hashes and counts are committed in
`results/summary.json`.

## Local reproduction

The XLSX adapters require the bundled spreadsheet runtime and the ignored
licensed/open source objects described by `acquisition_receipts.jsonl`. Python
steps use the Aleph environment.

```bash
PY=/Users/sw1/miniconda3/envs/aleph/bin/python
$PY build_tfm_manifest.py --input-dir ../../../../data/external_training/experiment_factory/datasets/pilot --output results/tfm_observations.jsonl
$PY build_nuclear_manifest.py --workbook-export ../../../../data/external_training/experiment_factory/outer_library/nuclear_workbook.json --output results/nuclear_observations.jsonl
$PY build_mechanotransduction_manifest.py --workbook-export ../../../../data/external_training/experiment_factory/outer_library/mechanotransduction_workbooks.json --output results/mechanotransduction_observations.jsonl
$PY build_mdck_calcium_manifest.py --source-zip ../../../../data/external_training/experiment_factory/outer_library/downloads/mdck_model_data_fit.zip --output results/mdck_calcium_observations.jsonl
$PY build_eahy_calcium_manifest.py --source-csv ../../../../data/external_training/experiment_factory/outer_library/downloads/piezo_nanoswitch/Figure_4_A_ii.csv --output results/eahy_calcium_observations.jsonl
$PY build_tendon_multimodal_manifest.py --workbook-export ../../../../data/external_training/experiment_factory/outer_library/tendon_source_data_1.json --output results/tendon_multimodal_observations.jsonl
$PY evaluate_tendon_multimodal.py --manifest results/tendon_multimodal_observations.jsonl --output results/tendon_multimodal_report.json
$PY build_gse226374_manifest.py --normalized-counts ../../../../data/external_training/experiment_factory/outer_library/downloads/gse226374/GSE226374_normalized_counts.txt.gz --output results/gse226374_observations.jsonl
$PY evaluate_fibrotic_state_transfer.py --tendon-report results/tendon_multimodal_report.json --external-manifest results/gse226374_observations.jsonl --output results/fibrotic_state_transfer_report.json
$PY evaluate_calcium_uncertainty.py --train-manifest results/mechanotransduction_observations.jsonl --external-manifest results/eahy_calcium_observations.jsonl --output results/calcium_uncertainty_report.json
$PYREADR_PY build_multisite_migration_manifest.py --source ../../../../data/external_training/experiment_factory/outer_library/downloads/multisite_live_cell/analysis_full/'Data processing and analysis-2D'/R/dat_v2.R --output results/multisite_migration_observations.jsonl
$PY evaluate_multisite_migration.py --manifest results/multisite_migration_observations.jsonl --output results/multisite_migration_report.json
$PY build_electrotaxis_piv_manifest.py --control ../../../../data/external_training/experiment_factory/outer_library/downloads/electrotaxis_piv/control_tissue_piv.mat --stimulated ../../../../data/external_training/experiment_factory/outer_library/downloads/electrotaxis_piv/stimulated_tissue_piv.mat --output results/electrotaxis_piv_observations.jsonl
$PY evaluate_electrotaxis_piv.py --manifest results/electrotaxis_piv_observations.jsonl --output results/electrotaxis_piv_report.json
$PY build_dryad_migration_manifest.py --source-zip ../../../../data/external_training/experiment_factory/outer_library/downloads/dryad_9jh6m/Shafqat-Abbasi_Source_Datasets.zip --output results/dryad_migration_observations.jsonl --tensor-output results/dryad_migration_rows.npz
$PY evaluate_dryad_migration.py --tensor results/dryad_migration_rows.npz --output results/dryad_migration_report.json
$NODE export_workbook.mjs ../../../../data/external_training/experiment_factory/outer_library/downloads/elife_71032/elife-71032-fig4-data1-v2.xlsx ../../../../data/external_training/experiment_factory/outer_library/elife_71032_figure4.json
$PY build_elife_71032_manifest.py --workbook-export ../../../../data/external_training/experiment_factory/outer_library/elife_71032_figure4.json --output results/elife_71032_observations.jsonl
$PY evaluate_rock_context_transfer.py --dryad-manifest results/dryad_migration_observations.jsonl --elife-manifest results/elife_71032_observations.jsonl --output results/rock_context_transfer_report.json
$NODE export_workbook.mjs <extracted Figure2/control_dVdt_dAdt.xlsx> ../../../../data/external_training/experiment_factory/outer_library/elife_72381_control_dVdt_dAdt.json
$NODE export_workbook.mjs <extracted Figure2/Y27_dVdt_dAdt.xlsx> ../../../../data/external_training/experiment_factory/outer_library/elife_72381_Y27_dVdt_dAdt.json
$NODE export_workbook.mjs <extracted Figure3/tether_force_all.xlsx> ../../../../data/external_training/experiment_factory/outer_library/elife_72381_tether_force_all.json
$PY build_elife_72381_mechano_osmotic_manifest.py --figure2-source-zip ../../../../data/external_training/experiment_factory/outer_library/downloads/elife_72381/elife-72381-fig2-data1-v2.zip --figure3-source-zip ../../../../data/external_training/experiment_factory/outer_library/downloads/elife_72381/elife-72381-fig3-data1-v2.zip --control-spreading-export ../../../../data/external_training/experiment_factory/outer_library/elife_72381_control_dVdt_dAdt.json --y27-spreading-export ../../../../data/external_training/experiment_factory/outer_library/elife_72381_Y27_dVdt_dAdt.json --tether-export ../../../../data/external_training/experiment_factory/outer_library/elife_72381_tether_force_all.json --output results/elife_72381_mechano_osmotic_observations.jsonl
$PY evaluate_mechano_osmotic_evidence.py --manifest results/elife_72381_mechano_osmotic_observations.jsonl --data-sum-export ../../../../data/external_training/experiment_factory/outer_library/elife_72381_data_sum.json --output results/mechano_osmotic_evidence_report.json
$PY download_zenodo_tfm_figure_tables.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/zenodo_7432971 --receipt results/zenodo_7432971_figure_tables_receipt.json
$NODE export_workbooks.mjs ../../../../data/corpus/external_training/experiment_factory/outer_library/zenodo_7432971_figure_tables.json ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/zenodo_7432971/Kreysing_McHugh_key_resources.xlsx ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/zenodo_7432971/Figure2_Data.xlsx ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/zenodo_7432971/Figure3_Data.xlsx ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/zenodo_7432971/FigureS2_Data.xlsx
$PY build_zenodo_tfm_tension_manifest.py --workbook-export ../../../../data/corpus/external_training/experiment_factory/outer_library/zenodo_7432971_figure_tables.json --receipt results/zenodo_7432971_figure_tables_receipt.json --output results/zenodo_7432971_tfm_tether_observations.jsonl --report results/zenodo_7432971_manifest_report.json
$PY evaluate_zenodo_tfm_tension.py --manifest results/zenodo_7432971_tfm_tether_observations.jsonl --output results/zenodo_7432971_tfm_tether_evidence.json
$PY download_zenodo_tfm_fret_source.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/zenodo_14692589 --receipt results/zenodo_14692589_acquisition.json
$PY build_zenodo_tfm_fret_manifest.py --source-zip ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/zenodo_14692589/Source-Data.zip --receipt results/zenodo_14692589_acquisition.json --output results/zenodo_14692589_tfm_fret_observations.jsonl --report results/zenodo_14692589_manifest_report.json
$PY evaluate_zenodo_tfm_fret.py --manifest results/zenodo_14692589_tfm_fret_observations.jsonl --manifest-report results/zenodo_14692589_manifest_report.json --reference-manifest results/zenodo_7432971_tfm_tether_observations.jsonl --output results/zenodo_14692589_tfm_fret_evidence.json
$PY download_figshare_exif_emt_subset.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/figshare_exif_emt --receipt results/figshare_exif_emt_acquisition.json --workers 8 --raw-workers 8 --field-ids 0000
$PY build_figshare_exif_emt_manifest.py --raw-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/figshare_exif_emt/raw_exif_emt --receipt results/figshare_exif_emt_acquisition.json --manifest results/figshare_exif_emt_observations.jsonl --tensor results/figshare_exif_emt_rows.npz --output results/figshare_exif_emt_manifest_report.json
$PY evaluate_figshare_exif_emt.py --tensor results/figshare_exif_emt_rows.npz --manifest-report results/figshare_exif_emt_manifest_report.json --output results/figshare_exif_emt_evaluation.json
$PY download_mtrack_emt_pilot.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/mtrack_emt_xy01 --receipt results/mtrack_emt_acquisition.json --workers 6
$PY build_mtrack_emt_manifest.py --source-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/mtrack_emt_xy01 --receipt results/mtrack_emt_acquisition.json --manifest results/mtrack_emt_observations.jsonl --tensor results/mtrack_emt_rows.npz --output results/mtrack_emt_manifest_report.json
$PY evaluate_mtrack_emt.py --tensor results/mtrack_emt_rows.npz --manifest-report results/mtrack_emt_manifest_report.json --exif-report results/figshare_exif_emt_evaluation.json --output results/mtrack_emt_evaluation.json
$PY download_lincs_mcf10a_phenotypes.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/lincs_mcf10a --receipt results/lincs_mcf10a_acquisition.json
$PY build_lincs_mcf10a_phenotypes.py --source-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/lincs_mcf10a --receipt results/lincs_mcf10a_acquisition.json --manifest results/lincs_mcf10a_observations.jsonl --tensor results/lincs_mcf10a_rows.npz --output results/lincs_mcf10a_manifest_report.json
$PY evaluate_lincs_mcf10a_phenotypes.py --tensor results/lincs_mcf10a_rows.npz --manifest-report results/lincs_mcf10a_manifest_report.json --output results/lincs_mcf10a_evaluation.json
$PY download_karacosta_emt_cytof.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/karacosta_emt_cytof --receipt results/karacosta_emt_cytof_acquisition.json
$PY build_karacosta_emt_cytof.py --source-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/karacosta_emt_cytof --receipt results/karacosta_emt_cytof_acquisition.json --manifest results/karacosta_emt_cytof_observations.jsonl --tensor results/karacosta_emt_cytof_rows.npz --output results/karacosta_emt_cytof_manifest_report.json
$PY evaluate_karacosta_emt_cytof.py --tensor results/karacosta_emt_cytof_rows.npz --manifest-report results/karacosta_emt_cytof_manifest_report.json --exif-report results/figshare_exif_emt_evaluation.json --output results/karacosta_emt_cytof_evaluation.json
$PY download_geo_emt_rna.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/geo_emt_rna --receipt results/geo_emt_rna_acquisition.json
$PY build_geo_emt_rna.py --source-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/geo_emt_rna --receipt results/geo_emt_rna_acquisition.json --manifest results/geo_emt_rna_observations.jsonl --tensor results/geo_emt_rna_rows.npz --output results/geo_emt_rna_manifest_report.json
$PY evaluate_geo_emt_rna.py --tensor results/geo_emt_rna_rows.npz --manifest-report results/geo_emt_rna_manifest_report.json --protein-report results/karacosta_emt_cytof_evaluation.json --output results/geo_emt_rna_evaluation.json
$PY download_gse325309_pristine_holdout.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/gse325309 --protocol gse325309_preregistered_holdout_protocol.json --receipt results/gse325309_pristine_acquisition.json
$PY build_gse325309_pristine_holdout.py --source ../../../../data/corpus/external_training/experiment_factory/outer_library/gse325309/GSE325309_processed_data.xlsx --receipt results/gse325309_pristine_acquisition.json --protocol gse325309_preregistered_holdout_protocol.json --manifest results/gse325309_pristine_observations.jsonl --tensor results/gse325309_pristine_rows.npz --output results/gse325309_pristine_manifest_report.json
$PY evaluate_gse325309_pristine_holdout.py --tensor results/gse325309_pristine_rows.npz --manifest-report results/gse325309_pristine_manifest_report.json --output results/gse325309_pristine_evaluation.json
$PY audit_observation_operators.py --registry observation_operator_registry.json --repo-root ../../../.. --output results/observation_operator_audit.json
$PY audit_dryad_supracontractility_sources.py --source-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/dryad_08kprr59c --output results/dryad_08kprr59c_source_audit.json
$NODE export_workbooks.mjs ../../../../data/corpus/external_training/experiment_factory/outer_library/dryad_08kprr59c_workbooks.json ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/dryad_08kprr59c/Data_Fig2.xlsx ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/dryad_08kprr59c/Data_Fig4.xlsx
$PY build_dryad_supracontractility_manifest.py --workbook-export ../../../../data/corpus/external_training/experiment_factory/outer_library/dryad_08kprr59c_workbooks.json --source-audit results/dryad_08kprr59c_source_audit.json --output results/dryad_08kprr59c_supracontractility_observations.jsonl --report results/dryad_08kprr59c_manifest_report.json
$PY evaluate_dryad_supracontractility.py --manifest results/dryad_08kprr59c_supracontractility_observations.jsonl --output results/dryad_08kprr59c_supracontractility_evidence.json
$PY build_cell_monolayer_velocity_manifest.py --source-dir ../../../../data/external_training/experiment_factory/outer_library/downloads/cell_monolayer_velocity --output results/cell_monolayer_velocity_observations.jsonl --tensor-output results/cell_monolayer_velocity_vectors.npz
$PY evaluate_cell_monolayer_velocity.py --manifest results/cell_monolayer_velocity_observations.jsonl --tensor results/cell_monolayer_velocity_vectors.npz --output results/cell_monolayer_velocity_report.json
$PY train_calcium_temporal_cnn.py --train-manifest results/mechanotransduction_observations.jsonl --external-manifest results/eahy_calcium_observations.jsonl --checkpoint results/calcium_temporal_cnn_checkpoint.npz --output results/calcium_temporal_cnn_report.json
$PY predict_calcium_temporal_cnn.py --checkpoint results/calcium_temporal_cnn_checkpoint.npz --manifest results/eahy_calcium_observations.jsonl --condition 'LooPINS (+GsMTx4)' --condition 'LooPINS (+Mag)' --condition 'CaPINS (+GsMTx4)' --condition 'CaPINS (+Mag)' --output results/calcium_temporal_cnn_predictions.jsonl
$PY build_bbbc_translocation_if_manifest.py --bbbc013-zip ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc013/BBBC013_v1_images_bmp.zip --bbbc013-platemap ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc013/BBBC013_v1_platemap_all.txt --bbbc014-zip ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc014/BBBC014_v1_images.zip --bbbc014-platemap ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc014/BBBC014_v1_platemap_all.txt --output results/bbbc_translocation_if_observations.jsonl --tensor-output results/bbbc_translocation_if_rows.npz
$PY evaluate_bbbc_translocation_if.py --tensor results/bbbc_translocation_if_rows.npz --output results/bbbc_translocation_if_report.json
$PY build_hpa_if_embedding_tensor.py --archive ../../../../data/external_training/experiment_factory/outer_library/downloads/hpa_v25_1/subcell_image_umap_features.tsv.zip --output results/hpa_if_embedding_rows.npz --receipt results/hpa_if_embedding_receipt.json
$PY evaluate_hpa_if_localization.py --tensor results/hpa_if_embedding_rows.npz --output results/hpa_if_localization_report.json
$PY audit_outer_library_coverage.py --results results --output results/adversarial_coverage_audit.json
$PY build_hpa_raw_image_subset.py --tensor results/hpa_if_embedding_rows.npz --manifest results/hpa_raw_if_subset_manifest.json --download-dir ../../../../data/external_training/experiment_factory/outer_library/downloads/hpa_v25_1/raw_if_subset --receipt results/hpa_raw_if_download_receipt.json --workers 12
$PY download_bbbc054.py --output-dir ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc054 --receipt results/bbbc054_acquisition_receipt.json --segments 16 --workers 24
$PY build_bbbc054_microglia_manifest.py --replicate1-zip ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc054/Replicate_1.zip --annotation ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc054/Replicate1annotation.csv --output results/bbbc054_microglia_observations.jsonl --tensor-output results/bbbc054_microglia_rows.npz
$PY evaluate_bbbc054_microglia.py --tensor results/bbbc054_microglia_rows.npz --output results/bbbc054_microglia_report.json
$PY download_pbmc_cross_lab.py --output-dir data/corpus/external_training/experiment_factory/outer_library/downloads/pbmc_cross_lab --protocol pbmc_cross_lab_preregistered_protocol.json --receipt results/pbmc_cross_lab_acquisition.json
$PY build_pbmc_cross_lab_tensor.py --data-dir data/corpus/external_training/experiment_factory/outer_library/downloads/pbmc_cross_lab --receipt results/pbmc_cross_lab_acquisition.json --protocol pbmc_cross_lab_preregistered_protocol.json --output results/pbmc_cross_lab_tensor.npz --report results/pbmc_cross_lab_tensor_report.json
$PY train_pbmc_cross_lab_model.py --tensor results/pbmc_cross_lab_tensor.npz --protocol pbmc_cross_lab_preregistered_protocol.json --model results/pbmc_cross_lab_model.npz --report results/pbmc_cross_lab_model_report.json
$PY evaluate_pbmc_ifn_direction.py --data-dir data/corpus/external_training/experiment_factory/outer_library/downloads/pbmc_cross_lab --protocol pbmc_cross_lab_preregistered_protocol.json --soffice "$SOFFICE" --output results/pbmc_ifn_direction_report.json
$PY download_pbmc_two_lab_recalibration.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/pbmc_two_lab --protocol pbmc_two_lab_recalibration_preregistered_protocol.json --receipt results/pbmc_two_lab_acquisition.json
$PY calibrate_pbmc_on_gse164378.py --data-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/pbmc_two_lab --receipt results/pbmc_two_lab_acquisition.json --protocol pbmc_two_lab_recalibration_preregistered_protocol.json --source-model results/pbmc_cross_lab_model.npz --output-model results/pbmc_gse164378_calibrated_model.npz --output-report results/pbmc_gse164378_calibration_report.json
$PY download_wilk_pbmc_confirmation.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/wilk_pbmc_confirmation --protocol pbmc_final_confirmation_preregistered_protocol.json --receipt results/wilk_pbmc_confirmation_acquisition.json
# evaluate_wilk_pbmc_confirmation.py additionally needs h5py; install it in an ignored task-local dependency directory and prepend that directory to PYTHONPATH.
$PY evaluate_wilk_pbmc_confirmation.py --h5ad ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/wilk_pbmc_confirmation/GSE150728_Wilk_CELLxGENE.h5ad --receipt results/wilk_pbmc_confirmation_acquisition.json --protocol pbmc_final_confirmation_preregistered_protocol.json --model results/pbmc_gse164378_calibrated_model.npz --calibration-report results/pbmc_gse164378_calibration_report.json --output-report results/wilk_pbmc_confirmation_report.json --output-tensor results/wilk_pbmc_confirmation_tensor.npz
$PY evaluate_multitask_outer_model.py --calcium-report results/calcium_baseline_report.json --piezo1-report results/piezo1_transfer_report.json --calcium-uncertainty-report results/calcium_uncertainty_report.json --mechanism-report results/mechanism_contrast_report.json --migration-report results/multisite_migration_report.json --dryad-migration-report results/dryad_migration_report.json --fibrotic-report results/fibrotic_state_transfer_report.json --mechano-osmotic-report results/mechano_osmotic_evidence_report.json --piv-report results/electrotaxis_piv_report.json --cell-monolayer-report results/cell_monolayer_velocity_report.json --temporal-cnn-report results/calcium_temporal_cnn_report.json --if-translocation-report results/bbbc_translocation_if_report.json --hpa-if-report results/hpa_if_localization_report.json --microglia-report results/bbbc054_microglia_report.json --cell-cycle-report results/bbbc048_mlp_report.json --mitochondrial-stress-report results/bbbc053_mito_stress_report.json --sciplex-design-report results/sciplex3_design_audit.json --apoptosis-report results/biad2515_apoptosis_report.json --sciplex-state-report results/sciplex3_state_model_report.json --senscout-report results/senscout_morphology_model_report.json --figshare-cell-death-report results/figshare_cell_death_model_report.json --external-senescence-report results/external_senescence_model_report.json --senescence-confirmation-report results/gse301164_postfreeze_confirmation.json --tfm-tether-report results/zenodo_7432971_tfm_tether_evidence.json --piv-contractility-report results/dryad_08kprr59c_supracontractility_evidence.json --nematic-tfm-report results/dataverse_data2772_cellular_nematic_evidence.json --opto-force-report results/dryad_sj3tx9683_force_propagation_evidence.json --tfm-fret-report results/zenodo_14692589_tfm_fret_evidence.json --exif-emt-report results/figshare_exif_emt_evaluation.json --mtrack-emt-report results/mtrack_emt_evaluation.json --lincs-mcf10a-report results/lincs_mcf10a_evaluation.json --karacosta-emt-report results/karacosta_emt_cytof_evaluation.json --geo-emt-rna-report results/geo_emt_rna_evaluation.json --idr0072-if-report results/idr0072_if_evaluation.json --output results/multitask_outer_model_report.json
# The multitask command above must also include: --gse325309-emt-report results/gse325309_pristine_evaluation.json --pbmc-cross-lab-report results/pbmc_cross_lab_model_report.json --pbmc-ifn-direction-report results/pbmc_ifn_direction_report.json --pbmc-final-confirmation-report results/wilk_pbmc_confirmation_report.json --lightmycells-report results/lightmycells_confirmation_evaluation.json
# The sweep-gate command must include: --observation-operator-report results/observation_operator_audit.json
$PY regenerate_summary.py --results results --output results/summary.json
$PY download_bbbc048.py --output-dir ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc048 --receipt results/bbbc048_acquisition_receipt.json
$PY build_bbbc048_cell_cycle_manifest.py --archive ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc048/nested/CellCycle.zip --ground-truth ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc048/Ground_truth.lst --output results/bbbc048_cell_cycle_observations.jsonl --receipt results/bbbc048_manifest_receipt.json
$SYSTEM_PY build_bbbc048_image_tensor.py --manifest results/bbbc048_cell_cycle_observations.jsonl --archive ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc048/nested/CellCycle.zip --output results/bbbc048_image_rows.npz --receipt results/bbbc048_image_tensor_receipt.json
$PY train_bbbc048_mlp.py --tensor results/bbbc048_image_rows.npz --checkpoint results/bbbc048_mlp_checkpoint.npz --output results/bbbc048_mlp_report.json
$PY download_sciplex3_metadata.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/sciplex3 --receipt results/sciplex3_metadata_receipt.json
$PY audit_sciplex3_design.py --pdata ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/sciplex3/GSM4150378_sciPlex3_pData.txt.gz --cell-annotations ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/sciplex3/GSM4150378_sciPlex3_A549_MCF7_K562_screen_cell.annotations.txt.gz --matrix-receipt results/sciplex3_matrix_receipt.json --output results/sciplex3_design_audit.json
$PY build_reactome_state_panel.py --gene-annotations ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/sciplex3/GSM4150378_sciPlex3_A549_MCF7_K562_screen_gene.annotations.txt.gz --output reactome_state_panel.json
$PY build_sciplex3_state_tensor.py --matrix ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/sciplex3/GSM4150378_sciPlex3_A549_MCF7_K562_screen_UMI.count.matrix.gz --pdata ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/sciplex3/GSM4150378_sciPlex3_pData.txt.gz --panel reactome_state_panel.json --tensor results/sciplex3_state_rows.npz --output results/sciplex3_state_tensor_report.json
$PY train_sciplex3_state_model.py --tensor results/sciplex3_state_rows.npz --checkpoint results/sciplex3_state_model.npz --output results/sciplex3_state_model_report.json
$PY download_bbbc053.py --output-dir ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc053 --receipt results/bbbc053_acquisition_receipt.json
$SYSTEM_PY build_bbbc053_mito_stress.py --archive ../../../../data/external_training/experiment_factory/outer_library/downloads/bbbc053/FCCP-20220127T153841Z-001.zip --manifest results/bbbc053_mito_stress_observations.jsonl --tensor results/bbbc053_mito_stress_rows.npz --output results/bbbc053_mito_stress_report.json
$PY download_biad2515_apoptosis.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/biad2515 --receipt results/biad2515_acquisition_receipt.json
$PY build_biad2515_apoptosis.py --source-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/biad2515 --manifest results/biad2515_apoptosis_observations.jsonl --tensor results/biad2515_apoptosis_rows.npz --output results/biad2515_apoptosis_tensor_report.json
$PY evaluate_biad2515_apoptosis.py --tensor results/biad2515_apoptosis_rows.npz --checkpoint results/biad2515_apoptosis_ridge.npz --output results/biad2515_apoptosis_report.json
$PY download_figshare_cell_death.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/figshare_28202864 --receipt results/figshare_cell_death_acquisition.json
$PY download_figshare_cell_death.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/figshare_28202864 --receipt results/figshare_cell_death_acquisition.json --features --segments 384 --workers 8
$PY audit_figshare_cell_death_design.py --source-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/figshare_28202864 --output results/figshare_cell_death_design_audit.json
$PY extract_figshare_cell_death_aggregates.py --archive ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/figshare_28202864/features.tar.gz --acquisition-receipt results/figshare_cell_death_acquisition.json --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/figshare_28202864/aggregates --receipt results/figshare_cell_death_aggregate_receipt.json
$PARQUET_PY build_figshare_cell_death_tensor.py --aggregate-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/figshare_28202864/aggregates --extraction-receipt results/figshare_cell_death_aggregate_receipt.json --design-audit results/figshare_cell_death_design_audit.json --output results/figshare_cell_death_rows.npz --report results/figshare_cell_death_tensor_report.json
$PY train_figshare_cell_death_model.py --tensor results/figshare_cell_death_rows.npz --checkpoint results/figshare_cell_death_model.npz --output results/figshare_cell_death_model_report.json
$PY audit_gse250041_senescence_design.py --source-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/gse250041 --output results/gse250041_senescence_design_audit.json
$PY download_external_senescence_panels.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/external_senescence --receipt results/external_senescence_acquisition.json
$PY build_external_senescence_tensor.py --source-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/external_senescence --soffice "$SOFFICE" --output results/external_senescence_rows.npz --report results/external_senescence_tensor_report.json
$PY train_external_senescence_model.py --tensor results/external_senescence_rows.npz --checkpoint results/external_senescence_model.npz --output results/external_senescence_model_report.json
$PY evaluate_gse301164_confirmation.py --counts ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/external_senescence/GSE301164_BJ_RNA_seq_raw_counts.txt.gz --checkpoint results/external_senescence_model.npz --output results/gse301164_postfreeze_confirmation.json
$PY download_dataverse_cellular_nematics.py --output-dir ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/dataverse_data2772 --receipt results/dataverse_data2772_acquisition.json
$PY build_dataverse_cellular_nematics_manifest.py --source-root ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/dataverse_data2772 --receipt results/dataverse_data2772_acquisition.json --output results/dataverse_data2772_cellular_nematic_observations.jsonl --report results/dataverse_data2772_manifest_report.json --tensor-output results/dataverse_data2772_piv_vectors.npz
$PY evaluate_dataverse_cellular_nematics.py --manifest results/dataverse_data2772_cellular_nematic_observations.jsonl --manifest-report results/dataverse_data2772_manifest_report.json --output results/dataverse_data2772_cellular_nematic_evidence.json
$PY audit_dryad_force_propagation_source.py --source-root ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/dryad_sj3tx9683 --output results/dryad_sj3tx9683_source_audit.json
$PY build_dryad_force_propagation_manifest.py --source-root ../../../../data/corpus/external_training/experiment_factory/outer_library/downloads/dryad_sj3tx9683 --source-audit results/dryad_sj3tx9683_source_audit.json --output results/dryad_sj3tx9683_force_propagation_observations.jsonl --report results/dryad_sj3tx9683_manifest_report.json
$PY evaluate_dryad_force_propagation.py --manifest results/dryad_sj3tx9683_force_propagation_observations.jsonl --manifest-report results/dryad_sj3tx9683_manifest_report.json --source-audit results/dryad_sj3tx9683_source_audit.json --output results/dryad_sj3tx9683_force_propagation_evidence.json
$PY evaluate_theory_runtime.py --registry theory_registry.json --output results/theory_runtime_report.json
$PY evaluate_external_inference.py --output results/external_inference_report.json
$PY evaluate_latent_promotion.py
```

`PYREADR_PY` denotes a Python environment containing `pyreadr`; all evaluation
steps remain dependency-light and run in the Aleph environment.
`$NODE` is the bundled workspace Node runtime used with `@oai/artifact-tool`.
`$SYSTEM_PY` denotes `/usr/bin/python3`, whose Pillow decoder is used only to
turn verified BBBC048 JPEG members into the content-hashed compact tensor.
`$PARQUET_PY` is an isolated Python runtime containing PyArrow for the one-time
verified Parquet-to-NPZ conversion. `$SOFFICE` is the bundled LibreOffice
binary used to read the original legacy GSE63577 OLE spreadsheet; the converted
header and counts sheet were rendered and visually checked before tensor use.

See `REPORT.md` for measured capability and failure boundaries.
