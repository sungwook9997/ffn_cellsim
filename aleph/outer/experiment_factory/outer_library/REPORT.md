# Outer Library implementation report — 2026-08-06

## Cortical-tension comparison preparation (2026-08-10)

The Outer Library now carries a structured, protocol-conditioned reference for
the engine's next cortical-tension comparison. The admitted numeric source is
Hosseini 2020 for MCF-7 suspended/rounded interphase cells at 37 C under dynamic
AFM parallel-plate confinement. The bundle records the total effective Laplace
surface-tension summary as 270 pN/um median and 180--400 pN/um empirical IQR
for 27 cells. Its audit resolves the record to verified `KB-6.1.7`, `SE181`,
and `source_audit.verdict=OK`. `VG-γ-MCF7-suspended` remains draft and is not
scored by this package.

The engine adapter is fail-closed. It exposes no engine magnitude or ratio
unless every record establishes full-native CUDA execution with all
compartments, physical accepted steps, time-step independence, a stationary
total method-of-planes tension, a non-blocked quantitative claim, exact cell/
state/protocol identity, a native-validated geometry-matched observation
operator, and between-seed scatter from three distinct records. The current
stored cortex record is retained only as a refusal control: its values are
redacted because it is force-accepted, quantitatively blocked, lacks the
comparison context, and the AFM/Laplace operator mapping is only
`candidate_native`. MDA-MB-231 and adherent MCF-7 are explicitly unavailable;
no cross-state or cross-cell-line value is imputed.

Artifacts and the reproducible command live under
`cortical_tension/README.md`, with the audited JSON and one-row machine-readable
CSV beside the other Outer Library results. This head has no Aleph parameter,
gate, or physics-mutation authority and required no GPU execution.

## Outcome

The previous literature-candidate system has been extended into a real public
experimental-data path. Thirty-eight datasets including the separate HPA embedding
source across thirty-seven conservative lab groups now produce 193,832 traceable
samples and 329,303 canonical observations,
plus compact raw-row and paired-image feature tensors. This is enough to run genuine
single-cell cell-state experiments; it is not yet enough to authorize an Aleph
parameter sweep.

The current objective was deliberately changed because the present Aleph engine
cannot yet execute the required parameter-sweep volume. This is more than a
change in execution order. The current deliverable is an independently useful
external biophysical inference network that can infer observable state, predict
held-out measurements, compare compatible mechanisms, calibrate uncertainty,
refuse out-of-domain inputs, and recommend a discriminating assay without
depending on an Aleph sweep. Its success or failure must be judged on those
independent endpoints now.

This change does not discard a future Aleph connection. If broad sweeps become
executable, the validated posterior and ambiguity set may later be tested as a
pre-sweep constraint. That is a conditional future interface, not current Aleph
parameter authority and not the reason the present model may claim success.
Matched observation operators, identifiability, and Aleph forward sensitivity
must pass before it may reduce candidate axes or bound sweep ranges. Until then,
it must not emit numeric Aleph parameters or silently remove a sweep axis. The
scope correction is recorded in `objective_scope_clarification_2026-08-07.json`;
the already-frozen training preregistration is retained unchanged as an audit
artifact.

The newest force branch adds two independent laboratories and two cell types.
In the 2026 Science/Dataverse NIH3T3 source, the Figure S9 tables reproduce the
published myosin-inhibition direction exactly at the field level: mean traction
is 136.8858 before and 30.1780 Pa after para-nitroblebbistatin; mean maximum
principal tension is 25.4025 before and 5.1826 mN/m after. All 12 fields move
in the same direction, but the source path exposes one date/culture, so this is
one biological group with 12 technical fields. The undocumented role of `t2`
is not guessed: `t0/t1` reproduce the paper histogram and `t2` is withheld.

The eLife/Dryad opto-MDCK source adds the inverse perturbation: global RhoA
activation has positive author-derived relative surface-tension increase in
47/68 x-axis and 46/68 y-axis objects. Its released aggregate files contain a
material inconsistency—`fullstim_data.csv` has 22/30/16 objects in AR1to1d,
AR1to2d, and AR2to1d, while `stats_fullstim.csv` reports 35/15 objects in
AR1to1d/AR1to1s. The latter is excluded, nine pickle members are never
executed, and the 68 rows stay one split group because sample/date IDs are not
released.

Together with the already replicated C2C12 PIV/IF response, these sources form
a three-lab directional mechanism-family constraint. They rank actomyosin,
stress-fibre, cortex, and adhesion candidates for a future Aleph sweep; they do
not distinguish those axes numerically. The sweep gate therefore still emits
`numeric_range: null` and requires matched PIV/TFM/MSM/contour observation
operators, aligned biological replication, dataset/lab uncertainty and OOD,
identifiability, and monotonic Aleph forward sensitivity.

## Pre-sweep theory and promotion layer

The first executable registry contains three dependency-light experts:
quasi-static membrane tether force, single-mode standard-linear-solid
relaxation, and isotropic 2D persistent-random-walk MSD. Machine checks measured
a tether inversion relative error of 5.10e-16, zero single-mode equation
residual, and the expected ballistic/diffusive persistent-walk limits. Invalid
domains return typed refusals.

For the worked AFM tether example, 18.697398 pN identifies the product
`kappa*sigma_app = 4.427644115 pN^2`. Fixing kappa at 0.082 pN*um conditionally
gives `sigma_app = 53.9956484 pN/um`, but does not identify its bilayer-tension
and membrane-cortex-attachment components. Consequently the current promotion
report is `external_training_or_validation_only`: independent-dataset and
independent-lab holdouts, calibrated uncertainty, ambiguity resolution, a
geometry-matched Aleph observation operator, and forward sensitivity are all
still false. No numeric sweep range is emitted.

The execution path toward the unchanged sweep objective is explicit:

```text
external experiment/theory inference (independently useful now)
  -> independent dataset + lab validation and calibrated uncertainty
  -> identifiability audit
  -> exact-unit, geometry-matched Aleph observation operator
  -> Aleph forward-sensitivity calibration
  -> candidate-axis reduction and bounded sweep proposal for physics review
  -> executable Aleph sweep when the engine cost gate is met
```

BBBC048 now supplies the previously missing cell-cycle axis: 32,266 Jurkat
cells with three channels and seven exact author phases. The archive and nested
image ZIP are CRC-clean. The label distribution is extremely imbalanced (only
15 anaphase, 27 telophase, and 68 metaphase cells), so it is representation and
cell-cycle supervision from one provider, not an external validation set.
The first raw-image head uses disjoint train/selection/calibration/test cells
and achieves same-provider test macro F1 0.614, balanced accuracy 0.651, and
accuracy 0.793. Its calibrated 90% conformal set attains 0.895 marginal
coverage with mean size 1.28 of seven. Metaphase has only 12 test cells,
anaphase four, and telophase six; conditional rare-phase claims are therefore
refused even where point recall looks favorable.

For broader perturbational state, sci-Plex3 metadata and pData are verified for
799,317 cells. Exactly 762,795 have complete cell-type/replicate/time/treatment/
dose context and 36,522 do not; the latter are barred from supervised training.
The complete cells form 4,896 mandatory context groups (median 131 cells,
range 2--4,448). Replicate holdout remains same-lab validation. The official
3,096,224,606-byte UMI matrix now passes its gzip stream, 1,007,419,688 sparse
line, and SHA-256 checks (`7d632716...5d00a`), so the head is training-ready.

The expression feature basis is no longer a manually selected marker list.
Nine official Reactome human pathways supply all reference gene products for
cell cycle, apoptosis, senescence, TGF-beta EMT, hypoxia, DNA repair, UPR,
innate immunity, and interferon signaling. Exact intersection with the
sci-Plex annotation yields 2,228 genes; every API payload is SHA-256 locked and
there are zero manual additions. The matrix builder streams sparse
triplets and produces a 4,896 x 2,228 context tensor, rather than materializing
the impossible 110,983 x 799,317 dense cell matrix. Drug identities, including
all doses and replicates, remain disjoint across model splits; vehicle is a
reference and never a shared supervised target.

Once the matrix receipt passes size, gzip-stream, sparse-line-count, and SHA-256
checks, the design audit becomes training-ready from that receipt rather than a
hand-edited flag. The multitask contract then accepts the trained state head
only if treatment overlap is zero and sealed compounds were absent from fit,
selection, and calibration. Even a successful same-provider test remains
training-only and cannot enter the Aleph evidence or parameter channels.

The first treatment-disjoint model is now complete, and it fails rather than
being cosmetically promoted. Fit-only cell-type/replicate/time residualization
improves a 768-gene, 96-hidden-unit multitask MLP over its raw-expression
baseline, but it still reaches only 0.198 macro F1 and 0.429 accuracy over 21
sealed compounds (12 evaluable pathway classes). Its nominal 90% conformal set
covers 0.952 only by returning 11.52 of 16 classes on average, and the
vehicle OOD diagnostic rejects 0/8 contexts. Replicate latent cosine similarity
is 0.746 across 2,448 matched contexts, but that is same-provider consistency.
The tensor is a valid training asset; this pathway head is an uninformative
baseline that needs a stronger perturbational objective and external study
before it can represent general cell state.

BBBC053 adds a distinct mitochondrial-stress imaging axis with 29 DMSO and 29
FCCP TOM20 fields from murine CAD cells. All 58 source images are 2048x2048,
16-bit, super-resolution TIFF and have been decoded into content-hashed compact
images plus 13 morphology features. It covers FCCP mitochondrial uncoupling,
not hypoxia, and no biological-replicate identifiers are published; it is
therefore training-only and does not close the independent hypoxia-state gap.

S-BIAD2515 adds an author-labelled apoptosis endpoint with a stronger split
hierarchy than the other new state sets. The official CC0 BioImage Archive
indexes 6,240 raw TIFFs (45,085,609,920 bytes) from four-channel HeLa live-cell
imaging followed by terminal Annexin-V. The author analysis commit contributes
2,336 live morphology features and 530 terminal features for 30 wells: ten
staurosporine doses with three wells per dose. Terminal features are excluded
from model inputs. One well at every dose is used for fitting, one for
calibration, and the authors' third well is sealed for testing. A fit-only
dual-ridge model reaches sealed same-provider R2 0.878, MAE 0.0983, and
positive-versus-other AUROC 0.88; it beats the dose-only MAE baseline 0.175.
The ten-well calibration interval covers 9/10 sealed wells. This establishes a
useful training endpoint, not external dataset/lab authority.

An independent 2025 Uppsala/SciLifeLab Cell Painting corpus (Figshare
28202864 v2, CC BY 4.0) is now verified at the design-file level. It contains
MCF7 cells exposed for 48 hours to compounds spanning apoptosis, autophagy,
ferroptosis, immunogenic cell death, necroptosis, and pyroptosis across eight
plates. The 21,271-site metadata table, 15,430-row QC table, and all three
author well-split files pass their published MD5 values. This is independent
provider evidence for a broader programmed-cell-death representation, but it
is not the same endpoint as S-BIAD2515 Annexin-V and is therefore never counted
as a direct external regression test. The 2.451 GB feature archive remains a
required, separately verified input before its mechanism head can be trained.
The published well splits are not suitable for unseen-drug generalization: 45
to 47 compounds cross train/validation/test roles depending on feature
extractor. The fixed Aleph audit instead assigns all 55 labelled compounds to
40 train, four selection, five calibration, and six sealed-test identities
with zero chemical overlap. The six test compounds cover every death class;
rare pyroptosis and necroptosis classes deliberately lack full conditional
calibration authority.

The next direct state axis is now source-verified before its matrices are
admitted. NIH/NIA GSE250041 supplies 8,664 untreated proliferating and 4,949
10-Gy IR/day-10 senescent WI-38 CITE-seq barcodes. Both conditions use the
same 36,601 gene-expression plus eight antibody-capture features. This is an
independent dataset and laboratory relative to sci-Plex, but each condition is
one pooled capture; single cells cannot manufacture biological replication.
Both count matrices now pass exact-size, gzip, MatrixMarket-shape, and SHA-256
verification. A compact 13,613 x 43 tensor retains 35 source-declared RNA
markers and all eight antibody-capture channels using log1p CPM/CLR
normalization. Because capture ID is perfectly confounded with state, it has
zero valid dataset or lab splits: this removes the raw-matrix acquisition gap,
but it does not remove the independent senescence-validation gap.

SenSCOUT adds a morphology-to-senescence path without treating author-derived
clusters as truth. Its CC0 Dryad archive is exact-size, SHA-256, CRC, and
path-safety verified. The model uses the authors' exact 87-feature curated
morphology panel and excludes biomarker intensities, UMAP, KMEANS, labels, and
file identifiers from its inputs. Bio1 alone supplies fit/selection/calibration
with image groups kept intact; Bio2 is a sealed same-provider repeat and Bio3
is a frozen shift diagnostic. On 377 confidently labelled Bio2 BGAL cells the
cell-level diagnostic balanced accuracy is 0.659, but across the 36 image-group
prevalences the Bio1-calibrated 90% interval covers only 0.444. The direct BGAL
head therefore fails uncertainty transfer. P16 morphology association is more
stable (Bio2/Bio3 Pearson r 0.770/0.808), while HMGB1 and LMNB1 shift sharply in
Bio3 and P21 has negative external-repeat R2. All heads remain training-only;
no same-provider result is counted as an independent dataset/lab holdout or an
Aleph sweep constraint.

| Selected modality | Observations |
|---|---:|
| Single-cell calcium time series | 97,722 |
| Immunofluorescence-derived | 7,438 |
| Multisite high-content live-cell features | 2,916 |
| Multiscale tracked-cell organizational features | 21,388 |
| Author-supervised migration-mode summaries | 948 |
| 3D CDM per-cell migration trajectory observables | 106 |
| AFM membrane tether force | 314 |
| FXm/RICM spreading and volume-rate observables | 630 |
| Particle-image-velocimetry-derived | 1,212 |
| PIV and live-actin-derived C2C12 | 1,435 |
| Fibronectin-orientation IF-derived | 720 |
| Traction-force-derived | 2,523 |
| FLIM-FRET vinculin-tension-sensor-derived | 266 |
| Paired focal-adhesion TFM/FLIM-FRET-derived | 492 |
| Pa/mN-per-m paired TFM/MSM-derived | 120 |
| TFM/MSM/actin contour-model-derived | 476 |
| Single-cell FRET time series | 1,452 |
| Counting-derived | 72 |
| Cantilever-post force-derived | 89 |
| Nuclear morphology-derived | 800 |
| RT-qPCR-derived | 210 |
| Bulk RNA-seq marker panel | 40 |
| IF pixel-derived | 17 |
| BBBC013 paired-channel FKHR-EGFP/DRAQ IF field features | 1,152 |
| BBBC014 paired-channel NF-kB FITC/DAPI IF field features | 1,152 |
| BBBC054 brightfield microglia cell patches | 58,186 |
| BBBC048 three-channel imaging-flow cell-cycle labels | 32,266 |
| BBBC053 TOM20 super-resolution mitochondrial-stress images | 58 |
| S-BIAD2515 live morphology to terminal Annexin-V endpoints | 30 |
| **All canonical observations (`results/summary.json`)** | **236,913** |

Sources include Figshare 16826740 and 27241950 plus the official open Zenodo
records [18495219](https://zenodo.org/records/18495219),
[8215150](https://zenodo.org/records/8215150),
[19890985](https://zenodo.org/records/19890985), the 2025 ACTG1/occludin image
set, the 2026 tendon mechanoculture source workbook, and the three-lab
[MULTIMOT 2D live-cell dataset](https://doi.org/10.17044/scilifelab.21407402).
Those sources are CC
BY 4.0; all downloads are content-hashed and recorded in
`acquisition_receipts.jsonl`.
The 1,610 Figshare 16826740 TFM observations no longer carry guessed or unknown
units: the linked paper's publisher Figure 2E axes identify projected area as
µm², major/minor lengths as µm, and total/x/y force as nN. Every normalized row
records the article-content hash and the exact Figure 2E graphic locator.
The independent fibrotic-state confirmation is the public NCBI GEO record
[GSE226374](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226374),
whose record carries no explicit reusable-content licence; its provenance and
licence state remain distinct from the CC BY sources.
The CC0 [Dryad 9jh6m source dataset](https://doi.org/10.5061/dryad.9jh6m),
linked to eLife DOI 10.7554/eLife.11384, contributes fibronectin, ROCK, and
talin perturbation series from H1299-PL cells. It overlaps the Strömblad/KI
research group represented in MULTIMOT and is therefore same-lab expansion,
not an independent fourth lab.
The independent 2022 [eLife 71032 study](https://doi.org/10.7554/eLife.71032)
adds per-cell NIH3T3 migration speed and oscillation period in 3D
cell-derived matrix. Its Figure 4 workbook is CC BY 4.0 and its published MD5
and local SHA-256 are verified. The authors do not expose independent
experiment identifiers for these cells, which limits its confirmatory status.
The independent 2022 [eLife 72381 study](https://doi.org/10.7554/eLife.72381)
adds HeLa Kyoto spreading/volume dynamics and single-cell AFM tether-force
measurements. Its Figure 2 and Figure 3 ZIPs are CC BY 4.0 and exact-file
SHA-256 receipts are committed. Unlike the 3D migration workbook, these source
tables retain experiment labels that reproduce the authors' N=3--9 hierarchy.
The official [BBBC013](https://bbbc.broadinstitute.org/BBBC013) and
[BBBC014](https://bbbc.broadinstitute.org/BBBC014) benchmark plates contribute
192 intact two-channel IF fields under CC BY 3.0. Their image ZIPs, text plate
maps, and BBBC014 XLS map are byte-counted, SHA-256 pinned, and ZIP-CRC checked.
They share Ilya Ravkin provider lineage and are therefore one conservative lab
group even though acquisition sites, instruments, assays, and cell types differ.
The official Human Protein Atlas v25.1 subcellular image-embedding archive adds
81,007 IF images covering 13,366 genes, 15,700 antibodies, 39 cell lines, and
35 author localization labels. The 632,515,574-byte CC BY 4.0 ZIP is CRC-clean
and pinned by SHA-256. It is represented separately from canonical experimental
observations because provider image features are training inputs, not physical
observables.

## Measured model results

The first HeLa three-condition calcium baseline uses 150 complete single-cell
curves with a 90/30/30 cell-grouped split and zero cell overlap.

| Test | Result |
|---|---:|
| Accuracy | 0.567 |
| Macro F1 | 0.553 |
| NLL | 1.133 |
| ECE (10-bin) | 0.216 |
| Semantic leave-one-condition-out OOD AUROC | 0.596 |
| Independent-lab MDCK-II OOD AUROC | 0.910 |

Thus geometric OOD detects the distinct MDCK-II dataset well, while unseen
conditions inside the HeLa domain and probability calibration remain weak.

The independent-lab task transfer trains a binary GsMTx4/Piezo1-inhibition
signature on 100 HeLa cells, then evaluates the entire independent EA.hy926 lab
dataset without mixing samples. Of 625 external cells, 584 map to the predeclared
`+Mag` versus `+GsMTx4` mechanism-level task.

| External-lab transfer | Result |
|---|---:|
| Samples | 584 |
| Accuracy | 0.545 |
| Balanced accuracy | 0.518 |
| Macro F1 | 0.508 |

This legacy broad-arm task is near chance. A second, predeclared aligned task
keeps only the author-described active LooPINS/CaPINS strata and uses a
baseline-, duration-, and amplitude-invariant representation. Feature form,
ridge penalty, and threshold are selected on a deterministic HeLa-only split;
EA.hy926 labels are never used for selection.

| Aligned external-lab task | Result |
|---|---:|
| Samples | 311 |
| AUROC | 0.762 |
| Balanced accuracy | 0.726 |
| Macro F1 | 0.550 |

Discrimination improved materially, but the task still fails because macro F1
is below 0.70 under the untouched, highly imbalanced external prevalence.

A separate uncertainty audit reserves 30 HeLa cells that are used neither for
model fitting nor penalty/temperature selection. Pooled split conformal fails
under external class imbalance (coverage 0.666). Class-conditional Mondrian
split conformal, still calibrated only on HeLa, restores untouched EA.hy926
coverage to 0.923 against the 0.90 target. It emits singleton sets for 51.4% of
cells at 0.850 singleton accuracy; class-wise coverage is 0.981 and 0.911. A
training-centroid dataset-shift score separates HeLa calibration from EA.hy926
at AUROC 0.732, with 37.9% of external cells beyond the HeLa 95th percentile.
The uncertainty holdout now passes, while the pooled failure and detected shift
remain in the report and the per-cell macro-F1 failure remains independent.

### Three-lab migration-state transfer and uncertainty

The official MULTIMOT author object contains 472,173 cell-time rows from three
labs, three people per lab, three independent experiments per person, three
technical replicates per condition, and matched control/ROCK-inhibition arms.
Cells and frames are deliberately not counted as independent samples. The
builder uses the paper's exact 18-variable complete-case panel and emits 2,916
technical-replicate medians: 18 features for each of 162 assay replicates. The
inference units are the 27 independent experiments nested within nine people;
confidence intervals resample the three person clusters in each lab.

The initial committed protocol learned a paired-control perturbation direction
in lab 1, selected a threshold in lab 2, and opened lab 3 once. It produced
AUROC 0.985, balanced accuracy 0.889, and macro F1 0.888, but left no disjoint
lab-2 samples for uncertainty calibration. The audited revision therefore
removes threshold fitting: paired-control centering declares zero as the
no-effect origin, lab 2 is used for validation and class-conditional conformal
calibration, and no lab-3 label enters numeric fitting.

Because lab 3 had already been opened under the initial revision, the revised
results below are a retrospective external-lab evaluation, not a new pristine
confirmation. Matched-control centering is part of the assay and the model is
not claimed to work without a control arm.

| Revised external lab 3 result | Value |
|---|---:|
| Technical-replicate samples | 54 |
| AUROC | 0.985 |
| Balanced accuracy | 0.741 |
| Macro F1 | 0.722 |
| Independent experiments | 9 |
| Person-cluster bootstrap 95% CI | [1.780, 2.600] |
| Mondrian conformal coverage | 1.000 |
| Mean prediction-set size | 1.352 |
| Singleton fraction / accuracy | 0.648 / 1.000 |

Raw features almost perfectly identify lab 3 versus lab 1 (OOD AUROC 0.971).
After matched-control centering the same diagnostic falls to 0.347, while all
nine external-lab independent experiment projections remain positive. The
operational external metrics and retrospective uncertainty targets pass, but
`pristine_sealed_confirmation=false` and
`publication_confirmatory_holdout_passed=false` until a new lab is acquired.
This is not a claim that ROCK inhibition identifies one physical parameter;
`aleph_parameter_proposal_eligible` remains false.

### Representative electrotaxis PIV

Both official Zenodo MAT files now pass their published byte counts and MD5
checks. The adapter produces eight spatial field summaries for each of 60
control and 60 stimulated frames (960 observations total). Frames retain one
shared, explicitly unknown technical-replicate identity rather than being
mislabelled as 120 replicates.

| Descriptive frame mean | Stimulated minus control |
|---|---:|
| Mean speed | +0.552 |
| Mean x velocity | +0.939 |
| RMS speed | +0.666 |
| Polar order | +0.125 |
| Nematic order | +0.031 |

The released files contain one representative tissue per condition and do not
declare the velocity unit. Spatial vectors and time frames are correlated
measurements, not biological replicates. Consequently the evaluator emits no
inferential CI, no external holdout, no observable constraint, and no Aleph
parameter proposal. These data are usable for representation pretraining and
descriptive operator checks, not for a treatment-effect claim.

### Same-lab multiscale migration representation transfer

The verified Dryad archive contains 22,647 repeated cell-time observations
from 367 tracked H1299-PL cells across 33 author-reported experimental dates.
It provides 60 image-derived behavioral and organizational features under
fibronectin, ROCK/Y27632, and talin-siRNA perturbations. The canonical manifest
emits 22,336 cell-level medians and author-mode fractions; a 7.2 MB compressed
tensor retains all 1,282,754 finite feature values without converting time
observations into independent biological replicates.

For the migration-mode task, the five behavioral features that define the
author Continuous/Discontinuous labels are excluded. A class-balanced ridge
model uses 55 organizational features. Fibronectin experimental dates are
fixed into disjoint train, validation, and conformal-calibration groups using
only the fibronectin source; talin labels do not choose the split, penalty, or
threshold. All talin dates and tracked cells are then opened as a retrospective
same-lab, different-perturbation test.

| Organizational-state transfer | Result |
|---|---:|
| Validation observation AUROC / balanced accuracy / macro F1 | 0.961 / 0.884 / 0.884 |
| Talin observation AUROC / balanced accuracy / macro F1 | 0.967 / 0.890 / 0.888 |
| Talin tracked-cell AUROC / balanced accuracy / macro F1 | 0.981 / 0.907 / 0.901 |
| Talin experimental-date bootstrap balanced-accuracy 95% CI | [0.854, 0.947] |
| Cell-level Mondrian coverage | 0.967 |
| Prediction-set mean size | 1.546 |
| Singleton fraction / accuracy | 0.454 / 0.928 |

This passes the declared same-lab cross-perturbation representation and
uncertainty checks. It is evidence that the organizational state is not merely
the five behavioral label inputs. It does not satisfy the missing requirement:
`independent_lab_holdout=false`, `pristine_lab4_confirmation=false`, and
`aleph_parameter_proposal_eligible=false` remain hard-coded outputs.

### Independent 3D migration context and ROCK direction reversal

The eLife source provides paired speed/period columns for 24 untreated NIH3T3
cells and 9 Y-27632 cells, plus three other actomyosin/lamellipodia
perturbations. Y-27632 changes persistent speed by -0.367 µm/min (Hedges
g=-0.937) and oscillation period by +9.036 min in this 3D CDM assay. The
cell-bootstrap speed interval is [-0.553, -0.194], but it remains descriptive
because independent experiment IDs are absent.

This is not interchangeable with the Dryad H1299-PL 2D result. Across four
Dryad dates containing both DMSO and Y-27632 arms, the mean date-level speed
effect is positive (+7,785.5 author source units), while its date bootstrap
interval [-1,463.4, 15,321.1] includes zero. The two observed directions are
opposite. The library therefore records a context-heterogeneity warning rather
than inventing a universal “ROCK inhibition decreases/increases speed” rule.
This independent dataset expands context-diversity training, but
`independent_lab_mechanism_holdout_passed=false` and
`observable_constraint_eligible=false`: it cannot validate the full MULTIMOT
18-feature signature or select one Aleph parameter.

### Experiment-grouped mechano-osmotic and membrane-force evidence

The eLife 72381 Figure 2 tables contribute 194 control cells from three
independent experiments and 121 Y-27632 cells from four independent
experiments. Each cell retains paired initial spreading-area rate and volume
flux; treatment confidence intervals resample experiment means, not cells.
Their units are recovered from the author `data_sum.xlsx` headers, specifically
Sheet1 B1 (`dA/dt i (µm2/min)`) and E1 (`dV/dt i (µm3/min)`); all 630 spreading
observations carry the Figure 2 archive hash and the appropriate cell locator.
The source control table contains 195 finite rows, but the published Figure 2D
summary declares n=194. The builder excludes the one exact negative-spreading
row whose removal reproduces both published control means, and records that
rule explicitly rather than silently changing the source.
For Y-27632, the 121-cell raw volume-flux mean exactly matches `data_sum.xlsx`,
but the raw area-rate mean is 34.853 versus the summary's 34.177 µm²/min
(difference +0.676). The manifest preserves the per-cell raw workbook and the
report records this mismatch; it does not overwrite measurements to force the
aggregate table to agree.

| Y-27632 minus control | Experiment-level estimate | 95% bootstrap CI |
|---|---:|---:|
| Initial spreading-area rate (µm²/min) | +25.360 | [14.166, 37.591] |
| Initial volume flux (µm³/min) | -11.837 | [-17.491, -7.304] |

Figure 3 contributes 314 AFM per-cell tether-force points with force reported
in pN. File dates and experiment suffixes reproduce every published N. For the
three dates containing matched vehicle and Y-27632 arms during the 30--90
minute spreading phase, all three experiment effects are positive: +8.866,
+9.605, and +8.638 pN. Their mean is +9.036 pN with experiment-bootstrap CI
[8.638, 9.605]. Under the paper's constant-bending-rigidity assumption, the
matched mean-force ratio corresponds descriptively to a 2.200
apparent-membrane-tension ratio.

The detected effect is transient rather than universal. At 4--5 hours in
spread cells, the matched experiment effect is -0.067 pN with CI
[-0.625, 0.700]; on round
micropatterned cells it is +0.249 pN with CI [-2.221, 2.058]. This temporal
structure is now part of the evidence constraint; CI overlap with zero is not
claimed as an equivalence test.

The gate deliberately rejects the shortcut “tether force = cortical tension.”
AFM tether force also depends on bending rigidity, membrane-cortex attachment,
and tether geometry, while Y-27632 changes several cellular mechanisms. The
new Helfrich surface-tension candidate therefore remains
`numeric_range=null`, `may_emit_sweep_axis=false`, and blocked until an
AFM-tether observation operator, joint rigidity treatment, time-stratified
forward sensitivity, and an independent-lab confirmation exist.

### Human dermal fibroblast collective velocity representation

Dryad dataset `10.25349/D96W58` contributes a CC0, author-labelled collective
motion source from human dermal fibroblasts on isotropic and nematic liquid
crystal elastomer substrates. Eleven source files were checked against the
record's official byte sizes and MD5 values and independently pinned by
SHA-256. The normalized result contains 3,154 observations in 832 technical
samples. A compact tensor retains all 7,645 Figure 5 vectors with columns
`vx`, `vy`, `px`, `py`, and velocity angle. One stationary vector has the
author spreadsheet value `#DIV/0!` for angle; it is preserved as `NaN` and is
accepted only because both velocity components are exactly zero.

The spatial fields reproduce the expected organizational contrast:

| Substrate | Technical vectors | Mean speed (source unit) | Polar order | Nematic order |
|---|---:|---:|---:|---:|
| Isotropic | 2,944 | 0.002336 | 0.03478 | 0.05485 |
| Nematic | 4,701 | 0.002015 | 0.02903 | 0.48994 |

The longitudinal 1 µM FAKi arm has lower timepoint-mean density than its
untreated nematic comparison (221.773 versus 382.401 cells/mm²), lower
cell-substrate order (0.2597 versus 0.4723), and slightly lower mean speed in
the undeclared source velocity unit (11.492 versus 12.245). These are
descriptions, not causal estimates: the arms differ strongly in density, the
release identifies no independent cultures, and spatial vectors/time points
are correlated technical measurements. Consequently every inferential CI is
`null`; external holdout, uncertainty holdout, observable constraint, and
Aleph parameter proposal are all `false`.

The sweep bridge records this source as `admitted_training_only`. All 7,645
vectors share mandatory split group `valentine-lab-ucsb`, so no spatial point,
frame, or time point can leak across folds. This adds useful active-nematic and
FAK-context representation data without manufacturing replication or a sweep
axis.

### Paired-channel IF translocation dataset holdout

BBBC013 supplies one 640x640 FKHR-EGFP signal image and one DRAQ DNA image for
each of 96 U2OS wells. Four plate rows cover each of the author-labelled
Wortmannin and LY294002 dose series. BBBC014 supplies one 1360x1024 FITC NF-kB
signal image and one DAPI image for each of 96 wells: four MCF7 rows and four
A549 rows across twelve author TNF-alpha concentrations. The adapter decodes
the 8-bit BMP payloads without Pillow, validates all 384 archive members, and
pairs both channels before any split or feature is emitted.

The sample unit is one intact two-channel well field. Otsu DNA foreground is
cleaned and component-filtered only to construct an image-derived nuclear mask;
no component is promoted to a labelled cell. Twelve observables capture signal
intensity, nuclear-mask morphology, nuclear/perinuclear signal enrichment, and
DNA-signal association. This yields 192 tensor rows and 2,304 traceable
observations. The same nuclear/perinuclear log-ratio was fixed using BBBC013
Wortmannin rows A--C; BBBC014 was then opened once as the dataset holdout.

| Context | Fields | Dose-vs-score Spearman | Extreme-dose AUROC | Technical-row bootstrap 95% CI |
|---|---:|---:|---:|---:|
| BBBC013 Wortmannin | 48 | 0.920 | 1.000 | [0.891, 0.941] |
| BBBC013 LY294002 compound transfer | 48 | 0.927 | 1.000 | [0.913, 0.946] |
| BBBC014 MCF7 sealed dataset holdout | 48 | 0.667 | 1.000 | [0.631, 0.746] |
| BBBC014 A549 sealed dataset holdout | 48 | 0.839 | 1.000 | [0.786, 0.872] |

The representation cleanly separates the predeclared high/low dose extremes
in both held-out cell types and the shrunk-Mahalanobis dataset-shift detector
separates BBBC014 at AUROC 1.000. It nevertheless fails the declared transfer
gate: MCF7 Spearman is below 0.70. A fixed-ridge split-conformal audit trained
on BBBC013 rows A--C and calibrated on row D covers 83.3% of BBBC014 against a
90% target and has normalized-dose interval width 0.739, above the 0.60 width
limit. Neither threshold is relaxed after seeing the holdout.

Therefore this path is admitted for grouped IF representation training only.
`observable_evidence_eligible=false`, `uncertainty_holdout_passed=false`, and
`independent_lab_holdout=false`. The four author replica rows are technical
plate replicas rather than reported biological replicates, and shared provider
lineage prevents the two datasets from being counted as independent labs.

### HPA subcellular IF localization representation

The HPA archive provides 1,024 neural image features per IF image. A fixed
every-eighth feature view retains 128 dimensions without inspecting labels or
holdout performance. Splitting individual images would leak the same target or
reagent, so the split unit is the connected component of the gene--antibody
bipartite graph. This yields 51,976 train, 12,229 validation, 8,026 calibration,
and 8,776 sealed grouped-test images across 12,540 components, with exactly zero
gene or antibody crossing a split.

| Sealed grouped HPA test | Result |
|---|---:|
| Evaluable author location labels | 33 |
| Macro AUROC | 0.950 |
| Macro average precision | 0.585 |
| Macro F1 | 0.527 |
| Micro F1 | 0.787 |
| Exact multilabel match | 0.548 |
| Marginal calibration coverage (nominal 90%) | 0.898 |

The high rank discrimination demonstrates useful IF representation signal, but
rare-label macro F1 and exact-set prediction remain modest. More importantly,
the grouped test is still one HPA provider and one annotation pipeline. The
calibration result is diagnostic only because repeated images within reagent
components are not independent experiments. This ninth head is training-only,
has no independent-lab status, and cannot select an Aleph parameter.

### Adversarial axis audit

The canonical corpus contains 329,303 observation rows and 193,832 unique
dataset/sample pairs: 1.70 rows per sample. Time points, channels, and derived
observables are therefore never advertised as independent training subjects.
The audit finds 29,055 observations (8.8%) with unknown, author-processed, or
otherwise unspecified units and 238,071 (72.3%) without a reported biological
replicate identifier. These rows can train representations but cannot silently
become inferential replicates.

The existing `cell_state` field is also not yet a controlled biological
ontology: its values mix perturbations, migration modes, assay phases, and
phenotypes. Covered task axes include mechanosensitive calcium, context-specific
ROCK/migration, fibrotic condition contrasts, mechano-osmotic response, IF
localization, cell cycle, mitochondrial stress, one single-replicate LPS
microglia series, the same-provider S-BIAD2515 apoptosis endpoint, cross-lab
directional senescence, and replicated C2C12 supracellular contractility/ECM
alignment. Absolute senescence state remains OOD-blocked rather than promoted
to a production classifier. EMT now has fixed-IF treatment coverage and a
distinct-lab live temporal diagnostic, but not a common-task controlled holdout.
Explicitly missing axes are hypoxia, differentiation/stemness, DNA damage, and
infection response. Missing orthogonal modalities include WB,
verified single-cell RNA-seq expression (metadata is present while the matrix
transfer remains incomplete), spatial transcriptomics,
proteomics/phosphoproteomics, metabolomics, conventional flow/CyTOF, and EM.

The acquisition queue therefore prioritizes a new untouched raw-IF provider for
post-revision confirmation, independent cell-cycle/apoptosis and replicate-complete controlled EMT panels, donor-aware
single-cell transcriptomics, replicate-complete raw TFM/PIV, and protein-level
confirmation. A deterministic HPA pilot has also selected 512 images covering
all 35 localization labels and 29 cell lines. It downloads blue nucleus, red
microtubule, green target-protein, and yellow ER channels separately, with at
most one selected image per gene--antibody component and the original grouped
split preserved.
All 2,048 channel JPEGs were then downloaded from the official HPA image host:
410,923,719 bytes total. Every payload passes JPEG boundary validation and is
individually size- and SHA-256-pinned in the generated receipt.

### Cross-provider raw-IF transfer result

The HPA raw-image checkpoint was frozen before the final IDR0072 pixel
evaluation. Source resolution, optical channel order, and field selection were
audited separately; one malformed three-image field was refused before
prediction, leaving 2,593 evaluated confocal fields across 18 landmark proteins,
81 plates, two species screens, and six supported localization classes. This is
a real external-provider transfer diagnostic, but schema and optics discovery
mean it is not a pristine final confirmation.

Every frozen endpoint failed:

| Endpoint | Result | Required | Status |
|---|---:|---:|---|
| image macro average precision | 0.3169 | >= 0.45 | **fail** |
| image macro F1 | 0.3057 | >= 0.45 | **fail** |
| multilabel ECE | 0.2370 | <= 0.20 | **fail** |
| positive conformal coverage | 0.8056 | >= 0.85 | **fail** |
| mean prediction-set size | 5.723/9 | <= 4.0 | **fail** |

Human-screen macro AP was 0.4911, but the larger mouse screen reached only
0.2944; this species/provider asymmetry cannot be averaged into acceptance.
Landmark-plus-plate aggregation raised macro AP to 0.3939 but macro F1 remained
0.2805, and plate fields are technical measurements rather than new biological
replicates. The 29th multitask head is therefore
`blocked_external_domain_shift`: it has no observable-evidence, uncertainty,
production, Aleph, or numeric-range authority. IDR0072 is now consumed and may
not be used to tune or validate a replacement. A revised domain-robust IF model
must use other providers for fitting/selection and reserve another untouched
provider for single-use confirmation.

### Light My Cells spatial observation branch

The independent-network objective also requires a visual encoder that preserves
spatial measurement content rather than reducing every IF field to a class.
The 2026 Light My Cells release (BioImage Archive `S-BIAD1047`, CC BY 4.0)
provides the appropriate task: registered BF/PC/DIC stacks paired with
author-labelled nuclear, mitochondrial, tubulin, or actin fluorescence.

The protocol was committed before bulk FileList access. A deterministic hash
assigns the 30 independent Study Components to 16 fit, four selection, four
calibration, and six single-use confirmation studies. The complete metadata
audit exactly reproduces the paper totals: 56,984 files, 2,574 acquisition
sets, 41,213 BF planes, 7,670 PC planes, 3,499 DIC planes, 2,533 nuclear
targets, 1,819 mitochondrial targets, 223 tubulin targets, and 27 actin targets.
The listed source volume is 109,287,061,290 bytes and no filename is unmatched.

The sealed confirmation partition contains 102 acquisition sets across six
studies, with 102 nuclear and 91 mitochondrial targets. Tubulin (23) and actin
(10) fall below the frozen 25-set support threshold and are descriptive only.
The planned
head predicts fluorescence maps and uncertainty from transmitted light, with
study-level spatial correlation, SSIM, MAE improvement, pixel-interval
coverage, and embedding-OOD endpoints. It is an independent visual-observation
head; passing it would not create general cell-state or Aleph parameter
authority.

The first OME-header audit found that 2,782 files omit `Plane.PositionZ` and one
target is not a valid TIFF. That diagnostic was committed before defining a
fallback. The frozen v2 resolver uses physical-position alignment when present,
the only plane for singleton inputs, and otherwise the lower median z plane;
the latter is explicitly marked as approximate rather than exact registration.
It recovers 4,375 of 4,376 possible fit/selection/calibration pairs (99.977%):
3,253 exact-position, 57 singleton, and 1,065 fallback pairs. The one malformed
Study_25 nuclear target remains refused. The resulting 7,933 unique payloads
total 12,678,835,256 bytes, 11.6% of the full source listing. Confirmation
headers and pixels remain untouched.

Full payload acquisition exposed one further provider inconsistency:
`Study_28/image_2492_Nucleus.ome.tiff` is listed as 6,108,577 bytes, while the
official URL serves a 4,096-byte metadata-only BigTIFF with no image pixels.
The source-integrity protocol was committed before filtering; the pair is
refused without replacement. The final non-confirmation asset contains 4,374
pairs and 7,932 unique TIFFs totaling 12,672,726,679 bytes, each locally
SHA-256-pinned. Robust per-image percentile normalization and BOX reduction
produce a finite `[4374,64,64]` float16 input/target tensor spanning 24 studies
and 2,472 acquisition sets. Two inputs and one target are constant after source
decoding and remain explicit audit cases rather than silently deleted.

The first conditional U-Net comparison uses only 3,744 fit pairs for gradient
updates, 61 pairs from four studies for architecture/epoch selection, and 569
pairs from four other studies for calibration. An adversarial audit corrected
the initial pair-pooled selection calculation before any checkpoint was
committed; the final selector gives every study equal weight. It chooses a
26,641-parameter, 16-base-channel model at epoch 9.

This v1 model is not ready for the sealed studies. On Nucleus/Mitochondria
selection aggregates it reaches Pearson 0.1984, global SSIM 0.0699, and an
8.59% macro-study MAE improvement over the fit-only target mean, missing all
three provisional readiness targets (0.35, 0.25, and 10%). The conservative
study-cluster intervals cover 95.09% of mitochondrial and 95.23% of nuclear
pixels, slightly exceeding the frozen 95% upper bound, while normalized widths
are a broad 0.494 and 0.531. The checkpoint is retained as a reproducible
development diagnostic, but all external-confirmation, uncertainty,
production, and Aleph authority remain false. Confirmation data stay sealed
for a stronger spatial model.

The preregistered v2 model addresses the blur failure with a two-level residual
U-Net and an endpoint-aligned L1/correlation/global-SSIM/gradient loss. It has
244,001 parameters and uses deterministic joint flips on fit rows only. Epoch
selection maximizes the weakest of the three normalized selection endpoints,
not any single favorable metric. Epoch 16 is selected with Pearson 0.5580,
global SSIM 0.4755, and 30.09% macro-study MAE improvement; the weakest endpoint
is still 1.594 times its required value. Nuclear and mitochondrial calibration
coverage is 94.38% and 94.10%, within the frozen 85–95% band, although interval
widths remain broad at 0.514 and 0.492. The exact checkpoint, calibration
thresholds, and fit-embedding OOD reference are now frozen. This satisfies the
precondition to open the six single-use confirmation studies, but by itself
does not grant external, uncertainty, production, or Aleph authority.

The confirmation partition was then opened exactly once against that frozen
checkpoint. Source-integrity checks admit 386 TIFF payloads totaling
1,536,818,955 bytes with no size inconsistency, yielding 226 pairs from 102
acquisition sets. The six studies span five experimenter groups, but three
groups overlap training and only CRSA Paris and IGDR Bretagne-Loire are new;
this is therefore a complete study holdout and only a two-group partial
laboratory holdout.

The final result fails the preregistered all-endpoints rule. Required-target
macro-study Pearson is 0.3292 against 0.35, while study-level fit-reference
embedding-OOD AUROC is 0.3854 against 0.70. Global SSIM 0.2746, relative MAE
improvement 11.92%, and pixel-interval coverage 92.31% pass their thresholds.
Mitochondrial Pearson/SSIM are 0.3466/0.3007 and nuclear values are
0.3119/0.2485; the two novel experimenter groups are weaker at 0.3030/0.2538
for mitochondria and 0.2616/0.1903 for nuclei. The 30th multitask head is thus
`blocked_external_spatial_transfer`. It has no evidence, uncertainty,
production, general-cell-state, Aleph, or numeric-sweep authority, and the
consumed confirmation set may not be used to tune its replacement.

A domain-generalized v3 replacement was therefore trained without reopening
that consumed partition. Its 173,937-parameter conditional U-Net combines
instance-normalized residual blocks, deterministic microscope-style views,
prediction consistency, and a study-adversarial embedding. Fit gradients still
use only the original 3,744 pairs. Epoch 19 is selected internally with Pearson
0.5066, global SSIM 0.4727, and 31.23% MAE improvement; mitochondrial and
nuclear calibration coverage are 93.65% and 93.16%. Internal unseen-study OOD
AUROC is only 0.4219, so v3 remains a training-only checkpoint.

An independent non-disease Allen Cell audit then froze 32 normal M0 WTC11
hiPSC TOMM20 cells from 32 FOVs, 15 plates, two instruments, and two workflows.
All 64 raw crop/FOV payloads were acquired with content length, ETag, and SHA-256
receipts. The preregistered exact DNA-control reconstruction failed: the
isotropically resampled FOV DNA differs from the official crop by 15.51 raw
uint16 units on average. Consequently the single-use confirmation is formally
`refused_geometry_exact_match` and cannot grant authority.

A labelled post-refusal diagnostic nevertheless answers whether v3 improved
cross-provider transfer. Against the same Allen fields, v2 reaches Pearson
0.1306, SSIM 0.1022, and -60.29% MAE improvement, whereas v3 reaches 0.3679,
0.3590, and +2.30%. Thus v3 materially improves spatial correlation and
structure but still misses the frozen 10% MAE-improvement requirement. Its
95.84% interval coverage exceeds the 95% upper bound and OOD AUROC remains only
0.4292. The images show the same limitation: v3 recovers coarse mitochondrial
regions better than v2 but not the fine TOMM20 network. The head remains
blocked for production, observable evidence, uncertainty, general cell-state,
Aleph, and numeric sweep use.

BBBC054 then fills the first immune-response training gap with 60 annotated
brightfield fields and 58,186 exact-cell coordinates from replicate 1. The
actual CSV contains 43,813 `round`, 11,143 `ramified`, and 3,230 `amoeboid`
labels. This corrects a source-page inconsistency that says `stratified`; no
label was rewritten. Replicates 2 and 3 contain images but no corresponding
manual labels, so they are not promoted to labelled external tests.

The whole-field split has 36/12/6/6 train/validation/calibration/test fields and
zero field overlap. A 12-feature brightfield patch LDA reaches only 0.502 test
balanced accuracy, 0.402 macro F1, and 0.611 multiclass Brier score. Nominal 90%
split-conformal coverage is 0.905 only by returning 2.45 of 3 classes on
average, so the uncertainty result is not useful. The classifier is rejected.
The descriptive trajectory remains informative: from fields 0--9 to 50--59,
`ramified` fraction falls by 0.214 (field-bootstrap 95% CI
[-0.262, -0.166]) while `round` rises by 0.225 ([0.175, 0.272]). Without a
matched control or another annotated replicate, this is training context, not
causal immune-state evidence.

### Independent TFM/tether pre-sweep constraint

Zenodo 7432971 adds a second laboratory and dataset for traction and membrane-
tether observables. The official archive is 4,942,831,343 bytes; its central
directory (3,408 members) and the CRCs of four selected author XLSX members were
verified by byte-range reads. This avoids downloading 4.94 GB while preserving
member-level hashes and provenance. The source contributes 1,024 normalized
observations from 380 cell-level samples: 648 TFM-derived observations and 376
optical-tweezer observations. Source workbooks were rendered and visually
checked; no formula errors were found.

NIH3T3 mean traction has median 9.812 Pa at 1 kPa substrate stiffness and
32.585 Pa at 10 kPa (30 cells per condition). The direction is positive in
both author file batches, but the 10 kPa medians differ sharply by batch; two
unresolved batches are not enough for inferential replication. Xenopus growth-
cone mean traction has medians 1.274 Pa at 100 Pa and 2.927 Pa at 300 Pa
(33/99 cells). Acquisition date is partly confounded with stiffness, with only
2021-06-03 shared between groups. Optical-tweezer tables do not expose
biological replicate identifiers. All cells are therefore descriptive
measurement units, never biological replicates.

The operator audit preserves two exact unit facts without overclaiming their
semantics. Steady tether force is pN in both the experiment and Aleph candidate
observable, but it jointly constrains bending rigidity, apparent tension, and
membrane-cortex adhesion. Likewise, `1 pN/um2 = 1 Pa`, so Aleph clutch traction
divided by contact area is dimensionally compatible with mean TFM stress, but
the experimental ROI/mask averaging operator is not yet proven equivalent.
Both mappings are blocked. This source contributes to the independent external
mechanical-state posterior now. It may additionally reduce a future Aleph sweep
only after operator and forward-sensitivity validation; that optional adapter
does not define the present objective.

### Independent TFM--FRET load-transfer constraint

Zenodo 14692589 adds an independent 2026 KU Leuven MEF source with 89
cell-level mean-traction measurements in Pa, 266 cell-level VinTS/TSMod
FLIM-FRET measurements, and 246 paired focal adhesions from three representative
cells. The nested GraphPad Prism archives are CRC-clean, and every admitted
value retains its exact nested table/cell locator. In total, 847 observations
are grouped into 358 cell-level samples.

Mean-traction median rises from 163.1344 Pa on 4.5 kPa PAA to 220.2214 Pa on
13 kPa PAA, while VinTS median FRET falls from 11.8803% to 10.7436%, consistent
with increasing vinculin tension. The traction direction agrees with the
independent NIH3T3 1-to-10 kPa dataset in the same Pa unit. This is a genuine
cross-lab same-observable direction test, but not a numeric transfer: cell type
and stiffness range differ, and neither released table provides enough mapped
biological-replicate identifiers.

The three representative focal-adhesion tables reproduce the authors' direct,
inverse, and uncorrelated categories. FRET-on-traction slopes are -0.52775,
+0.16964, and +0.04005, with correlations -0.6652, +0.2899, and +0.0745,
respectively. This directly prevents the model from equating bulk traction with
one molecular-tension value. Two source defects are preserved rather than
silently repaired: the Figure 2b stiff-substrate cell-area column contains 90
values while the caption declares 85, so cell area is withheld; Figure 4 Prism
columns say Pa although the paper and 0--1 range establish per-cell normalized
traction, so those rows use `author_normalized_fraction`.

The resulting focal-adhesion load-transfer request ranks clutch bond stiffness,
anchor stiffness, active tension, and motor stiffness for a future sweep. It
requires matched cell-mask TFM, focal-adhesion traction, relative vinculin-FRET,
and adhesion-geometry operators. `numeric_range` remains null until mapped
replication, calibrated uncertainty/OOD, and Aleph forward sensitivity pass.

### Replicated C2C12 PIV/IF pre-sweep constraint

Dryad 08kprr59c adds an independent CU Boulder/Sandia source for collective
myoblast-to-myotube dynamics on mechanically anisotropic LCNs. The selected
Fig. 2 and Fig. 4 XLSX files and README are CC0; their local SHA-256 digests
exactly match the Dryad API records. The workbooks were rendered before
normalization. Their formatting is semantically active: red cells are author-
excluded stage-motion or frame-link failures and green cells identify maximum
cell density for temporal alignment. An OOXML audit retained all nine green
cells and prevented all 112 red cells from entering the manifest.

The normalized source contains 2,155 observations from 16 author-independent
LCN substrates: 1,435 live-actin/PIV measurements and 720 fibronectin IF angular
frequencies. Time points and multiple observables remain nested technical
measurements, never independent samples. At the nine-replicate maximum-density
time point (36 h after confluence), anisotropic mLCNs versus isotropic iLCNs had
mean order 0.751 versus 0.484, spatial-disorder fraction 0.062 versus 0.196,
and velocity-correlation length 165.43 versus 130.58 µm.

In Fig. 4, two vehicle and two daily 10 µM blebbistatin substrates were
evaluated after the reported delayed divergence (t >= 42 h). Every control
replicate exceeded every treated replicate for the three declared directional
observables. The replicate-mean control-minus-blebbistatin contrasts were
0.2881 for orientation order, 393.19 µm for nematic-correlation length, and
54.67 µm for velocity-correlation length. Fibronectin angular distributions
were separately converted to nematic order, giving replicate values
0.849/0.595 for controls and 0.209/0.172 for treated substrates, a contrast of
0.5316. With only two substrates per condition, the smallest attainable
two-sided exact randomization p-value is 1/3; cluster-bootstrap intervals are
therefore exploratory and cannot establish calibrated uncertainty. Treated
speed is available for only one substrate and is explicitly non-inferential.

This source closes the previous “PIV has no biological replicate IDs or exact
units” acquisition gap, but it does not close the same-task external-lab or
Aleph-operator gaps. The directional latent is supracellular contractility and
reciprocal cell-ECM alignment. It can rank future simulations by predicted
order/correlation direction only after Aleph implements matching live-actin,
nematic-correlation, velocity-correlation, and fibronectin-orientation
operators. Blebbistatin is a multiscale nonmuscle-myosin-II perturbation, not a
measurement of one Aleph parameter, so no numeric sweep range is emitted.

### Task-scoped multimodal model runtime

The Outer Library now has an executable twenty-one-head model contract rather than
one global score that mixes incompatible labels. Every head declares its
modality, model kind, train/validation/calibration/test grouping, external
status, uncertainty status, pristine-confirmation status, and Aleph authority.
The runtime rejects overlapping migration labs, use of external calcium labels
for model selection, and any input report that claims direct Aleph parameter
authority.

| Head | Current status | External evidence | Uncertainty | Pristine |
|---|---|---:|---:|---:|
| Calcium cell state | blocked external classification | no | yes | yes |
| PIEZO1 calcium mechanism | validated external evidence only | yes | no | yes |
| Paired migration state | validated retrospective external evidence only | yes | yes | no |
| Migration-mode representation | same-lab training only | no | no | no |
| Fibrotic-state mechanism | validated external evidence only | yes | no | yes |
| Mechano-osmotic observable | operator-blocked evidence only | no external holdout | no | no |
| PIV/velocity-field representation | training only, no biological replication | no | no | no |
| IF translocation representation | dataset threshold/uncertainty blocked, training only | no | no | yes (dataset only) |
| HPA IF localization representation | gene/antibody-grouped same-provider training only | no | no | no |
| ExIF raw-IF EMT condition representation | single-plate treatment-label training only | no | no | no |
| M-TRACK live VIM-RFP EMT trajectory | one treated trajectory; intensity direction conflict | no | no | no |
| LPS microglia morphophenotype | single-replicate classifier rejected, training only | no | no | no |
| Cell-cycle image representation | same-provider rare-phase-limited training only | no | no | no |
| Mitochondrial-stress IF representation | no-replicate training only | no | no | no |
| sci-Plex perturbation-state representation | unseen-compound same-provider training only | no | no | no |
| Annexin-V endpoint representation | same-provider well holdout training only | no | no | no |
| SenSCOUT morphology senescence | same-provider biorepeat training only | no | no | no |
| Programmed-cell-death mechanism | unseen-compound same-provider rejected | no | no | no |
| Senescence transcriptomic direction | validated external direction only; absolute state OOD-blocked | yes, direction only | no | yes |
| TFM/tether pre-sweep constraint | operator- and replication-blocked training/descriptive only | no | no | no |
| C2C12 PIV/IF contractility constraint | replicated direction; operator and aligned external-lab blocked | no | no | no |
| Cross-cell-type actomyosin force constraint | three-lab direction; same-observable holdout blocked | no | no | no |
| Cross-lab TFM--FRET load-transfer constraint | same-unit direction reproduced; operator and replicate mapping blocked | no | no | no |

The calcium head illustrates why these permissions are separate. Its active-
nanoswitch external AUROC is 0.762 and class-conditional conformal coverage is
0.923, but external macro F1 is only 0.550; therefore it remains blocked as a
per-cell classifier. The paired migration head reaches external macro F1
0.722, balanced accuracy 0.741, AUROC 0.985, and conformal coverage 1.000, but
the test lab had already been opened in an earlier revision, so it is admitted
only as retrospective evidence. PIEZO1 and fibrotic heads are matched-condition
contrasts, not per-cell classifiers.

The programmed-cell-death head uses a verified 2,451,859,048-byte Figshare
feature archive and only its aggregate DeepProfiler profiles. The resulting
954-well, 672-feature tensor contains 55 compounds. A compound-disjoint
40/4/5/6 fit/selection/calibration/test model obtains balanced accuracy 0.500,
macro-F1 0.389, and nominal 90% conformal coverage 0.667 on six sealed
compounds. These results fail acceptance and the 954 wells are not counted as
independent training units.

The frozen 35-marker senescence head uses GSE63577 only for fitting, selection,
and calibration. GSE297406 gives balanced accuracy 0.750 and four of four
positive paired directions, but sign-test `p=0.0625`, conformal coverage 0.125,
and OOD rejection 1.0 prohibit absolute-state promotion. GSE301164 was acquired
after checkpoint freeze and never entered fitting or selection. It confirms all
three replicative-senescence and all three doxorubicin contrasts (`6/6`, combined
one-sided sign-test `p=0.015625`). Its OOD rejection is also 1.0, so the runtime
admits only the paired directional contrast. Absolute state classification and
every numeric Aleph sweep action remain refused pending cross-lab calibration,
an exact observation operator, and forward-sensitivity validation.

A real trainable temporal encoder was also tested rather than assuming that a
larger model would solve transfer. The selected network has eight learned
length-9 temporal filters, ReLU, global mean/max pooling, and a two-class
softmax head: 114 trainable parameters in total. Filter count and L2 penalty
were selected only on deterministic HeLa train/validation cells. The selected
model achieves HeLa validation macro F1 0.950 and AUROC 0.990, but retrospective
EA.hy926 macro F1 falls to 0.529 despite AUROC 0.778 and balanced accuracy
0.711. Nearest-centroid embedding OOD AUROC is only 0.620, so that detector also
fails to warn reliably about the task-transfer failure. Unlabelled label-shift
EM also fails: it estimates only 44.1% positive cells while the retrospective
author labels contain 83.3%, and slightly reduces macro F1. Thus the failure is
not repairable as a simple prevalence shift. The exact-gradient
L-BFGS checkpoint is byte-reproducible and hash-pinned, but
`architecture_accepted=false`; it has no sweep or uncertainty authority.
The checkpoint now carries the preprocessing dimensions, learned filters,
linear head, training-embedding normalization, class centroids, and class
names. A separate read-only inference command emits per-sample probabilities
and embedding distance while hard-coding `aleph_authority=none` and
`may_select_aleph_parameter=false`; inference does not change the rejected
validation status.

There is deliberately no shared latent space yet. The current sources do not
contain enough samples aligned across modalities and labs to train one without
manufacturing correspondence. One calcium temporal encoder has now been
trained and rejected externally; the next stage is replicated raw IF/TFM/PIV
encoders followed by task-scoped late fusion. That is the point at which the
Slurm-allocated 4090 becomes useful. The recorded report and checkpoint hashes
pin the exact runtime result in `results/summary.json`.

### Lab-invariant mechanism contrast

The same data were then evaluated with a predeclared within-study contrast,
rather than absolute per-cell curve classification. The discovery endpoint is
normalized calcium AUC under GsMTx4 minus the matched reference. In the external
lab, LooPINS and CaPINS are the two author-described PIEZO1-targeting active
nanoswitch strata; BINPs and no-particle arms remain controls and were not
relabelled as active actuators.

| Contrast | Mean difference | Hedges g | Bootstrap 95% CI |
|---|---:|---:|---:|
| HeLa discovery | -0.305 | -1.294 | [-0.399, -0.218] |
| EA.hy926 LooPINS | -0.286 | -1.629 | [-0.465, -0.137] |
| EA.hy926 CaPINS | -0.048 | -0.271 | [-0.103, 0.006] |
| EA.hy926 equal-stratum estimate | -0.167 | — | [-0.261, -0.087] |

The external direction is consistent in 2/2 active strata and the equal-stratum
bootstrap interval excludes zero. This passes the mechanism-level observable
holdout. It does **not** rescue the failed per-cell classifier: the CaPINS
stratum alone remains uncertain, biological replicate IDs are unavailable, and
the result is a condition contrast rather than an individual-cell state label.

### Tendon multimodal state evidence

The 2026 tendon-mechanoculture author workbook contributes 1,147 exact-cell
observations: 89 force, 48 IF, 800 nuclear-shape, and 210 RT-qPCR points. No
cross-panel pairing is invented. Rigid minus compliant primary contrasts agree
across three modalities:

| Observable | Mean difference | Bootstrap 95% CI |
|---|---:|---:|
| Traction force/cell (nN/cell) | +111.858 | [90.009, 133.878] |
| alpha-SMA relative IF | +0.854 | [0.464, 1.255] |
| Acta2 relative expression | +4.575 | [3.124, 6.043] |

F-actin changes by +0.070 with CI [-0.048, 0.207] and is not promoted to a
confirmatory constraint. Nuclear circularity changes by -0.096, but author
biological-replicate identities are unavailable, so it remains descriptive.

For fibrosis, force is higher in systemic sclerosis in both tendon fibroblasts
(+3.928 nN/cell, CI [1.285, 6.807]) and the held-out dermal cell type
(+10.623, CI [6.153, 14.697]). This is a successful cell-type holdout inside one
lab, explicitly not an independent-lab holdout.

An additional untouched lab/dataset, GEO GSE226374, supplies two biological
replicates in each of four human dermal fibroblast conditions. A predeclared
ACTA2/CCN2/COL1A1/FN1/THBS1 log-expression signature increases under TGF-beta1
by +0.872 (bootstrap CI [0.770, 0.974]) and is reversed by celastrol by +1.063
when expressed as TGF-beta1 minus TGF-beta1+celastrol (CI [1.007, 1.118]).
Together with the rat tendon alpha-SMA IF discovery contrast, this passes a
cross-species, cross-modality, independent-lab mechanism holdout. It remains a
condition-level state signature, not a single-cell classifier.

## Current objective and future sweep connection

The near-term target was changed because the current Aleph runtime cannot yet
afford broad parameter sweeps. Outer Library is therefore an independent
deliverable now: it must infer cell/assay context, rank compatible mechanisms,
predict unobserved measurements, expose ambiguity, and recommend the next
discriminating experiment without depending on Aleph execution. Its scientific
validity is evaluated on those terms, not deferred until the simulator is fast.

The architecture still preserves a conditional route to a future sweep. Once
Aleph can execute the required studies, validated external posteriors can be
tested as **pre-sweep constraints** by ranking mechanism families and
context-incompatible regions. They are not yet demonstrated sweep compression.
Numeric ranges remain blocked until a matched observation operator, independent
validation, identifiability, and monotonic forward sensitivity exist. Thus the
current objective is genuinely different because the intended sweep cannot yet
run; future sweep compatibility remains an explicit interface rather than a
present capability or present success criterion.

The observation-operator inventory now makes that boundary executable rather
than rhetorical. CDH1 RNA, VIM RNA, E-cadherin protein, and Vimentin protein are
all `unavailable`: Aleph has no callable that emits those assay quantities. Its
vimentin-labelled mechanical population is not a transcript/protein measurement
model. Tissue traction is only `candidate_native` because a per-site exterior
traction readout exists, while TFM geometry, aggregation, and unit validation do
not. Neither status can enter numeric promotion; only a separately validated
native or surrogate operator with pinned validation evidence can do so.

`sweep_gate.py` still returns `refused` for executable Aleph parameter proposals. Presence
of an external holdout is not enough: the per-cell task holdout must pass, macro
F1 must reach 0.70, ECE must be at most 0.15, and semantic OOD AUROC must reach
0.70. The passed mechanism contrast is emitted through a separate
`evidence_constraint_eligible` channel. Tendon evidence now creates five typed,
explicitly confounded candidates spanning SF/arc active tension, NMII force and
attachment kinetics, and both focal-adhesion half-stiffnesses. Every candidate
has `numeric_range=null` and `may_emit_sweep_axis=false`. The
`aleph.outer_library.forward_sensitivity.v1` evaluator requires at least four
axis levels, three converged repeats per level, exact parameter/observable units,
a monotonic response, and full bracketing of the external bootstrap interval
inside the simulated hull before it can emit a non-authoritative range. The
independent-lab state requirement now passes; monotonic Aleph sensitivity and a
unit-matched traction operator remain required. The gate cannot mutate physics or write decision authority
under any outcome.
The multisite migration result enters the evidence-only channel as a third
constraint alongside calcium and fibrotic-state evidence. External metrics and
retrospective conformal targets pass, while a separate pristine-confirmation
check fails. It also carries a pathway-level identifiability warning and cannot
select a parameter, emit a numeric range, or override the still-failed calcium
classifier/calibration checks.
The eLife mechano-osmotic result is a fourth evidence-only constraint. The
forward-sensitivity contract now accepts either traction (`nN_per_cell`) or
membrane-tether plateau force (`pN`) and preserves a typed operator refusal as
evidence instead of coercing it to NaN or a plausible number. Traction remains
blocked while its upstream load-path precondition fails. Tether is numerically
available, but the evaluator refuses every current one-axis card because one
tether force identifies only `sigma_bilayer + W`; Aleph has no executable `W`
axis. The gate therefore records a non-executable attachment-energy axis
request, not a widened cortical/surface-tension range, and cannot mutate physics.

### Raw IF EMT context pilot

The ExIF release (Figshare `10.6084/m9.figshare.26500210.v1`; Nature
Communications `10.1038/s41467-025-59592-7`) adds a provider-independent raw IF
path for A549 EMT context. The acquisition code verifies the official 194 GB
ZIP64 archive by its remote central directory and range-extracts only an admitted
pilot: one field from each of 48 wells and five aligned channels per field. The
plate ledger preserves eight row-specific marker panels, paired columns for
control/EGF/TGF-beta1, the 48-hour exposure, and field/well grouping.

This closes the *missing training coverage* entry for EMT, not the validation
entry. The labels are author treatment conditions (control, weak EGF induction,
strong TGF-beta1 induction), not per-cell EMT ground truth; all fields come from
one plate date and are technical measurements. The paired-column diagnostic
excludes the variable EMT marker from its input, and the marker-direction check
is descriptive only. The resulting model head is training-only and may later
condition an Aleph sweep on EMT context, but it cannot emit a parameter axis,
numeric range, uncertainty claim, OOD claim, or independent HPA localization
holdout.

### Live temporal EMT contradiction audit

The 2020 Science Advances M-TRACK release (`10.1126/sciadv.aba9319`) adds ten
predeclared fluorescence timepoints and their 16-bit segmentation masks from
the paper-linked public GUI example. They are all from one A549 VIM-RFP
TGF-beta 4 ng/mL `xy01` trajectory, so all ten remain one indivisible split
group; the release supplies neither an untreated control nor biological
replicate identifiers. The exact 20 files total 76,575,370 bytes and each has a
local SHA-256 receipt. The dependency-light TIFF and PNG decoders match system
Pillow on shape, label values, and intensity extrema.

The result is intentionally adverse rather than cosmetically harmonized. Mask
area rises from an early-three median of 34,768 to a late-three median of 63,167
pixels (Spearman rho 0.952), but masked VIM-RFP mean falls from 486.600 to
471.477 gray units (rho -0.915). ExIF fixed-antibody IF showed Vimentin rising
under TGF-beta relative to control. These assays and contrasts are not
numerically equivalent, and the M-TRACK sequence has no control, so the
disagreement is not evidence against either paper. It is evidence that raw
Vimentin intensity alone is currently unsafe as an Aleph sweep-conditioning
observable: photobleaching, cell-area dilution, reporter dynamics, and temporal
alignment must be represented first. The 23rd task head is therefore
training-only, emits no numeric range or sweep axis, and acts as a blocker
against premature intensity-only sweep reduction.

### Controlled MCF10A IF collection holdout

The 2022 Communications Biology LINCS MCF10A resource
(`10.1038/s42003-022-03975-9`) supplies an author-released 159-row IF phenotype
table and a 365-row sample-metadata table under CC BY 4.0. Exact file sizes and
SHA-256 digests are pinned. The paper defines C1 and C2 as separate OHSU
collection periods, each with at least three biological replicates. The build
therefore uses `collection + replicate` as the split group and treats repeated
source rows/WellIndex values as technical wells. It produces 1,263 observations
without inflating the seven released IF biological groups into 159 samples.

Using C1 only for fitting, a six-feature TGFB+EGF-versus-EGF nearest-centroid
classifier scores 8/8 at 24 h and 5/6 at 48 h on C2; AUROC is 1.0 at both times.
At 48 h every paired replicate in both collections shows the same three primary
directions: cytoplasmic KRT5 mean decreases (median C1 -54.4, C2 -107.2 author
intensity units), normalized second-neighbor distance increases (C1 +0.299,
C2 +0.245), and mean cells per cluster decreases (C1 -6.959, C2 -4.097).
Paired-group bootstrap intervals exclude zero in the expected direction.

This is the 24th task head and materially improves controlled EMT-state
training, but it is not promoted. C2 is a later collection at the same OHSU
site, not an independent laboratory; semantic OOD refusal remains unvalidated.
The explicit other-ligand diagnostic confirms why: at 24 h its distance score
has OOD AUROC 0.781 with only 75% in-distribution acceptance and 75% other-ligand
rejection; at 48 h AUROC rises to 0.944 but in-distribution acceptance collapses
to 33.3%. A high OOD rank score therefore does not override the failed refusal
operating point.
In addition, Methods reports 10 ng/mL companion EGF for the TGFB condition while
Supplementary Data 23 records 20 ng/mL. The discrepancy is unresolved rather
than silently normalized. The sweep gate therefore exposes only a blocked
request for KRT5/spatial/cluster observation operators, with no numeric range
or parameter-selection authority.

### Independent-lab HCC827 CyTOF EMT direction bridge

Karacosta et al. (`10.1038/s41467-019-13441-6`) add a Stanford HCC827
TGF-beta/withdrawal CyTOF time course under CC BY 4.0. The exact 79,066,708-byte
Source Data workbook is pinned at SHA-256
`bd6e239008b40cb91b722f16a43c81f43f043adc34f6807a74d094746bfaadfa`.
It contains 90,066 nontransformed single-cell events, 87,354 events in the
eight retained author CCAST clusters, and 29,066 transformed downsampled cells
across E1/E2/E3, pEMT1/2/3, M, and MET. The six clustering features are
E-cadherin/CD324, Vimentin, CD44, CD24, MUC1, and Twist. A streaming OOXML
reader processes these sheets without constructing a large generic workbook
object graph.

The progressive EMT direction is strong descriptively. From untreated 0 d to
10 d TGF-beta, raw median E-cadherin falls from 41.9962 to 1.5490 while raw
median Vimentin rises from 2.1362 to 2,877.2588. After ten days of withdrawal,
E-cadherin rises to 15.3103 and Vimentin falls to 604.6917, reproducing a
partial MET direction. The progressive E-cadherin-down/Vimentin-up pair agrees
with the independent UNSW A549 fixed-IF source across dataset, laboratory,
cell line, and modality. The gate can now reject an EMT-state candidate that
requires the opposite pair of directions.

This is the 25th task head, but not an inferential promotion. Cytobank requires
login and the anonymously downloadable Source Data exposes no biological-
replicate identifier; every raw cell is therefore assigned to one mandatory
experiment split group. The paper reports a similar independent replicate in
Supplementary Figure 6, but that figure is not treated as a released per-cell
replicate. Moreover, the CCAST state labels were constructed using the same six
features, so a classifier trained and tested on those features would be label
circular. Absolute IF-to-CyTOF intensity transfer, uncertainty calibration,
semantic OOD authority, and numeric Aleph ranges all remain refused. The new
capability is strictly a cross-lab direction filter that compresses a future
Aleph sweep.

### Replicated multi-study A549 RNA EMT constraint

The EMT direction bridge now also has five independent GEO laboratories with
biological replication: GSE17708 has 26 samples over an untreated/0.5--72 h
TGFB1 time course, GSE42373 has eight samples across 2D/3D control and combined
TGFB/TNFA treatment, and GSE125369 has two long-TGFB3 and two control RNA-seq
samples. GSE69667 adds a duplicate 0--96 h A549 TGFB1 RNA-seq time course, and
GSE49644 adds triplicate parental/three-week-TGFB endpoints in A549, HCC827,
and NCI-H358. Official matrices, TPM tables, the raw RSEM archive, family SOFT
metadata, and GPL570 annotations are byte- and SHA-256-pinned. GPL570 probe
selection is deterministic (`201131_s_at` CDH1; `201426_s_at` VIM);
GSE49644 uses its sole exact-symbol/Entrez VIM row because
the author-annotated release omits `201426_s_at`.

At the five compatible endpoints, all ten study-by-gene effects have the
expected EMT sign and complete between-group separation. GSE17708 72 h changes
are CDH1 -1.2818 and VIM +0.2645 log2 units; GSE42373 3D changes are -0.9125
and +0.3162; GSE125369 normalized RSEM changes are -24.4623 and +694.5755 on
its own within-study scale. GSE69667 changes are -19.195 TPM and +371.32 TPM;
GSE49644 A549 changes are -5.8120 and +0.7892 log2 units. The HCC827 and
NCI-H358 GSE49644 contrasts also agree in all four marker directions but do not
count as extra independent studies. These directions agree with the independent
A549 IF and HCC827 CyTOF protein observations. Absolute RNA platform scales and
RNA/protein magnitudes are explicitly non-transferable.

This is a stronger external mechanism constraint. Five successful retrospective
study signatures give a nominal exact
one-sided sign-test p=0.03125. However, the two additions were selected from
literature already reporting successful EMT, so study selection was not
outcome-blind and confirmatory direction inference remains blocked alongside
numeric uncertainty promotion.
GSE42373 2D treatment leaves CDH1 slightly increased while the 3D endpoint
decreases it; this culture dependence is retained as a sealed semantic-OOD
challenge. A sealed direction-Hamming diagnostic trained on the five endpoint
signatures refuses this context, but one OOD context cannot validate a general
refusal model. Consequently the layer
can reject or demote opposite-direction EMT candidates before a future Aleph
sweep, while `numeric_range` and `may_emit_sweep_axis` remain null/false.

The replicated panel is now the 26th task-scoped head in the generated
multi-head runtime. It is deliberately excluded from `admissible_head_ids` and
enters `sweep_gate.py` only as a non-authoritative direction constraint. The
gate records all 72 biological samples and five study/lab groups, returns a
blocked calibration request for CDH1/VIM-compatible Aleph observation
operators, and still reports overall `status: refused`. An adversarial contract
test that falsely promotes retrospective direction evidence into numeric uncertainty is rejected
by both the model evaluator and the sweep gate.

To obtain a genuinely outcome-blind confirmation, GSE325309 was frozen before
its processed expression workbook was downloaded or inspected. The committed
protocol fixes all six GSM identities, CDH1/VIM exact-symbol matching, no sample
exclusions, log2(FPKM+1), expected directions, and an exact 20-assignment
composite permutation test. Only after commit `7984e5b` was the 6,266,161-byte
workbook downloaded and pinned at SHA-256
`0d9d636dff0fba25794c41150bbe56be3e57079cdf5bf924548c160faa64d71e`.
Control-to-TGFB1 CDH1 FPKM falls from 17.07/17.61/18.60 to
0.850/0.865/1.017, while VIM rises from 211.41/217.32/206.59 to
410.70/416.34/413.44. The frozen composite is the most extreme of all 20
three-versus-three assignments, giving exact one-sided `p=0.05`; both genes
also have all nine pairwise differences in the expected direction. This is a
pristine external direction holdout, not numeric uncertainty, OOD validation,
or permission to emit an Aleph range.

## PBMC context compression — preregistered external-lab result

This addition is part of the redefined independent-network objective, adopted
because broad Aleph parameter sweeps cannot presently run. The context model
may infer cell identity, rank mechanism families, and recommend the next
discriminating assay on its own evidence. Its optional future Aleph adapter may
not delete a sweep axis, select an Aleph parameter, or emit a numeric range.

The protocol was frozen in commit `b49f9ed` before expression or cell labels
were opened. GSE96583 supplied 24,413 author-labelled singlet PBMC cells. The
released donor identifiers differed from the anticipated list, so the frozen
lexical fallback assigned six whole donors to fitting, donor 1256 to
calibration, and donor 1488 to the sealed internal test. Feature selection used
only fit donors. Of the top 512 non-mitochondrial/non-ribosomal variable genes,
508 aligned by exact HGNC symbol to GSE132044.

The original Broad Single Cell Portal SCP424 annotations were aligned to
22,594 expression rows across two PBMC experiments and multiple technologies.
The seven Zenodo benchmark label CSVs were rejected as classifier targets
because their row counts do not align one-to-one with the official GEO
method-specific expression axes. This prevents a plausible-looking but invalid
cross-source row join.

The preregistered multinomial linear model selected its L2 coefficient only by
six-donor leave-one-donor-out macro-F1. Its results are:

| Endpoint | Result | Frozen threshold | Status |
|---|---:|---:|---|
| sealed donor macro-F1 | 0.9236 | >= 0.70 | pass |
| pristine Broad external macro-F1 | 0.7198 | >= 0.70 | pass |
| pristine Broad external ECE | 0.0462 | <= 0.15 | pass |
| megakaryocyte semantic-OOD AUROC | 0.8567 | >= 0.70 | pass |
| pristine Broad 90% conformal coverage | 0.5025 | >= 0.90 | **fail** |

The failure is scientifically important: class ranking transfers reasonably,
but UCSF donor calibration does not transfer to the Broad laboratory and assay
methods. Consequently the 28th runtime head is
`blocked_external_uncertainty`; `general_cell_type_state_model_ready` remains
false and the sweep gate remains `refused`. The correct next step is a disjoint
cross-lab calibration source, not calibration on the already-opened Broad test
labels.

The IFN-state classifier scores sealed-donor macro-F1 0.9946 with 90% conformal
coverage 0.9268 inside GSE96583. Independently, the frozen 14-gene signature in
GSE72502 increases in every paired donor after IFN-alpha: mean log1p(RPKM)
differences are +3.1203, +3.2819, and +3.9100 for donors A, B, and C. This is a
direction-only cross-lab confirmation; numeric scale transfer is prohibited.

### PBMC disjoint recalibration and final confirmation

The failure above was not repaired using Broad labels. The classifier weights
and 508-feature order were frozen, then only one scalar temperature and six
class-conditional conformal thresholds were fitted on GSE164378 (NYGC/Satija;
151,533 retained cells, eight donors, 490 common features). Classifier parameter
bytes remained identical. The calibration-lab diagnostic reached macro-F1
0.8337, ECE 0.0467, and marginal coverage 0.9000.

Two intended confirmation sources were adversarially refused before prediction:
GSE222647 did not satisfy the frozen donor/condition header aliases, and
GSE158055 exposed `Patients` after the singular `patient` contract had been
frozen. Neither source was relabelled post hoc. The final protocol was committed
before opening the paper-linked Wilk GSE150728 CELLxGENE H5AD, which contains
44,721 cells and exactly six healthy donors. Its raw integer count matrix shared
469 of 508 frozen model genes.

The single-use Stanford confirmation produced:

| Endpoint | Result | Frozen threshold | Status |
|---|---:|---:|---|
| healthy-donor macro-F1 | 0.7458 | >= 0.70 | pass |
| ECE (15 bins) | 0.1117 | <= 0.15 | pass |
| semantic-OOD AUROC | 0.7250 | >= 0.70 | pass |
| marginal conformal coverage | 0.6937 | >= 0.90 | **fail** |
| mean prediction-set size | 1.2135 | <= 3.0 | pass |
| B / DC / monocyte conditional coverage | 0.6113 / 0.5665 / 0.2654 | >= 0.80 | **fail** |

Thus cell-type ranking and coarse OOD detection transfer across labs, but the
claimed uncertainty does not. The asset is permanently consumed, the PBMC head
remains `blocked_external_uncertainty`, and its predictions may not exclude a
candidate mechanism or Aleph sweep region. This is also a concrete result for
the current independent-network objective: classification is promising, while
uncertainty calibration is not yet scientifically deployable.

### Multiprovider IF confirmation-provider refusal

The HPA/OpenCell development encoder and imaging-quality OOD head passed their
development gates, but that does not constitute independent-provider
validation. IDR0006 was preregistered as a single-use nuclear-subset candidate
before opening its pinned annotation CSV. The metadata gate found 16,224 rows
and resolved the frozen ORF, gene, plate, and well aliases without ambiguity.
However, the released phenotype representation is a repeated
`Phenotype 1...8 Term Name` structure, whereas the frozen protocol allowed one
long phenotype column or the declared wide author-label columns. Neither
predeclared schema matched.

IDR0006 was therefore refused without post-hoc alias expansion. The audit
result contains zero parquet requests, zero channel-metadata requests, zero
downloaded image objects, zero accessed pixels, and zero model predictions.
Independent-provider validation, external uncertainty authority, production
eligibility, and Aleph parameter authority all remain false. A replacement
provider requires a new preregistration; the IDR0006 failure cannot be repaired
by changing the consumed gate.

### IDR0168 single-use confirmation — source-integrity refusal

IDR0168/S-BIAD1891 was then frozen as a genuinely independent-provider,
nuclear-localization-subset confirmation. The protocol fixed one lexical-first
field for every gene-by-cell-line group: 45 raw four-channel z-stacks covering
nine genes and five cell lines, including 25 fields whose exact gene was unseen
in the frozen HPA/OpenCell development tensors. It also fixed the native DAPI
and Alexa488 target channels, all localization and corruption-OOD thresholds,
and the rule that any decode mismatch refuses the whole provider without a
replacement field.

All 45 selected archives were acquired before decoding. The receipt contains
34,187,550,441 bytes, 45 object SHA-256 hashes, balanced counts of five fields
per gene and nine per cell line, and records zero decoded pixel archives and
zero model predictions. The receipt was committed separately as `d1edb2f`.

Native decoding then stopped on the eighth field,
`A375_8_RBM23_Z-stack_HF_01_2024-03-07_YunHao_01.29.06.zarr.tgz`. The official
929,941,454-byte object fails gzip CRC
(`0x1a76f459 != 0x6f6f9669`) and contains pyramid metadata for levels 1–3 but
no native `./0/0/.zarray`. Five independent 1 MiB HTTP range checks at the
start, quarter, midpoint, three-quarter point, and end match the acquired file
byte-for-byte, supporting an official-object integrity failure rather than an
unrecorded local substitution.

The frozen no-replacement rule was enforced: seven fields had been decoded in
memory before the failure, but no tensor was written and the model was never
called. The single-use candidate is permanently consumed with status
`failed_source_integrity_and_permanently_consumed`. It contributes no external
IF evidence, uncertainty, production, or Aleph authority. The multitask
contract retains `raw_IF_subcellular_localization_transfer` as blocked and now
records that a different untouched provider is required.

## What remains before the requested capability is achieved

1. Use the newly passing three-lab migration split and the 22,647-row Dryad
   organizational-state tensor as the first training targets for raw-image or
   trajectory encoders. The available GPU host has RTX 4090 capacity, but
   Project policy requires an explicit PI grant and a Slurm `gpu-submit`
   allocation. GPU is not needed for the present feature/runtime audit.
2. Acquire a fourth matched-control migration lab and keep it untouched for
   pristine confirmation of the fixed-origin conformal revision.
3. Replace the single-head calcium curve baseline with a modality encoder that separates
   biological mechanism from lab, cell line, acquisition duration, and stimulus.
   Domain-adversarial or contrastive alignment must be selected using training
   labs only, never the external test lab.
4. Add at least one more independently labelled Piezo1/GsMTx4 dataset so that
   one lab can tune domain alignment and another remains untouched for final
   testing.
5. Expand the one decoded mixed-culture IF pilot into author-labelled,
   cell-segmented image groups. The pilot extracts real pixel statistics but is
   deliberately not treated as a genotype classifier or independent holdout.
6. Acquire biologically replicated PIV/velocity fields with declared units;
   both the completed MDCK Zenodo pair and the hdF Dryad velocity source are
   representation-only and cannot serve as external treatment holdouts.
7. Build observable-to-Aleph parameter likelihoods only after the external task
   gate passes; until then the system may retrieve evidence but must not narrow a
   physics sweep.
8. Finish and checksum the Figshare 28202864 feature archive, train the
   well-disjoint six-class programmed-cell-death head, and use apoptosis only as
   cross-assay latent evidence rather than pretending its MoA labels are
   Annexin-V measurements.
9. Replace the failed scalar PBMC recalibration with a preregistered
   multi-source/domain-conditional uncertainty model. GSE164378 calibration and
   the consumed GSE150728 confirmation prove that one global temperature plus
   class thresholds is insufficient. Reserve a new lab before fitting and a
   second untouched lab before evaluation. Until 90% external coverage and 80%
   supported-class coverage are demonstrated, PBMC context predictions may
   rank but not exclude mechanisms or future sweep regions.
10. Replace the source-integrity-failed IDR0168 confirmation with a new
    preregistered provider. Do not reuse IDR0168 as pristine and do not select a
    different field from its already-opened candidate pool.

## Verification

`test_outer_library.py` currently collects 129 tests for content hashes,
pairing, provenance, authority isolation, sample grouping, time coordinates,
independent lab separation, mechanism contrasts, multitask split contracts,
and refusal behavior. The two new IDR0168 refusal tests, the regenerated summary
hash test, and the multitask/sweep non-authority test pass together (4/4).

The latest whole-file run is not reported as green. It exposed three unrelated
pre-existing numerical replay differences: approximately `3.5e-13` in the HPA
IF conformal radius, plus optimizer-result differences in the BBBC054 microglia
and calcium temporal-CNN reports. Their committed scientific results were not
rewritten. A fourth summary-hash mismatch caused by the intentional multitask
report update was regenerated and passes on retest.
