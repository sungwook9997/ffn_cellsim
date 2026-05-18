# V2 Layer 1 — Sanity-Gate Reference Anchors (DRAFT)

> ⚠️ **PI FRAMING CORRECTION — 2026-05-06 17:55 KST**
> (sister to commit `4ffcbbe`; claude-work `mcp_msg:2710`, codex
> routing `mcp_msg:2706` / `mcp_msg:2708`).
>
> v2 is a **cell-resolved spheroid simulator** (NOT a
> single-cell-isolation framework). 260313 is **v2's primary
> validation anchor at integrated spheroid stages** (NOT a
> reservoir, NOT outside scope). Layer 1 dynamics are
> **component building blocks** of the integrated spheroid
> simulation, validated downstream when V2-3 / V2-5 stages mature,
> **NOT** validated in isolation against single-cell PI data.
>
> This draft's technical anchor content (§1 anchor table, §3
> Magic-Number list, §4 measurement-protocol) is correct as a
> Layer 1 *component-level* contract. The §6 Q1-Q3 design-
> discussion questions remain **open** and are reframed as Layer 1
> component / V2-3 integration choices, not as evidence that v2
> lacks a PI validation anchor.
>
> Older "single-cell anchor" / "reservoir" / "Layer 1 calibration
> source" framing has been corrected throughout this edit pass.

**Status: DRAFT — not locked.** This document expands the Layer 1
sanity-gate anchor table from `docs/v2/10_dev_roadmap_v2.md:63-93` into
inline-citation form so that any future V2-2 Layer 1 component
dynamics module (e.g. active-contour boundary, FA per-state
machine, protrusion event scheduler) can satisfy CLAUDE.md's
Sanity Gate Protocol *before first execution*. This doc is a
**reference table only** — it does **not** select a boundary
representation (vertex / phase-field / hybrid), does **not** fix
any rate constant, and does **not** lock anything. Choices over
algorithm and parameter values belong to the `design-discussion`
room.

Cross-room dependencies (resolved before this draft is promoted to a
locked Sanity Gate spec):

- Active-contour boundary representation choice (vertex vs phase-field
  vs hybrid) — `design-discussion` room, V2-2 unit P1.
- Layer 1 component validation pathway — per the PI framing correction
  (commit `4ffcbbe`), 260313 is v2's primary validation anchor at
  integrated spheroid stages, **not** a single-cell-isolation source.
  Layer 1 component dynamics are validated downstream through
  integrated spheroid simulation when V2-3 / V2-5 stages mature.
  Whether any *single-cell* dataset is acquired separately for
  isolated component testing is a design-discussion choice (see §6
  Q2), not a precondition for v2's primary validation chain.

This file is named `_draft.md` deliberately so the workflow can
rename or relocate it on promotion without retraction overhead.
This draft **intentionally remains untracked at root `docs/`** until
promotion; it is not committed in the post-restructure layout.

**Path-stability note** (updated 2026-05-06 after Option A+B
restructure landing in commit `0533422`; references the original
restructure scope map `mcp_msg:2462` and correctness-review PASS
`mcp_msg:2487`): the Option A+B docs/tests restructure has
**already landed** in main; `docs/v2/*` and `tests/v2/*` are the
canonical post-merge paths. This draft was authored before the
merge with internal references against the pre-merge main-branch
layout (e.g. `docs/v2_*.md`, `docs/10_dev_roadmap_v2.md`,
`tests/test_v2_*.py`); those pre-merge references remain in the
body as historical-context until promotion. When this draft is
promoted, the following sweep is required as part of the
promotion checklist:

- Every `docs/v2_*.md` and `docs/10_dev_roadmap_v2.md` reference in
  §§1, 4, 5, 7 must be re-pointed to `docs/v2/...`.
- This file itself moves to `docs/v2/v2_layer_1_sanity_gate_anchors_draft.md`
  (or its post-promotion locked equivalent).
- `tests/test_v2_*.py` references (none currently in this draft, but
  any added by v1 must be tracked).
- References to `acs/v2/*` files (e.g. `acs/v2/single_cell.py`,
  `acs/v2/imaging_contract/area_extractor.py`, `acs/v2/metrics.py`,
  `acs/v2/measurement_boundary.py`) are **stable** — Option A+B does
  **not** move `acs/v2/` code.

**IF disclaimer** (per Codex `mcp_msg:2492`): all impact-factor
estimates marked `IF ≈ N` in §1 are **provisional** author estimates
based on recent journal-level averages. Promotion to a locked Sanity
Gate spec (§7 criterion (e)) requires explicit cross-check
against current Journal Citation Reports for every cited journal
before any anchor is cited as "IF ≥ 15" in a contract context. Until
that cross-check is recorded, no `IF ≥ 15` claim from this draft is
load-bearing for downstream lock specs.

## §0 Scope and non-scope

Scope (what this doc fixes):
- The set of literature anchors that back each Sanity Gate 1–6 check
  for any future Layer 1 component mechanobiology module (e.g. an
  active-contour boundary solver, FA state machine, or protrusion
  event scheduler considered as a building block of the integrated
  cell-resolved spheroid simulator).
- For each anchor: the exact claim it justifies *and* the explicit
  claim it does **not** justify.
- The Magic-Number risk list (parameters that must declare ranges +
  anchors before any run, per CLAUDE.md Magic-Number Block).
- The measurement-protocol contract for Layer 1 component observables
  (Hard Rule 11 compliance), where component observables are the
  per-module quantities (per-vertex force, per-FA traction,
  per-event protrusion rate) that feed into the integrated
  spheroid-stage validation against 260313 downstream.

Non-scope (out of this doc, must be resolved elsewhere):
- Algorithmic discretisation choices (boundary representation, time
  integrator, vertex spacing, smoothing kernel) — `design-discussion`.
- Numerical CFL / stability bounds — derive from our own force law +
  mobility *after* implementation; biology references **never** import
  numerical bounds.
- Whether any separate single-cell imaging dataset is acquired for
  isolated component testing (orthogonal to v2's primary integrated
  spheroid-stage validation against 260313 — see §6 Q2).
- Cell-cell junction or collective dynamics — Layer 2 / V2-3 territory
  (covered separately in `10_dev_roadmap_v2.md:111-132`).

## §1 Anchor table — expansion of `10_dev_roadmap_v2.md:70-79`

The table below preserves the roadmap row order. For each row: column
(a) re-states the sub-channel, (b) expands citations with author/year/
journal/IF estimate, (c) names the precise claim each anchor backs,
(d) names what each anchor explicitly does **not** justify, and (e)
maps to one or more of the six Sanity Gate checks (`SG1` dimensional,
`SG2` boundary cases, `SG3` conservation, `SG4` numerical sanity,
`SG5` sign/sense, `SG6` measurement-protocol consistency).

### 1.1 Active-matter / time-scale separation

- Anchor A: **Marchetti, Joanny, Ramaswamy, Liverpool, Prost, Rao,
  Aditi Simha 2013, *Rev Mod Phys*** ("Hydrodynamics of soft active
  matter"), IF ≈ 50.
- Anchor B: **Prost, Jülicher, Joanny 2015, *Nat Phys*** ("Active gel
  physics"), IF ≈ 22.
- Backed claim: separation of length/time scales for an active gel
  representation of a single cell — the cortex-cytoplasm composite is
  active-stress-driven, viscoelastic, and dissipative on a polarity
  axis. Justifies SG1 (dimensional analysis) only.
- Explicitly **not** backed: any specific value of `dt_crit`, mobility
  `μ`, viscoelastic relaxation `τ_R`, or contraction-rate `ζΔμ` — these
  are model-and-numerics-specific and must be derived from the actual
  force law after implementation.
- Sanity-Gate map: **SG1** (Re, Ca, De, Pe order-of-magnitude bands
  for the cortex-cytoplasm composite); **SG3** (active-stress channel
  is *not* conservative — SG3 must record this dissipation channel
  explicitly).

### 1.2 Stochastic edge / shape coupling

- Anchor A: **Keren, Pincus, Allen, Barnhart, Marriott, Mogilner,
  Theriot 2008, *Nature*** ("Mechanism of shape determination in
  motile cells"), IF ≈ 50.
- Anchor B: **Raynaud, Ambühl, Gabella, Bornert, Sbalzarini, Meister,
  Verkhovsky 2016, *Nat Phys*** ("Minimal model for spontaneous cell
  polarization and edge activity in oscillating, rotating and migrating
  cells"), IF ≈ 22.
- Backed claim: single-cell shape and migration timescale couple
  through stochastic edge switching; polarity persistence is a
  measurable autocorrelation / displacement-cosine quantity, not a
  fitted relaxation constant.
- Explicitly **not** backed: a specific value of `τ_p` (polarity
  memory) — Sanity Gate 6 requires reporting `C(Δt)` and
  displacement-polarity cosine over a **pre-registered** window before
  any optimisation.
- Sanity-Gate map: **SG6** (polarity-vector autocorrelation `C(Δt)`
  and displacement-polarity cosine over pre-registered window — the
  measurement protocol for polarity *itself*); **SG5** (edge activity
  fluctuations are sign-symmetric on average, polarity bias breaks the
  symmetry monotonically — sign check).

### 1.3 Cortex / surface tension

- Anchor A (single-cell): **Chugh, Clark, Smith, Cassani, Dierkes,
  Ragab, Roux, Charras, Salbreux, Paluch 2017, *Nat Cell Biol***
  ("Actin cortex architecture regulates cell surface tension"),
  IF ≈ 21.
- Anchor B (cell-cell, **moves to Layer 2**): Maître, Berthoumieux,
  Krens, Salbreux, Jülicher, Paluch, Heisenberg 2012, *Science*
  ("Adhesion functions in cell sorting by mechanically coupling the
  cortices of adhering cells"), IF ≈ 56. Listed for cross-reference
  only; this paper anchors interfacial tension at cell-cell contacts,
  not single-cell free-boundary cortex.
- Backed claim (Anchor A only): the actin cortex contributes a
  shape-restoring contractile surface tension at the free-cell
  boundary; the magnitude of the tension is set by myosin activity and
  cortex mesh architecture, both measurable.
- Explicitly **not** backed: any Col1-spreading-induced cortex
  modulation. No clean IF≥15 anchor exists for the latter — if a
  Col1-induced cortex term is introduced, it must be flagged
  exploratory in the parameter-role taxonomy.
- Sanity-Gate map: **SG5** (cortex tension sign — contractile,
  shape-restoring, opposes any expansion of the boundary); **SG3**
  (cortex performs negative work on isotropic expansion ⇒ acts as a
  restoring channel only, never an energy source).

### 1.4 Protrusion / lamellipodia events

- Anchor A: **Machacek, Hodgson, Welch, Elliott, Pertz, Nalbant, Abell,
  Johnson, Hahn, Danuser 2009, *Nature*** ("Coordination of Rho GTPase
  activities during cell protrusion"), IF ≈ 50.
- Anchor B: **Tkachenko, Sabouri-Ghomi, Pertz, Kim, Gutierrez, Machacek,
  Ginsberg, Danuser, Hahn 2011, *Nat Cell Biol*** ("Protein kinase A
  governs a RhoA-RhoGDI protrusion-retraction pacemaker in migrating
  cells"), IF ≈ 21.
- Backed claim: lamellipodia protrusion is event-like with a Rho/Rac/
  Cdc42 + PKA signalling cycle; protrusion edge-events are
  inter-arrival-distributed in a way *consistent* with an
  inhomogeneous Poisson hazard, but the Poisson form is a
  **simulation ansatz**, not a measured law.
- Explicitly **not** backed: any single value of `λ0` (baseline event
  rate). Validation is via inter-arrival distribution / coefficient
  of variation and spatial pair-correlation `g(r)` first peak — never
  by fitting `λ0` to PI's `A(t)`.
- Sanity-Gate map: **SG6** (the *measurement* of `λ` requires live
  boundary-tracking imaging at sub-minute cadence; if the available
  channel cannot resolve individual events, the gate must compare
  inter-arrival CV / `g(r)` and *not* a raw `λ`); **SG2** (boundary
  case `λ → 0` ⇒ static cell is well-defined; `λ → ∞` ⇒ saturate to
  continuous protrusion, must not NaN).

### 1.5 Focal adhesion state machine

This is a 5-anchor stack. Each anchor backs a **distinct** part of the
state machine; no single paper backs all rates.

- Anchor A — **primary measurement paper** (rigidity sensing / force
  fluctuations): **Plotnikov, Pasapera, Sabass, Waterman 2012,
  *Cell*** ("Force fluctuations within focal adhesions mediate
  ECM-rigidity sensing to guide directed cell migration"),
  IF ≈ 64. Provides traction-force-microscopy (TFM) measurements of
  intra-FA force fluctuations on substrates of graded rigidity.
- Anchor B — **primary measurement paper** (assembly / maturation
  lifetimes): **Choi, Vicente-Manzanares, Zareno, Whitmore, Mogilner,
  Horwitz 2008, *Nat Cell Biol*** ("Actin and α-actinin orchestrate
  the assembly and maturation of nascent adhesions in a myosin II
  motor-independent manner"), IF ≈ 21. Live-imaging measurement of
  nascent → mature transitions in fibroblast leading edges.
- Anchor C — **primary measurement paper** (state identity /
  molecular architecture): **Kanchanawong, Shtengel, Pasapera, Ramko,
  Davidson, Hess, Waterman 2010, *Nature*** ("Nanoscale architecture
  of integrin-based cell adhesions"), IF ≈ 50. iPALM
  super-resolution measurement establishing the layered nm-scale FA
  architecture (integrin / talin / actin strata).
- Anchor D — **review** (clutch organisation): **Case, Waterman 2015,
  *Nat Cell Biol*** ("Integration of actin dynamics and cell adhesion
  by a three-dimensional, mechanosensitive molecular clutch"),
  IF ≈ 21. *Review article* — cites primary clutch measurement
  papers; rate constants must be sourced from those primaries (e.g.
  Hu, Plotnikov, Wang papers cited within), not from this review.
- Anchor E — **primary measurement paper** (talin / clutch rigidity
  threshold): **Elosegui-Artola, Oria, Chen, Kosmalska,
  Pérez-González, Castro, Zhu, Trepat, Roca-Cusachs 2016, *Nat Cell
  Biol*** ("Mechanical regulation of a molecular clutch defines force
  transmission and transduction in response to matrix rigidity"),
  IF ≈ 21. Talin-unfolding rigidity-threshold measurement on PDMS
  gradient substrates.
- Backed claim (composite): five-state machine identity
  `unbound → nascent → mature → slipping → released` with
  rigidity-dependent transition rates, force-fluctuation-mediated
  rigidity sensing, and a talin-unfolding-based slip threshold.
- Explicitly **not** backed: any single condition-and-cell-type
  -independent rate constant. All transition rates must be expressed
  as ranges with cell-type and rigidity dependence flagged. Plotnikov
  alone does **not** give all five rates and must not be cited as
  such.
- Sanity-Gate map: **SG3** (FA state-machine population is conserved
  per cell; transitions move probability mass between states without
  loss); **SG5** (slipping-state force *opposes* protrusion direction
  ⇒ negative work; mature-state stiffens substrate coupling
  monotonically); **SG6** (FA gates require TFM / paxillin marker
  channels; without them, only an *aggregate* observable like total
  adhered area or aggregate clutch engagement can be reported, never
  per-state lifetimes).

### 1.6 Filopodial traction (introduce only if needed)

- Anchor: **Chan, Odde 2008, *Science*** ("Traction dynamics of
  filopodia on compliant substrates"), IF ≈ 47.
- Backed claim: filopodia generate substrate-stiffness-dependent
  traction in a stick-slip-like dynamics, distinct from lamellipodial
  traction.
- Explicitly **not** backed: a need for a filopodial channel in the
  base Layer 1 model. Add only if the imaging data shows
  filopodia-dominated edge events that the lamellipodial channel
  cannot reproduce.
- Sanity-Gate map: **SG5** (filopodial traction sign — stick-slip
  cycle with adhesion phase pulling toward substrate, slip phase
  releasing); **SG2** (boundary case: zero filopodia ⇒ degenerate to
  Anchor 1.4 lamellipodial channel cleanly).

### 1.7 Vertex / active-contour numerics

- Anchor: **none** — biology references do **not** justify
  discretisation choices.
- Note: Bi & Park 2014, *Phys Rev X*; Park, Atia, Mitchel, Fredberg,
  Butler 2015, *Nat Mater* are valid Layer 2 *tissue* references for
  collective jamming/shape-index gates — they are **insufficient** for
  Layer 1 single-cell active-contour numerics. Citing them as a
  numerical justification is a category error.
- Backed claim: nothing — discretisation is justified by our own
  convergence / symmetry / area-drift tests, not by literature.
- Sanity-Gate map: **SG1** (`dt_crit` derived from our force law and
  mobility, *not* imported); **SG4** (grid resolution `dx` ≪ smallest
  physical feature; vertex spacing convergence test); **SG3** (area
  drift / mass conservation under symmetry — must be measured per
  scheme, not assumed).

### 1.8 Input segmentation

- Anchor A — **primary segmentation method paper**: **Stringer, Wang,
  Michaelos, Pachitariu 2021, *Nature Methods*** ("Cellpose: a
  generalist algorithm for cellular segmentation"), IF ≈ 36
  (provisional). Cellpose provides a generalist deep-learning cell-
  segmentation algorithm trained on a heterogeneous cytoplasmic-
  fluorescence corpus. This anchor justifies *how* a per-frame
  boundary mask may be produced from a fluorescence/transmitted-light
  image, in single-cell or multi-cell input.
- Anchor B — **computer-vision foundation-model method paper, not a
  biology/layer anchor**: **Kirillov, Mintun, Ravi, Mao, Rolland,
  Gustafson, Xiao, Whitehead, Berg, Lo, Dollár, Girshick 2023**
  ("Segment Anything", arXiv:2304.02643 / Meta AI Research). SAM is a
  *computer-vision* foundation-model method paper without an IF≥15
  *biology* anchor; it is admissible only as **optional segmentation
  method provenance** (e.g. as one of several segmentation backends
  whose mask output feeds the §1.8 input pipeline). SAM does **not**
  back any biology, mechanobiology, or Layer 1 / Layer 2 claim, and
  must never be cited as a biology anchor in a Sanity Gate.
- Backed claim: input segmentation only — the boundary mask used as
  initial condition or per-frame observation. The two anchors are
  interchangeable for this purpose; the choice is a tooling decision,
  not a science decision.
- Explicitly **not** backed: solver vertex spacing, `dx`, smoothing
  bandwidth, or any internal numerical resolution choice. The mask is
  evidence about a single observed boundary, not a numerical scheme.
  Cellpose / SAM also do **not** justify any *measurement* of
  protrusion events, traction, FA dynamics, or any Layer 1+ biology
  observable — they only justify the mask.
- Sanity-Gate map: **SG6** (the *input* of any Layer 1 measurement is
  the segmentation mask — this anchor only justifies how that mask is
  produced, not the downstream metric).

## §2 Sanity Gate ↔ anchor coverage matrix

| Sanity Gate check | Anchors that back it | Anchors that explicitly do **not** back it |
|---|---|---|
| SG1 dimensional | 1.1 (active-matter scales), 1.7 note (`dt_crit` derived locally) | 1.4 (Poisson `λ0`), 1.5 (single rate constants) |
| SG2 boundary cases | 1.4 (`λ → 0`/`∞`), 1.6 (zero filopodia degenerate) | — |
| SG3 conservation | 1.1 (active-stress dissipation channel declared), 1.3 (cortex shape-restoring only), 1.5 (FA state population), 1.7 (own area-drift tests) | — |
| SG4 numerical sanity | 1.7 (own convergence / symmetry tests) | 1.7 (biology refs explicitly do **not** validate `dx` / `dt`) |
| SG5 sign/sense | 1.1 (active stress on polarity axis), 1.2 (polarity symmetry-breaking), 1.3 (cortex contractile), 1.5 (slip vs mature signs), 1.6 (filopodial stick-slip sign) | — |
| SG6 measurement protocol | 1.2 (polarity autocorrelation), 1.4 (inter-arrival CV / `g(r)` not raw `λ`), 1.5 (TFM / paxillin requirement), 1.8 (mask provenance) | 1.4 (raw `λ0` fitting forbidden) |

If a Layer 1 module's Sanity Gate analysis cites an anchor outside its
column-2 row above, that is a contract violation and must be flagged
in code review.

## §3 Magic-number risk list

(Verbatim from `10_dev_roadmap_v2.md:81-86`, repeated here so any
future module's Sanity Gate doc can cross-reference a single canonical
list.)

Parameters that must declare a numeric **range** plus a literature
anchor *before* the run, per CLAUDE.md Magic-Number Block:

- `λ0` — baseline protrusion event rate. Anchor 1.4 (Machacek /
  Tkachenko). Must declare range from inter-arrival distribution, not
  a single value.
- `τ_ref` — refractory period after a protrusion event. Anchor 1.4.
- `ξ_p` — protrusion angular width / correlation length. Anchor 1.4
  (`g(r)` first peak).
- `γ_cortex` — cortex tension. Anchor 1.3 (Chugh).
- area / perimeter stiffnesses — no clean biology anchor; declare
  exploratory and constrain by SG2 boundary cases (must keep cell
  area positive and finite).
- `k_on` / `k_off` for FA states — Anchor 1.5 (Choi / Plotnikov /
  Elosegui-Artola). Express as cell-type + rigidity-dependent ranges.
- FA slip threshold — Anchor 1.5 (Elosegui-Artola, talin unfolding).
- `τ_p` — polarity memory. Anchor 1.2 (Keren / Raynaud). Must be a
  *measured* autocorrelation, not a fitted constant.
- traction-coupling gain — **a derived model parameter**, not a
  directly-measured constant. It is *constrained* by Anchor 1.5
  primary measurements (Plotnikov 2012 for intra-FA force-fluctuation
  magnitudes; Elosegui-Artola 2016 for rigidity-threshold response),
  but the numeric gain value depends on the **chosen clutch / force-
  law form** and cannot be imported directly from either paper.
  Promotion of any specific gain value to a Magic-Number range
  therefore requires (i) the chosen clutch/force-law expression
  recorded explicitly, (ii) the constraints from Plotnikov 2012 +
  Elosegui-Artola 2016 propagated through that expression, and
  (iii) a sensitivity check showing the gain range is not tuned to a
  PI-data target. **Not** sourced from the Case & Waterman 2015
  review — that paper is a synthesis and any numeric value cited
  through it would be a Magic-Number loophole.

**Rule (CLAUDE.md Hard Rules)**: never tune any of the above to PI
A(t). Any parameter whose value is selected to make a specific
simulated number match a target *must* be flagged in the
Magic-Number Block check and surfaced to the PI.

## §4 Measurement-protocol notes (Hard Rule 11 compliance)

(Synthesised from `10_dev_roadmap_v2.md:88-93` plus the imaging
measurement-protocol gate seal `49e246c`.)

- Top-down projected polygon area — `SingleCellState.projected_area_um2()`
  shoelace, in `acs/v2/single_cell.py` — is the correct simulation
  analog for top-down segmentation observations. This is the Layer 1
  single-cell analog of the spheroid-stage measurement that
  `acs/v2/imaging_contract/area_extractor.py` reproduces from PI's
  saved-mask PNGs.
- Substrate-contact area is **not** a valid analog of a top-down
  projection observation. The April 2026 v13 episode (Stage 1a)
  recorded in CLAUDE.md Hard Rule 11 is the canonical anti-pattern;
  any new Layer 1 metric must explicitly state which experimental
  imaging modality it matches.
- Protrusion / FA gates require live boundary tracking + TFM /
  paxillin / fluorescent marker channels. If those channels are not
  available in the chosen single-cell dataset, the gate must compare
  *distributional* observables (inter-arrival CV, spatial `g(r)`),
  never a raw `λ`.

### 4.1 Existing measurement registry (`acs/v2/metrics.py`)

(Per Codex `mcp_msg:2492` review note 1 — listed for measurement
modality discipline only.)

The repository already exposes a **pure measurement registry** at
`acs/v2/metrics.py` (added in `75824af`, expanded in `6e48f18`),
keyed on `(metric_name, ArtifactKind)`. This file is **not** a
Layer 1 physics anchor and does **not** justify any of the rate
constants, sign senses, or stability bounds discussed in §1; it
implements the *measurement modality* side of Hard Rule 11 only.

Concretely, the registry's relevance to Sanity Gate 6 is that it
fixes the channel between an `ArtifactKind` (e.g. `BOUNDARY_CONTOURS`,
`SEGMENTATION_MASK`) and a unit-stamped scalar metric (`um2`, `um`,
…) so that any future Layer 1 module:

- declares **which artifact** it measures from (preventing the
  substrate-contact-area-vs-top-down-projection class of mismatch);
- reports in the units that match `MetricSpec.unit` strings in
  `acs/v2/data_contract.py`;
- defensively calls `MeasurementBoundary.validate()` at the
  measurement boundary (preventing direct-construction shortcuts that
  bypass geometry checks).

For Layer 1 work, this means:

- The top-down projected polygon area metric in §4 above maps onto
  the registry's `projected_area` family (boundary-keyed) and
  matches the top-down imaging modality used by the
  `imaging_measurement_protocol_gate`.
- Any *new* Layer 1 morphology metric (perimeter, roughness, etc.)
  added in support of single-cell active-contour validation must be
  registered through this registry, with its `ArtifactKind` declared.
  The registry is the **measurement-modality contract**, not a
  physics decision.
- The registry itself imposes no parameter range, no rate constant,
  and no boundary representation choice; those decisions belong to
  §1, §3, and the `design-discussion` room respectively.

## §5 Known limitations and explicit exclusions

- v2's primary validation anchor is 260313 at the integrated
  spheroid stage (PI framing correction `4ffcbbe`); Layer 1
  *component* parameters are derived from literature anchors at
  the component level and validated downstream through integrated
  spheroid simulation, **not** by direct fit to single-cell PI
  data. Whether any separate single-cell dataset is acquired for
  isolated component testing is a design-discussion choice (see §6
  Q2), independent of v2's primary integrated-validation chain.
  No fitted parameter on any 260313 row at the Layer 1 component
  layer (Hard Rule 1 + integrated-validation-only).
- The full V2-2 anchor stack assumes a *boundary-resolved* single-cell
  representation. If the eventual algorithm choice is phase-field
  rather than vertex/active-contour, anchors 1.5 (FA per-state) and
  1.4 (per-event protrusion) require a re-mapping to a continuum
  density observable — flagged as an open `design-discussion` item.
- Reviews vs primaries (per Codex `mcp_msg:2492` review note 2 +
  v1 §1.5 annotations): Anchor 1.1 (Marchetti 2013 *Rev Mod Phys*,
  Prost et al. 2015 *Nat Phys*) are review/perspective articles —
  they justify SG1 dimensional analysis only and are **not** sources
  for any rate constant. Anchor 1.5 Anchor D (Case & Waterman 2015
  *Nat Cell Biol*) is a review and is now explicitly *not* the
  source for the traction-coupling-gain range in §3 (corrected v1).
  All Magic-Number ranges in §3 must trail back to primary
  measurement papers — review citations alone are insufficient.
- IF estimates in §1 are **provisional** author estimates from recent
  journal-level averages; they have not been cross-checked against
  the current Journal Citation Reports in this draft. Final lock
  spec (§7 criterion (e)) must verify each IF against JCR
  before any anchor is cited as "IF ≥ 15" in a contract context.
  Until that verification is recorded, no `IF ≥ 15` claim from this
  draft is load-bearing for downstream lock specs.
- The roadmap-anchored measurement-protocol claim in §4 is operationally
  channeled through `acs/v2/metrics.py` (per §4.1), but that registry
  itself does **not** justify any physics or numerics — it is a
  measurement-modality contract only. Treating the registry as a
  Layer 1 anchor would be a category error.

## §6 Open design-discussion questions (out of scope here)

Do **not** answer in this draft. Per the PI framing correction
(commit `4ffcbbe`), all three questions are reframed as Layer 1
component / V2-3 integration choices, **not** as evidence that v2
lacks a PI validation anchor — v2's primary validation is at the
integrated spheroid stage against 260313, downstream of these
component-level decisions. Surfacing for the `design-discussion`
room:

1. Boundary representation — vertex, phase-field, or hybrid? This
   selection cascades into which anchors are *operationally*
   applicable to which Layer 1 component observables (see §5), and
   into how component dynamics integrate into the V2-3 cell-resolved
   spheroid simulation.
2. Optional single-cell dataset for isolated component testing —
   whether to acquire any separate single-cell dataset to test
   Layer 1 components in isolation, *in addition to* v2's primary
   integrated-validation against 260313. This is an optional
   isolation-test channel, not a precondition for the primary
   validation chain.
3. Whether to require **TFM** for any Layer 1 component dynamics
   that includes the FA state machine, or to drop FA from the
   *first* runnable Layer 1 component build and add it as a Layer 1+
   extension. Independent of v2's spheroid-stage validation against
   260313 (TFM availability gates component-level isolation tests,
   not v2's primary integrated validation).

## §7 Provenance and intended evolution

This draft was written 2026-05-06 KST in the `implementation-work`
room while Codex held the active claim on
`v2-repo-structure-cleanup` (Option A+B). It is intentionally
**not** locked. Promotion to `docs/v2/v2_layer_1_sanity_gate_anchors_locked.md`
(or whatever the post-restructure naming convention dictates)
requires:

- (a) `design-discussion` resolution of §6 questions
- (b) Codex review checkpoint with explicit `id=...` ack on the
  anchor table contract (per memory `reviewed_by_only_after_explicit_diff_ack`)
- (c) PI sign-off on the Magic-Number range declarations in §3
- (d) Verification that 260313 is **not** referenced as a direct
  Layer 1 component parameter-fitting source (per PI framing
  correction `4ffcbbe`: 260313 is v2's primary validation anchor
  at integrated spheroid stages, **not** a single-cell-isolation
  parameter-fitting source); Hard Rule 1 + this draft's §5 first
  bullet.
- (e) JCR cross-check of every `IF ≈ N` value in §1 against current
  Journal Citation Reports, with each (journal, IF, JCR-edition)
  triple recorded inline. Until (e) is satisfied, no anchor in §1
  may be cited as "IF ≥ 15" in a downstream lock spec.
- (f) Path-stability sweep across §§1, 4, 5, 7. The Option A+B
  restructure has already landed in commit `0533422`; remaining
  pre-merge references in the body must be re-pointed to the
  post-merge layout: `docs/v2_*.md` → `docs/v2/...`,
  `docs/10_dev_roadmap_v2.md` → `docs/v2/10_dev_roadmap_v2.md`,
  `tests/test_v2_*.py` → `tests/v2/...`. This draft file itself
  must be relocated from root `docs/v2_layer_1_sanity_gate_anchors_draft.md`
  to `docs/v2/v2_layer_1_sanity_gate_anchors_draft.md` (or its
  post-promotion locked equivalent).

Until (a)–(f) are met, treat this file as a working reference, not as
a contract.

### 7.1 v0 → v1 changelog (this round)

- §0: expanded path-stability note into a sweep checklist; added IF
  disclaimer block (per Codex `mcp_msg:2492` review note 2).
- §1.5: annotated each anchor in the FA 5-stack as primary
  measurement vs review; flagged Case & Waterman 2015 explicitly as
  a review whose rate constants must trail back to its primaries.
- §1.8: expanded both segmentation anchors (per Codex
  `mcp_msg:2497` required edit 1) — Cellpose given full citation
  coordinates and provisional IF; SAM explicitly classified as a
  computer-vision foundation-model method with no biology IF≥15
  anchor, admissible only as optional segmentation provenance and
  never as a biology layer anchor.
- §3: removed the bare "Case & Waterman" citation for
  traction-coupling gain; re-anchored to Plotnikov 2012 + Elosegui-
  Artola 2016 primaries (v1 first pass), then sharpened in v1
  second pass (per Codex `mcp_msg:2497` required edit 2) into
  "derived model parameter, constrained by primaries through a
  chosen clutch/force-law form, with explicit sensitivity-check
  requirement to avoid a Magic-Number loophole."
- §4: added §4.1 — measurement-registry implementation channel
  (`acs/v2/metrics.py`), framed as Hard Rule 11 modality discipline
  only, **not** a Layer 1 physics anchor (per Codex `mcp_msg:2492`
  review note 1).
- §5: added Codex `mcp_msg:2492` cross-references; explicit
  reviews-vs-primaries paragraph; promoted IF disclaimer to a
  load-bearing limitation.
- §7: added promotion criteria (e) JCR cross-check and (f)
  path-stability sweep. Also fixed two stale "§7 criterion d,
  expanded" references in §0 and §5 to point to current criterion
  (e) (per Codex `mcp_msg:2497` required edit 3).

The v1 changelog is preserved here for the duration of the
working-reference state. When promotion to a locked spec happens, the
changelog should either be moved to a separate revision-history page
or pruned to a single line that records the promotion event.
