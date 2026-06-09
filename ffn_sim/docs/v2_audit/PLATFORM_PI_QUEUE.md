# Platform PI Queue — 2026-06-09

## 🔀 SESSION SPLIT (2026-06-09, PI) — SF generation-limit, two follow-on sessions
The passive 2-layer basal contractile apparatus BUILD is COMPLETE (B1-B4 + S1 + bending,
all build-time gates PASS; commits `6f67b99`→`644f770` on `h7/compartment-platform`). The
ONLY remaining blocker to the live Kumar 10-30 nN gate is the **generation-limit** (the SF
instance of the cortical γ-floor: force-budget ~2000× under band at lit density). The CFL
wall is SOLVED (bending wired → bend-before-stretch ratio ~1.5e5 → AFINES soft-stretch →
explicit BAOAB feasible, no frozen integrator). PI split the generation-limit decision into
two follow-on sessions — **run each in its OWN git worktree** (1-session-1-worktree rule;
both touch basal_mesh.py / the force budget). (i) is decisive; run it first or in parallel.

### Session (i) — "CLOSE THE FLOOR": derive Route-B NMII force-scaling from a density datum
- **Mandate**: find the native NMII content/density of a SINGLE ventral stress fiber
  (minifilaments per cross-section / per µm; + per-minifilament stall force) — via
  deep-research (cross-verified, like the N_filaments/EA audit), since the data is likely
  old/sparse. DERIVE the Route-B mesoscale force-scale factor (`myosin.py`
  `mesoscale_force_scaling`: factor = native_NMII/effective_NMII, grid-invariant → NOT a
  magic number). Apply to `sf_myosin_`; test if the SF tension reaches Kumar 10-30 nN
  (force-budget `scripts/h7_basal_sf_force_budget.py` first; then optionally the soft-stretch
  explicit-BAOAB dynamic run on gbook). This ALSO fills the cortical γ-floor's missing density
  datum from the SF side.
- **Boot**: BASAL_MESH_DESIGN_2026-06-09.md, PLATFORM_AUTORUN_LOG_2026-06-09.md (Phase-3 table),
  the γ-floor docs (H7_ACTIVE_GAMMA_SYNTHESIS_2026-06-09 in the main repo + the
  `gamma-floor-layered-resolution` memory), `h7_basal_sf_force_budget.py`, `myosin.py`
  `mesoscale_force_scaling`, `basal_bending_report` (bend-before-stretch headroom).
- **Deliverable**: factor DERIVED + Kumar verdict (PASS/REFUTE with the density datum stated),
  or HALT to PI if no usable native-NMII-per-SF datum exists.

### Session (ii) — "ACCEPT THE FLOOR": generation-bound conclusion + γ-floor integration
- **Mandate**: document the SF tension as GENERATION-limited (force-budget ~2000× under Kumar
  at lit density), UNIFY with the cortical γ-floor conclusion (one generation-limit story —
  same ½·n·f·ℓ budget, same MCF7-density-datum gap), report honestly with the bounds. Consider
  the gate-reframe (active-fraction / band scope). Produce the SF-generation-limit synthesis
  doc; cross-link to the cortical γ-floor synthesis. This is the honest fallback if (i) finds
  no density datum that closes Kumar.
- **Boot**: same docs as (i) + the cortical γ-floor synthesis.
- **Deliverable**: SF generation-limit synthesis doc, integrated with the cortical γ-floor.

### Reconciliation
(i)'s outcome largely determines (ii): if a real native-NMII density closes Kumar via Route B,
the SF line is unblocked (and the cortical γ-floor gets its density datum); if not, (ii)'s
generation-bound conclusion stands (SF + cortical γ as one bounded generation-limit). Whichever
session lands the result updates this queue + the Dev-Logs board.

### ✅ Session (ii) LANDED 2026-06-09 (branch `h7/sf-generation-floor`) — generation-bound ACCEPTED
**Deliverable: `docs/v2_audit/SF_GENERATION_LIMIT_SYNTHESIS_2026-06-09.md`** (+
`scripts/h7_basal_sf_force_budget.py::decompose_generation_gap` + 5 consistency tests).
- **SF tension is GENERATION-limited, N_filaments-INDEPENDENT.** Raw budget ≈5 pN, ~2014×
  under Kumar. **Decomposed (labeled):** ~11× per-minifilament fidelity (brief 5 pN → lit 56 pN
  dipole = the SAME §9 cortical per-head/stiffness fix) × ~180×(floor)/~360×(centre)
  **cross-sectional NMII count per SF** = the Route-B native:effective factor = the MISSING
  MCF7 SF-NMII density datum (= session (i)'s target). Engagement realism ~9× WORSE.
- **UNIFIED with the cortical active-γ floor:** same ½·n·f·ℓ aggregate-motor budget, same f≈56 pN
  dipole, same missing-MCF7-density-datum gap (areal n₂D for cortex, cross-sectional N for SF),
  same architecture caveat — raw count is the WRONG lever; tension is set by ORGANISATION
  (cortex: actin overlap, Chugh/Truong Quang; SF: sarcomeric polarity, §11 SF-2c −61 pN
  slackening, Hotulainen-Lappalainen).
- **Bound honesty:** budget is an UPPER bound on a COHERENT bundle; the current random-polarity
  construction realises ~0 net traction even with the per-head fix → the **sarcomeric-polarity
  SF construction (§11 NEXT)** is the distinct organisation lever (separate from (i)'s count datum).
- **PI items surfaced (no gate-loosening):** (1) active-fraction scope (Kumar total vs
  blebb-sensitive; SF-2c shows ~+3476 pN passive pre-tension), (2) tool-role scope (fine-grained
  tool SUPPLIES ζΔμ ∝ ρ_NMII vs owns the absolute band), (3) density-conditioned provisional band,
  (4) cross-line per-head reconciliation (0.5 / 2 / 8.48 pN across SF-brief / cortical-dipole / §9).
- **Reconcile w/ (i):** if (i) lands a defensible MCF7 SF-NMII datum giving N_cross≈180–360, Route B
  closes Kumar AND fills the cortical density gap → this conclusion updates toward (i). Else it stands.

---


## ✅ ACTIVATED under the PI ownership grant (2026-06-09 "소유권 허용 / 모든 것 달려")
- **osmotic_regulation → LIVE** (gate PASS: τ_RVD=3.2s in band, RVD sign, no-contamination,
  off-identity). Post-build attach in manifest.py (no cell.py surgery).
- **microtubules → LIVE** (gate PASS: aster assembled, CFL passed at cytoplasm drag,
  γ_soft IDENTICAL OFF/ON = no contamination, L_p in band). Full snapshot-extension
  wiring in cell.py + manifest.py. **n_mt CAPPED at 7** (see new blocker below).
- **intermediate_filaments → LIVE** (gate PASS: cage assembled +480 part, γ_soft
  IDENTICAL OFF/ON, if_ denylisted). LINEAR small-strain path; the nonlinear
  strain-stiffening Table law stays PI-pending (NotImplementedError) → SECONDARY
  stiffening gate deferred. Fixed 2 crash-on-enable bugs (gsd-None velocity/angles
  guards in the IF extender — masked by the smoke, caught by the real build).
- **linc → LIVE** (gate PASS: bridges formed n=40, FORCE-FREE per-bond EXACT-r0,
  γ_soft IDENTICAL OFF/ON = no contamination, off-identity +0 particles, CFL ~4-order
  headroom). Option A: bonds nucleus_bead → nearest perinuclear if_bead (the LIVE IF
  cage); requires=('nucleus','intermediate_filaments'). Per-bond EXACT-r0 (FA-clutch
  convention) so each bridge is force-free at its as-built separation. Added
  `unique_acceptor` 1:1 matching (one nesprin per IF anchor) to bound per-if_bead bond
  degree — WITHOUT it the nlist per-particle exclusion cap overflows (same class as MTOC).
  **PI items to ratify (LINC is LIVE but these are pending):**
  - **k_linc = 1e-2 N/m** (route-A folded-rod pre-unfolding secant, Rief 1999; ORDER, conf
    MEDIUM). Module default stays None (un-configured enabled build raises); the config
    carries the candidate. PI to ratify route A vs B (WLC entropic ~6e-6) vs C (md.bond.Table
    nonlinear). Also the 2→8 pN f_rest re-anchor + Déjardin-vs-Arsenovic attribution + the 6
    SourceEvidence rows (Rief99, Déjardin20, Arsenovic16, Autore13, Crisp06, Lombardi11) before
    any deliverable cite.
  - ⚠ **Mesoscale-r0 fidelity caveat (FOLLOW-UP).** At the ×40 mesoscale the nucleus↔if_bead
    bead separation is ~0.5-1.5 µm (NOT the 50 nm real nesprin span), so r0 is an EFFECTIVE
    mesoscale coupling length and the 8 pN f_rest oracle is approximate. **Option C** (seed a
    dedicated R_nuc+50 nm perinuclear acceptor cap — a new particle layer) is the faithful-span
    alternative; surfaced per the "fidelity-impact → PI queue" rule. The gate's no-contamination
    + force-free + bridges-formed controls hold regardless; the [2,10] pN <T_linc> resting-tension
    oracle is the deferred physics gate (needs an equilibrated actomyosin run; raw full-cell run
    trips the BAOAB guard — equilibration prelude required, same as the other compartments).
  - ⚠ **GPU many-bond-types perf (FOLLOW-UP).** Per-bond EXACT-r0 mints one bond type per bridge
    (~n_IF_filaments); combined with IF's 16 crosslink bins HOOMD warns "many bond types perform
    poorly / shared-memory errors on the GPU." Fine for the CPU build-time gate; flag for the GPU
    production path (a coarser shared per-r0 bin needs ~hundreds of bins at k_linc=1e-2 to stay
    thermal, so exact-r0 is actually the leaner choice — a native bond plugin or a stiffness-aware
    binning is the real fix).
- **membrane_reservoir → LIVE** (gate PASS: mesh assembled n_tethers≈1969, cross-layer
  no self-pairs, FORCE-FREE at the offset, γ_soft IDENTICAL OFF/ON = no contamination,
  off-identity, bleb paths honestly blocked). STATIC mem_tether mesh on an OWN
  radially-offset mem_node layer (BLOCKER-1 fix); requires=('membrane_surface',).
  Degree-aware unique-acceptor selection (max_anchor_degree=5) bounds per-cortex bond
  degree to +1 (nlist exclusion-cap safe). Anchored defaults: W_MCA=1e-5 J/m² (KU-3.B1.4),
  k_tether=0.1 N/m (ERM order), max_tether_dist=200 nm (KU-3.17), membrane_offset=50 nm
  (¼ reach). **PI items still pending (the bleb/reservoir paths stay DISABLED until ratified):**
  - **σ_crit_bleb** (Tinevez 2009 critical cortical tension for bleb growth; MCF7 uncertain)
    is None → MembraneTetherUpdater (bleb nucleation) raises. Also the Bell-Evans rupture
    loop in `act()` is an unimplemented TODO (refused unconditionally). Anchor σ_crit_bleb
    AND implement the rupture loop before enabling the bleb-nucleation path.
  - **f_excess** (membrane reservoir excess-area fraction; Raucher-Sheetz 1999 / Figard 2014
    report a few %–tens of %; MCF7 uncertain) is None → `released_area()` raises; the
    tension-buffering reservoir-release path stays disabled.
  - ⚠ **membrane_offset / n_mem_nodes are mesoscale discretisation knobs** (geometry,
    grid-invariant; membrane_offset=¼ acceptor reach, n_mem_nodes=2000 shell nodes). Not
    physiological magic numbers — but flag if a specific MCF7 membrane area / node density
    is later required.
- **cortical_tension γ-denylist → registry-driven** (`NONCORTICAL_COMPARTMENT_PREFIXES`
  = `REGISTRY.gamma_denylist()` minus `cortex_*`). Unblocks the no-contamination control
  for every γ-contaminating compartment. ⚠️ sibling Gate-B owns this file → merge-coordinate.

## ✅ ventral_stress_fibers (③) → LIVE (passive) 2026-06-09 — FA-pair decision RESOLVED by PI
PI ratified **long-axis-aligned** basal FA pairing (2026-06-09). SF is now LIVE
(PASSIVE backbone): `select_aligned_fa_pairs` (basal subset → in-plane PCA principal
axis → low/high-projection pairing) yields aligned bundles (mean |cos|≈0.93, ~8 µm),
per-bundle EXACT-r0 backbone (force-free, max strain 1.6e-14). Gate PASS (5 controls:
bundles assembled n_SF=20, aligned, force-free, cortical γ_soft IDENTICAL OFF/ON =
no-contam, off-identity). N_filaments=20 (Cramer 1997) config candidate. **Still
DEFERRED to PI for the ACTIVE phase:** the Kumar 2006 10-30 nN single-SF tension
band needs NMII via the **sf_myosin_* prefix split** in `cortex/myosin.py` (else SF
motors contaminate cortical γ) + the equilibrated build (equilibration prelude). The
PASSIVE backbone is LIVE; the active contractile gate is the follow-on. Historical
detail of the now-resolved pairing decision:
- **FA-pair selection geometry.** `generate_stress_fiber_layout` pairs FA endpoints
  by **random shuffle** (`rng.permutation(M)`). At the real footprint the 2150
  integrins span z = −7.35 … 1.6 µm, so random pairing yields bundles of wildly
  different lengths spanning the whole cell — NOT physiological basal ventral stress
  fibers (which are basal-plane, aligned, ~10 µm in a 15 µm cell). A force-free
  build-time gate (assembled + no-contam) would PASS regardless, but the bundle
  GEOMETRY would be unphysical → marking SF LIVE on this would overstate fidelity.
  **Decision needed:** the vSF FA-pairing rule. Options: (a) nearest-basal-neighbour
  pairs at a physiological span (simplest defensible v1); (b) direction-aligned pairs
  along the cell long axis / traction field (more faithful, needs a polarity input);
  (c) keep random (rejected — unphysical). **Recommend (a)** as the v1 passive
  activation, (b) filed as the faithful follow-on. Implementing (a) is a small
  `stress_fibers.py` pairing-mode addition (no new physics).
- **Per-bundle k_actin / ell0.** μ_SF = N_filaments·EA_single is set (N_filaments=20
  candidate, Cramer 1997), but `k_actin` stays None and `ell0` varies per FA-pair
  distance, so `measure_sf_tension` needs the per-bundle registered r0 (the layout's
  `ell0_actin`), not a single k_actin — wire the measure to read `layout.ell0_actin`.
- Still also pending (from before): **N_filaments** PI-ratification (candidate 20),
  and the **sf_myosin_* prefix** split in `cortex/myosin.py` for the ACTIVE NMII
  phase (passive backbone needs neither; the Kumar 10-30 nN tension band is the
  active+equilibrated gate, deferred).
**Status:** RESOLVED — PI picked long-axis-aligned (b-style); SF is LIVE (passive).
Per-bundle EXACT-r0 backbone (force-free) implemented; measure_sf_tension prefix-
matches sf_actin_bond* for the deferred active gate.

### ⚙️ ACTIVE-NMII progress (Phase-3 ① in flight, 2026-06-09) + a NEW physics-design decision
**DONE (prerequisites, byte-identity-safe, committed):**
- **①a prefix-split** (`3c5795c`): `cortex/myosin.py` is now parameterized by a
  type-name `prefix` (default `"cortex_myosin_"`, byte-identical). An SF placement
  passes `prefix="sf_myosin_"` → its motor bonds fall under the existing `("sf_",)`
  γ-denylist, so SF NMII can never contaminate the cortical active-γ signal. This is
  the split the queue asked for. 42 myosin/SF tests + 11 new prefix tests pass.
- **①b actin-pool generalization** (`db00c20`): `MyosinStepUpdater` now takes an
  explicit `actin_pool_tags` (default `None`→`arange(n_cortex_actin)`, byte-identical),
  so the SAME Stam-Hocky/Hill machinery can bind the (non-contiguous) `sf_actin` chain
  instead of the cortex shell. 63 dynamics tests byte-identical + 2 new
  parity/confinement tests on a shifted pool.

**⛔ NEW PI DECISION NEEDED — SF contraction TOPOLOGY (blocks the live Kumar gate):**
The grip_walk NMII contraction (KU-3.5, PI-ratified 2026-05-31) is a **Stam-Hocky
BIPOLAR dipole**: `_bipolar_accepts` requires the +/− head-sets to grip **two
DIFFERENT (antiparallel) filaments** — that clause removes the zero-dipole degeneracy
in the cortex MESHWORK. But a ventral stress fiber is currently modeled as **ONE
`sf_actin` chain per bundle** (single filament), so the bipolar gate BLOCKS both
head-sets from gripping the same bundle → no contraction. The single-chain bundle
cannot host the antiparallel-filament dipole as written. Options:
- **(a) Antiparallel sub-chains** — split each bundle into 2 (or more) interdigitated
  antiparallel `sf_actin` sub-chains (true sarcomeric SF: NMII bridges antiparallel
  actin, Hotulainen-Lappalainen 2006). MOST FAITHFUL; a `stress_fibers.py` actin-
  scaffold topology change (doubles SF actin, re-does the backbone/anchor wiring) +
  the bipolar gate then works UNCHANGED. **Recommend** for fidelity.
- **(b) Same-filament center-walk** — a `cortex/myosin.py` SF-mode flag that lets the
  two head-sets grip the SAME chain at different positions and walk toward the
  minifilament centre (shorten the bundle). Smaller change, but a NEW contraction
  variant distinct from the ratified cortex bipolar gate → needs a modeling sign-off.
- **(c) Lumped cable tension** — REJECTED (violates the no-lumped-mechanism rule).
The PLACEMENT of `sf_myosin_*` minifilaments (where the beads sit, force-free) is the
same under (a)/(b); the actin SCAFFOLD differs, which is why I'm holding it for the
decision rather than building on an undecided topology. Also still pending: the Kumar
gate needs an **equilibration prelude** (raw full-cell SF run trips the BAOAB guard;
`cell/equilibration.equilibrate_cell` exists), `N_filaments` ratification (→ μ_SF),
and an explicit `k_actin` for `measure_sf_tension` (currently None → raises).

## ✅ cadherin_junction (④) → LIVE 2026-06-09 (FIRST multicell; two-cell doublet)
GATE-J PASS via the NEW two-cell assembler `cell/doublet.py::build_cell_doublet`
(two cortex shells facing across an interface; cadherins seeded on matched facing
caps; trans-dimers SEEDED pre-bound = the engaged-junction baseline; the catch-slip
binder maintains them). Build-time controls: 60 trans-dimers ALL A↔B (cross=60,
intra=0 — genuine two-cell junction), force-free (strain 1e-15), cadherin_ excluded
from cortical γ (denylisted), binder attached. Registry LIVE with `manifest_path=None`
(built by `build_cell_doublet`, NOT the single-cell loader). **Still DEFERRED to PI:**
- **Dynamic catch-slip maintenance** (bound-fraction φ ≈ k_on/(k_on+k_off) equilibrium)
  + the **Iturri ~6.5 nN ensemble de-adhesion** observable need an equilibrated run
  (the raw doublet run trips the BAOAB guard — equilibration prelude required; the
  binder's `n_sub` overflow guard was added defensively for that path).
- **Iturri-2020 SourceEvidence** (Cells 9(4):935) is UNREGISTERED in the Notion SoT —
  register before any deliverable cite (it sets the n_cad scale bridge + the ensemble
  de-adhesion report).
- **Scope:** GATE-J is a CORTEX doublet (two cortices + the junction). Per-cell full
  internal compartment stacks on a doublet (each cell with nucleus/MT/IF/membrane/SF)
  is a follow-on (run the single-cell extenders per cell of the doublet).

## ✅ junctional_actin (⑤) → LIVE 2026-06-09 (LAST; reserved STUB build implemented)
The reserved enabled+anchored build path is now IMPLEMENTED: built on the doublet
(`build_cell_doublet(with_junctional_actin=True)`) — one `junc_actin` head per
interface cadherin within the α-catenin reach of a SAME-cell cortex bead, a force-free
`junc_actin_anchor` (head↔cadherin) + per-r0-bin `junc_actin_couple_b{i}` (head↔cortex).
GATE PASS (belt assembled, same-cell coupling cross_cell=0, force-free, junc_actin_
γ-excluded, biphasic catch F*=6 pN). Catch-set are PI-CANDIDATES in the config
(Buckley x_catch/x_slip SOLID; k_catch0/k_slip0 ORDER; k_couple/k_anchor/k_on/
max_couple_dist/anchor_r0 DERIVED H.3); module default None → un-anchored build raises.
**DEFERRED to PI:** (1) ratify the catch-set / sanction the H.3 transfer + the
parallel-Pereverzev-vs-Buckley-sequential form distinction; (2) the DYNAMIC catch-slip
maintenance (JunctionalActinCouplingUpdater as a live Action); (3) ⚠ SPARSE-belt
fidelity — the single-particle cadherin is the ectodomain tip at the interface
(~0.5 µm from cortex), so only tips within the catch reach couple; a faithful DENSE
belt needs a cadherin-tail particle near the cortex (follow-on).

## 🏁 ACTIVATION BACKLOG COMPLETE (2026-06-09)
ALL 8 default-OFF compartments are LIVE (osmotic · MT · IF · LINC · membrane_reservoir ·
ventral_stress_fibers · cadherin_junction · junctional_actin). No EXPERIMENTAL/STUB
remains. Everything below + the per-compartment DEFERRED notes are active/dynamic-phase
work or PI parameter ratifications, NOT activations.

## ⚠️ NEW BLOCKER found during activation
- **MTOC single-hub degree vs HOOMD nlist exclusion cap (7).** The aster's single MTOC
  carries n_mt backbone bonds; with the full-cell LJ nlist on (`exclusions=("bond","1-3")`),
  n_mt>7 overflows the compile-time per-particle exclusion cap → `Too many bonds to process
  exclusions`. So a single-hub MTOC supports only n_mt≤7 (a sparse aster, now LIVE). A DENSER
  interphase aster (PI candidate 20, up to ~250) needs a **multi-bead MTOC core** in
  `cell/microtubules.py build_mt_topology` (distribute arms across ceil(n_mt/5) core beads) —
  a topology change that churns the 35 MT tests + 3 smokes. Filed for a focused follow-up.

---
# Original queue — 2026-06-09

Items the autonomous compartment-platform session (branch `h7/compartment-platform`)
could NOT do because they require **PI parameter ratification**, **shared-file edits**
(`cell/cell.py`, `cell/manifest.py`, `cortex/cortical_tension.py`, existing `cortex/*`
runtime, `integrator/`, `native/`), or a **gate-contract change** — all out of this
session's ownership. Each was SKIPPED here and recorded for PI. None block the
smoke-harness plumbing work (which used standalone harnesses + SMOKE-ONLY values).

The ENABLED-PATH smoke harnesses (`scripts/compartment_smoke/`) confirmed every
compartment's enabled build/step path runs (or is honestly blocked); the items below
are what stands between that plumbing and a real, PI-signed ACTIVATION.

## A. Shared-file edits required for activation (PI-gated; not owned here)

1. **cortical_tension γ-denylist extension (ALL γ-contaminating compartments).**
   `cortex/cortical_tension.py` `ADHESION_BOND_TYPE_PREFIXES` is hardcoded to FA
   prefixes only. Activating any of `ventral_stress_fibers (sf_, sf_myosin_)`, `linc
   (linc_)`, `intermediate_filaments (if_)`, `microtubules (mt_)`, `membrane_reservoir
   (mem_)`, `cadherin_junction (cadherin_)`, `junctional_actin (junc_actin)` needs its
   prefix added (cleanest: have the estimator read `REGISTRY.gamma_denylist()`), with a
   contamination regression test. **Shared file → PI sign-off.**

2. **stress_fibers BLOCKER #2 — `sf_myosin_*` prefix split.** SF NMII currently reuses
   the cortex `cortex_myosin_*` bond types (the active-γ signal). A LIVE SF build would
   contaminate cortical γ. Needs `cortex/myosin.py` parameterized to mint a distinct
   `sf_myosin_*` prefix. **Shared file (`cortex/myosin.py`) → PI sign-off.** Until then
   only the PASSIVE SF backbone is runnable (smoke#6 ran passive, NMII off).

3. **5-point wiring into the loader (every compartment).** `manifest.resolve_baseline`
   + `cell.py` (`_extend_snapshot_*`/`attach_*` + `CellBuildOptions` flag + `gamma_map`
   entry) + `configs/mcf7_baseline.yaml` `optional_subsystems` block + registry
   `status→LIVE`/`manifest_path`. `cell/cell.py` + `cell/manifest.py` are not owned.
   **PI-gated** (gate-contract: status flip).

## B. PI parameter ratifications (None-gated; smoke used SMOKE-ONLY placeholders)

These keep each enabled build raising until ratified (DETAILED_PLANS has candidates):
- **ventral_stress_fibers:** `N_filaments` (→ μ_SF); bundle bending EI (deferred).
- **linc:** `k_linc` (route A/B/C; smoke used 1e-2 N/m route-A candidate).
- **intermediate_filaments:** the nonlinear strain-stiffening Table curve (smoke ran the
  linear path only; nonlinear stays NotImplementedError).
- **microtubules:** `Y_stretch` (the dt lever) + `L_mt` host geometry; DI rates (default-OFF).
- **osmotic_regulation:** `Lp` (smoke used 1e-12 MCF7/AQP5 candidate; default 1e-13).
- **membrane_reservoir:** `σ_crit_bleb` + `f_excess` (bleb + reservoir-release paths blocked).
- **cadherin_junction:** k_trans/r_bind derived (OK to run); Iturri-2020 SourceEvidence
  registration before any deliverable cite.
- **junctional_actin:** the α-catenin catch set `k_catch0`/`k_slip0` (ORDER_ESTIMATE).

## C. STUB build-path implementations (not owned / PI-sanctioned release)

1. **junctional_actin HOOMD build path.** `extend_snapshot_with_junctional_actin` raises
   `NotImplementedError` EVEN when anchored (the explicit head-placement + per-r0 coupling
   topology is reserved for the PI-sanctioned release). Smoke#8 ran the SCALAR laws only
   (catch shape F*≈6 pN, engaged f≈0.91) + asserted the build barrier. Implementing the
   build path is a STUB→EXPERIMENTAL graduation; depends on `cadherin_junction` GATE-J
   passing first (junctional_actin requires cadherin_junction).

2. **membrane_reservoir bleb rupture updater.** `MembraneTetherUpdater` raises
   `NotImplementedError` unconditionally (the Bell-Evans `act()` loop is unwritten AND
   σ_crit_bleb is None). The static tether mesh is fully runnable (smoke#7); the rupture
   path needs the act() implementation + σ_crit_bleb ratification.

3. **microtubules dynamic instability updater** (`MTDynamicInstability`) raises
   `NotImplementedError` (plus-end snapshot-rebuild unwired); default-OFF, separate PI task.

## D. Two-cell / doublet build (cadherin activation)

`cadherin_junction` real activation (GATE-J) needs a `build_cell_doublet()` two-cell
assembler — `cell.py` is single-shell. Smoke#5 used a minimal facing-patch standalone
build to validate the binder; the production doublet is a `cell.py` change → PI-gated.

---
*Logged by the autonomous platform session; see `PLATFORM_AUTORUN_LOG_2026-06-09.md`.*
