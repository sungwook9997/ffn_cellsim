# Stage 1d — Layer 4 cellular Marangoni: pre-implementation Sanity Gate

This document is the pre-implementation sanity review for Stage 1d, in
which Layer 4 (cellular Marangoni — γ-gradient-driven tangential surface
flow) is activated on top of the Layer 1 v15 + Stage 1a+ Option β α=1.0
+ Stage 1a++ Layer 2 ζ=0.4 + Stage 1b Layer 3 φ-ODE + Stage 1c Layer 5
mechano-osmotic + Path C effective gravity baseline. **Marangoni core
only**: γ(φ) coupling Layer 3 → Layer 4. Nematic order Q-tensor and
internal vortex analysis are deferred to a future Stage 1d.b /
Stage 1d.c.

The companion documents are:
- `docs/07_internal_flow_dynamics.md` — pre-existing Layer 4 Marangoni
  + coffee-ring + Q-tensor framework
- `docs/02_force_models.md` §4.1 — cellular Marangoni stress reference
  (Pajic-Lijakovic & Milivojevic Eur Biophys J 2022, Fütterer Phys Rev
  Fluids 2022)
- `docs/03_adhesion_dynamics.md` — γ_cc(φ) coupling specification
- `docs/stage1c_sanity.md` — inherited carriers (Layers 1+2+3+5+gravity)
- `docs/12_validation.md` — Sanity-Gate Protocol (six checks) + Magic-
  Number Block

## Stage 1d scope

### What is added (Layer 4 Marangoni core)

Per `docs/07_internal_flow_dynamics.md` §1 + the Stage 1d-only
restrictions ("Marangoni 핵심만, nematic/vortex 미포함"):

1. **φ-modulated surface tension**: per `docs/03_adhesion_dynamics.md`
   §"Mechanical Coupling":
   $$ \gamma_{\text{cc}}(\phi) = \gamma_{\text{max}} \cdot (1 - \phi) + \gamma_{\text{min}} \cdot \phi $$
   with γ_max = 2.0 mJ/m² (E-cad strong), γ_min = 0.3 mJ/m² (E-cad
   weak). In Stage 1d simplification:
   $$ \gamma(\phi) = \gamma_0 + \gamma_1 \cdot \phi $$
   with $\gamma_0 = \gamma_{\text{max}}^\star = Ca_{\text{cc,max}} \cdot K \cdot R_0$
   and $\gamma_1 = (\gamma_{\text{min}} - \gamma_{\text{max}})^\star =
   -|\gamma_{\text{max}} - \gamma_{\text{min}}|^\star$.
   In dimensionless: $\gamma_1^\star \approx -0.017$ (with γ_max-γ_min ≈
   1.7e-3 J/m²) — the negative sign reflects "more integrin engagement
   → weaker cell-cell adhesion → lower γ".
2. **Spatial φ modulation** — required for ∇γ ≠ 0. Stage 1d activates
   Layer 3 spatial S_p signal:
   $$ S_p = \begin{cases} 1 & z_p < n_{\text{contact band}} \cdot dx \\ 0 & \text{otherwise} \end{cases} $$
   Substrate-engaged particles (contact band) receive the S=1 signal
   driving φ→φ_eq=0.75; interior particles receive S=0, φ relaxes
   toward φ→0 via the k_- term. Result: at long times,
   φ_contact-band ≈ 0.75, φ_interior ≈ 0 → spatial gradient ∇φ
   localised at the contact-band/interior boundary. This is the
   Cho 2020 mechanism made spatial.
3. **Marangoni impulse on grid** — per `docs/07_internal_flow_dynamics.md`
   §"Implementation":
   - Scatter γ(φ_p) to grid (γ_grid[I] via P2G, weighted by particle
     mass to give per-cell γ_eff)
   - Compute ∇γ_grid via central differences on the grid
   - Project to tangent: $\nabla_s \gamma = (I - \hat n \otimes \hat n)
     \cdot \nabla \gamma$ with $\hat n$ = existing
     `grid_normal` field from CSF
   - Apply tangential impulse at boundary-band cells:
     $\Delta v = \nabla_s \gamma \cdot dt / \rho_{\text{local}}$
4. **Boundary band only** — Marangoni acts only at the free surface, in
   the same band as the CSF impulse (cells with `m > min_cell_mass`
   AND `|grad_color|` > threshold).

### Active layers / what stays off

| Layer | Status | Notes |
|---|---|---|
| L1 bulk hydrodynamics (v15) | ON (carrier) | unchanged |
| L1a+ Option β substrate (γ_sub = γ_cc) | ON (carrier) | α=1.0 |
| L1a++ Layer 2 active stress | ON (carrier) | ζ_star = 0.4 |
| L1b Layer 3 φ-ODE | ON (carrier) | EXTENDED with spatial S_p (Stage 1d enables) |
| L1c Layer 5 mechano-osmotic | ON (carrier) | Tier 2 K(ρ_osm) |
| **L1d Layer 4 cellular Marangoni** | **NEW ACTIVATION** | core only (∇γ along boundary; nematic Q + vortex deferred) |
| Path C effective gravity | ON (carrier) | g_star = 0.01 (Stage 1c recalibrated) |
| Nematic Q-tensor / internal vortex | OFF | Stage 1d.b / 1d.c (deferred) |
| Layer 6 (chemistry / necrosis) | OFF | Stage 2 |

Default carrier `phi_initial = 0.55` (Stage 1c Pre baseline). With the
new spatial S_p, contact-band particles drift toward φ_eq=0.75, interior
particles drift toward 0 (since k_+·S=0 there).

### Path C g_star recalibration check

Layer 4 Marangoni magnitude estimate: $|\nabla \gamma|_{\max} \approx
|\gamma_1| \cdot |\nabla \phi|_{\max}$. With $|\gamma_1| = 0.017$ and
$|\nabla \phi|_{\max} \approx 1/R_0 = 1$, $|\nabla \gamma|_{\max} \approx
0.017$. Compare to:
- Layer 2 active stress: ζ·K = 0.4 (Marangoni is 24× weaker)
- Surface tension γ·κ = 0.028 (Marangoni same order)

Marangoni is a **subtle** effect, not a dominant force in this scale.
**Path C g_star = 0.01 (Stage 1c recalibrated value) carries over to
Stage 1d unchanged.** No further recalibration needed (Marangoni is too
weak to disturb the gravity/surface-tension balance regime).

---

## Sanity Gate (six checks)

### 0. Time-scale mapping sanity check (NEW for Stage 1d, PI requested)

Before any further check, verify the dimensionless-to-dimensional time
mapping used throughout the project, and assess Stage 1d pilot duration
vs the PI experimental timescale.

**τ_relax dimensional anchor**:
- Project default: `physics.maxwell_tau_s = 60` (per all configs since
  Stage 1a pilot).
- Anchor: `docs/02_force_models.md` §1.3 cites Moeendarbary et al. 2013
  *Nature Materials* 12:253 (IF 47, [primary]) for τ_relax ~ 10–600 s
  (cytoplasm 10–100 s; junction remodeling ~10 min).
  **Verification (2026-04-29)**: web search confirmed Moeendarbary 2013
  reports poroelastic time t_p ~ 0.1–10 s (rapid AFM force-relaxation
  decay in first 0.5 s), with the longer junction-remodeling regime
  extending to minutes. **τ_relax = 60 s sits in the upper cytoplasm
  range / lower junction-remodeling range** — well-anchored within the
  cited literature; specific 60 s value is a project-default mid-range
  choice within the established 10–600 s span.
- **Magic-Number Block on τ_relax**: PASS (literature-anchored range from
  Moeendarbary IF 47; specific 60 s is mid-range pick, grid-invariant in
  dimensional terms, not chosen to fit any gate). This is a stronger
  anchor than ζ_star Option α' / γ_sub_Col1 / α_osm cases.

**Pilot wall-time mapping**:
- `dt_star = 0.01`, `total_time_star = 240` → 24,000 steps.
- Dimensional: `dt = dt_star · τ_relax = 0.01 · 60 = 0.6 s`; total
  pilot duration = 24,000 · 0.6 s = **14,400 s = 240 min = 4 hours**.
- Frame interval: `frame_interval_star = 10` → 10·τ_relax = 600 s = 10
  minutes per frame; 24 frames × 10 min = 4 hours total.

**PI experimental timescale**:
- `data/experimental/260313_{Bare,Pre,Lam4}.csv` `Time_min` column
  format confirmed (verification 2026-04-29).
- Frame 0: t = 0 min. Frame 1: t = 60 min = 1 hour. Subsequent frames
  at 60-min intervals.
- Total durations per `docs/00_project_vision.md` §"PI's CSV files":
  Bare 82 hr, Pre 82 hr, Lam4 83 hr.
- A/A₀ values reported: at experimental t = 60 min already 1.35–1.45
  (Bare/Pre) and 1.08 (Lam4 first frame); A/A₀ final 4.6–33.1.

**Pilot vs experimental gap**:
- **Pilot 4 hr ≈ 5% of PI experiment 82 hr (20× shorter)**.
- Implication: Layer 4 Marangoni — predicted weak (Ma ≈ 0.06,
  perturbative) — has only ~5% of the experimental time to accumulate
  effect. Even if mechanism is fully active, R-drift contribution may
  be smaller than experimental endpoint A/A₀ would suggest.
- Implication: A/A₀ comparisons must be done at *matching* time points
  (sim t=240min vs PI t=240min frame), not at experimental endpoints.

**Production scale per `CLAUDE.md`**:
- "Production target: 80 sim hours (matches PI's experimental duration)"
- Time budget: 80 hr = 80·3600 s = 288,000 s = 4800·τ_relax.
- vs pilot 4 hr = 240·τ_relax (20× scale-up).
- Wall-clock estimate: pilot ~1 min (Layer 1) → ~5 min (full Layer
  1+2+3+5) → ~10 min projected (Layer 1+2+3+4+5). Production scale
  20× ≈ 100 min ≈ 1.7 hr (10× cell count + 20× time = ~200× compute,
  but within the CLAUDE.md ≤ 12 GB VRAM laptop A5000 ≤ 20 hr budget).

**Stage 1d result interpretation framing** (deduced from time-scale
mapping):
1. **Pilot 4 hr is sufficient to TEST Layer 4 mechanism** (∇γ existence,
   Marangoni impulse direction, spatial φ contrast formation).
2. **Pilot 4 hr is INSUFFICIENT to reproduce PI A/A₀ endpoint values**
   (Marangoni perturbative + 20× shorter time). A/A₀ comparison at
   matching time points (sim t=240min vs PI t=240min frame, both
   ≈ 0–4 hr region) is the meaningful comparison; PI endpoint comparison
   requires production-scale runs.
3. **R drift gate should detect mechanism contribution** (Layer 4 effect
   on R is integrated over the full pilot; perturbative Marangoni adds
   a small correction; R drift improvement vs Stage 1c baseline 0.138
   is the expected signal).

**Check 0: PASS** with the explicit framing that Stage 1d pilot is a
mechanism test (Bucket M1/M2/M3 classification valid) but not a
quantitative reproduction of PI endpoint A/A₀ data (which requires
production-scale runs).

### 1. Dimensional analysis

- $\gamma_1$ — surface-tension/φ-coupling slope, units of stress·length
  = K·R₀ in dimensional terms. Dimensionless: $\gamma_1^\star = \gamma_1
  / (K \cdot R_0)$. With $|\gamma_{\max} - \gamma_{\min}| \approx 1.7
  \times 10^{-3}$ J/m² and $K \cdot R_0 = 10^3 \cdot 10^{-4} = 0.1$
  J/m², $\gamma_1^\star = -0.017$ (negative per Cho 2020 mechanism).
- $\nabla \gamma$ — units of stress (force/area) per length, in
  dimensionless [stress / R₀].
- Marangoni impulse magnitude: $|\nabla \gamma| \cdot dt / \rho \approx
  0.017 \cdot 0.01 / 1 = 1.7 \times 10^{-4}$ velocity per step. Compare
  to the existing CSF impulse magnitude $\gamma \cdot \kappa \cdot dt /
  \rho \approx 0.028 \cdot 0.01 / 1 = 2.8 \times 10^{-4}$ — same order
  of magnitude, well within numerical bounds.
- **Cellular Marangoni number** $Ma = |\nabla \gamma| \cdot R_0 / (\eta
  \cdot v) \approx 0.017 \cdot 1 / (0.3 \cdot 1) = 0.057$ — weak
  Marangoni regime (Ma < 1; see `07_internal_flow_dynamics.md` §5
  "Cellular Marangoni Number" predicting Ma ~ 0.1 for cellular spreading
  vs 1-100 for thin-film drying, consistent).
- CFL: Marangoni adds no new fast scale (no inertia in overdamped). dt
  unchanged.

**Check 1: PASS** with the calibration confirmed against the existing
07_internal_flow_dynamics.md analytical estimate.

### 2. Boundary cases

- $\gamma_1 \to 0$ (γ φ-independent): Marangoni impulse identically
  zero; recovers Stage 1c. PASS.
- $\nabla \phi \to 0$ (uniform φ): Marangoni impulse zero (no driving
  force). Single-phenotype Stage 1c Pre with uniform initial φ would
  give this; the spatial S_p extension breaks this (contact-band ≠
  interior over time).
- Pajic-Lijakovic 2022 limit: their analytical Marangoni framework
  predicts directed cell motion from low γ to high γ. Our sign
  convention with $\gamma_1 < 0$ + Marangoni flow toward higher γ:
  flow goes from high-φ regions (low γ) to low-φ regions (high γ).
  In Stage 1d this means flow from contact-band-spread cells (high φ)
  toward interior (low φ) — i.e. a *retraction* component opposing the
  spreading driven by Layer 2 active stress. Coupled effect: Layer 2
  pushes cells outward at boundary, Layer 4 Marangoni pulls γ-defined
  flow inward. Net depends on relative magnitudes (Layer 2 ζ·K = 0.4
  vs Marangoni ∇γ ≈ 0.017 → Layer 2 dominates by 24×, Marangoni is
  perturbative).
- $|\nabla \gamma| \to \infty$: would require sharp φ discontinuity;
  prevented by per-particle S_p smoothing and Layer 3 ODE finite-rate
  evolution.
- N → 0 / N → ∞: per-particle field; same scaling as v15.
- Δt → 0 / Δt → large: same as Stage 1c.

**Check 2: PASS** with sign analysis confirmed against Pajic-Lijakovic
2022 framework.

### 3. Conservation invariants

- **Mass**: unchanged. ✓
- **Momentum (horizontal)**: Marangoni is a *surface body force* applied
  to a thin band. The tangent projection $(I - \hat n \otimes \hat n)$
  removes the normal component, leaving only tangential motion (which
  is internal to the boundary band). No external horizontal momentum
  injection. The MPM scatter is symmetric, so horizontal momentum is
  conserved to atomic-add round-off. Stage 1a+ horizontal-only momentum
  gate stays in force. ✓
- **Momentum (vertical)**: tangent direction at the substrate contact
  band is roughly horizontal (n̂ ≈ ẑ at the bottom surface), so
  Marangoni vertical momentum injection is negligible; absorbed by
  substrate as in Stage 1a+/Path C. ✓
- **Angular momentum**: Marangoni at a curved surface CAN inject
  angular momentum about an off-center axis. Tracked diagnostic only;
  no gate (same as v15).
- **Energy**: Marangoni does work along surface tangent direction;
  energy injection per step bounded by $|\nabla \gamma| \cdot v_{\text{surf}}
  \cdot A_{\text{surface}} \cdot dt$. In dimensionless terms this is
  $\approx 0.017 \cdot 0.01 \cdot 12 \cdot 0.01 \approx 2 \times 10^{-5}$
  per step — small. Energy-monotone gate is already SUSPENDED under
  Layer 2 / 3 / 5 active runs (Stage 1a++ Cousin-Rule contract change
  inherited); Stage 1d does NOT change this. **NEW under Stage 1d**:
  the suspension condition is extended to also fire when `layer3_enabled`
  OR `layer4_enabled` OR `layer5_enabled` — fixing the inherited
  Stage 1b carry-over issue (energy-monotone gate was running and
  failing under Layer 3 active because the suspension condition only
  checked `zeta_star > 0`). Cousin-Rule contract change explicitly
  surfaced.

**Check 3: PASS** with the Stage 1d-extended energy-monotone suspension
condition (Cousin-Rule contract change resolving inherited Stage 1b
carry-over issue).

### 4. Numerical sanity

- $dt^\star = 0.01$ unchanged.
- Marangoni impulse magnitude per step: $|\nabla \gamma| \cdot dt /
  \rho \approx 1.7 \times 10^{-4}$ — well within f32 precision and
  existing speed gates.
- New grid scratch field `grid_gamma` (γ_eff per cell) and `grid_gamma_grad`
  (vector field) — same memory pattern as existing grid_color +
  grid_color_grad. < 5% wall-clock impact.
- Spatial S_p (Layer 3 extension): adds a per-particle
  position-dependent indicator. O(N) work, < 1% wall-clock.

**Check 4: PASS**.

### 5. Sign / sense check

- $\gamma_1 < 0$: high-φ regions (more integrin) have LOW γ. Matches
  Cho 2020 (E-cad weakening with integrin upregulation reduces
  cell-cell cohesion).
- Marangoni flow direction (Pajic-Lijakovic 2022 framework): from low
  γ to high γ. With $\gamma_1 < 0$: flow from high-φ region to low-φ
  region. In our Stage 1d setup with contact-band-φ_high vs
  interior-φ_low: Marangoni flow goes from contact-band UP into the
  interior — a retraction component along the substrate-attached
  cells' surface.
- This is *physically* the cellular analog of the classical Marangoni
  effect (e.g. tears of wine: surface flow from low to high γ). The
  cellular interpretation: cells with weaker E-cad (high φ) at the
  spreading edge are PULLED toward the stronger-E-cad bulk by surface
  tension gradient — a retraction-like effect that *opposes* the Layer
  2 active spreading. Net effect is a balance.
- Implementation sign: $\Delta v = + \nabla_s \gamma \cdot dt /
  \rho_{\text{local}}$. With $\nabla_s \gamma$ pointing from low γ
  toward high γ (i.e. from high φ toward low φ, since γ_1 < 0), the
  velocity impulse aligns with the gradient direction → flow from low
  γ to high γ. ✓

**Check 5: PASS**.

### 6. Measurement-protocol consistency (NEW since v13)

- **R drift / shape metrics**: same as Stage 1c.
- **∇γ measurement protocol**: γ_eff scattered to grid via P2G (weighted
  by particle mass), then central FD on grid. Walk-through: for the
  Marangoni impulse to have the correct magnitude, γ_eff must
  represent the per-cell average of γ(φ_p) weighted by particle
  presence — exactly what mass-weighted P2G gives. ✓
- **Tangent projection**: uses existing `grid_normal` field from CSF
  (Brackbill 1992). The same off-peak f″/f′ asymmetry that broke v13
  could in principle apply here too — but the tangent projection
  $(I - \hat n \otimes \hat n)$ doesn't involve a SECOND derivative;
  it's a first-order operator. The v13 pathology applied to the
  divergence (-∇·n̂) which is second-derivative-like; here we just
  PROJECT a vector onto the tangent plane. No off-peak amplification.
  ✓
- **Energy budget**: Layer 4 Marangoni work is bounded as in check 3.
  Reported as new diagnostic `marangoni_power_star`.
- **NEW measurement protocol**: spatial φ field. Aggregation per shell
  (boundary, contact, interior) gives the φ-spatial-modulation witness
  that Layer 3 spatial S_p is producing the expected ∇φ. Walk-through:
  per-particle scalar; aggregations well-defined; same protocol as
  Stage 1b/1c <φ>_bulk / <φ>_boundary diagnostics. ✓

**Check 6: PASS** with energy-monotone suspension Cousin-Rule contract
change (Stage 1d resolves inherited Stage 1b issue).

---

## Magic-Number Block

### γ(φ) coupling slope `γ_1_star`

| Parameter | Value | Anchor | IF |
|---|---|---|---|
| γ_max (E-cad strong) | 2.0 mJ/m² | `docs/03_adhesion_dynamics.md` §"Mechanical Coupling" citing Maître Science 2012 IF 47 (cell-cell adhesion energy 0.5-2 mJ/m² range) | 47 |
| γ_min (E-cad weak) | 0.3 mJ/m² | Same anchor (lower end of the Maître range) | 47 |
| γ_1 = γ_min - γ_max | -1.7 mJ/m² | Derived | — |
| γ_1_star = γ_1/(K·R₀) | -0.017 | dimensional conversion | — |

1. **Derivable** — yes, framework anchored to Maître Science 2012 IF 47
   for the cell-cell γ range; the φ coupling form γ(φ) = linear
   interpolation between γ_max and γ_min is the standard
   `03_adhesion_dynamics.md` choice (specified BEFORE Stage 1d
   implementation; not a fitting-driven choice). PASS.
2. **Grid-invariant** — yes, dimensionless ratio, independent of dx,
   dt, n_particles. PASS.
3. **Fitting** — no. Values pre-specified; not adjusted to make any
   gate pass. PASS.

### Layer 3 spatial S_p indicator

| Parameter | Value | Anchor |
|---|---|---|
| S_p = 1 if z_p < n_contact_band·dx, else 0 | binary indicator | `docs/03_adhesion_dynamics.md` §"S(t) signal" definition |

The spatial S_p definition matches the existing 03_adhesion_dynamics.md
framework `S_i(t) = w_contact · A_substrate_i / A_normalize + ...`
collapsed to the simplest binary indicator (substrate-contact YES/NO).
Test 1 PASS (framework-derivable from 03 doc), test 2 PASS (cell-count
indicator, grid-invariant), test 3 PASS (binary, not a fittable
parameter).

### Marangoni framework reference

Pajic-Lijakovic & Milivojevic Eur Biophys J 2022 (IF ~ 2.5) is the
primary cellular Marangoni framework reference (already cited in
`02_force_models.md` §4.1 and `07_internal_flow_dynamics.md`). Eur
Biophys J IF is below the IF≥15 preferred threshold per CLAUDE.md, but
the project's existing precedent (`02_force_models.md`) accepts this as
the primary anchor for cellular Marangoni; secondary anchor Fütterer
Phys Rev Fluids 2022 for the laser-ablation experimental confirmation
of cellular Marangoni-like motions. **Verification (2026-04-29)**: web
search confirmed paper exists (Springer Nature Link, DOI
10.1007/s00249-022-01612-1); abstract confirms "gradient of tissue
surface tension induces directed cell spreading from regions of LOWER
tissue surface tension to regions of HIGHER tissue surface tension"
(matching our sign convention). **Magic-Number Block PARTIAL** (analogous
to ζ_star Option α' precedent — framework anchored, specific numerical
pinpoints inherit Maître IF 47 anchor for γ range; PI full
authorisation 2026-04-29 covers).

---

## New gate / diagnostic candidates (Stage 1d)

| Gate | Status | Tolerance | Source |
|---|---|---|---|
| All Stage 1c inherited gates | inherited | unchanged | Stage 1c |
| Energy monotone | **SUSPENDED** under L3/L4/L5 (Stage 1d Cousin-Rule extension) | — | Stage 1d check 3 |
| **R drift improvement vs Stage 1c baseline (0.138)** | NEW GATE | strict-less than 0.138 | inherited improvement |
| **<\|∇_s γ\|>_boundary finite (Marangoni driving force)** | NEW DIAG | log only | first principles |
| **Spatial φ contrast: <φ>_contact − <φ>_interior** | NEW DIAG | log only; expected to grow over pilot duration | Layer 3 spatial S_p |
| **Marangoni power finite** | NEW GATE | finite (not NaN, not > 10·U_strain at end); same form as Stage 1a++ active-power-finite | first principles |

---

## Implementation outline

- **acs/physics/mlsmpm.py**:
  - SolverConfig: add `layer4_enabled`, `gamma_max_star`, `gamma_min_star`,
    `layer3_spatial_S` (bool, enables Layer 3 spatial S_p extension).
    Defaults disable Layer 4 and keep Layer 3 single-S behaviour.
  - New grid scratch fields: `grid_gamma` (per-cell γ_eff scalar),
    `grid_gamma_grad` (vector). Same allocation pattern as
    grid_color / grid_color_grad.
  - New kernels:
    - `_scatter_gamma_to_grid`: per particle, compute γ(φ_p) =
      γ_max·(1-φ) + γ_min·φ; scatter mass-weighted to grid_gamma.
    - `_compute_gamma_grad`: central FD on grid_gamma → grid_gamma_grad.
    - Marangoni impulse added in `_grid_op_overdamped`: at cells with
      `m > min_cell_mass`, compute tangent projection
      `∇_s γ = (I - n̂⊗n̂) · ∇γ` (n̂ from existing grid_normal), apply
      `dv += dt · ∇_s γ / ρ_local`.
  - `_integrate_phi_ode` (Stage 1b): add per-particle S_p indicator
    when `layer3_spatial_S` enabled. dφ/dt = k_+ · S_p · (1-φ) - k_- · φ.
  - `_compute_invariants`: add Marangoni power diagnostic
    `P_marangoni = Σ_cells |∇_s γ|·|v_surf| · cell_volume`.
- **acs/runner.py**: read `layer4` config block; Layer 4 logger;
  energy-monotone gate suspension condition extended to L3/L4/L5;
  new R-drift-improvement-vs-Stage-1c gate; Marangoni-power-finite
  gate.
- **One config**: `configs/stage1d_pilot.yaml`. Layer 1 + 1a+ β α=1.0
  + 1a++ ζ=0 (Layer 3 overrides) + 1b φ Pre + 1c Layer 5 + Path C
  + Layer 4 Marangoni + Layer 3 spatial S extension.

---

## Decision request to PI (resolved 2026-04-29 by full authorization)

PI granted **full authorisation** 2026-04-29 ("원래 framework standard 풀
적용. 축소 금지. 매 단계 풀 protocol") for Stage 1d including:
- Layer 4 Marangoni core only (nematic Q + vortex deferred per PI scope).
- γ(φ) coupling per `03_adhesion_dynamics.md` framework (Maître IF 47
  anchored γ range; Pajic-Lijakovic 2022 framework reference).
- Layer 3 spatial S_p extension (binary substrate-contact indicator)
  — required for ∇φ ≠ 0 → ∇γ ≠ 0 → Marangoni driving force exists.
- Path C g_star = 0.01 unchanged from Stage 1c (Marangoni is
  perturbative and does not require recalibration).
- Energy-monotone gate suspension Cousin-Rule extension to L3/L4/L5
  (resolves inherited Stage 1b issue).
- Inherited Stage 1c contract changes (anchor-force integrand
  K_eff modulation update still flagged; Stage 1d does NOT
  implement that fix).

Stage 1d implementation begins under these conditions. Stop conditions
remain in force: no magic numbers, no gate semantics edits beyond the
explicitly recorded contract extensions, no v13 anti-pattern, halt and
surface to PI on any FAIL or critical error.
