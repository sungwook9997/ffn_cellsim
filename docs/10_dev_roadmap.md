# 10 — Development Roadmap (Bottom-Up)

This is the canonical "what to do next" document. **Claude Code reads this to determine the current milestone.** Each Stage has explicit entry conditions, deliverables, and validation tests. Do not advance to the next stage until validation passes.

---

## Stage 0 — Environment Setup
**Goal**: Reproducible, GPU-portable Python environment.

### Tasks
- [x] `requirements.txt` or `pyproject.toml` with pinned versions (Taichi >= 1.7, numpy, scipy, matplotlib, plotly, vispy, h5py, pyyaml, etc.)
- [x] Conda env `activecellsim` (or similar)
- [x] Verify Taichi GPU detection on RTX A5000 (`ti.init(arch=ti.gpu)`)
- [x] Verify GPU portability flag (env var `ACS_GPU_BACKEND` switches between cuda/vulkan/cpu)
- [x] `configs/` folder with `pilot.yaml`, `production.yaml`, `dev.yaml`
- [x] Logging setup (`logs/` dir, structured logging)
- [x] Pre-commit hook or linter (ruff/black) — optional but recommended

### Validation
- `python -c "import taichi as ti; ti.init(arch=ti.gpu); print(ti.cfg.arch)"` returns the GPU backend
- `pytest` runs (even with no tests yet — confirms framework)

**Stage 0 completed 2026-04-29.**
Verified on primary workstation: Python 3.11.15, Taichi 1.7.4, `Arch.cuda` detected on
RTX A5000 Laptop GPU (16 GB). `ACS_GPU_BACKEND=cpu` correctly falls back to `Arch.x64`.
12/12 pytest smoke tests pass; ruff clean. See `scripts/verify_env.py` for the
reproducible Stage 0 check. Note: CLAUDE.md states "RTX A5000 (24GB VRAM)" but the
actual primary GPU is the *Laptop* variant with 16GB — surfaced to PI for review.

---

## Stage 1a — Single Spheroid Equilibrium (Layer 1 minimal)
**Goal**: A single 3D spheroid sits in vacuum, holds together, reaches mechanical equilibrium — and we have empirical confidence the same code is safe to run at full production scale and duration.

### Numerical scheme decisions (locked at Stage 1a entry)
- **Solver baseline**: Taichi `examples/simulation/mpm99.py` extended to 3D (or `mpm3d.py` if present); APIC + quadratic kernel = MLS-MPM up to constants.
- **Constitutive integrator**: **exponential integrator** for the Maxwell deviatoric stress —
  `σ(t+Δt) = σ(t)·exp(-Δt/τ) + 2μ·τ·ε̇·(1 - exp(-Δt/τ))`. Closed-form, unconditionally stable in stiff τ, removes the need for any "reduced form" fallback.
- **Free-surface detection**: **density-based (option A)** — boundary particle iff local SPH-style kernel density `ρ̂_p < 0.6·ρ_bulk`. Reasoning: simplest, matches MLS-MPM particle data flow naturally, no extra reconstruction pass.
- **Surface tension imposition**: **Continuum Surface Force (CSF, Brackbill 1992) on the background grid**, with boundary particles as the only sources of the colour-function gradient. Standard, stable, and grid-resident which keeps it inside the existing P2G/G2P pipeline.

### Parameter anchoring policy (Stage 1a)
- **Stage 1a uses literature mid-range placeholders**, not formal IF≥15 anchors. Goal of Stage 1a is *equilibrium + numerical stability*, where relative ratios and stability dominate over absolute values.
- Placeholders (every site marked `# TODO[anchor]`):
  - Maxwell relaxation `τ = 60 s`
  - Cortical effective bulk `K = 1 kPa`
  - Cohesive surface energy `γ_cc = 1 mJ/m²`
  - Mass density `ρ = 1.05 g/cm³`
- **Formal IF≥15 anchoring runs as a separate gate between Stage 1a and Stage 1a+** (see `docs/02_force_models.md`, `docs/12_validation.md`). Stage 1a+ does not start until every `TODO[anchor]` is closed.

### Tasks
- [ ] Implement MLS-MPM 3D solver (`acs/physics/mlsmpm.py`) from the Taichi example baseline
- [ ] Maxwell viscoelastic constitutive via exponential integrator (`acs/physics/constitutive.py`)
- [ ] Density-based free-surface detection + CSF surface tension on grid (`acs/physics/cohesion.py`)
- [ ] Cortical effective bulk via volumetric Neo-Hookean term (in `constitutive.py`)
- [ ] Cell-equivalent material point initialization (random pack inside a sphere, ~5,000 points; pilot uses ~1,000)
- [ ] HDF5 frame export (`acs/io/hdf5_writer.py`) following `docs/13_data_schema.md`
- [ ] Conservation/diagnostics monitor (`acs/physics/invariants.py`) — mass, momentum, KE, strain energy each frame
- [ ] Vispy offline viewer → PNG sequence (`acs/viz/frame_export.py`); separated from sim process
- [ ] Shape-metrics post-processor (`acs/analysis/shape_metrics.py`) — sphericity (Wadell), R_g, R/R₀
- [ ] Unit tests: mass conservation, momentum conservation, monotone-energy decay, equilibrium sphericity

### Staircase burn-in protocol (mandatory)
**Pilot 30 min alone does not certify production-scale safety.** Run all four steps in order; each step is a PASS/FAIL gate; FAIL → halt and root-cause before advancing.

| Step                     | Config file                          | N points | Sim time | Wall-clock target | Purpose                                |
|--------------------------|--------------------------------------|----------|----------|-------------------|----------------------------------------|
| Stage 1a-pilot           | `configs/stage1a_pilot.yaml`         | 1,000    | 4 hr     | ≤ 30 min          | Code skeleton + baseline step rate     |
| Stage 1a-mid             | `configs/stage1a_mid.yaml`           | 1,000    | 24 hr    | ~1 hr             | Long-time numerical drift              |
| Stage 1a-prod-burnin     | `configs/stage1a_prod_burnin.yaml`   | 5,000    | 24 hr    | ~3–5 hr           | Scale + duration jointly stress-tested |
| Stage 1a-production      | `configs/stage1a_production.yaml`    | 5,000    | 80 hr    | 5–15 hr           | Full production-grade single-spheroid  |

Each step writes to `results/{run.name}/` (kept separate so a failing burn-in run never overwrites a passing pilot).

#### PASS/FAIL gate — measured at every step
1. **Conservation drift**
   - mass: `|Δm/m₀| < 1e-10` (exact in MLS-MPM; any drift is a code bug)
   - momentum (no external force): `|Δp/(m·v_rms)| < 1e-3` over the full run
   - energy: KE + strain + surface energy must be **monotonically non-increasing** within numerical noise (`+1e-3·E_max` tolerance for round-off)
2. **Step-time scaling vs N** — `t_step(N) / t_step(1000)` should be ≲ N/1000 (super-linear scaling = a quadratic neighbour search slipped in)
3. **GPU memory trajectory** — peak VRAM stable after the first 10 frames (no leak); production-scale step ≤ 12 GB on Laptop A5000
4. **Numerical stability** — zero NaN/Inf in positions, velocities, stress, deformation gradient; max velocity `< 10·v_rms` (no run-away particles)
5. **Equilibrium shape** (pilot+) — Wadell sphericity ψ ≥ 0.95 by `t = τ_relax`
6. **Geometric drift** — `|R(t)/R₀ − 1| < 0.05` after the initial relaxation transient. Detects volume drift bugs (J → ≠ 1 unphysical), constitutive sign errors, and CSF imbalance. The MLS-MPM scheme is incompressible only at K → ∞; a 5% radius envelope is a generous numerical margin against round-off.

A run-level gate report (`results/{run.name}/gate_report.md`) records every metric + PASS/FAIL verdict + any next-step blockers. The runner halts the staircase automatically on FAIL.

#### Pilot 30-min fallback ladder (apply only on the pilot step, in order)
1. Pin `device_memory_GB` so Taichi mallocs once
2. Drop background grid resolution from 128³ to 64³
3. Reduce N to 750 temporarily (revisit after Stage 1b enables further optimisation)
4. **If still > 30 min**: STOP and discuss with PI — the inner loop is too slow to iterate productively and we need to revisit the scheme before climbing the staircase.

---

## Stage 1a+ — Substrate Contact + Wetting Onset
**Goal**: Spheroid descends, contacts Col1-coated rigid substrate, partial wetting occurs.

### Tasks
- [ ] Gravity + buoyancy (medium density, cell density)
- [ ] Substrate plane (rigid, infinite, with adhesion energy density)
- [ ] Hertz contact mechanics at first touch
- [ ] Adhesion energy γ_substrate creating wetting force
- [ ] Spreading-coefficient-based asymptotic check

### Validation
- Initial contact area scales as Hertz prediction: a ~ (R·δ)^(1/2) where δ is indentation
- Long-term wetting follows Young-like equilibrium (deformed spheroid shape)
- Smooth transition (no instabilities) from free-floating to wetted

---

## Stage 1a++ — Layer 2: Boundary Biology
**Goal**: Lamellipodia, filopodia, discrete focal adhesions, leader cell formation.

> **Reopen the overdamped assumption here.** Stage 1a adopts overdamped dynamics
> (no inertia) on the basis that cellular spreading is Re ≪ 1 (Re ≈ 10⁻¹³, see
> `docs/12_validation.md`). Layer 2 introduces **lamellipodia events on the
> ms-second scale**, which can be fast enough that the overdamped approximation
> may break. Before activating Layer 2, re-derive Re and the local Stokes number
> for the lamellipodial event timescale and confirm overdamped is still valid;
> if not, switch to a semi-inertial scheme (e.g., add a small effective inertial
> term whose magnitude is justified by the lamellipodial mass × velocity scale).

### Tasks
- [ ] Free-surface (boundary material point) identification each step
- [ ] Stochastic lamellipodia events: rate depending on local stress, FA density, φ (later)
- [ ] Filopodia probing: directional, ECM-sensing
- [ ] Discrete FA patches: stochastic formation, force-dependent lifetime (catch/slip-bond)
- [ ] Leader cell stochastic polarization → enhanced active stress
- [ ] Visualization: lamellipodia/filopodia overlay on viewer

### Validation
- Leader cells emerge at ~1-5% of edge cells (literature consistent)
- Spreading rate increases vs Stage 1a+ baseline (active vs passive wetting)
- Edge protrusion events observable in frame sequence

---

## Stage 1b — Layer 3: Adhesion Network Dynamics (φ ODE)
**Goal**: ULA-like vs pV4D4-like phenotypic distinction emerges.

### Tasks
- [ ] Per-material-point φ state variable
- [ ] ODE `dφ/dt = k₊·S(t)·(1-φ) - k₋·φ` with literature-derived k₊, k₋
- [ ] φ modulates: γ_cc(φ), σ_a(φ), ρ_int(φ), K_cortex(φ) per `03_adhesion_dynamics.md`
- [ ] Two-anchor parameter set definitions (ULA-like, pV4D4-like)
- [ ] Cross-check predicted φ trajectory against Cho et al. 2020 Western blot timecourse
- [ ] Visualization: φ field colormap on material points

### Validation
- ULA-like simulation: φ stays low throughout (E-cad dominant maintained)
- pV4D4-like simulation: φ rises after ~24 sim hr (matches Cho et al. trend qualitatively)
- Two phenotypes show distinct spreading dynamics

---

## Stage 1c — Layer 5: Mechano-Osmotic Coupling (Tier 2)
**Goal**: Spreading-induced volume loss couples to mechanics.

### Tasks
- [ ] Per-material-point density ρ(t)
- [ ] ODE `dρ/dt = α·ε̇_spread - β·(ρ - ρ₀)` (literature-anchored: ~50% volume loss / hour for fast-spreading cells)
- [ ] ρ-modulated viscosity: η_eff(ρ) = η₀ (ρ/ρ₀)^n
- [ ] ρ-modulated cortex stiffness: K_cortex(ρ) = K₀ (ρ/ρ₀)
- [ ] Visualization: density field (drying-analog)

### Validation
- Cells in spreading region show density increase (mass conservation: water efflux)
- Cells in spheroid core (less spreading) maintain ρ ≈ ρ₀
- Mechanical stiffness gradient emerges naturally (stiffer at edges)

---

## Stage 1d — Layer 4: Internal Flow Dynamics
**Goal**: Cellular Marangoni, nematic order, vortex structure analysis.

### Tasks
- [ ] γ_eff(r,t) computed from local φ
- [ ] Tangential Marangoni stress at free surface
- [ ] Nematic order parameter Q from material-point polarizations
- [ ] Vorticity field on background grid
- [ ] Streamline / particle-tracer post-processing
- [ ] Coffee-ring metric: radial density profile
- [ ] Visualization: vector field, streamline, Q tensor ellipses

### Validation
- Internal toroidal flow detected (Phys Rev Fluids 2022 analog)
- Nematic order at edges (literature consistent)
- Coffee-ring index nonzero (cells slightly accumulate at edge)

---

## Stage 1e — Radial vs Full Comparison Framework
**Goal**: Two parallel simulations (Sim A = full 3D anisotropic, Sim B = derived radial reduced) with comparison metrics.

### Tasks
- [ ] Derive radial-reduced model via continuum mechanics (see `05_radial_approximation.md`)
- [ ] Implement Sim B as 1D radial ODE integrator (separate code path, same parameter library)
- [ ] Comparison layer: same initial conditions, same parameters, run both
- [ ] Divergence metrics: A/A₀(t) RMS deviation, circularity time series, stress anisotropy Δσ
- [ ] Phase diagram: parameter regions where radial valid vs invalid

### Validation
- For "ideal" parameters (high γ_cc, low σ_a, large R), Sim A ≈ Sim B (radial valid)
- For "active" parameters (high σ_a, anisotropic protrusions), Sim A diverges from Sim B (radial breaks down)
- Phase diagram boundary is smooth and physically interpretable

---

## Stage 1 Final — 5-Point Composition Sweep
**Goal**: Production runs across 5 mechanical-parameter-space points.

### Tasks
- [ ] Define 5 sweep points: 2 anchors (ULA-like, pV4D4-like) + 3 interpolated
- [ ] Batch runner with config-based parameterization
- [ ] 3 stochastic repeats per condition = 15 total runs
- [ ] Sequential GPU scheduling (A5000 first; spillover to RTX 4090×2 / Colab if backlog)
- [ ] Automated dashboard generation per run (matplotlib + plotly)
- [ ] Aggregate analysis across runs

### Validation
- All 15 runs complete with valid output
- Aggregate trends across the 5-point sweep are smooth (no pathological jumps)
- ULA-like and pV4D4-like predictions overlay onto PI's experimental Bare/Pre/Lam4 data within reasonable bounds (note: this is comparison, NOT fit validation — see `00_project_vision.md`)

---

## Stage 2 (Future, After Stage 1 Closeout)
- [ ] Layer 6: chemistry, O₂/nutrient diffusion, viability, necrotic core dynamics
- [ ] Lam-supplemented substrate condition (active in medium)
- [ ] Cell heterogeneity (size/stiffness distributions)
- [ ] Initial spheroid R > 250 μm regime (necrosis-relevant)

## Stage 3 (Future)
- [ ] Other epithelial cell lines (MDCK, etc.)
- [ ] ECM remodeling (collagen alignment, MMP degradation)
- [ ] Substrate stiffness sweep (PA gel comparison)

## Stage 4 (Future)
- [ ] Mesenchymal/EMT cell lines (MDA-MB-231) — model boundary test
- [ ] Chemical signaling (TGF-β coupling)
- [ ] Optional: nucleus mechanics

---

## Current Status (refreshed 2026-04-29, Option F Week 1)

**Active option**: Option F (engineering-Marangoni hybrid) per
`docs/codex_review_synthesis.md`.

### Stages completed
- **Stage 0** — completed 2026-04-29 (env, CUDA, configs, smoke tests).
- **Stage 1a** — completed v15 density-based volumetric stress
  σ_vol = K·(ρ_ref/ρ_kernel − 1)·I; commit 8e06d8e. Architectural R
  drift ceiling 0.244 documented in `docs/outcomes_v15.md`.
- **Stage 1a+** — completed Option β substrate CSF
  γ_sub_eff = γ_sub_star · ecm_strength; commit chain through Phase
  4-v2. R drift improvement vs v15 baseline documented in
  `docs/outcomes_stage1a_plus.md`.
- **Stage 1a++** — completed Layer 2 active boundary stress
  σ_act = -ζ·K·I; ζ=0.6 first PASS of R drift gate;
  `docs/outcomes_stage1a_plus_plus.md`.
- **Stage 1b** — completed Layer 3 φ-ODE (Cho 2020 k_+/k_-);
  `docs/outcomes_stage1b.md`.
- **Stage 1c** — completed Layer 5 Tier 2 mechano-osmotic
  (Guo 2017); `docs/outcomes_stage1c.md`.
- **Stage 1d** — completed Layer 4 Marangoni γ(φ);
  `docs/outcomes_stage1d.md`. Phase 4-v2 phenotype ordering reproduced
  (Bare<Pre<Lam4 in both R drift and A/A₀_topdown).
- **Stage 1e** — Sim A 1D radial reduction module shipped
  (`acs/analysis/radial_reduction.py`); `docs/outcomes_stage1e.md`.
  Sim A vs Sim B comparison **gated by F2 horizontal momentum drift
  HARD-BLOCKER** until anisotropy diagnosis is settled.
- **Stage 2** — Layer 6 chemistry/ECM degradation completed;
  `docs/outcomes_stage2.md`.
- **Production Lam4** — completed (commit b506b57); 5k × 80 hr;
  finding: **peak-and-decay**, top-down peak 1.570 at 7.75 hr, end
  1.409 at 80 hr; Bucket P3 strict (mechanism-missing). See
  `docs/production_lam4_finding.md`.

### Active work — Option F (Codex review + Marangoni)

**Week 1 (current)**: documentation pass.
- [x] `docs/marangoni_review.md` (commit d0a7d99)
- [x] `docs/codex_review_synthesis.md` (commit 6aa3b27)
- [x] `docs/gate_fail_taxonomy.md` (this commit)
- [x] `docs/parameter_registry.md` (this commit)
- [x] roadmap refresh (this commit)
- [ ] `docs/SESSION_HANDOFF.md` consolidation
- [ ] `docs/production_lam4_finding.md` labeling fix

**Week 2 (DONE this commit)**: solver investigation (no code changes).
- [x] `docs/horizontal_momentum_drift_investigation.md` — F2 root-cause
  is gate-normalization-too-tight in overdamped equilibrium + initial-
  pack asymmetry (~1.4% N⁻¹/² floor); absolute drift bounded ~1e-3 across
  all runs. **Reclassified ACCEPTED-LIMITATION** with Stage 1e
  seed-averaging requirement.
- [x] `docs/anchor_force_balance_investigation.md` — F3 root-cause is
  Adami-Hu-Adams §3 kernel truncation at substrate boundary making
  F_pressure_down formula's compressive assumption invalid (P becomes
  tensile, formula reports negative). NOT causal for asymptote
  (F_substrate is under-deflected). **Reclassified ACCEPTED-LIMITATION**.
  F4 same root cause, linked.

**Week 3+ (deferred until Week 2 outputs)**: PI re-decision between:
- Stage 1a++.b (discrete boundary events: lamellipodia / filopodia /
  leader cells) — Option β in marangoni_review.md
- Stage 1d.b (Marangoni Mechanism A / E / F) — Option α
- Or paper-as-is (Option G in codex_review_synthesis.md)

### Deferred items (Codex review items)
- Layer 3 audit (Codex item 4) — blocks Mechanism A/E/F upgrades
- Layer 5 audit (Codex item 5) — blocks Mechanism F coupling
- Test discipline expansion (Codex item 8) — 7 missing tests
- `mlsmpm.py` file split (Codex item 9) — prerequisite for new state vars

### Hard-blocker gate FAILs (publication-claim restricted)
Per `docs/gate_fail_taxonomy.md` (post-Week 2):
- ~~F2 horizontal momentum drift~~ — ACCEPTED-LIMITATION (Stage 1e
  seed-averaging required)
- ~~F3 anchor force balance~~ — ACCEPTED-LIMITATION (NOT causal,
  diagnostic formula bug only)
- ~~F4 contact-band ρ_kernel~~ — ACCEPTED-LIMITATION (linked to F3)
- F9 φ trajectory — REMAINING HARD-BLOCKER (Layer 3 audit pending)

When you (Claude Code) update this file after completing a stage,
mark the stage's checkboxes complete and add a "Stage X completed
YYYY-MM-DD" line.
