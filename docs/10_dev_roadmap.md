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
**Goal**: A single 3D spheroid sits in vacuum, holds together, reaches mechanical equilibrium.

### Tasks
- [ ] Implement MLS-MPM solver (use Taichi's official MLS-MPM example as starting point)
- [ ] Cell-equivalent material point initialization (random pack inside a sphere, ~5,000 points; pilot uses ~1,000)
- [ ] Cohesive force (surface-tension-like) keeping cells together
- [ ] Cortical tension (effective bulk modulus)
- [ ] Maxwell viscoelastic constitutive law
- [ ] HDF5 frame export
- [ ] Vispy real-time viewer (offline mode: writes PNG sequence)
- [ ] Unit tests: mass conservation, momentum conservation in absence of external force

### Validation
- Spheroid converges to spherical shape within τ_relax ~ minutes (sim time)
- Total energy decays monotonically (Maxwell relaxation)
- Mass conserved exactly (MLS-MPM is conservative)
- Pilot run completes in **< 30 minutes wall-clock**

### Pilot config (Stage 1a target)
- 1,000 material points
- 4 sim hours (after equilibrium; equilibrium phase ~30 sim min)
- Wall-clock target: ≤ 30 min on RTX A5000

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

## Current Status
**As of project initiation**: Stage 0 not yet started. Begin with Stage 0 environment setup, then proceed strictly in Stage order.

When you (Claude Code) update this file after completing a stage, mark the stage's checkboxes complete and add a "Stage X completed YYYY-MM-DD" line.
