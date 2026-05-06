# ActiveCellSim — Codex Project Context

## WHY (Project Vision)
This is an **independent first-principles computational framework** for simulating MCF7 spheroid spreading on Col1-coated substrate, built from mechanobiology fundamentals — NOT a tool to fit existing experimental data.

The PI's experimental research uses a radial-symmetric phenomenological model (`A/A0 = a + b/R + c/R²`) to analyze spheroid spreading. This simulation is built **independently from that model** to:
1. Generate new mechanobiological insight from a fluid-dynamics-inspired perspective
2. Quantify when/how the radial approximation is valid (radial vs full anisotropic comparison)
3. Test whether the PI's thin-film drying intuitions (Marangoni, coffee-ring, drying-induced concentration) systematically map to cellular spreading

If results converge with the PI's experimental data, that is **bonus validation** — but the simulation must stand on its own academic merit.

## WHAT (Architecture in 1 paragraph)
A 3D hybrid Eulerian-Lagrangian simulation using **Taichi MLS-MPM** with ~5,000 cell-equivalent material points carrying internal state. The model integrates 5 layers: (1) bulk active viscoelastic hydrodynamics, (2) boundary cell biology (lamellipodia/filopodia/focal adhesions), (3) adhesion network dynamics (φ ODE for E-cad ↔ Int-β1 transition), (4) internal flow dynamics (cellular Marangoni, nematic order, vortex), (5) mechano-osmotic coupling (phenomenological Tier 2). Layer 6 (chemistry/necrosis) reserved as Stage 2 hook. Comparison framework runs **two independent simulations** — full 3D anisotropic vs properly-derived radial-reduced ODE — to quantify the validity domain of radial approximations.

## HOW (Critical Technical Constraints)

### Stack
- **Language**: Python 3.11+
- **Simulation core**: Taichi (>= 1.7), prefer MLS-MPM examples as starting point
- **Hardware target**:
  - **Primary working environment**: NVIDIA RTX A5000 *Laptop* GPU, **16 GB VRAM** (Xeon W-11955M mobile workstation, Windows, headless SSH). Most pilots and many production runs land here. **16 GB is the baseline reference for all sizing.**
  - **Secondary / large-sweep environments**: 24 GB+ GPUs — RTX 4090 ×2 (lab workstation), Google Colab A100 (40 GB), or desktop A5000 (24 GB). Used for runs that exceed 16 GB or for parallel sweep throughput.
  - **GPU portability is mandatory** — no hard-coded device IDs or paths; the `ACS_GPU_BACKEND` env var (cuda/vulkan/opengl/metal/cpu/auto) and `acs.gpu.init_taichi` dispatcher route everything.
- **Execution**: Headless SSH, full automation (no mid-run user input)
- **Visualization**:
  - Real-time viewer: Vispy (lightweight, SSH-friendly via offline frame export)
  - Production render: Blender Cycles via `blender --background --python` (headless)
  - Analysis dashboard: matplotlib + plotly
- **Data format**: HDF5 for simulation snapshots (frames), JSON for metadata, CSV for metrics

### Performance Budget
Targets are written against the **16 GB Laptop A5000 baseline**. The same configs must also be runnable on 24 GB+ targets (typically faster, never blocked).
- **Pilot**: 1,000 material points × 4–8 sim hours → wall-clock **≤ 30 minutes** on Laptop A5000.
- **Production**: 5,000 material points × 80+ sim hours → wall-clock **~7–20 hours/run** on Laptop A5000 (Laptop GPU runs ~70% of desktop A5000). On 4090×2 / A100, expect closer to the original 5–15 hr.
- **VRAM ceiling**: pilot configs target ≤ 8 GB working set; production targets ≤ 12 GB so OS + frame buffers + Taichi runtime have headroom on the 16 GB device.
- **Stage 1a benchmark (mandatory)**: measure actual peak VRAM and wall-clock of the 5,000-point production config on Laptop A5000. If it exceeds 12 GB, split into a `production_24gb.yaml` configuration and reduce the 16 GB version to whatever fits.
- **Frame interval**: pilot 60 min, production 15 min (matches PI's experimental imaging).
- **Sweep**: 5 conditions × 3 repeats = 15 runs, sequential batch (Laptop A5000 by default; spillover to 4090×2 / Colab when backlog > 1 day).

### Validation Principle (Read Carefully)
- **Literature-first**: All physical parameters from peer-reviewed sources, prefer IF ≥ 15 (Nat Phys, Nat Mater, Nat Cell Biol, Cell, Science, PNAS, Nat Commun, eLife)
- **PI's data role**: Comparison overlay only. **Do NOT use for parameter fitting.** See `data/experimental/` for CSV files (Bare/Pre/Lam4 from 2026-03-13 experiment).
- **Cho et al. 2020 role**: Cross-check only. Use literature-derived ODE rates first; check against Cho's Western blot timecourse second.
- **If literature predictions disagree with experiments**: that is a finding, not a bug. Honesty over fit.

## Where to find what

Authoritative file→version map: `STRUCTURE.md` at repo root.

```
ActiveCellSim/
├── AGENTS.md                            # This file (always read at start)
├── STRUCTURE.md                         # Full repo file→v1/v2/v3/shared map
├── README.md                            # Human-facing project overview
├── docs/
│   ├── 00_project_vision.md             # Full WHY + framing + assumptions list
│   ├── 01_model_architecture.md         # 5-layer architecture detail
│   ├── 02_force_models.md               # Every force term with IF≥15 references
│   ├── 03_adhesion_dynamics.md          # φ ODE, E-cad ↔ Int-β1, Cho cross-check
│   ├── 04_simulation_setup.md           # Initial/boundary conditions, time step
│   ├── 05_radial_approximation.md       # Proper continuum derivation of A/A0 form
│   ├── 06_radial_full_comparison.md     # Sim A vs Sim B comparison protocol
│   ├── 07_internal_flow_dynamics.md     # Marangoni, nematic Q, vortex (drying-inspired)
│   ├── 08_mechano_osmotic.md            # Tier 2 phenomenological volume regulation
│   ├── 09_visualization.md              # Vispy viewer + Blender + dashboard specs
│   ├── 11_performance_protocol.md       # Benchmark, profiling, GPU portability
│   ├── 12_validation.md                 # Literature anchors, assumptions, limitations
│   ├── 13_data_schema.md                # HDF5 structure, metadata, naming
│   ├── references.bib                   # All citations
│   ├── v1/                              # v1 (frozen) docs: 10_dev_roadmap.md, outcomes_*, stage1*, path_c_*, production_lam4_*, parameter_registry, gate_fail_taxonomy, marangoni/layer3 reviews, …
│   └── v2/                              # v2 (active) docs: 00_project_vision_v2.md, 10_dev_roadmap_v2.md, all v2_*.md (briefs, locked specs, sanity gates)
├── data/
│   ├── experimental/                    # PI's CSVs (260313_Bare/Lam4/Pre.csv) — READ ONLY for overlay
│   └── literature/                      # Extracted parameter tables from papers
└── (code/, results/, etc. — Codex creates these)
```

**Progressive disclosure**: Start with this file + `docs/v2/10_dev_roadmap_v2.md` (active). The v1 roadmap is at `docs/v1/10_dev_roadmap.md` for historical reference. Pull other docs only when relevant to current task. Do NOT load all docs at once. Authoritative file→version map: `STRUCTURE.md` at repo root.

## Code Conventions
- Modular structure: `physics/`, `boundary/`, `adhesion/`, `viz/`, `io/`, `analysis/` separation
- Configuration via YAML. **VRAM-tier-segregated production configs** so the 16 GB baseline is never silently broken by a parameter intended for a larger GPU:
  - `configs/dev.yaml` — fast smoke iteration (≤ 2 GB)
  - `configs/pilot.yaml` — 16 GB Laptop A5000, Stage 1a target
  - `configs/production_16gb.yaml` — 16 GB baseline production config (default sweep)
  - `configs/production_24gb.yaml` — 24 GB+ target (RTX 4090, desktop A5000) — only created when a parameter set genuinely cannot fit in 16 GB
  - `configs/production.yaml` is reserved as a symlink/alias to whichever 16gb/24gb file is the current default; do not edit it directly
- Each Tier/Layer is independently togglable via config flag (for staged activation)
- Type hints + docstrings (Google style)
- Unit tests for physics modules (especially conservation laws, scaling tests)
- Logging with severity levels (DEBUG for dev, INFO for production runs)

## Active roadmap (v2 supersedes v1 for new work)

The original `docs/v1/10_dev_roadmap.md` (Stage 1a → 2 below) describes the
v1 spheroid-first continuum prototype, kept frozen for reference and
reproducibility of prior results. **All new development follows
`docs/v2/10_dev_roadmap_v2.md`** — image-constrained, cell-resolved,
single-cell first. See `docs/v2/00_project_vision_v2.md` for why we pivoted
and `docs/v1/v1_continuum_backup.md` for what v1 deliberately stays as.
The `acs/v2/` package is the v2 code home; v1 modules remain under
`acs/` (non-`v2/` subpaths) and are not deleted.

When in doubt about which roadmap to follow: v2.

## Development Order — v1 (frozen, historical reference)
**Strict order — do not skip stages:**
1. **Stage 1a**: Single deformable spheroid relaxation (free-floating, no substrate) → verify equilibrium shape, conservation laws
2. **Stage 1a+**: Add substrate contact (Hertz + adhesion) → verify wetting onset
3. **Stage 1a++**: Activate Layer 2 (boundary biology) → verify lamellipodia events, leader cell emergence
4. **Stage 1b**: Activate Layer 3 (φ adhesion dynamics) → verify ULA vs pV4D4 phenotype distinction
5. **Stage 1c**: Activate Layer 5 Tier 2 (mechano-osmotic) → verify spreading-induced volume loss
6. **Stage 1d**: Activate Layer 4 (Marangoni, nematic, vortex analysis) → verify internal flow patterns
7. **Stage 1e**: Build radial-reduced sim B and comparison framework
8. **Stage 2+**: Reserved for chemistry/necrosis, see `docs/v1/10_dev_roadmap.md`

Run validation tests at each stage end. Don't move forward with broken physics.

## Hard Rules
- **NEVER** fit parameters to the PI's experimental data. Literature-derived only.
- **NEVER** hard-code GPU device IDs or paths. Read from config.
- **ALWAYS** save full configuration + git commit hash with each run output.
- **ALWAYS** treat `data/experimental/*.csv` as read-only.
- If a physical parameter has no IF≥15 reference, flag it explicitly in code comments and `docs/12_validation.md`.
- **ALWAYS** run the Sanity Gate Protocol (below) on any new or modified physics/numerics module BEFORE the first execution.
- **NEVER** introduce an empirical scaling factor / ad-hoc tuning constant ("magic number") that fails any of the three tests in the Magic-Number Block (below). Such factors paper over a missing principled term and rot quickly under parameter sweeps. If a numerical correction is genuinely needed, derive it from literature or first principles (e.g., proper finite-difference curvature operator instead of a scalar κ-scale knob); if you cannot, halt and surface to the PI.
- **NEVER** modify a gate's tolerance, normalisation, or check window to make a failing run pass. Gates are validation contracts written before the run; if a gate is genuinely incorrect, surface to the PI for a contract change rather than editing it inline.
- **Rule 10 (Dimensional comparison verification)**: any order-of-magnitude comparison (e.g. "gravity ≈ surface tension") must be done with both quantities reduced to the SAME unit (force-per-area is the canonical comparison unit at a 2D interface; force-per-volume is the canonical bulk unit). Comparing per-volume to per-area without the characteristic-length conversion is a unit mismatch. The Path C `g_star = 0.1` episode (2026-04-29) is the canonical anti-pattern: the sanity-md proposed `g_star/κ ≈ 0.036` "same order as γ_star=0.01" — a per-volume-vs-per-area mismatch that hid a 7× over-anchoring; the corrected comparison `g_star · H` vs `γ · κ` gave 7×, predicting (and producing) the pancake regime. Before any dimensional balance claim, write the unit chain explicitly.
- **Rule 11 (Sim-experiment measurement matching)**: when comparing simulation outputs to PI experimental data (e.g. A/A₀ trajectories), the simulation measurement modality must match the experimental imaging modality. PI experimental A/A₀ is a top-down microscope projection (`Area_um2` column in `data/experimental/260313_*.csv`), NOT a substrate-contact area. Substrate-contact-area metrics depopulate on lift-off and give artefactual A/A₀ ≈ 0 even when the spheroid is intact. The correct simulation analog of PI's A/A₀ is the xy-plane convex hull of all particles (z-independent top-down projection), implemented in `acs/analysis/shape_metrics.py:top_down_projection_area`. The Stage 1a+ inherited "A/A₀ FAIL" pattern across Stage 1b/1c/1d/2 was a measurement-protocol mismatch caught at the Phase 1.2 outstanding-issue resolution (2026-04-29). Before adding any new comparison gate, document the experimental measurement modality and verify the simulation measures the same quantity.

## Magic-Number Block (mandatory check)

A "magic number" is any new empirical scaling factor, tuning constant, or correction multiplier introduced into a physics / numerics module. Before adding one, the author must answer all three of the following. **Any "yes" on test 3, or "no" on tests 1–2, blocks the change** — surface to the PI with the underlying scheme issue instead.

1. **Derivable**: can the value be derived from a literature reference (preferably IF ≥ 15) or from first principles (dimensional analysis, conservation law, asymptotic expansion)?
2. **Grid-invariant**: does the value remain valid as `dx`, `dt`, `grid_n`, or `n_particles` change? Or does each parameter sweep require re-tuning?
3. **Fitting**: was the value chosen to make a specific simulated number match a target (gate threshold, experimental datapoint, prior result)?

The April 2026 `csf_kappa_scale=0.05` episode is the canonical anti-pattern: it failed test 1 (no derivation), failed test 2 (would need re-calibration at every grid resolution), and was true on test 3 (chosen so the radius drift gate would pass). The principled fix was a proper finite-difference curvature operator (∇·n̂) — which has a literature reference, is grid-invariant up to discretisation error, and was not chosen to fit any target.

## Sanity Gate Protocol (mandatory before first execution of any physics/numerics code)

After writing or modifying a physics/numerics module, **before running it**, perform the following five checks. Record the results either as a docstring section in the module under the heading `Sanity Gate` or as a sibling `*_sanity.md` note when the analysis is too long for a docstring. If any check FAILs, halt and surface the failure to the PI with concrete options before any further code is written. The CFL violation discovered in the first `mlsmpm.py` draft (dt=0.02 s vs CFL limit ≈1 μs, off by 16,000×) is the canonical example of this protocol catching a fork-in-the-road decision early.

### 1. Dimensional analysis
- Every input and output annotated with units (μm, s, Pa, kg, …).
- Compute the characteristic scales: length L, time T, stress S, velocity V, mass M.
- Compute the dominant non-dimensional numbers: Reynolds (Re), capillary (Ca), Deborah (De), Péclet (Pe), Mach (Ma) as relevant.
- For explicit time-stepping schemes: compute the CFL / stability bound and verify `dt < dt_critical` with the configured parameters. If not, FAIL.
- Record the result as `# Dimensional check: PASS — Re=X, Ca=Y, De=Z, dt/dt_CFL=W` or `# Dimensional check: FAIL — <reason>`.

### 2. Boundary cases
- Walk through the parameter / discretisation extremes that the code must survive: N → 0, N → ∞, Δt → 0, Δt → large, γ → 0, K → ∞, R → 0, etc.
- Either prove (with a comment) that the limit is well-defined or add a guard. Silent NaN-on-extremes is a FAIL.

### 3. Conservation invariants
- State explicitly which quantities are conserved (mass, momentum, angular momentum, energy) and on which lines/kernels they are preserved.
- State the deliberate dissipation channels (viscous, drag, Maxwell relaxation) and verify they are *only* dissipative (no sign flip).
- Identify suspect leak points (e.g., reflective boundaries that may absorb momentum, atomic ops with non-deterministic order).

### 4. Numerical sanity
- Δt vs the slowest physically meaningful timescale (Δt ≪ τ_relax, etc.) and the fastest resolved timescale.
- Grid resolution vs the smallest physical feature you need to capture (`dx ≪ R₀`, `dx ≲ thickness of boundary layer`).
- Float precision: justify f32 vs f64 explicitly (typical: f32 for particle/grid hot fields, f64 for cumulative diagnostics).

### 5. Sign / sense check
- For every force / flux term, write one line on the *direction* it pushes (cohesion ⇒ attractive ⇒ negative work on expansion; pressure ⇒ repulsive ⇒ positive work on expansion; drag ⇒ opposes velocity; etc.).
- A force whose sign cannot be checked against intuition is itself a FAIL.

### 6. Measurement-protocol consistency
- An analytical proof of correctness is only meaningful if it covers the *measurement protocol* used to evaluate the gate or report the result. A point/peak analytical derivation (e.g., "exact at the gradient peak", "exact for a uniform field") is **insufficient** when the measurement is band-averaged, integrated over a region, sampled at multiple points, or otherwise evaluated *off* the analytical-proof point.
- Walk through where the gate or report value comes from. Identify every cell, particle, or sample the measurement averages or integrates over. For each, write down whether the analytical proof still holds — and if it does not, either (i) extend the proof to cover those samples, (ii) change the measurement to match the proof's domain, or (iii) record the off-proof contribution as a known systematic and bound it.
- The April 2026 v13 episode (Stage 1a `_build_curvature` Laplacian form) is the canonical anti-pattern: the analytical "κ = 2/R exact at the gradient peak" proof passed checks 1–5, but the band-averaged measurement included off-peak cells where `κ = f″/f′ + 2/r` picked up an asymmetric contribution `|f″/f′|_max ≈ 2/δ_smoothing ≈ 13`, giving a measured κ 5× worse than the prior scheme. The off-peak response was the missing protocol item.

### Failure handling
A FAIL halts further code work for the current module. Surface the issue to the PI with at least three concrete options (e.g., reduce Δt, switch to implicit, switch to overdamped). Wait for direction before proceeding. Never silently work around a Sanity Gate failure.

## Chat pane vs work pane discipline (acs-collab tmux)

The acs-collab tmux layout is **two panes per agent**: `codex-chat`
receives PI messages via the relay, `codex-work` runs long
implementation/review/test cycles. Claude mirrors this with `claude-chat` /
`claude-work`. The split exists so PI never has to wait for an
implementing agent — the chat pane stays free to acknowledge, route,
and report status.

**You can tell which pane you are in two ways**:
- The relay wrapper that delivers each PI message includes
  `pane=chat work_pane=<your_name>-work`. If you see `pane=chat`, you are
  the chat instance.
- Otherwise (e.g. a session directly briefed via `work_briefing.md`),
  inspect `tmux display -p '#S'` — `*-work` means work pane.

**Hard rules for the chat pane** (no exceptions without explicit PI ask):
1. **Acknowledge fast, finish nothing big.** Reply within one turn,
   ideally in 1–3 short sentences plus an MCP `send` and a `/agent_status`
   POST. Never start a multi-minute or multi-file edit yourself.
2. **Dispatch every non-trivial task to the work pane** via
   `/tmp/acs-collab/work_briefing.md`. Once dispatched, post a one-line
   "delegated to `<your_name>-work`" and stay idle for the next PI
   message.
3. **No commits, no `git add`, no destructive ops from chat.** The work
   pane handles those.
4. **No claim transfers without notifying PI.** If chat holds a claim,
   hand it off in the briefing and announce the transfer via MCP.

**Hard rules for the work pane**:
1. **Heartbeat aggressively.** `/agent_status` POST every milestone and
   at minimum every 1–2 minutes during long tasks. Use the pane-aware
   `agent` value (`codex-work`) so the sidebar 4-pane panel stays
   informative.
2. **Milestone reports via MCP `send`.** Short, structured: `[<pane> ·
   NN%] <one line>`. Batched.
3. **Open questions go back to PI through the work pane's MCP send**, not
   by hijacking the chat pane.

## When in doubt
- Read `docs/00_project_vision.md` for framing
- Read `docs/v2/10_dev_roadmap_v2.md` for what to do next (v2 active); `docs/v1/10_dev_roadmap.md` for v1 history
- Read `docs/12_validation.md` for assumption list and limitations
- Ask the PI before deviating from any principle in this file
