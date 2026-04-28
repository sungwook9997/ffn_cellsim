# ActiveCellSim — Claude Code Project Context

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
- **Hardware target**: NVIDIA RTX A5000 (primary), portable to RTX 4090×2, fallback Google Colab
- **GPU portability is mandatory** — no hard-coded device paths
- **Execution**: Headless SSH, full automation (no mid-run user input)
- **Visualization**:
  - Real-time viewer: Vispy (lightweight, SSH-friendly via offline frame export)
  - Production render: Blender Cycles via `blender --background --python` (headless)
  - Analysis dashboard: matplotlib + plotly
- **Data format**: HDF5 for simulation snapshots (frames), JSON for metadata, CSV for metrics

### Performance Budget
- **Pilot**: 1,000 material points × 4–8 sim hours → wall-clock **≤ 30 minutes**
- **Production**: 5,000 material points × 80+ sim hours → wall-clock **5–15 hours/run**
- **Frame interval**: pilot 60 min, production 15 min (matches PI's experimental imaging)
- **Sweep**: 5 conditions × 3 repeats = 15 runs, sequential batch

### Validation Principle (Read Carefully)
- **Literature-first**: All physical parameters from peer-reviewed sources, prefer IF ≥ 15 (Nat Phys, Nat Mater, Nat Cell Biol, Cell, Science, PNAS, Nat Commun, eLife)
- **PI's data role**: Comparison overlay only. **Do NOT use for parameter fitting.** See `data/experimental/` for CSV files (Bare/Pre/Lam4 from 2026-03-13 experiment).
- **Cho et al. 2020 role**: Cross-check only. Use literature-derived ODE rates first; check against Cho's Western blot timecourse second.
- **If literature predictions disagree with experiments**: that is a finding, not a bug. Honesty over fit.

## Where to find what

```
ActiveCellSim/
├── CLAUDE.md                            # This file (always read at start)
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
│   ├── 10_dev_roadmap.md                # Bottom-up Stage 1a → 1d → 2 milestones
│   ├── 11_performance_protocol.md       # Benchmark, profiling, GPU portability
│   ├── 12_validation.md                 # Literature anchors, assumptions, limitations
│   ├── 13_data_schema.md                # HDF5 structure, metadata, naming
│   └── references.bib                   # All citations
├── data/
│   ├── experimental/                    # PI's CSVs (260313_Bare/Lam4/Pre.csv) — READ ONLY for overlay
│   └── literature/                      # Extracted parameter tables from papers
└── (code/, results/, etc. — Claude Code creates these)
```

**Progressive disclosure**: Start with this file + `docs/10_dev_roadmap.md`. Pull other docs only when relevant to current task. Do NOT load all docs at once.

## Code Conventions
- Modular structure: `physics/`, `boundary/`, `adhesion/`, `viz/`, `io/`, `analysis/` separation
- Configuration via YAML (`configs/pilot.yaml`, `configs/production.yaml`)
- Each Tier/Layer is independently togglable via config flag (for staged activation)
- Type hints + docstrings (Google style)
- Unit tests for physics modules (especially conservation laws, scaling tests)
- Logging with severity levels (DEBUG for dev, INFO for production runs)

## Development Order (Bottom-up)
**Strict order — do not skip stages:**
1. **Stage 1a**: Single deformable spheroid relaxation (free-floating, no substrate) → verify equilibrium shape, conservation laws
2. **Stage 1a+**: Add substrate contact (Hertz + adhesion) → verify wetting onset
3. **Stage 1a++**: Activate Layer 2 (boundary biology) → verify lamellipodia events, leader cell emergence
4. **Stage 1b**: Activate Layer 3 (φ adhesion dynamics) → verify ULA vs pV4D4 phenotype distinction
5. **Stage 1c**: Activate Layer 5 Tier 2 (mechano-osmotic) → verify spreading-induced volume loss
6. **Stage 1d**: Activate Layer 4 (Marangoni, nematic, vortex analysis) → verify internal flow patterns
7. **Stage 1e**: Build radial-reduced sim B and comparison framework
8. **Stage 2+**: Reserved for chemistry/necrosis, see `docs/10_dev_roadmap.md`

Run validation tests at each stage end. Don't move forward with broken physics.

## Hard Rules
- **NEVER** fit parameters to the PI's experimental data. Literature-derived only.
- **NEVER** hard-code GPU device IDs or paths. Read from config.
- **ALWAYS** save full configuration + git commit hash with each run output.
- **ALWAYS** treat `data/experimental/*.csv` as read-only.
- If a physical parameter has no IF≥15 reference, flag it explicitly in code comments and `docs/12_validation.md`.

## When in doubt
- Read `docs/00_project_vision.md` for framing
- Read `docs/10_dev_roadmap.md` for what to do next
- Read `docs/12_validation.md` for assumption list and limitations
- Ask the PI before deviating from any principle in this file
