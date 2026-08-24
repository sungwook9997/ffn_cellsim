---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# Session handoff — ac/ resting baseline (2026-07-22 late) → next Claude session

Branch `codex/ff-ac-codex`, HEAD **`d81d368b`**, gbook (`~/ffn_ac_native`) synced to it. **You are the Lead**
(PI directive). Codex is reduced to a separate branch / debug-only; `codex/ff-ac-codex` is your integration line.

## The one-line state

The multi-session resting-convergence blocker is now fully diagnosed: it is **NOT a solver/preconditioner
problem** (the cortex/fibers are force-balanced at ~0.03 pN). It is that the **plasma membrane carries the
full turgor pressure (~40 pN/node) unheld**, because the **membrane↔cortex ERM load path does not transmit
force**. The immediate next step is a **geometry check**, not more code.

## What this session delivered (all committed)

- **FIRE static solver** (`2fac388f` kernels+CPU-validation, `83fde9b5` driver 3-way wiring, `c9732b07` CLI
  `--inner-solver fire`). Tangent-free minimizer of the SAME projected descent; no global tangent / linear
  solve / preconditioner. CPU-validated (Warp==NumPy <1e-9; converges 2e-13 vs plain descent 0.33 at cond~300).
  On gbook, plain explicit descent DIVERGES on the resting cell (47.4→47.8) while FIRE is stable. **Implicit is
  retired for statics.** Dynamic runs already default to explicit (bleb_growth).
- **Resting force-balance preload infra** (`14bb1c96`): `_preload_erm_resting_balance` (pre-stretch ERM to
  cancel the membrane node force, incremental fixed point `rest += f_in/k_erm`) + `_preload_cortex_pretension`
  (scale crosslink rest by `1-ε` for Laplace hoop tension). Both opt-in, default off. CLI `--preload-erm-balance`,
  `run_from_resting(..., cortex_prestrain=)`.
- **Pressure provenance correction** (`9c0f82c0`, params_i0b1.yaml): 40 Pa is a **HeLa proxy**
  (Fischer-Friedrich 2014), NOT MCF7 — **MCF7 value HELD**. It is `dP_hyd` (net hydrostatic ΔP the cortex
  bears), NOT osmotic. Keep three states: `dP_hyd` / `dPi_osm` / `p_excess`; water flux = `L_p(dPi_osm−dP_hyd)`.
  MCF7 dP set by a **PI-authored consistency gate** `dP=2(γ_cortex+γ_mem)/R` (MCF7 γ~0.41 mN/m → ~109 Pa,
  back-computed), NOT hand-swapped. **Still TODO in code**: rename the `PI_0_PA` constant + split
  `osmotic_difference=PI_0_PA` (driver.py ~976) into the two distinct states.
- **Radial ERM pairing** (`d81d368b`): `_pair_radial` + `CellConfig.erm_radial_pairing` (opt-in). Pairs each
  membrane node to its most-radial cortex node. **Did NOT fix the coupling** (see below).

## The diagnosis chain (evidence, gbook 8k & 20k filaments — identical, so it is the membrane not the cortex)

Per-compartment max|F| at the built resting geometry (`_accumulate_all`, one force eval):
- actin/cortex: **0.03 pN** (balanced) · myosin 0 · nucleus 0.16 · **membrane: max 47.4, mean 40.7 pN** (uniform
  outward = turgor pressure ΔP·A_node). Count-independent ⇒ it is the 642-node membrane, not the fibers.
- Physics: a curved membrane under turgor ΔP needs tension `γ=ΔP·R/2` (Young-Laplace) to balance it; `γ_mem`
  (lipid) is ~10 pN/µm, far too small — the **cortex** must bear it (cortical tension). A relaxed cortex +
  positive turgor has NO force-balanced equilibrium, so implicit-deflation AND FIRE both stall at the same 47.4.
- Attempted fixes and what they proved:
  - Cortex pre-tension: crosslink pre-strain explodes cortex force (0.03→**10,600** pN at ε=0.02) while the
    **membrane stays exactly 47.4/40.7** → cortex tension does not reach the membrane.
  - ERM preload: transfers only **~11 of 40 pN** to the cortex; membrane unchanged.
  - Radial pairing: ~11 pN still, membrane unchanged → pairing DIRECTION is not the issue.
  ⇒ The membrane and cortex are **mechanically decoupled**; the ERM linkage does not carry the load either way.

## ⭐ DO THIS FIRST (next session): verify the membrane↔cortex tether geometry

**Leading hypothesis:** the overlap-free cortex radial thickness (0.2 µm) spreads cortex nodes across R≈7.30–7.50,
reaching `R_mem = R_CELL_UM = 7.5` (cortex centre `R_CORTEX_UM = 7.40`). So membrane↔cortex ERM tethers are
**near-degenerate (length ~0, ill-defined direction)** — the preload math (`u = d/|d|`) and the force transmission
both break when |d|→0. Check before writing any fix:

1. Build a cell (8k, `overlap_free_cortex=True`), read `cell.membrane.erm_m_d/erm_c_d`, compute the tether
   length distribution `|pos[erm_c]-pos[erm_m]|` and the radial-vs-tangential split. Are lengths ~0? Is the
   direction radial?
2. Inspect the geometry: R_mem 7.5 vs cortex node radii (min/mean/max) with overlap-free thickness. Does the
   cortex top (7.50) collide with the membrane (7.50)?
3. If degenerate: the membrane must sit CLEANLY outside the cortex with a finite, consistent ERM gap (e.g. put
   the cortex shell fully inside R_mem, or set the ERM to a finite physiological rest length). Then re-test the
   ERM preload — it should transmit the full ~40 pN, and cortex pre-tension (calibrated to MCF7 γ, so
   `dP_hyd=2γ/R` is self-consistent) should hold the membrane → **FIRE converges the resting baseline**.

Only after force balance is reached (strict projected-force gate) is the baseline a valid production state.

## Then (already-built parts, just assembly)

full-native FIRE validation @70,686 → KERNEL_BOUND the ac/engine facades to native kernels (no dup) → compose
MT/IF/ECM into build_cell → bleb growth at membrane subdiv 6 → dynamic runs. FIRE, compartments, engine
scaffold, viz are all ready; the resting baseline is the ONLY closed gate.

## Operational notes
- gbook: `ssh gbook`, python `~/miniconda3/envs/ffn_sim/bin/python`, repo `~/ffn_ac_native`. **Syncthing: this
  Mac repo IS the source** — `rsync -az <file> gbook:~/ffn_ac_native/<file>` forces an immediate push (Syncthing
  won't revert; it converges to this repo). Wait on a run with `ssh gbook 'while kill -0 <pid> 2>/dev/null; do
  sleep 5; done; cat log'`. **Dev Mac is CPU-only Warp; `build_cell` is CUDA-gated (I0-A)** → full-cell runs are
  gbook-only. CPU here = kernel codegen + NumPy-reference gates only.
- Shared branch: always `git diff --stat` before commit (Codex bundling hazard). Commits end with
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. `ffn/foundation` push is PI-gated.
- Open PI gaps: **MCF7 ERM linker density** (KB has none; currently 1/node diagnostic) · **MCF7 γ_cortex source**
  for the pressure-consistency gate · ECM collagen gate HELD (30–100 vs 5–100 Pa).
- Respond in Korean (code/commit/docs in convention). Validate at FULL NATIVE. Codex = debug idea only.

## Key files
- `ac/cell/inner_mechanics.py` (FIRE kernels), `ac/cell/driver.py` (FIRE branch + `_preload_erm_resting_balance`
  + `_preload_cortex_pretension` + run_from_resting), `ac/cell/compartments.py` (`_pair_radial`,
  `build_membrane_compartment` ERM pairing, MembraneCompartment.accumulate), `ac/cell/assemble.py`
  (build_cell, R_CELL_UM/R_CORTEX_UM, CellConfig), `ff/cortex_assembly.py` (overlap-free thickness),
  `ac/cell/membrane_pressure.py` (turgor traction = (p_inside−p_ext)·A·n_out), `ac/fluid/params_i0b1.yaml`
  (dP_hyd provenance).
