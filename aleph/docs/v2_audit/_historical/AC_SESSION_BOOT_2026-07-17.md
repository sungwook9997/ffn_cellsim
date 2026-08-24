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

# SUPERSEDED — Active Cell engine handoff (2026-07-17; do not boot from this file)

> **Invalidated 2026-07-18.** This is a historical handoff only. Its `STABLE` verdict was based on a candidate
> residual (`2957 -> 1992`) that had not satisfied the registered displacement convergence criterion. The
> transactional runtime now correctly rejects that candidate, restores the authoritative state, and reports
> `STABLE=false`: 60 inner iterations, `9.946e-4 µm` candidate displacement versus `7.449e-9 µm` tolerance,
> committed residual restored to `2634.119 pN`, and committed time `0 s`. Boot from `AGENTS.md` and
> [`outputs/ac/foundation-hardening/REPORT.md`](../../outputs/ac/foundation-hardening/REPORT.md), not this file.

This document captures the state as it was understood on 2026-07-17. The historical branch was
`ac/new-engine` at `18b3da8`; its runtime verdicts are not current evidence.

## 0. Boot (minute 0)
1. Read `AGENTS.md`/`CLAUDE.md`, `NEW_ENGINE_BUILD_PLAN_2026-07-16.md` (the master I-spine), and this doc.
   `git log --oneline -15`; confirm branch `ac/new-engine`.
2. Runtime = Warp-CUDA only (I0-A). The **dev Mac cannot run CUDA** — author kernel SOURCE + run pure-numpy
   oracles locally; **all device runs are on gbook** (§2).
3. Confirm the I0-A contract is intact (Warp-only, no HOOMD, head-resolved NMII, 70,686 cortex, MCF7×collagen×α2β1).

## 1. What is DONE (this session)
- **I0-A ratified + I0-B closed per-increment** (fluid I0-B1 α=1/c_v/Π₀; others GAP'd — §4).
- **6 increments analytic-integrated** on `ac/new-engine`: fluid I1a/b/c, nucleus I2, solid I2b, motor I3
  (power-stroke fixed), weave I4 (dendritic Arp2/3, Gate-1 bit-identical), emergence I5-detector.
- **ALL 5 tracks' Warp kernels device-parity PROVEN** on the A5000 (NG-0, round-off vs numpy reference).
- **First runnable ASSEMBLED cell (P5)** — `ac/cell/{assemble,driver}.py` composes bending+link-spring+
  MyosinForce+StericForce+PressureCoupling under one outer-physical/inner-mechanical scheduler (lumped
  ff/ mechanisms retired by non-use). It was **incorrectly reported STABLE** at full native population;
  the 2026-07-18 transactional rerun invalidated that verdict (see the supersession notice above).
- ⭐ **NG-4 fluid-first payload PASSED** — `ac/cell/{fsi_coupling,ng4_payload,fsi_demo}.py`: div(v_s)→fluid
  coupling live (dp 427 Pa vs stub-control 6e-4 Pa); Skempton undrained/drained **2.02×** (balloon=1×);
  τ_p ∝ R²/c_v verified (ratio 2.00); fluid flows v_f 0.47 µm/s. **The cytosol is a genuine rate-dependent
  poroelastic medium, NOT a balloon** (FF's "개살구" replaced, demonstrated).
- **Active myosin activated** (bound fraction 0.14 emerges, 33 pN contractile) — honestly **density-floored**
  (~4 orders under physiological, §6.2, NOT closed by density).
- ⭐ **Assembly MILESTONE-2 landed** (`ac/cell/compartments.py` + assemble/driver wiring) — the composed cell
  now carries the **deformable-mesh nucleus** (`NucleusCompartment`: Helfrich bending + framework-#6 lamina
  3-regime areal tension + device-resident nucleoplasm volume ν→½ + LINC tether = the single nucleus↔cortex
  load path) and the **plasma-membrane Helfrich sheet** (`MembraneCompartment`: bending + γ_mem area tension,
  K_A reservoir OFF, + ERM tether riding on the cortex). `nucleus_shell_kernel` (radial bead ball) **RETIRED
  same-commit** by non-use. The nucleus's live oblate mesh is the fluid domain's **relative-no-flux inner
  inclusion** (`MeshNucleusMaskProvider`, refreshed each remap via `bind_live`), NOT a static sphere. The **I4
  walk_dir hand-off is wired** (`step_myosin_kinetics`: attach→`fill_walk_dir_kernel`→step/detach, so bound
  heads walk the REAL actin barbed-end polarity, not the minifilament-axis default). The **vacuous
  `NucleusNoFluxBC` probe is FIXED** (`p−p≡0` → `suppressed_flux` reads the real nucleus-neighbour cell:
  0 on a uniform field, >0 under a gradient = the mask is load-bearing, not vacuously satisfied).
  Historical native run: **the full 70,686 `STABLE` label is invalid**. Its residual decrease described an
  uncommitted candidate; the accepted state was never proven converged. Device NG-5 core classified
  **4,139 NUCLEUS cells**; `suppressed_flux` changed 0→**2.90**
  under an imposed gradient. 22 new host-oracle tests + full `tests/ac` **297 pass**. Adversarial 6-lens verify
  workflow: **0 confirmed runtime bugs** (2 findings triaged to hardening — provider live-refresh + whole-cell
  force-finiteness gate — both applied).
- **Visuals** (browser-verified): `outputs/ac/cell_t0/` (t0 geometry), `outputs/ac/cell_assembled/`
  (assembled + composed force), `outputs/ac/cell_dynamics/` (**flowing pore-pressure gradient + inward v_f** —
  the professor's target figure).

## 2. The gbook run loop (the ONLY device path)
- ssh alias **`gbook`** (Tailscale). Python `~/miniconda3/envs/ffn_sim/bin/python`. RTX A5000 `cuda:0`, Warp 1.14.
- Code lives at **`~/ffn_ac_native/`** (a SEPARATE dir, NOT the Syncthing folder — Syncthing↔git race, see
  `reference-syncthing-git-race`). Sync before a run: `rsync -az -e "ssh -o BatchMode=yes" aleph/ac/<pkg>/
  gbook:ffn_ac_native/aleph/ac/<pkg>/`. Run: `ssh gbook 'cd ~/ffn_ac_native && PYTHONPATH=. <py> -m
  aleph.components.incumbent.driver --from-resting --n-filaments 70686 ...'`.
- Native-gate runners already written: `ac/fluid/native_gates/ng0_parity.py`,
  `ac/{solid,nucleus,motor}/native_gates/ng0_*_parity.py`, and `ac/cell/{driver,ng4_payload,fsi_demo,
  myosin_activate,dump_coupled}.py`. Track worktrees `../ffn_ac-{fluid-spine,nucleus,solid-ev,motor-nmii,
  emergence,weave}` exist (all merged; reuse or `git worktree remove`).

## 3. NEXT work (priority order — continue the master I-spine)
0. ✅ **Assembly milestone-2 — DONE (2026-07-17).** Deformable-mesh nucleus + membrane wired into
   `ac/cell/{compartments,assemble,driver}.py`, `nucleus_shell_kernel` retired same-commit, I4 `walk_dir`
   hand-off wired (`step_myosin_kinetics`), vacuous no-flux probe fixed (`suppressed_flux`). The former full
   70,686 `STABLE` claim is invalidated; NG-5 nucleus-mask/no-flux teeth were device-verified; 297 tests +
   6-lens adversarial verify (0 confirmed bugs). See §1. **Remaining sub-items now open:** (a) live-mesh
   MEMBRANE fluid provider (nucleus one is live; membrane outer BC still the static sphere = exact at rest);
   (b) the nucleus/membrane **3-D HTML viz** (dump_state emits combined pos but not mesh faces — extend +
   render + browser-verify); (c) interior chromatin WLC net + nucleoplasm viscosity (deferred, named).
1. **Native gates NG-1..NG-9** on the assembled cell — NG-5 nucleus-no-flux core DONE (device: mask + teeth);
   still to run: **NG-1 FSI-OFF regression** (`--no-pressure`/`biot_fsi=None` bit-identical to Warp-FF),
   **NG-2 device conservation** (content balance == membrane flux, no nucleus term — pairs with the NG-5
   suppressed-flux teeth), **NG-6 zero-roundtrip** profiler gate; magnitude NG-3/4/7/8 need I0-B.
3. **Remaining increments** (master order): **I5-proof** (emergence — the falsifiable native run: does the
   isotropic seed CONDENSE bundles? flat → FINDING + SEEDED fallback, never re-tune), **I6** adhesion
   (needs I4 + ECM + adherent), **I7** nucleus-load module (cap→LINC→flatten + cell-type IF), **I8** actin-MT,
   **I9** full production.
   - ✅ **Milestone-2 viz DONE** (`b87d086`): `dump_state` emits nucleus/membrane mesh faces + LINC pairs;
     `ac_cell_assembled_viz` draws the composed NUCLEUS (amber) + MEMBRANE (cyan) + LINC (violet) on the real
     511k cell. Browser-verified (whole + cut-away, 0 JS errors); cut-away reveals the amber nucleus + violet
     LINC inside. `outputs/ac/cell_assembled/ac_cell_assembled.html`.
4. **Adherent cell** — currently SUSPENDED-only; the interesting physics (SF, cap, flatten, migration, real
   traction) needs I6 + ECM + FA. 
5. **σ_EV / cortex excluded volume** — build the cortex WITH excluded volume (pre-relax or sparser placement)
   so the load-bearing steric verdict is valid (currently force_cap'd, ~63k interpenetrating nodes).

## 4. OPEN PI decisions — I0-B GAPs (do NOT choose; ask PI, report-not-tune)
- **σ_EV** discretization (node-node-coarse vs **segment-segment-physical** — Lead recommends segment-segment
  per the fine-grained principle; the assembly demonstrated the node-node interpenetration at full population).
- **φ** (porosity — NOT the 0.7 water fraction; Darcy pore fraction), **c₀** (MCF7 free G-actin monomer).
- **lamin-A/C vs -B split, k_linc** (8 pN is a TENSION not a stiffness — HALT), **rupture threshold**,
  **η_nucleoplasm**, **oblate aspect** (I2/I7); knee_strain + E_nuc should be NULLED to GAP (B minor fix).
- **F_stall/N_side/v0/k_xb** (I3 motor), **Arp2/3 branch/capping rate + NPF density + fascin params** (I4).
- **membrane L_p** (near-seals at τ_p scale — a finding; whole-cell drainage ≫ internal poroelastic time).
- **D_c** draft-provisional (2-6 µm²/s); magnitude verdicts on NG-3/4/7/8 HELD until these close.

## 5. KNOWN STRUCTURAL limits (master §6 — report-not-tune, will NOT just go away)
- **magnitudes density-floored**: γ ~530× floor; single SF O(0.1-0.9) vs 5-6 nN; cap flatten ~0.2 vs ~2.5 µm;
  active myosin 33 pN ~4 orders under physiological. FINDINGS to PI; NEVER add heads/density to close a floor.
- **emergence may not condense** in the overdamped quasi-static solve (I5 falsifiable; the myosin-turnover→
  local-density SF-condensation trigger is still absent; flat → SEEDED-labeled fallback).
- **N-FIXED** (dormant-daughter activation; true growing-N nucleation/severing + the reaction network +
  thymosin = increment-2 deferrals).

## 6. Discipline (carried from CLAUDE.md, enforced this session)
Analytic-first (numpy oracle → Warp SOURCE → gbook gate); tracks never edit `ff/` (INTEGRATION.md patch-notes,
lead applies by-non-use); no I0-B magnitude tuned to an outcome; **validate at FULL native population** (the
coarse trap was demonstrated twice — 2,000-filament runs mislead); browser-verify every HTML viz; each
compartment at its physiological setpoint. Adversarially verify track deliverables before integration.

---

## Retired boot prompt

Do not paste the former prompt: it propagated the invalid `STABLE` claim and stale branch state. Follow the
repository `AGENTS.md`, current Git head, and the foundation-hardening report linked in the notice instead.

## Change log
- 2026-07-17: created at the session-end checkpoint. Fluid-first cytosol demonstrated (NG-4) + visualized;
  assembly (P5) runs at full population; next phase = nucleus/membrane wiring + native gates + I5-I9.
