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

# Phase C — session handoff + next-session boot (2026-06-22)

**This is the boot artifact for a FRESH session.** Read top-to-bottom, then start at §NEXT.
The previous session was very long; state lives here + in git + Notion Dev Logs.

## What this session did (all committed on `h7/compartment-platform`)

### 1. Full fine-grained physics backlog on the Warp DCM engine — DONE, GPU-verified
All ported off the lumped/proxy mechanisms onto GPU-resident Warp kernels (CLAUDE.md mechanistic rule):
- **A1/A2/A3 mesh integrity** (`ffbfe07`,`8a67d51`,`e875238`): host SWAP/SPLIT/COLLAPSE remesh (frozen
  node pool) · edge-edge contact (seg-seg, catches node-face pokethrough) · dt-CFL adaptive substep + diagnostic.
- **E1 cadherin catch-bond ensemble** (`aaf2101`): explicit trans-dimer bonds (Rakshit catch-slip),
  REPLACES the cohesion tent + M3 binary junction switch. De-cohesion = emergent bond rupture. ×40
  scale bridge: rest/capture at mesh scale, f0=29.2pN molecular.
- **C6 Pereverzev catch-slip ECM clutch** (`7d08d6f`): explicit basal-node↔substrate integrin clutch,
  REPLACES the wetting body-force proxy. k_on=1/s (KB-verified Kim2012/ChanOdde/Bangasser).
- **E2 nucleus** (`8fe5819`): deformable chromatin/lamin bilinear core (H.9 KU-3.B2).
- **B4 surface tension + global area** (`6391d21`); **B5 thin-plate bending** (`4e5202d`).
- **C7 division/proliferation** (`2cbd4cd`): parked cell pool → rim daughters. **C8 necrosis 3-zone**
  (`02b2d5c`): depth O2-proxy + division gating.
- **D10 FCC/Voronoi isotropic builder** (`13b0b11`); **D11 validation gates** (`572d12e`).
- hybrid: cadherin bond cadence 8→50 (`1ddb073`) — host sync was ~6%, not the bottleneck.

Driver: `aleph/warp_port/dcm_warp_decohesion.py` (all flags). Hosts:
`dcm_{cadherin,ecm_clutch,division,necrosis,lamellipodium,junction_switch}_host.py`. Kernels added to
`dcm_neighbor_warp.py` + `dcm_substrate_warp.py`. Viz: `scripts/warp_decoh_mesh_viz.py` (contact/junction
modes, per-frame topology, pen annotation).

### 2. Implicit/IMEX acceleration — coded + validated (CPU/standalone), GPU verify PENDING
`aleph/warp_port/dcm_warp_implicit.py` + driver `--integrator implicit`. Plan: `PHASE_C_IMPLICIT_IMEX_ACCEL_PLAN_2026-06-22.md`.
- I1/I2/I3 linearly-implicit overdamped + matrix-free CG (`7b81dd2`,`3d84c8d`,`69c4c9f`): stable @100×
  CFL where explicit blows up; accuracy ceiling ~100× dt (V/V0 conserved).
- I4 Newton (`f4a66c8`): ceiling 100×→1000× dt.
- all-device Warp CG (`bd6c26a`): == scipy to 3.6e-10, host-transfer-free.
- driver wiring (`19243b0`): IMEX (stiff contact/turgor/edges/well/nucleus/bending implicit w/ frozen-grid
  JVP; soft wetting/lamel/cadherin/clutch explicit in RHS). CPU: same equilibrium @10× dt, 10× fewer steps.
- contact/preconditioner investigation (`2b12a83`): plain CG is correct (3 iters near equilibrium); contact
  does NOT dominate conditioning (γ/dt regularises); a Hutchinson-Jacobi preconditioner HURT → removed.
- **Realised production speedup so far = 0** (no production has used `--integrator implicit` yet).

### 3. Code review — all 10 findings fixed (`44f714c`)
#1 implicit dormant-node freeze · #2 p_div 0.5→0.04 (ref-anchored) · #3 viz cof>=0 mask · #4 implicit
requires grid · #5 lamellipodium live-rim (division daughters) · #6 cadherin koff table 90→300pN · #7 CFL
retune n_sub double-count · #8 necrosis default PROLIFERATING · #9 ConvexHull QJ · #10 turgor_necrotic
doc + lamellipodium upload doubling.

## PRODUCTION RESULT (full fine-grained N=100, explicit) — the headline finding
`outputs/warp_decohesion/figs/n100_fg_{contact,junction}_{montage.png,surface.mp4}` + `n100_fg_viz.npz`.
- **A/A0 = 1.001 (FLAT)** across all 80k spread steps. maxZ 80µm held, V/V0=1.000.
- `[GATES] finite=PASS  G2_interpenetration=FAIL (pen peak 2.49)  volume_conservation=PASS`.
- Visual: compact fcc spheroid that does NOT spread once the wetting PROXY is removed.

## ⚠️ OPEN QUESTIONS (PI decision — NOT concluded; do not auto-fix)
1. **Why A/A0 flat?** Two non-exclusive reads, both = "the deleted SettlingForce/wetting PROXY was
   doing the spreading": (i) the fully-mechanistic stack's equilibrium is ~no-spread; (ii) real-time
   faithfulness — at physiological γ + v_front, 1.6 s of sim is far too short for minutes-scale
   spreading (lamellipodium advances ~0.16 µm in the window). **Decisive test = run LONGER physical
   time via the implicit accelerator** and see if A/A0 climbs (ii) or stays flat (i). Do NOT add a
   spreading driver to force A/A0 up — that is the magic-number/proxy sin PI explicitly rejected.
2. **pen=2.49 (G2 FAIL):** internal cell interpenetration the A-stack (remesh+edge-edge+cfl) does not
   drive to 0. Surface looks fine; overlap is internal. Needs investigation (contact stiffness / dt /
   active-set), not judged.
3. **Validation framing:** PI's A/A0=a+b/R+c/R² is likely a quasi-static EQUILIBRIUM observable → if so,
   timescale-matching isn't required, only reaching mechanical equilibrium (which the accelerator enables).
   Raw PI data = `~/ActiveCellSim/data/experimental/260313_{Bare,Pre,Lam4}.csv` (overlay-only, NO fitting).

## NEXT (start here)
1. `git -C /Users/sw1/ffn_cellsim-platform log --oneline -12` (confirm `44f714c` HEAD).
2. `conda activate ffn_sim`; gbook A5000 is the production GPU (rsync `ffn_sim`→`gbook:~/ffn_phase_c/`,
   `PYTHONPATH=. python -m aleph.warp_port.dcm_warp_decohesion --device cuda:0 ...`).
3. **GPU-verify the implicit driver path** on gbook (`--integrator implicit --accel-dt 8e-4`), then
   **implicit-vs-explicit wall-clock benchmark** (same physical time, real seconds) → the actual speedup.
4. **Decisive (i)/(ii) test:** with implicit, run a much LONGER physical-time N=100 spread → does A/A0
   climb? Co-track pen/V/V0/maxZ (never A/A0 alone — landmine rule).
5. Then: pen=2.49 investigation; optional I5 (implicit-Langevin, kT>0); fold Newton into device_cg for >100× dt.
6. `make kb-check` failed this session on a `python` PATH quirk (no KB changes made → no drift); fix the
   Makefile python or run the tag_kb checks with the conda python.
