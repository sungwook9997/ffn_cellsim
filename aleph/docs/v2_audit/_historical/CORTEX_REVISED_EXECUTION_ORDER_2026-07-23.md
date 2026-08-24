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

# Cortex — revised execution order (PI dynamic-cortex direction folded into the current work), 2026-07-23

PI direction added to the IN-PROGRESS percolation/resting work — **do not stop or revert it; do not treat it as
the final cortex.** This supersedes the framing in `CORTEX_EMERGENT_AND_REMAINING_PLAN_2026-07-23.md` (that plan's
emergent-binding piece is now one sub-item under Stage 3 of the broader dynamic-cortex program below).

## §0. Reframe (the key decision)

- The committed `density_per_fil=20` (`bc5ff3b0`) + `overlap_free_cortex` default (`81805535`) are a **STATIC
  CONNECTED CONTROL** (native-verified: single-spanning, rigid, 0.2µm, 0 interpenetrations). They are a diagnostic
  fix for the fragmentation artifact — **NOT the physiological dynamic cortex.**
- The final cortex = a **CellState-conditioned, finite-monomer, variable-length, variable-N, KMC-topology,
  accepted-step-committed DYNAMIC network.**
- `70,686 / 3µm / 20-xl-per-fil` are the **INITIAL active population at the first specified CellState**, NOT
  permanent runtime invariants. When length/count change, `xl/fil=20` is not conserved; the conserved quantity is
  **state-conditioned total crosslinker inventory / contour-length-based inventory**, with bound topology from KMC.

## §1. Organizing principle — CellState → conditional inventory/rate (do NOT collapse to a phenotype axis)

The cortex is NOT split into suspended/adherent 2-stage or one linear phenotype. Instead the whole cell's SLOW
state is a **multi-axis CellState** that sets the **rates + initial inventory of LOCAL stochastic events**:

`CellState = { cell_line/clone, EMT, cell_cycle_phase, adhesion/ECM, geometry/confinement, osmotic/volume, perturbation }`

- **Conditional, not additive/scalar.** EMT's effect can flip sign with cell-cycle/adhesion → treat as conditional
  combinations, never a summed scalar or independent addition.
- **Compute compromise:** for sec–min runs hold EMT + cell-cycle as a FIXED physiological context; run only the
  explicit molecular events at runtime (below).
- **Events change reaction RATE or explicit molecular STATE — never overwrite force/tension/filament-count.**
  Topology and forces emerge from mechanics. CellState never directly specifies force/tension.

Runtime explicit events: actin poly/depoly · nucleation · capping/uncapping · severing · crosslinker
bind/unbind/rebind · NMII bind/step/detach · ERM bind/unbind · local myosin/RhoA pulse · cortex rupture + bleb
reassembly.

## §2. Subsystem re-classification — FIRST PASS (Stage 1 refines with a full audit)

`BOUND` = wired + runs in the accepted-step path · `SOURCE-ONLY` = implemented/tested but not wired ·
`PLACEHOLDER` = stub · `GAP` = absent.

| subsystem | status | evidence / what's missing |
|---|---|---|
| monomer transport ↔ polymerization | **SOURCE-ONLY→GAP** | `ac/fluid/transport.MonomerField` exists (meant to RETIRE the infinite scalar `G_actin` in `ff/polymerization_warp.py`); the **conservation coupling** `A_total=A_free+A_polymer+A_sequestered` to poly/depoly is the GAP |
| variable filament length | **SOURCE-ONLY→GAP** | `length_dist='exponential'` implemented (`ff/cortex_assembly.py`) but unexposed; **runtime length change** = GAP |
| growing-N nucleation/severing | **PARTIAL→GAP** | `barbed_end_growth_kernel` / `actin_assembly_kernel` exist (`ff/motility_warp.py`); **severing + depoly-driven N change with persistent IDs** = GAP |
| topology-reforming crosslink CUDA KMC | **SOURCE-ONLY / PLACEHOLDER** | host oracle `ac/weave/crosslink_kmc.py` complete+tested; **device `crosslink_kmc_warp.py` neighbor-pass = PLACEHOLDER, unwired** |
| accepted-step transaction wiring | **BOUND (mech) → GAP (bio)** | outer accept/rollback exists (driver); **topology/monomer/RNG commit-only-on-accept** = GAP |
| device-resident population/ID ledger | **SOURCE-ONLY/BOUND → GAP** | `GlobalCellLedger` + unique actor/filament IDs exist (composition); **growing-N dormant/active device ledger** = GAP |

(Double-count note: the topology-reforming KMC (`kmc_reattach`) and `ff.motility_warp.xl_turnover_kernel`
(`r0_creep`) are mutually exclusive — `assert_single_channel`.)

## §3. Revised execution order

**Stage 0 — close the static-control percolation stage (CURRENT WORK; continue, do not revert).**
0.1 ✅ density=20 (`bc5ff3b0`) + thickness default (`81805535`) — native-verified single-spanning 3D shell.
0.2 Finish the resting-convergence probe (Stage-D / inner-solver re-tune for the now-stiff operator) on this
    static-control cortex; record honestly. **Claim only "static connected CONTROL converged"** (if it does) — never
    "dynamic cortex."

**Stage 1 — re-classify (full audit) + design CellState→conditional inventory/rate schema → REPORT TO PI (before
next coding).** Refine §2; specify the schema data model (how each CellState axis conditionally sets each event's
rate + initial inventory; EMT×cell-cycle×adhesion conditionality).

**Stage 2 — smallest vertical slice: `finite monomer + variable length + accepted rollback`.**
Wire `MonomerField ↔ polymerization` (conserve `A_total`; growth depletes local free monomer, shrink returns);
expose variable length; accepted-step transaction (topology/monomer/RNG commit ONLY on accepted outer step,
rejected = zero change). Built in `ac/weave` CUDA mechanics under `ac/engine` ownership — **NOT** in feature-frozen
`ac/cell`. Do NOT rewrite the whole cortex in one commit.

**Stage 3+ — build up the dynamic cortex.**
3a growing-N (nucleation/severing/complete-depoly; persistent global IDs + one-component ownership; dormant
   capacity pre-allocated, not counted active).
3b topology-reforming crosslink CUDA KMC (finish device neighbor-pass via `SegmentQuery`; wire runtime; double-count
   guard vs `xl_turnover_kernel`; native re-fragmentation gate: connectivity==1). ← the old "emergent binding" plan.
3c head-on-actin via the same KMC; capping/uncapping; NMII bind/step/detach; ERM; local myosin/RhoA pulse; cortex
   rupture + bleb reassembly — each accepted-step-committed rate/molecular-state changes.
3d CellState-conditioned inventory/rates driving all of the above.

## §4. Invariants + prohibitions (throughout)

- Persistent global filament IDs + exactly one component ownership; dormant capacity allocated but NOT counted active.
- Events change rate/molecular-state, **never** overwrite force/tension/count; topology+force emerge from mechanics.
- Topology/KMC commit ONLY on an accepted physical outer step; inner iteration ≠ biological time; **no "every N
  steps" magic number** — derive event probability/batching cadence from physical rate/hazard.
- Finite monomer conserved (`A_total`); a rejected outer step changes monomer/topology/RNG by nothing.
- **PROHIBITED:** tuning density/rate to hit a tension gate · CellState directly specifying tension/force ·
  resampling the filament population to a target distribution every N steps · overclaiming static crosslink
  placement as a dynamic cortex · host-authoritative monomer/topology/event queue · physics conclusions from a
  coarse population · gate loosening.
- Native full-70,686, physiological-baseline throughout; each increment adversarially validated before commit;
  shared-branch = explicit staging + md5, coordinate with Codex.
