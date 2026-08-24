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

# Dynamic cortex — execution plan (PI direction + Lead assessment folded in), 2026-07-23

Refines `CORTEX_REVISED_EXECUTION_ORDER_2026-07-23.md` with the *how* (the methodology the Lead's assessment
surfaced): the two non-negotiable gates, the shared transaction spine, the ensemble-validation shift, the
CellState sourcing discipline, and the KMC time-stepping decision. **The direction is sound; the risk is
engineering + validation + parameterization, not concept.** Build as validated increments; never rewrite the
cortex in one commit.

## §0. Target + framing

- **Target:** a **CellState-conditioned, finite-monomer, variable-length, variable-N, KMC-topology,
  accepted-step-committed DYNAMIC cortex** (PI 2026-07-23), built in `ac/weave` CUDA mechanics under `ac/engine`
  ownership — NOT in feature-frozen `ac/cell`.
- **Reframe (done, native-verified):** `density_per_fil=20` (`bc5ff3b0`) + `overlap_free_cortex` default
  (`81805535`) = a **STATIC CONNECTED CONTROL** (single-spanning, rigid, 0.2µm, 0 interpenetrations). NOT the
  final cortex.
- **CellState = { cell_line/clone, EMT, cell_cycle, adhesion/ECM, geometry/confinement, osmotic/volume,
  perturbation }** sets the **rates + initial inventory of local stochastic events, CONDITIONALLY** (EMT sign
  flips with cell-cycle/adhesion; never scalar-summed). Events change rate/molecular-state, **never** overwrite
  force/tension/count; topology + force emerge from mechanics.

## §1. ⭐ Two non-negotiable gates (before the dynamic build)

**GATE A — bank a STATIC resting baseline first.** A dynamic cortex (more DOF + stochasticity) is *harder* to
converge, not easier; if we build the dynamic machinery before one converged force-balanced cell exists, the
target keeps moving and we may never get a baseline. So: close the static-control resting gate → **one converged
native force-balanced cell**, frozen as the fixed reference. Claim only *"static connected control converged."*
This is itself the project's first converged native baseline — worth banking.

**GATE B — build the accepted-step TRANSACTION ABSTRACTION as shared infra before any biology.** monomer +
topology + population + IDs + RNG must commit atomically on an accepted outer step and roll back to the exact
pre-step state on reject, device-resident, zero host authority (I0-A). If each subsystem invents its own rollback
they will not compose. Build ONE transaction spine first; every dynamic subsystem binds to it.

## §2. Phases (each gated; native full-70,686; ensemble-validated where stochastic)

**Phase 0 — close the static control (CURRENT).**
0.1 ✅ density=20 + thickness default (native-verified single-spanning 3D shell).
0.2 Re-tune the inner solver for the now-STIFF operator (Stage-D showed overshoot: candidate 0.78→2.79 rejected)
    → reduce the ~0.78 membrane-turgor residual toward the 0.21 gate.
0.3 Add the resting cortical-tension source (`resting_bound_myosin`, now transmissible on the connected cortex;
    PI-GAP fraction/force) → **converged force-balanced cell = GATE A.** Honest claim: static control only.

**Phase 1 — the dynamic-cortex SPINE (infrastructure, no new biology).**
1a **Accepted-step transaction** (GATE B): a device transaction over the mutable state arrays; commit-on-accept /
   rollback-on-reject; profiler gate = zero authoritative GPU→CPU roundtrip in the loop.
1b **Device population/ID ledger:** persistent global filament IDs + exactly-one-component ownership; a
   dormant/active free-list so N can grow/shrink on device without host authority; dormant slots never counted
   active.
1c **Ensemble-validation harness:** the validation regime SHIFT — a stochastic cortex is validated by steady-state
   DISTRIBUTIONS across N native realizations (connectivity %, γ, filament count, length dist), not single-run
   asserts. Build this before the stochastic subsystems so "is it correct?" has an answer.
1d **KMC time-stepping decision (deliberate):** at native scale events/step reach millions. Choose + justify
   fixed-τ Poisson (`1−e^{−τr}`, current) vs τ-leaping vs exact Gillespie by accuracy-vs-cost; "derive cadence
   from hazard" ≠ a new magic `every-N-steps`. Document the accuracy bound of the chosen τ.

**Phase 2 — smallest vertical slice: `finite monomer + variable length + accepted rollback`.**
Wire `ac/fluid/transport.MonomerField ↔ polymerization` (retire the infinite scalar `G_actin`): conserve
`A_total = A_free + A_polymer + A_sequestered` (growth depletes LOCAL free monomer; shrink returns it). Expose
`length_dist` (variable length). All mutations go through the Phase-1 transaction. **Gate 2 (ensemble native):**
`A_total` conserved to machine precision across accepted steps; a rejected step changes monomer/topology/RNG by
*nothing*; length distribution matches the sourced form.

**Phase 3 — growing-N: nucleation / severing / complete-depoly.**
Filament count changes with persistent IDs + one-component ownership; **severing changes N only — total contour +
actin mass conserved** (mass returns to the MonomerField on depoly). Dormant capacity absorbs growth.

**Phase 4 — topology-reforming crosslink CUDA KMC + head-on-actin.**
Finish the `crosslink_kmc_warp` device neighbor pass (reuse `SegmentQuery` HashGrid); seed a free crosslinker-end
population; evolve each accepted step; feed bound ends to `link_spring`. **Double-count guard:** select
`kmc_reattach`, DISABLE `xl_turnover_kernel` (`assert_single_channel`). Head-on-actin binds by the SAME KMC (fixes
the 0.13µm splay / 17% capture cap). **Gate 4 (ensemble):** connectivity DISTRIBUTION stays percolated (never a
re-fragmentation the static fix `bc5ff3b0` removed); bound count steady; topology genuinely turns over.

**Phase 5 — full event set.** capping/uncapping · NMII bind/step/detach · ERM bind/unbind · local myosin/RhoA
pulse · cortex rupture + bleb reassembly — each accepted-step-committed rate/molecular-state changes, forces
emergent.

**Phase 6 — CellState conditioning (thin → thick, SOURCED-ONLY).** Start with ONE sourced CellState (the first
baseline). Add an axis / a conditional dependence **only when it has a sourced value** — never build unsourced
conditional branches (the schema must not become a numberless scaffold). EMT/cell-cycle stay fixed context for
sec–min runs; they enter as conditional *rate/inventory setters*, never as force.

## §3. Invariants + prohibitions (throughout)

Persistent global IDs + one-component ownership; dormant capacity not counted active. Events change
rate/molecular-state, never force/count. Topology/KMC commit ONLY on accepted outer step; inner iter ≠ bio time;
no `every-N-steps` magic (derive from hazard). Finite monomer conserved; rejected step = zero change. **PROHIBITED:**
density/rate tuning to a tension gate · CellState → force/tension directly · resampling the population to a target
distribution · overclaiming static placement as dynamic · host-authoritative monomer/topology/event queue ·
physics conclusions from a coarse population · gate loosening. Native full-70,686 + adversarial (ensemble)
validation before every commit; shared-branch = explicit staging + md5, coordinate with Codex.

## §4. Risk register (Lead assessment)

| risk | phase | mitigation |
|---|---|---|
| moving target — never get a baseline | GATE A | close the static resting baseline FIRST; freeze it as the reference |
| non-composing rollbacks | GATE B / 1a | one shared transaction spine before any subsystem |
| stochastic "correctness" undefined | 1c | ensemble/distribution validation harness before stochastic subsystems |
| CellState parameter explosion / unsourced branches | 6 | one sourced baseline; add axes only with sourced conditional dependence |
| KMC accuracy vs cost at native scale | 1d | deliberate time-stepping choice + documented accuracy bound |
| re-fragmentation when topology goes dynamic | 4 | ensemble connectivity gate = percolated distribution, not a single run |
| scope/velocity across a large surface | all | validated vertical slices, never a one-commit rewrite; bank each gate |

## §5. Immediate next action

Continue **Phase 0.2** (inner-solver re-tune on the connected+thick static control) → **0.3** (resting
bound-myosin) toward **GATE A**. In parallel, spec **Phase 1a (transaction)** + **1c (ensemble harness)** as the
spine. Report GATE A + the Phase-1 spine spec to PI before opening Phase 2.
