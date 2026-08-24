# INTEGRATION_ecm_world — far-field boundary KERNEL_BOUND increment

Owner file this session: `aleph/engine/ecm_world.py` + `aleph/tests/ac/engine/test_ecm_world.py`
only. No shared file was edited. Coordinates with `codex/ecm-gpu-topology` (which owns
`aleph/components/ecm/` device schema / Mikado topology / remodel) — the far-field boundary is a **separate
participant** per ECM_WORLD_PLAN §2.3 and is not touched by codex's commits, so there is no overlap.

## What landed (SEAMED → KERNEL_BOUND for the ECM boundary sub-runtime)

`FarFieldDirichletBoundaryRuntime` — a concrete `BoundaryAnchorDelegate` backed by three real Warp
kernels, replacing the pure-protocol/spy boundary seam:

- `_far_field_reaction_ledger_kernel` — **reaction-preserving**: records the constraint reaction
  `R = -f` on each pinned node into `reaction_d` *before* zeroing that node's force, so cell traction
  reaching the truncated domain is ledgered, not discarded (ECM_WORLD_PLAN §3).
- `_far_field_reaction_rollback_kernel` / `_far_field_reaction_commit_kernel` — accepted-step
  transaction: reaction + reaction-conjugate boundary work (`sum(R . delta_target)`, exactly 0 for a
  static baseline anchor) commit only under the device `accepted[0]` scalar; no private clock.

Ownership gates: all six ledger arrays are CUDA-resident, dtype/shape-checked, one device, distinct
storage; an empty (unanchored) pinned set is rejected. No magic numbers (kinematic pin + reaction
readout has no material constant; collagen constitutive law stays P2-gated and untouched).

Evidence: CPU-green structural gates pass on the Mac (fake-CUDA doubles); a `skipif(not CUDA)` numeric
gate (`_far_field_reaction_kernel_preserves_then_removes_pinned_force`) is authored for the GPU lane to
close the CUDA_UNIT force/sign/boundary check.

## Requested Lead actions (shared files — PI/Lead owned, not edited here)

1. **Typed global boundary ledger schema (ECM_WORLD_PLAN §11 gap #5).** `LedgerContributor` defines
   only methods. `FarFieldDirichletBoundaryRuntime.accumulate_ledger` forwards to a duck-typed
   `ledger.add_far_field_reaction(committed_reaction_d, boundary_work_d)` when present, else is a safe
   no-op. Add a `add_far_field_reaction(reaction_d: wp.array(vec3d), work_d: wp.array(float64))` slot to
   the typed device ledger so the global acceptance gate reads reaction/work under a stable name and can
   close `sum(reaction) + cell-traction-resultant = 0`.
2. **World dump_state / dispatch registration.** Register the concrete runtime instance so the
   `ecm_far_field_anchor` edge dispatches this boundary once per candidate (it already plugs into the
   existing `ECMBoundaryAnchorFacade`; only instance wiring + dump_state exposure of `reaction_d` /
   `boundary_work_d` is needed).

## Remaining gap to CUDA_UNIT

Run the authored `skipif` reaction gate on the A5000 (force/sign/zeroing), plus commit/rollback
bit-exact restoration and boundary-work closure on real CUDA. Continuum-impedance boundary mode and the
prescribed-motion boundary-work path remain future (baseline is `far_field_dirichlet`, static target).
