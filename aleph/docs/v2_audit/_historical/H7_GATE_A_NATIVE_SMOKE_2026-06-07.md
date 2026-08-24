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

# Gate-A native driver smoke — integration verified + the env-var + binder reality (2026-06-07)

Sub-session smoke of the Lead's `scripts/h7_gate_a_native.py` (820f630) on the full native
stack, to catch integration issues before the long Gate-A production launch.

## ✅ Verified
- The driver runs the native constrained integrator (`NativeConstrainedBaoabUpdater`) on
  `build_baseline_cell` (N=15,400): **stable, nonconv=0**, γ_soft + γ_rigid sampled, figure written.

## ⚠️ ACTION 1 — the native compartment forces are env-gated; production MUST set the flag
`build_baseline_cell` only uses the native `FFNRadialShellForce` when
**`FFN_GPU_DEVICE_COMPARTMENTS=1`** (cell.py:1402, a3a63ca). The driver's docstring says
"compartments via env" but does NOT set it itself. Smoke (n_fil=1000):
- env OFF (cpu compartment): **95 steps/s**
- env ON  (native compartment): **167 steps/s**
⇒ the Gate-A launch command MUST export `FFN_GPU_DEVICE_COMPARTMENTS=1` (and
`FFN_GPU_DEVICE_BAOAB=1`), else the run silently falls back to the 1.4 ms/force cpu
compartment path and loses most of the compartment speedup.

## ⚠️ ACTION 2 — the real driver is binder-bound; the 6×/4.6 d figure was optimistic
The everything-native ARM C (`h7_native_gate_a_verify`) measured 498 steps/s, but that timed a
short window in which the **active grip_walk myosin** (the production binder) barely fired. The
real driver, running the active myosin + xlink (each a global `set_snapshot` = 51.9 ms/firing),
is **binder-dominated at ~167 steps/s** with the native integrator + native compartment forces.
⇒ for a realistic Gate-A ETA, the **binder pool is not optional** — it is the dominant remaining
cost. The enabler (`FFNAttachmentSpringForce`, 3a26565) + the pattern (`h7_binder_pool_demo`,
5.3× end-to-end) are ready; wiring the myosin/xlink updaters to `set_attachments` (instead of
`set_snapshot` bond mutation) is the Lead's lane and the key to a tractable Gate-A wall-time.

## Net
Integration works; two production-launch must-dos surfaced (set the compartment env flag; wire
the binder pool). The native integrator + compartment forces are correct + stable in the real
driver; the binder set_snapshot is the live wall.
