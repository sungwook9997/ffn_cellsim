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

# Resting-baseline diagnosis (2026-07-23d) — the resting cortical-tension MECHANISM is BUILT, VERIFIED, and opt-in; physiological activation is a single surfaced PI-GAP

Branch `codex/ff-ac-codex`. This session executes the PI-approved next step from 23c: **build candidate 1 —
the resting bound-myosin setpoint** — as a fine-grained, faithful, additive/default-OFF mechanism. 23c proved
(decisively, at full native) that the resting membrane turgor can ONLY be balanced by CORTICAL TENSION
transmitted through the ERM, that passive crosslink pre-straining is non-conservative (§3), and that a
force-free ERM transmits nothing (§2). The missing pieces were (1) a resting cortical-tension SOURCE and (2) a
force-transmitting ERM coupling. **Both are now implemented and verified against a pure-NumPy oracle.** The
mechanism is ready for the PI to activate with the physiological value; that value is the only remaining GAP.

## TL;DR

- **Source (candidate 1) — resting bound-myosin setpoint.** A fraction of the EXPLICIT Stam-Hocky NMII heads are
  bound to the cortex actin at t0 and each carries an isometric per-head tension through the SAME power-stroke
  crossbridge every running motor uses (`crossbridge_segment_kernel`). No new force law, no lumped tension.
- **Transmission — radial, un-preloaded ERM.** The setpoint engages RADIAL ERM pairing (rest = the real
  membrane-cortex separation, NOT preloaded). The actively-tensioned cortex contracts, stretching the radial
  ERM, which then holds the membrane — the transmission 23c §2 found missing, WITHOUT the forbidden preload.
- **Verified (CPU oracle, runs on any host).** Composing the same conservative families the cell composes, the
  membrane residual drops from the full turgor floor to a minimum **exactly when the emergent cortical γ hits
  the Laplace target `ΔP·R/2 − γ_mem` (§5 of 23c)** — the mechanism, demonstrated, not asserted.
- **PI-GAP (the only thing left).** The resting **bound-head fraction** and the resting **per-head force** are
  absent from the Contract-Graph. Both are UNSET by default ⇒ the mechanism stays OFF (heads unbound at rest,
  bit-identical parity) and the build surfaces
  `PI-GAP: resting bound-myosin fraction + per-head force`.

**No native gate claim.** A TEST fraction proves the MECHANISM; the physiological fraction/force is a PI call,
so the native `< 0.21` gate is NOT closed tonight (it cannot be until the PI supplies the value).

## 1. What was built (files, all in ac/cell + ac/motor)

| file | role |
|---|---|
| `ac/motor/resting_setpoint.py` (new) | `RestingBoundMyosinSetpoint` (source-gated contract, no defaults) + `plan_resting_bound_heads` (pure-NumPy, deterministic seed planner) + `resting_tangential_load_reference` + `apply_resting_bound_heads` (device write). |
| `ac/cell/assemble.py` (edit) | `CellConfig.resting_bound_myosin_{fraction,force_pn,source,capture_um}` (default None/OFF). On activation: plan + seed the heads into the production `SegmentMotorRuntime` state after `enable_segment_runtime`, engage RADIAL ERM pairing, and record the seed ledger. OFF ⇒ untouched build + PI-GAP ledger message. |
| `ac/cell/resting_balance_oracle.py` (new) | pure-NumPy coupled two-shell oracle (turgor + membrane area-gradient + radial unilateral ERM + explicit active cortical dipoles) — the mechanism proof. |
| `ac/cell/driver.py` (edit) | `--resting-bound-myosin-{fraction,force,source,capture}` CLI so the PI can activate the setpoint on the resting driver. |
| tests | `tests/ac/motor/test_resting_setpoint.py` (6 host + 2 CUDA), `tests/ac/cell/test_resting_balance_oracle.py` (5 host). |

### The faithful realization (why it is not a lumped tension)

A seeded head is a real bound crossbridge in the production segment runtime
(`bound`/`seg_a`/`seg_b`/`bary_t`/`abscissa`/`walk_dir`); its force is the same `crossbridge_segment_kernel`
power-stroke spring. The resting per-head tension is realized THROUGH the power stroke: because the head is
bound at the NEAREST point on its actin segment, the head→anchor offset is perpendicular to the walk direction,
so setting `abscissa = f_head/k_xb + r0_xb − ⟨offset, ŵ⟩` makes the tangential crossbridge load equal `f_head`
**exactly** (verified: max |load − f_head| = 1.7e-13 pN in the host planner; the abscissa is a physical,
non-negative power-stroke displacement ≈ `f_head/k_xb`, not a bond pre-strain). Selection is RNG-free (the
`round(fraction·n_heads)` targets are strided over the eligible-within-capture heads); ineligible heads are
skipped and the shortfall recorded — the fraction is never faked up.

## 2. The mechanism, demonstrated (CPU oracle `coupled_resting_balance`, subdiv-3, 20k iters)

The membrane is HELD at its measured physiological radius (7.5 µm — the cell IS this size at rest); only the
cortex relaxes under its active tension. Sweeping the TEST cortical tension:

| per-dipole tension (pN) | emergent γ_cortex (pN/µm) | membrane residual (mean, pN) |
|---|---|---|
| 0 (**OFF**) | 0.00 | **40.74** (= the full turgor floor) |
| 40 | 62.4 | 22.71 |
| 60 | 93.6 | 13.70 |
| **88** | **137.3** (≈ target) | **2.17** (min — 95 % reduction) |
| 90 | 140.5 | 2.31 |
| 120 | 187.4 | 13.36 |

Laplace target `γ_cortex = ΔP·R/2 − γ_mem = 40·7.5/2 − 10 = 140 pN/µm`. **The membrane residual bottoms out
exactly when the emergent cortical γ reaches the derived Laplace target** — the 23c §5 prediction, now produced
by the composed mechanism (not fitted). OFF: γ_cortex = 0, ERM force-free (no preload), membrane at the turgor
floor. ON: the cortex contracts, the radial ERM stretches (mean tension ~40 pN), and the transmitted pull
cancels the turgor. The per-node turgor here (44 pN) is the coarse-mesh value; at native subdiv-6 (40,962
nodes) it equals the diagnosis's 0.644 pN — the reduction is **mesh-independent** (residual → 0 at the target),
so the coarse mesh is a sufficient mechanism proof. The residual does not reach exactly 0 because of the
coarse-mesh pentagonal-vertex non-uniformity (the max stays above the mean — the same mean-vs-max split the
native run shows, 0.6442 vs 0.7766).

## 3. Parity (OFF is bit-identical) + verification

- **OFF** (both PI-GAP fields None, the default): the seed block is skipped, radial pairing is left at the
  config value, and only three INFORMATIONAL ledger keys are added (`resting_bound_myosin_status =
  OFF_UNBOUND_AT_REST`, the PI-GAP message, `n_bound = 0`). The physics state is unchanged — the production
  segment motor is unbound at rest exactly as before.
- **All-or-none, source-gated:** a partial setpoint (fraction without force, or vice-versa) raises at config
  validation BEFORE any CUDA (dev-Mac testable); a complete setpoint requires an explicit provenance source.
- **Tests:** 11 new (6 host planner/contract + 5 host oracle) pass on the dev Mac; 2 CUDA tests (device apply +
  small-cell ON/OFF build) skip on the Mac and run on the gbook A5000. The full `aleph/tests/ac/` suite shows
  only the 2 PRE-EXISTING failures (`ac/engine/ledger.py` hard-coded `cuda:0`; `tests/ac/engine/
  test_surface_body.py` `device='cpu'`) — both in files this session did not touch and unrelated to this work.

## 4. The exact PI-GAP (what physiological activation needs)

Two magnitudes, both absent from the Contract-Graph, neither invented here:

1. **`resting_bound_myosin_fraction`** — the resting bound-head fraction (NMII duty ratio at rest).
2. **`resting_bound_myosin_force_pn`** — the resting isometric per-head tangential tension [pN].

Their product sets the resting cortical tension `γ_cortex` the ERM transmits; the Card-3 pressure-consistency
gate lands the required `γ_cortex = ΔP·R/2 − γ_mem` (140 pN/µm at 40 Pa; 260 at 72 Pa). When the PI supplies a
sourced fraction + per-head force (and a source label), activate with, e.g.:

```
python -m aleph.components.incumbent.driver --from-resting --n-filaments 70686 \
  --resting-bound-myosin-fraction <PI> --resting-bound-myosin-force <PI pN> \
  --resting-bound-myosin-source "<citation>" \
  --membrane-subdivisions 6 --inner-solver erm_gauss_seidel_tournament --n-inner 1600
```

The native `< 0.21` gate stays **OPEN** until a PI-supplied physiological fraction actually drives it there
(it will not tonight — the fraction is the GAP). Do NOT loosen the gate, hand-swap ΔP, or default the fraction.

## 5. Reproduce

```
# CPU oracle (dev Mac, no CUDA):
python -c "from aleph.components.incumbent.resting_balance_oracle import coupled_resting_balance as c; \
  print(c(dipole_fraction=0,per_dipole_tension_pn=0,subdivisions=3).membrane_residual_mean_pn, \
        c(dipole_fraction=1,per_dipole_tension_pn=88,subdivisions=3).membrane_residual_mean_pn)"
pytest aleph/tests/ac/cell/test_resting_balance_oracle.py aleph/tests/ac/motor/test_resting_setpoint.py -q
# device apply + small-cell ON/OFF build gates (gbook A5000): the two skipped CUDA tests above.
```
