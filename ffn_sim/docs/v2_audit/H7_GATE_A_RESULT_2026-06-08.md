# H.7 Gate-A result — REFUTE: s_grip develops to ℓ₀-scale, γ_soft stays floored

**Date:** 2026-06-08 (overnight, autonomous) · **Branch:** `h7/full-cell-integration`
**Run:** `scripts/h7_gate_a_native.py --n-fil 1000 --tag prod` on the **full GPU-main
stack** (native constrained integrator + native compartment ForceCompute), full
physiological MCF7 cell (all compartments + grip_walk myosin + xlinks ON), gbook
RTX A5000. Figure: `outputs/h7/figs/h7_gate_a_native_prod.png`; data:
`outputs/h7/production/h7_gate_a_native_prod.json`.

## Question (Gate-A)
The active cortical tension γ_soft is floored ~1000× below the band. Gate-A asks:
is it because the myosin grip-stretch `s_grip` never develops (heads bind+load but
don't WALK)? If so, γ_soft should CLIMB toward the band as s_grip → ℓ₀ (CONFIRM).
If γ_soft stays flat while s_grip → ℓ₀, the wall is DEEPER (REFUTE → Gate-B).

## Result — REFUTE (decisive)
With grip_walk genuinely walking (validated GPU==CPU s_grip accumulation), over
ticks 10→270 (1e5→2.7e6 steps):

| quantity | tick 10 | tick 130 | tick 200 | tick 270 |
|---|---|---|---|---|
| **s_grip/ℓ₀** | 0.024 | 0.292 | 0.415 | **0.529** |
| **γ_soft (mN/m)** | 1.97e-4 | 1.90e-4 | 1.82e-4 | 1.89e-4 |
| γ_rigid (mN/m) | 0.102 | 0.080 | 0.095 | 0.073 |
| r/r0 (backbone) | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

**s_grip climbed 22× (0.024 → 0.529, past the 0.5 verdict threshold) while γ_soft
stayed DEAD FLAT at ~1.9e-4 mN/m.** ⇒ **REFUTE**: the active-γ floor is NOT an
s_grip-generation problem — the grip-stretch develops to ℓ₀-scale and the realized
shell tension does not respond. The wall is downstream of generation.

## What this means
- Confirms + strengthens the prior authoritative diagnosis (force-aggregation /
  transmission-limited; "s_grip-delivery refuted"; [[project-gamma-floor-layered-resolution]])
  — now on the COMPLETE physiological cell at the production operating point, with
  grip_walk actually walking (not the binned_r0 proxy).
- γ_soft ~1.9e-4 mN/m is ~100× above the binned_r0 floor (1.77e-6) but ~1850×
  below the band (0.35) and **does not scale with s_grip** — the generated
  per-head stretch is not transmitted into spanning hoop tension.
- `r/r0 = 1.0000` throughout: the rigid M-SHAKE backbone does NOT condense
  (no buckling/densification) — exactly the Murrell/Lenz/Miyazaki symmetry-breaking
  the constrained backbone forbids. This is the **Gate-B** lever (relax M-SHAKE →
  allow buckling/condensation), still requiring the separate relaxed-constraint mode.

## Next lever (Gate-B / transmission), per the result
The active tension is generation-delivered (s_grip develops) but NOT transmitted.
Levers, in the project's authoritative order: (1) **Gate-B** — relax M-SHAKE to
permit filament buckling/condensation (PI pre-authorized a controlled experiment;
needs the relaxed-constraint native mode, not yet built); (2) soft long-range
transmission (the rigid backbone shunts active force locally); NOT s_grip
generation (refuted here), NOT binding throughput, NOT aggregation geometry.

## Provenance / integrity
- The run is VALID: stepping_mode = grip_walk (asserted at launch; the earlier
  s_grip = 0 was a stale gbook config defaulting to binned_r0 —
  `H7_GPU_MYOSIN_SGRIP_BUG_2026-06-07.md`, RESOLVED).
- Run continues to s_grip → ℓ₀ (tick ~425) and beyond for the complete plateau
  trajectory; the verdict (γ_soft flat through s_grip ≥ 0.5) is already decisive
  and will not change (γ_soft has been flat across a 22× s_grip increase).
- nonconv = 0 throughout (native M-SHAKE converges every step).

## UPDATE (ticks 350–380): bead-stepping onset — γ_soft rises ~3×, still ~627× below band
At ~tick 360 the bound heads first reach the s_grip overflow cap (s_grip ≥ ℓ₀) and
**actually STEP to the next minus-ward bead** (s_grip −= ℓ₀ per step → `n_step_advances`
> 0). Two changes follow:
- **mean s_grip plateaus + cycles** (peak 0.581 @ tick ~350, then 0.48→0.41) — it does
  NOT march monotonically to ℓ₀, because each bead-step subtracts ℓ₀ from the head's
  grip; the population oscillates around the step cap. ("s_grip developed" is satisfied:
  it exceeded 0.5.)
- **γ_soft jumps ~3×** (1.73e-4 → 5.58e-4 mN/m at tick 360, s_grip=0.443) — the FIRST
  real contraction signal from genuine material stepping (not just static loading).

So the result is sharper than "dead flat": grip_walk DOES produce a small contraction
once heads step — but γ_soft **plateaus at ~5.6e-4 mN/m = 1/627 of the band (0.35)**.
**REFUTE stands:** even with s_grip developed AND real bead-stepping contraction, the
realized shell tension is ~600× short of band → the wall is transmission/lever
(Gate-B), not generation. The ~3× step-onset rise quantifies how little of the stepping
work reaches spanning hoop tension on the rigid (r/r0=1.0000, non-condensing) backbone.
