# FF Stage 6Z — rigid AFM plate fix: cortex was floating through the soft plate; force was a ~10³× artifact

**Date:** 2026-07-02  **Engine:** FF (Warp, A5000)  **Branch:** dcm/main
**PI catch:** "누를 때 제대로 면으로 누른 거 맞냐? 면 위로 코르텍스가 떠오르는데" — the AFM plate wasn't
actually pressing the cell; the cortex floats above the plate. Also: do it NATIVE; show cytoplasm + membrane.

## One line

The virtual-AFM plate was a **soft penalty too weak to confine the (stiff) cortex** — at native, 15–75% of
cortex nodes floated ABOVE the plate and the reported plate force was a ~10³–10⁵× artifact of
Σ(k_plate·penetration). A **rigid plate (hard z-clamp)** confines the cell exactly (0% above) and gives
**realistic ~nN forces** — and the 6V confinement-independent γ≈0.15 mN/m **survives** the fix.

## The bug (soft penalty)

`plate_kernel` was a one-sided spring penalty k_plate·(|z|−half_gap). k_plate = 10× the turgor breathing mode
is far too soft to overcome the cortex's resistance to flattening, so the cell stayed ~spherical and poked
through the plates (measured on the committed 6W snaps):

| strain | half_gap | max\|z_cortex\| | penetration | cortex nodes ABOVE plate |
|---|---|---|---|---|
| 0.15 | 6.375 | 6.617 | +0.24 | 15% |
| 0.60 | 3.000 | 3.927 | +0.93 | 58% |
| 0.78 | 1.650 | 2.903 | +1.25 | **75%** |

Making k_plate stiffer doesn't help: dt_mu = 0.1/kmax ∝ 1/k_plate, so a stiffer plate just shrinks the step
(CFL) and relaxes less per step. The penalty approach is fundamentally CFL-limited for a stiff cortex.

## The fix — rigid plate (hard constraint, no CFL penalty)

`rigid_plate_step_kernel` (+ `rigid_plate=True` in `simulate_whole_cell_compression_on_device`): do the
overdamped step, then HARD-clamp |z−cz| ≤ half_gap (a geometric constraint, not a force → no CFL throttle).
The AFM reaction is the removed outward displacement / dt (= the force the plate exerts to hold each node),
summed on the top plate — the physical plate force by Newton's third law.

## Native result (N=70686, rigid plate, all compartments; cortex+turgor+nucleus+membrane)

| strain | cortex above plate | R_eq | **F_plate** | γ_apparent | V/V0 |
|---|---|---|---|---|---|
| 0.15 | **0.0%** | 7.51 | **3.3 nN** | 0.150 | 0.97 |
| 0.35 | **0.0%** | 7.51 | **18.7 nN** | 0.150 | 0.86 |
| 0.55 | **0.0%** | 7.84 | **38.8 nN** | 0.157 | 0.68 |

- **Exact confinement:** 0% of cortex nodes above the plate at every strain (the flat caps sit ON the plates;
  the cell is a proper pancake). The floating is gone.
- **Realistic force:** F_plate = 3–39 nN over 15–55% strain — the real AFM range. The soft plate had reported
  ~8500 nN (8.5e6 pN) at strain 0.1. **So FF_STAGE6V "gap #1: AFM force ~10³× too high" was substantially a
  SOFT-PLATE MEASUREMENT ARTIFACT** (Σ k_plate·penetration over deeply-penetrating nodes), not real cortex
  over-stiffness. A residual cortex resistance remains (at 55% strain F≈39 nN vs turgor·contact≈6 nN → ~6×,
  the real over-stiff-cortex tail), but it is ~10× not ~10³×.
- **Confinement-independence HOLDS:** γ_apparent ≈ 0.150 mN/m (0.150/0.150/0.157) across strain — the 6V/6W/6X
  conclusion survives the rigid plate. Under regulated ΔP the cell loses volume (V/V0 0.97→0.68) rather than
  bulging (R_eq ~7.5), so ½·ΔP·R_eq stays ~0.15. (Fischer-Friedrich interphase 0.17; confinement-independent.)

## Cytoplasm + membrane — implemented and ACTIVE (were just not drawn)

The PI asked "왜 cytoplasm이랑 plasma membrane은 안 보이냐 구현되어있는거 맞음?" — they were implemented but
not *drawn*. They are forces/media, not particle clouds:
- **Cytoplasm** = the turgor-pressurised incompressible interior (ΔP_turgor = 40 Pa held; V_cyto = V_hull −
  V_nuc tracked, V/V0 → 0.68) + the η = 65.9 Pa·s viscous drag (used in the DYNAMIC relax/BAOAB paths; the
  static AFM equilibrium is η-independent, so η sets the path not the equilibrium force).
- **Plasma membrane** = the γ_mem = 10 pN/µm surface tension (ΔP_mem = 2.7–3.0 Pa inward, active every step;
  6X).
Both are now drawn in the cross-section (cytoplasm fill + membrane envelope) with their live metrics.

## Sanity gate
- No magic number: rigid plate is a geometric constraint (no tuned stiffness); reaction is the physical
  removed-force. ✓
- Native: all numbers at N=70686. ✓
- No regression: `rigid_plate=False` default keeps the historical soft path; new option is additive; test
  asserts exact confinement. ✓
- Honest: the soft-plate artifact is stated; gap #1 corrected (force artifact, ~10× residual not ~10³×);
  confinement-independence re-confirmed with the rigid plate. ✓

## Figures
- `outputs/ff/figs/wholecell_compartments_native.png` — native rigid-plate cross-sections at strain
  0.15/0.35/0.55 showing ALL 4 compartments (cytoplasm fill, membrane envelope, cortex, nucleus) confined
  between the rigid plates (0% above), with live ΔP_turgor / ΔP_mem / V/V0.

## Files
- `ff/network_warp.py` — `rigid_plate_step_kernel` + `rigid_plate` option in
  `simulate_whole_cell_compression_on_device`.
- `tests/ff/test_compression.py` — `test_rigid_plate_confines_cortex_exactly`.

## Follow-ups (PI)
- Re-render the 6U/6W AFM figures (which used the soft plate) with the rigid plate; the soft-plate F_plate
  magnitudes in 6U/6V/6W/6X should be read as artifacts (the γ conclusions are unaffected — they come from
  turgor-Laplace / method-of-planes, not the plate force).
- The proper force-derived γ (Fischer-Friedrich contact-mechanics extraction) is now feasible with a real
  nN-scale F_plate — a cleaner γ than the ½·ΔP·R_eq proxy (6W caveat).

Related: [[project-cell-mechanics-extend]], FF_STAGE6X (membrane), FF_STAGE6W (nucleus), FF_STAGE6V (γ pressure-
borne — gap #1 corrected here), FF_STAGE6Y (arrangement).
