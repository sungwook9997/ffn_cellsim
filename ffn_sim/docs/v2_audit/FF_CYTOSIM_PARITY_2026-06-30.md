# FF ↔ Cytosim parity oracle — first cross-check (Stage 6f)

**Date:** 2026-06-30 · **Engine:** `ffn_sim/ff/` · **Branch:** `dcm/main` · **Status:** oracle BUILT + first finding → PI

## What changed

The runnable-Cytosim parity oracle is **no longer blocked**. Cytosim (Nédélec & Foethke's own C++
code, the independent reference) was cloned from gitlab and the headless `sim` binary built
(user-authorized; macOS arm64, Apple Accelerate BLAS/LAPACK, no GUI). This is the independent oracle
ENGINE.md's design requires — "the active layer is never validated against only itself".

Cytosim's config units are **s · µm · pN** — identical to the FF unit system (`ff.units`), so
quantities compare with no conversion. Harness: `ff/cytosim_parity.py` (degrades gracefully when no
binary; `tests/ff/test_cytosim_parity.py` skips without `CYTOSIM_SIM`, runs the real check with it).

## First cross-check — bending energy of a fixed circular arc

Same shape (a circular arc, R_arc=5µm, L=2µm, κ=20 pN·µm²) placed at identical model-points in both
engines; compare the static bending energy to each other and to the continuum κ·L/(2R²)=0.80 pN·µm.

| n_pts | FF | Cytosim | FF/analytic | Cytosim/analytic | FF·(n−1)/(n−2) |
|---|---|---|---|---|---|
| 5  | 0.5998 | 0.7997 | 0.750 | 0.9996 | 0.7997 |
| 9  | 0.6999 | 0.7999 | 0.875 | 0.9999 | 0.7999 |
| 17 | 0.7500 | 0.8000 | 0.937 | 1.0000 | 0.8000 |
| 33 | 0.7750 | 0.8000 | 0.969 | 1.0000 | 0.8000 |

**Finding.** Cytosim is **continuum-accurate at every resolution** (even n=5). FF's NF2007
interior-triple discrete energy **under-counts by exactly the end-factor (n−2)/(n−1)** — `FF·(n−1)/
(n−2)` recovers the continuum to 5 digits at all n, and FF→Cytosim as n→∞. At the cortex's n=7
beads/fiber the factor is 5/6 ⇒ **~17 % low** in the bending energy and force.

**Why.** Both are valid discretisations of the same continuum operator (κ/2)∫c²ds. FF is faithful to
the *paper's published* interior-vertex sum (`forces_warp`, NF2007 p9); Cytosim (the *code*)
end-corrects the quadrature so the coarse-mesh energy already equals the continuum. The analytic-only
checks we had (force=−∇E, the κq⁴ dispersion) are *self-consistent within FF* and therefore could not
see this — it took the independent oracle. This is exactly the value of building Cytosim.

## Impact + decision (PI)

- **γ-floor conclusion UNAFFECTED.** Bending is sub-dominant to the myosin prestress + turgor in γ
  (the actin axial tension that carries γ comes from the inextensibility constraint, not bending), so
  the 6d result (actomyosin γ ~10³× under band, matching the BAOAB-MD g_soft) stands.
- **Quantitative-fidelity choice — surfaced, NOT silently changed** (changing a force law is
  PI-gated; hard rule). Options:
  1. Keep FF faithful to the paper's interior sum; document the O(1/n) end-bias vs the oracle (it
     vanishes as fibers are refined; cortex fibers are coarse at 7 beads).
  2. Adopt Cytosim's end-corrected quadrature in the FF bending kernel (continuum-accurate at the
     cortex's small n) — a one-line-ish change to the energy/force assembly, re-validated against the
     oracle. Recommended if cortex-scale bending accuracy matters for a downstream observable.
- **Next parity targets** (same harness): single-fiber relaxation rate (μκq⁴ slowest mode),
  Euler buckling threshold π²κ/L², and the Hand (motor) stepping/detachment kinetics — to cross-check
  the §10.1 layer against Cytosim's `motor`/`couple` implementations.

## Reproduce

```
git clone https://gitlab.com/f-nedelec/cytosim && cd cytosim
cmake -B build -DCMAKE_BUILD_TYPE=Release && make -C build sim
CYTOSIM_SIM=$PWD/build/bin/sim python -m ffn_sim.ff.cytosim_parity
```
