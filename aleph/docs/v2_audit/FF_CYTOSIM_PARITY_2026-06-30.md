# FF ↔ Cytosim parity oracle — first cross-check (Stage 6f)

**Date:** 2026-06-30 · **Engine:** `aleph/laws/` · **Branch:** `dcm/main` · **Status:** oracle BUILT + first finding → PI

## What changed

The runnable-Cytosim parity oracle is **no longer blocked**. Cytosim (Nédélec & Foethke's own C++
code, the independent reference) was cloned from gitlab and the headless `sim` binary built
(user-authorized; macOS arm64, Apple Accelerate BLAS/LAPACK, no GUI). This is the independent oracle
ENGINE.md's design requires — "the active layer is never validated against only itself".

Cytosim's config units are **s · µm · pN** — identical to the FF unit system (`ff.units`), so
quantities compare with no conversion. Harness: `aleph/validation/cytosim_parity.py` (degrades gracefully when no
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

## Resolution (PI-authorized 2026-06-30) — end-correction ADOPTED

PI chose option 2: **adopt Cytosim's end-corrected quadrature** in the FF bending kernel. Grounded
directly from the oracle's source (`chain.cc::bendingEnergy0`):

```
e *= ( lsp + 1 ) / ( fnCut * lsp );   // lsp = nPoints-2, fnCut = segmentation
// "we only considered (nPoints-2) junctions, ... only a fraction of the total length"
```

Decoding it: Cytosim's energy = (κ/seg)·Σ(1−cosθ) × (n−1)/(n−2), and FF's `|m_{i-1}−2m_i+m_{i+1}|²`
form equals (κ/seg)·Σ(1−cosθ) exactly (since |Δ²|²=2seg²(1−cosθ)). So the ONLY difference is the
factor **g_f = (n−1)/(n−2) = p_f/(p_f−1)** (p_f = segments). It depends only on the point count
(topology), so multiplying each fiber's α by g_f scales energy AND force by the same constant ⇒
**force=−∇E preserved, q⁴ dispersion shape unchanged**. Implemented in `forces_warp._per_triple_alpha`
(`end_correction=True` default; `=False` recovers the paper-literal interior sum).

Result after adoption: **FF now matches Cytosim to FFc/cyto = 1.0000 at every resolution** (n=5→33),
including the cortex's coarse mesh. `tests/ff/test_cytosim_parity.py` asserts the match (and that the
raw interior sum still shows the (n−2)/(n−1) factor, documenting what the correction fixes).

- **γ-floor conclusion UNAFFECTED — re-verified.** Bending is sub-dominant to myosin/turgor in γ (the
  actin axial tension that carries γ comes from the inextensibility constraint, not bending). The
  grounded production γ_active went 1.52e-4 → **1.51e-4 mN/m** with the correction (still ~2300× under
  band, still on the BAOAB-MD g_soft). The conclusion is robust to the bending quadrature.
- **Next parity targets** (same harness): single-fiber relaxation rate (μκq⁴ slowest mode),
  Euler buckling threshold π²κ/L², and the Hand (motor) stepping/detachment kinetics — to cross-check
  the §10.1 layer against Cytosim's `motor`/`couple` implementations.

### Dynamic relaxation-rate parity — ATTEMPTED, DEFERRED (2026-06-30)

Tried the single-fiber relaxation-rate cross-check (per-frame energy via Cytosim's `repeat{ run;
report }` block — that mechanism works). It is **not a clean parity target**: the absolute decay rate
depends on Cytosim's drag MODEL (`drag_radius` / `drag_length` / `surface_effect` cylinder drag),
which differs from FF's `units.fiber_mobility` (μ=log(L_h/δ)/3πηL) by an O(1)–O(10) factor, and
Cytosim integrates *implicitly* (its own NF2007 Eq 2) so it relaxes stiff bending almost instantly —
even the convention-independent normalized E(t)/E0-vs-(t/τ) curves did not overlay cleanly in the
sampled window (Cytosim under-relaxed at matched wall-time; FF followed a clean exponential).
**Conclusion:** matching dynamic rates needs deep Cytosim drag-model calibration whose validation
value is low given the static energy + functional parity already pins the bending convention
(convention-independent). Deferred to a focused session that first calibrates the drag model. The
Hand/motor-kinetics parity (higher value for the γ result) similarly needs a dedicated setup.

## Reproduce

```
git clone https://gitlab.com/f-nedelec/cytosim && cd cytosim
cmake -B build -DCMAKE_BUILD_TYPE=Release && make -C build sim
CYTOSIM_SIM=$PWD/build/bin/sim python -m aleph.validation.cytosim_parity
```
