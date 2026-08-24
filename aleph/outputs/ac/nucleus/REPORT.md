# nucleus track — REPORT (I2 deformable-mesh nucleus, analytic foundation)

Track: **nucleus** (I2 → I7 later) · branch `ac/nucleus` · base `ac/new-engine`.
Status: I2 analytic oracle suite committed + visualized (**41 host-numpy tests**, dev Mac; 66/66 across
`tests/ac/`). Warp kernel SOURCE authored (`ac/nucleus/envelope.py`) but NOT run — the dev Mac has a
CPU-only Warp build; the lead runs the native CUDA gates on gbook (serial spine). Copies the fluid-spine
`REPORT.md` + `figs/` pattern (AC_PARALLEL_SESSIONS §1.9).

## What landed (NET-NEW under `ac/nucleus/`, no `ff/` edit)

The DCM/FF nucleus was a radial-spring **bead cloud** (one deformation mode, no bending, no rupture).
This rebuilds it as a genuinely **deformable triangulated shell** using the membrane Helfrich machinery:

- `geometry.py` — oblate-capable icosphere, dihedral hinges, divergence-theorem volume + exact gradient.
- `lamina_analytic.py` — dihedral bending (reuses the ff/membrane pattern, κ̃=8πκ/Σ_ref) + the
  **framework-#6 lamin split** (chromatin + lamin-B soft → lamin-A/C strain-stiffening past the knee) +
  Young-Laplace + EMERGENT rupture.
- `chromatin_analytic.py` — Marko–Siggia WLC internal polymer net (soft small-strain; finite-extensible).
- `linc_analytic.py` — nonlinear LINC tether + **pre-authored I7** capstan ∮T·κ ds + oblate-flatten.
- `mask_provider.py` — **`DeformableNucleusMaskProvider`** (§1.4 fluid interface) + static-sphere
  test-double + `NucleusBoundaryProvider` Protocol.
- `envelope.py` — the Warp kernel SOURCE (per-face lamin tension + rupture, incompressible volume,
  nucleoplasm viscosity, LINC tether, chromatin WLC), reusing `ff.membrane_surface.helfrich_bending_kernel`.
- `params_i0b2.yaml` — the I0-B2 ledger (real KB claim IDs; 5 GAP-PI unknowns).
- `INTEGRATION.md` — the deferred `ff/`/fluid patch-notes + the native-gate spec.

## Figures

Regenerate all: `PYTHONPATH=. python ffn_sim/scripts/ac_nucleus_vis.py` (pure numpy/matplotlib on the
committed oracles). Every figure overlays the closed-form oracle/band on the numeric result, annotates SI
units, and does not truncate axes (the visualization-integrity rules).

- **`figs/i2_bending.png`** — dihedral Helfrich bending. The raw icosphere dihedral sum Σ converges to the
  topological constant ≈7.50 and is **R-invariant** (R=1.5, 6.0 µm both land at 7.50), so κ̃=8πκ/Σ_ref
  makes the resting sphere reproduce the continuum Willmore energy **8πκ=2.0810 pN·µm exactly**
  (Magic-Number Block satisfied). The FD-gradient scatter sits on **y=x** (force = −∂E/∂x, the sign
  arbiter) with **Σf=1.6e-14 pN** (translation-invariant).
- **`figs/i2_lamin_knee.png`** — framework-#6 areal elasticity. σ(ε) is a continuous bilinear: soft
  chromatin+lamin-B below the knee (ε=0.10), strain-stiffening lamin-A/C above. The tangent areal modulus
  **jumps ×3.00 UP** at the knee (K_soft→K_laminAC) — the strain-stiffening crossover (KB-3.B2.2), a
  **REPORTED** ratio, not tuned.
- **`figs/i2_volume.png`** — incompressible nucleoplasm. Signed-tet mesh volume → continuum (4/3)πR³ with
  refinement; oblate flatten at **constant volume** gives a∝A^⅓ with V/V₀≡1 (the I7 Khatau flatten
  kinematics, pre-authored).
- **`figs/i2_rupture_capstan.png`** — EMERGENT rupture (σ→0 past the threshold ε=0.5 GAP-PI; clean on/off,
  threshold not tuned) + the **I7 capstan** ∮T·κ(n̂·d̂)ds reproducing the classical belt-over-cylinder
  **2T·sin(φ/2)** (net 2T at a full half-wrap) — tangential cap tension becomes a normal pressure, NOT a
  radial strut.

## Gates owned (CPU host numpy — 41 tests, all green)

| Gate | File | Result |
|---|---|---|
| sphere → 8πκ (Willmore), κ̃-calibrated | test_lamina_bending | exact at subdiv 2/3/4 |
| Σ_ref grid- & R-invariant (≈7.50) | test_lamina_bending | R-independent to 1e-6 |
| bending force = −∂E/∂x (FD sign arbiter) | test_lamina_bending | <1e-6·‖f‖; Σf<1e-10 |
| flat hinge → 0 force | test_lamina_bending | <1e-12 |
| lamin knee: continuity + tangent jump UP | test_lamina_area | ×ratio crossover (stiffening) |
| Young-Laplace ΔP=2σ/R | test_lamina_area | exact |
| volume → continuum + gradient FD + translation-inv | test_volume | <1e-6·‖g‖; rel 1e-10 |
| rupture EMERGES above / ABSENT below (on/off) | test_rupture | falsifiable step, any threshold |
| WLC force=dE/dx, monotone, diverges, small-strain k | test_chromatin | <1e-6·‖f‖ |
| LINC force=dE/dx, force-free rest, tension-only | test_linc_capstan | <1e-5·‖f‖ |
| capstan ∮T·κ ds = 2T (half-wrap) + ring cancels | test_linc_capstan | rel 2e-3 |
| oblate flatten conserves volume (a∝A^⅓) | test_linc_capstan | rel 1e-9 |
| DeformableNucleusMaskProvider ≡ oblate ellipsoid; Protocol; moving-mask | test_mask_provider | exact match |

## Native-gate viz spec (lead runs on gbook)

Interactive 3-D cell-morphology **HTML** rendering the **deformable oblate nucleus mesh** inside the
membrane-minus-nucleus fluid domain — real geometry, peel/slab/cut views, **full-res (no downsampling)**,
**browser-verified** (`browser_check.py`; a WebGL grep is meaningless — the render must be eyeballed).
Fields/geometry to NAME in the render: the envelope mesh colored by per-face areal strain (chromatin vs
lamin-A/C regime) + rupture flags; the LINC tether stubs; the reconciled I1a no-flux mask boundary. See
`INTEGRATION.md` §Native-gate spec for the physics acceptance bands.
