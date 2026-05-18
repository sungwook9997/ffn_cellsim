# Phase 1 — Unit 1.1 Progress Report (P1-corrected)

**Date**: 2026-05-18 (Day 1, second pass)
**Scope**: ECM Unit 1.1 — 2D Mikado collagen fiber network with **real segment-intersection cross-links**, **integrated XL mechanics**, energy/force only (no time integration, no cells).
**Repo HEAD**: `4ffcbbe` (untracked `acs_kb/`)
**Notion source**: KU-1.1, KU-1.2, KU-1.3, KU-1.7, KU-1.9, KU-1.22, KU-1.24, KU-1.25, KU-1.27, KU-1.28
**Mode**: production (`demo_mode = false`)
**Previous report**: superseded by P1 corrections below.

---

## P1 corrections — what changed since the first pass

Per PI, the first MVP was acknowledged as smoke-passing but **not scientifically complete**. The four required corrections, all implemented:

1. **Cross-link generation is now real 2D Mikado segment intersection** (`acs_kb/ecm/cross_links.py:_segment_intersections`): pairwise AABB pre-filter → Cramer's-rule line-line solve → in-segment (0 < t,u < 1) acceptance. The previous bead-neighbor proxy is gone. The Phase 1 generator still anchors each crossing to the **nearest bead** on each fiber (discrete approximation, documented), and periodic-image intersections are not enumerated (~5 % boundary fibers — to be added in Phase 2). Both approximations are flagged in the module docstring.
2. **Cross-link harmonic energy/forces are integrated into total network mechanics** (`fiber_mechanics.compute_energy / compute_forces` now accept `cross_links=…`). Newton 3rd law verified pairwise. FD vs analytic gradient agrees to 1.9e-7 relative including XL contribution.
3. **`⟨z⟩ is now emergent**, not constructed.** The generator no longer samples to a target; it returns every intersection it finds. The report distinguishes the emergent measurement from both the Mikado-theory prediction and the KU-1.3 biological reference.
4. **ξ (biological mesh) and ℓ_c (Mikado segment) are separated** in config (`biological_mesh` vs `target_segment_length`), in derivation (`mikado_ell_c_predicted`), and in the gate output. They are equated by modelling convention (Storm–MacKintosh) but tracked as distinct quantities.

## What is built

```
acs_kb/
├── common/
│   ├── derived_params.py     # KU-1.27 invert: N_f = π L_box²/(2 ℓ_c L_fiber)
│   └── sanity_gate.py        # ξ / ℓ_c / ⟨z⟩ distinct; biology gap logged
├── configs/
│   └── phase1_unit1.yaml     # primary scales only; everything else derived
├── ecm/
│   ├── fiber_network.py      # 2D Mikado + unwrapped endpoints for intersection
│   ├── fiber_mechanics.py    # WLC backbone + XL coupling, single API
│   ├── cross_links.py        # true 2D segment intersection + harmonic XL kernel
│   └── visualization.py
├── tests/test_ecm.py         # 7 tests, all KU-tagged
├── notebooks/01_first_fiber_network.py
└── outputs/phase1/
    ├── REPORT.md             # this file
    ├── 01_visual_iso.png
    ├── 02_emergent_z.png     # NEW: emergent vs theory vs biology
    ├── 03_order_check.png
    └── 01_summary.json
```

## KU validation numbers (P1-corrected)

| Quantity | Measured | Theory / reference | KU | Status |
|---|---|---|---|---|
| Mikado ℓ_c (Buffon, π/(2ρ_L)) | **2.04 μm** (mean over 10 seeds) | predicted 2.00 μm, biological ξ = 2 μm | KU-1.27 / KU-1.7 | PASS ([1, 4] μm band) |
| Crossings per fiber | **4.89 ± 0.03** | predicted 5.00 | KU-1.27 | PASS (Mikado theory) |
| Emergent ⟨z⟩ (10 seeds) | **2.574 ± 0.015** | predicted 2.60 (10-seed mean within prediction) | KU-1.3 sub-isostatic | PASS ([2.5, 3.9] band) |
| Biology gap (⟨z⟩ − KU-1.3 ref 3.2) | **−0.62** (logged) | reflects absence of 3D covalent / bundle XLs in 2D Mikado | KU-1.3 (info) | logged, not a fail |
| Nematic S at target 0.5 | 0.497 ± 0.016 | within 1/√(2N) = 0.016 | KU-1.9 | PASS |
| Backbone+XL FD vs analytic | **max rel err 1.9e-7** | < 1e-6 | KU-1.24, KU-1.28 | PASS |
| XL Newton 3rd law (Σ F_xl) | exact 0 | — | KU-1.28 | PASS |
| Sanity Gate (5 checks) | **5/5 PASS** | — | CLAUDE.md | PASS |

## Test suite

```
acs_kb/tests/test_ecm.py — 7 passed in 0.46 s
  ✓ test_generation_contract
  ✓ test_nematic_order_aligned
  ✓ test_emergent_z_matches_mikado_theory             # ⟨z⟩ within 10% of theory
  ✓ test_emergent_segment_length_matches_target       # ℓ_c within 10% of target
  ✓ test_analytic_forces_match_finite_difference      # WLC + XL combined
  ✓ test_xl_forces_pairwise_newton_third_law
  ✓ test_sanity_gate_passes_on_resolved_defaults      # incl. biology-gap log
```

## Performance

For `n_fibers = 3142` (15710 beads, 7687 cross-links):

| Step | Wall time |
|---|---|
| Fiber generation | 0.88 ms |
| Segment-intersection XL generation (O(N²) AABB+Cramer) | **102 ms** |
| Energy + forces (WLC + XL) | **8.2 ms** |
| **Total per-call** | **≈ 111 ms** |

Under the brief's `< 1 s` budget by ~9× — still ample headroom for Unit 1.2. The dominant cost is intersection finding (one-time at network construction). If Unit 1.2 needs faster reconstruction (e.g. plastic ECM remodelling), a grid-bucket / KD-tree variant of `_segment_intersections` would drop this to O(N).

## Sanity Gate output

```
Sanity Gate · phase1_ecm_network
  PASS  mikado_segment_length_in_KU_range  [KU-1.27 / KU-1.7]
        emergent ℓ_c=2.04e-06 m, predicted 2.00e-06, target ξ=2.00e-06, accepted [1.0e-06, 4.0e-06]
  PASS  z_in_sub_isostatic_range           [KU-1.3]
        emergent ⟨z⟩=2.579 (predicted 2.600, KU-1.3 ref 3.20), accepted [2.5, 3.9], z_c=2d=4
  PASS  z_within_prediction                [construction self-consistency]
        |⟨z⟩−predicted|=0.022, tol=10% (probes the periodic-boundary loss)
  PASS  fiber_count_under_cost_cap         [cost]
        n_fibers=3142, cap=20000; n_beads=15710; n_links=7687
  PASS  biology_gap_logged                 [info]
        ⟨z⟩_measured − ⟨z⟩_KU-1.3 = -0.621; 2D Mikado lacks 3D covalent and bundle cross-links present in collagen
```

## Why ⟨z⟩ ≠ KU-1.3 biological reference

The Mikado theory says crossings per fiber = 2 ρ_L L / π. With ℓ_c = ξ = 2 μm (i.e. ρ_L = π/4 μm⁻¹, KU-1.27 derivation), this gives 5 crossings per fiber. Each crossing creates one XL touching one bead on each of two fibers, so ⟨z⟩_XL = crossings/beads = 5/5 = 1. Total ⟨z⟩ = backbone (1.6) + 1.0 = **2.6**.

To hit KU-1.3's biological ⟨z⟩ = 3.2 with the **same** Mikado, we'd need ℓ_c ≈ 1.25 μm (≈ 5000 fibers), which is denser than the experimental porosity ξ ≈ 2 μm. The two KU values are *not simultaneously satisfiable* in 2D Mikado:

- Real 3D collagen has covalent cross-links (lysyl oxidase) and bundled-fiber connectivity not captured by 2D geometric crossings.
- We chose the Storm-MacKintosh modeling convention: match ξ first, accept the emergent ⟨z⟩ ≈ 2.6, document the gap.

The Sanity Gate's biology-gap line logs this explicitly. Phase 2 may close the gap either by promoting the network to 3D or by adding a chemistry-driven XL augmentation; both are deferred.

## Out of scope (per brief stop point)

- No time integration (Unit 1.2 — KU-1.26).
- No cells / FA / junctions (Units 2–4).
- No Bell-Evans XL dynamics (KU-1.15 deferred to Phase 2).
- No periodic-image intersection enumeration (Phase 2; accounts for ~1 % of the predicted-vs-measured ⟨z⟩ gap, the dominant of two error sources alongside Poisson noise).

## Suggested next step — Unit 1.2

1. Implement `EulerMaruyama` and `BAOAB` integrators behind an `Integrator` protocol (KU-1.26).
2. CFL Sanity Gate: `dt < γ_b/k_xl` (XL stiffness) and `dt < γ_b ℓ₀/(μ + κ/ℓ₀²)` (WLC).
3. Equipartition test: ⟨½ k_xl |Δr|²⟩ → ½ k_B T per cross-link.
4. Strain-stiffening validation (KU-1.4 + KU-1.30 #1, #2): Lees-Edwards shear ramp, fit G(γ) ∼ G₀/(1−γ/γ_c)², expect G₀ ≈ 36 Pa, γ_c ≈ 0.16.

**Awaiting PI review before proceeding to Unit 1.2.**
