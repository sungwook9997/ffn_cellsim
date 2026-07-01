# FF Stage 6Y — filament arrangement study: the γ-floor is arrangement-robust; structure is not

**Date:** 2026-07-02  **Engine:** FF (Warp, A5000)  **Branch:** dcm/main
**PI directive:** "다시 테스트해봐야지 … 크로스링크도 다시 검토 … 필라멘트 정렬이 어떻게 되어있는지가 제일 중요 …
여러가지 배열에 대해서 테스트해봐" — the filament ARRANGEMENT (alignment) is the key variable; re-review the
crosslinks; test across arrangements. **Follows** 6W/6X (whole-cell compartments).

## One line

The FF cortex had only ever been tested at ONE filament arrangement (isotropic random). Adding an alignment
control and sweeping arrangements shows: **arrangement is a real first-order variable for cortex STRUCTURE
(connectivity, structural tension) but the reported γ-FLOOR metric (γ_myo) is arrangement-robust** — the
γ-floor conclusion is not an artifact of one arrangement.

## What was built (additive; default = current isotropic, bit-identical → no regression)

`cortex_assembly`: an `orientation` ∈ {isotropic, circumferential, meridional} + `nematic_S` ∈ [0,1] control
on the fiber tangent field (`_tangent_field`), plus `nematic_order()` to measure alignment. Threaded through
`build_crosslinked_cortex` and `gamma_floor_run`. isotropic = the historical random tangent (test asserts
bit-identical). circumferential = tangent ∥ ê_φ (filaments wrap like latitude lines); meridional = tangent ∥
ê_θ (pole-to-pole).

## Crosslink re-review — the 1:1 default sits at the percolation edge

Connectivity = giant component of (intra-filament backbone ∪ crosslinks). At the production **1:1** crosslink
density (n_xl = n_fil), the network is **marginally connected AND strongly arrangement-sensitive**:

| arrangement | align→director | node-giant % (1:1) |
|---|---|---|
| isotropic | 0.64 | **78** |
| circumferential (aligned) | 1.00 | **39** |
| meridional (aligned) | 1.00 | 68 |

Crosslink-density sweep (giant %): 1:1 → {78/39/68}, **2:1 → {98/97/97}**, 4:1 → ~100 all. So ≥2:1 heals every
arrangement to a well-connected gel; 1:1 is right at the edge (only ~25% of nodes crosslinked, degree ~1.15).
→ **Surface to PI: is the 1:1 (γ-irrelevant, FF_STAGE6M) crosslink density the right connectivity choice, or
should production use ≥2:1?** (Real cortex is a well-connected gel.)

## γ across arrangements (native N=70686, f_myo=5 pN, 3 seeds, mean±std)

| condition | γ_myo | γ_active | γ_actin |
|---|---|---|---|
| isotropic 1:1 | 0.076±0.011 | 0.163±0.011 | 0.183±0.018 |
| circumferential 1:1 | 0.075±0.008 | **0.109±0.008** | **0.134±0.016** |
| isotropic 2:1 | 0.086±0.006 | 0.159±0.014 | 0.188±0.036 |
| circumferential 2:1 | 0.071±0.010 | 0.129±0.024 | 0.142±0.023 |

- **γ_myo (the reported γ-floor metric) is arrangement-ROBUST** — 0.07–0.09 pN/µm across isotropic/aligned/
  density, all differences within ~1–2σ. Since γ_myo ∝ myosin-dipole count × geometry, it does not depend on
  how the filaments are aligned → **the γ-floor conclusion (myosin ~10³× under the tension band, 6Q/6V) holds
  regardless of arrangement.** This is the key reassuring result.
- **γ_active and γ_actin (the structural frustration residual) ARE arrangement-sensitive** — aligned cortices
  run **~30% lower** than isotropic (γ_active 0.11 vs 0.16; γ_actin 0.13 vs 0.18 at 1:1), a real effect
  (~3–4σ over seeds). Physically: parallel/aligned filaments have less geodesic frustration on the curved
  shell than randomly-crossing isotropic ones. This is intrinsic (persists at 2:1), not just a connectivity
  artifact.

## Conclusions (for PI)

1. **Arrangement is a genuine first-order variable** for cortex structure — the PI intuition is confirmed.
   Fully-aligned filaments cross less → connectivity 78%→39% and structural tension ~30% lower.
2. **The 1:1 crosslink density is at the percolation edge** and is what makes the structure arrangement-
   sensitive; ≥2:1 gives a robust gel across all arrangements. → PI decision on the production density.
3. **The γ-floor result is arrangement-robust** (γ_myo invariant) → 6V/6Q stands; we were not fooled by
   testing one arrangement. The arrangement-sensitive quantity (γ_active) is the contaminated structural
   channel, not the floor metric.
4. Real cortex arrangement (isotropic vs partially-aligned) is a physiological input; if a specific MCF7
   cortex order parameter is known it should be a lit-anchored config value (surface to PI) rather than the
   implicit isotropic default.

## Sanity gate
- Additive, no regression: default isotropic bit-identical to the historical builder (test). ✓
- No magic number: arrangement is a modeling INPUT swept as a controlled variable (not tuned to a target). ✓
- Verified: arrangement γ effect confirmed over 3 seeds (3–4σ); connectivity is a direct graph computation. ✓
- Honest: γ_myo-robust vs γ_active-sensitive distinction stated; 1:1-vs-2:1 density surfaced to PI. ✓

## Figures
- `outputs/ff/figs/arrangement_study.png` — (top) 3D views of the three arrangements (isotropic tangle /
  circumferential latitude-wrap / meridional pole-to-pole); (bottom) connectivity vs crosslink density,
  γ channels across arrangements (seed error bars), and the conclusions.

## Files
- `ff/cortex_assembly.py` — `_tangent_field`, `_local_tangent_frame`, `nematic_order`, orientation/nematic_S
  in `build_cortex_network`.
- `ff/gamma_floor.py` — orientation pass-through in `build_crosslinked_cortex` + `gamma_floor_run`.
- `tests/ff/test_cortex_assembly.py` — orientation no-regression (bit-identical) + alignment order parameter.

Related: [[project-cell-mechanics-extend]], [[project-gamma-floor-likely-deficit]], FF_STAGE6Q (γ_myo floor),
FF_STAGE6M (1:1 crosslink density), FF_STAGE6V (tension pressure-borne).
