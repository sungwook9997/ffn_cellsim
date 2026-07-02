# DCM spheroid compaction — aggregate Foty-Steinberg σ is the driver, but native needs true IPC (2026-07-02)

> **⚠️ CORRECTION (native N=400):** the "SOLVED" headline (from N=100) was PREMATURE. At native N=400 the
> aggregate σ run **blows up** — V/V0=0.39 (cells collapse), pen=123 (massive interpenetration), maxZ explodes.
> The Rg −19% at N=400 is partly cell-collapse/overlap, NOT clean densification (porosity "0.96" = the collapse
> signature). **Verdict:** the aggregate σ IS the correct global compaction DRIVER (Rg drops, the mechanism the
> junction levers lacked), BUT it is not clean without a proper contact method — exactly §2e's conclusion that
> it needs **aggregate σ AND true IPC (log-barrier), built together.** The penalty contact fails under the
> compaction load at native density (pen 3.8 @N=100 → 123 @N=400). Next: build true log-barrier IPC. The N=100
> result below is real (V/V0 stable there) but is NOT a scale-valid claim on its own.

## (original N=100 finding — driver works, but see the correction above)


**Date:** 2026-07-02  **Engine:** DCM (Warp, A5000)  **Branch:** dcm/main
**PI directive (8h goal):** "dcm에서 현재 spheroid 조립 관련하여서 문제가 생기는 건데, 그것 해결할 방법 절대
멈추지 말고 계속해서 테스트해." **Follows** SPHEROID_STAGED_FORCES §2e (dynamic loose→compact FAILED with
junction levers; surfaced-to-PI). Executing my recommendation: build + test the missing aggregate σ.

## Result — dynamic loose→compact compaction now works

The missing driver (§2e: "needs AGGREGATE-level surface tension + IPC, built together") is the **aggregate-
level Foty-Steinberg liquid-drop surface tension**. Built as `dcm_aggregate_tension_warp.aggregate_laplace_
kernel`: inward Laplace ΔP=2σ/R_agg on the FREE (media-exposed) envelope faces toward the **global** aggregate
centroid — contracts a loose aggregate (which per-cell/differential γ, being local, cannot).

**Loose start (voronoi gap 2.4, N=100, bundle-10 stable, reach+contraction on), 8000 steps:**

| | ΔRg | porosity | asph |
|---|---|---|---|
| **baseline** (reach + contraction, no aggregate) | **+0.9%** | 0.75 → **0.85** (rises) | 0.027 → 0.029 |
| **+ aggregate σ = 5 mN/m** | **−25.6%** | 0.75 → **0.61** (drops) | 0.027 → **0.050** |

- The **baseline reproduces §2e's failure exactly**: Rg flat (+0.9%), porosity RISES (voids grow) — reach +
  contraction pull cells locally but do NOT globally densify.
- **Adding aggregate σ reverses it**: Rg drops 26%, porosity drops (void elimination), asphericity rises
  (faceting) — genuine dynamic loose→compact compaction, the first in the whole DCM series. This confirms the
  §2e diagnosis and supplies the missing global densification driver.

## Key debugging (all native/GPU, gbook) — the path to a stable test

1. **`--cad-rbind` is in METERS** (c_adh internal unit), not µm — passing 3.5 gave r_bind=3.5e6 µm =
   infinite-reach bond explosion → OOM. Fixed (×1e-6 in the driver).
2. **`--ipc` is broken in the `--conservative-contact` combo** — pen=73 (vs 3.8 without). Use
   `--conservative-contact --integrator implicit` (pen~3.8 = the known bundle level; true-IPC is a separate
   contact-fidelity item).
3. **Cadherin bundle-40 (~7 nN) DESTABILISES the loose start** (V/V0 diverges 0.4↔7.5 = numerical blowup).
   **bundle-10 (~1.7 nN, lit-range KB-4.11) is STABLE** (V/V0=0.995). The compaction test uses bundle-10.
4. σ is the lit Foty-Steinberg tissue tension (5 mN/m ∈ 1–20 band), swept as a controlled variable, NOT tuned
   to a compaction target.

## Honest caveats
- **Interpenetration** pen~3.8 (>0.3 gate) persists — the pre-existing strong-bond contact-fidelity issue
  (needs true log-barrier IPC; the `--ipc` flag is currently broken). The compaction (Rg/porosity) is measured
  despite this; a true-IPC fix is the follow-up for clean contact.
- Convex-hull porosity (1−ΣVcell/Vhull) is high in absolute terms (~0.6–0.85) because the hull is convex; the
  **trend** (aggregate drops vs baseline rises) is the robust comparison, not the absolute value.
- N=100 shown here; **native N=400 confirmation is running** (same config).

## Sanity gate
- No magic number: σ = lit Foty tissue tension (swept, not tuned); bundle-10 = lit KB-4.11 nN-range. ✓
- Physiological: aggregate σ is the emergent tissue surface tension (Foty-Steinberg DAH); mechanistic (envelope
  Laplace), not a lumped attractor. ✓
- Honest: baseline reproduces the known failure; the pen + porosity-absolute caveats are stated. ✓

## Deliverables (data + viz + HTML)
- **Data:** above table (compute_compaction.py: Rg/porosity/asph from frames; DCM pos in metres ×1e6→µm).
- **Viz:** `outputs/dcm_compaction_test.png` — Rg / porosity / asph trajectories, aggregate vs baseline.
- **Interactive HTML:** `outputs/dcm_compaction/dcm_agg_viewer.html` (+ `dcm_base_viewer.html`) — three.js
  animated playback (33 frames): rotate/zoom + play to WATCH the loose aggregate compact; cross-section modes.

## Files
- `dcm/dcm_aggregate_tension_warp.py` — `aggregate_laplace_kernel` + `aggregate_centroid_radius`.
- `dcm/dcm_warp_decohesion.py` — `--aggregate-tension --sigma-agg --agg-every` + `--cad-rbind` (µm→m fix).

Related: SPHEROID_STAGED_FORCES §2e (the diagnosis), CONFLUENT_INITIALIZER (the static shortcut), memory
[[project-dcm-compaction-turgor-blocked]], [[project-dcm-faceting-confluent-init]].
