---
kb_record:
  doc_id: FF_STAGE6O_NATIVE_GAMMA_2026-07-01
  title: ×40 mesoscale RETIRED — native ~38000-filament cortex; γ-floor magnitude updated (still floored)
  authoritative_as_of: 2026-07-01
  supersedes: []
  updates: [CORTICAL_TENSION_RECORD_2026-06-30, FF_STAGE6D_GAMMA_FLOOR_2026-06-29]
  note: "Updates the γ-floor MAGNITUDE to native scale (1.4e-4→6.2e-4 mN/m); the conclusion (force-magnitude floor) in the authoritative CORTICAL_TENSION_RECORD is unchanged + reflected there."
  status: current
---

# FF Stage 6O — native-scale cortex (×40 retired); γ-floor magnitude updated

**Date:** 2026-07-01 · **Engine:** `ffn_sim/ff/` · branch `dcm/main` · **PI directive: all work native**

## Why

The ×40 mesoscale (cortex N=1000 vs native ~38000) was ratified 2026-05-19 as a **CPU/HOOMD hardware
constraint**. FF is now GPU-native (the A5000 runs native ~38000 in ~5 s/realization, 236 MiB of 16 GB
— 1.4 % of memory). PI (2026-07-01): switch ALL work to native going forward. An A5000 N-convergence
test then showed the switch is not cosmetic — **γ_active is N-DEPENDENT**, and the ×40 was
**under-reporting it ~5×**.

## N-convergence (A5000, cortex γ-floor)

| N_fil | crosslink reach √(A/n) | γ_active [mN/m] | floor | settle | measure |
|---|---|---|---|---|---|
| 1000 (×40) | 1121 nm | 1.38e-4 | 2537× | 0.1 s | 0.0 s |
| 4000 | 560 nm | 2.96e-4 | 1181× | 0.0 s | 0.1 s |
| 10000 | 354 nm | 3.58e-4 | 977× | 0.1 s | 0.3 s |
| **38000 (native)** | 182 nm | **6.6e-4** | **~530×** | 0.3 s | 1.1 s |

γ rises with N and **converges at native** (settle 600 vs 2000 steps give identical γ). The ×40
crosslinker reach (1121 nm) was also wildly unphysical (α-actinin ε = 60 nm); native's 182 nm is far
closer. So native is more physical on both counts.

## Native γ-floor — re-verified, STILL floored, robust on every axis (A5000)

Production ensemble at native (N=38000 / n_xl=38000 / n_myo=3800, ratios preserved to isolate the
resolution change):

> **γ_active = 6.2 ± 0.13 × 10⁻⁴ mN/m, floor ~530–570× under the Salbreux band.**

- **Crosslink stiffness (6M) at native:** γ_myo is **bit-identical** at the corrected α-actinin
  stiffness (4.6e5) and the broken 0.1 (= 0.8366 pN/µm) — robust to 10⁶× crosslink stiffness, at native.
- **Motor turnover (6k) at native:** engaged fraction self-limits to 0.983 (= k_on/(k_on+p_off));
  floor holds.
- **Buckling / connectivity / finite-extensibility (6H):** N-independent contractile-network findings
  (screening, not amplification) — unchanged.

So the γ-floor conclusion **STANDS at native and is strengthened**: active actomyosin cortical tension
is force-magnitude-limited, robust to buckling/connectivity/extensibility/turnover/crosslink-stiffness.
The **magnitude is updated 1.4e-4 → ~6.2e-4 mN/m** (the ×40 under-reported it ~5×); the gap to band
narrows from ~2500× to ~530×, but the floor is real (not a coarse-graining artifact) — it persists at
the native, fully-resolved filament count.

## Open levers (unchanged by this switch)

The remaining γ-floor lever is the **experimental load-engaged motor-density datum** (Truong-Quang
penetration/overlap ≫ the only direct datum Nie 2015 0.625/µm²). NOTE: the production n_myo=3800
(ratio-preserved) gives a myosin areal density ~3 µm⁻² — this is a SEPARATE density question from the
resolution switch; reconciling it with Nie 2015 is the γ-floor's open lever, surfaced to PI.

## What changed (code)

`PROD_N_FIL/XL/MYO` → 38000/38000/3800; `CortexParams.n_filaments` → 38000; `CORTEX` ArchitectureSpec
→ 38000; `gamma_floor_production` gains n-overrides (default native, device path required at scale).
Tests use small-N fixtures (validate mechanism not scale; native magnitude validated on the A5000).
The other structures (filopodium/microvillus/SF cross-section/lamellipodium) were already at native
physiological counts — only the cortex was ×40-coarse-grained.
