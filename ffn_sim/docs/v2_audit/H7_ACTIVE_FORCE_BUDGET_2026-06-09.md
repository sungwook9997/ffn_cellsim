# H.7 active-myosin force-budget audit — the active-γ floor is GENERATION-bound at the literature parameters

**Date:** 2026-06-09 (autonomous). Branch `h7/full-cell-integration`. Author: Lead.
Frontier from `H7_GATE_B_RESULT_2026-06-09.md` §5 PRIMARY: Gate-A REFUTE (γ_soft≈3.06e-3
mN/m, 114× under band), Gate-B buckling REFUTE, active floor persists on the connected
mesh. The §5 question: *why does myosin contribute ~0 to cortical tension (active stress
~400× < turgor)?* Candidates: (1) force/density too low, (2) turgor too high, (3)
local→spanning conversion failing.

Tool: `ffn_sim/scripts/h7_active_force_budget.py`. It builds the **same** Gate-B operating
point (suspended/rounded, FA OFF, turgor ON, connected mesh, grip_walk; identical harness
to `h7_gate_b_probe._build_settled_cell`) and decomposes the floor in one measurement,
reading the live updater + the snapshot the integrator sees (so γ_soft == the gate's γ_soft).

## The decisive number — the analytic parameter envelope

The active-gel surface tension of an isotropic minifilament population on a thin shell is
the dipole 2D stress

> γ_active ≈ ½ · n_2D · f_minifil · ℓ_minifil

at **full engagement + full stall** — the *upper envelope* of what the configured
parameters can produce, independent of the sim. With the live config values (no tuning):

| quantity | value | provenance |
|---|---|---|
| n_2D (native areal density) | 0.60 /µm² | Nie 2015 Cytoskeleton, **HeLa interphase** (the only proxy; NO MCF7 datum; LOW-MED confidence) |
| f_minifil (2·28·F_stall) | 112 pN | Billington 2013 (56 heads) × Chugh 2017 (2 pN/head) |
| ℓ_minifil (backbone) | 301 nm | Billington 2013 EM contour |
| **γ_active envelope** | **0.0101 mN/m** | = ½·n2D·f·ℓ |
| Hosseini band_lo | 0.18 mN/m | contract §7 (MCF7 interphase IQR) |
| **gap** | **17.8× UNDER band_lo** | |

**Even at the theoretical maximum (every native minifilament engaged, every head at full
stall, perfectly isotropic), the configured parameters predict γ_active ≈ 0.01 mN/m — ~18×
below the band floor.** No engagement-throughput, transmission, buckling, or mesh fix can
cross this ceiling; it is set by the force budget alone.

## The simulation sits below the envelope, consistently

Smoke audit (n_filaments=160, CPU, loading phase s_grip=0; full data
`outputs/h7/production/h7_force_budget_smoke.json`):

| channel | value |
|---|---|
| engagement | 80→452 heads (1.4→8.1%, still climbing) of 5600 |
| per-head delivered T | 0.90 pN (mean attach r 213 nm; F_stall 8.48 pN scaled) |
| generated Σ\|T\| | 0.41 nN |
| γ_soft (MOP, active) | 1.75e-4 mN/m |
| γ_IK (virial cross-check) | 6.3e-5 mN/m (same order — confirms MOP) |
| γ_passive (turgor YL) | 0.499 mN/m |

Stacking with the prior full-scale Gate-A contraction plateau (γ_soft = 3.06e-3 mN/m, s_grip
developed): the picture is fully coherent —

```
band_lo                       0.18    mN/m   (target)
analytic envelope             0.010   mN/m   (18× under)  ← PARAMETER ceiling, irreducible
Gate-A contraction plateau    0.0031  mN/m   (60× under)  ← engagement+transmission losses
loading phase (this smoke)    0.00017 mN/m  (1000× under) ← no contraction yet
```

The contraction (grip_walk stepping) lifts the loading floor ~16× toward the envelope; the
remaining gap below the envelope is engagement (~8%, climbing) + geometric transmission. But
the **dominant, irreducible gap is the 18× parameter envelope** — Gate-A/Gate-B were probing
levers (s_grip generation, buckling, mesh connectivity) that all live *below* a ceiling that
is itself far under band.

## ⮕ FINAL synthesis is `H7_ACTIVE_GAMMA_SYNTHESIS_2026-06-09.md`

The interpretation evolved across three loops (recorded below for the trail). **Net result:
the active-γ floor is dominated by GENERATION (myosin force-budget / density datum ~10–20×
too low); network amplification is small (~1.3×, crosslink-bound) and cannot rescue it.**
Loop-2's "amplification could be 10–100×" was over-speculative and is corrected by the data
(g_actin/g_myo ≈ 28%). See the synthesis doc for the PI-facing version.

## Interpretation — loops 2–3 trail (superseded by the synthesis above)

⚠️ **Loop-1 said "generation-bound, 18× irreducible." Loop-2 over-corrected toward
"transmission-dominant via 10–100× amplification." Loop-3 + the crosslink-spacing bound show
the truth is in between and closer to loop-1: generation/density is dominant, amplification is
small.** The bond-type-resolved γ split:

**Bond-resolved active γ (loading phase, smoke):**

| channel | γ (mN/m) | trend as engagement grows 80→357 heads (4.5×) |
|---|---|---|
| total g_soft | 1.5e-4 | rises with engagement |
| ├ **myosin bonds** (direct dipole, ℓ≈301 nm) | 1.6e-4 | **rises** 1.15e-4→1.62e-4 (tracks engagement) |
| └ **actin-network bonds** (ℓ≈500 nm–3 µm) | 4.6e-5 | **FLAT** 3.8e-5→4.6e-5 (myosin-insensitive) |

So the active tension lives almost entirely in the **direct myosin dipole bonds**; it does
**NOT propagate into the actin network** (the actin-bond γ does not track myosin engagement).
This is decisive because a real cortex reaches band via **network prestress amplification**:
a local myosin contraction transmitted around a long actin/crosslink load path produces
tension ≈ f·L_path, and L_path/ℓ_minifilament can be ~10–100×. The MOP/virial of a *loaded*
actin network would far exceed the direct dipole. **That amplification is absent here** — the
contraction is shunted into the local myosin springs (Gate-A: actin backbone r/r0=1.0000;
γ_rigid also myosin-insensitive 0.078 vs 0.082).

The two walls:

- **WALL A — propagation/transmission (dominant, fixable lever).** Myosin tension does not
  load the actin network into spanning hoop tension. With prestress amplification the direct
  dipole (~0.01) could reach band; without it, γ stays at the direct-dipole scale. This is
  Gate-A's "transmission/lever" verdict, now pinned to *the actin network not carrying
  myosin-induced tension*. Candidate causes: isotropic local force-balance cancellation;
  the M-SHAKE-rigid actin not deforming to build prestress; weak myosin→actin geometric
  coupling.
- **WALL B — direct-dipole generation is itself sub-band (parameter).** Even the full-engage
  +full-stall direct-dipole envelope (½·n2D·f·ℓ, ℓ=minifilament) = 0.01 mN/m, 18× under band.
  So even if WALL A were removed at the *local* scale, the *direct* generation is short; band
  needs the longer-ℓ network path (WALL A) AND/OR a higher force budget (density/stall datum).

Candidate resolution:

- **(1) force/density too low — contributes (WALL B), but not the whole story.**
- **(2) turgor too high — NO.** Turgor (γ_passive 0.5, γ_rigid carries it) is physiological
  (Π₀=133 Pa) and reported separately (B3/B4); orthogonal to the active channel.
- **(3) local→spanning conversion failing — YES, dominant (WALL A).** The actin network does
  not carry myosin-induced tension; this is the missing prestress amplification.

**Caveat (must confirm):** the bond-resolved split above is the LOADING phase (s_grip=0). The
definitive test is whether the actin-network γ *rises as s_grip develops* (active stepping) —
a contraction-developed run is queued (gbook GPU). If actin-network γ stays flat through
s_grip→0.5, WALL A is confirmed as the dominant, mechanistically-located transmission wall.

## What this is NOT (discipline)

- NOT a license to bump the density/stall to pass. The band is LOCKED (contract) and the
  magic-number rule forbids tuning a parameter to land a gate. **This is surfaced to PI as a
  parameter/datum question, not acted on.**
- The analytic formula is a standard active-gel dipole estimate; it may miss a cooperative or
  effective-length factor (the contractile element's arm may be the actin mesh size, not the
  301 nm minifilament). That is the open reconciliation item below.

## Open / next (PI-surface candidates)

1. **MCF7 cortical myosin density.** The model uses a HeLa proxy (0.6/µm², explicitly LOW-MED
   confidence, no MCF7 datum). The band needs ~18×; is MCF7 NMII density/stall actually
   higher? → literature search (next loop) for an MCF7-anchored density or a direct
   active-stress measurement. If a real higher MCF7 value exists, that is a datum correction
   (PI-gated), not tuning.
2. **Active-gel formula effective length.** Reconcile ½·n2D·f·ℓ against the canonical
   Salbreux/Hannezo/Turlier myosin→cortical-tension relation; the effective contractile arm
   ℓ may be the actin filament/mesh length (µm), not the minifilament (301 nm) — a ~3–10×
   factor that would shrink the parameter gap.
3. **Band decomposition.** Does Hosseini 0.27 conflate passive (cortex elasticity + membrane)
   with active? The blebbistatin-sensitive (active) fraction sets the true γ_active target.

## Artifacts
- `ffn_sim/scripts/h7_active_force_budget.py` (reusable; `--device gpu` for full ×40 + a
  contraction-developed run).
- `ffn_sim/outputs/h7/production/h7_force_budget_smoke.json`,
  `ffn_sim/outputs/h7/figs/h7_force_budget_smoke.png`.
