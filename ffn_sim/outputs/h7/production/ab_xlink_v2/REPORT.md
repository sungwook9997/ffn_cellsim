# ab_xlink_v2 — crosslink transmission redesign (Slater 2021), CPU/GPU preliminary

**Date:** 2026-06-11  **Host:** gbook A5000.  **Scale:** reduced smoke (n_filaments=60,
n_cortex_actin=420, 100 motors, warmup 800, contract 8000, seed 1, connected mesh, grip_walk).
**Status:** PRELIMINARY — one axis VALID, one axis UNTESTED (harness gap found & flagged).

## Why
v1 ab_xlink varied crosslink **stiffness** (k_intra/k_attach 1e-3 vs 1e-7) → γ ~5% → "transmission
not the wall". Slater 2021 (Soft Matter 17:10274, same Taeyoon-Kim discrete family) says the actual
stress-transmission knobs are crosslink **density + force-dependent unbinding (Bell-Evans)**, not
stiffness. This re-tests on those knobs (brief: `docs/AB_XLINK_TRANSMISSION_REDESIGN_2026-06-11.md`).

## Results (γ_soft, ×1e3 mN/m; engaged 256/5600 heads, fixed across all points)

**Bell off-rate axis (VALID — `alpha_k_off0` is a per-crosslinker kinetic, count-independent):**
| k_off0 [1/s] | g_soft total | g_soft actin-bonds (transmission) | g_soft myosin-bonds |
|---|---|---|---|
| 0.066 (anchor, Ferrer2008) | 1.7445 | 0.8985 | 0.1626 |
| 0.0066 (10× slower) | 1.7565 (+0.7%) | **0.9389 (+4.5%)** | 0.1494 |
→ slower unbinding (more persistent connectivity) **does raise the transmission channel
(g_actin +4.5%)** — qualitatively consistent with Slater — but the **total γ moves only +0.7%**,
second-order vs the engagement floor.

**Density axis (NOT TESTED — harness gap):** sweeping `--xl-n` over **10 → 30 → 100 → 250 → 1000
→ 4000 (400×) gave γ_soft EXACTLY 1.7445 every time.** This is **not** density-independence — under
`connected_mesh=True` (this operating point), **`n_xl` is a SEEDING OUTPUT, not an input**
(`cell.py:745`, `connected_mesh.py:102`; `cell.py:931` uses `2·_cm_seed.n_xl`), so the `--xl-n`
config override is **silently ignored**. The "invariance" is the ignored override, not a result.
→ density remains **UNTESTED**. (Fix: vary connected-mesh structural density via `--cm-z-struct`,
or run `connected_mesh=False` with the dynamic `n_xl` path. `--xl-n` help now warns of this.)

## Verdict
- **Off-rate**: transmission channel is weakly off-rate-sensitive (+4.5% g_actin per 10× slower),
  but total γ barely moves (+0.7%) → does NOT challenge the magnitude floor.
- **Density**: OPEN — not yet validly tested (override ignored under connected_mesh).
- **Magnitude floor unchanged**: still **engagement/generation-limited** (256 of 5600 heads
  engaged; band closure needs ~26000) — loop24b "generation-limited, not transmission" narrative
  is **not overturned** by the valid (off-rate) axis. The Slater density hypothesis is neither
  confirmed nor refuted here.

## Next (PI-scoped)
1. **Working density knob**: sweep `--cm-z-struct` (connected-mesh structural connection density)
   and/or `connected_mesh=False` + dynamic `n_xl`, to actually vary crosslink density.
2. **Slater stress(r) protocol**: add the radial surface-crossing stress estimator (2nd, independent
   of γ_soft) per the brief §3.2.
3. Full-scale (n_filaments=120, longer contraction) confirmation once the density knob works.
This preliminary is reduced-scale + 1 seed; absolute γ (1.74e-3) is under-engaged vs the v1
plateau (1.37e-2) — trends only, not band evaluation.
