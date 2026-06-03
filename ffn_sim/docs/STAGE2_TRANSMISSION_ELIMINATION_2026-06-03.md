# STAGE-2 — Transmission/Turnover/Coherence Eliminated; γ-Floor Traces Upstream

**Date:** 2026-06-03 (autonomous block, PI asleep — git as safety net)
**Branch:** `test/mcf7-fullcell` (worktree)
**Owner:** Lead session
**Status:** STAGE-2 investigation closed with a redirect. Needs PI nod for the
proposed next diagnostics (one touches the stiffen-forbidden crosslinker `k`).

---

## TL;DR

The active cortical-tension floor (`g_soft ≈ 2e-4 mN/m`, ~1000× under the KU-3.5
band `[0.35, 0.65]`) is **NOT** a transmission problem. Five independent levers —
crosslink **percolation**, actin **turnover**, measurement **channel**, global
**coherence**, and **buckling** — were each tested in the full MCF7 constrained cell
and **none** move `g_soft` off the floor. The floor therefore traces **upstream** of
transmission, to force **generation / soft-coupling / force-budget** (STAGE-1
territory or model fidelity), not to the cortex network's connectivity or the
random-dipole √N cancellation we had hypothesised. The **structural/passive**
tension (`g_rigid`, native 0.57 mN/m, in-band) is real but myosin-independent —
not the KU-3.5 active target.

All runs: MCF7 full cell (cortex+membrane+nucleus+cytoplasm+enclosed-vol),
constrained M-SHAKE backbone, grip_walk motors, mesoscale `n_fil=1000` percolated
via the bind-scale recipe (or native `n_fil=38000`), CPU.

---

## The eliminations

| # | Hypothesis (γ-floor stack) | Test | Result |
|---|---|---|---|
| ④a | transmission via **percolation** | crosslink coordination z swept 0→3 (Ennomani gate) via `--bind-scale 6 --n-xl 1.5·n_fil --kon-scale 300` | `g_soft` floor **unchanged** (native z≈2 AND mesoscale z≈3 both ~2e-4) |
| ④b | **turnover** (sustained-stress, Hiraiwa-Salbreux) | implemented turnover-in-constrained (chain split, f_ss=0.591) and ran ON vs OFF | `g_soft` ON ≤ OFF at matched samples — **no lift** |
| — | measurement **channel** misattribution | motors ON (n=100) vs OFF (n=0) | `g_soft` ×340 motor-responsive (= the true active channel); `g_rigid` motor-**independent** (~0.09 both) = purely structural |
| ⑤ | global **coherence** (√N random-dipole cancellation) | `FFN_MYOSIN_ALIGN=meridional` — impose a perfectly coherent director (verified \|axis·meridian\|=1.000) | `g_soft` 2.32e-4 vs random 2.18e-4 — **no jump** |
| ⑥ | **buckling**-mediated contractility (Murrell-Gardel/Lenz) | cortex-angle θ readout under active load | θ≈179° (thermal, straight); **buckled<150° = 0.0%** — filaments do not buckle |

The coherence negative is the decisive surprise: a perfectly aligned bipolar
field still does not raise the active tension. So the force the motors produce is
not being converted into network tension **regardless of how it is arranged** —
the limit is at conversion/generation, not arrangement.

## What this isolates

Every rearrangement (more crosslinks, turnover flow, coherent alignment) funnels
the motor force through the same soft inter-filament coupling (crosslinker
`k=1e-7 N/m`, Furuike 2001 — a literature constant, **stiffening forbidden** by
the Magic-Number / no-tune hard rules) and the same grip-walk head springs. The
floor is set **before** transmission:

- the motors transport but do not build **sustained head-spring tension** (the
  active soft-bond tension `g_soft` saturates at ~2e-4 and the rise is a slow
  grip-walk buildup that extrapolates to ~40 M steps to reach band — i.e. a real
  floor, not a slow approach);
- inter-filament load that *is* produced is carried by soft springs that yield.

This is **STAGE-1 (force generation/ceiling) / model-fidelity**, not STAGE-2.

## What is NOT broken (positive results)

- **Structural/passive cortical tension reaches band.** `g_rigid` (rigid-bond
  Lagrange tension = construction prestress + enclosed-volume pressure +
  membrane) is **0.57 mN/m native (in band)**, myosin-independent. Real passive
  cortical tension — credit it as such, not as KU-3.5 active.
- **The mesoscale percolation recipe** (`bind×6` reach + `n_xl/n_fil=1.5` →
  z=2.96) is a reusable fast STAGE-2 tool: it reproduces the native Ennomani
  z-gate at `n_fil=1000` in ~12 min instead of a 3-4 h native run, by scaling the
  crosslinker partner-search radius (NOT `k`) — the geometric dual of the
  ratified ×40 areal coarse-graining.
- **Turnover-in-constrained is implemented and validated** (commit 7a4bb49):
  `ConstrainedLeimkuhlerMatthewsBAOAB.resync_chains()` (additive, default
  bit-identical) + `ActinTurnoverUpdater(baoab_action=)` chain-split, so a cofilin
  sever actually releases the M-SHAKE stretch constraint. cf→f_ss=0.591 verified.
  Ready for when turnover IS the lever; CPU-only (ragged) until GPU fragment
  padding.

## Proposed next diagnostics (need PI nod)

1. **Soft-coupling-ceiling probe** — sweep crosslinker `k × {10,100,1000}` as a
   DIAGNOSTIC ONLY (not adopted; same status as the meridional alignment probe).
   If `g_soft` rises with `k`, the soft-coupling ceiling is confirmed.
   ⚠️ Crosslinker stiffening is specifically forbidden by the hard rules, so even
   a labelled diagnostic is surfaced to PI rather than run unilaterally.
2. **Grip-walk head-tension instrumentation** — measure the head-spring extension
   distribution directly. Are the heads at rest (transport-only) or loaded but
   small? Decides generation-ceiling vs force-budget-magnitude.
3. ~~Buckling diagnostic~~ **DONE (negative)** — θ readout shows filaments stay
   straight (θ≈179°, buckled<150°=0.0%); the Euler threshold F_B≈6.9pN per 7-bead
   segment exceeds the distributed per-filament motor compression. Buckling
   contractility is not realized — a sixth eliminated lever.
4. **Force-budget re-examination** — is the per-motor stall force (after
   mesoscale_force_scaling) large enough to localize compression / build head
   tension? The 4b force-budget gap may be under-compensated.
5. **Model fidelity** — does the model miss a load path (filament-filament steric
   load-bearing under excluded volume; native filament density; a non-spring
   transmission element)? A focused question for the mechanism audit.

## Artifacts

- Driver flags added (`mcf7_fullcell_stage1.py`): `--bind-scale`, `--turnover-tau`.
- Diagnostics: `stage2_diagnostics.py` (`--bind-sweep`, `--bind-scale`).
- Figures: `outputs/h3/figs/fig_stage2_{percol_recipe,percolation_falsified,turnover_vs_off,conclusion}.png`.
- Logs: `outputs/h3/production/mcf7_{fullcell_NATIVE_percol,meso_PERCOL_bind6,meso_PERCOL_TURNOVER_t01,meso_PERCOL_MOTORS_OFF,meso_PERCOL_MERIDIONAL}.log`.
- Commits: `dadf50a` (bind_scale diag), `22fd5a2` (driver bind_scale), `7a4bb49`
  (turnover-in-constrained).
