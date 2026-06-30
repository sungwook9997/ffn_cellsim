---
kb_record:
  doc_id: FF_STAGE6N_CROSSLINK_TURNOVER_2026-06-30
  title: Robust stiff-crosslink relaxer (crosslink-turnover baseline) — lit-anchored link_k applied to production
  authoritative_as_of: 2026-06-30
  supersedes: []
  status: current
---

# FF Stage 6N — robust stiff-crosslink relaxer; the Ferrer link_k correction applied to production

**Date:** 2026-06-30 · **Engine:** `ffn_sim/ff/` · branch `dcm/main` · **PI-approved (둘 다)**

## Why

Stage 6M found that correcting the crosslinker stiffness to the lit-anchored Ferrer value (4.6e5/8.2e5
pN/µm) is physically right but, applied naively, makes the production cortex relaxation throttle (CFL)
and produce a SPURIOUS passive γ_xl (stiff crosslinks amplify the bending-settle residual stretch;
occasional over-stretch/blow-up). PI approved building the robust stiff-crosslink relaxer (item A) so
the correction can go into production. This is it.

## The robust relaxer — crosslink-turnover resting baseline

`gamma_floor.equilibrate(crosslink_turnover=True)` (default in `gamma_floor_run`):

1. Settle the resting SHELL **bending-only** (`relax_on_device(links=None)`) — the stiff crosslinks
   are EXCLUDED from the step, so they never throttle dt or over-stretch; the bending relaxes freely
   (E → 0, no blow-up, GPU-native).
2. Bind the crosslinks **FORCE-FREE at the settled geometry** (`xl_rest` ← settled lengths). This is
   the Hand §10.1 *attach* in the quasi-static / instantaneous-turnover limit: real crosslinks
   (α-actinin k_off ≈ 0.066/s, filamin catch-slip) rebind force-free at whatever the current cortex
   geometry is, so the RESTING baseline carries no pre-stress from the formation-vs-settle drift.

Mechanistic, not a tuning hack: the resting cortex has force-free crosslinks at ANY stiffness (with the
old soft 0.1 the drift stress was negligible anyway; the turnover baseline makes it exactly force-free,
which is what turnover physically does). The ACTIVE γ (myosin-transmitted, the production number) is
identical with or without it — myosin still stretches the (now force-free-rest) crosslinks under load,
and that response is force-magnitude-limited.

## Result — stiff crosslinks applied, production γ-floor clean + unchanged

`hand_kmc.ALPHA_ACTININ.link_k` = **4.6e5**, `FILAMIN.link_k` = **8.2e5 pN/µm** (Ferrer 2008,
PI-approved). Production γ-floor on the gbook A5000 (N=1000, n_real=12, stiff crosslinks +
crosslink_turnover, GPU-resident, 2.0 s):

> **γ_active = 1.459e-4 ± 0.21e-4 mN/m, floor 2399× under band** — clean (tight error bar, no blow-up),
> bit-consistent with the archived BAOAB-MD g_soft (~1.4e-4) and the soft-crosslink value (1.48e-4).

So the lit-anchored stiff crosslink stiffness is now in production with a clean, floored γ. Per-seed
diagnostics: γ_xl = 0.0 (force-free resting crosslinks), γ_myo ≈ 0.13–0.16 pN/µm (floored), all seeds
stable. This also **re-confirms, with the CORRECT stiffness, that the γ-floor is force-magnitude-limited
and robust to crosslink stiffness** (the 6M finding, now the production default rather than a workaround).

## Scope

- The γ-floor production path (resting settle + measure) uses `crosslink_turnover=True`.
- The DYNAMIC paths (`simulate_loaded_shell_on_device`, `simulate_turnover_on_device`) carry the stiff
  crosslinks as explicit springs → smaller dt (handled by the CFL in `relax_on_device`); they stay
  stable (stiff crosslinks make the shell MORE rigid) but a continuous crosslink-turnover during
  evolution is the faithful long-run treatment (extends the 6k myosin-turnover kernel to crosslinks).
- `kim_network.shear_modulus` uses the stiff k_xl via its own affine-shear solver (elastic axis).
- SE rows for the Ferrer stiffness datum + a `params_manifest` entry are the remaining KB-sync step
  (`SE_REGISTRATION_CANDIDATES_2026-06-30.md` §5); link_k is a code constant so kb-check shows no drift.

Tests: `test_gamma_floor_crosslink_turnover_robust_baseline` (γ_xl≈0, active floored, robust across
seeds); full FF suite 77 pass, kb-check green. Related: FF_STAGE6M (the finding),
[[project-gamma-floor-likely-deficit]] (crosslink stiffness closed as a root cause).
