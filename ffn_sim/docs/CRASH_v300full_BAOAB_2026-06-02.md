# Crash diagnostic — `v300full` BAOAB box-wrap blow-up (2026-06-02 ~06:01)

> **READ-ONLY analysis by the parallel prep session** (`AUTONOMOUS_LOG_2026-06-02.md`
> iter 13). I did **NOT** touch the integrator, the driver, the config, or relaunch —
> all are frozen / PI-gated / Lead-owned. This note is a grounded diagnostic + safe
> recovery options **for the Lead/PI to decide**, surfaced because the PI assigned
> crash-watch. Nothing here is applied.

## The failure (exact)

The Lead's `gpu_native38k_v300full` run (`scripts/h3_ku35_gripwalk_tier1.py`,
**v0_accel=300×, force_scaling ON, arm=grip_walk**, GPU native n_fil=38,000) died at
~06:01, **right after sample 11/40**, inside `sim.run()`:

```
FloatingPointError: Constrained BAOAB _wrap_into_box_xp: |fractional coord| exceeded
the int32-image guard (worst |round(f)|=3.255e+08 > 1e+08)
  integrator/constrained_baoab.py:949  act()  → step 4 (box wrap, GPU sibling)
  integrator/constrained_baoab.py:192  _wrap_into_box_xp()
```

## What the failure means (code-grounded, NOT speculation)

- The guard is the **int32-image overflow guard** (`INT32_GUARD = 1e8`,
  `constrained_baoab.py:182`) in the GPU-port box-wrap sibling of the frozen
  `_wrap_into_box`. It fires when a particle's fractional coordinate would need an
  image index > 1e8 — i.e. the particle is ~3.26×10⁸ box-lengths out. **The guard
  worked correctly**: it failed loud rather than letting the int32 image flag overflow
  / corrupt silently. No integrator bug; the guard is a safety net doing its job.
- Critically, the **SHAKE non-finite guard one step earlier did NOT fire**
  (`:938` checks `isfinite(projected)`). So the position was **finite but enormous** —
  not a NaN/Inf. A particle took a **single huge finite displacement**.
- The logged **constraint drift stayed ~e-15** through sample 11 (`_max_drift`,
  `:934`) — the rigid chains were healthy. So the blow-up was a **single over-forced /
  non-constrained particle**, not a slow constraint drift or a global instability.

## Why now — correlated with transport engaging (the likely cause)

| samples | adv (myosin steps) | drift | state |
|---|---|---|---|
| 1–8 | **0** (walk not engaged) | ~3.4e-15 | stable, γ on floor |
| 9 | **0 → 17424** (transport ENGAGES) | 3.39e-15 | first walk |
| 10–11 | 19073 → 19390 | 3.6–4.9e-15 | sustained |
| 11 → 12 | — | 💥 finite blow-up | crash |

The blow-up is a **few steps after the grip-walk transport engaged**. Grounding from
the driver (`h3_ku35_gripwalk_tier1.py:98-118`): `v0_accel` multiplies the literal
myosin `v0_per_head` (here **×300**), and the docstring itself flags this accelerant as
**diagnostic-only** — *"the v0-only accelerant speeds the myosin WALK but not …
BINDING … transport engages but coherence does not."* With `force_scaling`
(`mesoscale_force_scaling`) **ON** as well, once the accelerated walk actually moves a
grip point (sample 9), the **300× accelerated transport × mesoscale force-scaling**
delivers a force/displacement large enough to fling a head/actin bead out of the box in
one step → the int32 guard trips.

**This is consistent with the documented design caveat, not a new defect:** v0_accel is
an unphysical diagnostic accelerant; once transport engages, the accelerated force is
exactly the thing that can destabilize — most so with force_scaling ON. Per the Notion
board, **"force-scaling OFF = gold-standard native call"** — this v300full run is the
force-scaling-**ON** diagnostic arm, so the instability may be specific to the ON arm.

## Safe recovery options (for Lead/PI — all driver/config, NO integrator change)

Ranked least-invasive first. The integrator guard is correct and needs no change.

1. **Re-run as the force-scaling-OFF gold-standard.** The ratified gold-standard call is
   force-scaling OFF (Notion). This crash is in the ON diagnostic arm; the OFF arm
   avoids the mesoscale force multiplier that compounds the accelerated transport.
   *Lowest effort, matches the ratified config.*
2. **Lower `v0_accel`** (e.g. 100× or 50× — the driver's own example uses 50×). Less
   aggressive acceleration ⇒ smaller transport force at engage. Trades diagnostic speed
   for stability.
3. **`--couple-accel`** (scale binding `k_on` with `v0_accel`). The driver's documented
   fix for the walk-vs-binding timescale decoupling; keeps the network coherent in the
   accelerated window, which may reduce the over-forcing at engage. (Changes physics —
   PI call.)
4. **Smaller `dt` (`dt_factor`).** Generic stability lever; but the guard trips on force
   *magnitude*, so dt alone is a blunt tool unless paired with (2).

## Status / preserved data

- **Samples 1–11 + the 68 MB GSD are intact** (transport-engaged data through sample 11
  captured) — the run is recoverable, not lost.
- As of iter 14 (~06:21) the Lead had **not yet relaunched** (no new production artifact
  in 25 min; Lead HEAD 877cf40 unchanged). I am holding and watching; recovery is the
  Lead/PI decision. I have not SSHed gbook, not edited any file, not relaunched.
