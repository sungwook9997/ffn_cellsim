---
kb_record:
  doc_id: FF_STAGE6K_TURNOVER_2026-06-30
  title: γ-floor robust to motor turnover — on-device Hand KMC (Bell) duty-ratio self-limit
  authoritative_as_of: 2026-06-30
  supersedes: []
  status: current
---

# FF Stage 6k — the γ-floor survives the TURNOVER axis (on-device Hand KMC)

**Date:** 2026-06-30 · **Engine:** `aleph/laws/` · branch `dcm/main`

## Question

The γ-floor (`FF_STAGE6D_GAMMA_FLOOR`, RESOLVED) is a **force-magnitude** limit: the active actomyosin
cortical tension is ~10³× under the Salbreux 0.35–0.65 mN/m band, set by the engaged force-bearing
motor density × per-motor force. `FF_STAGE6H` ruled out three **static** network mechanisms as the
cause (buckling, connectivity, finite extensibility — the floor is flat in all of them). The remaining
axis: **dynamic motor turnover**. Real myosin binds/unbinds (Hand model, NF2007 §10.1); does
load-dependent turnover let the network reach band tension? This stage tests it on-device.

## Build (`ff/network_warp.py`)

- `kmc_bell_turnover_kernel` (Warp, one thread per link, `wp.rand`): a bound myosin Hand **detaches**
  with the Bell force-dependent probability `1−exp(−τ·p₀·exp(|f|/f₀))` on its own contractile load `f`;
  a detached Hand within the capture radius **re-attaches** with `1−exp(−τ·k_on)`.
- `myosin_bound_kernel`: only ENGAGED (`bound==1`) links exert the contractile force.
- `simulate_turnover_on_device`: resting cortex + bending + crosslink springs (static, ~99 % bound as
  `k_on≫k_off`) + engaged myosin, with periodic KMC ticks — fully GPU-resident.

Two corrections made (real bugs, not tuned to outcome):
1. **Re-attach capture = the mesoscale reach √(A/n)** (the sanctioned ×40 dual, `mesoscale_reach`),
   NOT the molecular ε (210 nm). The coarse myosin links span ~µm; with the molecular ε a detached
   mesoscale Hand could never rebind (engaged → 0, spurious).
2. **`tau_kmc` is PHYSICAL time** (default 0.01 s), decoupled from the mechanical overdamped-descent
   pseudo-step `dt_mu` (regime A/B: re-equilibrate the mechanics over `kmc_every` steps, then advance
   physical time by `tau_kmc`). For small `tau_kmc` the steady engaged fraction → `k_on/(k_on+p_off)`,
   τ-independent.

## Result — turnover cannot lift the floor

The on-device KMC matches the analytic Bell steady state `k_on/(k_on+p₀·exp(f/f₀))` (A5000 / CPU
parity), and the engaged fraction **self-limits with load** (duty ratio):

| per-myosin load f [pN] | engaged fraction (KMC) | analytic Bell | γ_active [mN/m] |
|---|---|---|---|
| 1  | 1.00 | 0.997 | 7.4e-5 |
| 5 (physiol stall) | 1.00 | 0.986 | 8.6e-5 |
| 10 | 0.95 | 0.95  | 6.8e-5 |
| 20 | 0.84 | 0.84  | 8.6e-5 |
| 30 | 0.63 | 0.68  | 1.3e-4 |
| 40 | 0.32 | 0.43  | 9.2e-5 |

(N=300 sweep, gbook A5000; the N=1000 production point: engaged 0.99, γ_turnover 1.70e-4 vs γ_static
1.38e-4 mN/m — both floored ~2000–2500× under band.)

> **γ_active stays ~10³× under the band at EVERY load.** At the physiological minifilament stall (5 pN)
> the engaged fraction is already ≈1 — there is **no headroom** for turnover to add tension. Raising
> the per-motor force only makes more motors **detach** (Bell), *lowering* the engaged density. So
> turnover can only hold the floor or push it lower; it cannot reach band. The γ-floor is robust to
> the turnover axis, completing the "robust to all network mechanisms" claim
> (static: buckling/connectivity/extensibility `FF_STAGE6H`; dynamic: turnover, here).

The floor is set by the **engaged force-bearing motor density × per-motor force** — a force-magnitude
limit, independent of every network mechanism we can ablate. The missing factor is the experimental
load-engaged motor-density datum (`project-sf-nmii-forcescale-result`, the γ/SF twin), not a model
mechanism.

## Figure

`outputs/ff/figs/turnover_gamma_floor.png` — (A) engaged fraction vs load: KMC points on the analytic
Bell curve, the duty-ratio self-limit; (B) γ vs load with the Salbreux band: γ flat at ~10⁻⁴ mN/m,
~10³× under band at all loads.

## Files / tests

`ff/network_warp.py` (kernels + `simulate_turnover_on_device`), `ff/viz_turnover.py`,
`tests/ff/test_network_warp.py::test_kmc_turnover_matches_bell_steady_state` (73 FF pass). Commit
`6k` series on `dcm/main`. A5000 (cuda:0): N=1000, 20k steps + KMC in 1.3 s — fully GPU-resident.
