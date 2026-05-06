# Stage 1a++ — bounded outcomes (Option II 3-run sweep, α' framing)

This document records the **PI-approved bounded-outcome decision tree** for
the Stage 1a++ Option II 3-run sweep (`ζ/K ∈ {0.1, 0.3, 1.0}`, no Pa
claim per Option α' resolution). It is written *before* the runs execute
so the result is read against a contract that exists prior to any data.

The companion documents are:
- `docs/v1/stage1a_plus_plus_layer2_sanity.md` — six-check Sanity Gate +
  Magic-Number Block + PI Option α' resolution (no Pa claim;
  K-anchored to Fischer-Friedrich Nat Cell Biol 2014 IF 30; Marchetti
  Rev Mod Phys 2013 IF 50 framework reference)
- `docs/v1/v15_deferred_diagnostics_plan.md` — Track B reference (trigger
  conditions defined here)
- `docs/v1/outcomes_stage1a_plus.md` — Stage 1a+ Option β baseline
  (R drift = 0.247 at α=1.0 carrier; this is the comparison baseline
  for Stage 1a++)
- `docs/v1/outcomes_v15.md` — v15 Layer 1 baseline (R drift = 0.244)

## Mechanism question (one)

**Does adding boundary-cell active stress (Layer 2 continuum entry,
σ_act = -ζ·K·I on `is_boundary[p]==1` particles) reduce the Stage 1a+
Option β α=1.0 baseline R drift of 0.247 (= the Layer 1 + substrate
ceiling)? If so, how does the response scale with ζ/K?**

The sweep is a **single response curve** — three (ζ/K, R drift) data
points giving the qualitative shape of the dependence. All other gates
are pass/fail correctness contracts and play no role in the bounded-
outcome bucketing.

## Bounded outcomes (4 sweep meta-buckets)

The buckets are defined on the *shape of the response curve* across the
three runs (ζ/K = 0.1, 0.3, 1.0). All R drift values use the same
measurement protocol as v15 / Stage 1a+ (`max_t |R(t)/R₀ − 1|` for
`t* ≥ 10`).

| Sweep meta-bucket | Definition | Action |
|---|---|---|
| **Bucket A — linear / monotone reduction** | R drift decreases monotonically with ζ/K (e.g. 0.22 → 0.18 → 0.10), with at least one ζ value giving R drift < 0.20 (i.e. < 20%, materially better than the 24.7% Stage 1a+ baseline) | STOP. Surface to PI as the *layer-by-layer narrative quantitative pillar 3* finding: Layer 2 boundary active stress contributes to spheroid integrity, with magnitude scaling with ζ/K. PI decides whether to (i) advance to Stage 1b (Layer 3 φ-ODE; introduces Bare/Pre/Lam4 phenotype mapping per `docs/SESSION_HANDOFF.md` §"Stage 1b framing memo"), (ii) refine the sweep, or (iii) activate Stage 1a++.b (lamellipodia / filopodia stochastic events). **No automatic Layer 3 entry** — Layer 3 mapping requires PI input (see "Layer 3 entry constraint" below). |
| **Bucket B — threshold / saturation** | R drift shows a clear non-monotone or saturating shape (e.g. 0.24 → 0.20 → 0.20, or 0.24 → 0.16 → 0.16); at least one ζ improves drift, but the response is not linear in ζ/K | STOP. Surface to PI. *Suggestive* of a regime change (active stress activates above some threshold, or saturates against another mechanism). Track B trigger may apply if the saturation level is still ≥ 22% (i.e. Layer 2 alone is bounded by the same architectural ceiling). PI decides among (i) refine sweep around the threshold, (ii) enter Track B selectively (per `docs/v1/v15_deferred_diagnostics_plan.md`), (iii) advance to Stage 1b regardless. |
| **Bucket C — small / no improvement** | All three ζ/K give R drift ≥ 22% (within 2.7% of the Stage 1a+ Option β α=1.0 baseline 0.247). Layer 2 continuum boundary active stress is also insufficient to break the architectural ceiling | STOP. **Track B fires** per `docs/v1/v15_deferred_diagnostics_plan.md` (Scenario A trigger). Three independent diagnostics — (i) Maxwell deviatoric residual / (ii) G2P kernel-density bias / (iii) small-strain effective bulk modulus — become live; PI selects which to execute first. The Layer 1 architectural-ceiling interpretation is at risk of being a numerical artifact and the layer-by-layer narrative needs verification before Stage 1b proceeds. |
| **Bucket D — divergence or NaN** | Any of the three runs halts on NaN/Inf, max-speed runaway, mass-conservation break, or any pre-existing scheme-correctness gate (calibration <J>=1, CSF localisation, anchor force balance, etc.) | STOP **immediately**. This is **Layer 2 implementation issue**, not a mechanism finding (per Scenario B in PI specification). Surface the failing gate + numerical state. Do **not** proceed to remaining sweep runs; do **not** trigger Track B (the bug is in Layer 2 code, not in v15). |

## Per-run gates (each of the 3 runs)

Same gate set as Stage 1a+ Option β plus the 3 PI-approved contract
changes (energy-monotone suspended, anchor-force-balance integrand
recomputed with σ_act,zz, contact-band ρ_kernel reinterpretation) and
the new gates introduced by the Layer 2 sanity-md:

- All Stage 1a+ gates (mass / momentum-horizontal / no-NaN / max-speed
  / VRAM / v15-inherited gates / contact-band ρ_kernel)
- **Anchor force balance** (recomputed): integrand now includes
  σ_act,zz at boundary-tagged contact-band particles; tolerance ≤ 0.20
  unchanged
- **R drift improvement vs Stage 1a+ Option β α=1.0 baseline** (new):
  `R_drift_1aplusplus < 0.247` (strict-less; this is *the* mechanism
  question gate per run)
- **Boundary-tag stability** (new): `|Δn_boundary / n_boundary|` per
  frame ≤ 0.10 (numerical hygiene; if the boundary tag flickers
  pathologically, the active-stress measurement is meaningless)
- **Active-work finite** (new): `W_active(t)` finite at end of run
  (not NaN, not > 10·U_strain at end). Diagnostic for the energy-
  monotone suspension contract change
- **Energy-monotone (KE + U_strain)** — **SUSPENDED** under Stage 1a++
  (active stress injects energy by construction; Cousin-Rule contract
  change recorded in sanity-md)

## Diagnostics (no gate, log-only)

- shell-averaged `<ρ_kernel>(r/R₀)` (carried over from v15)
- substrate contact area `A_contact` (carried over from Stage 1a+)
- apparent contact angle θ (carried over from Stage 1a+)
- per-frame `n_boundary`, `n_contact_band`
- per-frame cumulative active work `W_active(t)`
- per-frame boundary-tag flicker (`Δn_boundary` per frame)

## Files of record (per run)

```
results/stage1a_plus_plus_zeta_0p1/{gate_report.md, metrics.csv, ...}
results/stage1a_plus_plus_zeta_0p3/{gate_report.md, metrics.csv, ...}
results/stage1a_plus_plus_zeta_1p0/{gate_report.md, metrics.csv, ...}
```

Plus a sweep-level summary written to `results/stage1a_plus_plus_summary.md`
after the third run completes (response-curve table + meta-bucket
classification).

## Stop conditions (auto-progression rules per PI 2026-04-29 directive)

**Auto-progression IS allowed within Stage 1a++** (Layer 2 끝까지):
- Run 1 (ζ/K = 0.1) completes with no NaN/divergence → automatically
  proceed to Run 2 (ζ/K = 0.3).
- Run 2 completes with no NaN/divergence → automatically proceed to
  Run 3 (ζ/K = 1.0).
- Run 3 completes → automatically aggregate the response curve and
  classify the meta-bucket.

**Auto-STOP conditions (any of these halts the sweep immediately)**:
- NaN/Inf in any frame of any run (existing v15 NaN-detection gate
  catches this).
- Max-speed runaway (`max v* > max_speed_over_vrms · max(v_rms,
  V_FLOOR)`).
- Mass conservation violation (`|Δm/m₀| > 1e-10`).
- pytest `tests/test_physics_conservation.py` fails (run before
  starting the sweep).
- Active-work cumulative `W_active(t)` becomes unbounded (> 10·U_strain
  at any frame; new gate per Stage 1a++ contract).
- Boundary-tag flicker > 0.30 per frame (deemed degenerate
  measurement).

In any auto-STOP case: surface the failing gate + numerical state +
which run was running. Do **not** trigger Track B from a Bucket D
failure (that is a Layer 2 implementation issue, not a mechanism
finding).

**Auto-STOP at end of Stage 1a++** (regardless of result):
- After 3-run sweep completes (or auto-STOP fires), STOP.
- **Layer 3 (Stage 1b) auto-entry is FORBIDDEN.** PI input is required
  before Stage 1b begins because Layer 3 introduces the Bare/Pre/Lam4
  phenotype mapping (formation-environment → initial φ / γ_cc; see
  `docs/SESSION_HANDOFF.md` §"Stage 1b framing memo") which is a
  PI-experiment-data-direct touchpoint.

## Layer 3 entry constraint

Per PI 2026-04-29 directive, **Stage 1b auto-entry is forbidden**. The
reason: Layer 3 activates the φ-ODE for E-cad ↔ Int-β1 transition, and
the three Bare/Pre/Lam4 conditions enter as initial conditions on the φ
field and on initial γ_cc, mapped from the formation-environment
phenotype to starting-state parameters. This mapping is a *direct
touchpoint with PI's experimental data* (Cho et al. 2020 mechanism +
PI's own laminin-presentation extensions on pV4D4 hydrogel) and
requires PI's domain expertise. Stage 1b initial-condition mapping is
not derivable from purely mechanistic principles; PI input on the
formation-environment → initial-state mapping function is required.

After Stage 1a++ completes (regardless of meta-bucket), the runner
STOPs and surfaces the result to PI. PI then decides whether to (i)
advance to Stage 1b with their selected φ / γ_cc mapping, (ii) refine
the Stage 1a++ sweep, (iii) activate Stage 1a++.b (Layer 2 stochastic
events: lamellipodia / filopodia / leader cells / discrete focal
adhesions), or (iv) trigger Track B selectively.

## Stop conditions (general; unchanged)

- Magic number 도입 금지
- Gate semantic 변경 금지 (Cousin Rule); the three Stage 1a++ contract
  changes (energy-monotone suspend, anchor-force integrand recomp,
  contact-band ρ_kernel reinterpretation) are explicit pre-recorded
  exceptions per the sanity-md
- v13 anti-pattern (measurement-protocol consistency 위반) 시 즉시 STOP
- 한도 외 escalation 금지 — Stage 1b 자동 진입 금지 (위 명시)
