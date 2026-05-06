# Production Lam4 — bounded outcomes (5k × 80 sim hr, single pilot)

This document records the **PI-pre-authorized bounded-outcome decision tree**
for the production Lam4 single-phenotype run. Written before execution.

The companion documents are:
- `docs/v1/production_lam4_sanity.md` — six-check Sanity Gate (production-
  scale specific)
- `configs/production_lam4.yaml` — production config (Lam4 phenotype +
  zeta_max=0.6 + 5k × 80 hr)
- Phase 4-v2 commit `708b44d` — Lam4 pilot baseline (R drift 0.073,
  A/A₀_topdown 2.657 at pilot 4 hr)
- PI experimental data: `data/experimental/260313_Lam4.csv`
  (A/A₀_final range 8.0–33.1 over 83 hr)

## Mechanism question (one)

**Does the full Stage 2 framework (Layer 1+2+3+4+5+6 + Path C +
Phase 1 fixes + zeta_max=0.6) at production scale (5k × 80 hr =
~PI experimental duration) reproduce the PI Lam4 experimental
A/A₀_final range [8.0, 33.1]?**

The pilot 4-hr A/A₀_topdown value 2.657 was substantially below PI
endpoint range 8-33; the Phase 4-v2 finding showed ordering
reproduction. Production scale tests whether the gap was time-budget
(pilot 4 hr is 5% of PI 82 hr; trajectory may continue to spread)
OR mechanism-missing (the model has reached an asymptote near 2.7 and
will not reach PI range regardless of time).

The pilot is a **single simulation** producing R(t), A/A₀_topdown(t),
sphericity(t), φ(t), ρ_osm(t), mmp(t), ecm(t) trajectories at 320
frames over 80 sim hr.

## Bounded outcomes (4 buckets)

| Outcome | Definition | Action |
|---|---|---|
| **Bucket P1 — Time-budget was sufficient cause; framework reaches PI range** | A/A₀_topdown(end) ∈ [4.0, 33.1] (within or above PI Lam4 lower bound, with monotone-increasing trajectory) AND no Bucket-P4 failure | STOP. Surface to PI as **time-budget diagnostic positive**: the pilot result was a finite-time artifact, not a missing mechanism. PI decides whether to (i) launch 3-phenotype production (Bare/Pre/Lam4 ~6-9 hr wall-clock); (ii) consolidate the publication-strong narrative (26 simulations + 1 production = framework validated at experimental scale); (iii) extend to longer times for asymptotic behaviour. **No automatic next-stage entry.** |
| **Bucket P2 — Time + mechanism both contribute; partial reproduction** | A/A₀_topdown(end) ∈ [2.5, 4.0] (substantially above pilot 2.657 but below PI lower bound 8.0) | STOP. Surface to PI. The model continued to spread beyond pilot but plateaued before reaching PI range. PI decides among (i) Stage 1d.b (lamellipodia / filopodia stochastic events for additional spreading driver); (ii) phenotype-dependent ζ_max refinement; (iii) accept as partial finding with quantitative gap analysis. |
| **Bucket P3 — Mechanism-missing; framework asymptotes below PI range** | A/A₀_topdown(end) < 2.5 (essentially unchanged from pilot 2.657, OR even decreased) | STOP. Surface to PI. Stage 2 framework cannot reach PI experimental Lam4 range even at full experimental duration. PI decides among (i) Stage 1d.b (additional active mechanisms); (ii) Layer 3 φ-ODE refinement (per-particle spatial heterogeneity); (iii) accept as paper limitation with quantitative narrative ("our framework captures the qualitative phenotype ordering but the absolute spreading magnitude requires additional mechanisms beyond the 6-layer phenomenology"). |
| **Bucket P4 — Numerical / scaling failure** | NaN/Inf, VRAM overflow > 12 GB, max-speed runaway, mass drift > 1e-9 over 480k steps, ecm_strength bounds violation, runaway substrate degradation, or any pre-existing scheme-correctness gate fail by orders of magnitude | STOP **immediately**. Production-scale numerical issue (drift accumulation, VRAM, etc.). Diagnose before considering further runs. |

## Per-run gates (production-scale specific)

Per Production Lam4 sanity-md:
- All Phase 4-v2 inherited gates (mass / horizontal-momentum / no-NaN
  / max-speed / VRAM / v15-inherited gates / Layer 2 active power
  / boundary-tag stability / Layer 3 φ ∈ [0,1] / φ trajectory / Layer
  5 ρ_osm bounds / Layer 5 trajectory / Layer 6 ecm bounds / Layer 6
  mmp monotone).
- Energy-monotone SUSPENDED (Phase 1.1 Cousin-Rule extension covers
  L3/L4/L5/L6 active).
- **Mass drift gate relaxed to 1e-9** (480k-step f32 round-off
  budget); other tolerances unchanged.

## Diagnostics (no gate, log-only)

- All Phase 4-v2 diagnostics inherited.
- **NEW long-time diagnostic — A/A₀_topdown trajectory shape**: plot
  A/A₀_topdown(t) over 320 frames; classify trajectory shape
  (monotone-increasing / saturating / non-monotone) for Bucket
  classification.
- **NEW PI cross-comparison overlay**: at common time points (sim
  t = PI t = 60, 600, 1200, ..., 5040 min), compare
  A/A₀_topdown_sim(t) vs A/A₀_PI_Lam4(t). Direct quantitative overlay.

## Files of record

```
results/production_lam4/{gate_report.md, metrics.csv, shell_profile.csv,
  contact_metrics.csv, snapshots.h5, ...}
results/production_lam4_console.log
```

## Auto-progression rules (PI 2026-04-29 full authorization, production only)

**Auto-progression IS allowed within production Lam4 single run**:
- Sanity_md → outcomes → config → background pilot → analysis →
  bucket classification → commit. All within one production run.

**Auto-STOP conditions**:
- All inherited auto-STOP conditions.
- Bucket P4 trigger (NaN/Inf, VRAM > 12 GB, mass drift > 1e-9, etc.).
- Wall-clock > 5 hr (3× original estimate); PI re-evaluation needed
  if compute is taking longer than budgeted.

**Auto-STOP at end of production Lam4**:
- After production completes (or auto-STOP fires) + bucket
  classification + commit, STOP.
- **Other phenotype production (Bare/Pre) auto-entry FORBIDDEN**.
- Stage 1d.b / paper-draft auto-entry FORBIDDEN.
- Inherited fixes auto-entry FORBIDDEN.

## Stop conditions (general; unchanged)

- Magic number 도입 금지 (no new parameters in production).
- Gate semantic 변경 금지 beyond inherited contract changes + the
  one explicit relaxation (mass drift 1e-10 → 1e-9 for long-time
  budget).
- v13 anti-pattern (measurement-protocol consistency 위반) 시 즉시
  STOP.
- 한도 외 escalation 금지 — other-phenotype production / Stage 1d.b
  / paper draft / framework next-step 모두 자동 금지.
