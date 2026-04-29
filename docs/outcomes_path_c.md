# Path C — bounded outcomes (v15 baseline single pilot, default scope)

This document records the **PI-pre-authorized bounded-outcome decision tree**
for the Path C v15 baseline pilot. Written before runs execute.

The companion documents are:
- `docs/path_c_sanity.md` — six-check Sanity Gate + Magic-Number Block + PI
  full-authorisation record
- `docs/outcomes_v15.md` — original v15 baseline (R drift 0.244, free-floating)
- `docs/outcomes_stage1a_plus.md` — Stage 1a+ Option α/β baselines (R drift
  0.247 under mechanical-only / energetic-wetting substrate)
- `docs/outcomes_stage1b.md` — Stage 1b Bucket I result (3-phenotype ordering
  reproduced via R drift; A/A₀ measurement broken by Stage 1a+ inherited
  lift-off pattern)

## Mechanism question (one)

**Does adding effective gravity to the Layer 1 + substrate baseline (i)
prevent the inherited Stage 1a+ lift-off pattern, (ii) restore meaningful
A/A₀ measurement, and (iii) give a R drift that is materially different
from the v15 free-spheroid baseline 0.244 and the Stage 1a+ Option α
mechanical-only baseline 0.247?**

The pilot is a **single simulation** producing a **multi-metric snapshot**
(R drift, A/A₀ trajectory, contact-band ρ_kernel, anchor force balance,
sphericity) evaluated against the inherited v15 / 1a+ baselines.

## Bounded outcomes (4 buckets)

The buckets are defined on the *combination* of R drift and A/A₀
trajectory, since Path C is fundamentally a measurement-protocol fix
(restoring meaningful A/A₀) AND a mechanism addition (gravity prevents
lift-off).

| Outcome | Definition | Action |
|---|---|---|
| **Bucket P1 — Path C succeeds** | R drift < 0.247 (Stage 1a+ Option α baseline) AND A/A₀ trajectory passes the new "min ≥ 0.5" gate (no catastrophic lift-off) AND contact-band ρ_kernel ∈ [0.85, 1.15] for majority of post-transient frames | STOP. Surface to PI as **the framework fix**: Path C resolves the inherited Stage 1a+ measurement-protocol issue. PI decides whether to (i) re-baseline Stage 1a+ Option β / Stage 1a++ / Stage 1b under Path C (which would fully reconcile the A/A₀ ordering with PI experimental data, completing layer-by-layer narrative pillar 4 with proper measurement), (ii) advance to Stage 1c (Layer 5 mechano-osmotic), (iii) refine `g_star` via a small sweep, (iv) report as publication-strong framework pillar. **No automatic re-baseline / Stage 1c entry.** |
| **Bucket P2 — Path C partial success** | R drift < 0.247 AND A/A₀ trajectory passes ≥ 0.5 gate but contact-band ρ_kernel still fails OR sphericity degrades (gravity flattening too aggressive) | STOP. Surface to PI. Path C reduces lift-off but does not fully fix substrate anchoring. PI decides among (i) refine `g_star` (lower for less flattening, higher for more anchor), (ii) revisit substrate-CSF formulation, (iii) accept Path C as partial fix and proceed to re-baseline. |
| **Bucket P3 — Path C insufficient** | R drift ≥ 0.247 OR A/A₀ trajectory still fails (lift-off persists) at provisional `g_star = 0.1` | STOP. Surface to PI. The provisional `g_star` magnitude is too small to anchor; either sweep to higher `g_star` (PI sub-decision) or revisit the framework (Path C alone is insufficient; substrate-CSF + gravity together still cannot anchor). |
| **Bucket P4 — Critical failure** | NaN/Inf, max-speed runaway, mass-conservation break, energy-monotone (extended sum) FAIL by orders of magnitude, or `gravity_star` causes spheroid to disintegrate | STOP **immediately**. Implementation issue (gravity sign, energy-budget bug, etc.). Diagnose before proceeding. |

## Per-run gates

Per Path C sanity-md §"New gate / diagnostic candidates":
- Stage 1a+ Option α inherited: mass / horizontal momentum / no NaN /
  max-speed / VRAM / v15 inherited gates.
- **Energy monotone (extended sum)**: `KE + U_strain + U_surface_free
  + U_grav`, monotone-decay tolerance unchanged (gravity is conservative;
  spheroid settling does positive work on cells, but the integrated
  total energy decreases via overdamped drag).
- **Anchor force balance (recomputed integrand)**: `|F_sub - (∫σ_vol_zz
  dA + M_spheroid · g_star)| / |F_sub|` ≤ 0.20.
- **A/A₀ trajectory finite & ≥ 0.5**: NEW MEANINGFUL GATE under Path C
  (resolves Stage 1a+ inherited lift-off pattern).
- **R drift improvement vs Stage 1a+ Option α baseline**: strict-less
  than 0.247.

## Diagnostics (no gate, log-only)

- Per-frame `U_grav`, `total_energy = KE + U_strain + U_surface_free + U_grav`
- Per-frame contact-band particle count `n_contact_band_particles`
- Per-frame contact area `A_contact_xy_hull` and ratio `A_contact / A_0`
- Per-frame sedimentation depth `<z_p>` (spheroid centre-of-mass z position)
- Per-frame apparent contact angle θ (geometric, log only)

## Files of record

```
results/path_c_v15_baseline/{gate_report.md, metrics.csv, shell_profile.csv,
  contact_metrics.csv, snapshots.h5, ...}
```

## Auto-progression rules (PI 2026-04-29 default scope)

**Default scope per PI directive**: **v15 only re-baseline 후 STOP, 다른
stage 결정 대기**. Concretely:
- Run the Path C v15 baseline pilot once.
- Classify into bucket P1 / P2 / P3 / P4.
- STOP. Surface to PI. **No automatic re-baselining of Stage 1a+ /
  1a++ / 1b under Path C.**

**Auto-STOP conditions** (any of these halts execution immediately):
- All Stage 1a+ Option α auto-STOP conditions inherited.
- Bucket P4 trigger (NaN/Inf, max-speed runaway, mass-conservation
  break, energy-monotone fail by orders of magnitude).
- pytest fails before the pilot starts.
- `gravity_star` causes per-step velocity to exceed v_rms · 100 (factor
  10× safety margin above the existing max_speed_over_vrms gate).

## Stop conditions (general; unchanged)

- Magic number 도입 금지 (`g_star = 0.1` is explicitly framed as a
  body-force coefficient with PARTIAL Magic-Number Block per the
  dimensional-coefficient framing; PI authorisation suffices).
- Gate semantic 변경 금지 beyond the recorded Path C extensions
  (energy-monotone sum to include `U_grav`; anchor-force-balance
  integrand to include gravity term).
- v13 anti-pattern (measurement-protocol consistency 위반) 시 즉시 STOP.
- 한도 외 escalation 금지 — Stage 1a+/1a++/1b re-baseline 자동 진입 금지.
