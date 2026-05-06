# Stage 1c — bounded outcomes (Layer 5 mechano-osmotic, single pilot)

This document records the **PI-pre-authorized bounded-outcome decision tree**
for the Stage 1c Layer 5 mechano-osmotic Tier 2 pilot. Written before runs
execute.

The companion documents are:
- `docs/v1/stage1c_sanity.md` — six-check Sanity Gate + Magic-Number Block
- `docs/08_mechano_osmotic.md` — pre-existing Tier 2 framework
- `docs/v1/outcomes_stage1b.md` — Stage 1b Bucket I baseline (R drift
  0.158-0.185 across 3 phenotypes)
- `docs/v1/outcomes_path_c.md` — Path C v15 baseline (g_star=0.1 pancake;
  recalibrated to 0.01 for Stage 1c)
- `docs/12_validation.md` — Pillar 3 (PI experimental data is cross-check)

## Mechanism question (one)

**Does Layer 5 (mechano-osmotic Tier 2 K(ρ_osm) coupling) on top of the
full carrier (Layer 1 + 1a+ β α=1.0 + 1a++ ζ=0.4 + 1b φ Pre + Path C
g_star=0.01) give a meaningful contribution to spheroid stability /
spreading vs the Stage 1b baseline?**

The pilot is a **single simulation** producing a multi-metric snapshot
(R drift, A/A₀ trajectory, ρ_osm trajectory, K_eff bulk values, contact
band, sphericity). Comparison baselines:
- Stage 1b Pre carrier (without Layer 5, without Path C) R drift = 0.166
- Path C v15 baseline (Layer 1 + α + Path C g_star=0.1) R drift = 0.326
  (over-anchored pancake)

## Bounded outcomes (4 buckets)

| Outcome | Definition | Action |
|---|---|---|
| **Bucket S1 — Layer 5 + Path C g_star=0.01 succeeds** | R drift < 0.166 (Stage 1b Pre baseline) AND A/A₀ ≥ 0.5 throughout AND sphericity ≥ 0.85 (no pancake) AND ρ_osm trajectory monotone-rise above 1.0 (Guo mechanism reproduced) | STOP. Surface to PI as the **layer-by-layer narrative quantitative pillar 5** finding: Layer 5 Tier 2 + recalibrated Path C contributes to spheroid stability + spreading, with spreading-induced water efflux (Guo 2017 mechanism) reproducing the experimental volume-stiffness coupling. PI decides whether to (i) advance to Stage 1d (Layer 4 Marangoni / nematic / vortex), (ii) re-baseline Stage 1b 3-phenotype with Layer 5 + Path C (Bare/Pre/Lam4 comparison with proper A/A₀), (iii) refine Layer 5 (activate η or σ_active(ρ) couplings), (iv) report as publication-strong finding. **No automatic Stage 1d entry.** |
| **Bucket S2 — partial improvement / qualified success** | R drift < 0.247 (substrate-stage baseline) but ≥ 0.166 (Stage 1b baseline), OR sphericity moderately degraded (∈ [0.7, 0.85]), OR ρ_osm trajectory plateaus instead of rising | STOP. Surface to PI. Layer 5 contributes but does not fully break the Stage 1b ceiling, OR Path C g_star=0.01 is mildly over-anchor. PI decides among (i) refine g_star (smaller), (ii) refine Layer 5 ODE constants, (iii) activate additional Layer 5 couplings, (iv) accept as partial finding. |
| **Bucket S3 — no improvement / lift-off returns** | R drift ≥ 0.247 OR A/A₀ trajectory still fails (lift-off persists with g_star=0.01 too small) OR ρ_osm trajectory flat (water-flux mechanism inactive) | STOP. Surface to PI. Path C g_star=0.01 may be under-anchored; OR Layer 5 ODE rate constants need adjustment; OR Layer 5 K(ρ) coupling alone is insufficient. PI decides: (i) sweep g_star upward, (ii) sweep α_osm_star, (iii) activate σ_active(ρ_osm) coupling. |
| **Bucket S4 — divergence / NaN / pancake / runaway** | NaN/Inf, max-speed runaway, mass-conservation break, ρ_osm exits [0.5, 1.6] bounds, sphericity < 0.5 (strong pancake), or any pre-existing scheme-correctness gate fail by orders of magnitude | STOP **immediately**. Layer 5 implementation issue OR Path C g_star=0.01 still too large OR ODE stiffness violation (constructor invariant). Diagnose before proceeding. |

## Per-run gates

Per Stage 1c sanity-md §"Per-run gates":
- All Stage 1b inherited (mass / horizontal momentum / no NaN / max-speed
  / VRAM / v15 / Layer 2 active-power / boundary-tag stability / Layer 3
  φ ∈ [0,1] / φ trajectory / A/A₀ ≥ 0.5).
- Energy-monotone SUSPENDED (Stage 1a++ contract change).
- **ρ_osm ∈ [0.5, 1.6] per-particle invariant** (NEW): min(ρ_osm) ≥ 0.5,
  max(ρ_osm) ≤ 1.6 across all frames.
- **<ρ_osm> trajectory finite & non-pathological** (NEW): `<ρ_osm>(end)`
  finite (not NaN); `<ρ_osm>(end) ∈ [0.95, 1.6]` (allows slight relaxation
  below 1.0 but no catastrophic loss).

## Diagnostics (no gate, log-only)

- Per-frame `<ρ_osm>_bulk`, `<ρ_osm>_boundary`, `<ρ_osm>_contact`,
  `min(ρ_osm)`, `max(ρ_osm)`.
- Per-frame `<K_eff>_bulk = K · <ρ_osm>_bulk`.
- Per-frame mean strain rate `<-tr(C)>_bulk` (Guo mechanism driver).
- All Stage 1b diagnostics inherited.

## Files of record

```
results/stage1c_pilot/{gate_report.md, metrics.csv, shell_profile.csv,
  contact_metrics.csv, snapshots.h5, ...}
```

## Auto-progression rules (PI 2026-04-29 full authorization, Stage 1c only)

**Auto-progression IS allowed within Stage 1c**:
- Sanity_md → outcomes → implementation → pytest → pilot → analysis →
  commit. All within one stage, automatic.

**Auto-STOP conditions**:
- All Stage 1b auto-STOP conditions inherited.
- Bucket S4 trigger (NaN/Inf, runaway, ρ_osm bounds violation,
  catastrophic pancake).
- pytest fails before pilot.
- Layer 5 ODE introduces verifiable inconsistency (e.g. ρ_osm trajectory
  diverges from predicted Guo mechanism by orders of magnitude).

**Auto-STOP at end of Stage 1c**:
- After 1-pilot completes (or auto-STOP fires), STOP.
- Stage 1d (Layer 4 Marangoni) auto-entry is FORBIDDEN.
- Layer 5 refinement / phenotype re-baseline / publication-finding
  consolidation requires PI input.

## Stop conditions (general; unchanged)

- Magic number 도입 금지 (Layer 5 ODE constants are doc-attributed
  PI-authorized values per `08_mechano_osmotic.md`; not magic numbers).
- Gate semantic 변경 금지 beyond inherited Stage 1a++ contract changes
  + the new ρ_osm invariant gate (which is a state-variable bound, not
  a fitted physics gate).
- v13 anti-pattern (measurement-protocol consistency 위반) 시 즉시 STOP.
- 한도 외 escalation 금지 — Stage 1d 자동 진입 금지.
