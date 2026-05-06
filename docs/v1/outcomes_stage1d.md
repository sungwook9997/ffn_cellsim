# Stage 1d — bounded outcomes (Layer 4 Marangoni single pilot)

This document records the **PI-pre-authorized bounded-outcome decision tree**
for the Stage 1d Layer 4 cellular Marangoni pilot. Written before runs
execute.

The companion documents are:
- `docs/v1/stage1d_sanity.md` — six-check Sanity Gate + Magic-Number Block
- `docs/07_internal_flow_dynamics.md` — pre-existing Layer 4 framework
- `docs/v1/outcomes_stage1c.md` — Stage 1c S1-S2 baseline (R drift 0.138)
- `docs/12_validation.md` — Pillar 3 (PI experimental data is cross-check)

## Mechanism question (one)

**Does Layer 4 Marangoni (φ-modulated γ → ∇γ → tangential surface flow)
on top of the Stage 1c full carrier give a measurable contribution to
spheroid spreading dynamics, AND reproduce the cellular-Marangoni
direction predicted by Pajic-Lijakovic 2022 (low-γ → high-γ surface
flow)?**

The pilot is a **single simulation** producing per-frame R drift,
A/A₀, ∇γ field magnitude, spatial φ contrast, Marangoni power, plus
all inherited Stage 1c diagnostics.

## Bounded outcomes (4 buckets)

| Outcome | Definition | Action |
|---|---|---|
| **Bucket M1 — Layer 4 contributes measurably** | R drift < 0.138 (Stage 1c baseline) AND `<\|∇_s γ\|>_boundary > 0` (Marangoni active) AND `<φ>_contact > <φ>_interior` (spatial modulation working) AND no Bucket-M4 failure | STOP. Surface to PI as **layer-by-layer narrative quantitative pillar 6 (final)**: Layer 4 cellular Marangoni measurably contributes to spheroid spreading via φ-driven γ-gradient, completing the Stage 1 layer-by-layer mechanism table. PI decides whether to (i) advance to Stage 1d.b (nematic Q-tensor), (ii) re-baseline 3-phenotype with Layer 4, (iii) Stage 1e (radial-reduction comparison), (iv) consolidate the 16-simulation publication-strong narrative. **No automatic next-stage entry.** |
| **Bucket M2 — Layer 4 mechanism active but R drift unchanged** | `<\|∇_s γ\|>_boundary > 0` AND spatial φ contrast active, BUT R drift ≥ 0.138 within ±5% (Marangoni perturbative as predicted, doesn't materially shift R) | STOP. Surface to PI as **mechanism confirmation without dominant contribution**: Marangoni IS active per Pajic-Lijakovic 2022 mechanism but is too weak (Ma ≈ 0.06 < 1) to dominate R drift in this regime. This is consistent with the Stage 1d sanity-md prediction. PI decides whether to (i) report as expected-perturbative finding, (ii) increase γ_max-γ_min coupling (parameter sweep), (iii) test in different layer carrier (e.g. without Layer 2 to isolate). |
| **Bucket M3 — Layer 4 mechanism inactive (∇γ ≈ 0)** | `<\|∇_s γ\|>_boundary ≈ 0` OR spatial φ contrast not developing OR Marangoni impulse not influencing surface velocity | STOP. Surface to PI. Layer 3 spatial S_p extension OR γ(φ) scatter OR Marangoni impulse implementation has a bug; or the physics regime makes the effect negligible. PI decides among (i) debug spatial S_p, (ii) refine γ scatter / FD precision, (iii) report as "Marangoni below detectability in this regime". |
| **Bucket M4 — Critical failure** | NaN/Inf, max-speed runaway, mass-conservation break, ρ_osm bounds violation, φ ∉ [0,1], or spheroid disintegration | STOP **immediately**. Implementation issue (Marangoni sign wrong, ∇γ blowup, etc.). |

## Per-run gates

Per Stage 1d sanity-md:
- All Stage 1c inherited gates.
- Energy-monotone gate suspension condition EXTENDED to also fire when
  Layer 3 / 4 / 5 active (Cousin-Rule contract change resolves inherited
  Stage 1b issue).
- **R drift improvement vs Stage 1c baseline 0.138** (NEW GATE):
  strict-less than 0.138.
- **Marangoni power finite** (NEW GATE): same form as Stage 1a++
  active-power-finite (max |P_M| ≤ 10 · max(U_strain)).

## Diagnostics (no gate, log-only)

- Per-frame `<|∇_s γ|>_boundary`, `<|∇γ|>_grid_max`.
- Per-frame `<φ>_contact_band`, `<φ>_interior` (spatial contrast witness).
- Per-frame mean γ_eff at boundary band.
- All Stage 1c inherited diagnostics (ρ_osm, contact, shell, gravity,
  φ aggregates).

## Files of record

```
results/stage1d_pilot/{gate_report.md, metrics.csv, shell_profile.csv,
  contact_metrics.csv, snapshots.h5, ...}
```

## Auto-progression rules (PI 2026-04-29 full authorization, Stage 1d only)

**Auto-progression IS allowed within Stage 1d**:
- Sanity_md → outcomes → implementation → pytest → pilot → analysis →
  commit. All within one stage, automatic.

**Auto-STOP conditions**:
- All Stage 1c auto-STOP conditions inherited.
- Bucket M4 trigger.
- pytest fails before pilot.
- Marangoni power exceeds 10·U_strain (gate fail).

**Auto-STOP at end of Stage 1d**:
- After 1-pilot completes (or auto-STOP fires), STOP.
- Stage 1d.b (nematic Q-tensor) auto-entry FORBIDDEN.
- Stage 1e (radial reduction) auto-entry FORBIDDEN.
- Layer 4 refinement / phenotype re-baseline / publication consolidation
  requires PI input.

## Stop conditions (general; unchanged)

- Magic number 도입 금지.
- Gate semantic 변경 금지 beyond inherited Stage 1a++ contract changes
  + Stage 1d Cousin-Rule extensions explicitly recorded in sanity-md.
- v13 anti-pattern (measurement-protocol consistency 위반) 시 즉시 STOP.
- 한도 외 escalation 금지 — Stage 1d.b / 1e / Stage 2 자동 진입 금지.
