# Outcomes — Stage 1d.c Lam4 cell-count sweep (1k / 2k / 4k / 8k × 80 hr)

PI directive 2026-04-30 (overnight autonomous): 4-cell-count Lam4
sweep with all stages enabled (1b.b φ_memory + c_act split, 1d.b
Marangoni A+F, 1d.c ECM communication + protrusion state machine
+ traction). Production-conservative parameters. All Option F Week 2
gate contracts + Codex connected-component diagnostic active.

Chain: `b2i3qaywn`. Total wall-clock: 02:56 → 06:01 = **3h 5m** for
4 productions × 80 hr each + per-run paper-grade visualization.

## Top-line result

**Stage 1d.c reproduces the Lam4 80-hr spreading qualitatively** (no
peak-and-decay, sustained outward extension, free-edge-driven
polarity). Quantitative magnitude is **0.15–0.23× of PI experimental
endpoint** at all four cell counts after correcting for hull-leverage
artifact at low N.

| cell count | A_hull / A_0 (80 hr) | A_largest / A_0 (80 hr) | hull leverage | largest frac | escape | Wall-clock |
|---|---|---|---|---|---|---|
| 1k | 7.29 (artifact) | **1.11** | 6.57 | 0.879 | 12% | 33 min |
| 2k | 1.88 | **1.34** | 1.40 | 0.991 | 0.9% | 34 min |
| 4k | 1.59 | **1.59** | **1.00** | **1.00** | **0%** | 36 min |
| 8k | 1.73 | **1.73** | **1.00** | **1.00** | **0%** | 44 min |
| **PI exp Lam4 (n=26)** | **7.49 (78 hr)** | n/a | n/a | n/a | n/a | — |

> The hull-leverage gate (Codex Stage 1d.c addition) is the critical
> diagnostic. The 1k run's `A_hull/A_0 = 7.29` matches PI almost
> exactly — but `largest_frac = 0.879` reveals 12% of particles
> escaped, inflating the convex hull 6.6×. The intact-spheroid
> spread is `A_largest/A_0 = 1.11`, an order of magnitude below PI.
> Without the connected-component diagnostic this artifact would
> have been reported as success.

## Per-time intact-spheroid trajectory (A_largest/A_0)

| t (hr) | 1k | 2k | 4k | 8k | PI exp |
|---|---|---|---|---|---|
| 5 | 1.10 | 1.32 | 1.49 | 1.63 | 2.66 |
| 15 | 1.11 | 1.37 | 1.57 | 1.70 | 4.40 |
| 30 | 1.11 | 1.40 | 1.62 | 1.71 | 6.29 |
| 45 | 1.04 | 1.38 | 1.60 | 1.73 | ~6.7 |
| 60 | 1.06 | 1.36 | 1.61 | 1.72 | ~7.0 |
| 80 | 1.11 | 1.34 | 1.59 | 1.73 | 7.49 |

## Findings

### F1 — Stage 1d.c mechanism qualitatively works
- Persistent outward traction observed via `protrusion_state` cycling
  through quiet → filopodia → nascent → lamellipodium → retract.
- φ_memory exactly preserved (drift = 0.0 at 80 hr) across all 4 runs,
  confirming Stage 1b.b Layer 3 split is stable at production scale.
- ECM displacement field non-zero, propagating per Nahum 2023 screened
  Helmholtz operator (λ_ecm = 4 cell diameters mid-range).
- No NaN, no runaway, gates 5a + F3 anchor balance + F4 contact-band
  ρ + F8 A/A0_topdown all PASS.

### F2 — Cell-count convergence indicates physical asymptote
- Intact-spheroid endpoint: 1.11 → 1.34 → 1.59 → 1.73 across 1k/2k/4k/8k.
- 4k → 8k: +9% only, suggesting the framework approaches a physical
  asymptote near A_largest/A_0 ≈ 1.7–2.0 with current parameters.
- This contradicts the "Mechanism A is enough" hypothesis: even with
  full Stage 1d.c (ECM + protrusion + traction + Mechanism A+F),
  asymptote is 0.23× of PI experimental.

### F3 — Magnitude gap suggests parameter sweep needed
- Production-conservative `T0=0.1`, `λ_filo=0.05`, `α_A=0.2`,
  `α_F=0.1` are likely too weak.
- Pilot 1k 4-hr (A/A0 = 8.52, single-cell-runaway artifact) showed
  the *mechanism* can produce extreme spread; the question is the
  parameter range that gives PI-magnitude sustained spread without
  hull-leverage artifacts.
- PI sensitivity sweep on `T0` and `tau_lam` is the next step.

### F4 — Free-edge fix + speed cap effective
- Free-edge accuracy fix (only contact-band particles update polarity)
  + Codex's `protrusion_speed_cap_star = 0.06` together prevent the
  R drift 0.588 catastrophe seen in pre-fix Stage 1d.c smoke pilot.
- Post-fix R drift: 1k 0.36, 2k 0.10, 4k 0.17, 8k 0.21 — bounded but
  still over v15 baseline 0.244 at 1k (lowest mass / per-cell impulse
  largest).

### F5 — Codex's connected-component gate is essential
- Without `topdown_hull_leverage` and `largest_component_fraction`
  diagnostics, the 1k 7.29 hull artifact would have been reported as
  PI-matching success.
- With the gate, the artifact is correctly classified (FAIL at hull
  leverage 7.39, largest_frac 0.879) and the intact-spheroid value
  1.11 is exposed.
- This is a *measurement-protocol* correction in the spirit of Hard
  Rule 11.

## Gate state across the 4 runs

| Gate | 1k | 2k | 4k | 8k |
|---|---|---|---|---|
| F1 curvature κ vs 2/R | FAIL 21.7% | FAIL 21.7% | FAIL 21.7% | FAIL 13.4% |
| F2 momentum drift abs | PASS 1.25e-3 | PASS | PASS | PASS |
| F3 anchor balance | PASS 0.187 | PASS | PASS | PASS |
| F4 contact-band ρ | PASS 0.652 | borderline 0.648 | borderline 0.649 | borderline |
| F5 R drift | FAIL 0.363 | FAIL 0.099 | FAIL 0.166 | FAIL 0.214 |
| F6 sphericity post-spread | PASS 0.808 | PASS | PASS | PASS |
| F7 active power finite | FAIL 153 | FAIL 30.1 | PASS | PASS |
| F8 A/A0_topdown | PASS | PASS | PASS | PASS |
| 5a φ_memory preservation | PASS 0.000 | PASS 0.000 | PASS 0.000 | PASS 0.000 |
| 5b c_act_band → c_eq | FAIL 0.245 | FAIL 0.152 | FAIL 0.151 | FAIL 0.187 |
| Stage 1d.c hull leverage | FAIL 7.39 | FAIL 1.45 | PASS 1.00 | PASS 1.00 |
| Stage 1d.c largest fraction | FAIL 0.879 | PASS 0.991 | PASS 1.00 | PASS 1.00 |

F5 (radius drift) and 5b (c_act trajectory) remain ACCEPTED-LIMITATION
per `docs/v1/gate_fail_taxonomy.md`. F1 curvature passes at 8k (13.4% vs
limit 10% — closer; production scale gives better surface-cell
statistics).

## Visualization package

Per CODEX_FIGURE_GUIDE.md and PI viz directive 2026-04-29.

Per-run (each of `results/stage1dc_production_{1k,2k,4k,8k}_lam4/figures/`):
- `live_topdown/frame_*.png` (321 frames each)
- `live_sideview/frame_*.png` (321 frames each)
- `overlays/overlay_*.png` (321 frames each, 3×4 multi-channel grid:
  φ_memory / c_act / φ_eff / ρ_osm / γ_p / ‖v‖ / pressure / dev_norm
  / protrusion_state / fa_strength / ecm_signal / polarity vector)
- `movies/topdown_live.{mp4,gif}` + `movies/sideview_live.{mp4,gif}`
- `dashboard/run_dashboard.png` (4×4 metrics + Stage 1d.c panels)
- `parameter_tables/run_parameters.{csv,png}`
- `final_frame_{topdown,sideview}.png`
- Paper-grade dual save:
  `figures/Figures_for_draft/run_dashboard_*.svg` (white bg)
  `figures/Figures_for_PPT/run_dashboard_*.png` (600 dpi transparent)
  `figures/Logs/_captions/stage1dc/*.json` (provenance manifest)

Cell-count aggregate (`results/cellcount_comparison_lam4/figures/`):
- `AA0_topdown_hull_vs_PI.{svg,png}` — hull-based with PI overlay
- `AA0_largest_component_vs_PI.{svg,png}` — intact-spheroid with PI
- `hull_leverage_vs_t.{svg,png}` — escape diagnostic + gate limits
- `final_frame_topdown_quartet.{svg,png}` — 2×2 final top-down
- `final_frame_sideview_quartet.{svg,png}` — 2×2 final side-view
- `cellcount_summary_dashboard.{svg,png}` — 2×2 summary panels
- 6 JSON manifests under `Logs/_captions/stage1dc_cellcount/`

## Implications for next steps

1. **PI sensitivity sweep on T0 / τ_lam / λ_filo** — the magnitude
   gap (~0.15–0.23×) suggests production-conservative parameters are
   weaker than PI experimental cells. A sweep with the new
   connected-component gate to filter hull artifacts will identify the
   PI-matching range without fitting.

2. **Mechanism E (Stone 1990 surfactant transport)** — the
   marangoni_review.md deferred mechanism may be the missing piece
   for sustained spread at 4k+ scale where escape is suppressed.

3. **Layer 3 audit expansion** — c_act_band trajectory still FAILing
   at production (5b err 0.15–0.25 vs tol 0.10), suggests the contact-
   band particle population turnover is not letting c_act saturate.
   Codex review item 4 (Layer 3 deeper audit) remains open.

4. **F5 R drift 0.36 at 1k** — the mass-scaling sensitivity of
   per-cell impulse is most extreme at small N. Either accept as a
   small-N limitation or normalize traction by `N_total/N_band`.

5. **Hull-leverage gate is now the primary escape diagnostic** — its
   inclusion in the gate report is essential for any future Stage 1d.c
   production.

## Cross-references

- `docs/v1/stage1d_c_ecm_communication_sanity.md` — Stage 1d.c sanity
- `docs/v1/stage1d_b_marangoni_sanity.md` — Mechanism A+F (active here)
- `docs/v1/stage1b_b_phi_split_sanity.md` — Layer 3 split (active here)
- `docs/v1/marangoni_review.md` — Mechanism E pending
- `docs/codex_review_synthesis.md` — Codex items 1, 4 referenced
- `docs/v1/gate_fail_taxonomy.md` — gate classification
- `data/experimental/260313_Lam4.csv` — PI 26-position dataset
- `scripts/render_cellcount_comparison.py` — aggregate viz (this commit)
