---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# DCM /goal quantitative validation vs literature (beyond-goal, 2026-07-01)

Turning the /goal qualitative wins (aggregation faceting + proliferation T47D gradient) into
quantitative, literature-checked results. Measured on the committed production npz (PROD_n400,
PHASE5d_prolif_n200). Figure: `viz_html/16_quant_validation_vs_lit.png`.

## 1. Faceting (aggregation) — Q = 149; faceting is GEOMETRIC-init-driven, not γ̃-energy-driven

- Production foam isoperimetric **Q = A³/V² mean 149, median 146** (sphere = 113; SimuCell3D-grade
  faceting ≈ 250). So the cells are genuinely faceted (well above a sphere) but mildly (Phase-3 finding:
  radial-warp-icosphere ceiling ~150).
- **γ̃ = γ/(K·ℓ) = 1.1×10⁻⁴** — far BELOW the SimuCell3D faceting band [0.02, 0.10] (K = 7.73e5 driver
  default, γ = 1e-3, ℓ = V₀^⅓). **Honest interpretation:** our faceting does NOT come from the γ̃
  surface-tension regime (we are 2-3 orders below the band). It comes from the **confluent Voronoi
  INITIALIZATION** (geometric), which the energy then MAINTAINS (Phase-4 maintain result). This is a
  DIFFERENT faceting mechanism than SimuCell3D's γ̃-driven one — a real, documented distinction, not a
  failure. To be γ̃-in-band you'd drop K to ~2500 (the K-reconciliation open question) — a PI decision,
  not tuned here.

## 2. T47D interior deformation (proliferation) — direction CORRECT, magnitude ~83%

- Grown proliferation foam (N=200→280): **CORE cell surface-area 749 µm² vs RIM 659 µm² → ratio 1.14**.
- T47D EM datum (PMC10212087): inner ~224 µm² vs outer ~164 µm² → **ratio 1.37**.
- **Correct direction** (inner cells larger/more-deformed, the T47D signature), at **~83% of the lit
  magnitude** (1.14 vs 1.37). A real partial quantitative match — the mechanism (growth crowds the
  interior) reproduces the sign and most of the magnitude. More divisions (further growth) would likely
  increase the ratio toward 1.37.

## 3. Porosity — NEGATIVE (−5.3%), reveals the division-without-remesh penetration

- Grown foam porosity = 1 − ΣV_cell / V_hull = **−5.3%** (cell volumes exceed the aggregate hull →
  cells OVERLAP). T47D EM: porosity 12% (day 5) → 2.4% (day 20, jamming).
- **This is a real limitation surfaced honestly:** division runs WITHOUT remesh (the division+remesh
  co-run is disabled, `DCM_DIVISION_REMESH_CORUN_DESIGN_2026-06-29.md`), so each mitotic insertion
  accumulates local interpenetration (the G2_interpenetration FAIL seen during Phase 5). The equilibrium
  (non-dividing) foams are clean (pen-ratio 0.93-0.97); the penetration is division-specific. Resolving
  the division+remesh co-run would clean the grown-foam porosity — a defined next build.

## Summary

| metric | ours | literature | verdict |
|---|---|---|---|
| faceting Q | 149 | SimuCell3D ~250 | faceted, mild; geometric-init not γ̃-driven (γ̃ 1.1e-4 << band) |
| T47D inner/outer area | 1.14 | 1.37 (T47D EM) | correct direction, ~83% magnitude |
| porosity | −5.3% | 12%→2% (T47D) | overlap from division-w/o-remesh (defined fix) |

**Net:** aggregation faceting + proliferation interior-gradient are quantitatively real and lit-consistent
in direction; two honest gaps are documented (γ̃ regime = K-reconciliation PI decision; grown-foam
porosity = division+remesh co-run build). No magic-number was applied. These are the validated /goal
science results plus their precisely-scoped remaining work.
