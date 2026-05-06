# Stage 2 — bounded outcomes (Layer 6 chemistry/ECM single pilot, exploratory)

This document records the **PI-pre-authorized bounded-outcome decision tree**
for the Stage 2 Layer 6 chemistry/ECM remodeling pilot. Written before
runs execute. PI exploratory framing 2026-04-29: "Layer 6 돌리고 의미 확인".

The companion documents are:
- `docs/v1/stage2_sanity.md` — six-check Sanity Gate + Magic-Number Block
- `docs/v1/outcomes_stage1d.md` — Stage 1d Bucket M2 baseline (R drift 0.147)
- `docs/v1/outcomes_path_c.md` — Path C dimensional self-correction (relevant
  for substrate-weakening ↔ gravity rebalance question)

## Mechanism question (one)

**Does Layer 6 (MMP secretion → ECM degradation → time-varying γ_sub_eff)
on top of the Stage 1d full carrier introduce a measurable change in
spheroid spreading dynamics, AND reproduce the cancer-invasion framework
qualitative mechanism (substrate-engaged cells degrade substrate over time
→ adhesion weakens → spheroid less anchored)?**

The pilot is a **single simulation** producing per-frame R drift, A/A₀,
ecm_strength(t), mmp_total(t), plus all inherited Stage 1d diagnostics.

## Bounded outcomes (4 buckets)

| Outcome | Definition | Action |
|---|---|---|
| **Bucket Z1 — Layer 6 contributes meaningfully + cancer-invasion mechanism reproduced** | Stage 2 R drift differs from Stage 1d by ≥ 5% (either direction); ecm_strength(end) < 0.7 (substantial degradation); mmp_total(t) monotone-increasing throughout; A/A₀ trajectory reflects evolving substrate (e.g. mid-pilot spread acceleration as substrate weakens, OR late-pilot lift-off if degradation overshoots) | STOP. Surface to PI as **the project's 5+1 layer framework finale**: Layer 6 cancer-invasion mechanism reproduced quantitatively; Stage 1 layer-by-layer narrative complete with all 6 layers + Path C + Sim A/B comparison. PI decides whether to (i) Stage 2.b (de novo ECM secretion for full remodeling), (ii) production-scale runs (5k cells × 80 hr) for PI experimental endpoint comparison, (iii) consolidate the 18-simulation publication-strong narrative. **No automatic next-stage entry.** |
| **Bucket Z2 — Layer 6 mechanism active but R drift unchanged** | ecm_strength decreases as predicted, mmp_total accumulates, but R drift change < 5% (Layer 6 perturbative, doesn't dominate) | STOP. Surface to PI as **mechanism confirmation without dominant contribution**: Layer 6 is biologically active (MMP/ECM dynamics work), but in this regime + this duration, the mechanical effect on spreading is small. Consistent with Lu 2011 (ECM remodeling on 2-4 hr timescale; pilot at 4 hr captures only the initial degradation). PI decides among (i) extend simulation duration (production scale), (ii) increase α_MMP / β_deg rates for clearer effect, (iii) report as expected-perturbative finding. |
| **Bucket Z3 — Layer 6 mechanism inactive (no degradation)** | ecm_strength stays > 0.95 (no measurable degradation) OR mmp_total stays at 0 | STOP. Surface to PI. Layer 6 ODE implementation has a bug; OR the rates are too low; OR n_contact_band stays at 0 (mechanism never triggers — Stage 1a+ inherited substrate-anchoring issue persisted to this stage). PI decides among (i) debug Layer 6 ODE, (ii) increase rates, (iii) verify n_contact_band populated. |
| **Bucket Z4 — Critical failure** | NaN/Inf, max-speed runaway, mass-conservation break, ecm_strength bounds violation, runaway substrate degradation (ecm_strength → 0 negative or → +∞), or any pre-existing scheme-correctness gate fail by orders of magnitude | STOP **immediately**. Layer 6 implementation issue. Diagnose before proceeding. |

## Per-run gates

Per Stage 2 sanity-md:
- All Stage 1d inherited gates.
- Energy monotone SUSPENDED (Stage 1a++).
- **ecm_strength ∈ [ecm_strength_min, 1.0] invariant** (NEW): min ≥
  ecm_strength_min, max ≤ 1.0 across all frames.
- **mmp_total finite & non-decreasing** (NEW): mmp_total(end) finite
  (not NaN, not > 1e6); mmp_total(t+1) ≥ mmp_total(t) for all
  consecutive frames (within numerical tolerance 1e-9).

## Diagnostics (no gate, log-only)

- Per-frame `mmp_total`, `ecm_strength`.
- Per-frame γ_sub_eff (= γ_sub_star · ecm_strength).
- All Stage 1d inherited diagnostics (φ aggregates, ρ_osm, contact band,
  shell, gravity).

## Files of record

```
results/stage2_pilot/{gate_report.md, metrics.csv, shell_profile.csv,
  contact_metrics.csv, snapshots.h5, ...}
```

## Auto-progression rules (PI 2026-04-29 full authorization, Stage 2 only)

**Auto-progression IS allowed within Stage 2**:
- Sanity_md → outcomes → implementation → pytest → pilot → analysis →
  commit. All within one stage, automatic.

**Auto-STOP conditions**:
- All Stage 1d auto-STOP conditions inherited.
- Bucket Z4 trigger.
- pytest fails before pilot.
- ecm_strength bounds violation OR mmp_total monotonicity break.

**Auto-STOP at end of Stage 2**:
- After 1-pilot completes (or auto-STOP fires), STOP.
- Stage 2.b (de novo ECM secretion / drug effects) auto-entry FORBIDDEN.
- Inherited fixes (energy-monotone Cousin-Rule, anchor-force K_eff)
  auto-entry FORBIDDEN.
- Production-scale runs auto-entry FORBIDDEN.

## Stop conditions (general; unchanged)

- Magic number 도입 금지 (Layer 6 ODE constants are PI-authorized
  PARTIAL with explicit honest disclosure; not magic numbers per the
  verify pattern).
- Gate semantic 변경 금지 beyond inherited Stage 1a++ contract changes
  + the 2 new Stage 2 ODE invariant gates.
- v13 anti-pattern (measurement-protocol consistency 위반) 시 즉시 STOP.
- 한도 외 escalation 금지 — Stage 2.b / inherited fixes / production
  scale 자동 진입 금지.
