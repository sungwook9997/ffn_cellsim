# Stage 1b — bounded outcomes (3-phenotype pilot, Bare/Pre/Lam4)

This document records the **PI-pre-authorized bounded-outcome decision tree**
for the Stage 1b 3-phenotype pilot (Layer 3 φ-ODE active, Bare/Pre/Lam4
initial-condition mapping). Written before runs execute.

The companion documents are:
- `docs/stage1b_layer3_sanity.md` — six-check Sanity Gate + Magic-Number
  Block + PI authorization record
- `docs/03_adhesion_dynamics.md` — pre-existing φ-ODE specification
- `docs/SESSION_HANDOFF.md` §"Stage 1b framing memo" — phenotype-not-
  substrate framing
- `docs/outcomes_stage1a_plus_plus.md` — Stage 1a++ Bucket B baseline
- `docs/12_validation.md` — Pillar 3 (PI experimental data is cross-check,
  not fit target)

## Mechanism question (one)

**Does φ-modulated active stress (driven by formation-environment phenotype
per Cho 2020) reproduce the inter-condition R-drift / A/A₀ ordering observed
in the PI's experimental Bare → Pre → Lam4 progression?**

The pilot is a **3-run sweep** (one per phenotype), producing three
(R drift, A/A₀(t), φ trajectory) trajectories. The cross-comparison
question is qualitative ordering, NOT numerical fit (per
`docs/12_validation.md` Pillar 3).

## Bounded outcomes (4 sweep meta-buckets)

The buckets are defined on the *phenotype ordering pattern* across the
three runs. PI experimental ordering of A/A₀_final: Bare (4.6–13.4) <
Pre (4.8–26.3) < Lam4 (8.0–33.1) — i.e. **higher laminin presentation
during formation → more spreading at long times**.

| Sweep meta-bucket | Definition (3-phenotype outcome) | Action |
|---|---|---|
| **Bucket I — qualitative ordering reproduced** | Simulation A/A₀ ordering matches PI: A/A₀(Bare) < A/A₀(Pre) < A/A₀(Lam4) (or equivalently R drift Bare > Pre > Lam4 if active stress is the dominant driver and Lam4 spreads more); all three runs complete normally (no auto-STOP triggered) | STOP. Surface to PI as the **layer-by-layer narrative quantitative pillar 4** finding: phenotype-dependent active stress (Cho 2020 mechanism) reproduces the experimental ordering. PI decides whether to (i) refine Layer 3 (e.g. activate K_cortex(φ) or γ_cc(φ) couplings for richer comparison), (ii) advance to Stage 1c (Layer 5 mechano-osmotic), (iii) extend to longer simulation duration to compare with PI's 24–96 hr A/A₀ data, (iv) report as publication-strong finding. **No automatic Stage 1c entry.** |
| **Bucket II — partial ordering / weak phenotype effect** | One or two pairwise orderings match PI; the third is reversed or indistinguishable (e.g. Bare ≈ Pre < Lam4, or Bare < Lam4 < Pre) | STOP. Surface to PI. Suggestive: φ → ζ coupling is correct in direction but quantitatively imprecise. PI decides among (i) refine the φ_initial values, (ii) activate a second φ-coupling channel (e.g. γ_cc(φ)), (iii) refine the ζ_min, ζ_max bounds, (iv) accept as a partial finding. |
| **Bucket III — ordering reversed or no phenotype effect** | All three runs give similar R drift / A/A₀ within ~5%, OR the ordering is opposite to PI experiment | STOP. Surface to PI. Indicates that φ → ζ is NOT the dominant phenotype channel; other φ-couplings (γ_cc, K_cortex) may be needed, or the Cho 2020 framework needs reinterpretation. PI decides among (i) try different φ-coupling channel, (ii) activate multiple couplings simultaneously, (iii) revisit the φ_initial mapping, (iv) accept the negative finding (φ → ζ alone insufficient for phenotype reproduction). |
| **Bucket IV — divergence / NaN / scheme failure** | Any of the three runs hits NaN/Inf, max-speed runaway, mass-conservation break, φ outside [0,1] bounds, active-power-finite gate failure, or any pre-existing scheme-correctness gate | STOP **immediately**. This is **Layer 3 implementation issue**, not a phenotype mechanism finding. Surface the failing gate + state of the run. Do not proceed to remaining phenotype runs in the sweep. |

## Per-run gates

Per Stage 1b sanity-md §"Per-run gates", inheriting all Stage 1a++ gates
plus 3 new Stage 1b additions: φ ∈ [0,1] invariant, φ trajectory monotone
toward predicted φ_eq, A/A₀ trajectory finite & non-pathological. Energy-
monotone gate stays SUSPENDED (Stage 1a++ contract change inherited).

## Diagnostics (no gate, log-only)

- Per-frame `<φ>_bulk`, `<φ>_boundary`, `<φ>_contact`, `min(φ)`, `max(φ)`,
  `std(φ)`
- Per-frame `A/A₀(t)`
- Per-frame mean-field-equivalent ζ (for cross-comparison with Stage
  1a++ point-ζ runs)
- Per-frame contact-band ρ_kernel (Stage 1a+ inherited)
- Per-frame shell `<ρ_kernel>(r/R₀)` profile (v15 inherited)

## Files of record (per phenotype)

```
results/stage1b_pilot_bare/{gate_report.md, metrics.csv, ...}
results/stage1b_pilot_pre/{gate_report.md, metrics.csv, ...}
results/stage1b_pilot_lam4/{gate_report.md, metrics.csv, ...}
```

Plus a sweep-level summary written to `results/stage1b_summary.md` after
the third phenotype run completes (3-phenotype response table + meta-bucket
classification + PI experimental-data cross-comparison overlay).

## Auto-progression rules (per PI 2026-04-29 full authorization)

**Auto-progression IS allowed within Stage 1b**:
- Bare run completes with no Bucket-IV failure → automatically proceed
  to Pre run.
- Pre run completes with no Bucket-IV failure → automatically proceed
  to Lam4 run.
- Lam4 run completes → automatically aggregate the response curve and
  classify the meta-bucket I / II / III.

**Auto-STOP conditions** (any of these halts the sweep immediately):
- All Stage 1a++ auto-STOP conditions inherited.
- φ exits [0, 1] bounds at any frame.
- φ trajectory diverges (does not relax toward predicted φ_eq within
  reasonable margin).
- pytest fails before the sweep starts.

**Auto-STOP at end of Stage 1b** (regardless of bucket):
- After 3-phenotype pilot completes (or auto-STOP fires), STOP.
- Stage 1c (Layer 5 mechano-osmotic) auto-entry is FORBIDDEN.
- Layer 3 refinement / Stage 1a++.b activation requires PI input.

## Stop conditions (general; unchanged)

- Magic number 도입 금지 (φ_initial values are explicitly PI-authorized
  estimates per the sanity-md, not magic numbers per the project
  convention).
- Gate semantic 변경 금지 beyond the inherited Stage 1a++ contract
  changes.
- v13 anti-pattern (measurement-protocol consistency 위반) 시 즉시 STOP.
- 한도 외 escalation 금지 — Stage 1c 자동 진입 금지.
