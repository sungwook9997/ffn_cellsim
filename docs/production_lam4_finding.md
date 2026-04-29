# Production Lam4 finding — Bucket P3 (peak-and-decay; framework retracts in long-time vs PI [8, 33])

> **REVISION 2 2026-04-29 (Layer 3 audit)**: per PI directive
> `docs/layer3_phi_audit.md`, the peak-and-decay interpretation
> below is **provisionally contaminated** by artificial interior φ
> decay (the v11 spatial S=0 patch erases formation phenotype memory
> on the Cho 2020 transition timescale, contradicting the experimental
> setup where Bare/Pre/Lam4 is set at t=0). Final asymptote
> classification is deferred until Stage 1b.b Pilot 2 (5k Lam4 80hr
> with φ_memory + c_act split).
>
> **REVISION 1 2026-04-29 (Option F Week 1)**: original commit b506b57
> labeled the gate metric `contact_xy_hull/A₀` as "A/A₀_topdown" and
> reported "Bucket P2 plateau ≈ 2.65". This was a metric mislabel.
> The true `A_over_A0_topdown` (the publication-relevant top-down
> projection per Hard Rule 11) is **peak-and-decay**: peak 1.570 at
> 7.75 hr, end 1.409 at 80 hr. Reclassified to **Bucket P3 strict
> (mechanism-missing)**. See `docs/codex_review_synthesis.md` and
> `docs/marangoni_review.md` for the full reading.

Production Lam4 (commit `bkzdh05b1`, 5k cells × 80 sim hr = 480k steps,
wall-clock 23.46 min on Laptop A5000) completed. This document records
the finding and PI-diagnostic answer.

## TL;DR (revised)

**Pilot Phase 4-v2 Lam4 4 hr A_over_A0_topdown(end) = 2.657** (was
already past the spreading peak at this 4 hr endpoint — see "raw vs
peak" note below)

**Production Lam4 80 hr A_over_A0_topdown trajectory**:
- t* = 0 → 1.000
- t* = 465 (~7.75 hr) → **peak 1.570**
- t* = 4800 (~80 hr) → **end 1.409**

→ Framework exhibits **peak-and-decay** behavior. The continuum
reaches a maximum spreading at ~8 hr that is far below the PI Lam4
experimental range [8, 33], then *retracts* slowly to 1.409 by 80 hr.
The 5+1 layer continuum framework is **NOT time-budget-limited and is
NOT plateau** — it actively retracts in the long-time regime. Bucket
**P3 strict (mechanism-missing)** per `docs/production_lam4_outcomes.md`.

The previous "2.643" number is the **substrate `contact_xy_hull/A₀`
gate metric** (basal patch only, not all-particle top-down). It also
shows peak-and-decay (peak 2.643 at t* = 165 ≈ 2.75 hr, end 2.142),
but it is not the publication-relevant metric.

## Per-metric snapshot (Production Lam4 vs Pilot)

| Metric | Pilot Phase 4-v2 (4 hr) | Production (80 hr) | Δ |
|---|---|---|---|
| Wall-clock | ~1 min | 23.5 min | 23× (much faster than 2-3 hr estimate) |
| Steps | 24,000 | 480,000 | 20× (matches sim time scaling) |
| **A/A₀_topdown(end)** | 2.657 | **2.643** | **−0.5% (plateau)** |
| R drift | 0.073 | 0.151 | +106% |
| Sphericity (min) | 0.829 | 0.881 | +6% |
| φ_min / φ_max | 0.493 / 0.622 | 0.088 / 0.733 | spatial dispersion ↑ |
| <φ>(end) | 0.502 | 0.159 | interior S=0 → φ→0 |
| ρ_osm_max | 1.600 (saturate) | 1.224 | saturation released |
| ecm_strength(end) | 0.747 (mild) | **0.100 (floor hit)** | full degradation |
| mmp_total(end) | 1.36 | 152.99 | ~110× (proportional to cells × time) |
| n_boundary | 1000 (all) | 4978/5000 | similar |
| **Curvature κ vs 2/R** | 22% over (FAIL) | **9.5% (FIRST PASS)** | better statistics at 5k particles |
| A/A₀ ≥ 0.5 trajectory | 0.225 min (FAIL) | [1.0, 2.643] (PASS) | top-down measurement robust |

## Bucket P3 classification (revised)

`docs/production_lam4_outcomes.md` boundaries applied to the *correct*
top-down endpoint (1.409, not 2.643):
- P1 [4.0, 33.1]: time-budget cause; framework reaches PI range — **NO**
- P2 [2.5, 4.0]: partial reproduction — **NO** (1.409 < 2.5)
- P3 [< 2.5]: mechanism-missing strict — **YES** (1.409 < 2.5)
- P4: numerical failure — NO (no NaN, no runaway, all bounds preserved)

The framework **does not plateau** — it peaks at ~8 hr then retracts.
Phase 4-v2 4 hr endpoint 2.657 was on the metric variant that includes
the contact-patch envelope; the true top-down projection over 80 hr
shows monotone retraction after t* = 465. This is consistent with
incomplete Layer 4 Marangoni dynamics: the boundary φ saturates at
φ_eq ≈ 0.75 while interior φ decays to 0.088, generating a *contracting*
∇γ that the current implementation does not counter (no time-relaxation
γ as in Yadav 2022, no Stone 1990 surfactant transport, no Layer 5 ↔
Layer 4 coupling). See `docs/marangoni_review.md` Mechanisms A/E/F.

## Mechanism gap interpretation

The **A/A₀ asymptote gap** (~2.65 vs PI 8-33, factor of 3-12×) is
attributable to several missing mechanisms in the current 6-layer
continuum framework:

1. **ecm_strength floor hit (0.100)**: Layer 6 substrate degradation
   completed at production scale; further degradation impossible
   under the floor clamp. Yet additional spreading does not occur —
   indicating substrate weakening alone is not the primary spreading
   driver beyond the asymptote.

2. **φ interior decay** (φ_min → 0.088): Layer 3 spatial S_p extension
   drives interior particles toward φ→0 (S=0 there, k_- relaxation
   dominates). This produces a φ gradient driving Layer 4 Marangoni
   *retraction* (low-φ interior pulls high-φ boundary inward via
   γ-gradient flow). The retraction component is now significant
   (boundary φ stays at 0.733, near φ_eq=0.75), counteracting
   active-stress spreading.

3. **ρ_osm desaturation** (1.224 vs pilot 1.600 max): at long times,
   the spreading rate becomes slow enough that osmotic relaxation
   (β_osm) catches up with spreading-driven water efflux (α_osm·tr(ε̇)).
   Layer 5 stiffening contribution diminished compared to pilot.

4. **Active stress saturation**: with φ_eff stable near φ_eq = 0.75 at
   the boundary band, ζ_eff(φ) = 0.50 stable, but R drift increases
   from 0.073 → 0.151. Active stress is doing work at the boundary
   but the bulk does not continue to spread because the energy is
   dissipated through overdamped drag faster than active power injects.

5. **Continuum-only limitation (Stage 1a++.b deferred)**: the
   continuum-level Layer 2 boundary stress does not include the
   discrete stochastic events (lamellipodia, filopodia, leader cells,
   discrete focal adhesions) which produce intermittent extension
   pulses at the spreading edge. PI experimental MCF7 spheroids
   exhibit these discrete events; their absence in the continuum
   framework likely accounts for a substantial fraction of the
   A/A₀ asymptote gap.

## Curvature κ FIRST PASS

The 9.5% relative error (vs limit 10%) is the first PASS of the
curvature gate across all 27 simulations. Mechanism: with 5,000
particles spread over a substrate-spread shape, the colour-field
gradient is better resolved at the surface (838 surface cells vs 355
in pilot), reducing the f″/f′ off-peak amplification (the v13
anti-pattern that was structurally bounded at v12 by the ε² regulariser
in `_build_curvature` per `docs/stage1a_aha_div_sanity.md`).

## 27-simulation finalised narrative

| Stage | R drift | A/A₀_topdown(end) | Notes |
|---|---|---|---|
| (25 prior simulations) | 0.244 → 0.043 → 0.073 | various | layer-by-layer narrative |
| Phase 4-v2 Bare pilot 4hr | 0.131 | 1.746 | phenotype |
| Phase 4-v2 Pre pilot 4hr | 0.095 | 2.359 | phenotype |
| Phase 4-v2 Lam4 pilot 4hr | 0.073 | 2.657 | best phenotype |
| **Production Lam4 80hr (this run)** | **0.151** | **2.643** | **Bucket P2, plateau confirmed** |
| PI experimental Lam4 endpoint | n/a | 8.0–33.1 | 3-12× above framework asymptote |

## Implication for the publication narrative

The production result is **publication-relevant in two ways**:

1. **Quantitative validation of the layer-by-layer mechanism narrative**:
   Each layer (1-6) contributes a quantitative R drift / A/A₀
   reduction; the framework reproduces phenotype ordering (Bare < Pre
   < Lam4 in spreading) consistent with PI experimental data; the
   curvature operator is validated at production statistics; the
   measurement-protocol-matched A/A₀_topdown is robustly bounded.

2. **Quantitative negative finding on continuum-only framework**:
   Even at full PI experimental duration (80 hr) and production-scale
   particle count (5k), the 5+1 layer continuum framework asymptotes at
   A/A₀ ≈ 2.65, a factor of 3-12× below PI experimental [8, 33].
   This identifies the **continuum-only limitation** as the missing
   mechanism (most likely Stage 1a++.b stochastic boundary events:
   lamellipodia / filopodia / leader-cell discrete activations).

The publication can present the layer-by-layer narrative as the
**successful continuum result** AND the production asymptote as a
**quantified validity boundary** — analogous to the Stage 1e Sim A vs
Sim B comparison ("the radial reduction misses 3D physics; the 3D
continuum reaches an asymptote that misses discrete-events physics").

## Next-step PI decision options (no auto-entry)

1. **Stage 1a++.b** — activate Layer 2 stochastic events (lamellipodia
   /  filopodia / leader cells). Test whether discrete events close
   the asymptote gap.
2. **Multi-phenotype production** — Bare/Pre/Lam4 production scale to
   confirm asymptote behaviour is phenotype-independent (or
   phenotype-dependent gap).
3. **Paper draft** — consolidate 27-simulation narrative + production
   asymptote finding into manuscript outline.
4. **Sensitivity sweep** — vary φ_initial / ζ_max within Phase 2.2
   stable regime to map the asymptote sensitivity.

PI input required for next step. Auto-entry FORBIDDEN per Production
Lam4 outcome auto-STOP.
