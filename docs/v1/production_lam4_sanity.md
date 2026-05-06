# Production Lam4 — pre-execution Sanity Gate

This document is the pre-execution sanity review for the **first
production-scale run**: single Lam4 phenotype, 5,000 material points
× 80 sim hours (= 4,800 · τ_relax), full Stage 2 framework + Path C
+ Phase 4-v2 zeta_max=0.6 update + Phase 1 fixes. Estimated wall-clock
~2-3 hr on Laptop A5000 16 GB baseline. PI directive 2026-04-29:
A/A₀ absolute-value time-budget vs mechanism-missing diagnostic.

The companion documents are:
- `docs/v1/path_c_sanity.md` / `docs/v1/stage1c_sanity.md` /
  `docs/v1/stage1d_sanity.md` / `docs/v1/stage2_sanity.md` — inherited
  carriers' Sanity Gates
- `docs/v1/production_lam4_outcomes.md` — bounded-outcome decision tree
  (P1 time-budget / P2 partial / P3 mechanism-missing / P4 numerical)
- `configs/production_5k_80hr.yaml` — original production scaffold
  (used as Pre-baseline; Lam4 config derived with zeta_max=0.6 update)
- `docs/12_validation.md` — Sanity-Gate Protocol + Magic-Number Block

## Production scope (Lam4 single-phenotype)

### What is run

- **Single Lam4 phenotype**: φ_initial = 0.80 (Phase 4-v2 best phenotype:
  R drift 0.073, A/A₀_topdown 2.657 at pilot 4 hr).
- **5,000 material points** (5× pilot), **80 sim hours = 4800·τ_relax
  = 480,000 steps** at dt_star=0.01 (20× pilot duration).
- **Full Stage 2 framework**: Layer 1 (v15) + Stage 1a+ Option β α=1.0
  substrate + Stage 1a++ Layer 2 (overridden by Layer 3) + Stage 1b
  Layer 3 φ-ODE + Stage 1c Layer 5 mechano-osmotic + Stage 1d Layer 4
  Marangoni + Stage 2 Layer 6 chemistry/ECM + Path C effective gravity.
- **zeta_max = 0.6** (Phase 4-v2 update, Phase 2.2 first-PASS finding).
- **Phase 1 fixes inherited**: energy-monotone Cousin-Rule extension,
  anchor-force K_eff modulation, top-down A/A₀ measurement, Hard Rules
  10/11.

### Carrier parameter inheritance

All values from Phase 4-v2 (commit `708b44d`). No new physical
parameters introduced for production scale; only computational scale
parameters change:
- `n_material_points`: 1000 → 5000
- `total_sim_time_s`: 14,400 (4 hr) → 288,000 (80 hr)
- `total_time_star`: 240 → 4800
- `frame_interval_s`: 600 (10 min) → 900 (15 min — matches PI
  experimental imaging cadence per `CLAUDE.md`)
- `frame_interval_star`: 10 → 15
- `diagnostic_interval_steps`: 100 → 1000 (less-frequent diagnostics
  for 480k-step run; still ~480 diagnostic frames)
- `device_memory_GB`: 4 → 12 (A5000 16 GB ceiling)

### Time-scale mapping (per Stage 1d sanity-md check 0)

- τ_relax = 60 s (Moeendarbary 2013 Nat Mater IF 47 anchored).
- Pilot: 4 hr ≈ 5% of PI experiment 82 hr. PI A/A₀_final 8-33 for Lam4.
- Production: **80 hr ≈ 100% of PI experimental duration**. A/A₀
  trajectory directly comparable to PI `260313_Lam4.csv` time series.
- Sampling: 320 frames at 15 min/frame matches PI experimental imaging
  cadence.

---

## Sanity Gate (six checks, production-scale specific)

### 1. Dimensional analysis

All physical parameters inherited unchanged from Phase 4-v2; no new
non-dimensional numbers. Production-scale checks:

- **CFL stability** (Stage 1a v15 / Layer 5 / Layer 3 stiffness
  invariants) all hold at production scale (rates are physical
  constants per τ_relax; dt unchanged at 0.01·τ_relax).
- **Numerical noise accumulation** over 480,000 steps: f32 round-off
  ε ~ 10⁻⁷ per operation; momentum drift ~ √N · ε per step. At 480k
  steps + 5k particles, expected drift ~ √480k · 10⁻⁷ ≈ 10⁻⁵ —
  compatible with the 1e-3 momentum gate. Monitor.
- **VRAM ceiling (12 GB)**: pilot peaked at ~2 GB for 1k particles
  + 64³ grid. Linear in N: estimate 5×2 = 10 GB peak for 5k particles.
  Below 12 GB ceiling but tight margin. Monitor.

**Check 1: PASS** with explicit monitoring of drift accumulation +
VRAM peak.

### 2. Boundary cases

- **5k particles**: well above N ≥ 8 constructor minimum; well below
  the N → ∞ limit at 16 GB VRAM. Pilot validation (1k) extrapolates
  linearly per CLAUDE.md framework.
- **80 hr long-time run**: ~480k steps. Layer 6 ecm_strength clamp
  (≥ 0.1) prevents runaway full degradation at long times. Layer 5
  ρ_osm clamp ([0.5, 1.6]) prevents osmotic runaway. Layer 3 φ ∈ [0,1]
  preserved by ODE structure.
- **No new boundary cases** vs pilot.

**Check 2: PASS**.

### 3. Conservation invariants

Inherited from Stage 2 + Path C contracts:
- Mass exact (per-particle, no creation).
- Horizontal momentum gated; vertical absorbed by substrate.
- Energy-monotone gate SUSPENDED (Phase 1.1 Cousin-Rule extension).
- ecm_strength bounds invariant + mmp_total monotone-increasing
  (Stage 2 gates).
- ρ_osm bounds invariant (Stage 1c gate).
- φ ∈ [0,1] invariant (Stage 1b gate).

Long-time-specific monitoring:
- Mass drift watched per frame (gate tolerance 1e-10).
- Momentum drift accumulation reported (horizontal-only gate at 1e-3).

**Check 3: PASS**.

### 4. Numerical sanity

- dt_star = 0.01 unchanged; all stiffness invariants satisfied at
  pilot are preserved at production (rates are time-scale-invariant).
- Snapshot writing: 320 frames × 5k particles × ~50 bytes/particle ≈
  80 MB per HDF5; well within disk budget.
- 480k steps × ~3 ms/step (Stage 2 pilot mean) ≈ 24 minutes for
  pilot-equivalent compute. Production scale-up: 5× particles →
  ~5× ms/step → ~10-15 ms/step → 80-120 minutes. Combined with
  diagnostic / writer overhead: ~2-3 hr total. Matches PI estimate.

**Check 4: PASS**.

### 5. Sign / sense check

All inherited from Phase 4-v2. No new force terms; Lam4-specific
ζ_eff = 0.50 sits at Phase 2.2 ζ=0.50 stable-regime data point
(R drift 0.102 at pure Layer 2 baseline; with full Layer 4/5/6 + Path
C carrier reaches R drift 0.073 at pilot 4 hr). Sign chain:
- Layer 2 active stress (compressive on boundary cells) → outward push
- Layer 3 φ-ODE (substrate-engaged → Int-β1 → high ζ_eff)
- Layer 4 Marangoni (γ-gradient → tangential surface flow, retraction)
- Layer 5 K(ρ_osm) (volume contraction → water efflux → bulk stiffening)
- Layer 6 MMP (substrate degradation → reduced anchor over time)
- Path C gravity (constant downward body force → substrate anchoring)

**Check 5: PASS**.

### 6. Measurement-protocol consistency

Production-scale measurements match PI experimental modality (per Hard
Rule 11):
- **A/A₀_topdown** (PI Area_um2 column equivalent): xy-projection
  convex hull, z-independent. Direct comparison with PI's Lam4
  trajectory (A/A₀ 1 → 8-33 over 82 hr).
- **R**: effective_radius from particle second moment (3D shape
  measure). Cross-check with R_eff_xy = √(A_topdown / π) for
  consistency.
- **shell_profile.csv**: bulk transmission witness (carry-over from
  v15).
- **contact_metrics.csv**: substrate diagnostics (carry-over from
  Stage 1a+).

Long-time-specific protocols:
- Frame writing every 15 min sim time matches PI imaging cadence.
- 320 frames per run sufficient for trajectory shape + endpoint
  comparison.

**Check 6: PASS**.

---

## Magic-Number Block

**No new parameters introduced for production**. All carrier values
inherited from Phase 4-v2 (commit `708b44d`):
- γ_star, ξ_star (Layer 1 v15)
- gamma_sub_alpha (Stage 1a+ Option β)
- zeta_min, **zeta_max=0.6** (Phase 4-v2 update from Phase 2.2 finding)
- k_+, k_-, phi_initial=0.80 (Lam4 phenotype)
- alpha_osm, beta_osm, rho_osm bounds (Stage 1c)
- gamma_max, gamma_min (Stage 1d Marangoni)
- alpha_mmp, beta_deg, ecm bounds (Stage 2 Layer 6)
- gravity_star = 0.01 (Path C recalibrated)

All Magic-Number Block PARTIAL/PASS classifications inherited from
prior stages (ρ_floor, γ_sub_Col1, ζ_star Option α', α_osm,
gravity_star, α_MMP/β_deg, γ_max/γ_min — all framework-anchored with
honest disclosure per the verify pattern).

**Magic-Number Block PASS** (no new parameters; inherited).

---

## New gate / diagnostic candidates (production-scale specific)

| Gate | Status | Tolerance | Source |
|---|---|---|---|
| All Phase 4-v2 inherited gates | inherited | unchanged | Phase 4-v2 |
| **A/A₀_topdown(end) vs PI experimental Lam4 range [8, 33]** | NEW DIAGNOSTIC (no gate; bucketing per outcomes_lam4) | log only | PI cross-check (Pillar 3) |
| **Long-time mass drift** | NEW GATE | ≤ 1e-9 (relaxed from 1e-10 to allow 480k-step round-off) | first principles |
| **Long-time horizontal momentum drift** | inherited | unchanged (1e-3 with V_FLOOR) | Stage 1a+ |
| **VRAM peak** | inherited | ≤ 12 GB | CLAUDE.md ceiling |

---

## Implementation outline

- **NEW config**: `configs/production_lam4.yaml` derived from
  `configs/production_5k_80hr.yaml` with:
  - `run.name = production_lam4`
  - `run.output_dir = results/production_lam4`
  - `layer3.phi_initial = 0.80` (Lam4)
  - `layer3.zeta_max = 0.6` (Phase 4-v2 update)
  - `gate.mass_drift_rel_max = 1.0e-9` (relaxed for long-time)
  - All other parameters inherited from `production_5k_80hr.yaml`.
- **No code changes** to `acs/physics/mlsmpm.py` or `acs/runner.py` —
  Phase 1-5 fixes already in place.

---

## Decision request to PI (resolved 2026-04-29 by full authorization)

PI granted **full authorisation** 2026-04-29 ("Production scale 진입.
단일 phenotype Lam4 (best phenotype, ~2-3hr wall-clock background)")
for:
- Single Lam4 phenotype production run (no Bare/Pre concurrent).
- 5k particles × 80 sim hr scale (4800·τ_relax = 480k steps).
- zeta_max=0.6 update from Phase 4-v2.
- Frame interval 15 min (matches PI experimental imaging).
- Background execution + automatic 4-bucket result classification +
  commit + STOP.

Stop conditions remain in force: no magic numbers (no new parameters),
no gate semantics edits beyond inherited Phase 1.1 Cousin-Rule
extension, no v13 anti-pattern. Critical-error auto-STOP applies
(NaN/Inf, VRAM overflow, mass drift > 1e-9 across 480k steps,
horizontal momentum > gate). Other phenotype / stage / paper-draft
auto-entry FORBIDDEN; PI explicit decision required.
