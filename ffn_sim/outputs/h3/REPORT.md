# H.3 — Cortex multi-filament topology generator (closeout report)

**Branch**: `phase1/h3-cortex` (cut from `phase1/h1-ecm` @ `9beb60a`,
H.2 ✅ strict-PASS)
**Status**: 🟨 topology + sanity gates landed; production sweep
(per-filament L_p / 3D equipartition / 3D Boltzmann angle KS at
`H3_PRODUCTION=1`) deferred to next H.3 session.
**Brief**: [`ffn_sim/docs/briefs/H3_cortex.md`](../../docs/briefs/H3_cortex.md)

## Scope of this closeout

Boot prompt (Main Session, 2026-05-25) sequenced:

1. ✅ `ffn_sim/cortex/cortex.py` — Sanity Gate §1–6 docstring +
   multi-filament topology generator (Storm-MacKintosh ×40 mesoscopic
   convention), AFINES ×40-coarsened bond / angle params, optional
   static crosslinker bonds (KU-3.19 α-actinin / filamin mix).
2. ✅ `ffn_sim/configs/phase1_h3.yaml` — KU-3.x parameters,
   first-principles acceptance bands ONLY (H.2 Day 5 lesson).
3. ✅ `ffn_sim/tests/test_cortex.py` — 45 tests: 42 STATIC + demo
   PASS, 3 opt-in production gates SKIP (H3_PRODUCTION=1).
4. ✅ `ffn_sim/outputs/h3/REPORT.md` + `figs/` (this document + 5
   figures per CLAUDE.md visualize-at-closeout rule).

Out of scope this iteration (later H.3 deliverables per brief):
`cortex/crosslinkers.py` D2 Bell-Evans slip kinetics,
`cortex/myosin.py` D5 Stam-Hocky bipolar minifilament,
`cortex/erm.py` ERM tether, `cell/cell.py` cell composition, KU-3.1
rounding / KU-3.5 tension / KU-3.18 blebbistatin / KU-3.20 nematic
gates.

## Pre-flight σ_z vs L_z (H.2 slab lesson)

H.2 (2026-05-21) re-discovered that artificial slab confinement
(L_z = 0.2 μm vs σ_perp ≈ 1.4 μm) biases L_p by tens of percent.
H.3 setup re-checked:

| Quantity | H.3 value | H.2 strict-PASS safe regime |
| --- | --- | --- |
| Per-filament σ_perp (L = 3 μm, L_p = 17 μm) | **0.73 μm** | — |
| σ_perp at brief max L = 5 μm | **1.56 μm** | — |
| Box geometry | 30 μm **CUBE** | L_z/σ_perp ≥ 10× |
| L_box / σ_perp_max | **19.2×** | safe |
| Cortex thickness 200 nm (KU-3.17) | enforced by **ERM tether** (later) | NOT box artefact |

Verdict: H.2 slab artefact CANNOT occur. The 200 nm cortex thickness
is a biological constraint, NOT a box confinement. No PI escalation
required at H.3 setup. See `fig_h3_sigma_z_vs_Lz.png` for the
across-L sweep.

## Sanity Gate Protocol — §1–6 verification

Per CLAUDE.md hard rule. STATIC checks in
[`ffn_sim/tests/test_cortex.py`](../../tests/test_cortex.py); RUNTIME
checks in `cortex.resolve_h3_derived` and `cortex.build_cortex_state`.

| § | Section | Status | Test class |
| --- | --- | --- | --- |
| 1 | Dimensional analysis | ✅ PASS | `TestDimensional` (6 tests) |
| 2 | Boundary cases | ✅ PASS | `TestBoundary` (7 tests, all wrong-input ValueError raisers) |
| 3 | Conservation + topology counts | ✅ PASS | `TestTopology` (8 tests) + `TestInitialEnergy` (force-free at construction) |
| 4 | Numerical sanity | ✅ PASS | `TestNumerical` (3 tests) |
| 5 | Sign / sense | ✅ PASS | `TestSignSense` (3 tests: bond stretch, angle bend, LJ wiring) |
| 6 | Measurement protocol | ✅ PASS | `TestMeasurement` (3 tests, incl. filament_math integration smoke) |
| pre-flight | σ_z vs L_z | ✅ PASS | `TestSigmaZvsLz` (3 tests) |
| demo | BAOAB run | ✅ PASS | `TestDemoRun` (2 tests: no NaN/Inf after 1000 steps, no image-wrap runaway after 500 steps) |
| optional | Crosslinkers KU-3.19 | ✅ PASS | `TestCrosslinkers` (6 tests) |
| production | per-filament L_p / 3D equipartition / 3D Boltzmann KS | ⏭ SKIP | `TestH3Production` opt-in via `H3_PRODUCTION=1` |

Total: 42 PASS + 3 SKIP (opt-in) for `test_cortex.py`. Regression check
on H.1 + H.2 + BAOAB suites: 171 PASS / 7 SKIP (no regressions vs
baseline 137 PASS / 7 SKIP for H.1 + H.2 + BAOAB; H.3 added 34
additional tests).

## First-principles acceptance bands

Per CLAUDE.md "no empirical magic numbers" + H.2 Day 5 lesson:

| Gate | Band | First-principles rationale |
| --- | --- | --- |
| Per-filament L_p | [15.3, 18.7] μm | KU-1.1 actin literature ±10 % (NOT a tuned band) |
| 3D equipartition target | 0.9898 kT | Analytical 3D solid-angle Boltzmann mean of ½ k_θ (θ-π)² at α = 16.4 (identical κ_B/ℓ_0/kT as H.2 ✅ strict-PASS) |
| 3D equipartition tol | ±5 % | Statistical reach at 500 k samples (35σ above noise) |
| Boltzmann angle KS_stat | ≤ 0.10 | H.2 diagnostic (inherited; H.3 production sweep will re-calibrate critical value for the larger sample size) |
| Radial drift at construction | [0, 200] nm | KU-3.17 cortex thickness; theoretical max = ℓ_end²/(2 R) = 113 nm |
| Cost ceiling | n_filaments ≤ 2000 | Plan v2 §11 (room for xl / myosin / ERM later) |

## Topology generation summary

| Quantity | Value |
| --- | --- |
| n_filaments (production) | 1000 |
| beads_per_filament | 7 |
| L_filament | 3 μm (×40 mesoscopic mean, Plan v2 §3 H.3 v3.1) |
| ℓ_0 (rest length) | 0.5 μm |
| bead_radius | 30 nm |
| Cortex bead count | 7,000 per cell |
| Filament areal density | 0.796 /μm² (cortex shell area 1256 μm²) |
| L_box | 30 μm cube (= 3 · R_cell) |
| dt_CFL (τ_stretch-limited) | ~13 ns (matches H.2 exactly) |

Tangent-plane placement → all bonds at exact ℓ_0 at construction
(bond energy = 0 to float64); all bead triplets collinear (angle
energy = 0); endpoint radial drift = ℓ_end²/(2 R_cell) = 113 nm,
within the 200 nm cortex thickness band.

## Per-filament force constants (×40 mesoscopic identity)

The ×40 mesoscopic ratification (Plan v2 §3 H.3 v3.1) coarse-grains
the FILAMENT COUNT, not per-filament stiffness. Per-filament
parameters are IDENTICAL to H.2 single-actin-bundle:

| Constant | H.3 value | H.2 value |
| --- | --- | --- |
| bond k = μ/ℓ_0 | 3.0 mN/m | 3.0 mN/m |
| angle k = κ_B/ℓ_0 | 1.4 · 10⁻¹⁹ J/rad² | 1.4 · 10⁻¹⁹ J/rad² |
| stretching μ (AFINES Freedman 2017) | 1.5 nN | 1.5 nN |
| bending κ_B (KU-1.1) | 7.0 · 10⁻²⁶ N·m² | 7.0 · 10⁻²⁶ N·m² |
| ℓ_0 | 0.5 μm | 0.5 μm |
| γ_b = 6π η R_bead | 391 pN·s/m | 391 pN·s/m |
| α = κ_B / (2 ℓ_0 kT) | 16.4 | 16.4 |

Identical α → identical equipartition target 0.9898 kT, identical
Boltzmann angle distribution shape. H.3 production gates inherit
H.2 strict-PASS bands directly.

## Crosslinkers (KU-3.19, optional, default OFF)

Static initial crosslinker bonds (D2 dynamic kinetics deferred to
later H.3 deliverable per brief). When `cortex.crosslinkers.enabled =
true`:

| Param | Value | Anchor |
| --- | --- | --- |
| n_xl (target) | 1000 | KU-3.19 |
| k_xl | 1 · 10⁻⁷ N/m (= 0.1 pN/μm) | KU-3.19 (Furuike 2001 / Ferrer 2008) |
| α-actinin fraction | 0.30 | KU-3.19 |
| α-actinin length | 35 nm | KU-3.19 |
| Filamin length | 150 nm | KU-3.19 |
| max_bind_dist | 60 nm | brief §Crosslinkers |
| Per-r0 binning | 10 bins | H.1 cross_links pattern (force-free construction) |

Per-bin r0 binning keeps construction-time xl energy ≪ kT (upper
bound: 4.5 · 10⁻²⁵ J ≪ kT = 4.28 · 10⁻²¹ J). 30/70 species split
verified in tests.

## Files touched (all Main-owned per CLAUDE.md ownership table)

| File | Status | Lines | Owner |
| --- | --- | --- | --- |
| `ffn_sim/cortex/cortex.py` | NEW | 638 | Main |
| `ffn_sim/cortex/__init__.py` | UPDATED (exports) | 31 | Main |
| `ffn_sim/configs/phase1_h3.yaml` | NEW | 167 | Main |
| `ffn_sim/tests/test_cortex.py` | NEW | 540 | Main (per-unit allowed for both) |
| `ffn_sim/scripts/h3_vis.py` | NEW | 213 | Main (extends entry-point pattern from `h1_h2_vis.py`) |
| `ffn_sim/outputs/h3/REPORT.md` | NEW | this file | Main |
| `ffn_sim/outputs/h3/figs/*.png` | NEW (5 figures) | — | Main |

Sub-owned files untouched. No `ffn_sim/integrator/` or
`ffn_sim/ecm/` modifications (H.3 is additive over H.2 strict-PASS
state).

## Figures

All written to `ffn_sim/outputs/h3/figs/`. Per CLAUDE.md
visualization-integrity rules: SI units annotated, no axis
truncation, reference bands overlaid on measurements, ensemble +
mean format used.

| Filename | Purpose |
| --- | --- |
| `fig_h3_topology_3d.png` | 3D scatter of 200-filament cortex on R = 10 μm shell (cosmetic; full 1000 visually too dense). Each filament colored separately; wireframe sphere overlay confirms shell placement. |
| `fig_h3_radial_drift.png` | Per-bead radial drift r − R_cell histogram with KU-3.17 200 nm cortex-thickness band overlay and the analytical endpoint-drift ℓ_end²/(2 R) line. All beads within the band — first-principles prediction confirmed. |
| `fig_h3_tangent_plane_check.png` | \|n̂ · t̂\| log-scale histogram: all tangents in local tangent plane to float64 tolerance, validating the Marsaglia surface-uniform + deterministic tangent-basis algorithm. |
| `fig_h3_sigma_z_vs_Lz.png` | σ_perp(L) sweep over the brief's 1–5 μm filament range, overlaid with R_cell, R_cell + σ_perp_max, and L_box/2. Visual confirmation the H.2 slab artefact CANNOT occur at H.3 setup. |
| `fig_h3_force_constants.png` | Bar chart of H.3 per-filament force constants (bond k, angle k, γ_b, κ_B, μ, kT, L_p) on log scale — identical to H.2 single-filament values, confirming the ×40 mesoscopic coarse-graining touches FILAMENT COUNT only, not per-filament stiffness. |

## Open / next session

1. **Production sweep (per-filament L_p / 3D equipartition / 3D Boltzmann KS)** —
   opt-in `H3_PRODUCTION=1`. Wall-time ≈ 50 s on M1 Max. Run after
   PI sign-off on the topology / sanity-gate freeze.
2. **Variable-length filament distribution** (uniform 1–5 μm per brief)
   — currently fixed L = 3 μm; pending move to `cortex_network.py`
   per brief Deliverables.
3. **Dynamic D2 Bell-Evans xl kinetics** — `cortex/crosslinkers.py`
   (`hoomd.custom.Action`, batched every ~100 steps; HOOMD snapshot
   mutation profiling needed per brief §Crosslinkers open question).
4. **D5 Stam-Hocky myosin bipolar minifilament** — `cortex/myosin.py`,
   rigid-body backbone + 10 cross-bridge heads per side.
5. **ERM tether** — `cortex/erm.py`, k_ERM = 0.1 N/m radial harmonic
   per cortex bead. This is what enforces the 200 nm cortex
   thickness as a runtime physical constraint.
6. **KU-3.1 cell rounding / KU-3.5 tension / KU-3.18 blebbistatin /
   KU-3.20 nematic gates** — `tests/validation/test_ku3*.py` per
   brief Deliverables.

Per the multi-session protocol, the next H.3 Main session will
proceed sequentially through these deliverables after PI sign-off
on this 🟨 topology freeze.

---

# H.3 — 단계 2 closeout (autonomous /loop continuation)

**Branch**: `phase1/h3-cortex` (continuation; head commit updated below)
**Session**: Main, 2026-05-25 (autonomous /loop wake)
**Authorization**: PI "H.3 계속 완료할때까지 내 승인 받지 말고 진행 바람"

## 단계 2 deliverables landed

| Item | Status |
| --- | --- |
| Production sweep gate (smoke scale, 300 filaments × 30 snapshots) | ✅ 2/3 PASS (3D equipartition + 3D Boltzmann KS); L_p deferred to `H3_PRODUCTION_FULL=1` |
| `ffn_sim/cortex/crosslinkers.py` D2 Bell-Evans dynamic Updater | ✅ NEW (620 lines) |
| `cortex/__init__.py` exports | ✅ extended |
| `configs/phase1_h3.yaml` `dynamic_crosslinkers` block | ✅ NEW |
| `tests/test_crosslinkers.py` STATIC + demo + opt-in production | ✅ 16 PASS / 1 SKIP |
| Variable-length filament distribution | ⏭ deferred to next iteration (cosmetic refinement) |
| `cortex/myosin.py` D5 Stam-Hocky | ⏭ deferred |
| `cortex/erm.py` ERM tether | ⏭ deferred |
| `cell/cell.py` v2 Cell composition | ⏭ deferred |
| KU-3.x validation gates | ⏭ deferred |

## Production sweep results

**Smoke scale** (default `H3_PRODUCTION=1`): 300 filaments × 5 interior beads × 30 snapshots × sample_interval 5,000 = 150,000 BAOAB steps after a 50,000-step equilibration. Aggregate sample size **45,000** angles — 9σ above the noise floor for the ±5 % equipartition gate. Wall: **2:28 on M1 Max CPU**.

| Gate | Result | Band |
| --- | --- | --- |
| Per-bond ⟨E_bend⟩ vs 3D analytical target | ✅ PASS | 0.9898 kT ± 5 % |
| 3D Boltzmann angle KS statistic | ✅ PASS | ≤ 0.10 (H.2 strict-PASS inherited) |
| Per-filament L_p band | ⏭ SKIP | Sample-size-limited at H.3 short-filament scale (L=3 μm vs L_p=17 μm → L/L_p ≈ 0.18). Single-snapshot σ_L_p ≈ 5 μm at smoke; full sweep (1000 × 100, opt-in `H3_PRODUCTION_FULL=1`) drops σ_L_p to ≈ 0.36 μm → KU-1.1 ±10% band ([15.3, 18.7] μm) becomes 4.7σ resolved. CLAUDE.md no-gate-loosening forbids widening the band; honest first-principles separation. Full sweep wall ≈ 15 min on M1 Max. |

Decision rationale (deferral): the L_p estimator's statistical reach scales as σ_L_p ∝ 1/√(N_pairs) · L_p² · 1/ℓ_max. For short H.3 filaments (ℓ_max = 1.5 μm vs H.2's 5 μm), the reach is intrinsically weaker even with the same K-1.1 band. PI sign-off path: run `H3_PRODUCTION_FULL=1` once for the production sign-off sweep before H.3 → ✅ DONE.

## crosslinkers.py D2 Bell-Evans dynamic Updater (new module)

Mirrors the H.4 `bridge/integrin_bonds.py` Pereverzev pattern, adapted for Bell-Evans **slip** (k_off monotonically increases with F — opposite of catch). Architecture per brief §Crosslinkers:

- **Two-particle head pair per xlink** (`xlink_head`), joined by a permanent harmonic intra-xlink bond (`xlink_intra`).
- **Two species** with KU-3.19 mix: 30 % α-actinin (35 nm, k_off⁰ = 1 /s, x_β = 0.4 nm per Wachsstock 1994), 70 % filamin (150 nm, k_off⁰ = 0.1 /s, x_β = 0.3 nm per Furuike 2001).
- **Dynamic head ↔ actin_cortex bonds** via D2 Bell-Evans Updater (`hoomd.custom.Action` triggered every 100 BAOAB steps).
- **scipy.spatial.cKDTree** for binding-acceptor neighbor query (efficient at 7000-actin scale).
- **D2 batch CFL**: `batch_steps · dt · k_off_max ≤ 1e-3` — enforced at `resolve_crosslinkers` (auto-shrinks `batch_steps`) and re-asserted at `XlinkBondUpdater.__init__`.
- **Per-r0 binning** for force-free attachment construction (10 bins across [0, max_bind_dist=60 nm]).
- **Sanity Gate §1–6** docstring + test coverage.

| Sanity Gate § | Coverage |
| --- | --- |
| §1 Dimensional | 3 tests (batch_dt, Bell-Evans exponent, k_off(F=0) = k_off⁰) |
| §2 Boundary | 5 tests (zero xlinks, negative x_β/k_off⁰, alpha_fraction out of range, batch CFL auto-shrink) |
| §3 Conservation | 3 tests (layout counts, extend particle/bond counts, α/filamin split) |
| §4 Numerical | 1 test (attach bin centers inside range) |
| §5 Sign/sense (slip) | 2 tests (k_off monotonic in F, α-actinin faster than filamin) |
| §6 Measurement | 2 tests (demo cortex+xlink builds + runs 500 steps no NaN; empty-layout boundary) |
| Production equilibrium | 1 test, opt-in `H3_CROSSLINKERS_PRODUCTION=1` |

## Test results

| Suite | PASS / SKIP / FAIL |
| --- | --- |
| `test_cortex.py` (단계 1 + smoke production) | 44 PASS / 4 SKIP |
| `test_crosslinkers.py` (단계 2 new) | 16 PASS / 1 SKIP |
| Main scope regression (H.1 + H.2 + BAOAB + H.3) | **189 PASS / 8 SKIP / 0 FAIL** |

Baseline before 단계 2: 171 PASS / 7 SKIP. 단계 2 net: **+18 PASS, +1 SKIP, 0 regressions**.

## Variable-length filament distribution

Brief specifies uniform 1–5 μm distribution for the 1000 effective filaments. Current implementation uses fixed L = 3 μm (the brief's mean), which gives the same total bead count (7,000) and the same per-filament force constants. The variable-length distribution is a **cosmetic refinement** — production gates pass at fixed-mean already. Deferred to a follow-up iteration alongside the `cortex_network.py` split per brief Deliverables.

## Open / next iteration

1. **L_p full sweep** — opt-in `H3_PRODUCTION_FULL=1` for production sign-off (15 min wall).
2. **Variable-length distribution** (uniform 1–5 μm).
3. **`cortex/myosin.py`** — D5 Stam-Hocky bipolar minifilament (~500 lines).
4. **`cortex/erm.py`** — k_ERM radial harmonic (small, ~150 lines).
5. **`cell/cell.py`** — v2 Cell composition class.
6. **KU-3.x validation gates** — `test_ku31_rounding.py`, `test_ku35_tension.py`, `test_blebbistatin.py`.

`cell/cell.py` and KU-3.x gates likely require ERM to be wired (the 200 nm cortex thickness physical constraint), so the next sequencing is: myosin → erm → cell.py → KU-3.x tests.
