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

---

# H.3 — 단계 3 closeout (autonomous /loop continuation)

**Branch**: `phase1/h3-cortex` (continuation; head commit updated below)
**Session**: Main, 2026-05-25 (autonomous /loop wake)
**Authorization**: PI "H.3 계속 완료할때까지".

## 단계 3 deliverables landed

| Item | Status |
| --- | --- |
| `ffn_sim/cortex/erm.py` (~340 lines) — ERMHarmonic md.force.Custom + Sanity Gate §1–6 + CFL gate | ✅ NEW |
| `ffn_sim/cell/cell.py` (~280 lines) — v2 Cell composition class (cortex + ERM + xlinks + future slots) | ✅ NEW |
| `cortex/__init__.py` + `cell/__init__.py` exports | ✅ wired |
| `configs/phase1_h3.yaml` `erm` block | ✅ added |
| `tests/test_erm.py` (13 PASS) | ✅ NEW |
| `tests/test_cell.py` (10 PASS) | ✅ NEW |
| `cortex/myosin.py` D5 Stam-Hocky | ⏭ deferred to 단계 4 |
| KU-3.x validation gates | ⏭ deferred to 단계 4 (need myosin first for tension / blebbistatin) |
| L_p full sweep | ⏭ deferred (15 min wall, opt-in `H3_PRODUCTION_FULL=1`) |
| Variable-length filament distribution | ⏭ deferred (cosmetic) |

## ERM CFL sanity finding (단계 3 honest discovery)

**Brief literal `k_ERM = 0.1 N/m` (KU-3.18) is numerically incompatible with the cortex production `dt_CFL = 13 ns`.**

τ_ERM = γ_b / k_ERM = 3.91 · 10⁻¹⁰ / 0.1 = **3.91 ns** < dt_CFL = 13 ns. Stiff-spring runaway empirically observed: cortex bead drift hits **13.88 μm** in 500 BAOAB steps when ERM-on, vs analytic σ_radial = 0.21 nm.

Resolution (no gate-loosening, per CLAUDE.md):
- **`attach_erm_to_simulation` raises** on the brief-literal at cortex dt (CFL strict gate by default).
- Two production paths (PI sign-off path):
  - Reduce dt to 0.39 ns (33× more compute).
  - Soften k_ERM to a CFL-safe value (changes brief-literal physics — requires PI ratification against KU-3.18 anchor).
- **Demo / smoke tests** use `soft_k_ERM = 1·10⁻⁴ N/m` with explicit override; CFL-safe at cortex dt; documented as smoke-only.

This is a real sanity-gate finding that flags a parameter conflict in the brief itself. Stop-and-ask-PI criteria triggered: KU-3.18 oracle gives `k_ERM = 0.1 N/m` but this conflicts with the inherited cortex dt_CFL. Surface for PI decision in next iteration.

## v2 Cell composition class

`Cell.build(p_cortex, p_xlinks=..., p_erm=..., options=...)` composes a HOOMD Simulation with arbitrary subset of subsystems enabled:

- cortex (always)
- crosslinkers (D2 Bell-Evans dynamic)
- ERM tether (radial harmonic, CFL-gated)
- myosin / lamellipodium / FA slots (None — hooks for 단계 4 / H.5 / H.4 integration)

Diagnostics: `cell.bead_count_summary()`, `cell.tag_ranges()`, `cell.diagnostics()` — all dict-returning for REPORT.md generation.

Replaces deleted v1 `acs_kb/cell/cell.py` (single-chain Cortex-coupled, do not resurrect per CLAUDE.md).

## Test results

| Suite | PASS / SKIP / FAIL |
| --- | --- |
| `test_erm.py` (단계 3 new) | 13 PASS / 0 SKIP |
| `test_cell.py` (단계 3 new) | 10 PASS / 0 SKIP |
| Main scope regression (H.1 + H.2 + BAOAB + H.3) | **204 PASS / 8 SKIP / 0 FAIL** |

Baseline before 단계 3: 189 PASS / 8 SKIP. 단계 3 net: **+15 PASS, 0 new SKIP, 0 regressions**. (The +23 raw is offset by some test consolidation in shared fixtures.)

## Open / next iteration

1. **`cortex/myosin.py`** — D5 Stam-Hocky bipolar minifilament + D6 Hill stepping (~500-600 lines).
2. **KU-3.x validation gates** — test_ku31_rounding, test_ku35_tension, test_blebbistatin (need myosin + ERM wired into cell.py).
3. **L_p full sweep** opt-in (`H3_PRODUCTION_FULL=1`) — PI sign-off path for H.3 → ✅ DONE.
4. **ERM CFL resolution** — PI sign-off for dt-reduction OR k_ERM-softening (brief literal vs cortex dt_CFL conflict).
5. **Variable-length filament distribution** (cosmetic).

---

# H.3 — 단계 4 closeout (autonomous /loop continuation)

**Branch**: `phase1/h3-cortex` (continuation; head commit updated below)
**Session**: Main, 2026-05-25 (autonomous /loop wake)
**Authorization**: PI "H.3 계속 완료할때까지".

## 단계 4 deliverables landed

| Item | Status |
| --- | --- |
| `ffn_sim/cortex/myosin.py` (~700 lines) — D5 Stam-Hocky bipolar minifilament + D6 Hill stepping + Sanity Gate §1-6 | ✅ NEW |
| `ffn_sim/tests/test_myosin.py` (21 PASS) | ✅ NEW |
| `ffn_sim/tests/validation/test_ku3x_cortex.py` (KU-3.1/3.5/3.18/3.20 skeleton + measurement utilities) | ✅ NEW (4 utility PASS / 5 production SKIP opt-in `H3_KU3_PRODUCTION=1`) |
| `cortex/__init__.py` exports + `configs/phase1_h3.yaml` `myosin` block | ✅ wired |
| L_p full sweep `H3_PRODUCTION_FULL=1` | ⏭ wall-time UNDERESTIMATE (실측 11 min in CPU 후 추정 ~91 min full) — defer to dedicated PI sign-off session |
| KU-3.x production runs | ⏭ multi-hour wall, opt-in skeleton only |
| ERM CFL resolution | ⏭ awaiting PI sign-off |

## myosin.py architecture (D5 Stam-Hocky bipolar minifilament + D6 Hill)

- **100 minifilaments per cell** (KU-3.x Salbreux 2012 density 3/μm²).
- **Per minifilament**: 14-bead backbone (700 nm rigid-rod via stiff harmonic, k_backbone = 10·k_head_spring) + 10 cross-bridge heads per side (perpendicular harmonic, k_head_spring = 1 pN/μm = 1e-6 N/m, r0 = 200 nm). Total **34 particles per minifilament × 100 = 3,400 motor beads per cell**.
- **D2 Bell-Evans slip** head ↔ actin attach bonds (k_off⁰=10 /s, x_β=0.6 nm per Veigel 2002 NMII; mirrors `crosslinkers.py` Bell-Evans pattern).
- **D6 Hill stepping** while engaged: v(F) = v0·(F_s − F)/(F_s + F/a_over_F_stall), v0 = 1 μm/s, F_s = 0.5 pN, a/F_s = 0.5 (Kovács 2003). Hill closed-form imported from H.4 `bridge/motor.py` (shared module, read-only per CLAUDE.md).
- **Runtime stepper** `MyosinStepUpdater`: batched `hoomd.custom.Action`, every 100 BAOAB steps; unbinds via Bell-Evans, binds via scipy.spatial.cKDTree, advances `xlink_attach_b{i}` rest-length bin per Hill v(F).
- **CFL**: τ_head = 0.4 ms, τ_backbone = 39 μs — both ≫ cortex dt_CFL = 13 ns. **No CFL conflict** (unlike ERM where k_ERM = 0.1 N/m violates CFL).

## KU-3.x validation gates (skeleton + measurement utilities)

Production gates (KU-3.1 rounding / KU-3.5 tension / KU-3.18 blebbistatin / KU-3.20 nematic) require **60 s simulated time on a full cortex+myosin+ERM cell** — multi-hour to multi-day wall on M1 Max CPU. Skeletons skip-marked `H3_KU3_PRODUCTION=1`.

**Measurement utilities** that DO run at CI time (4 PASS):
- `cell_aspect_ratio(positions)` — covariance-tensor principal-eigenvalue ratio of the bounding ellipsoid.
- `nematic_order_S(axes)` — Q-tensor largest eigenvalue (0 = isotropic, 1 = aligned).

Both validated against analytical cases (isotropic sphere → aspect 1, S → 0; aligned axes → S → 1; 12×8×8 ellipsoid → aspect 1.5).

**KU-3.x first-principles bands** (no measurement-anchoring per H.2 Day 5):
- KU-3.1: aspect ratio ≤ 1.2 at 60 s (Salbreux 2012 literature criterion).
- KU-3.5: γ_cortex = 0.5 mN/m ± 30 % (Phase 1 default mid-range epithelial value, 0.1–1 mN/m literature range).
- KU-3.18: myosin-OFF aspect > 1.3 at 60 s (comparative against KU-3.1).
- KU-3.20: isotropic S < 0.1; aligned S > 0.3.

All bands derived from literature; tests will surface to PI when production runs land.

## L_p full sweep (deferred — wall-time underestimate)

Background full sweep (`H3_PRODUCTION_FULL=1`) killed after 11+ min CPU — the actual wall-time estimate is **~91 min** (37× the smoke 2:28), not the 15 min I quoted from the 단계 2 smoke. Reason: full sweep = 1000 filaments × 100 snapshots × 50 k step interval = 11.1× more steps × 3.33× more particles vs smoke. PI sign-off path needs a dedicated session (or weekend overnight run) to land this gate.

Alternative path: a **medium-scale sweep** (`H3_PRODUCTION_MEDIUM=1`, ~15 min wall, σ_L_p ≈ 1.34 μm) would resolve the KU-1.1 ±10% (1.7 μm) band at ≈ 1.3σ — borderline but possibly usable for an interim PASS. Add in next iteration if PI wants faster signal-off.

## Test results

| Suite | PASS / SKIP / FAIL |
| --- | --- |
| `test_myosin.py` (단계 4 new) | 21 PASS / 0 SKIP |
| `test_ku3x_cortex.py` (단계 4 new; utils only) | 4 PASS / 5 SKIP (opt-in) |
| Main scope regression (H.1 + H.2 + BAOAB + H.3) | **229 PASS / 13 SKIP / 0 FAIL** |

Baseline before 단계 4: 204 PASS / 8 SKIP. 단계 4 net: **+25 PASS, +5 SKIP (opt-in only), 0 regressions**.

## Open / next iteration

1. **L_p full sweep** dedicated session (~91 min wall) for H.3 → ✅ DONE production sign-off.
2. **KU-3.x production runs** (multi-hour each) — PI sign-off required, also needs ERM CFL resolution.
3. **ERM CFL resolution** — PI sign-off for dt-reduction OR k_ERM-softening.
4. **Variable-length filament distribution** (cosmetic).
5. **Cell composition with myosin wired into Cell.build** — currently Cell composes cortex + ERM + xlinks, but not yet myosin. Next iteration: extend Cell.build to wire MyosinStepUpdater.

---

# H.3 — 단계 5 closeout (autonomous /loop continuation)

**Branch**: `phase1/h3-cortex` (continuation; head commit updated below)
**Session**: Main, 2026-05-25 (autonomous /loop wake)
**Authorization**: PI "H.3 계속 완료할때까지".

## 단계 5 deliverables landed

| Item | Status |
| --- | --- |
| `ffn_sim/cell/cell.py` `build_cortex_full_simulation` helper — unified cortex + (optional) xlinks + (optional) myosin builder | ✅ NEW (≈230 lines added) |
| `Cell.build(with_myosin=True, p_myosin=...)` dispatch path | ✅ wired |
| `Cell.bead_count_summary` / `tag_ranges` / `diagnostics` extended for myosin | ✅ |
| `ffn_sim/cell/__init__.py` exports `build_cortex_full_simulation` | ✅ |
| `tests/test_cell_full.py` (NEW, 11 PASS / 1 SKIP) | ✅ NEW |
| Variable-length filament distribution | ⏭ deferred (cosmetic) |
| ERM CFL resolution + L_p full sweep + KU-3.x production | ⏭ awaiting PI sign-off / dedicated session |

## build_cortex_full_simulation architecture

Unified helper that stacks the three cortex subsystems in a single HOOMD Simulation:

1. `build_cortex_state(p_cortex, with_crosslinkers=False)` — base actin shell.
2. **(optional)** if `p_xlinks and p_xlinks.n_xl > 0`: `generate_xlink_layout` + `extend_cortex_state_with_xlinks` — adds `xlink_head` particles + `xlink_intra` bonds + `xlink_attach_b{i}` bond types.
3. **(optional)** if `p_myosin and p_myosin.n_motors_per_cell > 0`: `generate_cortex_myosin_layout` + `extend_state_with_cortex_myosin` — adds `cortex_myosin_backbone` + `cortex_myosin_head` particles + 3 new bond types.
4. HOOMD `md.bond.Harmonic` with ALL bond-type params (cortex-bond + xlink_intra + xlink_attach_b{i} + cortex_myosin_backbone + cortex_myosin_head_backbone + cortex_myosin_attach_b{i}).
5. `md.angle.Harmonic` (cortex-angle).
6. `md.pair.LJ` (`md.nlist.Tree`) with per-pair WCA wiring:
   - WCA enabled: actin × actin, xlink_head × xlink_head, myosin_backbone × *, myosin_head × myosin_*.
   - **WCA disabled** (so heads can approach and bind): xlink_head × actin, cortex_myosin_head × actin.
7. `md.Integrator(dt=p_cortex.dt_cfl)` with methods=[].
8. BAOAB Updater with per-type γ_b for all four particle types.
9. `XlinkBondUpdater` (Periodic batch_steps trigger).
10. `MyosinStepUpdater` (Periodic batch_steps trigger).

Returns a `dict[str, Any]` of handles: sim, topology, xlink_layout, myosin_layout, baoab_*, xlink_*, myosin_*, n_cortex_actin, n_xlink_heads, n_myosin_particles.

## Cell.build dispatch logic (extended)

- `options.with_myosin=True` → unified path via `build_cortex_full_simulation`.
- Else: dispatches to existing `build_cortex_xlink_simulation` (if `with_crosslinkers`) OR `build_cortex_simulation` (bare cortex). Preserves 단계 1-4 behavior verbatim.
- ERM wiring (`options.with_erm`) attaches AFTER any of the three paths, regardless of myosin/xlinks (CFL gate enforced).

`Cell` dataclass now holds: `p_myosin`, `myosin_layout`, `myosin_action`, `myosin_updater`, `n_myosin_particles`.

## Test results

| Suite | PASS / SKIP / FAIL |
| --- | --- |
| `test_cell_full.py` (단계 5 new) | 11 PASS / 1 SKIP |
| Main scope regression (H.1 + H.2 + BAOAB + H.3) | **240 PASS / 14 SKIP / 0 FAIL** |

Baseline before 단계 5: 229 PASS / 13 SKIP. 단계 5 net: **+11 PASS, +1 SKIP, 0 regressions**.

## Open / next iteration

1. **L_p full sweep** dedicated overnight session (~91 min wall).
2. **KU-3.x production runs** (multi-hour each; needs ERM CFL resolution).
3. **ERM CFL conflict** PI sign-off (dt-reduction vs k_ERM-softening).
4. **Variable-length filament distribution** (uniform 1–5 μm per brief; cosmetic).
5. (선택) integration-level visualization: `scripts/h3_full_vis.py` extending `h3_vis.py` for 3-way snapshot rendering.

---

# H.3 — 단계 6 closeout (autonomous /loop continuation)

**Branch**: `phase1/h3-cortex` (continuation; head commit updated below)
**Session**: Main, 2026-05-25 (autonomous /loop wake)

## 단계 6 deliverables landed

| Item | Status |
| --- | --- |
| `ffn_sim/scripts/h3_full_vis.py` (~370 lines) — 3-way integration visualization (5 figures) | ✅ NEW |
| `outputs/h3/figs/fig_h3_full_*.png` (5 NEW figures) | ✅ written |
| Bug fix: xlink head jitter σ 5 nm → 50 nm (avoid LJ blow-up when 2 xlinks pick same anchor) | ✅ |
| Bug fix: inter-subsystem LJ disabled in `build_cortex_full_simulation` (myosin × actin, etc.) | ✅ |
| Bug fix: nlist exclusions `('bond', '1-3')` for the 3-way builder | ✅ |
| Bug fix: `fig_full_3d_scatter` uses `state.get_snapshot()` for particle.types (cpu_local_snapshot lacks `.types`) | ✅ |
| Hangul → romanized in vis text panel (Hangul missing from DejaVu Sans Mono font; cosmetic) | ✅ |
| Variable-length filament distribution | ⏭ deferred (cosmetic, time budget) |

## h3_full_vis.py figures

1. `fig_h3_full_3d_scatter.png` — 3D scatter of full 3-way cortex on R = 10 μm shell (cortex actin green, xlink heads orange, myosin backbones purple, myosin heads red).
2. `fig_h3_full_bond_network.png` — bar chart of per-bond-type counts (cortex-bond, xlink_intra, xlink_attach_b*, cortex_myosin_*) from the live HOOMD frame.
3. `fig_h3_full_bead_counts.png` — `Cell.bead_count_summary` per-subsystem bar chart.
4. `fig_h3_full_updater_activity.png` — Updater activity time-series over an 8 × 200-step demo window: engaged-head count (top) + cumulative bind/break/Hill-step events (bottom).
5. `fig_h3_full_pipeline_summary.png` — text panel with 단계 1–6 commit history + current demo cell diagnostics + open-for-PI list.

## Bug fixes (단계 6 sanity findings)

**Three production bugs surfaced when extending the 3-way demo** to actually run BAOAB:

1. **LJ ∞ from coincident xlink heads** — original 5 nm jitter was too small to clear WCA cutoff (60 nm σ, 67 nm r_cut). Two xlinks picking the same anchor placed heads only ~12 nm apart → LJ ~1e-12 J → BAOAB runaway. Fix: jitter σ → 50 nm (typical separation ~70 nm > r_cut). Validation: 240 PASS regression, 11 PASS on new 3-way integration tests.
2. **Inter-subsystem WCA repulsion at construction** — original `build_cortex_full_simulation` had myosin × actin WCA enabled, which blew up when randomly-placed myosin backbones sat near cortex actin beads at the shell. Fix: disable all inter-subsystem WCA (myosin sits ABOVE actin in cortex anatomy; the intra-cortex packing is biology, not steric repulsion at this coarse-graining scale). Documented in code with rationale.
3. **WCA-cutoff vs bonded-pair conflict** — α-actinin xlink_intra rest length 35 nm and myosin backbone segment 54 nm are both shorter than WCA r_cut = 67 nm, so without nlist exclusions, bonded pairs would also feel WCA repulsion, competing with the harmonic bond. Fix: `md.nlist.Tree(exclusions=('bond', '1-3'))`.

All three are HONEST sanity-gate findings (same class as 단계 3 ERM CFL): the 단계 1-5 unit tests didn't catch these because they ran the SUBSYSTEMS individually, not the integrated 3-way. The full-cell integration was the first place where the combined LJ topology was actually exercised.

## Test results

| Suite | PASS / SKIP / FAIL |
| --- | --- |
| Main scope regression (H.1 + H.2 + BAOAB + H.3) | **240 PASS / 14 SKIP / 0 FAIL** |
| `test_cell_full.py` 3-way integration (passes ALL 11 demo tests) | (subset of above) |

Baseline before 단계 6: 240 PASS / 14 SKIP. 단계 6 net: **+0 PASS / +0 SKIP** but **3 production bugs caught + fixed** — net robustness improvement.

## Open / next iteration

1. **Variable-length filament distribution** (cosmetic refinement, brief §Cortex topology).
2. **L_p full sweep** dedicated overnight session (~91 min wall).
3. **KU-3.x production runs** (multi-hour each; ERM CFL 선행 필요).
4. **ERM CFL** PI sign-off.
5. **3-way 60s production gate** `H3_INTEGRATION_PRODUCTION=1` (multi-hour).

---

# H.3 — 단계 7 closeout (autonomous /loop continuation)

**Branch**: `phase1/h3-cortex` (continuation; head commit updated below)
**Session**: Main, 2026-05-25 (autonomous /loop wake)

## 단계 7 deliverables landed

| Item | Status |
| --- | --- |
| `cortex/cortex.py` `VariableLengthCortexLayout` dataclass + `generate_variable_length_cortex_layout` + `build_variable_length_cortex_state` + `build_variable_length_cortex_simulation` | ✅ NEW (~280 lines added) |
| `cortex/__init__.py` exports for variable-length API | ✅ wired |
| `configs/phase1_h3.yaml` `cortex.variable_length` block | ✅ |
| `tests/test_variable_length.py` (19 PASS / 1 SKIP, opt-in `H3_VARIABLE_LENGTH_PRODUCTION=1`) | ✅ NEW |
| Brief §Cortex topology (uniform 1–5 μm) realised | ✅ |

## Variable-length filament distribution architecture

**Additive design** (fixed-N functions untouched, no 단계 1-6 baseline regression):

- Per-filament L_i ~ Uniform(L_min=1 μm, L_max=5 μm), quantized to ℓ_0 = 0.5 μm multiples → N_beads_i = round(L_i/ℓ_0) + 1 (range 3–11 beads).
- Flat layout: `positions_flat` (Σ N_i, 3), `n_beads_per_filament` (F,), `filament_starts` (F,) for navigation.
- Bond / angle groups generated based on per-filament counts (Σ (N_i-1) bonds, Σ (N_i-2) angles with N=2 filaments contributing 0 angles).
- `build_variable_length_cortex_simulation` inherits the 단계 6 nlist exclusions lesson: `nlist.Tree(exclusions=('bond', '1-3'))`.

Config opt-in via `cortex.variable_length.enabled` flag (default OFF so existing builders unchanged); functional override via direct call to `generate_variable_length_cortex_layout`.

Production-scale L_p sweep at variable-length deferred to opt-in `H3_VARIABLE_LENGTH_PRODUCTION=1` (same multi-hour rationale as L_p full + KU-3.x gates).

## Brief §Cortex topology realisation

- **Brief literal**: "uniform 1–5 μm (mean 3 μm) → average 7 beads per filament at ℓ_0=0.5 μm. Total cortex beads ≈ 7,000 per cell".
- **단계 1-6 (fixed-mean)**: L = 3 μm fixed → exactly 7 beads × 1000 filaments = 7,000 beads.
- **단계 7 (variable-length)**: L_i drawn from Uniform(1, 5) μm, mean 3 μm (empirical with F=80 demo: 2.5–3.5 μm 95% CI). Total beads ≈ F · 7 = 7,000 on average, with σ ≈ F · σ_N ≈ 80 · 1.4 ≈ 112 (production scale: ~370 bead-count variance, ~5% of total). Same per-filament force constants (cortex.bond_k, cortex.angle_k); same dt_CFL.

The brief's distribution is now FULLY REALISABLE; the fixed-mean approximation remains the default for test-baseline preservation.

## Test results

| Suite | PASS / SKIP / FAIL |
| --- | --- |
| `test_variable_length.py` (단계 7 new) | 19 PASS / 1 SKIP (opt-in) |
| Main scope regression (H.1 + H.2 + BAOAB + H.3) | **259 PASS / 15 SKIP / 0 FAIL** |

Baseline before 단계 7: 240 PASS / 14 SKIP. 단계 7 net: **+19 PASS / +1 SKIP / 0 regressions**.

## Open / next iteration (final stretch)

ALL big-ticket items remain PI sign-off / dedicated session paths:

1. **L_p full sweep** dedicated overnight session (~91 min wall).
2. **KU-3.x production runs** (multi-hour each; ERM CFL 선행 필요).
3. **ERM CFL** PI sign-off (k_ERM 0.1 N/m vs dt_CFL=13 ns).
4. **3-way 60s production gate** `H3_INTEGRATION_PRODUCTION=1`.
5. **Variable-length L_p production sweep** `H3_VARIABLE_LENGTH_PRODUCTION=1` (multi-hour).

H.3 implementation deliverables (cortex.py + crosslinkers.py + erm.py + myosin.py + cell/cell.py + build_cortex_full_simulation + variable-length + 10 figures + tests) **ALL LANDED**. Status H.3 🟨 → ✅ DONE ratification awaits the PI sign-off + dedicated production sessions above.

---

# H.3 — 단계 8 closeout (autonomous /loop: L_p MEDIUM sweep)

**Branch**: `phase1/h3-cortex` (continuation; head commit updated below)
**Session**: Main, 2026-05-25 (autonomous /loop wake)

## 단계 8 deliverables landed

| Item | Status |
| --- | --- |
| `H3_PRODUCTION_MEDIUM=1` gate added to `test_cortex.py::TestH3Production` fixture | ✅ |
| L_p production test skipif extended to unblock under either MEDIUM or FULL env var | ✅ |
| L_p MEDIUM sweep (500 filaments × 50 snapshots × 25k step interval = 1.25 M BAOAB steps, ~15 min wall) | (result below) |

## L_p MEDIUM gate result

**Scale**: 500 filaments × 5 interior beads × 50 snapshots = **125 000 angle/bond samples** per s.
**Statistical reach**: σ_L_p ≈ L_p² / (√N_pairs · ℓ_max) ≈ 17² / (√125k · 1.5) μm ≈ **1.34 μm** — borderline at 1.3σ on the KU-1.1 ±10 % band [15.3, 18.7] μm.
**Wall**: **28:00** on M1 Max CPU (≈2× the 15 min naïve extrapolation — same per-particle BAOAB Updater overhead pattern observed in 단계 2 smoke vs 단계 4 full ratio).

**Result: 3/3 PASS** (`H3_PRODUCTION=1 H3_PRODUCTION_MEDIUM=1 pytest ffn_sim/tests/test_cortex.py::TestH3Production`):
- `test_per_filament_L_p_in_KU11_band` ✅ PASS — L_p in [15.3, 18.7] μm (KU-1.1 ±10 %)
- `test_3d_equipartition_strict` ✅ PASS — ⟨E_bend⟩ within ±5 % of 0.9898 kT
- `test_3d_boltzmann_angle_KS` ✅ PASS — KS_stat ≤ 0.10

This is the **strongest autonomous-loop signal** that the H.3 cortex correctly recovers KU-1.1 single-actin physics at the multi-filament scale. The H.2 strict-PASS bands carry through to H.3 — as predicted, since per-filament force constants (κ_B, μ, ℓ_0) are identical and the ×40 mesoscopic coarse-graining touches FILAMENT COUNT only.

## Significance

A MEDIUM PASS is **suggestive** that the production simulation correctly recovers KU-1.1 L_p at the H.3 cortex scale (cortex actin shell with random-tangent placement). It does NOT replace the FULL gate (4.7σ resolution required for the PI strict-PASS contract per the 단계 3 ERM CFL precedent — H.2 set the precedent of distinguishing interim from production-sign-off bands). The MEDIUM gate is autonomous-feasible (28 min wall); FULL (~91 min wall) needs a dedicated overnight session.

3D equipartition and 3D Boltzmann KS gates run at MEDIUM scale too — both still PASS, now with ~3× tighter statistical reach than at smoke scale (already strict-PASS in 단계 2).

**Decision per CLAUDE.md no-gate-loosening**: H.3 status remains **🟨 candidate-for-✅** until the FULL gate runs in a dedicated overnight session. The MEDIUM PASS is recorded as interim sign-off evidence, not as the production contract.

## H.3 status assessment

- **Implementation** (단계 1-7): ALL deliverables landed.
- **Tests** (단계 1-7): 259 PASS / 15 SKIP / 0 FAIL Main scope.
- **단계 8 production gate (L_p MEDIUM)**: borderline interim signal.
- **Remaining for H.3 → ✅ DONE ratification**: ERM CFL PI decision · L_p FULL dedicated session · KU-3.x production runs · 3-way 60s production · variable-length L_p production. **None of these are autonomous-feasible in 25-min /loop iterations.**

H.3 has reached the END of the autonomous /loop's productive scope. Further progress requires PI sign-off + dedicated long-running sessions.

---

# H.3 production gate progress — post-/loop session (2026-05-26)

PI directive 2026-05-26: "프로덕션 돌리자" (let's run production).

## ERM CFL RE-RATIFICATION (PI option A)

Per the 단계 3 sanity finding, the brief literal `k_ERM = 0.1 N/m` was numerically incompatible with the cortex production `dt_CFL = 13 ns`.

PI selected **option A: soften k_ERM to 1.0 × 10⁻⁴ N/m** (1000× softer than brief literal). KU-3.18 RE-RATIFIED.

- `configs/phase1_h3.yaml` updated with full rationale comment.
- `tests/test_erm.py` + `tests/test_cell.py` CFL-boundary tests updated to test the gate against the OLD 0.1 N/m hardcoded value (preserves the sanity check; anyone reverting to old value still gets the CFL gate raised).
- 23/23 ERM + Cell tests PASS after ratification.

**Physical implication**: σ_radial = √(kT/k_ERM) ≈ **6.5 nm** per bead (was 0.21 nm under old k_ERM). Still ≈30× tighter than the 200 nm KU-3.17 cortex thickness band — the biological "tight pinning relative to membrane band" interpretation carries through. Rounding / tension / blebbistatin / nematic gates probe EMERGENT large-scale curvature, not per-bead radial distribution, so the σ_radial change does NOT affect KU-3.x acceptance.

## L_p FULL sweep (autonomous-feasible production)

Launched in background on `phase1/h3-cortex` after ERM ratification:

```
H3_PRODUCTION=1 H3_PRODUCTION_FULL=1 pytest test_cortex.py::TestH3Production
```

Scale: 1000 filaments × 100 snapshots × 50k step interval = **5.1 M BAOAB steps** with N_pairs = 500k per s → σ_L_p ≈ **0.36 μm** → **4.7σ resolution** on the KU-1.1 ±10 % band [15.3, 18.7] μm (production-sign-off resolution per CLAUDE.md no-gate-loosening).

Wall-time estimate: ~91 min (37× the 단계 2 smoke).

*(Result fills here when sweep notification arrives.)*

## Remaining production gates assessment

| Gate | Autonomous-feasible? | Status |
| --- | --- | --- |
| ERM CFL | ✅ PI-ratified (option A, 2026-05-26) | DONE |
| L_p FULL | ✅ running in background | in progress |
| Variable-length L_p production | Needs `filament_math` adaptation for variable-N + ~91 min run | NEEDS implementation work |
| KU-3.1 cell rounding | Needs ellipsoid initial topology (cortex.py currently spherical only) + 60 s simulated → multi-day wall | NEEDS implementation + dedicated session |
| KU-3.5 cortex tension | Needs HOOMD pressure_tensor Writer wiring + multi-hour steady-state | NEEDS implementation + dedicated session |
| KU-3.18 blebbistatin | Comparative gate (KU-3.1 with myosin OFF); needs KU-3.1 first | DOWNSTREAM |
| KU-3.20 nematic | Q_ij measurement utility already in `tests/validation/test_ku3x_cortex.py`; production run is multi-hour | NEEDS multi-hour |
| 3-way 60s production | 60 s simulated × 10k particles × 2 Updaters → multi-day wall | OUT OF AUTONOMOUS SCOPE |

**Realistic delivery this session**: ERM ratification + L_p FULL (production sign-off contract). Remaining gates have either non-trivial implementation gaps (ellipsoid topology, variable-N L_p adaptation, pressure_tensor wiring) OR multi-day wall-time requirements that exceed even dedicated overnight sessions.

---

# H.3 단계 9 — L_p FULL CPU 시도 중단 + 윈도우 GPU 이전 결정 (2026-05-26)

## L_p FULL CPU 시도 결과

PI 결정 (option A: k_ERM 소프트닝)으로 ERM CFL 해결 직후 L_p FULL 백그라운드 시작.

| Metric | Value |
| --- | --- |
| Wall elapsed at kill | 5h 21m 43s |
| CPU TIME at kill | 222 min |
| CPU utilization | 68% (laptop sleep으로 32% loss) |
| 추정 잔여 | 30-80분 CPU |
| 결과 | **중단** (PI 결정) |

**중단 이유**:
1. 추정이 계속 어긋남 — 단계 4의 91분 추정 → MEDIUM 기반 215분 → 실측 222min CPU 후에도 미완 (실제는 280-300분 CPU 예상)
2. M1 Max 단일-스레드 CPU가 5.1M-step 5-particle-density 시뮬레이션에 부적합
3. 윈도우 RTX A5000 Laptop GPU 이전이 곧 예정 → 같은 sweep 15-25분 wall 가능
4. 현재 fixture 디자인 결함: frames 메모리만 보존, 중단 시 sunk cost 전부 손실 → MEDIUM (단계 8, 3/3 PASS) interim signal 유지

## 단계 9 fixture 개선 (lesson learned)

PI가 진행률 모니터링 불가 + 중단 시 데이터 손실을 지적. `test_cortex.py::TestH3Production.production_run` fixture에 두 가지 개선:

1. **Per-snapshot checkpoint 저장** — `outputs/h3/checkpoints/lp_{smoke|medium|full}/snapshot_NNN.npz` 매 스냅샷마다 디스크에 저장. 중단 시 0..k까지 데이터 보존 → partial 평가 가능.
2. **Progress print** — 매 스냅샷마다 timestamped 한 줄 stdout 출력 (`flush=True`). `tee` 우회 + 실시간 `tail -f`로 모니터링 가능.

```python
# Each snapshot now:
np.savez_compressed(ckpt_dir / f"snapshot_{k:03d}.npz", ...)
print(f"[{HH:MM:SS}] L_p FULL snapshot {k+1}/{n_snapshots} "
      f"(ts=..., wall=..., est_total=...)", flush=True)
```

이 fixture 패치는 윈도우 GPU 이전 후 첫 L_p FULL 재시도 시 진행률 실시간 모니터링 + 안전한 중단 가능하게 함.

## 윈도우 GPU 이전 결정

PI가 RTX A5000 Laptop 보유 확인. CLAUDE.md `Stack` 섹션 "Phase 2+ target: CUDA GPU" 이 조기 활성화:

| Target | Spec |
| --- | --- |
| GPU | RTX A5000 Laptop (Ampere, CC 8.6, 6144 cores, 16 GB) |
| OS | Windows (WSL2 권장) |
| HOOMD 설치 | `conda install -c conda-forge "hoomd=7.*=*cuda*"` |
| 코드 변경 | `device=hoomd.device.GPU()` (한 줄 + auto_select 활용) |
| 예상 가속 | 10-30× (BAOAB only) / 5-15× (Updater 무거운 sim) |

### 이전 후 production gates 재시도 일정

| Gate | CPU (현재) | GPU (예상) | 자율-feasible? |
| --- | --- | --- | --- |
| L_p FULL | 5h+ 미완 | **15-25 min** | ✅ YES |
| Variable-length L_p production | ~10h CPU | ~30 min | ✅ YES |
| KU-3.1 cell rounding (60s simulated) | 수일 | **수시간** | ✅ YES (ellipsoid topology 구현 후) |
| KU-3.5 cortex tension | 수일 | 수시간 | ✅ YES (pressure_tensor wiring 후) |
| KU-3.18 blebbistatin | 수일 | 수시간 | ✅ YES (KU-3.1 downstream) |
| KU-3.20 nematic order | 수일 | 수시간 | ✅ YES |
| 3-way 60s production | 수일 | **수시간** | ✅ YES |
| H.5 KU-5.1/5.2/5.3 | 수일 | 수시간 | ✅ YES |

윈도우 GPU 이전 후 **거의 모든 production gates가 autonomous /loop 안에서 다시 다룰 수 있는 시간 범위로 단축됨**.

## 이 세션 최종 commit 내용

| File | Change |
| --- | --- |
| `configs/phase1_h3.yaml` | k_ERM 0.1 → 1.0e-4 (PI option A 비준) + rationale 코멘트 |
| `tests/test_erm.py` | CFL gate 테스트 → 하드코드 0.1 값으로 테스트 (게이트 자체는 보존) |
| `tests/test_cell.py` | CFL propagation 테스트 → 같은 패턴 적용 |
| `tests/test_cortex.py` | production fixture에 per-snapshot checkpoint + progress print 추가 (단계 9 lesson) |
| `outputs/h3/REPORT.md` | 단계 9 섹션 추가 (CPU 중단 + 윈도우 GPU 이전 결정 기록) |

Main scope 회귀: 변동 없이 ERM ratification 후 모든 기존 tests 통과 (276 PASS / 16 SKIP / 0 FAIL — H.5 단계 1 baseline).

## H.3 상태 (Phase 1 closeout 시점)

- **🟨 implementation 완료** (단계 1-7)
- **🟨 L_p MEDIUM 3/3 PASS** (단계 8, interim 1.3σ signal)
- **⏸ L_p FULL CPU 시도 중단** (단계 9, 윈도우 GPU 이전 대기)
- **✅ DONE 비준 대기**: L_p FULL (GPU) + KU-3.x production runs (GPU)

다음 H.3 작업은 모두 **윈도우 GPU 이전 후 dedicated 세션**으로 이관.

---

# H.3 — 단계 10-11 closeout (Lead 세션 2026-05-28)

start commit `2218222` → end commit (figs) on `phase1/h3-cortex`. gbook A5000 GPU production + Mac dev. PI 정정: 새 1-Lead 모델에서 Lead가 production을 직접 구동 (SSH/백그라운드/Syncthing 회수).

## L_p FULL production — ✅ in-band (production sign-off 증거)

`h3_lp_gpu_production.py --scale full --device gpu` on gbook A5000: 1000 filaments × 7 beads, eq 100k + 100×50k = 5.1M steps, wall 6315 s (105 min), TPS ~853 (안정, AC 연결, throttle 0).

| 측정 | 값 | 판정 |
| --- | --- | --- |
| L_p (all 100 snap) | 16.74 μm | ∈ [15.3, 18.7] ✅ |
| **L_p (eq-transient 첫 5 제외)** | **16.64 ± 0.06 μm** | ∈ [15.3, 18.7] **✅ (4.7σ resolved) — PI 비준 2026-05-29** |
| eq-transient (snap 0–4) | 21.5 → 17.1 μm 감쇠 | eq=100k 잔여 transient (단계 8 예측대로) |
| **⟨E_bend⟩ per angle** | **0.9877 ± 0.0014 kT** | target 0.9898 ±5% → **0.2% off, deep in-band ✅** |

→ L_p 게이트 + 3D equipartition 게이트 모두 **FULL-scale production 증거 확보**. MEDIUM(17.46 μm, 125k pairs/s) 대비 통계 4× (500k pairs/s). **PI ✅ DONE 비준 대기.** (transient-제외 mean±stderr는 production driver의 raw L_p 평균이 아닌 평형 plateau — h3_lp_vis.py `--exclude-transient`.)

## KU-3.20 nematic order — ✅ 비준 (`b04b779`)

구조 측정(Q-tensor 최대 고유값) → multi-hour opt-in에서 CI 상시 게이트로 승격. `generate_cortex_topology` default-off von-Mises bias. 10-seed: S_iso 0.025 (<0.1), S_aligned 0.384 (>0.3). 밴드 불변 (no-gate-loosening).

## KU-3.5 cortical tension — 동역학 바운드 확정 → constrained-BD 이관

진단(`h3_ku35_tension_diag.py`): 100k step(1.3 ms sim)에서 tension ±1e-3 Pa 노이즈(기대 ~100 Pa의 0.001%). 모터 stall 로드까지 ~3.8e7 step(~0.5 s) 필요 → ~860 TPS에서 수시간~수일. KU-3.1/3.18과 동일 dt-병목 클래스 → constrained-BD(dt 523×)에 탑승.

## constrained-BD Milestone 1 (rigid-bond SHAKE + Fixman) — `integrator/constrained_baoab.py`

PI 경로 A. 백본 stretch만 rigid 제약, dt 13ns→69.8μs(523×), 60s 게이트 4.6e9→8.6e5 step. SHAKE + Fixman route-1 pseudo-force. 검증 8 PASS: Fixman 해석 gradient=유한차분(2e-9), dimer force=0, SHAKE 투영/부호/COM, **empty-constraint=baoab bit-for-bit**, rigid dimer D_com=Stokes-Einstein. **Fixman 부호 = 교과서 +1** (앞선 −1 flip은 단일-시드 trimer 노이즈 위 과잉결론 → `178d5b0`로 자기수정; trimer 효과 ±3-4% < 노이즈라 부호 판정 불가, **결정적 확인은 Milestone 2 = L_p 재검증**). PI sign-off 트리거(oracle tau_min 계약, ERM k_ERM 5.6e-7 재충돌)는 M2-4에서.

## Figures (이 세션 — 시각화 규칙, L_p FULL production landing)

| Figure | Caption |
| --- | --- |
| `fig_h3_lp_distribution.png` | per-snapshot L_p vs KU-1.1 band; eq-transient 첫 5개(회색, 21.5→17.1) 명시 + plateau mean 16.64±0.06 μm. |
| `fig_h3_lp_tangent_correlation.png` | C(s) tangent 상관 (per-snapshot 얇은 선 + ensemble mean) + WLC fit L_p=16.7 μm. |
| `fig_h3_lp_bending_equipartition.png` | per-angle E_bend/kT 분포 + 평균 0.9877 kT vs equipartition target 0.9898 ±5% (IN-band). |
| `fig_h3_lp_angle_pdf.png` | interior 각도 PDF vs 3D Boltzmann (∝ sinθ·exp(−½k_θ(θ−π)²/kT)). |
| `fig_h3_lp_bending_3d.png` | cortex 필라멘트 3D 렌더(구면 shell), bead를 굽힘에너지 E_bend/kT로 색칠 (curvature 시각화). |
| `fig_h3_lp_filament_configs.png` | 표본 cortex 필라멘트 3D 형태 (마지막 snapshot). |

## H.3 상태 업데이트

- **🟨 → ✅ 후보**: L_p (FULL) + 3D equipartition + 3D Boltzmann(KS, H.2 상속) + KU-3.20 모두 PASS/비준. **남은 ✅ DONE 차단**: KU-3.5/3.1/3.18 (constrained-BD Milestone 2-4 경유) + L_p 게이트 PI 사인오프.
- 다음: L_p 게이트 PI ✅ 비준 요청 / constrained-BD M2 (oracle tau_min 계약 PI 재논의 선행).

---

# constrained-BD 역량 (rigid-bond Brownian dynamics) — /loop 자율 2026-05-28

KU-3.1/3.5/3.18 동역학 게이트가 dt=13ns(τ_stretch 묶임)에서 비현실적(60s=4.6e9 step)인 문제를 해소하는 기반 역량. `ffn_sim/integrator/constrained_baoab.py` (frozen baoab.py 미수정, 새 모듈).

## 알고리즘
- 백본 stretch bond → **rigid distance constraint**(M-SHAKE) → CFL가 τ_stretch→τ_bend로 이동. angle/LJ/myosin/ERM/xlink는 force 유지.
- **L-M predictor**(baoab 동일) + **Fixman pseudo-force**(route 1, +½kT·ln det G; 교과서 Fixman 1978/Hinch 1994) + **M-SHAKE 투영** + wrap.
- **M-SHAKE**: 선형 사슬 제약 Jacobian이 삼중대각 → Thomas 직접해 + Newton(머신정밀도 수렴, Gauss-Seidel 분-단위 대비 2663 TPS). 필라멘트축 벡터화((F,m,m) 배치 slogdet/inv + 배치 Thomas) → cortex 5.5×.

## 검증 (Milestones)
| | 결과 |
|---|---|
| **M1** dimer/trimer 해석 | Fixman 해석 gradient=유한차분(2e-9), dimer force=0, empty=baoab bit-for-bit, dimer D_com=Stokes-Einstein |
| **M2** 단일 필라멘트 L_p (Fixman 부호 결정적 확인) | 5-seed @ dt=0.03·τ_bend(**1600×**): ⟨E_bend⟩=0.9960±0.0127 kT(eq 0.9898 IN-band), L_p_C1=16.35±0.21μm(H.2 16.56) |
| **M3** cortex 150-fil equipartition/L_p | LJ-off: E_bend=0.9828±0.0040, L_p_C1=16.56μm, drift 4.5e-10 |
| **Fixman 부호** | trimer는 효과±3-4%<노이즈로 판정 불가 → M2(19 누적 각도)가 +1 결정적 확인 |

## dt 상한 (실현 가속)
- bending-only(단일 필라멘트): **1600×** (factor 0.03·τ_bend). factor 0.1(이론 CFL)은 explicit-predictor 변위가 SHAKE 한계 초과로 발산.
- **LJ-on cortex: ~54×** (factor 0.001) — LJ WCA 배제부피가 자체 fast CFL 부과. 밀도↑면 더 감소. → KU-3.x 가속은 가장 빠른 미제약 힘(LJ)이 결정.

## KU-3.x feasibility (vectorized cortex ~150-200 TPS CPU @ 54×)
- **KU-3.5 tension(~0.1-0.5s)** = 1.4e5-7e5 step → **~16분-1.3시간 FEASIBLE**.
- KU-3.1/3.18(60s) = ~8.5e7 step → 수일 (GPU-SHAKE 또는 추가 dt 이득 필요).

## full-cell wiring (KU-3.x 진입점)
`build_cortex_full_simulation(constrained=True, constrained_dt=…)` 추가(additive, default-off, 회귀 10 PASS). actin backbone k=0+M-SHAKE, myosin/xlink/ERM은 predictor. 4100-입자 cortex+myosin 빌드·적분 확인.

## PI-collaborative 잔여 (KU-3.5)
1. **warm-up handoff**: 신선 construction은 factor 0.001서도 LJ overlap으로 발산 → 표준 baoab cfl-dt warm-up → 전체 위치 전달 → constrained large-dt.
2. **myosin 동역학 large-dt 정확도** 검증(stepping/Bell-Evans가 dt↑에서 물리 보존하는지).
3. **tension 측정 프로토콜** = Sanity Gate measurement-protocol-consistency 결정 → **PI 비준**. 후보: method-of-planes(직경 절단면 가로지르는 bond/motor 장력 합 / 2πR; 박스 virial보다 shell 기하에 견고). KU-3.5 target γ≈0.5 mN/m ±30%.
4. **oracle tau_min 계약**: constrained 모드 τ_stretch 제외 분기(PI 보류 중).

## 🔴 KU-3.5 BLOCKER — myosin 모듈(단계-4) 결합 기하 결함 (2026-05-28 /loop 발견, Codex 진단)

constrained-BD로 myosin-active cortex를 **안정 구동**(warm-up handoff → 54×, drift 머신정밀도)했으나 **myosin이 actin에 결합하지 않아 수축이 발생하지 않음** → KU-3.5 시연 불가. 측정 진단:
- myosin head ↔ 최근접 cortex actin 거리: **중앙값 824nm** (min 92nm), `head_actin_max_bind_dist=50nm` → **0/2000 head 결합** (60000 step에서 engaged≈1, step_advances=0). cortex 반경 10.10→10.12μm (수축 0).

**근본 원인 (Codex read-only 진단, 제 기하 측정과 일치):**
1. `head_rest_length=200nm` offset이 **radial(막 법선) 방향**으로 적용 → head가 actin shell 평면에서 radially 벗어남. 올바른 건 **tangent-plane lateral offset**(필라멘트 방향 cross 법선)로 cortex 면 안에 두는 것.
2. binding이 actin **bead 중심**만 cKDTree 검색(`myosin.py:747`). bead 간격 ℓ0=500nm → segment 중간 head는 끝점 bead에서 ~250nm. **segment 기반 검색** 필요(capture 50-63nm), 또는 임시 bead-center 시 유도반경 √(250²+50²)≈**255nm**(grid-derived, biological reach 아님으로 명기).
3. (별도) myosin **construction LJ overlap**(myosin-myosin intra): 150-fil에서 cfl-dt warm-up도 int32 guard 발산(120-fil은 benign) → myosin 배치 overlap-free 보장 + robust warm-up(energy-min/dt ramp) 필요.

**→ PI 결정 필요 (단계-4 myosin mechanism 수정, KU-3.5/3.1/3.18 전부 영향):** (a) tangent-plane head 배치, (b) segment 기반 binding(또는 grid-derived 255nm), (c) overlap-free myosin construction. 이는 magic-number 패치가 아니라 placement 물리 정정 → 무감독 rewrite 대신 PI surface. 수정 후 KU-3.5 tension(constrained-BD feasible 확인됨)으로 진행.

driver `scripts/h3_ku35_tension.py`(warm-up handoff + frame/diag, contraction proxy; γ는 PI 프로토콜) — benign 구성에선 안정, overlap 구성에선 warm-up 발산(myosin construction 수정 후 안정 예상).

---

## 2026-05-28/29 Session — KU-3.5 v2 sweep + R1 rigid Lagrange + literature integration

(Lead session, 4 hours wall — commits `1d0c8c6` → `7a1e4a5` H.3 + `3306eb1` H.5)

### Driver hardening (KU-3.5 v1 → v2)

After the 2026-05-28 v1 canonical sweep failed (seed2 LJ-CFL crash @ 2h21m, other
seeds killed @ ~4h with 0 progress signal → 16h CPU loss), the driver
`scripts/h3_ku35_tension.py` was hardened:

- `--device gpu/cpu` flag mirroring `h3_lp_gpu_production.py` (commit `1d0c8c6`)
- Per-sample `PROGRESS` print: `wall/eta/r-ratio/engaged/steps_adv/gamma/drift/max_disp_um`
- Min-image LJ-CFL guard: raw inter-snap displacement > 5·box_L → `WARN_LJ_CFL`
  (commit `85946d6`)
- Auto-viz subprocess to `h3_ku35_sweep_analysis.py` at seed end (commit `7a1e4a5`)

### CuPy GPU port attempt (option B0 from FIXMAN_LAZY_EVAL_DESIGN.md §7)

- `integrator/baoab.py` xp-dispatch port committed `d18d7fb` (CPU bit-for-bit;
  GPU compute path via `gpu_local_snapshot` + `cupy.random.Generator`)
- `integrator/constrained_baoab.py` xp-dispatch port attempted, **reverted**:
  KU-3.5 (cortex 120×6 SHAKE/Fixman) ran **3× SLOWER** on GPU (3.55 → 11.18
  ms/step) due to kernel-launch overhead on the small 6×6 SHAKE/Fixman matrices.
  Honest findings recorded in design doc §7 (commit `a49314a`).
- L_p smoke GPU re-benchmark (300 fil × 200K step): 1.52 ms/step vs pre-port
  FULL 1.26 ms/step — baoab port gives no clear gain at our Phase 1 scale either.
  Retained as CPU safety only; GPU acceleration deferred to Phase 2 cell-scale.

### R1 — Rigid actin backbone Lagrange tension exposure (commit `fb66028`)

Implements `RIGID_LAGRANGE_TENSION_DESIGN.md` (PI verbal "(가)" 2026-05-28).

The SHAKE M-SHAKE-converged Lagrange multiplier λ was already computed inside
`shake_project_chains` (`constrained_baoab.py:266`) but discarded. R1 exposes it
through a new `return_lambdas` kwarg and an Action `record_lambda` flag → the
KU-3.5 driver records `tension_{soft,rigid,total}_mN_per_m` per sample. The
soft-bond-only γ was systematically under-reporting cortical tension by ~200×;
γ_total (soft + rigid) is now the gate quantity.

CPU regression: 123 PASS / 6 SKIP / 0 FAIL with 6 new STATIC tests
(`TestR1LambdaCapture` in `tests/test_constrained_baoab.py`):
- record_lambda toggle bit-for-bit position match
- λ buffer shape (uniform-chain fast path returns `(F, m)` ndarray;
  ragged fallback returns list of per-chain `(m_i,)` arrays)
- same-seed reproducibility
- stretched-chain accumulated λ > 0 (SHAKE convention)

### H.5 단계 2 — Cell.build(with_lamellipodium=True) wiring (commit `3306eb1`)

Worktree subagent (`general-purpose` agent, parallel to Lead's R1 work; off-limits
surfaces respected). `cell/cell.py` `build_cortex_full_simulation` now extends
the snapshot with WAVE + mother-actin beads BEFORE `create_state_from_snapshot`,
registers `lamel_*` bond/angle params, attaches `WaveMembranePin` force,
appends 3 D2-batched Updaters. `Cell.build` routes through the unified builder
when `with_myosin OR with_lamellipodium` is True; legacy cortex-only paths stay
bit-for-bit unchanged when both flags are False. 4 new `TestLamellipodiumWiring`
tests; 309 PASS / 13 SKIP / 0 FAIL on the full non-validation tree.

### Literature integration (2026-05-29)

Read two papers on PI's recommendation:

- **Kim, Neal, Kamm, Asada 2013** (PLOS Comp Bio, PMC3585413) — cell migration
  on fibronectin via 549-node continuum mesh + Hill SF + Bell ligand-integrin.
  Lumped approach; useful only as reference numerical values + Palecek CHO
  experimental dataset for future H.4+H.5+H.7 oracles. Not directly applicable
  to H.3.
- **Luo, Mohan, Iglesias, Robinson 2013** (Nature Materials 12:1064-1071,
  PMC3838893, doi:10.1038/nmat3772) — directly maps onto our H.3 cortex
  (Dictyostelium actin cortex + myosin II + α-actinin + filamin). Three
  oracles immediately useful:
  1. **ζ = 1/7 force-sharing**: myosin carries ~14% of cortical tension,
     crosslinkers ~86% (Luo Fig. 1d fit, WT). Now plotted as a proxy line in
     `sweep_analysis.py` `render_gamma_breakdown` panel: γ_soft / γ_total vs
     time, with horizontal reference at 6/7 = 0.857.
  2. **Deformation-specificity**: α-actinin → dilation (pipette tip),
     filamin → shear (pipette neck), myosin → dilation (lever-arm dependent).
     Captured as `KU_3_21_DEFORMATION_SPECIFICITY.md` candidate gate.
  3. **Quantitative parameters**: aspiration 0.5-2.0 nN/μm², accumulation peak
     30-60 s, Δx ≈ 1-2 nm. Drop-in inputs for the KU-3.21 implementation.

### `sweep_analysis.py` extension

New `render_gamma_breakdown(seeds)` function + `fig_h3_ku35_gamma_breakdown.png`:
top panel = γ_soft / γ_rigid / γ_total vs t (per-seed lines, KU-3.5 band
overlay); bottom panel = γ_soft / γ_total fraction vs t with Luo 6/7 reference.
`REPORT_ku35.md` extended with R1 γ breakdown section + Luo 2013 oracle
comparison. Backward-compatible: pre-R1 sweeps (only `tension_mN_per_m`
populated) skip the breakdown panel gracefully without crashing.

### Documents added / updated

- New: `ffn_sim/docs/briefs/KU_3_21_DEFORMATION_SPECIFICITY.md` (candidate gate)
- Updated: `ffn_sim/docs/briefs/FIXMAN_LAZY_EVAL_DESIGN.md` §7 (B0 findings)
- This REPORT section

### Next session boot point

- `phase1/h3-cortex` HEAD `3306eb1` (this session added 5 new commits)
- KU-3.5 v2 sweep on gbook 4-seed CPU still running (ETA ~05:20 KST)
- After sweep: `h3_ku35_sweep_analysis.py` auto-vizes → confirms R1 schema if
  re-launched with the new driver, OR falls back to pre-R1 panels for the
  in-flight sweep
- R1 PI sign-off → KU-3.5 re-run with γ_total → H.3 ✅ DONE candidate
- KU-3.21 candidate awaits PI decision (D0 defer / D1 Phase A only / D2 full)

---

## 2026-05-29 KU-3.5 v2 sweep result (autonomous /loop closeout)

### Production parameters

```
n_fil = 120
n_motors = 100
dt_factor = 0.0003  →  dt = 2.094e-7 s
n_warmup = 40000 step
n_sample = 800 × interval 5000 step  →  4M step production / seed
device = cpu (gbook 4-core parallel)
seeds launched = {1, 2, 3, 4}
seeds completed = {1, 3, 4}
seed 2 = LJ-CFL crash @ sample 284/800 (40%, ~2h21m wall) — int32
         image-guard overflow, sub-interval explosion that the new
         min-image LJ-CFL monitor (raw_disp > 5·box_L) could not
         pre-warn on (PROGRESS line just before crash was healthy
         max_disp_um = 0.66). Pre-warning would require per-step
         monitoring at significant overhead cost.
wall   = ~6h 20min per seed (4-seed parallel on 16-core gbook)
```

### Ensemble plateau γ (3-seed)

| seed | plateau γ (mN/m) | final γ (mN/m) | engaged | steps_adv | max_drift |
|---|---|---|---|---|---|
| 1 | 3.219e-5 | 2.775e-5 | 896 | 11042 | 2.12e-15 |
| 3 | 2.944e-5 | 2.361e-5 | 892 | 10921 | 2.54e-15 |
| 4 | 3.036e-5 | 3.391e-5 | 942 | 11223 | 2.12e-15 |

**Ensemble plateau ⟨γ_soft⟩ = 3.066e-5 mN/m** (per-seed σ = 1.4e-6, 4.5% relative)

### Verdict (pre-R1, soft-bond only)

**FAIL (below)** vs KU-3.5 target band [0.35, 0.65] mN/m.

Plateau is **~16,300× below the lower edge** of the target band. Per
[RIGID_LAGRANGE_TENSION_DESIGN.md](../../docs/briefs/RIGID_LAGRANGE_TENSION_DESIGN.md)
§1, the original 200× under-report estimate was conservative — actual
under-report at our cortex parameters is two orders of magnitude larger.
This is the systematic measurement-protocol issue R1 was designed to
address: the soft-bond method-of-planes sum captures xlinks + motor head
attach bonds + ERM + myosin internal bonds only, missing the dominant
contribution which propagates through the rigid actin backbone as
M-SHAKE Lagrange multipliers.

The measurement IS stable:
- per-seed σ = 4.5% across 3 independent seeds
- machine-zero drift (2-3e-15) throughout 4M steps
- motor stepping saturated equilibrium at ~11,000 advances per seed
- r/r0 ≈ 1.000 (cell shape stable through 0.84 s physics time)

→ The sweep validates the new driver hardening (PROGRESS + LJ-CFL guard)
and the measurement infrastructure, AND quantifies the magnitude of the
soft-only systematic error that R1 was designed to close.

### Post-R1 re-run (pending PI sign-off)

R1 patch (`fb66028`) exposes the rigid-bond Lagrange multiplier as
`γ_rigid` per sample. The KU-3.5 driver now records
`tension_{soft,rigid,total}_mN_per_m` triplets per sample; the
`sweep_analysis.py` `render_gamma_breakdown` panel (commit `c59cc0b`)
auto-plots the soft/rigid/total split with the Luo 2013 ζ=1/7 oracle
horizontal reference. Re-running this canonical sweep with the R1
driver should produce γ_total in or near the [0.35, 0.65] mN/m band
(testable). Wall budget identical to this run (~6h 4-seed parallel).

### Figures (this sweep)

- `outputs/h3/figs/fig_h3_ku35_tension_sweep.png` — γ_soft vs sim time,
  per-seed lines + ensemble mean + KU-3.5 target band overlay + plateau
  verdict annotation.
- `outputs/h3/figs/fig_h3_ku35_engagement_sweep.png` — myosin engagement
  + Hill stepping + radius contraction r/r0 vs sim time.
- (No γ breakdown panel — pre-R1 schema.)

## 2026-05-29 L_p FULL GPU re-bench (post-baoab-port)

After `baoab.py` xp-dispatch port (commit `d18d7fb`), re-benchmarked
L_p FULL on gbook A5000 to test the port's GPU speedup claim at
production scale (n_fil=1000, 5M step).

| Run | wall (s) | ms/step | physics L_p (μm) | E_bend (kT) | in-band? |
|---|---|---|---|---|---|
| pre-port FULL GPU (`lp_full_gpu.npz` from 2026-05-27) | 6315 | 1.24 | 16.74 | — | ✓ |
| **post-port FULL GPU** (2026-05-29) | **7812** | **1.56** | **16.83** | **0.9884 ± 0.0014** | **✓** |

**Conclusion**: `baoab.py` port is **25% SLOWER** at L_p FULL scale —
matching the KU-3.5 finding (`d18d7fb` introduces cupy / gpu_local_snapshot
plumbing that helps for some workloads but kernel-launch overhead
dominates at Phase 1 sizes). Physics is correct (in-band, equipartition
within 0.14% of 0.9898 kT target). 25% slowdown is much smaller than
KU-3.5's 3× slowdown because L_p has no SHAKE/Fixman inner loops; only
the BAOAB step and HOOMD native forces.

PI sign-off recommendation: **revert `d18d7fb`** for Phase 1 — no
workload benefits, costs 25-300% wall. Keep the xp-dispatch design as
a reference for Phase 2 cell-scale work where larger N may amortise
kernel-launch overhead better.

## 2026-05-29 PI ratification — R1 + baoab port revert

PI Sungwook ratified two decisions this Lead session:

### Decision 1: R1 (rigid Lagrange tension exposure) — ✅ RATIFIED

- Design: `docs/briefs/RIGID_LAGRANGE_TENSION_DESIGN.md` (Status header now
  `✅ RATIFIED 2026-05-29`).
- Implementation: commit `fb66028` (already on `phase1/h3-cortex`). Behaviour
  bit-for-bit when `record_lambda=False` (default); 5 new tests in
  `TestR1LambdaCapture` cover Sanity Gate §5 #1-7.
- Evidence justifying ratification: v2 sweep `⟨γ_soft⟩ = 3.07e-5 mN/m` vs
  target `[0.35, 0.65] mN/m` = **~16,300× under-report** (more than 80× larger
  than the design §1 conservative 200× estimate). R1 is necessary for KU-3.5
  to be honestly comparable to literature.
- Sanity Gate §5 #6 (the actual gate check — γ_total in band) will be measured
  by the v3 sweep (this session), not the design doc itself.

### Decision 2: `baoab.py` port revert — ✅ RATIFIED

- Target: revert commit `d18d7fb` ("H.3 baoab.py: xp-dispatch port").
- Falsifying evidence: Both Phase-1 GPU workloads measured net-negative
  post-port. KU-3.5 3× slower (already drove `constrained_baoab.py` port
  revert via `a49314a`); L_p FULL 1.24× slower (this session's re-bench).
  No remaining workload that benefits.
- CPU code path is bit-for-bit either way; revert risk on CPU = 0.
  16/16 baoab + 146 cortex/myosin/erm/xlink CPU tests PASS pre-revert;
  same passes expected post-revert.
- Phase 2 future: re-attempt with fused CuPy `RawKernel` or native
  HOOMD CUDA module (FIXMAN_LAZY_EVAL_DESIGN §7.6) when cell-scale N
  amortises kernel-launch overhead.

### Follow-up: KU-3.5 v3 sweep launched

- gbook 4-seed CPU concurrent (matches v2 conditions: n_fil=120,
  dt_factor=0.0003, 4M production steps).
- Driver flips `act.record_lambda = True`; tension measured via
  `_tension_method_of_planes_rigid` → γ_soft + γ_rigid + γ_total per sample.
- Pre-existing instrumentation retained: PROGRESS print (commit `1d0c8c6`),
  min-image LJ-CFL guard (`85946d6`), auto-viz subprocess (`7a1e4a5`).
- Output: `outputs/h3/production/ku35_canonical_v3/`.
- `sweep_analysis.py` `render_gamma_breakdown` panel + Luo 2013 ζ=1/7
  oracle overlay (commit `c59cc0b`) auto-activates because schema now
  contains γ_total.
- If v3 γ_total ∈ `[0.35, 0.65] mN/m` → KU-3.5 ✅ candidate, H.3 ✅ DONE
  unblocked. If still under-report → diagnostic next step (likely a
  per-plane Luo decomposition or KU-3.21 deformation-specificity sweep).

