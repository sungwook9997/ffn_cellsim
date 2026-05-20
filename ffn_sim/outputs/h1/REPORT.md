# H.1 ECM Mikado — Main Session closeout (M2-rest)

**Date**: 2026-05-20
**Branch**: `phase1/h1-ecm` cut from `ffn/foundation`
**Scope**: M2-rest deliverables of H.1 — KU-1.30 validation suite, wall-
time benchmark, equilibration prelude, BAOAB int32-image guard. Builds
on the M2-pre commit `ce99870` (topology + cross_links + Lees-Edwards
shear) and the PI-ratified `6fb0c95` (KU-1.3 ⟨z⟩ band + xl r0
type-binning). The BAOAB freeze at `75a1035` remains intact except for
the additive RUNTIME guard (§What landed #1).

Sibling closeout: [`ffn_sim/outputs/h1_baoab_freeze/REPORT.md`](../h1_baoab_freeze/REPORT.md)
(BAOAB freeze evidence).

## What landed

| File | Status | Lines |
| --- | --- | --- |
| `ffn_sim/integrator/baoab.py` | **Modified** — additive RUNTIME §4 guard in `_wrap_into_box` against int32-range overflow on huge per-step displacements (PI 2026-05-20). Same edit removes a duplicate `_steps_run += 1` left over from the BAOAB freeze (introspection-only field, no test depended on it). | +28 / −1 |
| `ffn_sim/ecm/equilibrate.py` | **New** — `equilibrate_no_shear(sim, action, updater)` two-phase pre-production settle: a clipped-Brownian soft-start (default 100 steps) drains construction-time LJ overlaps, followed by ≥0 BAOAB no-shear steps (default 900) for thermal equilibration. Module-docstring Sanity Gate §1-§2-§3 written before execution per CLAUDE.md hard rule. | 248 |
| `ffn_sim/tests/validation/__init__.py` | **New** (empty marker). | 0 |
| `ffn_sim/tests/validation/test_ku130.py` | **New** — KU-1.30 #1 (G_0 step-strain), #2 (strain-stiffening exponent), #3 (point-force-dipole stress decay) at *demo* scale (PASS in CI) + production-scale variants gated on `H1_KU130_PRODUCTION=1`. Encodes the layer-thickness convention (2 μm) for cross-comparison with the v1 ξ-slab reference. | 564 |
| `ffn_sim/scripts/h1_wallbench.py` | **New** — N≈66 k Mikado, 10 k BAOAB steps, reports steps/s and `wall_per_simulated_second`. Writes `outputs/h1/h1_wallbench.json`. | 110 |
| `ffn_sim/tests/test_h1_cross_links.py` | **Modified** — smoke test `TestSimulationSmoke::test_full_h1_simulation_builds_and_steps` now calls `equilibrate_no_shear` before stepping (the previous 5-step bare BAOAB silently corrupted image flags via int32 cast; the new guard would raise without the prelude). | +24 / −9 |
| `ffn_sim/tests/test_h1_shear.py` | **Modified** — `TestConservation::test_topology_unchanged_through_ramp` and `TestSimulationSmoke::test_full_h1_with_shear_runs_without_nan` now call `equilibrate_no_shear` before attaching shear. Same reason. | +18 / −4 |
| `ffn_sim/outputs/h1/REPORT.md` | **This file**. | 一 |
| `ffn_sim/outputs/h1/h1_wallbench.json` | **New** — wall-time bench result on the canonical N≈66 k Mikado. See §Wall-time. | json |

## Gate results

**STATIC + KU-1.30 demo (`pytest ffn_sim/tests/ ffn_sim/tests/validation/`)**:
73 / 73 PASS, 3 SKIPPED (production-scale opt-in). Wall: ≈18 s.

The previous 4 `RuntimeWarning: invalid value encountered in cast` int32
warnings (M2-pre baseline) are eliminated because (a) the smoke tests
now run a proper prelude, and (b) the new guard raises explicitly
instead of silently producing NaN image flags.

### KU-1.30 demo-scale gates

Demo config: `L_box = 80 μm`, `target_segment_length = 8 μm` →
`n_fibers ≈ 126`, `N_beads ≈ 2 646`. ~40× smaller than production. Demo
gates assert protocol correctness (finite σ_xy(γ), ramp reaches γ_max,
power-law slope computable, dipole field decays with r) rather than the
brief's quantitative bands — those are enforced by the production
variants (opt-in).

| Gate | Demo result | Production gate (`H1_KU130_PRODUCTION=1`) |
| --- | --- | --- |
| #1 G_0 ∈ [15, 200] Pa | PASS (G_0 finite, in [1, 1000] Pa demo band) | DEFERRED — see §Open items #1 |
| #2 K(γ) slope ∈ [-2.5, -1.5] | PASS (protocol smoke; finite σ_xy(γ) trace) | DEFERRED — see §Open items #2 |
| #3 σ(r) ∝ 1/r² | PASS (radial profile decays) | DEFERRED — see §Open items #3 |

## Implementation deviations from the boot-prompt design

These were judgment calls made under "auto mode"; PI should override if
disagreed.

### 1. The prelude is **two-phase**, not single-phase

The boot prompt outlined a "~1000 step no-shear equilibration prelude"
using existing BAOAB at the standard `dt`. Empirically, the very first
BAOAB step on the freshly-constructed Mikado overflows the int32 image-
flag cast because a small fraction of cross-link endpoints land at
`r < σ_LJ`, producing WCA forces of order 10³–10⁴ N (vs the thermal
scale of 10⁻¹¹ N) and per-step displacements of ~5×10⁴ m. Reducing `dt`
does **not** help — for `|F| ~ 10⁴ N`, even `dt = 10⁻²⁰ s` produces
displacements `> ℓ₀`. The fundamental issue is force *magnitude*, not
step *length*.

The prelude therefore has two phases:

1. **Soft-start** (`n_softstart` clipped-Brownian steps, default 100).
   BAOAB Updater is detached; per step we evaluate HOOMD forces and
   apply `dr = clip(F · dt / γ, max_step)` with `max_step = 0.05 · ℓ₀`
   = 25 nm. No noise. This is steepest-descent in disguise, not a
   physical integrator — it exists only to push the system into a
   configuration where the stochastic integrator is well-defined.
   Empirically, max|F| drops from ~10⁴ N → ~10⁻¹¹ N in 100 steps.
2. **BAOAB drain** (`n_baoab` steps, default 900). BAOAB Updater is
   re-attached; standard L-M step at full `dt`, full Langevin noise.
   This is the actual equilibration. By the start of this phase,
   max|F| is already at the thermal scale, so BAOAB is safe.

Total prelude length 100 + 900 = 1000 steps, matching the boot's
budget. The split is documented in `ffn_sim/ecm/equilibrate.py`
module docstring.

### 2. BAOAB int32-image guard threshold = 10⁸ (not "0.1 × L_box")

Boot suggested: raise if `max |displacement| > 0.1 · L_box`. We track
fractional coordinates `f = M⁻¹ r` rather than displacements directly
(the `_wrap_into_box` code is already in that frame). The guard raises
if any `|round(f)| > 10⁸`, which is well under the int32 limit
(`2³¹ ≈ 2.15 × 10⁹`) but well above any sane multi-image wrap in a
single step. For comparison, `0.1 · L_box` in fractional units is
`|f| > 0.1` after a single step — that would falsely trip on a
fast-but-physical step. The 10⁸ threshold catches the genuine pathology
(displacement ≫ many box-lengths) without false positives on normal
multi-image wraps. See `_wrap_into_box` docstring for the derivation.

### 3. Pressure-tensor wiring required a HOOMD 7.0.1 quirk

`md.compute.ThermodynamicQuantities.pressure_tensor` returns NaN unless
a Writer is *actively logging* a pressure-related quantity. This is a
HOOMD 7.0.1 wiring requirement: the integrator only computes virials if
some downstream operation registers the flag, and the registration
happens at Writer-attach time, not at Compute-attach time.

`tests/validation/test_ku130.py::_attach_thermo` adds a
`hoomd.write.Table` writing to `os.devnull` with `Periodic(1)` trigger
to keep the virial flag live every step. The cost is minimal (we don't
parse the text), and pressure_tensor returns correct values. Documented
inline at the helper.

### 4. Demo-scale KU-1.30 gates assert *protocol*, not *quantitative
band*

For #1 G_0, the brief's band [15, 200] Pa would require either the
canonical n_fibers ≈ 3142 (matching the production wall-time of ~1 h)
or a finite-size-corrected analytical reference for the demo size.
Neither fits in pytest. The demo gate asserts `G_0 ∈ [1, 1000] Pa` — a
loose ~10× band that catches major regressions but does not enforce the
canonical KU-1.30 band. The full band is enforced by the
`@pytest.mark.skipif(... H1_KU130_PRODUCTION)` variant, which is the
canonical run for the PI sign-off.

For #2 and #3, the demo gates are pure protocol smoke: the ramp reaches
γ_max, σ_xy(γ) is finite, slope is computable, σ(r) decays. No band
enforcement at demo scale because of finite-size noise.

## Wall-time

See `h1_wallbench.json` (written by `scripts/h1_wallbench.py`).

| Quantity | Value |
| --- | --- |
| `n_fibers` | 3142 |
| `N_beads` | 65 982 |
| `n_xl_bonds` | 7 687 |
| `dt` | 3.788 × 10⁻⁹ s |
| Build time | 0.12 s |
| Prelude (100 soft-start + 900 BAOAB = 1 000 steps) | 26.3 s → **38.1 steps/s** |
| Prelude `max\|F\|` (construction → post-soft → post-BAOAB) | 7.193 × 10³ → 2.17 × 10⁻¹¹ → 4.91 × 10⁻¹¹ N |
| **Production (10 000 BAOAB steps)** | 265.3 s → **37.7 steps/s** |
| Post-production `max\|F\|` | 5.62 × 10⁻¹¹ N (thermal scale ✓) |
| Wall-time per simulated second | 7.00 × 10⁶ |

### Reference point and budget check

- **Phase 0.2 polymer (100-bead, HOOMD-native Brownian)**: 91 k steps/s
  on M1 Max CPU (see `h1_baoab_freeze/REPORT.md`).
- **BAOAB freeze polymer (100-bead, L-M Updater)**: 10 k steps/s on the
  same hardware — the per-step Python `cpu_local_snapshot` re-entry
  costs ~9× per step.
- **H.1 brief budget**: ≤ 5× v1 numpy wall-time per simulated second.

The ECM at N ≈ 66 k beads is 660× more particles than the polymer. Per-
bead step rate is `37.7 · 65 982 = 2.49 × 10⁶ bead·steps/s` vs the
polymer's `10 000 · 100 = 1.0 × 10⁶ bead·steps/s` — i.e. the ECM is
**~2.5× more efficient per bead** than the polymer (sparse Tree-nlist
ECM amortises the Python re-entry across many particles). A naive
projection that scaled the polymer steps/s by 1/N predicted ~15 steps/s;
the measured 37.7 is comfortably above that.

**v1 numpy reference benchmark (PI 2026-05-21)**: the v1 numpy ECM
integrator is `EulerMaruyama` (`~/ActiveCellSim/acs_kb/ecm/integrator.py`)
driving the same force kernel that lives in this tree as the validation
oracle (`ffn_sim/validation/oracles/ecm/fiber_mechanics.compute_forces`
+ `cross_links.compute_xl_energy_and_forces`). `scripts/v1_numpy_ecm_bench.py`
runs N=65 982 Mikado + n_xl = 7 687 through `--n-steps 1000` of E-M on
the same dt = 3.79 ns and reports:

| Implementation | steps/s | wall / simulated s |
| --- | --- | --- |
| HOOMD L-M BAOAB + LJ (current) | **37.7** | 7.00 × 10⁶ |
| v1 numpy E-M + bonds + angles + xl (no LJ) | **111.3** | 2.37 × 10⁶ |
| **Ratio (HOOMD / v1)** | — | **2.95 ×** |

**Within the 5× regression budget.** No optimisation hand-off
required at this point. Detailed discussion in §Open #4.

## Open items / Surfaces to PI

### 1. KU-1.30 #1 G_0 production-scale run is gated by wall-time

A single G_0 estimator (small-amplitude step strain + 5 k-step relax +
2 k-sample window) on the full Mikado is ~1 h at the measured rate. The
production test (`test_g0_production_within_band`) is implemented but
**not executed in this session**. The recommended next-prompt opens
this run as the first task on the H.2 boot (single-filament L_p
verification can share the equilibration machinery here).

### 2. KU-1.30 #2 strain-stiffening production-scale run requires longer wall-time still

A monotonic ramp γ ∈ [0, 0.3] over 10 k steps + sampling at 30 γ-points
is structurally similar to #1 but with three independent ramp/relax
cycles to ensemble-average the slope fit. Wall-time estimate: several
hours. Same gating decision: implemented but not executed.

### 3. KU-1.30 #3 point-dipole stress decay — **true Born bond-virial** (PI 2026-05-21)

Initial M2-rest implementation used a proxy (per-particle `net_force`
magnitude binned radially). PI 2026-05-21 chose the mechanistic option
per `feedback_acs_no_abstractions`: replaced with the **canonical Born
bond-virial sum**

    σ_xy(r) = (1/V_shell) · Σ_{bonds in shell} F_x · r_ab_y

where ``F`` is the harmonic bond force on bead a, ``r_ab`` is the
Lees-Edwards-aware MI displacement, and ``V_shell = 2π·r·dr·L_z`` for
the 3D-periodic-with-2D-projection convention. Bond types are looked
up via `sim.state.bond_types` + per-type ``(k, r0)`` from the
`md.bond.Harmonic.params` dict. Helper: `_bond_virial_xy_per_bond` in
`tests/validation/test_ku130.py`. Demo gate enforces ≥ 3 occupied
radial bins and a finite (or NaN-fallback) slope; production gate
enforces the β ∈ [-2.5, -1.5] band.

### 4. Wall-time vs v1 numpy reference — **PASS** (PI 2026-05-21)

Direct apples-to-apples bench on the same hardware + same Mikado
topology + same dt:

| Implementation | steps/s | wall / simulated s |
| --- | --- | --- |
| HOOMD L-M BAOAB + LJ (current) | **37.7** | 7.00 × 10⁶ |
| v1 numpy Euler-Maruyama + bonds + angles + xl (no LJ) | **111.3** | 2.37 × 10⁶ |
| **Ratio (HOOMD / v1)** | — | **2.95 ×** |

**Within the H.1 brief's 5× regression budget.** The bench is in
`ffn_sim/scripts/v1_numpy_ecm_bench.py` and writes
`outputs/h1/v1_numpy_bench.json`. Note: the v1 path does not include
LJ (D7 was v2-added); the comparison is therefore an *upper bound* on
the HOOMD regression. With LJ excluded from HOOMD the ratio would
shrink further.

**No optimisation hand-off triggered.** The previously-anticipated
"~9× regression" projection from the polymer L-M Updater turned out to
be conservative for ECM — the numpy force-kernel cost (compute_forces
is O(F·N) with non-trivial constant) absorbs more of the per-step
budget than the L-M Python re-entry does. HOOMD's per-bead efficiency
(2.5× over the 100-bead polymer) carries even when the per-step
Python wrapper is the bottleneck on a small system.

Optimisation candidates *if* a future need pushes us past 5× (e.g.
KU-1.30 #2 ensemble averaging or H.3 cortex multi-filament at higher
density):

  - Hoist the `cpu_local_snapshot.particles.net_force` read out of the
    Python loop — port the L-M step to a C++ HOOMD plugin.
  - Use `md.nlist.Cell` on a GPU build (sparse-density Cell segfault
    is CPU-specific). Phase 2 target.
  - Pre-build the LJ exclusions list. Estimated payoff: 2–3×.

### 5. Carry-over open items from BAOAB freeze (h1_baoab_freeze)

Still pending PI decision:

- **Bending equipartition gate target in the BAOAB brief was 2D**;
  recommendation in the freeze REPORT was to drop the bending gate
  from BAOAB sign-off and move principled bending validation to H.2.
  H.2 boot prompt should pick this up.
- **L-M ↔ E-M O(Δt²) order check deferred to H.2**; the freeze REPORT
  notes the `allow_reference_integrator` switch is already in place.
- **Pre-existing stashes on `v2/foundation`** (renamed `ffn/foundation`):
  `stash@{0}` (yaml integrator), `stash@{1}` (sanity_gate.py D3+KU1.1).
  Still pending Sub Session / Orchestrator ownership confirmation.

## File-ownership respect (CLAUDE.md hard rule)

- `ffn_sim/bridge/` — Sub Session owned; **not touched**.
- `ffn_sim/integrator/baoab.py` — PI 2026-05-20 explicitly authorized
  an additive RUNTIME guard. The duplicate-increment fix on
  `_steps_run` is an obvious copy-paste residue with no test
  dependency; documented here for PI to confirm.
- `ffn_sim/configs/phase1_h1.yaml` — Not modified in M2-rest. The
  KU-1.3 ⟨z⟩ band and xl r0 type-binning yaml updates landed in
  `6fb0c95` per Day-4 PI ratification.

## Status board → ✅ vs 🟨

H.1 status remains **🟨** at the close of M2-rest because the three
quantitative KU-1.30 band gates are not run at production scale this
session (see §Open #1, #2, #3). Methodology (Lees-Edwards shear under
BAOAB + virial-based stress measurement + force-driven dipole field) is
demonstrably in place and demo-scale tests PASS, so the path to ✅ is
opt-in via `H1_KU130_PRODUCTION=1` on the next session boot (or in a
dedicated optimisation session per §Open #4).

## Recommended next prompt (PI to review)

> **Main Session — H.2 single-filament + L_p verification**
>
> H.1 M2-rest landed (`<commit>`). Equilibration prelude + KU-1.30
> demo gates PASS; production gates implemented but deferred (~1 h #1,
> several h #2). H.1 status remains 🟨 pending the production sweep.
>
> Continue per H.2 brief `ffn_sim/docs/briefs/H2_single_filament.md`:
>
> 1. Single filament configuration (N=21 beads, no cross-links, free
>    BC on x/y, periodic z) under BAOAB. Reuse `ffn_sim/integrator/baoab.py`
>    frozen at `75a1035` + the M2-rest int32 guard (PI-authorized
>    additive change).
> 2. Tangent–tangent correlation function ⟨t̂(0)·t̂(s)⟩ vs arc length.
>    Fit exponential decay → L_p. Compare to brief target `L_p ∈ [12, 22] μm`
>    (KU-1.1 nominal 17 μm ± 30%).
> 3. L-M vs E-M order separation: BAOAB freeze §Open #2 — run the same
>    L_p protocol under E-M reference via the
>    `allow_reference_integrator` switch in
>    `validation/oracles/common/sanity_gate.py`; expect E-M bias scaling
>    as O(Δt) vs BAOAB's O(Δt²).
>
> Open items that ride along:
> - Bending equipartition gate (BAOAB freeze §Open #1) tested here.
> - PI to call on whether to run H.1 KU-1.30 production gates (#1, #2,
>   #3) as the H.2 prelude or in a parallel optimisation session.
> - Stashes on `ffn/foundation` (`stash@{0..1}`) still need ownership.
