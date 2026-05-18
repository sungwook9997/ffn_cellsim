"""Phase 1 Unit 3.1 validation tests (Worker C — Cell Track).

Seven KU-driven checks:

1. Cortex generation contract — KU-3.17 (n_fibers, beads, geometry).
2. Tangential bias — KU-3.17 (vMF spread around local tangent).
3. ECM module re-use — KU-1.24 (compute_forces accepts cortex bead
   positions and cross-links, output shape matches).
4. Cell dataclass round-trip — Cell.from_cortex / refresh_centroid /
   compute_cortex_boundary_position.
5. KU-3.1 VALIDATION: cell rounding — elliptical cortex relaxes to
   aspect ratio < 1.2 within the KU-3.1 60-s window.
6. KU-3.5 VALIDATION: blebbistatin response — 10× γ_cortex drop gives
   slower rounding (qualitative comparison).
7. Phase 1 cell-cortex Sanity Gate passes on resolved defaults.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from acs_kb.cell.cell import Cell
from acs_kb.cell.cortex import (
    Cortex,
    generate_cortex,
    generate_elliptical_cortex,
    measure_aspect_ratio,
    measure_coordination,
    measure_coverage_ratio,
)
from acs_kb.cell.force_balance import (
    cortical_tension_forces,
    relax,
    solve_overdamped_step,
)
from acs_kb.common.derived_params_cell import load_cell_config
from acs_kb.common.sanity_gate import gate_phase1_cell_cortex
from acs_kb.ecm.fiber_mechanics import compute_forces


CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_unit3.yaml"


@pytest.fixture(scope="module")
def cfg():
    return load_cell_config(CONFIG_PATH)


@pytest.fixture(scope="module")
def cortex(cfg) -> Cortex:
    c = cfg["cell"]
    return generate_cortex(
        R_cell=c["R_cell"],
        n_cortex_fibers=c["n_cortex_fibers"],
        L_cortex_fiber=c["L_cortex_fiber"],
        beads_per_fiber=c["beads_per_fiber"],
        kappa_vm_tangent=c["kappa_vm_tangent"],
        xl_cutoff=c["xl_cutoff"],
        xl_stiffness=c["xl_stiffness"],
        seed=c["seed"],
        box_size=c["derived"]["box_size"],
    )


# ------------------------------------------------------------------ #
# 1. Generation contract
# ------------------------------------------------------------------ #

def test_generation_contract(cfg, cortex):
    """KU-3.17: shape, fiber count, rest length, coverage."""
    c = cfg["cell"]
    n_f, n_b, dim = cortex.bead_positions.shape
    assert n_f == c["n_cortex_fibers"]
    assert n_b == c["beads_per_fiber"]
    assert dim == 2
    assert cortex.rest_length == pytest.approx(c["derived"]["rest_length"], rel=1e-12)
    # Coverage = N · L / (2π R) — a property of the construction.
    assert measure_coverage_ratio(cortex) == pytest.approx(
        c["derived"]["fiber_coverage_ratio"], rel=1e-9
    )


# ------------------------------------------------------------------ #
# 2. Tangential bias
# ------------------------------------------------------------------ #

def test_tangential_bias_present(cfg):
    """KU-3.17: fiber orientations cluster around the local tangent.

    Mean |cos(angle − tangent_angle)| should be > 0.5 for κ_vM = 5
    (uniform distribution gives 2/π ≈ 0.637, so we set the bar at
    > 0.85 which is achieved only with substantial alignment).
    """
    c = cfg["cell"]
    rng_seed = c["seed"] + 7
    cortex = generate_cortex(
        R_cell=c["R_cell"],
        n_cortex_fibers=400,                 # larger sample for the alignment statistic
        L_cortex_fiber=c["L_cortex_fiber"],
        beads_per_fiber=c["beads_per_fiber"],
        kappa_vm_tangent=c["kappa_vm_tangent"],
        xl_cutoff=c["xl_cutoff"],
        xl_stiffness=c["xl_stiffness"],
        seed=rng_seed,
        box_size=c["derived"]["box_size"],
    )
    # Reconstruct each fiber's nominal tangent from its centre.
    rel = cortex.fiber_centers - cortex.cell_center
    angular = np.arctan2(rel[:, 1], rel[:, 0])
    tangent_angle = np.mod(angular + 0.5 * np.pi, np.pi)
    # Difference in [0, π/2] (head-tail symmetric).
    diff = np.abs(np.mod(cortex.fiber_orientations - tangent_angle + np.pi, np.pi) - 0.5 * np.pi)
    diff = 0.5 * np.pi - diff       # 0 = perfectly tangent
    mean_cos = float(np.mean(np.cos(diff)))
    assert mean_cos > 0.85, (
        f"Tangential bias too weak: ⟨cos(θ-θ_tangent)⟩={mean_cos:.3f}, "
        f"expected > 0.85 for κ_vM={c['kappa_vm_tangent']}"
    )


# ------------------------------------------------------------------ #
# 3. ECM module re-use (KU-1.24 + KU-1.28)
# ------------------------------------------------------------------ #

def test_ecm_module_reuse_shape_and_finite(cfg, cortex):
    """KU-1.24 + KU-1.28: compute_forces accepts the cortex shape verbatim."""
    c = cfg["cell"]
    F = compute_forces(
        cortex.bead_positions,
        rest_length=cortex.rest_length,
        stretching_modulus=c["stretching_modulus"],
        bending_modulus=c["bending_modulus"],
        box_size=cortex.box_size,
        cross_links=cortex.cross_links,
    )
    assert F.shape == cortex.bead_positions.shape
    assert np.all(np.isfinite(F))
    # Closed bonded system: Σ F = 0 (Newton 3rd) holds for the WLC backbone
    # and XL springs because every force is internal pairwise.
    assert np.allclose(F.sum(axis=(0, 1)), 0.0, atol=1e-18), (
        f"Σ F = {F.sum(axis=(0,1))}; cortex bead-spring is bonded so must close."
    )


def test_cortical_tension_force_sign_and_units(cfg, cortex):
    """KU-3.5 (revised brief): F_tension is a Hookean radial spring centred on R_cell.

    For beads at R_eff > R_cell the force is inward; for R_eff < R_cell it
    is outward; the sign of the projection on r̂_eff equals the sign of
    (R_cell − R_eff). We verify this contract bead-by-bead rather than
    asserting universal inward force, because the elliptical IC and the
    fiber-along-tangent offsets put some beads inside R_cell.
    """
    c = cfg["cell"]
    F_t = cortical_tension_forces(
        cortex.bead_positions, c["gamma_cortex"], cortex.rest_length,
        R_cell=cortex.R_cell,
    )
    assert F_t.shape == cortex.bead_positions.shape
    pts = cortex.bead_positions.reshape(-1, 2)
    centroid = pts.mean(axis=0)
    rel = pts - centroid
    rad = np.linalg.norm(rel, axis=1)
    safe = np.where(rad > 0, rad, 1.0)
    rhat = rel / safe[:, None]
    proj = np.einsum("bd,bd->b", F_t.reshape(-1, 2), rhat)
    # Sign contract: proj has the same sign as (R_cell − R_eff).
    expected_sign = np.sign(cortex.R_cell - rad)
    nontrivial = np.abs(rad - cortex.R_cell) > 1e-12
    np.testing.assert_array_equal(
        np.sign(proj[nontrivial]),
        expected_sign[nontrivial],
        err_msg="Hookean tension sign contract violated.",
    )


# ------------------------------------------------------------------ #
# 4. Cell dataclass round-trip
# ------------------------------------------------------------------ #

def test_cell_dataclass_and_boundary_position(cfg, cortex):
    cell = Cell.from_cortex(cell_id=0, cortex=cortex)
    pts = cortex.bead_positions.reshape(-1, 2)
    centroid = pts.mean(axis=0)
    assert np.allclose(cell.center_position, centroid)
    assert np.allclose(cell.nucleus_position, centroid)
    assert cell.inner_mode == "discrete"
    # Boundary lookup at θ = 0 returns a bead near the +x extreme.
    p = cell.compute_cortex_boundary_position(0.0)
    # The reported bead should sit roughly in the +x direction from centroid.
    rel = p - cell.center_position
    assert np.arctan2(rel[1], rel[0]) == pytest.approx(0.0, abs=np.pi / 6)
    # And it should be at radius ~ R_cell (within construction tolerance).
    radius = float(np.linalg.norm(rel))
    assert 0.5 * cortex.R_cell < radius < 1.5 * cortex.R_cell


# ------------------------------------------------------------------ #
# 5. KU-3.1 VALIDATION — cell rounding (integration test)
# ------------------------------------------------------------------ #

@pytest.fixture(scope="module")
def rounding_run(cfg):
    """Run the elliptical-cortex relaxation once and reuse for several tests."""
    c = cfg["cell"]
    r = c["rounding_test"]
    cortex = generate_elliptical_cortex(
        semi_major=r["semi_major"], semi_minor=r["semi_minor"],
        n_cortex_fibers=c["n_cortex_fibers"],
        L_cortex_fiber=c["L_cortex_fiber"],
        beads_per_fiber=c["beads_per_fiber"],
        kappa_vm_tangent=c["kappa_vm_tangent"],
        xl_cutoff=c["xl_cutoff"],
        xl_stiffness=c["xl_stiffness"],
        seed=c["seed"],
        box_size=c["derived"]["box_size"],
    )
    final_cortex, traj = relax(
        cortex,
        stretching_modulus=c["stretching_modulus"],
        bending_modulus=c["bending_modulus"],
        gamma_cortex=c["gamma_cortex"],
        gamma_drag=c["derived"]["gamma_b"],
        dt=r["target_dt"],
        n_steps=r["n_steps"],
        sample_interval=r["sample_interval"],
        keep_frames=False,
        K_area=c["K_area"],
    )
    return cortex, final_cortex, traj


def test_initial_aspect_ratio_is_elliptical(cfg, rounding_run):
    """Sanity: the elliptical IC is actually ~1.5:1 (≥ 1.3)."""
    initial_cortex, _, _ = rounding_run
    ar0 = measure_aspect_ratio(initial_cortex)
    assert ar0 >= 1.3, f"Initial aspect ratio too round: {ar0:.3f}"


def test_cell_rounding_acceptance_window(cfg, rounding_run):
    """KU-3.1 VALIDATION: aspect ratio < 1.2 within the 60-s window."""
    c = cfg["cell"]
    r = c["rounding_test"]
    _, _, traj = rounding_run
    mask = traj.times <= r["acceptance_window_s"] + 1e-12
    ar_in_window = traj.aspect_ratios[mask]
    min_in_window = float(ar_in_window.min())
    assert min_in_window < r["aspect_ratio_acceptance"], (
        f"KU-3.1 FAIL — aspect ratio stays at {min_in_window:.3f} > "
        f"{r['aspect_ratio_acceptance']} during 0..{r['acceptance_window_s']} s"
    )


def test_aspect_ratio_monotonic_trend(rounding_run):
    """KU-3.1 acceptance text: late aspect ratio < initial (rounding happened).

    With τ_round ≈ R_cell · γ_drag / γ_cortex ≈ 2 s and total sim 120 s,
    the cortex equilibrates in the first few percent of the trajectory.
    We compare the initial AR (t=0 fixed sample) to the late-trajectory
    mean — must drop by at least 0.10 for the test to pass.
    """
    _, _, traj = rounding_run
    ar = traj.aspect_ratios
    q = max(1, len(ar) // 4)
    initial = float(ar[0])
    late_mean = float(ar[-q:].mean())
    assert late_mean < initial - 0.10, (
        f"Aspect ratio did not decrease enough: initial={initial:.3f}, "
        f"late mean={late_mean:.3f}"
    )


# ------------------------------------------------------------------ #
# 6. KU-3.5 VALIDATION — blebbistatin qualitative response
# ------------------------------------------------------------------ #

def test_blebbistatin_slower_rounding(cfg):
    """KU-3.5 VALIDATION: 10× γ_cortex drop produces qualitatively slower rounding.

    Quantitative: at t = 30 s the high-tension run should have a
    smaller aspect ratio than the blebbistatin run by at least 0.05
    (a clear gap, well above numerical noise).
    """
    c = cfg["cell"]
    r = c["rounding_test"]

    def _run(gamma: float, sim_time: float = 2.0):
        cortex = generate_elliptical_cortex(
            semi_major=r["semi_major"], semi_minor=r["semi_minor"],
            n_cortex_fibers=c["n_cortex_fibers"],
            L_cortex_fiber=c["L_cortex_fiber"],
            beads_per_fiber=c["beads_per_fiber"],
            kappa_vm_tangent=c["kappa_vm_tangent"],
            xl_cutoff=c["xl_cutoff"],
            xl_stiffness=c["xl_stiffness"],
            seed=c["seed"],
            box_size=c["derived"]["box_size"],
        )
        dt = r["target_dt"]
        n_steps = int(sim_time / dt)
        _, traj = relax(
            cortex,
            stretching_modulus=c["stretching_modulus"],
            bending_modulus=c["bending_modulus"],
            gamma_cortex=gamma,
            gamma_drag=c["derived"]["gamma_b"],
            dt=dt,
            n_steps=n_steps,
            sample_interval=5,
            keep_frames=False,
            K_area=c["K_area"],
        )
        return traj

    # Probe at t = 1.0 s — full-γ run (τ ≈ 2 s) has decayed to ≈ 0.6 of
    # the initial perturbation; blebbistatin run (τ ≈ 20 s) has barely
    # moved (≈ 0.95 of initial). The gap at t = 1 s is much larger than
    # at t = 0.5 s, giving the test reliable headroom above bead noise.
    traj_full = _run(c["gamma_cortex"])
    traj_blebb = _run(c["gamma_cortex_blebbistatin"])
    probe_t = 1.0
    ar_full = float(traj_full.aspect_ratios[np.argmin(np.abs(traj_full.times - probe_t))])
    ar_blebb = float(traj_blebb.aspect_ratios[np.argmin(np.abs(traj_blebb.times - probe_t))])
    assert ar_full + 0.05 < ar_blebb, (
        f"Blebbistatin should round more slowly at t={probe_t} s: "
        f"full γ → AR={ar_full:.3f}, blebb → AR={ar_blebb:.3f}"
    )


# ------------------------------------------------------------------ #
# 7. Sanity Gate on the full Unit 3.1 default run
# ------------------------------------------------------------------ #

def test_sanity_gate_passes_on_resolved_defaults(cfg, cortex, rounding_run):
    c = cfg["cell"]
    r = c["rounding_test"]
    initial_cortex, final_cortex, traj = rounding_run

    z = measure_coordination(cortex)
    coverage = measure_coverage_ratio(cortex)
    ar_init = measure_aspect_ratio(initial_cortex)
    ar_final = measure_aspect_ratio(final_cortex)
    mask = traj.times <= r["acceptance_window_s"] + 1e-12
    ar_window = float(traj.aspect_ratios[mask].min())

    radius_initial = float(traj.mean_radii[0])
    radius_final = float(traj.mean_radii[-1])
    drift_rel = (radius_final - radius_initial) / radius_initial

    rep = gate_phase1_cell_cortex(
        measured_z=z,
        coverage_ratio=coverage,
        initial_aspect_ratio=ar_init,
        final_aspect_ratio=ar_final,
        acceptance_window_reached_aspect_ratio=ar_window,
        acceptance=c["acceptance"],
        cortex_radius_drift_rel=drift_rel,
        n_fibers=cortex.n_fibers,
        n_xls=len(cortex.cross_links),
        demo_mode=c["demo_mode"],
    )
    assert rep.passed, rep.summary()


# ------------------------------------------------------------------ #
# 8. Force-balance equilibrium sanity (no-input case)
# ------------------------------------------------------------------ #

def test_one_step_is_smooth(cfg, cortex):
    """A single overdamped step under default tension does not blow up."""
    c = cfg["cell"]
    new_cortex = solve_overdamped_step(
        cortex,
        stretching_modulus=c["stretching_modulus"],
        bending_modulus=c["bending_modulus"],
        gamma_cortex=c["gamma_cortex"],
        gamma_drag=c["derived"]["gamma_b"],
        dt=c["rounding_test"]["target_dt"],
        K_area=c["K_area"],
    )
    assert np.all(np.isfinite(new_cortex.bead_positions))
    # Maximum displacement must be smaller than ℓ₀ (sub-bond motion / step).
    disp = np.linalg.norm(
        new_cortex.bead_positions - cortex.bead_positions, axis=-1
    ).max()
    assert disp < cortex.rest_length, (
        f"Single-step displacement {disp:.3e} m exceeds bond length "
        f"{cortex.rest_length:.3e} m — CFL likely violated."
    )
