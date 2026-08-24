"""Structure gates for the 3-model ERM membrane--cortex connector spectrum.

Mac-side: host acceptance oracles + the source-gated per-model parameter contract + spectrum-continuity
algebra.  The CUDA force/commit kernels are exercised under ``importorskip`` and skip without a GPU (I0-A) —
native validation belongs to the GATE A session that owns the A5000.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.incumbent.erm_spectrum import (
    OFF_RATE_BELL_SLIP,
    OFF_RATE_PEREVERZEV_CATCH_SLIP,
    ERMSpectrumModel,
    ERMSpectrumParams,
)
from aleph.components.incumbent.erm_spectrum_analytic import (
    erm_bell_slip_off_rate,
    erm_bond_density_rhs,
    erm_bond_density_steady_state,
    erm_bond_density_tick,
    erm_catch_slip_off_rate,
    erm_catch_slip_peak_force,
    erm_elastic_friction_axial_force,
    erm_rigid_weld_axial_force,
)

# ---------------------------------------------------------------------------
# model 1 -- RIGID_WELD (완전결합)
# ---------------------------------------------------------------------------


def test_rigid_weld_is_bilateral_and_zero_at_rest() -> None:
    length = np.array([0.8, 1.0, 1.2], dtype=np.float64)
    f = erm_rigid_weld_axial_force(length, k_weld=10.0, rest_length=1.0)
    assert f[0] == pytest.approx(-2.0)  # compression pushes apart (bilateral)
    assert f[1] == pytest.approx(0.0)
    assert f[2] == pytest.approx(+2.0)  # extension pulls together


def test_rigid_weld_stiffening_drives_displacement_to_zero() -> None:
    # A fixed tensile load L-rest carried by an ever-stiffer weld demands an ever-smaller extension:
    # for a target force, extension = F/k_weld -> 0 as k_weld -> inf (the Δu = 0 co-motion limit).
    target = 5.0
    exts = [target / k for k in (1e2, 1e4, 1e6)]
    assert exts[0] > exts[1] > exts[2]
    for k, ext in zip((1e2, 1e4, 1e6), exts):
        assert erm_rigid_weld_axial_force(1.0 + ext, k_weld=k, rest_length=1.0) == pytest.approx(target)


# ---------------------------------------------------------------------------
# model 2 -- ELASTIC_FRICTION (탄성·마찰)
# ---------------------------------------------------------------------------


def test_elastic_spring_is_unilateral_and_density_scaled() -> None:
    length = np.array([0.9, 1.0, 1.2], dtype=np.float64)
    zero_rate = np.zeros_like(length)
    f = erm_elastic_friction_axial_force(length, zero_rate, k_b=10.0, rho_b=2.0, xi=1.0, rest_length=1.0)
    np.testing.assert_array_equal(f[:2], 0.0)          # compression + rest: no tension (ezrin cannot push)
    assert f[2] == pytest.approx(2.0 * 10.0 * 0.2)     # rho_b · k_b · extension


def test_friction_opposes_sliding_in_either_direction() -> None:
    # At rest length the spring is zero, so the axial force is pure friction xi·du_dot.
    f_out = erm_elastic_friction_axial_force(1.0, +0.5, k_b=10.0, rho_b=1.0, xi=4.0, rest_length=1.0)
    f_in = erm_elastic_friction_axial_force(1.0, -0.5, k_b=10.0, rho_b=1.0, xi=4.0, rest_length=1.0)
    assert float(f_out) == pytest.approx(+2.0)
    assert float(f_in) == pytest.approx(-2.0)


def test_elastic_friction_superposes_spring_and_dashpot() -> None:
    f = erm_elastic_friction_axial_force(1.2, 0.5, k_b=10.0, rho_b=1.0, xi=4.0, rest_length=1.0)
    assert float(f) == pytest.approx(1.0 * 10.0 * 0.2 + 4.0 * 0.5)


# ---------------------------------------------------------------------------
# off-rate laws (Bell slip vs Pereverzev catch/slip)
# ---------------------------------------------------------------------------


def test_bell_slip_is_monotone_increasing() -> None:
    loads = np.array([0.0, 1.0, 2.0, 4.0], dtype=np.float64)
    rate = erm_bell_slip_off_rate(loads, k_off0=0.5, bell_force=3.0)
    assert rate[0] == pytest.approx(0.5)
    assert np.all(np.diff(rate) > 0.0)


def test_catch_slip_falls_then_rises_with_a_peak() -> None:
    kw = dict(k_catch0=0.4, x_catch_um=0.61e-3, k_slip0=0.5, x_slip_um=0.14e-3, kT=4.28e-3)
    f_star = erm_catch_slip_peak_force(**kw)
    assert f_star > 0.0
    lo = float(erm_catch_slip_off_rate(f_star - 0.5 * f_star, **kw))
    at = float(erm_catch_slip_off_rate(f_star, **kw))
    hi = float(erm_catch_slip_off_rate(f_star + 0.5 * f_star, **kw))
    assert at < lo and at < hi  # F* is the off-rate minimum (catch -> slip crossover)


def test_catch_slip_peak_requires_catch_dominance() -> None:
    with pytest.raises(ValueError, match="slip-only"):
        erm_catch_slip_peak_force(k_catch0=0.1, x_catch_um=0.1e-3, k_slip0=0.5,
                                  x_slip_um=0.14e-3, kT=4.28e-3)


# ---------------------------------------------------------------------------
# model 3 -- DYNAMIC_CLUTCH (동적결합) bond-density kinetics
# ---------------------------------------------------------------------------


def test_bond_density_rhs_vanishes_at_steady_state() -> None:
    k_off = np.array([0.2, 1.0, 5.0], dtype=np.float64)
    rho_ss = erm_bond_density_steady_state(k_on=1.0, rho_max=1.0, k_off=k_off)
    rhs = erm_bond_density_rhs(rho_ss, k_on=1.0, rho_max=1.0, k_off=k_off)
    np.testing.assert_allclose(rhs, 0.0, atol=1e-15)


def test_steady_state_brackets_the_occupancy() -> None:
    # rho_ss -> rho_max as k_off -> 0 (all bound); rho_ss -> 0 as k_off -> inf (all stripped).
    assert float(erm_bond_density_steady_state(k_on=1.0, rho_max=1.0, k_off=1e-9)) == pytest.approx(1.0, abs=1e-6)
    assert float(erm_bond_density_steady_state(k_on=1.0, rho_max=1.0, k_off=1e9)) == pytest.approx(0.0, abs=1e-6)


def test_exact_tick_relaxes_monotonically_to_steady_state() -> None:
    rho = 0.1
    ss = float(erm_bond_density_steady_state(k_on=1.0, rho_max=1.0, k_off=0.5))
    prev = rho
    for _ in range(200):
        rho = float(erm_bond_density_tick(rho, k_on=1.0, rho_max=1.0, k_off=0.5, tau=0.05))
        assert prev <= rho <= ss + 1e-12  # monotone up, never overshoots the steady state
        prev = rho
    assert rho == pytest.approx(ss, abs=1e-6)


def test_exact_tick_limits() -> None:
    assert float(erm_bond_density_tick(0.3, k_on=1.0, rho_max=1.0, k_off=0.5, tau=0.0)) == pytest.approx(0.3)
    huge = float(erm_bond_density_tick(0.3, k_on=1.0, rho_max=1.0, k_off=0.5, tau=1e6))
    ss = float(erm_bond_density_steady_state(k_on=1.0, rho_max=1.0, k_off=0.5))
    assert huge == pytest.approx(ss)


# ---------------------------------------------------------------------------
# spectrum continuity — the models are limits of one another
# ---------------------------------------------------------------------------


def test_dynamic_clutch_reduces_to_elastic_spring_when_fully_bound() -> None:
    # k_off -> 0 => rho -> rho_max; the model-3 force form is then model-2 with rho_b = rho_max.
    rho_ss = float(erm_bond_density_steady_state(k_on=1.0, rho_max=1.5, k_off=1e-9))
    assert rho_ss == pytest.approx(1.5, abs=1e-6)
    f3 = erm_elastic_friction_axial_force(1.2, 0.0, k_b=10.0, rho_b=rho_ss, xi=1.0, rest_length=1.0)
    f2 = erm_elastic_friction_axial_force(1.2, 0.0, k_b=10.0, rho_b=1.5, xi=1.0, rest_length=1.0)
    assert float(f3) == pytest.approx(float(f2), rel=1e-5)


def test_elastic_spring_reduces_to_weld_tensile_branch() -> None:
    # rho_b·k_b == k_weld, no friction, extension => the unilateral spring equals the weld's tensile branch.
    f2 = erm_elastic_friction_axial_force(1.2, 0.0, k_b=5.0, rho_b=2.0, xi=1.0, rest_length=1.0)
    f1 = erm_rigid_weld_axial_force(1.2, k_weld=10.0, rest_length=1.0)
    assert float(f2) == pytest.approx(float(f1))


# ---------------------------------------------------------------------------
# source-gated per-model parameter contract
# ---------------------------------------------------------------------------


def test_params_reject_missing_source() -> None:
    with pytest.raises(ValueError, match="provenance"):
        ERMSpectrumParams(model=ERMSpectrumModel.RIGID_WELD, source="  ", k_weld_pn_um=10.0)


def test_rigid_weld_requires_only_stiffness() -> None:
    p = ERMSpectrumParams(model=ERMSpectrumModel.RIGID_WELD, source="TEST", k_weld_pn_um=10.0)
    assert p.model is ERMSpectrumModel.RIGID_WELD
    with pytest.raises(ValueError, match="k_weld_pn_um"):
        ERMSpectrumParams(model=ERMSpectrumModel.RIGID_WELD, source="TEST")


def test_elastic_friction_requires_all_clutch_fields() -> None:
    ERMSpectrumParams(model=ERMSpectrumModel.ELASTIC_FRICTION, source="TEST",
                      k_b_pn_um=10.0, rho_b=1.0, xi_pn_s_um=2.0)
    with pytest.raises(ValueError, match="xi_pn_s_um"):
        ERMSpectrumParams(model=ERMSpectrumModel.ELASTIC_FRICTION, source="TEST",
                          k_b_pn_um=10.0, rho_b=1.0)


def test_dynamic_clutch_bell_and_pereverzev_contracts() -> None:
    ERMSpectrumParams(model=ERMSpectrumModel.DYNAMIC_CLUTCH, source="TEST",
                      k_b_pn_um=10.0, rho_b=0.5, xi_pn_s_um=2.0,
                      k_on_s=1.0, rho_max=1.0, off_rate_mode=OFF_RATE_BELL_SLIP,
                      k_off0_s=0.5, bell_force_pn=3.0)
    ERMSpectrumParams(model=ERMSpectrumModel.DYNAMIC_CLUTCH, source="TEST",
                      k_b_pn_um=10.0, rho_b=0.5, xi_pn_s_um=2.0,
                      k_on_s=1.0, rho_max=1.0, off_rate_mode=OFF_RATE_PEREVERZEV_CATCH_SLIP,
                      k_catch0_s=0.4, x_catch_um=0.61e-3, k_slip0_s=0.5, x_slip_um=0.14e-3, kT_pn_um=4.28e-3)
    # Bell law selected but slip datum absent -> rejected.
    with pytest.raises(ValueError, match="k_off0_s"):
        ERMSpectrumParams(model=ERMSpectrumModel.DYNAMIC_CLUTCH, source="TEST",
                          k_b_pn_um=10.0, rho_b=0.5, xi_pn_s_um=2.0,
                          k_on_s=1.0, rho_max=1.0, off_rate_mode=OFF_RATE_BELL_SLIP, bell_force_pn=3.0)


def test_dynamic_clutch_rejects_overfull_and_bad_mode() -> None:
    with pytest.raises(ValueError, match="rho_b must not exceed rho_max"):
        ERMSpectrumParams(model=ERMSpectrumModel.DYNAMIC_CLUTCH, source="TEST",
                          k_b_pn_um=10.0, rho_b=2.0, xi_pn_s_um=2.0,
                          k_on_s=1.0, rho_max=1.0, off_rate_mode=OFF_RATE_BELL_SLIP,
                          k_off0_s=0.5, bell_force_pn=3.0)
    with pytest.raises(ValueError, match="off_rate_mode"):
        ERMSpectrumParams(model=ERMSpectrumModel.DYNAMIC_CLUTCH, source="TEST",
                          k_b_pn_um=10.0, rho_b=0.5, xi_pn_s_um=2.0,
                          k_on_s=1.0, rho_max=1.0, off_rate_mode=7,
                          k_off0_s=0.5, bell_force_pn=3.0)


# ---------------------------------------------------------------------------
# CUDA kernels (I0-A: native gate belongs to the GATE A A5000 session)
# ---------------------------------------------------------------------------


def _cuda_device():
    wp = pytest.importorskip("warp")
    wp.init()
    device = next((d for d in wp.get_devices() if d.is_cuda), None)
    if device is None:
        pytest.skip("I0-A: ERM spectrum kernel gate requires a CUDA GPU")
    return wp, str(device)


def test_cuda_weld_force_is_bilateral_and_newton_third() -> None:
    wp, device = _cuda_device()
    from aleph.components.incumbent.erm_spectrum import erm_weld_force_kernel

    pos = np.array([[1.2, 0, 0], [0.8, 0, 0], [1.0, 0, 0],
                    [0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=np.float64)
    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    membrane_d = wp.array(np.arange(3, dtype=np.int32), dtype=wp.int32, device=device)
    cortex_d = wp.array(np.arange(3, 6, dtype=np.int32), dtype=wp.int32, device=device)
    bound_d = wp.ones(3, dtype=wp.int32, device=device)
    rest_d = wp.array(np.ones(3, dtype=np.float64), dtype=wp.float64, device=device)
    force_d = wp.zeros(6, dtype=wp.vec3d, device=device)
    wp.launch(erm_weld_force_kernel, dim=3,
              inputs=[pos_d, membrane_d, cortex_d, bound_d, wp.float64(10.0), rest_d],
              outputs=[force_d], device=device)
    f = force_d.numpy()
    np.testing.assert_allclose(f[0], [-2.0, 0, 0], atol=1e-14)  # extended: membrane pulled in
    np.testing.assert_allclose(f[1], [+2.0, 0, 0], atol=1e-14)  # compressed: membrane pushed out (bilateral)
    np.testing.assert_allclose(f[2], 0.0, atol=1e-14)
    np.testing.assert_allclose(f.sum(axis=0), 0.0, atol=1e-14)


def test_cuda_clutch_force_spring_plus_friction() -> None:
    wp, device = _cuda_device()
    from aleph.components.incumbent.erm_spectrum import erm_clutch_force_kernel

    pos = np.array([[1.2, 0, 0], [1.0, 0, 0],
                    [0, 0, 0], [0, 0, 0]], dtype=np.float64)
    vel = np.array([[0.5, 0, 0], [0.5, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=np.float64)
    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    vel_d = wp.array(vel, dtype=wp.vec3d, device=device)
    membrane_d = wp.array(np.arange(2, dtype=np.int32), dtype=wp.int32, device=device)
    cortex_d = wp.array(np.arange(2, 4, dtype=np.int32), dtype=wp.int32, device=device)
    bound_d = wp.ones(2, dtype=wp.int32, device=device)
    rho_d = wp.array(np.ones(2, dtype=np.float64), dtype=wp.float64, device=device)
    rest_d = wp.array(np.ones(2, dtype=np.float64), dtype=wp.float64, device=device)
    force_d = wp.zeros(4, dtype=wp.vec3d, device=device)
    wp.launch(erm_clutch_force_kernel, dim=2,
              inputs=[pos_d, vel_d, membrane_d, cortex_d, bound_d, wp.float64(10.0), rho_d,
                      wp.float64(4.0), rest_d],
              outputs=[force_d], device=device)
    f = force_d.numpy()
    # pair 0: spring 1·10·0.2 = 2.0, friction 4·0.5 = 2.0 -> axial 4.0 (membrane pulled in)
    np.testing.assert_allclose(f[0], [-4.0, 0, 0], atol=1e-13)
    # pair 1: at rest length spring 0, friction only 4·0.5 = 2.0
    np.testing.assert_allclose(f[1], [-2.0, 0, 0], atol=1e-13)
    np.testing.assert_allclose(f.sum(axis=0), 0.0, atol=1e-13)


def test_cuda_bond_density_commit_is_accepted_gated_and_matches_oracle() -> None:
    wp, device = _cuda_device()
    from aleph.components.incumbent.erm_spectrum import erm_bond_density_commit_kernel

    k_b, rest, ext = 10.0, 1.0, 0.2
    pos = np.array([[rest + ext, 0, 0], [0, 0, 0]], dtype=np.float64)
    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    membrane_d = wp.array(np.array([0], dtype=np.int32), dtype=wp.int32, device=device)
    cortex_d = wp.array(np.array([1], dtype=np.int32), dtype=wp.int32, device=device)
    bound_d = wp.ones(1, dtype=wp.int32, device=device)
    rest_d = wp.array(np.array([rest], dtype=np.float64), dtype=wp.float64, device=device)
    rho_d = wp.array(np.array([0.1], dtype=np.float64), dtype=wp.float64, device=device)
    accepted_d = wp.zeros(1, dtype=wp.int32, device=device)

    k_on, rho_max, k_off0, bell_force, tau = 1.0, 1.0, 0.5, 3.0, 0.2
    inputs = [pos_d, membrane_d, cortex_d, bound_d, wp.float64(k_b), rest_d, rho_d,
              wp.float64(k_on), wp.float64(rho_max), wp.int32(OFF_RATE_BELL_SLIP),
              wp.float64(k_off0), wp.float64(bell_force),
              wp.float64(1.0), wp.float64(1.0), wp.float64(1.0), wp.float64(1.0), wp.float64(1.0),
              wp.float64(tau), accepted_d]

    wp.launch(erm_bond_density_commit_kernel, dim=1, inputs=inputs, device=device)
    assert float(rho_d.numpy()[0]) == pytest.approx(0.1)  # rejected transaction: exact no-op

    accepted_d.fill_(1)
    wp.launch(erm_bond_density_commit_kernel, dim=1, inputs=inputs, device=device)
    f_bond = k_b * ext
    k_off = float(erm_bell_slip_off_rate(f_bond, k_off0=k_off0, bell_force=bell_force))
    expected = float(erm_bond_density_tick(0.1, k_on=k_on, rho_max=rho_max, k_off=k_off, tau=tau))
    assert float(rho_d.numpy()[0]) == pytest.approx(expected, rel=1e-12)
