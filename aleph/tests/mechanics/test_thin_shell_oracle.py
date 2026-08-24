"""Self-consistency + known-physics tests for the pressurized thin-shell S1/S2 oracle."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.validation.oracles.mechanics import thin_shell as ts

R = 7.5e-6        # m
H = 0.2e-6        # m (thin cortex, h/R ~ 0.027)
E = 1.0e3         # Pa
NU = 0.45
GAMMA = 1.0e-4    # N/m
P = 40.0          # Pa (resting turgor)


def test_laplace_and_inverse_are_consistent():
    dP = float(ts.laplace_pressure(GAMMA, R))
    assert dP == pytest.approx(2.0 * GAMMA / R)
    assert ts.tension_from_laplace(dP, R) == pytest.approx(GAMMA)


def test_cortical_tension_equals_membrane_stress_times_thickness():
    tension = ts.cortical_tension(P, R)
    sigma = float(ts.thin_wall_membrane_stress(P, R, H))
    assert tension == pytest.approx(P * R / 2.0)
    assert tension == pytest.approx(sigma * H)


def test_cortical_tension_matches_laplace_gamma():
    # A shell at pressure p has tension pR/2, which is exactly the Laplace gamma.
    tension = ts.cortical_tension(P, R)
    assert float(ts.laplace_pressure(tension, R)) == pytest.approx(P)


def test_hoop_strain_formula_and_volume_is_three_times():
    eps = float(ts.thin_shell_hoop_strain(P, R, H, E, NU))
    assert eps == pytest.approx(P * R * (1.0 - NU) / (2.0 * H * E))
    assert float(ts.thin_shell_radial_expansion(P, R, H, E, NU)) == pytest.approx(R * eps)
    assert float(ts.thin_shell_volume_change(P, R, H, E, NU)) == pytest.approx(3.0 * eps)


def test_zero_pressure_zero_response():
    assert float(ts.thin_wall_membrane_stress(0.0, R, H)) == pytest.approx(0.0)
    assert float(ts.thin_shell_hoop_strain(0.0, R, H, E, NU)) == pytest.approx(0.0)


def test_membrane_stress_scales_linearly_with_pressure():
    p = np.array([10.0, 20.0, 40.0])
    sigma = ts.thin_wall_membrane_stress(p, R, H)
    assert sigma[1] / sigma[0] == pytest.approx(2.0)
    assert sigma[2] / sigma[0] == pytest.approx(4.0)


def test_reissner_stiffness_formula_and_scaling():
    k = ts.reissner_point_stiffness(E, H, R, NU)
    assert k == pytest.approx(4.0 * E * H * H / (R * np.sqrt(3.0 * (1.0 - NU**2))))
    # Doubling thickness quadruples stiffness (h^2); doubling E doubles it.
    assert ts.reissner_point_stiffness(E, 2 * H, R, NU) == pytest.approx(4.0 * k)
    assert ts.reissner_point_stiffness(2 * E, H, R, NU) == pytest.approx(2.0 * k)


def test_thin_wall_validity_window():
    assert ts.valid_thin_wall(H, R)["valid"] is True
    assert ts.valid_thin_wall(2.0e-6, R)["valid"] is False   # h/R ~ 0.27


def test_guards():
    with pytest.raises(ValueError):
        ts.laplace_pressure(-1.0, R)
    with pytest.raises(ValueError):
        ts.thin_shell_hoop_strain(P, R, H, E, 1.0)
    with pytest.raises(ValueError):
        ts.thin_wall_membrane_stress(-1.0, R, H)
