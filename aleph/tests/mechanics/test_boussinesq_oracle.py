"""Self-consistency + known-physics tests for the Boussinesq S2 spatial-decay oracle."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.validation.oracles.mechanics import boussinesq as bq

P = 1.0e-9        # N
E = 1.0e3         # Pa
NU = 0.45


def test_surface_settlement_formula_and_inverse_r_decay():
    r = np.array([1.0, 2.0, 4.0]) * 1e-6
    u = bq.surface_settlement(P, E, NU, r)
    assert u[0] == pytest.approx(P * (1.0 - NU**2) / (np.pi * E * r[0]))
    # u_z ~ 1/r : doubling r halves the settlement.
    assert u[1] / u[0] == pytest.approx(0.5)
    assert u[2] / u[0] == pytest.approx(0.25)


def test_axial_stress_formula_and_inverse_square_decay():
    z = np.array([1.0, 2.0, 4.0]) * 1e-6
    s = bq.sigma_zz_axis(P, z)
    assert s[0] == pytest.approx(3.0 * P / (2.0 * np.pi * z[0] ** 2))
    # sigma_zz ~ 1/z^2 : doubling z quarters the stress.
    assert s[1] / s[0] == pytest.approx(0.25)
    assert s[2] / s[0] == pytest.approx(1.0 / 16.0)


def test_general_field_reduces_to_axis_at_r_zero():
    z = np.array([1.0, 2.0, 3.0]) * 1e-6
    on_axis = bq.sigma_zz(P, 0.0, z)
    np.testing.assert_allclose(on_axis, bq.sigma_zz_axis(P, z), rtol=1e-12)


def test_decay_exponent_recovered_from_field_loglog_fit():
    # Along a fixed ray, |sigma_zz| ~ 1/rho^2. Fit the exponent from a log-log slope.
    z = np.geomspace(1e-6, 1e-4, 40)
    r = z                                   # 45-degree ray
    rho = np.sqrt(r**2 + z**2)
    s = bq.sigma_zz(P, r, z)
    slope = np.polyfit(np.log(rho), np.log(s), 1)[0]
    assert slope == pytest.approx(-bq.STRESS_DECAY_EXPONENT, abs=1e-9)


def test_exponent_constants():
    assert bq.STRESS_DECAY_EXPONENT == 2.0
    assert bq.DISPLACEMENT_DECAY_EXPONENT == 1.0


def test_guards():
    with pytest.raises(ValueError):
        bq.surface_settlement(P, E, NU, 0.0)
    with pytest.raises(ValueError):
        bq.sigma_zz_axis(P, 0.0)
    with pytest.raises(ValueError):
        bq.sigma_zz(P, -1e-6, 1e-6)
