"""Self-consistency + known-physics tests for the linear viscoelastic S1/M2 oracle.

Validates ``aleph.validation.oracles.mechanics.viscoelastic_relaxation`` against
the analytic identities of the Maxwell, Kelvin-Voigt and SLS (Zener) models. Pure
Python, no GPU.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.validation.oracles.mechanics import viscoelastic_relaxation as ve

G0 = 1.0e3      # Pa
TAU = 5.0       # s
G_INF = 4.0e2   # Pa (equilibrium)
G1 = 6.0e2      # Pa (relaxing)


# --------------------------------------------------------------------------- #
# Maxwell
# --------------------------------------------------------------------------- #
def test_maxwell_relaxation_boundaries():
    assert ve.maxwell_relaxation(0.0, G0, TAU) == pytest.approx(G0)
    assert ve.maxwell_relaxation(TAU, G0, TAU) == pytest.approx(G0 / np.e)
    assert ve.maxwell_relaxation(1e6, G0, TAU) == pytest.approx(0.0, abs=1e-6)


def test_maxwell_creep_is_linear_fluid():
    # J(t) = 1/G0 + t/eta, eta = G0*tau -> slope 1/eta, intercept 1/G0.
    t = np.linspace(0.0, 100.0, 50)
    J = ve.maxwell_creep_compliance(t, G0, TAU)
    assert J[0] == pytest.approx(1.0 / G0)
    slope = np.polyfit(t, J, 1)[0]
    assert slope == pytest.approx(1.0 / (G0 * TAU))


# --------------------------------------------------------------------------- #
# Kelvin-Voigt
# --------------------------------------------------------------------------- #
def test_kelvin_voigt_creep_boundaries():
    E = 800.0
    assert ve.kelvin_voigt_creep(0.0, E, TAU) == pytest.approx(0.0)
    assert ve.kelvin_voigt_creep(TAU, E, TAU) == pytest.approx((1.0 / E) * (1.0 - 1.0 / np.e))
    assert ve.kelvin_voigt_creep(1e6, E, TAU) == pytest.approx(1.0 / E)


def test_kelvin_voigt_creep_monotone_increasing():
    t = np.linspace(0.0, 50.0, 100)
    J = ve.kelvin_voigt_creep(t, 800.0, TAU)
    assert np.all(np.diff(J) >= 0.0)


# --------------------------------------------------------------------------- #
# SLS (Zener)
# --------------------------------------------------------------------------- #
def test_sls_relaxation_boundaries():
    assert ve.sls_relaxation(0.0, G_INF, G1, TAU) == pytest.approx(G_INF + G1)  # glassy
    assert ve.sls_relaxation(TAU, G_INF, G1, TAU) == pytest.approx(G_INF + G1 / np.e)
    assert ve.sls_relaxation(1e6, G_INF, G1, TAU) == pytest.approx(G_INF)       # equilibrium


def test_sls_relaxation_monotone_non_increasing():
    t = np.linspace(0.0, 60.0, 200)
    G = ve.sls_relaxation(t, G_INF, G1, TAU)
    assert np.all(np.diff(G) <= 0.0)


def test_sls_creep_boundaries_and_glassy_equilibrium_consistency():
    J0 = ve.sls_creep_compliance(0.0, G_INF, G1, TAU)
    Jinf = ve.sls_creep_compliance(1e7, G_INF, G1, TAU)
    # J(0) = 1/G(0) glassy; J(inf) = 1/G(inf) equilibrium.
    assert J0 == pytest.approx(1.0 / (G_INF + G1))
    assert Jinf == pytest.approx(1.0 / G_INF)


def test_sls_creep_monotone_increasing():
    t = np.linspace(0.0, 200.0, 300)
    J = ve.sls_creep_compliance(t, G_INF, G1, TAU)
    assert np.all(np.diff(J) >= 0.0)


def test_retardation_time_exceeds_relaxation_time():
    tau_eps = ve.sls_retardation_time(G_INF, G1, TAU)
    assert tau_eps > TAU
    assert tau_eps == pytest.approx(TAU * (G_INF + G1) / G_INF)


def test_retardation_equals_relaxation_when_no_relaxing_arm():
    assert ve.sls_retardation_time(G_INF, 0.0, TAU) == pytest.approx(TAU)


# --------------------------------------------------------------------------- #
# Recovery from a measured curve
# --------------------------------------------------------------------------- #
def test_relaxation_time_recovered_from_ideal_sls_curve():
    t = np.linspace(0.0, 8.0 * TAU, 400)
    G = ve.sls_relaxation(t, G_INF, G1, TAU)
    rec = ve.relaxation_time_from_curve(t, G)
    assert rec["G_inf"] == pytest.approx(G_INF, rel=1e-3)
    assert rec["tau"] == pytest.approx(TAU, rel=1e-2)
    assert rec["G1"] == pytest.approx(G1, rel=2e-2)
    assert rec["G0"] == pytest.approx(G_INF + G1, rel=1e-2)
    assert rec["r_squared"] > 0.999


def test_recovery_rejects_non_decaying_curve():
    t = np.linspace(0.0, 10.0, 20)
    G = np.linspace(100.0, 200.0, 20)     # increasing
    with pytest.raises(ValueError):
        ve.relaxation_time_from_curve(t, G)


# --------------------------------------------------------------------------- #
# Guards
# --------------------------------------------------------------------------- #
def test_negative_time_raises():
    with pytest.raises(ValueError):
        ve.sls_relaxation(-1.0, G_INF, G1, TAU)


def test_nonpositive_params_raise():
    with pytest.raises(ValueError):
        ve.maxwell_relaxation(1.0, 0.0, TAU)
    with pytest.raises(ValueError):
        ve.sls_relaxation(1.0, G_INF, -1.0, TAU)
