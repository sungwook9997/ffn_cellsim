"""Self-test of the LINC tether + pre-authored I7 load-path oracles (pure NumPy — no Warp/CUDA).

LINC: nonlinear tension-only nesprin link, force-free at rest, f=dE/dx. I7 (pre-authored): the capstan
∮T·κ ds tangential-tension→normal-pressure relation reproduces the classical belt-over-cylinder net
force 2T; oblate flatten conserves nucleoplasm volume (a∝A^{1/3}).
"""
from __future__ import annotations

import numpy as np
import pytest

from aleph.components.nucleus.linc_analytic import (
    capstan_line_integral,
    capstan_normal_pressure,
    linc_tether_energy,
    linc_tether_force,
    oblate_aspect_from_flatten,
    oblate_equatorial_radius_at_constant_volume,
)


# ---- LINC tether ----------------------------------------------------------------------------------
def test_linc_force_matches_fd_energy() -> None:
    """LINC tension f = dE/dΔ (the sign arbiter), for the cubic-stiffening tether."""
    x = np.linspace(1.0, 1.6, 50)
    f_ana = linc_tether_force(x, rest=1.0, k_linc=50.0, stiffening=2.0)
    eps = 1e-7
    f_fd = (linc_tether_energy(x + eps, 1.0, 50.0, 2.0) - linc_tether_energy(x - eps, 1.0, 50.0, 2.0)) / (2 * eps)
    assert np.max(np.abs(f_ana - f_fd)) < 1e-5 * np.max(np.abs(f_ana))


def test_linc_force_free_at_rest_and_tension_only() -> None:
    """f(rest)=0; a compressed tether (L<rest) transmits NO force (a tether pulls, never pushes)."""
    assert float(linc_tether_force(np.array([1.0]), 1.0, 50.0, 2.0)[0]) == pytest.approx(0.0, abs=1e-12)
    assert float(linc_tether_force(np.array([0.8]), 1.0, 50.0, 2.0)[0]) == pytest.approx(0.0, abs=1e-12)


def test_linc_monotone_and_stiffening() -> None:
    """Tension is monotone in extension and the cubic term stiffens it above the linear reference."""
    x = np.linspace(1.0, 1.5, 100)
    f = linc_tether_force(x, 1.0, 50.0, 2.0)
    assert np.all(np.diff(f) >= 0.0)
    f_lin = linc_tether_force(x, 1.0, 50.0, 0.0)
    assert f[-1] > f_lin[-1]                                     # stiffening raises the large-Δ tension


# ---- capstan (I7 pre-authored) --------------------------------------------------------------------
def test_capstan_normal_pressure() -> None:
    """p⊥ = T·κ (tangential tension over curvature κ=1/R becomes a normal pressure)."""
    assert capstan_normal_pressure(6.0, 1.0 / 3.0) == pytest.approx(2.0)
    assert np.allclose(capstan_normal_pressure(np.array([2.0, 4.0]), np.array([0.5, 0.25])), [1.0, 1.0])


def test_capstan_semicircle_net_force_is_2T() -> None:
    """The classical belt-over-cylinder result: a cap fiber of tension T over a half-wrap presses the
    nucleus with net force 2T along the wrap bisector — ∮T·κ(n̂·d̂)ds reproduces it."""
    T, R = 5.0, 3.0
    theta = np.linspace(0.0, np.pi, 4000)
    ds = np.full_like(theta, R * (theta[1] - theta[0]))
    kappa = np.full_like(theta, 1.0 / R)
    inward = -np.column_stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)])  # toward centre
    net = capstan_line_integral(T, kappa, ds, inward, direction=(0.0, -1.0, 0.0))
    assert net == pytest.approx(2.0 * T, rel=2e-3)


def test_capstan_full_ring_cancels() -> None:
    """A full symmetric wrap → net directional force 0 (inward normals cancel — a hydrostatic squeeze,
    no net translation; the residual pressure is carried by the volume constraint)."""
    T, R = 5.0, 3.0
    theta = np.linspace(0.0, 2 * np.pi, 4000, endpoint=False)
    ds = np.full_like(theta, R * (2 * np.pi / theta.size))
    kappa = np.full_like(theta, 1.0 / R)
    inward = -np.column_stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)])
    net = capstan_line_integral(T, kappa, ds, inward, direction=(0.0, -1.0, 0.0))
    assert abs(net) < 1e-6


# ---- oblate flatten at constant volume (I7 pre-authored) ------------------------------------------
def test_oblate_flatten_conserves_volume() -> None:
    """As the nucleus flattens (aspect↑) at fixed nucleoplasm volume, the equatorial radius grows as
    a ∝ A^{1/3} and the reconstructed volume is invariant."""
    v0 = 4.0 / 3.0 * np.pi * 3.0 ** 3                            # resting sphere R=3
    a1 = oblate_equatorial_radius_at_constant_volume(v0, 1.0)
    assert a1 == pytest.approx(3.0, rel=1e-9)                    # aspect 1 → sphere radius
    for aspect in (1.5, 2.0, 3.0):
        a = oblate_equatorial_radius_at_constant_volume(v0, aspect)
        c = a / aspect
        assert 4.0 / 3.0 * np.pi * a * a * c == pytest.approx(v0, rel=1e-9)
        assert a == pytest.approx(3.0 * aspect ** (1.0 / 3.0), rel=1e-9)


def test_oblate_aspect_from_flatten() -> None:
    """Polar-height ratio → aspect A=h^{-3/2} (0.5 height → A≈2.83), volume conserved throughout."""
    assert oblate_aspect_from_flatten(1.0) == pytest.approx(1.0)
    assert oblate_aspect_from_flatten(0.5) == pytest.approx(2.0 ** 1.5, rel=1e-9)
