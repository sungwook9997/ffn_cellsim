"""Self-test of the diffusion Green's-function oracle (pure NumPy — no Warp/CUDA)."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.fluid.greens_analytic import (
    continuous_point_source_3d,
    diffusion_greens_function,
    mean_square_radius,
)

C_V = 50.0  # um^2/s (KB-3.B3.2 anchor)


def test_integral_conserved_1d() -> None:
    """Space integral of the 1-D kernel equals the injected content (mass conservation)."""
    t = 0.5
    sigma = np.sqrt(2.0 * C_V * t)
    x = np.linspace(-8.0 * sigma, 8.0 * sigma, 20_001)
    g = diffusion_greens_function(np.abs(x), t, C_V, dim=1, source_strength=1.0)
    assert np.trapezoid(g, x) == pytest.approx(1.0, abs=1e-4)


def test_integral_conserved_3d() -> None:
    """Radial integral of the 3-D kernel (4 pi r^2 weight) equals the injected content."""
    t = 0.5
    sigma = np.sqrt(2.0 * C_V * t)
    r = np.linspace(0.0, 10.0 * sigma, 40_001)
    g = diffusion_greens_function(r, t, C_V, dim=3, source_strength=1.0)
    integral = np.trapezoid(g * 4.0 * np.pi * r**2, r)
    assert integral == pytest.approx(1.0, abs=1e-4)


def test_mean_square_radius_matches_diffusion_law() -> None:
    """Numerical <r^2> of the 3-D kernel equals the exact 2 d c_v t = 6 c_v t."""
    t = 0.5
    sigma = np.sqrt(2.0 * C_V * t)
    r = np.linspace(0.0, 12.0 * sigma, 60_001)
    g = diffusion_greens_function(r, t, C_V, dim=3)
    weight = g * 4.0 * np.pi * r**2
    r2_numeric = np.trapezoid(r**2 * weight, r) / np.trapezoid(weight, r)
    assert r2_numeric == pytest.approx(mean_square_radius(t, C_V, 3), rel=1e-3)
    assert mean_square_radius(t, C_V, 3) == pytest.approx(6.0 * C_V * t)


def test_kernel_satisfies_diffusion_pde_1d() -> None:
    """G_1d satisfies dp/dt = c_v d2p/dx2 (finite-difference check away from the peak)."""
    t, dt, dx = 0.5, 1e-4, 0.02
    x = np.linspace(2.0, 12.0, 400)  # off-peak so the Gaussian is smooth on this dx
    g0 = diffusion_greens_function(np.abs(x), t, C_V, 1)
    gp = diffusion_greens_function(np.abs(x), t + dt, C_V, 1)
    gm = diffusion_greens_function(np.abs(x), t - dt, C_V, 1)
    dgdt = (gp - gm) / (2.0 * dt)
    xw = np.array([x - dx, x, x + dx])
    gw = diffusion_greens_function(np.abs(xw), t, C_V, 1)
    d2gdx2 = (gw[0] - 2.0 * gw[1] + gw[2]) / dx**2
    # compare where the field is non-negligible (relative residual small)
    mask = g0 > 1e-4 * g0.max()
    resid = np.abs(dgdt[mask] - C_V * d2gdx2[mask]) / (np.abs(dgdt[mask]) + 1e-30)
    assert np.max(resid) < 1e-3


def test_continuous_source_approaches_steady_poisson() -> None:
    """Continuous 3-D source -> Q/(4 pi c_v r) as t -> inf (erfc -> 1)."""
    r = np.array([1.0, 2.0, 5.0])
    steady = 1.0 / (4.0 * np.pi * C_V * r)
    late = continuous_point_source_3d(r, t=1e8, c_v=C_V, rate=1.0)
    assert np.allclose(late, steady, rtol=1e-4)
    # and it is strictly below the steady value at finite time
    early = continuous_point_source_3d(r, t=0.01, c_v=C_V, rate=1.0)
    assert np.all(early < steady)
