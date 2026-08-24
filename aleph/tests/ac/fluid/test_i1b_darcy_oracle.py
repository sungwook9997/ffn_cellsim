"""I1b Darcy / pore-fluid-velocity gates — the fluid genuinely FLOWS (q, v_f), verified analytically.

Pure NumPy (no Warp/CUDA). Gates the kinematics of the mixed q-p form: Darcy flux linear in Delta p and
in mobility; v_f - v_s = q/phi exact for any porosity (phi swept — I0-B1b GAP); hydrostatic column carries
no discharge; steady point-source radial conservation; and no relative flow at drained equilibrium (the
bridge to the I1a solver's settled state).

Uses the I0-B1 closure c_v = mobility/S = 50 um^2/s. phi is UNSOURCED (I0-B1b GAP -> PI) and appears only
in v_f; the analytic gate is phi-agnostic.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.fluid.darcy_analytic import (
    darcy_flux,
    darcy_slab_flux,
    hydrostatic_gradient,
    pore_fluid_velocity,
    steady_point_source_flux,
)

C_V = 50.0          # um^2/s
MOBILITY = 5.0e-3   # k/mu, um^2/(Pa*s)
STORAGE_S = 1.0e-4  # 1/Pa


def test_slab_flux_linear_in_delta_p() -> None:
    """Doubling the pressure drop doubles the Darcy discharge (linearity in Delta p)."""
    q1 = darcy_slab_flux(10.0, length=5.0, mobility=MOBILITY)
    q2 = darcy_slab_flux(20.0, length=5.0, mobility=MOBILITY)
    assert q2 == pytest.approx(2.0 * q1)
    assert q1 == pytest.approx(MOBILITY * 10.0 / 5.0)


def test_slab_flux_linear_in_mobility() -> None:
    """Discharge is linear in mobility k/mu (permeability), at fixed pressure drop."""
    q1 = darcy_slab_flux(10.0, length=5.0, mobility=MOBILITY)
    q2 = darcy_slab_flux(10.0, length=5.0, mobility=3.0 * MOBILITY)
    assert q2 == pytest.approx(3.0 * q1)


def test_uniform_gradient_slab_has_uniform_flux() -> None:
    """A linear pressure profile gives spatially-uniform q (div q = 0 -> steady state, no source)."""
    x = np.linspace(0.0, 5.0, 51)
    p = 40.0 - 2.0 * x  # constant gradient dp/dx = -2 Pa/um
    grad_p = np.gradient(p, x)
    q = darcy_flux(grad_p, MOBILITY)
    assert np.allclose(q, q[0], atol=1e-9)
    assert q[0] == pytest.approx(-MOBILITY * (-2.0))


@pytest.mark.parametrize("phi", [0.5, 0.65, 0.7, 0.8, 0.99])
def test_pore_velocity_reconstruction_exact(phi: float) -> None:
    """v_f - v_s == q/phi to round-off, for a sweep of the (GAP) porosity."""
    rng = np.random.default_rng(0)
    v_s = rng.standard_normal(3)
    q = rng.standard_normal(3) * 1e-3
    v_f = pore_fluid_velocity(v_s, q, phi)
    assert np.allclose(v_f - v_s, q / phi, atol=1e-15)


def test_hydrostatic_column_carries_no_discharge() -> None:
    """grad p == rho_f b (hydrostatic equilibrium) -> q == 0 (no gravity-driven seepage)."""
    rho_f = 1.0e3            # kg/m^3 (water)
    g = 9.81                 # m/s^2
    grad_hydro = hydrostatic_gradient(rho_f, g)
    q = darcy_flux(grad_hydro, MOBILITY, rho_f_b=grad_hydro)
    assert q == pytest.approx(0.0, abs=1e-18)


def test_gravity_shifts_flux_by_mobility_rho_b() -> None:
    """With no pressure gradient, gravity drives q = +mobility * rho_f b (buoyant seepage)."""
    rho_f_b = 3.0
    q = darcy_flux(0.0, MOBILITY, rho_f_b=rho_f_b)
    assert q == pytest.approx(MOBILITY * rho_f_b)


def test_point_source_radial_flux_conserved() -> None:
    """Steady 3-D point source: 4 pi r^2 q(r) is independent of r (radial conservation = S*rate)."""
    r = np.array([1.0, 2.5, 5.0, 9.0])
    rate = 2.0
    q = steady_point_source_flux(r, rate, C_V, MOBILITY)
    shell_flux = 4.0 * np.pi * r**2 * q
    assert np.allclose(shell_flux, shell_flux[0], rtol=1e-12)
    assert shell_flux[0] == pytest.approx(STORAGE_S * rate * (MOBILITY / (C_V * STORAGE_S)))  # = S*rate


def test_no_relative_flow_at_equilibrium() -> None:
    """A uniform (drained-equilibrium) pressure field carries no discharge -> v_f == v_s (no streaming)."""
    grad_p = np.zeros(3)
    q = darcy_flux(grad_p, MOBILITY)
    v_s = np.array([0.1, -0.2, 0.05])
    v_f = pore_fluid_velocity(v_s, q, phi=0.7)
    assert np.allclose(v_f, v_s, atol=1e-15)
