"""Darcy discharge + pore-fluid velocity — the closed-form ground truth for the I1b mixed q-p form.

Host-side acceptance oracles (pure NumPy, NO Warp/simulation runtime). I1a solves pressure/mass; I1b
makes the fluid genuinely FLOW by reconstructing the Darcy discharge and the absolute pore-fluid velocity
(NEW_ENGINE_BUILD_PLAN §1 P1, §3 I1b):

    q   = phi (v_f - v_s) = -(k/mu) (grad p - rho_f b)      relative discharge (per unit bulk area)
    v_f = v_s + q / phi                                     absolute pore-fluid velocity

The payload of the flow is streaming / advection / rate-dependence (I1c rides ``v_f``), NOT dominating the
static balance — at drained steady state the pressure deviation is small (~Pa). These oracles gate the
KINEMATICS: Darcy flux is LINEAR in the pressure drop and in the mobility ``k/mu``; the velocity
reconstruction ``v_f - v_s = q/phi`` is exact for ANY porosity ``phi`` (so the analytic gate is
phi-agnostic — ``phi`` is an I0-B1b GAP, swept as an oracle parameter like ``alpha`` was in I1a); a
hydrostatic column (``grad p = rho_f b``) carries ZERO discharge; and for a steady point source the flux
through any enclosing sphere equals the source rate (radial conservation), the ``q``-form shadow of the
Green oracle.

Sanity Gate (self-tested in tests/ac/fluid/test_i1b_darcy_oracle.py):
  * q linear in Delta p (double the drop -> double the flux) and in mobility;
  * a uniform-gradient slab has spatially-uniform q (div q = 0, no source) -> steady state;
  * v_f - v_s == q/phi to round-off, for a sweep of phi;
  * grad p == rho_f b (hydrostatic) -> q == 0;
  * steady 3-D point source: 4 pi r^2 q(r) == rate for all r (radial conservation).

References: Darcy 1856; Biot 1941; Bear 1972, Dynamics of Fluids in Porous Media. Mobility/permeability
anchor: params_i0b1.yaml (k, mu_pore, c_v=k*M/mu, KB-3.B3.2). Porosity phi: I0-B1b (GAP -> PI).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "darcy_flux",
    "darcy_slab_flux",
    "pore_fluid_velocity",
    "hydrostatic_gradient",
    "steady_point_source_flux",
]


def darcy_flux(
    grad_p: npt.ArrayLike, mobility: float, rho_f_b: npt.ArrayLike = 0.0
) -> npt.NDArray[np.float64]:
    """Darcy discharge ``q = -mobility (grad p - rho_f b)`` (relative to the solid, per unit bulk area).

    Args:
        grad_p: Pressure gradient component(s) ``dp/dx_i`` (Pa/length).
        mobility: Darcy mobility ``k/mu`` (length^2 / (Pa*time)).
        rho_f_b: Fluid body-force density component(s) ``rho_f * b`` (Pa/length); 0 = no gravity.

    Returns:
        Discharge ``q`` (length/time), same shape as ``grad_p``.
    """
    if mobility <= 0.0:
        raise ValueError("mobility k/mu must be positive")
    return -mobility * (np.asarray(grad_p, dtype=np.float64) - np.asarray(rho_f_b, dtype=np.float64))


def darcy_slab_flux(delta_p: float, length: float, mobility: float) -> float:
    """Steady 1-D slab discharge ``q = mobility * Delta p / L`` (uniform, div q = 0).

    Args:
        delta_p: Pressure drop across the slab ``p(0) - p(L)`` (Pa).
        length: Slab thickness ``L`` (length).
        mobility: Darcy mobility ``k/mu``.

    Returns:
        The (spatially uniform) discharge magnitude through the slab.
    """
    if length <= 0.0:
        raise ValueError("length must be positive")
    return mobility * delta_p / length


def pore_fluid_velocity(
    v_s: npt.ArrayLike, q: npt.ArrayLike, phi: float
) -> npt.NDArray[np.float64]:
    """Absolute pore-fluid velocity ``v_f = v_s + q/phi``.

    Args:
        v_s: Solid-skeleton velocity component(s) (length/time).
        q: Darcy discharge component(s) (length/time).
        phi: Porosity (fluid volume fraction), in (0, 1]. I0-B1b GAP — swept as an oracle parameter.

    Returns:
        Pore-fluid velocity ``v_f``, same shape as the inputs.
    """
    if not (0.0 < phi <= 1.0):
        raise ValueError("phi must be in (0, 1]")
    return np.asarray(v_s, dtype=np.float64) + np.asarray(q, dtype=np.float64) / phi


def hydrostatic_gradient(rho_f: float, gravity: float) -> float:
    """Hydrostatic pressure gradient ``rho_f * g`` at which the Darcy discharge vanishes.

    Args:
        rho_f: Pore-fluid mass density (mass/length^3).
        gravity: Gravitational acceleration magnitude (length/time^2).

    Returns:
        ``rho_f * g`` (Pa/length). Feeding this as both ``grad_p`` and ``rho_f_b`` gives ``q = 0``.
    """
    return rho_f * gravity


def steady_point_source_flux(r: npt.ArrayLike, rate: float, c_v: float, mobility: float) -> npt.NDArray[np.float64]:
    """Radial discharge of the steady 3-D point source, ``q(r) = rate / (4 pi r^2 * (mobility/... ))``.

    For the steady Poisson limit ``p(r) = rate/(4 pi c_v r)`` (``greens_analytic.continuous_point_source_3d``),
    Darcy gives ``q(r) = -mobility dp/dr = mobility * rate / (4 pi c_v r^2)``. Since ``c_v = mobility/S`` and
    the injected content rate maps through storage, the flux through any sphere ``4 pi r^2 q(r) =
    mobility*rate/c_v = S*rate`` is INDEPENDENT of ``r`` — radial conservation.

    Args:
        r: Radius/radii from the source (> 0).
        rate: Source content rate (content/time).
        c_v: Consolidation coefficient (= mobility/S), length^2/time.
        mobility: Darcy mobility ``k/mu``.

    Returns:
        Radial discharge ``q(r)`` (length/time).
    """
    rr = np.asarray(r, dtype=np.float64)
    if np.any(rr <= 0.0):
        raise ValueError("r must be positive")
    return mobility * rate / (4.0 * np.pi * c_v * rr**2)
