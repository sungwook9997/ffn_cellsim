"""2D Mikado fiber network generation.

Implements KU-1.22 (Phase 1 reference parameters) and KU-1.27
(Mikado network construction). Anisotropic orientation distribution
follows KU-1.9 (nematic order parameter S).

All quantities are SI (m, radians, K).

Sanity Gate
-----------
1. Dimensional analysis: L_box [m], L_fiber [m], rest_length = L_fiber
   / (N−1) [m]. No time integration here → no CFL check.
2. Boundary cases:
   - beads_per_fiber < 2 → ValueError.
   - S_order = 0 → uniform angles (von Mises bypass).
   - S_order ≥ 1 − ε → ValueError (perfect alignment singular).
   - Periodic wrap via np.mod ensures positions ∈ [0, L_box) for any
     finite center / orientation.
3. Conservation invariants: stateless — generates positions only.
   Bond lengths are uniform ℓ₀ at construction (within float rounding).
4. Numerical sanity: float64; Bessel ratio uses exponentially scaled
   `ive` to stay numerically stable up to κ ~ 1e3.
5. Sign / sense: nematic angle distribution P(θ) ∝ exp(κ cos 2θ)
   peaks at θ = 0 (mod π), as expected for the head-tail symmetric
   2D nematic.
6. Measurement protocol: `measure_nematic_order` averages cos 2(θ−θ₀)
   over all fibers — the exact estimator of the 2D nematic S
   (KU-1.9); ⟨S_measured⟩ → S_target as n_fibers → ∞ with statistical
   fluctuation σ_S ≈ 1/√(2N).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.special import ive
from scipy.optimize import brentq


@dataclass(slots=True)
class FiberNetwork:
    """Discrete 2D fiber network state.

    Attributes
    ----------
    bead_positions : np.ndarray, shape (n_fibers, n_beads, 2)
        Bead coordinates in meters, wrapped to [0, L_box).
    fiber_centers : np.ndarray, shape (n_fibers, 2)
        Unwrapped fiber centers (drawn uniformly in [0, L_box)²).
    fiber_endpoints : np.ndarray, shape (n_fibers, 2, 2)
        Unwrapped (p1, p2) endpoints of each rod for segment-intersection
        finding. Components may sit slightly outside [0, L_box).
    fiber_orientations : np.ndarray, shape (n_fibers,)
        Per-fiber director angle θ in radians.
    fiber_ids : np.ndarray, shape (n_fibers,)
        Integer identifiers.
    box_size : float
        Periodic box side length L_box in meters.
    rest_length : float
        Bond rest length ℓ₀ = L_fiber / (n_beads - 1) in meters.
    fiber_length : float
        Contour length L_fiber in meters (needed to map an intersection
        parameter s ∈ [-L_f/2, +L_f/2] back to a bead index).
    params : dict
        Snapshot of the input parameters (KU-1.22 defaults).
    """

    bead_positions: np.ndarray
    fiber_centers: np.ndarray
    fiber_endpoints: np.ndarray
    fiber_orientations: np.ndarray
    fiber_ids: np.ndarray
    box_size: float
    rest_length: float
    fiber_length: float
    params: dict[str, Any] = field(default_factory=dict)


def _kappa_for_nematic_order(S: float, tol: float = 1e-9) -> float:
    """Solve I₁(κ)/I₀(κ) = S for the von Mises concentration κ.

    The 2D nematic order parameter for a distribution
    P(θ) ∝ exp(κ cos 2θ) is S = I₁(κ) / I₀(κ).
    Returns κ ≥ 0. For S = 0 returns 0 exactly.
    """
    if S <= tol:
        return 0.0
    if S >= 1.0 - tol:
        raise ValueError("S must be < 1 (perfect alignment is singular).")

    def residual(k: float) -> float:
        # ive(n, k) = iv(n, k) * exp(-k); the exp(-k) factors cancel
        # in the ratio, giving a numerically stable form for large k.
        return ive(1, k) / ive(0, k) - S

    return brentq(residual, 1e-6, 1e3, xtol=tol)


def _sample_orientations(
    n: int, S: float, theta0: float, rng: np.random.Generator
) -> np.ndarray:
    """Sample n fiber director angles from a 2D nematic distribution.

    The distribution is symmetric under θ → θ + π (head-tail) and
    centred on `theta0`. For S = 0 it reduces to a uniform draw
    on [0, π).
    """
    if S <= 0.0:
        return rng.uniform(0.0, np.pi, size=n)

    kappa = _kappa_for_nematic_order(S)
    # Sample φ = 2θ ∈ [-π, π) from von Mises(μ=0, κ), then θ = (φ/2) + θ0.
    phi = rng.vonmises(mu=0.0, kappa=kappa, size=n)
    return np.mod(0.5 * phi + theta0, np.pi)


def generate_2d_fiber_network(
    L_box: float,
    n_fibers: int,
    L_fiber: float,
    beads_per_fiber: int = 5,
    S_order: float = 0.0,
    theta0: float = 0.0,
    seed: int | None = None,
    params: dict[str, Any] | None = None,
) -> FiberNetwork:
    """Generate a 2D Mikado-style fiber network with periodic BC (KU-1.27).

    Each fiber is a straight rod of length `L_fiber` with center drawn
    uniformly in the periodic box and orientation drawn from a 2D
    nematic distribution with order parameter `S_order` about `theta0`
    (KU-1.9). The fiber is discretized into `beads_per_fiber` beads
    with equal spacing ℓ₀ = L_fiber / (beads_per_fiber − 1).

    Parameters
    ----------
    L_box : float
        Box side length (meters).
    n_fibers : int
        Number of fibers.
    L_fiber : float
        Contour length of each fiber (meters).
    beads_per_fiber : int, default 5
        Number of beads per fiber (KU-1.22).
    S_order : float, default 0.0
        Nematic order parameter; 0 = isotropic, ~0.5 = aligned.
    theta0 : float, default 0.0
        Director angle for the aligned case (radians).
    seed : int or None
        RNG seed for reproducibility.
    params : dict or None
        Optional metadata to attach for downstream modules.

    Returns
    -------
    FiberNetwork
    """
    if beads_per_fiber < 2:
        raise ValueError("beads_per_fiber must be ≥ 2.")
    if L_fiber <= 0 or L_box <= 0:
        raise ValueError("Lengths must be positive.")

    rng = np.random.default_rng(seed)

    centers = rng.uniform(0.0, L_box, size=(n_fibers, 2))
    angles = _sample_orientations(n_fibers, S_order, theta0, rng)

    # Discretize each fiber: beads at offsets s ∈ [-L/2, +L/2] along the director.
    s = np.linspace(-0.5 * L_fiber, 0.5 * L_fiber, beads_per_fiber)
    directors = np.stack([np.cos(angles), np.sin(angles)], axis=1)  # (n, 2)
    offsets = s[None, :, None] * directors[:, None, :]                # (n, beads, 2)
    bead_positions_unwrapped = centers[:, None, :] + offsets
    bead_positions = np.mod(bead_positions_unwrapped, L_box)          # periodic wrap

    # Unwrapped segment endpoints (used by segment-intersection finder).
    endpoint_offsets = np.array([-0.5 * L_fiber, +0.5 * L_fiber])      # (2,)
    fiber_endpoints = (
        centers[:, None, :] + endpoint_offsets[None, :, None] * directors[:, None, :]
    )

    rest_length = L_fiber / (beads_per_fiber - 1)

    return FiberNetwork(
        bead_positions=bead_positions.astype(np.float64, copy=False),
        fiber_centers=centers.astype(np.float64, copy=False),
        fiber_endpoints=fiber_endpoints.astype(np.float64, copy=False),
        fiber_orientations=angles.astype(np.float64, copy=False),
        fiber_ids=np.arange(n_fibers, dtype=np.int64),
        box_size=float(L_box),
        rest_length=float(rest_length),
        fiber_length=float(L_fiber),
        params=dict(params or {}),
    )


def measure_nematic_order(angles: np.ndarray, theta0: float = 0.0) -> float:
    """Return the 2D nematic order parameter S = ⟨cos 2(θ − θ₀)⟩ (KU-1.9)."""
    return float(np.mean(np.cos(2.0 * (angles - theta0))))
