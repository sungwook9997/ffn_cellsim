"""Acceptance algebra for the Active Cell analytic projected implicit operator.

The runtime solves ``(a I + P K P) dx = P F`` on Warp CUDA.  ``K`` is assembled matrix-free from the exact
NF2007 bending Hessian and current central-spring tangents; ``P`` is the existing inextensibility projector.
This module provides small dense NumPy oracles only.

Sanity Gate:
    * Dimensions: ``K`` and ``a`` are pN/um, ``F`` is pN, and ``dx`` is um.
    * Symmetry: bending and central-spring tangents are symmetric; ``P K P`` is symmetric for orthogonal P.
    * Positivity: tension-side transverse tangents are nonnegative; compression-side transverse curvature is
      omitted from the SPD solve rather than smuggled in as a negative stiffness.
    * Boundary cases: ``a>0`` removes rigid singularities; zero force gives exactly zero displacement.
    * Fixed point: multiplying the force by an SPD inverse changes the convergence path, not the root ``P F=0``.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "bending_stiffness",
    "central_spring_tangent",
    "projected_implicit_solve",
    "rigid_strain_coarse_basis",
    "two_level_coarse_correction",
]


def bending_stiffness(
    n_nodes: int, triples: npt.ArrayLike, alpha: npt.ArrayLike,
) -> npt.NDArray[np.float64]:
    """Return the dense ``3N x 3N`` NF2007 bending stiffness."""
    tri = np.asarray(triples, dtype=np.int64)
    coeff = np.asarray(alpha, dtype=np.float64)
    if n_nodes < 0 or tri.ndim != 2 or tri.shape[1:] != (3,) or coeff.shape != (tri.shape[0],):
        raise ValueError("invalid bending topology")
    if np.any(tri < 0) or np.any(tri >= n_nodes):
        raise ValueError("triple index out of range")
    if np.any(~np.isfinite(coeff)) or np.any(coeff < 0.0):
        raise ValueError("alpha must be finite and nonnegative")
    scalar = np.zeros((n_nodes, n_nodes), dtype=np.float64)
    stencil = np.array([1.0, -2.0, 1.0])
    for nodes, value in zip(tri, coeff, strict=True):
        scalar[np.ix_(nodes, nodes)] += value * np.outer(stencil, stencil)
    return np.kron(scalar, np.eye(3, dtype=np.float64))


def central_spring_tangent(
    delta: npt.ArrayLike, *, stiffness: float, rest_length: float, tension_only: bool = False,
) -> npt.NDArray[np.float64]:
    """Return the SPD part of a current central Hookean spring tangent."""
    vector = np.asarray(delta, dtype=np.float64)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError("delta must be a finite three-vector")
    if not np.isfinite(stiffness) or stiffness < 0.0:
        raise ValueError("stiffness must be finite and nonnegative")
    if not np.isfinite(rest_length) or rest_length < 0.0:
        raise ValueError("rest_length must be finite and nonnegative")
    length = float(np.linalg.norm(vector))
    extension = length - rest_length
    if length == 0.0 or (tension_only and extension <= 0.0):
        return np.zeros((3, 3), dtype=np.float64)
    unit = vector / length
    axial = stiffness * np.outer(unit, unit)
    transverse = max(stiffness * extension / length, 0.0) * (
        np.eye(3, dtype=np.float64) - np.outer(unit, unit))
    return np.asarray(axial + transverse, dtype=np.float64)


def projected_implicit_solve(
    force: npt.ArrayLike,
    stiffness: npt.ArrayLike,
    projector: npt.ArrayLike,
    *,
    regularization: float,
) -> npt.NDArray[np.float64]:
    """Solve the dense acceptance system ``(aI + PKP) dx = PF``."""
    f = np.asarray(force, dtype=np.float64).reshape(-1)
    k = np.asarray(stiffness, dtype=np.float64)
    p = np.asarray(projector, dtype=np.float64)
    if k.shape != (f.size, f.size) or p.shape != k.shape:
        raise ValueError("stiffness and projector must match force size")
    if not np.isfinite(regularization) or regularization <= 0.0:
        raise ValueError("regularization must be finite and positive")
    operator = regularization * np.eye(f.size) + p @ k @ p
    return np.linalg.solve(operator, p @ f)


def rigid_strain_coarse_basis(
    positions: npt.ArrayLike, *, n_modes: int = 12,
) -> npt.NDArray[np.float64]:
    r"""Return an orthonormal ``(K, N, 3)`` rigid-body + constant-strain coarse space.

    The slow modes that neither the diagonal preconditioner nor the per-fiber blocks can see are the *global*
    smooth displacement fields: a whole-cell breathing (l=0), rigid translation/rotation (l=1), and ellipsoidal
    shape change (l=2). Expressed in the natural analytic basis these are exactly the classical near-kernel of a
    3-D elastic operator — the 3 rigid translations, 3 rigid rotations, and 6 constant symmetric strains
    ``u(x) = E (x - c)`` — and ``span{translation, isotropic strain (l=0), deviatoric strain (l=2), rotation}``
    is precisely the ``l <= 2`` spherical space (1 + 3 + 3 + 5 = 12 modes). The basis is derived from geometry
    alone (no relaxation parameter, no tuned constant); orthonormalizing the columns only conditions the coarse
    Galerkin matrix and never changes the span. Deflating this space is what removes the global-radial residual
    component that per-fiber and node-local preconditioning leave behind.

    Args:
        positions: ``(N, 3)`` reference nodal coordinates [um].
        n_modes: number of leading modes to keep (12 = full ``l <= 2``; ordered translations, rotations,
            strains). A degenerate (e.g. collinear) cloud that cannot support ``n_modes`` independent columns
            raises rather than silently returning a rank-deficient space.

    Returns:
        ``(K, N, 3)`` array whose ``K`` node-fields are L2-orthonormal over the flattened ``3N`` displacement
        space (``sum_i basis[j, i] . basis[k, i] == delta_jk``).
    """
    pos = np.asarray(positions, dtype=np.float64)
    if pos.ndim != 2 or pos.shape[1] != 3 or not np.all(np.isfinite(pos)):
        raise ValueError("positions must be a finite (N, 3) array")
    n = pos.shape[0]
    if not 1 <= n_modes <= 12:
        raise ValueError("n_modes must be between 1 and 12")
    if n == 0:
        raise ValueError("positions must contain at least one node")
    relative = pos - pos.mean(axis=0)

    columns: list[npt.NDArray[np.float64]] = []
    identity = np.eye(3, dtype=np.float64)
    # l = 1 solid translations.
    for axis in range(3):
        field = np.zeros((n, 3), dtype=np.float64)
        field[:, axis] = 1.0
        columns.append(field)
    # l = 1 tangential rotations u = e_axis x (x - c).
    for axis in range(3):
        columns.append(np.cross(identity[axis], relative))
    # l = 0 isotropic + l = 2 deviatoric constant strains u = E (x - c) over the 6 symmetric E.
    strain_pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    for a, b in strain_pairs:
        strain = np.zeros((3, 3), dtype=np.float64)
        strain[a, b] = 1.0
        strain[b, a] = 1.0
        columns.append(relative @ strain)

    matrix = np.stack(columns[:n_modes], axis=0).reshape(n_modes, 3 * n).T  # (3N, K)
    q, r = np.linalg.qr(matrix)
    rank = int(np.sum(np.abs(np.diag(r)) > 1.0e-12 * max(1.0, float(np.abs(r).max()))))
    if rank < n_modes:
        raise ValueError(
            f"geometry supports only {rank} independent coarse modes, {n_modes} requested")
    # Fix the QR sign gauge so the basis is deterministic across platforms.
    q = q * np.sign(np.diag(r))
    return np.ascontiguousarray(q.T.reshape(n_modes, n, 3))


def two_level_coarse_correction(
    operator: npt.ArrayLike, basis: npt.ArrayLike, residual: npt.ArrayLike,
) -> npt.NDArray[np.float64]:
    r"""Return the exact additive Galerkin coarse correction ``B (B^T A B)^-1 B^T r``.

    This is the acceptance reference for the device two-level preconditioner. ``A`` is the SPD projected
    operator ``a I + P K P``; ``B`` are the orthonormal coarse columns from :func:`rigid_strain_coarse_basis`.
    Adding ``B A_c^{-1} B^T`` (rank ``K``, SPSD) to the existing SPD fine preconditioner keeps the composite SPD,
    so PCG stays valid — the correction only accelerates the global smooth modes, it cannot move the fixed point
    or relax any residual gate.

    Args:
        operator: dense ``(3N, 3N)`` SPD operator ``A``.
        basis: ``(K, N, 3)`` orthonormal coarse space.
        residual: ``(N, 3)`` (or flat ``3N``) residual to precondition.

    Returns:
        ``(N, 3)`` coarse correction ``z_c`` with ``A_c z_c`` coefficients solved exactly.
    """
    a = np.asarray(operator, dtype=np.float64)
    b = np.asarray(basis, dtype=np.float64)
    r = np.asarray(residual, dtype=np.float64).reshape(-1)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("operator must be square")
    if b.ndim != 3 or b.shape[2] != 3 or b.shape[1] * 3 != a.shape[0]:
        raise ValueError("basis must be (K, N, 3) matching the operator size")
    if r.size != a.shape[0]:
        raise ValueError("residual must match the operator size")
    flat = b.reshape(b.shape[0], -1).T  # (3N, K)
    coarse_operator = flat.T @ a @ flat
    coefficients = np.linalg.solve(coarse_operator, flat.T @ r)
    return (flat @ coefficients).reshape(-1, 3)
