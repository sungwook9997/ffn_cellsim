r"""Explicit assembly of a matrix-free tangent operator, and its conservativity check.

WHY THIS EXISTS.  Every inner solve in this engine applies ``K = ∇²U`` as a *matrix-free* sequence of
Warp launches (:meth:`aleph.engine.sf_implicit.SFImplicitCG._stiffness`,
``ac/cell/implicit_mechanics``).  A matrix-free operator can be *used* without ever being *looked at* —
so the project has run PCG against it for months while knowing none of its spectral properties: not its
largest eigenvalue (the CFL bound is a Gershgorin OVER-estimate, never the real λ_max), not its null
space, not its condition number, and not whether it is symmetric at all.

Column probing turns the operator into an explicit matrix without touching the physics: applying it to
each Cartesian unit vector ``e_j`` returns column ``j`` exactly, because the operator is linear in its
input vector by construction (it is a tangent).  That costs ``3·n_nodes`` applications, so it is a
SLICE-SCALE instrument — at 494,802 nodes it would need 1.5M applications and is not the right tool.
The honest statement of its scope lives in :func:`assemble_dense_operator`'s ``Raises``.

THE SYMMETRY RESIDUAL IS A PHYSICS TEST, NOT A NUMERICS ONE.  A force field is conservative — i.e. it
is ``−∇U`` for SOME potential ``U`` — if and only if its Jacobian is symmetric.  So ``‖K − Kᵀ‖``
measures, with no ``U`` in hand and no energy function written anywhere, whether the slice's assembled
forces admit an energy at all.  A term that is double-counted asymmetrically, a scatter that violates
Newton's third law, or a tangent that does not match its own force all show up here as a nonzero
antisymmetric part.  This is the magnitude-free double-count guard the execution plan asks for
(§2 consequence 3), available before any energy ledger exists.

Its floor is DERIVED, never chosen: the same Higham summation bound the accepted-step balance gate uses
(:func:`aleph.engine.ledger.assemble_balance_tolerance_ratio`, PI decision D8), scaled by the
operator's own largest entry.  Reusing that function rather than writing a second tolerance is
deliberate — two independently-derived floors would drift.

engine units: length µm, force pN, stiffness pN/µm.

Sanity Gate:
    * dimensional: a stiffness column is ``[pN/µm]`` because the probe vector is a unit displacement
      ``[µm]`` and the operator returns a force ``[pN]``.  The symmetry residual is reported both
      absolutely ``[pN/µm]`` and as a dimensionless ratio.
    * boundary: a zero operator has zero norm; the relative asymmetry of a zero operator is defined as
      ``0.0`` rather than ``nan``, and is reported alongside the absolute norm so the caller can see it
      is a degenerate case.
    * conservation/invariant: symmetry of ``K`` is exactly the integrability condition for ``F = −∇U``.
      An asymmetric ``K`` means no potential exists, which is reported, never repaired.
    * numerical: the acceptance floor is the float64 order-independent accumulation bound for the
      declared contribution count — atomics make the summation order nondeterministic, so no
      order-dependent (and no bit-identity) criterion may be used.
    * sign-sense: not applicable; the probe computes no force of its own.
    * measurement-protocol: the probe records the exact ``regularization`` it was applied at.  A
      regularized operator ``aI + K`` has every eigenvalue shifted by ``a``, so a spectrum measured at
      ``a > 0`` is NOT the spectrum of ``K`` and must not be reported as one.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from aleph.engine.ledger import assemble_balance_tolerance_ratio

__all__ = [
    "NodeMatvec",
    "SymmetryReport",
    "assemble_dense_operator",
    "symmetry_report",
]

#: A matrix-free tangent application: ``(n, 3) -> (n, 3)``, linear in its argument.
#:
#: The driver supplies a closure that uploads the argument, launches the real device kernels, and
#: downloads the result.  Keeping the probe's contract at this boundary is what lets it be exercised
#: against a pure-NumPy reference operator with no device at all.
NodeMatvec = Callable[[np.ndarray], np.ndarray]


def assemble_dense_operator(
    matvec: NodeMatvec,
    n_nodes: int,
    *,
    max_dof: int = 8192,
) -> np.ndarray:
    """Return the explicit ``(3n, 3n)`` matrix of a matrix-free node operator by column probing.

    Column ``3i + c`` is the operator applied to the unit displacement of node ``i`` along axis ``c``.
    That is exact — not an approximation — for any operator that is linear in its input vector, which
    every tangent in this engine is.

    Args:
        matvec: The matrix-free application, ``(n, 3) -> (n, 3)``.
        n_nodes: Node count ``n``; the returned matrix is ``(3n, 3n)``.
        max_dof: Refusal threshold on ``3n``.  Column probing is ``O(n)`` operator applications and
            ``O(n²)`` memory, so it is a slice instrument; a full-population request is a wiring
            mistake and is rejected loudly rather than started and left to exhaust memory hours later.

    Returns:
        The dense operator, float64, shape ``(3n, 3n)``.

    Raises:
        ValueError: If ``n_nodes`` is not positive, if ``3·n_nodes`` exceeds ``max_dof``, or if the
            supplied ``matvec`` returns something other than a finite ``(n, 3)`` array.
    """
    n = int(n_nodes)
    if n <= 0:
        raise ValueError(f"n_nodes must be positive; got {n_nodes!r}")
    n_dof = 3 * n
    if n_dof > int(max_dof):
        raise ValueError(
            f"column probing needs {n_dof} operator applications and a {n_dof}x{n_dof} float64 matrix "
            f"({n_dof * n_dof * 8 / 1e9:.2f} GB); it is a SLICE instrument and this population exceeds "
            f"max_dof={max_dof}. Use an iterative (Lanczos) spectrum at full population instead"
        )

    matrix = np.zeros((n_dof, n_dof), np.float64)
    probe = np.zeros((n, 3), np.float64)
    for column in range(n_dof):
        node, axis = divmod(column, 3)
        probe[node, axis] = 1.0
        response = np.asarray(matvec(probe), np.float64)
        probe[node, axis] = 0.0
        if response.shape != (n, 3):
            raise ValueError(
                f"matvec returned shape {response.shape}, expected {(n, 3)} at column {column}"
            )
        if not np.all(np.isfinite(response)):
            raise ValueError(f"matvec returned a non-finite response at column {column}")
        matrix[:, column] = response.reshape(-1)
    return matrix


@dataclass(frozen=True, slots=True)
class SymmetryReport:
    """Whether an assembled tangent admits a potential, with a derived acceptance floor.

    Attributes:
        n_dof: Width of the operator.
        max_abs_entry: ``max|K_ij|`` [pN/µm] — the magnitude the floor is scaled by.
        frobenius_norm: ``‖K‖_F`` [pN/µm].
        asymmetry_frobenius: ``‖K − Kᵀ‖_F`` [pN/µm].
        max_abs_asymmetry: ``max|K_ij − K_ji|`` [pN/µm] — the entrywise worst case, which a Frobenius
            norm can hide when a single pair is badly wrong among many correct ones.
        relative_asymmetry: ``‖K − Kᵀ‖_F / ‖K‖_F``, dimensionless; ``0.0`` for a zero operator.
        round_off_floor: ``γ_n · max|K_ij|`` [pN/µm] — the largest entrywise asymmetry attributable to
            float64 accumulation alone.
        contribution_count: The ``n_terms`` the floor was derived at, recorded so the floor can be
            re-derived from the artifact.
        conservative: Whether ``max_abs_asymmetry`` sits at or below ``round_off_floor``, i.e. whether
            a potential ``U`` exists to within accumulation error.
    """

    n_dof: int
    max_abs_entry: float
    frobenius_norm: float
    asymmetry_frobenius: float
    max_abs_asymmetry: float
    relative_asymmetry: float
    round_off_floor: float
    contribution_count: int
    conservative: bool


def symmetry_report(matrix: np.ndarray, *, contribution_count: int) -> SymmetryReport:
    """Measure the antisymmetric part of an assembled tangent against a derived round-off floor.

    Args:
        matrix: The assembled operator, square and float64-convertible.
        contribution_count: How many accumulated terms land in the worst-case matrix entry — the
            ``n_terms`` of :func:`~aleph.engine.ledger.assemble_balance_tolerance_ratio`.  A COUNT,
            never a magnitude: a floor that took a force scale would relocate the magic number rather
            than remove it (PI D8).  Callers pass the operator's total launched-contribution count,
            which is an over-estimate of the per-entry count and therefore a conservative floor.

    Returns:
        The :class:`SymmetryReport`.

    Raises:
        ValueError: If ``matrix`` is not square, is not finite, or ``contribution_count`` is below 1.
    """
    k = np.asarray(matrix, np.float64)
    if k.ndim != 2 or k.shape[0] != k.shape[1]:
        raise ValueError(f"operator must be square; got shape {k.shape}")
    if not np.all(np.isfinite(k)):
        raise ValueError("operator contains non-finite entries")

    skew = k - k.T
    max_abs_entry = float(np.max(np.abs(k))) if k.size else 0.0
    frobenius = float(np.linalg.norm(k, "fro"))
    asym_frobenius = float(np.linalg.norm(skew, "fro"))
    max_abs_asym = float(np.max(np.abs(skew))) if skew.size else 0.0
    floor = assemble_balance_tolerance_ratio(int(contribution_count)) * max_abs_entry
    return SymmetryReport(
        n_dof=int(k.shape[0]),
        max_abs_entry=max_abs_entry,
        frobenius_norm=frobenius,
        asymmetry_frobenius=asym_frobenius,
        max_abs_asymmetry=max_abs_asym,
        relative_asymmetry=(asym_frobenius / frobenius) if frobenius > 0.0 else 0.0,
        round_off_floor=float(floor),
        contribution_count=int(contribution_count),
        conservative=bool(max_abs_asym <= floor),
    )
