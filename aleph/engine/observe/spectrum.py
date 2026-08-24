r"""The eigen-spectrum of ``K = ∇²U`` as a FIELD observable, plus what it decides.

WHY THE SPECTRUM AND NOT A RESIDUAL.  The engine's convergence verdicts are currently forces compared
against a force floor, and a force residual cannot distinguish "the solver is grinding down a genuinely
slow mode" from "the solver is stuck at a metastable barrier" — the two produce the same plateau.  The
spectrum separates them, and it does so with no threshold and no fitted constant:

* ``λ_max`` is the TRUE largest eigenvalue.  Every explicit step in this repo is sized by a Gershgorin
  bound instead, which is an over-estimate by an unmeasured factor — so the CFL step has always been
  conservative by an unknown amount, and ``2/λ_max`` is the real stability limit.
* the near-null space is the rigid-body / floppy content.  Its DIMENSION is a structural fact about the
  slice (an unbound minifilament is a free rigid body; an unanchored bundle has six zero modes), so
  counting it against the structurally expected count is a build check that no force gate performs.
* ``λ_max/λ_min`` is the time-scale separation, i.e. how many decades of dynamics an explicit integrator
  is forced to resolve to follow the slowest one.  That number IS the case for an implicit solve, stated
  without appeal to anyone's experience.
* the eigenvectors localize or they do not.  A mode occupying few nodes is a local rearrangement; a mode
  spanning the slice is collective.  Measuring the spatial extent of the softest modes is the direct
  measurement of the stress decay length, and therefore of whether the network is a "tensor network" in
  the sense the design conversation asked about.

The eigenvalue array is kept and reported IN FULL, not summarised — it is a field observable per the
plan's §3 data-representation rule, and the eigenvectors are this project's POD basis.

THE ZERO TOLERANCE IS DERIVED.  A symmetric eigensolve is backward stable: the computed eigenvalues are
the exact eigenvalues of ``K + E`` with ``‖E‖₂ ≤ p(n)·ε·‖K‖₂`` (Demmel, *Applied Numerical Linear
Algebra*, §5.2), so no eigenvalue is resolvable below that.  We take ``p(n) = n``, the conventional
conservative surrogate for the modest polynomial.  Nothing here is chosen to make a mode count come out
right — the floor is a property of float64 and the operator norm, and it is reported so a reader can
re-derive it.

engine units: length µm, force pN, stiffness pN/µm; a mobility step is µm/pN.

Sanity Gate:
    * dimensional: eigenvalues of a stiffness are ``[pN/µm]``.  ``2/λ_max`` is ``[µm/pN]``, which is the
      mobility step ``dt/γ`` the explicit path already uses — the two are directly comparable.
    * boundary: an operator with no eigenvalue above the round-off floor yields
      ``lambda_min_nonzero = None`` and a ``None`` condition number rather than a division by zero.
    * conservation/invariant: the count of modes below the floor is compared against a caller-declared
      structurally expected count; a mismatch is REPORTED, never absorbed by widening the floor.
    * numerical: eigenvalues come from ``eigh`` on the SYMMETRISED operator ``(K + Kᵀ)/2``, and the
      asymmetry that symmetrisation discarded must be reported separately
      (:func:`~aleph.engine.observe.operator_probe.symmetry_report`) — silently symmetrising an
      operator that is NOT symmetric would hide exactly the defect that check exists to find.
    * sign-sense: a positive-definite ``K`` has all eigenvalues positive; a negative eigenvalue means the
      configuration is not a minimum of ``U`` (a saddle), which is reported as a count, not clipped.
    * measurement-protocol: the spectrum must be measured at zero regularization.  ``aI + K`` shifts
      every eigenvalue by ``a``, so a regularized probe reports ``λ + a`` and would fabricate a null
      space of size zero.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

__all__ = [
    "LanczosSpectrum",
    "ModeLocalization",
    "StiffnessSpectrum",
    "lanczos_extremal_spectrum",
    "lanczos_ritz",
    "mode_localization",
    "stiffness_spectrum",
]

#: float64 unit round-off, used for the backward-stability floor.  A property of the format, not a knob.
_EPS64: float = float(np.finfo(np.float64).eps)


@dataclass(frozen=True, slots=True)
class StiffnessSpectrum:
    """The full eigen-spectrum of one assembled tangent, with the quantities it decides.

    Attributes:
        eigenvalues: All ``3n`` eigenvalues [pN/µm], ascending.  Reported in full — this is a field.
        n_dof: Operator width.
        operator_norm: ``λ_max`` in absolute value [pN/µm], the scale the floor is derived from.
        zero_tolerance: ``n·ε·‖K‖₂`` [pN/µm] — below this an eigenvalue is not resolvable.
        n_zero_modes: Count of eigenvalues whose magnitude is at or below ``zero_tolerance``.
        n_negative_modes: Count of eigenvalues below ``−zero_tolerance`` (saddle directions).
        rigid_body_zero_modes: The caller's structural LOWER bound on the null-space dimension — six
            per connected component carrying no Dirichlet anchor.  It is a lower bound and not an
            equality, because a component can carry internal floppy modes on top of its rigid-body
            content (a head held to a backbone by one spring is free to swing about it), and those are
            a real property of the built model rather than an error.  ``None`` if none was declared.
        null_space_at_least_rigid_body: Whether the measured count reaches that bound.  FALSE is a
            hard defect — fewer zero modes than free rigid bodies means the operator is constraining
            a motion nothing physically constrains.
        n_excess_zero_modes: Measured minus rigid-body count: the INTERNAL floppy content of the
            built model, reported as a measurement.
        lambda_max: Largest eigenvalue [pN/µm].
        lambda_min_nonzero: Smallest resolvable positive eigenvalue [pN/µm], or ``None``.
        condition_number: ``λ_max / λ_min_nonzero``, dimensionless, or ``None``.
        timescale_separation_decades: ``log10`` of the condition number, or ``None``.
        stable_explicit_mobility_step: ``2/λ_max`` [µm/pN] — the exact forward-Euler stability limit of
            the overdamped update ``x ← x + (dt/γ)·F``, against which the run's Gershgorin CFL step is a
            conservative under-estimate by a factor this measurement finally makes visible.
        gershgorin_bound: ``max_i Σ_j |K_ij|`` [pN/µm] — the Gershgorin row-sum bound computed FROM THE
            SAME assembled operator.  Every explicit step in this repo is sized by a hand-assembled
            version of this bound; computing it from the operator rather than re-deriving it is what
            makes the comparison below a measurement of the bound's conservatism rather than of two
            different formulas.
        gershgorin_over_lambda_max: How many times larger the bound is than the truth, dimensionless.
            An explicit step sized by the bound is exactly this factor smaller than it needs to be.
    """

    eigenvalues: np.ndarray
    n_dof: int
    operator_norm: float
    zero_tolerance: float
    n_zero_modes: int
    n_negative_modes: int
    rigid_body_zero_modes: int | None
    null_space_at_least_rigid_body: bool | None
    n_excess_zero_modes: int | None
    lambda_max: float
    lambda_min_nonzero: float | None
    condition_number: float | None
    timescale_separation_decades: float | None
    stable_explicit_mobility_step: float | None
    gershgorin_bound: float
    gershgorin_over_lambda_max: float | None
    eigenvectors: np.ndarray | None = field(default=None, repr=False)

    def as_artifact_fields(self, *, n_report: int = 64) -> dict[str, Any]:
        """Return the JSON-able record, keeping the lowest ``n_report`` eigenvalues verbatim.

        The full array stays on the object (and belongs in the run's ``npz``); what goes into the JSON
        record is the low end, because that is the end that decides convergence and time-scale
        separation, plus every scalar this spectrum determines.

        Args:
            n_report: How many of the lowest eigenvalues to carry into the JSON record.

        Returns:
            A dict of plain Python types.
        """
        low = self.eigenvalues[: max(0, int(n_report))]
        return {
            "n_dof": self.n_dof,
            "operator_norm_pN_per_um": self.operator_norm,
            "zero_tolerance_pN_per_um": self.zero_tolerance,
            "n_zero_modes": self.n_zero_modes,
            "n_negative_modes": self.n_negative_modes,
            "rigid_body_zero_modes": self.rigid_body_zero_modes,
            "null_space_at_least_rigid_body": self.null_space_at_least_rigid_body,
            "n_excess_zero_modes": self.n_excess_zero_modes,
            "lambda_max_pN_per_um": self.lambda_max,
            "lambda_min_nonzero_pN_per_um": self.lambda_min_nonzero,
            "condition_number": self.condition_number,
            "timescale_separation_decades": self.timescale_separation_decades,
            "stable_explicit_mobility_step_um_per_pN": self.stable_explicit_mobility_step,
            "gershgorin_bound_pN_per_um": self.gershgorin_bound,
            "gershgorin_over_lambda_max": self.gershgorin_over_lambda_max,
            "lowest_eigenvalues_pN_per_um": [float(v) for v in low],
        }


def stiffness_spectrum(
    matrix: np.ndarray,
    *,
    rigid_body_zero_modes: int | None = None,
    keep_vectors: bool = True,
) -> StiffnessSpectrum:
    """Diagonalise an assembled tangent and derive what its spectrum decides.

    Args:
        matrix: The assembled ``(3n, 3n)`` operator [pN/µm], measured at ZERO regularization.  It is
            symmetrised as ``(K + Kᵀ)/2`` before diagonalisation; report the discarded antisymmetric
            part separately rather than relying on this step to be harmless.
        rigid_body_zero_modes: The structural LOWER bound on the null space — six per connected
            component carrying no Dirichlet anchor.  Supplied as a CHECK on the operator; any excess
            over it is the built model's internal floppy content and is reported, not absorbed.
        keep_vectors: Whether to retain the eigenvectors (the POD basis).  They are ``(3n, 3n)``
            float64, so a caller that only wants scalars can decline them.

    Returns:
        The :class:`StiffnessSpectrum`.

    Raises:
        ValueError: If ``matrix`` is not square or not finite.
    """
    k = np.asarray(matrix, np.float64)
    if k.ndim != 2 or k.shape[0] != k.shape[1]:
        raise ValueError(f"operator must be square; got shape {k.shape}")
    if not np.all(np.isfinite(k)):
        raise ValueError("operator contains non-finite entries")

    symmetric = 0.5 * (k + k.T)
    gershgorin = float(np.max(np.sum(np.abs(symmetric), axis=1))) if k.size else 0.0
    values, vectors = np.linalg.eigh(symmetric)
    n_dof = int(k.shape[0])
    norm = float(np.max(np.abs(values))) if values.size else 0.0
    tolerance = float(n_dof) * _EPS64 * norm

    n_zero = int(np.count_nonzero(np.abs(values) <= tolerance))
    n_negative = int(np.count_nonzero(values < -tolerance))
    positive = values[values > tolerance]
    lambda_max = float(values[-1]) if values.size else 0.0
    lambda_min_nonzero = float(positive[0]) if positive.size else None

    condition: float | None = None
    decades: float | None = None
    if lambda_min_nonzero is not None and lambda_min_nonzero > 0.0 and lambda_max > 0.0:
        condition = lambda_max / lambda_min_nonzero
        decades = float(np.log10(condition))

    stable_step = (2.0 / lambda_max) if lambda_max > 0.0 else None
    at_least_rigid: bool | None = None
    excess: int | None = None
    if rigid_body_zero_modes is not None:
        at_least_rigid = bool(n_zero >= int(rigid_body_zero_modes))
        excess = n_zero - int(rigid_body_zero_modes)

    return StiffnessSpectrum(
        eigenvalues=values,
        n_dof=n_dof,
        operator_norm=norm,
        zero_tolerance=tolerance,
        n_zero_modes=n_zero,
        n_negative_modes=n_negative,
        rigid_body_zero_modes=None if rigid_body_zero_modes is None else int(rigid_body_zero_modes),
        null_space_at_least_rigid_body=at_least_rigid,
        n_excess_zero_modes=excess,
        lambda_max=lambda_max,
        lambda_min_nonzero=lambda_min_nonzero,
        condition_number=condition,
        timescale_separation_decades=decades,
        stable_explicit_mobility_step=stable_step,
        gershgorin_bound=gershgorin,
        gershgorin_over_lambda_max=(gershgorin / lambda_max) if lambda_max > 0.0 else None,
        eigenvectors=vectors if keep_vectors else None,
    )


@dataclass(frozen=True, slots=True)
class ModeLocalization:
    """How far one eigenmode reaches through the slice — the stress decay length, measured.

    Attributes:
        index: Index into the ascending eigenvalue array.
        eigenvalue: The mode's eigenvalue [pN/µm].
        participation_ratio: ``(Σ w_i)² / Σ w_i²`` over per-node energies ``w_i = |u_i|²``, normalised
            so ``Σ w_i = 1``; equals the number of nodes the mode effectively occupies.  ``1`` is a
            single-node mode, ``n`` is a perfectly delocalised one.
        participation_fraction: ``participation_ratio / n_nodes``, dimensionless.
        gyration_radius_um: The energy-weighted RMS distance from the mode's own energy centroid [µm] —
            the mode's spatial extent in physical length, which is what "decay length" means.
    """

    index: int
    eigenvalue: float
    participation_ratio: float
    participation_fraction: float
    gyration_radius_um: float


def mode_localization(
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    positions_um: np.ndarray,
    *,
    indices: list[int] | tuple[int, ...] | np.ndarray,
) -> list[ModeLocalization]:
    """Measure the spatial extent of selected eigenmodes over the slice's node positions.

    Args:
        eigenvalues: Ascending eigenvalues [pN/µm], shape ``(3n,)``.
        eigenvectors: Column-wise eigenvectors, shape ``(3n, 3n)``.
        positions_um: Node positions [µm], shape ``(n, 3)``, in the SAME index order the operator uses.
        indices: Which modes to measure.

    Returns:
        One :class:`ModeLocalization` per requested index, in the order requested.

    Raises:
        ValueError: On any shape mismatch, or an index outside the spectrum.
    """
    values = np.asarray(eigenvalues, np.float64).reshape(-1)
    vectors = np.asarray(eigenvectors, np.float64)
    positions = np.asarray(positions_um, np.float64)
    n_nodes = int(positions.shape[0])
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError(f"positions_um must be (n, 3); got {positions.shape}")
    if vectors.shape != (3 * n_nodes, values.size) or values.size != 3 * n_nodes:
        raise ValueError(
            f"eigen-decomposition shape {vectors.shape} / {values.size} does not match "
            f"{n_nodes} nodes"
        )

    out: list[ModeLocalization] = []
    for raw in indices:
        index = int(raw)
        if not (0 <= index < values.size):
            raise ValueError(f"mode index {index} outside the spectrum of size {values.size}")
        mode = vectors[:, index].reshape(n_nodes, 3)
        weight = np.einsum("ij,ij->i", mode, mode)
        total = float(weight.sum())
        if total <= 0.0:
            raise ValueError(f"mode {index} has zero norm")
        weight = weight / total
        participation = 1.0 / float(np.sum(weight * weight))
        centroid = weight @ positions
        offset = positions - centroid
        gyration = float(np.sqrt(np.sum(weight * np.einsum("ij,ij->i", offset, offset))))
        out.append(
            ModeLocalization(
                index=index,
                eigenvalue=float(values[index]),
                participation_ratio=participation,
                participation_fraction=participation / float(n_nodes),
                gyration_radius_um=gyration,
            )
        )
    return out


# ── Iterative (Lanczos) spectrum: the only path to the operator at native population ──────────────────
#
# WHY THIS EXISTS.  ``assemble_dense_operator`` refuses above 8,192 DOF by design — column probing is
# ``O(n)`` operator applications and ``O(n²)`` memory, so at the native cell's ~1.5M DOF it is not merely
# slow but structurally impossible (a dense float64 operator would be ~19 PB).  Every spectral statement
# so far therefore came from a slice, and a slice's spectrum is a property of the slice: the largest
# eigenvalue follows the stiffest bond present, the smallest is a GLOBAL mode and a small system
# UNDERSTATES the separation, and the null-space fraction has the slice's own composition in its
# denominator.  Lanczos removes that limit: it touches the operator only through matrix-vector products,
# which the engine already performs on device every inner iteration.
#
# HOW THE SMALL END IS REACHED WITHOUT A SOLVE.  Lanczos converges first at the extremes of the spectrum,
# so ``λ_max`` is cheap and ``λ_min`` is not — the usual remedy, shift-and-invert, needs a linear solve
# per iteration and would make the measurement depend on the very solver it is meant to diagnose.  We
# use spectral FOLDING instead: run the same recurrence on ``σI − K`` with ``σ`` any upper bound on the
# spectrum (the Gershgorin row-sum bound the explicit step is already sized by).  That operator's LARGEST
# eigenvalues are ``σ − λ_min``, i.e. K's smallest, and it costs one extra axpy per product.
#
# WHAT IS REPORTED, AND WHAT IS NOT.  Every Ritz value carries its own rigorous residual bound
# ``|β_m · s_{m,i}|``: for a symmetric operator, some true eigenvalue lies within that distance of the
# Ritz value (Parlett, *The Symmetric Eigenvalue Problem*, §13.2).  So a converged extreme announces
# itself and an unconverged interior one does too, and nothing needs a chosen tolerance to be believed.
# What Lanczos does NOT give is the full eigenvalue field — the interior stays unresolved at any
# affordable iteration count — so it answers "what are the extremes, and how far apart", never "how many
# modes are floppy", which remains a dense-probe question and therefore a per-object one.


@dataclass(frozen=True, slots=True)
class LanczosSpectrum:
    """Extremal Ritz values of a matrix-free symmetric operator, each with a rigorous error bound.

    Attributes:
        ritz_values: The Ritz values [pN/µm], ascending.  These approximate eigenvalues of the operator
            the recurrence was run on — which, for the folded pass, is not ``K``; use
            :func:`lanczos_extremal_spectrum` if you want both ends of ``K`` itself.
        residual_bounds: ``|β_m · s_{m,i}|`` per Ritz value [pN/µm].  A true eigenvalue lies within this
            distance of the corresponding Ritz value.  Small at the extremes, large in the interior —
            that asymmetry is the method, not a defect, and it is reported rather than hidden.
        n_iterations: Krylov dimension actually built (may be below the request if the recurrence broke
            down on an invariant subspace, which is exact convergence and is recorded as such).
        n_dof: Operator width.
        breakdown: Whether the recurrence terminated early on an invariant subspace.
        reorthogonalized: Whether full reorthogonalization was applied.  Without it, Lanczos loses
            orthogonality once a Ritz value converges and produces spurious COPIES of converged
            eigenvalues — which would be read as degeneracy, i.e. as physics.
        start_seed: The seed of the random start vector, so the run is reproducible.
    """

    ritz_values: np.ndarray = field(repr=False)
    residual_bounds: np.ndarray = field(repr=False)
    n_iterations: int
    n_dof: int
    breakdown: bool
    reorthogonalized: bool
    start_seed: int

    @property
    def extreme_high(self) -> tuple[float, float]:
        """The largest Ritz value and its residual bound [pN/µm]."""
        return float(self.ritz_values[-1]), float(self.residual_bounds[-1])

    @property
    def extreme_low(self) -> tuple[float, float]:
        """The smallest Ritz value and its residual bound [pN/µm]."""
        return float(self.ritz_values[0]), float(self.residual_bounds[0])

    def as_artifact_fields(self) -> dict[str, Any]:
        """Return the JSON-able record; the Ritz field is kept in full, never summarised."""
        return {
            "ritz_values_pN_per_um": self.ritz_values.tolist(),
            "residual_bounds_pN_per_um": self.residual_bounds.tolist(),
            "n_iterations": self.n_iterations,
            "n_dof": self.n_dof,
            "breakdown_on_invariant_subspace": self.breakdown,
            "reorthogonalized": self.reorthogonalized,
            "start_seed": self.start_seed,
        }


def lanczos_ritz(
    matvec: Any,
    n_nodes: int,
    *,
    n_iterations: int,
    seed: int = 0,
    reorthogonalize: bool = True,
    start_vector: np.ndarray | None = None,
) -> LanczosSpectrum:
    """Run Lanczos on a matrix-free node operator and return its Ritz values with error bounds.

    Args:
        matvec: The matrix-free application, ``(n, 3) -> (n, 3)``, linear and symmetric.  Symmetry is
            the method's premise, not something it checks — pair it with
            :func:`~aleph.engine.observe.operator_probe.symmetry_report` on a slice, where the
            operator can be assembled, before trusting it at a population where it cannot.
        n_nodes: Node count ``n``; the operator is ``3n × 3n``.
        n_iterations: Krylov dimension.  Cost is one operator application per iteration; with
            reorthogonalization the memory is ``3n × n_iterations`` float64, which the caller must budget
            (that product, not the operator, is what bounds this method at native population).
        seed: Seed of the random start vector.  A random start is used deliberately: a structured one
            can be orthogonal to the very eigenvector being sought, and the method would then converge to
            a spectrum missing it, silently.  Ignored when ``start_vector`` is supplied.
        reorthogonalize: Apply full reorthogonalization against the stored basis.
        start_vector: An explicit start, shape ``(n, 3)`` or ``(3n,)``, normalized here.  It exists for
            one purpose: when the operator is a PROJECTED one (``P K P``, the operator this engine's CG
            actually inverts), the constraint directions are an exact null space that a random start
            would put into the Krylov space, and the recurrence would then report the projector's
            zeros as the operator's small end.  Starting inside ``range(P)`` — which ``P K P`` leaves
            invariant — measures the RESTRICTION instead, i.e. the spectrum the solve really sees.
            A random start is still the right default everywhere else, for the reason under ``seed``.

    Returns:
        The :class:`LanczosSpectrum`.

    Raises:
        ValueError: If the sizes are not positive, if ``n_iterations`` exceeds the operator width, if
            ``start_vector`` has the wrong shape or vanishing norm, or if ``matvec`` returns a
            non-finite or mis-shaped array.
    """
    n = int(n_nodes)
    if n <= 0:
        raise ValueError(f"n_nodes must be positive; got {n_nodes!r}")
    n_dof = 3 * n
    steps = int(n_iterations)
    if steps < 1:
        raise ValueError(f"n_iterations must be at least 1; got {n_iterations!r}")
    if steps > n_dof:
        raise ValueError(
            f"n_iterations={steps} exceeds the operator width {n_dof}; a Krylov space cannot be larger "
            "than the space it lives in"
        )

    if start_vector is None:
        rng = np.random.default_rng(int(seed))
        q = rng.standard_normal(n_dof)
    else:
        q = np.asarray(start_vector, np.float64).reshape(-1).copy()
        if q.size != n_dof:
            raise ValueError(f"start_vector has {q.size} entries; the operator is {n_dof}-wide")
        if not np.all(np.isfinite(q)):
            raise ValueError("start_vector contains non-finite entries")
    start_norm = float(np.linalg.norm(q))
    if not (start_norm > 0.0):
        raise ValueError(
            "start_vector has zero norm; there is no direction to build a Krylov space from"
        )
    q /= start_norm

    basis = np.zeros((steps, n_dof), np.float64) if reorthogonalize else None
    alpha = np.zeros(steps, np.float64)
    beta = np.zeros(steps, np.float64)
    previous = np.zeros(n_dof, np.float64)
    beta_previous = 0.0
    breakdown = False
    built = steps

    # The breakdown floor is the float64 round-off of the vector norm, not a chosen tolerance: below it
    # the next basis vector is numerical noise rather than a direction of the operator.
    breakdown_floor = float(np.sqrt(n_dof)) * _EPS64

    for index in range(steps):
        if basis is not None:
            basis[index] = q
        response = np.asarray(matvec(q.reshape(n, 3)), np.float64)
        if response.shape != (n, 3):
            raise ValueError(f"matvec returned shape {response.shape}, expected {(n, 3)}")
        if not np.all(np.isfinite(response)):
            raise ValueError(f"matvec returned a non-finite response at iteration {index}")
        w = response.reshape(-1) - beta_previous * previous
        alpha[index] = float(np.dot(q, w))
        w -= alpha[index] * q
        if basis is not None:
            # Twice is enough (Kahan): one pass restores orthogonality to round-off, the second one
            # guards the case where the first pass itself cancelled catastrophically.
            for _ in range(2):
                w -= basis[: index + 1].T @ (basis[: index + 1] @ w)
        beta_value = float(np.linalg.norm(w))
        if index + 1 < steps:
            beta[index] = beta_value
        if beta_value <= breakdown_floor:
            breakdown = True
            built = index + 1
            break
        previous = q
        q = w / beta_value
        beta_previous = beta_value

    tridiagonal = np.diag(alpha[:built])
    if built > 1:
        off = beta[: built - 1]
        tridiagonal += np.diag(off, 1) + np.diag(off, -1)
    values, vectors = np.linalg.eigh(tridiagonal)
    last_beta = 0.0 if breakdown or built < steps else float(beta[built - 1])
    if built >= 1 and not breakdown and built == steps:
        # beta[built-1] was deliberately not stored (it is the residual norm, not a matrix entry), so
        # recompute it from the final w — carried in `beta_value` above.
        last_beta = float(beta_value)
    bounds = np.abs(last_beta * vectors[-1, :])

    return LanczosSpectrum(
        ritz_values=values,
        residual_bounds=bounds,
        n_iterations=built,
        n_dof=n_dof,
        breakdown=breakdown,
        reorthogonalized=bool(reorthogonalize),
        start_seed=int(seed),
    )


def lanczos_extremal_spectrum(
    matvec: Any,
    n_nodes: int,
    *,
    upper_bound: float,
    n_iterations: int,
    seed: int = 0,
    reorthogonalize: bool = True,
    projector: Any | None = None,
) -> dict[str, Any]:
    """Measure BOTH ends of a matrix-free symmetric spectrum: ``λ_max`` directly, ``λ_min`` by folding.

    The second pass runs the identical recurrence on ``σI − K``.  Its largest Ritz values are
    ``σ − λ_min``, so the small end converges as fast as a large end does — with no linear solve, hence
    with no dependence on the solver being diagnosed.

    Args:
        matvec: The matrix-free application, ``(n, 3) -> (n, 3)``.
        n_nodes: Node count.
        upper_bound: ``σ``, any upper bound on the spectrum [pN/µm] — the Gershgorin row-sum bound the
            explicit step is already sized by is the natural choice, and comparing it against the
            measured ``λ_max`` is itself the measurement of how conservative that step is.
        n_iterations: Krylov dimension for each pass.
        seed: Start-vector seed; the folded pass uses ``seed + 1`` so the two passes are independent.
        reorthogonalize: Full reorthogonalization (see :func:`lanczos_ritz`).
        projector: An orthogonal projector ``P``, ``(n, 3) -> (n, 3)``, when ``matvec`` is the PROJECTED
            operator ``P K P`` — the operator this engine's CG inverts.  Supplying it changes two
            things, and both are required for the answer to be about ``K`` rather than about the
            constraints: both passes start from ``P(random)`` so the Krylov space lives in
            ``range(P)``, and the folded operator becomes ``σ·P(v) − matvec(v)`` rather than
            ``σ·v − matvec(v)``.  Without the second substitution the identity term would carry the
            constraint directions at eigenvalue ``σ`` — the folded operator's LARGEST — so the folded
            pass would converge to the projector's null space and report ``λ_min = 0`` for any
            operator whatsoever.  ``None`` (default) measures the unprojected operator.

    Returns:
        A dict carrying both passes' records, the recovered ``λ_max``/``λ_min`` with their bounds, the
        ``upper_bound`` and its ratio to the measured ``λ_max``, and — when both ends are resolved and
        positive — the condition number and time-scale separation in decades.  With a ``projector``
        supplied, every one of those quantities describes the RESTRICTION of the operator to
        ``range(P)``, which the record states in its ``subspace`` key.

    Raises:
        ValueError: If ``upper_bound`` is not finite, or is below the measured ``λ_max`` (which makes the
            folded operator's extreme the WRONG end, so the result would be silently meaningless).
    """
    if not np.isfinite(upper_bound):
        raise ValueError(f"upper_bound must be finite; got {upper_bound!r}")

    n = int(n_nodes)
    start: np.ndarray | None = None
    folded_start: np.ndarray | None = None
    if projector is not None:
        if n <= 0:
            raise ValueError(f"n_nodes must be positive; got {n_nodes!r}")
        start = np.asarray(
            projector(np.random.default_rng(int(seed)).standard_normal((n, 3))), np.float64)
        folded_start = np.asarray(
            projector(np.random.default_rng(int(seed) + 1).standard_normal((n, 3))), np.float64)

    high = lanczos_ritz(matvec, n_nodes, n_iterations=n_iterations, seed=seed,
                        reorthogonalize=reorthogonalize, start_vector=start)
    lambda_max, lambda_max_bound = high.extreme_high
    if float(upper_bound) < lambda_max:
        raise ValueError(
            f"upper_bound={upper_bound!r} is below the measured lambda_max={lambda_max!r}; folding "
            "against it would return the wrong end of the spectrum instead of failing"
        )

    def folded(vector: np.ndarray) -> np.ndarray:
        identity_term = (np.asarray(vector, np.float64) if projector is None
                         else np.asarray(projector(vector), np.float64))
        return float(upper_bound) * identity_term - np.asarray(matvec(vector), np.float64)

    low_pass = lanczos_ritz(folded, n, n_iterations=n_iterations, seed=int(seed) + 1,
                            reorthogonalize=reorthogonalize, start_vector=folded_start)
    folded_high, folded_bound = low_pass.extreme_high
    lambda_min = float(upper_bound) - folded_high

    # A Ritz value is RESOLVED when its own rigorous bound is small compared to itself.  This matters far
    # more at the small end than the large one and the asymmetry is structural, not a tuning failure: the
    # folded operator's top eigenvalues are separated by ``(λ₂ − λ₁)/(σ − λ₁)``, so the more
    # ill-conditioned ``K`` is, the closer together the folded extremes sit and the slower the small end
    # converges.  Measured on a well-conditioned reference the large end is exact to round-off while the
    # small end can still be tens of percent out — so the verdict is reported per end, and a run that
    # quotes an unresolved ``λ_min`` (or a condition number built from one) is quoting the start vector.
    lambda_max_resolved = bool(lambda_max > 0.0 and lambda_max_bound < 0.01 * abs(lambda_max))
    lambda_min_resolved = bool(lambda_min > 0.0 and folded_bound < 0.01 * abs(lambda_min))

    record: dict[str, Any] = {
        "lambda_max_pN_per_um": lambda_max,
        "lambda_max_residual_bound": lambda_max_bound,
        "lambda_max_resolved": lambda_max_resolved,
        "lambda_min_pN_per_um": lambda_min,
        "lambda_min_residual_bound": folded_bound,
        "lambda_min_resolved": lambda_min_resolved,
        "resolution_criterion": "a Ritz value is called resolved when its rigorous residual bound "
                                "|beta_m·s_mi| is below 1% of the value itself; the bound is Parlett's, "
                                "so this is a statement about where a true eigenvalue must lie, not a "
                                "convergence heuristic",
        "upper_bound_pN_per_um": float(upper_bound),
        "upper_bound_over_lambda_max": (
            float(upper_bound) / lambda_max if lambda_max > 0.0 else None),
        "stable_explicit_mobility_step_um_per_pN": (
            2.0 / lambda_max if lambda_max > 0.0 else None),
        "folding": "lambda_min measured as sigma − lambda_max(sigma·I − K); no linear solve, so the "
                   "measurement does not depend on the solver it diagnoses",
        "subspace": (
            "the whole space: no projector was supplied, so these are eigenvalues of the operator as "
            "given" if projector is None else
            "range(P): both passes start inside the projector's range and the folded operator uses "
            "sigma·P in place of sigma·I, so every value here is an eigenvalue of the RESTRICTION of "
            "the projected operator — the constraint null space is excluded by construction rather "
            "than by discarding small values afterwards"
        ),
        "direct_pass": high.as_artifact_fields(),
        "folded_pass": low_pass.as_artifact_fields(),
    }
    if lambda_max > 0.0 and lambda_min > 0.0 and lambda_max_resolved and lambda_min_resolved:
        record["condition_number"] = lambda_max / lambda_min
        record["timescale_separation_decades"] = float(np.log10(lambda_max / lambda_min))
    elif lambda_max > 0.0 and lambda_min > 0.0:
        record["condition_number"] = None
        record["timescale_separation_decades"] = None
        record["condition_number_note"] = (
            "one end is not resolved at this Krylov dimension, so a condition number formed from it "
            "would be a property of the iteration count rather than of the operator; the two Ritz values "
            "and their bounds are reported instead"
        )
    else:
        record["condition_number"] = None
        record["timescale_separation_decades"] = None
        record["condition_number_note"] = (
            "lambda_min is at or below zero within the Krylov space reached, so the operator's small end "
            "is a null/near-null direction here and a condition number is not defined from these passes"
        )
    return record
