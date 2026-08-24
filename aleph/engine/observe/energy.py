r"""Energy accounting that needs no energy function: closed-loop work, path work, and the balance.

WHY A LOOP INTEGRAL COMES FIRST.  The execution plan makes energy the engine's common language, and the
strongest thing it asks for — "the balance ``dU/dt + dissipation + Σ(event jumps) = P_active`` closes" —
presumes a potential ``U`` that this engine has never written down: forces are assembled kernel by
kernel, and no module returns an energy.  Writing one per force family would be a second implementation
of the same physics, free to drift from the forces it claims to differentiate, which is the
"gate a broken build cannot fail" antipattern.

The loop integral avoids that entirely.  For any conservative field ``F = −∇U`` the work around ANY
closed path is exactly zero, whatever ``U`` is.  So::

    W_loop = ∮ F · dx

is a direct, magnitude-free test that an energy EXISTS for the forces the slice actually launches — and
it uses only force evaluations, the same launches the run already performs.  A nonzero ``W_loop`` means
the assembled force field creates or destroys energy on a cycle, which is what every double-count and
every broken adjoint pair does.  All six of the PI framework's documented traps are "counting the same
physics twice", and a term counted twice injects energy twice.

This is the same defect that :func:`~aleph.engine.observe.operator_probe.symmetry_report` detects
through ``‖K − Kᵀ‖``, reached from the opposite side, and the pair is worth having: the symmetry test
uses the TANGENT and is therefore local and exact but blind to any force whose tangent is not launched;
the loop test uses the FORCE at finite amplitude and so probes the assembled field itself, including
anything the tangent omits.  Disagreement between them localises the defect to the tangent.

WHAT ``EnergyBalance`` IS FOR.  Once a lane does own an energy and a dissipation channel, the closure
gate is the residual of the balance, reported as a dimensionless ratio against the largest single term.
The dataclass is the shape of that record, and it deliberately carries no threshold: what counts as
closed is a gate contract, not a property of this module.

HOW THE BALANCE IS FORMED WITHOUT WRITING A SECOND IMPLEMENTATION OF THE PHYSICS.  Implementing a
potential per force family would be that second implementation — free to drift from the forces it claims
to differentiate — so :func:`step_energy_balance` forms every term out of the SAME force launches the
step already performs:

``ΔU``
    ``−∫F_passive·dx`` along the step's own trajectory.  Legitimate because the passive field has been
    MEASURED conservative (loop residual on the quadrature curve), so this integral is a state function;
    the same integral along the straight chord between the endpoints is computed as well, and the two
    agreeing is the path-independence check that licenses calling it a potential difference at all.
``W_active``
    the same walk with the force evaluated TWICE per midpoint — heads bound and heads detached — and
    differenced.  The difference is exactly the crossbridge, i.e. the term the conservativity measurement
    identified as the non-conservative content, so the active channel is defined by the control that
    attributed it rather than by a separate model of it.
``dissipated``
    ``Σ γ|Δx|²/dt = Σ|Δx|²/μ`` over the trajectory at the run's own overdamped mobility ``μ = dt/γ``,
    accumulated ON DEVICE per inner iteration by :data:`accumulate_squared_displacement_kernel` so it is
    exact rather than subsampled (a line integral survives subsampling; a sum of squared displacements
    does not).  The "time" underneath it is the inner relaxation's own pseudo-time at the run's mobility:
    the balance is an identity of the trajectory the solver actually walked, not a statement about the
    outer physical clock.
``event_jump``
    zero BY DERIVATION on a lane whose only binding-dependent force is the crossbridge: the KMC commits
    at fixed configuration and the passive potential does not depend on the binding state at all, so no
    energy discontinuity enters ``U``.  It stays an explicit term because a lane with binding-dependent
    PASSIVE bonds (a crosslink turnover) does have one, and it must then be supplied rather than assumed.

The residual is therefore not a tautology: ``dissipated`` is a KINEMATIC measurement of the displacements
and the two work terms are FIELD integrals, so the balance closing says the discrete update is consistent
with the forces it was built from, and its failure to close localises a leak, a double count, or a
mis-accounted event.

WHAT THE RESIDUAL IS WHEN IT IS NOT ZERO, AND HOW IT SAYS SO.  It does not vanish even for a perfect
accounting, and the reason is worth stating exactly: an explicit overdamped step takes
``Δx = μ·F(x_old)``, so ``Σ|Δx|²/μ = Σ F(x_old)·Δx`` is a LEFT-endpoint quadrature of the same line
integral the work terms take by the midpoint rule, and the two differ by ``O(μ)`` per step — the
integrator's own first-order consistency error, not an accounting failure.  Refining the trajectory
sampling does not remove it (measured: it converges to a fixed value as the polyline approaches the
true trajectory); halving the MOBILITY does, exactly linearly.  A genuine leak — a force channel the
solver used and the ledger's evaluation omits, a double count, a one-sided adjoint — is a work integral
in its own right and is independent of ``μ``.  :func:`balance_convergence` compares the measured exponent
against those two structural predictions, so the residual identifies itself and no tolerance is chosen
anywhere.  This is the same argument :func:`loop_convergence` makes one level down.

engine units: length µm, force pN, work pN·µm (= 1e-18 J).

Sanity Gate:
    * dimensional: ``F [pN] · dx [µm]`` -> ``[pN·µm]``.  Closure ratios are dimensionless.
    * boundary: a loop of zero amplitude has zero work and an undefined ratio, so amplitude must be
      positive; a force field that is identically zero yields ``W_loop = 0`` with a zero scale, and the
      ratio is reported as ``0.0`` rather than ``nan``.
    * conservation/invariant: the loop work of a conservative field is identically zero at any
      amplitude and any segment count — that invariance is the test, and both are recorded so the
      caller can vary them.
    * numerical: the midpoint rule is second-order and EXACT for a linear force field (quadratic
      ``U``), so any residual at small amplitude is field asymmetry rather than quadrature error; the
      round-off floor is the same order-independent float64 accumulation bound the balance gate uses
      (PI D8), never a chosen epsilon.
    * sign-sense: work done BY the field along the path is positive; a conservative field returns
      exactly what it took, hence zero around a cycle.
    * measurement-protocol: the force callable must be a pure function of position with all discrete
      state (bound heads, contacts) FROZEN.  A KMC event fired mid-loop makes the path integral measure
      the event, not the field, and the freeze is the caller's responsibility.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
import warp as wp

from aleph.engine.ledger import assemble_balance_tolerance_ratio

__all__ = [
    "BalanceConvergence",
    "ClosedLoopWork",
    "EnergyBalance",
    "LoopConvergence",
    "StepEnergyLedger",
    "accumulate_squared_displacement_kernel",
    "balance_convergence",
    "closed_loop_work",
    "dissipated_work",
    "loop_convergence",
    "paired_path_work",
    "path_work",
    "step_energy_balance",
]

#: The amplitude exponent of the midpoint rule's truncation error around a closed circuit: per segment
#: the error is ``O(h³·f'')`` and ``4N`` segments of length ``a/N`` give ``O(a³/N²)``.
QUADRATURE_AMPLITUDE_EXPONENT: float = 3.0

#: The amplitude exponent of a genuine circulation: ``∮F·dx`` over a circuit of edge ``a`` encloses an
#: area ``a²``, so a field with nonzero curl scales as ``a²`` — and, unlike quadrature error, does not
#: depend on the segment count at all.
CIRCULATION_AMPLITUDE_EXPONENT: float = 2.0

#: The segment exponent of the midpoint rule's truncation error at fixed amplitude: ``O(N^-2)``.
QUADRATURE_SEGMENT_EXPONENT: float = 2.0

#: The segment exponent of a genuine circulation at fixed amplitude: none — it is exactly N-independent.
CIRCULATION_SEGMENT_EXPONENT: float = 0.0

#: The mobility exponent of a first-order overdamped integrator's own consistency error in the energy
#: balance.  ``Δx = μ·F(x_old)`` makes the kinematic dissipation ``Σ|Δx|²/μ`` a LEFT-endpoint quadrature
#: of ``∫F·dx`` while the ledger's work terms use the midpoint rule, and the two differ by ``O(μ)`` per
#: step.  So the balance of an explicit step closes only as fast as the step is consistent — which is a
#: property of the integrator, measurable, and not a defect in the accounting.
INTEGRATOR_CONSISTENCY_MOBILITY_EXPONENT: float = 1.0

#: The mobility exponent of a GENUINE accounting failure — a force channel the solver used but the
#: ledger's evaluation omits, a double count, or a broken adjoint pair.  Such a residual is a work
#: integral in its own right, so refining the mobility leaves it where it is.  The two exponents are
#: what let a residual say which of the two it is, with no tolerance chosen anywhere.
LEAK_MOBILITY_EXPONENT: float = 0.0

#: A force evaluation at a configuration: ``(n, 3) [µm] -> (n, 3) [pN]``, with discrete state frozen.
ForceFn = Callable[[np.ndarray], np.ndarray]


def path_work(force: ForceFn, path: np.ndarray) -> tuple[float, float]:
    """Return ``(∫F·dx, Σ|F·dx|)`` along a polyline in configuration space by the midpoint rule.

    Args:
        force: The force evaluation, ``(n, 3) -> (n, 3)``.
        path: Configurations, shape ``(m, n, 3)`` [µm], visited in order.  ``m >= 2``.

    Returns:
        The signed work [pN·µm] and the total unsigned work exchanged [pN·µm]; the latter is the scale
        a closure ratio is taken against, so the ratio compares the residual to the traffic rather than
        to an invented reference.

    Raises:
        ValueError: If the path has fewer than two configurations, is not ``(m, n, 3)``, or the force
            returns a non-finite or wrongly-shaped value.
    """
    points = np.asarray(path, np.float64)
    if points.ndim != 3 or points.shape[2] != 3 or points.shape[0] < 2:
        raise ValueError(f"path must be (m>=2, n, 3); got shape {points.shape}")

    signed = 0.0
    unsigned = 0.0
    for index in range(points.shape[0] - 1):
        start = points[index]
        end = points[index + 1]
        midpoint = 0.5 * (start + end)
        value = np.asarray(force(midpoint), np.float64)
        if value.shape != start.shape:
            raise ValueError(f"force returned shape {value.shape}, expected {start.shape}")
        if not np.all(np.isfinite(value)):
            raise ValueError(f"force returned a non-finite value on segment {index}")
        increment = float(np.sum(value * (end - start)))
        signed += increment
        unsigned += abs(increment)
    return signed, unsigned


def paired_path_work(
    force_a: ForceFn, force_b: ForceFn, path: np.ndarray
) -> tuple[float, float, float, float]:
    """Return ``(∫F_a·dx, ∫F_b·dx, Σ|F_a·dx|, Σ|F_b·dx|)`` for TWO fields over one shared polyline.

    Two separate :func:`path_work` walks would evaluate the two fields at the same midpoints only by
    construction; doing one walk makes it structural.  That matters because the whole point of the pair
    is their DIFFERENCE (total minus passive = the active channel), and a difference of integrals taken
    at different quadrature points measures the quadrature as much as the physics.

    Args:
        force_a: First force evaluation, ``(n, 3) -> (n, 3)``.
        force_b: Second force evaluation at the SAME configurations, ``(n, 3) -> (n, 3)``.
        path: Configurations, shape ``(m, n, 3)`` [µm], visited in order.  ``m >= 2``.

    Returns:
        The two signed works and the two unsigned totals, all [pN·µm].

    Raises:
        ValueError: If the path is malformed or either force returns a non-finite or mis-shaped value.
    """
    points = np.asarray(path, np.float64)
    if points.ndim != 3 or points.shape[2] != 3 or points.shape[0] < 2:
        raise ValueError(f"path must be (m>=2, n, 3); got shape {points.shape}")

    signed_a = signed_b = unsigned_a = unsigned_b = 0.0
    for index in range(points.shape[0] - 1):
        start = points[index]
        end = points[index + 1]
        midpoint = 0.5 * (start + end)
        delta = end - start
        for force, label in ((force_a, "force_a"), (force_b, "force_b")):
            value = np.asarray(force(midpoint), np.float64)
            if value.shape != start.shape:
                raise ValueError(f"{label} returned shape {value.shape}, expected {start.shape}")
            if not np.all(np.isfinite(value)):
                raise ValueError(f"{label} returned a non-finite value on segment {index}")
            increment = float(np.sum(value * delta))
            if label == "force_a":
                signed_a += increment
                unsigned_a += abs(increment)
            else:
                signed_b += increment
                unsigned_b += abs(increment)
    return signed_a, signed_b, unsigned_a, unsigned_b


@wp.kernel
def accumulate_squared_displacement_kernel(
    position: wp.array(dtype=wp.vec3d),
    previous: wp.array(dtype=wp.vec3d),
    squared_sum: wp.array(dtype=wp.float64),
) -> None:
    """Accumulate ``Σ|x − x_prev|²`` [µm²] over one array, for the overdamped dissipation channel.

    Launched once per inner iteration on each relaxed array, with ``previous`` holding the configuration
    before that iteration's update.  Solver-agnostic on purpose: it reads the displacement the solver
    actually produced rather than assuming ``Δx = μF``, so the explicit and the implicit inner paths are
    measured by the same instrument and their answers may legitimately differ.

    ``squared_sum`` is a ``(1,)`` float64 accumulator; the caller zeroes it at the start of a step and
    reads it out of the loop.  No host read happens here.
    """
    i = wp.tid()
    delta = position[i] - previous[i]
    wp.atomic_add(squared_sum, 0, wp.dot(delta, delta))


def dissipated_work(squared_displacement_um2: float, mobility_um_per_pN: float) -> float:
    """Return the overdamped drag dissipation ``Σ|Δx|²/μ`` [pN·µm] for a measured displacement sum.

    Overdamped mechanics has ``γ ẋ = F`` and dissipation rate ``γ|ẋ|²``, so over a discrete update of
    displacement ``Δx`` in a sub-step ``dt`` the energy lost to drag is ``γ|Δx|²/dt = |Δx|²/μ`` with the
    mobility ``μ = dt/γ`` [µm/pN] the run is actually integrating with.  No constant is introduced: ``μ``
    is the run's own CFL-derived (or IMEX) mobility step.

    Args:
        squared_displacement_um2: ``Σ|Δx|²`` [µm²] accumulated over the step's inner iterations.
        mobility_um_per_pN: The mobility step ``dt/γ`` [µm/pN].  Must be positive-finite.

    Returns:
        The dissipated work [pN·µm], nonnegative by construction.

    Raises:
        ValueError: If the mobility is not positive-finite or the displacement sum is negative.
    """
    if not (mobility_um_per_pN > 0.0 and np.isfinite(mobility_um_per_pN)):
        raise ValueError(f"mobility must be positive-finite; got {mobility_um_per_pN!r}")
    if squared_displacement_um2 < 0.0 or not np.isfinite(squared_displacement_um2):
        raise ValueError(
            f"squared displacement must be nonnegative-finite; got {squared_displacement_um2!r}")
    return float(squared_displacement_um2) / float(mobility_um_per_pN)


@dataclass(frozen=True, slots=True)
class ClosedLoopWork:
    """The work a force field does around one closed rectangular circuit.

    Attributes:
        work: ``∮F·dx`` [pN·µm].  Exactly zero for a conservative field.
        exchanged: ``Σ|F·dx|`` around the circuit [pN·µm] — the scale the residual is judged against.
        closure_ratio: ``|work| / exchanged``, dimensionless; ``0.0`` when nothing was exchanged.
        round_off_floor: The float64 accumulation bound for the segment count, times ``exchanged``.
        at_round_off_floor: Whether ``|work|`` sits at or below that floor.  This is NOT by itself a
            conservativity verdict: a nonlinear force field has midpoint-rule truncation error far
            above the round-off floor even when it is perfectly conservative, so a ``False`` here means
            only "there is a residual to explain".  What explains it is
            :func:`loop_convergence`, which identifies the residual by its scaling.
        amplitude_um: The circuit's edge amplitude [µm].
        n_segments: Segments per edge (four edges).
    """

    work: float
    exchanged: float
    closure_ratio: float
    round_off_floor: float
    at_round_off_floor: bool
    amplitude_um: float
    n_segments: int

    def as_artifact_fields(self) -> dict[str, Any]:
        """Return the JSON-able record."""
        return {
            "loop_work_pN_um": self.work,
            "exchanged_work_pN_um": self.exchanged,
            "closure_ratio": self.closure_ratio,
            "round_off_floor": self.round_off_floor,
            "at_round_off_floor": self.at_round_off_floor,
            "amplitude_um": self.amplitude_um,
            "n_segments_per_edge": self.n_segments,
        }


def closed_loop_work(
    force: ForceFn,
    origin: np.ndarray,
    direction_a: np.ndarray,
    direction_b: np.ndarray,
    *,
    amplitude_um: float,
    n_segments: int = 8,
) -> ClosedLoopWork:
    """Integrate ``∮F·dx`` around the rectangle spanned by two configuration-space directions.

    The circuit is ``x0 → x0+a·u → x0+a·u+a·v → x0+a·v → x0`` with each edge split into ``n_segments``
    midpoint-rule segments.  Both directions are normalised, so ``amplitude_um`` is the edge length in
    configuration space and the result is comparable across different direction choices.

    Args:
        force: The force evaluation with all discrete state frozen, ``(n, 3) -> (n, 3)``.
        origin: The configuration the circuit starts and ends at, ``(n, 3)`` [µm].
        direction_a: First edge direction, ``(n, 3)``; normalised internally.
        direction_b: Second edge direction, ``(n, 3)``; normalised internally.
        amplitude_um: Edge length [µm].  Must be positive.  Small amplitudes make the midpoint rule
            exact to higher order, so a residual that does NOT shrink with amplitude is field
            asymmetry rather than quadrature error — vary it to tell the two apart.
        n_segments: Segments per edge.

    Returns:
        The :class:`ClosedLoopWork`.

    Raises:
        ValueError: If the amplitude or segment count is not positive, the shapes disagree, or either
            direction has zero norm.
    """
    x0 = np.asarray(origin, np.float64)
    if x0.ndim != 2 or x0.shape[1] != 3:
        raise ValueError(f"origin must be (n, 3); got {x0.shape}")
    if not (amplitude_um > 0.0 and np.isfinite(amplitude_um)):
        raise ValueError(f"amplitude_um must be positive-finite; got {amplitude_um!r}")
    if int(n_segments) < 1:
        raise ValueError(f"n_segments must be at least 1; got {n_segments!r}")

    def unit(vector: np.ndarray, label: str) -> np.ndarray:
        array = np.asarray(vector, np.float64)
        if array.shape != x0.shape:
            raise ValueError(f"{label} must match origin shape {x0.shape}; got {array.shape}")
        norm = float(np.linalg.norm(array))
        if norm <= 0.0:
            raise ValueError(f"{label} has zero norm, so it spans no circuit")
        return array / norm

    u = unit(direction_a, "direction_a") * float(amplitude_um)
    v = unit(direction_b, "direction_b") * float(amplitude_um)

    corners = (x0, x0 + u, x0 + u + v, x0 + v, x0)
    segments = int(n_segments)
    path = []
    for index in range(4):
        start, end = corners[index], corners[index + 1]
        for k in range(segments):
            path.append(start + (end - start) * (k / segments))
    path.append(x0)

    signed, unsigned = path_work(force, np.asarray(path, np.float64))
    floor = assemble_balance_tolerance_ratio(4 * segments) * unsigned
    return ClosedLoopWork(
        work=signed,
        exchanged=unsigned,
        closure_ratio=(abs(signed) / unsigned) if unsigned > 0.0 else 0.0,
        round_off_floor=float(floor),
        at_round_off_floor=bool(abs(signed) <= floor),
        amplitude_um=float(amplitude_um),
        n_segments=segments,
    )


@dataclass(frozen=True, slots=True)
class LoopConvergence:
    """Which of two structurally predicted power laws the measured loop work follows.

    A nonzero ``∮F·dx`` has exactly two possible origins, and they scale differently, so the residual's
    OWN scaling identifies it — no tolerance is chosen anywhere:

    ============================  =====================  ==========================
    origin                        amplitude ``a``        segments ``N``
    ============================  =====================  ==========================
    midpoint-rule truncation      ``a³``                 ``N⁻²``
    genuine circulation (curl)    ``a²``                 independent of ``N``
    ============================  =====================  ==========================

    So a field whose residual halves-cubed with amplitude AND falls fourfold with segment count is
    conservative and what was measured was the quadrature; a residual that scales as ``a²`` and does
    not move with ``N`` is a real non-conservative force, i.e. a double count or a broken adjoint pair.

    The verdict compares the two MEASURED exponents against those two PREDICTED pairs and reports
    which is nearer.  That is a comparison between structural predictions, not a threshold.

    Attributes:
        amplitude_exponent: ``log2(|W(a)| / |W(a/2)|)``, or ``None`` if either work is zero.
        segment_exponent: ``log2(|W(N)| / |W(2N)|)`` at fixed amplitude, or ``None``.
        amplitude_distance_to_quadrature: ``|amplitude_exponent − 3|``.
        amplitude_distance_to_circulation: ``|amplitude_exponent − 2|``.
        segment_distance_to_quadrature: ``|segment_exponent − 2|``.
        segment_distance_to_circulation: ``|segment_exponent − 0|``.
        consistent_with_quadrature: Whether BOTH exponents sit nearer their quadrature prediction.
        consistent_with_circulation: Whether BOTH sit nearer their circulation prediction.
        verdict: ``"quadrature"`` / ``"circulation"`` / ``"ambiguous"`` — ambiguous when the two
            exponents disagree with each other, which is a real outcome and is not resolved by picking
            one.
    """

    amplitude_exponent: float | None
    segment_exponent: float | None
    amplitude_distance_to_quadrature: float | None
    amplitude_distance_to_circulation: float | None
    segment_distance_to_quadrature: float | None
    segment_distance_to_circulation: float | None
    consistent_with_quadrature: bool
    consistent_with_circulation: bool
    verdict: str

    def as_artifact_fields(self) -> dict[str, Any]:
        """Return the JSON-able record."""
        return {
            "amplitude_exponent": self.amplitude_exponent,
            "segment_exponent": self.segment_exponent,
            "predicted_quadrature": {
                "amplitude": QUADRATURE_AMPLITUDE_EXPONENT,
                "segments": QUADRATURE_SEGMENT_EXPONENT,
            },
            "predicted_circulation": {
                "amplitude": CIRCULATION_AMPLITUDE_EXPONENT,
                "segments": CIRCULATION_SEGMENT_EXPONENT,
            },
            "amplitude_distance_to_quadrature": self.amplitude_distance_to_quadrature,
            "amplitude_distance_to_circulation": self.amplitude_distance_to_circulation,
            "segment_distance_to_quadrature": self.segment_distance_to_quadrature,
            "segment_distance_to_circulation": self.segment_distance_to_circulation,
            "consistent_with_quadrature": self.consistent_with_quadrature,
            "consistent_with_circulation": self.consistent_with_circulation,
            "verdict": self.verdict,
        }


def loop_convergence(
    *,
    base: ClosedLoopWork,
    half_amplitude: ClosedLoopWork,
    refined_segments: ClosedLoopWork,
) -> LoopConvergence:
    """Identify a loop-work residual as quadrature error or as a genuine circulation, by its scaling.

    Args:
        base: The reference circuit, at amplitude ``a`` with ``N`` segments per edge.
        half_amplitude: The SAME circuit at amplitude ``a/2``, same ``N``.
        refined_segments: The SAME circuit at amplitude ``a``, with more segments per edge.

    Returns:
        The :class:`LoopConvergence`.

    Raises:
        ValueError: If the three circuits are not the intended refinements of one another — a scaling
            argument over circuits that differ in more than one respect measures nothing.
    """
    if half_amplitude.n_segments != base.n_segments:
        raise ValueError("the half-amplitude circuit must keep the base segment count")
    if not np.isclose(half_amplitude.amplitude_um, 0.5 * base.amplitude_um, rtol=1e-12):
        raise ValueError("the half-amplitude circuit must be at exactly half the base amplitude")
    if not np.isclose(refined_segments.amplitude_um, base.amplitude_um, rtol=1e-12):
        raise ValueError("the refined circuit must keep the base amplitude")
    if refined_segments.n_segments <= base.n_segments:
        raise ValueError("the refined circuit must use MORE segments than the base circuit")

    def exponent(coarse: float, fine: float, ratio: float) -> float | None:
        if coarse == 0.0 or fine == 0.0:
            return None
        return float(np.log(abs(coarse) / abs(fine)) / np.log(ratio))

    amplitude_exponent = exponent(base.work, half_amplitude.work, 2.0)
    segment_exponent = exponent(
        base.work, refined_segments.work, refined_segments.n_segments / base.n_segments)

    def distance(value: float | None, prediction: float) -> float | None:
        return None if value is None else abs(value - prediction)

    a_quad = distance(amplitude_exponent, QUADRATURE_AMPLITUDE_EXPONENT)
    a_circ = distance(amplitude_exponent, CIRCULATION_AMPLITUDE_EXPONENT)
    s_quad = distance(segment_exponent, QUADRATURE_SEGMENT_EXPONENT)
    s_circ = distance(segment_exponent, CIRCULATION_SEGMENT_EXPONENT)

    quadrature = bool(
        a_quad is not None and a_circ is not None and s_quad is not None and s_circ is not None
        and a_quad < a_circ and s_quad < s_circ
    )
    circulation = bool(
        a_quad is not None and a_circ is not None and s_quad is not None and s_circ is not None
        and a_circ < a_quad and s_circ < s_quad
    )
    verdict = "quadrature" if quadrature else ("circulation" if circulation else "ambiguous")
    return LoopConvergence(
        amplitude_exponent=amplitude_exponent,
        segment_exponent=segment_exponent,
        amplitude_distance_to_quadrature=a_quad,
        amplitude_distance_to_circulation=a_circ,
        segment_distance_to_quadrature=s_quad,
        segment_distance_to_circulation=s_circ,
        consistent_with_quadrature=quadrature,
        consistent_with_circulation=circulation,
        verdict=verdict,
    )


@dataclass(frozen=True, slots=True)
class EnergyBalance:
    """One accepted step's energy balance, as a residual against the largest single term.

    The balance is ``ΔU + dissipated + event_jump = active_input``.  Every term is supplied by the
    caller, in ``[pN·µm]``; this class only forms the residual and the dimensionless ratio, and it
    carries NO threshold — what counts as closed belongs to a gate contract, not to an observer.

    Attributes:
        delta_potential: ``ΔU`` over the step [pN·µm].
        dissipated: Work lost to drag over the step [pN·µm], sign-positive when energy leaves.
        event_jump: Energy discontinuity injected by accepted discrete events [pN·µm].  An ACTIVE
            system's total energy does not decrease, so monotonicity is not the gate — closure is.
        active_input: Work injected by the motors [pN·µm].
        residual: ``ΔU + dissipated + event_jump − active_input``.
        largest_term: ``max`` absolute term, the scale the residual is judged against.
        closure_ratio: ``|residual| / largest_term``, dimensionless; ``0.0`` if every term is zero.
    """

    delta_potential: float
    dissipated: float
    event_jump: float
    active_input: float
    residual: float
    largest_term: float
    closure_ratio: float

    @classmethod
    def from_terms(
        cls,
        *,
        delta_potential: float,
        dissipated: float,
        event_jump: float,
        active_input: float,
    ) -> EnergyBalance:
        """Form the balance from its four terms.

        Args:
            delta_potential: ``ΔU`` [pN·µm].
            dissipated: Drag dissipation [pN·µm].
            event_jump: Accepted-event energy jump [pN·µm].
            active_input: Motor work input [pN·µm].

        Returns:
            The :class:`EnergyBalance`.

        Raises:
            ValueError: If any term is not finite.
        """
        terms = (delta_potential, dissipated, event_jump, active_input)
        if not all(np.isfinite(term) for term in terms):
            raise ValueError(f"every energy term must be finite; got {terms!r}")
        residual = delta_potential + dissipated + event_jump - active_input
        largest = float(max(abs(term) for term in terms))
        return cls(
            delta_potential=float(delta_potential),
            dissipated=float(dissipated),
            event_jump=float(event_jump),
            active_input=float(active_input),
            residual=float(residual),
            largest_term=largest,
            closure_ratio=(abs(residual) / largest) if largest > 0.0 else 0.0,
        )


@dataclass(frozen=True, slots=True)
class StepEnergyLedger:
    """One accepted step's energy balance, with the evidence that each term means what it is called.

    Attributes:
        balance: The four-term :class:`EnergyBalance` and its dimensionless closure ratio.
        delta_potential_chord: ``ΔU`` recomputed along the straight chord between the step's endpoints,
            at the same quadrature effort.  Equal to :attr:`EnergyBalance.delta_potential` for a
            conservative passive field; this is the path-independence check that licenses calling the
            trajectory integral a potential difference.
        path_independence_ratio: ``|ΔU_path − ΔU_chord| / max(|ΔU_path|, |ΔU_chord|)``, dimensionless.
        exchanged_total: ``Σ|F_total·dx|`` along the trajectory [pN·µm] — the traffic the residual is
            small or large COMPARED TO, so nothing is judged against an invented reference.
        exchanged_passive: ``Σ|F_passive·dx|`` along the same walk [pN·µm].
        round_off_floor: The float64 accumulation bound for the walk's term count, times the traffic.
            A residual at this floor is arithmetic; a residual above it has a physical or numerical
            explanation, which is the point of reporting it rather than a verdict.
        at_round_off_floor: Whether ``|residual|`` sits at or below that floor.  NOT a closure verdict:
            midpoint quadrature over a subsampled trajectory has truncation error far above round-off
            even when the accounting is perfect, so what a run must show is that the residual FALLS with
            the sampling refinement — the same self-identifying argument the loop test uses.
        n_path_segments: Trajectory segments integrated (the subsampled polyline).
        path_length_um: ``Σ|Δx|`` along the integrated polyline [µm].
        squared_displacement_um2: The exact ``Σ|Δx|²`` [µm²] the dissipation was formed from, measured on
            every inner iteration rather than on the subsampled polyline.
        mobility_um_per_pN: The mobility the dissipation used [µm/pN].
    """

    balance: EnergyBalance
    delta_potential_chord: float
    path_independence_ratio: float
    exchanged_total: float
    exchanged_passive: float
    round_off_floor: float
    at_round_off_floor: bool
    n_path_segments: int
    path_length_um: float
    squared_displacement_um2: float
    mobility_um_per_pN: float

    def as_artifact_fields(self) -> dict[str, Any]:
        """Return the JSON-able record."""
        return {
            "delta_potential_pN_um": self.balance.delta_potential,
            "dissipated_pN_um": self.balance.dissipated,
            "event_jump_pN_um": self.balance.event_jump,
            "active_input_pN_um": self.balance.active_input,
            "residual_pN_um": self.balance.residual,
            "largest_term_pN_um": self.balance.largest_term,
            "closure_ratio": self.balance.closure_ratio,
            "delta_potential_chord_pN_um": self.delta_potential_chord,
            "path_independence_ratio": self.path_independence_ratio,
            "exchanged_total_pN_um": self.exchanged_total,
            "exchanged_passive_pN_um": self.exchanged_passive,
            "round_off_floor_pN_um": self.round_off_floor,
            "at_round_off_floor": self.at_round_off_floor,
            "n_path_segments": self.n_path_segments,
            "path_length_um": self.path_length_um,
            "squared_displacement_um2": self.squared_displacement_um2,
            "mobility_um_per_pN": self.mobility_um_per_pN,
        }


def step_energy_balance(
    *,
    path: np.ndarray,
    force_total: ForceFn,
    force_passive: ForceFn,
    squared_displacement_um2: float,
    mobility_um_per_pN: float,
    event_jump: float = 0.0,
) -> StepEnergyLedger:
    """Close ``ΔU + dissipated + event_jump = W_active`` over one step, from force launches alone.

    Every term comes from the launches the step already performs — see this module's header for why
    that, rather than a per-family potential, is what makes the balance a check instead of a second
    implementation of the same physics.

    Args:
        path: The step's own trajectory, shape ``(m, n, 3)`` [µm], subsampled from the inner iterations
            and INCLUDING both endpoints.  Refining the subsampling is how the residual is shown to be
            quadrature rather than a leak.
        force_total: The composed force with the binding state frozen, ``(n, 3) -> (n, 3)`` [pN].
        force_passive: The SAME evaluation with every active element disengaged (for a motor lane: every
            head detached, the binding SoA restored bit-for-bit afterwards).  Exactly one thing differs
            between the two callables.
        squared_displacement_um2: ``Σ|Δx|²`` [µm²] over EVERY inner iteration of the step, measured on
            device with :data:`accumulate_squared_displacement_kernel`.
        mobility_um_per_pN: The run's mobility step ``dt/γ`` [µm/pN].
        event_jump: Energy injected by accepted discrete events at fixed configuration [pN·µm].  Zero by
            derivation when the only binding-dependent force is the active one — pass it explicitly
            anyway on a lane where a passive bond can appear or vanish.

    Returns:
        The :class:`StepEnergyLedger`.

    Raises:
        ValueError: If the path is malformed, the mobility is not positive-finite, or a force evaluation
            returns a non-finite or mis-shaped value.
    """
    points = np.asarray(path, np.float64)
    if points.ndim != 3 or points.shape[2] != 3 or points.shape[0] < 2:
        raise ValueError(f"path must be (m>=2, n, 3); got shape {points.shape}")

    work_total, work_passive, exchanged_total, exchanged_passive = paired_path_work(
        force_total, force_passive, points)

    n_segments = int(points.shape[0] - 1)
    chord = np.stack([
        points[0] + (points[-1] - points[0]) * (index / n_segments)
        for index in range(n_segments + 1)
    ], axis=0)
    chord_work, _ = path_work(force_passive, chord)

    delta_potential = -work_passive
    active_input = work_total - work_passive
    dissipated = dissipated_work(squared_displacement_um2, mobility_um_per_pN)

    balance = EnergyBalance.from_terms(
        delta_potential=delta_potential, dissipated=dissipated,
        event_jump=float(event_jump), active_input=active_input,
    )

    chord_potential = -chord_work
    scale = max(abs(delta_potential), abs(chord_potential))
    floor = assemble_balance_tolerance_ratio(2 * n_segments) * (exchanged_total + exchanged_passive)
    return StepEnergyLedger(
        balance=balance,
        delta_potential_chord=float(chord_potential),
        path_independence_ratio=(
            float(abs(delta_potential - chord_potential) / scale) if scale > 0.0 else 0.0),
        exchanged_total=float(exchanged_total),
        exchanged_passive=float(exchanged_passive),
        round_off_floor=float(floor),
        at_round_off_floor=bool(abs(balance.residual) <= floor),
        n_path_segments=n_segments,
        path_length_um=float(np.sum(np.linalg.norm(np.diff(points, axis=0), axis=(1, 2)))),
        squared_displacement_um2=float(squared_displacement_um2),
        mobility_um_per_pN=float(mobility_um_per_pN),
    )


@dataclass(frozen=True, slots=True)
class BalanceConvergence:
    """Whether a non-zero energy-balance residual is the integrator's consistency error or a real leak.

    The same argument :class:`LoopConvergence` uses, one level up.  A residual has exactly two structural
    origins and they respond differently to the mobility ``μ`` the step was integrated at:

    ==================================  ==========================
    origin                              mobility ``μ``
    ==================================  ==========================
    first-order integrator consistency  ``μ¹``
    genuine leak / double count         independent of ``μ``
    ==================================  ==========================

    So halving the mobility (at the same relaxation, i.e. twice the iterations) halves an integrator
    residual and leaves a leak untouched.  The verdict compares the MEASURED exponent against those two
    PREDICTED ones and reports which is nearer — a comparison between structural predictions, not a
    threshold.

    Attributes:
        mobility_ratio: ``μ_base / μ_refined``, greater than one.
        residual_exponent: ``log(|R_base| / |R_refined|) / log(mobility_ratio)``, or ``None`` if either
            residual is exactly zero.
        distance_to_integrator_consistency: ``|residual_exponent − 1|``.
        distance_to_leak: ``|residual_exponent − 0|``.
        consistent_with_integrator: Whether the exponent sits nearer the first-order prediction.
        consistent_with_leak: Whether it sits nearer the leak prediction.
        verdict: ``"integrator_consistency"`` / ``"leak"`` / ``"ambiguous"``.
        base_closure_ratio: The coarser step's ``|residual| / largest_term``.
        refined_closure_ratio: The finer step's.
    """

    mobility_ratio: float
    residual_exponent: float | None
    distance_to_integrator_consistency: float | None
    distance_to_leak: float | None
    consistent_with_integrator: bool
    consistent_with_leak: bool
    verdict: str
    base_closure_ratio: float
    refined_closure_ratio: float

    def as_artifact_fields(self) -> dict[str, Any]:
        """Return the JSON-able record."""
        return {
            "mobility_ratio": self.mobility_ratio,
            "residual_exponent": self.residual_exponent,
            "predicted_integrator_consistency": INTEGRATOR_CONSISTENCY_MOBILITY_EXPONENT,
            "predicted_leak": LEAK_MOBILITY_EXPONENT,
            "distance_to_integrator_consistency": self.distance_to_integrator_consistency,
            "distance_to_leak": self.distance_to_leak,
            "consistent_with_integrator_consistency": self.consistent_with_integrator,
            "consistent_with_leak": self.consistent_with_leak,
            "verdict": self.verdict,
            "base_closure_ratio": self.base_closure_ratio,
            "refined_closure_ratio": self.refined_closure_ratio,
        }


def balance_convergence(
    *, base: StepEnergyLedger, refined: StepEnergyLedger
) -> BalanceConvergence:
    """Identify an energy-balance residual as integrator consistency or as a leak, by its own scaling.

    Args:
        base: The ledger of a step integrated at mobility ``μ``.
        refined: The ledger of the SAME relaxation integrated at a SMALLER mobility (more iterations, so
            the trajectory reaches the same place).  Only the mobility may differ; comparing two
            different relaxations measures the relaxation.

    Returns:
        The :class:`BalanceConvergence`.

    Raises:
        ValueError: If the refined ledger does not use a strictly smaller mobility.
    """
    if not (refined.mobility_um_per_pN < base.mobility_um_per_pN):
        raise ValueError(
            "the refined ledger must use a STRICTLY smaller mobility than the base ledger; "
            f"got {refined.mobility_um_per_pN!r} against {base.mobility_um_per_pN!r}"
        )
    ratio = float(base.mobility_um_per_pN / refined.mobility_um_per_pN)

    coarse = base.balance.residual
    fine = refined.balance.residual
    exponent: float | None = None
    if coarse != 0.0 and fine != 0.0:
        exponent = float(np.log(abs(coarse) / abs(fine)) / np.log(ratio))

    def distance(prediction: float) -> float | None:
        return None if exponent is None else abs(exponent - prediction)

    to_integrator = distance(INTEGRATOR_CONSISTENCY_MOBILITY_EXPONENT)
    to_leak = distance(LEAK_MOBILITY_EXPONENT)
    integrator = bool(to_integrator is not None and to_leak is not None and to_integrator < to_leak)
    leak = bool(to_integrator is not None and to_leak is not None and to_leak < to_integrator)
    return BalanceConvergence(
        mobility_ratio=ratio,
        residual_exponent=exponent,
        distance_to_integrator_consistency=to_integrator,
        distance_to_leak=to_leak,
        consistent_with_integrator=integrator,
        consistent_with_leak=leak,
        verdict=("integrator_consistency" if integrator else ("leak" if leak else "ambiguous")),
        base_closure_ratio=base.balance.closure_ratio,
        refined_closure_ratio=refined.balance.closure_ratio,
    )
