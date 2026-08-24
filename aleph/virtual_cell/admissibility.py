"""CPU sandbox oracle for accepted-step admissibility — the conjunction the engine does not make.

WHY THIS EXISTS.  The canonical whole-cell transaction
(:class:`aleph.engine.transaction.CellTransaction`) takes ONE device predicate and hands the same
flag to ``rollback``/``commit_irreversible`` and to the event clock.  When the caller does not supply
``accepted_d`` that predicate falls back to the global force-balance flag
(``transaction.py:154-161`` -> ``ledger.balance_ok_d``).  That flag is an **adjoint force-closure**
check: ``|Σreaction + Σtraction|² ≤ tol²`` on two independently accumulated halves of one Newton pair
(``ledger.py:242-258``).  The engine says so itself — ``ledger.py:419-421`` records that the two
resultants are equal and opposite *"whatever the configuration and however unconverged the solve.
Which is the point: the gate then tests the ADJOINT WIRING … not the convergence"*, and
``sf_motor_slice.py:65`` repeats it.  Newton's third law is therefore satisfied by construction on a
correctly wired slice, so the gate is blind to the inner mechanical residual, and an outer KMC/clock
step can commit on top of a mechanical solve that did not converge.

This module is the ORACLE for the missing conjunction.  It is an **unverified observer**: it imports
no engine runtime, imports no Warp, touches no gate, and can never emit
:attr:`~aleph.virtual_cell.contracts.EvidenceSource.NATIVE_ACCEPTED`.  Adopting its predicate inside
the engine is a PI gate-contract decision; see
``aleph/docs/v2_audit/PI_DECISION_CARD_ADMISSIBILITY_2026-07-29.md``.

ACCEPTANCE IS NOT STATIONARITY.  Nothing here asks a trajectory to settle.  ``inner_converged`` means
*this step's* mechanical solve reached the tolerance it was launched with; it does not mean forces are
going to zero, and a step inside a limit cycle, a NESS, a periodic or a chaotic regime is perfectly
admissible.  Regime classification is a different question and lives in another module.  No clause
here compares a state to its predecessor.

Sanity Gate:
  * dimensional: forces [pN]; ``inner_residual_pn`` is a force residual [pN] judged against the
    tolerance the solve was actually launched with, carried on the record rather than re-derived.
    Iteration counts and step indices are dimensionless; the clock digest covers a time [s].
  * boundary / numerical: adjoint closure is judged **scale-safely** — an absolute pN threshold is
    wrong across connectors whose force traffic differs by orders of magnitude.  The measure is
    ``|Σf_a + Σf_b| ≤ rel_tol · (Σ|f_a,i| + Σ|f_b,i|) + abs_floor``.  The scale is the sum of the
    contributions' MAGNITUDES, not the magnitude of their sum: a bipolar minifilament is a force
    dipole whose resultant is near zero while its traffic is not, so bracketing by the resultant would
    put the tolerance orders of magnitude below the arithmetic error it exists to admit.  This mirrors
    the engine's own D8 derivation (``ledger.py:215-235``).  The absolute floor covers the genuinely
    near-zero pair, where a relative test degenerates.  Neither number has a default: both are
    required caller arguments and choosing them is a PI gate decision.
  * numerical (NaN): a NaN comparison is False in both directions, so a naive ``residual > tol``
    REJECT-test silently accepts NaN.  Finiteness is therefore clause ZERO and short-circuits: when it
    fails, every downstream numeric clause is reported ``NOT_EVALUATED`` and no ordered comparison is
    performed at all.
  * conservation / measurement protocol: a rejected step must restore state, RNG, topology AND clock
    **bit-exactly**.  Verification is over the actual bytes (``ndarray.tobytes()``) via SHA-256, never
    a float comparison — ``a == b`` is False for two NaNs that are byte-identical and True for ``+0.0``
    vs ``-0.0``, and both of those are rollback bugs a float compare would misreport.  Restoring 3 of
    4 channels is a FAILURE and the audit names the channel that leaked.  A rollback check over a
    candidate that never mutated the channel is vacuous and is flagged as such.
  * ownership: this module owns no physical state, advances no clock, and adjudicates only records
    handed to it.  It is CPU-only and asserts nothing about the device.
"""

from __future__ import annotations

import hashlib
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
from types import MappingProxyType
from typing import Any

import numpy as np

from aleph.virtual_cell.contracts import EvidenceSource
from aleph.virtual_cell.step_receipt import predicate_vector_sha256

__all__ = [
    "ADMISSIBILITY_CLAUSE_ORDER",
    "EXPANSION_ROUTES",
    "AdmissibilityClause",
    "AdmissibilityVerdict",
    "AttemptedStepRecord",
    "ClauseOutcome",
    "ClauseStatus",
    "ConnectorForcePair",
    "RollbackAudit",
    "RollbackChannel",
    "StateDigestSet",
    "adjoint_closure_only",
    "digest_array",
    "digest_named_arrays",
    "evaluate_admissibility",
    "verify_rollback",
]

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_OBSERVER_AUTHORITY = "unverified-observer"


class AdmissibilityClause(StrEnum):
    """Named clauses of the conjunctive predicate, plus the transaction-integrity clause."""

    STATE_FINITE = "state_finite"
    ADJOINT_CLOSURE = "adjoint_closure"
    INNER_CONVERGED = "inner_converged"
    ITERATION_BUDGET = "iteration_budget"
    ROLLBACK_BIT_EXACT = "rollback_bit_exact"


class ClauseStatus(StrEnum):
    """Outcome of one clause; ``NOT_EVALUATED`` is not a pass and never counts toward acceptance."""

    PASS = "pass"
    FAIL = "fail"
    NOT_EVALUATED = "not-evaluated"


class RollbackChannel(StrEnum):
    """The four channels a rejected step must restore bit-exactly."""

    STATE = "state"
    RNG = "rng"
    TOPOLOGY = "topology"
    CLOCK = "clock"


#: The conjunction, in evaluation order.  ``STATE_FINITE`` is first and short-circuits the rest so a
#: NaN never reaches an ordered comparison.  ``ROLLBACK_BIT_EXACT`` is deliberately NOT in this tuple:
#: it audits the transaction's handling of a REJECTED step, not the step's admissibility.
ADMISSIBILITY_CLAUSE_ORDER: tuple[AdmissibilityClause, ...] = (
    AdmissibilityClause.STATE_FINITE,
    AdmissibilityClause.ADJOINT_CLOSURE,
    AdmissibilityClause.INNER_CONVERGED,
    AdmissibilityClause.ITERATION_BUDGET,
)

#: Where a failing clause routes when it becomes a programme wall, in the shared failure-routing
#: vocabulary.  A failure is never a reason to delete scope.
EXPANSION_ROUTES: Mapping[AdmissibilityClause, str] = MappingProxyType(
    {
        AdmissibilityClause.STATE_FINITE: (
            "native cost failure -> solver / preconditioner / adaptive design "
            "(a nonfinite candidate is a timestep / integrator / precision design problem)"
        ),
        AdmissibilityClause.ADJOINT_CLOSURE: (
            "archetype failure -> missing component / missing law project "
            "(a one-sided or sign-flipped scatter means a connector law is missing or miswired)"
        ),
        AdmissibilityClause.INNER_CONVERGED: (
            "native cost failure -> solver / preconditioner / adaptive design"
        ),
        AdmissibilityClause.ITERATION_BUDGET: (
            "native cost failure -> solver / preconditioner / adaptive design"
        ),
        AdmissibilityClause.ROLLBACK_BIT_EXACT: (
            "archetype failure -> missing component / missing law project "
            "(a leaked channel is state outside the transaction's coverage manifest)"
        ),
    }
)


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _digest(value: object, *, what: str) -> str:
    normalized = _nonempty(value, what=what)
    if not _SHA256.fullmatch(normalized):
        raise ValueError(f"{what} must be a lowercase SHA-256 digest")
    return normalized


def _nonnegative_int(value: object, *, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{what} must be a nonnegative integer")
    return int(value)


def _real_scalar(value: object, *, what: str) -> float:
    """Coerce to float WITHOUT rejecting NaN/inf — finiteness is a clause, not a constructor check."""
    if isinstance(value, bool) or not isinstance(value, (int, float, np.floating, np.integer)):
        raise TypeError(f"{what} must be a real scalar")
    return float(value)


def _readonly_float_array(values: Any, *, what: str) -> np.ndarray:
    array = np.array(values, dtype=np.float64, copy=True)
    if array.ndim == 1:
        array = array.reshape(1, -1)
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError(f"{what} must be a (3,) or (N, 3) force array [pN]")
    if array.shape[0] == 0:
        raise ValueError(f"{what} must not be empty")
    array.setflags(write=False)
    return array


def digest_array(values: Any) -> str:
    """SHA-256 over an array's dtype, shape and exact bytes.

    Bit-exact by construction: two arrays digest equal only if their bytes are identical, so ``+0.0``
    and ``-0.0`` differ (as they must for a rollback audit) and two identically-encoded NaNs match.

    Args:
        values: Anything :func:`numpy.ascontiguousarray` accepts and that is not an object array.

    Returns:
        The lowercase hex SHA-256 digest.

    Raises:
        TypeError: If the value resolves to a ``dtype=object`` array, whose bytes are pointers.
    """
    array = np.ascontiguousarray(values)
    if array.dtype == np.dtype(object):
        raise TypeError("object arrays have no byte identity; digest concrete numeric arrays")
    hasher = hashlib.sha256()
    hasher.update(array.dtype.str.encode())
    hasher.update(repr(array.shape).encode())
    hasher.update(array.tobytes())
    return hasher.hexdigest()


def digest_named_arrays(arrays: Mapping[str, Any]) -> str:
    """SHA-256 over a named channel of arrays, order-independent but name-sensitive.

    Args:
        arrays: Non-empty mapping of channel-field name to array-like value.

    Returns:
        The lowercase hex SHA-256 digest binding every field name to its byte digest.

    Raises:
        ValueError: If the mapping is empty or a name is blank or collides after normalization.
    """
    if not isinstance(arrays, Mapping) or not arrays:
        raise ValueError("named arrays must be a non-empty mapping")
    normalized: dict[str, str] = {}
    for raw_name, value in arrays.items():
        name = _nonempty(raw_name, what="array field name")
        if name in normalized:
            raise ValueError("array field names must not collide after whitespace normalization")
        normalized[name] = digest_array(value)
    hasher = hashlib.sha256()
    for name in sorted(normalized):
        hasher.update(name.encode())
        hasher.update(b"\x00")
        hasher.update(normalized[name].encode())
    return hasher.hexdigest()


@dataclass(frozen=True, slots=True)
class StateDigestSet:
    """Byte digests of the four channels an accepted-step transaction must own."""

    state: str
    rng: str
    topology: str
    clock: str

    def __post_init__(self) -> None:
        for field_name in ("state", "rng", "topology", "clock"):
            object.__setattr__(
                self, field_name, _digest(getattr(self, field_name), what=field_name)
            )

    @classmethod
    def from_channels(
        cls,
        *,
        state: Mapping[str, Any],
        rng: Mapping[str, Any],
        topology: Mapping[str, Any],
        clock: Mapping[str, Any],
    ) -> StateDigestSet:
        """Digest the four channels from their raw arrays.

        Args:
            state: Physical state arrays (positions, velocities, binding SoA, fields).
            rng: RNG stream arrays (counters, keys, epoch).
            topology: Connectivity arrays (edge lists, free lists, active-ID sets).
            clock: Clock arrays (accepted step index, accepted physical time [s]).

        Returns:
            The digest set for those bytes.
        """
        return cls(
            state=digest_named_arrays(state),
            rng=digest_named_arrays(rng),
            topology=digest_named_arrays(topology),
            clock=digest_named_arrays(clock),
        )

    def digest_for(self, channel: RollbackChannel) -> str:
        """Return the digest of one channel.

        Args:
            channel: Which of the four channels to read.

        Returns:
            The channel's lowercase hex SHA-256 digest.
        """
        return str(getattr(self, channel.value))


@dataclass(frozen=True, slots=True)
class ConnectorForcePair:
    """One connector's adjoint force pair — the two sides of a Newton's-third-law transfer.

    ``force_a_pn`` is what the connector scattered onto component A and ``force_b_pn`` what it
    scattered onto component B, each as the connector's OWN contribution (never a copy computed from
    the other side, which would make closure a tautology).  Both accept ``(3,)`` or ``(N, 3)`` [pN].
    """

    connector_id: str
    force_a_pn: np.ndarray
    force_b_pn: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "connector_id", _nonempty(self.connector_id, what="connector_id"))
        object.__setattr__(
            self, "force_a_pn", _readonly_float_array(self.force_a_pn, what="force_a_pn")
        )
        object.__setattr__(
            self, "force_b_pn", _readonly_float_array(self.force_b_pn, what="force_b_pn")
        )

    @property
    def finite(self) -> bool:
        """Whether every scattered force component on both sides is finite."""
        return bool(np.all(np.isfinite(self.force_a_pn)) and np.all(np.isfinite(self.force_b_pn)))

    @property
    def residual_pn(self) -> float:
        """``|Σf_a + Σf_b|`` [pN] — zero in exact arithmetic on a correctly wired connector."""
        resultant = self.force_a_pn.sum(axis=0) + self.force_b_pn.sum(axis=0)
        return float(np.linalg.norm(resultant))

    @property
    def magnitude_scale_pn(self) -> float:
        """``Σ|f_a,i| + Σ|f_b,i|`` [pN] — the connector's own force traffic, the Higham-bound scale."""
        magnitudes = float(np.linalg.norm(self.force_a_pn, axis=1).sum())
        magnitudes += float(np.linalg.norm(self.force_b_pn, axis=1).sum())
        return magnitudes


@dataclass(frozen=True, slots=True)
class AttemptedStepRecord:
    """One attempted physical step as an immutable, adjudicable record.

    Carries what an acceptance decision needs and nothing the decision may not see: the connectors'
    adjoint residuals, the inner mechanical residual AND the tolerance it was judged against, the
    iteration count against its budget, the candidate state for a finiteness check, the pre/post
    channel digests, and — when the transaction rolled back — the digests measured AFTER the rollback.

    Non-finite ``inner_residual_pn`` is accepted by the constructor on purpose: rejecting it here
    would hide the NaN-comparison defect this oracle exists to catch.
    """

    attempted_step_index: int
    attempted_dt_s: float
    connector_pairs: tuple[ConnectorForcePair, ...]
    inner_residual_pn: float
    inner_tolerance_pn: float
    inner_iterations: int
    inner_iteration_budget: int
    candidate_state: Mapping[str, np.ndarray]
    digests_before: StateDigestSet
    digests_after: StateDigestSet
    rollback_digests: StateDigestSet | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "attempted_step_index",
            _nonnegative_int(self.attempted_step_index, what="attempted_step_index"),
        )
        dt_s = _real_scalar(self.attempted_dt_s, what="attempted_dt_s")
        if not math.isfinite(dt_s) or dt_s <= 0.0:
            raise ValueError("attempted_dt_s must be finite and positive")
        object.__setattr__(self, "attempted_dt_s", dt_s)

        pairs = tuple(self.connector_pairs)
        if not pairs:
            raise ValueError("an attempted step must carry at least one connector force pair")
        if any(not isinstance(pair, ConnectorForcePair) for pair in pairs):
            raise TypeError("connector_pairs must contain ConnectorForcePair values")
        identifiers = [pair.connector_id for pair in pairs]
        if len(set(identifiers)) != len(identifiers):
            raise ValueError("connector ids must be unique within one attempted step")
        object.__setattr__(self, "connector_pairs", pairs)

        object.__setattr__(
            self,
            "inner_residual_pn",
            _real_scalar(self.inner_residual_pn, what="inner_residual_pn"),
        )
        tolerance = _real_scalar(self.inner_tolerance_pn, what="inner_tolerance_pn")
        if not math.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError("inner_tolerance_pn must be finite and positive")
        object.__setattr__(self, "inner_tolerance_pn", tolerance)

        object.__setattr__(
            self,
            "inner_iterations",
            _nonnegative_int(self.inner_iterations, what="inner_iterations"),
        )
        budget = _nonnegative_int(self.inner_iteration_budget, what="inner_iteration_budget")
        if budget == 0:
            raise ValueError("inner_iteration_budget must be positive")
        object.__setattr__(self, "inner_iteration_budget", budget)

        if not isinstance(self.candidate_state, Mapping) or not self.candidate_state:
            raise ValueError("candidate_state must be a non-empty mapping of arrays")
        frozen_state: dict[str, np.ndarray] = {}
        for raw_name, values in self.candidate_state.items():
            name = _nonempty(raw_name, what="candidate_state field")
            if name in frozen_state:
                raise ValueError("candidate_state fields must not collide after normalization")
            array = np.array(values, copy=True)
            if array.dtype == np.dtype(object):
                raise TypeError("candidate_state arrays must be numeric, not object arrays")
            array.setflags(write=False)
            frozen_state[name] = array
        object.__setattr__(self, "candidate_state", MappingProxyType(frozen_state))

        for field_name in ("digests_before", "digests_after"):
            if not isinstance(getattr(self, field_name), StateDigestSet):
                raise TypeError(f"{field_name} must be a StateDigestSet")
        if self.rollback_digests is not None and not isinstance(
            self.rollback_digests, StateDigestSet
        ):
            raise TypeError("rollback_digests must be a StateDigestSet or None")

    @property
    def state_is_finite(self) -> bool:
        """Whether every candidate-state array holds only finite values."""
        for array in self.candidate_state.values():
            if not np.issubdtype(array.dtype, np.number):
                continue
            if np.issubdtype(array.dtype, np.integer):
                continue
            if not bool(np.all(np.isfinite(array))):
                return False
        return True

    @property
    def nonfinite_sources(self) -> tuple[str, ...]:
        """Names of every nonfinite source found, in stable order (state fields, then residual)."""
        sources: list[str] = []
        for name in sorted(self.candidate_state):
            array = self.candidate_state[name]
            if np.issubdtype(array.dtype, np.floating) and not bool(np.all(np.isfinite(array))):
                sources.append(f"candidate_state.{name}")
        for pair in self.connector_pairs:
            if not pair.finite:
                sources.append(f"connector.{pair.connector_id}")
        if not math.isfinite(self.inner_residual_pn):
            sources.append("inner_residual_pn")
        return tuple(sources)


@dataclass(frozen=True, slots=True)
class ClauseOutcome:
    """One clause's result with the number it measured and the threshold it was judged against."""

    clause: AdmissibilityClause
    status: ClauseStatus
    detail: str
    measured: float | None = None
    threshold: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.clause, AdmissibilityClause):
            raise TypeError("clause must be an AdmissibilityClause")
        if not isinstance(self.status, ClauseStatus):
            raise TypeError("status must be a ClauseStatus")
        object.__setattr__(self, "detail", _nonempty(self.detail, what="clause detail"))

    @property
    def passed(self) -> bool:
        """True only for :attr:`ClauseStatus.PASS`; ``NOT_EVALUATED`` is never a pass."""
        return self.status is ClauseStatus.PASS


@dataclass(frozen=True, slots=True)
class RollbackAudit:
    """Bit-exactness audit of a rejected step's restoration, per channel."""

    restored: tuple[RollbackChannel, ...]
    leaked: tuple[RollbackChannel, ...]
    mutated_by_candidate: tuple[RollbackChannel, ...]

    def __post_init__(self) -> None:
        for field_name in ("restored", "leaked", "mutated_by_candidate"):
            channels = tuple(getattr(self, field_name))
            if any(not isinstance(channel, RollbackChannel) for channel in channels):
                raise TypeError(f"{field_name} must contain RollbackChannel values")
            object.__setattr__(self, field_name, channels)
        if set(self.restored) & set(self.leaked):
            raise ValueError("a channel cannot be both restored and leaked")
        if set(self.restored) | set(self.leaked) != set(RollbackChannel):
            raise ValueError("a rollback audit must classify every channel exactly once")

    @property
    def bit_exact(self) -> bool:
        """True only when all four channels were restored bit-exactly."""
        return not self.leaked

    @property
    def vacuous(self) -> bool:
        """True when the candidate mutated nothing, so restoration proves nothing."""
        return not self.mutated_by_candidate

    @property
    def detail(self) -> str:
        """Human-readable summary naming the channels that leaked."""
        if self.leaked:
            names = ", ".join(channel.value for channel in self.leaked)
            return f"rollback leaked {len(self.leaked)}/4 channel(s): {names}"
        if self.vacuous:
            return "rollback bit-exact on all 4 channels, but the candidate mutated none (vacuous)"
        mutated = ", ".join(channel.value for channel in self.mutated_by_candidate)
        return f"rollback bit-exact on all 4 channels (candidate had mutated: {mutated})"


@dataclass(frozen=True, slots=True)
class AdmissibilityVerdict:
    """Typed, per-clause verdict for one attempted step.  Never a bare bool.

    ``admissible`` is the conjunction of :data:`ADMISSIBILITY_CLAUSE_ORDER`.  ``transaction_sound``
    additionally requires that a rejected step was rolled back bit-exactly on all four channels.
    """

    attempted_step_index: int
    admissible: bool
    clauses: tuple[ClauseOutcome, ...]
    rejected_by: tuple[AdmissibilityClause, ...]
    rollback: RollbackAudit | None
    expansions: tuple[str, ...]
    authority: str = _OBSERVER_AUTHORITY

    def __post_init__(self) -> None:
        if self.authority != _OBSERVER_AUTHORITY:
            raise ValueError("the admissibility oracle authority is fixed to unverified-observer")
        clauses = tuple(self.clauses)
        if any(not isinstance(clause, ClauseOutcome) for clause in clauses):
            raise TypeError("clauses must contain ClauseOutcome values")
        object.__setattr__(self, "clauses", clauses)
        object.__setattr__(self, "rejected_by", tuple(self.rejected_by))
        object.__setattr__(self, "expansions", tuple(self.expansions))
        if self.admissible and self.rejected_by:
            raise ValueError("an admissible verdict cannot name a rejecting clause")

    @property
    def evidence_source(self) -> EvidenceSource:
        """This oracle produces analytic-oracle evidence and can never produce native evidence."""
        return EvidenceSource.ANALYTIC_ORACLE

    @property
    def transaction_sound(self) -> bool:
        """Admissible, or rejected AND rolled back bit-exactly on all four channels."""
        if self.admissible:
            return True
        return self.rollback is not None and self.rollback.bit_exact

    def outcome(self, clause: AdmissibilityClause) -> ClauseOutcome:
        """Return one clause's outcome.

        Args:
            clause: The clause to look up.

        Returns:
            The recorded :class:`ClauseOutcome`.

        Raises:
            KeyError: If the clause was not part of this adjudication.
        """
        for recorded in self.clauses:
            if recorded.clause is clause:
                return recorded
        raise KeyError(f"clause {clause.value} was not adjudicated")

    def predicate_results(self) -> Mapping[str, bool]:
        """Return the clause vector in the shape a ``StepReceipt`` binds.

        ``NOT_EVALUATED`` maps to ``False``: a clause that was skipped never granted acceptance.

        Returns:
            A read-only mapping of clause name to boolean pass.
        """
        return MappingProxyType({c.clause.value: c.passed for c in self.clauses})

    def predicate_results_sha256(self) -> str:
        """Bind the clause vector with the receipt layer's own hash.

        Returns:
            The lowercase hex SHA-256 digest of the normalized clause vector.
        """
        return predicate_vector_sha256(dict(self.predicate_results()))


def _closure_tolerance_pn(scale_pn: float, *, rel_tol: float, abs_floor_pn: float) -> float:
    return rel_tol * scale_pn + abs_floor_pn


def _validate_tolerances(rel_tol: object, abs_floor_pn: object) -> tuple[float, float]:
    relative = _real_scalar(rel_tol, what="adjoint_rel_tol")
    floor = _real_scalar(abs_floor_pn, what="adjoint_abs_floor_pn")
    if not math.isfinite(relative) or relative < 0.0 or relative >= 1.0:
        raise ValueError("adjoint_rel_tol must be finite and in [0, 1)")
    if not math.isfinite(floor) or floor < 0.0:
        raise ValueError("adjoint_abs_floor_pn must be finite and nonnegative")
    if relative == 0.0 and floor == 0.0:
        raise ValueError("adjoint closure needs a nonzero relative tolerance or absolute floor")
    return relative, floor


def adjoint_closure_only(
    record: AttemptedStepRecord,
    *,
    adjoint_rel_tol: float,
    adjoint_abs_floor_pn: float,
) -> bool:
    """Model the INCUMBENT force-only acceptance predicate: global adjoint closure, nothing else.

    This reproduces what ``GlobalCellLedger.assemble_balance`` computes (``ledger.py:242-258``): the
    two channels are accumulated GLOBALLY — every connector's A-side into the reaction resultant and
    every B-side into the traction resultant — and the gate compares ``|Σreaction + Σtraction|``
    against a tolerance derived from the ledger's own accumulated magnitudes.  It exists here so a
    test can demonstrate, on the same record, what the incumbent accepts and what
    :func:`evaluate_admissibility` rejects.  It is a MODEL of the engine predicate, not a call into it.

    Two weaknesses are intentional and faithful to the incumbent.  It is blind to the inner mechanical
    residual entirely.  And because the accumulation is global, two connectors with equal and opposite
    one-sided errors cancel and the gate passes — a per-connector test is strictly stronger.

    Args:
        record: The attempted step to judge.
        adjoint_rel_tol: Relative closure tolerance (dimensionless, in ``[0, 1)``).  No default: this
            is a PI gate decision.
        adjoint_abs_floor_pn: Absolute closure floor [pN] for the near-zero case.  No default.

    Returns:
        Whether the incumbent force-only gate would accept.  A nonfinite residual compares False and
        therefore rejects — but only if the nonfinite value reaches the resultant; a NaN elsewhere in
        the candidate state is invisible to this gate.

    Raises:
        TypeError: If ``record`` is not an :class:`AttemptedStepRecord`.
        ValueError: If the tolerances are out of range or both zero.
    """
    if not isinstance(record, AttemptedStepRecord):
        raise TypeError("record must be an AttemptedStepRecord")
    relative, floor = _validate_tolerances(adjoint_rel_tol, adjoint_abs_floor_pn)
    reaction = np.zeros(3, dtype=np.float64)
    traction = np.zeros(3, dtype=np.float64)
    scale_pn = 0.0
    for pair in record.connector_pairs:
        reaction = reaction + pair.force_a_pn.sum(axis=0)
        traction = traction + pair.force_b_pn.sum(axis=0)
        scale_pn += pair.magnitude_scale_pn
    residual_pn = float(np.linalg.norm(reaction + traction))
    tolerance_pn = _closure_tolerance_pn(scale_pn, rel_tol=relative, abs_floor_pn=floor)
    return bool(residual_pn <= tolerance_pn)


def verify_rollback(record: AttemptedStepRecord) -> RollbackAudit:
    """Verify bit-exact restoration of state, RNG, topology and clock after a rejected step.

    Comparison is over SHA-256 digests of the raw bytes, never over floats.  A channel counts as
    restored only when its post-rollback digest equals its pre-step digest exactly; restoring three of
    four is a failure and the returned audit names the channel that leaked.

    The audit also records which channels the candidate had actually mutated before rollback: if none
    did, the restoration check is vacuous and :attr:`RollbackAudit.vacuous` says so.

    Args:
        record: An attempted step carrying ``rollback_digests``.

    Returns:
        The per-channel :class:`RollbackAudit`.

    Raises:
        ValueError: If the record carries no post-rollback digests.
    """
    if record.rollback_digests is None:
        raise ValueError("rollback verification needs post-rollback digests on the record")
    restored: list[RollbackChannel] = []
    leaked: list[RollbackChannel] = []
    mutated: list[RollbackChannel] = []
    for channel in RollbackChannel:
        before = record.digests_before.digest_for(channel)
        after = record.digests_after.digest_for(channel)
        rolled_back = record.rollback_digests.digest_for(channel)
        if after != before:
            mutated.append(channel)
        if rolled_back == before:
            restored.append(channel)
        else:
            leaked.append(channel)
    return RollbackAudit(
        restored=tuple(restored),
        leaked=tuple(leaked),
        mutated_by_candidate=tuple(mutated),
    )


def _finiteness_clause(record: AttemptedStepRecord) -> ClauseOutcome:
    sources = record.nonfinite_sources
    if sources:
        return ClauseOutcome(
            clause=AdmissibilityClause.STATE_FINITE,
            status=ClauseStatus.FAIL,
            detail=f"nonfinite value in: {', '.join(sources)}",
        )
    return ClauseOutcome(
        clause=AdmissibilityClause.STATE_FINITE,
        status=ClauseStatus.PASS,
        detail="candidate state, connector forces and inner residual are all finite",
    )


def _adjoint_clause(
    record: AttemptedStepRecord, *, rel_tol: float, abs_floor_pn: float
) -> ClauseOutcome:
    worst_pair: ConnectorForcePair | None = None
    worst_excess = -math.inf
    worst_tolerance = 0.0
    for pair in record.connector_pairs:
        tolerance_pn = _closure_tolerance_pn(
            pair.magnitude_scale_pn, rel_tol=rel_tol, abs_floor_pn=abs_floor_pn
        )
        excess = pair.residual_pn - tolerance_pn
        if excess > worst_excess:
            worst_excess = excess
            worst_pair = pair
            worst_tolerance = tolerance_pn
    assert worst_pair is not None  # a record always carries at least one pair
    if worst_excess > 0.0:
        return ClauseOutcome(
            clause=AdmissibilityClause.ADJOINT_CLOSURE,
            status=ClauseStatus.FAIL,
            detail=(
                f"connector {worst_pair.connector_id!r} force pair does not close: "
                f"|Σf_a + Σf_b| = {worst_pair.residual_pn:.6g} pN over a traffic scale of "
                f"{worst_pair.magnitude_scale_pn:.6g} pN"
            ),
            measured=worst_pair.residual_pn,
            threshold=worst_tolerance,
        )
    return ClauseOutcome(
        clause=AdmissibilityClause.ADJOINT_CLOSURE,
        status=ClauseStatus.PASS,
        detail=(
            f"every connector force pair closes; worst is {worst_pair.connector_id!r} at "
            f"{worst_pair.residual_pn:.6g} pN"
        ),
        measured=worst_pair.residual_pn,
        threshold=worst_tolerance,
    )


def _convergence_clause(record: AttemptedStepRecord) -> ClauseOutcome:
    converged = record.inner_residual_pn <= record.inner_tolerance_pn
    return ClauseOutcome(
        clause=AdmissibilityClause.INNER_CONVERGED,
        status=ClauseStatus.PASS if converged else ClauseStatus.FAIL,
        detail=(
            f"inner mechanical residual {record.inner_residual_pn:.6g} pN vs the tolerance "
            f"{record.inner_tolerance_pn:.6g} pN the solve was judged against"
        ),
        measured=record.inner_residual_pn,
        threshold=record.inner_tolerance_pn,
    )


def _budget_clause(record: AttemptedStepRecord) -> ClauseOutcome:
    within = record.inner_iterations <= record.inner_iteration_budget
    return ClauseOutcome(
        clause=AdmissibilityClause.ITERATION_BUDGET,
        status=ClauseStatus.PASS if within else ClauseStatus.FAIL,
        detail=(
            f"inner solve used {record.inner_iterations} of "
            f"{record.inner_iteration_budget} permitted iterations"
        ),
        measured=float(record.inner_iterations),
        threshold=float(record.inner_iteration_budget),
    )


def _skipped(clause: AdmissibilityClause) -> ClauseOutcome:
    return ClauseOutcome(
        clause=clause,
        status=ClauseStatus.NOT_EVALUATED,
        detail=(
            "not evaluated: the finiteness clause failed first, and an ordered comparison against a "
            "nonfinite value is False in both directions"
        ),
    )


def evaluate_admissibility(
    record: AttemptedStepRecord,
    *,
    adjoint_rel_tol: float,
    adjoint_abs_floor_pn: float,
) -> AdmissibilityVerdict:
    """Adjudicate one attempted step under the CONJUNCTIVE admissibility predicate.

    The conjunction is, in order: finiteness, per-connector adjoint closure, inner mechanical
    convergence against the tolerance the solve was launched with, and iteration count within budget.
    Finiteness is evaluated FIRST and short-circuits — every later clause is reported
    ``NOT_EVALUATED`` rather than compared, because ``nan <= tol`` and ``nan > tol`` are both False and
    a naive predicate therefore accepts NaN.

    Adjoint closure is judged per connector and scale-safely; see the module Sanity Gate for why an
    absolute threshold is wrong and why the scale is the sum of contribution magnitudes rather than
    the magnitude of their sum.

    This is an oracle.  A verdict here changes no engine behaviour, and its ``authority`` is fixed to
    ``unverified-observer``.

    Args:
        record: The attempted step to adjudicate.
        adjoint_rel_tol: Relative adjoint-closure tolerance (dimensionless, in ``[0, 1)``).  There is
            deliberately NO default — choosing it is a PI gate-contract decision, and a value picked
            to make a step pass would violate the no-tuning rule.
        adjoint_abs_floor_pn: Absolute adjoint-closure floor [pN] for pairs whose traffic is near
            zero, where a purely relative test degenerates.  Also no default, also a PI decision.

    Returns:
        A typed :class:`AdmissibilityVerdict` with a per-clause breakdown, the clauses that rejected,
        the rollback audit when the record carries post-rollback digests, and the expansions each
        failure opens.

    Raises:
        TypeError: If ``record`` is not an :class:`AttemptedStepRecord`.
        ValueError: If the tolerances are out of range or both zero.
    """
    if not isinstance(record, AttemptedStepRecord):
        raise TypeError("record must be an AttemptedStepRecord")
    relative, floor = _validate_tolerances(adjoint_rel_tol, adjoint_abs_floor_pn)

    finiteness = _finiteness_clause(record)
    if finiteness.passed:
        clauses = (
            finiteness,
            _adjoint_clause(record, rel_tol=relative, abs_floor_pn=floor),
            _convergence_clause(record),
            _budget_clause(record),
        )
    else:
        clauses = (
            finiteness,
            _skipped(AdmissibilityClause.ADJOINT_CLOSURE),
            _skipped(AdmissibilityClause.INNER_CONVERGED),
            _skipped(AdmissibilityClause.ITERATION_BUDGET),
        )

    rejected = tuple(clause.clause for clause in clauses if not clause.passed)
    admissible = not rejected

    rollback: RollbackAudit | None = None
    if record.rollback_digests is not None:
        rollback = verify_rollback(record)
        clauses = clauses + (
            ClauseOutcome(
                clause=AdmissibilityClause.ROLLBACK_BIT_EXACT,
                status=ClauseStatus.PASS if rollback.bit_exact else ClauseStatus.FAIL,
                detail=rollback.detail,
                measured=float(len(rollback.restored)),
                threshold=float(len(RollbackChannel)),
            ),
        )

    expansion_keys = list(rejected)
    if rollback is not None and not rollback.bit_exact:
        expansion_keys.append(AdmissibilityClause.ROLLBACK_BIT_EXACT)
    expansions = tuple(dict.fromkeys(EXPANSION_ROUTES[key] for key in expansion_keys))

    return AdmissibilityVerdict(
        attempted_step_index=record.attempted_step_index,
        admissible=admissible,
        clauses=clauses,
        rejected_by=rejected,
        rollback=rollback,
        expansions=expansions,
    )
