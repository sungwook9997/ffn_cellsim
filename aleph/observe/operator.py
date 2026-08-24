"""The observation seam: accepted state in, an *apparent* quantity out, with the pipeline visible.

The shape this module enforces
------------------------------
An experiment does not hand anyone a modulus. It hands over a voltage; a calibration turns that into
a deflection; a contact model turns *that* into a modulus. Each step is a modelling choice, and the
number at the end is **apparent** — it is what this analysis pipeline reports, not what the material
is.

A simulator that computes the true quantity directly and compares it against a published one has
skipped every step where the disagreement actually lives. It will then attribute the mismatch to
physics. So :class:`ObservationOperator` keeps the two halves apart by type:

    raw_observable(state, context) -> RawObservable  | Refusal
    apparent_quantity(raw)         -> ApparentQuantity | Refusal

and :class:`ApparentQuantity` carries ``analysis_chain``, the ordered list of what was done to the
raw observable to get it. A quantity that cannot say how it was produced cannot be compared to one
that can.

Refusal is a return value, not an exception and not a NaN
---------------------------------------------------------
:class:`Refusal` is a frozen dataclass and is **not** a subclass of :class:`BaseException`. A test
asserts this so it can never quietly become one. Three failure modes are being ruled out at once:

1. **Exception-as-control-flow.** A caller who forgets the ``try`` gets a crash in production; one
   who remembers writes a broad ``except`` that swallows genuine bugs alongside the refusal.
2. **A silent NaN.** It propagates, surfaces far from its origin, and carries no reason. Worse, a
   NaN that lands in a reduction turns a whole array into a NaN, so the eventual complaint is about
   the wrong quantity.
3. **``None``.** Indistinguishable from "not implemented yet".

A refusal carries a typed :class:`RefusalCode`, a sentence a human can act on, and a
machine-readable ``detail`` mapping. It is a first-class result: refusing to answer is a scientific
finding, and it is the finding this project most wants to be able to record.

Defensive posture toward other lanes
------------------------------------
Nothing here imports :mod:`aleph.state` or :mod:`aleph.units`. The state is reached through the
:class:`AcceptedState` Protocol — owner names and arrays — so this lane cannot be broken by another
lane's API settling, and neither lane has to wait for the other. Units are carried as documented
strings on the results; when L1's dimension algebra settles those fields are the attachment point.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any, Protocol, runtime_checkable

import numpy as np

from aleph.observe.manifest import ProtocolManifest

__all__ = [
    "AcceptedState",
    "ApparentQuantity",
    "ApplicabilityDomain",
    "ObservationContext",
    "ObservationOperator",
    "RawObservable",
    "Refusal",
    "RefusalCode",
    "is_refusal",
]


class RefusalCode(StrEnum):
    """Why an operator declined to produce a number.

    Each member is a *different finding*, and they are kept apart for the same reason
    :class:`~aleph.observe.manifest.Declaration` keeps its two members apart: collapsing them loses
    the information that decides what to do next.
    """

    #: A required input quantity is outside the range over which this operator is claimed to work.
    #: The operator may well return a number there; it just would not be entitled to.
    OUT_OF_APPLICABILITY_DOMAIN = "OUT_OF_APPLICABILITY_DOMAIN"

    #: The state does not carry an owner this operator needs. A structural mismatch, not a data one.
    MISSING_STATE_OWNER = "MISSING_STATE_OWNER"

    #: There is data, and there is not enough of it to support any claim. Distinct from
    #: ``OUT_OF_APPLICABILITY_DOMAIN``: "we looked and could not tell" versus "we should not look".
    INSUFFICIENT_SAMPLES = "INSUFFICIENT_SAMPLES"

    #: A NaN or an infinity reached the operator. Never absorbed, never propagated.
    NON_FINITE_INPUT = "NON_FINITE_INPUT"

    #: The question is well posed but beyond what the discretisation can resolve — e.g. a spectrum
    #: requested past the mesh Nyquist limit. Answering would return an alias wearing the label of
    #: the thing that was asked for.
    BEYOND_RESOLUTION_LIMIT = "BEYOND_RESOLUTION_LIMIT"

    #: This operator does not report that quantity. A refusal, not a failure.
    UNSUPPORTED_QUANTITY = "UNSUPPORTED_QUANTITY"

    #: Inputs are mutually inconsistent — shapes, lengths, or units that cannot both be right.
    INCONSISTENT_INPUT = "INCONSISTENT_INPUT"

    #: A **known, open defect** on the path this operator reads makes its input untrustworthy. The
    #: arithmetic would succeed and the number would be wrong in a way that looks right.
    #:
    #: Kept apart from ``OUT_OF_APPLICABILITY_DOMAIN`` because the two decide different next
    #: actions. Out-of-domain says *ask a different question, or move the input*. This one says
    #: *nothing you can do to the input will help; another lane owns the fix, and when it lands
    #: this operator will start answering without being edited*. Collapsing them would put a
    #: blocked-on-someone-else finding in the column a caller reads as their own mistake.
    #:
    #: A refusal carrying this code must name, in ``detail``, the precondition that failed and the
    #: lane that owns it. An unroutable one is a log line.
    UPSTREAM_DEFECT_OPEN = "UPSTREAM_DEFECT_OPEN"


@dataclass(frozen=True, slots=True)
class Refusal:
    """A typed decline. Deliberately **not** an exception and never a NaN.

    Attributes:
        code: The machine-readable reason.
        reason: One sentence a human can act on. Should name the offending quantity and its value.
        operator_id: Which operator declined.
        detail: Machine-readable specifics — the value, the bound it violated, the owner that was
            missing. A refusal a caller cannot programmatically inspect is a log line.
    """

    code: RefusalCode
    reason: str
    operator_id: str
    detail: Mapping[str, Any] = field(default_factory=dict)

    def __str__(self) -> str:
        return f"[{self.code.value}] {self.operator_id}: {self.reason}"

    def as_dict(self) -> dict[str, Any]:
        """JSON-able form, for a run record."""
        return {
            "refused": True,
            "code": self.code.value,
            "reason": self.reason,
            "operator_id": self.operator_id,
            "detail": dict(self.detail),
        }


def is_refusal(value: object) -> bool:
    """``True`` when a result is a :class:`Refusal`.

    A one-line helper exists so that call sites read as a branch on the result rather than as an
    ``isinstance`` check someone might be tempted to replace with a truthiness test. A ``Refusal``
    is truthy.
    """
    return isinstance(value, Refusal)


@dataclass(frozen=True, slots=True)
class ApplicabilityDomain:
    """Where an operator is entitled to report a number.

    An operator will usually *return* something outside its domain — floating point is obliging.
    The domain is the statement of where that something means anything, and it must be declared
    before the measurement, not inferred from where the answer looked reasonable.

    Attributes:
        description: What the domain is, in words, for a reader of a run record.
        bounds: ``{quantity_name: (low, high)}``, closed intervals, in the operator's declared
            units. Use ``-inf``/``inf`` for a one-sided bound.
        required_state_owners: Owner names the operator needs present in the accepted state.
    """

    description: str
    bounds: Mapping[str, tuple[float, float]] = field(default_factory=dict)
    required_state_owners: tuple[str, ...] = ()

    def check(self, operator_id: str, **values: float) -> Refusal | None:
        """Check named quantities against the declared bounds.

        A quantity with no declared bound is *ignored*, not rejected: the domain says where the
        operator is entitled to speak, and silence about a quantity is not a constraint on it.
        Passing a non-finite value is always a refusal, whether or not it is bounded.

        Args:
            operator_id: Recorded on any refusal produced.
            **values: Quantities to check, by the names used in :attr:`bounds`.

        Returns:
            A :class:`Refusal` for the first violation found, or ``None`` if every checked quantity
            is inside its bounds.
        """
        for name, raw in values.items():
            value = float(raw)
            if not math.isfinite(value):
                return Refusal(
                    code=RefusalCode.NON_FINITE_INPUT,
                    reason=f"{name} is {value!r}, which is not a finite number",
                    operator_id=operator_id,
                    detail={"quantity": name, "value": value},
                )
            bound = self.bounds.get(name)
            if bound is None:
                continue
            low, high = float(bound[0]), float(bound[1])
            if not (low <= value <= high):
                return Refusal(
                    code=RefusalCode.OUT_OF_APPLICABILITY_DOMAIN,
                    reason=(
                        f"{name} = {value:.6g} is outside the declared applicability domain "
                        f"[{low:.6g}, {high:.6g}]. The operator would return a number here and "
                        "would not be entitled to it."
                    ),
                    operator_id=operator_id,
                    detail={
                        "quantity": name,
                        "value": value,
                        "low": low,
                        "high": high,
                        "domain": self.description,
                    },
                )
        return None

    def check_owners(self, operator_id: str, available: frozenset[str]) -> Refusal | None:
        """Refuse before any arithmetic if a required state owner is absent.

        This check runs first, deliberately. An operator that computes for a while and *then*
        discovers the state was wrong has already spent effort deciding what to report, and the
        temptation at that point is to report the part that worked.

        Args:
            operator_id: Recorded on any refusal produced.
            available: Owner names present in the accepted state.

        Returns:
            A :class:`Refusal` naming every missing owner, or ``None``.
        """
        missing = tuple(sorted(set(self.required_state_owners) - set(available)))
        if missing:
            return Refusal(
                code=RefusalCode.MISSING_STATE_OWNER,
                reason=(
                    f"the accepted state does not carry required owner(s) {list(missing)}; "
                    f"present: {sorted(available)}"
                ),
                operator_id=operator_id,
                detail={"missing": list(missing), "available": sorted(available)},
            )
        return None


@runtime_checkable
class AcceptedState(Protocol):
    """The minimum an operator needs from a committed simulation state.

    Deliberately tiny, and deliberately not :mod:`aleph.state`. An observer needs to know which
    owners exist and to read their arrays; it has no business knowing how they are stored, and
    coupling it to L2's evolving schema would make the two lanes block each other.

    Implementations must expose:

    ``owners``
        Owner names present. An owner is a subsystem that holds state — ``"membrane"``,
        ``"cortex"``.

    ``array(owner, name)``
        A read-only view of one named array belonging to one owner.
    """

    @property
    def owners(self) -> frozenset[str]: ...

    def array(self, owner: str, name: str) -> np.ndarray: ...


@dataclass(frozen=True, slots=True)
class ObservationContext:
    """What the observer is told about the run, separate from the state itself.

    Kept apart from :class:`AcceptedState` because these are facts about the *measurement*, not
    about the system, and conflating them is how a sampling cadence ends up inferred from a solver
    step size.

    Attributes:
        sample_interval_s: Physical time between consecutive recorded frames [s]. This is the
            *sampling* interval, not the integrator step: with a recording stride the two differ,
            and handing the integrator step to a correlation-time estimator understates ``tau_int``
            by exactly that stride — an error in the flattering direction.
        temperature_k: Temperature [K].
        extras: Anything else a specific operator declares it needs. Named so that a reader can see
            it was passed in rather than assumed.
    """

    sample_interval_s: float
    temperature_k: float
    extras: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class RawObservable:
    """What the instrument produces, before any interpretation.

    This is the analogue of the voltage trace: the quantity the measurement literally yields, with
    its own units and its own sample accounting. It is *not* the answer, and keeping it as a
    separate type is what stops it being mistaken for one.

    Attributes:
        name: What was measured.
        operator_id: Which operator produced it.
        values: The raw array, in :attr:`units`.
        units: Units of ``values``, as text — e.g. ``"um^2"``.
        n_samples: Frames that entered each entry of ``values``. Same shape as ``values``, or a
            scalar count.
        manifest_hash: Digest of the protocol manifest in force. Travels with the data from here on.
        support: Coordinates the values are indexed by, e.g. ``{"q_inv_um": array}``. Empty for a
            scalar-support observable.
        detail: Anything else the operator wants on the record.
    """

    name: str
    operator_id: str
    values: np.ndarray
    units: str
    n_samples: np.ndarray | int
    manifest_hash: str
    support: Mapping[str, np.ndarray] = field(default_factory=dict, repr=False)
    detail: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class ApparentQuantity:
    """A number the analysis pipeline reports — *apparent*, and saying so in its own name.

    The word carries the load. A modulus recovered by fitting a Hertz model to a force curve is the
    modulus *of the Hertz model that best fit that curve*, under that contact assumption, over that
    fit window. Calling it "the modulus" asserts something the measurement never established.

    Attributes:
        name: The quantity.
        value: Its value, in :attr:`units`.
        units: Units, as text.
        standard_error: Error on ``value``, or ``None`` when the pipeline cannot support one. Never
            zero as a stand-in for "not computed".
        n_effective: Independent samples behind ``value``. Not the frame count — see
            :mod:`aleph.observe.stationarity` for why those differ by ``sqrt(2*tau_int/dt)``.
        manifest_hash: Digest of the protocol manifest in force.
        analysis_chain: Ordered list of the transformations applied to the raw observable. This is
            the field that makes the quantity comparable: two numbers with different chains are
            different quantities regardless of sharing a name.
        identifiability: Directions in parameter space this quantity actually constrains. A quantity
            that constrains nothing is not evidence, however precise it is.
        detail: Anything else on the record.
    """

    name: str
    operator_id: str
    value: float
    units: str
    standard_error: float | None
    n_effective: float
    manifest_hash: str
    analysis_chain: tuple[str, ...]
    identifiability: tuple[str, ...] = ()
    detail: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.analysis_chain:
            raise ValueError(
                "an apparent quantity must state the analysis chain that produced it; an empty "
                "chain claims the number came from nowhere, which is the one thing it cannot have "
                "done"
            )
        if self.standard_error == 0.0:
            raise ValueError(
                f"{self.name!r} carries standard_error = 0.0, which claims infinite precision. "
                "This field is documented as never zero as a stand-in: a reduction that found no "
                "scatter must report None and say why, so the reader sees 'nothing varied' rather "
                "than 'exact'. Every operator here reaches this by the same route -- a constant "
                "input series -- so the guard is shared rather than repeated per operator"
            )

    def as_dict(self) -> dict[str, Any]:
        """JSON-able form, for a run record."""
        return {
            "refused": False,
            "name": self.name,
            "operator_id": self.operator_id,
            "value": self.value,
            "units": self.units,
            "standard_error": self.standard_error,
            "n_effective": self.n_effective,
            "manifest_hash": self.manifest_hash,
            "analysis_chain": list(self.analysis_chain),
            "identifiability": list(self.identifiability),
            "detail": dict(self.detail),
        }


@runtime_checkable
class ObservationOperator(Protocol):
    """Everything an observation must declare about itself before it may report a number.

    The declarations are the point. An operator that cannot say which owners it reads, which
    parameter directions it constrains, or where it stops being valid is not an instrument, it is a
    function that returns a float.
    """

    @property
    def operator_id(self) -> str:
        """Stable identifier, including a version. Recorded on every result and every refusal."""
        ...

    @property
    def manifest(self) -> ProtocolManifest:
        """The protocol in force. Its hash travels with every result this operator produces."""
        ...

    @property
    def required_state_owners(self) -> tuple[str, ...]:
        """Owners this operator reads. Checked before any arithmetic happens."""
        ...

    @property
    def identifiability_directions(self) -> tuple[str, ...]:
        """Which directions in parameter space this observation can constrain.

        Stated so that a degeneracy is a declared property rather than a discovery made after a
        posterior fails to contract. For the passive spectrum the honest statement is that it
        constrains ``kappa`` and ``sigma`` only in the combination the accessible ``q`` band
        resolves, and not at all along an overall amplitude rescale.
        """
        ...

    @property
    def applicability(self) -> ApplicabilityDomain:
        """Where this operator is entitled to report a number."""
        ...

    def raw_observable(
        self, state: AcceptedState, context: ObservationContext
    ) -> RawObservable | Refusal:
        """Produce what the instrument would produce, with no interpretation applied."""
        ...

    def apparent_quantity(self, raw: RawObservable) -> ApparentQuantity | Refusal:
        """Run the raw observable through the analysis a real experiment would use."""
        ...


def observe(
    operator: ObservationOperator, state: AcceptedState, context: ObservationContext
) -> ApparentQuantity | Refusal:
    """Run both halves in order, short-circuiting on the first refusal.

    A convenience, and a narrow one on purpose: it does not hide the two-step structure, it just
    saves the caller from writing the refusal branch twice. Any operator that wants to report only
    the raw observable should call :meth:`ObservationOperator.raw_observable` directly.
    """
    raw = operator.raw_observable(state, context)
    if isinstance(raw, Refusal):
        return raw
    return operator.apparent_quantity(raw)


@dataclass(frozen=True, slots=True)
class ArrayState:
    """A minimal :class:`AcceptedState` over plain arrays.

    Exists so that this lane can be tested, and so that another lane can hand over arrays without
    either side importing the other. When L2's accepted-state object lands it will satisfy the
    Protocol structurally and this class stays useful for tests.
    """

    arrays: Mapping[str, Mapping[str, np.ndarray]]

    @property
    def owners(self) -> frozenset[str]:
        return frozenset(self.arrays)

    def array(self, owner: str, name: str) -> np.ndarray:
        try:
            block = self.arrays[owner]
        except KeyError as exc:  # pragma: no cover - guarded by check_owners at every call site
            raise KeyError(f"no such state owner: {owner!r}") from exc
        try:
            return block[name]
        except KeyError as exc:
            raise KeyError(f"owner {owner!r} carries no array named {name!r}") from exc


def first_non_finite(arrays: Sequence[np.ndarray]) -> int | None:
    """Index of the first array containing a non-finite entry, or ``None``.

    Used by operators to refuse before computing rather than after. Checking up front costs one
    pass and buys the guarantee that no NaN was ever fed to a reduction.
    """
    for index, item in enumerate(arrays):
        data = np.asarray(item)
        if data.size and not np.all(np.isfinite(data)):
            return index
    return None
