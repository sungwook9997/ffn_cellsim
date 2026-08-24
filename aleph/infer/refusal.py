"""The vocabulary the inference layer uses to say "no", in a form a caller can act on.

An inference layer that can only return numbers will return a number for every question, including
the questions its data cannot answer.  The usual disguise is a wide error bar: the machinery reports
``kappa = 3 +/- 400`` and the reader takes away "3".  The honest answer in that case is not a number
at all, it is *this observation cannot see this quantity*, and the two statements are different in
kind, not in degree.  So refusal here is a first-class return type rather than an exception someone
remembers to catch, and it is a type with a required remedy.

Three properties are enforced rather than encouraged.

**A refusal names the offending input.**  "Out of distribution" with nothing attached cannot be
reproduced or argued with.  ``offending_input`` carries the thing that triggered it.

**A refusal names what would lift it.**  ``lift_requires`` may not be empty.  A refusal without a
remedy is a dead end, and dead ends get worked around by deleting the check.  If the remedy is
"nothing you can do with this dataset", that must be *said*, because then the right move is to run a
different simulation --- which is what :class:`NativeQueryRequest` exists to express.

**A refusal is both a value and an exception.**  Some call sites want to branch on a result; others
sit under an optimizer and need to abort.  Every :class:`Refusal` carries the exception type that
transports it, so the two paths never drift apart: ``refusal.as_error()`` and
``except UnidentifiableDirectionError as e: e.refusal`` are inverses.

:class:`NativeQueryRequest` is the structured version of "I cannot answer this; run this specific
simulation instead".  It is deliberately concrete --- which observable, at which parameters, how
many samples, and what result would lift the refusal --- because a request that cannot be executed
without a conversation is a complaint, not a request.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, ClassVar, Mapping

__all__ = [
    "InsufficientData",
    "InsufficientDataError",
    "ModelInadequate",
    "ModelInadequateError",
    "NativeQueryRequest",
    "OutOfDistribution",
    "OutOfDistributionError",
    "ProvenanceIncomplete",
    "ProvenanceIncompleteError",
    "Refusal",
    "RefusalError",
    "Unidentifiable",
    "UnidentifiableDirectionError",
    "is_refusal",
]


class RefusalError(Exception):
    """Transport for a :class:`Refusal` through a call stack that cannot return one.

    The refusal itself is the payload and stays reachable as :attr:`refusal`; the message is a
    rendering of it.  Callers should read the payload, not parse the string.
    """

    def __init__(self, refusal: "Refusal") -> None:
        super().__init__(refusal.message())
        self.refusal = refusal


class OutOfDistributionError(RefusalError):
    """The input lies outside the region where the model or surrogate was ever checked."""


class UnidentifiableDirectionError(RefusalError):
    """The requested quantity lies in, or too near, the null space of the observation."""


class InsufficientDataError(RefusalError):
    """There is not enough data for the requested inference to be meaningful."""


class ModelInadequateError(RefusalError):
    """The model cannot represent the data at all; a parameter estimate would be meaningless."""


class ProvenanceIncompleteError(RefusalError):
    """The evidence lineage is missing, so no number derived from it may be quoted."""


def _short(value: Any, limit: int = 160) -> str:
    """A repr short enough to live in an exception message without burying it."""
    try:
        text = repr(value)
    except Exception:  # pragma: no cover - a repr that raises is still worth reporting around
        text = f"<unreprable {type(value).__name__}>"
    if len(text) > limit:
        text = text[: limit - 3] + "..."
    return text


@dataclass(frozen=True)
class NativeQueryRequest:
    """A specific simulation that would answer the question the current data cannot.

    Attributes:
        question: The thing that could not be answered, in one sentence.
        observable: Which observation operator to run.  A name, not a description.
        parameters: Where in parameter space to run it.
        n_samples: How many independent samples are needed.  Must be positive --- "more data"
            without a number cannot be scheduled.
        rationale: Why *this* run, rather than reprocessing what already exists.
        acceptance: What result would lift the refusal.  Stated before the run, so the answer
            cannot be graded after the fact.
    """

    question: str
    observable: str
    parameters: Mapping[str, float]
    n_samples: int
    rationale: str
    acceptance: str

    def __post_init__(self) -> None:
        for name in ("question", "observable", "rationale", "acceptance"):
            if not str(getattr(self, name)).strip():
                raise ValueError(
                    f"NativeQueryRequest.{name} may not be empty: a request that cannot be "
                    "executed without a follow-up conversation is a complaint, not a request."
                )
        if self.n_samples <= 0:
            raise ValueError(
                f"NativeQueryRequest.n_samples must be positive, got {self.n_samples}: "
                "'more data' with no number attached cannot be scheduled."
            )
        object.__setattr__(self, "parameters", dict(self.parameters))

    def render(self) -> str:
        params = ", ".join(f"{k}={v!r}" for k, v in sorted(self.parameters.items()))
        return (
            f"NATIVE QUERY: run observable {self.observable!r} at ({params}) "
            f"with n_samples={self.n_samples}. "
            f"Question: {self.question} Rationale: {self.rationale} Accept if: {self.acceptance}"
        )


@dataclass(frozen=True)
class Refusal:
    """Base class for every typed refusal in the inference layer.

    Attributes:
        reason: Why the answer is being withheld.
        offending_input: The input that triggered it, kept so the refusal is reproducible.
        lift_requires: What would have to change for the answer to become available.  Required.
        detail: Machine-readable context; whatever the specific check computed.
        native_query: An optional concrete simulation that would resolve the refusal.
    """

    reason: str
    offending_input: Any
    lift_requires: str
    detail: Mapping[str, Any] = field(default_factory=dict)
    native_query: NativeQueryRequest | None = None

    #: The exception subclass that carries this refusal.  Set by each concrete refusal.
    error_type: ClassVar[type[RefusalError]] = RefusalError
    #: Stable wire label.  Written into artifacts; never renumber, only add.
    kind: ClassVar[str] = "refusal"

    def __post_init__(self) -> None:
        if not str(self.reason).strip():
            raise ValueError("A refusal must state a reason.")
        if not str(self.lift_requires).strip():
            raise ValueError(
                "A refusal must state what would lift it. A refusal with no remedy is a dead end, "
                "and dead ends get worked around by deleting the check. If nothing about this "
                "dataset can lift it, say that and attach a NativeQueryRequest."
            )
        object.__setattr__(self, "detail", dict(self.detail))

    def message(self) -> str:
        text = (
            f"[{self.kind}] {self.reason} "
            f"| offending input: {_short(self.offending_input)} "
            f"| to lift: {self.lift_requires}"
        )
        if self.native_query is not None:
            text += " | " + self.native_query.render()
        return text

    def as_error(self) -> RefusalError:
        """The exception that transports this refusal, for call sites that cannot return a value."""
        return self.error_type(self)

    def raise_(self) -> "None":
        raise self.as_error()


@dataclass(frozen=True)
class OutOfDistribution(Refusal):
    """The input is outside the region where the forward model or surrogate was ever checked.

    Extrapolation is not an error, it is an unbounded one; the point of this refusal is that the
    error bar beside an extrapolated estimate describes the noise, not the extrapolation.
    """

    statistic: float | None = None
    threshold: float | None = None

    error_type: ClassVar[type[RefusalError]] = OutOfDistributionError
    kind: ClassVar[str] = "out_of_distribution"


@dataclass(frozen=True)
class Unidentifiable(Refusal):
    """The requested quantity lies in, or too near, the null space of this observation.

    ``direction`` is the offending parameter combination, in the coordinates the analysis ran in.
    Reporting it matters: it converts "we cannot measure kappa" into "we can only measure
    kappa*sigma", which is a usable scientific statement and often the actual result.
    """

    direction: tuple[float, ...] | None = None
    parameter_names: tuple[str, ...] = ()

    error_type: ClassVar[type[RefusalError]] = UnidentifiableDirectionError
    kind: ClassVar[str] = "unidentifiable"


@dataclass(frozen=True)
class InsufficientData(Refusal):
    """Fewer independent observations than the requested inference needs."""

    have: int | None = None
    need: int | None = None

    error_type: ClassVar[type[RefusalError]] = InsufficientDataError
    kind: ClassVar[str] = "insufficient_data"


@dataclass(frozen=True)
class ModelInadequate(Refusal):
    """The model cannot represent the data, so its best-fit parameters describe nothing.

    This is the refusal that an a-priori noise variance makes possible.  Fit the variance instead
    and the reduced chi-squared is one by construction, so this check can never fire.
    """

    statistic: float | None = None
    p_value: float | None = None

    error_type: ClassVar[type[RefusalError]] = ModelInadequateError
    kind: ClassVar[str] = "model_inadequate"


@dataclass(frozen=True)
class ProvenanceIncomplete(Refusal):
    """The chain from artifact to number is broken, so the number is an assertion."""

    missing_fields: tuple[str, ...] = ()

    error_type: ClassVar[type[RefusalError]] = ProvenanceIncompleteError
    kind: ClassVar[str] = "provenance_incomplete"


def is_refusal(value: Any) -> bool:
    """True when ``value`` is a refusal, for call sites that return ``Result | Refusal``."""
    return isinstance(value, Refusal)
