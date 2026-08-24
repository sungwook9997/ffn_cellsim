"""Two independent axes for grading a result: how far it ran, and whether its number may be quoted.

Aleph separates these because collapsing them is how a project talks itself into publishing.  A
single ladder forces every reader to guess which question a label answers, and the guess is always
generous: "it ran on the GPU at native population" gets read as "so the number is good".  It is not.
Running is a *structural* fact about code paths.  Quotability is a *quantitative* fact about whether
a magnitude survived its own error budget.  A run can be maximal on the first axis and refused on
the second, and that combination is a perfectly good outcome, not a shortfall.

The two types here are therefore deliberately non-comparable.  :class:`EvidenceRung` is ordered.
:class:`QuantitativeStatus` is ordered.  There is no function anywhere in Aleph that maps one to the
other, and :func:`classify_claim` is written so that adding one later would be an obvious edit rather
than an accident.

The single enforcement point is :func:`classify_claim`.  It is the only supported way to build an
:class:`EvidenceLabel`, and it refuses to mint ``QUOTABLE`` unless the caller names the artifact
digest that carries the supporting evidence.  A number without a retrievable artifact behind it is an
assertion, and Aleph records assertions as ``ASSUMED``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum, IntEnum, unique
from typing import Sequence

__all__ = [
    "DIGEST_PATTERN",
    "EvidenceLabel",
    "EvidenceRung",
    "QuantitativeStatus",
    "RUNG_ORDER",
    "classify_claim",
    "is_artifact_digest",
    "rung_at_least",
]


@unique
class EvidenceRung(IntEnum):
    """How far a result actually executed, weakest to strongest.  Structural only.

    The integer values are sparse on purpose.  A rung that has to be inserted later (say, a CPU
    native-population slice) gets a value between two existing ones without renumbering anything that
    is already written into artifacts on disk, and without a side table that can drift out of step
    with the enum.  The ordering *is* the type: ``EvidenceRung.CPU_UNIT < EvidenceRung.GPU_UNIT``
    is true because of how the members are declared, not because of a dict somebody maintains.

    Serialize with :attr:`wire` (the member name) and read back with :meth:`from_wire`.  The integer
    is an ordering device and must never reach an artifact, because reordering the ladder would then
    silently reinterpret every record ever written.

    Members:
        ASSUMED: A human asserted it.  Nothing executed.  This is the floor and the default, and it
            is where every claim starts.
        ANALYTIC_ORACLE: A closed form was evaluated.  Real arithmetic happened, but no Aleph
            simulation ran, so this grades the oracle rather than the runtime.
        CPU_UNIT: Aleph runtime code executed on CPU at unit scale.
        GPU_UNIT: The same executed on the ratified production backend at unit scale.
        GPU_NATIVE_SLICE: Executed at native population and geometry, but the candidate step was not
            put through the acceptance transaction — so the state it produced is not accepted state.
        GPU_NATIVE_ACCEPTED: Executed at native population and the acceptance transaction accepted
            the step.  This is the strongest structural fact Aleph can establish, and it still says
            nothing about whether the resulting magnitude may be quoted.
    """

    ASSUMED = 0
    ANALYTIC_ORACLE = 10
    CPU_UNIT = 20
    GPU_UNIT = 30
    GPU_NATIVE_SLICE = 40
    GPU_NATIVE_ACCEPTED = 50

    @property
    def wire(self) -> str:
        """The stable string written into artifacts."""
        return self.name

    @classmethod
    def from_wire(cls, wire: str) -> EvidenceRung:
        """Read a rung back from its artifact spelling.

        Args:
            wire: The member name as written by :attr:`wire`.

        Returns:
            The matching rung.

        Raises:
            ValueError: If the spelling names no declared rung.  Unknown rungs are refused rather
                than mapped to the floor, because a record written by a newer schema must not be
                quietly downgraded into something a reader will accept.
        """
        try:
            return cls[wire]
        except KeyError:
            raise ValueError(
                f"{wire!r} is not a declared evidence rung; "
                f"expected one of {[member.name for member in cls]}"
            ) from None


@unique
class QuantitativeStatus(Enum):
    """Whether the magnitude a result produced may leave the artifact.  Quantitative only.

    Ordered weakest to strongest, but the ordering is *within this axis only*.  Nothing on
    :class:`EvidenceRung` moves a claim along this axis, and nothing here moves a claim along that
    one.

    Members:
        BLOCKED: The number may not be quoted anywhere, for any purpose.  This is the default, and
            it is the correct state for a run whose mechanism is demonstrated but whose error budget
            is unestablished, whose parameters are unsourced, or whose solver did not converge.
        PROVISIONAL: The number may be quoted *inside Aleph*, marked as provisional, to steer the
            next experiment.  It may not leave the project and may not enter a manuscript.
        QUOTABLE: The number may be quoted externally.  Reaching this requires naming the artifact
            digest that supports it — see :func:`classify_claim`.
    """

    BLOCKED = 0
    PROVISIONAL = 1
    QUOTABLE = 2

    @property
    def wire(self) -> str:
        """The stable string written into artifacts."""
        return self.name

    @classmethod
    def from_wire(cls, wire: str) -> QuantitativeStatus:
        """Read a status back from its artifact spelling.

        Args:
            wire: The member name as written by :attr:`wire`.

        Returns:
            The matching status.

        Raises:
            ValueError: If the spelling names no declared status.
        """
        try:
            return cls[wire]
        except KeyError:
            raise ValueError(
                f"{wire!r} is not a declared quantitative status; "
                f"expected one of {[member.name for member in cls]}"
            ) from None


#: The rungs in ladder order, weakest first.  Derived from the enum so it cannot disagree with it.
RUNG_ORDER: tuple[EvidenceRung, ...] = tuple(sorted(EvidenceRung))

#: The shape an artifact digest must have to count as *naming* an artifact: an algorithm label, a
#: colon, and at least 32 hex characters.  The pattern exists so that a caller cannot satisfy the
#: quotability requirement with ``"see the spreadsheet"`` or ``"TODO"``.  It checks form, not
#: existence — proving the artifact is on disk is the reader's job in ``aleph/artifacts``.
DIGEST_PATTERN: re.Pattern[str] = re.compile(r"^(sha256|blake2b):[0-9a-f]{32,128}$")


def is_artifact_digest(candidate: object) -> bool:
    """Return whether ``candidate`` has the form of an artifact digest.

    Args:
        candidate: Anything; non-strings are simply not digests.

    Returns:
        ``True`` when the value matches :data:`DIGEST_PATTERN`.
    """
    return isinstance(candidate, str) and DIGEST_PATTERN.match(candidate) is not None


def rung_at_least(rung: EvidenceRung, floor: EvidenceRung) -> bool:
    """Return whether ``rung`` sits at or above ``floor`` on the structural ladder.

    Provided so that callers comparing rungs do so through a named function rather than by reaching
    for the integer values, which are an implementation detail of the ordering.

    Args:
        rung: The rung to test.
        floor: The rung to test it against.

    Returns:
        ``True`` when ``rung`` is at least as strong as ``floor``.

    Raises:
        TypeError: If either argument is not an :class:`EvidenceRung`.
    """
    if not isinstance(rung, EvidenceRung) or not isinstance(floor, EvidenceRung):
        raise TypeError("rung comparison takes EvidenceRung members, not strings or integers")
    return int(rung) >= int(floor)


@dataclass(frozen=True, slots=True)
class EvidenceLabel:
    """The pair of axis values a claim carries, plus the evidence for each.

    Build these with :func:`classify_claim`.  The constructor is intentionally permissive about
    field combinations because the *policy* lives in one place; a dataclass that re-implemented the
    policy in ``__post_init__`` would be a second copy of it to keep in step.  What the constructor
    does enforce is that the fields are well typed and that the recorded outcome is internally
    consistent with itself.

    Attributes:
        rung: How far the result executed.
        status: Whether its magnitude may be quoted.
        basis: What was actually measured or executed to earn the rung.  Required and non-empty:
            a rung with no stated basis records an intention rather than an observation.
        supporting_artifact_digest: The digest of the artifact carrying the supporting evidence.
            Required whenever ``status`` is ``QUOTABLE``; optional otherwise.
        blocking_reasons: Every reason the quantitative axis was held below what the caller asked
            for.  Empty when nothing was held back.
        requested_status: What the caller asked for, kept so a demotion is visible in the record
            rather than looking like the caller's own modest choice.
    """

    rung: EvidenceRung
    status: QuantitativeStatus
    basis: str
    supporting_artifact_digest: str | None = None
    blocking_reasons: tuple[str, ...] = ()
    requested_status: QuantitativeStatus | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.rung, EvidenceRung):
            raise TypeError("rung must be an EvidenceRung, never a free string")
        if not isinstance(self.status, QuantitativeStatus):
            raise TypeError("status must be a QuantitativeStatus, never a free string")
        if not self.basis.strip():
            raise ValueError(
                f"rung {self.rung.wire} was asserted with no basis; name what was executed or "
                "measured to earn it"
            )
        if self.status is QuantitativeStatus.QUOTABLE and not is_artifact_digest(
            self.supporting_artifact_digest
        ):
            raise ValueError(
                "a QUOTABLE label must name the artifact digest that supports it; got "
                f"{self.supporting_artifact_digest!r}"
            )
        if self.status is QuantitativeStatus.QUOTABLE and self.blocking_reasons:
            raise ValueError(
                "a label cannot be QUOTABLE while carrying blocking reasons: "
                f"{list(self.blocking_reasons)}"
            )

    @property
    def was_demoted(self) -> bool:
        """Whether the recorded status is weaker than what the caller asked for."""
        return (
            self.requested_status is not None
            and self.status.value < self.requested_status.value
        )

    def as_artifact_fields(self) -> dict[str, object]:
        """Return the flat fields a record writes.  Both axes are always present.

        Returns:
            A JSON-able mapping.  Both axis keys appear unconditionally, so a reader never has to
            infer a missing quantitative status from a present rung.
        """
        fields: dict[str, object] = {
            "evidence_rung": self.rung.wire,
            "quantitative_status": self.status.wire,
            "evidence_basis": self.basis,
            "supporting_artifact_digest": self.supporting_artifact_digest,
        }
        if self.blocking_reasons:
            fields["quantitative_blocking_reasons"] = list(self.blocking_reasons)
        if self.was_demoted:
            assert self.requested_status is not None  # narrowed by was_demoted
            fields["requested_quantitative_status"] = self.requested_status.wire
        return fields


def classify_claim(
    *,
    rung: EvidenceRung,
    basis: str,
    requested_status: QuantitativeStatus = QuantitativeStatus.BLOCKED,
    supporting_artifact_digest: str | None = None,
    blocking_reasons: Sequence[str] = (),
) -> EvidenceLabel:
    """Grade a claim on both axes, demoting the quantitative one rather than refusing outright.

    This is the single place where Aleph's quotability policy is written down.  Three rules, and
    the first is the one the whole module exists for.

    1. **Structural coverage never implies quantitative validity.**  ``rung`` is read here only to
       apply rule 3.  It cannot raise ``requested_status``, and there is no code path in which a
       stronger rung yields a stronger status than a weaker rung would have.  Reaching
       ``GPU_NATIVE_ACCEPTED`` grants nothing on the quantitative axis.

    2. **Any blocking reason wins.**  If the caller supplies even one reason the number is not
       trustworthy, the result is ``BLOCKED``, whatever was asked for and whatever the rung.  The
       reasons are carried into the label so the record says *why*, and so a later reader can check
       whether the reason has since been retired.

    3. **QUOTABLE requires a named artifact, and requires something to have executed.**  Without a
       well-formed ``supporting_artifact_digest`` the request is demoted to ``PROVISIONAL`` with a
       blocking reason — Aleph does not raise here, because refusing to construct the label would
       tempt callers to hand-build one.  Demoting produces a record that is honest *and* keeps the
       provenance of the attempt.  Separately, ``ASSUMED`` cannot carry a quotable number; that is
       definitional rather than a gate, since nothing was executed to produce one.

    Args:
        rung: The structural rung the result earned.
        basis: What was executed or measured to earn ``rung``.  Must be non-empty.
        requested_status: The quantitative status the caller believes is warranted.  Defaults to
            ``BLOCKED`` so that silence never produces a quotable number.
        supporting_artifact_digest: Digest of the artifact carrying the supporting evidence, in the
            form checked by :data:`DIGEST_PATTERN`.
        blocking_reasons: Reasons the magnitude should not be quoted, each a short sentence.

    Returns:
        An :class:`EvidenceLabel` whose ``status`` is at most ``requested_status``.

    Raises:
        TypeError: If ``rung`` or ``requested_status`` is not the right enum type.
        ValueError: If ``basis`` is empty, or if a supplied blocking reason is empty.
    """
    if not isinstance(rung, EvidenceRung):
        raise TypeError("rung must be an EvidenceRung, never a free string")
    if not isinstance(requested_status, QuantitativeStatus):
        raise TypeError("requested_status must be a QuantitativeStatus, never a free string")
    if not basis.strip():
        raise ValueError("classify_claim needs a basis naming what was executed or measured")

    reasons: list[str] = []
    for reason in blocking_reasons:
        if not str(reason).strip():
            raise ValueError("a blocking reason must say something; got an empty string")
        reasons.append(str(reason).strip())

    status = requested_status

    # Rule 2 — a stated doubt outranks everything, including a maximal rung.
    if reasons:
        status = QuantitativeStatus.BLOCKED

    # Rule 3 — quotability is bought with an artifact, not with a rung.
    if status is QuantitativeStatus.QUOTABLE and not is_artifact_digest(
        supporting_artifact_digest
    ):
        reasons.append(
            "QUOTABLE was requested without naming a supporting artifact digest, so the magnitude "
            "cannot be retrieved or re-checked by a reader; demoted to PROVISIONAL"
        )
        status = QuantitativeStatus.PROVISIONAL

    if status is QuantitativeStatus.QUOTABLE and rung is EvidenceRung.ASSUMED:
        reasons.append(
            "the rung is ASSUMED, meaning nothing executed to produce this magnitude; an artifact "
            "digest cannot substitute for a measurement that was never made"
        )
        status = QuantitativeStatus.BLOCKED

    # A demotion path may have appended reasons after clearing a QUOTABLE request; the invariant the
    # label enforces is that QUOTABLE and blocking reasons never coexist.
    if status is QuantitativeStatus.QUOTABLE:
        reasons = []

    return EvidenceLabel(
        rung=rung,
        status=status,
        basis=basis.strip(),
        supporting_artifact_digest=supporting_artifact_digest,
        blocking_reasons=tuple(reasons),
        requested_status=requested_status,
    )
