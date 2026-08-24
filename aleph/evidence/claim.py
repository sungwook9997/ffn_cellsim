"""Quantitative claims and an append-only register in which retraction is a record, not a deletion.

A claim is one number Aleph is prepared to say out loud: a statement, a value, a unit, a physical
dimension, where on the two axes of :mod:`aleph.evidence.ladder` it sits, and which artifact backs
it.  A register is the ordered log of everything that has ever been claimed.

The design rule is that **the log only grows**.  Retracting a claim appends a retraction record;
nothing is removed and nothing is edited in place.  This is not tidiness for its own sake.  A
project that deletes withdrawn numbers loses the only evidence that it withdrew them, and the same
number quietly comes back six weeks later because nobody can point at the reason it left.  Aleph
keeps the reason attached to the identifier forever.

The second rule follows from the first: **a retracted claim cannot be re-promoted to QUOTABLE on the
evidence that was already retracted.**  If artifact ``sha256:aaa…`` supported a claim that was
withdrawn, then ``sha256:aaa…`` is spent.  Re-quoting the number requires a *new* artifact — a new
run, a new oracle, a new measurement.  :meth:`ClaimRegister.register` enforces this by refusing any
QUOTABLE claim whose supporting digest appears among the register's retracted digests.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field, replace
from enum import Enum, unique
from typing import Iterator, Mapping, Sequence

from aleph.evidence.ladder import (
    EvidenceLabel,
    EvidenceRung,
    QuantitativeStatus,
    is_artifact_digest,
)

__all__ = [
    "CLAIM_ID_PATTERN",
    "ClaimRegister",
    "ClaimRegisterError",
    "QuantitativeClaim",
    "RegisterEntry",
    "RegisterEvent",
    "Retraction",
]

#: A claim identifier is a short, stable, lower-kebab token.  Constrained so that identifiers can be
#: used as filenames, dictionary keys, and cross-references in a manuscript without escaping.
CLAIM_ID_PATTERN: re.Pattern[str] = re.compile(r"^[a-z0-9][a-z0-9._-]{2,63}$")


class ClaimRegisterError(RuntimeError):
    """Raised when an operation would break the register's append-only or re-promotion rules."""


@unique
class RegisterEvent(Enum):
    """What kind of thing one entry in the register log is.

    Members:
        ASSERTED: A claim entered the register for the first time.
        RETRACTED: A previously asserted claim was withdrawn, with a reason.
        SUPERSEDED: A retracted claim was replaced by a named successor claim.  Recorded as its own
            event so that "withdrawn and abandoned" and "withdrawn and replaced" are distinguishable
            in the log without reading the claims.
    """

    ASSERTED = "ASSERTED"
    RETRACTED = "RETRACTED"
    SUPERSEDED = "SUPERSEDED"


@dataclass(frozen=True, slots=True)
class Retraction:
    """The record of a claim being withdrawn.  Attached to the claim, never replacing it.

    Attributes:
        reason: Why the number is no longer sayable.  Required and non-empty — a retraction with no
            reason is indistinguishable from a mistake, and gives a later reader no way to tell
            whether the problem has since been fixed.
        retracted_by: Who or what performed the retraction (a person, a gate, a review).  Required,
            because "the number was withdrawn" and "the PI withdrew the number" carry different
            weight and the register should not blur them.
        superseded_by: The identifier of the claim that replaces this one, if any.  ``None`` means
            the number was withdrawn and nothing stands in its place, which is a legitimate and
            common outcome.
        retracted_at: An ISO-8601 timestamp, recorded as a string so the register is exactly what
            gets serialized with no timezone reinterpretation on the way through.
    """

    reason: str
    retracted_by: str
    superseded_by: str | None = None
    retracted_at: str | None = None

    def __post_init__(self) -> None:
        if not self.reason.strip():
            raise ValueError("a retraction must state a reason")
        if not self.retracted_by.strip():
            raise ValueError("a retraction must name who or what retracted the claim")
        if self.superseded_by is not None and not CLAIM_ID_PATTERN.match(self.superseded_by):
            raise ValueError(
                f"superseded_by must be a claim identifier matching {CLAIM_ID_PATTERN.pattern!r}; "
                f"got {self.superseded_by!r}"
            )

    def as_json_obj(self) -> dict[str, object]:
        """Return the JSON-able record of this retraction."""
        return {
            "reason": self.reason,
            "retracted_by": self.retracted_by,
            "superseded_by": self.superseded_by,
            "retracted_at": self.retracted_at,
        }


@dataclass(frozen=True, slots=True)
class QuantitativeClaim:
    """One number Aleph is prepared to state, with everything a reader needs to check it.

    A claim carries its own evidence label rather than a bare rung, so the two axes travel together
    and neither can be read without the other.

    Attributes:
        claim_id: Stable identifier, matching :data:`CLAIM_ID_PATTERN`.  Never reused, even after
            retraction; a superseding claim gets a new identifier.
        statement: The sentence a human would read, in Aleph's own words.
        value: The magnitude.  Must be finite.
        unit: The unit the magnitude is expressed in, spelled exactly as the record will carry it
            (for example ``"J"``, ``"N/m"``, ``"m^-1"``).  Required.
        dimension: The physical dimension the unit belongs to (for example ``"energy"``,
            ``"surface_tension"``).  Required and separate from ``unit``, because a unit typo and a
            dimensional error are different mistakes and only the second is fatal.  Lane L1 owns the
            dimensional vocabulary; this field records which term of it applies.
        label: The two-axis grading from :func:`aleph.evidence.ladder.classify_claim`.
        retraction: ``None`` while the claim stands; a :class:`Retraction` once withdrawn.
    """

    claim_id: str
    statement: str
    value: float
    unit: str
    dimension: str
    label: EvidenceLabel
    retraction: Retraction | None = None

    def __post_init__(self) -> None:
        if not CLAIM_ID_PATTERN.match(self.claim_id):
            raise ValueError(
                f"claim_id must match {CLAIM_ID_PATTERN.pattern!r}; got {self.claim_id!r}"
            )
        if not self.statement.strip():
            raise ValueError(f"claim {self.claim_id!r} needs a statement a human can read")
        if not isinstance(self.value, (int, float)) or isinstance(self.value, bool):
            raise TypeError(f"claim {self.claim_id!r} value must be a real number")
        if self.value != self.value or self.value in (float("inf"), float("-inf")):
            raise ValueError(f"claim {self.claim_id!r} value must be finite; got {self.value!r}")
        if not self.unit.strip():
            raise ValueError(
                f"claim {self.claim_id!r} needs a unit; a bare magnitude is not a physical claim"
            )
        if not self.dimension.strip():
            raise ValueError(
                f"claim {self.claim_id!r} needs a dimension, separately from its unit, so that a "
                "dimensional error is visible without parsing the unit string"
            )
        if not isinstance(self.label, EvidenceLabel):
            raise TypeError(f"claim {self.claim_id!r} label must be an EvidenceLabel")

    @property
    def is_retracted(self) -> bool:
        """Whether this claim has been withdrawn."""
        return self.retraction is not None

    @property
    def is_quotable(self) -> bool:
        """Whether this claim may be quoted externally right now.

        A retracted claim is never quotable, whatever its label says.  The label is the grade the
        evidence earned; the retraction is a later fact that overrides it.
        """
        return (
            not self.is_retracted
            and self.label.status is QuantitativeStatus.QUOTABLE
        )

    @property
    def supporting_artifact_digest(self) -> str | None:
        """The digest of the artifact backing this claim, if one was named."""
        return self.label.supporting_artifact_digest

    @property
    def rung(self) -> EvidenceRung:
        """The structural rung this claim's evidence reached."""
        return self.label.rung

    @property
    def status(self) -> QuantitativeStatus:
        """The quantitative status the evidence earned, before any retraction is applied."""
        return self.label.status

    def with_retraction(self, retraction: Retraction) -> QuantitativeClaim:
        """Return a copy of this claim carrying ``retraction``.

        The original object is unchanged, which is what lets the register keep both the asserted and
        the retracted form of the same identifier in its log.

        Args:
            retraction: The withdrawal record.

        Returns:
            A new claim with :attr:`retraction` set.

        Raises:
            ClaimRegisterError: If the claim is already retracted.  Re-retracting would overwrite
                the first reason, which is the one deletion this module exists to prevent.
        """
        if self.retraction is not None:
            raise ClaimRegisterError(
                f"claim {self.claim_id!r} is already retracted for: {self.retraction.reason!r}. "
                "Retracting again would overwrite the original reason"
            )
        return replace(self, retraction=retraction)

    def as_json_obj(self) -> dict[str, object]:
        """Return the JSON-able record of this claim, both axes and the retraction included."""
        record: dict[str, object] = {
            "claim_id": self.claim_id,
            "statement": self.statement,
            "value": float(self.value),
            "unit": self.unit,
            "dimension": self.dimension,
            "retracted": self.is_retracted,
            "quotable_now": self.is_quotable,
        }
        record.update(self.label.as_artifact_fields())
        record["retraction"] = None if self.retraction is None else self.retraction.as_json_obj()
        return record


@dataclass(frozen=True, slots=True)
class RegisterEntry:
    """One immutable line of the register log.

    Attributes:
        sequence: Position in the log, starting at zero.  Gaps are impossible by construction.
        event: What happened.
        claim: The full claim as it stood *after* the event.  Storing the whole claim rather than a
            delta means any line of the log can be read in isolation.
        note: Optional free text about this particular log line.
    """

    sequence: int
    event: RegisterEvent
    claim: QuantitativeClaim
    note: str | None = None

    def as_json_obj(self) -> dict[str, object]:
        """Return the JSON-able record of this log line."""
        return {
            "sequence": self.sequence,
            "event": self.event.value,
            "note": self.note,
            "claim": self.claim.as_json_obj(),
        }


@dataclass(slots=True)
class ClaimRegister:
    """An append-only log of claims, retractions, and supersessions.

    Nothing in the public surface removes or mutates an entry.  :meth:`history` hands back a tuple,
    :meth:`current` reconstructs the present state of an identifier by reading the log forward, and
    the log itself is private.

    The register also holds the spent-evidence rule.  Every digest that ever supported a retracted
    claim is remembered in :attr:`retracted_digests`, and :meth:`register` refuses to mint a new
    QUOTABLE claim on any of them.  This is the guard that stops a withdrawn number from returning
    under a new identifier with the same evidence behind it.
    """

    _log: list[RegisterEntry] = field(default_factory=list, repr=False)

    def __len__(self) -> int:
        """Return the number of log lines, not the number of live claims."""
        return len(self._log)

    def __iter__(self) -> Iterator[RegisterEntry]:
        """Iterate the log oldest first."""
        return iter(tuple(self._log))

    def history(self) -> tuple[RegisterEntry, ...]:
        """Return every log line, oldest first.

        Returns:
            An immutable snapshot.  Callers cannot append through it.
        """
        return tuple(self._log)

    def history_for(self, claim_id: str) -> tuple[RegisterEntry, ...]:
        """Return every log line touching one identifier, oldest first.

        Args:
            claim_id: The identifier to filter on.

        Returns:
            The matching log lines, possibly empty.
        """
        return tuple(entry for entry in self._log if entry.claim.claim_id == claim_id)

    def current(self, claim_id: str) -> QuantitativeClaim | None:
        """Return the present state of one identifier, or ``None`` if it was never asserted.

        Args:
            claim_id: The identifier to resolve.

        Returns:
            The claim as of the last log line that touched it.
        """
        found: QuantitativeClaim | None = None
        for entry in self._log:
            if entry.claim.claim_id == claim_id:
                found = entry.claim
        return found

    def live_claims(self) -> tuple[QuantitativeClaim, ...]:
        """Return every identifier's current state, excluding retracted ones.

        Returns:
            The standing claims, in the order they were first asserted.
        """
        seen: dict[str, QuantitativeClaim] = {}
        for entry in self._log:
            seen[entry.claim.claim_id] = entry.claim
        return tuple(claim for claim in seen.values() if not claim.is_retracted)

    def quotable_claims(self) -> tuple[QuantitativeClaim, ...]:
        """Return every claim that may currently be quoted externally."""
        return tuple(claim for claim in self.live_claims() if claim.is_quotable)

    @property
    def retracted_digests(self) -> frozenset[str]:
        """Every artifact digest that has ever supported a retracted claim.

        These digests are *spent*: they may still be cited as context, but they can no longer buy a
        QUOTABLE status for any claim in this register.
        """
        spent = {
            entry.claim.supporting_artifact_digest
            for entry in self._log
            if entry.claim.is_retracted and entry.claim.supporting_artifact_digest is not None
        }
        return frozenset(spent)

    def register(self, claim: QuantitativeClaim, *, note: str | None = None) -> RegisterEntry:
        """Append a newly asserted claim.

        Args:
            claim: The claim to assert.  Must not already carry a retraction, and its identifier
                must be unused.
            note: Optional free text for this log line.

        Returns:
            The appended entry.

        Raises:
            ClaimRegisterError: If the identifier is already in the log, if the claim arrives
                pre-retracted, or if it is QUOTABLE on evidence that a retraction already spent.
        """
        if not isinstance(claim, QuantitativeClaim):
            raise TypeError("register takes a QuantitativeClaim")
        if self.current(claim.claim_id) is not None:
            raise ClaimRegisterError(
                f"claim identifier {claim.claim_id!r} is already in the register; identifiers are "
                "never reused, so a replacement claim needs a new one"
            )
        if claim.is_retracted:
            raise ClaimRegisterError(
                f"claim {claim.claim_id!r} arrived already retracted; assert it first so the log "
                "records what was claimed before it records that it was withdrawn"
            )
        self._reject_spent_evidence(claim)
        return self._append(RegisterEvent.ASSERTED, claim, note)

    def retract(
        self,
        claim_id: str,
        *,
        reason: str,
        retracted_by: str,
        superseded_by: str | None = None,
        retracted_at: str | None = None,
        note: str | None = None,
    ) -> RegisterEntry:
        """Withdraw a claim by appending a retraction record.

        Nothing is deleted.  The asserted form of the claim stays in the log at its original
        sequence number, and a new line is appended carrying the same identifier with the retraction
        attached.

        Args:
            claim_id: The identifier to withdraw.
            reason: Why the number may no longer be stated.
            retracted_by: Who or what withdrew it.
            superseded_by: Identifier of the replacement claim, if one exists.  It does not have to
                be registered yet, so that a supersession can be recorded in the order it happened.
            retracted_at: ISO-8601 timestamp string.
            note: Optional free text for this log line.

        Returns:
            The appended entry.

        Raises:
            ClaimRegisterError: If the identifier is unknown or already retracted.
        """
        standing = self.current(claim_id)
        if standing is None:
            raise ClaimRegisterError(
                f"cannot retract {claim_id!r}: it was never asserted in this register"
            )
        retraction = Retraction(
            reason=reason,
            retracted_by=retracted_by,
            superseded_by=superseded_by,
            retracted_at=retracted_at,
        )
        withdrawn = standing.with_retraction(retraction)
        event = (
            RegisterEvent.SUPERSEDED if superseded_by is not None else RegisterEvent.RETRACTED
        )
        return self._append(event, withdrawn, note)

    def as_json_obj(self) -> dict[str, object]:
        """Return the whole log as a JSON-able object.

        Returns:
            The log plus a summary of which identifiers currently stand and which digests are spent.
        """
        return {
            "schema": "aleph-claim-register@1",
            "entries": [entry.as_json_obj() for entry in self._log],
            "live_claim_ids": [claim.claim_id for claim in self.live_claims()],
            "quotable_claim_ids": [claim.claim_id for claim in self.quotable_claims()],
            "spent_artifact_digests": sorted(self.retracted_digests),
        }

    def _append(
        self, event: RegisterEvent, claim: QuantitativeClaim, note: str | None
    ) -> RegisterEntry:
        """Append one line to the log and return it."""
        entry = RegisterEntry(
            sequence=len(self._log), event=event, claim=claim, note=note
        )
        self._log.append(entry)
        return entry

    def _reject_spent_evidence(self, claim: QuantitativeClaim) -> None:
        """Refuse a QUOTABLE claim standing on evidence that a retraction already spent.

        Args:
            claim: The claim about to be appended.

        Raises:
            ClaimRegisterError: If the claim is QUOTABLE and its supporting digest already backed a
                retracted claim.
        """
        if claim.status is not QuantitativeStatus.QUOTABLE:
            return
        digest = claim.supporting_artifact_digest
        if digest is None or digest not in self.retracted_digests:
            return
        prior = [
            entry.claim.claim_id
            for entry in self._log
            if entry.claim.is_retracted and entry.claim.supporting_artifact_digest == digest
        ]
        raise ClaimRegisterError(
            f"claim {claim.claim_id!r} cannot be QUOTABLE on artifact {digest}: that artifact "
            f"already supported the retracted claim(s) {sorted(set(prior))}. Re-quoting this "
            "number requires new evidence, not a new identifier on the old evidence"
        )


def build_register(claims: Sequence[QuantitativeClaim] | Mapping[str, QuantitativeClaim]) -> (
    ClaimRegister
):
    """Build a register by asserting each claim in order.

    Args:
        claims: Claims to assert, either as a sequence or as an identifier-keyed mapping whose
            insertion order is used.

    Returns:
        The populated register.

    Raises:
        ClaimRegisterError: Propagated from :meth:`ClaimRegister.register`.
    """
    register = ClaimRegister()
    iterable = claims.values() if isinstance(claims, Mapping) else claims
    for claim in iterable:
        register.register(claim)
    return register
