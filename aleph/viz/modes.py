"""The three viewer modes, and the contract each one is allowed to make.

Manuscript §11 gives the viewer three modes with three different relationships to truth, and the
whole point of this module is that the relationship is a *type* rather than a flag someone remembers
to check.

* ``LIVE`` is best-effort. Frames may drop, the consumer may die, and none of it is evidence. A live
  frame does not know how many frames were lost between it and the previous one at the moment it is
  built, so a rate computed from live frames is not a rate.
* ``EVIDENCE_CAPTURE`` publishes at a cadence declared *before* the run, and every frame carries the
  identifiers a later reader needs in order to find the state it depicts.
* ``REPLAY`` re-reads what evidence capture wrote. It creates no new evidence and it must not be
  allowed to: a replay that could mint frames would be a way to manufacture accepted state that no
  transaction ever accepted.

The predicates below are deliberately not derived from one another. ``publishes_evidence`` and
``is_evidence_grade`` differ exactly on ``REPLAY``, and collapsing them into one property is how a
replay session ends up authorised to write.
"""

from __future__ import annotations

from enum import StrEnum, unique

__all__ = [
    "ModeContractError",
    "ViewerMode",
]


class ModeContractError(RuntimeError):
    """An operation was attempted in a mode whose contract does not permit it."""


@unique
class ViewerMode(StrEnum):
    """Which of §11's three contracts a viewer session is operating under."""

    LIVE = "live"
    """Best-effort observation. Frames may drop; nothing here is evidence."""

    EVIDENCE_CAPTURE = "evidence_capture"
    """Predeclared cadence, full identifiers, exact cost. The only mode that mints evidence."""

    REPLAY = "replay"
    """Re-reads a capture. Evidence-grade because the capture was; may not publish."""

    @property
    def may_drop_frames(self) -> bool:
        """Whether losing a frame is an acceptable outcome rather than a capture failure."""
        return self is ViewerMode.LIVE

    @property
    def publishes_evidence(self) -> bool:
        """Whether this mode is allowed to create a new accepted frame."""
        return self is ViewerMode.EVIDENCE_CAPTURE

    @property
    def is_evidence_grade(self) -> bool:
        """Whether frames read in this mode may be cited.

        ``REPLAY`` is grade-bearing and non-publishing at the same time, which is why this is a
        separate question from :attr:`publishes_evidence`.
        """
        return self is not ViewerMode.LIVE

    @property
    def requires_predeclared_cadence(self) -> bool:
        """Whether a cadence must have been fixed before the first accepted step."""
        return self is ViewerMode.EVIDENCE_CAPTURE

    @property
    def not_evidence_reason(self) -> str | None:
        """Why this mode's output may not be cited, or ``None`` when it may."""
        if self is ViewerMode.LIVE:
            return (
                "LIVE is best-effort: frames are dropped whenever the consumer falls behind, so a "
                "live frame cannot support a rate, a count over time, or any claim that depends on "
                "having seen every accepted step"
            )
        return None

    def require_publishing(self, what: str) -> None:
        """Raise unless this mode may mint a new accepted frame.

        Args:
            what: The operation being attempted, named in the error.

        Raises:
            ModeContractError: When the mode is not ``EVIDENCE_CAPTURE``.
        """
        if not self.publishes_evidence:
            raise ModeContractError(
                f"{what} is not permitted in mode {self.value!r}: only "
                f"{ViewerMode.EVIDENCE_CAPTURE.value!r} may publish an accepted frame. "
                + (self.not_evidence_reason or "REPLAY re-reads evidence; it does not create it.")
            )
