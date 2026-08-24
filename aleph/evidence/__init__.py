"""Provenance, claims, and retractions — Aleph's record of what it is entitled to say.

Two ideas live here and they are kept apart on purpose:

* :mod:`aleph.evidence.ladder` grades a result on two orthogonal axes — how far it executed
  (:class:`~aleph.evidence.ladder.EvidenceRung`) and whether its magnitude may be quoted
  (:class:`~aleph.evidence.ladder.QuantitativeStatus`).  Nothing maps one to the other.
* :mod:`aleph.evidence.claim` holds the claims themselves in an append-only register where
  retraction is a record rather than a deletion.
"""

from aleph.evidence.claim import (
    ClaimRegister,
    ClaimRegisterError,
    QuantitativeClaim,
    RegisterEntry,
    RegisterEvent,
    Retraction,
)
from aleph.evidence.ladder import (
    EvidenceLabel,
    EvidenceRung,
    QuantitativeStatus,
    classify_claim,
    rung_at_least,
)

__all__ = [
    "ClaimRegister",
    "ClaimRegisterError",
    "EvidenceLabel",
    "EvidenceRung",
    "QuantitativeClaim",
    "QuantitativeStatus",
    "RegisterEntry",
    "RegisterEvent",
    "Retraction",
    "classify_claim",
    "rung_at_least",
]
