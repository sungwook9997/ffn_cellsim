"""Deterministic identifiers for proposed external experiment annotations."""

from __future__ import annotations

import hashlib


def _identifier(namespace: str, *parts: object) -> str:
    material = "\x1f".join(str(part).strip() for part in parts)
    return f"aleph:{namespace}:{hashlib.sha256(material.encode('utf-8')).hexdigest()[:16]}"


def experiment_id(source_family_id: str, evidence_locator_key: str, local_ordinal: int) -> str:
    return _identifier("experiment", source_family_id.lower(), evidence_locator_key, local_ordinal)


def observation_id(
    experiment_record_id: str,
    observation_type: str,
    evidence_locator_key: str,
    local_ordinal: int,
) -> str:
    return _identifier(
        "observation",
        experiment_record_id,
        observation_type,
        evidence_locator_key,
        local_ordinal,
    )
