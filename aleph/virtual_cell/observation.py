"""Composite validation for synthetic observation artifacts."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

from aleph.virtual_cell.contracts import (
    ArtifactEnvelope,
    CellStateManifest,
    EvidenceSource,
    RepresentationKind,
)
from aleph.virtual_cell.sidecar import ArraySidecarDescriptor
from aleph.virtual_cell.synthetic_microscopy import SyntheticObservationProvenance

__all__ = ["SyntheticObservationArtifact"]


@dataclass(frozen=True, slots=True)
class SyntheticObservationArtifact:
    """One hash-connected source/manifest/config/sidecar unit; still synthetic evidence only."""

    envelope: ArtifactEnvelope
    sidecar: ArraySidecarDescriptor
    provenance: SyntheticObservationProvenance
    manifest: CellStateManifest

    def __post_init__(self) -> None:
        if self.envelope.evidence_source is not EvidenceSource.SYNTHETIC:
            raise ValueError("synthetic observation artifact must retain SYNTHETIC evidence")
        if self.envelope.representation.kind is not RepresentationKind.SYNTHETIC_IMAGE:
            raise ValueError("synthetic observation artifact needs SYNTHETIC_IMAGE representation")
        if self.envelope.cell_state_manifest_hash != self.manifest.manifest_hash:
            raise ValueError("envelope manifest hash does not match supplied manifest")
        if self.envelope.observation_provenance_sha256 != self.provenance.record_hash:
            raise ValueError("envelope does not bind the supplied observation provenance")
        if self.provenance.observed_sidecar_sha256 != self.sidecar.sha256:
            raise ValueError("observation provenance does not bind the supplied sidecar")
        representation_axes = self.envelope.representation.axes
        sidecar_axes = tuple(axis.name for axis in self.sidecar.axes)
        if representation_axes != sidecar_axes:
            raise ValueError("envelope representation axes do not match sidecar axes")

    @property
    def artifact_hash(self) -> str:
        payload = json.dumps(
            {
                "artifact_id": self.envelope.artifact_id,
                "manifest_hash": self.manifest.manifest_hash,
                "sidecar_sha256": self.sidecar.sha256,
                "provenance_sha256": self.provenance.record_hash,
                "axes": [axis.name for axis in self.sidecar.axes],
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
        return hashlib.sha256(payload).hexdigest()
