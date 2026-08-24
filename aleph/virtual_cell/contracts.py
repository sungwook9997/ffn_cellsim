"""Representation-neutral contracts shared by the virtual-cell research lanes.

The explicit Warp runtime has a deliberately narrow ``ac.engine.CellState``: six slow-context axes
condition rates and inventories but own no physical state.  That contract remains useful inside the
reference engine, but it cannot be the project-wide cell vocabulary: it fixes one axis set and its
first factory is one cell line.  :class:`CellStateManifest` wraps rather than replaces that runtime
object and can describe arbitrary species, lineages, environments, observations, and missing laws.

The other types prevent the three forms of evidence laundering that become possible once tensor,
surrogate, and image lanes exist:

* an approximate representation must carry its measured error and fallback;
* a surrogate or synthetic image cannot claim native-accepted provenance;
* a sandbox experiment must predeclare controls, a falsifier, held-out interventions, and where a
  failure expands the programme.

No class in this module imports Warp, changes a physics gate, or grants a quantitative claim.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from types import MappingProxyType
from typing import Any

__all__ = [
    "ApproximationErrorReport",
    "ArtifactEnvelope",
    "CellStateManifest",
    "EvidenceSource",
    "RepresentationDescriptor",
    "RepresentationKind",
    "SandboxExperimentCard",
    "UncertaintySource",
    "manifest_from_legacy_cell_state",
]

type JSONScalar = str | int | float | bool | None
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


class RepresentationKind(StrEnum):
    """How an artifact represents cell state or an observation."""

    EXPLICIT_STATE = "explicit-state"
    EMPIRICAL_ENSEMBLE = "empirical-ensemble"
    LOCAL_GAUSSIAN = "local-gaussian"
    GAUSSIAN_MIXTURE = "gaussian-mixture"
    TENSOR_TRAIN = "tensor-train"
    LATENT_DENSITY = "latent-density"
    SURROGATE = "surrogate"
    SYNTHETIC_IMAGE = "synthetic-image"
    SCALAR_PROJECTION = "scalar-projection"


class EvidenceSource(StrEnum):
    """Where an artifact's numbers came from, independent of representation."""

    NATIVE_ACCEPTED = "native-accepted"
    NATIVE_ATTEMPTED = "native-attempted"
    ANALYTIC_ORACLE = "analytic-oracle"
    REDUCED_MODEL = "reduced-model"
    SURROGATE = "surrogate"
    SYNTHETIC = "synthetic"
    PUBLIC_EXPERIMENT = "public-experiment"


class UncertaintySource(StrEnum):
    """Mechanistically distinct uncertainty sources that must not be merged silently."""

    THERMAL = "thermal"
    STOCHASTIC_EVENT = "stochastic-event"
    STRUCTURAL = "structural"
    BIOLOGICAL = "biological"
    MEASUREMENT = "measurement"
    PARAMETER = "parameter"
    NUMERICAL = "numerical"
    REPRESENTATION = "representation"


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _unique_strings(
    values: Sequence[str], *, what: str, allow_empty: bool = True
) -> tuple[str, ...]:
    if isinstance(values, str) or not isinstance(values, Sequence):
        raise TypeError(f"{what} must be a sequence of strings")
    normalized = tuple(_nonempty(value, what=f"{what} entry") for value in values)
    if not allow_empty and not normalized:
        raise ValueError(f"{what} must not be empty")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{what} must not contain duplicates")
    return normalized


def _frozen_scalar_mapping(
    values: Mapping[str, JSONScalar], *, what: str
) -> Mapping[str, JSONScalar]:
    if not isinstance(values, Mapping):
        raise TypeError(f"{what} must be a mapping")
    frozen: dict[str, JSONScalar] = {}
    for key, value in values.items():
        name = _nonempty(key, what=f"{what} key")
        if name in frozen:
            raise ValueError(f"{what} keys must not collide after whitespace normalization")
        if not isinstance(value, (str, int, float, bool, type(None))):
            raise TypeError(f"{what}[{name!r}] must be a JSON scalar")
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError(f"{what}[{name!r}] must be finite")
        frozen[name] = value
    return MappingProxyType(frozen)


@dataclass(frozen=True, slots=True)
class CellStateManifest:
    """A cell-archetype/state/environment manifest above the narrow runtime ``CellState``.

    ``missing_components`` and ``missing_laws`` are first-class outputs.  A new cell archetype can
    therefore enter the registry before the current mechanochemical grammar supports it; the honest
    result is a missing-law verdict rather than silently forcing it into an existing cell-line model.
    """

    manifest_id: str
    species: str
    lineage: str
    identity: Mapping[str, JSONScalar]
    biological_state: Mapping[str, JSONScalar]
    environment: Mapping[str, JSONScalar]
    components: tuple[str, ...]
    connectors: tuple[str, ...]
    observations: tuple[str, ...] = ()
    missing_components: tuple[str, ...] = ()
    missing_laws: tuple[str, ...] = ()
    provenance: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "manifest_id", _nonempty(self.manifest_id, what="manifest_id"))
        object.__setattr__(self, "species", _nonempty(self.species, what="species"))
        object.__setattr__(self, "lineage", _nonempty(self.lineage, what="lineage"))
        object.__setattr__(self, "identity", _frozen_scalar_mapping(self.identity, what="identity"))
        object.__setattr__(
            self,
            "biological_state",
            _frozen_scalar_mapping(self.biological_state, what="biological_state"),
        )
        object.__setattr__(
            self, "environment", _frozen_scalar_mapping(self.environment, what="environment")
        )
        object.__setattr__(
            self,
            "components",
            _unique_strings(self.components, what="components", allow_empty=False),
        )
        object.__setattr__(self, "connectors", _unique_strings(self.connectors, what="connectors"))
        object.__setattr__(
            self, "observations", _unique_strings(self.observations, what="observations")
        )
        object.__setattr__(
            self,
            "missing_components",
            _unique_strings(self.missing_components, what="missing_components"),
        )
        object.__setattr__(
            self, "missing_laws", _unique_strings(self.missing_laws, what="missing_laws")
        )
        object.__setattr__(self, "provenance", _unique_strings(self.provenance, what="provenance"))
        overlap = set(self.components) & set(self.missing_components)
        if overlap:
            raise ValueError(
                "a component cannot be both present and missing: " + ", ".join(sorted(overlap))
            )

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical JSON-able content whose hash identifies this manifest."""
        return {
            "manifest_id": self.manifest_id,
            "species": self.species,
            "lineage": self.lineage,
            "identity": dict(sorted(self.identity.items())),
            "biological_state": dict(sorted(self.biological_state.items())),
            "environment": dict(sorted(self.environment.items())),
            "components": sorted(self.components),
            "connectors": sorted(self.connectors),
            "observations": sorted(self.observations),
            "missing_components": sorted(self.missing_components),
            "missing_laws": sorted(self.missing_laws),
            "provenance": sorted(self.provenance),
        }

    @property
    def manifest_hash(self) -> str:
        """SHA-256 of canonical manifest content; independent of mapping insertion order."""
        payload = json.dumps(
            self.to_dict(), sort_keys=True, separators=(",", ":"), ensure_ascii=False
        ).encode()
        return hashlib.sha256(payload).hexdigest()


_APPROXIMATE_REPRESENTATIONS = frozenset(
    {
        RepresentationKind.LOCAL_GAUSSIAN,
        RepresentationKind.GAUSSIAN_MIXTURE,
        RepresentationKind.TENSOR_TRAIN,
        RepresentationKind.LATENT_DENSITY,
        RepresentationKind.SURROGATE,
    }
)


@dataclass(frozen=True, slots=True)
class ApproximationErrorReport:
    """Content bindings and measured scalar error for one reduced representation."""

    report_sha256: str
    source_blob_sha256: str
    reduced_blob_sha256: str
    metric_id: str
    requested_tolerance: float
    measured_error: float
    projection_report_sha256: str | None = None
    metric_unit: str = "1"
    relative_metric: bool = True

    def __post_init__(self) -> None:
        for field_name in (
            "report_sha256",
            "source_blob_sha256",
            "reduced_blob_sha256",
        ):
            if not _SHA256.fullmatch(getattr(self, field_name)):
                raise ValueError(f"{field_name} must be a lowercase SHA-256 digest")
        if self.projection_report_sha256 is not None and not _SHA256.fullmatch(
            self.projection_report_sha256
        ):
            raise ValueError("projection_report_sha256 must be a lowercase SHA-256 digest")
        object.__setattr__(self, "metric_id", _nonempty(self.metric_id, what="metric_id"))
        object.__setattr__(
            self,
            "metric_unit",
            _nonempty(self.metric_unit, what="metric_unit"),
        )
        if not isinstance(self.relative_metric, bool):
            raise TypeError("relative_metric must be bool")
        for field_name in ("requested_tolerance", "measured_error"):
            value = getattr(self, field_name)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0.0
            ):
                raise ValueError(f"{field_name} must be finite and nonnegative")
        if self.relative_metric and self.requested_tolerance >= 1.0:
            raise ValueError("relative requested_tolerance must be less than one")
        if self.measured_error > self.requested_tolerance:
            raise ValueError("measured_error cannot exceed requested_tolerance")


@dataclass(frozen=True, slots=True)
class RepresentationDescriptor:
    """The representation, uncertainty sources, approximation error, and fallback of an artifact."""

    kind: RepresentationKind
    uncertainty_sources: tuple[UncertaintySource, ...]
    axes: tuple[str, ...]
    approximation_error_bound: float | None = None
    error_metric: str | None = None
    error_report: ApproximationErrorReport | None = None
    fallback: str | None = None
    trust_region: Mapping[str, JSONScalar] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.kind, RepresentationKind):
            raise TypeError("kind must be a RepresentationKind")
        sources = tuple(self.uncertainty_sources)
        if any(not isinstance(source, UncertaintySource) for source in sources):
            raise TypeError("uncertainty_sources must contain UncertaintySource values")
        if len(set(sources)) != len(sources):
            raise ValueError("uncertainty_sources must not contain duplicates")
        object.__setattr__(self, "uncertainty_sources", sources)
        object.__setattr__(
            self, "axes", _unique_strings(self.axes, what="representation axes", allow_empty=False)
        )
        object.__setattr__(
            self, "trust_region", _frozen_scalar_mapping(self.trust_region, what="trust_region")
        )

        if self.kind in _APPROXIMATE_REPRESENTATIONS:
            if (
                self.approximation_error_bound is None
                or not math.isfinite(self.approximation_error_bound)
                or self.approximation_error_bound < 0.0
            ):
                raise ValueError(
                    f"{self.kind.value} needs a finite nonnegative approximation_error_bound"
                )
            _nonempty(self.error_metric, what=f"{self.kind.value} error_metric")
            if not isinstance(self.error_report, ApproximationErrorReport):
                raise ValueError(f"{self.kind.value} needs a measured error_report")
            if self.error_report.metric_id != self.error_metric:
                raise ValueError("error_report metric_id must match error_metric")
            if self.error_report.measured_error != self.approximation_error_bound:
                raise ValueError("approximation_error_bound must equal error_report measured_error")
            _nonempty(self.fallback, what=f"{self.kind.value} fallback")
        elif self.approximation_error_bound is not None:
            if (
                not math.isfinite(self.approximation_error_bound)
                or self.approximation_error_bound < 0
            ):
                raise ValueError("approximation_error_bound must be finite and nonnegative")

        if self.kind is RepresentationKind.SURROGATE and not self.trust_region:
            raise ValueError("a surrogate representation needs an explicit non-empty trust_region")


@dataclass(frozen=True, slots=True)
class ArtifactEnvelope:
    """Cross-lane provenance wrapper; it grants no evidence rung or quantitative claim."""

    artifact_id: str
    cell_state_manifest_hash: str
    evidence_source: EvidenceSource
    representation: RepresentationDescriptor
    source_artifact_ids: tuple[str, ...] = ()
    projection_hash: str | None = None
    observation_provenance_sha256: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "artifact_id", _nonempty(self.artifact_id, what="artifact_id"))
        object.__setattr__(
            self,
            "cell_state_manifest_hash",
            _nonempty(self.cell_state_manifest_hash, what="cell_state_manifest_hash"),
        )
        if not _SHA256.fullmatch(self.cell_state_manifest_hash):
            raise ValueError("cell_state_manifest_hash must be a lowercase SHA-256 digest")
        if not isinstance(self.evidence_source, EvidenceSource):
            raise TypeError("evidence_source must be an EvidenceSource")
        if not isinstance(self.representation, RepresentationDescriptor):
            raise TypeError("representation must be a RepresentationDescriptor")
        object.__setattr__(
            self,
            "source_artifact_ids",
            _unique_strings(self.source_artifact_ids, what="source_artifact_ids"),
        )

        if self.evidence_source in {
            EvidenceSource.NATIVE_ACCEPTED,
            EvidenceSource.NATIVE_ATTEMPTED,
        }:
            raise ValueError(
                "generic ArtifactEnvelope cannot grant native evidence authority; use a future "
                "runtime-attested adapter after external adjudication"
            )
        if (
            self.representation.kind
            in _APPROXIMATE_REPRESENTATIONS
            | {RepresentationKind.SYNTHETIC_IMAGE, RepresentationKind.SCALAR_PROJECTION}
            and not self.source_artifact_ids
        ):
            raise ValueError(
                "derived, approximate, and synthetic artifacts require source_artifact_ids ancestry"
            )
        if (
            self.representation.kind is RepresentationKind.SURROGATE
            and self.evidence_source is not EvidenceSource.SURROGATE
        ):
            raise ValueError("a surrogate representation must retain SURROGATE evidence provenance")
        if (
            self.representation.kind is RepresentationKind.SYNTHETIC_IMAGE
            and self.evidence_source is not EvidenceSource.SYNTHETIC
        ):
            raise ValueError(
                "a synthetic-image representation must retain SYNTHETIC evidence provenance"
            )
        if (
            self.representation.kind
            in _APPROXIMATE_REPRESENTATIONS - {RepresentationKind.SURROGATE}
            and self.evidence_source is EvidenceSource.NATIVE_ACCEPTED
        ):
            raise ValueError(
                "an approximate representation is REDUCED_MODEL evidence; point back to the "
                "native artifact through source_artifact_ids instead of laundering its provenance"
            )
        if self.representation.kind is RepresentationKind.SCALAR_PROJECTION:
            if self.projection_hash is None or not _SHA256.fullmatch(self.projection_hash):
                raise ValueError("scalar projection_hash must be a lowercase SHA-256 digest")
            if not self.source_artifact_ids:
                raise ValueError("a scalar projection needs at least one source artifact")
        if self.representation.kind is RepresentationKind.SYNTHETIC_IMAGE and (
            self.observation_provenance_sha256 is None
            or not _SHA256.fullmatch(self.observation_provenance_sha256)
        ):
            raise ValueError("a synthetic image needs an observation_provenance_sha256 binding")


@dataclass(frozen=True, slots=True)
class SandboxExperimentCard:
    """Predeclared simulation experiment with controls, falsifier, and expansion routes."""

    experiment_id: str
    question: str
    manifest_hashes: tuple[str, ...]
    intervention: str
    observables: tuple[str, ...]
    positive_controls: tuple[str, ...]
    negative_controls: tuple[str, ...]
    adversarial_controls: tuple[str, ...]
    held_out_interventions: tuple[str, ...]
    falsifier: str
    failure_expansions: tuple[str, ...]
    budget_class: str
    wet_lab: bool = False
    information_gain_rationale: str | None = None

    def __post_init__(self) -> None:
        for field_name in (
            "experiment_id",
            "question",
            "intervention",
            "falsifier",
            "budget_class",
        ):
            object.__setattr__(
                self, field_name, _nonempty(getattr(self, field_name), what=field_name)
            )
        for field_name in (
            "manifest_hashes",
            "observables",
            "positive_controls",
            "negative_controls",
            "adversarial_controls",
            "held_out_interventions",
            "failure_expansions",
        ):
            object.__setattr__(
                self,
                field_name,
                _unique_strings(getattr(self, field_name), what=field_name, allow_empty=False),
            )
        if self.wet_lab:
            _nonempty(
                self.information_gain_rationale,
                what="wet-lab information_gain_rationale",
            )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able pre-registration record."""
        return {
            "experiment_id": self.experiment_id,
            "question": self.question,
            "manifest_hashes": list(self.manifest_hashes),
            "intervention": self.intervention,
            "observables": list(self.observables),
            "positive_controls": list(self.positive_controls),
            "negative_controls": list(self.negative_controls),
            "adversarial_controls": list(self.adversarial_controls),
            "held_out_interventions": list(self.held_out_interventions),
            "falsifier": self.falsifier,
            "failure_expansions": list(self.failure_expansions),
            "budget_class": self.budget_class,
            "wet_lab": self.wet_lab,
            "information_gain_rationale": self.information_gain_rationale,
        }

    @property
    def card_hash(self) -> str:
        """SHA-256 of the canonical pre-registration record."""
        payload = json.dumps(
            self.to_dict(), sort_keys=True, separators=(",", ":"), ensure_ascii=False
        ).encode()
        return hashlib.sha256(payload).hexdigest()


def manifest_from_legacy_cell_state(
    state: object,
    *,
    manifest_id: str,
    species: str,
    components: Sequence[str],
    connectors: Sequence[str],
    observations: Sequence[str] = (),
    provenance: Sequence[str] = (),
) -> CellStateManifest:
    """Wrap a six-axis runtime ``CellState`` without making those axes universal.

    The adapter is structural and intentionally imports no ``aleph.ac`` type.  This lets the broad
    virtual-cell layer consume the incumbent runtime today while future runtimes expose different
    state axes.
    """
    axes = getattr(state, "axes", None)
    if axes is None:
        raise TypeError("legacy cell state must expose an axes sequence")
    values = tuple(axes)
    if len(values) != 6 or any(not isinstance(value, str) or not value.strip() for value in values):
        raise ValueError("legacy cell state axes must be six non-empty strings")
    lineage, emt, cycle, adhesion, geometry, osmotic_state = values
    return CellStateManifest(
        manifest_id=manifest_id,
        species=species,
        lineage=lineage,
        identity={"legacy_runtime_schema": "ac.engine.CellState@6-axis"},
        biological_state={"emt": emt, "cycle": cycle},
        environment={
            "adhesion": adhesion,
            "geometry": geometry,
            "osmotic_state": osmotic_state,
        },
        components=tuple(components),
        connectors=tuple(connectors),
        observations=tuple(observations),
        provenance=tuple(provenance),
    )
