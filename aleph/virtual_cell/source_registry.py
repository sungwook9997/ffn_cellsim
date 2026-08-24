"""Provenance DAG for evidence sources, with structural exclusion of self-scoring evidence.

WHY this module exists.  Master plan section 10.2 attack A6: if an encoder is trained on synthetic
images and then SCORED on synthetic images from the same renderer family, it learns renderer
artifacts and self-scores perfectly.  The plan's own verdict is that such a test is a structure test
and not a biological validation.  A convention cannot enforce that, because the laundering happens
one hop at a time: a renderer output becomes a training corpus, the training corpus becomes a fitted
statistic, and the fitted statistic gets quoted as if it described cells.  Each hop looks innocent.
So the exclusion is TRANSITIVE and structural here: a source is citable as biological evidence only
if EVERY node in its ancestry is empirically sourced and in good standing.  Failing that, the query
returns the offending ancestor and the path to it, not a bare boolean.

WHY a DAG, and why cycles must be refused.  Evidence provenance is acyclic by definition: a claim
cannot be an ancestor of its own support.  A cycle in a provenance graph is exactly the shape of
circular justification, and it silently defeats the transitive exclusion above because a naive walk
either loops forever or terminates early and declares the corpus clean.  :meth:`SourceRegistry.freeze`
refuses cycles and reports the offending cycle as a PATH so the circular argument can be read off.

WHY ``CALIBRATION_ONLY`` is separately gated.  A calibration set is real biological data, so it
passes the corpus query; but a model tuned on it cannot be validated against it.  That restriction is
also transitive: deriving a "held-out" split from a calibration set leaks the calibration.

Provenance kinds are the master plan's section 2.6 set.  Everything not empirically sourced is
excluded from the biological-truth corpus by a WHITELIST of origins, so a new origin added later
fails closed rather than being silently admitted.

Public microscopy datasets: this module defines the ADAPTER CONTRACT that a dataset must satisfy to
be usable as an observation-validation source.  It downloads nothing, and it claims no integration.
Every real dataset in :data:`BLOCKED_PUBLIC_MICROSCOPY_DATASETS` is recorded as
``blocked / not yet integrated``; the only dataset that satisfies the contract is a worked in-memory
FIXTURE, which is itself excluded from the biological corpus because a fixture is not evidence.

Sanity Gate (self-tested in tests/ac/virtual_cell/test_optics_sources.py):
  * boundary: self-loop, two-cycle, and long-cycle provenance graphs are all refused, each with the
    cycle reported as a path; a dangling parent reference is refused by name.
  * conservation-invariant (the one that matters): the citable-as-biological-evidence set is closed
    under ancestry.  If a node is excluded, every descendant is excluded, at any depth, with the
    blocking ancestor and path reported.  A source two hops downstream of a synthetic renderer output
    must be excluded — that is the A6 failure.
  * measurement protocol: held-out validation admissibility requires an empirical ancestry, a
    declared split identifier, and no ``CALIBRATION_ONLY`` anywhere in the ancestry.
  * dimensional/typing: dataset contracts declare pixel size, axial step, exposure, bit depth, and a
    channel-to-marker map; a contract claiming ``INTEGRATED`` while missing any requirement is refused
    at construction.
  * fail closed: origins are admitted by whitelist, and the default status of a real public dataset is
    ``BLOCKED``.

Nothing here grants an evidence rung.  It only refuses ones that were never earned.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from types import MappingProxyType

__all__ = [
    "BLOCKED_PUBLIC_MICROSCOPY_DATASETS",
    "DECLARED_EVIDENCE_GAPS",
    "REQUIRED_DATASET_REQUIREMENTS",
    "BiologicalEvidenceExclusion",
    "DatasetIntegrationStatus",
    "EvidenceGapVerdict",
    "EvidenceGraph",
    "ExclusionReason",
    "HeldOutValidationRefused",
    "MicroscopyDatasetContract",
    "ProvenanceKind",
    "SourceCycleError",
    "SourceOrigin",
    "SourceRecord",
    "SourceRegistry",
    "SourceStatus",
    "UnknownParentError",
    "register_public_dataset_contracts",
    "worked_fixture_dataset",
]

_IDENTIFIER = re.compile(r"^[a-z0-9][a-z0-9._:-]*$")


class ProvenanceKind(StrEnum):
    """Master plan section 2.6 provenance kinds for a parameter, claim, or corpus."""

    SOURCED = "SOURCED"
    DERIVED = "DERIVED"
    COMPUTED_QM_MD = "COMPUTED_QM_MD"
    INFERRED_POSTERIOR = "INFERRED_POSTERIOR"
    CALIBRATION_ONLY = "CALIBRATION_ONLY"
    CONVENIENCE = "CONVENIENCE"
    PI_GAP = "PI_GAP"


class SourceOrigin(StrEnum):
    """Where the numbers physically came from, independent of provenance kind."""

    LITERATURE = "literature"
    PUBLIC_EXPERIMENT = "public-experiment"
    IN_HOUSE_EXPERIMENT = "in-house-experiment"
    SYNTHETIC_RENDER = "synthetic-render"
    SIMULATION = "simulation"
    COMPUTATION = "computation"
    EXPERT_JUDGEMENT = "expert-judgement"


# Whitelist, deliberately.  An origin added later is NOT biological evidence until someone argues it
# into this set, so the failure mode of forgetting to update this module is over-exclusion.
_EMPIRICAL_ORIGINS = frozenset(
    {
        SourceOrigin.LITERATURE,
        SourceOrigin.PUBLIC_EXPERIMENT,
        SourceOrigin.IN_HOUSE_EXPERIMENT,
    }
)

_PLACEHOLDER_PROVENANCE = frozenset({ProvenanceKind.CONVENIENCE, ProvenanceKind.PI_GAP})


class SourceStatus(StrEnum):
    """Standing of a source in the corpus."""

    ACTIVE = "active"
    REJECTED = "rejected"
    RETRACTED = "retracted"
    BLOCKED = "blocked"
    FIXTURE = "fixture"


class ExclusionReason(StrEnum):
    """Why a source may not be cited as biological evidence."""

    NON_EMPIRICAL_ORIGIN = "non-empirical-origin"
    REJECTED = "rejected"
    RETRACTED = "retracted"
    BLOCKED = "blocked"
    FIXTURE = "contract-fixture"
    PLACEHOLDER_PROVENANCE = "placeholder-provenance"
    DERIVED_FROM_EXCLUDED = "derived-from-excluded"


class DatasetIntegrationStatus(StrEnum):
    """How far a public dataset has actually been wired in.  Nothing here downloads data."""

    BLOCKED_NOT_INTEGRATED = "blocked/not-yet-integrated"
    IN_MEMORY_FIXTURE = "in-memory-fixture"
    INTEGRATED = "integrated"


class SourceCycleError(ValueError):
    """A provenance graph contained a cycle; the cycle is reported as a path."""

    def __init__(self, cycle: Sequence[str]) -> None:
        self.cycle: tuple[str, ...] = tuple(cycle)
        super().__init__("provenance cycle refused: " + " -> ".join(self.cycle))


class UnknownParentError(ValueError):
    """A source named a parent that was never registered."""

    def __init__(self, source_id: str, parent_id: str) -> None:
        self.source_id = source_id
        self.parent_id = parent_id
        super().__init__(f"source {source_id!r} names unregistered parent {parent_id!r}")


class HeldOutValidationRefused(ValueError):
    """A source may not be used as held-out validation; the reason is machine-readable."""

    def __init__(self, source_id: str, reason: str) -> None:
        self.source_id = source_id
        self.reason = reason
        super().__init__(f"{source_id!r} is not admissible as held-out validation: {reason}")


def _identifier(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    normalized = value.strip()
    if not _IDENTIFIER.fullmatch(normalized):
        raise ValueError(
            f"{what} must match {_IDENTIFIER.pattern} (lowercase, no spaces): {value!r}"
        )
    return normalized


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _optional_nonempty(value: object, *, what: str) -> str | None:
    if value is None:
        return None
    return _nonempty(value, what=what)


def _optional_positive(value: object, *, what: str) -> float | None:
    if value is None:
        return None
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0.0
    ):
        raise ValueError(f"{what} must be positive-finite when supplied")
    return float(value)


def _unique_identifiers(values: object, *, what: str) -> tuple[str, ...]:
    if isinstance(values, str) or not isinstance(values, Sequence):
        raise TypeError(f"{what} must be a sequence of identifiers")
    normalized = tuple(_identifier(value, what=f"{what} entry") for value in values)
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{what} must not contain duplicates")
    return normalized


@dataclass(frozen=True, slots=True)
class EvidenceGapVerdict:
    """A wall this lane hit, typed so it routes to an expansion instead of shrinking scope."""

    gap_id: str
    failure: str
    status: str
    routes_to: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "gap_id", _identifier(self.gap_id, what="gap_id"))
        for name in ("failure", "routes_to"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        if self.status not in {"blocked", "oracle-only", "pi-decision"}:
            raise ValueError("status must be 'blocked', 'oracle-only', or 'pi-decision'")


DECLARED_EVIDENCE_GAPS: tuple[EvidenceGapVerdict, ...] = (
    EvidenceGapVerdict(
        gap_id="public-microscopy-not-integrated",
        failure=(
            "no public microscopy dataset is wired in, so every image comparison available today is "
            "renderer-internal and cannot falsify the observation model"
        ),
        status="blocked",
        routes_to=(
            "sim-data mismatch -> observation-model extension: land a dataset adapter that satisfies "
            "MicroscopyDatasetContract, then split real-image domain shift from model inadequacy"
        ),
    ),
    EvidenceGapVerdict(
        gap_id="single-renderer-family",
        failure=(
            "one renderer family means a passing image test is a structure test; A6 self-scoring is "
            "not excluded by accuracy, only by provenance"
        ),
        status="blocked",
        routes_to=(
            "image nonidentifiability -> renderer holdout plus a second independent rendering engine, "
            "and adversarial artifact controls"
        ),
    ),
    EvidenceGapVerdict(
        gap_id="calibration-holdout-partition",
        failure=(
            "which real sources are calibration and which are held out is a programme decision, not "
            "a property this module can infer"
        ),
        status="pi-decision",
        routes_to=(
            "PI declares held_out_split_id per source; the registry only enforces that a calibration "
            "ancestry can never appear on the held-out side"
        ),
    ),
)


@dataclass(frozen=True, slots=True)
class SourceRecord:
    """One node of the provenance DAG.

    ``parents`` names the sources this one was derived from.  Edges point child -> parent, and the
    graph must be acyclic.
    """

    source_id: str
    provenance: ProvenanceKind
    origin: SourceOrigin
    status: SourceStatus
    description: str
    parents: tuple[str, ...] = ()
    citation: str | None = None
    held_out_split_id: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "source_id", _identifier(self.source_id, what="source_id"))
        for name, enum_type in (
            ("provenance", ProvenanceKind),
            ("origin", SourceOrigin),
            ("status", SourceStatus),
        ):
            if not isinstance(getattr(self, name), enum_type):
                raise TypeError(f"{name} must be a {enum_type.__name__}")
        object.__setattr__(self, "description", _nonempty(self.description, what="description"))
        object.__setattr__(self, "parents", _unique_identifiers(self.parents, what="parents"))
        object.__setattr__(self, "citation", _optional_nonempty(self.citation, what="citation"))
        object.__setattr__(
            self,
            "held_out_split_id",
            None
            if self.held_out_split_id is None
            else _identifier(self.held_out_split_id, what="held_out_split_id"),
        )
        if self.provenance is ProvenanceKind.SOURCED and self.parents:
            raise ValueError(
                f"{self.source_id!r}: SOURCED is a primary source and must not name parents"
            )
        if self.provenance is ProvenanceKind.PI_GAP and self.parents:
            raise ValueError(
                f"{self.source_id!r}: PI_GAP is a declared absence of evidence and has no parents"
            )
        if (
            self.provenance in {ProvenanceKind.DERIVED, ProvenanceKind.INFERRED_POSTERIOR}
            and not self.parents
        ):
            raise ValueError(
                f"{self.source_id!r}: {self.provenance.value} must name at least one parent"
            )
        if self.origin in _EMPIRICAL_ORIGINS and self.citation is None:
            raise ValueError(
                f"{self.source_id!r}: an empirical origin must carry a citation or accession"
            )

    def to_dict(self) -> dict[str, object]:
        """Return the canonical JSON-able record content."""
        return {
            "source_id": self.source_id,
            "provenance": self.provenance.value,
            "origin": self.origin.value,
            "status": self.status.value,
            "description": self.description,
            "parents": sorted(self.parents),
            "citation": self.citation,
            "held_out_split_id": self.held_out_split_id,
        }


@dataclass(frozen=True, slots=True)
class BiologicalEvidenceExclusion:
    """Why one source may not be cited as biological evidence, and which ancestor blocks it."""

    source_id: str
    reason: ExclusionReason
    blocking_source_id: str
    path: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "source_id", _identifier(self.source_id, what="source_id"))
        object.__setattr__(
            self,
            "blocking_source_id",
            _identifier(self.blocking_source_id, what="blocking_source_id"),
        )
        if not isinstance(self.reason, ExclusionReason):
            raise TypeError("reason must be an ExclusionReason")
        path = tuple(self.path)
        if not path or path[0] != self.source_id or path[-1] != self.blocking_source_id:
            raise ValueError("path must run from source_id to blocking_source_id inclusive")
        object.__setattr__(self, "path", path)

    def __str__(self) -> str:
        return f"{self.source_id}: {self.reason.value} via " + " -> ".join(self.path)


def _self_exclusion_reason(record: SourceRecord) -> ExclusionReason | None:
    """Reason this record is excluded on its own account, ignoring ancestry."""
    status_reasons = {
        SourceStatus.REJECTED: ExclusionReason.REJECTED,
        SourceStatus.RETRACTED: ExclusionReason.RETRACTED,
        SourceStatus.BLOCKED: ExclusionReason.BLOCKED,
        SourceStatus.FIXTURE: ExclusionReason.FIXTURE,
    }
    if record.status in status_reasons:
        return status_reasons[record.status]
    if record.origin not in _EMPIRICAL_ORIGINS:
        return ExclusionReason.NON_EMPIRICAL_ORIGIN
    if record.provenance in _PLACEHOLDER_PROVENANCE:
        return ExclusionReason.PLACEHOLDER_PROVENANCE
    return None


def _canonical_cycle(path: Sequence[str]) -> tuple[str, ...]:
    """Rotate a closed walk so it starts at its lexicographically smallest node."""
    nodes = list(path[:-1])
    pivot = min(range(len(nodes)), key=lambda index: nodes[index])
    rotated = nodes[pivot:] + nodes[:pivot]
    return (*rotated, rotated[0])


def _topological_order(records: Mapping[str, SourceRecord]) -> tuple[str, ...]:
    """Return parents-before-children order, refusing dangling parents and cycles."""
    white, grey, black = 0, 1, 2
    colour = dict.fromkeys(records, white)
    order: list[str] = []
    for start in sorted(records):
        if colour[start] != white:
            continue
        colour[start] = grey
        path = [start]
        stack: list[tuple[str, Iterable[str]]] = [(start, iter(records[start].parents))]
        while stack:
            node, parent_iter = stack[-1]
            descended = False
            for parent in parent_iter:
                if parent not in records:
                    raise UnknownParentError(node, parent)
                if colour[parent] == grey:
                    raise SourceCycleError(_canonical_cycle([*path[path.index(parent) :], parent]))
                if colour[parent] == white:
                    colour[parent] = grey
                    path.append(parent)
                    stack.append((parent, iter(records[parent].parents)))
                    descended = True
                    break
            if not descended:
                colour[node] = black
                order.append(node)
                stack.pop()
                path.pop()
    return tuple(order)


@dataclass(frozen=True, slots=True)
class EvidenceGraph:
    """A validated, acyclic provenance graph that answers citation-admissibility questions."""

    records: Mapping[str, SourceRecord]
    topological_order: tuple[str, ...]
    _exclusions: Mapping[str, BiologicalEvidenceExclusion]
    _calibration_ancestry: Mapping[str, tuple[str, ...]]

    def __contains__(self, source_id: object) -> bool:
        return source_id in self.records

    def __len__(self) -> int:
        return len(self.records)

    def record(self, source_id: str) -> SourceRecord:
        """Return one record.

        Raises:
            KeyError: If the source is not registered.
        """
        return self.records[source_id]

    def ancestors(self, source_id: str) -> frozenset[str]:
        """Return every source reachable through parent edges, excluding ``source_id`` itself."""
        if source_id not in self.records:
            raise KeyError(source_id)
        found: set[str] = set()
        frontier = list(self.records[source_id].parents)
        while frontier:
            current = frontier.pop()
            if current in found:
                continue
            found.add(current)
            frontier.extend(self.records[current].parents)
        return frozenset(found)

    def exclusion(self, source_id: str) -> BiologicalEvidenceExclusion | None:
        """Return why the source is excluded from the biological corpus, or ``None`` if citable."""
        if source_id not in self.records:
            raise KeyError(source_id)
        return self._exclusions.get(source_id)

    def citable_as_biological_evidence(self) -> tuple[str, ...]:
        """Return the sorted sources that may be cited as biological evidence.

        A source qualifies only if it and every ancestor is empirically sourced, in good standing,
        and not a placeholder.  This is what structurally refuses the A6 synthetic self-scoring loop.
        """
        return tuple(
            source_id for source_id in sorted(self.records) if source_id not in self._exclusions
        )

    def excluded_from_biological_evidence(self) -> Mapping[str, BiologicalEvidenceExclusion]:
        """Return every exclusion, keyed by source id."""
        return self._exclusions

    def calibration_ancestry(self, source_id: str) -> tuple[str, ...]:
        """Return the ``CALIBRATION_ONLY`` sources in this source's closure, itself included."""
        if source_id not in self.records:
            raise KeyError(source_id)
        return self._calibration_ancestry[source_id]

    def assert_admissible_as_held_out_validation(self, source_id: str) -> None:
        """Refuse a source that cannot serve as held-out validation.

        Args:
            source_id: Source to check.

        Raises:
            KeyError: If the source is not registered.
            HeldOutValidationRefused: If the source is not citable as biological evidence, carries a
                ``CALIBRATION_ONLY`` source anywhere in its closure, or declares no split identifier.
        """
        if source_id not in self.records:
            raise KeyError(source_id)
        excluded = self._exclusions.get(source_id)
        if excluded is not None:
            raise HeldOutValidationRefused(source_id, str(excluded))
        calibration = self._calibration_ancestry[source_id]
        if calibration:
            raise HeldOutValidationRefused(
                source_id,
                "CALIBRATION_ONLY may not be used as held-out validation; calibration closure: "
                + ", ".join(calibration),
            )
        if self.records[source_id].held_out_split_id is None:
            raise HeldOutValidationRefused(
                source_id, "no held_out_split_id declared, so the split cannot be audited"
            )

    def admissible_held_out_validation_ids(self) -> tuple[str, ...]:
        """Return the sorted sources that pass every held-out validation requirement."""
        admissible: list[str] = []
        for source_id in sorted(self.records):
            try:
                self.assert_admissible_as_held_out_validation(source_id)
            except HeldOutValidationRefused:
                continue
            admissible.append(source_id)
        return tuple(admissible)

    @property
    def graph_hash(self) -> str:
        """SHA-256 over the canonical content of every record, independent of insertion order."""
        payload = json.dumps(
            [self.records[source_id].to_dict() for source_id in sorted(self.records)],
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        ).encode()
        return hashlib.sha256(payload).hexdigest()


class SourceRegistry:
    """Mutable builder for a provenance DAG; :meth:`freeze` validates and returns the graph.

    Parents may be named before they are registered, so a cyclic or dangling graph can be built and
    then REFUSED with a readable diagnosis rather than being impossible to express.
    """

    __slots__ = ("_records",)

    def __init__(self) -> None:
        self._records: dict[str, SourceRecord] = {}

    def register(self, record: SourceRecord) -> SourceRegistry:
        """Add one record.

        Args:
            record: Source record to add.

        Returns:
            This registry, for chaining.

        Raises:
            TypeError: If ``record`` is not a ``SourceRecord``.
            ValueError: If the source id is already registered.
        """
        if not isinstance(record, SourceRecord):
            raise TypeError("record must be a SourceRecord")
        if record.source_id in self._records:
            raise ValueError(f"source {record.source_id!r} is already registered")
        self._records[record.source_id] = record
        return self

    def register_all(self, records: Iterable[SourceRecord]) -> SourceRegistry:
        """Add several records in order."""
        for record in records:
            self.register(record)
        return self

    def __contains__(self, source_id: object) -> bool:
        return source_id in self._records

    def __len__(self) -> int:
        return len(self._records)

    def freeze(self) -> EvidenceGraph:
        """Validate the graph and return an immutable queryable view.

        Returns:
            The validated :class:`EvidenceGraph`.

        Raises:
            UnknownParentError: If any source names an unregistered parent.
            SourceCycleError: If the parent edges contain a cycle, including a self-loop.
        """
        records = MappingProxyType(dict(self._records))
        order = _topological_order(records)

        exclusions: dict[str, BiologicalEvidenceExclusion] = {}
        calibration: dict[str, tuple[str, ...]] = {}
        for source_id in order:
            record = records[source_id]
            own_calibration = set()
            if record.provenance is ProvenanceKind.CALIBRATION_ONLY:
                own_calibration.add(source_id)
            for parent in record.parents:
                own_calibration.update(calibration[parent])
            calibration[source_id] = tuple(sorted(own_calibration))

            reason = _self_exclusion_reason(record)
            if reason is not None:
                exclusions[source_id] = BiologicalEvidenceExclusion(
                    source_id=source_id,
                    reason=reason,
                    blocking_source_id=source_id,
                    path=(source_id,),
                )
                continue
            for parent in record.parents:
                inherited = exclusions.get(parent)
                if inherited is None:
                    continue
                exclusions[source_id] = BiologicalEvidenceExclusion(
                    source_id=source_id,
                    reason=ExclusionReason.DERIVED_FROM_EXCLUDED,
                    blocking_source_id=inherited.blocking_source_id,
                    path=(source_id, *inherited.path),
                )
                break

        return EvidenceGraph(
            records=records,
            topological_order=order,
            _exclusions=MappingProxyType(exclusions),
            _calibration_ancestry=MappingProxyType(calibration),
        )


# --------------------------------------------------------------------------------------------------
# Public microscopy dataset adapter contract.  Nothing below downloads or integrates any dataset.
# --------------------------------------------------------------------------------------------------

REQUIRED_DATASET_REQUIREMENTS: tuple[str, ...] = (
    "citation",
    "licence",
    "accession",
    "modality",
    "channels",
    "channel_marker_map",
    "pixel_size_um",
    "exposure_s",
    "bit_depth",
    "psf_calibration",
    "background_calibration",
    "species",
    "cell_line",
    "biological_state",
)


@dataclass(frozen=True, slots=True)
class MicroscopyDatasetContract:
    """What a public microscopy dataset must supply to serve as an observation-validation source.

    The contract is deliberately about the OBSERVATION MODEL, not about pixels: without pixel size,
    exposure, bit depth, a channel-to-marker map, and a PSF/background calibration, a real image
    cannot be compared to a rendered one at all, and any agreement or disagreement would be
    uninterpretable.  ``axial_step_um`` is optional because a single-plane dataset is still usable.
    """

    dataset_id: str
    integration_status: DatasetIntegrationStatus
    modality: str | None = None
    citation: str | None = None
    licence: str | None = None
    accession: str | None = None
    channels: tuple[str, ...] = ()
    channel_marker_map: Mapping[str, str] = MappingProxyType({})
    pixel_size_um: float | None = None
    axial_step_um: float | None = None
    exposure_s: float | None = None
    bit_depth: int | None = None
    psf_calibration: str | None = None
    background_calibration: str | None = None
    species: str | None = None
    cell_line: str | None = None
    biological_state: str | None = None
    perturbations: tuple[str, ...] = ()
    held_out_split_id: str | None = None
    notes: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "dataset_id", _identifier(self.dataset_id, what="dataset_id"))
        if not isinstance(self.integration_status, DatasetIntegrationStatus):
            raise TypeError("integration_status must be a DatasetIntegrationStatus")
        for name in (
            "modality",
            "citation",
            "licence",
            "accession",
            "psf_calibration",
            "background_calibration",
            "species",
            "cell_line",
            "biological_state",
            "notes",
        ):
            object.__setattr__(self, name, _optional_nonempty(getattr(self, name), what=name))
        object.__setattr__(self, "channels", _unique_identifiers(self.channels, what="channels"))
        object.__setattr__(
            self, "perturbations", _unique_identifiers(self.perturbations, what="perturbations")
        )
        if not isinstance(self.channel_marker_map, Mapping):
            raise TypeError("channel_marker_map must be a mapping")
        marker_map = {
            _identifier(key, what="channel_marker_map key"): _nonempty(
                value, what="channel_marker_map value"
            )
            for key, value in self.channel_marker_map.items()
        }
        if marker_map and set(marker_map) != set(self.channels):
            raise ValueError("channel_marker_map must cover exactly the declared channels")
        object.__setattr__(self, "channel_marker_map", MappingProxyType(marker_map))
        for name in ("pixel_size_um", "axial_step_um", "exposure_s"):
            object.__setattr__(self, name, _optional_positive(getattr(self, name), what=name))
        if self.bit_depth is not None and (
            isinstance(self.bit_depth, bool)
            or not isinstance(self.bit_depth, int)
            or not 1 <= self.bit_depth <= 32
        ):
            raise ValueError("bit_depth must be an integer in [1, 32] when supplied")
        object.__setattr__(
            self,
            "held_out_split_id",
            None
            if self.held_out_split_id is None
            else _identifier(self.held_out_split_id, what="held_out_split_id"),
        )
        if (
            self.integration_status is DatasetIntegrationStatus.INTEGRATED
            and self.missing_requirements()
        ):
            raise ValueError(
                f"{self.dataset_id!r} claims INTEGRATED but is missing: "
                + ", ".join(self.missing_requirements())
            )

    def missing_requirements(self) -> tuple[str, ...]:
        """Return the contract requirements this dataset has not supplied, in declared order."""
        missing = []
        for requirement in REQUIRED_DATASET_REQUIREMENTS:
            value = getattr(self, requirement)
            if value is None or (isinstance(value, (tuple, Mapping)) and len(value) == 0):
                missing.append(requirement)
        return tuple(missing)

    def satisfies_contract(self) -> bool:
        """Return whether every required field is supplied.  This is NOT a claim of integration."""
        return not self.missing_requirements()

    def to_source_record(self) -> SourceRecord:
        """Return the provenance record this dataset would occupy.

        A blocked dataset becomes a ``BLOCKED`` source and an in-memory fixture becomes a ``FIXTURE``
        source; both are structurally excluded from the biological corpus, so nothing downstream of
        them can be cited as biological evidence until the dataset is really integrated.
        """
        status = {
            DatasetIntegrationStatus.BLOCKED_NOT_INTEGRATED: SourceStatus.BLOCKED,
            DatasetIntegrationStatus.IN_MEMORY_FIXTURE: SourceStatus.FIXTURE,
            DatasetIntegrationStatus.INTEGRATED: SourceStatus.ACTIVE,
        }[self.integration_status]
        citation = self.citation or f"declared, not integrated: {self.dataset_id}"
        return SourceRecord(
            source_id=self.dataset_id,
            provenance=ProvenanceKind.SOURCED,
            origin=SourceOrigin.PUBLIC_EXPERIMENT,
            status=status,
            description=self.notes or f"public microscopy dataset {self.dataset_id}",
            citation=citation,
            held_out_split_id=self.held_out_split_id,
        )


BLOCKED_PUBLIC_MICROSCOPY_DATASETS: Mapping[str, MicroscopyDatasetContract] = MappingProxyType(
    {
        contract.dataset_id: contract
        for contract in (
            MicroscopyDatasetContract(
                dataset_id="jump-cp-cpg0016",
                integration_status=DatasetIntegrationStatus.BLOCKED_NOT_INTEGRATED,
                modality="widefield-fluorescence-cell-painting",
                accession="cpg0016-jump",
                notes=(
                    "JUMP Cell Painting chemical/genetic perturbation morphology; declared as a "
                    "future observation-validation source, nothing downloaded or integrated"
                ),
            ),
            MicroscopyDatasetContract(
                dataset_id="allen-cell-hipsc-single-cell-image",
                integration_status=DatasetIntegrationStatus.BLOCKED_NOT_INTEGRATED,
                modality="spinning-disk-confocal-fluorescence",
                accession="allencell-hipsc-single-cell-image-dataset",
                notes=(
                    "Allen Institute hiPSC single-cell images for intracellular organization and "
                    "shape variation; blocked, and note the confocal axial gate this lane does not "
                    "implement"
                ),
            ),
            MicroscopyDatasetContract(
                dataset_id="idr-live-cell-timelapse",
                integration_status=DatasetIntegrationStatus.BLOCKED_NOT_INTEGRATED,
                modality="live-cell-timelapse-fluorescence",
                accession="idr-studies",
                notes=(
                    "Image Data Resource live-cell timelapse studies for dynamics and cell-to-cell "
                    "variation; blocked, not integrated"
                ),
            ),
        )
    }
)


def worked_fixture_dataset() -> tuple[MicroscopyDatasetContract, Mapping[str, object]]:
    """Return one in-memory dataset that actually satisfies the adapter contract, plus its payload.

    This exists to show what a compliant adapter must hand over, and to give the contract test
    something real to check.  It is a FIXTURE, not data: the payload is a deterministic constant, and
    :meth:`MicroscopyDatasetContract.to_source_record` gives it ``SourceStatus.FIXTURE`` so it and
    everything derived from it are excluded from the biological corpus.

    Returns:
        The contract and an in-memory payload mapping channel name to a nested list of rows.
    """
    contract = MicroscopyDatasetContract(
        dataset_id="fixture-two-channel-widefield",
        integration_status=DatasetIntegrationStatus.IN_MEMORY_FIXTURE,
        modality="widefield-fluorescence",
        citation="in-memory contract fixture; not a publication and not data",
        licence="not-applicable-fixture",
        accession="fixture:two-channel-widefield:v1",
        channels=("dna", "actin"),
        channel_marker_map={"dna": "Hoechst 33342", "actin": "phalloidin-568"},
        pixel_size_um=0.325,
        axial_step_um=0.29,
        exposure_s=0.05,
        bit_depth=16,
        psf_calibration="declared Gaussian sigma 0.21 um from a bead measurement protocol",
        background_calibration="declared constant 100 ADU offset with a flat-field reference",
        species="fixture-species",
        cell_line="fixture-line",
        biological_state="fixture-interphase-adherent",
        perturbations=("vehicle-control",),
        held_out_split_id="fixture-holdout-a",
        notes="worked contract example; deterministic constant payload, never biological evidence",
    )
    payload: Mapping[str, object] = MappingProxyType(
        {
            "shape": (4, 4),
            "dtype": "<u2",
            "frames": MappingProxyType(
                {
                    "dna": tuple(
                        tuple(100 + 3 * row + column for column in range(4)) for row in range(4)
                    ),
                    "actin": tuple(
                        tuple(100 + 2 * row * column for column in range(4)) for row in range(4)
                    ),
                }
            ),
        }
    )
    return contract, payload


def register_public_dataset_contracts(registry: SourceRegistry) -> SourceRegistry:
    """Register every declared public dataset plus the worked fixture into ``registry``.

    Every real dataset lands as ``BLOCKED`` and the fixture lands as ``FIXTURE``, so the registry's
    biological-evidence query is empty until a dataset is genuinely integrated.

    Args:
        registry: Registry to populate.

    Returns:
        The same registry, for chaining.
    """
    if not isinstance(registry, SourceRegistry):
        raise TypeError("registry must be a SourceRegistry")
    for contract in BLOCKED_PUBLIC_MICROSCOPY_DATASETS.values():
        registry.register(contract.to_source_record())
    fixture_contract, _ = worked_fixture_dataset()
    registry.register(fixture_contract.to_source_record())
    return registry
