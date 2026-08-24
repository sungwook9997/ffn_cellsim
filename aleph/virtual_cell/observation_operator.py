"""The one-way map from mechanistic state to an observable — master plan section 6.1.

The plan declares two central objects and the adversarial audit (`audit_checklist` item A13) found
mechanically that neither existed: eleven lanes of good parts, and nothing that makes them a circuit.
This module is one half of the connective tissue.  It is deliberately thin — it computes almost
nothing new — because its job is to give every lane a common shape to plug into, not to add a
twelfth part.

Why a Protocol and not a base class
-----------------------------------
An observation operator is a structural contract: *given a mechanistic state, produce an observable
in declared units, or REFUSE*.  Nothing about that requires a shared implementation, and demanding
inheritance would force a future runtime's renderer to subclass this package to be usable.
:class:`ObservationOperator` is therefore a ``runtime_checkable`` :class:`typing.Protocol`, and the
concrete adapters below satisfy it without being related to each other.

Why an operator may refuse
--------------------------
The failure this prevents is the one the project has actually paid for: an observable computed from a
state that did not contain what the operator needed, returning a number that is arithmetically
well-defined and physically meaningless.  Every operator declares an :class:`InputRequirement` — named
fields and their declared unit or kind — and :meth:`ObservationOperator.applicability` compares the
state against it BEFORE any arithmetic.  A refusal is a typed :class:`ApplicabilityVerdict` naming the
missing fields, the unit mismatches, and the expansion the failure routes to.  A refusal never carries
a number, for the same reason ``surrogate.SurrogateResponse.value`` raises: a caller who receives a
number will use it.

What is implemented, and what is NAMED as missing
-------------------------------------------------
Implemented, as adapters over code that already exists and is already tested:

* :class:`FluorescenceOperator` over :func:`~aleph.virtual_cell.optics.render_emitters` and
  :func:`~aleph.virtual_cell.optics.render_widefield_stack`;
* :class:`SummaryProjectionOperator` over
  :class:`~aleph.virtual_cell.telemetry.ScalarProjectionSpec` — the legacy scalar-gate path, which
  is versioned and content-hashed rather than being a free function.

NOT implemented, and named rather than stubbed: brightfield, traction, AFM force-indentation and
sequencing.  Each is a member of :data:`UNIMPLEMENTED_OPERATORS` with a typed
:class:`ObservationCapabilityVerdict` saying what is physically missing and which expansion it opens
— the same discipline :data:`~aleph.virtual_cell.optics.UNIMPLEMENTED_AXIAL_MODELS` already applies
to the confocal pinhole.  :class:`UnimplementedOperator` satisfies the protocol so the gap is
enumerable and can be carried in a posterior's operator table, and its :meth:`observe` raises
unconditionally so it can never contribute a number.  Naming what is missing is the deliverable; a
stub that returns a plausible array is the failure.

That refusal is enforced at the REGISTRATION boundary and not only inside the class.  An external
review (2026-07-29, finding 3) subclassed :class:`SummaryProjectionOperator`, overrode ``modality`` to
``TRACTION`` and ``operator_id`` to ``unimplemented::traction``, satisfied the structural protocol,
was carried in a posterior's operator table, and returned the scalar ``5.0`` for a modality the
capability table calls ``BLOCKED``.  :func:`validate_operator_registration` is the answer: any table
that admits an operator must call it, and it refuses an operator claiming a blocked modality or the
reserved ``unimplemented::`` id namespace unless it is the genuine
:class:`UnimplementedOperator` carrying that modality's own capability verdict.

Provenance
----------
Every result carries an :class:`~aleph.virtual_cell.contracts.ArtifactEnvelope`, and every operator
can mint a :class:`~aleph.virtual_cell.source_registry.SourceRecord` whose origin is
``SYNTHETIC_RENDER`` or ``SIMULATION`` — both outside the empirical-origin whitelist.  The registry's
transitive exclusion then does the rest: anything derived from a synthetic observation is excluded
from the biological-truth corpus at every hop, which is master-plan verdict A6 expressed as a graph
rather than as a warning.

Sanity Gate (recorded before first execution):

* **Dimensional** — a fluorescence result is EXPECTED PHOTONS per pixel (the unit
  :mod:`aleph.virtual_cell.optics` produces and
  :func:`~aleph.virtual_cell.synthetic_microscopy.apply_camera_model` consumes); a summary
  projection carries the unit declared on its :class:`ScalarProjectionSpec`.  Units are strings on
  both the requirement and the result, and a state whose declared field unit differs from the
  requirement is REFUSED rather than silently reinterpreted.
* **Boundary** — an empty emitter sequence renders the background frame and does not raise; a missing
  field, a wrong-kind field element, a wrong array rank and a unit mismatch all refuse; an
  unimplemented modality refuses before touching the state at all.
* **Conservation / invariant** — this module performs no photon arithmetic of its own.  The operator
  output is required to be BIT-IDENTICAL to calling the underlying renderer directly with the same
  arguments, which is what makes optics' photon-budget gate still apply here; that identity is
  tested rather than asserted.
* **Measurement protocol** — an observation's identity is
  ``sha256(operator_id, modality, renderer/projection identity, state id, state provenance)``, and the
  state's provenance digest covers the array CONTENT (canonical dtype, shape and the SHA-256 of the
  C-order bytes, in :class:`~aleph.virtual_cell.sidecar.ArraySidecarDescriptor`'s vocabulary), so two
  results are comparable only if the operator, its configuration AND the observed numbers all match.
  The projection path additionally binds ``ScalarProjectionSpec.definition_sha256``, so relaxing a
  projection's tolerance mints a different projection instead of masquerading as the historical one.
* **Numerical** — no new floating-point path exists here; the closed forms, quadrature and truncation
  policy are optics', and the reduction whitelist is telemetry's.

Nothing in this module imports Warp, runs on a GPU, or grants any evidence rung.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from types import MappingProxyType
from typing import Any, Protocol, runtime_checkable

import numpy as np

from aleph.virtual_cell.contracts import (
    ArtifactEnvelope,
    EvidenceSource,
    RepresentationDescriptor,
    RepresentationKind,
    UncertaintySource,
)
from aleph.virtual_cell.optics import (
    LineEmitter,
    OpticsConfig,
    PointEmitter,
    TriangleEmitter,
    render_emitters,
    render_widefield_stack,
)
from aleph.virtual_cell.source_registry import (
    ProvenanceKind,
    SourceOrigin,
    SourceRecord,
    SourceStatus,
)
from aleph.virtual_cell.telemetry import ScalarProjectionSpec

__all__ = [
    "BLOCKED_MODALITIES",
    "RESERVED_OPERATOR_ID_PREFIX",
    "UNIMPLEMENTED_OPERATORS",
    "ApplicabilityVerdict",
    "CapabilityStatus",
    "FluorescenceOperator",
    "InputRequirement",
    "MechanisticState",
    "ObservationCapabilityVerdict",
    "ObservationModality",
    "ObservationOperator",
    "ObservationRefused",
    "ObservationResult",
    "OperatorRegistrationRefused",
    "SummaryProjectionOperator",
    "UnimplementedOperator",
    "unimplemented_operator",
    "validate_operator_registration",
]

_EMITTER_TYPES = (PointEmitter, LineEmitter, TriangleEmitter)

#: Declared unit/kind string for a field holding :data:`~aleph.virtual_cell.optics.Emitter` values.
#: Emitters are not scalars, so what a requirement pins is their KIND, not a physical unit.
EMITTER_KIND = "optics-emitter-sequence"

#: Unit of every rendered fluorescence frame.  Optics normalises the PSF before the frame crop, so
#: this is a photon COUNT per pixel and never a density.
PHOTON_UNIT = "expected-photons/pixel"


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _sha256_json(payload: Any) -> str:
    """Return the SHA-256 of a canonical JSON encoding of ``payload``."""
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode()).hexdigest()


#: Memory order the content address is taken in.  Fixed, so the same numbers in a Fortran-ordered and
#: a C-ordered buffer hash the same; the vocabulary is
#: :class:`~aleph.virtual_cell.sidecar.ArraySidecarDescriptor`'s, reused rather than reinvented so an
#: array has ONE content address in this package.
_ARRAY_ORDER = "C"

#: A default ``object.__repr__`` carries the instance's memory address, which is neither content nor
#: stable across interpreters.  A field whose repr looks like that cannot be content-addressed.
_ADDRESS_REPR = re.compile(r"^<.* object at 0x[0-9a-f]+>$")


def _stable_repr(value: object, *, what: str) -> str:
    text = repr(value)
    if _ADDRESS_REPR.match(text):
        raise ValueError(
            f"state field {what!r} holds a {type(value).__qualname__} whose repr is its memory "
            "address, so it has no content address; give it a value-based __repr__ (a frozen "
            "dataclass has one), or present it as an array"
        )
    return text


def _content_token(value: object, *, what: str) -> Any:
    """Return a JSON-able CONTENT token for one state field value.

    An array contributes its canonical dtype, its shape and the SHA-256 of its bytes — never its
    identity — so two states holding different numbers can never present the same digest.

    Raises:
        ValueError: If the value cannot be content-addressed (object-dtype array, or an object whose
            repr is its memory address).  Failing closed is the point: a value that cannot be hashed
            must not be silently omitted from the hash.
    """
    if isinstance(value, (np.ndarray, np.generic)):
        array = np.asarray(value)
        if array.dtype.hasobject:
            raise ValueError(
                f"state field {what!r} is an object-dtype array; its buffer holds pointers, so "
                "hashing it would address the pointers rather than the numbers"
            )
        return {
            "kind": "ndarray",
            "dtype": array.dtype.str,
            "shape": list(array.shape),
            "array_order": _ARRAY_ORDER,
            "sha256": hashlib.sha256(
                np.ascontiguousarray(array).tobytes(order=_ARRAY_ORDER)
            ).hexdigest(),
        }
    if value is None or isinstance(value, (bool, int, float, str, bytes)):
        return {"kind": "scalar", "repr": repr(value)}
    if isinstance(value, Mapping):
        return {
            "kind": "mapping",
            "items": [
                [str(key), _content_token(item, what=f"{what}[{key!r}]")]
                for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
            ],
        }
    if isinstance(value, Sequence):
        return {
            "kind": "sequence",
            "items": [
                _content_token(item, what=f"{what}[{index}]") for index, item in enumerate(value)
            ],
        }
    return {
        "kind": "opaque",
        "type": type(value).__qualname__,
        "repr": _stable_repr(value, what=what),
    }


def _own_content(value: object) -> Any:
    """Return a value the state OWNS: arrays copied and frozen, sequences and mappings immutable.

    Without this the state keeps the caller's object, so a caller can mutate the array after
    construction and make two different observations share one identity (external review 2026-07-29,
    finding 2 — the same counterexample, reached by a second route).
    """
    if isinstance(value, np.ndarray):
        array = np.array(value, order=_ARRAY_ORDER, copy=True)
        array.setflags(write=False)
        return array
    if isinstance(value, Mapping):
        return MappingProxyType({key: _own_content(item) for key, item in value.items()})
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        return value
    return tuple(_own_content(item) for item in value)


def _unique_strings(
    values: Sequence[str], *, what: str, allow_empty: bool = False
) -> tuple[str, ...]:
    if isinstance(values, str) or not isinstance(values, Sequence):
        raise TypeError(f"{what} must be a sequence of strings")
    normalized = tuple(_nonempty(value, what=f"{what} entry") for value in values)
    if not allow_empty and not normalized:
        raise ValueError(f"{what} must not be empty")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{what} must not contain duplicates")
    return normalized


class ObservationModality(StrEnum):
    """The observation modalities master plan section 6.1 enumerates.

    A modality is a statement about the physical measurement, not about the code: two operators of the
    same modality may disagree completely, and an operator whose modality is listed in
    :data:`UNIMPLEMENTED_OPERATORS` exists only to refuse.
    """

    FLUORESCENCE = "fluorescence"
    BRIGHTFIELD = "brightfield"
    TRACTION = "traction"
    FORCE_INDENTATION = "force-indentation"
    SEQUENCING = "sequencing"
    SUMMARY_PROJECTION = "summary-projection"


class CapabilityStatus(StrEnum):
    """Why a declared capability is unavailable.

    The three values are exactly the vocabulary
    :class:`~aleph.virtual_cell.optics.OpticsCapabilityVerdict` already validates against, kept
    identical on purpose so a capability gap reads the same wherever it is raised; the test suite pins
    the two vocabularies together rather than trusting this comment.
    """

    BLOCKED = "blocked"
    ORACLE_ONLY = "oracle-only"
    PI_DECISION = "pi-decision"


@dataclass(frozen=True, slots=True)
class ObservationCapabilityVerdict:
    """A named, unimplemented observation capability and the expansion it opens.

    ``why_not_approximated`` is required and is the whole point of the type: an observation operator
    that could be *approximated* by an existing one would not be a gap, and writing down why the
    approximation is wrong is what stops the approximation being quietly adopted later.
    """

    modality: ObservationModality
    capability: str
    status: CapabilityStatus
    why_not_approximated: str
    routes_to: str

    def __post_init__(self) -> None:
        if not isinstance(self.modality, ObservationModality):
            raise TypeError("modality must be an ObservationModality")
        if not isinstance(self.status, CapabilityStatus):
            raise TypeError("status must be a CapabilityStatus")
        for name in ("capability", "why_not_approximated", "routes_to"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of this capability gap."""
        return {
            "modality": self.modality.value,
            "capability": self.capability,
            "status": self.status.value,
            "why_not_approximated": self.why_not_approximated,
            "routes_to": self.routes_to,
        }


UNIMPLEMENTED_OPERATORS: Mapping[ObservationModality, ObservationCapabilityVerdict] = (
    MappingProxyType(
        {
            ObservationModality.BRIGHTFIELD: ObservationCapabilityVerdict(
                modality=ObservationModality.BRIGHTFIELD,
                capability="label-free transmitted-light contrast",
                status=CapabilityStatus.BLOCKED,
                why_not_approximated=(
                    "brightfield contrast comes from the refractive-index/phase structure of the "
                    "sample modulating a transmitted field; the emitter model carries photon RATES "
                    "and no index map, so inverting a fluorescence render produces a picture with "
                    "the right sign and the wrong physics"
                ),
                routes_to="observation-model extension: refractive-index field plus transmitted-field "
                "propagation (partially coherent transfer function)",
            ),
            ObservationModality.TRACTION: ObservationCapabilityVerdict(
                modality=ObservationModality.TRACTION,
                capability="substrate traction map",
                status=CapabilityStatus.BLOCKED,
                why_not_approximated=(
                    "a traction map is the substrate's surface stress field, obtained from a "
                    "substrate constitutive law and the adjoint of its Green's function; the census "
                    "carries per-clutch connector forces, and summing them onto a grid is a force "
                    "histogram, not a traction field, because it omits the substrate's own mechanics"
                ),
                routes_to="connector-law extension: substrate Green's-function operator with its "
                "declared adjoint (Newton's 3rd law), then a traction reconstruction operator",
            ),
            ObservationModality.FORCE_INDENTATION: ObservationCapabilityVerdict(
                modality=ObservationModality.FORCE_INDENTATION,
                capability="AFM force-indentation protocol",
                status=CapabilityStatus.BLOCKED,
                why_not_approximated=(
                    "the observable is produced by a moving contact boundary driven on the "
                    "instrument's own ramp schedule; a Hertz or Bell-Evans closed form is the "
                    "acceptance ORACLE for such a curve and using it AS the operator would reinstate "
                    "exactly the lumped-model runtime the architecture principle inverts"
                ),
                routes_to="reference-dynamics extension: indenter contact boundary condition plus a "
                "protocol-resolved ramp, with the closed form retained only as the oracle",
            ),
            ObservationModality.SEQUENCING: ObservationCapabilityVerdict(
                modality=ObservationModality.SEQUENCING,
                capability="sequencing composition/state likelihood",
                status=CapabilityStatus.BLOCKED,
                why_not_approximated=(
                    "the likelihood needs per-species molecular counts and a capture/library/depth "
                    "model; the mechanical state carries densities and mechanical inventories, so a "
                    "count drawn from a density would be a convenience number wearing a likelihood"
                ),
                routes_to="knowledge-layer extension: composition state variable plus a "
                "capture/library sampling likelihood, kept separate from the mechanical posterior",
            ),
        }
    )
)


#: Modalities the capability table declares unavailable.  Derived FROM the table rather than typed out
#: again, so a modality cannot become blocked in one place and open in another.
BLOCKED_MODALITIES: frozenset[ObservationModality] = frozenset(UNIMPLEMENTED_OPERATORS)

#: Reserved operator-id namespace.  Only the declared :class:`UnimplementedOperator` refusals may use
#: it; :func:`validate_operator_registration` refuses anything else that does.
RESERVED_OPERATOR_ID_PREFIX = "unimplemented::"


@dataclass(frozen=True, slots=True)
class InputRequirement:
    """What an operator needs to be able to observe at all.

    Args:
        fields: Names the state must carry.  An EMPTY requirement is legal and means the operator
            needs nothing — the only operator for which that is true is
            :class:`UnimplementedOperator`, whose gap is in the operator rather than in the state.
        field_units: Declared unit — or, for a non-scalar field, KIND — per required field.  Every
            required field must appear here: an unstated unit is how a length gets rendered as a rate.
        description: One line saying what the operator does with the fields.
    """

    fields: tuple[str, ...]
    field_units: Mapping[str, str]
    description: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "fields", _unique_strings(self.fields, what="required fields", allow_empty=True)
        )
        if not isinstance(self.field_units, Mapping):
            raise TypeError("field_units must be a mapping")
        units = {
            _nonempty(key, what="field_units key"): _nonempty(value, what=f"field_units[{key!r}]")
            for key, value in self.field_units.items()
        }
        missing = set(self.fields) - set(units)
        if missing:
            raise ValueError(
                "every required field needs a declared unit or kind; missing: "
                + ", ".join(sorted(missing))
            )
        extra = set(units) - set(self.fields)
        if extra:
            raise ValueError(
                "field_units declares units for fields that are not required: "
                + ", ".join(sorted(extra))
            )
        object.__setattr__(self, "field_units", MappingProxyType(units))
        object.__setattr__(self, "description", _nonempty(self.description, what="description"))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of this requirement."""
        return {
            "fields": list(self.fields),
            "field_units": dict(sorted(self.field_units.items())),
            "description": self.description,
        }


@dataclass(frozen=True, slots=True)
class MechanisticState:
    """The state side of the one-way map: named fields, their declared units, and provenance.

    This is deliberately NOT a physics object.  It is the narrow interface an operator sees, so a
    Warp-resident engine state, an analytic fixture and a reduced representation can all present
    themselves to the same operator without this module knowing which it is holding.  ``source_id``
    is the provenance-graph node the observation will be derived FROM, which is what makes the
    transitive synthetic exclusion work without the operator knowing anything about the corpus.

    The state OWNS its content.  Arrays are copied and frozen at construction and sequences and
    mappings are made immutable, and :attr:`content_sha256` is computed once, there.  Both halves come
    from an external review (2026-07-29, finding 2): with ``state_id``, ``manifest_hash``,
    ``source_id`` and the field names held equal, ``tension=[1, 2, 3]`` and ``tension=[1, 2, 30]``
    observed to ``3.0`` and ``30.0`` and presented IDENTICAL ``provenance_digest``, envelope
    ``artifact_id`` and ``observation_provenance_sha256``; and because the mapping was copied only
    shallowly, mutating the caller's ndarray after construction produced the same collision a second
    way.  Two different observations sharing one artifact identity is not a documentation problem, so
    the schema now enforces what the docstring used to ask of the caller.
    """

    state_id: str
    manifest_hash: str
    fields: Mapping[str, Any]
    field_units: Mapping[str, str]
    source_id: str
    content_sha256: str = field(init=False, default="", repr=False, compare=False)

    def __post_init__(self) -> None:
        for name in ("state_id", "manifest_hash", "source_id"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        if len(self.manifest_hash) != 64 or any(
            character not in "0123456789abcdef" for character in self.manifest_hash
        ):
            raise ValueError("manifest_hash must be a lowercase SHA-256 digest")
        if not isinstance(self.fields, Mapping):
            raise TypeError("fields must be a mapping")
        values = {
            _nonempty(key, what="field name"): _own_content(value)
            for key, value in self.fields.items()
        }
        object.__setattr__(self, "fields", MappingProxyType(values))
        if not isinstance(self.field_units, Mapping):
            raise TypeError("field_units must be a mapping")
        units = {
            _nonempty(key, what="field_units key"): _nonempty(value, what=f"field_units[{key!r}]")
            for key, value in self.field_units.items()
        }
        undeclared = set(values) - set(units)
        if undeclared:
            raise ValueError(
                "every state field needs a declared unit or kind; missing: "
                + ", ".join(sorted(undeclared))
            )
        object.__setattr__(self, "field_units", MappingProxyType(units))
        object.__setattr__(
            self,
            "content_sha256",
            _sha256_json(
                {name: _content_token(values[name], what=name) for name in sorted(values)}
            ),
        )

    @property
    def provenance_digest(self) -> str:
        """SHA-256 over the state's identity, provenance AND the content of its fields.

        The content half is :attr:`content_sha256`, computed once at construction: a
        :class:`MechanisticState` is assembled on the host BETWEEN accepted steps, so its fields are
        already host-resident by the time this type holds them and the digest adds no readback that
        building the state did not already perform.  It is never taken inside the physical-time loop.
        """
        return _sha256_json(
            {
                "state_id": self.state_id,
                "manifest_hash": self.manifest_hash,
                "source_id": self.source_id,
                "fields": sorted(self.fields),
                "field_units": dict(sorted(self.field_units.items())),
                "content_sha256": self.content_sha256,
            }
        )


@dataclass(frozen=True, slots=True)
class ApplicabilityVerdict:
    """Whether an operator can observe a state, and — when it cannot — exactly what is missing.

    A refusal names ``routes_to`` because the master plan's failure routing is expansionary: an
    operator that cannot observe a state has found either a missing state field or a missing
    observation model, and both open work rather than closing scope.
    """

    operator_id: str
    state_id: str
    applicable: bool
    reason: str
    missing_fields: tuple[str, ...] = ()
    unit_mismatches: tuple[str, ...] = ()
    routes_to: str | None = None

    def __post_init__(self) -> None:
        for name in ("operator_id", "state_id", "reason"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        if not isinstance(self.applicable, bool):
            raise TypeError("applicable must be bool")
        object.__setattr__(
            self,
            "missing_fields",
            _unique_strings(self.missing_fields, what="missing_fields", allow_empty=True),
        )
        object.__setattr__(
            self,
            "unit_mismatches",
            _unique_strings(self.unit_mismatches, what="unit_mismatches", allow_empty=True),
        )
        if self.routes_to is not None:
            object.__setattr__(self, "routes_to", _nonempty(self.routes_to, what="routes_to"))
        if self.applicable:
            if self.missing_fields or self.unit_mismatches:
                raise ValueError("an applicable verdict cannot also name a defect")
            if self.routes_to is not None:
                raise ValueError("an applicable verdict routes to no expansion")
        elif self.routes_to is None:
            raise ValueError(
                "a refusal must name the expansion it opens; a bare refusal closes scope, which is "
                "the failure routing this project forbids"
            )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of this verdict."""
        return {
            "operator_id": self.operator_id,
            "state_id": self.state_id,
            "applicable": self.applicable,
            "reason": self.reason,
            "missing_fields": list(self.missing_fields),
            "unit_mismatches": list(self.unit_mismatches),
            "routes_to": self.routes_to,
        }


class ObservationRefused(RuntimeError):
    """Raised when an operator is asked to observe a state it cannot observe.

    It is an error rather than a sentinel value for the reason the surrogate harness gives: a caller
    who receives a number will use it, so a refusal must not be able to look like a result.
    """

    def __init__(self, verdict: ApplicabilityVerdict) -> None:
        if not isinstance(verdict, ApplicabilityVerdict):
            raise TypeError("verdict must be an ApplicabilityVerdict")
        self.verdict = verdict
        detail = verdict.reason
        if verdict.missing_fields:
            detail += "; missing fields: " + ", ".join(verdict.missing_fields)
        if verdict.unit_mismatches:
            detail += "; unit mismatches: " + ", ".join(verdict.unit_mismatches)
        super().__init__(
            f"{verdict.operator_id!r} refused state {verdict.state_id!r}: {detail}; "
            f"routes to {verdict.routes_to}"
        )


@dataclass(frozen=True, slots=True)
class ObservationResult:
    """One observation, its units, and the envelope that keeps it out of the biological corpus."""

    operator_id: str
    values: np.ndarray
    units: str
    envelope: ArtifactEnvelope
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "operator_id", _nonempty(self.operator_id, what="operator_id"))
        array = np.asarray(self.values, dtype=np.float64)
        if not np.all(np.isfinite(array)):
            raise ValueError(
                "an observation must be finite; a NaN observable is not an observation"
            )
        array = array.copy()
        array.setflags(write=False)
        object.__setattr__(self, "values", array)
        object.__setattr__(self, "units", _nonempty(self.units, what="units"))
        if not isinstance(self.envelope, ArtifactEnvelope):
            raise TypeError("envelope must be an ArtifactEnvelope")
        if not isinstance(self.metadata, Mapping):
            raise TypeError("metadata must be a mapping")
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    @property
    def modality(self) -> ObservationModality:
        """Modality recorded in the metadata, as an enum member."""
        return ObservationModality(self.metadata["modality"])

    @property
    def scalar(self) -> float:
        """The single value of a scalar observation.

        Raises:
            ValueError: If the observation is not a single value.  A mean over an image is a
                different observable from the image, and this refuses to conflate them.
        """
        if self.values.size != 1:
            raise ValueError(
                f"observation {self.operator_id!r} has {self.values.size} values; reducing it here "
                "would be an undeclared projection"
            )
        return float(self.values.reshape(-1)[0])


@runtime_checkable
class ObservationOperator(Protocol):
    """The one-way map from mechanistic state to an observable.

    An implementation must declare (a) what state it requires, (b) what units it emits, and (c) its
    own applicability, and must REFUSE rather than return a number for a state it cannot observe.  The
    map is one-way by construction: nothing in this protocol inverts an observation back into state.
    Inversion is a separate, harder problem and lives in
    :mod:`aleph.virtual_cell.sandbox_inverse`, where its non-identifiability is measured.
    """

    @property
    def operator_id(self) -> str:
        """Stable identifier, recorded on every artifact this operator produces."""
        ...

    @property
    def modality(self) -> ObservationModality:
        """Which physical measurement this operator claims to model."""
        ...

    @property
    def input_requirement(self) -> InputRequirement:
        """The state fields and declared units/kinds this operator needs."""
        ...

    @property
    def output_units(self) -> str:
        """Unit of the values this operator emits."""
        ...

    def applicability(self, state: MechanisticState) -> ApplicabilityVerdict:
        """Judge whether this operator can observe ``state``, without observing it."""
        ...

    def observe(self, state: MechanisticState) -> ObservationResult:
        """Observe ``state``, or raise :class:`ObservationRefused`."""
        ...

    def provenance_record(
        self, *, source_id: str, parent_source_ids: Sequence[str], description: str
    ) -> SourceRecord:
        """Mint the provenance node an observation from this operator should be registered as."""
        ...


class _RequirementChecked:
    """Applicability derived from the declared requirement.  Not part of the public protocol.

    Shared so that "missing field" and "unit mismatch" mean the same thing for every operator; an
    operator with extra conditions checks them after calling this.
    """

    __slots__ = ()

    def _base_verdict(self, state: MechanisticState, *, routes_to: str) -> ApplicabilityVerdict:
        if not isinstance(state, MechanisticState):
            raise TypeError("state must be a MechanisticState")
        requirement = self.input_requirement  # type: ignore[attr-defined]
        operator_id = self.operator_id  # type: ignore[attr-defined]
        missing = tuple(name for name in requirement.fields if name not in state.fields)
        mismatches = tuple(
            f"{name}: state declares {state.field_units[name]!r}, operator requires "
            f"{requirement.field_units[name]!r}"
            for name in requirement.fields
            if name in state.fields and state.field_units.get(name) != requirement.field_units[name]
        )
        if missing or mismatches:
            return ApplicabilityVerdict(
                operator_id=operator_id,
                state_id=state.state_id,
                applicable=False,
                reason="the state does not carry what this operator declared it needs",
                missing_fields=missing,
                unit_mismatches=mismatches,
                routes_to=routes_to,
            )
        return ApplicabilityVerdict(
            operator_id=operator_id,
            state_id=state.state_id,
            applicable=True,
            reason=f"state carries {', '.join(requirement.fields)} in the declared units",
        )


@dataclass(frozen=True, slots=True)
class FluorescenceOperator(_RequirementChecked):
    """Fluorescence renderer as an observation operator, over the existing optics module.

    ``focal_planes_um`` selects between the two shapes optics already produces: ``None`` gives the
    single widefield z-projection at the configuration's own focal plane, and a non-empty tuple gives
    the through-focus stack.  Nothing is re-derived here; the operator's contribution is the declared
    requirement, the refusal, and the provenance binding.

    Args:
        operator_id: Stable identifier.
        channel: Which fluorescence channel this operator renders.  Two channels of the same sample
            are two operators, because they see different labels and therefore different states.
        optics: The optics configuration.  Constructing one with an unimplemented axial model already
            raises inside optics, so a confocal or light-sheet operator cannot be built at all.
        emitter_field: Name of the state field holding the emitters.
    """

    operator_id: str
    channel: str
    optics: OpticsConfig
    emitter_field: str
    focal_planes_um: tuple[float, ...] | None = None

    def __post_init__(self) -> None:
        for name in ("operator_id", "channel", "emitter_field"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        if not isinstance(self.optics, OpticsConfig):
            raise TypeError("optics must be an OpticsConfig")
        if self.focal_planes_um is not None:
            planes = tuple(float(value) for value in self.focal_planes_um)
            if not planes:
                raise ValueError(
                    "focal_planes_um must be None (single projection) or a non-empty tuple; an "
                    "empty stack is not a z stack"
                )
            if any(not np.isfinite(value) for value in planes):
                raise ValueError("focal_planes_um entries must be finite")
            object.__setattr__(self, "focal_planes_um", planes)

    @property
    def modality(self) -> ObservationModality:
        """Always :attr:`ObservationModality.FLUORESCENCE`; a property, so there is no field to set."""
        return ObservationModality.FLUORESCENCE

    @property
    def output_units(self) -> str:
        """Expected photons per pixel, at unit detection gain."""
        return PHOTON_UNIT

    @property
    def input_requirement(self) -> InputRequirement:
        """The emitter sequence, declared as a kind rather than a unit."""
        return InputRequirement(
            fields=(self.emitter_field,),
            field_units={self.emitter_field: EMITTER_KIND},
            description=(
                f"render channel {self.channel!r} from point/line/triangle emitters through "
                f"{self.optics.renderer_id}"
            ),
        )

    @property
    def axes(self) -> tuple[str, ...]:
        """Axis names of the rendered array, in order."""
        if self.focal_planes_um is None:
            return ("y_px", "x_px")
        return ("z_plane", "y_px", "x_px")

    def applicability(self, state: MechanisticState) -> ApplicabilityVerdict:
        """Refuse a state with no emitter field, a wrong-kind field, or a non-emitter element."""
        routes_to = (
            "state-model extension: the observed component must expose a geometry-to-emitter map "
            "(see archetypes/grammar for which component owns the geometry)"
        )
        verdict = self._base_verdict(state, routes_to=routes_to)
        if not verdict.applicable:
            return verdict
        emitters = state.fields[self.emitter_field]
        if isinstance(emitters, _EMITTER_TYPES) or not isinstance(emitters, Sequence):
            return ApplicabilityVerdict(
                operator_id=self.operator_id,
                state_id=state.state_id,
                applicable=False,
                reason=(
                    f"field {self.emitter_field!r} must be a SEQUENCE of emitters; a single emitter "
                    f"or a scalar is not a labelled structure"
                ),
                routes_to=routes_to,
            )
        bad = tuple(
            f"{self.emitter_field}[{index}]: {type(entry).__name__}"
            for index, entry in enumerate(emitters)
            if not isinstance(entry, _EMITTER_TYPES)
        )
        if bad:
            return ApplicabilityVerdict(
                operator_id=self.operator_id,
                state_id=state.state_id,
                applicable=False,
                reason="the emitter sequence holds values this renderer has no forward model for",
                unit_mismatches=bad,
                routes_to=routes_to,
            )
        return verdict

    def provenance_digest(self, state: MechanisticState) -> str:
        """SHA-256 binding this operator's configuration to the state it observed."""
        return _sha256_json(
            {
                "operator_id": self.operator_id,
                "modality": self.modality.value,
                "channel": self.channel,
                "renderer_id": self.optics.renderer_id,
                "optics_config_sha256": self.optics.config_hash,
                "focal_planes_um": list(self.focal_planes_um or ()),
                "state_provenance_sha256": state.provenance_digest,
            }
        )

    def observe(self, state: MechanisticState) -> ObservationResult:
        """Render the state's emitters, or refuse.

        Args:
            state: The state to observe.

        Returns:
            An :class:`ObservationResult` whose values are expected photons per pixel and whose
            envelope is SYNTHETIC by construction.

        Raises:
            ObservationRefused: If :meth:`applicability` refuses the state.
        """
        verdict = self.applicability(state)
        if not verdict.applicable:
            raise ObservationRefused(verdict)
        emitters = list(state.fields[self.emitter_field])
        if self.focal_planes_um is None:
            values = render_emitters(emitters, self.optics)
        else:
            values = render_widefield_stack(emitters, self.optics, self.focal_planes_um)
        provenance = self.provenance_digest(state)
        envelope = ArtifactEnvelope(
            artifact_id=f"{self.operator_id}@{provenance[:16]}",
            cell_state_manifest_hash=state.manifest_hash,
            evidence_source=EvidenceSource.SYNTHETIC,
            representation=RepresentationDescriptor(
                kind=RepresentationKind.SYNTHETIC_IMAGE,
                uncertainty_sources=(UncertaintySource.MEASUREMENT,),
                axes=self.axes,
            ),
            source_artifact_ids=(state.state_id,),
            observation_provenance_sha256=provenance,
        )
        return ObservationResult(
            operator_id=self.operator_id,
            values=values,
            units=self.output_units,
            envelope=envelope,
            metadata={
                "modality": self.modality.value,
                "channel": self.channel,
                "renderer_id": self.optics.renderer_id,
                "n_emitters": len(emitters),
                "n_planes": 1 if self.focal_planes_um is None else len(self.focal_planes_um),
            },
        )

    def provenance_record(
        self, *, source_id: str, parent_source_ids: Sequence[str], description: str
    ) -> SourceRecord:
        """Mint a ``DERIVED`` / ``SYNTHETIC_RENDER`` node for a rendered observation.

        ``SYNTHETIC_RENDER`` is outside the registry's empirical-origin whitelist, so this node and
        everything downstream of it is excluded from the biological corpus transitively.  That is the
        structural form of master-plan verdict A6.
        """
        return SourceRecord(
            source_id=_nonempty(source_id, what="source_id"),
            provenance=ProvenanceKind.DERIVED,
            origin=SourceOrigin.SYNTHETIC_RENDER,
            status=SourceStatus.ACTIVE,
            description=_nonempty(description, what="description"),
            parents=tuple(parent_source_ids),
        )


@dataclass(frozen=True, slots=True)
class SummaryProjectionOperator(_RequirementChecked):
    """The legacy scalar-gate path, as an observation operator over a versioned projection.

    Master plan section 6.1 lists ``summary projection -> legacy scalar observable`` beside the image
    and force modalities on purpose: the scalar a gate reads IS an observation, and treating it as
    something more direct is how a headline number came to be produced by a free function nobody could
    re-run.  :class:`~aleph.virtual_cell.telemetry.ScalarProjectionSpec` already fixed that by
    storing the projection as content-hashed DATA; this operator is the adapter that puts it on the
    same interface as the renderer.

    Args:
        operator_id: Stable identifier.
        spec: The stored projection, including the tolerance that is part of its hash.
        array_field: Name of the state field holding the array.
        axis_names: Ordered axis names of that array.
        array_unit: Declared unit of the array's entries.
    """

    operator_id: str
    spec: ScalarProjectionSpec
    array_field: str
    axis_names: tuple[str, ...]
    array_unit: str

    def __post_init__(self) -> None:
        for name in ("operator_id", "array_field", "array_unit"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        if not isinstance(self.spec, ScalarProjectionSpec):
            raise TypeError("spec must be a ScalarProjectionSpec")
        object.__setattr__(self, "axis_names", _unique_strings(self.axis_names, what="axis_names"))
        unknown = set(self.spec.selection) - set(self.axis_names)
        if unknown:
            raise ValueError(
                "the projection selects axes this operator does not declare: "
                + ", ".join(sorted(unknown))
            )

    @property
    def modality(self) -> ObservationModality:
        """Always :attr:`ObservationModality.SUMMARY_PROJECTION`; a property, with no field behind it."""
        return ObservationModality.SUMMARY_PROJECTION

    @property
    def output_units(self) -> str:
        """The unit declared on the stored projection."""
        return self.spec.unit

    @property
    def input_requirement(self) -> InputRequirement:
        """The named array, in its declared unit."""
        return InputRequirement(
            fields=(self.array_field,),
            field_units={self.array_field: self.array_unit},
            description=(
                f"apply projection {self.spec.name}@{self.spec.version} "
                f"({self.spec.reduction.value}) over axes {', '.join(self.axis_names)}"
            ),
        )

    def applicability(self, state: MechanisticState) -> ApplicabilityVerdict:
        """Refuse a missing array, a unit mismatch, or an array whose rank is not the declared one."""
        routes_to = (
            "telemetry extension: the projection's axes must exist in the archive; see "
            "telemetry.build_telemetry_archive for the axis contract"
        )
        verdict = self._base_verdict(state, routes_to=routes_to)
        if not verdict.applicable:
            return verdict
        array = np.asarray(state.fields[self.array_field], dtype=np.float64)
        if array.ndim != len(self.axis_names):
            return ApplicabilityVerdict(
                operator_id=self.operator_id,
                state_id=state.state_id,
                applicable=False,
                reason=(
                    f"array rank {array.ndim} does not match the {len(self.axis_names)} declared "
                    f"axes ({', '.join(self.axis_names)}); a projection applied to the wrong rank "
                    "reduces the wrong axis and still returns a number"
                ),
                routes_to=routes_to,
            )
        return verdict

    def provenance_digest(self, state: MechanisticState) -> str:
        """SHA-256 binding this projection definition to the state it was applied to."""
        return _sha256_json(
            {
                "operator_id": self.operator_id,
                "modality": self.modality.value,
                "projection_sha256": self.spec.definition_sha256,
                "axis_names": list(self.axis_names),
                "array_field": self.array_field,
                "state_provenance_sha256": state.provenance_digest,
            }
        )

    def observe(self, state: MechanisticState) -> ObservationResult:
        """Project the state's array to one scalar, or refuse.

        Args:
            state: The state to observe.

        Returns:
            An :class:`ObservationResult` holding a single value in the projection's declared unit.

        Raises:
            ObservationRefused: If :meth:`applicability` refuses the state.
        """
        verdict = self.applicability(state)
        if not verdict.applicable:
            raise ObservationRefused(verdict)
        value = self.spec.evaluate(state.fields[self.array_field], self.axis_names)
        provenance = self.provenance_digest(state)
        envelope = ArtifactEnvelope(
            artifact_id=f"{self.operator_id}@{provenance[:16]}",
            cell_state_manifest_hash=state.manifest_hash,
            evidence_source=EvidenceSource.REDUCED_MODEL,
            representation=RepresentationDescriptor(
                kind=RepresentationKind.SCALAR_PROJECTION,
                uncertainty_sources=(UncertaintySource.REPRESENTATION,),
                axes=self.axis_names,
            ),
            source_artifact_ids=(state.state_id,),
            projection_hash=self.spec.definition_sha256,
            observation_provenance_sha256=provenance,
        )
        return ObservationResult(
            operator_id=self.operator_id,
            values=np.asarray(value, dtype=np.float64),
            units=self.output_units,
            envelope=envelope,
            metadata={
                "modality": self.modality.value,
                "projection": f"{self.spec.name}@{self.spec.version}",
                "projection_sha256": self.spec.definition_sha256,
                "reduction": self.spec.reduction.value,
                "tolerance": self.spec.tolerance,
            },
        )

    def provenance_record(
        self, *, source_id: str, parent_source_ids: Sequence[str], description: str
    ) -> SourceRecord:
        """Mint a ``DERIVED`` / ``SIMULATION`` node for a projected scalar.

        ``SIMULATION`` is likewise outside the empirical whitelist: a number a simulation produced is
        not evidence about a cell, however many gates it passed.
        """
        return SourceRecord(
            source_id=_nonempty(source_id, what="source_id"),
            provenance=ProvenanceKind.DERIVED,
            origin=SourceOrigin.SIMULATION,
            status=SourceStatus.ACTIVE,
            description=_nonempty(description, what="description"),
            parents=tuple(parent_source_ids),
        )


@dataclass(frozen=True, slots=True)
class UnimplementedOperator(_RequirementChecked):
    """A NAMED observation operator that does not exist yet, and refuses unconditionally.

    It satisfies :class:`ObservationOperator` so the gap is a first-class, enumerable member of a
    posterior's operator table rather than an absence nobody can query — the same reason
    :class:`~aleph.virtual_cell.archetypes.MissingLawVerdict` is an object and not a missing entry.
    :meth:`observe` raises whatever it is given, so it can never contribute a number to anything.
    """

    operator_id: str
    verdict: ObservationCapabilityVerdict

    def __post_init__(self) -> None:
        object.__setattr__(self, "operator_id", _nonempty(self.operator_id, what="operator_id"))
        if not isinstance(self.verdict, ObservationCapabilityVerdict):
            raise TypeError("verdict must be an ObservationCapabilityVerdict")

    @property
    def modality(self) -> ObservationModality:
        """The modality this operator would model if it existed."""
        return self.verdict.modality

    @property
    def output_units(self) -> str:
        """Deliberately not a physical unit: this operator emits nothing."""
        return "none"

    @property
    def input_requirement(self) -> InputRequirement:
        """An empty requirement — the gap is in the operator, not in the state."""
        return InputRequirement(
            fields=(),
            field_units={},
            description=f"{self.verdict.capability} is {self.verdict.status.value}",
        )

    def _refusal(self, state: MechanisticState) -> ApplicabilityVerdict:
        if not isinstance(state, MechanisticState):
            raise TypeError("state must be a MechanisticState")
        return ApplicabilityVerdict(
            operator_id=self.operator_id,
            state_id=state.state_id,
            applicable=False,
            reason=f"{self.verdict.capability} is {self.verdict.status.value}: "
            f"{self.verdict.why_not_approximated}",
            routes_to=self.verdict.routes_to,
        )

    def applicability(self, state: MechanisticState) -> ApplicabilityVerdict:
        """Always inapplicable, for every state, with the capability gap as the reason."""
        return self._refusal(state)

    def observe(self, state: MechanisticState) -> ObservationResult:
        """Always raise.

        Raises:
            ObservationRefused: Unconditionally.  There is no state for which this returns a value.
        """
        raise ObservationRefused(self._refusal(state))

    def provenance_record(
        self, *, source_id: str, parent_source_ids: Sequence[str], description: str
    ) -> SourceRecord:
        """Refuse: an operator that produces no observation has no provenance node to mint.

        Raises:
            ObservationRefused: Always.
        """
        raise ObservationRefused(
            ApplicabilityVerdict(
                operator_id=self.operator_id,
                state_id=_nonempty(source_id, what="source_id"),
                applicable=False,
                reason=f"{self.verdict.capability} produces no artifact, so it has no provenance",
                routes_to=self.verdict.routes_to,
            )
        )


def unimplemented_operator(modality: ObservationModality) -> UnimplementedOperator:
    """Return the declared, refusing operator for an unimplemented modality.

    Args:
        modality: The modality to look up.

    Returns:
        The :class:`UnimplementedOperator` carrying that modality's capability verdict.

    Raises:
        KeyError: If the modality is implemented, which is the honest answer: asking for the
            unimplemented form of a modality that exists is a caller error, not a gap.
    """
    if not isinstance(modality, ObservationModality):
        raise TypeError("modality must be an ObservationModality")
    if modality not in UNIMPLEMENTED_OPERATORS:
        raise KeyError(
            f"{modality.value} is implemented; see FluorescenceOperator / SummaryProjectionOperator"
        )
    return UnimplementedOperator(
        operator_id=reserved_operator_id(modality), verdict=UNIMPLEMENTED_OPERATORS[modality]
    )


def reserved_operator_id(modality: ObservationModality) -> str:
    """Return the one operator id a blocked modality is allowed to carry."""
    return f"{RESERVED_OPERATOR_ID_PREFIX}{modality.value}"


class OperatorRegistrationRefused(TypeError):
    """Raised when an operator table is asked to admit an operator it must not admit.

    A :class:`TypeError` because the object is the wrong KIND of operator for the claim it makes, and
    because the protocol check this replaces already raised one.
    """


def validate_operator_registration(operator: object) -> ObservationOperator:
    """Judge whether ``operator`` may be admitted to an operator table, at the REGISTRY boundary.

    The protocol is structural, so satisfying it is not evidence about behaviour: an external review
    (2026-07-29, finding 3) subclassed :class:`SummaryProjectionOperator`, overrode ``modality`` to
    ``TRACTION`` and ``operator_id`` to ``unimplemented::traction``, passed ``isinstance`` against
    :class:`ObservationOperator`, was carried in a :class:`~aleph.virtual_cell.cell_state_posterior.\
CellStatePosterior` and returned the ``reduced-model`` scalar ``5.0`` for a modality
    :data:`UNIMPLEMENTED_OPERATORS` calls ``BLOCKED``.  ``UnimplementedOperator.observe`` raising
    unconditionally could not prevent that, because the forged operator was never an
    ``UnimplementedOperator``: the guarantee was class-local and the hole was at registration.

    Two things are therefore checked HERE, where a table admits an operator:

    * an operator whose modality is blocked must be the genuine :class:`UnimplementedOperator`
      carrying that modality's own :class:`ObservationCapabilityVerdict` — not a subclass of it and
      not something else wearing the modality;
    * the ``unimplemented::`` id namespace is reserved, so an operator using it must be the genuine
      article registered under exactly its own reserved id.

    Args:
        operator: The candidate operator.

    Returns:
        The same object, typed as an :class:`ObservationOperator`.

    Raises:
        OperatorRegistrationRefused: If the object does not satisfy the protocol, claims a blocked
            modality without being the declared refusal, or claims the reserved id namespace.
    """
    if not isinstance(operator, ObservationOperator):
        raise OperatorRegistrationRefused(
            f"{type(operator).__name__} does not satisfy the ObservationOperator protocol"
        )
    modality = operator.modality
    if not isinstance(modality, ObservationModality):
        raise OperatorRegistrationRefused(
            f"{operator.operator_id!r} declares a modality that is not an ObservationModality"
        )
    genuine = type(operator) is UnimplementedOperator
    operator_id = operator.operator_id
    if modality in UNIMPLEMENTED_OPERATORS:
        expected = UNIMPLEMENTED_OPERATORS[modality]
        if not genuine or operator.verdict != expected:  # type: ignore[attr-defined]
            raise OperatorRegistrationRefused(
                f"{operator_id!r} claims modality {modality.value!r}, which the capability table "
                f"declares {expected.status.value}: {expected.capability}. A blocked modality may "
                "only be carried by the declared UnimplementedOperator holding that modality's own "
                "capability verdict, because that is the only operator that cannot return a number. "
                f"Registering a {type(operator).__name__} here would turn a named gap back into a "
                f"silent estimate; route to: {expected.routes_to}"
            )
        if operator_id != reserved_operator_id(modality):
            raise OperatorRegistrationRefused(
                f"the declared refusal for {modality.value!r} is registered as "
                f"{reserved_operator_id(modality)!r}, not {operator_id!r}; a second id for one gap "
                "makes the gap countable twice"
            )
    elif operator_id.startswith(RESERVED_OPERATOR_ID_PREFIX):
        raise OperatorRegistrationRefused(
            f"the {RESERVED_OPERATOR_ID_PREFIX!r} id namespace is reserved for the declared "
            f"UnimplementedOperator refusals; {operator_id!r} is a "
            f"{type(operator).__name__} for modality {modality.value!r} and would inherit the "
            "credibility of a name that promises to return nothing"
        )
    return operator
