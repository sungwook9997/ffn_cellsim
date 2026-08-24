"""What must be declared before a measurement means anything.

A number without its protocol is not a measurement, it is a rumour. Two indentation moduli taken at
different loading rates on the same cell are not two estimates of one quantity — they are two
different quantities, and averaging them is nonsense dressed as statistics. The manuscript's §8.1
lists what has to be on the record before a reported number can be compared to anything, and
:class:`ProtocolManifest` is that list made into a type.

Thirty-four fields in six families. None of them has a default.

The distinction this module exists for
--------------------------------------
Most of the geometry-and-control family does not apply to a passive fluctuation measurement. There
is no probe, so there is no probe radius. That is a **property of the protocol**, and it is
categorically different from an indentation run where a probe certainly existed and nobody wrote
down its radius.

The first is a complete record. The second is a hole in one.

    Declaration.NOT_APPLICABLE   the protocol has no such quantity     -- an assertion
    Declaration.NOT_RECORDED     it existed and was not captured       -- a confession

If both were ``None``, or if the second could simply be omitted from the mapping, the difference
would vanish at exactly the moment someone needs it: when they are deciding whether two datasets can
be pooled. So neither is ``None``, and no field may be omitted.

Why no field has a default
--------------------------
Adding a field to :class:`ProtocolManifest` breaks every construction site in the repository. That
is the intended cost. A new declarable is a new question every existing protocol has to answer, and
a default answers it for them — wrongly, silently, and in the direction that makes the record look
complete. ``tests/observe/test_manifest.py::test_no_manifest_field_has_a_default`` enforces this
structurally so it survives the next person to add a field.

There is deliberately no generic "fill the remainder with NOT_APPLICABLE" helper. That helper is
precisely how an unrecorded field becomes an inapplicable one without anybody deciding that it is.
:func:`passive_fluctuation_manifest` spells out all thirty-four values in source, in one place, and
can be audited by reading it.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
from dataclasses import dataclass, fields
from enum import StrEnum
from typing import Final

__all__ = [
    "ANALYSIS_FIELDS",
    "ENVIRONMENT_FIELDS",
    "EVIDENCE_FIELDS",
    "FIELD_FAMILIES",
    "GEOMETRY_AND_CONTROL_FIELDS",
    "MANIFEST_SCHEMA_TAG",
    "OBSERVATION_FIELDS",
    "TIME_AND_AMPLITUDE_FIELDS",
    "Declaration",
    "ManifestValue",
    "ProtocolManifest",
    "passive_fluctuation_manifest",
]


class Declaration(StrEnum):
    """The two ways a field can carry no concrete value, kept apart on purpose.

    ``NOT_APPLICABLE``
        The protocol has no such quantity. A passive fluctuation measurement has no probe, so
        ``probe_geometry`` is not missing information — there is nothing there to miss. Declaring
        this is an *assertion about the protocol* and it makes the record complete.

    ``NOT_RECORDED``
        The quantity existed and nobody captured it. This is a defect in the record, and the whole
        point of naming it is that a reader can see the defect instead of inferring completeness
        from a field that merely looks filled in.

    Collapsing these two into ``None`` destroys the only distinction that decides whether two
    datasets may be pooled.
    """

    NOT_APPLICABLE = "NOT_APPLICABLE"
    NOT_RECORDED = "NOT_RECORDED"


#: What a manifest field may hold. Note the absence of ``None``: there is no third way to say
#: nothing, because the two above are the only two that mean anything different.
ManifestValue = str | int | float | bool | Declaration

#: Prefixed into the hashed payload so a future change to the field set cannot collide with the
#: present one. Bump this when fields are added, removed or renamed.
MANIFEST_SCHEMA_TAG: Final[str] = "aleph.observe.manifest/v1"


GEOMETRY_AND_CONTROL_FIELDS: Final[tuple[str, ...]] = (
    "probe_geometry",
    "probe_location",
    "probe_direction",
    "control_mode",
    "contact_model",
    "adhesion_state",
)

TIME_AND_AMPLITUDE_FIELDS: Final[tuple[str, ...]] = (
    "rate",
    "frequency",
    "dwell",
    "loading_history",
    "max_strain",
    "preconditioning",
    "sampling_cadence",
)

ENVIRONMENT_FIELDS: Final[tuple[str, ...]] = (
    "temperature",
    "medium",
    "osmotic_state",
    "substrate",
    "confinement",
    "pharmacology",
)

OBSERVATION_FIELDS: Final[tuple[str, ...]] = (
    "raw_output_schema",
    "raw_output_units",
    "calibration",
)

ANALYSIS_FIELDS: Final[tuple[str, ...]] = (
    "segmentation",
    "inverse_model",
    "regularization",
    "fit_window",
    "noise_model",
    "resolution_limit",
)

EVIDENCE_FIELDS: Final[tuple[str, ...]] = (
    "specimen_identity",
    "batch",
    "citation_audit",
    "applicability_statement",
    "evidence_role",
    "operator_version",
)

FIELD_FAMILIES: Final[dict[str, tuple[str, ...]]] = {
    "geometry_and_control": GEOMETRY_AND_CONTROL_FIELDS,
    "time_and_amplitude": TIME_AND_AMPLITUDE_FIELDS,
    "environment": ENVIRONMENT_FIELDS,
    "observation": OBSERVATION_FIELDS,
    "analysis": ANALYSIS_FIELDS,
    "evidence": EVIDENCE_FIELDS,
}


@dataclass(frozen=True, slots=True)
class ProtocolManifest:
    """Everything that must be declared before a measurement can be compared to anything.

    Every field is required and every field is a :data:`ManifestValue`. A field that does not apply
    to this protocol takes :attr:`Declaration.NOT_APPLICABLE`; a field that applies and was not
    captured takes :attr:`Declaration.NOT_RECORDED`. Those are different claims and they hash
    differently.

    Fields carrying a physical quantity carry it as text *with its unit written in*, e.g.
    ``"296.0 K"`` rather than ``296.0``. A bare float in a manifest is a unit waiting to be guessed.
    """

    # -- geometry and control -----------------------------------------------------------------
    probe_geometry: ManifestValue
    probe_location: ManifestValue
    probe_direction: ManifestValue
    control_mode: ManifestValue
    contact_model: ManifestValue
    adhesion_state: ManifestValue

    # -- time and amplitude -------------------------------------------------------------------
    rate: ManifestValue
    frequency: ManifestValue
    dwell: ManifestValue
    loading_history: ManifestValue
    max_strain: ManifestValue
    preconditioning: ManifestValue
    sampling_cadence: ManifestValue

    # -- environment --------------------------------------------------------------------------
    temperature: ManifestValue
    medium: ManifestValue
    osmotic_state: ManifestValue
    substrate: ManifestValue
    confinement: ManifestValue
    pharmacology: ManifestValue

    # -- observation --------------------------------------------------------------------------
    raw_output_schema: ManifestValue
    raw_output_units: ManifestValue
    calibration: ManifestValue

    # -- analysis -----------------------------------------------------------------------------
    segmentation: ManifestValue
    inverse_model: ManifestValue
    regularization: ManifestValue
    fit_window: ManifestValue
    noise_model: ManifestValue
    resolution_limit: ManifestValue

    # -- evidence -----------------------------------------------------------------------------
    specimen_identity: ManifestValue
    batch: ManifestValue
    citation_audit: ManifestValue
    applicability_statement: ManifestValue
    evidence_role: ManifestValue
    operator_version: ManifestValue

    def __post_init__(self) -> None:
        bad: list[str] = []
        for field in fields(self):
            value = getattr(self, field.name)
            if value is None or not isinstance(value, str | int | float | bool):
                # ``Declaration`` is a ``StrEnum``, so it passes the ``str`` arm above.
                bad.append(f"{field.name}={value!r}")
        if bad:
            raise TypeError(
                "every manifest field must carry a concrete value or a Declaration; "
                "None is not a permitted value because it cannot distinguish "
                "NOT_APPLICABLE from NOT_RECORDED. Offenders: " + ", ".join(bad)
            )

    # -- views --------------------------------------------------------------------------------

    def as_mapping(self) -> dict[str, ManifestValue]:
        """The full field mapping, field order preserved. Every field is present, always."""
        return {f.name: getattr(self, f.name) for f in fields(self)}

    def by_family(self) -> dict[str, dict[str, ManifestValue]]:
        """The same content grouped into the six §8.1 families, for a human reading a record."""
        return {
            family: {name: getattr(self, name) for name in names}
            for family, names in FIELD_FAMILIES.items()
        }

    def not_recorded_fields(self) -> tuple[str, ...]:
        """Fields that apply to this protocol and were not captured — the holes in the record.

        This is the list a reviewer should look at first. It is deliberately not merged with
        :meth:`not_applicable_fields`.
        """
        return tuple(
            f.name for f in fields(self) if getattr(self, f.name) is Declaration.NOT_RECORDED
        )

    def not_applicable_fields(self) -> tuple[str, ...]:
        """Fields this protocol asserts have no referent. Not a defect; a description."""
        return tuple(
            f.name for f in fields(self) if getattr(self, f.name) is Declaration.NOT_APPLICABLE
        )

    def is_complete(self) -> bool:
        """True when nothing that applies was left uncaptured.

        A complete manifest may still be mostly ``NOT_APPLICABLE`` — the passive fluctuation
        protocol is, and that is the correct description of it.
        """
        return not self.not_recorded_fields()

    # -- identity -----------------------------------------------------------------------------

    def canonical_payload(self) -> str:
        """The exact string that is hashed. Exposed so a digest can be reproduced by hand.

        Sorted keys and fixed separators, so the digest depends on content and never on field
        declaration order or on JSON formatting.
        """
        body = {
            name: (value.value if isinstance(value, Declaration) else value)
            for name, value in self.as_mapping().items()
        }
        return json.dumps(
            {"schema": MANIFEST_SCHEMA_TAG, "fields": body},
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        )

    def manifest_hash(self) -> str:
        """SHA-256 of :meth:`canonical_payload`, hex.

        This is the join key between a reported number and the protocol that produced it. Its one
        required property is that manifests differing in *any* field — including differing only in
        ``NOT_APPLICABLE`` versus ``NOT_RECORDED`` — get different digests.
        """
        return hashlib.sha256(self.canonical_payload().encode("utf-8")).hexdigest()

    @property
    def short_hash(self) -> str:
        """First sixteen hex characters, for logs and filenames. Never the join key."""
        return self.manifest_hash()[:16]

    def amended(self, **changes: ManifestValue) -> ProtocolManifest:
        """A copy with named fields replaced — and therefore a different manifest hash.

        Amending is the honest way to record that a protocol changed. Mutating in place is not
        available: the manifest is frozen, because a number that travels with a hash must travel
        with the content that produced the hash.
        """
        unknown = set(changes) - {f.name for f in fields(self)}
        if unknown:
            raise TypeError(f"not manifest fields: {sorted(unknown)}")
        return dataclasses.replace(self, **changes)


def passive_fluctuation_manifest(
    *,
    temperature: ManifestValue,
    medium: ManifestValue,
    sampling_cadence: ManifestValue,
    fit_window: ManifestValue,
    resolution_limit: ManifestValue,
    specimen_identity: ManifestValue,
    batch: ManifestValue,
    citation_audit: ManifestValue,
    applicability_statement: ManifestValue,
    evidence_role: ManifestValue,
    operator_version: ManifestValue,
    osmotic_state: ManifestValue = Declaration.NOT_RECORDED,
    substrate: ManifestValue = Declaration.NOT_RECORDED,
    confinement: ManifestValue = Declaration.NOT_RECORDED,
    pharmacology: ManifestValue = Declaration.NOT_RECORDED,
) -> ProtocolManifest:
    """The manifest for `ALEPH-DQ-102`, the passive membrane thermal fluctuation spectrum.

    Every one of the thirty-four fields is written out below. That verbosity is the feature: this
    function is the audit trail for the claim "these fields do not apply to this protocol", and it
    can be checked by reading it against the manuscript's §8.1 list.

    The keyword arguments split into two groups, and the split is itself a claim:

    * **Required** — quantities that genuinely exist for a passive measurement and that a caller
      must supply. There is no default for ``temperature``, because a thermal fluctuation spectrum
      without a temperature is not underspecified, it is meaningless.
    * **Defaulting to** :attr:`Declaration.NOT_RECORDED` — environmental facts that exist for any
      real specimen but are routinely not captured. They default to the *confession*, never to the
      assertion, so a caller who does not think about substrate ends up with an honest record
      rather than a false claim that there was no substrate.

    Everything not in either group is :attr:`Declaration.NOT_APPLICABLE`, and each block below says
    why.

    Args:
        temperature: With units, e.g. ``"296.0 K"``. Sets ``kB*T``; the whole spectrum scales with it.
        medium: Bathing solution. Sets the viscosity and hence the relaxation times, which is what
            makes consecutive frames correlated (see :mod:`aleph.observe.stationarity`).
        sampling_cadence: Frame interval, with units. The one time-family field that applies: a
            passive measurement imposes no rate, but it does *sample* at one.
        fit_window: The wavevector band the reported spectrum is claimed over. Declared before the
            fit, because choosing it after seeing the spectrum is how a power law gets manufactured.
        resolution_limit: The mesh Nyquist statement, with units — see
            :mod:`aleph.observe.fluctuation` §2.3.
        specimen_identity: What was measured.
        batch: Preparation batch.
        citation_audit: Where the parameters came from, or that they came from nowhere.
        applicability_statement: The domain over which the result is claimed to hold.
        evidence_role: What this measurement is being used *for* — e.g. calibration versus test.
            The same number is not equally admissible in both roles.
        operator_version: Version of the code that produced it.
        osmotic_state: Defaults to ``NOT_RECORDED``.
        substrate: Defaults to ``NOT_RECORDED``.
        confinement: Defaults to ``NOT_RECORDED``.
        pharmacology: Defaults to ``NOT_RECORDED``.

    Returns:
        The frozen :class:`ProtocolManifest`.
    """
    na = Declaration.NOT_APPLICABLE
    return ProtocolManifest(
        # Geometry and control: there is no probe. Nothing touches the membrane, nothing is
        # displaced, nothing is held. Every field in this family is an assertion of absence, and
        # that is exactly why ALEPH-DQ-102 chose this modality first — it removes the three largest
        # sources of protocol ambiguity (probe geometry, contact model, adhesion) before the first
        # number is ever reported.
        probe_geometry=na,
        probe_location=na,
        probe_direction=na,
        control_mode=na,
        contact_model=na,
        adhesion_state=na,
        # Time and amplitude: a passive measurement imposes no loading. There is no rate, no drive
        # frequency, no dwell, no history, no strain amplitude, and nothing to precondition. The
        # sole exception is the cadence at which the observer samples, which is a property of the
        # observation and not of any forcing.
        rate=na,
        frequency=na,
        dwell=na,
        loading_history=na,
        max_strain=na,
        preconditioning=na,
        sampling_cadence=sampling_cadence,
        # Environment: all of it applies. Temperature and medium are required because they set the
        # fluctuation amplitude and the correlation time respectively; the rest default to the
        # confession rather than the assertion.
        temperature=temperature,
        medium=medium,
        osmotic_state=osmotic_state,
        substrate=substrate,
        confinement=confinement,
        pharmacology=pharmacology,
        # Observation: fixed by this operator, so they are stated rather than passed in. A caller
        # who needs different ones is running a different protocol and should say so with
        # ``amended()``, which changes the hash.
        raw_output_schema="mode power spectrum P(q) with per-mode sample counts",
        raw_output_units="P in um^2, q in 1/um, amplitudes in um, areas in um^2",
        calibration=na,  # no transducer; displacements are the state, in the state's own units
        # Analysis: there is no image to segment and no inverse problem to regularize. The spectrum
        # is a projection and a variance, both of which are exact operations on the samples. The
        # noise model is stated because the reported error bars depend on it.
        segmentation=na,
        inverse_model=na,
        regularization=na,
        fit_window=fit_window,
        noise_model="finite-sampling scatter of a variance estimate, with the sample count corrected "
        "to n_eff = T/(2*tau_int) because consecutive frames are correlated",
        resolution_limit=resolution_limit,
        # Evidence: all of it applies, always, to every measurement.
        specimen_identity=specimen_identity,
        batch=batch,
        citation_audit=citation_audit,
        applicability_statement=applicability_statement,
        evidence_role=evidence_role,
        operator_version=operator_version,
    )
