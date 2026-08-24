"""An external dataset, recorded so that a later reader can fetch the same bytes and judge them.

The experimental reader lane (`aleph.learn.lanes.MTG_PN_EXP`) is blocked on four things: dataset
identity, protocol, calibration, and split leakage.  Its gate's own rationale says why they are one
question rather than four — together they "decide whether the lane is measuring biology or measuring
its own batch structure."  All four are booleans on `GateEvidence` with nothing behind them, so
asserting any of them today would be a claim with no artifact.  This is the artifact.

One manifest per external source answers all four for that source.  Nothing here trains anything or
touches the seal; it records what a source is, in a form a guard can check.

What is deliberately strict
---------------------------
**No field has a default.**  A dataset whose temperature was never reported must say
``Declaration.NOT_RECORDED``, which hashes differently from ``NOT_APPLICABLE`` and differently again
from a plausible number somebody filled in.  This is the rule `ExperimentContextManifest` already
applies to temperature, and for the same reason: a default laundered into a manifest is an assumption
that becomes invisible one step downstream.

**Units are carried twice, as published and as converted.**  The conversion is where errors hide, and
keeping the published form lets a reader recompute it.  This project has already paid for that once:
handing the vesicle oracle ``4πR²`` where it wanted ``R²`` produces a *perfect* fit, χ² of 7e-23, with
κ wrong by exactly ``1/(4π)``.  A wrong unit conversion does not look wrong.

**Calibration has three states, not two.**  Resolved, free parameter, or not recorded.  The middle
one is legitimate and it is not free: measured on the closed-vesicle oracle, the passive fluctuation
spectrum has Jacobian rank 2 over ``(κ, σ, R, T)``, and an unknown amplitude calibration is an exact
degeneracy — ``(A, κ, σ) → (fA, fκ, fσ)`` leaves every prediction unchanged.  So declaring a
calibration free **consumes one of the two dimensions the observation has**, and
:meth:`ExternalDatasetManifest.dimensions_consumed` makes that cost countable instead of leaving it
to be discovered during inference.

**Split keys are the four §8.4 names.**  `LaneSpec.fit` takes `split_declaration` as a required
argument because the rule is a hard constraint: entire donor, lab, batch and perturbation held out,
random neighbouring-frame splits prohibited.  A manifest that cannot name those keys cannot be split
correctly, and this repository has already met the failure — the vesicle identifiability result
records that its frames are treated as independent while real frames correlate over the mode
relaxation time, longest exactly where κ information is scarcest.

Direction of dependency
-----------------------
This module does **not** import `aleph.learn`.  `aleph/learn/evaluation.py` already imports
`aleph.evidence.claim`, so the reverse would be circular; more to the point, a record of what a
dataset is should not know that a neural reader exists.  :meth:`ExternalDatasetManifest.gate_readiness`
returns a plain verdict and a caller maps it onto whatever gate it is answering.
"""

from __future__ import annotations

import importlib
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Final

from aleph.observe.manifest import Declaration

__all__ = [
    "AxisRef",
    "Calibration",
    "CalibrationStatus",
    "DatasetManifestError",
    "ExternalDatasetManifest",
    "GateReadiness",
    "SPLIT_KEYS",
]

#: The four keys §8.4 requires an experimental split to hold out whole.
SPLIT_KEYS: Final[tuple[str, ...]] = ("donor", "lab", "batch", "perturbation")

_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_ISO_UTC = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")


class DatasetManifestError(ValueError):
    """A manifest is incomplete or inconsistent, and nothing derived from it may be quoted."""


class CalibrationStatus(StrEnum):
    """How raw instrument units reach physical units, and what each answer costs."""

    RESOLVED = "resolved"
    """A factor and its uncertainty are known. Costs no dimension."""

    FREE_PARAMETER = "free_parameter"
    """Legitimate, and it consumes one dimension of what the observation can distinguish."""

    NOT_RECORDED = "not_recorded"
    """Nobody captured it. Not the same claim as declaring it free, and it hashes differently."""


@dataclass(frozen=True, slots=True, kw_only=True)
class Calibration:
    """The conversion from recorded units to physical units, and its evidence.

    Attributes:
        status: Which of the three answers this dataset gives.
        factor: Multiplicative factor into the unit named by the manifest's conversion, or None.
        relative_uncertainty: Fractional uncertainty on the factor, or None.
        note: Where the factor came from, or why it is free or unrecorded.
    """

    status: CalibrationStatus
    factor: float | None = None
    relative_uncertainty: float | None = None
    note: str = ""

    def __post_init__(self) -> None:
        if self.status is CalibrationStatus.RESOLVED:
            if self.factor is None or self.relative_uncertainty is None:
                raise DatasetManifestError(
                    "A resolved calibration needs both a factor and its relative uncertainty. A "
                    "factor without an uncertainty is a point estimate presented as exact, which is "
                    "the shape of error this manifest exists to prevent."
                )
            if not (self.factor > 0.0):
                raise DatasetManifestError(f"calibration factor must be positive; got {self.factor}")
            if not (0.0 <= self.relative_uncertainty < 1.0):
                raise DatasetManifestError(
                    "relative_uncertainty must lie in [0, 1); a fractional uncertainty at or above "
                    f"one means the factor's sign is undetermined. Got {self.relative_uncertainty}"
                )
        elif self.factor is not None or self.relative_uncertainty is not None:
            raise DatasetManifestError(
                f"A {self.status} calibration carries no factor. Supplying one says the calibration "
                "is known while declaring that it is not."
            )
        if not self.note.strip():
            raise DatasetManifestError(
                "Calibration.note must say where the factor came from, or why it is free or "
                "unrecorded. An unexplained calibration cannot be checked by anyone."
            )

    @property
    def dimensions_consumed(self) -> int:
        """One when the calibration is a free parameter, zero otherwise."""
        return 1 if self.status is CalibrationStatus.FREE_PARAMETER else 0


@dataclass(frozen=True, slots=True, kw_only=True)
class AxisRef:
    """One axis of the registry, named so it can be resolved rather than trusted.

    Attributes:
        carrier: Dotted path of the dataclass the axis lives on.
        field: Field name on it.
    """

    carrier: str
    field: str

    def __post_init__(self) -> None:
        module_path, _, type_name = self.carrier.rpartition(".")
        if not module_path or not type_name:
            raise DatasetManifestError(f"AxisRef.carrier must be a dotted path; got {self.carrier!r}")
        try:
            module = importlib.import_module(module_path)
        except ModuleNotFoundError as missing:
            raise DatasetManifestError(
                f"AxisRef names {self.carrier!r} but {module_path!r} does not import: {missing}"
            ) from missing
        carrier_type = getattr(module, type_name, None)
        if carrier_type is None:
            raise DatasetManifestError(f"{module_path!r} has no {type_name!r}")
        import dataclasses as _dc

        if not _dc.is_dataclass(carrier_type):
            raise DatasetManifestError(f"{self.carrier} is not a dataclass and carries no axes")
        names = {f.name for f in _dc.fields(carrier_type)}
        if self.field not in names:
            raise DatasetManifestError(
                f"{self.carrier} has no field {self.field!r}. An axis name that does not resolve "
                "makes a dataset silently inform nothing, which reads as a dataset that informs "
                "everything it was filed under."
            )

    def __str__(self) -> str:
        return f"{self.carrier}.{self.field}"


@dataclass(frozen=True, slots=True, kw_only=True)
class GateReadiness:
    """Which of the experimental lane's four preconditions this manifest answers.

    Deliberately not a `GateEvidence`: that type lives in the sealed reader package, which already
    imports this one's neighbours. A caller maps this onto whatever gate it is answering.
    """

    dataset_identity: bool
    protocol: bool
    calibration: bool
    split_leakage: bool
    unmet: tuple[str, ...]

    @property
    def complete(self) -> bool:
        return not self.unmet


@dataclass(frozen=True, slots=True, kw_only=True)
class ExternalDatasetManifest:
    """One external dataset, identified, converted, calibrated and keyed for splitting.

    Attributes:
        dataset_id: Stable local name. Appears in every record derived from this source.
        accession: DOI, repository accession or URL. What a reader fetches.
        retrieved_utc: When it was fetched, ``YYYY-MM-DDTHH:MM:SSZ``.
        content_digest: SHA-256 of the retrieved payload, lowercase hex.
        axes_informed: Registry axes this source carries evidence about.
        published_units: Axis field name to the unit **as the source printed it**.
        unit_conversion: Axis field name to the conversion into the ``pN-µm-s-K`` system.
        protocol_manifest_digest: Digest of the filled `ProtocolManifest`, or a `Declaration`.
        calibration: How recorded units reach physical units.
        split_keys: The four §8.4 keys, each a value or a `Declaration`.
        note: Anything else worth carrying.
    """

    dataset_id: str
    accession: str
    retrieved_utc: str
    content_digest: str
    axes_informed: tuple[AxisRef, ...]
    published_units: Mapping[str, str]
    unit_conversion: Mapping[str, str]
    protocol_manifest_digest: str | Declaration
    calibration: Calibration
    split_keys: Mapping[str, str | Declaration]
    note: str = field(default="")

    def __post_init__(self) -> None:
        for name in ("dataset_id", "accession", "retrieved_utc", "content_digest"):
            if not str(getattr(self, name)).strip():
                raise DatasetManifestError(
                    f"ExternalDatasetManifest.{name} must be a non-empty string. A source that "
                    "cannot be fetched again is an assertion, not evidence."
                )
        if not _ISO_UTC.match(self.retrieved_utc):
            raise DatasetManifestError(
                f"retrieved_utc must be YYYY-MM-DDTHH:MM:SSZ; got {self.retrieved_utc!r}"
            )
        if not _DIGEST.match(self.content_digest):
            raise DatasetManifestError(
                f"content_digest must be 64 lowercase hex characters; got {self.content_digest!r}"
            )
        if not self.axes_informed:
            raise DatasetManifestError(
                f"{self.dataset_id}: a dataset informing no axis has nothing to contribute and "
                "should not be filed."
            )

        informed = {ref.field for ref in self.axes_informed}
        for label, table in (("published_units", self.published_units),
                             ("unit_conversion", self.unit_conversion)):
            covered = set(table)
            if covered != informed:
                missing = sorted(informed - covered)
                extra = sorted(covered - informed)
                raise DatasetManifestError(
                    f"{self.dataset_id}: {label} must name exactly the informed axes. "
                    f"Missing {missing}; unexpected {extra}. A unit table that does not line up "
                    "with the axes leaves a conversion to be guessed at use."
                )
            for axis, value in table.items():
                if not str(value).strip():
                    raise DatasetManifestError(f"{self.dataset_id}: {label}[{axis!r}] is blank")

        declared_keys = set(self.split_keys)
        if declared_keys != set(SPLIT_KEYS):
            raise DatasetManifestError(
                f"{self.dataset_id}: split_keys must name exactly {list(SPLIT_KEYS)}; got "
                f"{sorted(declared_keys)}. §8.4 holds out entire donor, lab, batch and "
                "perturbation, so a key that is absent cannot be held out."
            )
        for key, value in self.split_keys.items():
            if isinstance(value, Declaration):
                continue
            if not str(value).strip():
                raise DatasetManifestError(
                    f"{self.dataset_id}: split key {key!r} is blank. Say NOT_RECORDED or "
                    "NOT_APPLICABLE — the two are different claims and a blank is neither."
                )

    @property
    def dimensions_consumed(self) -> int:
        """Dimensions this dataset's unresolved calibration costs the inference."""
        return self.calibration.dimensions_consumed

    def gate_readiness(self) -> GateReadiness:
        """Which of the experimental lane's four preconditions this manifest answers.

        **Identity is never answered here.** An audit demonstrated why: a `content_digest` of
        ``"a" * 64`` satisfies the format check and is the digest of nothing, so an earlier version
        of this method returning ``dataset_identity=True`` was trusting the manifest author's string
        rather than verifying anything. Identity needs a fetch of the payload and a comparison
        against this digest, which is a driver's job and not a dataclass's; until that verification
        exists, this reports ``False`` and lists it unmet.
        """
        unmet: list[str] = ["dataset_identity"]

        protocol = not isinstance(self.protocol_manifest_digest, Declaration)
        if protocol and not _DIGEST.match(str(self.protocol_manifest_digest)):
            protocol = False
        if not protocol:
            unmet.append("protocol_resolved")

        calibration = self.calibration.status is not CalibrationStatus.NOT_RECORDED
        if not calibration:
            unmet.append("calibration_resolved")

        split = all(not isinstance(value, Declaration) for value in self.split_keys.values())
        if not split:
            unmet.append("split_leakage_resolved")

        return GateReadiness(
            dataset_identity=False,
            protocol=protocol,
            calibration=calibration,
            split_leakage=split,
            unmet=tuple(unmet),
        )
