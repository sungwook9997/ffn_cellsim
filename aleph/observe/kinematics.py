"""Spreading rate, volume flux and PIV speed — and the unit discipline they are exported under.

What this module exports
------------------------
Three kinematic observables that a migration or spreading assay reports, all of them time
derivatives of something the state already carries:

* **projected-area rate** ``d(A)/dt`` [um^2/s] — how fast a cell spreads
* **volume flux** ``d(V)/dt`` [um^3/s] — signed, because it is *negative* during Y-27632 spreading
  and taking a magnitude would erase the finding
* **PIV field speed** — the mean speed of a velocity field [um/s]

Why the unit discipline is a type here and a convention nowhere
----------------------------------------------------------------
The external records these are compared against are incomplete in a specific, recoverable way. The
eLife 72381 report gives ``initial_spreading_area_rate = 34.18`` and ``initial_volume_flux =
-16.23`` and never states what they are in. Codex's ``build_tfm_manifest.py`` is explicit about the
same problem in its own data — it tags every field ``source_unit_unknown`` and says it "never infers
missing units".

A number with no unit is still usable, and it is usable for exactly two things: its **sign** and its
**ratio to another number from the same series**. Both survive multiplication by an unknown positive
constant. A magnitude does not.

So this module refuses to compare magnitudes against a unit-unknown external contrast, by name and
as a return value, and offers :func:`direction_agreement` instead. That is not caution. If the
comparison were merely discouraged, the first person under deadline would do it, the arithmetic
would succeed, the answer would be plausible, and the invented unit would be undiscoverable
afterwards.

And :class:`DirectionAgreement` says on its own face that it is **not an acceptance criterion**,
because the second failure mode is a direction check quietly hardening into a gate.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any, Final

import numpy as np

from aleph.observe.manifest import Declaration, ManifestValue, ProtocolManifest
from aleph.observe.operator import (
    AcceptedState,
    ApparentQuantity,
    ApplicabilityDomain,
    ObservationContext,
    RawObservable,
    Refusal,
    RefusalCode,
)
from aleph.observe.stationarity import (
    effective_sample_size,
    integrated_autocorrelation_time,
    standard_error_of_correlated_mean,
)

__all__ = [
    "MIN_RATE_FRAMES",
    "DirectionAgreement",
    "ExternalContrast",
    "KinematicRateOperator",
    "MagnitudeAgreement",
    "PivFieldSpeedOperator",
    "UnitProvenance",
    "direction_agreement",
    "kinematic_manifest",
    "magnitude_agreement",
    "piv_field_speed_operator",
    "spreading_area_rate_operator",
    "volume_flux_operator",
]

#: Three frames is the floor for a central difference with one interior point. Below it there is no
#: derivative to average and the "rate" would be a single forward difference wearing a mean's name.
MIN_RATE_FRAMES: Final[int] = 3


class UnitProvenance(StrEnum):
    """How well an external series' unit is known. There is deliberately no third value.

    ``SOURCE_UNIT_UNKNOWN``
        The source published a number and no unit. This is the state of the eLife 72381 spreading
        and volume-flux tables and of every traction field in Codex's TFM manifest. It is not a
        defect in the import — it is a faithful record of a defect in the publication — and it must
        survive into every comparison made downstream.

    Any other value is the unit itself, as text, exactly as the operators carry theirs.
    """

    SOURCE_UNIT_UNKNOWN = "source_unit_unknown"


@dataclass(frozen=True, slots=True)
class ExternalContrast:
    """A treated-minus-control effect from outside this repository, with what qualifies it.

    Attributes:
        name: The external observable's own name.
        mean_difference: Treated minus control, in :attr:`units`.
        ci95: The published 95% interval, in :attr:`units`.
        units: The unit as text, or :attr:`UnitProvenance.SOURCE_UNIT_UNKNOWN`.
        inference_unit: What the interval is over — "cell", "author_independent_experiment". Carried
            because an interval over experiments and an interval over cells differ by roughly an
            order of magnitude and the difference is invisible in the number itself.
        n_inference_units: How many of them. Three matched experiments is not three hundred cells.
        source: Accession or DOI.
    """

    name: str
    mean_difference: float
    ci95: tuple[float, float]
    units: str | UnitProvenance
    inference_unit: str
    n_inference_units: int
    source: str

    @property
    def sign(self) -> int:
        """``+1``/``-1``/``0``. Zero when the interval covers zero — no direction was
        established."""
        low, high = float(self.ci95[0]), float(self.ci95[1])
        if low > 0.0:
            return 1
        if high < 0.0:
            return -1
        return 0

    @property
    def unit_is_known(self) -> bool:
        return self.units != UnitProvenance.SOURCE_UNIT_UNKNOWN


@dataclass(frozen=True, slots=True)
class DirectionAgreement:
    """Whether Aleph moved the same way the experiment did. Never whether it moved far enough.

    Attributes:
        agrees: Signs match.
        aleph_sign: ``+1``/``-1``/``0``.
        external_sign: The external contrast's established sign.
        external: The contrast compared against, on the record.
        basis: Why a sign is admissible where a magnitude is not.
        is_acceptance_criterion: Always ``False``, and a field rather than a docstring so that a
            caller that promotes it has to write the promotion down.
        admissible_use: What this verdict may be used for.
    """

    agrees: bool
    aleph_sign: int
    external_sign: int
    external: ExternalContrast
    basis: str
    is_acceptance_criterion: bool = False
    admissible_use: str = (
        "regression testing and qualitative direction only. This verdict is not an acceptance "
        "criterion and must not be used as a gate: it establishes that a sign was reproduced, and "
        "a sign is one bit"
    )


@dataclass(frozen=True, slots=True)
class MagnitudeAgreement:
    """Whether Aleph's effect lands inside the external interval. Only ever reachable in known
    units."""

    within_ci: bool
    aleph_effect: float
    units: str
    external: ExternalContrast
    z_like: float
    detail: Mapping[str, Any] = field(default_factory=dict)


def direction_agreement(
    *, aleph_effect: float, external: ExternalContrast
) -> DirectionAgreement | Refusal:
    """Compare signs. Admissible against a unit-unknown external series; a magnitude is not.

    Refuses when the external interval covers zero: there is then no external direction to agree
    with, and reporting agreement against a sign the data does not establish would manufacture a
    confirmation out of noise.
    """
    external_sign = external.sign
    if external_sign == 0:
        return Refusal(
            code=RefusalCode.UNSUPPORTED_QUANTITY,
            reason=(
                f"the external 95% interval for {external.name!r} is {list(external.ci95)} and "
                "covers zero, so the source establishes no direction. Agreement with a sign the "
                "data does not carry would be a confirmation manufactured out of noise"
            ),
            operator_id="aleph.observe.kinematics.direction_agreement/v1",
            detail={
                "external_name": external.name,
                "ci95": list(external.ci95),
                "source": external.source,
            },
        )

    value = float(aleph_effect)
    if not math.isfinite(value):
        return Refusal(
            code=RefusalCode.NON_FINITE_INPUT,
            reason=f"aleph_effect is {value!r}, which has no sign",
            operator_id="aleph.observe.kinematics.direction_agreement/v1",
            detail={"aleph_effect": value},
        )

    aleph_sign = 0 if value == 0.0 else (1 if value > 0.0 else -1)
    return DirectionAgreement(
        agrees=aleph_sign == external_sign,
        aleph_sign=aleph_sign,
        external_sign=external_sign,
        external=external,
        basis=(
            "a sign is invariant under multiplication by an unknown POSITIVE constant, so it "
            "survives an unrecorded unit. A magnitude does not, which is why magnitude_agreement "
            "refuses the same comparison"
        ),
    )


def magnitude_agreement(
    *, aleph_effect: float, aleph_units: str, external: ExternalContrast
) -> MagnitudeAgreement | Refusal:
    """Compare magnitudes — and refuse unless both sides state the same unit.

    Two refusals, and they are different findings:

    * the external series has **no** unit — nothing can be done here, the record is incomplete;
    * the two sides state **different** units — no conversion is attempted, deliberately. A unit
      conversion buried in a comparison function is a factor nobody will ever see again.
    """
    operator_id = "aleph.observe.kinematics.magnitude_agreement/v1"

    if not external.unit_is_known:
        return Refusal(
            code=RefusalCode.INCONSISTENT_INPUT,
            reason=(
                f"the external series {external.name!r} from {external.source} carries "
                f"units={UnitProvenance.SOURCE_UNIT_UNKNOWN.value!r}, and Aleph reports "
                f"{aleph_units!r}. Comparing the magnitudes would invent the missing unit, and the "
                "invention would be undiscoverable afterwards because the arithmetic succeeds. Use "
                "direction_agreement(), which needs only the sign"
            ),
            operator_id=operator_id,
            detail={
                # Named so the two refusals below can be told apart programmatically. They are
                # different findings: an absent unit is a defect in the SOURCE that no caller can
                # repair, a mismatch is a defect in the CALL that the caller can. A test that could
                # not distinguish them let a mutation through on an earlier draft of this module.
                "failure": "external_unit_absent",
                "external_name": external.name,
                "external_units": UnitProvenance.SOURCE_UNIT_UNKNOWN.value,
                "aleph_units": aleph_units,
                "source": external.source,
                "repairable_by_caller": False,
                "admissible_alternative": "aleph.observe.kinematics.direction_agreement",
            },
        )

    if str(external.units) != str(aleph_units):
        return Refusal(
            code=RefusalCode.INCONSISTENT_INPUT,
            reason=(
                f"Aleph reports {aleph_units!r} and {external.name!r} is in "
                f"{str(external.units)!r}. No conversion is attempted here on purpose: a factor "
                "applied inside a comparison is a factor nobody sees again"
            ),
            operator_id=operator_id,
            detail={
                "failure": "unit_mismatch",
                "aleph_units": aleph_units,
                "external_units": str(external.units),
                "repairable_by_caller": True,
            },
        )

    low, high = float(external.ci95[0]), float(external.ci95[1])
    half_width = 0.5 * (high - low)
    centre = float(external.mean_difference)
    return MagnitudeAgreement(
        within_ci=low <= float(aleph_effect) <= high,
        aleph_effect=float(aleph_effect),
        units=str(aleph_units),
        external=external,
        z_like=(float(aleph_effect) - centre) / half_width if half_width > 0.0 else math.inf,
        detail={
            "ci95": [low, high],
            "inference_unit": external.inference_unit,
            "n_inference_units": external.n_inference_units,
            "z_like_is_in_units_of": "half-widths of the external 95% interval, not sigma",
        },
    )


# -- the manifest --------------------------------------------------------------------------------


def kinematic_manifest(
    *,
    temperature: ManifestValue,
    medium: ManifestValue,
    substrate: ManifestValue,
    sampling_cadence: ManifestValue,
    fit_window: ManifestValue,
    resolution_limit: ManifestValue,
    specimen_identity: ManifestValue,
    batch: ManifestValue,
    citation_audit: ManifestValue,
    applicability_statement: ManifestValue,
    evidence_role: ManifestValue,
    operator_version: ManifestValue,
    segmentation: ManifestValue = (
        "NOT_RECORDED: a projected area, a volume and a velocity field all come out of an image, "
        "and where the edge was drawn sets the number. A rate reported without its segmentation "
        "cannot be pooled with one reported under a different threshold"
    ),
    calibration: ManifestValue = Declaration.NOT_RECORDED,
    osmotic_state: ManifestValue = Declaration.NOT_RECORDED,
    confinement: ManifestValue = Declaration.NOT_RECORDED,
    pharmacology: ManifestValue = Declaration.NOT_RECORDED,
    adhesion_state: ManifestValue = Declaration.NOT_RECORDED,
    loading_history: ManifestValue = Declaration.NOT_RECORDED,
    preconditioning: ManifestValue = Declaration.NOT_RECORDED,
) -> ProtocolManifest:
    """The manifest for a free-spreading or migration kinematic observable. All thirty-four fields.

    ``segmentation`` defaults to a spelled-out ``NOT_RECORDED`` sentence rather than the bare enum,
    for the same reason ``regularization`` does in
    :func:`~aleph.observe.traction.tfm_traction_manifest`:
    it is the single largest source of between-study disagreement for this modality and deserves to
    be legible rather than terse in the record.
    """
    na = Declaration.NOT_APPLICABLE
    return ProtocolManifest(
        # Geometry and control: nothing touches the cell. It spreads or migrates on its own.
        probe_geometry=na,
        probe_location=na,
        probe_direction=na,
        control_mode=na,
        contact_model=na,
        adhesion_state=adhesion_state,
        # Time and amplitude: unforced, so no rate, frequency, dwell, strain or preconditioning is
        # imposed. Only the observation cadence applies — the same shape as the passive protocol.
        rate=na,
        frequency=na,
        dwell=na,
        loading_history=loading_history,
        max_strain=na,
        preconditioning=preconditioning,
        sampling_cadence=sampling_cadence,
        # Environment: the substrate is required. A spreading rate on glass and one on a 1 kPa gel
        # are different quantities.
        temperature=temperature,
        medium=medium,
        osmotic_state=osmotic_state,
        substrate=substrate,
        confinement=confinement,
        pharmacology=pharmacology,
        raw_output_schema="per-frame scalar series and its central-difference time derivative",
        raw_output_units="declared per operator: um^2/s, um^3/s or um/s. Always per SECOND, never "
        "per frame — a rate per frame is the one error this module's rescale control exists for",
        calibration=calibration,
        segmentation=segmentation,
        inverse_model=(
            "none: a central-difference time derivative of a directly observed series. There is no "
            "ill-posed inversion, which is what separates this modality from TFM"
        ),
        regularization=na,
        fit_window=fit_window,
        noise_model=(
            "frame-to-frame scatter of the per-frame rate, with the sample count corrected to "
            "n_eff = N/(2*tau_samples); consecutive frames of a spreading cell are strongly "
            "correlated and sigma/sqrt(N) is never reported"
        ),
        resolution_limit=resolution_limit,
        specimen_identity=specimen_identity,
        batch=batch,
        citation_audit=citation_audit,
        applicability_statement=applicability_statement,
        evidence_role=evidence_role,
        operator_version=operator_version,
    )


# -- the operators ---------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class KinematicRateOperator:
    """Time derivative of a scalar per-frame series. Satisfies ``ObservationOperator``.

    Args:
        manifest: The protocol in force.
        quantity: What is being differentiated, for the result's name.
        owner: State owner. Default ``"cell"``.
        array_name: The per-frame scalar series.
        value_units: Units of the series itself, e.g. ``"um^2"``. The reported rate is this per
            second, and the ``/s`` is appended here rather than passed in so it cannot be forgotten.
    """

    manifest: ProtocolManifest
    quantity: str
    array_name: str
    value_units: str
    owner: str = "cell"

    @property
    def operator_id(self) -> str:
        return f"aleph.observe.kinematics.KinematicRateOperator[{self.quantity}]/v1"

    @property
    def required_state_owners(self) -> tuple[str, ...]:
        return (self.owner,)

    @property
    def rate_units(self) -> str:
        return f"{self.value_units}/s"

    @property
    def identifiability_directions(self) -> tuple[str, ...]:
        return (
            f"the sign and the time course of d({self.quantity})/dt: constrained",
            "NOT constrained: any single mechanism. A spreading rate integrates protrusion, "
            "adhesion turnover, membrane reservoir and osmotic flux at once",
            "NOT constrained in absolute magnitude against the external record, whose own units "
            "are source_unit_unknown — see magnitude_agreement, which refuses that comparison",
        )

    @property
    def applicability(self) -> ApplicabilityDomain:
        return ApplicabilityDomain(
            description=(
                "an unforced single cell observed at a fixed cadence over a declared window, with "
                "the segmentation that produced the series on the record"
            ),
            bounds={"sample_interval_s": (0.0, math.inf)},
            required_state_owners=(self.owner,),
        )

    def raw_observable(
        self, state: AcceptedState, context: ObservationContext
    ) -> RawObservable | Refusal:
        """Central differences in the interior, one-sided at the ends. Per SECOND, never per
        frame."""
        domain = self.applicability

        owner_refusal = domain.check_owners(self.operator_id, state.owners)
        if owner_refusal is not None:
            return owner_refusal

        if not (context.sample_interval_s > 0.0):
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    f"sample_interval_s = {context.sample_interval_s!r}; a rate divided by a "
                    "non-positive interval is not a rate, and defaulting it would silently report "
                    "the derivative per frame"
                ),
                operator_id=self.operator_id,
                detail={"sample_interval_s": context.sample_interval_s},
            )

        try:
            series = np.asarray(state.array(self.owner, self.array_name), dtype=np.float64)
        except KeyError as exc:
            return Refusal(
                code=RefusalCode.MISSING_STATE_OWNER,
                reason=str(exc),
                operator_id=self.operator_id,
                detail={"owner": self.owner, "array": self.array_name},
            )

        if series.ndim != 1:
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=f"expected a 1-D per-frame series; got shape {series.shape}",
                operator_id=self.operator_id,
                detail={"shape": list(series.shape)},
            )

        if series.size < MIN_RATE_FRAMES:
            return Refusal(
                code=RefusalCode.INSUFFICIENT_SAMPLES,
                reason=(
                    f"{series.size} frame(s); a rate needs at least {MIN_RATE_FRAMES}. Two frames "
                    "give one forward difference, which is a difference and not a mean rate"
                ),
                operator_id=self.operator_id,
                detail={"n_frames": int(series.size), "minimum": MIN_RATE_FRAMES},
            )

        if not np.all(np.isfinite(series)):
            n_bad = int(np.count_nonzero(~np.isfinite(series)))
            return Refusal(
                code=RefusalCode.NON_FINITE_INPUT,
                reason=(
                    f"{n_bad} non-finite sample(s) in {self.owner}.{self.array_name}. Refused "
                    "rather than masked: a differenced NaN contaminates both its neighbours"
                ),
                operator_id=self.operator_id,
                detail={"n_non_finite": n_bad},
            )

        dt = float(context.sample_interval_s)
        rate = np.gradient(series, dt)

        return RawObservable(
            name=f"{self.quantity}_rate",
            operator_id=self.operator_id,
            values=rate,
            units=self.rate_units,
            n_samples=int(series.size),
            manifest_hash=self.manifest.manifest_hash(),
            support={
                "time_s": np.arange(series.size, dtype=np.float64) * dt,
                "series": series,
            },
            detail={
                "sample_interval_s": dt,
                "series_units": self.value_units,
                "first_value": float(series[0]),
                "last_value": float(series[-1]),
            },
        )

    def apparent_quantity(self, raw: RawObservable) -> ApparentQuantity | Refusal:
        """Mean rate over the declared window, with the correlated error bar. Signed, never
        absolute."""
        rate = np.asarray(raw.values, dtype=np.float64)
        dt = float(raw.detail.get("sample_interval_s", 0.0))
        if not (dt > 0.0):
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason="the raw observable carries no positive sample_interval_s",
                operator_id=self.operator_id,
                detail={"sample_interval_s": dt},
            )

        tau = integrated_autocorrelation_time(rate, dt)
        n = int(rate.size)
        std = float(np.std(rate))

        return ApparentQuantity(
            name=f"apparent_{self.quantity}_rate",
            operator_id=self.operator_id,
            value=float(np.mean(rate)),
            units=self.rate_units,
            standard_error=standard_error_of_correlated_mean(std, n, tau.tau_samples),
            n_effective=effective_sample_size(n, tau.tau_samples),
            manifest_hash=raw.manifest_hash,
            analysis_chain=(
                f"central-difference time derivative of the {self.quantity} series, per second",
                "mean over the declared window, SIGNED (no magnitude is taken: the volume flux is "
                "negative during Y-27632 spreading and a magnitude would erase that)",
                "standard error of the correlated mean with a Sokal-windowed tau_int",
            ),
            identifiability=self.identifiability_directions,
            detail={
                "n_frames": n,
                "tau_samples": float(tau.tau_samples),
                "rate_std": std,
                "external_comparison": (
                    "direction only unless the external series states a unit; see "
                    "aleph.observe.kinematics.magnitude_agreement, which refuses otherwise"
                ),
            },
        )


def spreading_area_rate_operator(manifest: ProtocolManifest, owner: str = "cell"):
    """``d(projected area)/dt`` [um^2/s]."""
    return KinematicRateOperator(
        manifest=manifest,
        quantity="projected_area",
        array_name="projected_area_um2",
        value_units="um^2",
        owner=owner,
    )


def volume_flux_operator(manifest: ProtocolManifest, owner: str = "cell"):
    """``d(volume)/dt`` [um^3/s]. Signed: it is negative during Y-27632 spreading."""
    return KinematicRateOperator(
        manifest=manifest,
        quantity="volume",
        array_name="volume_um3",
        value_units="um^3",
        owner=owner,
    )


@dataclass(frozen=True, slots=True)
class PivFieldSpeedOperator:
    """Mean speed of a PIV velocity field. Satisfies ``ObservationOperator``.

    Reads ``(frames, nodes, dimensions)`` velocities and reports the mean speed per frame. Kept
    separate from :class:`KinematicRateOperator` because it differentiates nothing — the field is
    already a velocity — and pretending otherwise would put a spurious ``dt`` in the chain.
    """

    manifest: ProtocolManifest
    owner: str = "field"
    array_name: str = "velocity_um_per_s"

    @property
    def operator_id(self) -> str:
        return "aleph.observe.kinematics.PivFieldSpeedOperator/v1"

    @property
    def required_state_owners(self) -> tuple[str, ...]:
        return (self.owner,)

    @property
    def identifiability_directions(self) -> tuple[str, ...]:
        return (
            "the mean speed of the field: constrained",
            "NOT constrained: direction or coherence. A rotating field and a translating field of "
            "the same speed are the same number here",
            "NOT constrained in absolute magnitude against a PIV record whose velocity units are "
            "not stated — the brief for this lane names this case specifically",
        )

    @property
    def applicability(self) -> ApplicabilityDomain:
        return ApplicabilityDomain(
            description="a reconstructed velocity field with its units on the record",
            bounds={"sample_interval_s": (0.0, math.inf)},
            required_state_owners=(self.owner,),
        )

    def raw_observable(
        self, state: AcceptedState, context: ObservationContext
    ) -> RawObservable | Refusal:
        owner_refusal = self.applicability.check_owners(self.operator_id, state.owners)
        if owner_refusal is not None:
            return owner_refusal

        try:
            velocity = np.asarray(state.array(self.owner, self.array_name), dtype=np.float64)
        except KeyError as exc:
            return Refusal(
                code=RefusalCode.MISSING_STATE_OWNER,
                reason=str(exc),
                operator_id=self.operator_id,
                detail={"owner": self.owner, "array": self.array_name},
            )

        if velocity.ndim != 3:
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    f"expected (frames, nodes, dimensions); got shape {velocity.shape}. A 2-D "
                    "array is ambiguous between one frame of vectors and many frames of scalars"
                ),
                operator_id=self.operator_id,
                detail={"shape": list(velocity.shape)},
            )

        if not np.all(np.isfinite(velocity)):
            n_bad = int(np.count_nonzero(~np.isfinite(velocity)))
            return Refusal(
                code=RefusalCode.NON_FINITE_INPUT,
                reason=f"{n_bad} non-finite velocity component(s); refused rather than masked",
                operator_id=self.operator_id,
                detail={"n_non_finite": n_bad},
            )

        speed_per_frame = np.mean(np.linalg.norm(velocity, axis=-1), axis=-1)
        return RawObservable(
            name="piv_mean_field_speed",
            operator_id=self.operator_id,
            values=speed_per_frame,
            units="um/s",
            n_samples=int(velocity.shape[1]),
            manifest_hash=self.manifest.manifest_hash(),
            support={
                "time_s": np.arange(velocity.shape[0], dtype=np.float64)
                * float(context.sample_interval_s)
            },
            detail={
                "sample_interval_s": float(context.sample_interval_s),
                "n_frames": int(velocity.shape[0]),
                "n_nodes": int(velocity.shape[1]),
                "n_dimensions": int(velocity.shape[2]),
            },
        )

    def apparent_quantity(self, raw: RawObservable) -> ApparentQuantity | Refusal:
        speed = np.asarray(raw.values, dtype=np.float64)
        dt = float(raw.detail.get("sample_interval_s", 0.0))
        if not (dt > 0.0):
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason="the raw observable carries no positive sample_interval_s",
                operator_id=self.operator_id,
                detail={"sample_interval_s": dt},
            )
        tau = integrated_autocorrelation_time(speed, dt)
        n = int(speed.size)
        return ApparentQuantity(
            name="apparent_piv_mean_field_speed",
            operator_id=self.operator_id,
            value=float(np.mean(speed)),
            units="um/s",
            standard_error=standard_error_of_correlated_mean(
                float(np.std(speed)), n, tau.tau_samples
            ),
            n_effective=effective_sample_size(n, tau.tau_samples),
            manifest_hash=raw.manifest_hash,
            analysis_chain=(
                "per-node speed as the Euclidean norm of the velocity vector",
                "mean over nodes, per frame",
                "mean over frames, with the correlated standard error",
            ),
            identifiability=self.identifiability_directions,
            detail={"n_frames": n, "tau_samples": float(tau.tau_samples)},
        )


def piv_field_speed_operator(manifest: ProtocolManifest, owner: str = "field"):
    """Mean speed of a PIV velocity field [um/s]."""
    return PivFieldSpeedOperator(manifest=manifest, owner=owner)
