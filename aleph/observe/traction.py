"""Traction force per cell, from TFM — an operator built to refuse, and to stop refusing by itself.

Why this is written now, while it cannot answer
-----------------------------------------------
Traction is the observable a mechanics engine is most obviously judged by, and it is the one an
external sweep will reach for first. It is also the one this engine cannot currently supply, because
two defects are open on the load path it reads:

**(a) the NMII actin-side force is discarded.** ``aleph/scenarios/whole_cell_couplings.py:735`` is
``return builder(owner).sites``: the ``ProtrusionActinTarget`` wrapper is unwrapped and thrown away,
so ``deliver_forces_pn`` is never called. Measured at 0 calls in 12 steps.

**(b) the delivery probe reports a false positive.** ``aleph/scenarios/whole_cell.py:1643``
``_probe_delivery`` scores a connector row as *delivering* when any ledger pair carries a non-zero
force. That is a statement about bookkeeping, not about motion. Four NMII rows book 1.800 pN and
displace their owner by exactly ``0.000000e+00`` um under a 1000x head-strain amplification, and the
probe scores all four as delivering.

A traction number computed over either defect looks like traction and is not. If a sweep starts on
it, everything downstream is trained on a broken scale — and the failure is silent, because the
number has the right units and a plausible magnitude.

So: the operator is built, its arithmetic is under control, and it **returns a refusal** until the
path is sound. The refusal is not a hand-set flag. It is a runtime check against audit arrays the
state carries, which means that when the owning lanes fix the engine, this operator begins reporting
without a line of it being edited. That is the difference between a guard and a comment.

How the precondition is checkable without importing the engine
--------------------------------------------------------------
:mod:`aleph.observe` does not import :mod:`aleph.scenarios` and must not start. It does not need to:
both defects have the same *state-visible* signature — **a row that books force and moves nobody** —
and that is two arrays and a comparison. The audit therefore lives under its own state owner,
:data:`LOAD_PATH_AUDIT_OWNER`, kept separate from the physics arrays so that nobody can mistake it
for an input to the number.

On the external record, and why it is not an acceptance criterion
-----------------------------------------------------------------
The external TFM tables reach this repository through Codex's ``build_tfm_manifest.py``, which tags
every traction field ``source_unit_unknown`` and states that it "never infers missing units". The
external traction magnitudes are therefore usable as **directions and ratios** and are not pN. This
operator reports pN and Pa; comparing the two as absolute quantities would invent the missing unit.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Final

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

__all__ = [
    "DEFAULT_BOOKED_FORCE_FLOOR_PN",
    "LOAD_PATH_AUDIT_ARRAYS",
    "LOAD_PATH_AUDIT_OWNER",
    "LoadPathPrecondition",
    "TractionPerCellOperator",
    "tfm_traction_manifest",
]

_OPERATOR_ID: Final[str] = "aleph.observe.traction.TractionPerCellOperator/v1"

#: State owner carrying the load-path audit. Deliberately **not** the physics owner: an audit array
#: that lived beside ``traction_force_pn`` would eventually be summed into something.
LOAD_PATH_AUDIT_OWNER: Final[str] = "load_path_audit"

#: The three arrays the audit must supply, one entry per connector row feeding the traction.
#:
#: ``booked_force_pn``         what the ledger recorded for the row [pN]
#: ``owner_displacement_um``   how far the row's owner actually moved [um]
#: ``deliver_call_count``      how many times the row's delivery entry point was called
LOAD_PATH_AUDIT_ARRAYS: Final[tuple[str, ...]] = (
    "booked_force_pn",
    "owner_displacement_um",
    "deliver_call_count",
)

#: Below this, a booked force is round-off and carries no claim that anything should have moved.
#: Set well under the 1.800 pN the four NMII rows book, so the open defect is caught, and well above
#: the 1e-30-to-1e-14 round-off band the scenario census measured on three build-time rows.
DEFAULT_BOOKED_FORCE_FLOOR_PN: Final[float] = 1.0e-9

#: TFM needs a substrate that visibly deforms. The band is a *declaration of where this operator is
#: entitled to speak*, not a claim about where traction microscopy stops working: gels in the
#: literature run from a few hundred Pa to tens of kPa, and glass at ~69 GPa is outside it by six
#: orders of magnitude.
DEFAULT_MIN_SUBSTRATE_PA: Final[float] = 1.0e2
DEFAULT_MAX_SUBSTRATE_PA: Final[float] = 1.0e5


@dataclass(frozen=True, slots=True)
class LoadPathPrecondition:
    """One named, runtime-checkable condition the load path must satisfy before traction means
    anything.

    Attributes:
        name: Stable identifier, reported in the refusal's ``detail``.
        why: One sentence a human can act on, naming the defect.
        owning_lanes: Which session rows own the fix. A refusal that cannot be routed is a log line.
        check: ``(audit_arrays, row_names) -> tuple of offending row labels``. Empty means it
        passes.
    """

    name: str
    why: str
    owning_lanes: tuple[str, ...]
    check: Callable[[Mapping[str, np.ndarray], Sequence[str]], tuple[str, ...]]


def _rows_that_never_deliver(
    audit: Mapping[str, np.ndarray], row_names: Sequence[str]
) -> tuple[str, ...]:
    counts = np.asarray(audit["deliver_call_count"], dtype=np.float64)
    return tuple(row_names[i] for i in np.flatnonzero(counts <= 0.0))


def _rows_that_book_force_and_move_nobody(
    audit: Mapping[str, np.ndarray], row_names: Sequence[str]
) -> tuple[str, ...]:
    booked = np.asarray(audit["booked_force_pn"], dtype=np.float64)
    moved = np.asarray(audit["owner_displacement_um"], dtype=np.float64)
    offending = np.flatnonzero((booked > DEFAULT_BOOKED_FORCE_FLOOR_PN) & (moved <= 0.0))
    return tuple(row_names[i] for i in offending)


#: The two open defects, as checks. Ordered: the delivery-count check runs first because a row that
#: never delivers cannot meaningfully be asked whether its booked force moved anyone.
_DEFAULT_PRECONDITIONS: Final[tuple[LoadPathPrecondition, ...]] = (
    LoadPathPrecondition(
        name="every_wired_row_delivers_at_least_once",
        why=(
            "a wired row whose delivery entry point is never called contributes exactly zero to "
            "the traction while appearing in the census as present. aleph/scenarios/"
            "whole_cell_couplings.py:735 returns builder(owner).sites, unwrapping and discarding "
            "the ProtrusionActinTarget whose deliver_forces_pn would have been that call: measured "
            "at 0 calls in 12 steps"
        ),
        owning_lanes=("2729da45 (S-STROKE)", "46143f30 (Lane W2)"),
        check=_rows_that_never_deliver,
    ),
    LoadPathPrecondition(
        name="booked_force_moves_its_owner",
        why=(
            "a row that books force in the ledger and displaces its owner by exactly zero has a "
            "load path that is disconnected downstream of the ledger. aleph/scenarios/"
            "whole_cell.py:1643 _probe_delivery scores such a row as *delivering* because it reads "
            "the ledger pair and not the motion: four NMII rows book 1.800 pN and move their owner "
            "0.000000e+00 um under a 1000x head-strain amplification"
        ),
        owning_lanes=("2729da45 (S-STROKE)", "46143f30 (Lane W2)"),
        check=_rows_that_book_force_and_move_nobody,
    ),
)


def tfm_traction_manifest(
    *,
    temperature: ManifestValue,
    medium: ManifestValue,
    substrate: ManifestValue,
    control_mode: ManifestValue,
    sampling_cadence: ManifestValue,
    fit_window: ManifestValue,
    resolution_limit: ManifestValue,
    specimen_identity: ManifestValue,
    batch: ManifestValue,
    citation_audit: ManifestValue,
    applicability_statement: ManifestValue,
    evidence_role: ManifestValue,
    operator_version: ManifestValue,
    regularization: ManifestValue = (
        "NOT_RECORDED: TFM is ill-posed and every reported traction depends on the regulariser and "
        "its parameter. A record without one cannot be compared to a record with one"
    ),
    segmentation: ManifestValue = Declaration.NOT_RECORDED,
    calibration: ManifestValue = Declaration.NOT_RECORDED,
    adhesion_state: ManifestValue = Declaration.NOT_RECORDED,
    loading_history: ManifestValue = Declaration.NOT_RECORDED,
    preconditioning: ManifestValue = Declaration.NOT_RECORDED,
    osmotic_state: ManifestValue = Declaration.NOT_RECORDED,
    confinement: ManifestValue = Declaration.NOT_RECORDED,
    pharmacology: ManifestValue = Declaration.NOT_RECORDED,
) -> ProtocolManifest:
    """The manifest for a traction force microscopy measurement. All thirty-four fields.

    Two families separate this protocol sharply from the other two in :mod:`aleph.observe`:

    * **Analysis.** TFM is an *ill-posed inverse problem*. Substrate bead displacements are
      inverted to a traction field, and every published traction depends on the regulariser and its
      parameter. Both ``inverse_model`` and ``regularization`` therefore apply, and asserting
      ``NOT_APPLICABLE`` for them — as the passive fluctuation protocol correctly does — would be a
      false record. ``segmentation`` applies too: the cell footprint has to be cut out of an image.
    * **Geometry and control.** There is no probe, so the probe fields are genuine
      ``NOT_APPLICABLE`` assertions, but the *substrate* is the instrument and its modulus is not
      optional.

    Args:
        temperature: With units.
        medium: Bathing solution.
        substrate: The gel and its modulus, with units, and its coating. This is the instrument.
        control_mode: How the cell was held. A free-contraction traction and a stretched-substrate
            traction are different quantities.
        sampling_cadence: Frame interval, with units.
        fit_window: The footprint the traction is summed over.
        resolution_limit: Bead density and the traction below which nothing is resolved.
        specimen_identity: What was measured.
        batch: Preparation batch.
        citation_audit: Where the parameters came from, or that they came from nowhere.
        applicability_statement: The domain over which the result is claimed.
        evidence_role: What the number is used for.
        operator_version: Version of the code that produced it.
        regularization: Defaults to a spelled-out ``NOT_RECORDED`` sentence rather than the bare
            enum, because for TFM specifically an unrecorded regulariser is the single largest
            source of between-study disagreement and deserves to be legible in the record.
        segmentation: Defaults to ``NOT_RECORDED``.
        calibration: Bead tracking and the substrate modulus measurement. Defaults to
            ``NOT_RECORDED``.
        adhesion_state: Defaults to ``NOT_RECORDED``.
        loading_history: Defaults to ``NOT_RECORDED``.
        preconditioning: Defaults to ``NOT_RECORDED``.
        osmotic_state: Defaults to ``NOT_RECORDED``.
        confinement: Defaults to ``NOT_RECORDED``.
        pharmacology: Defaults to ``NOT_RECORDED``.

    Returns:
        The frozen :class:`~aleph.observe.manifest.ProtocolManifest`.
    """
    na = Declaration.NOT_APPLICABLE
    return ProtocolManifest(
        # Geometry and control: nothing indents the cell. The substrate is the transducer, and it is
        # recorded under `substrate` in the environment family where it belongs.
        probe_geometry=na,
        probe_location=na,
        probe_direction=na,
        control_mode=control_mode,
        contact_model=(
            "cell-substrate traction recovered from substrate displacement on a linearly elastic "
            "half-space; no indenter contact model is involved"
        ),
        adhesion_state=adhesion_state,
        # Time and amplitude: a free-contraction traction imposes no rate, frequency, dwell or
        # strain amplitude. The cell sets its own. Only the observation cadence applies.
        rate=na,
        frequency=na,
        dwell=na,
        loading_history=loading_history,
        max_strain=na,
        preconditioning=preconditioning,
        sampling_cadence=sampling_cadence,
        # Environment: all of it applies, and `substrate` is required rather than defaulted, because
        # a traction without a substrate modulus is not underspecified, it is meaningless.
        temperature=temperature,
        medium=medium,
        osmotic_state=osmotic_state,
        substrate=substrate,
        confinement=confinement,
        pharmacology=pharmacology,
        # Observation.
        raw_output_schema="per-adhesion-site traction force magnitude and site area",
        raw_output_units="force in pN, area in um^2, stress in Pa (1 pN/um^2 = 1 Pa exactly)",
        calibration=calibration,
        # Analysis: this is where TFM differs most from the passive protocol.
        segmentation=segmentation,
        inverse_model=(
            "substrate displacement -> traction, an ill-posed inversion. This operator reads a "
            "traction field that has ALREADY been inverted; it does not perform the inversion, and "
            "the regulariser that produced its input is recorded in `regularization`"
        ),
        regularization=regularization,
        fit_window=fit_window,
        noise_model=(
            "per-site traction uncertainty is NOT_RECORDED by this operator: the total is a sum of "
            "magnitudes over sites whose individual errors are correlated through the shared "
            "inversion, and a per-site independent error model would understate the total's error"
        ),
        resolution_limit=resolution_limit,
        # Evidence.
        specimen_identity=specimen_identity,
        batch=batch,
        citation_audit=citation_audit,
        applicability_statement=applicability_statement,
        evidence_role=evidence_role,
        operator_version=operator_version,
    )


@dataclass(frozen=True, slots=True)
class TractionPerCellOperator:
    """Total traction force a cell exerts on its substrate. Satisfies ``ObservationOperator``.

    Refuses with :attr:`~aleph.observe.operator.RefusalCode.UPSTREAM_DEFECT_OPEN` while any
    :class:`LoadPathPrecondition` fails, and begins reporting — with no edit to this file — as soon
    as the audit arrays say the path is sound.

    Args:
        manifest: The protocol in force. Build it with :func:`tfm_traction_manifest`.
        adhesion_owner: State owner holding the traction field. Default ``"adhesion"``.
        traction_array: Per-site traction force magnitude [pN]. Default ``"traction_force_pn"``.
        area_array: Per-site contact area [um^2]. Default ``"site_area_um2"``.
        audit_owner: State owner holding the load-path audit. Default :data:`LOAD_PATH_AUDIT_OWNER`.
        load_path_preconditions: The checks. Exposed as a field, and readable before any call, so a
            caller can see what would block it without first triggering a refusal.
        min_substrate_pa: Lower applicability bound on the substrate modulus.
        max_substrate_pa: Upper applicability bound.
    """

    manifest: ProtocolManifest
    adhesion_owner: str = "adhesion"
    traction_array: str = "traction_force_pn"
    area_array: str = "site_area_um2"
    audit_owner: str = LOAD_PATH_AUDIT_OWNER
    load_path_preconditions: tuple[LoadPathPrecondition, ...] = _DEFAULT_PRECONDITIONS
    min_substrate_pa: float = DEFAULT_MIN_SUBSTRATE_PA
    max_substrate_pa: float = DEFAULT_MAX_SUBSTRATE_PA

    # -- declarations ---------------------------------------------------------------------------

    @property
    def operator_id(self) -> str:
        return _OPERATOR_ID

    @property
    def required_state_owners(self) -> tuple[str, ...]:
        return (self.adhesion_owner, self.audit_owner)

    @property
    def identifiability_directions(self) -> tuple[str, ...]:
        """A whole-cell integral constrains a product, and it is worth writing down which one."""
        return (
            "the product (number of engaged adhesions) x (mean force per engaged adhesion): "
            "constrained, as their sum",
            "NOT constrained: that product's two factors separately. Halving the engaged count and "
            "doubling the per-adhesion force is exactly invariant, and those are different cells",
            "NOT constrained: which generator supplied the force. Motor stroke, passive crosslink "
            "tension and substrate pre-strain all arrive at the same total",
            "NOT constrained: anything, in absolute units, against the external TFM record — its "
            "own manifest tags every traction field source_unit_unknown, so that comparison is "
            "available as a ratio and a direction only",
        )

    @property
    def applicability(self) -> ApplicabilityDomain:
        return ApplicabilityDomain(
            description=(
                "a single isolated cell adhered to a linearly elastic compliant substrate whose "
                "modulus is on the record, with a load path that has been audited as delivering"
            ),
            bounds={
                "substrate_youngs_modulus_pa": (self.min_substrate_pa, self.max_substrate_pa),
            },
            required_state_owners=(self.adhesion_owner, self.audit_owner),
        )

    # -- the load-path audit ----------------------------------------------------------------------

    def audit_load_path(
        self, state: AcceptedState, row_names: Sequence[str]
    ) -> Refusal | None:
        """Run every precondition in order. ``None`` when the path is sound.

        Public because the audit is a finding in its own right: a caller may want to record *why*
        traction is unavailable without asking for traction.
        """
        audit: dict[str, np.ndarray] = {}
        for name in LOAD_PATH_AUDIT_ARRAYS:
            try:
                audit[name] = np.asarray(state.array(self.audit_owner, name), dtype=np.float64)
            except KeyError as exc:
                return Refusal(
                    code=RefusalCode.MISSING_STATE_OWNER,
                    reason=(
                        f"{exc}. The load-path audit is not optional: an absent audit is not the "
                        "same as a clean one, and defaulting it to 'sound' is exactly the hazard "
                        "this operator exists to refuse"
                    ),
                    operator_id=self.operator_id,
                    detail={"owner": self.audit_owner, "array": name},
                )

        not_reached = [p.name for p in self.load_path_preconditions]
        for precondition in self.load_path_preconditions:
            not_reached.remove(precondition.name)
            offending = precondition.check(audit, row_names)
            if not offending:
                continue
            booked = audit["booked_force_pn"]
            moved = audit["owner_displacement_um"]
            stalled = booked[moved <= 0.0]
            return Refusal(
                code=RefusalCode.UPSTREAM_DEFECT_OPEN,
                reason=(
                    f"load-path precondition {precondition.name!r} fails on "
                    f"{len(offending)} row(s): {precondition.why}. This operator would return a "
                    "number and it would look like traction; the fix is not in aleph/observe"
                ),
                operator_id=self.operator_id,
                detail={
                    "failed_precondition": precondition.name,
                    "preconditions_not_reached": not_reached,
                    "n_offending_rows": len(offending),
                    "offending_rows": list(offending),
                    "owning_lanes": list(precondition.owning_lanes),
                    "max_booked_force_pn_with_zero_displacement": (
                        float(np.max(stalled)) if stalled.size else 0.0
                    ),
                },
            )
        return None

    # -- the two halves, kept apart ---------------------------------------------------------------

    def raw_observable(
        self, state: AcceptedState, context: ObservationContext
    ) -> RawObservable | Refusal:
        """Read the per-site traction field, **after** the load path has passed its audit.

        The audit runs before any arithmetic, for the reason
        :meth:`~aleph.observe.operator.ApplicabilityDomain.check_owners` gives about owners: an
        operator that computes first and discovers the problem afterwards has already decided what
        to report, and the temptation at that point is to report the part that worked.
        """
        domain = self.applicability

        owner_refusal = domain.check_owners(self.operator_id, state.owners)
        if owner_refusal is not None:
            return owner_refusal

        modulus = context.extras.get("substrate_youngs_modulus_pa")
        if modulus is None:
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    "context.extras carries no 'substrate_youngs_modulus_pa'. It is not defaulted: "
                    "a traction is recovered *through* the substrate's compliance, so a default "
                    "would silently choose the instrument"
                ),
                operator_id=self.operator_id,
                detail={
                    "missing_extra": "substrate_youngs_modulus_pa",
                    "extras": sorted(context.extras),
                },
            )

        bounds_refusal = domain.check(
            self.operator_id, substrate_youngs_modulus_pa=float(modulus)
        )
        if bounds_refusal is not None:
            return bounds_refusal

        try:
            traction = np.asarray(
                state.array(self.adhesion_owner, self.traction_array), dtype=np.float64
            )
            area = np.asarray(state.array(self.adhesion_owner, self.area_array), dtype=np.float64)
        except KeyError as exc:
            return Refusal(
                code=RefusalCode.MISSING_STATE_OWNER,
                reason=str(exc),
                operator_id=self.operator_id,
                detail={"owner": self.adhesion_owner},
            )

        if traction.shape != area.shape or traction.ndim != 1:
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    f"traction {traction.shape} and site area {area.shape} must be the same 1-D "
                    "per-site shape; a traction with no area cannot become a stress"
                ),
                operator_id=self.operator_id,
                detail={"traction_shape": list(traction.shape), "area_shape": list(area.shape)},
            )

        n_sites = int(traction.size)
        raw_names = context.extras.get("load_path_row_names") or ()
        row_names = (
            tuple(str(name) for name in raw_names)
            if len(raw_names) == n_sites
            else tuple(f"row[{i}]" for i in range(n_sites))
        )

        audit_lengths = {
            name: int(np.asarray(state.array(self.audit_owner, name)).size)
            for name in LOAD_PATH_AUDIT_ARRAYS
            if self.audit_owner in state.owners
            and _has_array(state, self.audit_owner, name)
        }
        if audit_lengths and set(audit_lengths.values()) != {n_sites}:
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    f"the load-path audit covers {sorted(set(audit_lengths.values()))} row(s) and "
                    f"the traction field has {n_sites} site(s); an audit that does not cover the "
                    "rows being read is auditing something else"
                ),
                operator_id=self.operator_id,
                detail={"n_sites": n_sites, "audit_lengths": audit_lengths},
            )

        audit_refusal = self.audit_load_path(state, row_names)
        if audit_refusal is not None:
            return audit_refusal

        if not np.all(np.isfinite(traction)) or not np.all(np.isfinite(area)):
            n_bad = int(
                np.count_nonzero(~np.isfinite(traction)) + np.count_nonzero(~np.isfinite(area))
            )
            return Refusal(
                code=RefusalCode.NON_FINITE_INPUT,
                reason=(
                    f"{n_bad} non-finite entr(y/ies) in the traction field. Refused rather than "
                    "masked: a single infinity in a sum of magnitudes destroys the total"
                ),
                operator_id=self.operator_id,
                detail={"n_non_finite": n_bad},
            )

        if np.any(area <= 0.0):
            worst = float(np.min(area))
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    f"site area <= 0 (worst {worst:.6g} um^2); a traction stress would divide by it"
                ),
                operator_id=self.operator_id,
                detail={"min_site_area_um2": worst},
            )

        magnitude = np.abs(traction)
        return RawObservable(
            name="per_adhesion_traction_force",
            operator_id=self.operator_id,
            values=magnitude,
            units="pN",
            n_samples=n_sites,
            manifest_hash=self.manifest.manifest_hash(),
            support={"site_area_um2": area},
            detail={
                "n_sites": n_sites,
                "substrate_youngs_modulus_pa": float(modulus),
                "row_names": list(row_names),
                "load_path_audited": True,
                "preconditions_passed": [p.name for p in self.load_path_preconditions],
            },
        )

    def apparent_quantity(self, raw: RawObservable) -> ApparentQuantity | Refusal:
        """Sum the magnitudes over the footprint, and report the RMS stress alongside.

        ``standard_error`` is ``None``, not zero. The per-site tractions of a TFM inversion are
        correlated through the shared regulariser, so an independent-error propagation would be
        wrong in the flattering direction, and zero would be a claim of exactness. ``None`` is the
        honest statement that this pipeline cannot support an error bar on the total.
        """
        magnitude = np.asarray(raw.values, dtype=np.float64)
        area = np.asarray(raw.support["site_area_um2"], dtype=np.float64)
        total_area = float(np.sum(area))
        if not (total_area > 0.0):
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason="total footprint area is zero; there is no cell to report a traction for",
                operator_id=self.operator_id,
                detail={"total_area_um2": total_area},
            )

        stress = magnitude / area  # pN/um^2, which is Pa exactly
        rms_stress = float(math.sqrt(float(np.sum(stress**2 * area)) / total_area))

        return ApparentQuantity(
            name="apparent_total_traction_force_per_cell",
            operator_id=self.operator_id,
            value=float(np.sum(magnitude)),
            units="pN",
            standard_error=None,
            n_effective=float(magnitude.size),
            manifest_hash=raw.manifest_hash,
            analysis_chain=(
                "audit the load path: every wired row delivers, every booked force moves its owner",
                "take the magnitude of each per-site traction force",
                "sum the magnitudes over the declared footprint",
                "report the area-weighted RMS traction stress alongside, in Pa",
            ),
            identifiability=self.identifiability_directions,
            detail={
                "n_sites": int(magnitude.size),
                "total_footprint_area_um2": total_area,
                "rms_traction_stress_pa": rms_stress,
                "max_site_traction_pn": float(np.max(magnitude)),
                "mean_site_traction_pn": float(np.mean(magnitude)),
                "standard_error_is_none_because": (
                    "per-site tractions from one TFM inversion share a regulariser and are "
                    "therefore correlated; an independent-error sum would understate the total"
                ),
                "external_comparison_is_ratio_only_because": (
                    "the external TFM record tags every traction field source_unit_unknown"
                ),
            },
        )


def _has_array(state: AcceptedState, owner: str, name: str) -> bool:
    """``True`` when the owner carries the named array. Used to size the audit before running it."""
    try:
        state.array(owner, name)
    except KeyError:
        return False
    return True
