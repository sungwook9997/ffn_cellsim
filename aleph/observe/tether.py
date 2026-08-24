"""The AFM membrane-tether forward operator: what an instrument would read off this cell.

Why this operator, and why it is not a socket for a cortical-tension number
--------------------------------------------------------------------------
Pulling a membrane tether with an AFM is the cheapest experiment that interrogates the
**membrane-cortex composite**. The measured plateau force does not report the bilayer's own tension
and it does not report the cortex's tension. It reports the sum:

    sigma_app = sigma_bilayer + W

where ``W`` is the membrane-cortex attachment energy per unit area. ``W`` has units of energy per
area, which in this repository's units is ``pN.um/um^2 = pN/um`` — numerically the same units as a
tension, and that coincidence is precisely why the two get conflated.

So an operator that accepted one lumped "tension" would already have had the answer added in by its
caller. This one takes the two apart, composes them itself, and therefore *can* be asked the
question this lane exists to ask: does Aleph's membrane-cortex coupling, as built, reproduce a real
tether force?

The physics, and where it comes from
------------------------------------
A cylindrical tether of radius ``R`` drawn from a reservoir at apparent tension ``sigma_app`` costs
bending energy and area energy per unit length:

    E/L = pi*kappa/R + 2*pi*R*sigma_app

Minimising over ``R`` gives the two closed forms this module implements:

    R   = sqrt(kappa / (2*sigma_app))
    f   = 2*pi*sqrt(2*kappa*sigma_app)          [and f == E/L at the optimum, which is asserted]

Both are standard (Derenyi/Julicher/Prost; Hochmuth). Note what they imply: ``f^2`` is *linear* in
``sigma_app``, so a tether force is a square-root probe of tension. Two-fold changes in apparent
tension move the force by only 41%, which is the single most important thing to know before
comparing a simulated force to a measured one.

What is NOT modelled, and is refused rather than approximated
-------------------------------------------------------------
* **Finite retraction speed.** A real pull at speed ``v`` reads ``f(v) = f_0 + (viscous term)``.
  This operator computes ``f_0`` only, and refuses above a declared quasi-static speed rather than
  returning the plateau under a label that claims more than it is.
* **The force barrier at tether formation.** The transient peak before the plateau is a separate
  observable with a separate protocol. This operator reports the plateau.
* **Membrane reservoir depletion.** A long pull draws area and raises ``sigma_app``. Modelling it
  needs a reservoir compliance this repository has not sourced, so the operator reports the plateau
  frame by frame from the state it is given and lets the state carry whatever drift is real.

Refusal is a return value here, as everywhere in :mod:`aleph.observe`. In particular a non-positive
apparent tension is a refusal and never a NaN: ``sqrt`` of it would produce one, and a NaN in a
reduction turns a whole array into a complaint about the wrong quantity.
"""

from __future__ import annotations

import math
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
from aleph.observe.stationarity import (
    effective_sample_size,
    integrated_autocorrelation_time,
    standard_error_of_correlated_mean,
)

__all__ = [
    "DEFAULT_MAX_QUASI_STATIC_SPEED_UM_PER_S",
    "MIN_PLATEAU_FRAMES",
    "AfmTetherForceOperator",
    "afm_tether_manifest",
    "apparent_tension_from_force_pn_per_um",
    "apparent_tension_pn_per_um",
    "equilibrium_tether_force_pn",
    "equilibrium_tether_radius_um",
]

_OPERATOR_ID: Final[str] = "aleph.observe.tether.AfmTetherForceOperator/v1"

#: A plateau mean over fewer frames than this is not a measurement. Four is the floor
#: :mod:`aleph.observe.stationarity` needs to estimate a correlation time at all; eight is the floor
#: at which the estimate is worth reporting. Below it the operator refuses rather than reporting a
#: mean whose error bar it cannot support.
MIN_PLATEAU_FRAMES: Final[int] = 8

#: Declared upper bound on retraction speed for the quasi-static plateau [um/s]. Above it the
#: measured force carries a viscous contribution this operator does not model, so it refuses. This
#: is a **declaration about where the operator is entitled to speak**, made before any measurement,
#: not a claim about where membrane friction becomes physically important.
DEFAULT_MAX_QUASI_STATIC_SPEED_UM_PER_S: Final[float] = 1.0


# -- closed forms ------------------------------------------------------------------------------


def apparent_tension_pn_per_um(
    bilayer_tension_pn_per_um: float, membrane_cortex_attachment_pn_per_um: float
) -> float:
    """``sigma_app = sigma_bilayer + W`` [pN/um].

    Kept as a named function with two arguments rather than inlined, because the whole design claim
    of this module is that these are two quantities and not one.

    Args:
        bilayer_tension_pn_per_um: In-plane tension of the lipid bilayer alone [pN/um].
        membrane_cortex_attachment_pn_per_um: Membrane-cortex attachment energy per unit area
            [pN.um/um^2 = pN/um]. Zero for a membrane that is not attached to a cortex.

    Returns:
        The apparent tension a tether pull actually probes.
    """
    return float(bilayer_tension_pn_per_um) + float(membrane_cortex_attachment_pn_per_um)


def equilibrium_tether_force_pn(
    bending_rigidity_pn_um: float, apparent_tension_pn_per_um_value: float
) -> float:
    """``f = 2*pi*sqrt(2*kappa*sigma_app)`` [pN].

    Args:
        bending_rigidity_pn_um: ``kappa`` [pN.um].
        apparent_tension_pn_per_um_value: ``sigma_app`` [pN/um]. Must be positive; this function is
            the arithmetic and the caller is responsible for having refused a non-positive value
            (see :meth:`AfmTetherForceOperator.raw_observable`).

    Returns:
        The equilibrium plateau force [pN].
    """
    return 2.0 * math.pi * math.sqrt(
        2.0 * float(bending_rigidity_pn_um) * float(apparent_tension_pn_per_um_value)
    )


def equilibrium_tether_radius_um(
    bending_rigidity_pn_um: float, apparent_tension_pn_per_um_value: float
) -> float:
    """``R = sqrt(kappa/(2*sigma_app))`` [um].

    Reported alongside the force because it is the operator's own falsifiability check: a recovered
    radius outside roughly 10-50 nm says the ``(kappa, sigma_app)`` pair is not a cell membrane,
    whatever the force came out as.
    """
    return math.sqrt(
        float(bending_rigidity_pn_um) / (2.0 * float(apparent_tension_pn_per_um_value))
    )


def apparent_tension_from_force_pn_per_um(
    force_pn: float, bending_rigidity_pn_um: float
) -> float:
    """``sigma_app = f^2 / (8*pi^2*kappa)`` [pN/um] — the inversion, with its assumption named.

    This is the step a published tether-force paper takes to report an "apparent membrane tension",
    and it is **only** valid at constant ``kappa``. Any comparison of two tensions obtained this way
    is a comparison of two forces squared, with a bending rigidity assumed equal on both sides.

    It is exposed as a separate function, and not folded into the operator, so that a caller cannot
    take the step without writing its name.
    """
    return (float(force_pn) ** 2) / (
        8.0 * math.pi * math.pi * float(bending_rigidity_pn_um)
    )


# -- the manifest ------------------------------------------------------------------------------


def afm_tether_manifest(
    *,
    temperature: ManifestValue,
    medium: ManifestValue,
    probe_geometry: ManifestValue,
    control_mode: ManifestValue,
    rate: ManifestValue,
    sampling_cadence: ManifestValue,
    fit_window: ManifestValue,
    resolution_limit: ManifestValue,
    specimen_identity: ManifestValue,
    batch: ManifestValue,
    citation_audit: ManifestValue,
    applicability_statement: ManifestValue,
    evidence_role: ManifestValue,
    operator_version: ManifestValue,
    probe_location: ManifestValue = Declaration.NOT_RECORDED,
    probe_direction: ManifestValue = Declaration.NOT_RECORDED,
    adhesion_state: ManifestValue = Declaration.NOT_RECORDED,
    dwell: ManifestValue = Declaration.NOT_RECORDED,
    loading_history: ManifestValue = Declaration.NOT_RECORDED,
    preconditioning: ManifestValue = Declaration.NOT_RECORDED,
    calibration: ManifestValue = Declaration.NOT_RECORDED,
    osmotic_state: ManifestValue = Declaration.NOT_RECORDED,
    substrate: ManifestValue = Declaration.NOT_RECORDED,
    confinement: ManifestValue = Declaration.NOT_RECORDED,
    pharmacology: ManifestValue = Declaration.NOT_RECORDED,
) -> ProtocolManifest:
    """The manifest for an AFM membrane-tether pull. All thirty-four fields, written out.

    Read this next to :func:`~aleph.observe.manifest.passive_fluctuation_manifest` — the two are
    nearly complementary, and that is the argument for keeping the manifest a type rather than a
    convention. The passive protocol asserts ``NOT_APPLICABLE`` across the whole
    geometry-and-control family because nothing touches the membrane. Here something does, and
    asserting the same thing would be a false record rather than a terse one.

    Only **two** fields are asserted ``NOT_APPLICABLE`` here, and each is argued below. Everything
    else is either required of the caller or defaults to :attr:`Declaration.NOT_RECORDED` — the
    confession, never the assertion. ``pharmacology`` in particular defaults to the confession,
    which matters because the external record this operator is compared against is a drug
    perturbation.

    Args:
        temperature: With units, e.g. ``"310.0 K"``.
        medium: Bathing solution. Sets the viscosity, hence the speed at which the quasi-static
            assumption fails.
        probe_geometry: Cantilever and bead. There *is* a probe.
        control_mode: Constant-speed retraction, force clamp, or otherwise. The plateau force means
            different things under different control.
        rate: Retraction speed, with units.
        sampling_cadence: Frame interval, with units.
        fit_window: The extension window the plateau is claimed over, declared before the fit.
        resolution_limit: The force noise floor, with units.
        specimen_identity: What was measured.
        batch: Preparation batch.
        citation_audit: Where the parameters came from, or that they came from nowhere.
        applicability_statement: The domain over which the result is claimed.
        evidence_role: What the number is used for.
        operator_version: Version of the code that produced it.
        probe_location: Where on the cell. Defaults to ``NOT_RECORDED``.
        probe_direction: Pull direction. Defaults to ``NOT_RECORDED``.
        adhesion_state: Bead-membrane attachment chemistry. Defaults to ``NOT_RECORDED``.
        dwell: Contact dwell before retraction. Defaults to ``NOT_RECORDED``.
        loading_history: Prior pulls on the same cell. Defaults to ``NOT_RECORDED``.
        preconditioning: Defaults to ``NOT_RECORDED``.
        calibration: Cantilever spring constant and optical lever sensitivity. Defaults to
            ``NOT_RECORDED``. A forward-operator run has no transducer and should say so in words
            here rather than take the ``NOT_APPLICABLE`` assertion, because the record it will be
            compared against does have one.
        osmotic_state: Defaults to ``NOT_RECORDED``.
        substrate: Defaults to ``NOT_RECORDED``.
        confinement: Defaults to ``NOT_RECORDED``.
        pharmacology: Defaults to ``NOT_RECORDED``.

    Returns:
        The frozen :class:`~aleph.observe.manifest.ProtocolManifest`.
    """
    na = Declaration.NOT_APPLICABLE
    return ProtocolManifest(
        # Geometry and control: every field applies. A probe exists, it is somewhere, it pulls in
        # some direction, under some control, through some contact, against some adhesion.
        probe_geometry=probe_geometry,
        probe_location=probe_location,
        probe_direction=probe_direction,
        control_mode=control_mode,
        contact_model=(
            "cylindrical membrane tether at its energy minimum: "
            "E/L = pi*kappa/R + 2*pi*R*sigma_app, "
            "giving R = sqrt(kappa/2*sigma_app) and f = 2*pi*sqrt(2*kappa*sigma_app). Fixed by "
            "this "
            "operator; a caller needing a different one is running a different protocol and should "
            "say so with amended(), which changes the hash"
        ),
        adhesion_state=adhesion_state,
        # Time and amplitude.
        rate=rate,
        # NOT_APPLICABLE, argued: a constant-speed retraction imposes no oscillation. There is no
        # drive frequency to record. This is an assertion about the protocol and it is correct.
        frequency=na,
        dwell=dwell,
        loading_history=loading_history,
        # NOT_APPLICABLE, argued: a tether at its plateau is a fluid structure held at constant
        # force and constant radius. It has an extension, which `fit_window` records, but it has no
        # strain amplitude — there is no reference configuration being strained.
        max_strain=na,
        preconditioning=preconditioning,
        sampling_cadence=sampling_cadence,
        # Environment: all of it applies.
        temperature=temperature,
        medium=medium,
        osmotic_state=osmotic_state,
        substrate=substrate,
        confinement=confinement,
        pharmacology=pharmacology,
        # Observation: fixed by this operator.
        raw_output_schema=(
            "per-frame equilibrium plateau force, with the tether radius on the record"
        ),
        raw_output_units="force in pN, tether radius in um, tension in pN/um, kappa in pN.um",
        calibration=calibration,
        # Analysis. Neither an image nor an ill-posed inverse problem is involved: the plateau force
        # follows in closed form from (kappa, sigma_app), so there is nothing to segment and nothing
        # to regularise. Both are assertions and both are correct for THIS analysis; a pipeline that
        # fitted a plateau out of a noisy trace would have to amend them.
        segmentation=na,
        inverse_model=(
            "closed-form cylindrical tether energy minimum; no fit. The reverse step "
            "sigma_app = f^2/(8*pi^2*kappa) is a SEPARATE function and assumes kappa constant"
        ),
        regularization=na,
        fit_window=fit_window,
        noise_model=(
            "frame-to-frame scatter of the plateau force, with the sample count corrected to "
            "n_eff = N/(2*tau_samples) because consecutive frames of a stepped simulation are "
            "correlated; sigma/sqrt(N) is never reported"
        ),
        resolution_limit=resolution_limit,
        # Evidence: all of it applies, always.
        specimen_identity=specimen_identity,
        batch=batch,
        citation_audit=citation_audit,
        applicability_statement=applicability_statement,
        evidence_role=evidence_role,
        operator_version=operator_version,
    )


# -- the operator ------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AfmTetherForceOperator:
    """Forward model of an AFM membrane-tether pull. Satisfies ``ObservationOperator``.

    Reads three per-frame series and composes them. The state contract is deliberately in *named
    arrays* rather than in an imported material object, because :mod:`aleph.observe` does not import
    :mod:`aleph.vertical` and must not start: an observer that holds a reference to the owner it
    observes will eventually be tempted to ask it for the answer.

    Args:
        manifest: The protocol in force. Build it with :func:`afm_tether_manifest`.
        membrane_owner: State owner holding the bilayer series. Default ``"membrane"``.
        attachment_owner: State owner holding the membrane-cortex attachment energy series. Default
            ``"cortex"``. A cell with no cortex should supply a zero series under this name rather
            than omit the owner, so that "not attached" is a measured zero and not a missing input.
        tension_array: Default ``"tension_pn_per_um"``, shape ``(T,)``.
        bending_rigidity_array: Default ``"bending_rigidity_pn_um"``, shape ``(T,)``.
        attachment_array: Default ``"membrane_attachment_energy_pn_per_um"``, shape ``(T,)``.
        min_temperature_k: Lower applicability bound.
        max_temperature_k: Upper applicability bound.
        max_quasi_static_speed_um_per_s: Declared upper bound on retraction speed.
    """

    manifest: ProtocolManifest
    membrane_owner: str = "membrane"
    attachment_owner: str = "cortex"
    tension_array: str = "tension_pn_per_um"
    bending_rigidity_array: str = "bending_rigidity_pn_um"
    attachment_array: str = "membrane_attachment_energy_pn_per_um"
    min_temperature_k: float = 250.0
    max_temperature_k: float = 330.0
    max_quasi_static_speed_um_per_s: float = DEFAULT_MAX_QUASI_STATIC_SPEED_UM_PER_S

    # -- declarations ---------------------------------------------------------------------------

    @property
    def operator_id(self) -> str:
        return _OPERATOR_ID

    @property
    def required_state_owners(self) -> tuple[str, ...]:
        return (self.membrane_owner, self.attachment_owner)

    @property
    def identifiability_directions(self) -> tuple[str, ...]:
        """What one tether force can constrain — and the degeneracy that is built into it.

        The honest statement is short and unflattering. A tether force is a square-root probe of a
        *sum*. It cannot separate the two terms of that sum, and it responds to either of them with
        half the leverage the number's precision suggests.
        """
        return (
            "sigma_app = sigma_bilayer + W: constrained, as f^2/(8*pi^2*kappa)",
            "kappa: constrained only jointly with sigma_app, through the product kappa*sigma_app",
            "NOT constrained: the split of sigma_app into sigma_bilayer + W. One force cannot "
            "separate a sum; a second observable at a different kappa, or an independent bilayer "
            "tension, is required",
            "NOT constrained: kappa and sigma_app separately under (kappa, sigma_app) -> "
            "(kappa*g, sigma_app/g), which leaves the force exactly invariant and moves the tether "
            "radius by g — so the RADIUS is the observable that breaks this one",
        )

    @property
    def applicability(self) -> ApplicabilityDomain:
        return ApplicabilityDomain(
            description=(
                "quasi-static plateau of a single cylindrical membrane tether pulled from a cell "
                "whose membrane-cortex attachment energy is on the record; no reservoir depletion "
                "modelled, no tether-formation barrier reported"
            ),
            bounds={
                "temperature_k": (self.min_temperature_k, self.max_temperature_k),
                "retraction_speed_um_per_s": (0.0, self.max_quasi_static_speed_um_per_s),
                "sample_interval_s": (0.0, math.inf),
            },
            required_state_owners=(self.membrane_owner, self.attachment_owner),
        )

    # -- the forward map, exposed so a control can call it without the wrapper --------------------

    def apparent_tension_series(self, state: AcceptedState) -> np.ndarray:
        """Per-frame ``sigma_app = sigma_bilayer + W`` [pN/um]. Refusals must already have run.

        This is the **only** place the two terms are composed. A second copy of this sum elsewhere
        in the module is how the composition gets silently dropped from one code path while the
        other keeps a control green — which is exactly what a mutation check on an earlier draft of
        this file demonstrated.
        """
        sigma = np.asarray(state.array(self.membrane_owner, self.tension_array), dtype=np.float64)
        attach = np.asarray(
            state.array(self.attachment_owner, self.attachment_array), dtype=np.float64
        )
        return sigma + attach

    def force_series(self, state: AcceptedState) -> np.ndarray:
        """Per-frame equilibrium plateau force [pN]. Refusals must already have run.

        The **single** implementation of the forward map. :meth:`raw_observable` calls this rather
        than repeating the arithmetic, and a control pins the two together with ``np.array_equal``.
        """
        kappa = np.asarray(
            state.array(self.membrane_owner, self.bending_rigidity_array), dtype=np.float64
        )
        return 2.0 * np.pi * np.sqrt(2.0 * kappa * self.apparent_tension_series(state))

    def radius_series(self, state: AcceptedState) -> np.ndarray:
        """Per-frame equilibrium tether radius [um]. Refusals must already have run."""
        kappa = np.asarray(
            state.array(self.membrane_owner, self.bending_rigidity_array), dtype=np.float64
        )
        return np.sqrt(kappa / (2.0 * self.apparent_tension_series(state)))

    # -- the two halves, kept apart --------------------------------------------------------------

    def raw_observable(
        self, state: AcceptedState, context: ObservationContext
    ) -> RawObservable | Refusal:
        """Produce the per-frame plateau force the instrument would read, with no interpretation.

        Refuses, in this order and before any arithmetic that could produce a number: missing owner,
        missing or non-quasi-static retraction speed, out-of-domain context, missing array, ragged
        shapes, too few frames, non-finite samples, and finally a non-positive apparent tension.

        The order is not cosmetic. The finiteness check must precede the positivity check, or a NaN
        would be reported as an out-of-domain tension and the caller would go looking for a physics
        problem that is actually a state problem.
        """
        domain = self.applicability

        owner_refusal = domain.check_owners(self.operator_id, state.owners)
        if owner_refusal is not None:
            return owner_refusal

        speed = context.extras.get("retraction_speed_um_per_s")
        if speed is None:
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    "context.extras carries no 'retraction_speed_um_per_s'. It is not defaulted: "
                    "supplying one would assert on the caller's behalf that the pull was "
                    "quasi-static, which is the single assumption this operator's plateau rests on"
                ),
                operator_id=self.operator_id,
                detail={
                    "missing_extra": "retraction_speed_um_per_s",
                    "extras": sorted(context.extras),
                },
            )

        bounds_refusal = domain.check(
            self.operator_id,
            temperature_k=context.temperature_k,
            retraction_speed_um_per_s=float(speed),
            sample_interval_s=context.sample_interval_s,
        )
        if bounds_refusal is not None:
            return bounds_refusal

        if not (context.sample_interval_s > 0.0):
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    f"sample_interval_s = {context.sample_interval_s!r}; the correlation "
                    "correction "
                    "on the plateau mean cannot be expressed without a positive sampling interval"
                ),
                operator_id=self.operator_id,
                detail={"sample_interval_s": context.sample_interval_s},
            )

        wanted = (
            (self.membrane_owner, self.tension_array),
            (self.membrane_owner, self.bending_rigidity_array),
            (self.attachment_owner, self.attachment_array),
        )
        series: list[np.ndarray] = []
        for owner, name in wanted:
            try:
                series.append(np.asarray(state.array(owner, name), dtype=np.float64))
            except KeyError as exc:
                return Refusal(
                    code=RefusalCode.MISSING_STATE_OWNER,
                    reason=str(exc),
                    operator_id=self.operator_id,
                    detail={"owner": owner, "array": name},
                )

        shapes = [tuple(a.shape) for a in series]
        if any(a.ndim != 1 for a in series):
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    f"every input must be a 1-D per-frame series; got shapes {shapes}. A scalar is "
                    "one frame and carries no scatter to put an error bar on"
                ),
                operator_id=self.operator_id,
                detail={"shapes": [list(s) for s in shapes]},
            )
        if len(set(shapes)) != 1:
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    f"the three per-frame series have different lengths {shapes}; they would have "
                    "to be aligned by an assumption this operator is not entitled to make"
                ),
                operator_id=self.operator_id,
                detail={"shapes": [list(s) for s in shapes]},
            )

        n_frames = int(shapes[0][0])
        if n_frames < MIN_PLATEAU_FRAMES:
            return Refusal(
                code=RefusalCode.INSUFFICIENT_SAMPLES,
                reason=(
                    f"{n_frames} frame(s); a plateau force is a mean over a window and needs at "
                    f"least {MIN_PLATEAU_FRAMES} to carry an error bar that means anything"
                ),
                operator_id=self.operator_id,
                detail={"n_frames": n_frames, "minimum": MIN_PLATEAU_FRAMES},
            )

        for (owner, name), data in zip(wanted, series, strict=True):
            if not np.all(np.isfinite(data)):
                n_bad = int(np.count_nonzero(~np.isfinite(data)))
                return Refusal(
                    code=RefusalCode.NON_FINITE_INPUT,
                    reason=(
                        f"{n_bad} non-finite sample(s) in {owner}.{name}. Refused rather than "
                        "masked: dropping frames changes the sampling interval the correlation "
                        "time is expressed in"
                    ),
                    operator_id=self.operator_id,
                    detail={"owner": owner, "array": name, "n_non_finite": n_bad},
                )

        tension, kappa, attach = series
        sigma_app = self.apparent_tension_series(state)

        if not np.all(sigma_app > 0.0):
            worst = float(np.min(sigma_app))
            n_bad = int(np.count_nonzero(sigma_app <= 0.0))
            return Refusal(
                code=RefusalCode.OUT_OF_APPLICABILITY_DOMAIN,
                reason=(
                    f"{n_bad} frame(s) have apparent tension sigma_app = sigma_bilayer + W <= 0 "
                    f"(worst {worst:.6g} pN/um). There is no tether: with no positive apparent "
                    "tension nothing sets a finite radius, and the closed form would return a NaN"
                ),
                operator_id=self.operator_id,
                detail={
                    "n_non_positive": n_bad,
                    "min_apparent_tension_pn_per_um": worst,
                    "quantity": "apparent_tension_pn_per_um",
                },
            )

        if not np.all(kappa > 0.0):
            worst = float(np.min(kappa))
            return Refusal(
                code=RefusalCode.OUT_OF_APPLICABILITY_DOMAIN,
                reason=(
                    f"bending rigidity kappa <= 0 (worst {worst:.6g} pN.um). A membrane with no "
                    "bending cost has no tether radius; the minimisation has no interior solution"
                ),
                operator_id=self.operator_id,
                detail={"min_bending_rigidity_pn_um": worst, "quantity": "bending_rigidity_pn_um"},
            )

        force = self.force_series(state)
        radius = self.radius_series(state)

        return RawObservable(
            name="membrane_tether_plateau_force",
            operator_id=self.operator_id,
            values=force,
            units="pN",
            n_samples=n_frames,
            manifest_hash=self.manifest.manifest_hash(),
            support={
                "time_s": np.arange(n_frames, dtype=np.float64) * float(context.sample_interval_s),
                "tether_radius_um": radius,
                "apparent_tension_pn_per_um": sigma_app,
            },
            detail={
                "sample_interval_s": float(context.sample_interval_s),
                "retraction_speed_um_per_s": float(speed),
                "mean_bending_rigidity_pn_um": float(np.mean(kappa)),
                "mean_bilayer_tension_pn_per_um": float(np.mean(tension)),
                "mean_attachment_energy_pn_per_um": float(np.mean(attach)),
                "mean_tether_radius_um": float(np.mean(radius)),
            },
        )

    def apparent_quantity(self, raw: RawObservable) -> ApparentQuantity | Refusal:
        """Reduce the plateau to one number the way a real analysis would, and say how.

        The error bar is the correlated one. On a stepped simulation consecutive frames are
        correlated by construction, so ``sigma/sqrt(N)`` would understate the error in the
        flattering direction — see :mod:`aleph.observe.stationarity`.
        """
        force = np.asarray(raw.values, dtype=np.float64)
        n_frames = int(force.size)
        dt = float(raw.detail.get("sample_interval_s", 0.0))
        if not (dt > 0.0):
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    "the raw observable carries no positive sample_interval_s, so no correlation "
                    "time can be formed and no honest error bar can be reported"
                ),
                operator_id=self.operator_id,
                detail={"sample_interval_s": dt},
            )

        tau = integrated_autocorrelation_time(force, dt)
        std = float(np.std(force))
        n_eff = effective_sample_size(n_frames, tau.tau_samples)
        # A plateau with no scatter has no error bar. `sigma*sqrt(2*tau/N)` returns 0.0 here, and a
        # zero standard error reads as infinite precision rather than as "nothing varied" -- which
        # is what `ApparentQuantity.standard_error` means by "never zero as a stand-in". Today's
        # runtime reaches this branch by construction: `gamma_mem` is the constant Raucher-Sheetz
        # plateau (`compartments.py:579`, K_A upturn default-OFF), `kappa` is a resolved constant,
        # and `GAMMA_MCA_PN_UM` is a literal -- so the composed force is an echo of three inputs.
        # Reported as None with the reason on the record, not refused: the VALUE is still the
        # forward map the state implies, and refusing would hide that the state carried no dynamics.
        constant_plateau = std == 0.0
        sem = (
            None
            if constant_plateau
            else standard_error_of_correlated_mean(std, n_frames, tau.tau_samples)
        )

        radius = np.asarray(raw.support["tether_radius_um"], dtype=np.float64)

        return ApparentQuantity(
            name="apparent_membrane_tether_force",
            operator_id=self.operator_id,
            value=float(np.mean(force)),
            units="pN",
            standard_error=sem,
            n_effective=n_eff,
            manifest_hash=raw.manifest_hash,
            analysis_chain=(
                "compose sigma_app = sigma_bilayer + W per frame",
                "closed-form cylindrical tether plateau f = 2*pi*sqrt(2*kappa*sigma_app) per frame",
                "mean over the declared plateau window",
                "standard error of the correlated mean, sigma*sqrt(2*tau_samples/N), with "
                "tau_samples from a Sokal-windowed integrated autocorrelation time",
            ),
            identifiability=self.identifiability_directions,
            detail={
                "n_frames": n_frames,
                "constant_plateau": constant_plateau,
                "standard_error_note": (
                    "zero frame-to-frame scatter: the plateau force is a restatement of constant "
                    "(kappa, sigma_bilayer, W) inputs, not a measurement of a varying state"
                    if constant_plateau
                    else "standard error of the correlated mean"
                ),
                "tau_samples": float(tau.tau_samples),
                "tau_int_s": float(tau.tau_int_s),
                "tau_window_closed": bool(tau.window_closed),
                "tau_note": tau.note,
                "plateau_std_pn": std,
                "mean_tether_radius_um": float(np.mean(radius)),
                "mean_apparent_tension_pn_per_um": float(
                    np.mean(np.asarray(raw.support["apparent_tension_pn_per_um"], dtype=np.float64))
                ),
                "mean_bending_rigidity_pn_um": raw.detail.get("mean_bending_rigidity_pn_um"),
                "mean_bilayer_tension_pn_per_um": raw.detail.get("mean_bilayer_tension_pn_per_um"),
                "mean_attachment_energy_pn_per_um": raw.detail.get(
                    "mean_attachment_energy_pn_per_um"
                ),
            },
        )
