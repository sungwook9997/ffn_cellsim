"""Open-loop ECM dynamics preflight for prescribed inputs.

This module covers the separated-dynamics-mode ECM updates that take
caller-supplied prescribed input fields and integrate them into
:class:`acs.v2.ecm_substrate.ECMSubstrateState` **without closing the
loop back to cell motion, FA state, or any biology-derived input**.
Closed-loop ECM activation is gated behind the six-item gate in
``docs/v2/v2_phase1_forward_roadmap.md`` "Closed-Loop ECM Gate".

Four preflight functions live here, one per ECM field family:

- ``accumulate_prescribed_traction(ecm, traction_density_nN_per_um2,
  dt_s)`` — ECM-OL-1 (commit ``29a360f``): integrates a prescribed
  non-negative scalar traction-density field into the cumulative
  ``accumulated_traction_nNs_per_um2`` history. Monotone
  non-decreasing by construction (no relaxation in preflight).
- ``apply_prescribed_stiffness_rate(ecm, stiffness_rate_kpa_per_s,
  dt_s)`` — 6.4-open-A (commit ``7dc1767``): applies a prescribed
  signed stiffness rate field into ``stiffness_kpa``. Signed
  rates are allowed (open-loop preflight permits stiffening *and*
  softening), but the **post-state** ``stiffness_kpa`` must remain
  non-negative per the ECM schema; an update that would push any
  cell below zero raises ``stiffness_negative_post_update`` instead
  of clamping. No upper cap is invented.
- ``apply_prescribed_density_rate(ecm, ligand_density_rate_per_s,
  fiber_density_rate_per_s, dt_s)`` — 6.4-open-B (commit ``e4a71c9``):
  applies prescribed signed density-rate fields into
  ``ligand_density`` and ``fiber_density`` simultaneously. Both
  density fields are dimensionless and bounded to ``[0, 1]`` by
  the ECM schema; signed rates are allowed (positive thickens,
  negative thins), but post-state values must stay in ``[0, 1]``.
  Lower-bound violations raise
  ``ligand_density_negative_post_update`` or
  ``fiber_density_negative_post_update``; upper-bound violations
  raise ``ligand_density_exceeds_one_post_update`` or
  ``fiber_density_exceeds_one_post_update``. No clamping. The two
  rate fields are processed in strict cross-isolation: the ligand
  update only consumes ``ligand_density_rate_per_s`` and the
  fiber update only consumes ``fiber_density_rate_per_s``.
- ``apply_prescribed_orientation_rate(ecm, orientation_rate_per_s,
  dt_s)`` — 6.4-open-C (commit ``8b2c1ab``): applies a prescribed signed
  ``(nx, ny, 2, 2)`` orientation-tensor rate field into
  ``orientation_tensor``. Signed rates are allowed (positive grows
  the local component, negative shrinks it), but the rate field
  must already be component-wise symmetric to within the ECM
  schema's ``_ORIENTATION_SYMMETRY_TOL`` (no auto-symmetrization),
  and the post-state must stay component-wise within
  ``[-_ORIENTATION_BOUND, +_ORIENTATION_BOUND]`` (no clamping).
  Failure kinds: ``orientation_rate_shape``,
  ``orientation_rate_non_finite``,
  ``orientation_rate_asymmetric``,
  ``orientation_negative_bound_post_update``,
  ``orientation_exceeds_bound_post_update``. The ECM schema
  constants are imported from :mod:`acs.v2.ecm_substrate` so a
  future tightening (e.g. a smaller symmetry tolerance) propagates
  here without a separate magic number.

Sanity Gate (open-loop preflight, applies to every function in this
module):

    1. Dimensional analysis. Each function performs exactly one Rule
       10 unit-chain reduction and stores the result in the matching
       ``ECMSubstrateState`` field:
         - ``accumulate_prescribed_traction``:
           ``[nN/μm²] · [s] = [nN·s/μm²]`` →
           ``accumulated_traction_nNs_per_um2``.
         - ``apply_prescribed_stiffness_rate``:
           ``[kPa/s] · [s] = [kPa]`` → ``stiffness_kpa``.
         - ``apply_prescribed_density_rate``:
           ``[1/s] · [s] = [dimensionless density]`` →
           ``ligand_density`` and ``fiber_density``.
         - ``apply_prescribed_orientation_rate``:
           ``[1/s] · [s] = [dimensionless orientation tensor]`` →
           ``orientation_tensor``.
       No kPa↔nN/μm² conversion, no force-per-volume comparison, no
       CFL/Re/Ca/De/Pe/Ma. Closed-loop ECM gate item 6 (the kPa↔
       nN/μm² unit-chain proof) is the closed-loop gate's
       responsibility, not this preflight's.
    2. Boundary cases. Each function rejects: shape mismatch with
       the ECM grid (``traction_shape``, ``stiffness_rate_shape``,
       ``ligand_density_rate_shape``, ``fiber_density_rate_shape``,
       ``orientation_rate_shape``); non-finite or bool inputs
       (``*_non_finite`` and ``dt_invalid``); ``dt_s < 0`` or
       non-finite (``dt_invalid``); per-field schema violations on
       the *post-state*:
         - ``accumulate_prescribed_traction``: negative input
           rejected (``traction_negative``); the cumulative field
           is monotone non-decreasing by construction.
         - ``apply_prescribed_stiffness_rate``: signed input
           accepted; post-state ``stiffness_kpa < 0`` rejected
           (``stiffness_negative_post_update``).
         - ``apply_prescribed_density_rate``: signed input
           accepted per density field; post-state must stay in
           ``[0, 1]`` per field, with split lower/upper failure
           kinds (``ligand_density_negative_post_update``,
           ``ligand_density_exceeds_one_post_update``,
           ``fiber_density_negative_post_update``,
           ``fiber_density_exceeds_one_post_update``).
         - ``apply_prescribed_orientation_rate``: input must be
           component-wise symmetric within
           ``_ORIENTATION_SYMMETRY_TOL``
           (``orientation_rate_asymmetric``); post-state must stay
           component-wise in ``[-_ORIENTATION_BOUND,
           +_ORIENTATION_BOUND]`` with split lower/upper failure
           kinds (``orientation_negative_bound_post_update``,
           ``orientation_exceeds_bound_post_update``).
       ``dt_s == 0`` and zero input fields are exact no-ops
       returning a fresh state. No clamping in any function.
    3. Conservation / provenance invariants. Input ``ecm`` is never
       mutated. Every non-target ECM field is copied verbatim; the
       returned state runs through ``ECMSubstrateState.validate()``
       so any breakage of schema invariants surfaces here, not at
       downstream consumers. ``accumulate_prescribed_traction``
       additionally guarantees monotonic non-decrease of the
       cumulative field by rejecting negative input.
    4. Numerical sanity. float64 throughout. There is no explicit
       stability bound because every preflight function in this
       module is a linear open-loop recorder, not a feedback update.
       Grid resolution is inherited from ``ECMSubstrateState``.
    5. Sign / sense check.
         - traction: positive prescribed traction adds to the
           cumulative exposure; negative input rejected at validation
           because preflight has no relaxation mechanism.
         - stiffness rate: positive rate stiffens, negative rate
           softens, zero rate / zero dt no-ops. The function never
           infers stiffening from cell forces; the rate field is
           caller-supplied (closed-loop is the only place where
           biology-derived rates are wired in).
         - density rates: positive rates thicken (toward 1),
           negative rates thin (toward 0). ligand and fiber
           density rates are processed in strict cross-isolation —
           ligand update never reads fiber rate and vice versa.
         - orientation rate: positive components grow the local
           orientation tensor entry, negative components shrink it.
           No symmetrization or normalization is performed — the
           caller is responsible for supplying a symmetric rate
           tensor. The schema's symmetry and bound checks are
           reused (not redefined here) so the open-loop preflight
           cannot drift away from the schema contract.
    6. Measurement-protocol consistency (Hard Rule 11). The grid
       layout (origin, spacing, shape) is preserved verbatim so the
       same coordinates that downstream visualization or comparison
       reads still resolve to the same physical cell. ``frame_dump``
       round-trip already re-validates these invariants per the ECM
       schema; this module ensures it never produces an invalid
       output.

Magic-Number Block: no tunable constants in this module. The
schema constants ``_ORIENTATION_SYMMETRY_TOL`` (numerical tie-break
for symmetry) and ``_ORIENTATION_BOUND`` (dimensionless cap on the
orientation tensor's component magnitudes) are imported from
:mod:`acs.v2.ecm_substrate` and reused by
``apply_prescribed_orientation_rate`` rather than duplicated. They
remain documented as numerical tie-breaks at their definition site,
not physics tunables.
"""

from __future__ import annotations

import numpy as np

from acs.v2.ecm_substrate import (
    _ORIENTATION_BOUND,
    _ORIENTATION_SYMMETRY_TOL,
    ECMSubstrateState,
)


class ECMOpenLoopError(ValueError):
    """Validation error carrying a machine-readable failure kind."""

    def __init__(self, failure_kind: str, message: str) -> None:
        self.failure_kind = failure_kind
        super().__init__(message)


def accumulate_prescribed_traction(
    ecm: ECMSubstrateState,
    traction_density_nN_per_um2,
    dt_s: float,
) -> ECMSubstrateState:
    """Accumulate prescribed scalar traction density into ECM history.

    Args:
        ecm: Existing ECM substrate state. It is validated and never mutated.
        traction_density_nN_per_um2: Array with shape ``ecm.grid_shape`` and
            units [nN/um^2]. Values must be finite and non-negative.
        dt_s: Non-negative scalar timestep [s]. ``0`` is allowed and produces
            an exact no-op on accumulated history.

    Returns:
        A new :class:`ECMSubstrateState` with identical stiffness, ligand,
        fiber, orientation, source, origin, and spacing, and with
        ``accumulated_traction_nNs_per_um2`` incremented by
        ``traction_density_nN_per_um2 * dt_s``.

    Raises:
        ECMOpenLoopError: If the timestep or traction field violates the
            open-loop contract.
        ValueError: If ``ecm`` itself fails its schema validation.
    """

    ecm.validate()
    dt = _validate_dt(dt_s)
    traction = _validate_traction_density(
        traction_density_nN_per_um2, expected_shape=ecm.grid_shape
    )

    accumulated = (
        np.asarray(ecm.accumulated_traction_nNs_per_um2, dtype=np.float64)
        + traction * dt
    )
    next_ecm = ECMSubstrateState(
        origin_um_xy=tuple(ecm.origin_um_xy),
        spacing_um=float(ecm.spacing_um),
        stiffness_kpa=np.array(ecm.stiffness_kpa, dtype=np.float64, copy=True),
        ligand_density=np.array(ecm.ligand_density, dtype=np.float64, copy=True),
        fiber_density=np.array(ecm.fiber_density, dtype=np.float64, copy=True),
        orientation_tensor=np.array(ecm.orientation_tensor, dtype=np.float64, copy=True),
        accumulated_traction_nNs_per_um2=accumulated,
        source=ecm.source,
    )
    next_ecm.validate()
    return next_ecm


def _validate_dt(dt_s: float) -> float:
    if isinstance(dt_s, bool):
        raise ECMOpenLoopError("dt_invalid", "dt_s must be a finite non-negative scalar")
    # Reject numpy bool scalars and 0-d bool arrays. ``float(np.True_) -> 1.0``
    # and ``float(np.array(True)) -> 1.0`` would otherwise sneak past the
    # python-bool guard above and silently coerce to dt_s == 1.0.
    if isinstance(dt_s, np.bool_):
        raise ECMOpenLoopError("dt_invalid", "dt_s must be a finite non-negative scalar")
    if isinstance(dt_s, np.ndarray):
        if dt_s.dtype == bool or dt_s.shape != ():
            raise ECMOpenLoopError(
                "dt_invalid", "dt_s must be a finite non-negative scalar"
            )
    try:
        dt = float(dt_s)
    except (TypeError, ValueError):
        raise ECMOpenLoopError(
            "dt_invalid", "dt_s must be a finite non-negative scalar"
        ) from None
    if not np.isfinite(dt) or dt < 0.0:
        raise ECMOpenLoopError("dt_invalid", "dt_s must be finite and non-negative")
    return dt


def _validate_traction_density(value, *, expected_shape: tuple[int, int]) -> np.ndarray:
    traction = np.asarray(value, dtype=np.float64)
    if traction.shape != expected_shape:
        raise ECMOpenLoopError(
            "traction_shape",
            f"traction_density_nN_per_um2 must have shape {expected_shape!r}, "
            f"got {traction.shape!r}",
        )
    if not np.isfinite(traction).all():
        raise ECMOpenLoopError(
            "traction_non_finite",
            "traction_density_nN_per_um2 must contain only finite values",
        )
    if (traction < 0.0).any():
        raise ECMOpenLoopError(
            "traction_negative",
            "traction_density_nN_per_um2 must be non-negative at every grid cell",
        )
    return traction


def _validate_stiffness_rate(value, *, expected_shape: tuple[int, int]) -> np.ndarray:
    """Signed stiffness-rate field validator. Negative entries are
    permitted (softening); only the *post-state* non-negativity is
    enforced by the caller."""

    rate = np.asarray(value, dtype=np.float64)
    if rate.shape != expected_shape:
        raise ECMOpenLoopError(
            "stiffness_rate_shape",
            f"stiffness_rate_kpa_per_s must have shape {expected_shape!r}, "
            f"got {rate.shape!r}",
        )
    if not np.isfinite(rate).all():
        raise ECMOpenLoopError(
            "stiffness_rate_non_finite",
            "stiffness_rate_kpa_per_s must contain only finite values",
        )
    return rate


def apply_prescribed_stiffness_rate(
    ecm: ECMSubstrateState,
    stiffness_rate_kpa_per_s,
    dt_s: float,
) -> ECMSubstrateState:
    """Apply a prescribed signed stiffness-rate field to ``ecm.stiffness_kpa``.

    The update is the simplest linear open-loop recorder:
    ``stiffness_new = stiffness_old + stiffness_rate_kpa_per_s · dt_s``
    with the unit chain ``[kPa/s] · [s] = [kPa]``. Negative rate
    entries are allowed (softening); the post-state
    ``stiffness_kpa`` must remain non-negative per the ECM schema.
    A would-be-negative cell raises
    ``stiffness_negative_post_update`` rather than being clamped to
    zero (no magic clipping).

    All non-stiffness ECM fields (``ligand_density``, ``fiber_density``,
    ``orientation_tensor``, ``accumulated_traction_nNs_per_um2``) are
    copied verbatim so the only field this function ever changes is
    ``stiffness_kpa``. The grid layout (``origin_um_xy``,
    ``spacing_um``, shape) is preserved verbatim. The input ``ecm``
    is never mutated; the returned state runs through
    ``ECMSubstrateState.validate()``.

    ``dt_s == 0`` and a zero rate field are both exact no-ops
    returning a fresh state with bit-equivalent ``stiffness_kpa``.

    Args:
        ecm: Existing ECM substrate state. Validated and never mutated.
        stiffness_rate_kpa_per_s: Array with shape ``ecm.grid_shape``
            and units [kPa/s]. Values must be finite. Negative entries
            permitted (softening).
        dt_s: Non-negative scalar timestep [s]. ``0`` is allowed and
            produces an exact no-op on stiffness.

    Raises:
        ECMOpenLoopError: ``dt_invalid``, ``stiffness_rate_shape``,
            ``stiffness_rate_non_finite``, or
            ``stiffness_negative_post_update`` per the gate §2 boundary
            table.
    """

    ecm.validate()
    dt = _validate_dt(dt_s)
    rate = _validate_stiffness_rate(
        stiffness_rate_kpa_per_s, expected_shape=ecm.grid_shape
    )

    new_stiffness = (
        np.asarray(ecm.stiffness_kpa, dtype=np.float64) + rate * dt
    )
    if (new_stiffness < 0.0).any():
        # Surface with input context, not the schema's generic
        # validation message. No clamping (gate §2 + §5).
        min_value = float(new_stiffness.min())
        raise ECMOpenLoopError(
            "stiffness_negative_post_update",
            f"prescribed stiffness rate would push stiffness_kpa below zero "
            f"(min post-update value = {min_value!r}); open-loop preflight "
            f"does not clamp",
        )

    next_ecm = ECMSubstrateState(
        origin_um_xy=tuple(ecm.origin_um_xy),
        spacing_um=float(ecm.spacing_um),
        stiffness_kpa=new_stiffness,
        ligand_density=np.array(ecm.ligand_density, dtype=np.float64, copy=True),
        fiber_density=np.array(ecm.fiber_density, dtype=np.float64, copy=True),
        orientation_tensor=np.array(ecm.orientation_tensor, dtype=np.float64, copy=True),
        accumulated_traction_nNs_per_um2=np.array(
            ecm.accumulated_traction_nNs_per_um2, dtype=np.float64, copy=True
        ),
        source=ecm.source,
    )
    next_ecm.validate()
    return next_ecm


def _validate_density_rate(
    value, *, expected_shape: tuple[int, int], field_name: str
) -> np.ndarray:
    """Signed density-rate field validator. Negative entries are
    permitted (thinning); only the *post-state* [0, 1] bounds are
    enforced by the caller."""

    rate = np.asarray(value, dtype=np.float64)
    if rate.shape != expected_shape:
        raise ECMOpenLoopError(
            f"{field_name}_rate_shape",
            f"{field_name}_density_rate_per_s must have shape {expected_shape!r}, "
            f"got {rate.shape!r}",
        )
    if not np.isfinite(rate).all():
        raise ECMOpenLoopError(
            f"{field_name}_rate_non_finite",
            f"{field_name}_density_rate_per_s must contain only finite values",
        )
    return rate


def apply_prescribed_density_rate(
    ecm: ECMSubstrateState,
    ligand_density_rate_per_s,
    fiber_density_rate_per_s,
    dt_s: float,
) -> ECMSubstrateState:
    """Apply prescribed signed density-rate fields to ``ligand_density``
    and ``fiber_density`` simultaneously.

    The update is the simplest linear open-loop recorder for each
    field independently:

    - ``ligand_new = ligand_old + ligand_density_rate_per_s · dt_s``
    - ``fiber_new = fiber_old + fiber_density_rate_per_s · dt_s``

    The unit chain is ``[1/s] · [s] = [dimensionless]``. Negative rate
    entries are allowed (thinning); the post-state densities must
    remain in ``[0, 1]`` per the ECM schema. Lower-bound violations
    raise ``{ligand,fiber}_density_negative_post_update``;
    upper-bound violations raise
    ``{ligand,fiber}_density_exceeds_one_post_update``. No clamping
    in either direction (gate §2 + §5).

    Strict cross-isolation: the ligand update consumes only
    ``ligand_density_rate_per_s`` and the fiber update consumes only
    ``fiber_density_rate_per_s``. The two rate fields are validated
    and applied independently.

    All non-density ECM fields (``stiffness_kpa``,
    ``orientation_tensor``, ``accumulated_traction_nNs_per_um2``,
    ``origin_um_xy``, ``spacing_um``, ``source``) are copied
    verbatim. The input ``ecm`` is never mutated; the returned state
    runs through ``ECMSubstrateState.validate()``.

    ``dt_s == 0`` and zero rate fields are exact no-ops returning a
    fresh state.

    Args:
        ecm: Existing ECM substrate state. Validated and never mutated.
        ligand_density_rate_per_s: Array with shape ``ecm.grid_shape``
            and units [1/s]. Values must be finite. Signed entries
            permitted.
        fiber_density_rate_per_s: Same shape/units as the ligand
            rate. Validated independently of the ligand rate.
        dt_s: Non-negative scalar timestep [s]. ``0`` is allowed and
            produces an exact no-op.

    Raises:
        ECMOpenLoopError: per the gate §2 boundary table.
    """

    ecm.validate()
    dt = _validate_dt(dt_s)
    ligand_rate = _validate_density_rate(
        ligand_density_rate_per_s,
        expected_shape=ecm.grid_shape,
        field_name="ligand_density",
    )
    fiber_rate = _validate_density_rate(
        fiber_density_rate_per_s,
        expected_shape=ecm.grid_shape,
        field_name="fiber_density",
    )

    new_ligand = (
        np.asarray(ecm.ligand_density, dtype=np.float64) + ligand_rate * dt
    )
    if (new_ligand < 0.0).any():
        min_value = float(new_ligand.min())
        raise ECMOpenLoopError(
            "ligand_density_negative_post_update",
            f"prescribed ligand-density rate would push ligand_density below "
            f"zero (min post-update value = {min_value!r}); open-loop "
            f"preflight does not clamp",
        )
    if (new_ligand > 1.0).any():
        max_value = float(new_ligand.max())
        raise ECMOpenLoopError(
            "ligand_density_exceeds_one_post_update",
            f"prescribed ligand-density rate would push ligand_density above "
            f"one (max post-update value = {max_value!r}); open-loop "
            f"preflight does not clamp",
        )

    new_fiber = (
        np.asarray(ecm.fiber_density, dtype=np.float64) + fiber_rate * dt
    )
    if (new_fiber < 0.0).any():
        min_value = float(new_fiber.min())
        raise ECMOpenLoopError(
            "fiber_density_negative_post_update",
            f"prescribed fiber-density rate would push fiber_density below "
            f"zero (min post-update value = {min_value!r}); open-loop "
            f"preflight does not clamp",
        )
    if (new_fiber > 1.0).any():
        max_value = float(new_fiber.max())
        raise ECMOpenLoopError(
            "fiber_density_exceeds_one_post_update",
            f"prescribed fiber-density rate would push fiber_density above "
            f"one (max post-update value = {max_value!r}); open-loop "
            f"preflight does not clamp",
        )

    next_ecm = ECMSubstrateState(
        origin_um_xy=tuple(ecm.origin_um_xy),
        spacing_um=float(ecm.spacing_um),
        stiffness_kpa=np.array(ecm.stiffness_kpa, dtype=np.float64, copy=True),
        ligand_density=new_ligand,
        fiber_density=new_fiber,
        orientation_tensor=np.array(ecm.orientation_tensor, dtype=np.float64, copy=True),
        accumulated_traction_nNs_per_um2=np.array(
            ecm.accumulated_traction_nNs_per_um2, dtype=np.float64, copy=True
        ),
        source=ecm.source,
    )
    next_ecm.validate()
    return next_ecm


def _validate_orientation_rate(
    value, *, expected_grid_shape: tuple[int, int]
) -> np.ndarray:
    """Signed orientation-tensor rate validator. Negative entries
    permitted; only the *post-state* component magnitudes are
    enforced by the caller. Symmetry is required at the input rate
    so the post-state symmetric invariant is preserved without
    auto-symmetrization."""

    expected_shape = (*expected_grid_shape, 2, 2)
    rate = np.asarray(value, dtype=np.float64)
    if rate.shape != expected_shape:
        raise ECMOpenLoopError(
            "orientation_rate_shape",
            f"orientation_rate_per_s must have shape {expected_shape!r}, "
            f"got {rate.shape!r}",
        )
    if not np.isfinite(rate).all():
        raise ECMOpenLoopError(
            "orientation_rate_non_finite",
            "orientation_rate_per_s must contain only finite values",
        )
    asymmetry = float(
        np.max(np.abs(rate - np.swapaxes(rate, -1, -2)))
    )
    if asymmetry > _ORIENTATION_SYMMETRY_TOL:
        raise ECMOpenLoopError(
            "orientation_rate_asymmetric",
            f"orientation_rate_per_s must be symmetric in its 2x2 blocks "
            f"to within {_ORIENTATION_SYMMETRY_TOL!r}; max |R - R.T| = "
            f"{asymmetry!r} (no auto-symmetrization in open-loop preflight)",
        )
    return rate


def apply_prescribed_orientation_rate(
    ecm: ECMSubstrateState,
    orientation_rate_per_s,
    dt_s: float,
) -> ECMSubstrateState:
    """Apply a prescribed signed orientation-tensor rate field to
    ``ecm.orientation_tensor``.

    The update is the simplest linear open-loop recorder:
    ``orientation_new = orientation_old + orientation_rate_per_s · dt_s``
    with the unit chain ``[1/s] · [s] = [dimensionless]``. Signed rate
    entries are allowed (positive component grows local entry,
    negative shrinks). The rate field must already be component-wise
    symmetric to within ``_ORIENTATION_SYMMETRY_TOL``; the function
    performs **no auto-symmetrization** so that a caller bug does not
    silently drift away from the schema's symmetric invariant.

    The post-state must satisfy
    ``|orientation_new[..., i, j]| <= _ORIENTATION_BOUND``
    component-wise, again per the ECM schema. Lower- and upper-bound
    violations raise distinct failure kinds
    (``orientation_negative_bound_post_update`` for components below
    ``-_ORIENTATION_BOUND``; ``orientation_exceeds_bound_post_update``
    for components above ``+_ORIENTATION_BOUND``). Exact ``±1.0``
    boundary values are accepted; no clamping in either direction.

    All non-orientation ECM fields (``stiffness_kpa``,
    ``ligand_density``, ``fiber_density``,
    ``accumulated_traction_nNs_per_um2``, ``origin_um_xy``,
    ``spacing_um``, ``source``) are copied verbatim; the input
    ``ecm`` is never mutated; the returned state runs through
    ``ECMSubstrateState.validate()``.

    ``dt_s == 0`` and a zero rate field are exact no-ops.

    Args:
        ecm: Existing ECM substrate state. Validated and never mutated.
        orientation_rate_per_s: Array with shape
            ``(*ecm.grid_shape, 2, 2)`` and units [1/s]. Values must
            be finite and component-wise symmetric.
        dt_s: Non-negative scalar timestep [s]. ``0`` is allowed and
            produces an exact no-op.

    Raises:
        ECMOpenLoopError: per the gate §2 boundary table.
    """

    ecm.validate()
    dt = _validate_dt(dt_s)
    rate = _validate_orientation_rate(
        orientation_rate_per_s, expected_grid_shape=ecm.grid_shape
    )

    new_orientation = (
        np.asarray(ecm.orientation_tensor, dtype=np.float64) + rate * dt
    )
    component_max = float(np.max(new_orientation))
    component_min = float(np.min(new_orientation))
    if component_min < -_ORIENTATION_BOUND:
        raise ECMOpenLoopError(
            "orientation_negative_bound_post_update",
            f"prescribed orientation rate would push a component below "
            f"-{_ORIENTATION_BOUND} (min post-update value = "
            f"{component_min!r}); open-loop preflight does not clamp",
        )
    if component_max > _ORIENTATION_BOUND:
        raise ECMOpenLoopError(
            "orientation_exceeds_bound_post_update",
            f"prescribed orientation rate would push a component above "
            f"+{_ORIENTATION_BOUND} (max post-update value = "
            f"{component_max!r}); open-loop preflight does not clamp",
        )

    next_ecm = ECMSubstrateState(
        origin_um_xy=tuple(ecm.origin_um_xy),
        spacing_um=float(ecm.spacing_um),
        stiffness_kpa=np.array(ecm.stiffness_kpa, dtype=np.float64, copy=True),
        ligand_density=np.array(ecm.ligand_density, dtype=np.float64, copy=True),
        fiber_density=np.array(ecm.fiber_density, dtype=np.float64, copy=True),
        orientation_tensor=new_orientation,
        accumulated_traction_nNs_per_um2=np.array(
            ecm.accumulated_traction_nNs_per_um2, dtype=np.float64, copy=True
        ),
        source=ecm.source,
    )
    next_ecm.validate()
    return next_ecm

