"""Open-loop ECM dynamics preflight for prescribed inputs.

This module covers the separated-dynamics-mode ECM updates that take
caller-supplied prescribed input fields and integrate them into
:class:`acs.v2.ecm_substrate.ECMSubstrateState` **without closing the
loop back to cell motion, FA state, or any biology-derived input**.
Closed-loop ECM activation is gated behind the six-item gate in
``docs/v2_phase1_forward_roadmap.md`` "Closed-Loop ECM Gate".

Three preflight functions live here, one per ECM field family:

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
  fiber_density_rate_per_s, dt_s)`` — 6.4-open-B (this commit):
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
       No kPa↔nN/μm² conversion, no force-per-volume comparison, no
       CFL/Re/Ca/De/Pe/Ma. Closed-loop ECM gate item 6 (the kPa↔
       nN/μm² unit-chain proof) is the closed-loop gate's
       responsibility, not this preflight's.
    2. Boundary cases. Each function rejects: shape mismatch with
       ECM grid; non-finite or bool inputs; ``dt_s < 0`` or
       non-finite; per-field sign violations on the *post-state*
       (negative traction history, negative stiffness). ``dt_s == 0``
       and zero input fields are exact no-ops returning a fresh
       state. ``apply_prescribed_stiffness_rate`` accepts negative
       rate entries (softening); only the post-state non-negativity
       is enforced.
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
    6. Measurement-protocol consistency (Hard Rule 11). The grid
       layout (origin, spacing, shape) is preserved verbatim so the
       same coordinates that downstream visualization or comparison
       reads still resolve to the same physical cell. ``frame_dump``
       round-trip already re-validates these invariants per the ECM
       schema; this module ensures it never produces an invalid
       output.

Magic-Number Block: no tunable constants. The schema bound
``_ORIENTATION_BOUND = 1.0`` for the orientation tensor lives in
``acs.v2.ecm_substrate`` and is not relevant to this preflight.
"""

from __future__ import annotations

import numpy as np

from acs.v2.ecm_substrate import ECMSubstrateState


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

