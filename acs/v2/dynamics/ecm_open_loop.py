"""Open-loop ECM dynamics preflight for prescribed traction history.

This module implements the smallest separated ECM dynamics mode after P1
alpha: accumulate a caller-supplied, non-negative scalar traction density
field into :class:`acs.v2.ecm_substrate.ECMSubstrateState` without closing
the loop back to cell motion, FA state, stiffness, density, or fiber
orientation.

Sanity Gate:
    1. Dimensional analysis: input ``traction_density_nN_per_um2`` is
       [nN/um^2], ``dt_s`` is [s], and the stored history
       ``accumulated_traction_nNs_per_um2`` is [nN*s/um^2]. The only update is
       ``H_new = H_old + traction_density * dt_s``. No kPa conversion, force
       per-volume comparison, CFL, Re, Ca, De, Pe, or Ma term exists in this
       open-loop accumulator.
    2. Boundary cases: ``dt_s = 0`` and zero traction are exact no-ops;
       negative, non-finite, bool, or non-scalar ``dt_s`` is rejected;
       traction fields with wrong shape, negative entries, or non-finite
       entries are rejected before any state is created.
    3. Conservation invariants: no mass, momentum, or energy state is stored
       or changed here. All ECM fields except accumulated traction are copied
       unchanged. The only intentional monotone quantity is scalar traction
       history, which can only increase because the prescribed density and
       timestep are non-negative.
    4. Numerical sanity: float64 is used for accumulation and returned arrays.
       There is no explicit stability bound because this is a linear
       open-loop recorder, not a feedback update. Grid resolution is inherited
       from ``ECMSubstrateState`` and the input traction must match that grid.
    5. Sign/sense check: positive prescribed traction increases stored ECM
       traction exposure; zero prescribed traction leaves it unchanged. The
       function deliberately does not infer traction direction or substrate
       reaction forces.
    6. Measurement-protocol consistency: this preflight records scalar
       substrate traction exposure per ECM grid cell. It is not a projected
       area metric, not a stiffness/remodeling readout, and not a closed-loop
       force applied to cell boundaries.

Magic-Number Block: no tunable constants are introduced. N/A.
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

