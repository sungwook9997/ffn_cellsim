"""Static-boundary focal-adhesion traction preflight (6.3a).

Implements the locked 6.3a dynamics from
``docs/v2/v2_focal_adhesion_dynamics_sanity_gate.md``: deterministic
explicit-rate maturation/binding updates plus per-FA cell-on-substrate
traction algebra on a frozen-centroid static contour. Returns updated
adhesions, per-FA forces (cell-on-substrate and substrate-on-cell), and
radial/tangential decompositions relative to the supplied centroid.

The unit deliberately does not call any active-contour ``step()``, does
not update the ECM, does not introduce stochastic nucleation or
lifetime samplers, does not auto-transition state labels, and does not
clip computed traction. Closed-loop FA / ECM coupling is staged behind
its own Sanity Gate.

Sanity Gate scope (dynamics/focal_adhesion.py):

- §1 dimensional: traction in nN; rates in 1/s; ``dt_fa_s`` in s; the
  per-FA update reductions are spelled out in
  ``docs/v2/v2_focal_adhesion_dynamics_sanity_gate.md`` §1.
- §2 boundary: every failure_kind from the gate doc §2 is raised before
  any state mutation; bool/non-finite/negative/wrong-shape inputs all
  reject.
- §3 conservation: per-FA action-reaction is exact (substrate reaction
  is the algebraic negation of the cell-on-substrate traction);
  aggregate Newton-3 sum is therefore (0, 0) by construction.
- §4 numerical: ``dt_fa_s · max(supplied rates) ≤ _DT_RATE_SAFETY_MARGIN``
  (= 0.5, the same named module-level constant the active-contour gate
  uses; rationale identical and reused, not re-derived).
- §5 sign: cell-on-substrate traction is inward (toward centroid) under
  the default radial axis; substrate-on-cell reaction is outward;
  per-state rate policy follows the §5 table (slipping never increases
  ``bound_fraction``; unbound/released traction forced to (0, 0)).
- §6 measurement-protocol: radial/tangential decomposition uses
  ``n_hat_radial_outward`` so the inward default traction has
  ``radial_components < 0``; the round-trip identity
  ``radial · n_hat_radial_outward + tangential · n_hat_tangential_ccw
  == traction_force_nN_xy`` holds to float64 round-off.

Magic-Number Block: this module declares no tunable numeric. The reused
0.5 stability margin is the only numerical safety factor and is
imported from ``acs.v2.active_contour`` rather than redefined.

B1 typed-schema migration scope (per
``docs/v2/v2_focal_adhesion_dynamics_result_typed_locked.md`` +
``docs/v2/v2_focal_adhesion_dynamics_result_typed_sanity_gate.md``,
commits ``ed5c0ca`` / ``7924730`` / ``77d4e91`` / ``88eaa4b`` /
``3b10df6``): :class:`FocalAdhesionDynamicsResult` is
``frozen=True, slots=True`` and its ``diagnostics`` field is the
typed :class:`FocalAdhesionDynamicsDiagnostics` dataclass (4
fields: ``n_adhesions``, ``aggregate_cell_force_nN_xy``,
``aggregate_substrate_reaction_nN_xy``,
``max_traction_magnitude_nN``). The migration is **shape-only**;
no behavior change to 6.3a force / rate algebra.

B1 forbidden (text-level guard layered on top of the runtime
tests in ``tests/v2/test_v2_focal_adhesion_dynamics.py``):

- No 6.3a force / rate algebra change (typing migration only).
- No HB#4 ``diagnostics_dict`` migration here (HB#4 surface is
  grandfathered; separate unit).
- No nested mapping deep-freeze (``multiplier_histogram`` /
  ``max_effective_rate_per_name`` stay ``Mapping[...]`` in the
  6.3b extension; separate unit).
- No ``__getitem__`` shim, deprecation hybrid, or any other
  backward-compat dict access on :class:`FocalAdhesionDynamicsResult`
  or :class:`FocalAdhesionDynamicsDiagnostics`. Tests
  ``test_focal_adhesion_diagnostics_no_dict_access`` enforce.
- No preemptive ``# type: ignore[assignment]`` on any subclass
  diagnostics override (deferred to a future strict-type-check
  unit per locked Y4).
- No ``Generic[T_Diag]`` parameterization (overengineered for
  B1).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import ArrayLike

from acs.v2.active_contour import _DT_RATE_SAFETY_MARGIN
from acs.v2.focal_adhesion import FocalAdhesionState


class FocalAdhesionDynamicsError(ValueError):
    """Validation error carrying a machine-readable failure kind."""

    def __init__(self, failure_kind: str, message: str) -> None:
        self.failure_kind = failure_kind
        super().__init__(message)


@dataclass(frozen=True)
class FocalAdhesionDynamicsParameters:
    """Caller-supplied 6.3a static FA dynamics parameters.

    No project defaults: every value must be explicit. Optional rates
    that are absent mean "no rate-driven update for that field" (per
    the locked §5 state table). ``max_traction_nN`` is validation-only
    and never clamps; computed traction exceeding it raises
    ``max_traction_exceeded``.
    """

    dt_fa_s: float
    traction_scale_nN: float
    k_maturity_per_s: Optional[float] = None
    k_bind_per_s: Optional[float] = None
    k_unbind_per_s: Optional[float] = None
    max_traction_nN: Optional[float] = None

    def validate(self) -> "FocalAdhesionDynamicsParameters":
        if isinstance(self.dt_fa_s, bool) or not isinstance(self.dt_fa_s, (int, float)):
            raise FocalAdhesionDynamicsError("dt_invalid", "dt_fa_s must be a finite non-negative scalar")
        dt = float(self.dt_fa_s)
        if not math.isfinite(dt) or dt < 0.0:
            raise FocalAdhesionDynamicsError("dt_invalid", "dt_fa_s must be finite and non-negative")

        if (
            isinstance(self.traction_scale_nN, bool)
            or not isinstance(self.traction_scale_nN, (int, float))
        ):
            raise FocalAdhesionDynamicsError(
                "traction_scale_invalid",
                "traction_scale_nN must be a finite non-negative scalar",
            )
        scale = float(self.traction_scale_nN)
        if not math.isfinite(scale) or scale < 0.0:
            raise FocalAdhesionDynamicsError(
                "traction_scale_invalid", "traction_scale_nN must be finite and non-negative"
            )

        for name, value in (
            ("k_maturity_per_s", self.k_maturity_per_s),
            ("k_bind_per_s", self.k_bind_per_s),
            ("k_unbind_per_s", self.k_unbind_per_s),
        ):
            if value is None:
                continue
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise FocalAdhesionDynamicsError(
                    "rate_invalid", f"{name} must be a finite non-negative scalar"
                )
            v = float(value)
            if not math.isfinite(v) or v < 0.0:
                raise FocalAdhesionDynamicsError(
                    "rate_invalid", f"{name} must be finite and non-negative"
                )

        if self.max_traction_nN is not None:
            if (
                isinstance(self.max_traction_nN, bool)
                or not isinstance(self.max_traction_nN, (int, float))
            ):
                raise FocalAdhesionDynamicsError(
                    "max_traction_invalid",
                    "max_traction_nN must be a finite non-negative scalar (or None)",
                )
            mt = float(self.max_traction_nN)
            if not math.isfinite(mt) or mt < 0.0:
                raise FocalAdhesionDynamicsError(
                    "max_traction_invalid",
                    "max_traction_nN must be finite and non-negative when provided",
                )

        rates = [
            r for r in (self.k_maturity_per_s, self.k_bind_per_s, self.k_unbind_per_s)
            if r is not None
        ]
        max_rate = float(max(rates)) if rates else 0.0
        if dt * max_rate > _DT_RATE_SAFETY_MARGIN:
            raise FocalAdhesionDynamicsError(
                "dt_rate_violation",
                f"dt_fa_s · max(rate) = {dt * max_rate!r} exceeds safety margin "
                f"{_DT_RATE_SAFETY_MARGIN!r}; reduce dt_fa_s or supplied rates "
                f"(auto-shrink is not allowed)",
            )
        return self


@dataclass(frozen=True, slots=True)
class FocalAdhesionDynamicsDiagnostics:
    """Typed diagnostics for one 6.3a static FA dynamics step.

    All fields are populated by ``step_focal_adhesions_static``;
    units carry per locked plan
    ``docs/v2/v2_focal_adhesion_dynamics_result_typed_locked.md`` §3.

    Attributes:
        n_adhesions: count of FAs in the step (dimensionless).
        aggregate_cell_force_nN_xy: 2-vector, each component in
            ``nN``, summed across FAs.
        aggregate_substrate_reaction_nN_xy: 2-vector, each
            component in ``nN``, Newton-3 partner of the cell
            force aggregate.
        max_traction_magnitude_nN: max L2 norm of per-FA cell
            force in ``nN``; ``0.0`` when ``n_adhesions == 0``.
    """

    n_adhesions: int
    aggregate_cell_force_nN_xy: tuple[float, float]
    aggregate_substrate_reaction_nN_xy: tuple[float, float]
    max_traction_magnitude_nN: float


@dataclass(frozen=True, slots=True)
class FocalAdhesionDynamicsResult:
    """Outputs of one 6.3a static FA dynamics step."""

    updated_adhesions: tuple[FocalAdhesionState, ...]
    cell_force_nN_xy: np.ndarray
    substrate_reaction_nN_xy: np.ndarray
    radial_components: np.ndarray
    tangential_components: np.ndarray
    diagnostics: FocalAdhesionDynamicsDiagnostics


def _validate_centroid(centroid_um_xy: ArrayLike) -> np.ndarray:
    arr = np.asarray(centroid_um_xy, dtype=np.float64)
    if arr.shape != (2,):
        raise FocalAdhesionDynamicsError(
            "centroid_shape", f"centroid_um_xy must have shape (2,), got {arr.shape!r}"
        )
    if not np.isfinite(arr).all():
        raise FocalAdhesionDynamicsError(
            "non_finite_centroid", "centroid_um_xy must contain only finite values"
        )
    return arr


def _validate_axis(axis_xy: ArrayLike) -> np.ndarray:
    arr = np.asarray(axis_xy, dtype=np.float64)
    if arr.shape != (2,):
        raise FocalAdhesionDynamicsError(
            "axis_shape", f"traction_axis_xy must have shape (2,), got {arr.shape!r}"
        )
    if not np.isfinite(arr).all():
        raise FocalAdhesionDynamicsError(
            "axis_shape", "traction_axis_xy must contain only finite values"
        )
    if math.hypot(*arr) == 0.0:
        raise FocalAdhesionDynamicsError(
            "axis_zero_magnitude", "traction_axis_xy must have non-zero magnitude"
        )
    return arr


def _radial_unit(position: np.ndarray, centroid: np.ndarray) -> np.ndarray:
    """Outward radial unit vector at `position` w.r.t. `centroid`. Raises
    `radial_axis_undefined` when position == centroid."""

    delta = position - centroid
    norm = math.hypot(*delta)
    if norm == 0.0:
        raise FocalAdhesionDynamicsError(
            "radial_axis_undefined",
            "FA position equals centroid; radial decomposition is undefined",
        )
    return delta / norm


def compute_radial_tangential_decomposition(
    traction_xy: ArrayLike,
    position_xy: ArrayLike,
    centroid_xy: ArrayLike,
) -> tuple[float, float]:
    """Return ``(radial, tangential)`` components in nN.

    `radial` = traction · n_hat_radial_outward
    `tangential` = traction · n_hat_tangential_ccw
    where `n_hat_tangential_ccw = (-n_y, n_x)` (90° CCW rotation of
    the outward radial). Round-trip
    ``radial · n_hat_radial_outward + tangential · n_hat_tangential_ccw``
    equals the input traction to float64 round-off.
    """

    traction = np.asarray(traction_xy, dtype=np.float64)
    if traction.shape != (2,):
        raise FocalAdhesionDynamicsError(
            "traction_shape", f"traction_xy must have shape (2,), got {traction.shape!r}"
        )
    position = np.asarray(position_xy, dtype=np.float64)
    centroid = _validate_centroid(centroid_xy)
    n_radial = _radial_unit(position, centroid)
    n_tangential = np.array([-n_radial[1], n_radial[0]], dtype=np.float64)
    return float(traction @ n_radial), float(traction @ n_tangential)


def step_focal_adhesions_static(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    centroid_um_xy: ArrayLike,
    params: FocalAdhesionDynamicsParameters,
    *,
    traction_axis_xy_per_fa: Optional[list[Optional[ArrayLike]]] = None,
) -> FocalAdhesionDynamicsResult:
    """Advance a list of adhesions by one 6.3a static FA dynamics step.

    The cell contour/centroid is a frozen geometric reference: this
    function does not call any active-contour ``step()`` and does not
    update the ECM. Per-FA per-step rate updates follow the §5 state
    table from the gate doc; absent rates mean "no change" for the
    relevant field.

    Per-FA traction direction defaults to the inward radial axis
    ``(centroid - position) / ||centroid - position||``. A per-FA
    override may be supplied via ``traction_axis_xy_per_fa`` (same
    length as ``adhesions``; entries that are ``None`` fall back to
    the default radial axis). The supplied axis is normalised
    internally; its sign is the caller's; the magnitude is always
    ``traction_scale_nN · maturity · bound_fraction``.
    """

    params.validate()
    centroid = _validate_centroid(centroid_um_xy)

    adhesions_tuple = tuple(adhesions)
    n = len(adhesions_tuple)
    if traction_axis_xy_per_fa is not None and len(traction_axis_xy_per_fa) != n:
        raise FocalAdhesionDynamicsError(
            "axis_shape",
            f"traction_axis_xy_per_fa length {len(traction_axis_xy_per_fa)} != "
            f"adhesions length {n}",
        )

    cell_force = np.zeros((n, 2), dtype=np.float64)
    substrate_reaction = np.zeros((n, 2), dtype=np.float64)
    radial_components = np.zeros((n,), dtype=np.float64)
    tangential_components = np.zeros((n,), dtype=np.float64)
    updated: list[FocalAdhesionState] = []

    for i, fa in enumerate(adhesions_tuple):
        fa.validate()
        position = np.asarray(fa.position_um_xy, dtype=np.float64)
        if not np.isfinite(position).all():
            raise FocalAdhesionDynamicsError(
                "non_finite_position", "FA position_um_xy must contain only finite values"
            )

        # Always require radial unit even when an explicit axis is
        # supplied, because the §6 decomposition uses it. This is the
        # gate's `radial_axis_undefined` guard for position == centroid.
        n_radial_outward = _radial_unit(position, centroid)

        if traction_axis_xy_per_fa is None or traction_axis_xy_per_fa[i] is None:
            # Default radial axis: inward (cell-on-substrate toward centroid).
            unit_axis = -n_radial_outward
        else:
            axis = _validate_axis(traction_axis_xy_per_fa[i])
            unit_axis = axis / math.hypot(*axis)

        new_state = fa.state
        new_maturity = fa.maturity
        new_bound = fa.bound_fraction

        # State-specific rate policy table (gate doc §5).
        if fa.state == "unbound":
            new_bound = 0.0  # forced
        elif fa.state == "released":
            pass  # no rate updates; traction forced to zero below
        elif fa.state == "nascent":
            if params.k_maturity_per_s is not None:
                new_maturity = fa.maturity + params.dt_fa_s * params.k_maturity_per_s * (
                    1.0 - fa.maturity
                )
            if params.k_bind_per_s is not None:
                new_bound = fa.bound_fraction + params.dt_fa_s * params.k_bind_per_s * (
                    1.0 - fa.bound_fraction
                )
        elif fa.state == "mature":
            if params.k_maturity_per_s is not None:
                new_maturity = fa.maturity + params.dt_fa_s * params.k_maturity_per_s * (
                    1.0 - fa.maturity
                )
            if params.k_bind_per_s is not None:
                new_bound = fa.bound_fraction + params.dt_fa_s * params.k_bind_per_s * (
                    1.0 - fa.bound_fraction
                )
        elif fa.state == "slipping":
            if params.k_unbind_per_s is not None:
                # Toward zero: bound_new = bound + dt · k · (0 - bound)
                #            = bound - dt · k · bound
                #            = bound · (1 - dt · k)
                new_bound = fa.bound_fraction * (
                    1.0 - params.dt_fa_s * params.k_unbind_per_s
                )

        # Clamp to [0, 1] — the rate updates above are linearly stable
        # under the §4 dt·rate ≤ 0.5 guard, so this is a defensive
        # numerical-floor clamp (not a magic cap).
        new_maturity = float(np.clip(new_maturity, 0.0, 1.0))
        new_bound = float(np.clip(new_bound, 0.0, 1.0))

        # Compute traction for the updated state.
        if fa.state in ("unbound", "released"):
            traction = np.zeros((2,), dtype=np.float64)
        else:
            magnitude = params.traction_scale_nN * new_maturity * new_bound
            if (
                params.max_traction_nN is not None
                and magnitude > params.max_traction_nN
            ):
                raise FocalAdhesionDynamicsError(
                    "max_traction_exceeded",
                    f"FA {fa.adhesion_id!r} computed traction magnitude {magnitude!r} "
                    f"exceeds max_traction_nN={params.max_traction_nN!r}",
                )
            traction = magnitude * unit_axis

        cell_force[i] = traction
        substrate_reaction[i] = -traction
        # Decomposition uses outward radial convention so default inward
        # traction has radial_components < 0.
        n_tangential_ccw = np.array(
            [-n_radial_outward[1], n_radial_outward[0]], dtype=np.float64
        )
        radial_components[i] = float(traction @ n_radial_outward)
        tangential_components[i] = float(traction @ n_tangential_ccw)

        updated.append(
            FocalAdhesionState(
                adhesion_id=fa.adhesion_id,
                cell_id=fa.cell_id,
                position_um_xy=fa.position_um_xy,
                age_s=fa.age_s + params.dt_fa_s,
                maturity=new_maturity,
                bound_fraction=new_bound,
                state=new_state,
                traction_force_nN_xy=(float(traction[0]), float(traction[1])),
                linked_protrusion_id=fa.linked_protrusion_id,
                source=fa.source,
            )
        )

    cell_force_sum = cell_force.sum(axis=0).tolist()
    substrate_reaction_sum = substrate_reaction.sum(axis=0).tolist()
    diagnostics = FocalAdhesionDynamicsDiagnostics(
        n_adhesions=n,
        aggregate_cell_force_nN_xy=(
            float(cell_force_sum[0]),
            float(cell_force_sum[1]),
        ),
        aggregate_substrate_reaction_nN_xy=(
            float(substrate_reaction_sum[0]),
            float(substrate_reaction_sum[1]),
        ),
        max_traction_magnitude_nN=float(
            np.linalg.norm(cell_force, axis=1).max() if n else 0.0
        ),
    )

    return FocalAdhesionDynamicsResult(
        updated_adhesions=tuple(updated),
        cell_force_nN_xy=cell_force,
        substrate_reaction_nN_xy=substrate_reaction,
        radial_components=radial_components,
        tangential_components=tangential_components,
        diagnostics=diagnostics,
    )
