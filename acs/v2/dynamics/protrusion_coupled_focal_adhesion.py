"""Protrusion-coupled focal-adhesion dynamics (6.3b).

Implements the locked 6.3b dynamics from
``docs/v2_63b_protrusion_coupling_locked.md`` and
``docs/v2_63b_protrusion_coupled_fa_sanity_gate.md``: a deterministic
wrapper around 6.3a static FA traction preflight that computes per-FA
effective maturation/binding rates from a caller-supplied protrusion
registry plus a typed multiplier table, then delegates per-FA to
:func:`acs.v2.dynamics.focal_adhesion.step_focal_adhesions_static`.

The wrapper deliberately does **not** reimplement 6.3a rate or
traction algebra inline. It only changes which effective-rate
parameters reach 6.3a and aggregates per-FA results back in input
order. The 5 physical fields of the result match 6.3a's contract
verbatim; the 6.3b-specific diagnostics extension is captured by
the locked B1 subclass result type
:class:`ProtrusionCoupledDynamicsResult` (IS-A
:class:`acs.v2.dynamics.focal_adhesion.FocalAdhesionDynamicsResult`),
whose ``diagnostics`` field is the typed
:class:`ProtrusionCoupledDynamicsDiagnostics` (IS-A
:class:`acs.v2.dynamics.focal_adhesion.FocalAdhesionDynamicsDiagnostics`)
dataclass with 4 6.3a-shared fields + 4 6.3b extension fields
(``linked_missing``, ``reciprocal_missing``,
``multiplier_histogram``, ``max_effective_rate_per_name``).

Sanity Gate scope (dynamics/protrusion_coupled_focal_adhesion.py):

- §1 dimensional: multiplier dimensionless, effective rate 1/s,
  ``dt_fa_s · max(effective_rate)`` dimensionless, traction nN
  (6.3a delegate). Spelled out in
  ``docs/v2_63b_protrusion_coupled_fa_sanity_gate.md`` §1.
- §2 boundary: every failure_kind from the gate doc §2 is raised
  before any state mutation. Typed multiplier-key validation runs
  upfront via :meth:`ProtrusionStateMultipliers.validate`; missing
  linked id raises ``linked_protrusion_missing`` per FA, before
  stepping. Reciprocal mismatch is recorded in diagnostics, not
  raised.
- §3 conservation: per-FA Newton-3 inherits exactly from 6.3a's
  delegate; aggregate Newton-3 is the sum of per-FA cancellations.
  Determinism enforced by zero RNG and pure dict.get lookups.
- §4 numerical: ``dt_fa_s · max(effective_rate over all FAs and all
  rate names) <= _DT_RATE_SAFETY_MARGIN`` (= 0.5, the same named
  module-level constant 6.3a / P1 alpha use; rationale identical
  and reused, not re-derived).
- §5 sign: multiplier ``>= 0`` enforced at validate; multiplier
  ``1.0`` produces base 6.3a behavior; multiplier ``0`` zeroes the
  effective rate but **does not** change the FA state label
  (6.3a's no-auto-state-transition contract is preserved).
- §6 measurement-protocol: 6.3a output contract is unchanged
  (cell-on-substrate inward radial-to-centroid traction by default,
  decomposition reconstructs the original); the 6.3b diagnostics
  extension lives on the typed
  :class:`ProtrusionCoupledDynamicsDiagnostics` subclass, every
  base-class field is inherited and accessible via attribute
  access, and ``__getitem__`` access is rejected per B1 (see B1
  forbidden block below).

Magic-Number Block: this module declares no tunable numeric. The
reused 0.5 stability margin is the only numerical safety factor and
is imported from :mod:`acs.v2.active_contour` rather than redefined.
The reference biology multiplier table from the lock §3 is
documentation only; runtime requires explicit caller-supplied
multipliers.

B1 typed-schema migration scope (per
``docs/v2_focal_adhesion_dynamics_result_typed_locked.md`` +
``docs/v2_focal_adhesion_dynamics_result_typed_sanity_gate.md``,
commits ``ed5c0ca`` / ``7924730`` / ``77d4e91`` / ``88eaa4b`` /
``3b10df6``): :func:`step_protrusion_coupled_focal_adhesions`
returns :class:`ProtrusionCoupledDynamicsResult` (IS-A
:class:`acs.v2.dynamics.focal_adhesion.FocalAdhesionDynamicsResult`),
``frozen=True, slots=True``; its ``diagnostics`` field is the
typed :class:`ProtrusionCoupledDynamicsDiagnostics` (IS-A
:class:`acs.v2.dynamics.focal_adhesion.FocalAdhesionDynamicsDiagnostics`),
8 fields total (4 inherited base + 4 6.3b extension). The
migration is **shape-only**; no behavior change to 6.3a/6.3b
force / rate algebra.

B1 forbidden (text-level guard layered on top of the runtime
tests in ``tests/test_v2_protrusion_coupled_focal_adhesion.py``):

- No 6.3a/6.3b force / rate algebra change (typing migration
  only).
- No HB#4 ``diagnostics_dict`` migration here (HB#4 surface is
  grandfathered; separate unit).
- No nested mapping deep-freeze: ``multiplier_histogram`` is
  ``Mapping[str, Counter]``, ``max_effective_rate_per_name`` is
  ``Mapping[str, float]`` (separate unit).
- No ``__getitem__`` shim, deprecation hybrid, or any other
  backward-compat dict access on
  :class:`ProtrusionCoupledDynamicsResult` or
  :class:`ProtrusionCoupledDynamicsDiagnostics`. Tests
  ``test_protrusion_coupled_diagnostics_no_dict_access`` enforce.
- No preemptive ``# type: ignore[assignment]`` on the subclass
  diagnostics field override (locked Y4 plain override; deferred
  to a future strict-type-check unit).
- No ``Generic[T_Diag]`` parameterization (overengineered for
  B1).
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from typing import Mapping, Optional

import numpy as np
from numpy.typing import ArrayLike

from acs.v2.active_contour import _DT_RATE_SAFETY_MARGIN
from acs.v2.dynamics.focal_adhesion import (
    FocalAdhesionDynamicsDiagnostics,
    FocalAdhesionDynamicsError,
    FocalAdhesionDynamicsParameters,
    FocalAdhesionDynamicsResult,
    step_focal_adhesions_static,
)
from acs.v2.focal_adhesion import FocalAdhesionState
from acs.v2.protrusion import ProtrusionEvent, ProtrusionState

_ALLOWED_STATE_KEYS: frozenset[str] = frozenset(
    {"growing", "stalled", "retracting", "ended"}
)
_ALLOWED_RATE_KEYS: frozenset[str] = frozenset(
    {"k_maturity_per_s", "k_bind_per_s", "k_unbind_per_s"}
)


@dataclass(frozen=True, slots=True)
class ProtrusionCoupledDynamicsDiagnostics(FocalAdhesionDynamicsDiagnostics):
    """Typed diagnostics for one 6.3b protrusion-coupled FA step.

    IS-A :class:`FocalAdhesionDynamicsDiagnostics` — the four base
    fields (``n_adhesions``, ``aggregate_cell_force_nN_xy``,
    ``aggregate_substrate_reaction_nN_xy``,
    ``max_traction_magnitude_nN``) inherit from the 6.3a base. The
    four extension fields below are 6.3b-specific.

    Attributes:
        linked_missing: count of FAs whose linked-protrusion lookup
            silently missed (currently always 0; reserved for
            future linkage-resolution paths).
        reciprocal_missing: count of FAs with resolved
            ``linked_protrusion_id`` whose linked protrusion did
            NOT list the FA in ``associated_adhesion_ids``.
        multiplier_histogram: per-rate-name :class:`Counter` of
            effective multiplier values rounded to 12 decimals
            for float64-round-off-stable comparison.
        max_effective_rate_per_name: per-rate-name maximum
            effective rate across all FAs, in ``[1/s]``.
    """

    linked_missing: int
    reciprocal_missing: int
    multiplier_histogram: Mapping[str, Counter]
    max_effective_rate_per_name: Mapping[str, float]


@dataclass(frozen=True, slots=True)
class ProtrusionCoupledDynamicsResult(FocalAdhesionDynamicsResult):
    """Outputs of one 6.3b protrusion-coupled FA step.

    IS-A :class:`FocalAdhesionDynamicsResult` — the five physical
    fields are inherited unchanged. Only the ``diagnostics`` field
    is overridden to the 6.3b-specific
    :class:`ProtrusionCoupledDynamicsDiagnostics` subclass per
    locked Y4 (plain override, no
    ``# type: ignore[assignment]``).
    """

    diagnostics: ProtrusionCoupledDynamicsDiagnostics


@dataclass(frozen=True)
class ProtrusionStateMultipliers:
    """Caller-supplied per-state multiplier table for 6.3b.

    No project default. Unknown state or rate keys raise at
    :meth:`validate`. Missing known keys return neutral 1.0
    multiplier (absent = neutral, but spelling enforced — prevents
    typo silent-neutral bugs).

    Allowed state keys (exactly): ``"growing"``, ``"stalled"``,
    ``"retracting"``, ``"ended"`` (the four ``ProtrusionState``
    enum values).

    Allowed rate keys (exactly): ``"k_maturity_per_s"``,
    ``"k_bind_per_s"``, ``"k_unbind_per_s"`` (the three
    :class:`FocalAdhesionDynamicsParameters` rate field names).

    Multiplier values must be finite, non-bool, ``>= 0``. Negative
    multipliers, NaN, ``+/-inf``, and Python ``bool`` are rejected.
    """

    table: Mapping[str, Mapping[str, float]]

    def validate(self) -> "ProtrusionStateMultipliers":
        for state_key, rate_subtable in self.table.items():
            if state_key not in _ALLOWED_STATE_KEYS:
                raise FocalAdhesionDynamicsError(
                    "multiplier_table_unknown_key",
                    f"unknown protrusion state key {state_key!r}; "
                    f"allowed: {sorted(_ALLOWED_STATE_KEYS)}",
                )
            for rate_key, value in rate_subtable.items():
                if rate_key not in _ALLOWED_RATE_KEYS:
                    raise FocalAdhesionDynamicsError(
                        "multiplier_table_unknown_key",
                        f"unknown rate key {rate_key!r} under state "
                        f"{state_key!r}; allowed: "
                        f"{sorted(_ALLOWED_RATE_KEYS)}",
                    )
                if isinstance(value, bool) or not isinstance(value, (int, float)):
                    raise FocalAdhesionDynamicsError(
                        "multiplier_table_unknown_key",
                        f"multiplier value for state {state_key!r} rate "
                        f"{rate_key!r} must be a finite non-negative "
                        f"non-bool scalar, got {type(value).__name__}",
                    )
                v = float(value)
                if not math.isfinite(v):
                    raise FocalAdhesionDynamicsError(
                        "multiplier_table_unknown_key",
                        f"multiplier value for state {state_key!r} rate "
                        f"{rate_key!r} must be finite, got {value!r}",
                    )
                if v < 0.0:
                    raise FocalAdhesionDynamicsError(
                        "multiplier_table_unknown_key",
                        f"multiplier value for state {state_key!r} rate "
                        f"{rate_key!r} must be non-negative, got {value!r}",
                    )
        return self

    def lookup(self, state: str, rate_name: str) -> float:
        """Return the multiplier for (state, rate_name); ``1.0`` if
        either key is absent. Caller is expected to call
        :meth:`validate` before any lookup so unknown spellings have
        already raised."""

        return float(self.table.get(state, {}).get(rate_name, 1.0))


def _effective_rate(
    base_rate: Optional[float], multiplier: float
) -> Optional[float]:
    """Return ``base_rate * multiplier`` or ``None`` when the base
    rate is absent (preserving 6.3a's "absent rate ⇒ no update for
    that field" contract)."""

    if base_rate is None:
        return None
    return float(base_rate) * float(multiplier)


def _build_effective_params(
    base_params: FocalAdhesionDynamicsParameters,
    multipliers: ProtrusionStateMultipliers,
    protrusion_state: str,
) -> FocalAdhesionDynamicsParameters:
    """Construct a per-FA :class:`FocalAdhesionDynamicsParameters`
    with effective rates, leaving ``traction_scale_nN`` /
    ``max_traction_nN`` / ``dt_fa_s`` unchanged."""

    return FocalAdhesionDynamicsParameters(
        dt_fa_s=base_params.dt_fa_s,
        traction_scale_nN=base_params.traction_scale_nN,
        k_maturity_per_s=_effective_rate(
            base_params.k_maturity_per_s,
            multipliers.lookup(protrusion_state, "k_maturity_per_s"),
        ),
        k_bind_per_s=_effective_rate(
            base_params.k_bind_per_s,
            multipliers.lookup(protrusion_state, "k_bind_per_s"),
        ),
        k_unbind_per_s=_effective_rate(
            base_params.k_unbind_per_s,
            multipliers.lookup(protrusion_state, "k_unbind_per_s"),
        ),
        max_traction_nN=base_params.max_traction_nN,
    )


def _max_effective_rate_per_name(
    effective_params_per_fa: list[FocalAdhesionDynamicsParameters],
) -> dict[str, float]:
    """Per-rate-name maximum effective rate across all FAs. Absent
    rate (None) is treated as 0.0 (no contribution to the gate)."""

    def _value(p: FocalAdhesionDynamicsParameters, name: str) -> float:
        v = getattr(p, name)
        return float(v) if v is not None else 0.0

    return {
        name: max(
            (_value(p, name) for p in effective_params_per_fa), default=0.0
        )
        for name in ("k_maturity_per_s", "k_bind_per_s", "k_unbind_per_s")
    }


def step_protrusion_coupled_focal_adhesions(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    centroid_um_xy: ArrayLike,
    params: FocalAdhesionDynamicsParameters,
    protrusion_registry: Mapping[str, ProtrusionEvent],
    multipliers: ProtrusionStateMultipliers,
    *,
    traction_axis_xy_per_fa: Optional[list[Optional[ArrayLike]]] = None,
) -> ProtrusionCoupledDynamicsResult:
    """Advance a list of adhesions by one 6.3b protrusion-coupled FA
    dynamics step.

    Per-FA the wrapper looks up the FA's ``linked_protrusion_id`` in
    ``protrusion_registry`` and reads the ``state`` of the linked
    protrusion. It then constructs per-FA effective
    :class:`FocalAdhesionDynamicsParameters` whose maturation /
    binding / unbinding rates are
    ``base_rate · multipliers.lookup(state, rate_name)``, and
    delegates to
    :func:`acs.v2.dynamics.focal_adhesion.step_focal_adhesions_static`.
    The traction algebra and state-table semantics are 6.3a's
    unchanged contract.

    Failure modes (raised before any stepping):
    - ``multiplier_table_unknown_key`` on unknown keys / bool /
      non-finite / negative multiplier values (caller should call
      ``multipliers.validate()`` upfront; this wrapper does so before
      anything else).
    - ``linked_protrusion_missing`` when an FA's
      ``linked_protrusion_id`` is set but is not present in
      ``protrusion_registry``.
    - ``linked_protrusion_invalid_state`` when a linked
      ``ProtrusionEvent.state`` value is not one of the typed
      ``ProtrusionState`` enum members. Detected/manual records with
      a misspelled state must fail loudly here rather than silently
      becoming neutral 1.0 multiplier.
    - ``dt_rate_violation`` when ``dt_fa_s · max(effective_rate over
      all FAs and rate names) > _DT_RATE_SAFETY_MARGIN``. Auto-shrink
      is forbidden — the caller must reduce ``dt_fa_s`` or the
      multiplier table.
    - All inherited 6.3a failures (``rate_invalid``,
      ``traction_scale_invalid``, ``dt_invalid``, ``centroid_shape``,
      ``radial_axis_undefined``, etc.) propagate from the per-FA
      delegate without reinterpretation.

    FAs whose ``linked_protrusion_id is None`` are treated as
    "no protrusion influence" — multiplier ``1.0`` for every rate,
    equivalent to base 6.3a behavior for that FA.

    Diagnostics (additive — every 6.3a diagnostic key is preserved):
    - ``linked_missing``: count of FAs whose linked id was missing
      from the registry. Always ``0`` on success because the wrapper
      raises before stepping; kept for symmetry with
      ``reciprocal_missing``.
    - ``reciprocal_missing``: count of FAs whose
      ``linked_protrusion_id`` resolved but whose linked protrusion
      did not list the FA in ``associated_adhesion_ids``.
    - ``multiplier_histogram``: per-rate-name :class:`Counter` of
      effective multiplier values (rounded to 12 decimals to make
      the histogram comparable across float64 round-off).
    - ``max_effective_rate_per_name``: per-rate-name maximum
      effective rate across all FAs, in ``[1/s]``.

    Output array order matches input adhesion order.
    """

    multipliers.validate()

    # Base-params validation BEFORE the empty-FA early return path so an
    # invalid scalar (bool dt_fa_s, negative traction_scale_nN, etc.)
    # cannot silently succeed on `adhesions=()`. We split the validation
    # into two passes: scalar/max-traction first with rates washed to
    # None (so the 6.3a dt-rate gate is vacuous), then rate-value
    # validation with dt washed to 0 (so the dt-rate gate is again
    # vacuous). 6.3b runs its own stronger dt-rate gate over the
    # effective rates further down.
    FocalAdhesionDynamicsParameters(
        dt_fa_s=params.dt_fa_s,
        traction_scale_nN=params.traction_scale_nN,
        k_maturity_per_s=None,
        k_bind_per_s=None,
        k_unbind_per_s=None,
        max_traction_nN=params.max_traction_nN,
    ).validate()
    FocalAdhesionDynamicsParameters(
        dt_fa_s=0.0,
        traction_scale_nN=params.traction_scale_nN,
        k_maturity_per_s=params.k_maturity_per_s,
        k_bind_per_s=params.k_bind_per_s,
        k_unbind_per_s=params.k_unbind_per_s,
        max_traction_nN=params.max_traction_nN,
    ).validate()

    adhesions_tuple = tuple(adhesions)
    n = len(adhesions_tuple)

    if traction_axis_xy_per_fa is not None and len(traction_axis_xy_per_fa) != n:
        raise FocalAdhesionDynamicsError(
            "axis_shape",
            f"traction_axis_xy_per_fa length {len(traction_axis_xy_per_fa)} "
            f"!= adhesions length {n}",
        )

    # Per-FA: resolve linked protrusion (raise on missing), record
    # reciprocal mismatch as diagnostic, build effective params.
    effective_params_per_fa: list[FocalAdhesionDynamicsParameters] = []
    multiplier_values_per_name: dict[str, list[float]] = {
        "k_maturity_per_s": [],
        "k_bind_per_s": [],
        "k_unbind_per_s": [],
    }
    reciprocal_missing = 0

    for fa in adhesions_tuple:
        linked_id = fa.linked_protrusion_id
        if linked_id is None:
            protrusion_state: str = ""  # sentinel: no lookup attempted
        else:
            protrusion = protrusion_registry.get(linked_id)
            if protrusion is None:
                raise FocalAdhesionDynamicsError(
                    "linked_protrusion_missing",
                    f"FA {fa.adhesion_id!r} declares "
                    f"linked_protrusion_id={linked_id!r} but the "
                    f"protrusion_registry does not contain that id",
                )
            protrusion_state = protrusion.state
            # Hard-fail on a runtime protrusion.state outside the typed
            # ProtrusionState enum: otherwise the multiplier lookup
            # would silently fall through to neutral 1.0 and a
            # detected/manual record with a misspelled state would be
            # treated as base 6.3a behavior. The typed-state contract
            # (lock §2) demands a loud failure here.
            if protrusion_state not in _ALLOWED_STATE_KEYS:
                raise FocalAdhesionDynamicsError(
                    "linked_protrusion_invalid_state",
                    f"FA {fa.adhesion_id!r} links to protrusion "
                    f"{linked_id!r} whose state {protrusion_state!r} is "
                    f"not in the typed ProtrusionState enum "
                    f"{sorted(_ALLOWED_STATE_KEYS)}",
                )
            if fa.adhesion_id not in protrusion.associated_adhesion_ids:
                reciprocal_missing += 1

        if linked_id is None:
            effective = params  # neutral identity — keep base params
            for name in multiplier_values_per_name:
                multiplier_values_per_name[name].append(1.0)
        else:
            effective = _build_effective_params(
                params, multipliers, protrusion_state
            )
            for name in multiplier_values_per_name:
                multiplier_values_per_name[name].append(
                    multipliers.lookup(protrusion_state, name)
                )
        effective_params_per_fa.append(effective)

    # Pre-step dt-rate gate over all FAs and all rate names. This is a
    # stronger guarantee than per-FA delegation's check because the
    # 6.3a per-FA validate() only sees its own params; a heterogeneous
    # FA list could otherwise hide one violating FA behind passing
    # ones.
    max_per_name = _max_effective_rate_per_name(effective_params_per_fa)
    overall_max_effective_rate = max(max_per_name.values(), default=0.0)
    dt_s = float(params.dt_fa_s)
    if dt_s * overall_max_effective_rate > _DT_RATE_SAFETY_MARGIN:
        raise FocalAdhesionDynamicsError(
            "dt_rate_violation",
            f"dt_fa_s · max(effective_rate over all FAs and rate names) = "
            f"{dt_s * overall_max_effective_rate!r} exceeds safety margin "
            f"{_DT_RATE_SAFETY_MARGIN!r}; reduce dt_fa_s or supplied "
            f"multiplier table (auto-shrink is not allowed)",
        )

    # Per-FA delegation. Each call delegates a single-FA list to 6.3a;
    # this preserves input order trivially and avoids any inline
    # re-implementation of 6.3a's algebra. The grouping optimization
    # (group FAs by identical effective param tuple, call 6.3a once
    # per group) is permitted by the lock as behavior-equivalent but
    # is not implemented in this first version.
    cell_force = np.zeros((n, 2), dtype=np.float64)
    substrate_reaction = np.zeros((n, 2), dtype=np.float64)
    radial_components = np.zeros((n,), dtype=np.float64)
    tangential_components = np.zeros((n,), dtype=np.float64)
    updated: list[FocalAdhesionState] = []

    for i, (fa, effective) in enumerate(
        zip(adhesions_tuple, effective_params_per_fa)
    ):
        per_axis = (
            [traction_axis_xy_per_fa[i]]
            if traction_axis_xy_per_fa is not None
            else None
        )
        single = step_focal_adhesions_static(
            (fa,),
            centroid_um_xy,
            effective,
            traction_axis_xy_per_fa=per_axis,
        )
        cell_force[i] = single.cell_force_nN_xy[0]
        substrate_reaction[i] = single.substrate_reaction_nN_xy[0]
        radial_components[i] = float(single.radial_components[0])
        tangential_components[i] = float(single.tangential_components[0])
        updated.append(single.updated_adhesions[0])

    multiplier_histogram = {
        name: Counter(round(v, 12) for v in values)
        for name, values in multiplier_values_per_name.items()
    }

    cell_force_sum = cell_force.sum(axis=0).tolist()
    substrate_reaction_sum = substrate_reaction.sum(axis=0).tolist()
    diagnostics = ProtrusionCoupledDynamicsDiagnostics(
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
        linked_missing=0,
        reciprocal_missing=reciprocal_missing,
        multiplier_histogram=multiplier_histogram,
        max_effective_rate_per_name=max_per_name,
    )

    return ProtrusionCoupledDynamicsResult(
        updated_adhesions=tuple(updated),
        cell_force_nN_xy=cell_force,
        substrate_reaction_nN_xy=substrate_reaction,
        radial_components=radial_components,
        tangential_components=tangential_components,
        diagnostics=diagnostics,
    )
