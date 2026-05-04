"""Phase D no-op closed-loop scaffolding integrator.

Implements the locked Phase D no-op design from
``docs/v2_phase_d_no_op_scaffolding_locked.md`` (commit ``18b5430`` /
``7786b20``) and the pre-execution Sanity Gate
``docs/v2_phase_d_no_op_scaffolding_sanity_gate.md`` (commit
``42b8e34``): a pure composition layer that wires Hard Blocker #3
(FA→ECM scatter primitive) and Hard Blocker #4 (ECM→FA bias
primitive) into a single integrator step entry point with three
pure functions and two typed dataclasses.

The module deliberately does **not**:

- mutate the input ECM or FA states (pure composition);
- introduce any new validation kinds beyond what HB#3 and HB#4
  raise (delegate fully);
- introduce a loose ``diagnostics_dict`` field on
  :class:`FAToECMResponseResult` (rejected in lock round 1 C1 —
  typed-contract regression);
- accept optional Phase E parameters such as ``response_law``,
  ``mapping_law``, or ``bias_params`` (silent-activation guard);
- mutate the ECM in place or deep-copy it as a pretend-update;
- scalarize the traction density at any wrapper layer (Hard Rule
  11 — Phase E response law owns the reduction);
- consume an ``effective_stiffness`` helper or any geometry-field
  shortcut (locked phased plan §3 side-door guard);
- import any Phase E module (structural absence guard);
- expose any ``step_*_active`` symbol — Phase E active behavior
  lands as a separate module under its own lock.

Sanity Gate scope (acs/v2/dynamics/closed_loop_phase_d.py):

- §1 dimensional: wrapper inherits HB#3 + HB#4 unit chains; no
  new chain at the composition layer.
- §2 boundary: empty FA list returns the locked invariants
  (zero scatter + identity ECM + ``(0, 3)`` multipliers + non-
  ``None`` empty :class:`ECMSampledAtFAs`); FA position
  validation order delegates fully to HB#3 / HB#4 sister-pattern
  (``len(position) != 2`` → unpack → ``non_finite_fa_position``
  → ``fa.validate()`` → out-of-grid).
- §3 conservation: pure composition preserves HB#3 component-
  wise scatter conservation and HB#4 structural multiplier
  invariant; triple-identity ``result.updated_ecm is
  fa_to_ecm.updated_ecm is ecm`` invariant in Phase D.
- §4 numerical: float64 throughout; per-call work = HB#3 work +
  HB#4 work + O(1) dataclass instantiation; no new tolerance
  constant.
- §5 sign: wrapper performs no arithmetic; HB#3 scatter sign and
  HB#4 multiplier sign inherited.
- §6 measurement-protocol: three-layer Hard Rule 11 guard —
  locked text (this docstring) + runtime meta-test
  ``test_phase_d_step_does_not_satisfy_phase_e_response_law`` +
  structural meta-tests
  ``test_phase_d_step_no_phase_e_kwargs`` /
  ``test_phase_d_step_no_diagnostics_dict_field``.

Magic-Number Block: zero new tunables at the wrapper layer.
Inherits HB#3 / HB#4 ``_BOUNDARY_TOL_RELATIVE = 1e-12`` only.

Phase E activation is BLOCKED on all 5 Hard Blockers + the
effective_stiffness law decision. Phase E lands as a separate
module (e.g., ``closed_loop_phase_e.py``) with its own Sanity
Gate and lock; nothing in this module authorizes any active
response or active bias law.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from acs.v2.dynamics.ecm_to_fa_bias import (
    ECMToFABiasResult,
    compute_ecm_to_fa_bias_neutral,
)
from acs.v2.dynamics.fa_to_ecm_scattering import (
    scatter_fa_traction_to_ecm_bilinear,
)
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState


@dataclass(frozen=True, slots=True)
class FAToECMResponseResult:
    """Phase D FA→ECM leg output.

    Attributes:
        traction_density_xy: ``(nx, ny, 2)`` float64 array in
            ``nN/μm²``, the bilinear-scattered per-cell vector
            traction density from
            :func:`acs.v2.dynamics.fa_to_ecm_scattering.scatter_fa_traction_to_ecm_bilinear`.
        updated_ecm: the ECM state after the Phase D FA→ECM leg.
            In Phase D the leg performs no ECM update, so this
            **is** the input ``ecm`` (object identity).
    """

    traction_density_xy: np.ndarray
    updated_ecm: ECMSubstrateState


@dataclass(frozen=True, slots=True)
class PhaseDNoOpStepResult:
    """Phase D no-op closed-loop step composite result.

    Attributes:
        fa_to_ecm: :class:`FAToECMResponseResult` from the FA→ECM
            leg.
        ecm_to_fa: :class:`ECMToFABiasResult` from the ECM→FA leg
            (neutral all-1.0 multipliers in Phase D).
        updated_ecm: the ECM state after the Phase D step. Triple-
            identity invariant in Phase D:
            ``result.updated_ecm is fa_to_ecm.updated_ecm is ecm``.
    """

    fa_to_ecm: FAToECMResponseResult
    ecm_to_fa: ECMToFABiasResult
    updated_ecm: ECMSubstrateState


def step_fa_to_ecm_response(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
) -> FAToECMResponseResult:
    """Phase D FA→ECM leg: scatter traction + identity ECM response.

    Calls
    :func:`acs.v2.dynamics.fa_to_ecm_scattering.scatter_fa_traction_to_ecm_bilinear`
    for the scatter array and returns :class:`FAToECMResponseResult`
    with the scatter array plus the input ``ecm`` as
    ``updated_ecm`` (object identity in Phase D).

    Failure: any FA validation failure from HB#3 propagates with
    the HB#3 ``failure_kind`` (``"fa_position_outside_ecm_grid"``,
    ``"non_finite_fa_position"``, ``"non_finite_fa_traction"``)
    via :class:`FAToECMScatteringError`.

    Args:
        adhesions: tuple or list of :class:`FocalAdhesionState`.
        ecm: :class:`ECMSubstrateState`. Validated by HB#3 at
            entry.

    Returns:
        :class:`FAToECMResponseResult` with the scatter array and
        the input ``ecm`` (object identity).
    """

    traction_density_xy = scatter_fa_traction_to_ecm_bilinear(adhesions, ecm)
    return FAToECMResponseResult(
        traction_density_xy=traction_density_xy,
        updated_ecm=ecm,
    )


def step_ecm_to_fa_bias(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
) -> ECMToFABiasResult:
    """Phase D ECM→FA leg: neutral bias wrapper.

    Calls
    :func:`acs.v2.dynamics.ecm_to_fa_bias.compute_ecm_to_fa_bias_neutral`
    directly. The wrapper exists so Phase E can add
    ``step_ecm_to_fa_bias_active`` as a separate function (Phase D
    / Phase E function naming separation per HB#4 lock §0
    silent-activation guard).

    Failure: any FA validation failure from HB#4 propagates with
    the HB#4 ``failure_kind`` (``"fa_bias_position_outside_ecm_grid"``,
    ``"non_finite_fa_position"``) via :class:`FAToECMBiasError`.

    Args:
        adhesions: tuple or list of :class:`FocalAdhesionState`.
        ecm: :class:`ECMSubstrateState`. Validated by HB#4 at
            entry.

    Returns:
        :class:`ECMToFABiasResult` with all-1.0 multipliers
        (shape ``(N_FA, 3)``) and non-``None`` sampled
        diagnostics.
    """

    return compute_ecm_to_fa_bias_neutral(adhesions, ecm)


def step_phase_d_no_op(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
) -> PhaseDNoOpStepResult:
    """Phase D no-op closed-loop step: composition wrapper.

    Calls :func:`step_fa_to_ecm_response` first; only if FA→ECM
    succeeds, then calls :func:`step_ecm_to_fa_bias`. Returns
    :class:`PhaseDNoOpStepResult` with both leg results plus
    ``updated_ecm`` (= input ``ecm`` in Phase D, object identity).

    The deterministic call order means that for an FA outside the
    inclusive ECM footprint, the wrapper always raises
    :class:`FAToECMScatteringError` with
    ``failure_kind="fa_position_outside_ecm_grid"`` (HB#3) — the
    HB#4 ``"fa_bias_position_outside_ecm_grid"`` failure kind is
    only reachable via a mocked/stubbed HB#4 failure scenario,
    not natural geometry (per locked test 9 amendment ``7786b20``).

    Failure: surfaces the FIRST leg's failure (HB#3 first, then
    HB#4 only if HB#3 succeeded). Wrapper introduces no new
    failure_kinds.

    Args:
        adhesions: tuple or list of :class:`FocalAdhesionState`.
        ecm: :class:`ECMSubstrateState`. Validated by both legs at
            entry.

    Returns:
        :class:`PhaseDNoOpStepResult` with both leg results and
        ``updated_ecm`` (object identity).
    """

    fa_to_ecm = step_fa_to_ecm_response(adhesions, ecm)
    ecm_to_fa = step_ecm_to_fa_bias(adhesions, ecm)
    return PhaseDNoOpStepResult(
        fa_to_ecm=fa_to_ecm,
        ecm_to_fa=ecm_to_fa,
        updated_ecm=ecm,
    )
