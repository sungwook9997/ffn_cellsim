"""V2 Phase E v1 closed-loop ECM gate composition (ECM-side evidence only).

Implements the locked Phase E v1 composition design from
``docs/v2_phase_e_composition_locked.md`` (commit ``c8b9550``)
and the impl-work Sanity Gate
``docs/v2_phase_e_composition_sanity_gate.md`` (commit
``b0b5baa`` + amendments ``a3bc8a6`` + ``e5a8186``; Codex Sanity
Gate PASS at ``id=1552``).

Phase E v1 wording boundary (paraphrased intent guard per Codex
``id=1546`` blocker resolution + ``id=1549`` example refinement):
this composition provides ECM-side evidence only; complete
bidirectional closure remains future work requiring Phase E v2,
the Item 5 sweep harness, and PI approval.

This is the **first integrated** closed-loop ECM gate composition
combining all 5 upstream hard blockers + Phase D infrastructure
into a single per-step function. Together with HB#3 (FA→ECM
scattering), HB#4 (ECM→FA bias neutral), HB#1+#2 (active
orientation update), HB#5 (Lyapunov-like metric), and Phase D
no-op orchestrator, all 5 upstream Phase E blockers are now
implemented.

Phase E v1 composition order (locked Y5):

1. HB#3 scatter: FA traction → ECM grid traction density.
2. HB#1+#2 orientation update: instantaneous traction stimulus
   drives convex orientation update; returns fresh ECM with no
   aliasing.
3. HB#4 neutral bias readout: post-HB#1+#2 ECM (forward-compat
   for v2 active variant).
4. HB#5 Lyapunov-like metric: post-HB#1+#2 ECM + the original
   HB#3 scatter output (NOT recomputed; call-order spy test
   enforces).

Identity invariant (Y4): ``result.updated_ecm is
result.orientation_response.updated_ecm`` — object identity, no
copy, no mutation in this composition wrapper.

The module deliberately does **not**:

- introduce v2 placeholder, enum, hook, strategy, or optional
  bias law parameter (Y3 — silent-activation risk; v2 gets a
  separate function lock when HB#4-active exists);
- mutate ECM in place, mutate adhesions, or copy any sub-result
  (Y4 — composition wrapper is a thin orchestrator);
- introduce any wrapper failure kind: no PhaseE*Error class, no
  failure_kind assignment in the module body (Y17; sub-call
  errors propagate verbatim);
- reference the deferred ECM→FA mechanosensing helper anywhere
  in the source body (Y8 + Y13 — string + AST guard via
  meta-test; the lock doc may name the deferred helper, but
  the .py file cannot);
- redefine TRACTION_REF_NN_PER_UM2 or K_ORIENT_PER_S locally
  (Y14 — must be imported from ecm_constitutive_response;
  meta-test enforces object identity);
- reorder composition (Y5 — locked HB#3 → HB#1+#2 → HB#4 →
  HB#5; call-order spy test enforces);
- recompute scatter inside the HB#5 call (Y5 — the original
  traction_density_xy from step 1 is passed to HB#5).

Failure-kind discipline diverges from HB#1+#2: this composition
wrapper has NO domain-specific failure modes; sub-call errors
(HB#3 / HB#1+#2 / HB#4 / HB#5) propagate as-is with their
original types and messages. Stop-before-next-leg deterministic
order ensures no partial-update contamination.

Sanity Gate scope (acs/v2/dynamics/closed_loop_phase_e.py):

- §1 dimensional: composition wrapper performs no unit
  transformation; 5 sub-call unit chains preserved by the
  pseudocode (Hard Rule 10 inline-derived in the locked §1 +
  Sanity Gate §1).
- §2 boundary: 5 boundary classes (empty FA per Y7 cite source
  line 170 / dt=0 valid no-op / dt invalid HB#1+#2 propagation
  / single-cell / stop-before-next-leg per Y6) all locked +
  test-covered.
- §3 conservation: composition stateless; identity invariant
  ``updated_ecm is orientation_response.updated_ecm`` is the
  wrapper's structural invariant; stop-before-next-leg ensures
  no partial-update contamination.
- §4 numerical: float64 throughout (inherited); no new
  tolerance per Y14 (constants imported, not redefined).
- §5 sign: no composition-introduced sign; sub-call sign
  conventions inherited.
- §6 measurement-protocol (Hard Rule 11, central anchor):
  ECM-side wording boundary (Y1) + 4 source-level meta-tests
  (Y2 + Y13 + Y14 + Y17) + Step 6 sister-gate-mirror across 4
  layers + wording-boundary.

Magic-Number Block: zero new tunables at composition layer per
Y14. Both ``TRACTION_REF_NN_PER_UM2`` and ``K_ORIENT_PER_S``
re-imported from ``ecm_constitutive_response``, NOT redefined;
meta-test enforces object identity
``ce.TRACTION_REF_NN_PER_UM2 is cr.TRACTION_REF_NN_PER_UM2``.

Phase E v2 (HB#4-active variant) and the Item 5 sweep harness
are separate later cycles; this module landing does NOT
authorize Phase E v2 composition or full bidirectional closure
itself.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState
from acs.v2.dynamics.ecm_constitutive_response import (
    ECMOrientationResponseResult,
    K_ORIENT_PER_S,
    TRACTION_REF_NN_PER_UM2,
    step_ecm_orientation_response,
)
from acs.v2.dynamics.ecm_lyapunov_metric import (
    ECMOrientationLyapunovMetricResult,
    compute_ecm_orientation_lyapunov_metric,
)
from acs.v2.dynamics.ecm_to_fa_bias import (
    ECMToFABiasResult,
    compute_ecm_to_fa_bias_active,
    compute_ecm_to_fa_bias_neutral,
    _validate_k_active as _validate_hb4_active_k_active,
)
from acs.v2.dynamics.fa_to_ecm_scattering import (
    scatter_fa_traction_to_ecm_bilinear,
)


@dataclass(frozen=True, slots=True)
class PhaseEStepResult:
    """Outputs of one Phase E v1 closed-loop ECM gate composition step.

    Identity invariant (Y4): ``updated_ecm is
    orientation_response.updated_ecm`` — object identity, no
    copy, no mutation in the composition wrapper.

    Attributes:
        traction_density_xy: ``(nx, ny, 2)`` float64 array in
            ``nN/μm²`` from HB#3 scatter
            (:func:`acs.v2.dynamics.fa_to_ecm_scattering.scatter_fa_traction_to_ecm_bilinear`).
        orientation_response: HB#1+#2 active orientation update
            output
            (:class:`acs.v2.dynamics.ecm_constitutive_response.ECMOrientationResponseResult`)
            with updated ECM + 13-field diagnostics.
        ecm_to_fa_bias: HB#4 neutral bias readout
            (:class:`acs.v2.dynamics.ecm_to_fa_bias.ECMToFABiasResult`)
            with all-1.0 multipliers (Phase D default; Phase E
            v1 hard-wires neutral). Read on the post-HB#1+#2
            ECM for forward-compat with future Phase E v2
            (HB#4-active variant).
        lyapunov_metric: HB#5 Lyapunov-like metric
            (:class:`acs.v2.dynamics.ecm_lyapunov_metric.ECMOrientationLyapunovMetricResult`)
            providing ECM-side evidence for the closed-loop ECM
            gate. Computed on the post-HB#1+#2 ECM + the
            original HB#3 scatter output (NOT a recomputed
            scatter — call-order spy test enforces).
        updated_ecm: the ECM state after the Phase E v1 step.
            Object-identity equal to ``orientation_response.updated_ecm``;
            this composition wrapper does NOT add another copy
            layer.
    """

    traction_density_xy: np.ndarray
    orientation_response: ECMOrientationResponseResult
    ecm_to_fa_bias: ECMToFABiasResult
    lyapunov_metric: ECMOrientationLyapunovMetricResult
    updated_ecm: ECMSubstrateState


def step_closed_loop_phase_e_v1(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
    dt_s: float,
    *,
    traction_ref_nN_per_um2: float = TRACTION_REF_NN_PER_UM2,
    k_orient_per_s: float = K_ORIENT_PER_S,
) -> PhaseEStepResult:
    """One Phase E v1 ECM-side closed-loop ECM gate composition step.

    This composition provides ECM-side evidence only; complete
    bidirectional closure remains future work requiring Phase E
    v2, the Item 5 sweep harness, and PI approval.

    Composes (locked Y5 deterministic order):

    1. HB#3 scatter: FA traction → ECM grid traction density.
    2. HB#1+#2 orientation update: instantaneous traction
       stimulus drives convex orientation update; returns fresh
       ECM with no aliasing per HB#1+#2 Y12.
    3. HB#4 neutral bias readout: on the post-HB#1+#2 ECM
       (forward-compat for v2 active variant).
    4. HB#5 Lyapunov-like metric: on the post-HB#1+#2 ECM +
       the original HB#3 scatter output (NOT recomputed;
       call-order spy test enforces).

    No wrapper failure kinds (Y17): sub-call errors (HB#3 /
    HB#1+#2 / HB#4 / HB#5) propagate verbatim. Stop-before-
    next-leg deterministic order ensures no partial-update
    contamination.

    Identity invariant (Y4): ``result.updated_ecm is
    result.orientation_response.updated_ecm`` — object
    identity, no copy, no mutation in the wrapper.

    Args:
        adhesions: tuple or list of
            :class:`acs.v2.focal_adhesion.FocalAdhesionState`.
            Empty list valid (per HB#3 contract); produces
            all-zero traction → identity orientation update →
            empty bias `(0, 3)` → V_active=0.
        ecm: input
            :class:`acs.v2.ecm_substrate.ECMSubstrateState`.
            Validated by HB#1+#2 entry (Y1 silent-heal-path
            closure inherited per HB#1+#2 ``id=1486``).
        dt_s: timestep in seconds. ``0.0`` is valid no-op
            (HB#1+#2 Y15); negative / boolean / non-finite
            raise from HB#1+#2.
        traction_ref_nN_per_um2: literature-pinned reference
            traction scale (default
            :data:`acs.v2.dynamics.ecm_constitutive_response.TRACTION_REF_NN_PER_UM2`).
            Object-identity re-imported per Y14; meta-test
            enforces no local redefinition.
        k_orient_per_s: alignment rate constant (default
            :data:`acs.v2.dynamics.ecm_constitutive_response.K_ORIENT_PER_S`).
            Object-identity re-imported per Y14.

    Returns:
        :class:`PhaseEStepResult` with the 5 sub-results and
        ``updated_ecm`` object-identical to
        ``orientation_response.updated_ecm`` (Y4 locked
        invariant). Note: Phase E v1 has single-identity
        only; ``result.updated_ecm is ecm`` is **False**
        (HB#1+#2 returns a fresh ECM per its Y12 no-aliasing
        contract).
    """

    traction_density_xy = scatter_fa_traction_to_ecm_bilinear(adhesions, ecm)
    orientation_response = step_ecm_orientation_response(
        ecm,
        traction_density_xy,
        dt_s,
        traction_ref_nN_per_um2=traction_ref_nN_per_um2,
        k_orient_per_s=k_orient_per_s,
    )
    ecm_to_fa_bias = compute_ecm_to_fa_bias_neutral(
        adhesions, orientation_response.updated_ecm
    )
    lyapunov_metric = compute_ecm_orientation_lyapunov_metric(
        orientation_response.updated_ecm, traction_density_xy
    )
    return PhaseEStepResult(
        traction_density_xy=traction_density_xy,
        orientation_response=orientation_response,
        ecm_to_fa_bias=ecm_to_fa_bias,
        lyapunov_metric=lyapunov_metric,
        updated_ecm=orientation_response.updated_ecm,
    )


def step_closed_loop_phase_e_v2(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
    dt_s: float,
    *,
    k_active: float,
    traction_ref_nN_per_um2: float = TRACTION_REF_NN_PER_UM2,
    k_orient_per_s: float = K_ORIENT_PER_S,
) -> PhaseEStepResult:
    """One Phase E v2 ECM-side composition step (active HB#4 in the loop).

    Phase E v2 step 2 achieves per-step **composition-structure
    closure** by placing HB#4-active non-neutral ECM->FA bias inside
    the FA->ECM->FA loop. It does **not** by itself establish Items 1-4
    multi-step satisfaction over time; that requires a sweep harness
    or simulation-engine trajectory (locked Q4 = (c) hybrid wording
    boundary).

    Composes (sister with v1 with HB#4-active swap at step 3):

    1. HB#3 scatter: FA traction -> ECM grid traction density.
    2. HB#1+#2 orientation update: instantaneous traction stimulus
       drives convex orientation update; returns fresh ECM with no
       aliasing per HB#1+#2 Y12.
    3. HB#4-**active** bias readout: deviatoric Rayleigh score on
       post-HB#1+#2 ECM (non-neutral multipliers per HB#4-active Y1).
    4. HB#5 Lyapunov-like metric: on the post-HB#1+#2 ECM + the
       original HB#3 scatter object (NOT recomputed; call-order spy
       test enforces).

    Validation order (locked Codex C3 / HB#4-active Y12 sister):
      a. ``ecm.validate()``
      b. ``_validate_hb4_active_k_active(k_active)`` — cheap parameter
         failure first (BEFORE HB#3 scatter O(N_FA) work)
      c. sub-call sequence (HB#3 -> HB#1+#2 -> HB#4-active -> HB#5)

    No wrapper failure kinds (Y17 inheritance from v1): sub-call errors
    (HB#3 / HB#1+#2 / HB#4-active / HB#5) propagate verbatim. ``k_active``
    validation surfaces the imported HB#4-active ``k_active_invalid``
    error kind from the shared validator (no v2-local error kind
    invented).

    Identity invariant Y4 (inherited from v1): ``result.updated_ecm is
    result.orientation_response.updated_ecm`` — object identity, no
    copy, no mutation in the wrapper.

    Args:
        adhesions: tuple or list of
            :class:`acs.v2.focal_adhesion.FocalAdhesionState`. Empty
            list valid (per HB#3 contract); produces all-zero traction
            -> identity orientation update -> empty bias ``(0, 3)`` ->
            V_active=0.
        ecm: input :class:`acs.v2.ecm_substrate.ECMSubstrateState`.
            Validated at function entry.
        dt_s: timestep in seconds. ``0.0`` is valid no-op (HB#1+#2
            Y15 + locked anti-collapse fixture); negative / boolean /
            non-finite raise from HB#1+#2.
        k_active: dimensionless coupling strength for HB#4-active
            deviatoric Rayleigh law. **REQUIRED, no default** (Q3 = (b)
            sister with HB#4-active Y3). Must be finite, positive, and
            produce a finite ``exp(k_active*sqrt(2))`` upper bound.
        traction_ref_nN_per_um2: literature-pinned reference traction
            scale (default
            :data:`acs.v2.dynamics.ecm_constitutive_response.TRACTION_REF_NN_PER_UM2`).
        k_orient_per_s: alignment rate constant (default
            :data:`acs.v2.dynamics.ecm_constitutive_response.K_ORIENT_PER_S`).

    Returns:
        :class:`PhaseEStepResult` (Q1 = (a) reused; HB#4-active Y15
        sister) with the 5 sub-results and ``updated_ecm`` object-
        identical to ``orientation_response.updated_ecm`` (Y4 invariant).

    Raises:
        FAToECMBiasError: reports the HB#4-active ``k_active_invalid``
            error kind if ``k_active`` is invalid (raised BEFORE any
            sub-call fires).
        Errors from sub-calls (HB#3 / HB#1+#2 / HB#4-active / HB#5)
            propagate verbatim — no wrapper-level translation.
    """

    ecm.validate()
    _validate_hb4_active_k_active(k_active)

    traction_density_xy = scatter_fa_traction_to_ecm_bilinear(adhesions, ecm)
    orientation_response = step_ecm_orientation_response(
        ecm,
        traction_density_xy,
        dt_s,
        traction_ref_nN_per_um2=traction_ref_nN_per_um2,
        k_orient_per_s=k_orient_per_s,
    )
    ecm_to_fa_bias = compute_ecm_to_fa_bias_active(
        adhesions,
        orientation_response.updated_ecm,
        k_active=k_active,
    )
    lyapunov_metric = compute_ecm_orientation_lyapunov_metric(
        orientation_response.updated_ecm, traction_density_xy
    )
    return PhaseEStepResult(
        traction_density_xy=traction_density_xy,
        orientation_response=orientation_response,
        ecm_to_fa_bias=ecm_to_fa_bias,
        lyapunov_metric=lyapunov_metric,
        updated_ecm=orientation_response.updated_ecm,
    )
