"""V2 Adhesion-Network Layer dynamics: cadherin/integrin program ODE step.

Locked at ``docs/v2_adhesion_network_layer_locked.md`` (A-tier full
debate, Y1-Y17 trace, PI ``id=2003`` directive integration).

PI ``id=2003`` directive (verbatim wording boundary, Y16 anchor):

    This framework absorbs PI ``id=2003`` directive: integrin-β1 and
    E-cadherin are NOT a strict-opposition or conserved-binary-switch
    relationship. The two programs are partly competing axes, not
    exclusive. Coexistence (both high, mixed transitional state, both
    low) MUST remain representable. Reciprocal cross terms in the rate
    law are signal-dependent coupling, NOT enforced antagonism.

ODE form (Y3 reciprocal coupling, Y9 exact affine analytical update):

    dI/dt = k_int_on  · S_fa       · (1 - I) - k_int_off · S_junction · I
    dC/dt = k_cad_on  · S_junction · (1 - C) - k_cad_off · S_fa       · C

Exact analytical solution per axis (Y9/Y17, sister with HB#1+#2 expm1
pattern at ``docs/v2_hard_blocker_1_2_constitutive_direction_locked.md``):

    a = k_on  · signal_on
    b = k_off · signal_off
    if a + b == 0:    x_new = x         (exact no-op for zero combined rate)
    else:             x_ss  = a / (a + b)
                      x_new = x_ss + (x - x_ss) · exp(-(a + b) · dt_s)

This first unit is **state-evolving diagnostic-only** (Y5):

- State evolves over time
- No mutation of FA rates, traction, junction forces, or contour motility
- Future feedback into FA rates is a separate A-tier unit
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from acs.v2.adhesion_network_state import (
    AdhesionNetworkRates,
    AdhesionNetworkSignals,
    AdhesionNetworkStepDiagnostics,
    AdhesionNetworkStepResult,
    CellAdhesionNetworkState,
)
from acs.v2.focal_adhesion import FocalAdhesionState
from acs.v2.junction import JunctionState


def _validate_dt_s(dt_s: float) -> None:
    """Y14 ``dt_s`` validation. ``dt_s == 0`` is a valid no-op
    (analytical solution gives ``x_new = x`` exactly)."""
    if isinstance(dt_s, bool):
        raise ValueError(
            f"dt_s must be float (not bool — Python bool subset int trap), "
            f"got {type(dt_s).__name__}"
        )
    if not np.isfinite(dt_s):
        raise ValueError(f"dt_s must be finite, got {dt_s!r}")
    if float(dt_s) < 0.0:
        raise ValueError(f"dt_s must be non-negative, got {dt_s!r}")


def _exact_affine_step_with_diagnostics(
    x: float,
    signal_on: float,
    signal_off: float,
    k_on: float,
    k_off: float,
    dt_s: float,
) -> tuple[float, Optional[float], float]:
    """Y9/Y17 exact analytical solution.

    Solves ``dx/dt = k_on*signal_on*(1 - x) - k_off*signal_off*x``.

    Returns ``(x_new, x_steady_state_or_None, relaxation_rate)``. The
    steady state is ``None`` when ``a + b == 0`` (both signals or
    both rates absent), and ``x_new == x`` exactly in that branch.
    """
    a = float(k_on) * float(signal_on)
    b = float(k_off) * float(signal_off)
    rate = a + b
    if rate == 0.0:
        return float(x), None, 0.0
    x_ss = a / rate
    x_new = x_ss + (float(x) - x_ss) * float(np.exp(-rate * float(dt_s)))
    return float(x_new), float(x_ss), float(rate)


def step_adhesion_network_state(
    state: CellAdhesionNetworkState,
    signals: AdhesionNetworkSignals,
    rates: AdhesionNetworkRates,
    dt_s: float,
) -> AdhesionNetworkStepResult:
    """One adhesion-network state step using exact affine update (Y9).

    PI ``id=2003`` directive: signal-dependent coupling, NOT strict
    opposition. Coexistence (both programs high under symmetric high
    signals) is representable; test 13 enforces.

    Args:
        state: current :class:`CellAdhesionNetworkState`.
        signals: realized engagement signals (Y10 typed contract; built
            via :func:`compute_fa_engagement_signal` /
            :func:`compute_junction_engagement_signal`).
        rates: rate constants (Y6 — caller-supplied, no defaults).
        dt_s: timestep in seconds. ``dt_s == 0`` is a valid exact
            no-op (Y14).

    Returns:
        :class:`AdhesionNetworkStepResult` with the new state and typed
        diagnostics (steady states, relaxation rates, echoed signals,
        echoed ``dt_s``).

    Raises:
        ValueError: on invalid ``dt_s`` (Y14) or invalid input dataclass
            fields (validated at construction).
    """
    _validate_dt_s(dt_s)

    integrin_new, integrin_steady_state, integrin_rate = (
        _exact_affine_step_with_diagnostics(
            x=state.integrin_program,
            signal_on=signals.fa_engagement_signal,
            signal_off=signals.junction_engagement_signal,
            k_on=rates.k_int_on,
            k_off=rates.k_int_off,
            dt_s=dt_s,
        )
    )

    cadherin_new, cadherin_steady_state, cadherin_rate = (
        _exact_affine_step_with_diagnostics(
            x=state.cadherin_program,
            signal_on=signals.junction_engagement_signal,
            signal_off=signals.fa_engagement_signal,
            k_on=rates.k_cad_on,
            k_off=rates.k_cad_off,
            dt_s=dt_s,
        )
    )

    new_state = CellAdhesionNetworkState(
        cell_id=state.cell_id,
        cadherin_program=cadherin_new,
        integrin_program=integrin_new,
    )
    diagnostics = AdhesionNetworkStepDiagnostics(
        integrin_steady_state=integrin_steady_state,
        cadherin_steady_state=cadherin_steady_state,
        integrin_relaxation_rate_per_s=integrin_rate,
        cadherin_relaxation_rate_per_s=cadherin_rate,
        fa_engagement_signal=signals.fa_engagement_signal,
        junction_engagement_signal=signals.junction_engagement_signal,
        dt_s=float(dt_s),
    )
    return AdhesionNetworkStepResult(state=new_state, diagnostics=diagnostics)


def compute_fa_engagement_signal(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
) -> float:
    """Y7 normalized FA signal.

    Returns the mean of ``bound_fraction * maturity`` over all FAs
    (branch-defined normalization, no magic blend weights). Empty FA
    list returns 0.0 (no integrin engagement signal).

    Each input FA is validated via ``fa.validate()`` before its
    contribution is computed (Codex ``id=2079`` P0 BLOCKER fix). This
    preserves the "normalized signal in [0, 1]" contract at the helper
    boundary; invalid FA schema (e.g., ``maturity > 1``) surfaces as
    a plain :class:`ValueError` from the FA schema validator instead
    of producing an out-of-range signal value.
    """
    if not adhesions:
        return 0.0
    for fa in adhesions:
        fa.validate()
    products = np.array(
        [float(fa.bound_fraction) * float(fa.maturity) for fa in adhesions],
        dtype=np.float64,
    )
    return float(products.mean())


def compute_junction_engagement_signal(
    junctions: tuple[JunctionState, ...] | list[JunctionState],
) -> float:
    """Y7 normalized junction signal, Y13 fail-closed on ``cadherin_proxy``.

    Returns the mean of ``cadherin_proxy`` over all junctions. Empty
    junction list (single-cell context) returns 0.0 (no cadherin
    engagement signal). If any junction has ``cadherin_proxy is None``,
    raises :class:`ValueError` — silent zero-coercion is forbidden;
    caller must populate ``cadherin_proxy`` or filter junctions
    before passing.

    Each input Junction is validated via ``j.validate()`` before its
    proxy is consumed (Codex ``id=2079`` P0 BLOCKER fix). This
    preserves the "normalized signal in [0, 1]" contract at the helper
    boundary; invalid junction schema (e.g., ``cadherin_proxy > 1``)
    surfaces as a plain :class:`ValueError` from the junction schema
    validator before the helper returns.
    """
    if not junctions:
        return 0.0
    proxies: list[float] = []
    for j in junctions:
        j.validate()
        if j.cadherin_proxy is None:
            raise ValueError(
                f"Junction with cadherin_proxy=None cannot contribute to "
                f"engagement signal. Caller must populate cadherin_proxy "
                f"or filter junctions before passing."
            )
        proxies.append(float(j.cadherin_proxy))
    return float(np.array(proxies, dtype=np.float64).mean())
