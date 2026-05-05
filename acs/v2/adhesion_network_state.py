"""V2 Adhesion-Network Layer state contract: cadherin/integrin program dynamics.

This is the v2 reformulation of the v1 Layer 3 adhesion-network state
contract. The v2 roadmap's `Layer 3` (per ``docs/00_project_vision_v2.md``)
refers to the ECM fiber network — a separate unit. This unit is the
adhesion-network layer per the original CLAUDE.md Layer 3 definition,
locked at ``docs/v2_adhesion_network_layer_locked.md``.

PI ``id=2003`` directive (verbatim wording boundary, Y16 anchor):

    This framework absorbs PI ``id=2003`` directive: integrin-β1 and
    E-cadherin are NOT a strict-opposition or conserved-binary-switch
    relationship. The two programs are partly competing axes, not
    exclusive. Coexistence (both high, mixed transitional state, both
    low) MUST remain representable. Reciprocal cross terms in the rate
    law are signal-dependent coupling, NOT enforced antagonism.

This module provides the typed state, signals, rates, diagnostics, and
result dataclasses. The ODE step function lives in
:mod:`acs.v2.dynamics.adhesion_network_dynamics`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass(frozen=True, slots=True)
class CellAdhesionNetworkState:
    """Cell-level adhesion-program state (v2 first lock).

    Upstream cell-program signal favoring cadherin (cell-cell) vs
    integrin (cell-substrate) adhesion mode. **Distinct from FA
    mechanical state (FocalAdhesionState.bound_fraction/maturity) and
    junction mechanical state (JunctionState.cadherin_proxy)** — those
    are downstream effectors / input signals, NOT the adhesion-network
    program itself (Y2).

    PI ``id=2003`` directive: NOT strict opposition; NOT conserved sum
    (Y16). Coexistence (both high, both low, mixed) is representable.
    """

    cell_id: str
    cadherin_program: float
    integrin_program: float

    def __post_init__(self) -> None:
        if not isinstance(self.cell_id, str) or not self.cell_id:
            raise ValueError("cell_id must be non-empty")
        for name, v in (
            ("cadherin_program", self.cadherin_program),
            ("integrin_program", self.integrin_program),
        ):
            if isinstance(v, bool):
                raise ValueError(
                    f"{name} must be float (not bool — Python bool subset int trap), "
                    f"got {type(v).__name__}"
                )
            if not np.isfinite(v):
                raise ValueError(f"{name} must be finite, got {v!r}")
            if not (0.0 <= float(v) <= 1.0):
                raise ValueError(f"{name} must be in [0, 1], got {v!r}")

    @property
    def phi_integrin_dominance(self) -> Optional[float]:
        """Derived diagnostic ratio. ``None`` when both programs are 0 (Y11).

        NOT a stored field; computed property only. Branch-defined, no
        epsilon hiding. Returns ``None`` for both-zero (undefined
        dominance), explicit ratio otherwise.
        """
        total = float(self.integrin_program) + float(self.cadherin_program)
        if total == 0.0:
            return None
        return float(self.integrin_program) / total


@dataclass(frozen=True, slots=True)
class AdhesionNetworkSignals:
    """Realized engagement signals input to the ODE (Y10 + Y14).

    Built from FA / Junction objects via separate signal builders
    :func:`acs.v2.dynamics.adhesion_network_dynamics.compute_fa_engagement_signal`
    and
    :func:`acs.v2.dynamics.adhesion_network_dynamics.compute_junction_engagement_signal`.
    Caller cannot pass raw counts ambiguously; signals are pre-normalized
    to [0, 1].
    """

    fa_engagement_signal: float
    junction_engagement_signal: float

    def __post_init__(self) -> None:
        for name, v in (
            ("fa_engagement_signal", self.fa_engagement_signal),
            ("junction_engagement_signal", self.junction_engagement_signal),
        ):
            if isinstance(v, bool):
                raise ValueError(
                    f"{name} must be float (not bool), got {type(v).__name__}"
                )
            if not np.isfinite(v):
                raise ValueError(f"{name} must be finite, got {v!r}")
            if not (0.0 <= float(v) <= 1.0):
                raise ValueError(f"{name} must be in [0, 1], got {v!r}")


@dataclass(frozen=True, slots=True)
class AdhesionNetworkRates:
    """All rates required, no defaults (Y6 sister with HB#4-active Y3
    ``k_active`` discipline). Production must carry external provenance
    (Cho 2020 + general epithelial spheroid literature)."""

    k_int_on: float
    k_int_off: float
    k_cad_on: float
    k_cad_off: float

    def __post_init__(self) -> None:
        for name, v in (
            ("k_int_on", self.k_int_on),
            ("k_int_off", self.k_int_off),
            ("k_cad_on", self.k_cad_on),
            ("k_cad_off", self.k_cad_off),
        ):
            if isinstance(v, bool):
                raise ValueError(
                    f"{name} must be float (not bool — Python bool subset int trap), "
                    f"got {type(v).__name__}"
                )
            if not np.isfinite(v):
                raise ValueError(f"{name} must be finite, got {v!r}")
            if not float(v) > 0.0:
                raise ValueError(f"{name} must be strictly positive, got {v!r}")


@dataclass(frozen=True, slots=True)
class AdhesionNetworkStepDiagnostics:
    """Typed diagnostics for one adhesion-network step (Y12).

    ``integrin_steady_state`` and ``cadherin_steady_state`` are
    ``None`` when the corresponding axis has zero combined rate
    (``a + b == 0``); the analytical solution gives ``x_new = x``
    (exact no-op) in that case (Y9). Echoed signals + ``dt_s`` per
    Codex ``id=2034`` FYI (run traceability).
    """

    integrin_steady_state: Optional[float]
    cadherin_steady_state: Optional[float]
    integrin_relaxation_rate_per_s: float
    cadherin_relaxation_rate_per_s: float
    fa_engagement_signal: float
    junction_engagement_signal: float
    dt_s: float


@dataclass(frozen=True, slots=True)
class AdhesionNetworkStepResult:
    """Single-step result (Y12)."""

    state: CellAdhesionNetworkState
    diagnostics: AdhesionNetworkStepDiagnostics
