"""Focal adhesion state schema.

Canonical Phase 1 ``FocalAdhesionState`` per
``docs/v2_phase1_plan_consolidated.md`` §5.5. Adhesions live on the
substrate plane in v2, so position is two-dimensional (xy in
micrometers). Per Plan §5.5, the third (z) coordinate is implied zero
and intentionally not part of the schema; legacy code that used a 3D
position must drop the z component when migrating.

Dynamics (state-machine transition rates, traction generation, slip /
release rules) are explicitly out of P0 scope and gated behind their
own Sanity Gate before execution. This schema only validates state
identity and unit invariants.

Legacy compatibility:

- ``acs.v2.single_cell`` re-imports this class so existing imports keep
  working.
- ``maturity`` and ``bound_fraction`` are preserved as numeric scalar
  state in ``[0, 1]`` per the original v2 schema.

Sanity Gate scope: non-physics schema module. The full physics 6-item
gate is N/A; this module only owns the boundary-case checks (state
literal, position finite, age non-negative, [0,1] bounds), units
(``traction_force_nN_xy`` in nN), and the measurement-protocol
consistency that ``state="unbound"`` and ``bound_fraction == 0`` are
the only legitimate ``traction_force_nN_xy == (0, 0)`` configuration
under the no-traction-without-attachment rule (Plan §6.3).

Magic-Number Block: this module declares no tunable numeric. N/A.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, Optional

FocalAdhesionStateLabel = Literal["unbound", "nascent", "mature", "slipping", "released"]
FocalAdhesionSource = Literal["simulated", "marker_detected", "manual"]


@dataclass(frozen=True)
class FocalAdhesionState:
    """Schema for one focal adhesion observation or simulation patch."""

    adhesion_id: str
    cell_id: str
    position_um_xy: tuple[float, float]
    age_s: float
    maturity: float
    bound_fraction: float
    state: FocalAdhesionStateLabel = "nascent"
    traction_force_nN_xy: tuple[float, float] = (0.0, 0.0)
    """Per-FA cell-on-substrate traction vector in nN. The substrate-on-cell
    reaction is the negation; Newton 3 holds per FA. Convention locked in
    `docs/v2_focal_adhesion_dynamics_sanity_gate.md` §6.

    State conventions (locked in the same gate, §5 state table):

    - ``state == "unbound"``: ``bound_fraction`` is forced to 0; traction
      is forced to (0, 0). The schema validator enforces this on
      construction.
    - ``state == "released"``: terminal absorbing state; the 6.3a dynamics
      step forces traction to (0, 0) regardless of stored
      ``bound_fraction``. Schema permits ``bound_fraction > 0`` to record
      a "what was bound right before release" history value; downstream
      consumers must read ``traction_force_nN_xy``, not ``bound_fraction``,
      for the active force.
    - ``state == "slipping"``: traction follows §5 magnitude formula but
      ``bound_fraction`` cannot increase (decreases only under explicit
      ``k_unbind_per_s``).
    - ``state in {"nascent", "mature"}``: traction follows §5 formula;
      both ``maturity`` and ``bound_fraction`` can increase under explicit
      caller-supplied rates.
    """
    linked_protrusion_id: Optional[str] = None
    source: FocalAdhesionSource = "simulated"

    def validate(self) -> None:
        if not self.adhesion_id.strip():
            raise ValueError("adhesion_id must be non-empty")
        if not self.cell_id.strip():
            raise ValueError("cell_id must be non-empty")
        if len(self.position_um_xy) != 2 or not all(
            math.isfinite(v) for v in self.position_um_xy
        ):
            raise ValueError("position_um_xy must contain exactly 2 finite values")
        if not math.isfinite(self.age_s) or self.age_s < 0.0:
            raise ValueError("age_s must be finite and non-negative")
        if not math.isfinite(self.maturity) or not (0.0 <= self.maturity <= 1.0):
            raise ValueError("maturity must be in [0, 1]")
        if not math.isfinite(self.bound_fraction) or not (0.0 <= self.bound_fraction <= 1.0):
            raise ValueError("bound_fraction must be in [0, 1]")
        if self.state not in {"unbound", "nascent", "mature", "slipping", "released"}:
            raise ValueError(
                "state must be unbound, nascent, mature, slipping, or released"
            )
        if len(self.traction_force_nN_xy) != 2 or not all(
            math.isfinite(v) for v in self.traction_force_nN_xy
        ):
            raise ValueError("traction_force_nN_xy must contain exactly 2 finite values")

        traction_magnitude = math.hypot(*self.traction_force_nN_xy)
        if self.state == "unbound" and traction_magnitude > 0.0:
            raise ValueError(
                "traction_force_nN_xy must be (0, 0) when state is 'unbound' "
                "(no traction without attachment, Plan §6.3)"
            )
        if self.state == "unbound" and self.bound_fraction > 0.0:
            raise ValueError(
                "bound_fraction must be 0 when state is 'unbound'"
            )

        if self.linked_protrusion_id is not None and not self.linked_protrusion_id.strip():
            raise ValueError("linked_protrusion_id must be non-empty when provided")
        if self.source not in {"simulated", "marker_detected", "manual"}:
            raise ValueError("source must be simulated, marker_detected, or manual")
