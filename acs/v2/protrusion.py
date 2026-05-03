"""Protrusion (lamellipodia / filopodia) event schema.

This module defines the canonical Phase 1 ``ProtrusionEvent`` schema per
``docs/v2_phase1_plan_consolidated.md`` §5.4. The event records *what*
happens at the cell boundary; dynamics (event hazard, growth/retraction
rules, RNG reproducibility) are explicitly out of P0 scope and must
write their own Sanity Gate before execution.

Legacy compatibility:

- ``acs.v2.single_cell`` re-imports ``ProtrusionEvent`` from this module so
  existing imports `from acs.v2.single_cell import ProtrusionEvent` keep
  working.
- ``lifetime_s`` and ``confidence`` are preserved as optional legacy
  fields. When ``end_time_s`` and ``lifetime_s`` are both set, validation
  enforces consistency.

Sanity Gate scope: schema-only module with no time stepping or numerics.
The full physics 6-item gate is N/A; this module only owns the
boundary-case checks (state literal, event_type literal, length/width
non-negative, finite times) and the measurement-protocol consistency
that ``ProtrusionEvent.boundary_angle_rad`` is in radians and positions
are in micrometers.

Magic-Number Block: this module declares no tunable numeric. N/A.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, Optional

# Numerical validation tolerance for the legacy ``lifetime_s`` <-> ``end_time_s
# - start_time_s`` consistency check. This is a numerical tie-break that
# absorbs float64 round-off, NOT a physical or fitted tunable.
_LIFETIME_CONSISTENCY_TOL_S = 1e-9

ProtrusionEventType = Literal["lamellipodium", "filopodium"]
ProtrusionState = Literal["growing", "stalled", "retracting", "ended"]
ProtrusionSource = Literal["simulated", "detected", "manual"]


@dataclass(frozen=True)
class ProtrusionEvent:
    """Schema for a single boundary protrusion (lamellipodium / filopodium).

    Attributes:
        cell_id: Identifier of the parent cell.
        start_time_s: Event start time in seconds, finite and non-negative.
        boundary_angle_rad: Angular position of the event root on the
            boundary, in radians.
        length_um: Current event length in micrometers.
        event_id: Optional event identifier; when set, must be non-empty.
        end_time_s: Optional event end time in seconds. Must be finite,
            non-negative, and ``>= start_time_s`` when present.
        event_type: ``lamellipodium`` for broad edge events,
            ``filopodium`` for narrow probing events.
        state: Phase of the event lifecycle.
        root_position_um_xy: Optional xy position of the event root in
            micrometers, world y-up.
        tip_position_um_xy: Optional xy position of the event tip.
        direction_um_xy: Optional unit-ish direction vector in xy.
        width_um: Optional event width in micrometers (lamellipodia).
        width_rad: Optional event angular width in radians (filopodia).
        force_candidate_nN: Optional candidate tip force in nanonewtons.
        associated_adhesion_ids: Optional ids of focal adhesions linked to
            this event.
        source: Where the event came from.
        lifetime_s: Optional legacy duration field. When set together
            with ``end_time_s`` it must agree with
            ``end_time_s - start_time_s`` to within 1e-9 s.
        confidence: Optional legacy detection confidence in ``[0, 1]``.
    """

    cell_id: str
    start_time_s: float
    boundary_angle_rad: float
    length_um: float
    event_id: Optional[str] = None
    end_time_s: Optional[float] = None
    event_type: ProtrusionEventType = "lamellipodium"
    state: ProtrusionState = "growing"
    root_position_um_xy: Optional[tuple[float, float]] = None
    tip_position_um_xy: Optional[tuple[float, float]] = None
    direction_um_xy: Optional[tuple[float, float]] = None
    width_um: Optional[float] = None
    width_rad: Optional[float] = None
    force_candidate_nN: Optional[float] = None
    associated_adhesion_ids: tuple[str, ...] = ()
    source: ProtrusionSource = "simulated"
    lifetime_s: Optional[float] = None
    confidence: Optional[float] = None

    def validate(self) -> None:
        if not self.cell_id.strip():
            raise ValueError("cell_id must be non-empty")
        if not math.isfinite(self.start_time_s) or self.start_time_s < 0.0:
            raise ValueError("start_time_s must be finite and non-negative")
        if not math.isfinite(self.boundary_angle_rad):
            raise ValueError("boundary_angle_rad must be finite")
        if not math.isfinite(self.length_um) or self.length_um < 0.0:
            raise ValueError("length_um must be finite and non-negative")

        if self.event_id is not None and not self.event_id.strip():
            raise ValueError("event_id must be non-empty when provided")

        if self.end_time_s is not None:
            if not math.isfinite(self.end_time_s) or self.end_time_s < 0.0:
                raise ValueError("end_time_s must be finite and non-negative when provided")
            if self.end_time_s < self.start_time_s:
                raise ValueError("end_time_s must be >= start_time_s")

        if self.event_type not in {"lamellipodium", "filopodium"}:
            raise ValueError("event_type must be lamellipodium or filopodium")
        if self.state not in {"growing", "stalled", "retracting", "ended"}:
            raise ValueError("state must be growing, stalled, retracting, or ended")
        if self.source not in {"simulated", "detected", "manual"}:
            raise ValueError("source must be simulated, detected, or manual")

        for label, xy in (
            ("root_position_um_xy", self.root_position_um_xy),
            ("tip_position_um_xy", self.tip_position_um_xy),
            ("direction_um_xy", self.direction_um_xy),
        ):
            if xy is None:
                continue
            if len(xy) != 2 or not all(math.isfinite(v) for v in xy):
                raise ValueError(f"{label} must contain exactly 2 finite values")

        if self.direction_um_xy is not None:
            if math.hypot(*self.direction_um_xy) == 0.0:
                raise ValueError(
                    "direction_um_xy must have non-zero magnitude when provided"
                )

        if self.width_um is not None:
            if not math.isfinite(self.width_um) or self.width_um < 0.0:
                raise ValueError("width_um must be finite and non-negative when provided")
        if self.width_rad is not None:
            if not math.isfinite(self.width_rad) or self.width_rad < 0.0:
                raise ValueError("width_rad must be finite and non-negative when provided")
        if self.force_candidate_nN is not None:
            if not math.isfinite(self.force_candidate_nN) or self.force_candidate_nN < 0.0:
                raise ValueError(
                    "force_candidate_nN must be finite and non-negative when provided"
                )

        seen_ids: set[str] = set()
        for fa_id in self.associated_adhesion_ids:
            if not fa_id.strip():
                raise ValueError("associated_adhesion_ids must be non-empty strings")
            if fa_id in seen_ids:
                raise ValueError(
                    f"associated_adhesion_ids must be unique, got duplicate {fa_id!r}"
                )
            seen_ids.add(fa_id)

        if self.lifetime_s is not None:
            if not math.isfinite(self.lifetime_s) or self.lifetime_s < 0.0:
                raise ValueError("lifetime_s must be finite and non-negative when provided")
            if self.end_time_s is not None:
                expected = self.end_time_s - self.start_time_s
                if abs(self.lifetime_s - expected) > _LIFETIME_CONSISTENCY_TOL_S:
                    raise ValueError(
                        f"lifetime_s {self.lifetime_s!r} disagrees with "
                        f"end_time_s - start_time_s = {expected!r}"
                    )

        if self.confidence is not None:
            if not math.isfinite(self.confidence) or not (0.0 <= self.confidence <= 1.0):
                raise ValueError("confidence must be in [0, 1] when provided")
