"""Single-cell state schemas for v2 mechanobiology.

This module intentionally defines state and validation only. Dynamical update
rules belong in later physics modules after the data contract and sanity checks
are agreed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from acs.v2.measurement_boundary import MeasurementBoundary


@dataclass(frozen=True)
class ProtrusionEvent:
    """A lamellipodia/filopodia-like boundary event extracted or simulated."""

    time_s: float
    cell_id: str
    event_type: str
    boundary_angle_rad: float
    length_um: float
    lifetime_s: float
    confidence: float = 1.0

    def validate(self) -> None:
        if self.time_s < 0.0:
            raise ValueError("time_s must be non-negative")
        if not self.cell_id.strip():
            raise ValueError("cell_id must be non-empty")
        if self.event_type not in {"lamellipodium", "filopodium", "retraction"}:
            raise ValueError(
                "event_type must be lamellipodium, filopodium, or retraction"
            )
        if self.length_um < 0.0:
            raise ValueError("length_um must be non-negative")
        if self.lifetime_s <= 0.0:
            raise ValueError("lifetime_s must be positive")
        if not (0.0 <= self.confidence <= 1.0):
            raise ValueError("confidence must be in [0, 1]")


@dataclass(frozen=True)
class FocalAdhesionState:
    """State of one focal adhesion observation or simulated adhesion patch."""

    adhesion_id: str
    cell_id: str
    position_um_xyz: tuple[float, float, float]
    age_s: float
    maturity: float
    bound_fraction: float

    def validate(self) -> None:
        if not self.adhesion_id.strip():
            raise ValueError("adhesion_id must be non-empty")
        if not self.cell_id.strip():
            raise ValueError("cell_id must be non-empty")
        if len(self.position_um_xyz) != 3:
            raise ValueError("position_um_xyz must contain exactly 3 values")
        if self.age_s < 0.0:
            raise ValueError("age_s must be non-negative")
        if not (0.0 <= self.maturity <= 1.0):
            raise ValueError("maturity must be in [0, 1]")
        if not (0.0 <= self.bound_fraction <= 1.0):
            raise ValueError("bound_fraction must be in [0, 1]")


@dataclass
class SingleCellState:
    """Geometry and mechanobiology state for one cell at one timepoint."""

    cell_id: str
    time_s: float
    measurement_boundary: MeasurementBoundary
    height_um: Optional[float] = None
    polarity_xy: Optional[tuple[float, float]] = None
    protrusions: list[ProtrusionEvent] = field(default_factory=list)
    adhesions: list[FocalAdhesionState] = field(default_factory=list)

    def validate(self) -> None:
        if not self.cell_id.strip():
            raise ValueError("cell_id must be non-empty")
        if self.time_s < 0.0:
            raise ValueError("time_s must be non-negative")
        self.measurement_boundary.validate()
        if self.height_um is not None and self.height_um <= 0.0:
            raise ValueError("height_um must be positive when provided")
        if self.polarity_xy is not None:
            p = np.asarray(self.polarity_xy, dtype=float)
            if p.shape != (2,) or not np.isfinite(p).all():
                raise ValueError("polarity_xy must contain two finite values")
            norm = float(np.linalg.norm(p))
            if norm <= 0.0:
                raise ValueError("polarity_xy must be non-zero when provided")
        for event in self.protrusions:
            event.validate()
            if event.cell_id != self.cell_id:
                raise ValueError("all protrusion events must match SingleCellState.cell_id")
        for adhesion in self.adhesions:
            adhesion.validate()
            if adhesion.cell_id != self.cell_id:
                raise ValueError("all adhesions must match SingleCellState.cell_id")

    def projected_area_um2(self) -> float:
        """Return the measurement-boundary projected area in um2."""

        self.validate()
        return self.measurement_boundary.projected_area_um2()
