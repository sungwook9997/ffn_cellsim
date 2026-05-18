"""Cell dataclass (KU-3.16 integration target; week-6 freeze interface).

Phase 1 Unit 3.1 ships the minimal :class:`Cell` shape that Worker D
(junction) and Worker B (FA traction) will consume in subsequent
units. The dataclass owns the cortex, a polarity vector (for Unit
3.2 motility), optional focal-adhesion list (Worker B Unit 2.x
output), and a switch between discrete- and surface-tension cortex
representations (KU-3.1 outer-vs-inner cortex distinction).

Interface contract — **week 6 freeze**
--------------------------------------
The fields and the :meth:`compute_cortex_boundary_position` signature
below are frozen at week 6. Any change after that requires an
emergency Notion sync because Worker D imports this class to build
two-cell E-cadherin junctions and Worker B injects FocalAdhesion
items into ``focal_adhesions`` during Unit 2.2.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

from acs_kb.cell.cortex import Cortex

if TYPE_CHECKING:
    # Worker B's FocalAdhesion (Unit 2.x) — only typed, never imported at
    # runtime, so Unit 3.1 can ship without depending on Unit 2.x being
    # merged. Replace with the concrete import after Worker B's week 5
    # freeze.
    from typing import Protocol

    class FocalAdhesion(Protocol):  # noqa: D401 — protocol stub.
        ...


@dataclass(slots=True)
class Cell:
    """Phase 1 single-cell aggregate state.

    Attributes
    ----------
    id : int
        Unique cell identifier within a simulation run.
    cortex : Cortex
        Discrete bead-spring cortex (Phase 1 Unit 3.1).
    center_position : np.ndarray, shape (2,)
        Centroid of the cortex bead cloud (m). Computed from the
        cortex centroid; kept as a field so external modules (Worker
        D) can read the current value without recomputing.
    polarity : np.ndarray, shape (2,)
        Unit polarity vector (Unit 3.2 motility); defaults to +x.
    focal_adhesions : list
        Worker B's :class:`FocalAdhesion` instances (Unit 2.x). Empty
        in Phase 1 Unit 3.1.
    inner_mode : Literal["discrete", "surface_tension"]
        KU-3.1 outer-vs-inner cortex switch. Phase 1 Unit 3.1 only
        constructs ``"discrete"``; ``"surface_tension"`` is reserved
        for Phase 1 Unit 5 (spheroid interior) and triggers the
        continuum cortex code path when implemented.
    surface_tension : float
        γ used by the ``"surface_tension"`` mode (N/m); default
        KU-3.5 mid-range value.
    nucleus_position : np.ndarray, shape (2,)
        Nuclear position (m); Phase 1 places it at the centroid.
    """

    id: int
    cortex: Cortex
    center_position: np.ndarray
    polarity: np.ndarray = field(
        default_factory=lambda: np.array([1.0, 0.0], dtype=np.float64)
    )
    focal_adhesions: list[Any] = field(default_factory=list)
    inner_mode: Literal["discrete", "surface_tension"] = "discrete"
    surface_tension: float = 0.5e-3
    nucleus_position: np.ndarray = field(
        default_factory=lambda: np.zeros(2, dtype=np.float64)
    )

    @classmethod
    def from_cortex(cls, cell_id: int, cortex: Cortex, **kwargs: Any) -> "Cell":
        """Build a :class:`Cell` from a freshly generated cortex.

        The centroid is recomputed from the cortex bead cloud (not the
        construction centre), so an elliptical cortex starts with its
        true geometric centre rather than the cortex's nominal
        ``cell_center``.
        """
        pts = cortex.bead_positions.reshape(-1, 2)
        centroid = pts.mean(axis=0)
        return cls(
            id=int(cell_id),
            cortex=cortex,
            center_position=centroid.astype(np.float64, copy=False),
            nucleus_position=centroid.astype(np.float64, copy=False),
            **kwargs,
        )

    def refresh_centroid(self) -> None:
        """Recompute ``center_position`` from the current cortex bead cloud."""
        pts = self.cortex.bead_positions.reshape(-1, 2)
        self.center_position = pts.mean(axis=0).astype(np.float64, copy=False)

    def compute_cortex_boundary_position(self, angle: float) -> np.ndarray:
        """Return the cortex bead nearest to the ray at ``angle`` from centre.

        Used by Worker D's junction code to find contact points
        between two cells. The angle is measured CCW from +x with
        respect to ``center_position``. The returned point is the
        bead whose direction-from-centre is closest to the requested
        angle (the cortex is a discrete bead cloud, so we report the
        bead position rather than a synthetic boundary point).
        """
        pts = self.cortex.bead_positions.reshape(-1, 2)
        rel = pts - self.center_position
        bead_angles = np.arctan2(rel[:, 1], rel[:, 0])
        # Wrap angular distance into [-π, +π].
        dtheta = np.angle(np.exp(1j * (bead_angles - angle)))
        idx = int(np.argmin(np.abs(dtheta)))
        return pts[idx].astype(np.float64, copy=False)
