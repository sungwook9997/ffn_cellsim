"""Bridge value objects (week-5 interface contract).

The :class:`FocalAdhesion` dataclass is the shared data contract between
Worker B (this Unit, motor-clutch) and Worker C (Unit 3.x, cell motility):
Worker C imports it to populate ``Cell.focal_adhesions``. Field names and
types here are **frozen at week 5**; any change must be coordinated with
Worker C via Notion alert before merge.

All quantities are SI (metres, seconds, Newtons). 1-D loading axis is
assumed for Phase 1 Unit 2.1; the 2-D extension enters at Unit 2.2.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(slots=True)
class FocalAdhesion:
    """Single focal adhesion state (Phase 1 Unit 2.1).

    Attributes
    ----------
    position : np.ndarray, shape (2,)
        FA centre in the lab frame, in metres.
    n_clutches_total : int
        Total integrin–ligand bonds carried by this FA (engaged + free).
        Default 50 follows KU-2.18.
    clutches_engaged : np.ndarray, shape (n_clutches_total,), dtype=bool
        ``True`` for clutches currently bound to the substrate.
    clutch_forces : np.ndarray, shape (n_clutches_total,), dtype=float64
        Per-clutch axial force in Newtons; disengaged entries are 0.
    actin_position : float
        1-D actin-bundle position along the loading axis, in metres.
        Increases as myosin pulls retrograde flow into the FA.
    area : float
        FA projected area in m² (1 μm² Phase 1 default, KU-2.2).
    vinculin_count : int
        Bound vinculin molecules. **Held at 0 in Unit 2.1**; populated
        from Unit 2.2 onward (KU-2.7).
    talin_unfolded_domains : int
        Number of unfolded talin R-domains (0–13). **Held at 0 in Unit
        2.1**; populated from Unit 2.2 onward (KU-2.6).
    age : float
        Seconds since this FA was instantiated.
    """

    position: np.ndarray
    n_clutches_total: int
    clutches_engaged: np.ndarray
    clutch_forces: np.ndarray
    actin_position: float = 0.0
    area: float = 1.0e-12
    vinculin_count: int = 0
    talin_unfolded_domains: int = 0
    age: float = 0.0


def make_focal_adhesion(
    position: np.ndarray,
    n_clutches_total: int = 50,
    area: float = 1.0e-12,
) -> FocalAdhesion:
    """Construct a fresh FA with all clutches disengaged (KU-2.18 defaults)."""
    pos = np.asarray(position, dtype=np.float64).reshape(2)
    return FocalAdhesion(
        position=pos,
        n_clutches_total=int(n_clutches_total),
        clutches_engaged=np.zeros(n_clutches_total, dtype=bool),
        clutch_forces=np.zeros(n_clutches_total, dtype=np.float64),
        actin_position=0.0,
        area=float(area),
        vinculin_count=0,
        talin_unfolded_domains=0,
        age=0.0,
    )
