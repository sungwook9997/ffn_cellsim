"""Cell-cell junction state schema.

Canonical Phase 1 ``JunctionState`` per
``docs/v2_phase1_plan_consolidated.md`` §5.6. A junction lives between
exactly two distinct cells; cells in contact remain separate polygons
and never merge geometry. Junction dynamics (creation, aging,
maturation, contact-edge protrusion suppression, contact inhibition
signal) are explicitly out of P0 scope.

Sanity Gate scope: non-physics schema module. Full 6-item gate N/A;
this module owns the boundary-case checks (distinct cell ids, no
self-junction, contact_length non-negative, age non-negative, bounded
maturity / scalar proxies, valid state literal) and the
measurement-protocol consistency that the cell-id pair is canonicalized
(``cell_id_a`` < ``cell_id_b`` lexicographically) so one junction is
not double-counted in symmetric neighbor graphs.

Magic-Number Block: no tunable numeric. N/A.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal, Optional

JunctionStateLabel = Literal[
    "free",
    "contacting",
    "nascent",
    "cortex_coupled",
    "remodeling",
    "separating",
]


@dataclass(frozen=True)
class JunctionState:
    """Schema for one cell-cell junction between two distinct cells.

    Cell-id pair is stored canonicalized (``cell_id_a < cell_id_b``).
    Use :meth:`from_cell_pair` to construct from an unordered pair.
    """

    junction_id: str
    cell_id_a: str
    cell_id_b: str
    contact_length_um: float
    age_s: float
    maturity: float
    state: JunctionStateLabel = "contacting"
    contact_edge_indices_a: tuple[int, ...] = field(default_factory=tuple)
    contact_edge_indices_b: tuple[int, ...] = field(default_factory=tuple)
    cadherin_proxy: Optional[float] = None
    tension_proxy_nN: Optional[float] = None
    contact_inhibition_signal: Optional[float] = None

    @classmethod
    def from_cell_pair(
        cls,
        junction_id: str,
        cell_id_pair: tuple[str, str],
        contact_length_um: float,
        age_s: float,
        maturity: float,
        **kwargs,
    ) -> "JunctionState":
        if len(cell_id_pair) != 2:
            raise ValueError("cell_id_pair must contain exactly 2 ids")
        a, b = sorted(cell_id_pair)
        return cls(
            junction_id=junction_id,
            cell_id_a=a,
            cell_id_b=b,
            contact_length_um=contact_length_um,
            age_s=age_s,
            maturity=maturity,
            **kwargs,
        )

    def validate(self) -> None:
        if not self.junction_id.strip():
            raise ValueError("junction_id must be non-empty")
        if not self.cell_id_a.strip() or not self.cell_id_b.strip():
            raise ValueError("cell_id_a and cell_id_b must be non-empty")
        if self.cell_id_a == self.cell_id_b:
            raise ValueError("a junction cannot link a cell to itself")
        if self.cell_id_a >= self.cell_id_b:
            raise ValueError(
                "cell_id_a must be lexicographically less than cell_id_b "
                "(canonicalize via JunctionState.from_cell_pair)"
            )
        if not math.isfinite(self.contact_length_um) or self.contact_length_um < 0.0:
            raise ValueError("contact_length_um must be finite and non-negative")
        if not math.isfinite(self.age_s) or self.age_s < 0.0:
            raise ValueError("age_s must be finite and non-negative")
        if not math.isfinite(self.maturity) or not (0.0 <= self.maturity <= 1.0):
            raise ValueError("maturity must be in [0, 1]")
        if self.state not in {
            "free",
            "contacting",
            "nascent",
            "cortex_coupled",
            "remodeling",
            "separating",
        }:
            raise ValueError(
                "state must be free, contacting, nascent, cortex_coupled, "
                "remodeling, or separating"
            )

        for label, indices in (
            ("contact_edge_indices_a", self.contact_edge_indices_a),
            ("contact_edge_indices_b", self.contact_edge_indices_b),
        ):
            for idx in indices:
                if not isinstance(idx, int) or isinstance(idx, bool):
                    raise ValueError(f"{label} must contain integers, got {idx!r}")
                if idx < 0:
                    raise ValueError(f"{label} indices must be non-negative")
            if len(set(indices)) != len(indices):
                raise ValueError(f"{label} must contain unique indices")

        for label, value in (
            ("cadherin_proxy", self.cadherin_proxy),
            ("contact_inhibition_signal", self.contact_inhibition_signal),
        ):
            if value is None:
                continue
            if not math.isfinite(value) or not (0.0 <= value <= 1.0):
                raise ValueError(f"{label} must be in [0, 1] when provided")

        if self.tension_proxy_nN is not None:
            if not math.isfinite(self.tension_proxy_nN):
                raise ValueError("tension_proxy_nN must be finite when provided")
