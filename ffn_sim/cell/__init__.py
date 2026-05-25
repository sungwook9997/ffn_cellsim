"""H.3 Cell composition module (v2 — replaces deleted v1 acs_kb/cell/cell.py).

Public API:
- :class:`Cell` — composition over cortex + ERM + crosslinkers (+ future
  myosin / lamellipodium / FA slots).
- :class:`CellBuildOptions` — per-build feature flags DTO.
"""

from ffn_sim.cell.cell import Cell, CellBuildOptions

__all__ = ["Cell", "CellBuildOptions"]
