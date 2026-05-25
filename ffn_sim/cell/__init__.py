"""H.3 Cell composition module (v2 — replaces deleted v1 acs_kb/cell/cell.py).

Public API:
- :class:`Cell` — composition over cortex + ERM + crosslinkers + myosin
  (+ future lamellipodium / FA slots).
- :class:`CellBuildOptions` — per-build feature flags DTO.
- :func:`build_cortex_full_simulation` — unified cortex + xlinks +
  myosin builder, used internally by :meth:`Cell.build` when
  ``with_myosin=True``.
"""

from ffn_sim.cell.cell import (
    Cell,
    CellBuildOptions,
    build_cortex_full_simulation,
)

__all__ = ["Cell", "CellBuildOptions", "build_cortex_full_simulation"]
