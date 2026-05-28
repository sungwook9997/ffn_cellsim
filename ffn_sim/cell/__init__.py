"""H.3 + H.5 Cell composition module (v2 — replaces deleted v1 acs_kb/cell/cell.py).

Public API:
- :class:`Cell` — composition over cortex + ERM + crosslinkers + myosin
  + lamellipodium (+ future FA slot).
- :class:`CellBuildOptions` — per-build feature flags DTO.
- :func:`build_cortex_full_simulation` — unified cortex + xlinks +
  myosin builder, used internally by :meth:`Cell.build` when
  ``with_myosin=True``.
- :class:`ResolvedH5` / :func:`resolve_h5_lamellipodium` — H.5
  lamellipodium config.
- :func:`build_lamellipodium_simulation` — lamellipodium-only sim
  (WAVE plane + Arp2/3 branching + elongation + capping); for KU-5.x
  validation gates.
- :func:`extend_cortex_snapshot_with_lamellipodium` /
  :func:`attach_lamellipodium_to_simulation` — composition helpers used
  by :meth:`Cell.build` when ``with_lamellipodium=True``.
"""

from ffn_sim.cell.cell import (
    Cell,
    CellBuildOptions,
    build_cortex_full_simulation,
)
from ffn_sim.cell.lamellipodium import (
    ArpBranchingUpdater,
    BarbedEndElongationUpdater,
    CappingUpdater,
    LamellipodiumLayout,
    LamellipodiumState,
    ResolvedH5,
    WaveMembranePin,
    attach_lamellipodium_to_simulation,
    build_lamellipodium_simulation,
    extend_cortex_snapshot_with_lamellipodium,
    generate_lamellipodium_layout,
    resolve_h5_lamellipodium,
)

__all__ = [
    "ArpBranchingUpdater",
    "BarbedEndElongationUpdater",
    "CappingUpdater",
    "Cell",
    "CellBuildOptions",
    "LamellipodiumLayout",
    "LamellipodiumState",
    "ResolvedH5",
    "WaveMembranePin",
    "attach_lamellipodium_to_simulation",
    "build_cortex_full_simulation",
    "build_lamellipodium_simulation",
    "extend_cortex_snapshot_with_lamellipodium",
    "generate_lamellipodium_layout",
    "resolve_h5_lamellipodium",
]
