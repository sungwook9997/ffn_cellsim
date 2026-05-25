"""H.3 Cortex multi-filament network module.

Public API:
- :func:`resolve_h3_derived` — resolve cortex config YAML to ``ResolvedH3``.
- :func:`generate_cortex_topology` — sphere-shell topology generator.
- :func:`generate_static_xl_bonds` — optional static KU-3.19 crosslinker bonds.
- :func:`build_cortex_state` — HOOMD GSD frame builder.
- :func:`build_cortex_simulation` — full HOOMD Simulation with BAOAB.
- :func:`resolve_crosslinkers` — D2 Bell-Evans dynamic xlink config.
- :class:`XlinkBondUpdater` — D2 Bell-Evans dynamic xlink_head ↔ actin_cortex.
- :func:`build_cortex_xlink_simulation` — cortex + dynamic xlink wired sim.
"""

from ffn_sim.cortex.cortex import (
    CortexTopology,
    CrosslinkerBonds,
    ResolvedH3,
    build_cortex_simulation,
    build_cortex_state,
    generate_cortex_topology,
    generate_static_xl_bonds,
    resolve_h3_derived,
    xl_bin_rest_lengths,
    xl_bin_type_names,
)
from ffn_sim.cortex.crosslinkers import (
    ResolvedCrosslinkers,
    XlinkBondUpdater,
    XlinkLayout,
    build_cortex_xlink_simulation,
    extend_cortex_state_with_xlinks,
    generate_xlink_layout,
    make_xlink_updater,
    resolve_crosslinkers,
    xlink_attach_bin_names,
    xlink_attach_bin_rest_lengths,
)

__all__ = [
    "CortexTopology",
    "CrosslinkerBonds",
    "ResolvedCrosslinkers",
    "ResolvedH3",
    "XlinkBondUpdater",
    "XlinkLayout",
    "build_cortex_simulation",
    "build_cortex_state",
    "build_cortex_xlink_simulation",
    "extend_cortex_state_with_xlinks",
    "generate_cortex_topology",
    "generate_static_xl_bonds",
    "generate_xlink_layout",
    "make_xlink_updater",
    "resolve_crosslinkers",
    "resolve_h3_derived",
    "xl_bin_rest_lengths",
    "xl_bin_type_names",
    "xlink_attach_bin_names",
    "xlink_attach_bin_rest_lengths",
]
