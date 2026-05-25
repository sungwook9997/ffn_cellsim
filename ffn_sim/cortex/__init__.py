"""H.3 Cortex multi-filament network module.

Public API:
- :func:`resolve_h3_derived` — resolve config YAML to ``ResolvedH3``.
- :func:`generate_cortex_topology` — sphere-shell topology generator.
- :func:`generate_static_xl_bonds` — optional KU-3.19 crosslinker bonds.
- :func:`build_cortex_state` — HOOMD GSD frame builder.
- :func:`build_cortex_simulation` — full HOOMD Simulation with BAOAB.
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

__all__ = [
    "CortexTopology",
    "CrosslinkerBonds",
    "ResolvedH3",
    "build_cortex_simulation",
    "build_cortex_state",
    "generate_cortex_topology",
    "generate_static_xl_bonds",
    "resolve_h3_derived",
    "xl_bin_rest_lengths",
    "xl_bin_type_names",
]
