"""The arena engine: one fixed-capacity world, claimed as contiguous ID ranges.

A separate implementation line from ``aleph/engine/``, run BESIDE the incumbent and never replacing
it, with node-by-node force parity as the A/B. Nothing here computes physics yet; :mod:`.arena` is the
allocation and bookkeeping layer everything else will claim out of.
"""

from aleph.world.arena import Claim, Kind, WorldArena
from aleph.world.bond import BondCount, BondFamily, SourceClass, build_bond_family
from aleph.world.floorplan import render_floorplan, write_floorplan
from aleph.world.strand import Strand, build_strand
from aleph.world.surface import Surface, build_surface, enclosed_volume, icosphere
from aleph.world.viewer import render_cell_view, write_cell_view

__all__ = [
    "Claim", "Kind", "WorldArena",
    "Strand", "build_strand",
    "BondCount", "BondFamily", "SourceClass", "build_bond_family",
    "render_floorplan", "write_floorplan",
    "Surface", "build_surface", "icosphere", "enclosed_volume",
    "render_cell_view", "write_cell_view",
]
