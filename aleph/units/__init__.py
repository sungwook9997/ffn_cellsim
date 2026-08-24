"""Cited acceptance bands. `ALEPH-PORT-4001`.

This package holds the provenance side of numbers: a band, its unit, and where the bound comes
from. It is **not** the unit system — that is `aleph/laws/units.py`, the pN·µm·s
nondimensionalisation, and it stays there.

Layer position: the bottom. Nothing here may import from `ff/`, `ac/`, or anything above.
"""

from __future__ import annotations

from aleph.units.band import LiteratureBand, OutOfBandError, band_guard
from aleph.units.provenance import PROVENANCE

__all__ = [
    "PROVENANCE",
    "LiteratureBand",
    "OutOfBandError",
    "band_guard",
]
