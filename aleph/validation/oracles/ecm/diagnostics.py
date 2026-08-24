"""Diagnostic helpers for Unit 1.2 dynamics validation."""

from __future__ import annotations

import numpy as np

from aleph.validation.oracles.ecm.cross_links import CrossLink


def link_extension_energies(
    bead_positions: np.ndarray,
    cross_links: list[CrossLink],
    box_size: float,
) -> np.ndarray:
    """Return per-link harmonic energy ½ k_xl (|d|−r_0)² in joules."""
    if not cross_links:
        return np.zeros(0)
    fa = np.fromiter((xl.fiber_a for xl in cross_links), dtype=np.int64,
                     count=len(cross_links))
    ba = np.fromiter((xl.bead_a for xl in cross_links), dtype=np.int64,
                     count=len(cross_links))
    fb = np.fromiter((xl.fiber_b for xl in cross_links), dtype=np.int64,
                     count=len(cross_links))
    bb = np.fromiter((xl.bead_b for xl in cross_links), dtype=np.int64,
                     count=len(cross_links))
    r0 = np.fromiter((xl.rest_length for xl in cross_links), dtype=np.float64,
                     count=len(cross_links))
    k = np.fromiter((xl.stiffness for xl in cross_links), dtype=np.float64,
                    count=len(cross_links))
    d = bead_positions[fb, bb, :] - bead_positions[fa, ba, :]
    d -= box_size * np.round(d / box_size)
    dist = np.linalg.norm(d, axis=1)
    return 0.5 * k * (dist - r0) ** 2
