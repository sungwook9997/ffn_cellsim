"""H.7 — generate_lamellipodium_layout geometry dispatch (flat / basal_ring /
polarized_patch).

The per-geometry placement physics is unit-tested in test_lamellipodium_basal_ring
and test_lamellipodium_polarized_patch. THIS file pins the integration: the
`geometry=` dispatch on generate_lamellipodium_layout routes to the right
builder, returns a drop-in LamellipodiumLayout, populates `mother_tangents` for
the single-cell geometries (None for flat), and leaves the flat path unchanged.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.cell.lamellipodium import (
    LamellipodiumLayout,
    generate_lamellipodium_layout,
    resolve_h5_lamellipodium,
)

_H5_CFG = Path(__file__).resolve().parents[1] / "configs" / "phase1_h5.yaml"
_R_CELL = 7.5e-6
_L_BOX = 3.0 * _R_CELL  # box_factor (3) * R_cell, matches the cortex box


@pytest.fixture(scope="module")
def p_h5():
    cfg = yaml.safe_load(open(_H5_CFG))
    # phase1_h5.yaml's Y_max / wave_area were authored for a 30 um box; for the
    # R_cell=7.5 um cell (L_box=22.5 um) let them auto-derive from the box
    # (Y_max=0.45*L_box, wave_area=L_box^2) so the flat resolve is valid. This
    # box-vs-Y_max mismatch is exactly the geometry problem the H.7 single-cell
    # geometries fix; here we only need a valid ResolvedH5 for the dispatch test.
    cfg["lamellipodium"].pop("Y_max", None)
    cfg["lamellipodium"].pop("wave_area", None)
    return resolve_h5_lamellipodium(cfg, L_box=_L_BOX, dt=1.0e-8)


def test_flat_plane_default_unchanged(p_h5):
    """Default geometry stays flat_plane: WAVE on the y=Y_max plane, no
    mother_tangents (seed loops then use the legacy -y)."""
    lay = generate_lamellipodium_layout(p_h5, wave_tag_start=0)
    assert isinstance(lay, LamellipodiumLayout)
    assert lay.mother_tangents is None
    # flat WAVE plane: all y == Y_max
    np.testing.assert_allclose(lay.wave_positions[:, 1], p_h5.Y_max)


def test_basal_ring_dispatch_sets_radial_tangents(p_h5):
    """geometry='basal_ring' routes to the ring builder, returns a drop-in
    layout with outward-radial basal tangents."""
    lay = generate_lamellipodium_layout(
        p_h5, wave_tag_start=0, geometry="basal_ring", R_cell=_R_CELL
    )
    assert isinstance(lay, LamellipodiumLayout)
    assert lay.mother_tangents is not None
    assert lay.mother_tangents.shape == (p_h5.n_WAVE, 3)
    assert np.isfinite(lay.mother_tangents).all()
    # outward-radial in the basal plane: ~unit, small z-component
    norms = np.linalg.norm(lay.mother_tangents, axis=1)
    np.testing.assert_allclose(norms, 1.0, atol=1e-6)
    assert np.abs(lay.mother_tangents[:, 2]).max() < 0.2


def test_polarized_patch_dispatch_sets_forward_tangents(p_h5):
    """geometry='polarized_patch' routes to the patch builder; all mother
    tangents point forward along the (default +x) polarization."""
    lay = generate_lamellipodium_layout(
        p_h5, wave_tag_start=0, geometry="polarized_patch", R_cell=_R_CELL
    )
    assert isinstance(lay, LamellipodiumLayout)
    assert lay.mother_tangents is not None
    assert lay.mother_tangents.shape == (p_h5.n_WAVE, 3)
    # forward growth: in-plane component aligns with +x (mean dot > 0)
    inplane = lay.mother_tangents[:, :2]
    inplane = inplane / np.linalg.norm(inplane, axis=1, keepdims=True)
    assert float(inplane[:, 0].mean()) > 0.5


def test_non_flat_requires_R_cell(p_h5):
    for geom in ("basal_ring", "polarized_patch"):
        with pytest.raises(ValueError, match="R_cell"):
            generate_lamellipodium_layout(p_h5, wave_tag_start=0, geometry=geom)


def test_unknown_geometry_raises(p_h5):
    with pytest.raises(ValueError, match="unknown geometry"):
        generate_lamellipodium_layout(
            p_h5, wave_tag_start=0, geometry="banana", R_cell=_R_CELL
        )
