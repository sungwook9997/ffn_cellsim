"""Self-test of the §1.4 nucleus relative-no-flux boundary provider (pure NumPy — no Warp/CUDA).

The frozen fluid-track interface: both providers satisfy the ``NucleusBoundaryProvider`` Protocol; the
deformable winding-number occupancy matches the analytic oblate ellipsoid; the boundary faces carry
outward unit normals; and ``update`` makes the mask track the MOVING mesh (not a static R_nuc).
"""
from __future__ import annotations

import numpy as np
import pytest

from aleph.components.nucleus.geometry import build_oblate_mesh
from aleph.components.nucleus.mask_provider import (
    DeformableNucleusMaskProvider,
    NucleusBoundaryProvider,
    StaticSphereMaskProvider,
    inside_ellipsoid,
    solid_angle_winding_number,
)

GRID = np.mgrid[-5:5:12j, -5:5:12j, -5:5:12j].reshape(3, -1).T


def test_both_providers_satisfy_protocol() -> None:
    """The static placeholder (I1a) and the deformable mesh (I2) are interchangeable at the contract."""
    assert isinstance(StaticSphereMaskProvider((0, 0, 0), 3.0), NucleusBoundaryProvider)
    v, f = build_oblate_mesh(3.0, 2.0, 3)
    assert isinstance(DeformableNucleusMaskProvider(v, f), NucleusBoundaryProvider)


def test_static_sphere_matches_analytic() -> None:
    """The placeholder mask equals ‖x−c‖≤R exactly (mirrors what fluid-spine ships)."""
    sp = StaticSphereMaskProvider((0.0, 0.0, 0.0), 3.0)
    ref = np.linalg.norm(GRID, axis=1) <= 3.0
    assert np.array_equal(sp.inside_mask(GRID), ref)


def test_deformable_winding_matches_oblate_ellipsoid() -> None:
    """Winding-number occupancy of the oblate MESH == the analytic ellipsoid occupancy (robust
    point-in-closed-mesh; exact for the smooth baseline)."""
    v, f = build_oblate_mesh(3.0, 2.0, 3)                        # a=3, c=1.5
    prov = DeformableNucleusMaskProvider(v, f)
    m_mesh = prov.inside_mask(GRID)
    m_ellip = inside_ellipsoid(GRID, (0, 0, 0), 3.0, 1.5)
    assert np.array_equal(m_mesh, m_ellip)


def test_point_in_mesh_poles_and_centre() -> None:
    """Occupancy is correct at the centre, far field, and just inside/outside the flattened pole."""
    v, f = build_oblate_mesh(3.0, 2.0, 3)                        # polar semi-axis c=1.5
    prov = DeformableNucleusMaskProvider(v, f)
    assert prov.inside_mask(np.array([[0.0, 0.0, 0.0]]))[0]
    assert not prov.inside_mask(np.array([[10.0, 0.0, 0.0]]))[0]
    assert prov.inside_mask(np.array([[0.0, 0.0, 1.4]]))[0]      # inside c=1.5
    assert not prov.inside_mask(np.array([[0.0, 0.0, 1.6]]))[0]  # outside c=1.5


def test_boundary_normals_outward_and_unit() -> None:
    """Boundary faces carry outward unit normals (x·n̂>0 about the centred nucleus) for the flux BC."""
    v, f = build_oblate_mesh(3.0, 2.0, 3)
    cen, nrm = DeformableNucleusMaskProvider(v, f).boundary()
    assert np.allclose(np.linalg.norm(nrm, axis=1), 1.0)
    assert np.all(np.einsum("ij,ij->i", cen, nrm) > 0.0)        # outward


def test_update_tracks_moving_mesh() -> None:
    """The relative-no-flux mask follows the MOVING envelope — a point outside the resting mesh becomes
    inside after the mesh inflates (the whole point of a deformable, not static, boundary)."""
    v, f = build_oblate_mesh(3.0, 2.0, 3)
    prov = DeformableNucleusMaskProvider(v, f)
    pt = np.array([[0.0, 0.0, 2.0]])                             # outside resting pole c=1.5
    assert not prov.inside_mask(pt)[0]
    prov.update(v * 1.5)                                         # inflate → pole now c=2.25
    assert prov.inside_mask(pt)[0]


def test_inward_wound_mesh_rejected() -> None:
    """A negative-volume (inward-wound) mesh is rejected — normals must point outward for the BC."""
    v, f = build_oblate_mesh(3.0, 1.0, 2)
    with pytest.raises(ValueError):
        DeformableNucleusMaskProvider(v, f[:, ::-1])            # reversed winding → inward normals
