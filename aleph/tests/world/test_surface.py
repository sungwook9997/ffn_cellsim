"""Surfaces — the two invariants a closed mesh must carry, and the volume that must not be analytic."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.world.arena import Kind, WorldArena
from aleph.world.surface import build_surface, enclosed_volume, icosphere
from aleph.world import surface as surf_mod


@pytest.fixture()
def arena() -> WorldArena:
    return WorldArena(capacity={Kind.NODE: 20_000, Kind.FACE: 40_000, Kind.ANGLE4: 60_000})


def test_self_check_passes() -> None:
    surf_mod._demo()


@pytest.mark.parametrize("n,verts,faces", [(0, 12, 20), (1, 42, 80), (2, 162, 320), (3, 642, 1280)])
def test_icosphere_counts_follow_the_subdivision_law(n: int, verts: int, faces: int) -> None:
    v, f = icosphere(n)
    assert v.shape == (verts, 3) and f.shape == (faces, 3)
    assert np.allclose(np.linalg.norm(v, axis=1), 1.0), "vertices lie on the unit sphere"


def test_closed_mesh_has_three_halves_F_hinges(arena: WorldArena) -> None:
    """The hinge count, not the face count, sets the number of bending terms."""
    v, f = icosphere(2)
    s = build_surface(arena, "membrane", vertices=v, faces=f, radius_um=7.5)
    assert s.closed
    assert s.hinges.count == 3 * s.faces.count // 2 == 480


def test_volume_is_discrete_and_below_the_analytic_sphere(arena: WorldArena) -> None:
    """Substituting the analytic value produced a phantom force at t0 once already: a mesh inscribed in
    a sphere encloses LESS than the sphere."""
    v, f = icosphere(2)
    s = build_surface(arena, "membrane", vertices=v, faces=f, radius_um=7.5)
    analytic = 4.0 / 3.0 * np.pi * 7.5 ** 3
    assert 0.0 < s.volume0_um3 < analytic
    assert (analytic - s.volume0_um3) / analytic == pytest.approx(0.0339, abs=5e-4)


def test_refinement_shrinks_the_discrete_gap(arena: WorldArena) -> None:
    """Which is exactly why the substitution looks like noise instead of an error."""
    gaps = []
    for n in (1, 2, 3):
        v, f = icosphere(n)
        s = build_surface(arena, f"probe{n}", vertices=v, faces=f, radius_um=1.0)
        gaps.append((4.0 / 3.0 * np.pi - s.volume0_um3) / (4.0 / 3.0 * np.pi))
    assert gaps[0] > gaps[1] > gaps[2] > 0.0


def test_inward_winding_is_refused(arena: WorldArena) -> None:
    """An inverted face flips the sign of its own pressure traction."""
    v, f = icosphere(1)
    with pytest.raises(ValueError, match="wound INWARD"):
        build_surface(arena, "flipped", vertices=v, faces=f[:, ::-1], radius_um=1.0)


def test_enclosed_volume_sign_tracks_winding() -> None:
    v, f = icosphere(1)
    assert enclosed_volume(v, f) > 0.0
    assert enclosed_volume(v, f[:, ::-1]) == pytest.approx(-enclosed_volume(v, f))


def test_radius_refuses_to_default(arena: WorldArena) -> None:
    v, f = icosphere(1)
    with pytest.raises(ValueError, match="no default"):
        build_surface(arena, "sized", vertices=v, faces=f)


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan")])
def test_radius_must_be_finite_and_positive(arena: WorldArena, bad: float) -> None:
    v, f = icosphere(1)
    with pytest.raises(ValueError, match="finite and positive"):
        build_surface(arena, "sized", vertices=v, faces=f, radius_um=bad)


def test_area_scales_with_radius_squared(arena: WorldArena) -> None:
    v, f = icosphere(2)
    one = build_surface(arena, "a", vertices=v, faces=f, radius_um=1.0)
    ten = build_surface(arena, "b", vertices=v, faces=f, radius_um=10.0)
    assert ten.area0_total_um2 / one.area0_total_um2 == pytest.approx(100.0)


def test_too_few_faces_cannot_enclose(arena: WorldArena) -> None:
    with pytest.raises(ValueError, match="at least 4 faces"):
        build_surface(arena, "flat", vertices=np.eye(3), faces=np.array([[0, 1, 2]]), radius_um=1.0)


def test_out_of_range_face_index_is_refused(arena: WorldArena) -> None:
    v, f = icosphere(1)
    bad = f.copy()
    bad[0, 0] = v.shape[0] + 3
    with pytest.raises(ValueError, match="face indexes vertex"):
        build_surface(arena, "bad", vertices=v, faces=bad, radius_um=1.0)


def test_non_manifold_edge_is_refused(arena: WorldArena) -> None:
    """A dihedral over an edge used by three faces would be ambiguous, not merely wrong."""
    v, f = icosphere(1)
    extra = np.vstack([f, f[0][None, :]])
    with pytest.raises(ValueError, match="non-manifold"):
        build_surface(arena, "nm", vertices=v, faces=extra, radius_um=1.0)


def test_topology_is_global_and_the_arena_stays_partitioned(arena: WorldArena) -> None:
    v, f = icosphere(1)
    build_surface(arena, "membrane", vertices=v, faces=f, radius_um=7.5)
    s2 = build_surface(arena, "nucleus", vertices=v, faces=f, radius_um=3.0)
    assert int(s2.face_idx.min()) >= s2.nodes.lo
    assert int(s2.hinge_idx.max()) < s2.nodes.hi
    arena.assert_partitioned()
    assert arena.census()["populations"]["nucleus"]["face"] == f.shape[0] == 80


def test_centre_translates_without_changing_area_or_volume(arena: WorldArena) -> None:
    v, f = icosphere(2)
    at_origin = build_surface(arena, "a", vertices=v, faces=f, radius_um=3.0)
    offset = build_surface(arena, "b", vertices=v, faces=f, radius_um=3.0, centre=(10.0, -4.0, 2.5))
    assert offset.volume0_um3 == pytest.approx(at_origin.volume0_um3)
    assert offset.area0_total_um2 == pytest.approx(at_origin.area0_total_um2)
