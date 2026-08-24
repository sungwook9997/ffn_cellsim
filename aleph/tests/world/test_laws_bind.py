"""Binding ``laws/`` kernels to arena surfaces — the adapters, and the claim each one rests on.

No device.  Everything here is the HOST half of the bind: the hinge column permutation, the int32
narrowing, the range guard and the closure precondition.  The device half is
``aleph/scripts/arena_law_parity.py``, which is a CUDA run.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.laws.membrane_surface import build_membrane_hinges, membrane_bending_energy
from aleph.world import laws_bind
from aleph.world.arena import Kind, WorldArena
from aleph.world.laws_bind import HINGE_PERMUTATION, _as_int32_topology, upload_surface_topology
from aleph.world.surface import build_surface, icosphere


@pytest.fixture()
def surface():
    v, f = icosphere(2)
    arena = WorldArena(capacity={Kind.NODE: 400, Kind.FACE: 400, Kind.ANGLE4: 600})
    arena.claim("pad", Kind.NODE, 7)
    return arena, build_surface(arena, "membrane", vertices=v, faces=f, radius_um=5.0)


def test_self_check_passes() -> None:
    laws_bind._demo()


def test_hinge_permutation_is_energy_neutral(surface) -> None:
    """The claim the whole surface bind rests on: ``world``'s hinge tuple needs no orientation recovery.

    ``world.surface`` discards which face carries the directed half-edge; the kernel's docstring says
    ``T1`` must contain ``i -> j``.  Swapping ``i`` and ``j`` reverses BOTH discrete normals, so
    ``n1 . n2`` is unchanged — and the mesh is PERTURBED first, because on a resting icosphere every
    dihedral is nearly flat and that is exactly where a sign error hides.
    """
    _, s = surface
    rng = np.random.default_rng(7)
    pos = s.position + rng.normal(0.0, 0.05, s.position.shape)
    ours = np.asarray(s.hinge_idx, np.int64)[:, list(HINGE_PERMUTATION)] - s.nodes.lo
    theirs = build_membrane_hinges(s.face_idx) - s.nodes.lo
    e_ours = membrane_bending_energy(pos, ours, 1.0)
    e_theirs = membrane_bending_energy(pos, theirs, 1.0)
    assert e_theirs > 1.0, "a perturbed sphere must carry real bending energy or the test is vacuous"
    assert abs(e_ours - e_theirs) <= 1e-12 * abs(e_theirs)


def test_the_unpermuted_order_is_caught(surface) -> None:
    """The same check must SEPARATE a wrong order, or it certifies nothing about the right one."""
    _, s = surface
    rng = np.random.default_rng(7)
    pos = s.position + rng.normal(0.0, 0.05, s.position.shape)
    raw = np.asarray(s.hinge_idx, np.int64) - s.nodes.lo
    theirs = build_membrane_hinges(s.face_idx) - s.nodes.lo
    assert abs(membrane_bending_energy(pos, raw, 1.0)
               - membrane_bending_energy(pos, theirs, 1.0)) > 1.0


def test_topology_reaching_outside_its_claim_is_refused(surface) -> None:
    """One shared allocation makes a silent cross-population coupling possible; private arrays did not."""
    _, s = surface
    with pytest.raises(ValueError, match="outside its own claim"):
        _as_int32_topology("face_idx", np.asarray(s.face_idx, np.int64) - s.nodes.lo, s)
    with pytest.raises(ValueError, match="outside its own claim"):
        _as_int32_topology("face_idx", np.asarray(s.face_idx, np.int64) + s.nodes.count, s)


def test_bookkeeping_only_arena_refuses_upload(surface) -> None:
    """``device=None`` is the CPU-importable mode; there is no CPU simulation path and never may be."""
    arena, s = surface
    with pytest.raises(RuntimeError, match="no CPU simulation path"):
        upload_surface_topology(arena, s)
