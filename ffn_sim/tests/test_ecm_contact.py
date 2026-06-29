"""Sanity-gate tests for ``bridge/ecm_contact.py`` (cell↔ECM contact manifold).

Pure-numpy core (``bin_ecm_contact``), the tag-ordered sim reader, and the
FA-traction cross-consistency gate (FA-engaged patches ⊆ contact patches).
"""

from __future__ import annotations

from contextlib import contextmanager

import numpy as np
import pytest

from ffn_sim.archive.hoomd_legacy.bridge.ecm_contact import (
    TYPE_LIGAND,
    assert_traction_within_contact,
    bin_ecm_contact,
    measure_ecm_contact_manifold,
)
from ffn_sim.common.surface_manifold import SurfaceManifold

R_CELL = 7.5e-6


@pytest.fixture(scope="module")
def manifold() -> SurfaceManifold:
    return SurfaceManifold.icosphere(subdivisions=3, radius=R_CELL)


# ---------------------------------------------------------------------------
# Pure-core binning
# ---------------------------------------------------------------------------
class TestBinCore:
    def test_south_cap_ligands_in_contact(self, manifold):
        # Place one ligand just OUTWARD (toward substrate) of each south-cap
        # patch centroid: gap = 10 nm ≪ contact_gap → those patches in contact.
        south = manifold.tri_centroids[:, 2] < -0.5 * R_CELL
        south_ids = np.flatnonzero(south)
        gap = 10.0e-9
        lig = manifold.tri_centroids[south_ids] + gap * manifold.tri_normals[south_ids]
        out = bin_ecm_contact(manifold, lig, contact_gap=50.0e-9)
        # Every south-cap patch we seeded is in contact.
        assert np.all(out["contact_mask"][south_ids])
        # No NORTH patch is in contact (no ligands up there).
        north = manifold.tri_centroids[:, 2] > 0.5 * R_CELL
        assert not np.any(out["contact_mask"][np.flatnonzero(north)])
        # Partition: every ligand counted once.
        assert out["partition"]["passed"]
        assert out["n_ligands"] == lig.shape[0]
        assert out["patch_n_ligands"].sum() == lig.shape[0]

    def test_far_ligands_not_in_contact(self, manifold):
        # Ligands a full patch-width away from any centroid (gap ≫ contact_gap)
        # are homed but NOT in contact.
        south_id = int(np.argmin(manifold.tri_centroids[:, 2]))
        far = manifold.tri_centroids[south_id] + 2.0e-6 * manifold.tri_normals[south_id]
        out = bin_ecm_contact(manifold, far[None, :], contact_gap=50.0e-9)
        assert out["n_ligands"] == 1
        assert out["n_ligands_in_contact"] == 0
        assert out["n_contact_patches"] == 0
        assert out["contact_area_m2"] == 0.0

    def test_contact_gap_boundary_inclusive(self, manifold):
        sid = int(np.argmin(manifold.tri_centroids[:, 2]))
        gap = 30.0e-9
        lig = manifold.tri_centroids[sid] + gap * manifold.tri_normals[sid]
        out = bin_ecm_contact(manifold, lig[None, :], contact_gap=gap)
        # gap == contact_gap → in contact (≤).
        assert out["ligand_in_contact"][0]
        assert out["contact_mask"][sid]

    def test_signed_normal_offset_sense(self, manifold):
        # A ligand placed +n̂ (outward) of the centroid → positive normal offset.
        sid = int(np.argmin(manifold.tri_centroids[:, 2]))
        lig = manifold.tri_centroids[sid] + 12.0e-9 * manifold.tri_normals[sid]
        out = bin_ecm_contact(manifold, lig[None, :], contact_gap=50.0e-9)
        assert out["patch_normal_offset_m"][sid] == pytest.approx(12.0e-9, rel=1e-3)

    def test_contact_area_equals_sum_of_patch_areas(self, manifold):
        south = manifold.tri_centroids[:, 2] < -0.5 * R_CELL
        ids = np.flatnonzero(south)
        lig = manifold.tri_centroids[ids] + 10.0e-9 * manifold.tri_normals[ids]
        out = bin_ecm_contact(manifold, lig, contact_gap=50.0e-9)
        expect = float(np.sum(manifold.tri_areas[out["contact_mask"]]))
        assert out["contact_area_m2"] == pytest.approx(expect, rel=1e-12)
        assert 0.0 < out["contact_area_fraction"] < 1.0

    def test_min_gap_is_nearest_ligand(self, manifold):
        # Two ligands on the same patch: patch_min_gap is the NEARER one.
        sid = int(np.argmin(manifold.tri_centroids[:, 2]))
        n = manifold.tri_normals[sid]
        c = manifold.tri_centroids[sid]
        lig = np.stack([c + 40.0e-9 * n, c + 8.0e-9 * n])
        out = bin_ecm_contact(manifold, lig, contact_gap=50.0e-9)
        assert out["patch_n_ligands"][sid] == 2
        assert out["patch_min_gap_m"][sid] == pytest.approx(8.0e-9, rel=1e-3)

    def test_empty_ligands(self, manifold):
        out = bin_ecm_contact(
            manifold, np.empty((0, 3)), contact_gap=50.0e-9
        )
        assert out["n_ligands"] == 0
        assert out["n_contact_patches"] == 0
        assert out["contact_area_m2"] == 0.0
        assert out["partition"]["passed"]
        assert out["patch_min_gap_m"].shape == (manifold.n_tri,)
        assert np.all(np.isinf(out["patch_min_gap_m"]))

    def test_bad_contact_gap_raises(self, manifold):
        with pytest.raises(ValueError):
            bin_ecm_contact(manifold, np.empty((0, 3)), contact_gap=0.0)
        with pytest.raises(ValueError):
            bin_ecm_contact(manifold, np.empty((0, 3)), contact_gap=-1.0)


# ---------------------------------------------------------------------------
# Tag-ordered sim reader (fake HOOMD sim)
# ---------------------------------------------------------------------------
class _FakeParticles:
    def __init__(self, position, tag, typeid):
        self.position = position
        self.tag = tag
        self.typeid = typeid


class _FakeSnap:
    def __init__(self, particles):
        self.particles = particles


class _FakeState:
    def __init__(self, particle_types, position, tag, typeid):
        self.particle_types = particle_types
        self._snap = _FakeSnap(_FakeParticles(position, tag, typeid))

    @property
    @contextmanager
    def cpu_local_snapshot(self):
        yield self._snap


class _FakeSim:
    def __init__(self, state):
        self.state = state


class TestSimReader:
    def test_reads_tag_ordered_ligands(self, manifold):
        # Two ligand particles + one non-ligand, in a SCRAMBLED row order with
        # tags that re-index them; the reader must recover the ligands by type
        # in tag order.
        sid = int(np.argmin(manifold.tri_centroids[:, 2]))
        n = manifold.tri_normals[sid]
        c = manifold.tri_centroids[sid]
        lig0 = c + 10.0e-9 * n
        lig1 = c + 12.0e-9 * n
        other = np.array([0.0, 0.0, 0.0])
        # Row order: [other, lig1, lig0]; tags map row→tag so tag order is
        # [lig0(tag0), lig1(tag1), other(tag2)].
        position = np.stack([other, lig1, lig0])
        tag = np.array([2, 1, 0])
        typeid = np.array([0, 1, 1])  # 0=cortex(other), 1=ligand
        state = _FakeState(["cortex", TYPE_LIGAND], position, tag, typeid)
        sim = _FakeSim(state)
        out = measure_ecm_contact_manifold(sim, manifold, contact_gap=50.0e-9)
        assert out["n_ligands"] == 2
        assert out["contact_mask"][sid]
        # Both ligands homed on the south patch, within contact.
        assert out["n_ligands_in_contact"] == 2

    def test_missing_ligand_type_is_empty(self, manifold):
        state = _FakeState(
            ["cortex"], np.zeros((1, 3)), np.array([0]), np.array([0])
        )
        out = measure_ecm_contact_manifold(_FakeSim(state), manifold, contact_gap=1e-7)
        assert out["n_ligands"] == 0
        assert out["n_contact_patches"] == 0


# ---------------------------------------------------------------------------
# Cross-consistency gate: FA-engaged patches ⊆ contact patches
# ---------------------------------------------------------------------------
class TestTractionWithinContact:
    def _contact_with(self, manifold, contact_ids):
        mask = np.zeros(manifold.n_tri, dtype=bool)
        mask[contact_ids] = True
        return {"contact_mask": mask}

    def test_subset_passes(self, manifold):
        contact = self._contact_with(manifold, [3, 7, 11, 42])
        fa = np.zeros(manifold.n_tri, dtype=bool)
        fa[[7, 42]] = True  # FA ⊂ contact
        gate = assert_traction_within_contact(contact, {"basal_patch_mask": fa})
        assert gate["passed"]
        assert gate["n_fa_patches"] == 2
        assert gate["n_fa_outside_contact"] == 0

    def test_fa_outside_contact_fails(self, manifold):
        contact = self._contact_with(manifold, [3, 7])
        fa = np.zeros(manifold.n_tri, dtype=bool)
        fa[[7, 99]] = True  # 99 ∉ contact
        gate = assert_traction_within_contact(contact, {"basal_patch_mask": fa})
        assert not gate["passed"]
        assert gate["n_fa_outside_contact"] == 1
        assert 99 in gate["fa_outside_ids"].tolist()

    def test_mismatched_manifolds_raise(self, manifold):
        contact = {"contact_mask": np.zeros(manifold.n_tri, dtype=bool)}
        fa = np.zeros(manifold.n_tri + 1, dtype=bool)
        with pytest.raises(ValueError):
            assert_traction_within_contact(contact, {"basal_patch_mask": fa})


# ---------------------------------------------------------------------------
# Type literal agrees with the FA single-source-of-truth
# ---------------------------------------------------------------------------
def test_ligand_type_matches_fa():
    from ffn_sim.archive.hoomd_legacy.bridge.fa import TYPE_LIGAND as FA_TYPE_LIGAND

    assert TYPE_LIGAND == FA_TYPE_LIGAND
