"""LAYOUT-level sanity-gate tests for the H.7 basal-ring lamellipodium geometry.

Covers ``ffn_sim/cell/lamellipodium_basal_ring.py`` Sanity Gate §1-6 (the
fast, static, no-full-HOOMD-run subset; the dynamic footprint-advance §6 check
belongs to the integrated membrane-ON run downstream):

- §1 Dimensional analysis / derivation (TestRingDerivation)
- §2 Boundary cases                     (TestBoundary)
- §3 Geometry invariants                (TestRingPlacement)
- §5 Sign / sense (OUTWARD-radial)      (TestRadialSignSense)

These are pure-NumPy geometry checks (no ``sim.run``), so the suite is fast.

Design: ``ffn_sim/docs/briefs/H7_LAMELLIPODIUM_SPHERICAL_INTEGRATION.md`` §3.
"""

from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.archive.hoomd_legacy.cell.lamellipodium import (
    LamellipodiumLayout,
    resolve_h5_lamellipodium,
)
from ffn_sim.archive.hoomd_legacy.cell.lamellipodium_basal_ring import (
    BasalRingGeometry,
    BasalRingLayout,
    generate_basal_ring_lamellipodium_layout,
    resolve_basal_ring_geometry,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h5.yaml"
)

# Physiological MCF7 single cell (CLAUDE.md platform conventions): R_cell =
# 7.5 µm (Wagner 2011), box_factor = 3 → L_box = 22.5 µm. We resolve the H.5
# rates against this box so the basal ring lives in the same box the full-cell
# cortex uses.
R_CELL = 7.5e-6
BOX_FACTOR = 3.0
L_BOX = BOX_FACTOR * R_CELL
DT = 13.0e-9


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _resolved(n_WAVE: int = 60, *, L_box: float = L_BOX) -> object:
    """Resolve H.5 params for the MCF7 single-cell box (R_cell = 7.5 µm).

    The flat-model ``Y_max`` is irrelevant to the basal-ring path, but the
    shared :func:`resolve_h5_lamellipodium` still validates it against the box
    half-edge. ``phase1_h5.yaml`` pins ``Y_max = 13.5 µm`` for its 30 µm box; we
    rescale it to ``0.45·L_box`` for the smaller MCF7 box so the resolve passes
    (this value is never read by the basal-ring layout).
    """
    cfg = deepcopy(_load_cfg())
    cfg["lamellipodium"]["n_WAVE"] = n_WAVE
    cfg["lamellipodium"]["Y_max"] = 0.45 * L_box
    return resolve_h5_lamellipodium(cfg, L_box=L_box, dt=DT)


@pytest.fixture(scope="module")
def p_lamel():
    return _resolved(n_WAVE=60)


# ---------------------------------------------------------------------------
# §1 — derivation of ring radius + basal z (no magic numbers)
# ---------------------------------------------------------------------------
class TestRingDerivation:
    def test_spherical_cap_base_radius_formula(self, p_lamel):
        """Ring radius = sqrt(cap_depth·(2 R_cell − cap_depth)) (exact cap base)."""
        cap_depth = 0.8e-6
        geom = resolve_basal_ring_geometry(
            p_lamel, R_cell=R_CELL, cap_depth=cap_depth,
        )
        expected = math.sqrt(cap_depth * (2.0 * R_CELL - cap_depth))
        assert geom.ring_radius == pytest.approx(expected, rel=1e-12)
        # Units: sqrt(m·m) = m, and the value is sub-R_cell for a shallow cap.
        assert 0.0 < geom.ring_radius < R_CELL

    def test_basal_plane_is_south_pole_plus_rest_length(self, p_lamel):
        """z_basal = −R_cell + ℓ₀ (south pole + one backbone rest-length)."""
        geom = resolve_basal_ring_geometry(p_lamel, R_cell=R_CELL)
        assert geom.z_south == pytest.approx(-R_CELL, rel=1e-12)
        assert geom.basal_z_offset == pytest.approx(p_lamel.rest_length, rel=1e-12)
        assert geom.z_basal == pytest.approx(
            -R_CELL + p_lamel.rest_length, rel=1e-12
        )
        # Lamella sits just ABOVE the south pole, below the equator.
        assert -R_CELL < geom.z_basal < 0.0

    def test_default_cap_depth_is_rest_length(self, p_lamel):
        """Default cap_depth = ℓ₀ → ring radius = sqrt(ℓ₀·(2R−ℓ₀))."""
        geom = resolve_basal_ring_geometry(p_lamel, R_cell=R_CELL)
        h = p_lamel.rest_length
        assert geom.cap_depth == pytest.approx(h, rel=1e-12)
        assert geom.ring_radius == pytest.approx(
            math.sqrt(h * (2.0 * R_CELL - h)), rel=1e-12
        )

    def test_contact_radius_frac_override_bypasses_cap(self, p_lamel):
        """Explicit fraction → ring_radius = frac·R_cell (cap formula bypassed)."""
        geom = resolve_basal_ring_geometry(
            p_lamel, R_cell=R_CELL, contact_radius_frac=0.6,
        )
        assert geom.ring_radius == pytest.approx(0.6 * R_CELL, rel=1e-12)
        assert geom.contact_radius_frac == pytest.approx(0.6, rel=1e-12)

    def test_all_geometry_fields_finite(self, p_lamel):
        geom = resolve_basal_ring_geometry(p_lamel, R_cell=R_CELL)
        for name in ("R_cell", "z_south", "cap_depth", "ring_radius",
                     "z_basal", "basal_z_offset"):
            val = getattr(geom, name)
            assert math.isfinite(val), f"{name} not finite: {val!r}"


# ---------------------------------------------------------------------------
# §2 — boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_zero_wave_empty_layout(self, p_lamel):
        """n_WAVE == 0 → empty (no ring), still a well-formed LamellipodiumLayout."""
        p0 = _resolved(n_WAVE=0)
        out = generate_basal_ring_lamellipodium_layout(
            p0, wave_tag_start=1000, R_cell=R_CELL,
        )
        assert isinstance(out, BasalRingLayout)
        assert isinstance(out.layout, LamellipodiumLayout)
        assert out.layout.wave_positions.shape == (0, 3)
        assert out.layout.mother_seed_positions.shape == (0, 3)
        assert out.layout.wave_to_mother_bond_pairs.shape == (0, 2)
        assert out.mother_tangents.shape == (0, 3)
        assert out.azimuths.shape == (0,)

    @pytest.mark.parametrize("n_WAVE", [1, 3])
    def test_small_n_wave_well_formed(self, n_WAVE):
        """Degenerate small ring is still on-ring + outward-radial."""
        p = _resolved(n_WAVE=n_WAVE)
        out = generate_basal_ring_lamellipodium_layout(
            p, wave_tag_start=0, R_cell=R_CELL,
        )
        assert out.layout.wave_positions.shape == (n_WAVE, 3)
        assert out.mother_tangents.shape == (n_WAVE, 3)
        # Each tangent is a unit vector in the basal plane.
        norms = np.linalg.norm(out.mother_tangents, axis=1)
        assert np.allclose(norms, 1.0, atol=1e-12)
        assert np.allclose(out.mother_tangents[:, 2], 0.0, atol=1e-12)

    def test_bad_R_cell_raises(self, p_lamel):
        for bad in (0.0, -1.0e-6, float("nan"), float("inf")):
            with pytest.raises(ValueError):
                resolve_basal_ring_geometry(p_lamel, R_cell=bad)

    def test_bad_contact_radius_frac_raises(self, p_lamel):
        for bad in (0.0, 1.0, 1.5, -0.2):
            with pytest.raises(ValueError):
                resolve_basal_ring_geometry(
                    p_lamel, R_cell=R_CELL, contact_radius_frac=bad,
                )

    def test_cap_depth_clamped_to_sphere(self, p_lamel):
        """cap_depth > 2·R_cell clamps to the whole sphere (base radius → 0)."""
        geom = resolve_basal_ring_geometry(
            p_lamel, R_cell=R_CELL, cap_depth=10.0 * R_CELL,
        )
        # h clamped to 2R → sqrt(2R·(2R − 2R)) = 0; real, not NaN.
        assert math.isfinite(geom.ring_radius)
        assert geom.ring_radius == pytest.approx(0.0, abs=1e-15)

    def test_ring_too_big_for_box_raises(self):
        """A ring that exceeds the box half-edge surfaces (no silent clip)."""
        # L_box = R_cell (box_factor 1) → half-edge = R_cell/2 < the |z_basal| ≈
        # R_cell extent of the south-pole lamella, so the ring cannot fit.
        p_tiny = _resolved(n_WAVE=10, L_box=R_CELL)
        with pytest.raises(ValueError):
            resolve_basal_ring_geometry(p_tiny, R_cell=R_CELL)


# ---------------------------------------------------------------------------
# §3 — geometry invariants: WAVE on ring, mother one ℓ₀ inward, inside box
# ---------------------------------------------------------------------------
class TestRingPlacement:
    def test_wave_beads_on_basal_ring(self, p_lamel):
        """All WAVE beads at radius r_contact in (x,y), z = z_basal."""
        out = generate_basal_ring_lamellipodium_layout(
            p_lamel, wave_tag_start=0, R_cell=R_CELL,
        )
        geom = out.geometry
        wp = out.layout.wave_positions
        r_xy = np.hypot(wp[:, 0], wp[:, 1])
        assert np.allclose(r_xy, geom.ring_radius, rtol=1e-10, atol=1e-12)
        assert np.allclose(wp[:, 2], geom.z_basal, rtol=0.0, atol=1e-12)

    def test_mothers_one_rest_length_inward(self, p_lamel):
        """Each mother seed sits one ℓ₀ inward of its WAVE along −t̂ (radial)."""
        out = generate_basal_ring_lamellipodium_layout(
            p_lamel, wave_tag_start=0, R_cell=R_CELL,
        )
        wp = out.layout.wave_positions
        mp = out.layout.mother_seed_positions
        # r_mother = r_wave − ℓ₀·t̂; with t̂ purely radial-in-plane, the in-plane
        # radius drops by exactly ℓ₀ and z is unchanged.
        r_wave_xy = np.hypot(wp[:, 0], wp[:, 1])
        r_mother_xy = np.hypot(mp[:, 0], mp[:, 1])
        assert np.allclose(
            r_wave_xy - r_mother_xy, p_lamel.rest_length, rtol=1e-9, atol=1e-12
        )
        assert np.allclose(mp[:, 2], wp[:, 2], atol=1e-12)
        # Displacement WAVE→mother has magnitude exactly ℓ₀.
        disp = wp - mp
        assert np.allclose(
            np.linalg.norm(disp, axis=1), p_lamel.rest_length, rtol=1e-10
        )

    def test_n_wave_count_consistent(self, p_lamel):
        out = generate_basal_ring_lamellipodium_layout(
            p_lamel, wave_tag_start=0, R_cell=R_CELL,
        )
        n = p_lamel.n_WAVE
        assert out.layout.wave_positions.shape == (n, 3)
        assert out.layout.mother_seed_positions.shape == (n, 3)
        assert out.layout.wave_to_mother_bond_pairs.shape == (n, 2)
        assert out.mother_tangents.shape == (n, 3)
        assert out.azimuths.shape == (n,)

    def test_tag_bookkeeping_matches_flat_schema(self, p_lamel):
        """Anchor pairs = (wave_tag, mother_tag) with the flat tag layout."""
        start = 4242
        out = generate_basal_ring_lamellipodium_layout(
            p_lamel, wave_tag_start=start, R_cell=R_CELL,
        )
        n = p_lamel.n_WAVE
        lay = out.layout
        assert lay.wave_tag_start == start
        assert lay.mother_tag_start == start + n
        pairs = lay.wave_to_mother_bond_pairs
        assert np.array_equal(pairs[:, 0], np.arange(n) + start)
        assert np.array_equal(pairs[:, 1], np.arange(n) + start + n)

    def test_all_positions_inside_box_half_edge(self, p_lamel):
        out = generate_basal_ring_lamellipodium_layout(
            p_lamel, wave_tag_start=0, R_cell=R_CELL,
        )
        half = 0.5 * p_lamel.L_box
        for arr in (out.layout.wave_positions, out.layout.mother_seed_positions):
            assert np.all(np.abs(arr) < half), "position outside box half-edge"


# ---------------------------------------------------------------------------
# §5 — sign / sense: per-WAVE tangent is OUTWARD-radial in the basal plane
# ---------------------------------------------------------------------------
class TestRadialSignSense:
    def test_tangents_are_outward_radial(self, p_lamel):
        """t̂ · r̂_out > 0 (here ≈ 1) and t_z ≈ 0 for every WAVE."""
        out = generate_basal_ring_lamellipodium_layout(
            p_lamel, wave_tag_start=0, R_cell=R_CELL,
        )
        wp = out.layout.wave_positions
        tang = out.mother_tangents
        # Outward radial unit at each WAVE (in the basal plane).
        r_xy = wp[:, :2]
        r_out = r_xy / np.linalg.norm(r_xy, axis=1, keepdims=True)
        dot = tang[:, 0] * r_out[:, 0] + tang[:, 1] * r_out[:, 1]
        # SIGN/SENSE GATE: strictly outward (not inward, not tangential).
        assert np.all(dot > 0.0), "a mother tangent points inward/tangential"
        # By construction t̂ IS the outward radial in-plane unit → dot ≈ 1.
        assert np.allclose(dot, 1.0, atol=1e-9)
        # Growth stays in the basal plane.
        assert np.allclose(tang[:, 2], 0.0, atol=1e-12)

    def test_tangents_are_unit_vectors(self, p_lamel):
        out = generate_basal_ring_lamellipodium_layout(
            p_lamel, wave_tag_start=0, R_cell=R_CELL,
        )
        norms = np.linalg.norm(out.mother_tangents, axis=1)
        assert np.allclose(norms, 1.0, atol=1e-12)

    def test_tangents_azimuthally_distributed(self):
        """Many WAVE → tangents spread around the ring (mean ≈ 0, isotropic)."""
        p = _resolved(n_WAVE=400)
        out = generate_basal_ring_lamellipodium_layout(
            p, wave_tag_start=0, R_cell=R_CELL,
        )
        tang = out.mother_tangents
        mean_dir = tang.mean(axis=0)
        # A full uniform ring has zero net tangent (no preferred direction).
        # 400 random azimuths → |mean| is O(1/sqrt(N)) ≈ 0.05; bound loosely.
        assert np.linalg.norm(mean_dir[:2]) < 0.2, (
            f"tangents not azimuthally balanced: mean={mean_dir!r}"
        )
        # Azimuths must actually span the circle (not clustered).
        phi = out.azimuths
        assert phi.min() < 0.5 and phi.max() > 2.0 * math.pi - 0.5

    def test_mother_to_wave_points_outward(self, p_lamel):
        """The mother→WAVE displacement is the outward growth direction."""
        out = generate_basal_ring_lamellipodium_layout(
            p_lamel, wave_tag_start=0, R_cell=R_CELL,
        )
        wp = out.layout.wave_positions
        mp = out.layout.mother_seed_positions
        grow_dir = wp - mp  # the barbed end advances from mother toward WAVE
        grow_unit = grow_dir / np.linalg.norm(grow_dir, axis=1, keepdims=True)
        # Equals the stored outward-radial tangent.
        assert np.allclose(grow_unit, out.mother_tangents, atol=1e-9)
        # And points outward (mother is closer to centre than WAVE).
        r_wave = np.linalg.norm(wp[:, :2], axis=1)
        r_mother = np.linalg.norm(mp[:, :2], axis=1)
        assert np.all(r_wave > r_mother)


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
def test_layout_is_seed_reproducible(p_lamel):
    rng1 = np.random.default_rng(123)
    rng2 = np.random.default_rng(123)
    a = generate_basal_ring_lamellipodium_layout(
        p_lamel, wave_tag_start=0, R_cell=R_CELL, rng=rng1,
    )
    b = generate_basal_ring_lamellipodium_layout(
        p_lamel, wave_tag_start=0, R_cell=R_CELL, rng=rng2,
    )
    assert np.array_equal(a.layout.wave_positions, b.layout.wave_positions)
    assert np.array_equal(a.mother_tangents, b.mother_tangents)
    assert np.array_equal(a.azimuths, b.azimuths)
