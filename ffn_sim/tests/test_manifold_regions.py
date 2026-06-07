"""Sanity-gate tests for ``cortex/manifold_regions.py`` (H.7 spatial substrate).

Covers the manifold-region masks' Sanity Gate (CLAUDE.md): dimensional / unit
checks, boundary (empty region), resolution invariance of the region AREA
FRACTION (the manifold must not let mesh resolution control the physical region
size), sign/sense of the per-patch in-plane growth direction, and — the wiring
proof — CROSS-CONSISTENCY with the lamellipodium WAVE placement (every WAVE
bead's home patch lies inside the corresponding mask).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.cortex.manifold_regions import (
    ManifoldRegion,
    basal_ring_region,
    patch_polar_azimuth,
    polarized_patch_region,
)
from ffn_sim.cortex.surface_manifold import SurfaceManifold

# Physiological MCF7 single cell (CLAUDE.md platform conventions).
R_CELL = 7.5e-6
REST_LENGTH = 0.5e-6  # ℓ₀ backbone bond (H.5 default)


@pytest.fixture(scope="module")
def manifold() -> SurfaceManifold:
    """Level-3 icosphere (1280 faces) at the MCF7 cortex radius."""
    return SurfaceManifold.icosphere(subdivisions=3, radius=R_CELL)


# ---------------------------------------------------------------------------
# Per-patch angles
# ---------------------------------------------------------------------------
class TestPolarAzimuth:
    def test_south_pole_is_theta_zero(self, manifold):
        theta, _ = patch_polar_azimuth(manifold.tri_centroids)
        # The single most-southern patch (max −z) is near θ = 0.
        south = int(np.argmin(manifold.tri_centroids[:, 2]))
        assert theta[south] < 0.3  # within a patch of the pole
        # θ is bounded to [0, π].
        assert theta.min() >= 0.0 and theta.max() <= math.pi + 1e-9

    def test_azimuth_range(self, manifold):
        _, phi = patch_polar_azimuth(manifold.tri_centroids)
        assert phi.min() >= -math.pi - 1e-9 and phi.max() <= math.pi + 1e-9


# ---------------------------------------------------------------------------
# basal_ring collar
# ---------------------------------------------------------------------------
class TestBasalRing:
    def test_collar_nonempty_and_straddles_contact_latitude(self, manifold):
        reg = basal_ring_region(
            manifold, R_cell=R_CELL, rest_length=REST_LENGTH,
            collar_half_angle=0.2,
        )
        assert isinstance(reg, ManifoldRegion)
        assert reg.name == "basal_ring"
        assert reg.patch_ids.size > 0
        theta_c = reg.meta["theta_c_rad"]
        # Every selected patch sits within the collar of the contact latitude.
        assert np.all(np.abs(reg.theta - theta_c) <= 0.2 + 1e-9)
        # The collar spans the FULL azimuth (isotropic ring): φ covers > 5 rad.
        assert reg.phi.max() - reg.phi.min() > 5.0

    def test_area_fraction_in_unit_interval(self, manifold):
        reg = basal_ring_region(
            manifold, R_cell=R_CELL, rest_length=REST_LENGTH, collar_half_angle=0.2,
        )
        assert 0.0 < reg.area_fraction < 1.0
        # area == Σ selected patch areas.
        assert reg.area == pytest.approx(
            float(np.sum(manifold.tri_areas[reg.mask])), rel=1e-12
        )

    def test_in_plane_dir_is_tangent_and_outward(self, manifold):
        reg = basal_ring_region(
            manifold, R_cell=R_CELL, rest_length=REST_LENGTH, collar_half_angle=0.2,
        )
        # Unit vectors.
        assert np.allclose(np.linalg.norm(reg.in_plane_dir, axis=1), 1.0, atol=1e-9)
        # Tangent to the patch: n̂ · dir ≈ 0.
        dots = np.einsum("ij,ij->i", reg.in_plane_dir, reg.normals)
        assert np.all(np.abs(dots) < 1e-9)
        # Outward-radial: positive projection on the xy-radial direction.
        r_xy = reg.centroids.copy()
        r_xy[:, 2] = 0.0
        r_xy /= np.linalg.norm(r_xy, axis=1, keepdims=True)
        radial_dot = np.einsum("ij,ij->i", reg.in_plane_dir, r_xy)
        assert np.all(radial_dot > 0.0)

    def test_contact_radius_frac_path(self, manifold):
        frac = 0.5
        reg = basal_ring_region(
            manifold, R_cell=R_CELL, rest_length=REST_LENGTH,
            contact_radius_frac=frac, collar_half_angle=0.2,
        )
        # θ_c = asin(frac); sin θ_c == frac.
        assert math.sin(reg.meta["theta_c_rad"]) == pytest.approx(frac, rel=1e-12)

    def test_bad_inputs_raise(self, manifold):
        with pytest.raises(ValueError):
            basal_ring_region(manifold, R_cell=-1.0, rest_length=REST_LENGTH)
        with pytest.raises(ValueError):
            basal_ring_region(manifold, R_cell=R_CELL, rest_length=0.0)
        with pytest.raises(ValueError):
            basal_ring_region(
                manifold, R_cell=R_CELL, rest_length=REST_LENGTH,
                contact_radius_frac=1.5,
            )


# ---------------------------------------------------------------------------
# polarized_patch
# ---------------------------------------------------------------------------
class TestPolarizedPatch:
    def test_patch_nonempty_and_localized(self, manifold):
        reg = polarized_patch_region(
            manifold, R_cell=R_CELL, half_angle_azimuth=math.pi / 6.0,
        )
        assert reg.name == "polarized_patch"
        assert reg.patch_ids.size > 0
        # A localized cap (not the whole shell): a small minority of patches.
        assert reg.patch_ids.size < 0.25 * manifold.n_tri
        assert reg.area_fraction < 0.2
        # Every selected patch lies within the derived geodesic cap radius of the
        # leading-edge centre direction (definitional, pole-safe).
        center = np.asarray(reg.meta["cap_center_dir"])
        cap_r = reg.meta["cap_radius_rad"]
        dirs = reg.centroids / np.linalg.norm(reg.centroids, axis=1, keepdims=True)
        geo = np.arccos(np.clip(dirs @ center, -1.0, 1.0))
        assert np.all(geo <= cap_r + 1e-9)

    def test_in_plane_dir_is_tangent_and_forward(self, manifold):
        p_hat = np.array([1.0, 0.0, 0.0])
        reg = polarized_patch_region(
            manifold, R_cell=R_CELL, polarization=p_hat,
            half_angle_azimuth=math.pi / 6.0,
        )
        assert np.allclose(np.linalg.norm(reg.in_plane_dir, axis=1), 1.0, atol=1e-9)
        dots = np.einsum("ij,ij->i", reg.in_plane_dir, reg.normals)
        assert np.all(np.abs(dots) < 1e-9)  # tangent
        # Forward: positive projection on p̂.
        fwd = reg.in_plane_dir @ p_hat
        assert np.all(fwd > 0.0)

    def test_patch_rotates_with_polarization(self, manifold):
        reg_x = polarized_patch_region(
            manifold, R_cell=R_CELL, polarization=(1.0, 0.0, 0.0),
        )
        reg_y = polarized_patch_region(
            manifold, R_cell=R_CELL, polarization=(0.0, 1.0, 0.0),
        )
        # Rigid 90° rotation about the box axis → φ₀ shifts by π/2 and the
        # patch's mean (area-weighted) azimuth follows. (The two caps still share
        # the south-pole vicinity, so they are NOT disjoint — but the bulk moves.)
        assert reg_x.meta["phi0_rad"] == pytest.approx(0.0, abs=1e-9)
        assert reg_y.meta["phi0_rad"] == pytest.approx(math.pi / 2.0, abs=1e-9)

        def _mean_azimuth(reg):
            w = manifold.tri_areas[reg.mask]
            cx = np.sum(w * np.cos(reg.phi)) / np.sum(w)
            cy = np.sum(w * np.sin(reg.phi)) / np.sum(w)
            return math.atan2(cy, cx)

        d = (_mean_azimuth(reg_y) - _mean_azimuth(reg_x)) % (2.0 * math.pi)
        assert d == pytest.approx(math.pi / 2.0, abs=0.3)
        # The masks are clearly different (the rotation moved most of the patch).
        assert np.count_nonzero(reg_x.mask ^ reg_y.mask) > np.count_nonzero(
            reg_x.mask & reg_y.mask
        )

    def test_vertical_polarization_raises(self, manifold):
        with pytest.raises(ValueError):
            polarized_patch_region(
                manifold, R_cell=R_CELL, polarization=(0.0, 0.0, 1.0),
            )

    def test_ring_azimuth_raises(self, manifold):
        with pytest.raises(ValueError):
            polarized_patch_region(
                manifold, R_cell=R_CELL, half_angle_azimuth=math.pi + 0.1,
            )


# ---------------------------------------------------------------------------
# Boundary: empty region
# ---------------------------------------------------------------------------
class TestEmptyRegion:
    def test_tiny_collar_far_from_any_patch_is_empty(self, manifold):
        # A vanishingly thin collar at a latitude between patch centroids can
        # select nothing — the region must be well-formed and empty, not raise.
        reg = basal_ring_region(
            manifold, R_cell=R_CELL, rest_length=REST_LENGTH,
            contact_radius_frac=0.999_999, collar_half_angle=1e-6,
        )
        assert reg.patch_ids.size == 0
        assert reg.area == 0.0
        assert reg.area_fraction == 0.0
        assert reg.in_plane_dir.shape == (0, 3)


# ---------------------------------------------------------------------------
# Resolution invariance — mesh resolution must NOT control the region size
# ---------------------------------------------------------------------------
class TestResolutionInvariance:
    @pytest.mark.parametrize(
        "region_fn",
        [
            lambda m: basal_ring_region(
                m, R_cell=R_CELL, rest_length=REST_LENGTH, collar_half_angle=0.25,
            ),
            lambda m: polarized_patch_region(
                m, R_cell=R_CELL, half_angle_azimuth=math.pi / 6.0,
            ),
        ],
        ids=["basal_ring", "polarized_patch"],
    )
    def test_area_fraction_converges_across_grid(self, region_fn):
        fracs = []
        for s in (2, 3, 4):  # 320 / 1280 / 5120 faces — canonical grid
            m = SurfaceManifold.icosphere(subdivisions=s, radius=R_CELL)
            fracs.append(region_fn(m).area_fraction)
        # The finite-width angular region's area fraction must converge: the
        # coarse→fine change must shrink (resolution refining, not drifting).
        d_coarse = abs(fracs[1] - fracs[0])
        d_fine = abs(fracs[2] - fracs[1])
        assert d_fine <= d_coarse + 1e-3
        # And the absolute fraction must be stable to within a few % across grid.
        assert max(fracs) - min(fracs) < 0.05


# ---------------------------------------------------------------------------
# Cross-consistency with the lamellipodium WAVE placement (the WIRING proof)
# ---------------------------------------------------------------------------
class TestLamellipodiumCoRegistration:
    """Every WAVE bead's home patch must lie inside the matching region mask.

    This proves the manifold mask and the lamellipodium particle placement are
    the SAME region in two representations — the point of wiring the regions onto
    the spatial substrate.
    """

    @staticmethod
    def _resolved_h5(n_WAVE: int):
        from ffn_sim.cell.lamellipodium import resolve_h5_lamellipodium

        cfg_path = (
            Path(__file__).resolve().parents[1] / "configs" / "phase1_h5.yaml"
        )
        with open(cfg_path) as f:
            cfg = yaml.safe_load(f)
        L_box = 3.0 * R_CELL
        cfg["lamellipodium"]["n_WAVE"] = n_WAVE
        cfg["lamellipodium"]["Y_max"] = 0.45 * L_box
        return resolve_h5_lamellipodium(cfg, L_box=L_box, dt=13.0e-9)

    def test_basal_ring_wave_home_patches_in_mask(self, manifold):
        from ffn_sim.cell.lamellipodium_basal_ring import (
            generate_basal_ring_lamellipodium_layout,
        )

        p = self._resolved_h5(n_WAVE=120)
        wrapped = generate_basal_ring_lamellipodium_layout(
            p, wave_tag_start=0, R_cell=R_CELL,
        )
        wave_pos = wrapped.layout.wave_positions
        home = manifold.nearest_patch(wave_pos)
        # Collar wide enough to contain the home patches: ≥ a few patch widths
        # (the ring is thin; one patch is ~mean_edge/R wide). This is the
        # co-registration claim — the mask CONTAINS the WAVE footprint.
        collar = 3.0 * manifold.mean_edge_length / R_CELL
        reg = basal_ring_region(
            manifold, R_cell=R_CELL, rest_length=p.rest_length,
            cap_depth=wrapped.geometry.cap_depth, collar_half_angle=collar,
        )
        assert np.all(reg.mask[home]), (
            "some basal-ring WAVE home patches fall outside the collar mask"
        )

    def test_polarized_patch_wave_home_patches_in_mask(self, manifold):
        from ffn_sim.cell.lamellipodium_polarized_patch import (
            generate_polarized_patch_layout,
        )

        p = self._resolved_h5(n_WAVE=120)
        # Generate WAVE with a SMALLER patch than the mask so boundary beads can't
        # map to a centroid just outside the (identical-angle) mask.
        wrapped = generate_polarized_patch_layout(
            p, wave_tag_start=0, R_cell=R_CELL,
            polarization=(1.0, 0.0, 0.0),
            half_angle_azimuth=math.pi / 9.0, half_angle_linear=math.pi / 9.0,
        )
        wave_pos = wrapped.layout.wave_positions
        home = manifold.nearest_patch(wave_pos)
        reg = polarized_patch_region(
            manifold, R_cell=R_CELL, polarization=(1.0, 0.0, 0.0),
            half_angle_azimuth=math.pi / 6.0, half_angle_linear=math.pi / 6.0,
        )
        assert np.all(reg.mask[home]), (
            "some polarized-patch WAVE home patches fall outside the patch mask"
        )
