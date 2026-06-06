"""Sanity-Gate tests for the polarized leading-edge patch lamellipodium.

H.7 approach (b) — single migrating cell front (``geometry="polarized_patch"``).
These are LAYOUT-level static checks (fast, no full HOOMD run): they assert the
WAVE/NPF reservoir is a LOCALIZED patch (not an azimuthal ring), that every
mother barbed-end tangent points FORWARD (aligned with the polarization ``p̂``),
that counts / containment / basal placement are correct, and that changing
``p̂`` rigidly rotates the patch.

Sanity-Gate axes covered (CLAUDE.md hard rule):
- §2 boundary: ``n_WAVE == 0`` empty layout; degenerate (vertical) ``p̂``
  raises; out-of-range half-angle raises.
- §3 invariants: ``n_WAVE`` WAVE beads, ``n_WAVE`` mothers, ``n_WAVE`` anchor
  bonds; all positions inside the box half-edge.
- §5 sign / sense: per-WAVE tangents have positive ``p̂``-projection (forward
  fan); mean tangent ≈ ``p̂``; WAVE azimuthal spread is BOUNDED (a patch, not a
  ring); patch sits at the basal contact (south cap, ``z`` near ``−R_cell``);
  changing ``p̂`` rotates the patch.
"""

from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.cell.lamellipodium import (
    LamellipodiumLayout,
    resolve_h5_lamellipodium,
)  # resolve_h5_lamellipodium builds a valid ResolvedH5 in the fixture below.
from ffn_sim.cell.lamellipodium_polarized_patch import (
    PolarizedPatchLayout,
    generate_polarized_patch_layout,
)


# ---------------------------------------------------------------------------
# Fixtures — reuse the KU-5.x H.5 config for a valid ResolvedH5, with the box
# made self-consistent with the cortex sphere (L_box = box_factor · R_cell).
# ---------------------------------------------------------------------------
CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h5.yaml"
)

R_CELL = 7.5e-6           # m   MCF7 cortex radius (cortex default)
BOX_FACTOR = 3.0          # cortex.py default box_factor
L_BOX = BOX_FACTOR * R_CELL
DT = 13.0e-9              # s   host integrator dt (matches test_lamellipodium)


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _resolved(n_WAVE: int = 60, *, L_box: float = L_BOX):
    cfg = deepcopy(_load_cfg())
    cfg["lamellipodium"]["n_WAVE"] = n_WAVE
    # Y_max / wave_area are FLAT-PLANE parameters the polarized patch does not
    # use; scale Y_max with the (cortex-consistent) box so the flat-plane
    # validator (Y_max < L_box/2) passes for this box size.
    cfg["lamellipodium"]["Y_max"] = 0.45 * L_box
    return resolve_h5_lamellipodium(cfg, L_box=L_box, dt=DT)


@pytest.fixture(scope="module")
def p_lamel():
    return _resolved(n_WAVE=60)


def _make(
    p,
    *,
    polarization=None,
    half_angle_azimuth: float = math.pi / 6.0,
    half_angle_linear: float | None = None,
    fan_spread: float | None = None,
    seed: int = 0,
) -> PolarizedPatchLayout:
    return generate_polarized_patch_layout(
        p,
        wave_tag_start=0,
        R_cell=R_CELL,
        polarization=polarization,
        half_angle_azimuth=half_angle_azimuth,
        half_angle_linear=half_angle_linear,
        fan_spread=fan_spread,
        rng=np.random.default_rng(seed),
    )


def _azimuths(positions: np.ndarray) -> np.ndarray:
    """Azimuth (atan2(y, x)) of each row about the box axis."""
    return np.arctan2(positions[:, 1], positions[:, 0])


def _angular_spread(azimuths: np.ndarray) -> float:
    """Peak-to-peak angular spread on the circle (handles wrap-around).

    Returns the smallest arc (in radians) that contains all the azimuths.
    """
    a = np.sort(np.mod(azimuths, 2.0 * math.pi))
    if a.size < 2:
        return 0.0
    gaps = np.diff(np.concatenate([a, a[:1] + 2.0 * math.pi]))
    largest_gap = float(gaps.max())
    return 2.0 * math.pi - largest_gap


# ---------------------------------------------------------------------------
# Return-schema compatibility with LamellipodiumLayout
# ---------------------------------------------------------------------------
class TestSchema:
    def test_returns_lamellipodium_layout(self, p_lamel):
        out = _make(p_lamel)
        assert isinstance(out, PolarizedPatchLayout)
        assert isinstance(out.layout, LamellipodiumLayout)
        lay = out.layout
        # Same fields / shapes the flat generator returns.
        assert lay.wave_positions.shape == (p_lamel.n_WAVE, 3)
        assert lay.mother_seed_positions.shape == (p_lamel.n_WAVE, 3)
        assert lay.wave_to_mother_bond_pairs.shape == (p_lamel.n_WAVE, 2)
        assert lay.wave_tag_start == 0
        assert lay.mother_tag_start == p_lamel.n_WAVE
        assert out.mother_tangents.shape == (p_lamel.n_WAVE, 3)

    def test_anchor_bond_pairs_are_wave_to_mother(self, p_lamel):
        lay = _make(p_lamel).layout
        n = p_lamel.n_WAVE
        # WAVE tags [0, n); mother tags [n, 2n); each anchor pairs them 1:1.
        np.testing.assert_array_equal(
            lay.wave_to_mother_bond_pairs[:, 0], np.arange(n)
        )
        np.testing.assert_array_equal(
            lay.wave_to_mother_bond_pairs[:, 1], np.arange(n) + n
        )

    def test_wave_tag_start_offset_respected(self, p_lamel):
        out = generate_polarized_patch_layout(
            p_lamel, wave_tag_start=137, R_cell=R_CELL,
            rng=np.random.default_rng(0),
        )
        lay = out.layout
        assert lay.wave_tag_start == 137
        assert lay.mother_tag_start == 137 + p_lamel.n_WAVE
        assert int(lay.wave_to_mother_bond_pairs[0, 0]) == 137


# ---------------------------------------------------------------------------
# §3 invariants — counts + containment
# ---------------------------------------------------------------------------
class TestCountsAndContainment:
    def test_counts_consistent(self, p_lamel):
        lay = _make(p_lamel).layout
        n = p_lamel.n_WAVE
        assert lay.wave_positions.shape[0] == n
        assert lay.mother_seed_positions.shape[0] == n
        assert lay.wave_to_mother_bond_pairs.shape[0] == n

    def test_all_positions_inside_box_half_edge(self, p_lamel):
        out = _make(p_lamel)
        half = 0.5 * p_lamel.L_box
        for arr in (out.layout.wave_positions, out.layout.mother_seed_positions):
            assert np.all(np.abs(arr) < half)

    def test_wave_beads_on_cortex_shell(self, p_lamel):
        # WAVE beads are scattered on the spherical shell |r| = R_cell (the
        # cortex surface), so the patch is co-located with the cortex south cap.
        out = _make(p_lamel)
        radii = np.linalg.norm(out.layout.wave_positions, axis=1)
        np.testing.assert_allclose(radii, R_CELL, rtol=1e-9)

    def test_mother_offset_one_rest_length_from_wave(self, p_lamel):
        # Mother seed is one rest-length from its WAVE (behind, along tangent).
        out = _make(p_lamel)
        d = np.linalg.norm(
            out.layout.wave_positions - out.layout.mother_seed_positions, axis=1
        )
        np.testing.assert_allclose(d, p_lamel.rest_length, rtol=1e-9)


# ---------------------------------------------------------------------------
# §5 sign / sense — LOCALIZED patch (not a ring)
# ---------------------------------------------------------------------------
class TestLocalizedPatch:
    def test_azimuthal_spread_is_bounded_not_2pi(self, p_lamel):
        # The defining property of a PATCH (vs a ring): WAVE azimuths cluster
        # within a bounded arc about p̂, NOT spread over the full 2π.
        half_az = math.pi / 6.0
        out = _make(p_lamel, half_angle_azimuth=half_az)
        az = _azimuths(out.layout.wave_positions)
        spread = _angular_spread(az)
        # Bounded by the patch's full azimuthal extent (2·half_angle), with a
        # small allowance for the random scatter not perfectly hitting the
        # edges; HARD requirement is simply spread < π (a patch, not a ring).
        assert spread <= 2.0 * half_az + 1.0e-9
        assert spread < math.pi

    def test_wider_half_angle_widens_patch(self, p_lamel):
        narrow = _make(p_lamel, half_angle_azimuth=math.pi / 9.0, seed=1)
        wide = _make(p_lamel, half_angle_azimuth=math.pi / 3.0, seed=1)
        s_narrow = _angular_spread(_azimuths(narrow.layout.wave_positions))
        s_wide = _angular_spread(_azimuths(wide.layout.wave_positions))
        assert s_wide > s_narrow

    def test_patch_centered_on_polarization_azimuth(self, p_lamel):
        # Default p̂ = +x̂ → patch centred on azimuth 0.
        out = _make(p_lamel)
        az = _azimuths(out.layout.wave_positions)
        # Circular mean azimuth ≈ 0 (within the patch half-extent).
        mean_az = math.atan2(float(np.mean(np.sin(az))), float(np.mean(np.cos(az))))
        assert abs(mean_az) < math.pi / 6.0


# ---------------------------------------------------------------------------
# §5 sign / sense — basal (south-cap) placement
# ---------------------------------------------------------------------------
class TestBasalPlacement:
    def test_patch_near_south_cap(self, p_lamel):
        # The lamella sits at the basal contact (south cap): all WAVE beads in
        # the southern hemisphere (z < 0), near the south pole z ≈ −R_cell.
        out = _make(p_lamel)
        z = out.layout.wave_positions[:, 2]
        assert np.all(z < 0.0)
        # The patch spans the polar band θ_c ± Δθ_lin about the contact circle,
        # so z lies between the lower (θ_c − Δθ_lin) and upper (θ_c + Δθ_lin)
        # cap edges — all still in the basal/south region (z < 0).
        theta_c = math.acos(-out.z_basal / R_CELL)
        z_lower = -R_CELL * math.cos(max(theta_c - out.half_angle_linear, 0.0))
        z_upper = -R_CELL * math.cos(
            min(theta_c + out.half_angle_linear, math.pi)
        )
        assert np.all(z >= z_lower - 1.0e-12)
        assert np.all(z <= z_upper + 1.0e-12)

    def test_basal_plane_offset_derived_from_geometry(self, p_lamel):
        out = _make(p_lamel)
        # z_basal = −R_cell + R_cell·(1 − cos Δθ_lin) = −R·cos(Δθ_lin).
        expected_z = -R_CELL * math.cos(out.half_angle_linear)
        assert math.isclose(out.z_basal, expected_z, rel_tol=1e-12)
        # Contact-circle radius matches the basal plane cutting the shell.
        expected_r = math.sqrt(R_CELL**2 - out.z_basal**2)
        assert math.isclose(out.r_contact, expected_r, rel_tol=1e-12)

    def test_patch_center_is_forwardmost_contact_point(self, p_lamel):
        out = _make(p_lamel)
        # Default p̂ = +x̂: centre = (r_contact, 0, z_basal).
        np.testing.assert_allclose(
            out.patch_center,
            np.array([out.r_contact, 0.0, out.z_basal]),
            atol=1e-15,
        )


# ---------------------------------------------------------------------------
# §5 sign / sense — FORWARD growth (barbed-end tangents aligned with p̂)
# ---------------------------------------------------------------------------
class TestForwardGrowth:
    def test_all_tangents_point_forward(self, p_lamel):
        out = _make(p_lamel)
        p_hat = out.polarization
        dots = out.mother_tangents @ p_hat
        # Every mother barbed end grows FORWARD: positive p̂-projection.
        assert np.all(dots > 0.0), (
            f"min tangent·p̂ = {float(dots.min()):.4f} (must be > 0)"
        )

    def test_mean_tangent_aligns_with_polarization(self, p_lamel):
        # On the curved south-cap shell the local surface normal points mostly
        # DOWN (−z), so a forward (p̂) growth direction kept tangent to the
        # surface necessarily acquires a +z component (the protrusion advances
        # forward AND rides up the shell from the contact toward the equator —
        # the physical advance of the contact line). The meaningful alignment
        # check is therefore the IN-PLANE (basal x,y) component, plus the
        # absence of any lateral (cross-p̂) bias.
        out = _make(p_lamel)
        p_hat = out.polarization
        t_in = out.mother_tangents.copy()
        t_in[:, 2] = 0.0
        t_in /= np.linalg.norm(t_in, axis=1, keepdims=True)
        # Every basal-projected tangent points forward; the mean is ≈ p̂.
        assert np.all(t_in @ p_hat > 0.0)
        mean_in = t_in.mean(axis=0)
        mean_in /= np.linalg.norm(mean_in)
        assert float(mean_in @ p_hat) > math.cos(math.pi / 12.0)
        # No spurious lateral drift: mean basal lateral component ≈ 0.
        lateral = np.array([-p_hat[1], p_hat[0], 0.0])
        assert abs(float(out.mother_tangents.mean(axis=0) @ lateral)) < 0.1

    def test_tangents_are_unit_vectors(self, p_lamel):
        out = _make(p_lamel)
        norms = np.linalg.norm(out.mother_tangents, axis=1)
        np.testing.assert_allclose(norms, 1.0, rtol=1e-9)

    def test_tangents_tangent_to_surface(self, p_lamel):
        # Forward growth is tangent to the cell surface (protrusion advances the
        # front, it does not dive into / lift off the substrate): t̂ · n̂ ≈ 0.
        out = _make(p_lamel)
        normals = out.layout.wave_positions / np.linalg.norm(
            out.layout.wave_positions, axis=1, keepdims=True
        )
        radial = np.sum(out.mother_tangents * normals, axis=1)
        np.testing.assert_allclose(radial, 0.0, atol=1e-9)

    def test_zero_fan_is_deterministic_and_forward(self, p_lamel):
        # With fan_spread = 0 the tangent is p̂ projected onto each bead's local
        # tangent plane (no random tilt): forward, reproducible (no RNG), and
        # the patch-averaged lateral bias vanishes (the patch is symmetric
        # about p̂). NOTE: per-bead the projected tangent is NOT laterally pure
        # — a bead off the p̂ meridian projects p̂ into a tilted tangent plane —
        # which is why the lateral check is on the MEAN, not per-bead.
        out_a = _make(p_lamel, fan_spread=0.0, seed=7)
        out_b = _make(p_lamel, fan_spread=0.0, seed=7)
        np.testing.assert_array_equal(out_a.mother_tangents, out_b.mother_tangents)
        p_hat = out_a.polarization
        assert np.all(out_a.mother_tangents @ p_hat > 0.0)
        lateral = np.array([-p_hat[1], p_hat[0], 0.0])
        assert abs(float(out_a.mother_tangents.mean(axis=0) @ lateral)) < 0.05


# ---------------------------------------------------------------------------
# §5 sign / sense — changing p̂ rotates the patch
# ---------------------------------------------------------------------------
class TestRotation:
    def test_changing_polarization_rotates_patch(self, p_lamel):
        out_x = _make(p_lamel, polarization=(1.0, 0.0, 0.0), seed=3)
        out_y = _make(p_lamel, polarization=(0.0, 1.0, 0.0), seed=3)
        # Patch centre azimuth rotates from 0 (+x) to π/2 (+y).
        az_x = math.atan2(out_x.patch_center[1], out_x.patch_center[0])
        az_y = math.atan2(out_y.patch_center[1], out_y.patch_center[0])
        assert math.isclose(az_x, 0.0, abs_tol=1e-9)
        assert math.isclose(az_y, 0.5 * math.pi, abs_tol=1e-9)
        # Mean WAVE azimuth also rotates by ≈ π/2.
        maz_x = math.atan2(
            float(np.mean(np.sin(_azimuths(out_x.layout.wave_positions)))),
            float(np.mean(np.cos(_azimuths(out_x.layout.wave_positions)))),
        )
        maz_y = math.atan2(
            float(np.mean(np.sin(_azimuths(out_y.layout.wave_positions)))),
            float(np.mean(np.cos(_azimuths(out_y.layout.wave_positions)))),
        )
        delta = (maz_y - maz_x) % (2.0 * math.pi)
        assert abs(delta - 0.5 * math.pi) < math.pi / 9.0

    def test_diagonal_polarization_dropped_to_inplane(self, p_lamel):
        # A p̂ with a z-component is projected to the basal (x, y) plane.
        out = _make(p_lamel, polarization=(1.0, 0.0, 5.0))
        np.testing.assert_allclose(
            out.polarization, np.array([1.0, 0.0, 0.0]), atol=1e-15
        )


# ---------------------------------------------------------------------------
# §2 boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_zero_wave_empty_layout(self):
        p = _resolved(n_WAVE=0)
        out = _make(p)
        assert out.layout.wave_positions.shape == (0, 3)
        assert out.layout.mother_seed_positions.shape == (0, 3)
        assert out.layout.wave_to_mother_bond_pairs.shape == (0, 2)
        assert out.mother_tangents.shape == (0, 3)
        # Geometry fields still resolved (so callers can introspect).
        assert out.z_basal < 0.0
        assert out.r_contact > 0.0

    def test_vertical_polarization_raises(self, p_lamel):
        with pytest.raises(ValueError, match="in-plane"):
            _make(p_lamel, polarization=(0.0, 0.0, 1.0))

    def test_out_of_range_azimuth_half_angle_raises(self, p_lamel):
        with pytest.raises(ValueError, match="half_angle_azimuth"):
            _make(p_lamel, half_angle_azimuth=4.0)  # > π
        with pytest.raises(ValueError, match="half_angle_azimuth"):
            _make(p_lamel, half_angle_azimuth=0.0)

    def test_out_of_range_linear_half_angle_raises(self, p_lamel):
        with pytest.raises(ValueError, match="half_angle_linear"):
            _make(p_lamel, half_angle_linear=2.0)  # > π/2

    def test_out_of_range_fan_spread_raises(self, p_lamel):
        with pytest.raises(ValueError, match="fan_spread"):
            _make(p_lamel, fan_spread=2.0)  # ≥ π/2

    def test_patch_outside_box_raises(self):
        # A tiny box (relative to R_cell) makes the shell patch fall outside the
        # box half-edge → containment guard fires (no silent gate-loosening).
        # L_box = 1.2·R_cell so |coord| up to R_cell > L_box/2 = 0.6·R_cell.
        p_small = _resolved(n_WAVE=20, L_box=1.2 * R_CELL)
        with pytest.raises(ValueError, match="box half-edge"):
            generate_polarized_patch_layout(
                p_small, wave_tag_start=0, R_cell=R_CELL,
                rng=np.random.default_rng(0),
            )
