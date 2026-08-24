"""Oracle-B bleb-growth analysis gate — pure-NumPy discrete reference (CPU; no CUDA device required).

The device growth driver (:mod:`aleph.components.incumbent.bleb_growth`) computes its emergent readouts — cap membership,
patch/control outward displacement, the ΔP-set Laplace critical radius, and the final bulge/grows verdict — on
host reads taken AFTER each accepted outer step (a legitimate post-step diagnostic, not an in-hot-loop D2H).
Those readouts go through four pure functions. This gate exercises exactly those functions on analytic
geometry with hand-computed answers, so the analysis pipeline is validated on the dev Mac while the full native
growth run stays gated on (b) resting convergence and the gbook A5000.

This is ANALYSIS-ONLY validation. It makes no claim about the physical bleb outcome; that requires the native
Warp-CUDA run once resting converges.
"""

from __future__ import annotations

import numpy as np

from aleph.components.incumbent.bleb_growth import (
    bulge_grows,
    cap_membership,
    laplace_state,
    patch_outward_um,
)
from aleph.components.incumbent.bleb_perturbation import cap_cos_half_angle

_Z = np.array([0.0, 0.0, 1.0])


def _sphere(n_theta: int = 200, n_phi: int = 200, radius: float = 7.5, centroid=(1.0, -2.0, 0.5)):
    """An AREA-uniform membrane-like point cloud on a sphere about an off-origin centroid.

    Sampling uniformly in ``cos θ`` (not θ) gives equal area per node, so the cap membership fraction converges
    to the solid-angle fraction ``(1 − cos θ½)/2`` rather than the pole-over-weighted θ-uniform value.
    """
    cos_theta = np.linspace(-1.0 + 1e-3, 1.0 - 1e-3, n_theta)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    cc, pp = np.meshgrid(cos_theta, phi, indexing="ij")
    ss = np.sqrt(1.0 - cc**2)
    unit = np.stack([ss * np.cos(pp), ss * np.sin(pp), cc], axis=-1).reshape(-1, 3)
    return unit * radius + np.asarray(centroid), np.asarray(centroid, float), unit.reshape(-1, 3)


# ── cap membership matches the analytic solid-angle fraction and the device kernel's ≥ test ────────────

def test_cap_membership_matches_solid_angle_fraction() -> None:
    pos, centroid, _ = _sphere()
    half_deg = 30.0
    in_cap, radius = cap_membership(pos, centroid, _Z, cap_cos_half_angle(half_deg))
    # cap solid-angle fraction of a sphere = (1 − cos θ½) / 2.
    expected = (1.0 - np.cos(np.deg2rad(half_deg))) / 2.0
    assert abs(in_cap.mean() - expected) < 5e-3
    assert np.allclose(radius, 7.5, atol=1e-9)          # every node sits at R on the sphere


def test_cap_membership_boundary_is_inclusive() -> None:
    r = np.deg2rad(30.0)
    on_edge = np.array([[np.sin(r), 0.0, np.cos(r)]])   # exactly on the 30° rim, unit sphere at origin
    just_out = np.array([[np.sin(r + 1e-3), 0.0, np.cos(r + 1e-3)]])
    cos_half = cap_cos_half_angle(30.0)
    assert cap_membership(on_edge, np.zeros(3), _Z, cos_half)[0][0]
    assert not cap_membership(just_out, np.zeros(3), _Z, cos_half)[0][0]


# ── outward displacement recovers a known radial bulge over the cap, control unchanged ─────────────────

def test_patch_outward_recovers_known_bulge() -> None:
    pos, centroid, unit = _sphere()
    cos_half = cap_cos_half_angle(30.0)
    in_cap, base_r = cap_membership(pos, centroid, _Z, cos_half)
    base_patch_r = float(base_r[in_cap].mean())
    base_ctrl_r = float(base_r[~in_cap].mean())

    bulge_um = 0.8
    moved = pos.copy()
    moved[in_cap] += unit[in_cap] * bulge_um             # push only the cap radially outward by 0.8 µm
    r_now = np.linalg.norm(moved - centroid, axis=1)

    assert abs(patch_outward_um(r_now, base_patch_r, in_cap) - bulge_um) < 1e-6
    assert abs(patch_outward_um(r_now, base_ctrl_r, ~in_cap) - 0.0) < 1e-9


def test_patch_outward_empty_mask_is_zero() -> None:
    r_now = np.linspace(7.0, 8.0, 50)
    assert patch_outward_um(r_now, 7.5, np.zeros(50, bool)) == 0.0


# ── Laplace state anchors to the committed γ_mem=10 pN/µm, ΔP=40 Pa -> a_crit=0.5 µm ───────────────────

def test_laplace_state_matches_committed_nucleation_anchor() -> None:
    # Matches aleph/outputs/ac/bleb/bleb_70686.json: gamma_mem_pN_um=10, dP=40 -> a_crit_um=0.5.
    a_crit, dp_laplace, expands = laplace_state(gamma_mem_pn_um=10.0, dp_pa=40.0, a_patch_um=3.17)
    assert abs(a_crit - 0.5) < 1e-9
    assert abs(dp_laplace - (2.0 * 10.0 / 3.17)) < 1e-9
    assert expands is True                               # 3.17 µm cap >> 0.5 µm critical radius -> bulges


def test_laplace_state_non_positive_pressure_is_infinite_and_stable() -> None:
    a_crit, _, expands = laplace_state(10.0, 0.0, 3.0)
    assert a_crit == float("inf") and expands is False
    a_crit2, dp2, expands2 = laplace_state(10.0, 40.0, 0.0)
    assert dp2 == float("inf") and expands2 is False and np.isfinite(a_crit2)


# ── final bulge/grows verdict from the post-de-adhesion series ─────────────────────────────────────────

def test_bulge_grows_true_when_patch_expands_beyond_control() -> None:
    patch = np.array([0.0, 0.1, 0.3, 0.6])              # cap expands after de-adhesion (index 0)
    ctrl = np.array([0.0, 0.0, 0.01, 0.0])
    bulge, grows = bulge_grows(patch, ctrl)
    assert abs(bulge - 0.6) < 1e-9 and grows is True


def test_bulge_grows_false_when_patch_does_not_expand() -> None:
    patch = np.array([0.5, 0.4, 0.3, 0.2])              # patch retreats -> not a growing bleb
    ctrl = np.array([0.0, 0.0, 0.0, 0.0])
    _, grows = bulge_grows(patch, ctrl)
    assert grows is False
    # empty series is a well-defined no-growth
    assert bulge_grows([], []) == (0.0, False)
