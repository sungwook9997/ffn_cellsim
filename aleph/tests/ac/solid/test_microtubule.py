"""ac/ microtubule compartment — analytic buckling/persistence gates (CPU) + bending-force sanity (CUDA).

The MT physics is the reused, ff-verified `ff.microtubule` (Sanity Gate VERIFIED 2026-07-07). These tests
certify the ac/ compartment WRAPPER preserves that ground truth: the analytic Euler buckling load and
persistence length (pure, run on the dev Mac) and — on CUDA — that the reused Cytosim bending force is zero
for a straight tube and restores a bent one toward straight.
"""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.components.solid.microtubule import (
    KAPPA_MT,
    euler_buckling_load,
    persistence_length_um,
)

_CUDA_DEVICE = next((d for d in wp.get_devices() if d.is_cuda), None)


# ── CPU: analytic ground-truth gates (Gittes 1993 / Nédélec 2007; KAPPA_MT = 20 pN·µm²) ──────────────

def test_kappa_mt_in_sourced_band() -> None:
    """KAPPA_MT is the sourced MT flexural rigidity (Gittes 1993 EI≈22, Nédélec 2007 =20 pN·µm²)."""
    assert 18.0 <= KAPPA_MT <= 24.0, "MT EI must sit in the 20-22 pN·µm² sourced band"


def test_euler_buckling_load_matches_analytic() -> None:
    """F_crit = π²·EI/L²: at κ=20, L=5 µm → 7.90 pN (the ff-documented value)."""
    assert euler_buckling_load(20.0, 5.0) == pytest.approx(np.pi**2 * 20.0 / 25.0, rel=1e-12)
    assert euler_buckling_load(20.0, 5.0) == pytest.approx(7.896, abs=1e-2)
    # a longer strut buckles at a lower load (∝ 1/L²)
    assert euler_buckling_load(20.0, 10.0) == pytest.approx(euler_buckling_load(20.0, 5.0) / 4.0, rel=1e-12)


def test_persistence_length_in_mt_band() -> None:
    """L_p = EI/k_BT must land in the MT band ~1–8 mm (i.e. 1000–8000 µm)."""
    lp = persistence_length_um(KAPPA_MT)
    assert 1000.0 <= lp <= 8000.0, f"MT L_p={lp:.0f} µm outside the 1–8 mm band"


# ── CUDA: the reused Cytosim bending force through the ac/ compartment wrapper ────────────────────────

@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
def test_straight_aster_has_zero_bending_force_bent_restores() -> None:
    from aleph.components.solid.microtubule import build_microtubule_compartment
    dev = str(_CUDA_DEVICE)
    mt = build_microtubule_compartment(centre=(0.0, 0.0, 0.0), n_mt=12, reach_R_um=7.4, seg_um=0.5,
                                        node_off=0, device=dev)
    with wp.ScopedDevice(dev):
        pos = wp.array(mt.pos0, dtype=wp.vec3d, device=dev)
        f = wp.zeros(mt.n_nodes, dtype=wp.vec3d, device=dev)
        mt.accumulate(pos, f)
        f_straight = np.linalg.norm(f.numpy(), axis=1).max()
        # radial arms are straight ⇒ zero discrete curvature ⇒ ~zero bending force
        assert f_straight < 1e-6, f"straight MT aster must have ~zero bending force, got {f_straight:.3g} pN"

        # kink one interior node of the first arm off-axis → the reused kernel must produce a restoring force
        p = mt.pos0.copy()
        mid = mt.n_nodes // (2 * mt.n_mt) + 1
        p[mid] = p[mid] + np.array([0.3, 0.0, 0.0])
        posb = wp.array(p, dtype=wp.vec3d, device=dev)
        fb = wp.zeros(mt.n_nodes, dtype=wp.vec3d, device=dev)
        mt.accumulate(posb, fb)
        fb_h = fb.numpy()
        assert np.linalg.norm(fb_h, axis=1).max() > 1e-3, "a bent MT must feel a restoring bending force"
        # the kinked node is pushed back toward the arm axis (−x), i.e. opposes the +x kink
        assert fb_h[mid][0] < 0.0, "bending force must restore the kinked node toward straight"


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
def test_compartment_builds_at_native_scale() -> None:
    from aleph.components.solid.microtubule import build_microtubule_compartment
    dev = str(_CUDA_DEVICE)
    mt = build_microtubule_compartment(n_mt=300, reach_R_um=7.4, device=dev)   # proxy epithelial MT count
    assert mt.n_mt == 300 and mt.n_nodes > 0 and mt.n_tri > 0
    assert mt.euler_buckling_load_pn(5.0) == pytest.approx(7.896, abs=1e-2)
