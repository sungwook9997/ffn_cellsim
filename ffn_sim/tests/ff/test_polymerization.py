"""Actin barbed-end polymerization ratchet (ff/polymerization_warp) — piece 1/5 of active movement.

Validates the kernel against the ANALYTIC Mogilner-Oster force-velocity ground truth (oracle-is-crosscheck
rule) + the lit-anchored constants (Pollard 1986 kinetics; kBT/δ e-fold force), + load-sharing + elongation.
"""

import numpy as np
import pytest
import warp as wp

from ffn_sim.ff.polymerization_warp import (
    polymerization_kernel,
    ratchet_velocity_np,
    resolve_polymerization,
)


def test_resolve_lit_constants():
    """Pollard 1986 kinetics → v0=(k_on·M−k_off)·δ_full; e-fold force kBT/δ ≈ 3.2 pN (Mogilner-Oster)."""
    sp = resolve_polymerization(G_actin_uM=20.0)
    assert sp.v0_um_s == pytest.approx((11.6 * 20.0 - 1.4) * 2.7e-3)     # 0.6226 µm/s
    assert sp.delta_um == pytest.approx(1.35e-3)                          # half-monomer
    assert sp.efold_force_pN == pytest.approx(4.28e-3 / 1.35e-3, rel=1e-9)  # ~3.17 pN
    with pytest.raises(ValueError):
        resolve_polymerization(G_actin_uM=0.05)                           # below the physiological band


def _run_vf(loads, sp, device="cpu"):
    pos = wp.array(np.array([[0, 0, 0], [0.5, 0, 0]], dtype=np.float64), dtype=wp.vec3d, device=device)
    out = []
    for f in loads:
        frc = wp.array(np.array([[0, 0, 0], [-f, 0, 0]], dtype=np.float64), dtype=wp.vec3d, device=device)
        vo = wp.zeros(1, dtype=wp.float64, device=device)
        sr = wp.array(np.array([0.5]), dtype=wp.float64, device=device)
        wp.launch(polymerization_kernel, dim=1, inputs=[
            pos, frc, wp.array([1], dtype=wp.int32, device=device), wp.array([0], dtype=wp.int32, device=device),
            wp.array([0], dtype=wp.int32, device=device), wp.array([1.0], dtype=wp.float64, device=device),
            sr, wp.float64(sp.v0_um_s), wp.float64(sp.delta_um), wp.float64(sp.kBT_pN_um),
            wp.float64(1e-3), vo], device=device)
        out.append((float(vo.numpy()[0]), float(sr.numpy()[0])))
    return out


def test_kernel_matches_analytic_force_velocity():
    """The kernel's v(f) IS the analytic ratchet v0·exp(−fδ/kBT) — exact (the always-on ground truth)."""
    sp = resolve_polymerization(G_actin_uM=20.0)
    loads = np.array([0.0, 1.0, 3.17, 7.0, 15.0])
    vk = np.array([r[0] for r in _run_vf(loads, sp)])
    assert np.allclose(vk, ratchet_velocity_np(loads, sp), rtol=1e-9)
    assert vk[0] == pytest.approx(sp.v0_um_s)                             # unloaded → v0
    assert vk[1] < vk[0] and vk[-1] < vk[1]                               # monotone-decreasing with load
    assert ratchet_velocity_np(sp.efold_force_pN, sp) == pytest.approx(sp.v0_um_s / np.e, rel=1e-9)  # e-fold


def test_load_sharing_and_elongation():
    """Load shares as f/N across barbed ends; the barbed segment elongates by v(f)·dt."""
    sp = resolve_polymerization(G_actin_uM=20.0)
    d = "cpu"
    pos = wp.array(np.array([[0, 0, 0], [0.5, 0, 0]], dtype=np.float64), dtype=wp.vec3d, device=d)
    frc = wp.array(np.array([[0, 0, 0], [-10.0, 0, 0]], dtype=np.float64), dtype=wp.vec3d, device=d)
    vo = wp.zeros(1, dtype=wp.float64, device=d); sr = wp.array(np.array([0.5]), dtype=wp.float64, device=d)
    wp.launch(polymerization_kernel, dim=1, inputs=[
        pos, frc, wp.array([1], dtype=wp.int32, device=d), wp.array([0], dtype=wp.int32, device=d),
        wp.array([0], dtype=wp.int32, device=d), wp.array([10.0], dtype=wp.float64, device=d),
        sr, wp.float64(sp.v0_um_s), wp.float64(sp.delta_um), wp.float64(sp.kBT_pN_um), wp.float64(1e-3), vo], device=d)
    v = float(vo.numpy()[0])
    assert v == pytest.approx(ratchet_velocity_np(1.0, sp), rel=1e-9)     # f=10/N=10 ≡ f=1
    assert sr.numpy()[0] == pytest.approx(0.5 + v * 1e-3, rel=1e-9)       # elongated by v·dt
