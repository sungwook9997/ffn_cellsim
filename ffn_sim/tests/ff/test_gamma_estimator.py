"""γ estimator + actin axial tension (ff/gamma_estimator, ff/constraints) — Stage 6d.

Anchors: the constraint-multiplier axial tension is +T under a pure stretch, −T under compression,
0 unloaded, and uniform along the fiber; the method-of-planes kernel reproduces a hand-computed γ
on a controlled single-plane bond set, scales linearly with tension, and ignores non-crossing
elements; unit conversions to mN/m and N/m are correct against the Salbreux band.
"""

import numpy as np
import pytest

from ffn_sim.ff.constraints import segment_axial_tension
from ffn_sim.ff.fiber_network import build_fiber_network
from ffn_sim.ff.gamma_estimator import (
    SALBREUX_BAND_PN_UM,
    gamma_to_mN_per_m,
    gamma_to_N_per_m,
    method_of_planes_gamma,
)

SEG = 0.5


def _straight(n=6):
    nodes = np.stack([np.arange(n) * SEG, np.zeros(n), np.zeros(n)], axis=1).astype(float)
    return build_fiber_network([nodes], kappa=0.073)


def test_axial_tension_stretch_compress_zero():
    net = _straight()
    z = segment_axial_tension(net, np.zeros((6, 3)))
    assert np.allclose(z, 0.0)
    F = np.zeros((6, 3)); F[0, 0] = -3.0; F[-1, 0] = 3.0      # pull apart
    t = segment_axial_tension(net, F)
    assert np.allclose(t, 3.0)                                 # uniform +T tension
    F2 = -F
    assert np.allclose(segment_axial_tension(net, F2), -3.0)   # compression


def test_axial_tension_scales_with_load():
    net = _straight()
    F = np.zeros((6, 3)); F[0, 0] = -1.0; F[-1, 0] = 1.0
    t1 = segment_axial_tension(net, F)
    t2 = segment_axial_tension(net, 2.5 * F)
    assert np.allclose(t2, 2.5 * t1)


def test_gamma_empty_is_zero():
    assert method_of_planes_gamma(np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0), 1.0) == 0.0


def test_gamma_single_plane_handcalc():
    """One x̂-normal plane (n_planes=1), one bond along x crossing x=0: γ = T·|û·n̂|/(2πR)."""
    R = 1.0
    rA = np.array([[-1.0, 0.0, 0.0]]); rB = np.array([[1.0, 0.0, 0.0]])
    T = np.array([5.0])
    g = method_of_planes_gamma(rA, rB, T, R, n_planes=1, centre=np.zeros(3))
    assert g == pytest.approx(5.0 / (2 * np.pi * R))


def test_gamma_scales_linearly_with_tension():
    R = 10.0
    th = np.linspace(0, 2 * np.pi, 40, endpoint=False)
    pts = np.stack([R * np.cos(th), R * np.sin(th), np.zeros(40)], axis=1)
    rA, rB = pts, np.roll(pts, -1, axis=0)
    g1 = method_of_planes_gamma(rA, rB, np.full(40, 1.0), R, n_planes=100, centre=np.zeros(3))
    g3 = method_of_planes_gamma(rA, rB, np.full(40, 3.0), R, n_planes=100, centre=np.zeros(3))
    assert g3 == pytest.approx(3.0 * g1, rel=1e-9)
    assert g1 > 0.0


def test_gamma_ignores_noncrossing_elements():
    """An element entirely on one side of every plane through the centre contributes nothing."""
    R = 1.0
    rA = np.array([[-1.0, 0.0, 0.0], [0.4, 0.4, 0.0]])     # 2nd bond tiny, same side, off-centre
    rB = np.array([[1.0, 0.0, 0.0], [0.45, 0.45, 0.0]])
    T = np.array([5.0, 100.0])
    g = method_of_planes_gamma(rA, rB, T, R, n_planes=1, centre=np.zeros(3))
    # plane normal x̂: bond 2 has both endpoints at x>0 → no crossing, ignored despite huge T
    assert g == pytest.approx(5.0 / (2 * np.pi * R))


def test_unit_conversions_and_band():
    assert gamma_to_mN_per_m(1.0) == pytest.approx(1e-3)
    assert gamma_to_N_per_m(1.0) == pytest.approx(1e-6)
    # Salbreux band 0.35–0.65 pN/µm = 0.35–0.65 mN/m
    lo, hi = SALBREUX_BAND_PN_UM
    assert gamma_to_mN_per_m(lo) == pytest.approx(0.35e-3)
    assert gamma_to_mN_per_m(hi) == pytest.approx(0.65e-3)
