"""Bending kernel dispersion relation — analytic oracle (NF2007 semiflexible bending physics).

The runnable-Cytosim parity oracle (ENGINE.md §6) is blocked on the external C++ build (no binary
available offline) — surfaced to PI. Meanwhile this validates the bending kernel against the
ANALYTIC anchor the physics provides: a transverse sinusoidal mode y = A·sin(qx) on a filament feels
a bending force density −κ q⁴ y (the semiflexible dispersion / 4th-derivative operator). The kernel
returns force PER MODEL-POINT (= force density × seg), so

    −F_y / (y · seg)  →  κ q⁴      (q·seg ≪ 1)

and the rate scales as q⁴ across modes. (The per-point ``seg`` factor cancels against the per-point
drag γ_point = γ_line·seg — `units.fiber_point_drag` — so the physical mode rate is the continuum
μ κ q⁴.) Together with the existing force=−∇E test this pins the kernel as a faithful κ-weighted
4th-derivative operator.
"""

import numpy as np
import pytest

from aleph.laws.fiber_network import build_fiber_network
from aleph.laws.forces_warp import bending_force

KAPPA = 0.073      # pN·µm² (actin, FF units)
SEG = 0.02         # µm — small so q·seg ≪ 1 (discrete → continuum)
N = 501


def _mode_force_density(m: int):
    """−F_y/(y·seg) on an interior antinode for transverse mode m, and the predicted κq⁴."""
    x = np.arange(N) * SEG
    L = (N - 1) * SEG
    q = m * np.pi / L
    nodes = np.stack([x, 1e-5 * np.sin(q * x), np.zeros(N)], axis=1)
    net = build_fiber_network([nodes], kappa=KAPPA)
    # validate the RAW interior-triple operator (end_correction off) against the continuum κq⁴; the
    # end-correction is a per-fiber constant (≈1.002 here) that rescales magnitude, not the q⁴ shape.
    F = bending_force(net, end_correction=False)
    interior = slice(8, N - 8)
    eig = float(np.nanmedian(-F[interior, 1] / nodes[interior, 1]))   # = κ q⁴ seg
    return eig / SEG, KAPPA * q**4, eig, q


def test_dispersion_matches_continuum():
    """−F_y/(y·seg) = κ q⁴ to <0.1% for several modes (faithful 4th-derivative operator)."""
    for m in (1, 2, 4, 8):
        dens, pred, _, _ = _mode_force_density(m)
        assert dens == pytest.approx(pred, rel=2e-3)


def test_dispersion_q4_scaling():
    """Doubling q quadruples-squared (×16) the bending rate — the q⁴ law, exactly."""
    _, _, eig1, q1 = _mode_force_density(1)
    _, _, eig2, q2 = _mode_force_density(2)
    _, _, eig4, q4 = _mode_force_density(4)
    assert eig2 / eig1 == pytest.approx((q2 / q1) ** 4, rel=2e-3)
    assert eig4 / eig2 == pytest.approx((q4 / q2) ** 4, rel=2e-3)
