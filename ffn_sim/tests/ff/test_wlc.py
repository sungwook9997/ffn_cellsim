"""Unit tests for the thermal WLC force-extension law (ffn_sim.ff.wlc) — the single source of the WLC law."""
import numpy as np
from ffn_sim.ff import wlc
from ffn_sim.ff.units import KBT

LP, EA = 17.0, 8639.0                       # collagen persistence + axial stiffness
XMAX = wlc.wlc_x_crossover(EA, LP, KBT)


def test_x_crossover_value():
    assert abs(XMAX - 0.99756) < 1e-4       # f'_WLC(x_max)=EA, fixed by (EA,Lp,kBT)


def test_small_x_linearizes_to_entropic_spring():
    x = np.array([1e-4, 1e-3])
    F = wlc.wlc_tension_np(x * 10.0, np.full(2, 10.0), LP, EA, KBT, XMAX)
    assert np.allclose(F / x, 3 * KBT / (2 * LP), rtol=1e-2)


def test_wall_is_C1_continuous():
    Lc, eps = 1.0, 1e-6
    Fm = wlc.wlc_tension_np([XMAX * Lc], [Lc], LP, EA, KBT, XMAX)[0]
    below = wlc.wlc_tension_np([(XMAX - eps) * Lc], [Lc], LP, EA, KBT, XMAX)[0]
    above = wlc.wlc_tension_np([(XMAX + eps) * Lc], [Lc], LP, EA, KBT, XMAX)[0]
    assert abs((Fm - below) / eps - (above - Fm) / eps) / EA < 1e-2   # slope match ≈ EA/Lc


def test_energy_is_integral_of_tension():
    L, Lc, eps = np.array([0.5, 0.8, 0.95]), np.ones(3), 1e-6
    dU = (wlc.wlc_energy_np(L + eps, Lc, LP, EA, KBT, XMAX) - wlc.wlc_energy_np(L, Lc, LP, EA, KBT, XMAX)) / eps
    assert np.allclose(dU, wlc.wlc_tension_np(L, Lc, LP, EA, KBT, XMAX), rtol=1e-3)


def test_mesh_scales_as_inverse_sqrt_conc():
    # xi = (V/L_total)^1/2 ; doubling L_total (∝c) should drop xi by sqrt(2)
    box_lo, box_hi = [0, 0, 0], [10, 10, 10]
    xi1 = wlc.geometric_mesh_um(np.full(100, 0.5), box_lo, box_hi, dim=3)
    xi2 = wlc.geometric_mesh_um(np.full(200, 0.5), box_lo, box_hi, dim=3)
    assert abs(xi1 / xi2 - np.sqrt(2.0)) < 1e-6


def test_baseline_operating_point():
    xi = 2.0
    Lc = wlc.segment_contour_lengths(np.array([0.5]), LP, xi)
    x0 = 0.5 / Lc[0]
    assert abs(x0 - 1.0 / (1.0 + xi / (6 * LP))) < 1e-9   # physiological resting extension
    f0 = wlc.wlc_tension_np([0.5], Lc, LP, EA, KBT, XMAX)[0]
    assert 0.05 < f0 < 1.0                                # small balanced pre-tension (not floppy/pre-stressed)
