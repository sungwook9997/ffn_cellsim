"""Cytosim parity oracle (ff/cytosim_parity) — Stage 6f.

Validates the FF bending kernel against the INDEPENDENT C++ Cytosim on the bending energy of a fixed
circular arc. SKIPPED when no ``sim`` binary is available (Cytosim is not a repo dependency), so the
suite stays green everywhere; runs the real cross-check when ``CYTOSIM_SIM`` points at a built
binary. Anchors the finding (FF_CYTOSIM_PARITY_2026-06-30): Cytosim is continuum-accurate at every
resolution; FF's interior-triple energy under-counts by exactly (n−2)/(n−1) and converges to it.
"""

import numpy as np
import pytest

from ffn_sim.ff.cytosim_parity import bending_energy_parity, find_cytosim_sim

pytestmark = pytest.mark.skipif(find_cytosim_sim() is None,
                                reason="Cytosim 'sim' binary not available (set CYTOSIM_SIM)")


def test_cytosim_is_continuum_accurate():
    """The oracle matches the analytic κL/2R² at every resolution (incl. coarse n=5)."""
    for r in bending_energy_parity(n_list=(5, 9, 17)):
        assert r["cyto_over_analytic"] == pytest.approx(1.0, abs=2e-3)


def test_ff_end_factor_and_convergence():
    """FF under-counts by the end-factor (n−2)/(n−1): FF·(n−1)/(n−2) recovers the continuum, and
    FF/analytic rises monotonically toward 1 as n grows."""
    rows = bending_energy_parity(n_list=(5, 9, 17, 33))
    for r in rows:
        assert r["ff_end_corrected"] == pytest.approx(r["analytic"], rel=3e-3)
        assert r["ff_over_analytic"] == pytest.approx((r["n"] - 2) / (r["n"] - 1), rel=3e-3)
    ratios = [r["ff_over_analytic"] for r in rows]
    assert ratios == sorted(ratios)                 # monotone convergence to 1
    assert ratios[-1] > ratios[0]


def test_ff_approaches_cytosim_with_resolution():
    """FF→Cytosim as the fiber is refined (the two discretisations agree in the continuum limit)."""
    rows = bending_energy_parity(n_list=(5, 33))
    err_coarse = abs(rows[0]["ff"] - rows[0]["cytosim"])
    err_fine = abs(rows[1]["ff"] - rows[1]["cytosim"])
    assert err_fine < err_coarse
