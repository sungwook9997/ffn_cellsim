"""Self-test of the chromatin WLC polymer-net oracle (pure NumPy — no Warp/CUDA).

The internal-polymer-net primitive: Marko–Siggia WLC tension is monotone, force = dE/dx (FD sign
arbiter), diverges as x→L (finite extensibility → never interpenetrates), and the small-extension limit
is the entropic spring constant k=3kT/(2ℓ_p L).
"""
from __future__ import annotations

import numpy as np
import pytest

from aleph.components.nucleus.chromatin_analytic import (
    KBT_PN_UM_310K,
    wlc_energy,
    wlc_small_strain_stiffness,
    wlc_tension,
)

L, LP = 1.0, 0.05  # µm (contour, persistence) — structural test values


def test_wlc_force_matches_fd_energy() -> None:
    """f(x) = dE/dx to tight tolerance (the sign arbiter): the tension is the energy gradient."""
    x = np.linspace(0.01, 0.9, 60) * L
    f_ana = wlc_tension(x, L, LP)
    eps = 1e-7
    f_fd = (wlc_energy(x + eps, L, LP) - wlc_energy(x - eps, L, LP)) / (2 * eps)
    assert np.max(np.abs(f_ana - f_fd)) < 1e-6 * np.max(np.abs(f_ana))


def test_wlc_monotone_and_diverges() -> None:
    """Tension is monotone increasing and blows up as x→L (finite extensibility)."""
    x = np.linspace(0.01, 0.99, 200) * L
    f = wlc_tension(x, L, LP)
    assert np.all(np.diff(f) > 0.0)
    assert f[-1] > 50.0 * f[0]                                   # steep near the contour length


def test_wlc_small_strain_stiffness() -> None:
    """Small-extension limit f/x → k = 3kT/(2 ℓ_p L) (the chromatin soft stiffness feeding k_chrom)."""
    k_ref = wlc_small_strain_stiffness(L, LP)
    x_small = 1e-4 * L
    k_num = float(wlc_tension(np.array([x_small]), L, LP)[0]) / x_small
    assert k_num == pytest.approx(k_ref, rel=1e-3)
    assert k_ref == pytest.approx(3.0 * KBT_PN_UM_310K / (2.0 * LP * L))


def test_wlc_force_free_at_zero_extension() -> None:
    """f(0)=0 and E(0)=0 (force-free, zero-energy reference)."""
    assert float(wlc_tension(np.array([0.0]), L, LP)[0]) == pytest.approx(0.0, abs=1e-12)
    assert float(wlc_energy(np.array([0.0]), L, LP)[0]) == pytest.approx(0.0, abs=1e-12)
