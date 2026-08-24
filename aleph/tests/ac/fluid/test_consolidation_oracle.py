"""Self-test of the Terzaghi consolidation oracle (pure NumPy — runs anywhere, no Warp/CUDA).

Validates the analytic ground truth itself against textbook values before it is used to gate the
Warp-CUDA I1a solver. This is acceptance-layer math, not simulation, so it is Warp-only-contract
exempt by construction (no simulation kernel here).
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.fluid.consolidation_analytic import (
    terzaghi_degree_of_consolidation,
    terzaghi_excess_pressure,
    time_factor,
)

# Terzaghi's classic U(T_v) table (Terzaghi 1943; any soil-mechanics text).
TEXTBOOK_TV_U = [
    (0.008, 0.10),
    (0.031, 0.20),
    (0.071, 0.30),
    (0.126, 0.40),
    (0.197, 0.50),
    (0.287, 0.60),
    (0.403, 0.70),
    (0.567, 0.80),
    (0.848, 0.90),
]


@pytest.mark.parametrize("t_v, u_expected", TEXTBOOK_TV_U)
def test_degree_of_consolidation_matches_textbook(t_v: float, u_expected: float) -> None:
    """U(T_v) reproduces the tabulated Terzaghi degree of consolidation to <1% absolute."""
    u = float(terzaghi_degree_of_consolidation(t_v))
    assert u == pytest.approx(u_expected, abs=0.01)


def test_degree_boundary_and_monotonicity() -> None:
    """U(0)=0, U(inf)->1, and U is strictly increasing in T_v."""
    tv = np.linspace(0.0, 3.0, 300)
    u = terzaghi_degree_of_consolidation(tv)
    assert u[0] == pytest.approx(0.0, abs=1e-9)
    assert float(terzaghi_degree_of_consolidation(50.0)) == pytest.approx(1.0, abs=1e-6)
    assert np.all(np.diff(u) > 0.0)


def test_half_consolidation_time_factor() -> None:
    """The classic T_v ~= 0.197 at U = 0.5 (and the U<0.6 approximation T_v = pi/4 U^2)."""
    # invert U(T_v) numerically over a fine grid
    tv = np.linspace(1e-4, 1.0, 200_000)
    u = terzaghi_degree_of_consolidation(tv)
    tv_at_half = float(np.interp(0.5, u, tv))
    assert tv_at_half == pytest.approx(0.197, abs=0.005)
    assert tv_at_half == pytest.approx(np.pi / 4 * 0.5**2, abs=0.01)  # small-U closed form


def test_pressure_field_symmetric_and_decays() -> None:
    """Excess pressure is symmetric about the mid-plane (z/H=1) and dissipates to 0 as T_v grows."""
    z = np.linspace(0.0, 2.0, 41)  # z/H across a doubly-drained layer
    u_early = terzaghi_excess_pressure(z, t_v=0.05, u0=1.0)
    # drained faces are ~0, mid-plane is the maximum, field is mirror-symmetric
    assert u_early[0] == pytest.approx(0.0, abs=1e-6)
    assert u_early[-1] == pytest.approx(0.0, abs=1e-6)
    assert np.allclose(u_early, u_early[::-1], atol=1e-9)
    assert np.argmax(u_early) == len(z) // 2
    u_late = terzaghi_excess_pressure(z, t_v=3.0, u0=1.0)
    assert np.max(u_late) < 1e-3


def test_series_truncation_converged() -> None:
    """Adding series terms past the default does not move U or the field past tolerance."""
    for tv in (0.02, 0.2, 1.0):
        u_a = float(terzaghi_degree_of_consolidation(tv, n_terms=200))
        u_b = float(terzaghi_degree_of_consolidation(tv, n_terms=800))
        assert u_a == pytest.approx(u_b, abs=1e-9)


def test_time_factor_roundtrip() -> None:
    """T_v = c_v t / H^2 with the I0-B1 anchor: tau_p (U~0.9) lands near ~1 s for a 7.5 um cell."""
    c_v = 50.0  # um^2/s (KB-3.B3.2 midpoint)
    h = 7.5     # um drainage path ~ cell radius
    # T_v ~ 0.848 at U=0.9 -> t = 0.848 H^2 / c_v
    t90 = 0.848 * h * h / c_v
    tv = float(time_factor(c_v, t90, h))
    assert tv == pytest.approx(0.848, rel=1e-6)
    assert 0.5 < t90 < 2.0  # ~1 s poroelastic relaxation, consistent with Moeendarbary/Charras
