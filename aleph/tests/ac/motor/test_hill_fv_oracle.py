"""Self-test of the Hill-1938 force-velocity oracle (pure NumPy — no Warp/CUDA).

Validates the analytic ground truth before it gates the Warp per-head stepping kernel. Acceptance-layer
math, Warp-only-contract exempt by construction.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.motor.hill_fv_analytic import (
    hill_force,
    hill_velocity,
    linear_velocity,
    power_per_head,
    step_work,
)

# illustrative (GAP/provisional) params — the SHAPE gates are magnitude-independent (params_i0b3.yaml)
V0 = 0.2          # um/s (GAP: 0.12 vs 0.2)
FS = 2.0          # pN   (GAP: 0.5 vs 2.0)


@pytest.mark.parametrize("kappa", [0.25, 0.5, 1.0, 10.0, 1e6])
def test_anchors_hold_for_all_kappa(kappa: float) -> None:
    """v(0)=v0 (unloaded) and v(F_s)=0 (stall) hold for every Hill curvature kappa."""
    assert float(hill_velocity(0.0, V0, FS, kappa)) == pytest.approx(V0, rel=1e-12)
    assert float(hill_velocity(FS, V0, FS, kappa)) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("kappa", [0.25, 0.5, 1.0, 5.0])
def test_forward_inverse_roundtrip(kappa: float) -> None:
    """hill_force is the exact inverse of hill_velocity on [0, F_s] x [0, v0]."""
    f = np.linspace(0.0, FS, 50)
    v = hill_velocity(f, V0, FS, kappa, clamp=False)
    f_rt = hill_force(v, V0, FS, kappa)
    assert np.allclose(f_rt, f, atol=1e-12)


def test_hill_1938_reference_point() -> None:
    """kappa=0.25 (muscle) gives v/v0 = 1/6 at F = F_s/2 — the Hill hyperbola, not the linear 1/2."""
    v = float(hill_velocity(FS / 2.0, V0, FS, kappa=0.25))
    assert v / V0 == pytest.approx(1.0 / 6.0, rel=1e-12)
    # the linear (kappa->inf) law would give 1/2 there — a clear, distinct value
    v_lin = float(linear_velocity(FS / 2.0, V0, FS))
    assert v_lin / V0 == pytest.approx(0.5, rel=1e-12)


def test_linear_is_kappa_to_infinity_limit() -> None:
    """kappa -> inf reproduces the PI-2026-07-07 linear law v0(1-F/F_s) to tolerance."""
    f = np.linspace(0.0, FS, 41)
    assert np.allclose(hill_velocity(f, V0, FS, kappa=1e8), linear_velocity(f, V0, FS), atol=1e-6)


def test_velocity_monotone_decreasing_and_convex_in_kappa() -> None:
    """v strictly decreases with load; and for interior F, v increases with kappa toward the linear chord."""
    f = np.linspace(0.0, FS, 60)
    v = hill_velocity(f, V0, FS, kappa=0.5, clamp=False)
    assert np.all(np.diff(v) < 0.0)  # strictly decreasing in F
    # finite-kappa Hill lies BELOW the linear chord; larger kappa -> closer to it (from below)
    fmid = FS / 2.0
    v_lo = float(hill_velocity(fmid, V0, FS, kappa=0.25))
    v_mid = float(hill_velocity(fmid, V0, FS, kappa=1.0))
    v_hi = float(hill_velocity(fmid, V0, FS, kappa=100.0))
    v_lin = float(linear_velocity(fmid, V0, FS))
    assert v_lo < v_mid < v_hi < v_lin


def test_clamp_prevents_active_lengthening() -> None:
    """Past stall (F > F_s) the clamped velocity is 0 (no active lengthening); unclamped goes negative."""
    assert float(hill_velocity(2.0 * FS, V0, FS, kappa=0.5, clamp=True)) == 0.0
    assert float(hill_velocity(2.0 * FS, V0, FS, kappa=0.5, clamp=False)) < 0.0


def test_force_free_at_v0() -> None:
    """Unloaded (v = v0) => F = 0: the minifilament exerts ~0 contractile force at free-running speed."""
    assert float(hill_force(V0, V0, FS, kappa=0.5)) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("kappa", [0.25, 0.5, 1.0, 100.0])
def test_power_sign_and_single_max(kappa: float) -> None:
    """P(F) = F v(F) >= 0 on [0, F_s], P(0) = P(F_s) = 0, and there is a single interior maximum."""
    f = np.linspace(0.0, FS, 501)
    p = power_per_head(f, V0, FS, kappa)
    assert np.all(p >= -1e-15)
    assert p[0] == pytest.approx(0.0, abs=1e-15)
    assert p[-1] == pytest.approx(0.0, abs=1e-12)
    # exactly one interior sign change of the discrete derivative (unimodal)
    dp = np.diff(p)
    sign_changes = int(np.sum(np.diff(np.sign(dp)) != 0))
    assert sign_changes == 1


def test_step_work_sign_and_efficiency() -> None:
    """Work per step W = F d_step >= 0 below stall and eta = W/dG_ATP < 1 (2nd law), for reference magnitudes."""
    d_step = 0.008   # um (~8 nm working stroke; I0-B3 GAP reference)
    dG_atp = 0.09    # pN*um (~20 kBT; I0-B3 GAP reference)
    f = np.linspace(0.0, FS, 50)
    w = step_work(f, d_step)
    assert np.all(w >= 0.0)
    assert np.all(w < dG_atp)  # eta = W/dG < 1 across the whole sub-stall range at these reference values
