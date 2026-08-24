"""Self-test of the elastic Brownian-ratchet barbed-end force-velocity oracle (pure NumPy — no Warp/CUDA)."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.weave.branch_angle import KT_310K_PN_UM
from aleph.components.weave.membrane_ratchet import (
    ACTIN_MONOMER_DELTA_UM,
    ensemble_protrusion_velocity,
    free_polymerization_velocity,
    per_filament_load,
    ratchet_velocity,
    stall_force,
)

# Representative I0-B4 GAP magnitudes used ONLY as oracle sweep variables (never a sourced value).
K_ON_C = 30.0   # 1/s
K_OFF = 3.0     # 1/s


def test_free_velocity_is_delta_times_net_rate() -> None:
    """V(0) = delta*(k_on_c - k_off): the load-free elongation speed."""
    v0 = free_polymerization_velocity(K_ON_C, K_OFF)
    assert v0 == pytest.approx(ACTIN_MONOMER_DELTA_UM * (K_ON_C - K_OFF))
    assert ratchet_velocity(0.0, K_ON_C, K_OFF) == pytest.approx(v0)


def test_force_velocity_is_strictly_decreasing_in_load() -> None:
    """Raising the opposing load monotonically slows protrusion (a force-velocity relation, never increasing)."""
    loads = np.linspace(0.0, 20.0, 50)
    v = ratchet_velocity(loads, K_ON_C, K_OFF)
    assert np.all(np.diff(v) < 0.0)


def test_stall_force_zeroes_the_velocity() -> None:
    """V(f_s) = 0 at f_s = (kT/delta)*ln(k_on_c/k_off); the closed form matches the ratchet root."""
    fs = stall_force(K_ON_C, K_OFF)
    assert fs == pytest.approx((KT_310K_PN_UM / ACTIN_MONOMER_DELTA_UM) * np.log(K_ON_C / K_OFF))
    assert float(ratchet_velocity(fs, K_ON_C, K_OFF)) == pytest.approx(0.0, abs=1e-12)
    assert float(ratchet_velocity(fs * 0.99, K_ON_C, K_OFF)) > 0.0     # below stall: protruding
    assert float(ratchet_velocity(fs * 1.01, K_ON_C, K_OFF)) < 0.0     # above stall: receding


def test_super_stall_floor_is_minus_delta_koff() -> None:
    """As f -> inf the ratchet argument -> 0, so V -> -delta*k_off (super-stall depolymerization, not clamped)."""
    v_huge = float(ratchet_velocity(1.0e6, K_ON_C, K_OFF))
    assert v_huge == pytest.approx(-ACTIN_MONOMER_DELTA_UM * K_OFF, rel=1e-9)


def test_stall_force_undefined_for_a_receding_end() -> None:
    """k_on_c <= k_off has no protruding regime: stall force is raised, never floored to a convenient value."""
    with pytest.raises(ValueError):
        stall_force(2.0, 5.0)
    with pytest.raises(ValueError):
        stall_force(5.0, 5.0)
    # a receding end still has a (negative) free velocity — reported, not clamped
    assert free_polymerization_velocity(2.0, 5.0) < 0.0


def test_ensemble_velocity_increases_with_load_sharing() -> None:
    """More load-sharing barbed ends drop the per-filament load, so protrusion speeds up (Mogilner-Oster)."""
    total = 40.0
    v = [ensemble_protrusion_velocity(total, n, K_ON_C, K_OFF) for n in range(1, 40)]
    assert np.all(np.diff(v) > 0.0)                                    # monotone increasing in N
    assert per_filament_load(total, 8) == pytest.approx(total / 8)
    # the single-end load is heaviest; the many-end load approaches free protrusion
    assert v[0] < free_polymerization_velocity(K_ON_C, K_OFF)
    assert v[-1] == pytest.approx(
        float(ratchet_velocity(total / 39.0, K_ON_C, K_OFF)))


def test_rejects_negative_load_and_bad_rates() -> None:
    """Dimensional / sign guards: negative load, negative rates, non-positive delta/kT are rejected."""
    with pytest.raises(ValueError):
        ratchet_velocity(-1.0, K_ON_C, K_OFF)
    with pytest.raises(ValueError):
        ratchet_velocity(1.0, -1.0, K_OFF)
    with pytest.raises(ValueError):
        ratchet_velocity(1.0, K_ON_C, K_OFF, delta_um=0.0)
    with pytest.raises(ValueError):
        per_filament_load(10.0, 0)
