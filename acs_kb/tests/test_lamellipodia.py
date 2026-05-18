"""Phase 1 Unit 3.1 — Unit 3.2 lamellipodia structural-prep tests.

Five KU-driven contract checks. The acceptance physics check
(free-protrusion v_p ≈ 30 nm/s) runs against the KU-3.6 closed-form
under the zero-load mock; the time-stepping body is stubbed (Unit 3.2
will implement) and the corresponding test asserts the stub raises.
"""

from __future__ import annotations

import numpy as np
import pytest

from acs_kb.cell.lamellipodia import (
    Lamellipodium,
    apply_tip_load_from_clutch,
    advance_lamellipodium,
    brownian_ratchet_velocity,
    zero_load_mock,
    DELTA_ACTIN_MONOMER_M,
)


# ---------------------------------------------------------------- #
# 1. KU-3.6 acceptance: free protrusion v_p ≈ 30 nm/s              #
# ---------------------------------------------------------------- #

def test_free_protrusion_velocity_matches_KU36_validation():
    """KU-3.6 VALIDATION: F=0, c_G physiological ⇒ v_p ≈ 30 nm/s ± 20 %."""
    v_p_m_per_s = brownian_ratchet_velocity(
        tip_force_N=0.0,
        g_actin_concentration_uM=1.0,
    )
    v_p_nm_per_s = v_p_m_per_s * 1e9
    assert 24.0 < v_p_nm_per_s < 36.0, (
        f"KU-3.6 free protrusion: expected ≈ 30 nm/s, got {v_p_nm_per_s:.2f}"
    )


def test_velocity_decreases_with_tip_force():
    """KU-3.6 monotonic load suppression: v_p strictly decreasing in F."""
    Fs = [0.0, 1e-12, 3e-12, 5e-12, 1e-11]
    vs = [brownian_ratchet_velocity(tip_force_N=F, g_actin_concentration_uM=1.0)
          for F in Fs]
    diffs = np.diff(vs)
    assert np.all(diffs < 0), (
        f"v_p must decrease monotonically with F; got v_p={vs} diffs={diffs}"
    )


def test_high_load_drives_depolymerization():
    """KU-3.6 stall: F·δ ≫ k_BT collapses on-rate; v_p → -δ k_off (m/s)."""
    v_p = brownian_ratchet_velocity(
        tip_force_N=1e-9,                     # 1 nN per filament, way past stall
        g_actin_concentration_uM=1.0,
    )
    expected_off = -DELTA_ACTIN_MONOMER_M * 1.0    # -δ · k_off  (m/s)
    assert v_p == pytest.approx(expected_off, rel=1e-3), (
        f"high-load limit should give v_p = -δ k_off = {expected_off}, got {v_p}"
    )


def test_input_validation():
    with pytest.raises(ValueError):
        brownian_ratchet_velocity(tip_force_N=0.0, delta_m=0.0)
    with pytest.raises(ValueError):
        brownian_ratchet_velocity(tip_force_N=0.0, g_actin_concentration_uM=0.0)


# ---------------------------------------------------------------- #
# 2. Lamellipodium dataclass contract                              #
# ---------------------------------------------------------------- #

def test_lamellipodium_construction_at_cortex_bead():
    lam = Lamellipodium.at_cortex_bead(
        lamellipodium_id=0,
        cell_id=0,
        bead_position=np.array([1e-5, 0.0]),
        outward_normal=np.array([1.0, 0.0]),
        n_filaments_at_tip=8,
        width_um=5.0,
    )
    assert lam.id == 0 and lam.cell_id == 0
    np.testing.assert_allclose(lam.base_position, lam.tip_position)
    assert lam.protrusion_length_m() == 0.0
    assert lam.polarity[0] == pytest.approx(1.0)
    assert lam.polarity[1] == pytest.approx(0.0)


def test_apply_tip_load_from_zero_load_mock():
    lam = Lamellipodium.at_cortex_bead(
        lamellipodium_id=1, cell_id=0,
        bead_position=np.array([1e-5, 0.0]),
        outward_normal=np.array([1.0, 0.0]),
    )
    lam_after = apply_tip_load_from_clutch(lam, zero_load_mock)
    assert lam_after is lam
    assert lam.last_tip_force_N == 0.0


# ---------------------------------------------------------------- #
# 3. Time advance is explicitly stubbed                            #
# ---------------------------------------------------------------- #

def test_advance_is_stub_until_unit_32():
    lam = Lamellipodium.at_cortex_bead(
        lamellipodium_id=0, cell_id=0,
        bead_position=np.array([1e-5, 0.0]),
        outward_normal=np.array([1.0, 0.0]),
    )
    with pytest.raises(NotImplementedError):
        advance_lamellipodium(lam, sampler=zero_load_mock, dt=0.01)
