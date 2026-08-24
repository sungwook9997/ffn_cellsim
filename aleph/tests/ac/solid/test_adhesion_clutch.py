"""ac/ ECM adhesion-clutch — α2β1 slip off-rate (CPU) + two-sided traction (CUDA).

Certifies the ac/ wrapper of the ff-audited clutch: the sourced α2β1-collagen slip off-rate (Attwood 2013)
is monotone in load, and on CUDA the engaged clutch pulls a basal actin node toward its collagen node while
an unbound clutch and a rest-length clutch exert nothing.
"""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.components.solid.adhesion_clutch import (
    A2B1_KOFF0_PER_S,
    CLUTCH_K_INT_PN_UM,
    attwood_slip_off_rate,
)

_CUDA_DEVICE = next((d for d in wp.get_devices() if d.is_cuda), None)


# ── CPU: sourced clutch params (Chan-Odde κ_c; Attwood α2β1 slip) ────────────────────────────────────

def test_clutch_stiffness_in_sourced_band() -> None:
    """k_int = 0.8 pN/nm = 800 pN/µm (Chan-Odde; literature band 0.5-5 pN/nm)."""
    assert 500.0 <= CLUTCH_K_INT_PN_UM <= 5000.0


def test_a2b1_slip_off_rate_is_bell_monotone() -> None:
    """k_off(0) = k_off0 = 0.44 s⁻¹ (Attwood); load increases the off-rate (slip bond)."""
    assert attwood_slip_off_rate(0.0) == pytest.approx(A2B1_KOFF0_PER_S, rel=1e-9)
    rates = attwood_slip_off_rate([0.0, 3.0, 6.0, 12.0])
    assert np.all(np.diff(rates) > 0.0), "slip off-rate must increase with force"
    # f_β = k_BT/x_β ≈ 5.9 pN → at F≈5.9 pN the rate ≈ e·k_off0
    assert attwood_slip_off_rate(5.9) == pytest.approx(A2B1_KOFF0_PER_S * np.e, rel=0.05)


# ── CUDA: two-sided clutch traction through the ac/ wrapper ───────────────────────────────────────────

@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
def test_engaged_clutch_pulls_actin_toward_ecm() -> None:
    from aleph.components.solid.adhesion_clutch import build_adhesion_clutch_compartment
    dev = str(_CUDA_DEVICE)
    # basal actin nodes just above a planar collagen substrate at z=0
    grid = np.stack(np.meshgrid(np.linspace(-2, 2, 6), np.linspace(-2, 2, 6)), -1).reshape(-1, 2)
    ecm_pos = np.column_stack([grid, np.zeros(len(grid))])                # collagen nodes at z=0
    actin_pos = np.column_stack([grid, np.full(len(grid), 0.05)])         # actin 0.05 µm above → within capture
    actin_gidx = np.arange(len(actin_pos))
    clutch = build_adhesion_clutch_compartment(actin_pos=actin_pos, actin_global_idx=actin_gidx,
                                               ecm_pos=ecm_pos, capture_um=0.1, rest_um=0.0, device=dev)
    assert clutch.n_clutch == len(actin_pos)

    n = len(actin_pos)
    with wp.ScopedDevice(dev):
        pos = wp.array(actin_pos, dtype=wp.vec3d, device=dev)
        f = wp.zeros(n, dtype=wp.vec3d, device=dev)
        clutch.accumulate(pos, f)
        fh = f.numpy()
        # engaged clutches (rest=0) pull each actin node toward its z=0 collagen node → net −z traction
        assert np.all(fh[:, 2] < 0.0), "engaged clutch must pull the basal actin toward the ECM (−z)"
        assert np.linalg.norm(fh, axis=1).max() > 1e-3

        # unbind every clutch → zero force
        clutch.ecm_node_d.assign(np.full(n, -1, np.int32))
        f2 = wp.zeros(n, dtype=wp.vec3d, device=dev)
        clutch.accumulate(pos, f2)
        assert np.abs(f2.numpy()).max() < 1e-12, "an unbound clutch must exert nothing"


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
def test_clutch_out_of_reach_is_unbound() -> None:
    from aleph.components.solid.adhesion_clutch import build_adhesion_clutch_compartment
    dev = str(_CUDA_DEVICE)
    ecm_pos = np.array([[0.0, 0.0, 0.0]])
    actin_pos = np.array([[0.0, 0.0, 5.0]])          # 5 µm away, capture 0.1 → no bond
    clutch = build_adhesion_clutch_compartment(actin_pos=actin_pos, actin_global_idx=np.array([0]),
                                               ecm_pos=ecm_pos, capture_um=0.1, device=dev)
    assert int(clutch.ecm_node_d.numpy()[0]) < 0
