"""Compartment/network parity gate: Warp harmonic bond vs HOOMD md.bond.Harmonic.

Reference fixture (``warp_port/fixtures/network_bond_ref.npz``) is the per-bead net
force + energy from HOOMD's native ``md.bond.Harmonic`` (see
``generate_network_fixture.py``) — the cortex axial-spring network force. The Warp
kernel is graded against THAT committed HOOMD output (guard-rail 2).

Tolerance-class: forces are accumulated per bond via ``wp.atomic_add`` (shared
beads), so order can differ from HOOMD's C++ accumulation; in a chain each bead has
≤2 bonds so the effect is negligible (achieved ~5e-16). Gated at rel < 1e-12.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

FIX = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "warp_port", "fixtures")


def _load(name: str = "network_bond_ref") -> dict:
    path = os.path.join(FIX, f"{name}.npz")
    if not os.path.exists(path):
        pytest.skip(
            f"fixture {path} missing — run the matching "
            "ffn_sim/warp_port/fixtures/generate_*_fixture.py"
        )
    return dict(np.load(path))


def test_harmonic_bond_warp_parity():
    """Warp harmonic-bond force + energy match HOOMD md.bond.Harmonic."""
    pytest.importorskip("warp")
    from ffn_sim.warp_port.network_warp import run_harmonic_bond_warp

    fx = _load()
    got = run_harmonic_bond_warp(
        pos=fx["pos"], bonds=fx["bonds"], k=float(fx["k"]), r0=float(fx["r0"]),
        box_L=(float(fx["L"]),) * 3, device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    Uscale = float(np.abs(fx["ref_energy"]).max()) + 1e-30
    dF = float(np.abs(got["force"] - fx["ref_force"]).max())
    dU = float(np.abs(got["energy"] - fx["ref_energy"]).max())
    assert dF / Fscale < 1e-12, f"bond force rel error too large: {dF/Fscale:.3e}"
    assert dU / Uscale < 1e-12, f"bond energy rel error too large: {dU/Uscale:.3e}"


def test_harmonic_angle_warp_parity():
    """Warp harmonic-angle (bending) force + energy match HOOMD md.angle.Harmonic."""
    pytest.importorskip("warp")
    from ffn_sim.warp_port.network_warp import run_harmonic_angle_warp

    fx = _load("network_angle_ref")
    got = run_harmonic_angle_warp(
        pos=fx["pos"], angles=fx["angles"], k=float(fx["k"]), t0=float(fx["t0"]),
        box_L=(float(fx["L"]),) * 3, device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    Uscale = float(np.abs(fx["ref_energy"]).max()) + 1e-30
    dF = float(np.abs(got["force"] - fx["ref_force"]).max())
    dU = float(np.abs(got["energy"] - fx["ref_energy"]).max())
    assert dF / Fscale < 1e-12, f"angle force rel error too large: {dF/Fscale:.3e}"
    assert dU / Uscale < 1e-12, f"angle energy rel error too large: {dU/Uscale:.3e}"
