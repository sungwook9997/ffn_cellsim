"""D0 sanity-gate pytest wrapper (Lane D). CUDA-only — skipped when no CUDA device is present.

These drive the same gate functions the CUDA driver runs, on a fresh rig each, and assert PASS.
Run on the A5000 (gbook):  pytest aleph/tests/ac/lane_d/test_d0_gates.py
"""
from __future__ import annotations

import pytest
import warp as wp

wp.init()
_CUDA = next((d for d in wp.get_devices() if d.is_cuda), None)
pytestmark = pytest.mark.skipif(_CUDA is None, reason="Lane D D0 gates require a CUDA device (A5000+)")


@pytest.fixture(scope="module")
def dev():
    return str(_CUDA)


def _rig(dev):
    from aleph.archive.lane_d.membrane_cytosol_rig import MembraneCytosolRig
    return MembraneCytosolRig(device=dev)


def test_g01_g02_traction_sign(dev):
    from aleph.archive.lane_d import gates_d0
    r = gates_d0.gate_traction_sign(_rig(dev))
    assert r["passed"], r


def test_g03_net_force_zero(dev):
    from aleph.archive.lane_d import gates_d0
    r = gates_d0.gate_net_force_zero(_rig(dev))
    assert r["passed"], r
    assert r["n_unresolved_faces"] == 0
    assert r["relative_net_force"] < 1e-3


def test_g04_adjoint_work(dev):
    from aleph.archive.lane_d import gates_d0
    r = gates_d0.gate_adjoint_work(_rig(dev))
    assert r["passed"], r


def test_g05_no_flux_mass_conservation(dev):
    from aleph.archive.lane_d import gates_d0
    r = gates_d0.gate_no_flux_mass_conservation(_rig(dev))
    assert r["passed"], r


def test_g06_reject_bit_exact(dev):
    from aleph.archive.lane_d import gates_d0
    r = gates_d0.gate_reject_bit_exact(_rig(dev))
    assert r["passed"], r
    assert r["pressure_bit_exact"] and r["membrane_bit_exact"]
