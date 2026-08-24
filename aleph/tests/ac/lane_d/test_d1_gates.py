"""D1 sanity-gate pytest wrapper (Lane D). CUDA-only — skipped when no CUDA device is present.

Run on the A5000 (gbook):  pytest aleph/tests/ac/lane_d/test_d1_gates.py
"""
from __future__ import annotations

import pytest
import warp as wp

wp.init()
_CUDA = next((d for d in wp.get_devices() if d.is_cuda), None)
pytestmark = pytest.mark.skipif(_CUDA is None, reason="Lane D D1 gates require a CUDA device (A5000+)")


@pytest.fixture(scope="module")
def dev():
    return str(_CUDA)


def _rig(dev):
    from aleph.archive.lane_d.force_path_rig import ForcePathRig
    return ForcePathRig(device=dev).relax_and_return()


def test_g11_force_ledger(dev):
    from aleph.archive.lane_d import gates_d1
    r = gates_d1.gate_force_ledger(_rig(dev))
    assert r["passed"], r
    assert r["net_internal_force_rel"] < 1e-9
    assert r["nucleus_reaction_pn"] > 1.0


def test_g12_mt_compression_sign(dev):
    from aleph.archive.lane_d import gates_d1
    r = gates_d1.gate_mt_compression_sign(dev)
    assert r["passed"], r
    assert r["restoring_force_pn"] > 0


def test_g13_if_tension_sign(dev):
    from aleph.archive.lane_d import gates_d1
    r = gates_d1.gate_if_tension_sign(_rig(dev))
    assert r["passed"], r


def test_g14_linc_equal_opposite(dev):
    from aleph.archive.lane_d import gates_d1
    r = gates_d1.gate_linc_equal_opposite(dev)
    assert r["passed"], r
    assert r["equal_opposite_residual_pn"] < 1e-9
    assert r["compression_force_pn"] < 1e-12   # tension-only


def test_g15_linc_detach_path(dev):
    from aleph.archive.lane_d import gates_d1
    r = gates_d1.gate_linc_detach_path(_rig(dev))
    assert r["passed"], r
    assert r["others_bit_exact"]


def test_g16_generation_validation(dev):
    from aleph.archive.lane_d import gates_d1
    r = gates_d1.gate_generation_validation(_rig(dev))
    assert r["passed"], r


def test_g17_reject_topology_bit_exact(dev):
    from aleph.archive.lane_d import gates_d1
    r = gates_d1.gate_reject_topology_bit_exact(_rig(dev))
    assert r["passed"], r


def test_g18_work_conservation(dev):
    from aleph.archive.lane_d import gates_d1
    r = gates_d1.gate_work_conservation(_rig(dev))
    assert r["passed"], r
    assert r["actuator_work_pn_um"] > 0
