"""Optimization #8 gate: end-to-end differentiability of the K-step DCM loop.

Asserts that reverse-mode autodiff of a scalar loss (final mean radius) w.r.t. the
physical parameters (turgor dP0, edge stiffness k_edge), taken THROUGH the whole
K-step simulation loop, matches a central finite-difference — i.e. the entire
simulation is differentiable, not just a single force evaluation (B4).
"""

from __future__ import annotations

import pytest


def test_diff_loop_grads_match_finite_difference():
    pytest.importorskip("warp")
    from ffn_sim.warp_port.dcm_warp_diff import run_diff

    r = run_diff(subdiv=1, K=40, dt=2.0e-8, device="cpu")
    assert r["rel_err_dP0"] < 1e-6, f"dP0 grad vs FD: {r['rel_err_dP0']:.3e}"
    assert r["rel_err_kedge"] < 1e-6, f"k_edge grad vs FD: {r['rel_err_kedge']:.3e}"
    # gradients must be non-trivial (a zero grad would pass vacuously)
    assert abs(r["grad_dP0_autodiff"]) > 0.0
    assert abs(r["grad_kedge_autodiff"]) > 0.0
    # physical sanity: more turgor -> larger cell; stiffer cortex -> smaller cell
    assert r["grad_dP0_autodiff"] > 0.0, "d(radius)/d(dP0) should be positive"
    assert r["grad_kedge_autodiff"] < 0.0, "d(radius)/d(k_edge) should be negative"
