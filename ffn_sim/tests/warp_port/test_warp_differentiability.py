"""B4 gate: Warp reverse-mode autodiff through a force law is clean + usable.

The Warp spike's prize (plan §B4): the native CUDA compartment plugin is not
differentiable. This asserts a clean reverse-mode gradient of the membrane
force-law loss w.r.t. the bare tension γ_mem, cross-checked against a closed-form
analytic gradient (machine precision) and a central finite-difference (numeric).
"""

from __future__ import annotations

import pytest


def test_warp_autodiff_matches_analytic():
    """Reverse-mode grad must match the closed-form analytic grad to ~machine eps."""
    pytest.importorskip("warp")
    from ffn_sim.warp_port.differentiability_b4 import run_check

    r = run_check(device="cpu")
    assert r["rel_err_vs_analytic"] < 1e-9, (
        f"autodiff vs analytic rel error too large: {r['rel_err_vs_analytic']:.3e}"
    )


def test_warp_autodiff_matches_finite_difference():
    """Reverse-mode grad must match a central finite-difference to FD accuracy."""
    pytest.importorskip("warp")
    from ffn_sim.warp_port.differentiability_b4 import run_check

    r = run_check(device="cpu")
    assert r["rel_err_vs_finite_diff"] < 1e-5, (
        f"autodiff vs finite-diff rel error too large: {r['rel_err_vs_finite_diff']:.3e}"
    )
    # the gradient must be non-trivial (a zero gradient would pass vacuously)
    assert abs(r["grad_autodiff"]) > 0.0, "autodiff gradient is exactly zero"
