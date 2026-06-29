"""Optimization #6 gate: lamellipodium tether wired into the device-resident loop.

Asserts that the committed lamellipodial traction-tether kernel composes with
turgor + cortex edge springs in one GPU-resident per-step loop (device gather of
leading-node positions + accumulating tether), and produces outward spreading on a
single cell on a substrate (the footprint grows toward the clutch-anchored actin).
"""

from __future__ import annotations

import pytest


def test_lamellipodium_crawl_spreads():
    pytest.importorskip("warp")
    from ffn_sim.dcm.dcm_warp_crawl import run_crawl

    r = run_crawl(subdiv=2, steps=2000, device="cpu")
    assert r["finite"], "crawl loop went non-finite"
    assert r["n_leading"] > 0, "no leading nodes selected"
    # the tether traction must pull the basal rim outward -> footprint grows
    assert r["footprint_ratio"] > 1.05, (
        f"lamellipodium did not spread the cell: footprint_ratio={r['footprint_ratio']:.3f}"
    )
