"""DCM substrate parity gate: Warp z-well + in-plane wetting vs committed HOOMD ref.

Reference fixture (``warp_port/fixtures/dcm_substrate_ref.npz``) holds the per-node net
forces from the REAL ``cell.dcm_gpu_forces.DcmSubstrateForceGPU`` (z-well + rigid-dish
floor) and ``DcmSubstrateWettingGPU`` (in-plane wetting), run on CPU on one icosahedral
cell straddling the substrate (see ``generate_dcm_substrate_fixture.py``). The Warp
kernels are graded against THOSE committed HOOMD outputs (guard-rail 2).

Gates:
* z-well (own-row write, no reduction): rel force < 1e-12 (machine-eps).
* wetting ``reduce='host'`` (per-face force law in Warp, numpy scatter): rel < 1e-12.
* wetting ``reduce='warp'`` (atomic scatter): rel < 1e-8 (reduction-order diagnostic).
"""

from __future__ import annotations

import os

import numpy as np
import pytest

FIX = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "warp_port", "fixtures")


def _load() -> dict:
    path = os.path.join(FIX, "dcm_substrate_ref.npz")
    if not os.path.exists(path):
        pytest.skip(
            f"fixture {path} missing — run "
            "ffn_sim/warp_port/fixtures/generate_dcm_substrate_fixture.py"
        )
    return dict(np.load(path))


def _rel(got: np.ndarray, ref: np.ndarray) -> float:
    return float(np.abs(got - ref).max()) / (float(np.abs(ref).max()) + 1e-30)


def test_dcm_substrate_well_parity():
    """z-well (plane-well + rigid floor), own-row write: rel error < 1e-12."""
    pytest.importorskip("warp")
    from ffn_sim.dcm.dcm_substrate_warp import run_dcm_substrate_well_warp

    fx = _load()
    got = run_dcm_substrate_well_warp(
        pos=fx["pos"], z0=float(fx["z0"]), W_cs=float(fx["W_cs"]),
        adh_range=float(fx["adh_range"]), k_sub=None, k_floor=float(fx["k_floor"]),
        device="cpu",
    )
    assert np.abs(fx["ref_force_well"]).max() > 0.0, "fixture has no z-well force"
    rel = _rel(got["force"], fx["ref_force_well"])
    assert rel < 1e-12, f"z-well force rel error too large: {rel:.3e}"


def test_dcm_substrate_wetting_parity_host():
    """Wetting per-face force law (host scatter): rel error < 1e-12."""
    pytest.importorskip("warp")
    from ffn_sim.dcm.dcm_substrate_warp import run_dcm_substrate_wetting_warp

    fx = _load()
    got = run_dcm_substrate_wetting_warp(
        pos=fx["pos"], faces=fx["faces"], z0=float(fx["z0"]),
        W_cs_Jm2=float(fx["W_cs_Jm2"]), adh_range=float(fx["adh_range"]),
        force_cap=float(fx["force_cap"]), reduce="host", device="cpu",
    )
    assert np.abs(fx["ref_force_wetting"]).max() > 0.0, "fixture has no wetting force"
    rel = _rel(got["force"], fx["ref_force_wetting"])
    assert rel < 1e-12, f"wetting host force rel error too large: {rel:.3e}"


def test_dcm_substrate_wetting_parity_warp_reduce():
    """Wetting full Warp port (atomic scatter): rel error < 1e-8."""
    pytest.importorskip("warp")
    from ffn_sim.dcm.dcm_substrate_warp import run_dcm_substrate_wetting_warp

    fx = _load()
    got = run_dcm_substrate_wetting_warp(
        pos=fx["pos"], faces=fx["faces"], z0=float(fx["z0"]),
        W_cs_Jm2=float(fx["W_cs_Jm2"]), adh_range=float(fx["adh_range"]),
        force_cap=float(fx["force_cap"]), reduce="warp", device="cpu",
    )
    rel = _rel(got["force"], fx["ref_force_wetting"])
    assert rel < 1e-8, f"wetting warp-reduce force rel error too large: {rel:.3e}"
