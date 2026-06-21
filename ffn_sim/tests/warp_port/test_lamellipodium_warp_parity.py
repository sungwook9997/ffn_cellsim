"""Lamellipodium parity gate: Warp traction-tether vs committed HOOMD reference.

Reference fixture (``warp_port/fixtures/lamellipodium_ref.npz``) is the per-node net
force from the REAL ``cell.dcm_lamellipodium.LamellipodialTractionTether`` run in a
HOOMD CPU sim on a 2-cell basal config (membrane rings + outward actin rings; see
``generate_lamellipodium_fixture.py``). The Warp port (host per-cell geometry +
leading-node selection, Warp per-node nearest-actin search + capped force) is graded
against THAT committed HOOMD output (guard-rail 2).
"""

from __future__ import annotations

import os

import numpy as np
import pytest

FIX = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "warp_port", "fixtures")


def _load() -> dict:
    path = os.path.join(FIX, "lamellipodium_ref.npz")
    if not os.path.exists(path):
        pytest.skip(
            f"fixture {path} missing — run "
            "ffn_sim/warp_port/fixtures/generate_lamellipodium_fixture.py"
        )
    return dict(np.load(path))


def test_lamellipodium_tether_warp_parity():
    """Warp lamellipodial traction-tether force matches HOOMD LamellipodialTractionTether."""
    pytest.importorskip("warp")
    from ffn_sim.warp_port.lamellipodium_warp import run_lamellipodium_tether_warp

    fx = _load()
    got = run_lamellipodium_tether_warp(
        pos=fx["pos"], typeid=fx["typeid"], cell_of_node=fx["cell_of_node"],
        rim_cells=fx["rim_cells"], actin_typeid=int(fx["actin_typeid"]),
        mem_typeid=int(fx["mem_typeid"]), z_basal=float(fx["z_basal"]),
        basal_band=float(fx["basal_band"]), k_tether=float(fx["k_tether"]),
        force_cap=float(fx["force_cap"]), tether_radius=float(fx["tether_radius"]),
        lead_frac=float(fx["lead_frac"]), device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    dF = float(np.abs(got["force"] - fx["ref_force"]).max())
    assert dF / Fscale < 1e-12, f"lamellipodium force rel error too large: {dF/Fscale:.3e}"
    assert np.abs(fx["ref_force"]).max() > 0.0, "fixture has no tether force"
