"""Fixman parity gate: Warp metric pseudo-force vs committed HOOMD reference.

Reference fixture (``warp_port/fixtures/fixman_ref.npz``) is the force + U_F from
the COMMITTED Python ``integrator.constrained_baoab.fixman_logdet_and_force`` on
the native-parity-test config (see ``generate_fixman_fixture.py``). The Warp
kernel is graded against THAT committed output (guard-rail 2).

Parity is tolerance-based by construction: the Python reference uses LAPACK
``slogdet``/``inv`` while the Warp kernel uses the tridiagonal continuant
recurrence — they differ in the last ULPs. Gate = the native plugin's own
contract (force rel < 1e-7, U_F rel < 1e-9, no non-positive metric). Achieved
here: force rel ~6e-16, U_F rel ~0 — far inside the bounds.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

FIX = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "warp_port", "fixtures")


def _load() -> dict:
    path = os.path.join(FIX, "fixman_ref.npz")
    if not os.path.exists(path):
        pytest.skip(
            f"fixture {path} missing — run "
            "ffn_sim/warp_port/fixtures/generate_fixman_fixture.py"
        )
    return dict(np.load(path))


def test_fixman_warp_parity():
    """Fixman force + U_F match the LAPACK reference; metric stays positive."""
    pytest.importorskip("warp")
    from ffn_sim.warp_port.fixman_warp import run_fixman_warp

    fx = _load()
    got = run_fixman_warp(
        pos=fx["pos"], chains=fx["chains"], inv_gamma=fx["inv_gamma"],
        kT=float(fx["kT"]), box_L=(float(fx["L"]),) * 3, device="cpu",
    )
    assert got["bad_sign"] == 0, f"{got['bad_sign']} chains had non-positive metric"

    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    Uscale = abs(float(fx["ref_U"])) + 1e-30
    dF = float(np.abs(got["force"] - fx["ref_force"]).max())
    dU = abs(got["U_F"] - float(fx["ref_U"]))
    assert dF / Fscale < 1e-7, f"Fixman force rel error too large: {dF/Fscale:.3e}"
    assert dU / Uscale < 1e-9, f"Fixman U_F rel error too large: {dU/Uscale:.3e}"
