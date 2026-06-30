"""DCM turgor parity gate: Warp exact-volume osmotic force vs committed HOOMD ref.

Reference fixture (``dcm/fixtures/dcm_turgor_ref.npz``) is the per-node net
force from the REAL ``cell.dcm.DcmTurgorForce`` run in a HOOMD CPU sim on a 2-cell
icosahedral mesh (see ``generate_dcm_turgor_fixture.py``). The Warp kernel is graded
against THAT committed HOOMD output (guard-rail 2).

Gate: ``reduce='host'`` (exact-volume reduced in numpy as the reference; force law
isolated) rel force < 1e-12. ``reduce='warp'`` (atomic-add volume) rel < 1e-8
(reduction-order diagnostic). Both achieve machine-eps here.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

FIX = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "dcm", "fixtures")


def _load() -> dict:
    path = os.path.join(FIX, "dcm_turgor_ref.npz")
    if not os.path.exists(path):
        pytest.skip(
            f"fixture {path} missing — run "
            "ffn_sim/dcm/fixtures/generate_dcm_turgor_fixture.py"
        )
    return dict(np.load(path))


def _run(fx: dict, reduce: str) -> dict:
    from ffn_sim.dcm.dcm_turgor_warp import run_dcm_turgor_warp

    return run_dcm_turgor_warp(
        pos=fx["pos"], faces=fx["faces"], face_cell=fx["face_cell"],
        n_cells=int(fx["n_cells"]), V0=float(fx["V0"]),
        turgor_dP0=float(fx["turgor_dP0"]), K_vol=float(fx["K_vol"]),
        reduce=reduce, device="cpu",
    )


def test_dcm_turgor_force_law_parity_host():
    """Force law (host-reduced exact volume): rel error < 1e-12."""
    pytest.importorskip("warp")
    fx = _load()
    got = _run(fx, "host")
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    dF = float(np.abs(got["force"] - fx["ref_force"]).max())
    assert dF / Fscale < 1e-12, f"turgor host force rel error too large: {dF/Fscale:.3e}"
    assert np.abs(fx["ref_force"]).max() > 0.0, "fixture has no turgor force"


def test_dcm_turgor_full_port_parity_warp_reduce():
    """Full Warp port (atomic-add volume reduction): rel error < 1e-8."""
    pytest.importorskip("warp")
    fx = _load()
    got = _run(fx, "warp")
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    dF = float(np.abs(got["force"] - fx["ref_force"]).max())
    assert dF / Fscale < 1e-8, f"turgor warp-reduce force rel error too large: {dF/Fscale:.3e}"
