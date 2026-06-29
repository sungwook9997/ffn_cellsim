"""DCM node-face contact parity gate: Warp vs committed numpy GROUND TRUTH.

Reference fixture (``dcm/fixtures/dcm_contact_ref.npz``) is the per-node
contact force from the brute-force numpy ground truth
``cell.dcm_face_contact.node_face_contact_forces`` on two interpenetrating
icosahedral cells (see ``generate_dcm_contact_fixture.py``) — exercising both the
repulsion and bilinear-tent adhesion branches + the per-cell adhesion multiplier.
The Warp kernel is graded against THAT committed output (guard-rail 2).

NOTE: the node-face CONTACT FORCE has no dynamic topology (fixed mesh arrays), so it
ports straight to Warp — the no-dynamic-topology constraint bites the REMESH, a
separate later piece. Forces are atomic-accumulated; parity is tolerance-class.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

FIX = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "dcm", "fixtures")


def _load() -> dict:
    path = os.path.join(FIX, "dcm_contact_ref.npz")
    if not os.path.exists(path):
        pytest.skip(
            f"fixture {path} missing — run "
            "ffn_sim/dcm/fixtures/generate_dcm_contact_fixture.py"
        )
    return dict(np.load(path))


def test_dcm_node_face_contact_warp_parity():
    """Warp node-face contact force matches the numpy brute-force ground truth."""
    pytest.importorskip("warp")
    from ffn_sim.dcm.dcm_contact_warp import run_node_face_contact_warp

    fx = _load()
    got = run_node_face_contact_warp(
        pos=fx["pos"], cell_of_node=fx["cell_of_node"], faces=fx["faces"],
        face_cell=fx["face_cell"], rep_strength=float(fx["rep_strength"]),
        adh_strength=float(fx["adh_strength"]), c_rep=float(fx["c_rep"]),
        c_adh=float(fx["c_adh"]), cad_mult=fx["cad_mult"], device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    dF = float(np.abs(got["force"] - fx["ref_force"]).max())
    assert dF / Fscale < 1e-12, f"contact force rel error too large: {dF/Fscale:.3e}"
    # the fixture must actually exercise contact (else parity is vacuous)
    assert np.abs(fx["ref_force"]).max() > 0.0, "fixture has no contact force"
