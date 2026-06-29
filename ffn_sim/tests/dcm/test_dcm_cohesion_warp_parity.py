"""DCM cohesion parity gate: Warp node-node cohesion vs committed HOOMD reference.

Reference fixture (``dcm/fixtures/dcm_cohesion_ref.npz``) is the per-node net
force from the REAL ``cell.dcm_contact.DcmTentContact`` (inter-cell repulsion +
adhesive tent) run in a HOOMD CPU sim on two overlapping node clouds spanning both
branches (see ``generate_dcm_cohesion_fixture.py``). The Warp kernel is graded
against THAT committed HOOMD output (guard-rail 2).

Atomic-free Newton-3 (each node computes its own pair forces; the symmetric pair is
recomputed by the other node with exact-opposite sign), so the only parity gap is
per-node summation order vs the reference's np.add.at — gated rel < 1e-12.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

FIX = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "dcm", "fixtures")


def _load() -> dict:
    path = os.path.join(FIX, "dcm_cohesion_ref.npz")
    if not os.path.exists(path):
        pytest.skip(
            f"fixture {path} missing — run "
            "ffn_sim/dcm/fixtures/generate_dcm_cohesion_fixture.py"
        )
    return dict(np.load(path))


def test_dcm_cohesion_warp_parity():
    """Warp node-node cohesion force matches HOOMD DcmTentContact."""
    pytest.importorskip("warp")
    from ffn_sim.dcm.dcm_cohesion_warp import run_dcm_cohesion_warp

    fx = _load()
    got = run_dcm_cohesion_warp(
        pos=fx["pos"], cell_of_node=fx["cell_of_node"],
        r_contact=float(fx["r_contact"]), c_adh=float(fx["c_adh"]),
        rep_strength=float(fx["rep"]), adh_strength=float(fx["adh"]),
        patch_area=float(fx["patch_area"]), force_cap=float(fx["force_cap"]),
        cad_mult=fx["cad_mult"], device="cpu",
    )
    Fscale = float(np.abs(fx["ref_force"]).max()) + 1e-30
    dF = float(np.abs(got["force"] - fx["ref_force"]).max())
    assert dF / Fscale < 1e-12, f"cohesion force rel error too large: {dF/Fscale:.3e}"
    assert np.abs(fx["ref_force"]).max() > 0.0, "fixture has no cohesion force"
