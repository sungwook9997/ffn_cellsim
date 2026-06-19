"""M-SHAKE parity gate: Warp Matrix-SHAKE vs committed HOOMD reference.

Reference fixture (``warp_port/fixtures/shake_ref.npz``) is the projected positions
+ accumulated Lagrange multipliers from the COMMITTED Python reference
``integrator.constrained_baoab.shake_project_chains`` (rigid, return_lambdas=True)
on the native-parity-test config (see ``generate_shake_fixture.py``). The Warp
kernel is graded against THAT committed output (guard-rail 2).

Why not strict bit-for-bit (unlike B1/B2): SHAKE is an ITERATIVE Newton solver, so
tiny float-op differences (Warp codegen may fuse multiply-add; numpy does not)
converge to a *slightly different but equally valid* point on the constraint
manifold. Both satisfy |bond|=ℓ₀ to ~tol. The gate is therefore the native
plugin's own contract: relative position/λ agreement plus constraint satisfaction
and zero nonconvergence. Achieved here: pos rel ~7e-12, λ rel ~4e-10 — both well
under the plan's 1e-9 and the native plugin's 1e-8.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

FIX = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "warp_port", "fixtures")


def _load() -> dict:
    path = os.path.join(FIX, "shake_ref.npz")
    if not os.path.exists(path):
        pytest.skip(
            f"fixture {path} missing — run "
            "ffn_sim/warp_port/fixtures/generate_shake_fixture.py"
        )
    return dict(np.load(path))


def _run(fx: dict) -> dict:
    from ffn_sim.warp_port.shake_warp import run_shake_warp

    return run_shake_warp(
        pred_pos=fx["pred"], ref_pos=fx["ref"], chains=fx["chains"],
        inv_mass=fx["inv_mass"], rest_length=float(fx["r0"]),
        box_L=(float(fx["L"]),) * 3, tol=float(fx["tol"]),
        max_iter=int(fx["max_iter"]), device="cpu",
    )


def test_shake_warp_parity():
    """Projected positions + λ match the reference; bonds satisfied; all converge."""
    pytest.importorskip("warp")
    fx = _load()
    got = _run(fx)

    assert got["nonconverged"] == 0, f"{got['nonconverged']} chains failed to converge"

    pscale = float(np.abs(fx["ref_proj"]).max())
    lscale = float(np.abs(fx["ref_lam"]).max()) + 1e-30
    dpos = float(np.abs(got["pos"] - fx["ref_proj"]).max())
    dlam = float(np.abs(got["lam"] - fx["ref_lam"]).max())
    assert dpos / pscale < 1e-9, f"position rel error too large: {dpos/pscale:.3e}"
    assert dlam / lscale < 1e-9, f"lambda rel error too large: {dlam/lscale:.3e}"

    # the Warp projection must itself satisfy the rigid-bond constraint to ~tol
    s = got["pos"][fx["chains"][:, :-1]] - got["pos"][fx["chains"][:, 1:]]
    bond = np.sqrt((s * s).sum(axis=2))
    bond_rel = float(np.abs(bond - float(fx["r0"])).max() / float(fx["r0"]))
    assert bond_rel < 1e-8, f"constraint not satisfied: max bond rel-err {bond_rel:.3e}"
