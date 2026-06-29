"""B2 parity gate: Warp radial-shell compartment force vs committed HOOMD reference.

Reference fixtures (``warp_port/fixtures/radial_ref_{nucleus,membrane,turgor}.npz``)
are per-bead force + energy from the three REAL production compartment forces
(``NucleusConfinement`` / ``MembraneSurfaceTension`` / ``EnclosedVolumePressure``)
run in a HOOMD CPU sim (see ``generate_radial_shell_fixture.py``). The Warp kernel
is graded against THAT committed numpy output (guard-rail 2), never a fresh oracle.

Gate (per WARP_SPIKE_THEN_DECOHESION_PLAN_2026-06-19.md §B2, "1e-12..1e-15-class"):
  * reduce="host" (force LAW isolated; centroid/R_mean computed identically to the
    reference): rel force + energy error < 1e-12 — the committed bit-parity claim.
  * reduce="warp" (full port; centroid/R_mean reduced in Warp via atomic_add): rel
    error < 1e-8 — the native plugin's own reduction-order tolerance, since atomic
    summation order differs from numpy's pairwise sum ("modulo atomic-reduction
    order", as the native parity script states). Reported as a diagnostic.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

FIX = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "warp_port", "fixtures")

LAWS = ("nucleus", "membrane", "turgor")


def _load(name: str) -> dict:
    path = os.path.join(FIX, f"radial_ref_{name}.npz")
    if not os.path.exists(path):
        pytest.skip(
            f"fixture {path} missing — run "
            "ffn_sim/warp_port/fixtures/generate_radial_shell_fixture.py"
        )
    return dict(np.load(path))


def _run(fx: dict, reduce: str) -> dict:
    from ffn_sim.dcm.radial_shell_warp import run_radial_shell_warp

    return run_radial_shell_warp(
        pos=fx["pos"], tag=fx["tag"],
        tag_range=(int(fx["t0"]), int(fx["t1"])),
        law=int(fx["law"]), R0=float(fx["R0"]), pa=float(fx["pa"]),
        pb=float(fx["pb"]), pc=float(fx["pc"]), pd=float(fx["pd"]),
        reduce=reduce, device="cpu",
    )


def _rel(got: dict, fx: dict) -> tuple[float, float]:
    dF = float(np.abs(got["force"] - fx["ref_force"]).max())
    sF = float(np.abs(fx["ref_force"]).max()) + 1e-300
    dU = float(np.abs(got["energy"] - fx["ref_energy"]).max())
    sU = float(np.abs(fx["ref_energy"]).max()) + 1e-300
    return dF / sF, dU / sU


@pytest.mark.parametrize("name", LAWS)
def test_radial_shell_force_law_parity_host(name):
    """Force LAW (host-reduced centroid/R_mean): rel error < 1e-12."""
    pytest.importorskip("warp")
    fx = _load(name)
    rF, rU = _rel(_run(fx, "host"), fx)
    assert rF < 1e-12, f"{name} host force rel error too large: {rF:.3e}"
    assert rU < 1e-12, f"{name} host energy rel error too large: {rU:.3e}"


@pytest.mark.parametrize("name", LAWS)
def test_radial_shell_full_port_parity_warp_reduce(name):
    """Full Warp port (atomic-add reduction): rel error < 1e-8 (reduction-order)."""
    pytest.importorskip("warp")
    fx = _load(name)
    rF, rU = _rel(_run(fx, "warp"), fx)
    assert rF < 1e-8, f"{name} warp-reduce force rel error too large: {rF:.3e}"
    assert rU < 1e-8, f"{name} warp-reduce energy rel error too large: {rU:.3e}"


def test_fixtures_nonzero_force():
    """Guard: every fixture must carry a non-trivial force (else parity is vacuous)."""
    for name in LAWS:
        fx = _load(name)
        assert np.abs(fx["ref_force"]).max() > 0.0, f"{name} fixture force is all-zero"
