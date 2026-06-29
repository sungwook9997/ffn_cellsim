"""B1 bit-parity gate: Warp overdamped L-M BAOAB vs the committed HOOMD reference.

The reference fixtures (``warp_port/fixtures/baoab_ref_{kt0,ktpos}.npz``) are the
trajectory endpoint of the FROZEN ``integrator.baoab.LeimkuhlerMatthewsBAOAB``
Action run inside a real HOOMD CPU simulation (see ``generate_baoab_fixture.py``).
This test grades the Warp kernel against THAT committed output — never against a
fresh self-authored oracle (guard-rail 2 of the Warp spike plan).

Gate (per WARP_SPIKE_THEN_DECOHESION_PLAN_2026-06-19.md §B1):
  * kT=0  → deterministic, bit-for-bit  (max abs pos diff < 1e-9, image exact).
  * kT>0  → identical injected noise, drift < 1e-7 (float-op ordering only).

Parity runs on the Warp CPU backend (the Mac has no CUDA build); the gbook A5000
is only for the optional speed benchmark.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

FIX = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "warp_port", "fixtures")


def _load(label: str) -> dict:
    path = os.path.join(FIX, f"baoab_ref_{label}.npz")
    if not os.path.exists(path):
        pytest.skip(
            f"fixture {path} missing — run "
            "ffn_sim/warp_port/fixtures/generate_baoab_fixture.py"
        )
    return dict(np.load(path))


def _run_warp(fx: dict) -> dict:
    # imported here so a missing warp install skips rather than errors at collect
    from ffn_sim.dcm.baoab_warp import run_baoab_warp

    return run_baoab_warp(
        pos0=fx["pos0"],
        image0=fx["image0"],
        force=fx["force"],
        gamma=fx["gamma"],
        kT=float(fx["kT"]),
        dt=float(fx["dt"]),
        noise=fx["noise"],
        box=tuple(float(x) for x in fx["box"]),
        device="cpu",
    )


def test_baoab_warp_parity_kt0_bitforbit():
    """kT=0: pure (F/γ)·Δt drift — Warp must match HOOMD bit-for-bit (< 1e-9)."""
    pytest.importorskip("warp")
    fx = _load("kt0")
    assert float(fx["kT"]) == 0.0
    got = _run_warp(fx)

    dpos = np.abs(got["pos"] - fx["ref_pos"]).max()
    dimg = np.abs(got["image"] - fx["ref_image"]).max()
    assert dpos < 1e-9, f"kT=0 position parity failed: max|Δpos|={dpos:.3e}"
    assert dimg == 0, f"kT=0 image-flag parity failed: max|Δimage|={dimg}"


def test_baoab_warp_parity_ktpos_drift():
    """kT>0: identical injected noise — drift from float-op ordering < 1e-7."""
    pytest.importorskip("warp")
    fx = _load("ktpos")
    assert float(fx["kT"]) > 0.0
    got = _run_warp(fx)

    dpos = np.abs(got["pos"] - fx["ref_pos"]).max()
    dimg = np.abs(got["image"] - fx["ref_image"]).max()
    dprv = np.abs(got["prv"] - fx["ref_prv"]).max()
    assert dpos < 1e-7, f"kT>0 position drift too large: max|Δpos|={dpos:.3e}"
    # Image flags are integers; identical noise should give identical wrap counts
    # unless a particle sits within drift of an exact half-box face. Allow ≤1.
    assert dimg <= 1, f"kT>0 image-flag drift too large: max|Δimage|={dimg}"
    assert dprv < 1e-12, f"kT>0 prv (W_n store) parity failed: max|Δprv|={dprv:.3e}"


def test_fixture_genuinely_wraps():
    """Guard: the fixture must actually cross image boundaries (else the wrap
    path is untested). Both fixtures should show non-zero image flags."""
    for label in ("kt0", "ktpos"):
        fx = _load(label)
        assert np.abs(fx["ref_image"]).max() >= 1, (
            f"fixture {label} never wraps — wrap path not exercised"
        )
