"""Emit the committed Warp-vs-HOOMD parity-result JSON the results-gate checks.

The pytest gate (``tests/warp_port/test_*.py``) is the real verification — it
asserts each Warp piece matches the committed HOOMD reference fixture within the
plan's thresholds. THIS produces the lightweight, committed JSON summary that
``outputs/tag_kb/verify_runs.py`` (the results-integrity gate) reads on disk, so
a verified Warp piece is enforced by the same gate that guards every other
headline result (the C12 lesson: claim a committed lightweight artifact, never
an uncommitted bulk file).

Reproducible on CPU (Warp CPU backend): re-run to regenerate the JSON.

    python ffn_sim/warp_port/parity_report.py
"""

from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")
OUT = os.path.join(FIX, "warp_parity_results.json")


def _baoab_parity() -> dict:
    """Run the B1 BAOAB Warp kernel against both committed fixtures."""
    from ffn_sim.warp_port.baoab_warp import run_baoab_warp

    results = {}
    for label in ("kt0", "ktpos"):
        path = os.path.join(FIX, f"baoab_ref_{label}.npz")
        fx = dict(np.load(path))
        got = run_baoab_warp(
            pos0=fx["pos0"], image0=fx["image0"], force=fx["force"],
            gamma=fx["gamma"], kT=float(fx["kT"]), dt=float(fx["dt"]),
            noise=fx["noise"], box=tuple(float(x) for x in fx["box"]),
            device="cpu",
        )
        results[label] = {
            "kT": float(fx["kT"]),
            "N": int(fx["pos0"].shape[0]),
            "K": int(fx["K"]),
            "max_abs_pos_diff": float(np.abs(got["pos"] - fx["ref_pos"]).max()),
            "max_abs_image_diff": int(np.abs(got["image"] - fx["ref_image"]).max()),
            "max_abs_prv_diff": float(np.abs(got["prv"] - fx["ref_prv"]).max()),
            "max_abs_image_flag": int(np.abs(fx["ref_image"]).max()),  # wrap exercised
        }
    return results


def main() -> None:
    report = {
        "_about": (
            "Warp-CPU vs HOOMD-numpy bit-parity per ported piece. Reference = "
            "committed HOOMD fixtures (frozen Action / native plugin). Gate "
            "thresholds: kT=0 < 1e-9 (bit-for-bit), kT>0 < 1e-7 (drift)."
        ),
        "B1_baoab": _baoab_parity(),
    }
    with open(OUT, "w") as f:
        json.dump(report, f, indent=2)
    print(f"wrote {OUT}")
    print(json.dumps(report["B1_baoab"], indent=2))


if __name__ == "__main__":
    main()
