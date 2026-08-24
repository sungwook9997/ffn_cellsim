"""Connected contractile network (ff/network_contractility) — Stage 6h decisive test.

Anchors: the patch builds a percolating network with sensible coordination z; dilution lowers z;
the contractile network develops a positive (contractile), finite σ that is deep UNDER the cortical
band (350-650 pN/µm) at physiological per-link force — the floor survives in the connected,
buckling-capable inextensible network. (Small/fast; the full f_act × z sweeps run separately.)
"""

import numpy as np
import pytest

from aleph.laws.network_contractility import build_patch, measure_sigma


def test_patch_builds_and_percolates():
    p = build_patch(nx=7, ny=7, keep_frac=1.0, rng=np.random.default_rng(0))
    assert p is not None
    assert 4.0 < p.z_mean <= 6.0                  # full triangular lattice interior z≈6
    assert p.pin.sum() >= 3                        # boundary pinned
    assert len(p.ep) > 10                          # many contractile fibers


def test_dilution_lowers_z():
    z_full = build_patch(nx=9, ny=9, keep_frac=1.0, rng=np.random.default_rng(1)).z_mean
    z_dil = build_patch(nx=9, ny=9, keep_frac=0.6, rng=np.random.default_rng(1)).z_mean
    assert z_dil < z_full                          # diluting edges reduces coordination


def test_contractile_sigma_positive_and_floored():
    """At physiological per-link force the connected contractile network develops a positive σ that
    is well under the 350 pN/µm band (the floor survives in the buckling-capable network)."""
    p = build_patch(nx=8, ny=8, keep_frac=1.0, rng=np.random.default_rng(2))
    out = measure_sigma(p, f_act=5.0, n_steps=2500)
    assert out is not None
    assert out["sigma_pN_um"] > 0.0               # contractile
    assert np.isfinite(out["sigma_pN_um"])
    assert out["sigma_pN_um"] < 350.0 / 5.0       # ≥5× under the band low edge (deeply floored)
