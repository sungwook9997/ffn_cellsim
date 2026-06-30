"""GPU-native FF network force kernels (ff/network_warp) — parity vs the numpy reference.

The FF network forces must run on-device (Warp, so device='cuda' on gbook A5000 works without a
separate "GPU port"). These assert the Warp kernels reproduce the original numpy forces
(gamma_floor._link_spring_force / _myosin_force) to float64 precision — so the GPU path is correct.
"""

import numpy as np
import pytest

from ffn_sim.ff.gamma_floor import _link_spring_force, _myosin_force
from ffn_sim.ff.network_warp import link_spring_force_np, myosin_force_np


def _rand_net(n=40, nlink=60, seed=0):
    rng = np.random.default_rng(seed)
    pos = rng.standard_normal((n, 3))
    i = rng.integers(0, n, nlink)
    j = (i + rng.integers(1, n, nlink)) % n
    links = np.stack([i, j], axis=1).astype(np.int64)
    return pos, links


def test_link_spring_parity():
    pos, links = _rand_net()
    k = np.full(links.shape[0], 0.3)
    r0 = np.full(links.shape[0], 0.5)
    ref = np.zeros_like(pos)
    _link_spring_force(pos, links[:, 0], links[:, 1], k, r0, ref)
    warp = link_spring_force_np(pos, links, k, r0)
    assert np.allclose(warp, ref, atol=1e-12, rtol=1e-9)


def test_myosin_parity():
    pos, links = _rand_net(seed=2)
    ref = np.zeros_like(pos)
    _myosin_force(pos, links[:, 0], links[:, 1], 5.0, ref)
    warp = myosin_force_np(pos, links, 5.0)
    assert np.allclose(warp, ref, atol=1e-12, rtol=1e-9)


def test_zero_links_safe():
    pos = np.random.default_rng(0).standard_normal((10, 3))
    empty = np.zeros((0, 2), np.int64)
    assert np.allclose(link_spring_force_np(pos, empty, np.zeros(0), np.zeros(0)), 0.0)
    assert np.allclose(myosin_force_np(pos, empty, 5.0), 0.0)
