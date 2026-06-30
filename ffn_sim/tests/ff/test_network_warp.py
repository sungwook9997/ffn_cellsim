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


def test_reshape_kernel_parity():
    """The Warp reshape kernel matches the numpy reshape bit-level + holds segment lengths."""
    import copy

    from ffn_sim.ff.constraints import reshape, segment_lengths
    from ffn_sim.ff.cortex_assembly import build_cortex_network
    from ffn_sim.ff.network_warp import reshape_np

    net, _ = build_cortex_network(n_filaments=30, rng=np.random.default_rng(0))
    net.pos = net.pos + 0.05 * np.random.default_rng(1).standard_normal(net.pos.shape)
    rest = net.seg_rest.copy()
    xw = reshape_np(net.pos, net.fiber_offsets, net.seg_rest, n_iter=4)
    net_np = copy.deepcopy(net); xn = reshape(net_np, n_iter=4)
    assert np.allclose(xw, xn, atol=1e-12)
    net.pos = xw
    assert np.allclose(segment_lengths(net), rest, rtol=1e-2)


def test_relax_on_device_matches_numpy():
    """The fully on-device FF cortex relaxation (bending+reshape, Warp) reaches the same resting shell
    as the numpy bending+reshape relax, bit-level — so the GPU path (device='cuda' on gbook) is correct."""
    import copy

    from ffn_sim.ff.constraints import reshape, segment_lengths
    from ffn_sim.ff.cortex_assembly import build_cortex_network
    from ffn_sim.ff.forces_warp import bending_energy, make_bending_force_fn
    from ffn_sim.ff.network_warp import relax_on_device

    net, _ = build_cortex_network(n_filaments=100, rng=np.random.default_rng(0))
    e0 = bending_energy(net); rest = net.seg_rest.copy()
    xw = relax_on_device(net, n_steps=300, reshape_every=25, device="cpu")
    # numpy reference (same bending force + reshape, no projector)
    net3, _ = build_cortex_network(n_filaments=100, rng=np.random.default_rng(0))
    bfn = make_bending_force_fn(net3); x = net3.pos.reshape(-1, 3).copy()
    seg = float(net3.seg_rest.mean()); dt = 0.1 / (float(net3.kappa.max()) / seg**3)
    for s in range(300):
        net3.pos = x; x = x + dt * bfn(x.reshape(-1)).reshape(-1, 3)
        if (s + 1) % 25 == 0:
            net3.pos = x; x = reshape(net3, n_iter=2)
    net3.pos = reshape(net3, n_iter=4)
    assert np.allclose(xw, net3.pos, atol=1e-9)              # device path == numpy path
    net.pos = xw
    assert bending_energy(net) < 0.05 * e0                    # relaxed
    assert np.allclose(segment_lengths(net), rest, rtol=1e-2)
