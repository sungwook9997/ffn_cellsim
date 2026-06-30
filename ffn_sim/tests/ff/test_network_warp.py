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


def test_branch_angle_kernel_parity_and_restoring():
    """The Arp2/3 angle-harmonic branch kernel matches the numpy reference bit-level AND restores the
    rest branch angle (a perturbed 120° branch relaxes to the 70° rest angle) — the mechanistic
    angle-harmonic branch (CLAUDE.md: NOT a rigid constraint)."""
    from ffn_sim.ff.network_warp import _branch_angle_force_numpy_ref, branch_angle_force_np

    rng = np.random.default_rng(0)
    pos = rng.standard_normal((12, 3))
    triples = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]], np.int64)
    th0 = np.deg2rad(70.0); k = 5.0
    fw = branch_angle_force_np(pos, triples, th0, k)
    fn = _branch_angle_force_numpy_ref(pos, triples, th0, k)
    assert np.allclose(fw, fn, atol=1e-12)

    # restoring: relax a single perturbed branch (120°) toward the 70° rest angle
    p = np.array([[1.0, 0, 0], [0, 0, 0],
                  [np.cos(np.deg2rad(120)), np.sin(np.deg2rad(120)), 0.0]])
    tri = np.array([[0, 1, 2]], np.int64)

    def angle(P):
        r1 = P[0] - P[1]; r2 = P[2] - P[1]
        return np.rad2deg(np.arccos(np.clip(np.dot(r1, r2) / (np.linalg.norm(r1) * np.linalg.norm(r2)), -1, 1)))

    x = p.copy()
    for _ in range(2000):
        f = _branch_angle_force_numpy_ref(x, tri, th0, k); f[1] = 0.0
        x = x + 0.01 * f
    assert abs(angle(x) - 70.0) < 1.0


def test_on_device_turgor_reduction_parity():
    """The on-device centroid + mean-radius reductions match the host turgor_pressure (so the GPU
    state-dependent turgor refresh is correct)."""
    import warp as wp

    from ffn_sim.ff.gamma_floor import CortexParams, build_crosslinked_cortex, turgor_pressure
    from ffn_sim.ff.network_warp import _csum, _rsum

    cx = build_crosslinked_cortex(CortexParams(), n_filaments=100, n_xl=100, n_myo=10,
                                  rng=np.random.default_rng(0))
    _, R_host = turgor_pressure(cx, cx.net.pos)
    N = cx.net.n_nodes
    pos_d = wp.array(np.ascontiguousarray(cx.net.pos), dtype=wp.vec3d, device="cpu")
    cacc = wp.zeros(1, dtype=wp.vec3d, device="cpu")
    wp.launch(_csum, dim=N, inputs=[pos_d, cacc], device="cpu")
    c = cacc.numpy()[0] / N
    centre = wp.vec3d(float(c[0]), float(c[1]), float(c[2]))
    racc = wp.zeros(1, dtype=wp.float64, device="cpu")
    wp.launch(_rsum, dim=N, inputs=[pos_d, centre, racc], device="cpu")
    R_dev = float(racc.numpy()[0]) / N
    assert abs(R_host - R_dev) < 1e-10


def test_loaded_shell_volume_stable_correct_cfl():
    """The on-device loaded shell (bending + crosslink + myosin + state-dependent turgor) is
    VOLUME-STABLE when the CFL includes the stiff turgor breathing-mode: V/V0 stays ≈1 (the tiny
    expansion self-relieves the resting osmotic excess to ΔP≈0) — NOT the CFL-artifact runaway that
    omitting the turgor stiffness produces."""
    from ffn_sim.ff.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
    from ffn_sim.ff.network_warp import simulate_loaded_shell_on_device

    cx = build_crosslinked_cortex(CortexParams(), n_filaments=120, n_xl=120, n_myo=12,
                                  rng=np.random.default_rng(0))
    _, traj = simulate_loaded_shell_on_device(cx, NMIIA_MINIFIL_STALL_PN, n_steps=8000,
                                              record_every=4000, device="cpu")
    vv0 = [t["V_over_V0"] for t in traj]
    assert all(0.95 < v < 1.05 for v in vv0)                  # stable, no blow-up
    assert abs(vv0[-1] - 1.0) < 0.01


def test_kmc_turnover_matches_bell_steady_state():
    """On-device Hand KMC turnover (Bell detach + Poisson re-attach) drives the engaged myosin
    fraction to the analytic steady state k_on/(k_on+p_off(f)) — and it self-limits under load
    (the engaged force-bearing density is what sets the γ-floor; turnover cannot raise it)."""
    from ffn_sim.ff.gamma_floor import (
        NMIIA_MINIFIL_STALL_PN,
        CortexParams,
        build_crosslinked_cortex,
        equilibrate,
    )
    from ffn_sim.ff.hand_kmc import NMIIA_MYOSIN, bell_off_rate
    from ffn_sim.ff.network_warp import simulate_turnover_on_device

    def steady(f):
        return NMIIA_MYOSIN.k_on / (NMIIA_MYOSIN.k_on + bell_off_rate(f, NMIIA_MYOSIN.p0, NMIIA_MYOSIN.f0))

    # physiological stall: stays ~fully engaged (k_on ≫ p_off) → floor unchanged
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=200, n_xl=200, n_myo=100,
                                  rng=np.random.default_rng(0))
    equilibrate(cx, 0.0, n_steps=200, method="device", device="cpu")
    _, bmask, _ = simulate_turnover_on_device(cx, NMIIA_MINIFIL_STALL_PN, n_steps=10000,
                                              kmc_every=50, tau_kmc=0.01, device="cpu")
    assert abs(bmask.mean() - steady(NMIIA_MINIFIL_STALL_PN)) < 0.05
    assert bmask.mean() > 0.9                                  # ~fully engaged at physiological load

    # high load: engaged fraction self-limits (Bell) well below 1 — turnover only LOWERS the density
    cx2 = build_crosslinked_cortex(CortexParams(), n_filaments=200, n_xl=200, n_myo=100,
                                   rng=np.random.default_rng(0))
    equilibrate(cx2, 0.0, n_steps=200, method="device", device="cpu")
    _, bmask2, _ = simulate_turnover_on_device(cx2, 30.0, n_steps=10000, kmc_every=50,
                                               tau_kmc=0.01, device="cpu")
    assert abs(bmask2.mean() - steady(30.0)) < 0.08
    assert bmask2.mean() < bmask.mean()                        # higher load ⇒ fewer engaged
