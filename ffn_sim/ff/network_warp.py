"""GPU-native FF network force kernels (Warp) — the cortex/network forces, device-agnostic.

The whole point of Warp (vs Cytosim's C++ CPU core) is GPU-native + differentiable execution: the
same kernel runs on CPU (Mac dev) or CUDA (gbook A5000) by passing ``device``. The FF bending force
(``forces_warp``) and the implicit CG (``dcm.dcm_warp_implicit.device_cg``) are already Warp; this
module adds the remaining FF network forces — crosslinker/actin axial springs, myosin contractile
links, and the radial osmotic turgor — as Warp kernels (one thread per element, ``wp.atomic_add``,
float64), matching the DCM ``network_warp`` pattern. Together they let the FF cortex relaxation run
fully on-device (no numpy round-trip), so gbook's GPU is used via ``device="cuda"`` — NOT a separate
"GPU port".

Validation: each kernel is bit-parity-checked against the original numpy force in
``gamma_floor`` / ``network_contractility`` (`tests/ff/test_network_warp.py`).
"""

from __future__ import annotations

import numpy as np
import warp as wp

wp.init()


@wp.kernel
def link_spring_kernel(
    pos: wp.array(dtype=wp.vec3d),              # (N,) node positions
    links: wp.array(dtype=wp.int32, ndim=2),    # (L, 2) node-index pairs (i, j)
    k_arr: wp.array(dtype=wp.float64),          # (L,) stiffness [pN/µm]
    r0_arr: wp.array(dtype=wp.float64),         # (L,) rest length [µm]
    force: wp.array(dtype=wp.vec3d),            # (N,) out (atomic accumulate)
):
    """Hookean link force k·(L−r0)·û on i (toward j when stretched), −on j. Crosslinker + actin axial."""
    t = wp.tid()
    i = links[t, 0]
    j = links[t, 1]
    d = pos[j] - pos[i]
    L = wp.length(d)
    if L > wp.float64(1e-12):
        f = (k_arr[t] * (L - r0_arr[t]) / L) * d
        wp.atomic_add(force, i, f)
        wp.atomic_add(force, j, -f)


@wp.kernel
def myosin_kernel(
    pos: wp.array(dtype=wp.vec3d),
    links: wp.array(dtype=wp.int32, ndim=2),    # (M, 2) myosin link node pairs
    f_myo: wp.float64,                          # constant contractile force [pN]
    force: wp.array(dtype=wp.vec3d),
):
    """Myosin contractile force: constant magnitude f_myo pulling i,j together (toward each other)."""
    t = wp.tid()
    i = links[t, 0]
    j = links[t, 1]
    d = pos[j] - pos[i]
    L = wp.length(d)
    if L > wp.float64(1e-12):
        f = (f_myo / L) * d
        wp.atomic_add(force, i, f)              # on i toward j (shorten)
        wp.atomic_add(force, j, -f)


@wp.kernel
def turgor_kernel(
    pos: wp.array(dtype=wp.vec3d),
    centre: wp.vec3d,
    dP_area: wp.float64,                         # ΔP · (area/node) [pN]
    force: wp.array(dtype=wp.vec3d),
):
    """Radial outward osmotic turgor force ΔP·(area/node)·r̂ per node (Young-Laplace shell)."""
    i = wp.tid()
    d = pos[i] - centre
    r = wp.length(d)
    if r > wp.float64(1e-12):
        wp.atomic_add(force, i, (dP_area / r) * d)


@wp.kernel
def axpy_kernel(x: wp.array(dtype=wp.vec3d), step: wp.float64, F: wp.array(dtype=wp.vec3d)):
    """Explicit overdamped step x += step·F (step = dt/γ)."""
    i = wp.tid()
    x[i] = x[i] + step * F[i]


@wp.kernel
def reshape_kernel(
    pos: wp.array(dtype=wp.vec3d),
    fiber_off: wp.array(dtype=wp.int32),        # (F+1,) node range [off[f], off[f+1]) per fiber
    seg_off: wp.array(dtype=wp.int32),          # (F+1,) segment range start per fiber
    seg_rest: wp.array(dtype=wp.float64),       # (S,) segment rest lengths
    n_iter: wp.int32,
):
    """NF2007 §5.3 reshape (one thread per fiber): restore |m_{k+1}−m_k|=seg_rest, conserving COG.

    Fibers own disjoint node ranges → no inter-thread races. Sequential over segments per the paper;
    ``n_iter`` passes converge the residual coupling."""
    f = wp.tid()
    a = fiber_off[f]
    b = fiber_off[f + 1]
    p = b - a - 1                                # segments
    s0 = seg_off[f]
    for _it in range(n_iter):
        for k in range(p):
            g = pos[a + k + 1] - pos[a + k]
            d = wp.length(g)
            if d > wp.float64(1e-300):
                u = g / d
                e = d - seg_rest[s0 + k]
                n_a = wp.float64(k + 1)
                n_b = wp.float64(p - k)
                tot = n_a + n_b
                d_a = e * n_b / tot
                d_b = e * n_a / tot
                for m in range(a, a + k + 1):
                    pos[m] = pos[m] + d_a * u
                for m in range(a + k + 1, b):
                    pos[m] = pos[m] - d_b * u


def reshape_np(pos, fiber_off, seg_rest, n_iter=4, device="cpu"):
    """Convenience: Warp reshape on a numpy (N,3) pos → reshaped numpy (for parity / numpy callers)."""
    fiber_off = np.ascontiguousarray(fiber_off, np.int32)
    seg_per = np.diff(fiber_off) - 1
    seg_off = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    pos_d = wp.array(np.ascontiguousarray(pos, np.float64), dtype=wp.vec3d, device=device)
    wp.launch(reshape_kernel, dim=fiber_off.shape[0] - 1,
              inputs=[pos_d, wp.array(fiber_off, dtype=wp.int32, device=device),
                      wp.array(seg_off, dtype=wp.int32, device=device),
                      wp.array(np.ascontiguousarray(seg_rest, np.float64), dtype=wp.float64, device=device),
                      wp.int32(n_iter)], device=device)
    wp.synchronize_device(device)
    return pos_d.numpy().astype(np.float64)


def link_spring_force_np(pos, links, k, r0, device="cpu"):
    """Convenience: Warp link-spring force as a numpy (N,3) array (for parity tests / numpy callers)."""
    N = pos.shape[0]
    pos_d = wp.array(np.ascontiguousarray(pos, dtype=np.float64), dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    if links.shape[0]:
        wp.launch(link_spring_kernel, dim=links.shape[0],
                  inputs=[pos_d, wp.array(np.ascontiguousarray(links, np.int32), dtype=wp.int32, device=device),
                          wp.array(np.ascontiguousarray(k, np.float64), dtype=wp.float64, device=device),
                          wp.array(np.ascontiguousarray(r0, np.float64), dtype=wp.float64, device=device),
                          force_d], device=device)
        wp.synchronize_device(device)
    return force_d.numpy().astype(np.float64)


def myosin_force_np(pos, links, f_myo, device="cpu"):
    """Convenience: Warp myosin contractile force as a numpy (N,3) array."""
    N = pos.shape[0]
    pos_d = wp.array(np.ascontiguousarray(pos, dtype=np.float64), dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    if links.shape[0]:
        wp.launch(myosin_kernel, dim=links.shape[0],
                  inputs=[pos_d, wp.array(np.ascontiguousarray(links, np.int32), dtype=wp.int32, device=device),
                          wp.float64(f_myo), force_d], device=device)
        wp.synchronize_device(device)
    return force_d.numpy().astype(np.float64)


@wp.kernel
def _zero(f: wp.array(dtype=wp.vec3d)):
    f[wp.tid()] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))


@wp.kernel
def _csum(pos: wp.array(dtype=wp.vec3d), acc: wp.array(dtype=wp.vec3d)):
    """Sum all node positions into acc[0] (→ centroid after /N). On-device reduction."""
    wp.atomic_add(acc, 0, pos[wp.tid()])


@wp.kernel
def _rsum(pos: wp.array(dtype=wp.vec3d), centre: wp.vec3d, acc: wp.array(dtype=wp.float64)):
    """Sum |pos − centre| into acc[0] (→ mean shell radius after /N). On-device reduction."""
    wp.atomic_add(acc, 0, wp.length(pos[wp.tid()] - centre))


def relax_on_device(net, *, links=None, k_xl=None, xl_rest=None, myo_links=None, f_myo=0.0,
                    n_steps=600, reshape_every=25, dt_mu=0.0, n_reshape_iter=2, device="cpu"):
    """Fully on-device FF cortex relaxation: bending (+ optional crosslink/myosin) forces + reshape,
    all Warp — NO per-step numpy round-trip. Runs on CPU (Mac) or CUDA (gbook A5000) via ``device``.

    Returns the relaxed positions (N,3) numpy. Inextensibility via periodic reshape (the robust path;
    no singular-prone projector). ``links``/``k_xl``/``xl_rest`` optional crosslinker springs;
    ``myo_links``/``f_myo`` optional myosin. Default (none) = resting-shell bending settle.
    """
    from ffn_sim.ff.forces_warp import _per_triple_alpha
    N = net.n_nodes
    tri = np.ascontiguousarray(net.bend_triples, np.int32)
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    fiber_off = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg_per = np.diff(fiber_off) - 1
    seg_off = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    seg = float(net.seg_rest.mean()) if net.seg_rest.size else 1.0
    if dt_mu <= 0.0:
        kmax = float(net.kappa.max()) / seg**3
        if k_xl is not None and len(k_xl):
            kmax = max(kmax, float(np.max(k_xl)))
        dt_mu = 0.1 / kmax

    from ffn_sim.ff.forces_warp import cytosim_bending_kernel
    d = device
    pos_d = wp.array(np.ascontiguousarray(net.pos, np.float64), dtype=wp.vec3d, device=d)
    f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    tri_d = wp.array(tri, dtype=wp.int32, device=d)
    alpha_d = wp.array(alpha, dtype=wp.float64, device=d)
    foff_d = wp.array(fiber_off, dtype=wp.int32, device=d)
    soff_d = wp.array(seg_off, dtype=wp.int32, device=d)
    srest_d = wp.array(np.ascontiguousarray(net.seg_rest, np.float64), dtype=wp.float64, device=d)
    has_xl = links is not None and len(links)
    if has_xl:
        xl_d = wp.array(np.ascontiguousarray(links, np.int32), dtype=wp.int32, device=d)
        kxl_d = wp.array(np.ascontiguousarray(k_xl, np.float64), dtype=wp.float64, device=d)
        r0_d = wp.array(np.ascontiguousarray(xl_rest, np.float64), dtype=wp.float64, device=d)
    has_myo = myo_links is not None and len(myo_links) and f_myo != 0.0
    if has_myo:
        myo_d = wp.array(np.ascontiguousarray(myo_links, np.int32), dtype=wp.int32, device=d)
    nT = tri.shape[0]
    for step in range(n_steps):
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        if nT:
            wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        if has_xl:
            wp.launch(link_spring_kernel, dim=len(links), inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            wp.launch(myosin_kernel, dim=len(myo_links), inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(axpy_kernel, dim=N, inputs=[pos_d, wp.float64(dt_mu), f_d], device=d)
        if (step + 1) % reshape_every == 0:
            wp.launch(reshape_kernel, dim=fiber_off.shape[0] - 1,
                      inputs=[pos_d, foff_d, soff_d, srest_d, wp.int32(n_reshape_iter)], device=d)
    wp.launch(reshape_kernel, dim=fiber_off.shape[0] - 1,
              inputs=[pos_d, foff_d, soff_d, srest_d, wp.int32(4)], device=d)
    wp.synchronize_device(d)
    return pos_d.numpy().astype(np.float64)


def simulate_loaded_shell_on_device(cortex, f_myo, *, n_steps=4000, reshape_every=25,
                                    turgor_every=20, record_every=100, dt_mu=0.0,
                                    n_reshape_iter=2, device="cpu"):
    """Evolve the ACTIVELY-LOADED cortex shell fully on-device (Warp): bending + crosslink springs +
    myosin contraction + STATE-DEPENDENT osmotic turgor + reshape, on the GPU (gbook A5000 via
    ``device="cuda:0"``). The position array stays device-resident; only ~4 scalars (centroid +
    mean-radius reductions) cross the bus every ``turgor_every`` steps to refresh ΔP via the Guo-2017
    closure (``gamma_floor.turgor_pressure``).

    This is the dynamic complement to ``relax_on_device``: where that settles the turgor-free resting
    shell, this runs the loaded shell — whose trajectory (V/V0, ΔP, R_mean) is the γ-floor signature
    (the floored actomyosin network has no static equilibrium against physiological turgor; it
    inflates/collapses rather than holding a tension at band). Returns (pos (N,3), traj) where traj is
    a list of dicts {step, R_mean, dP, V_over_V0}.
    """
    from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
    from ffn_sim.ff.gamma_floor import (
        TURGOR_DP0, TURGOR_PI_IN0, VMIN_FRAC,
    )
    N = cortex.net.n_nodes
    net = cortex.net
    tri = np.ascontiguousarray(net.bend_triples, np.int32)
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    fiber_off = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg_per = np.diff(fiber_off) - 1
    seg_off = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    seg = float(net.seg_rest.mean()) if net.seg_rest.size else 1.0
    # Guo closure params (host, computed once)
    R0 = cortex.R0_mean if cortex.R0_mean > 0.0 else cortex.R_um
    V0 = (4.0 / 3.0) * np.pi * R0**3
    vmin = VMIN_FRAC * V0
    if dt_mu <= 0.0:
        kmax = float(net.kappa.max()) / seg**3
        if cortex.xl_i.size:
            kmax = max(kmax, float(cortex.xl_k.max()))
        # CFL must include the STIFF turgor breathing-mode (else the shell numerically blows up):
        # K_vol = Π_in0/(1−vmin_frac); the collective radial stiffness per node for a uniform δr is
        # k_turgor = (K_vol/V0)·(4πR0²)²/N (δV=4πR²δr, δP=−(K_vol/V)δV, force=ΔP·area). This dominates.
        K_vol = TURGOR_PI_IN0 / (1.0 - VMIN_FRAC)
        k_turgor = (K_vol / V0) * (4.0 * np.pi * R0**2) ** 2 / N
        kmax = max(kmax, k_turgor)
        dt_mu = 0.1 / kmax

    d = device
    pos_d = wp.array(np.ascontiguousarray(net.pos, np.float64), dtype=wp.vec3d, device=d)
    f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    cacc = wp.zeros(1, dtype=wp.vec3d, device=d)
    racc = wp.zeros(1, dtype=wp.float64, device=d)
    tri_d = wp.array(tri, dtype=wp.int32, device=d)
    alpha_d = wp.array(alpha, dtype=wp.float64, device=d)
    foff_d = wp.array(fiber_off, dtype=wp.int32, device=d)
    soff_d = wp.array(seg_off, dtype=wp.int32, device=d)
    srest_d = wp.array(np.ascontiguousarray(net.seg_rest, np.float64), dtype=wp.float64, device=d)
    has_xl = bool(cortex.xl_i.size)
    if has_xl:
        xl_d = wp.array(np.ascontiguousarray(np.stack([cortex.xl_i, cortex.xl_j], 1), np.int32),
                        dtype=wp.int32, device=d)
        kxl_d = wp.array(np.ascontiguousarray(cortex.xl_k, np.float64), dtype=wp.float64, device=d)
        r0_d = wp.array(np.ascontiguousarray(cortex.xl_rest, np.float64), dtype=wp.float64, device=d)
    has_myo = bool(cortex.myo_i.size) and f_myo != 0.0
    if has_myo:
        myo_d = wp.array(np.ascontiguousarray(np.stack([cortex.myo_i, cortex.myo_j], 1), np.int32),
                         dtype=wp.int32, device=d)
    nT = tri.shape[0]
    centre = wp.vec3d(0.0, 0.0, 0.0)
    dP_area = 0.0
    traj = []

    def _refresh_turgor():
        nonlocal centre, dP_area
        cacc.zero_(); wp.launch(_csum, dim=N, inputs=[pos_d, cacc], device=d)
        c = cacc.numpy()[0] / N
        centre = wp.vec3d(float(c[0]), float(c[1]), float(c[2]))
        racc.zero_(); wp.launch(_rsum, dim=N, inputs=[pos_d, centre, racc], device=d)
        R_mean = float(racc.numpy()[0]) / N
        V = (4.0 / 3.0) * np.pi * R_mean**3
        dP = max(TURGOR_PI_IN0 * (V0 - vmin) / max(V - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0), 0.0)
        dP_area = dP * 4.0 * np.pi * R_mean**2 / N
        return R_mean, dP, V / V0

    for step in range(n_steps):
        if step % turgor_every == 0:
            R_mean, dP, vv0 = _refresh_turgor()
            if step % record_every == 0:
                traj.append({"step": step, "R_mean": R_mean, "dP": dP, "V_over_V0": vv0})
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        if nT:
            wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        if has_xl:
            wp.launch(link_spring_kernel, dim=cortex.xl_i.size, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            wp.launch(myosin_kernel, dim=cortex.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(turgor_kernel, dim=N, inputs=[pos_d, centre, wp.float64(dP_area), f_d], device=d)
        wp.launch(axpy_kernel, dim=N, inputs=[pos_d, wp.float64(dt_mu), f_d], device=d)
        if (step + 1) % reshape_every == 0:
            wp.launch(reshape_kernel, dim=fiber_off.shape[0] - 1,
                      inputs=[pos_d, foff_d, soff_d, srest_d, wp.int32(n_reshape_iter)], device=d)
    R_mean, dP, vv0 = _refresh_turgor()
    traj.append({"step": n_steps, "R_mean": R_mean, "dP": dP, "V_over_V0": vv0})
    wp.synchronize_device(d)
    net.pos = pos_d.numpy().astype(np.float64)
    return net.pos, traj
