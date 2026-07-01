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
def myosin_bound_kernel(
    pos: wp.array(dtype=wp.vec3d),
    links: wp.array(dtype=wp.int32, ndim=2),
    bound: wp.array(dtype=wp.int32),            # (M,) 1=engaged, 0=detached
    f_myo: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """Myosin contractile force, but only ENGAGED (bound==1) links pull — the rest are detached
    (Hand turnover). The engaged fraction is what the γ-floor identifies as the limiting density."""
    t = wp.tid()
    if bound[t] == 0:
        return
    i = links[t, 0]
    j = links[t, 1]
    d = pos[j] - pos[i]
    L = wp.length(d)
    if L > wp.float64(1e-12):
        f = (f_myo / L) * d
        wp.atomic_add(force, i, f)
        wp.atomic_add(force, j, -f)


@wp.kernel
def kmc_bell_turnover_kernel(
    pos: wp.array(dtype=wp.vec3d),
    links: wp.array(dtype=wp.int32, ndim=2),
    bound: wp.array(dtype=wp.int32),            # (M,) in/out engaged state
    f_load: wp.float64,                         # contractile load per engaged hand [pN] (=f_myo)
    p0: wp.float64,                             # Bell zero-force off-rate [1/s]
    f0: wp.float64,                             # Bell characteristic force [pN]
    cap_um: wp.float64,                         # capture radius for re-attachment [µm]
    k_on: wp.float64,                           # on-rate [1/s]
    tau: wp.float64,                            # KMC tick [s]
    seed: wp.int32,                             # per-call RNG seed (step-dependent)
):
    """NF2007 §10.1 Hand turnover (one thread per link): a bound hand DETACHES with the Bell
    force-dependent probability 1−exp(−τ·p₀·exp(|f|/f₀)); a detached hand within ``cap_um`` RE-ATTACHES
    with 1−exp(−τ·k_on). The load on an engaged myosin is its contractile force f_load (Bell sees the
    motor's own force). This makes the ENGAGED fraction self-limiting under load (the γ-floor mechanism)."""
    t = wp.tid()
    rstate = wp.rand_init(seed, t)
    i = links[t, 0]
    j = links[t, 1]
    L = wp.length(pos[j] - pos[i])
    if bound[t] == 1:
        p_off = p0 * wp.exp(wp.abs(f_load) / f0)        # Bell slip on the motor's own load
        p_det = wp.float64(1.0) - wp.exp(-tau * p_off)
        if wp.float64(wp.randf(rstate)) < p_det:
            bound[t] = wp.int32(0)
    else:
        if L < cap_um:
            p_att = wp.float64(1.0) - wp.exp(-tau * k_on)
            if wp.float64(wp.randf(rstate)) < p_att:
                bound[t] = wp.int32(1)


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
def branch_angle_kernel(
    pos: wp.array(dtype=wp.vec3d),
    triples: wp.array(dtype=wp.int32, ndim=2),   # (B, 3) [i, j(center/branch node), k]
    theta0: wp.float64,                          # rest branch angle [rad] (Arp2/3 ≈ 70° = 1.222 rad)
    k_angle: wp.float64,                         # angular stiffness [pN·µm/rad²]
    force: wp.array(dtype=wp.vec3d),
):
    """Harmonic angle force U = ½·k_angle·(θ−θ₀)² at each branch triple (one thread per branch).

    The mechanistic Arp2/3 branch (CLAUDE.md: angle-harmonic with thermal fluctuation, NOT a rigid
    72° constraint): θ is the angle at the branch node j between (i−j) and (k−j). Standard angle force
    F_i = (k(θ−θ₀)/sinθ)·∂cosθ/∂r_i, with F_j = −(F_i+F_k)."""
    t = wp.tid()
    i = triples[t, 0]
    j = triples[t, 1]
    k = triples[t, 2]
    r1 = pos[i] - pos[j]
    r2 = pos[k] - pos[j]
    d1 = wp.length(r1)
    d2 = wp.length(r2)
    if d1 < wp.float64(1e-12) or d2 < wp.float64(1e-12):
        return
    c = wp.dot(r1, r2) / (d1 * d2)
    c = wp.clamp(c, wp.float64(-1.0), wp.float64(1.0))
    s = wp.sqrt(wp.max(wp.float64(1.0) - c * c, wp.float64(1e-12)))
    theta = wp.acos(c)
    pref = k_angle * (theta - theta0) / s
    fi = pref * (r2 / (d1 * d2) - c * r1 / (d1 * d1))
    fk = pref * (r1 / (d1 * d2) - c * r2 / (d2 * d2))
    wp.atomic_add(force, i, fi)
    wp.atomic_add(force, k, fk)
    wp.atomic_add(force, j, -(fi + fk))


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


def branch_angle_force_np(pos, triples, theta0, k_angle, device="cpu"):
    """Convenience: Warp Arp2/3 branch-angle force as a numpy (N,3) array."""
    N = pos.shape[0]
    pos_d = wp.array(np.ascontiguousarray(pos, np.float64), dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    if triples.shape[0]:
        wp.launch(branch_angle_kernel, dim=triples.shape[0],
                  inputs=[pos_d, wp.array(np.ascontiguousarray(triples, np.int32), dtype=wp.int32, device=device),
                          wp.float64(theta0), wp.float64(k_angle), force_d], device=device)
        wp.synchronize_device(device)
    return force_d.numpy().astype(np.float64)


def _branch_angle_force_numpy_ref(pos, triples, theta0, k_angle):
    """Pure-numpy reference for the harmonic branch-angle force (parity oracle for the Warp kernel)."""
    f = np.zeros_like(pos)
    for (i, j, k) in triples:
        r1 = pos[i] - pos[j]; r2 = pos[k] - pos[j]
        d1 = np.linalg.norm(r1); d2 = np.linalg.norm(r2)
        if d1 < 1e-12 or d2 < 1e-12:
            continue
        c = np.clip(np.dot(r1, r2) / (d1 * d2), -1.0, 1.0)
        s = np.sqrt(max(1.0 - c * c, 1e-12))
        theta = np.arccos(c)
        pref = k_angle * (theta - theta0) / s
        fi = pref * (r2 / (d1 * d2) - c * r1 / (d1 * d1))
        fk = pref * (r1 / (d1 * d2) - c * r2 / (d2 * d2))
        f[i] += fi; f[k] += fk; f[j] -= (fi + fk)
    return f


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


@wp.kernel
def plate_kernel(
    pos: wp.array(dtype=wp.vec3d),
    cz: wp.float64,                 # centroid z (plates are symmetric about it)
    half_gap: wp.float64,           # plate half-separation h/2 [µm]
    k_plate: wp.float64,            # plate stiffness [pN/µm]
    force: wp.array(dtype=wp.vec3d),
    reaction: wp.array(dtype=wp.float64),   # [0]=top-plate force (AFM force), [1]=contact node count (top)
):
    """Two rigid parallel plates at z = cz ± half_gap; stiff one-sided penalty on nodes beyond them
    (the virtual-AFM / Fischer-Friedrich parallel-plate compression). Accumulates the TOP-plate reaction
    force reaction[0] = Σ k_plate·penetration (= the force the cell exerts on the plate, what AFM reads)."""
    i = wp.tid()
    z = pos[i][2] - cz
    if z > half_gap:
        pen = k_plate * (z - half_gap)
        wp.atomic_add(force, i, wp.vec3d(wp.float64(0.0), wp.float64(0.0), -pen))
        wp.atomic_add(reaction, 0, pen)
        wp.atomic_add(reaction, 1, wp.float64(1.0))
    elif z < -half_gap:
        pen = k_plate * (-half_gap - z)
        wp.atomic_add(force, i, wp.vec3d(wp.float64(0.0), wp.float64(0.0), pen))


@wp.kernel
def nucleus_shell_kernel(
    pos: wp.array(dtype=wp.vec3d),
    off: wp.int32,                  # nucleus beads occupy indices [off, off+n_nuc)
    centre: wp.vec3d,               # cell centroid (nucleus concentric with the cell)
    R0: wp.float64,                 # rest nuclear radius R_nuc [µm]
    k_chrom: wp.float64,            # inner (chromatin) slope [pN/µm]
    k_lamin: wp.float64,            # extra outer (lamina) slope past the knee [pN/µm]
    d_knee: wp.float64,             # lamin-engagement knee displacement [µm]
    F_knee: wp.float64,             # force at the knee [pN]
    force: wp.array(dtype=wp.vec3d),
):
    """Bilinear radial-shell (nucleus law 0) restoring force per nucleus bead — VERBATIM law-0 math from
    ``common/compartments.rsf_apply`` (chromatin slope inside the knee, +lamin slope past it: strain-
    stiffening). Accumulating (atomic_add) so it composes with the plate penalty on the same beads."""
    j = wp.tid()
    i = off + j
    p = pos[i]
    dx = p[0] - centre[0]; dy = p[1] - centre[1]; dz = p[2] - centre[2]
    r = wp.sqrt(dx * dx + dy * dy + dz * dz)
    rs = r
    if r <= wp.float64(0.0):
        rs = wp.float64(1.0)
    nx = dx / rs; ny = dy / rs; nz = dz / rs
    if r <= wp.float64(0.0):
        nx = wp.float64(0.0); ny = wp.float64(0.0); nz = wp.float64(0.0)
    dd = r - R0
    ad = wp.abs(dd)
    sgn = wp.float64(0.0)
    if dd > wp.float64(0.0):
        sgn = wp.float64(1.0)
    if dd < wp.float64(0.0):
        sgn = wp.float64(-1.0)
    Fmag = wp.float64(0.0)
    if ad <= d_knee:
        Fmag = -k_chrom * dd
    else:
        e = ad - d_knee
        Fmag = -sgn * (F_knee + (k_chrom + k_lamin) * e)
    wp.atomic_add(force, i, wp.vec3d(Fmag * nx, Fmag * ny, Fmag * nz))


def _seed_nucleus_cloud(centre, R_nuc, n_beads, rng):
    """Fibonacci-sphere nucleus bead cloud of radius R_nuc about ``centre`` (near-uniform, deterministic
    up to a small rng jitter so no two beads coincide). Returns (n_beads, 3) float64."""
    k = np.arange(n_beads) + 0.5
    phi = np.arccos(1.0 - 2.0 * k / n_beads)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * k
    u = np.stack([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)], axis=1)
    jitter = 1e-3 * R_nuc * rng.standard_normal((n_beads, 3))
    return (np.asarray(centre, dtype=np.float64) + R_nuc * u + jitter).astype(np.float64)


def simulate_whole_cell_compression_on_device(cortex, f_myo, *, strain=0.0, nucleus=None,
                                              n_steps=4000, reshape_every=25, turgor_every=20, dt_mu=0.0,
                                              n_reshape_iter=2, k_plate=None, pressure_setpoint=None,
                                              nucleus_seed=0, device="cpu"):
    """WHOLE-CELL virtual parallel-plate (AFM) compression: the cortex shell + turgor of
    :func:`simulate_compressed_shell_on_device` PLUS a mechanistic stiff nucleus (shared
    ``common/compartments`` law-0 radial shell). Two physical cortex↔nucleus couplings are modelled:

      1. **incompressible-cytoplasm displacement** — the nucleus occupies volume, so the compressible
         cytoplasm is V_cyto = V_hull − V_nuc; compressing the cell squeezes a smaller cytoplasm →
         ΔP rises faster (the nucleus is felt hydrostatically before any contact), and
      2. **direct compression** — once the plate half-gap drops below R_nuc (strain > 1−R_nuc/R_cell,
         ≈67% for R_nuc/R_cell≈1/3) the plates press on the nucleus beads directly and its stiff
         strain-stiffening law dominates the force response.

    ``nucleus`` is a ``common.compartments.ResolvedNucleus`` (or None → cortex-only, reproducing the
    non-nucleus AFM path). Returns (pos_all, metrics). To measure the nucleus contribution, run twice
    (nucleus=None vs a ResolvedNucleus) and diff F_plate. Validated on the A5000; the CPU path is for
    small-N interface tests."""
    from scipy.spatial import ConvexHull

    from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
    from ffn_sim.ff.gamma_floor import TURGOR_DP0, TURGOR_PI_IN0, VMIN_FRAC

    net = cortex.net
    Nc = net.n_nodes
    tri = np.ascontiguousarray(net.bend_triples, np.int32)
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    fiber_off = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg_per = np.diff(fiber_off) - 1
    seg_off = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    seg = float(net.seg_rest.mean()) if net.seg_rest.size else 1.0
    R0 = cortex.R0_mean if cortex.R0_mean > 0.0 else cortex.R_um

    # ---- seed the nucleus bead cloud concentric with the resting cortex ----
    c0 = net.pos.mean(axis=0)
    if nucleus is not None:
        rng = np.random.default_rng(nucleus_seed)
        nuc_pos = _seed_nucleus_cloud(c0, nucleus.R_nuc_um, nucleus.n_beads, rng)
        V_nuc = (4.0 / 3.0) * np.pi * nucleus.R_nuc_um ** 3      # incompressible nucleus displacement [µm³]
    else:
        nuc_pos = np.zeros((0, 3), dtype=np.float64)
        V_nuc = 0.0
    n_nuc = nuc_pos.shape[0]
    N = Nc + n_nuc
    pos0 = np.concatenate([net.pos, nuc_pos], axis=0)

    # V0 = RESTING compressible-cytoplasm volume (hull of the cortex minus the nucleus it displaces), so
    # V_cyto/V0 = 1 at strain 0 → ΔP = dP0 (same convention as simulate_compressed_shell_on_device).
    V0 = float(ConvexHull(net.pos).volume) - V_nuc
    vmin = VMIN_FRAC * V0
    half_gap = R0 * (1.0 - strain)
    K_vol = TURGOR_PI_IN0 / (1.0 - VMIN_FRAC)
    if k_plate is None:
        k_plate = 10.0 * (K_vol / V0) * (4.0 * np.pi * R0**2) ** 2 / Nc
    if dt_mu <= 0.0:
        kmax = float(net.kappa.max()) / seg**3
        if cortex.xl_i.size:
            kmax = max(kmax, float(cortex.xl_k.max()))
        kmax = max(kmax, (K_vol / V0) * (4.0 * np.pi * R0**2) ** 2 / Nc, k_plate)
        if nucleus is not None:
            kmax = max(kmax, nucleus.k_chrom + nucleus.k_lamin)     # nucleus stiffness sets the CFL too
        dt_mu = 0.1 / kmax

    d = device
    pos_d = wp.array(np.ascontiguousarray(pos0, np.float64), dtype=wp.vec3d, device=d)
    f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    react_d = wp.zeros(3, dtype=wp.float64, device=d)          # [0]=top force, [1]=cortex contacts, [2]=unused
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
    cz = 0.0
    dP_area = 0.0

    def _refresh_turgor():
        nonlocal centre, cz, dP_area
        p = pos_d.numpy().astype(np.float64)
        pcx = p[:Nc]                                  # cortex nodes define the cell surface + centroid
        c = pcx.mean(axis=0)
        centre = wp.vec3d(float(c[0]), float(c[1]), float(c[2]))
        cz = float(c[2])
        try:
            hull = ConvexHull(pcx)
            V, area = float(hull.volume), float(hull.area)
        except Exception:
            R_mean = float(np.linalg.norm(pcx - c, axis=1).mean()); V = (4/3)*np.pi*R_mean**3; area = 4*np.pi*R_mean**2
        V_cyto = V - V_nuc                            # incompressible-nucleus displacement coupling
        if pressure_setpoint is not None:
            dP = float(pressure_setpoint)
        else:
            dP = max(TURGOR_PI_IN0 * (V0 - vmin) / max(V_cyto - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0), 0.0)
        dP_area = dP * area / Nc
        return V, V_cyto, area, dP

    for step in range(n_steps):
        if step % turgor_every == 0:
            _refresh_turgor()
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        if nT:
            wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        if has_xl:
            wp.launch(link_spring_kernel, dim=cortex.xl_i.size, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            wp.launch(myosin_kernel, dim=cortex.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_area), f_d], device=d)  # cortex only
        if n_nuc:
            wp.launch(nucleus_shell_kernel, dim=n_nuc,
                      inputs=[pos_d, wp.int32(Nc), centre, wp.float64(nucleus.R_nuc_um),
                              wp.float64(nucleus.k_chrom), wp.float64(nucleus.k_lamin),
                              wp.float64(nucleus.d_knee_um), wp.float64(nucleus.F_knee_pN), f_d], device=d)
        wp.launch(plate_kernel, dim=N,     # plate acts on ALL nodes (cortex + nucleus)
                  inputs=[pos_d, wp.float64(cz), wp.float64(half_gap), wp.float64(k_plate), f_d, react_d], device=d)
        wp.launch(axpy_kernel, dim=N, inputs=[pos_d, wp.float64(dt_mu), f_d], device=d)
        if (step + 1) % reshape_every == 0:
            wp.launch(reshape_kernel, dim=fiber_off.shape[0] - 1,
                      inputs=[pos_d, foff_d, soff_d, srest_d, wp.int32(n_reshape_iter)], device=d)

    V, V_cyto, area, dP = _refresh_turgor()
    react_d.zero_()
    wp.launch(plate_kernel, dim=N,
              inputs=[pos_d, wp.float64(cz), wp.float64(half_gap), wp.float64(k_plate), f_d, react_d], device=d)
    wp.synchronize_device(d)
    react = react_d.numpy()
    F_plate = float(react[0])
    pos_all = pos_d.numpy().astype(np.float64)
    net.pos = pos_all[:Nc]
    p = net.pos; c = p.mean(axis=0)
    rxy = np.hypot(p[:, 0] - c[0], p[:, 1] - c[1])
    R_eq = float(rxy.max())
    zc = p[:, 2] - c[2]
    contact = np.abs(np.abs(zc) - half_gap) < 0.15 * seg + 0.05
    contact_radius = float(rxy[contact].max()) if contact.any() else 0.0
    gamma_laplace = 0.5 * dP * R_eq
    dP_force = F_plate / (np.pi * contact_radius**2) if contact_radius > 1e-6 else 0.0
    # nucleus geometry: equatorial radius + whether the plate reached the nucleus (direct compression)
    if n_nuc:
        pn = pos_all[Nc:]
        R_nuc_eq = float(np.hypot(pn[:, 0] - c[0], pn[:, 1] - c[1]).max())
        nucleus_contact = bool((np.abs(pn[:, 2] - cz) > half_gap).any())
    else:
        R_nuc_eq = 0.0
        nucleus_contact = False
    metrics = {"strain": strain, "half_gap": half_gap, "F_plate_pN": F_plate,
               "R_eq_um": R_eq, "contact_radius_um": contact_radius, "dP_turgor_Pa": dP,
               "dP_from_force_Pa": dP_force, "V_over_V0": V_cyto / V0, "V_hull_um3": V, "V_nuc_um3": V_nuc,
               "gamma_apparent_pN_um": gamma_laplace, "gamma_apparent_mN_m": gamma_laplace * 1e-3,
               "R_nuc_eq_um": R_nuc_eq, "nucleus_contact": nucleus_contact, "n_nuc": n_nuc,
               "k_plate": k_plate}
    return pos_all, metrics


def relax_on_device(net, *, links=None, k_xl=None, xl_rest=None, myo_links=None, f_myo=0.0,
                    branch_triples=None, branch_theta0=0.0, branch_k=0.0,
                    n_steps=600, reshape_every=25, dt_mu=0.0, n_reshape_iter=2, device="cpu"):
    """Fully on-device FF cortex relaxation: bending (+ optional crosslink/myosin/Arp2/3-branch) forces +
    reshape, all Warp — NO per-step numpy round-trip. Runs on CPU (Mac) or CUDA (gbook A5000) via ``device``.

    Returns the relaxed positions (N,3) numpy. Inextensibility via periodic reshape (the robust path;
    no singular-prone projector). ``links``/``k_xl``/``xl_rest`` optional crosslinker springs;
    ``myo_links``/``f_myo`` optional myosin; ``branch_triples``/``branch_theta0``/``branch_k`` optional
    angle-harmonic Arp2/3 branches (lamellipodium). Default (none) = resting-shell bending settle.
    """
    from ffn_sim.ff.forces_warp import _per_triple_alpha
    N = net.n_nodes
    tri = np.ascontiguousarray(net.bend_triples, np.int32)
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    fiber_off = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg_per = np.diff(fiber_off) - 1
    seg_off = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    seg = float(net.seg_rest.mean()) if net.seg_rest.size else 1.0
    has_branch = branch_triples is not None and len(branch_triples) and branch_k > 0.0
    if dt_mu <= 0.0:
        kmax = float(net.kappa.max()) / seg**3
        if k_xl is not None and len(k_xl):
            kmax = max(kmax, float(np.max(k_xl)))
        if has_branch:                                         # CFL: branch angular spring k_eff ≈ k_angle/ℓ²
            kmax = max(kmax, branch_k / seg**2)
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
    if has_branch:
        br_d = wp.array(np.ascontiguousarray(branch_triples, np.int32), dtype=wp.int32, device=d)
    nT = tri.shape[0]
    for step in range(n_steps):
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        if nT:
            wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        if has_xl:
            wp.launch(link_spring_kernel, dim=len(links), inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            wp.launch(myosin_kernel, dim=len(myo_links), inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        if has_branch:
            wp.launch(branch_angle_kernel, dim=len(branch_triples),
                      inputs=[pos_d, br_d, wp.float64(branch_theta0), wp.float64(branch_k), f_d], device=d)
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


def simulate_compressed_shell_on_device(cortex, f_myo, *, strain=0.0, n_steps=4000, reshape_every=25,
                                        turgor_every=20, dt_mu=0.0, n_reshape_iter=2, k_plate=None,
                                        pressure_setpoint=None, device="cpu"):
    """VIRTUAL PARALLEL-PLATE (AFM / Fischer-Friedrich) compression of the loaded cortex shell.

    Replicates the EXPERIMENTAL cortical-tension PROTOCOL (measurement-protocol-consistency sanity gate):
    the band papers (Chugh, Fischer-Friedrich, Hosseini) measure γ by DEFORMING the cell (AFM /
    parallel-plate / aspiration), not on a resting cortex. Here two rigid plates at z = centroid ±
    R0·(1−strain) compress the turgor-pressurised cortex; we read the plate reaction force and extract the
    apparent surface tension the SAME way the experiment does (liquid-drop Young-Laplace), so it can be
    compared apples-to-apples with the band — unlike the resting method-of-planes γ_myo.

    The compressed shape is NON-spherical, so the turgor closure uses the true enclosed volume (convex
    hull) instead of the spherical R_mean³. Returns (pos, metrics) where metrics has the plate force,
    contact/equatorial geometry, ΔP, V/V0, and the extracted γ_apparent [pN/µm].
    """
    from scipy.spatial import ConvexHull

    from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
    from ffn_sim.ff.gamma_floor import TURGOR_DP0, TURGOR_PI_IN0, VMIN_FRAC

    net = cortex.net
    N = net.n_nodes
    tri = np.ascontiguousarray(net.bend_triples, np.int32)
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    fiber_off = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg_per = np.diff(fiber_off) - 1
    seg_off = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    seg = float(net.seg_rest.mean()) if net.seg_rest.size else 1.0
    R0 = cortex.R0_mean if cortex.R0_mean > 0.0 else cortex.R_um
    # V0 MUST use the SAME volume convention as the running V (convex hull of the actual points), NOT the
    # sphere (4/3)πR0³: a ~few-% V0-vs-hull mismatch × the huge osmotic modulus (K_vol~7e5 Pa) makes ΔP a
    # pure volume-reference ARTIFACT (frustrated cortex → hull<sphere → ΔP huge; relaxed → hull>sphere → ΔP=0).
    # Reference to the RESTING (uncompressed) hull so V/V0=1 at strain 0 → ΔP=dP0.
    V0 = float(ConvexHull(net.pos).volume)
    vmin = VMIN_FRAC * V0
    half_gap = R0 * (1.0 - strain)                       # plate half-separation
    K_vol = TURGOR_PI_IN0 / (1.0 - VMIN_FRAC)
    if k_plate is None:
        k_plate = 10.0 * (K_vol / V0) * (4.0 * np.pi * R0**2) ** 2 / N   # 10× the turgor breathing mode → rigid plate
    if dt_mu <= 0.0:
        kmax = float(net.kappa.max()) / seg**3
        if cortex.xl_i.size:
            kmax = max(kmax, float(cortex.xl_k.max()))
        kmax = max(kmax, (K_vol / V0) * (4.0 * np.pi * R0**2) ** 2 / N, k_plate)
        dt_mu = 0.1 / kmax

    d = device
    pos_d = wp.array(np.ascontiguousarray(net.pos, np.float64), dtype=wp.vec3d, device=d)
    f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    react_d = wp.zeros(2, dtype=wp.float64, device=d)
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
    cz = 0.0
    dP_area = 0.0

    def _refresh_turgor():
        nonlocal centre, cz, dP_area
        p = pos_d.numpy().astype(np.float64)
        c = p.mean(axis=0)
        centre = wp.vec3d(float(c[0]), float(c[1]), float(c[2]))
        cz = float(c[2])
        try:
            hull = ConvexHull(p)
            V, area = float(hull.volume), float(hull.area)
        except Exception:
            R_mean = float(np.linalg.norm(p - c, axis=1).mean()); V = (4/3)*np.pi*R_mean**3; area = 4*np.pi*R_mean**2
        if pressure_setpoint is not None:
            # VOLUME REGULATION (water flux / osmoregulation) — the cell holds ΔP at a setpoint by letting
            # water/ions cross under excess pressure (Kedem-Katchalsky pump-leak). This is the perfect-
            # regulation limit (fast Lp); the real cell (finite Lp) is between this and the fixed-osmolyte
            # closure. A regulated cell behaves as a liquid drop of confinement-INDEPENDENT surface tension
            # (γ_apparent = ΔP_setpoint·R_eq/2), matching Fischer-Friedrich's confinement-independent γ.
            dP = float(pressure_setpoint)
        else:
            dP = max(TURGOR_PI_IN0 * (V0 - vmin) / max(V - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0), 0.0)
        dP_area = dP * area / N
        return V, area, dP

    for step in range(n_steps):
        if step % turgor_every == 0:
            _refresh_turgor()
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        if nT:
            wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        if has_xl:
            wp.launch(link_spring_kernel, dim=cortex.xl_i.size, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            wp.launch(myosin_kernel, dim=cortex.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(turgor_kernel, dim=N, inputs=[pos_d, centre, wp.float64(dP_area), f_d], device=d)
        wp.launch(plate_kernel, dim=N,
                  inputs=[pos_d, wp.float64(cz), wp.float64(half_gap), wp.float64(k_plate), f_d, react_d], device=d)
        wp.launch(axpy_kernel, dim=N, inputs=[pos_d, wp.float64(dt_mu), f_d], device=d)
        if (step + 1) % reshape_every == 0:
            wp.launch(reshape_kernel, dim=fiber_off.shape[0] - 1,
                      inputs=[pos_d, foff_d, soff_d, srest_d, wp.int32(n_reshape_iter)], device=d)

    # final measurement: plate reaction force + shape, on the equilibrated compressed shell
    V, area, dP = _refresh_turgor()
    react_d.zero_()
    wp.launch(plate_kernel, dim=N,
              inputs=[pos_d, wp.float64(cz), wp.float64(half_gap), wp.float64(k_plate), f_d, react_d], device=d)
    wp.synchronize_device(d)
    react = react_d.numpy()
    F_plate = float(react[0])                       # top-plate force = AFM force [pN]
    n_contact = float(react[1])
    net.pos = pos_d.numpy().astype(np.float64)
    p = net.pos; c = p.mean(axis=0)
    rxy = np.hypot(p[:, 0] - c[0], p[:, 1] - c[1])
    R_eq = float(rxy.max())                          # equatorial (bulge) radius [µm]
    zc = p[:, 2] - c[2]
    contact = np.abs(np.abs(zc) - half_gap) < 0.15 * seg + 0.05
    contact_radius = float(rxy[contact].max()) if contact.any() else 0.0
    gamma_laplace = 0.5 * dP * R_eq                  # liquid-drop Young-Laplace apparent tension [pN/µm]
    dP_force = F_plate / (np.pi * contact_radius**2) if contact_radius > 1e-6 else 0.0  # pressure from AFM force
    metrics = {"strain": strain, "half_gap": half_gap, "F_plate_pN": F_plate, "n_contact": n_contact,
               "R_eq_um": R_eq, "contact_radius_um": contact_radius, "dP_turgor_Pa": dP,
               "dP_from_force_Pa": dP_force, "V_over_V0": V / V0,
               "gamma_apparent_pN_um": gamma_laplace, "gamma_apparent_mN_m": gamma_laplace * 1e-3,
               "k_plate": k_plate}
    return net.pos, metrics


def simulate_turnover_on_device(cortex, f_myo, *, n_steps=8000, reshape_every=25, kmc_every=50,
                                tau_kmc=0.01, record_every=1000, dt_mu=0.0, n_reshape_iter=2,
                                seed0=12345, device="cpu"):
    """Evolve the resting cortex with MYOSIN HAND TURNOVER on-device (Warp): bending + crosslink
    springs (static, ~99 % bound since k_on≫k_off) + ENGAGED myosin (bound-aware) + reshape, with
    periodic Bell detach/re-attach KMC (``kmc_bell_turnover_kernel``) on the myosin population.

    This closes the γ-floor robustness claim on the TURNOVER axis (the static buckling/connectivity/
    extensibility ablations are in FF_STAGE6H): does load-dependent motor turnover lift the floor?
    With Bell detachment the engaged myosin fraction self-limits under load → the floor cannot rise
    (the engaged force-bearing density, not the network mechanism, sets γ). Returns (pos, bound_mask,
    traj) where traj records {step, bound_frac} and bound_mask is the final engaged state.

    ``tau_kmc`` is the PHYSICAL KMC tick [s] (default 0.01) — decoupled from the mechanical overdamped
    descent pseudo-step ``dt_mu`` (regime A/B: the mechanics re-equilibrates over ``kmc_every`` descent
    steps between ticks, then one Hand KMC tick advances physical time by ``tau_kmc``). For small
    ``tau_kmc`` the steady engaged fraction → k_on/(k_on+p_off), τ-independent."""
    from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
    from ffn_sim.ff.gamma_floor import mesoscale_reach
    from ffn_sim.ff.hand_kmc import NMIIA_MYOSIN

    net = cortex.net
    N = net.n_nodes
    tri = np.ascontiguousarray(net.bend_triples, np.int32)
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    fiber_off = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg_per = np.diff(fiber_off) - 1
    seg_off = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    seg = float(net.seg_rest.mean()) if net.seg_rest.size else 1.0
    if dt_mu <= 0.0:
        kmax = float(net.kappa.max()) / seg**3
        if cortex.xl_i.size:
            kmax = max(kmax, float(cortex.xl_k.max()))
        dt_mu = 0.1 / kmax

    d = device
    pos_d = wp.array(np.ascontiguousarray(net.pos, np.float64), dtype=wp.vec3d, device=d)
    f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
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
    M = int(cortex.myo_i.size)
    has_myo = M > 0 and f_myo != 0.0
    if has_myo:
        myo_d = wp.array(np.ascontiguousarray(np.stack([cortex.myo_i, cortex.myo_j], 1), np.int32),
                         dtype=wp.int32, device=d)
        bound_d = wp.array(np.ones(M, np.int32), dtype=wp.int32, device=d)  # start engaged
    nT = tri.shape[0]
    p0 = float(NMIIA_MYOSIN.p0); f0 = float(NMIIA_MYOSIN.f0); k_on = float(NMIIA_MYOSIN.k_on)
    # Re-attach capture = the MESOSCALE reach √(A/n) the links were formed at (sanctioned ×40 dual,
    # mesoscale_reach) — NOT the molecular ε (210 nm), which is far below the coarse link length so a
    # detached mesoscale hand could never rebind. Slack ×1.5 lets a fluctuating partner re-bind.
    cap = 1.5 * mesoscale_reach(cortex.R_um, net.n_fibers)
    traj = []

    def _bf(step):
        bf = float(bound_d.numpy().mean()) if has_myo else 0.0
        traj.append({"step": step, "bound_frac": bf})
        return bf

    for step in range(n_steps):
        if has_myo and step % kmc_every == 0:
            wp.launch(kmc_bell_turnover_kernel, dim=M,
                      inputs=[pos_d, myo_d, bound_d, wp.float64(f_myo), wp.float64(p0),
                              wp.float64(f0), wp.float64(cap), wp.float64(k_on),
                              wp.float64(tau_kmc), wp.int32(seed0 + step)], device=d)
            if step % record_every == 0:
                _bf(step)
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        if nT:
            wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        if has_xl:
            wp.launch(link_spring_kernel, dim=cortex.xl_i.size, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            wp.launch(myosin_bound_kernel, dim=M, inputs=[pos_d, myo_d, bound_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(axpy_kernel, dim=N, inputs=[pos_d, wp.float64(dt_mu), f_d], device=d)
        if (step + 1) % reshape_every == 0:
            wp.launch(reshape_kernel, dim=fiber_off.shape[0] - 1,
                      inputs=[pos_d, foff_d, soff_d, srest_d, wp.int32(n_reshape_iter)], device=d)
    _bf(n_steps)
    wp.synchronize_device(d)
    net.pos = pos_d.numpy().astype(np.float64)
    bound_mask = bound_d.numpy().astype(bool) if has_myo else np.zeros(0, bool)
    return net.pos, bound_mask, traj
