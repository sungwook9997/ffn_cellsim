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
def wlc_spring_kernel(
    pos: wp.array(dtype=wp.vec3d),              # (N,) node positions
    seg: wp.array(dtype=wp.int32, ndim=2),      # (S, 2) fiber-segment node pairs
    Lc: wp.array(dtype=wp.float64),             # (S,) per-segment contour length [µm]
    Lp: wp.float64,                             # persistence length [µm]
    EA: wp.float64,                             # axial stiffness [pN] (enthalpic backbone wall)
    kBT: wp.float64,                            # thermal energy [pN·µm]
    x_max: wp.float64,                          # entropic→enthalpic crossover (f'_WLC=EA)
    force: wp.array(dtype=wp.vec3d),            # (N,) out
):
    """Thermal extensible-WLC segment tension (mirrors ff.wlc.wlc_tension_np exactly): Marko-Siggia entropic
    F=(kBT/Lp)[1/(4(1−x)²)−1/4+x] for x=L/Lc ≤ x_max, else the finite EA wall F(x_max)+(EA/Lc)(L−x_max·Lc)
    (C¹-continuous, no 1/(1−x)² blow-up). FIBER SEGMENTS ONLY — crosslinks keep link_spring_kernel."""
    t = wp.tid()
    i = seg[t, 0]
    j = seg[t, 1]
    d = pos[j] - pos[i]
    L = wp.length(d)
    if L > wp.float64(1e-12):
        lc = Lc[t]
        x = L / lc
        one = wp.float64(1.0)
        qtr = wp.float64(0.25)
        if x <= x_max:
            om = one - x
            F = (kBT / Lp) * (one / (wp.float64(4.0) * om * om) - qtr + x)
        else:
            omm = one - x_max
            F_max = (kBT / Lp) * (one / (wp.float64(4.0) * omm * omm) - qtr + x_max)
            F = F_max + (EA / lc) * (L - x_max * lc)
        f = (F / L) * d
        wp.atomic_add(force, i, f)
        wp.atomic_add(force, j, -f)


@wp.kernel
def xl_turnover_kernel(
    pos: wp.array(dtype=wp.vec3d),              # (N,) node positions
    links: wp.array(dtype=wp.int32, ndim=2),    # (L, 2) crosslinker node pairs
    k_arr: wp.array(dtype=wp.float64),          # (L,) stiffness [pN/µm]
    r0_arr: wp.array(dtype=wp.float64),         # (L,) rest length [µm] — UPDATED IN PLACE
    koff0_dt: wp.float64,                        # k_off0 · dt_real per turnover event [dimensionless]
    x_beta_over_kT: wp.float64,                  # Bell-slip force scale x_β/kT [1/pN]
):
    """Viscoelastic crosslinker turnover (Bell slip, NF2007 §10.1 / Ferrer 2008 α-actinin): each event a
    fraction 1−exp(−k_off(F)·dt) of the crosslinker ensemble unbinds and rebinds FORCE-FREE at the current
    geometry, so the rest length creeps toward the current length (stress relaxation). Force-dependent
    off-rate k_off(F)=k_off0·exp(|F|·x_β/kT). dt≪1/k_off → elastic (stiff, fast AFM); dt≫1/k_off → relaxed
    (cortex remodels, creep). Makes each crosslinker a Maxwell element with relaxation time 1/k_off(F)."""
    t = wp.tid()
    i = links[t, 0]
    j = links[t, 1]
    L = wp.length(pos[j] - pos[i])
    F = k_arr[t] * (L - r0_arr[t])                                  # signed crosslinker force [pN]
    koff_dt = koff0_dt * wp.exp(wp.abs(F) * x_beta_over_kT)         # Bell slip, per-event
    frac = wp.float64(1.0) - wp.exp(-koff_dt)                       # ensemble fraction turned over
    r0_arr[t] = r0_arr[t] + frac * (L - r0_arr[t])                  # rebind force-free → rest length → current


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
def soft_contact_kernel(
    pos: wp.array(dtype=wp.vec3d),
    pairs: wp.array(dtype=wp.int32, ndim=2),    # (C, 2) contact pairs (i on structure A, j on structure B)
    r_contact: wp.float64,                       # contact radius [µm] (excluded-volume shell)
    k_contact: wp.float64,                       # contact stiffness [pN/µm]
    force: wp.array(dtype=wp.vec3d),
):
    """One-sided REPULSIVE (excluded-volume) contact: if |pos_j − pos_i| < r_contact, push the two nodes
    apart with k·(r_contact − L)·û. No attraction. Models steric contact (a growing filopodium tip pushing
    into an ECM fiber) — the ECM then supplies the protrusion's load EMERGENTLY, not as an imposed number."""
    t = wp.tid()
    i = pairs[t, 0]
    j = pairs[t, 1]
    d = pos[j] - pos[i]
    L = wp.length(d)
    if L > wp.float64(1e-12) and L < r_contact:
        f = (k_contact * (r_contact - L) / L) * d
        wp.atomic_add(force, i, -f)              # i pushed away from j
        wp.atomic_add(force, j, f)               # j pushed away from i


@wp.kernel
def substrate_plane_kernel(
    pos: wp.array(dtype=wp.vec3d),
    z_sub: wp.float64,                          # substrate plane height [µm]
    k_plane: wp.float64,                        # plane stiffness [pN/µm]
    force: wp.array(dtype=wp.vec3d),
):
    """One-sided excluded-volume substrate: a node below the plane (z<z_sub) is pushed UP with
    k_plane·(z_sub−z). The cell rests ON the substrate (does not sink through it); adhesion is the separate
    FA-clutch layer. Node index range is set by the launch dim (cortex nodes only)."""
    i = wp.tid()
    if pos[i][2] < z_sub:
        wp.atomic_add(force, i, wp.vec3d(0.0, 0.0, k_plane * (z_sub - pos[i][2])))


@wp.kernel
def freeze_kernel(force: wp.array(dtype=wp.vec3d), pinned: wp.array(dtype=wp.int32)):
    """Zero the net force on pinned (Dirichlet-BC) nodes so the integrator leaves them fixed — e.g. the far
    boundary of the ECM Mikado network (embedded in bulk matrix) or a clamped filopodium base."""
    t = wp.tid()
    if pinned[t] == 1:
        force[t] = wp.vec3d(0.0, 0.0, 0.0)


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
def _centroid_sum_kernel(
    pos: wp.array(dtype=wp.vec3d),
    off: wp.int32,                # sum pos[off : off+dim]
    out: wp.array(dtype=wp.float64),   # length 3 accumulator
):
    """GPU centroid reduction: Σ pos[off+i] → out[0..2] (divide by n on host)."""
    i = wp.tid()
    p = pos[off + i]
    wp.atomic_add(out, 0, p[0])
    wp.atomic_add(out, 1, p[1])
    wp.atomic_add(out, 2, p[2])


@wp.kernel
def _mesh_vol_area_kernel(
    pos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),   # (nf,3) OUTWARD-oriented triangle indices (local to the body)
    off: wp.int32,                              # node-index offset (0 cortex, Ne nucleus)
    out: wp.array(dtype=wp.float64),            # length 2: out[0]=6·V (divergence), out[1]=2·A
):
    """GPU enclosed-volume + area via the divergence theorem over a FIXED outward-oriented triangulation —
    replaces the per-refresh scipy ConvexHull + full pos→CPU sync (GPU-main). Exact for the deforming mesh."""
    t = wp.tid()
    a = pos[off + faces[t, 0]]
    b = pos[off + faces[t, 1]]
    c = pos[off + faces[t, 2]]
    wp.atomic_add(out, 0, wp.dot(a, wp.cross(b, c)))          # 6·(signed volume about the origin)
    wp.atomic_add(out, 1, wp.length(wp.cross(b - a, c - a)))  # 2·(triangle area)


@wp.kernel
def patch_force_kernel(
    pos: wp.array(dtype=wp.vec3d),
    centre: wp.vec3d,
    axis: wp.vec3d,               # unit patch axis (e.g. +z = top pole)
    cos_thresh: wp.float64,       # node in patch if dot(r_hat, axis) >= cos_thresh (spherical cap)
    f_node: wp.float64,           # per-node force along axis [pN]: >0 outward (tension, S2-2B), <0 inward (compression, S2-2A)
    force: wp.array(dtype=wp.vec3d),
):
    """S2 localized patch load: force f_node·axis on surface nodes inside the cap around `axis`."""
    i = wp.tid()
    d = pos[i] - centre
    r = wp.length(d)
    if r > wp.float64(1e-12):
        if (wp.dot(d, axis) / r) >= cos_thresh:
            wp.atomic_add(force, i, f_node * axis)


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
def add_force_kernel(F: wp.array(dtype=wp.vec3d), add: wp.array(dtype=wp.vec3d)):
    """Accumulate a precomputed per-node force ``add`` into the force accumulator ``F`` (F[i] += add[i]).

    Used by the FSI path to inject the cached −∇p Darcy back-reaction (interpolated from the Biot
    pore-pressure field once per physical step) into the inner mechanical-equilibrium relaxation."""
    i = wp.tid()
    F[i] = F[i] + add[i]


@wp.kernel
def fsi_compression_source_kernel(
    pos: wp.array(dtype=wp.vec3d),
    pos_prev: wp.array(dtype=wp.vec3d),
    centre: wp.vec3d,
    M_over_R: wp.float64,       # Biot coupling M/R0 [Pa/µm] — pressure per unit inward displacement
    src_out: wp.array(dtype=wp.vec3d),   # (Δp_i, 1.0, 0): value + partition-of-unity weight for the IBM deposit
):
    """Per-cortex-node pore-pressure INCREMENT this physical step Δp_i = M·(local inward strain δ_in/R0).

    A node moving INWARD (toward the cell centroid) compresses the local cytoplasm → deposits positive
    excess pore pressure there (Terzaghi/Biot: ∂p = +M·ε_compress, ε≈δ_in/R). Outward motion is clamped
    to 0 (the field relaxes by diffusion, not by a negative source). Component 0 carries Δp_i, component 1
    carries the weight 1.0 — so ONE Peskin ``ibm_spread`` gives both Σδ·Δp (numerator) and Σδ (denominator)
    and the grid deposit is the mass-matrix-normalized FIELD Δp(x)=Σδ·Δp/Σδ [Pa], not a raw density."""
    i = wp.tid()
    d = pos[i] - pos_prev[i]
    rvec = pos[i] - centre
    rn = wp.length(rvec)
    if rn < wp.float64(1e-9):
        src_out[i] = wp.vec3d(wp.float64(0.0), wp.float64(1.0), wp.float64(0.0))
        return
    rhat = rvec / rn
    inward = -wp.dot(d, rhat)                 # >0 ⇒ moved inward this physical step (compression)
    dp = M_over_R * wp.max(inward, wp.float64(0.0))
    src_out[i] = wp.vec3d(dp, wp.float64(1.0), wp.float64(0.0))


@wp.kernel
def fsi_pressure_force_kernel(
    pos: wp.array(dtype=wp.vec3d),
    centre: wp.vec3d,
    p_node: wp.array(dtype=wp.vec3d),   # per-node interpolated excess pore pressure (Pa) in component 0
    area_per_node: wp.float64,          # A_cortex / Nc [µm²]
    force: wp.array(dtype=wp.vec3d),
):
    """Spatially-varying pore-pressure back-reaction on the cortex SHELL: outward normal force
    p_i·(area/node)·r̂ per cortex node (the same Young-Laplace shell form as ``turgor_kernel``, but the
    pressure is the LOCAL interpolated Biot excess p_i instead of the uniform ΔP). p_i is a ZERO-MEAN
    spatial deviation ⇒ this adds only the spatial STRUCTURE (the uniform turgor stays in its own channel;
    no double-count). Units: Pa·µm² = pN."""
    i = wp.tid()
    d = pos[i] - centre
    r = wp.length(d)
    if r > wp.float64(1e-12):
        dP_area = p_node[i][0] * area_per_node
        wp.atomic_add(force, i, (dP_area / r) * d)


@wp.kernel
def fsi_nucleus_pressure_force_kernel(
    pos: wp.array(dtype=wp.vec3d),
    off: wp.int32,                      # nucleus beads live at [off, off+n_nuc)
    centre_nuc: wp.vec3d,               # nucleus centroid (its own body)
    p_node: wp.array(dtype=wp.vec3d),   # per-node interpolated cytoplasm excess pore pressure (Pa) in comp 0
    area_per_bead: wp.float64,          # A_nuc / n_nuc [µm²]
    force: wp.array(dtype=wp.vec3d),
):
    """Cytoplasm pore-pressure back-reaction on the NUCLEUS shell — the fluid→interior channel that
    unfreezes the nucleus. The nucleus is an INCLUSION (cytoplasm is EXTERIOR), so the pressure pushes its
    surface INWARD: −p_i·(A_nuc/n_nuc)·r̂_nuc per bead (note the MINUS — opposite the cortex, which the fluid
    fills from inside and pushes OUTWARD). p_i is the local ZERO-MEAN cytoplasm deviation, so the nucleus
    feels only the spatial (quadrupole) variation → it FLATTENS (poles-in / equator-out) at constant V_nuc
    (nucleus_pressure_kernel conserves volume). Physically DISJOINT from nucleus_shell_kernel (own
    chromatin/lamina elasticity) and nucleus_pressure_kernel (nucleoplasm interior) — no double-count.
    Momentum: a thermodynamic pressure pair with the cortex outward force, net ≈ 0 for the symmetric
    compression + zero-mean field (NOT strict Newtonian — that needs the u-p velocity field). Units: Pa·µm²=pN."""
    j = wp.tid()
    i = off + j
    d = pos[i] - centre_nuc
    r = wp.length(d)
    if r > wp.float64(1e-12):
        dP_area = p_node[i][0] * area_per_bead
        wp.atomic_add(force, i, (-dP_area / r) * d)             # INWARD (inclusion): cytoplasm squeezes the nucleus


@wp.kernel
def fsi_fiber_gradp_force_kernel(
    off: wp.int32,                      # fiber (MT) beads live at [off, off+n_fib)
    gnode: wp.array(dtype=wp.vec3d),    # per-node interpolated ∇p (Pa/µm)
    vseg: wp.array(dtype=wp.float64),   # displaced tube volume per bead π·a²·ℓ [µm³] (0 at the MTOC)
    force: wp.array(dtype=wp.vec3d),
):
    """Pore-pressure-gradient (buoyancy) body force on a 1-D immersed fiber (MT): −V_seg·∇p per bead. A
    fiber has NO enclosed surface, so it feels the fluid as the −∇p force on its displaced tube volume
    (the Darcy quasi-static limit where the fluid velocity is slaved to −∇p — no independent velocity DOF
    needed). HONEST: with a=12.5nm the displaced volume is tiny → this force is ~1e-3 pN, ~1000× below the
    MT bending/contact forces, so it UNFREEZES the MT numerically but is NOT its dominant coupling (that is
    the strut tip↔cortex contact + IF/LINC). Do NOT inflate a/V_seg to force motion (no-param-tuning). A
    genuine along-fiber viscous drag f=(η_f/k)(v_s−v_f) needs the u-p velocity extension (PI-authored KU)."""
    t = wp.tid()
    i = off + t
    wp.atomic_add(force, i, -vseg[t] * gnode[i])


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


@wp.kernel
def nucleus_pressure_kernel(
    pos: wp.array(dtype=wp.vec3d),
    off: wp.int32,                                # nucleus beads live at [off, off+n_nuc)
    centre_nuc: wp.vec3d,                         # nucleus centroid (its own body)
    dP_area: wp.float64,                          # nucleoplasm ΔP·(area/bead) [pN]
    force: wp.array(dtype=wp.vec3d),
):
    """Nucleus volume-conservation (incompressibility): radial outward pressure ΔP·(area/bead)·r̂ about the
    NUCLEUS centroid, ΔP = K_vol·(V0_nuc−V_nuc)/V0_nuc. When the plate flattens the nucleus (V<V0), ΔP>0
    pushes beads out — the free equator bulges, the clamped poles add to the plate reaction (stiffer, ν→½)."""
    j = wp.tid()
    i = off + j
    d = pos[i] - centre_nuc
    r = wp.length(d)
    if r > wp.float64(1e-12):
        wp.atomic_add(force, i, (dP_area / r) * d)


def _oriented_hull_faces(pts):
    """ConvexHull triangulation of ``pts`` (N,3), re-wound so every triangle normal points OUTWARD.
    Used once at build to give the GPU divergence-volume kernel a fixed, consistently-oriented mesh
    (ConvexHull.simplices come in arbitrary winding → signed volumes would cancel). Returns (nf,3) int32."""
    from scipy.spatial import ConvexHull
    h = ConvexHull(pts)
    faces = h.simplices.astype(np.int32).copy()
    a, b, c = pts[faces[:, 0]], pts[faces[:, 1]], pts[faces[:, 2]]
    flip = np.einsum("ij,ij->i", np.cross(b - a, c - a), h.equations[:, :3]) < 0.0
    faces[flip, 1], faces[flip, 2] = faces[flip, 2], faces[flip, 1].copy()
    return faces


def _seed_nucleus_cloud(centre, R_nuc, n_beads, rng):
    """Fibonacci-sphere nucleus bead cloud of radius R_nuc about ``centre`` (near-uniform, deterministic
    up to a small rng jitter so no two beads coincide). Returns (n_beads, 3) float64."""
    k = np.arange(n_beads) + 0.5
    phi = np.arccos(1.0 - 2.0 * k / n_beads)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * k
    u = np.stack([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)], axis=1)
    jitter = 1e-3 * R_nuc * rng.standard_normal((n_beads, 3))
    return (np.asarray(centre, dtype=np.float64) + R_nuc * u + jitter).astype(np.float64)


def simulate_whole_cell_compression_on_device(cortex, f_myo, *, strain=0.0, nucleus=None, membrane=None,
                                              microtubule=None, k_mt_contact=None,
                                              n_steps=4000, reshape_every=25, turgor_every=20, dt_mu=0.0,
                                              n_reshape_iter=2, k_plate=None, pressure_setpoint=None,
                                              Lp_um_s_Pa=None, K_drained_Pa=None, load_time_s=None,
                                              xl_koff_per_s=None, xl_x_beta_um=4.0e-4,
                                              v_press_um_s=None, dwell_steps=0, eta_bulk_Pa_s=None,
                                              rigid_plate=False, nucleus_seed=0, trace_every=0,
                                              local_load=None, biot_fsi=None, if_cage=None,
                                              membrane_surface=None, device="cpu"):
    """WHOLE-CELL virtual parallel-plate (AFM) compression: the cortex shell + turgor of
    :func:`simulate_compressed_shell_on_device` PLUS a mechanistic stiff nucleus (shared
    ``common/compartments`` law-0 radial shell) PLUS an optional plasma-membrane surface (law-1). Two
    physical cortex↔nucleus couplings are modelled:

      1. **incompressible-cytoplasm displacement** — the nucleus occupies volume, so the compressible
         cytoplasm is V_cyto = V_hull − V_nuc; compressing the cell squeezes a smaller cytoplasm →
         ΔP rises faster (the nucleus is felt hydrostatically before any contact), and
      2. **direct compression** — once the plate half-gap drops below R_nuc (strain > 1−R_nuc/R_cell,
         ≈67% for R_nuc/R_cell≈1/3) the plates press on the nucleus beads directly and its stiff
         strain-stiffening law dominates the force response.

    ``nucleus`` is a ``common.compartments.ResolvedNucleus`` (or None → no nucleus). ``membrane`` is a
    ``common.compartments.ResolvedMembrane`` (or None → no membrane): a distinct plasma-membrane surface
    tension applied INWARD (Young-Laplace ΔP_mem=2γ_mem/R) on the cortex nodes — an ADDITIVE lipid channel the
    FF cortex lacks (the cortex has no area-elastic term; reshape fixes segment length only). Wired as the
    **reservoir-buffered CONSTANT baseline tension γ_mem** (physiological state, Raucher-Sheetz plateau); the
    steep K_A elastic upturn past reservoir capacity is deferred (PI-blocked f_excess). γ_mem is a separate
    force channel; ΔP·R/2 is NOT re-reported as passive tension. Returns (pos_all, metrics). Validated on the
    A5000; the CPU path is for small-N interface tests.

    ``biot_fsi`` (default None → OFF, bit-identical monolithic loop = the regression gate) turns on the
    spatially-resolved Biot pore-pressure FSI (FSI_WIRING_DESIGN_2026-07-16.md). A dict with optional keys
    ``D_um2_s`` (Moeendarbary poroelastic diffusivity, def 50), ``M_biot_Pa`` (Biot coupling modulus, def
    K_drained or 300), ``dx_um`` (fluid grid, def 0.6), ``inner_steps`` (mech descents/physical step, def 40),
    ``n_phys`` (physical steps over the ramp, def 300), ``field_out_npz`` (path to save the pore-pressure
    field for the CFD-cell HTML). ON: an OUTER physical-time loop advances a Biot field (τ_p~L²/D) sourced by
    local cortex compression; an INNER mechanical relaxation equilibrates the skeleton each physical step
    (γ_solid drag ⇒ quasi-static, so 6πηR is RETIRED as the clock and the effective η emerges from D+M). The
    back-reaction is the spatially-varying pore pressure on the cortex SHELL (local p_i·area·r̂); the field
    carries the ZERO-MEAN spatial deviation so the uniform turgor stays in its own channel (no double-count).
    Requires ``v_press_um_s`` (the physical clock). Same (pos_all, metrics) return; metrics carry ``fsi_*``."""
    from scipy.spatial import ConvexHull

    from aleph.laws.forces_warp import _per_triple_alpha, cytosim_bending_kernel
    from aleph.laws.turgor_constants import TURGOR_DP0, TURGOR_PI_IN0, VMIN_FRAC

    net = cortex.net                                  # cortex net (surface/turgor/hull always use THIS)
    Nc = net.n_nodes
    # ---- Thread-C: MERGE the MT aster into the elastic network (additive; microtubule=None → bit-identical) ----
    # enet = [cortex(Nc) ; MT_arms(Nmt)] with the MTOC appended as node Ne-1; the elastic arrays (bending
    # triples, per-fiber κ, fiber offsets, reshape, seg_rest) come from enet so MT bending/inextensibility ride
    # the same loop, while all CORTEX surface physics (turgor/membrane/hull/centroid) stay on net[:Nc].
    from aleph.laws.microtubule import merge_aster_into_cortex
    _khub = float(cortex.xl_k.max()) if cortex.xl_i.size else 1.0e4
    _m = merge_aster_into_cortex(net, microtubule, k_hub_pn_um=_khub)
    enet = _m["net"]; Ne = _m["Ne"]                   # Ne == Nc when MT off
    mtoc_pos = _m["mtoc_pos"]                          # (1,3) MTOC, or (0,3) when off
    hub_i = _m["hub_i"].astype(np.int64); hub_j = _m["hub_j"].astype(np.int64)
    hub_k = np.asarray(_m["hub_k"], np.float64); hub_rest = np.asarray(_m["hub_rest"], np.float64)
    tri = np.ascontiguousarray(enet.bend_triples, np.int32)
    alpha = np.ascontiguousarray(_per_triple_alpha(enet), np.float64)
    fiber_off = np.ascontiguousarray(enet.fiber_offsets, np.int32)
    seg_per = np.diff(fiber_off) - 1
    seg_off = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    seg = float(enet.seg_rest.mean()) if enet.seg_rest.size else 1.0
    R0 = cortex.R0_mean if cortex.R0_mean > 0.0 else cortex.R_um

    # ---- seed the nucleus bead cloud concentric with the resting cortex ----
    c0 = net.pos.mean(axis=0)
    _K_vol_nuc = float(getattr(nucleus, "K_vol_Pa", 0.0) or 0.0) if nucleus is not None else 0.0
    _nuc_incompressible = _K_vol_nuc > 0.0
    # FSI needs the nucleus centroid/area/hull even when the nucleus is not incompressible (fluid→nucleus
    # coupling, FSI_ALL_COMPARTMENT_COUPLING_2026-07-16). Gate those on (_nuc_incompressible OR _nuc_fsi).
    _nuc_fsi = (biot_fsi is not None) and (nucleus is not None)
    if nucleus is not None:
        rng = np.random.default_rng(nucleus_seed)
        nuc_pos = _seed_nucleus_cloud(c0, nucleus.R_nuc_um, nucleus.n_beads, rng)
        V_nuc = (4.0 / 3.0) * np.pi * nucleus.R_nuc_um ** 3      # incompressible nucleus displacement [µm³]
        V0_nuc = float(ConvexHull(nuc_pos).volume) if _nuc_incompressible else 0.0   # rest nucleus hull vol
    else:
        nuc_pos = np.zeros((0, 3), dtype=np.float64)
        V_nuc = 0.0
        V0_nuc = 0.0
    n_nuc = nuc_pos.shape[0]
    N = Ne + n_nuc                                     # cortex + MT_arms + MTOC + nucleus (== Nc+n_nuc when MT off)
    pos0 = np.concatenate([enet.pos, mtoc_pos, nuc_pos], axis=0)   # [cortex ; MT_arms ; MTOC ; nucleus]
    # ---- IF cage (additive; if_cage=None → bit-identical). Radial-spoke perinuclear cage = the nucleus↔cortex
    # load path (INTERNAL_DISPLACEMENT_DIAGNOSIS): LINEAR-first, all bonds Hookean → ride link_spring_kernel
    # (NO new kernel). IF beads are a TAIL block (NOT in enet — reshape would kill their extensibility). ----
    _if_built = None
    n_if = 0
    if if_cage is not None:
        if nucleus is None or n_nuc == 0:
            raise ValueError("if_cage requires a nucleus (the cage links the nucleus surface to the cortex)")
        from aleph.laws.intermediate_filaments import resolve_intermediate_filaments, build_if_cage
        _if_res = resolve_intermediate_filaments(
            n_fil=int(if_cage.get("n_fil", 400)),
            E_if_Pa=float(if_cage.get("E_if_Pa", 6.0e6)),
            l_seg_um=float(if_cage.get("l_seg_um", 0.5)),
            ratio_xl=float(if_cage.get("ratio_xl", 0.1)))
        _if_built = build_if_cage(c0, nucleus.R_nuc_um, R0, nuc_pos, net.pos, _if_res,
                                  nuc_offset=Ne, cortex_offset=0, if_offset=Ne + n_nuc,
                                  crosslink=bool(if_cage.get("crosslink", True)))
        n_if = _if_built.n_if
        pos0 = np.concatenate([pos0, _if_built.pos], axis=0)      # [cortex ; MT ; MTOC ; nucleus ; IF]
        N = N + n_if
    # ---- Membrane surface — Template 2 real lipid sheet (additive; membrane_surface=None → bit-identical).
    # Replaces the LUMPED 2γ/R tension (`membrane`) with a distinct triangulated bilayer carrying area tension
    # (buffered-plateau γ_mem) + **Helfrich bending** (the defining mechanic the lumped model lacked) + ERM
    # tethers (two-way, rupture→bleb EMERGENT). The sheet is a TAIL block (NOT in enet — reshape/bending are
    # for fibers). Global membrane indices are [mem_off : mem_off+n_mem]. MEMBRANE_HELFRICH_DESIGN_2026-07-16. ----
    _mem_built = None
    n_mem = 0
    mem_off = N
    if membrane_surface is not None:
        from aleph.laws.membrane_surface import (build_membrane_mesh, build_membrane_hinges,
            calibrate_kappa_tilde, resolve_membrane_full, erm_rupture_force, GAMMA_MCA_PN_UM,
            membrane_area_kernel, helfrich_bending_kernel, erm_tether_kernel, membrane_bending_energy)
        from scipy.spatial import cKDTree
        _ms = membrane_surface if isinstance(membrane_surface, dict) else {}
        _mem_sub = int(_ms.get("subdiv", 4))
        _mem_gap = float(_ms.get("gap_um", 0.15))                 # membrane sits just outside the resting cortex
        _mem_res = resolve_membrane_full()                       # ResolvedMembrane (γ_mem, K_A, κ, τ_lysis)
        # Sheet radius = mean cortex radius + a small gap so the ERM tethers are SHORT (the membrane rides on the
        # cortex — real ezrin ~20 nm, coarse-mesh ~gap). Enclosing to the OUTERMOST cortex node instead would make
        # the tethers ~1 µm (unphysiological, sheet floating off the cortex). The resting cortex is mildly non-
        # spherical, so a few equatorial nodes can sit at/just outside the sheet; the ERM pin the membrane to the
        # cortex, and full leak-proof CONTAINMENT (per-node nearest-membrane radius) is a deferred follow-up
        # (design §5 — the scalar-envelope containment kernel is broken under non-spherical load; not wired here).
        R_mem = R0 + _mem_gap
        _mm = build_membrane_mesh(R_mem, subdivisions=_mem_sub, centre=tuple(c0))
        _mem_hinges = build_membrane_hinges(_mm.faces)
        # κ̃ calibrated to the RESTING sheet so the discrete sphere reproduces the continuum 8πκ_m (the SN 2/√3
        # factor is wrong for the icosphere — grid-invariant, measured; MEMBRANE_HELFRICH_DESIGN §calibration).
        _kappa_tilde = calibrate_kappa_tilde(_mem_res.kappa, _mm.verts, _mem_hinges)
        n_mem = _mm.verts.shape[0]
        mem_off = N                                              # global offset BEFORE appending the sheet
        _mem_faces_g = np.ascontiguousarray((_mm.faces + mem_off), np.int32)      # face indices → global pos_d
        _mem_hinges_g = np.ascontiguousarray((_mem_hinges + mem_off).astype(np.int32))
        # ERM tethers: each membrane node ↔ nearest cortex node, rest = the formation length (force-free at build)
        _cn = cKDTree(net.pos).query(_mm.verts)[1].astype(np.int64)
        _erm_m = np.ascontiguousarray((mem_off + np.arange(n_mem)).astype(np.int32))
        _erm_c = np.ascontiguousarray(_cn.astype(np.int32))      # cortex nodes live at [0:Nc]
        _erm_rest = np.ascontiguousarray(np.linalg.norm(_mm.verts - net.pos[_cn], axis=1).astype(np.float64))
        _k_erm = float(_ms.get("k_erm_pN_um", 50.0))             # ⚠ CFL-convenience (not single-molecule ezrin) → PI
        _f_rupt = erm_rupture_force(_mem_res.kappa, _mem_res.gamma_mem, GAMMA_MCA_PN_UM)  # KB-3.B1.4 derived (~11 pN)
        _sigma_mem = float(_mem_res.gamma_mem)                   # buffered-plateau tension (K_A off default → const)
        _mem_edge = float(np.linalg.norm(_mm.verts[_mm.edges[:, 0]] - _mm.verts[_mm.edges[:, 1]], axis=1).mean())
        pos0 = np.concatenate([pos0, _mm.verts], axis=0)         # [... ; nucleus ; IF ; membrane]
        N = N + n_mem
        _mem_built = True
    # MT tip nodes (outermost bead of each MT arm, merged indexing) for the tip↔cortex contact
    _mt_tips = np.array([int(fiber_off[f + 1]) - 1 for f in range(fiber_off.shape[0] - 1)
                         if int(fiber_off[f]) >= Nc], dtype=np.int64) if microtubule is not None else np.zeros(0, np.int64)

    # V0 = RESTING compressible-cytoplasm volume (hull of the cortex minus the nucleus it displaces), so
    # V_cyto/V0 = 1 at strain 0 → ΔP = dP0 (same convention as simulate_compressed_shell_on_device).
    _cortex_faces = _oriented_hull_faces(net.pos)           # FIXED outward mesh → GPU divergence volume/area
    V0 = float(ConvexHull(net.pos).volume) - V_nuc          # (divergence on _cortex_faces ≡ ConvexHull, validated)
    vmin = VMIN_FRAC * V0
    half_gap = R0 * (1.0 - strain)
    K_vol = TURGOR_PI_IN0 / (1.0 - VMIN_FRAC)
    A0_mem = float(ConvexHull(net.pos).area)                # resting cortex area → membrane baseline = γ_mem
    _nuc_faces = _oriented_hull_faces(nuc_pos) if ((_nuc_incompressible or _nuc_fsi) and n_nuc > 0) else np.zeros((0, 3), np.int32)
    # ---- Piece #2: biphasic poroelastic cytoplasm (additive; Lp/K_drained/load_time all None → bit-identical) ----
    # (A) DRAINAGE (Kedem-Katchalsky, σ≈1): the osmotic reference (water) volume V0_eff drains toward the geometric
    #     cytoplasm V_cyto so the van't Hoff turgor relaxes to dP0. Relaxation is the MEMBRANE clock
    #     τ_osm=(V0_eff−vmin)/(Lp·A·Π_in) (~30–200 s), NOT the poroelastic τ_p(~1 s); integrated by analytic
    #     exp-relaxation over the physical load time (unconditionally stable at any Lp). Lp default 1.6e-8 µm/(s·Pa)
    #     (COS-7 exosmotic/efflux, PMC3161049); AFM at v≫R/τ_osm stays ~undrained (correct fast-ramp behaviour).
    # (B) DRAINED SOLID (Terzaghi effective stress): dP_solid=K_drained·max(0,(V0−V_cyto)/V0) on the FIXED solid rest
    #     volume V0 (NOT the draining V0_eff — else it vanishes at the drained limit where it must carry the load).
    #     K_drained≈300 Pa (Moeendarbary 2013 soft-epithelial branch; ν≈0.25–0.3). Full anchors: FF_RESULTS_LOG.
    V0_eff = V0
    _K_drained = float(K_drained_Pa) if K_drained_Pa is not None else 0.0
    dP_osm = dP_solid = 0.0
    tau_osm = float("inf")
    _KBT_PN_UM = 4.142e-3                              # kB·T at 300 K [pN·µm]
    _xbeta_kT = float(xl_x_beta_um) / _KBT_PN_UM
    if k_plate is None:
        k_plate = 10.0 * (K_vol / V0) * (4.0 * np.pi * R0**2) ** 2 / Nc
    if dt_mu <= 0.0:
        kmax = float(net.kappa.max()) / seg**3
        if cortex.xl_i.size:
            kmax = max(kmax, float(cortex.xl_k.max()))
        kmax = max(kmax, (K_vol / V0) * (4.0 * np.pi * R0**2) ** 2 / Nc, k_plate)
        if nucleus is not None:
            kmax = max(kmax, nucleus.k_chrom + nucleus.k_lamin)     # nucleus stiffness sets the CFL too
            if _nuc_incompressible and V0_nuc > 0.0:               # nucleoplasm volume-pressure (same form as cortex turgor)
                kmax = max(kmax, (_K_vol_nuc / V0_nuc) * (4.0 * np.pi * nucleus.R_nuc_um ** 2) ** 2 / max(n_nuc, 1))
        if _if_built is not None and _if_built.bond_k.size:        # IF backbone/LINC/anchor stiffness sets the CFL
            kmax = max(kmax, float(_if_built.bond_k.max()))
        # LUMPED membrane (buffered-plateau γ_mem) is a soft ~few-Pa inward tension → no CFL term; the K_A
        # elastic upturn is deferred (reservoir PI-blocked). The REAL sheet (membrane_surface) adds the ERM
        # tether spring (k_erm) + Helfrich bending (~κ̃/ℓ³, soft at subdiv-4 per design §3) to the CFL.
        if _mem_built:
            kmax = max(kmax, _k_erm, _kappa_tilde / max(_mem_edge, 1e-6) ** 3)
        dt_mu = 0.1 / kmax

    # ---- Physical press-speed ramp (η=65.9 Pa·s cytoplasm-viscosity-limited deformation) ----
    # Each numerical step (kernel step = dt_mu, CFL-stable) equals dt_real = dt_mu/μ SECONDS, with μ the NF2007
    # cytoplasm mobility (μ=log(L/δ)/(3πηL)). The plate then advances v_press·dt_real per step, so relaxation
    # completeness is set by PHYSICS (press-speed vs the η-limited deformation speed), not an arbitrary n_steps —
    # the fix for the convergence artifact (a fast/instant press reads a NON-equilibrium transient; PI 2026-07-08).
    # v_press_um_s=None → instant strain applied at step 0 (legacy quasi-static; converge via n_steps).
    _ramp_on = v_press_um_s is not None
    _dt_real = _N_ramp = 0.0
    if _ramp_on:
        if eta_bulk_Pa_s is not None:
            # BULK-cytoplasm-calibrated per-node drag: the whole-cell Stokes drag 6π·η·R distributed over the Nc
            # cortex nodes (accounts for the hydrodynamic screening of the dense cortex that the single-fiber
            # NF2007 mobility ignores). Grid-intensive: Σγ = 6π·η·R, independent of Nc. η=65.9 Pa·s (Dessard 2024).
            # Validated against the analytic Newtonian ground truth σ_visc = η·ε̇ (ff_drag_calibration).
            gamma_node = 6.0 * np.pi * float(eta_bulk_Pa_s) * R0 / max(Nc, 1)   # pN·s/µm
            _mu_phys = 1.0 / gamma_node
        else:
            from aleph.laws.units import fiber_mobility
            _mu_phys = fiber_mobility(max(seg, 0.05))                  # µm/(pN·s) at η=65.9 (NF2007 §5.2)
        _dt_real = dt_mu / _mu_phys                                # physical seconds per numerical step
        _t_ramp = (R0 * strain) / float(v_press_um_s)              # time to indent one side by R0·strain
        _N_ramp = max(1, int(np.ceil(_t_ramp / _dt_real)))
        n_steps = _N_ramp + max(0, int(dwell_steps))
        load_time_s = n_steps * _dt_real                          # total real time → drives drainage + turnover

    # ---- timing for drainage (A) + crosslink turnover, using the (possibly ramp-updated) n_steps / load_time_s ----
    _drain_on = (Lp_um_s_Pa is not None and load_time_s is not None and pressure_setpoint is None)
    _n_refresh = max(1, n_steps // max(turgor_every, 1))
    _dt_refresh = (float(load_time_s) / _n_refresh) if _drain_on else 0.0
    _turnover_on = (xl_koff_per_s is not None and load_time_s is not None)
    _n_turn = max(1, n_steps // max(reshape_every, 1))
    _koff0_dt = (float(xl_koff_per_s) * float(load_time_s) / _n_turn) if _turnover_on else 0.0

    d = device
    pos_d = wp.array(np.ascontiguousarray(pos0, np.float64), dtype=wp.vec3d, device=d)
    f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    react_d = wp.zeros(3, dtype=wp.float64, device=d)          # [0]=top force, [1]=cortex contacts, [2]=unused
    # GPU turgor-geometry: enclosed volume/area/centroid via divergence over an OUTWARD-oriented mesh (GPU
    # kernels, no per-refresh full pos→CPU). The mesh is RE-HULLED every `_rehull_every` refreshes (not every
    # step) so it tracks the deforming surface — a fixed rest mesh drifts (V error compounds via turgor as the
    # shape leaves the rest sphere). This keeps ConvexHull cost ~`_rehull_every`× lower while staying accurate.
    cortex_faces_d = wp.array(np.ascontiguousarray(_cortex_faces, np.int32), dtype=wp.int32, device=d)
    nuc_faces_d = (wp.array(np.ascontiguousarray(_nuc_faces, np.int32), dtype=wp.int32, device=d)
                   if _nuc_faces.shape[0] > 0 else None)
    _cen_d = wp.zeros(3, dtype=wp.float64, device=d)           # centroid reduction accumulator
    _va_d = wp.zeros(2, dtype=wp.float64, device=d)            # [6V, 2A] accumulator (reused cortex+nucleus)
    _rehull_every = 5                                          # refreshes between mesh re-hulls (5·turgor_every steps)
    _refresh_i = [0]                                           # refresh counter (mutable, for the re-hull cadence)
    tri_d = wp.array(tri, dtype=wp.int32, device=d)
    alpha_d = wp.array(alpha, dtype=wp.float64, device=d)
    foff_d = wp.array(fiber_off, dtype=wp.int32, device=d)
    soff_d = wp.array(seg_off, dtype=wp.int32, device=d)
    srest_d = wp.array(np.ascontiguousarray(enet.seg_rest, np.float64), dtype=wp.float64, device=d)  # merged (incl MT arms)
    # crosslinks = cortex α-actinin/filamin + MTOC↔arm-base hub (both enter the same link_spring)
    xl_i_all = np.concatenate([cortex.xl_i, hub_i]).astype(np.int32)
    xl_j_all = np.concatenate([cortex.xl_j, hub_j]).astype(np.int32)
    xl_k_all = np.concatenate([cortex.xl_k, hub_k]).astype(np.float64)
    xl_r_all = np.concatenate([cortex.xl_rest, hub_rest]).astype(np.float64)
    if _if_built is not None:      # IF cage backbone + LINC + cortex-anchor + crosslinks (all Hookean → link_spring)
        xl_i_all = np.concatenate([xl_i_all, _if_built.bond_i]).astype(np.int32)
        xl_j_all = np.concatenate([xl_j_all, _if_built.bond_j]).astype(np.int32)
        xl_k_all = np.concatenate([xl_k_all, _if_built.bond_k]).astype(np.float64)
        xl_r_all = np.concatenate([xl_r_all, _if_built.bond_rest]).astype(np.float64)
    n_xl = int(xl_i_all.size)
    has_xl = n_xl > 0
    if has_xl:
        xl_d = wp.array(np.ascontiguousarray(np.stack([xl_i_all, xl_j_all], 1), np.int32), dtype=wp.int32, device=d)
        kxl_d = wp.array(np.ascontiguousarray(xl_k_all), dtype=wp.float64, device=d)
        r0_d = wp.array(np.ascontiguousarray(xl_r_all), dtype=wp.float64, device=d)
    # MT tip↔cortex soft contact: each MT tip ↔ nearest cortex node (excluded-volume prop under confinement)
    mtp_d = None
    if _mt_tips.size:
        from scipy.spatial import cKDTree
        _tip_cortex = cKDTree(net.pos).query(enet.pos[_mt_tips])[1].astype(np.int64)
        mtp_d = wp.array(np.ascontiguousarray(np.stack([_mt_tips, _tip_cortex], 1).astype(np.int32)),
                         dtype=wp.int32, ndim=2, device=d)
        k_mt_c = float(k_mt_contact) if k_mt_contact is not None else _khub
        r_mt_c = 0.5                                  # contact shell [µm]: engages when cortex reaches within 0.5µm of a tip
    has_myo = bool(cortex.myo_i.size) and f_myo != 0.0
    if has_myo:
        myo_d = wp.array(np.ascontiguousarray(np.stack([cortex.myo_i, cortex.myo_j], 1), np.int32),
                         dtype=wp.int32, device=d)
    # ---- Membrane-surface device arrays + force-launch closure (real sheet: area tension + Helfrich bending +
    # ERM tether). n_mem==0 → the closure is never invoked (lumped-tension path unchanged, bit-identical). ----
    if _mem_built:
        mem_faces_d = wp.array(_mem_faces_g, dtype=wp.int32, ndim=2, device=d)
        mem_hinges_d = wp.array(_mem_hinges_g, dtype=wp.int32, ndim=2, device=d)
        erm_m_d = wp.array(_erm_m, dtype=wp.int32, device=d)
        erm_c_d = wp.array(_erm_c, dtype=wp.int32, device=d)
        erm_bound_d = wp.array(np.ones(n_mem, np.int32), dtype=wp.int32, device=d)
        erm_rest_d = wp.array(_erm_rest, dtype=wp.float64, device=d)

    def _membrane_forces():
        """Accumulate the real-membrane forces into f_d at the sheet's global indices: in-plane area tension
        (buffered-plateau γ_mem) + Helfrich bending (calibrated κ̃) + ERM tether (two-way; tension>f_rupt →
        UNBIND = a bleb, emergent). All operate on the GLOBAL pos_d/f_d so they compose additively."""
        wp.launch(membrane_area_kernel, dim=mem_faces_d.shape[0],
                  inputs=[pos_d, mem_faces_d, wp.float64(_sigma_mem), f_d], device=d)
        wp.launch(helfrich_bending_kernel, dim=mem_hinges_d.shape[0],
                  inputs=[pos_d, mem_hinges_d, wp.float64(_kappa_tilde), f_d], device=d)
        wp.launch(erm_tether_kernel, dim=n_mem,
                  inputs=[pos_d, pos_d, erm_m_d, erm_c_d, erm_bound_d, wp.float64(_k_erm), erm_rest_d,
                          wp.float64(_f_rupt), f_d, f_d], device=d)

    def _mem_metrics(pos_all):
        """Membrane-surface diagnostics for the return metrics ({} when no sheet → metrics bit-identical to the
        no-membrane run). Reports the sheet area/areal-strain, the bending energy (vs the 8πκ sphere reference),
        and the ERM bleb count (ruptured tethers = emergent local detachments)."""
        if not _mem_built:
            return {}
        mp = pos_all[mem_off:mem_off + n_mem]
        fc = _mm.faces
        area_mem = float(0.5 * np.linalg.norm(np.cross(mp[fc[:, 1]] - mp[fc[:, 0]],
                                                       mp[fc[:, 2]] - mp[fc[:, 0]]), axis=1).sum())
        E_bend = membrane_bending_energy(mp, _mem_hinges, _kappa_tilde)
        bound = int(erm_bound_d.numpy().sum())
        return {"membrane_surface_on": True, "n_mem": int(n_mem), "mem_subdiv": int(_mem_sub),
                "R_mem_um": float(R_mem), "mem_area_um2": area_mem, "mem_A0_um2": float(_mm.A0),
                "mem_areal_strain": float((area_mem - _mm.A0) / _mm.A0),
                "mem_bend_energy_pN_um": float(E_bend),
                "mem_bend_over_8pikappa": float(E_bend / (8.0 * np.pi * _mem_res.kappa)),
                "kappa_tilde_pN_um": float(_kappa_tilde), "gamma_mem_sheet_pN_um": float(_sigma_mem),
                "erm_bound": bound, "erm_ruptured_bleb": int(n_mem - bound),
                "erm_f_rupt_pN": float(_f_rupt), "k_erm_pN_um": float(_k_erm)}

    nT = tri.shape[0]
    centre = wp.vec3d(0.0, 0.0, 0.0)
    cz = 0.0
    dP_area = 0.0
    dP_mem_area = 0.0
    gamma_mem_tot = 0.0
    dP_mem = 0.0
    # seed the nucleus centroid to its BUILD centroid (not the world origin) so the first FSI source/interp
    # before the first _refresh_turgor references the right body; area_nuc seeded to the sphere area.
    centre_nuc = wp.vec3d(float(c0[0]), float(c0[1]), float(c0[2])) if n_nuc else wp.vec3d(0.0, 0.0, 0.0)
    dP_nuc_area = 0.0                                # nucleoplasm ΔP·(area/bead)
    area_nuc = float(4.0 * np.pi * nucleus.R_nuc_um ** 2) if nucleus is not None else 0.0   # live nucleus surface area (FSI back-force)
    V_nuc_cur = V_nuc                                # current nucleus hull volume (tracked when incompressible)

    def _refresh_turgor():
        nonlocal centre, cz, dP_area, dP_mem_area, gamma_mem_tot, dP_mem, V0_eff, dP_osm, dP_solid, tau_osm
        nonlocal area_nuc
        nonlocal centre_nuc, dP_nuc_area, V_nuc_cur, cortex_faces_d, nuc_faces_d
        # periodic re-hull: refresh the mesh to the current (deformed) surface so the divergence volume stays
        # accurate. Only here do we touch the CPU (ConvexHull) + a full pos copy — every _rehull_every refreshes.
        if _refresh_i[0] % _rehull_every == 0:
            _pr = pos_d.numpy()
            cortex_faces_d = wp.array(np.ascontiguousarray(_oriented_hull_faces(_pr[:Nc]), np.int32),
                                      dtype=wp.int32, device=d)
            if (_nuc_incompressible or _nuc_fsi) and n_nuc > 0:
                cf = _oriented_hull_faces(_pr[Ne:Ne + n_nuc])   # nucleus-local indices
                nuc_faces_d = wp.array(np.ascontiguousarray(cf, np.int32), dtype=wp.int32, device=d)
        _refresh_i[0] += 1
        # GPU-resident cortex centroid + enclosed volume/area (divergence over the current oriented mesh).
        _cen_d.zero_()
        wp.launch(_centroid_sum_kernel, dim=Nc, inputs=[pos_d, wp.int32(0), _cen_d], device=d)
        _va_d.zero_()
        wp.launch(_mesh_vol_area_kernel, dim=cortex_faces_d.shape[0], inputs=[pos_d, cortex_faces_d, wp.int32(0), _va_d], device=d)
        wp.synchronize_device(d)
        cc = _cen_d.numpy() / float(Nc)
        centre = wp.vec3d(float(cc[0]), float(cc[1]), float(cc[2]))
        cz = float(cc[2])
        _va = _va_d.numpy()
        V = abs(float(_va[0])) / 6.0
        area = float(_va[1]) / 2.0
        R_mean = (3.0 * V / (4.0 * np.pi)) ** (1.0 / 3.0)     # equivalent-sphere radius (membrane Laplace 2γ/R)
        # nucleus volume conservation (incompressibility): ΔP_nuc = K_vol·(V0_nuc−V_nuc_cur)/V0_nuc (two-sided
        # bulk modulus), radially outward about the nucleus centroid so a plate-flattened nucleus bulges. GPU too.
        # centroid + area needed by BOTH the incompressibility law AND the FSI fluid→nucleus force → gate on OR.
        if (_nuc_incompressible or _nuc_fsi) and n_nuc > 0 and nuc_faces_d is not None:
            _cen_d.zero_()
            wp.launch(_centroid_sum_kernel, dim=n_nuc, inputs=[pos_d, wp.int32(Ne), _cen_d], device=d)
            _va_d.zero_()
            wp.launch(_mesh_vol_area_kernel, dim=nuc_faces_d.shape[0], inputs=[pos_d, nuc_faces_d, wp.int32(Ne), _va_d], device=d)
            wp.synchronize_device(d)
            cn = _cen_d.numpy() / float(n_nuc)
            centre_nuc = wp.vec3d(float(cn[0]), float(cn[1]), float(cn[2]))
            _van = _va_d.numpy()
            V_nuc_cur = abs(float(_van[0])) / 6.0
            area_nuc = float(_van[1]) / 2.0
            # nucleoplasm incompressibility pressure — ONLY when incompressible (V0_nuc>0; else 0/0 divide).
            if _nuc_incompressible and V0_nuc > 0.0:
                dP_nuc = _K_vol_nuc * (V0_nuc - V_nuc_cur) / V0_nuc
                dP_nuc_area = dP_nuc * area_nuc / n_nuc
        V_cyto = V - V_nuc                            # incompressible-nucleus displacement coupling
        # (B) drained solid skeleton (Terzaghi effective stress), referenced to the FIXED rest volume V0 — always
        #     present (total = pore pressure + effective stress); K_drained=0 → bit-identical to the pure-turgor prior.
        dP_solid = _K_drained * max(0.0, (V0 - V_cyto) / V0)
        if pressure_setpoint is not None:
            dP_osm = float(pressure_setpoint)          # perfect-water-flux (drained) pore-pressure limit
        else:
            # (A) van't Hoff osmotic turgor on the (draining) water reference V0_eff
            dP_osm = max(TURGOR_PI_IN0 * (V0_eff - vmin) / max(V_cyto - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0), 0.0)
            # drain the osmotic reference toward V_cyto over the physical step time (analytic exp-relaxation)
            if _drain_on:
                tau_osm = max(V0_eff - vmin, 1e-6 * V0) / (Lp_um_s_Pa * area * TURGOR_PI_IN0)
                V0_eff += (V_cyto - V0_eff) * (1.0 - np.exp(-_dt_refresh / tau_osm))
                V0_eff = min(max(V0_eff, vmin + 1e-6 * V0), V0)
        dP = dP_osm + dP_solid
        dP_area = dP * area / Nc
        if membrane is not None:
            # Piece #4: RESERVOIR PLATEAU + K_A UPTURN. The area reservoir (folds/microvilli/caveolae, Raucher &
            # Sheetz 1999) unfolds at ~constant baseline tension γ_mem until the APPARENT area exceeds
            # A0·(1+f_excess); past that the bilayer itself stretches → steep K_A elastic upturn (Rawicz 2000),
            # capped at the lysis tension. f_excess is a controlled variable (0.25 default; no MCF7 datum — see
            # FF_RESULTS_LOG PI-flag). f_excess=0 recovers the pure buffered plateau (bit-identical to prior).
            f_exc = float(getattr(membrane, "f_excess", 0.0) or 0.0)
            A_thresh = A0_mem * (1.0 + f_exc)
            if f_exc > 0.0 and area > A_thresh:
                eps_bil = (area - A_thresh) / A_thresh
                gamma_mem_tot = min(membrane.gamma_mem + membrane.K_A * eps_bil, membrane.tau_lysis)
            else:
                gamma_mem_tot = min(membrane.gamma_mem, membrane.tau_lysis)   # buffered plateau = γ_mem
            dP_mem = 2.0 * gamma_mem_tot / max(R_mean, 1e-9)
            dP_mem_area = -dP_mem * area / Nc          # inward (negative → turgor_kernel pushes toward centre)
        return V, V_cyto, area, dP

    # S2 localized patch load (default None → unused; plate path unchanged). Force f_node·axis on a cap.
    _ll = local_load is not None
    if _ll:
        _lax = np.asarray(local_load.get("axis", (0.0, 0.0, 1.0)), np.float64)
        _lax = _lax / (np.linalg.norm(_lax) + 1e-30)
        _ll_axis = wp.vec3d(float(_lax[0]), float(_lax[1]), float(_lax[2]))
        _ll_cos = wp.float64(float(local_load.get("cos_thresh", 0.9)))
        _ll_f = wp.float64(float(local_load.get("f_node_pN", 0.0)))
        _pc0 = pos0[:Nc] - pos0[:Nc].mean(0)
        _pmask = (_pc0 @ _lax) / (np.linalg.norm(_pc0, axis=1) + 1e-30) >= float(local_load.get("cos_thresh", 0.9))
        # rest patch position measured relative to the cortex CENTROID (_pc0 is already centred). The patch
        # load is an unbalanced net force → the free cell drifts rigidly; centroid-referencing every frame
        # removes that drift so patch_deflection is the true local indentation, not indentation + COM drift.
        _patch_rest = float((_pc0[_pmask] @ _lax).mean()) if _pmask.any() else 0.0

    # ================= FSI-ON path (biot_fsi set) — spatial Biot pore-pressure field, two-way ==================
    # ADDITIVE default-OFF: biot_fsi=None → this block is skipped and the monolithic loop below runs UNCHANGED
    # (bit-identical = the 28-part + native-cell regression gate). Architecture: FSI_WIRING_DESIGN_2026-07-16.md
    # — OUTER physical-time loop advances the Biot pore-pressure field on τ_p (Darcy redistribution); INNER
    # mechanical relaxation brings the skeleton to equilibrium each physical step (γ_solid drag ⇒ quasi-static,
    # so the 6πηR whole-cell drag is RETIRED as the physical clock — the poroelastic apparent viscosity now
    # EMERGES from the Biot D + modulus M via the spatially-varying pore-pressure back-reaction). The cortex is
    # a SHELL, so the back-reaction is a local outward normal p_i·area·r̂ (spatially-varying turgor), and the
    # field carries the ZERO-MEAN spatial deviation → the uniform turgor stays in its own channel (no η/turgor
    # double-count). See units.fiber_point_drag_solvent for the retired-6πηR γ_solid statement.
    if biot_fsi is not None:
        from aleph.laws.biot_fluid_warp import BiotField, ibm_spread_kernel, ibm_interp_kernel, biot_gradient_kernel
        if v_press_um_s is None:
            raise ValueError("biot_fsi requires v_press_um_s (the physical press clock; 6πηR is retired)")
        D_fsi = float(biot_fsi.get("D_um2_s", 50.0))          # Moeendarbary poroelastic diffusivity D (KB-3.B3.2)
        dx_fsi = float(biot_fsi.get("dx_um", 0.6))            # coarse fluid grid (P7.F.7 multiscale, dx ≪ R)
        # Biot coupling modulus M — the pore pressure generated per unit volumetric strain (p ≈ M·α·ε,
        # undrained). This is NOT K_drained (the skeleton stiffness already carried by dP_solid); the earlier
        # default conflated them. Physically M ≈ K_undrained − K_drained: Moeendarbary 2013 drained modulus
        # 0.4–0.9 kPa (HeLa 900, HT1080/MDCK 400 Pa, KB-3.B3.2) × the cell undrained/drained ratio ~2–3 →
        # K_u ~ 1–2.7 kPa → M ~ 0.6–1.8 kPa. Default 1 kPa (KB-DRAFT Biot modulus, FSI_ALL_COMPARTMENT_COUPLING).
        # Do NOT tune to a nucleus-response band (no-param-tuning). NB even at physical M the fluid-mediated
        # interior deformation is small (the angular pore-pressure deviation is ~Pa; the dominant nucleus
        # coupling is DIRECT — IF cage / LINC, not the fluid).
        M_fsi = float(biot_fsi.get("M_biot_Pa", 1000.0))     # Pa — Biot coupling modulus (undrained-derived)
        inner = int(biot_fsi.get("inner_steps", 40))         # inner mechanical descents per physical step
        n_phys = int(biot_fsi.get("n_phys", 300))            # outer physical steps over the ramp
        p_clip = float(biot_fsi.get("p_clip_Pa", 5.0e3))     # numerical guard on |excess p| (Pa)
        t_ramp_fsi = (R0 * strain) / float(v_press_um_s)     # physical ramp time (drag-independent)
        dt_phys = t_ramp_fsi / max(n_phys, 1)                # physical seconds per outer step (the retired clock)
        # coarse fluid grid centred on the cell (P7.F.7: dx ≪ R, grid spans the cell)
        c_cell = pos0[:Nc].mean(0)
        half_fsi = R0 * 1.6
        n_grid = int(2.0 * half_fsi / dx_fsi) + 1
        origin_fsi = wp.vec3d(float(c_cell[0] - half_fsi), float(c_cell[1] - half_fsi), float(c_cell[2] - half_fsi))
        bf = BiotField(n_grid, dx_fsi, D_fsi, device=d)
        # interior mask (grid cells within R0 of the cell centre → zero-mean is taken over the cell, not the box)
        _gx = (np.arange(n_grid) + 0.5) * dx_fsi
        _GX, _GY, _GZ = np.meshgrid(_gx, _gx, _gx, indexing="ij")
        _int_mask = ((_GX - half_fsi) ** 2 + (_GY - half_fsi) ** 2 + (_GZ - half_fsi) ** 2) <= (R0 * 1.05) ** 2
        # Biot sub-cycle count to respect the 3-D CFL dt ≤ dx²/6D over one physical step
        _biot_sub = max(1, int(np.ceil(dt_phys / bf.dt_max)))
        _biot_dt = dt_phys / _biot_sub
        src_d = wp.zeros(Nc, dtype=wp.vec3d, device=d)                 # per-node (Δp, weight) deposit
        gridsrc_d = wp.zeros((n_grid, n_grid, n_grid), dtype=wp.vec3d, device=d)
        pgrid_d = wp.zeros((n_grid, n_grid, n_grid), dtype=wp.vec3d, device=d)   # p packed in comp0 for interp
        pnode_d = wp.zeros(N, dtype=wp.vec3d, device=d)               # per-node interpolated excess p (comp0)
        pos_prev_d = wp.clone(pos_d)
        p_np = np.zeros((n_grid, n_grid, n_grid), np.float64)         # host-side excess pore-pressure field
        final_gap = R0 * (1.0 - strain)
        _fsi_trace = []
        _p_max_hist = []
        area = float(ConvexHull(net.pos).area)                       # seed area (updated each step by refresh)
        # ---- FSI coupling to the INTERIOR (nucleus 2-D shell + MT 1-D fiber) — FSI_ALL_COMPARTMENT_COUPLING_2026-07-16 ----
        _dt_refresh = dt_phys                                         # drainage integrates on the PHYSICAL clock (NOT _drain_on=False)
        _fsi_nuc = n_nuc > 0                                          # fluid → nucleus (the real interior-response channel)
        _rr = np.sqrt((_GX - half_fsi) ** 2 + (_GY - half_fsi) ** 2 + (_GZ - half_fsi) ** 2)   # grid |x−cell_centre|
        # MT (1-D fiber): −V_seg·∇p buoyancy — NEGLIGIBLE (~1e-3 pN) but wired for dimensional completeness (honest;
        # the MT's real fluid response needs the u-p velocity drag = PI-authored KU; its dominant coupling is the strut).
        _n_fib = Ne - Nc
        _mt_rad = float(getattr(microtubule, "radius_um", 0.0125)) if microtubule is not None else 0.0125
        gradp_d = wp.zeros((n_grid, n_grid, n_grid), dtype=wp.vec3d, device=d) if _n_fib > 0 else None
        gnode_d = wp.zeros(N, dtype=wp.vec3d, device=d) if _n_fib > 0 else None
        vseg_d = None
        if _n_fib > 0:
            _vseg = np.full(_n_fib, np.pi * (_mt_rad ** 2) * float(seg), np.float64)   # π·a²·ℓ displaced tube vol [µm³]
            _vseg[-1] = 0.0                                           # MTOC (last node in [Nc:Ne]) is the hub, no tube
            vseg_d = wp.array(np.ascontiguousarray(_vseg), dtype=wp.float64, device=d)
        _fsi_net_hist = []                                            # net fluid-coupling force magnitude (momentum monitor)
        for pstep in range(n_phys):
            _Vc, _Vcyto, area, _dPc = _refresh_turgor()              # centre, dP_area (uniform turgor), area
            # --- source: deposit the compression increment onto the field (mass-matrix-normalized IBM spread) ---
            wp.launch(fsi_compression_source_kernel, dim=Nc,
                      inputs=[pos_d, pos_prev_d, centre, wp.float64(M_fsi / max(R0, 1e-9)), src_d], device=d)
            gridsrc_d.zero_()
            wp.launch(ibm_spread_kernel, dim=Nc,
                      inputs=[pos_d, src_d, gridsrc_d, origin_fsi, wp.float64(dx_fsi), wp.int32(n_grid)], device=d)
            wp.synchronize_device(d)
            gs = gridsrc_d.numpy()                                    # (n,n,n,3): [...,0]=Σδ·Δp, [...,1]=Σδ
            num = gs[..., 0]; den = gs[..., 1]
            deposit = np.where(den > 1e-9, num / np.maximum(den, 1e-9), 0.0)
            p_np = p_np + deposit
            # --- drainage/redistribution: diffuse the field on D over one physical step (τ_p ~ L²/D) ---
            bf.set_field(np.ascontiguousarray(p_np))
            for _ in range(_biot_sub):
                bf.step(_biot_dt)
            p_np = bf.get_field()
            # zero-mean over the CYTOSOL (the MEAN pore pressure is the turgor, kept in its own channel; the Biot
            # field carries only the spatial DEVIATION → no double-count). When a nucleus is present, EXCLUDE its
            # interior (a HOLE in the cytoplasm whose field is unphysical + would bias the deviation the nucleus
            # surface feels), re-centred on the live nucleus centroid. n_nuc=0 → the cell ball, cortex-FSI bit-identical.
            _mask = _int_mask
            if _fsi_nuc:
                _cnx = float(centre_nuc[0]) - float(origin_fsi[0]); _cny = float(centre_nuc[1]) - float(origin_fsi[1])
                _cnz = float(centre_nuc[2]) - float(origin_fsi[2])
                _solid = ((_GX - _cnx) ** 2 + (_GY - _cny) ** 2 + (_GZ - _cnz) ** 2) <= (nucleus.R_nuc_um ** 2)
                _mask = _int_mask & ~_solid
            if _mask.any():
                p_np = p_np - float(p_np[_mask].mean())
            np.clip(p_np, -p_clip, p_clip, out=p_np)
            _p_max_hist.append(float(np.abs(p_np[_mask]).max()) if _mask.any() else 0.0)
            # --- back-reaction: interp p → per-node excess pressure over ALL nodes (dim=N) so the cortex, nucleus
            #     and MT each read their own local p; each back-reaction kernel indexes its slice (IBM kernels unchanged). ---
            pgrid_np = np.zeros((n_grid, n_grid, n_grid, 3), np.float64); pgrid_np[..., 0] = p_np
            pgrid_d = wp.array(np.ascontiguousarray(pgrid_np), dtype=wp.vec3d, device=d)
            pnode_d.zero_()
            wp.launch(ibm_interp_kernel, dim=N,
                      inputs=[pos_d, pgrid_d, pnode_d, origin_fsi, wp.float64(dx_fsi), wp.int32(n_grid)], device=d)
            area_per_node = area / max(Nc, 1)                         # A_cortex/Nc from the last refresh [µm²]
            area_per_bead_nuc = (area_nuc / max(n_nuc, 1)) if _fsi_nuc else 0.0   # A_nuc/n_nuc for the fluid→nucleus force
            # MT −∇p: compute ∇p on the grid + interp to every node (MT kernel indexes [Nc:Ne]); negligible force.
            if _n_fib > 0:
                wp.launch(biot_gradient_kernel, dim=(n_grid, n_grid, n_grid),
                          inputs=[pgrid_d, gradp_d, wp.float64(1.0 / (2.0 * dx_fsi)), wp.int32(n_grid)], device=d)
                gnode_d.zero_()
                wp.launch(ibm_interp_kernel, dim=N,
                          inputs=[pos_d, gradp_d, gnode_d, origin_fsi, wp.float64(dx_fsi), wp.int32(n_grid)], device=d)
            wp.copy(pos_prev_d, pos_d)                                # snapshot for next step's compression source
            # plate advances v_press·dt_phys per physical step (physical ramp)
            hg = max(final_gap, R0 - float(v_press_um_s) * (pstep + 1) * dt_phys)
            # --- INNER mechanical relaxation to equilibrium (quasi-static; γ_solid ⇒ effectively instantaneous) ---
            for _in in range(inner):
                wp.launch(_zero, dim=N, inputs=[f_d], device=d)
                if nT:
                    wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
                if has_xl:
                    wp.launch(link_spring_kernel, dim=n_xl, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
                if mtp_d is not None:
                    wp.launch(soft_contact_kernel, dim=mtp_d.shape[0], inputs=[pos_d, mtp_d, wp.float64(r_mt_c), wp.float64(k_mt_c), f_d], device=d)
                if has_myo:
                    wp.launch(myosin_kernel, dim=cortex.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
                wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_area), f_d], device=d)  # uniform turgor
                wp.launch(fsi_pressure_force_kernel, dim=Nc,
                          inputs=[pos_d, centre, pnode_d, wp.float64(area_per_node), f_d], device=d)  # spatial Biot deviation (cortex, OUTWARD)
                if _n_fib > 0:                                        # MT 1-D fiber: −V_seg·∇p buoyancy (negligible, honest)
                    wp.launch(fsi_fiber_gradp_force_kernel, dim=_n_fib,
                              inputs=[wp.int32(Nc), gnode_d, vseg_d, f_d], device=d)
                if membrane is not None and n_mem == 0:       # lumped 2γ/R only when no real sheet (no double-count)
                    wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_mem_area), f_d], device=d)
                if n_mem:                                     # real membrane sheet REPLACES the lumped tension
                    _membrane_forces()
                if n_nuc:
                    wp.launch(nucleus_shell_kernel, dim=n_nuc,
                              inputs=[pos_d, wp.int32(Ne), centre, wp.float64(nucleus.R_nuc_um),
                                      wp.float64(nucleus.k_chrom), wp.float64(nucleus.k_lamin),
                                      wp.float64(nucleus.d_knee_um), wp.float64(nucleus.F_knee_pN), f_d], device=d)
                    if _nuc_incompressible:
                        wp.launch(nucleus_pressure_kernel, dim=n_nuc,
                                  inputs=[pos_d, wp.int32(Ne), centre_nuc, wp.float64(dP_nuc_area), f_d], device=d)
                    # fluid → nucleus: local cytoplasm pore pressure squeezes the nucleus shell INWARD (the real
                    # interior-response channel; disjoint from shell-elasticity + nucleoplasm-incompressibility above).
                    wp.launch(fsi_nucleus_pressure_force_kernel, dim=n_nuc,
                              inputs=[pos_d, wp.int32(Ne), centre_nuc, pnode_d, wp.float64(area_per_bead_nuc), f_d], device=d)
                if rigid_plate:
                    react_d.zero_()
                    wp.launch(rigid_plate_step_kernel, dim=N,
                              inputs=[pos_d, wp.float64(dt_mu), f_d, wp.float64(cz), wp.float64(hg), react_d], device=d)
                else:
                    react_d.zero_()
                    wp.launch(plate_kernel, dim=N, inputs=[pos_d, wp.float64(cz), wp.float64(hg),
                              wp.float64(k_plate), f_d, react_d], device=d)
                    wp.launch(axpy_kernel, dim=N, inputs=[pos_d, wp.float64(dt_mu), f_d], device=d)
            # inextensibility projection once per physical step (after the inner relaxation)
            wp.launch(reshape_kernel, dim=fiber_off.shape[0] - 1,
                      inputs=[pos_d, foff_d, soff_d, srest_d, wp.int32(n_reshape_iter)], device=d)
            if trace_every and (pstep % max(1, n_phys // 40) == 0 or pstep == n_phys - 1):
                wp.synchronize_device(d)
                _fsi_trace.append({"pstep": pstep, "t_s": (pstep + 1) * dt_phys,
                                   "F_plate_pN": float(react_d.numpy()[0]),
                                   "p_excess_max_Pa": _p_max_hist[-1], "half_gap_um": hg})
        # ---- FSI final read + metrics (self-contained; the OFF monolithic loop below is not entered) ----
        V, V_cyto, area, dP = _refresh_turgor()
        wp.synchronize_device(d)
        react = react_d.numpy(); F_plate = float(react[0])
        pos_all = pos_d.numpy().astype(np.float64)
        net.pos = pos_all[:Nc]
        ppos = net.pos; cc = ppos.mean(axis=0)
        rxy = np.hypot(ppos[:, 0] - cc[0], ppos[:, 1] - cc[1]); R_eq = float(rxy.max())
        zc = ppos[:, 2] - cc[2]
        contact = np.abs(np.abs(zc) - final_gap) < 0.15 * seg + 0.05
        contact_radius = float(rxy[contact].max()) if contact.any() else 0.0
        tau_p = (R0 ** 2) / D_fsi                              # poroelastic redistribution time τ_p ~ L²/D
        metrics = {"strain": strain, "half_gap": final_gap, "F_plate_pN": F_plate, "trace": _fsi_trace,
                   "dt_s": float(dt_phys), "R_eq_um": R_eq, "contact_radius_um": contact_radius,
                   "dP_turgor_Pa": dP, "V_over_V0": V_cyto / V0, "V_hull_um3": V, "V_nuc_um3": V_nuc,
                   "gamma_apparent_pN_um": 0.5 * dP * R_eq, "gamma_apparent_mN_m": 0.5 * dP * R_eq * 1e-3,
                   "Nc": Nc, "Ne": Ne, "n_nuc": n_nuc, "has_mt": microtubule is not None,
                   "has_membrane": membrane is not None, "area_um2": area,
                   # FSI diagnostics
                   "fsi_on": True, "fsi_D_um2_s": D_fsi, "fsi_M_biot_Pa": M_fsi, "fsi_dx_um": dx_fsi,
                   "fsi_grid": f"{n_grid}³", "fsi_dt_phys_s": float(dt_phys), "fsi_n_phys": n_phys,
                   "fsi_inner": inner, "fsi_tau_p_s": float(tau_p), "fsi_biot_sub": _biot_sub,
                   "fsi_p_excess_max_Pa": float(np.abs(p_np[_int_mask]).max()) if _int_mask.any() else 0.0,
                   "fsi_p_excess_std_Pa": float(p_np[_int_mask].std()) if _int_mask.any() else 0.0,
                   "fsi_t_ramp_s": float(t_ramp_fsi), "v_press_um_s": v_press_um_s,
                   "fsi_gamma_solid_note": "6πηR retired; effective η emerges from Biot D+M (units.fiber_point_drag_solvent)"}
        # persist the pore-pressure field for the CFD-cell HTML (spatial, verifiable) — npz side channel so
        # the returned metrics stay JSON-serializable (no raw array in the dict).
        metrics["fsi_p_field"] = {"origin": [float(origin_fsi[0]), float(origin_fsi[1]), float(origin_fsi[2])],
                                  "dx": dx_fsi, "n": n_grid}
        _fld_out = biot_fsi.get("field_out_npz")
        if _fld_out:
            np.savez_compressed(_fld_out, p_field=p_np.astype(np.float32),
                                origin=np.array([float(origin_fsi[0]), float(origin_fsi[1]), float(origin_fsi[2])]),
                                dx=dx_fsi, n=n_grid, cortex_pos=pos_all[:Nc].astype(np.float32))
            metrics["fsi_p_field"]["npz"] = str(_fld_out)
        metrics.update(_mem_metrics(pos_all))
        return pos_all, metrics

    trace = []                                        # F(t) saturation trace (populated when trace_every>0)
    for step in range(n_steps):
        if step % turgor_every == 0:
            _refresh_turgor()
        if trace_every:
            react_d.zero_()                           # read-only diagnostic: react_d[0] = THIS step's plate force
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        if nT:
            wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        if has_xl:
            wp.launch(link_spring_kernel, dim=n_xl, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if mtp_d is not None:                 # MT tip ↔ cortex excluded-volume prop (Thread-C, engages under confinement)
            wp.launch(soft_contact_kernel, dim=mtp_d.shape[0], inputs=[pos_d, mtp_d, wp.float64(r_mt_c), wp.float64(k_mt_c), f_d], device=d)
        if has_myo:
            wp.launch(myosin_kernel, dim=cortex.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_area), f_d], device=d)  # cortex only
        if membrane is not None and n_mem == 0:
            # LUMPED membrane inward surface tension = turgor_kernel with a NEGATIVE (inward) dP·area/N (cortex
            # nodes only) — used ONLY when no real sheet is wired (the sheet replaces it; no double-count).
            wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_mem_area), f_d], device=d)
        if n_mem:                                     # real membrane sheet: area tension + Helfrich bending + ERM
            _membrane_forces()
        if n_nuc:
            wp.launch(nucleus_shell_kernel, dim=n_nuc,
                      inputs=[pos_d, wp.int32(Ne), centre, wp.float64(nucleus.R_nuc_um),
                              wp.float64(nucleus.k_chrom), wp.float64(nucleus.k_lamin),
                              wp.float64(nucleus.d_knee_um), wp.float64(nucleus.F_knee_pN), f_d], device=d)
            if _nuc_incompressible:                 # nucleoplasm volume-conservation pressure (ν→½)
                wp.launch(nucleus_pressure_kernel, dim=n_nuc,
                          inputs=[pos_d, wp.int32(Ne), centre_nuc, wp.float64(dP_nuc_area), f_d], device=d)
        # physical ramp: plate advances v_press·dt_real per step (half_gap → final over _N_ramp steps), else instant
        hg = half_gap if not _ramp_on else max(half_gap, R0 - float(v_press_um_s) * (step + 1) * _dt_real)
        if _ll:                             # S2 localized patch load (no plate): patch force + integrate
            # patch force on the CORTEX SURFACE only (dim=Nc) — the cap is a surface indentation; without
            # dim=Nc, interior MT-arm/MTOC/nucleus nodes inside the angular cap would also be pushed (they
            # are cortex-coupled → contaminates the readout in the native+full config). Cf. turgor dim=Nc.
            wp.launch(patch_force_kernel, dim=Nc, inputs=[pos_d, centre, _ll_axis, _ll_cos, _ll_f, f_d], device=d)
            wp.launch(axpy_kernel, dim=N, inputs=[pos_d, wp.float64(dt_mu), f_d], device=d)
        elif rigid_plate:                   # RIGID: step + hard-clamp |z|≤half_gap (exact confinement)
            wp.launch(rigid_plate_step_kernel, dim=N,
                      inputs=[pos_d, wp.float64(dt_mu), f_d, wp.float64(cz), wp.float64(hg), react_d], device=d)
        else:                               # SOFT penalty (historical): may under-confine a stiff cortex
            wp.launch(plate_kernel, dim=N, inputs=[pos_d, wp.float64(cz), wp.float64(hg),
                      wp.float64(k_plate), f_d, react_d], device=d)
            wp.launch(axpy_kernel, dim=N, inputs=[pos_d, wp.float64(dt_mu), f_d], device=d)
        if trace_every and step > 0 and (step % trace_every == 0 or step == n_steps - 1):
            wp.synchronize_device(d)  # skip step 0 (unphysical instant-clamp contact transient)
            if _ll:                   # localized load: trace the patch deflection(t) (saturation / convergence)
                _pp = pos_d.numpy()[:Nc]
                _pc = _pp - _pp.mean(0)       # centroid-relative → removes rigid COM drift from the free cell
                _pd = float((_pc[_pmask] @ _lax).mean()) - _patch_rest if _pmask.any() else 0.0
                trace.append({"step": step, "patch_deflection_um": _pd})
            else:
                trace.append({"step": step, "F_plate_pN": float(react_d.numpy()[0])})
        if (step + 1) % reshape_every == 0:
            wp.launch(reshape_kernel, dim=fiber_off.shape[0] - 1,
                      inputs=[pos_d, foff_d, soff_d, srest_d, wp.int32(n_reshape_iter)], device=d)
            if _turnover_on and has_xl:      # crosslink turnover: rest lengths relax toward current (Bell slip)
                wp.launch(xl_turnover_kernel, dim=n_xl,
                          inputs=[pos_d, xl_d, kxl_d, r0_d, wp.float64(_koff0_dt), wp.float64(_xbeta_kT)], device=d)

    V, V_cyto, area, dP = _refresh_turgor()
    react_d.zero_()
    if (not _ll) and rigid_plate:            # one more step to read the steady reaction (removed-disp / dt)
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        if nT:
            wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        if has_xl:
            wp.launch(link_spring_kernel, dim=n_xl, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if mtp_d is not None:
            wp.launch(soft_contact_kernel, dim=mtp_d.shape[0], inputs=[pos_d, mtp_d, wp.float64(r_mt_c), wp.float64(k_mt_c), f_d], device=d)
        if has_myo:
            wp.launch(myosin_kernel, dim=cortex.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_area), f_d], device=d)
        if membrane is not None and n_mem == 0:
            wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_mem_area), f_d], device=d)
        if n_mem:
            _membrane_forces()
        if n_nuc:
            wp.launch(nucleus_shell_kernel, dim=n_nuc,
                      inputs=[pos_d, wp.int32(Ne), centre, wp.float64(nucleus.R_nuc_um),
                              wp.float64(nucleus.k_chrom), wp.float64(nucleus.k_lamin),
                              wp.float64(nucleus.d_knee_um), wp.float64(nucleus.F_knee_pN), f_d], device=d)
            if _nuc_incompressible:                 # nucleoplasm volume-conservation pressure (ν→½)
                wp.launch(nucleus_pressure_kernel, dim=n_nuc,
                          inputs=[pos_d, wp.int32(Ne), centre_nuc, wp.float64(dP_nuc_area), f_d], device=d)
        wp.launch(rigid_plate_step_kernel, dim=N,
                  inputs=[pos_d, wp.float64(dt_mu), f_d, wp.float64(cz), wp.float64(half_gap), react_d], device=d)
    elif not _ll:
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
        pn = pos_all[Ne:Ne + n_nuc]                  # nucleus beads (the IF tail block, if any, is AFTER them)
        R_nuc_eq = float(np.hypot(pn[:, 0] - c[0], pn[:, 1] - c[1]).max())
        R_nuc_z = float(np.abs(pn[:, 2] - cz).max())     # polar half-height (flattens under the plate)
        nucleus_contact = bool((np.abs(pn[:, 2] - cz) > half_gap).any())
    else:
        R_nuc_eq = R_nuc_z = 0.0
        nucleus_contact = False
    # physical seconds per numerical step (overdamped: dt_phys = dt_mu * gamma_node = dt_mu / mu_phys).
    # dt_mu itself has units µm/pN (mobility descent step), NOT seconds — so trace/study time axes must use
    # this, computed even when the press-ramp is off (v_press=None).
    if _dt_real > 0.0:
        _dt_phys = _dt_real
    elif eta_bulk_Pa_s is not None:
        _dt_phys = dt_mu * (6.0 * np.pi * float(eta_bulk_Pa_s) * R0 / max(Nc, 1))
    else:
        _dt_phys = dt_mu               # no viscosity given → falls back to the numerical step (not physical s)
    metrics = {"strain": strain, "half_gap": half_gap, "F_plate_pN": F_plate, "trace": trace,
               "dt_s": float(_dt_phys),  # PHYSICAL seconds/step; trace physical time = step * dt_s
               "dt_mu_um_per_pN": float(dt_mu),  # the raw mobility descent step (µm/pN), for reference
               "R_eq_um": R_eq, "contact_radius_um": contact_radius, "dP_turgor_Pa": dP,
               "dP_from_force_Pa": dP_force, "V_over_V0": V_cyto / V0, "V_hull_um3": V, "V_nuc_um3": V_nuc,
               "gamma_apparent_pN_um": gamma_laplace, "gamma_apparent_mN_m": gamma_laplace * 1e-3,
               "R_nuc_eq_um": R_nuc_eq, "R_nuc_z_um": R_nuc_z, "nucleus_contact": nucleus_contact, "n_nuc": n_nuc,
               "K_vol_nuc_Pa": _K_vol_nuc, "V0_nuc_um3": V0_nuc, "V_nuc_cur_um3": float(V_nuc_cur),
               "nuc_vol_conserved": float(V_nuc_cur / V0_nuc) if V0_nuc > 0 else 1.0,
               "has_membrane": membrane is not None, "gamma_mem_channel_pN_um": gamma_mem_tot,
               "dP_mem_Pa": dP_mem, "area_um2": area, "A0_mem_um2": A0_mem,
               "k_plate": k_plate, "has_mt": microtubule is not None, "Nc": Nc, "Ne": Ne,
               "n_mt": int(_m.get("n_mt", 0)), "mtoc_idx": int(_m.get("mtoc_idx", -1)),
               # Piece #2 biphasic cytoplasm channels (dP_osm + dP_solid = dP_turgor when no setpoint)
               "dP_osm_Pa": dP_osm, "dP_solid_Pa": dP_solid, "K_drained_Pa": _K_drained,
               # fraction of the AVAILABLE drainage done: 0=undrained (V0_eff=V0), 1=fully drained (V0_eff=V_cyto)
               "drained_frac": float(np.clip((V0 - V0_eff) / max(V0 - V_cyto, 1e-9 * V0), 0.0, 1.0)),
               "V0_eff_over_V0": float(V0_eff / V0), "tau_osm_s": float(tau_osm), "load_time_s": load_time_s,
               "Lp_um_s_Pa": Lp_um_s_Pa, "f_excess": float(getattr(membrane, "f_excess", 0.0) or 0.0) if membrane is not None else 0.0,
               "xl_turnover_on": bool(_turnover_on), "xl_koff_per_s": xl_koff_per_s,
               "v_press_um_s": v_press_um_s, "dt_real_s": float(_dt_real), "n_ramp": int(_N_ramp),
               "t_ramp_s": float(_N_ramp * _dt_real) if _ramp_on else 0.0}
    if _ll:
        _pcn = pos_all[:Nc] - pos_all[:Nc].mean(0)             # centroid-relative → drift-free (see _patch_rest)
        _pnow = float((_pcn[_pmask] @ _lax).mean()) if _pmask.any() else 0.0
        metrics["patch_deflection_um"] = _pnow - _patch_rest   # signed: <0 inward (compression), >0 outward (tension)
        metrics["n_patch_nodes"] = int(_pmask.sum())
        metrics["f_node_pN"] = float(local_load.get("f_node_pN", 0.0))
    metrics.update(_mem_metrics(pos_all))
    return pos_all, metrics


@wp.kernel
def rigid_plate_step_kernel(
    pos: wp.array(dtype=wp.vec3d),
    step: wp.float64,                    # overdamped step dt/γ
    F: wp.array(dtype=wp.vec3d),
    cz: wp.float64,
    half_gap: wp.float64,
    reaction: wp.array(dtype=wp.float64),   # [0]=top-plate force (AFM), [1]=top contact count
):
    """RIGID parallel-plate = overdamped step x+=step·F, then HARD-clamp |z−cz| ≤ half_gap (a geometric
    CONSTRAINT, not a stiff penalty → no CFL throttle). The clamp removes the outward motion the plate
    resists; the AFM reaction is that removed displacement / step (= the force holding the node), summed on
    the top plate. This confines the cell EXACTLY to the gap (no penetration), unlike the soft penalty."""
    i = wp.tid()
    p = pos[i] + step * F[i]
    z = p[2] - cz
    if z > half_gap:
        wp.atomic_add(reaction, 0, (z - half_gap) / step)      # force = removed_disp / step
        wp.atomic_add(reaction, 1, wp.float64(1.0))
        p = wp.vec3d(p[0], p[1], cz + half_gap)
    if z < -half_gap:
        p = wp.vec3d(p[0], p[1], cz - half_gap)
    pos[i] = p


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
    from aleph.laws.forces_warp import _per_triple_alpha
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

    from aleph.laws.forces_warp import cytosim_bending_kernel
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
    from aleph.laws.forces_warp import _per_triple_alpha, cytosim_bending_kernel
    from aleph.laws.turgor_constants import (
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

    from aleph.laws.forces_warp import _per_triple_alpha, cytosim_bending_kernel
    from aleph.laws.turgor_constants import TURGOR_DP0, TURGOR_PI_IN0, VMIN_FRAC

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
    from aleph.laws.forces_warp import _per_triple_alpha, cytosim_bending_kernel
    from aleph.laws.fiber_network import mesoscale_reach
    from aleph.laws.hand_kmc import NMIIA_MYOSIN

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
