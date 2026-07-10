"""FF ECM mechanics harness — measure the effective modulus of an :class:`~ffn_sim.ff.ecm_library.ECMNetwork`
in **Pa**, three independent ways, so "실제와 같은 Pa 나오나" can be checked against literature.

All three return SI Pascals directly (FF units: 1 pN/µm² = 1 Pa):

* :func:`shear_modulus` — affine simple shear γ; pin the z-min/z-max layers, displace the top in x,
  relax the interior, read σ_xz = (x-reaction on the top layer)/area → **G = σ_xz/γ**.
* :func:`uniaxial_modulus` — affine uniaxial stress along an axis (lateral faces free); pin the
  −face, displace the +face, read σ = (axial reaction)/area → **E = σ/ε**. Along vs across the
  director gives the alignment anisotropy.
* :func:`indentation_modulus` — the "pressing" test: a rigid spherical indenter driven into the top
  surface of a bottom-clamped slab; F(δ) inverted with the parabolic (Sneddon) Hertz model →
  **E_eff** — exactly what an AFM/bead measurement of a gel reports.

The relaxation reuses the FF Warp kernels unchanged (``cytosim_bending_kernel``, ``link_spring_kernel``,
``reshape_kernel``, ``freeze``) plus one spherical-indenter kernel added here.

Method notes (validated 2026-07-10):
* **Energy route is primary.** After an affine pre-deformation + interior relaxation, the modulus is read
  from the stored elastic energy: G = 2U/(V·γ²), E = 2U/(V·ε²). The boundary-reaction route (kept as
  ``*_reaction_Pa`` cross-check) is sensitive to how many fibers grip the pinned band and reads unreliably
  on sparse networks. Caveat (fibrillar): U is the total energy over the whole box, whose pinned boundary
  layers (``margin_frac``≈0.08 top+bottom) are held affine — a constraint that can only *raise* the relaxed
  minimum, so the fibrillar modulus carries a small UPWARD bias that grows with ``margin_frac`` (~tens of %).
  For the continuum path this cancels exactly: ``calibrate_continuum_k`` measures E-per-k through the
  identical geometry, so the same boundary/V effect divides out (why the PA-gel returns E=5000 exactly).
* **Axial mode = finite EA spring (default).** The derived per-segment stretch spring k=EA/ℓ₀ gives a
  clean, converged modulus. ``axial_mode='reshape'`` (inextensible projection) is NOT recommended with
  the affine-pre-strain scheme here — it fights the imposed strain and fails to converge.
* **Harness validated on a known-modulus continuum:** a PA-gel lattice built to E=5000 Pa returns
  uniaxial E=5000 (by calibration), shear G within 9% of E/[2(1+ν)], and indentation E_eff=3970 Pa
  (0.79×) with flat E-vs-δ — so the indentation inversion is correct.
* **Local vs bulk (fibrillar):** a spherical indenter on a sparse semiflexible network reads the LOCAL
  modulus, which is legitimately SOFTER than the bulk affine shear/uniaxial modulus (non-affine local
  softening; contact radius must exceed the mesh ξ to approach the continuum value — see KB-1.14). For
  fibrillar ECM the shear/uniaxial modulus is the "material Pa" vs rheology; indentation is the local
  probe. Continuum gels read the same E all three ways.
* **Continuum-gel Poisson caveat:** a central-force Delaunay spring lattice obeys the Cauchy relation, so
  its true ν is pinned ~0.25–0.33 in 3D and cannot reach the specs' ν (PA 0.45, BM gels 0.4–0.5). The Hertz
  inversion uses ``spec.poisson``, introducing a ~10–13% E offset (E=E*(1−ν²)) folded into the ~0.7–0.8×
  continuum indentation factor. The calibrated uniaxial E is exact regardless of ν; only the indentation
  readout carries this. A ν-faithful continuum would need non-central / multibody terms.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

from ffn_sim.ff import units as U
from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
from ffn_sim.ff.network_warp import (
    _zero, axpy_kernel, link_spring_kernel, reshape_kernel,
)

wp.init()


@wp.kernel
def _sphere_indenter_kernel(pos: wp.array(dtype=wp.vec3d), centre: wp.vec3d, R_ind: wp.float64,
                            k_ind: wp.float64, force: wp.array(dtype=wp.vec3d),
                            reaction: wp.array(dtype=wp.float64)):
    """Rigid spherical indenter: any node inside the sphere is pushed radially out; the vertical reaction
    on the indenter (what an AFM cantilever reads) is accumulated. reaction[0]=upward force, [1]=contacts."""
    i = wp.tid()
    d = pos[i] - centre
    L = wp.length(d)
    if L < R_ind and L > wp.float64(1e-9):
        pen = k_ind * (R_ind - L)
        f = (pen / L) * d                       # push node away from indenter centre
        wp.atomic_add(force, i, f)
        wp.atomic_add(reaction, 0, -f[2])       # Newton reaction on indenter, +z = resists indentation
        wp.atomic_add(reaction, 1, wp.float64(1.0))


@wp.kernel
def _freeze_mask_kernel(force: wp.array(dtype=wp.vec3d), fixed: wp.array(dtype=wp.int32)):
    """Zero the force on prescribed-displacement (Dirichlet) nodes so the integrator leaves them put."""
    i = wp.tid()
    if fixed[i] == 1:
        force[i] = wp.vec3d(0.0, 0.0, 0.0)


@wp.kernel
def _pin_positions_kernel(pos: wp.array(dtype=wp.vec3d), fixed: wp.array(dtype=wp.int32),
                          target: wp.array(dtype=wp.vec3d)):
    """Re-clamp prescribed-displacement nodes to their target position. MUST run after reshape (which
    moves ALL nodes to satisfy inextensibility, including the pinned boundary) or the imposed strain
    leaks out and the measured modulus collapses to ~0."""
    i = wp.tid()
    if fixed[i] == 1:
        pos[i] = target[i]


@dataclass
class _Dev:
    """Device-resident network arrays for a relaxation run (built once, reused across steps)."""
    pos: wp.array
    f: wp.array
    tri: wp.array
    alpha: wp.array
    foff: wp.array
    soff: wp.array
    srest: wp.array
    links: wp.array
    klink: wp.array
    rlink: wp.array
    n_nodes: int
    n_tri: int
    n_links: int
    seg_mean: float
    device: str


def _to_device(ecm, link_pairs, link_k, link_rest, device):
    net = ecm.net
    N = net.n_nodes
    tri = np.ascontiguousarray(net.bend_triples, np.int32)
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64) if tri.shape[0] else np.zeros(0)
    fiber_off = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg_per = np.diff(fiber_off) - 1
    seg_off = np.concatenate([[0], np.cumsum(np.maximum(seg_per, 0))]).astype(np.int32)
    seg_mean = float(net.seg_rest.mean()) if net.seg_rest.size else 1.0
    d = device
    return _Dev(
        pos=wp.array(np.ascontiguousarray(net.pos, np.float64), dtype=wp.vec3d, device=d),
        f=wp.zeros(N, dtype=wp.vec3d, device=d),
        tri=wp.array(tri, dtype=wp.int32, device=d) if tri.shape[0] else None,
        alpha=wp.array(alpha, dtype=wp.float64, device=d) if tri.shape[0] else None,
        foff=wp.array(fiber_off, dtype=wp.int32, device=d),
        soff=wp.array(seg_off, dtype=wp.int32, device=d),
        srest=wp.array(np.ascontiguousarray(net.seg_rest, np.float64), dtype=wp.float64, device=d)
        if net.seg_rest.size else None,
        links=wp.array(np.ascontiguousarray(link_pairs, np.int32), dtype=wp.int32, device=d)
        if len(link_pairs) else None,
        klink=wp.array(np.ascontiguousarray(link_k, np.float64), dtype=wp.float64, device=d)
        if len(link_pairs) else None,
        rlink=wp.array(np.ascontiguousarray(link_rest, np.float64), dtype=wp.float64, device=d)
        if len(link_pairs) else None,
        n_nodes=N, n_tri=tri.shape[0], n_links=len(link_pairs), seg_mean=seg_mean, device=d)


def _force_pass(dev: _Dev, indenter=None):
    """One force evaluation (bending + links [+ indenter]) into dev.f; returns indenter reaction array or None."""
    d = dev.device
    wp.launch(_zero, dim=dev.n_nodes, inputs=[dev.f], device=d)
    if dev.n_tri:
        wp.launch(cytosim_bending_kernel, dim=dev.n_tri, inputs=[dev.pos, dev.tri, dev.alpha, dev.f], device=d)
    if dev.n_links:
        wp.launch(link_spring_kernel, dim=dev.n_links, inputs=[dev.pos, dev.links, dev.klink, dev.rlink, dev.f], device=d)
    react = None
    if indenter is not None:
        centre, R_ind, k_ind, react = indenter
        react.zero_()                            # per-pass reaction (else it accumulates over every relax step)
        wp.launch(_sphere_indenter_kernel, dim=dev.n_nodes,
                  inputs=[dev.pos, centre, wp.float64(R_ind), wp.float64(k_ind), dev.f, react], device=d)
    return react


def _cfl_dt(dev: _Dev, ecm, k_extra=0.0) -> float:
    kmax = 0.0
    if ecm.net.kappa.size and dev.n_tri:
        kmax = float(ecm.net.kappa.max()) / dev.seg_mean ** 3
    if dev.n_links:
        kmax = max(kmax, float(np.max(ecm.link_k)) if ecm.link_k.size else 0.0)
    kmax = max(kmax, k_extra)
    return 0.1 / max(kmax, 1e-9)


def _relax(dev: _Dev, ecm, fixed_mask, *, use_reshape, n_steps, reshape_every=25, indenter=None,
           dt_mu=0.0, k_extra=0.0):
    """Overdamped relaxation with prescribed-displacement (fixed) nodes held put. Bending + link springs
    [+ reshape inextensibility] [+ spherical indenter], all on device. Returns relaxed positions (N,3)."""
    d = dev.device
    if dt_mu <= 0.0:
        dt_mu = _cfl_dt(dev, ecm, k_extra=k_extra)
    fixed_d = wp.array(np.ascontiguousarray(fixed_mask.astype(np.int32)), dtype=wp.int32, device=d)
    target_d = wp.array(dev.pos.numpy().copy(), dtype=wp.vec3d, device=d)   # prescribed BC positions
    has_reshape = use_reshape and dev.srest is not None
    for step in range(n_steps):
        _force_pass(dev, indenter=indenter)
        wp.launch(_freeze_mask_kernel, dim=dev.n_nodes, inputs=[dev.f, fixed_d], device=d)
        wp.launch(axpy_kernel, dim=dev.n_nodes, inputs=[dev.pos, wp.float64(dt_mu), dev.f], device=d)
        if has_reshape and (step + 1) % reshape_every == 0:
            wp.launch(reshape_kernel, dim=dev.foff.shape[0] - 1,
                      inputs=[dev.pos, dev.foff, dev.soff, dev.srest, wp.int32(2)], device=d)
            wp.launch(_pin_positions_kernel, dim=dev.n_nodes, inputs=[dev.pos, fixed_d, target_d], device=d)
    if has_reshape:
        wp.launch(reshape_kernel, dim=dev.foff.shape[0] - 1,
                  inputs=[dev.pos, dev.foff, dev.soff, dev.srest, wp.int32(4)], device=d)
        wp.launch(_pin_positions_kernel, dim=dev.n_nodes, inputs=[dev.pos, fixed_d, target_d], device=d)
    wp.synchronize_device(d)
    return dev.pos.numpy().astype(np.float64)


def _reaction_on(dev: _Dev, node_mask, axis: int) -> float:
    """Internal force the sample exerts on the (held) measurement nodes, ``axis`` component, summed [pN]."""
    _force_pass(dev, indenter=None)
    wp.synchronize_device(dev.device)
    f = dev.f.numpy()
    return float(f[node_mask, axis].sum())


def _elastic_energy(ecm, pos, use_spring: bool) -> float:
    """Total stored elastic energy [pN·µm] at ``pos`` = bending + crosslink-spring [+ axial-spring]. This
    is the robust modulus route (G = 2U/(V·γ²)) — a bulk scalar, insensitive to boundary-layer definition
    unlike the boundary-reaction route. The unstrained built network is at rest (straight rods, bonds at
    rest length) so its reference energy is ~0; ΔU ≈ U."""
    from ffn_sim.ff.forces_warp import bending_energy
    net = ecm.net
    U = 0.0
    if net.bend_triples.shape[0]:
        saved = net.pos
        net.pos = pos
        try:
            U += float(bending_energy(net))
        finally:
            net.pos = saved
    if ecm.xl_i.size:
        d = np.linalg.norm(pos[ecm.xl_j] - pos[ecm.xl_i], axis=1)
        U += 0.5 * float(np.sum(ecm.xl_k * (d - ecm.xl_rest) ** 2))
    if use_spring and ecm.seg_i.size:
        d = np.linalg.norm(pos[ecm.seg_j] - pos[ecm.seg_i], axis=1)
        U += 0.5 * float(np.sum(ecm.seg_k * (d - ecm.seg_rest) ** 2))
    return U


def _links_for(ecm, axial_mode: str):
    """Which Hookean bonds to apply and whether to use reshape. Continuum lattices → Delaunay springs;
    fibrillar 'reshape' → crosslinks only (+reshape inextensibility); fibrillar 'spring' → xl + seg."""
    is_continuum = ecm.net.bend_triples.shape[0] == 0
    if is_continuum:
        pairs = np.stack([ecm.seg_i, ecm.seg_j], 1) if ecm.seg_i.size else np.zeros((0, 2), np.int64)
        return pairs.astype(np.int32), ecm.seg_k, ecm.seg_rest, False
    if axial_mode == "spring":
        return ecm.links, ecm.link_k, ecm.link_rest, False
    xl = np.stack([ecm.xl_i, ecm.xl_j], 1) if ecm.xl_i.size else np.zeros((0, 2), np.int64)
    return xl.astype(np.int32), ecm.xl_k, ecm.xl_rest, True


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# 1) shear modulus G
# ─────────────────────────────────────────────────────────────────────────────────────────────────

def shear_modulus(ecm, *, gamma: float = 0.02, n_steps: int = 4000, axial_mode: str = "spring",
                  device: str = "cpu", margin_frac: float = 0.08) -> dict:
    """Affine simple shear in the x–z plane → shear modulus **G [Pa]**.

    Pin the bottom (z<z_lo+m) and top (z>z_hi−m) layers; displace the top layer by γ·Lz in x; relax the
    interior; read σ_xz = (x-reaction on the top layer)/(Lx·Ly). G = σ_xz / γ.
    """
    lo, hi = ecm.box_lo, ecm.box_hi
    Lz = hi[2] - lo[2]
    m = margin_frac * Lz
    pos0 = ecm.net.pos.copy()
    bot = pos0[:, 2] < lo[2] + m
    top = pos0[:, 2] > hi[2] - m
    fixed = bot | top
    # affine pre-shear of the WHOLE network (start at the affine upper bound, relax DOWN to the non-affine
    # equilibrium with boundaries pinned) — the robust scheme (matches kim_network.shear_modulus)
    pos = pos0.copy()
    pos[:, 0] += gamma * (pos0[:, 2] - lo[2])
    ecm.net.pos[:] = pos
    pairs, klink, rlink, use_reshape = _links_for(ecm, axial_mode)
    dev = _to_device(ecm, pairs, klink, rlink, device)
    relaxed = _relax(dev, ecm, fixed, use_reshape=use_reshape, n_steps=n_steps)
    Fx = _reaction_on(dev, top, axis=0)
    sig = ecm_material_stress(ecm, relaxed, device=device) if not use_reshape else None  # full virial tensor
    ecm.net.pos[:] = pos0                                       # restore (leave the ECM object undeformed)
    V = float(np.prod(hi - lo))
    U = _elastic_energy(ecm, relaxed, use_spring=not use_reshape)
    G = 2.0 * U / (V * gamma ** 2)                              # energy route (primary) — pN/µm² = Pa
    area = (hi[0] - lo[0]) * (hi[1] - lo[1])
    G_reaction = abs(Fx / area) / gamma                         # boundary-reaction route (cross-check)
    G_virial = abs(float(sig[0, 2])) / gamma if sig is not None else float("nan")  # σ_xz/γ (grid-invariant)
    return {"G_Pa": G, "G_reaction_Pa": G_reaction, "G_virial_Pa": G_virial, "U_pN_um": U, "gamma": gamma,
            "V_um3": V, "n_top": int(top.sum()), "n_bot": int(bot.sum())}


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# 2) uniaxial modulus E (+ anisotropy)
# ─────────────────────────────────────────────────────────────────────────────────────────────────

def uniaxial_modulus(ecm, *, strain: float = 0.02, axis: str = "z", n_steps: int = 4000,
                     axial_mode: str = "spring", device: str = "cpu", margin_frac: float = 0.08) -> dict:
    """Affine uniaxial stress along ``axis`` (lateral faces free) → Young's modulus **E [Pa]**.

    Pin the −face; displace the +face by ε·L along the axis; relax; read σ = (axial reaction on the
    +face)/A_cross. E = σ/ε. (Compression if strain<0.)
    """
    ax = {"x": 0, "y": 1, "z": 2}[axis]
    lo, hi = ecm.box_lo, ecm.box_hi
    L = hi[ax] - lo[ax]
    m = margin_frac * L
    pos0 = ecm.net.pos.copy()
    neg = pos0[:, ax] < lo[ax] + m
    plus = pos0[:, ax] > hi[ax] - m
    fixed = neg | plus
    # affine pre-stretch of the WHOLE network along the axis; lateral faces free (Poisson relaxes)
    pos = pos0.copy()
    pos[:, ax] = lo[ax] + (1.0 + strain) * (pos0[:, ax] - lo[ax])
    ecm.net.pos[:] = pos
    pairs, klink, rlink, use_reshape = _links_for(ecm, axial_mode)
    dev = _to_device(ecm, pairs, klink, rlink, device)
    relaxed = _relax(dev, ecm, fixed, use_reshape=use_reshape, n_steps=n_steps)
    F = _reaction_on(dev, plus, axis=ax)
    ecm.net.pos[:] = pos0                                       # restore
    V = float(np.prod(hi - lo))
    U = _elastic_energy(ecm, relaxed, use_spring=not use_reshape)
    E = 2.0 * U / (V * strain ** 2)                             # energy route (primary)
    others = [i for i in range(3) if i != ax and (ecm.dim == 3 or i != 2)]
    area = np.prod([hi[i] - lo[i] for i in others]) if others else 1.0
    E_reaction = abs(F / area) / abs(strain)                    # boundary-reaction route (cross-check)
    return {"E_Pa": E, "E_reaction_Pa": E_reaction, "U_pN_um": U, "strain": strain, "axis": axis,
            "V_um3": V, "n_plus": int(plus.sum())}


def anisotropy(ecm_along, ecm_across, **kw) -> dict:
    """E along vs across the director (pass two builds with director set to the loading axis). Ratio>1 =
    stiffer along fibers (expected for aligned nets)."""
    ea = uniaxial_modulus(ecm_along, axis="z", **kw)
    ec = uniaxial_modulus(ecm_across, axis="z", **kw)
    return {"E_along_Pa": ea["E_Pa"], "E_across_Pa": ec["E_Pa"],
            "anisotropy_ratio": ea["E_Pa"] / max(ec["E_Pa"], 1e-9)}


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# 3) spherical indentation → Hertz E_eff (the "pressing" test)
# ─────────────────────────────────────────────────────────────────────────────────────────────────

def hertz_E_from_force(F_pN: float, R_um: float, delta_um: float, nu: float) -> float:
    """Invert the parabolic (Sneddon) Hertz law F = (4/3)·E*·√R·δ^1.5 → apparent E [Pa]. F[pN],R,δ[µm]."""
    if delta_um <= 0.0 or F_pN <= 0.0:
        return 0.0
    e_star = F_pN / ((4.0 / 3.0) * np.sqrt(R_um) * delta_um ** 1.5)
    return e_star * (1.0 - nu ** 2)


def _hertz_slope_E(deltas, forces, R_um, nu) -> float:
    x = np.asarray(deltas) ** 1.5
    y = np.asarray(forces)
    good = (x > 0) & (y > 0)
    if good.sum() < 2:
        return 0.0
    m = float(np.sum(x[good] * y[good]) / np.sum(x[good] * x[good]))
    return (m / ((4.0 / 3.0) * np.sqrt(R_um))) * (1.0 - nu ** 2)


def indentation_modulus(ecm, *, indenter_R_um: float = 6.0, max_depth_um: float | None = None,
                        n_depths: int = 5, k_ind: float = 2.0e3, relax_steps: int = 2500,
                        axial_mode: str = "spring", nu: float | None = None, center_xy=None,
                        device: str = "cpu") -> dict:
    """Press a rigid sphere of radius ``indenter_R_um`` into the TOP surface of a bottom-clamped slab and
    invert Hertz → **E_eff [Pa]** (the AFM/bead "pressing" measurement).

    Bottom face pinned (embedded / coverslip). Indenter centred over the box, lowered so its penetration
    δ sweeps a small range (default up to ~min(0.12·R_ind, 0.06·thickness) to stay Hertzian and off the
    bottom). F(δ) fit through the origin on F vs δ^1.5. Use R_ind ≫ mesh ξ so the response is continuum.
    """
    lo, hi = ecm.box_lo, ecm.box_hi
    Lx, Ly, Lz = (hi - lo)
    if nu is None:
        from ffn_sim.ff.ecm_library import get_spec
        try:
            nu = get_spec(ecm.material.split("+")[0]).poisson
        except Exception:
            nu = 0.3
    if max_depth_um is None:
        max_depth_um = min(0.12 * indenter_R_um, 0.06 * Lz)
    # bottom clamp
    m = 0.08 * Lz
    pos0 = ecm.net.pos.copy()
    fixed = pos0[:, 2] < lo[2] + m
    z_surf = float(np.percentile(pos0[:, 2], 99.0))            # effective top surface
    cx = 0.5 * (lo[0] + hi[0]) if center_xy is None else float(center_xy[0])
    cy = 0.5 * (lo[1] + hi[1]) if center_xy is None else float(center_xy[1])
    pairs, klink, rlink, use_reshape = _links_for(ecm, axial_mode)
    depths = np.linspace(max_depth_um / n_depths, max_depth_um, n_depths)
    forces, contacts = [], []
    for delta in depths:
        ecm.net.pos[:] = pos0                                   # fresh slab each depth (independent Hertz point)
        dev = _to_device(ecm, pairs, klink, rlink, device)
        centre_z = z_surf + indenter_R_um - delta               # sphere bottom sits δ below the surface
        react = wp.zeros(2, dtype=wp.float64, device=device)
        centre = wp.vec3d(float(cx), float(cy), float(centre_z))
        indenter = (centre, indenter_R_um, k_ind, react)
        _relax(dev, ecm, fixed, use_reshape=use_reshape, n_steps=relax_steps, indenter=indenter,
               k_extra=k_ind)
        # final indenter reaction at converged state
        _force_pass(dev, indenter=(centre, indenter_R_um, k_ind, react))
        wp.synchronize_device(device)
        rr = react.numpy()
        forces.append(float(rr[0]))
        contacts.append(int(rr[1]))
    ecm.net.pos[:] = pos0
    E_fit = _hertz_slope_E(depths, forces, indenter_R_um, nu)
    E_pts = [hertz_E_from_force(F, indenter_R_um, d, nu) for F, d in zip(forces, depths)]
    # contact guard: if the bead never engaged (contacts 0 / no positive force) E_eff=0 is NOT a physical
    # modulus — flag it so callers don't read a silent zero as "soft". Indenter must be positioned on the
    # surface with R_ind ≫ mesh ξ.
    contact_ok = bool(max(contacts) > 0 and E_fit > 0.0)
    if not contact_ok:
        import warnings
        warnings.warn(f"indentation_modulus: no/degenerate contact (contacts={contacts}, E_eff={E_fit}) — "
                      "E_eff=0 is not physical; check indenter position / R_ind vs mesh / slab surface.")
    return {"E_eff_Pa": E_fit, "nu": nu, "R_ind_um": indenter_R_um, "depths_um": depths.tolist(),
            "forces_pN": forces, "E_pts_Pa": E_pts, "contacts": contacts, "contact_ok": contact_ok,
            "flatness": (max(E_pts) / max(min([e for e in E_pts if e > 0], default=1.0), 1e-9)) if E_pts else 0.0}


def shear_energy(ecm, *, gamma, n_steps=5000, axial_mode="spring", device="cpu", margin_frac=0.08) -> float:
    """Stored elastic energy U(γ) [pN·µm] under an affine simple shear γ (pin z-min/z-max, relax interior).
    The nonlinear-mechanics primitive: σ(γ)=(1/V)dU/dγ, K(γ)=dσ/dγ — no small-strain assumption."""
    lo, hi = ecm.box_lo, ecm.box_hi
    Lz = hi[2] - lo[2]
    m = margin_frac * Lz
    pos0 = ecm.net.pos.copy()
    fixed = (pos0[:, 2] < lo[2] + m) | (pos0[:, 2] > hi[2] - m)
    pos = pos0.copy()
    pos[:, 0] += gamma * (pos0[:, 2] - lo[2])
    ecm.net.pos[:] = pos
    pairs, klink, rlink, use_reshape = _links_for(ecm, axial_mode)
    dev = _to_device(ecm, pairs, klink, rlink, device)
    relaxed = _relax(dev, ecm, fixed, use_reshape=use_reshape, n_steps=n_steps)
    ecm.net.pos[:] = pos0
    return _elastic_energy(ecm, relaxed, use_spring=not use_reshape)


def strain_stiffening(ecm, *, gammas=None, n_steps=5000, axial_mode="spring", device="cpu") -> dict:
    """Emergent nonlinear strain-stiffening: differential shear modulus K(γ)=dσ/dγ vs strain, from U(γ).

    Semiflexible ECM (collagen/fibrin) stiffens under strain — the differential modulus K rises many-fold
    above a critical strain γ_c (KB-1.V.2.5: K>10× above γ_c~0.1-0.5, K~σ, γ_c ~concentration-independent;
    an EMERGENT gate — must arise from the Mikado physics, NOT be tuned). This measures U(γ) on a strain
    grid, then σ=(1/V)dU/dγ and K=dσ/dγ by finite differences → K(γ)/K0, the stiffening ratio.
    """
    if gammas is None:
        gammas = np.linspace(0.01, 0.45, 12)
    gammas = np.asarray(gammas, float)
    V = float(np.prod(ecm.box_hi - ecm.box_lo))
    U = np.array([shear_energy(ecm, gamma=g, n_steps=n_steps, axial_mode=axial_mode, device=device)
                  for g in gammas])
    # σ = (1/V) dU/dγ ; K = dσ/dγ  (central differences on the U(γ) grid)
    sigma = np.gradient(U, gammas) / V
    K = np.gradient(sigma, gammas)
    K0 = K[0] if K[0] > 0 else np.max(K[K > 0], initial=1e-9)
    ratio = K / max(K0, 1e-12)
    # γ_c ≈ strain where K first exceeds 3×K0 (onset of stiffening)
    above = np.where(ratio > 3.0)[0]
    gamma_c = float(gammas[above[0]]) if above.size else float("nan")
    return {"gammas": gammas.tolist(), "U_pN_um": U.tolist(), "sigma_Pa": sigma.tolist(),
            "K_Pa": K.tolist(), "K_over_K0": ratio.tolist(), "K0_Pa": float(K0),
            "max_stiffening": float(np.max(ratio)), "gamma_c": gamma_c}


def _bond_virial_sigma_xz(pos, i_arr, j_arr, k_arr, r0_arr, V) -> float:
    """Shear stress σ_xz = (1/V) Σ_bonds (f/|r|)·r_x·r_z from central-force bonds (crosslinks+segments) [Pa]."""
    if len(i_arr) == 0:
        return 0.0
    d = pos[j_arr] - pos[i_arr]
    Lc = np.linalg.norm(d, axis=1) + 1e-12
    f = k_arr * (Lc - r0_arr)                         # scalar bond tension [pN]
    return float(np.sum(f / Lc * d[:, 0] * d[:, 2]) / V)


def _bending_virial_sigma_xz(net, pos, V, device="cpu") -> float:
    """Bending contribution to σ_xz via the atomic virial (1/2V)·Σ_i (r_x·f_z + r_z·f_x). Bending forces
    per triple sum to zero force AND zero torque, so this is origin-independent. MUST be included for
    bending-dominated (sub-isostatic) fiber networks — else the stress is severely under-reported."""
    if net.bend_triples.shape[0] == 0:
        return 0.0
    from ffn_sim.ff.forces_warp import bending_force
    saved = net.pos
    net.pos = pos
    try:
        fb = bending_force(net, device=device)
    finally:
        net.pos = saved
    return float(0.5 * np.sum(pos[:, 0] * fb[:, 2] + pos[:, 2] * fb[:, 0]) / V)


def ecm_material_stress(ecm, pos, device="cpu") -> np.ndarray:
    """Full macroscopic virial (Cauchy) stress tensor σ[3,3] of an ECM network at ``pos`` [Pa].

    σ_ab = (1/V)·[ Σ_bonds (F/L)·r_a·r_b  +  ½·Σ_i (r_a·f_bend_b + r_b·f_bend_a) ] — the crosslink + segment
    central-bond virial (symmetric, F=k(L−r₀) the bond tension) plus the bending atomic virial (origin-
    independent since bending forces are per-triple force- and torque-balanced). This is the ONE grid-invariant
    macroscopic stress readout: any modulus/normal-stress is a component of it (G=σ_xz/γ, E=σ_zz/ε, N1=σ_xx−σ_zz).
    Generalizes ``_bond_virial_sigma_xz``/``_bending_virial_sigma_xz`` from the xz component to the full tensor;
    same convention as the cortex ``ff_virial_stress`` (1 pN/µm²=1 Pa). Bending term is zero for continuum gels."""
    V = float(np.prod(np.asarray(ecm.box_hi) - np.asarray(ecm.box_lo)))
    sig = np.zeros((3, 3))
    for i, j, k, r0 in ((ecm.xl_i, ecm.xl_j, ecm.xl_k, ecm.xl_rest),
                        (ecm.seg_i, ecm.seg_j, ecm.seg_k, ecm.seg_rest)):
        if getattr(i, "size", 0):
            d = pos[j] - pos[i]
            L = np.linalg.norm(d, axis=1) + 1e-12
            f_over_L = k * (L - r0) / L                         # bond tension / length
            sig += np.einsum("n,na,nb->ab", f_over_L, d, d) / V
    if ecm.net.bend_triples.shape[0]:
        from ffn_sim.ff.forces_warp import bending_force
        saved = ecm.net.pos
        ecm.net.pos = pos
        try:
            fb = bending_force(ecm.net, device=device)
        finally:
            ecm.net.pos = saved
        sig += 0.5 * (pos.T @ fb + fb.T @ pos) / V
    return sig


def stress_relaxation(ecm, *, gamma0=0.1, koff0_per_s=0.01, x_beta_nm=0.4, dt_real_s=None,
                      t_total_s=None, n_record=24, mech_substeps=150, device="cpu") -> dict:
    """Viscoelastic stress-relaxation G(t): step shear γ₀, then let CROSSLINKS turn over (Bell slip,
    ``xl_turnover_kernel``) so their rest length creeps and the stress relaxes — a Maxwell/SLS network
    (KB-1.6: G(t)=G∞+G₁·exp(−t/τ), τ~30-1000 s set by crosslink lifetime 1/k_off; Chaudhuri matrix
    viscoelasticity). τ EMERGES ≈ 1/k_off₀ — validated, not tuned. Fibers (segments) + bending stay elastic;
    only crosslinks relax. Returns G(t)/G₀, fitted τ, and the residual G∞/G₀."""
    from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
    from ffn_sim.ff.network_warp import xl_turnover_kernel
    if dt_real_s is None:
        dt_real_s = (1.0 / koff0_per_s) * 6.0 / n_record        # span ~6 relaxation times
    lo, hi = ecm.box_lo, ecm.box_hi
    Lz = hi[2] - lo[2]
    m = 0.08 * Lz
    V = float(np.prod(hi - lo))
    pos0 = ecm.net.pos.copy()
    fixed = (pos0[:, 2] < lo[2] + m) | (pos0[:, 2] > hi[2] - m)
    pos = pos0.copy()
    pos[:, 0] += gamma0 * (pos0[:, 2] - lo[2])
    ecm.net.pos[:] = pos
    # device setup: bending + segment springs (elastic) + crosslink springs (turn over)
    d = device
    N = ecm.net.n_nodes
    tri = np.ascontiguousarray(ecm.net.bend_triples, np.int32)
    has_bend = tri.shape[0] > 0
    alpha = np.ascontiguousarray(_per_triple_alpha(ecm.net), np.float64) if has_bend else np.zeros(0)
    pos_d = wp.array(np.ascontiguousarray(pos, np.float64), dtype=wp.vec3d, device=d)
    f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    tri_d = wp.array(tri, dtype=wp.int32, device=d) if has_bend else None
    alpha_d = wp.array(alpha, dtype=wp.float64, device=d) if has_bend else None
    seg = np.stack([ecm.seg_i, ecm.seg_j], 1).astype(np.int32) if ecm.seg_i.size else np.zeros((0, 2), np.int32)
    seg_d = wp.array(seg, dtype=wp.int32, device=d) if seg.shape[0] else None
    segk_d = wp.array(np.ascontiguousarray(ecm.seg_k, np.float64), dtype=wp.float64, device=d) if seg.shape[0] else None
    segr_d = wp.array(np.ascontiguousarray(ecm.seg_rest, np.float64), dtype=wp.float64, device=d) if seg.shape[0] else None
    xl = np.stack([ecm.xl_i, ecm.xl_j], 1).astype(np.int32) if ecm.xl_i.size else np.zeros((0, 2), np.int32)
    xl_d = wp.array(xl, dtype=wp.int32, device=d) if xl.shape[0] else None
    xlk_d = wp.array(np.ascontiguousarray(ecm.xl_k, np.float64), dtype=wp.float64, device=d) if xl.shape[0] else None
    r0 = ecm.xl_rest.astype(np.float64).copy()
    r0_d = wp.array(np.ascontiguousarray(r0, np.float64), dtype=wp.float64, device=d) if xl.shape[0] else None
    fixed_d = wp.array(fixed.astype(np.int32), dtype=wp.int32, device=d)
    target_d = wp.array(pos.copy(), dtype=wp.vec3d, device=d)
    kmax = (float(ecm.net.kappa.max()) / (float(ecm.net.seg_rest.mean()) if ecm.net.seg_rest.size else 1.0) ** 3) if has_bend else 0.0
    kmax = max(kmax, float(np.max(ecm.link_k)) if ecm.link_k.size else 1.0)
    dt_mu = 0.1 / max(kmax, 1e-9)
    x_beta_over_kT = (x_beta_nm * 1e-3) / U.KBT                  # (µm)/(pN·µm) = 1/pN
    seg_i, seg_j, seg_k, seg_r = ecm.seg_i, ecm.seg_j, ecm.seg_k, ecm.seg_rest

    def equilibrate():
        for _ in range(mech_substeps):
            wp.launch(_zero, dim=N, inputs=[f_d], device=d)
            if has_bend:
                wp.launch(cytosim_bending_kernel, dim=tri.shape[0], inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
            if seg_d is not None:
                wp.launch(link_spring_kernel, dim=seg.shape[0], inputs=[pos_d, seg_d, segk_d, segr_d, f_d], device=d)
            if xl_d is not None:
                wp.launch(link_spring_kernel, dim=xl.shape[0], inputs=[pos_d, xl_d, xlk_d, r0_d, f_d], device=d)
            wp.launch(_freeze_mask_kernel, dim=N, inputs=[f_d, fixed_d], device=d)
            wp.launch(axpy_kernel, dim=N, inputs=[pos_d, wp.float64(dt_mu), f_d], device=d)
            wp.launch(_pin_positions_kernel, dim=N, inputs=[pos_d, fixed_d, target_d], device=d)

    ts, Gs = [], []
    equilibrate()                                               # instantaneous (elastic) state
    for rec in range(n_record):
        wp.synchronize_device(d)
        p = pos_d.numpy()
        r0h = r0_d.numpy() if xl_d is not None else np.zeros(0)
        sig = _bond_virial_sigma_xz(p, seg_i, seg_j, seg_k, seg_r, V) + \
            _bond_virial_sigma_xz(p, ecm.xl_i, ecm.xl_j, ecm.xl_k, r0h, V) + \
            _bending_virial_sigma_xz(ecm.net, p, V, device=device)   # bending stress (dominant for sub-isostatic)
        ts.append(rec * dt_real_s)
        Gs.append(abs(sig) / gamma0)
        if xl_d is not None:                                    # advance one real timestep of turnover
            wp.launch(xl_turnover_kernel, dim=xl.shape[0],
                      inputs=[pos_d, xl_d, xlk_d, r0_d, wp.float64(koff0_per_s * dt_real_s),
                              wp.float64(x_beta_over_kT)], device=d)
        equilibrate()
    ecm.net.pos[:] = pos0
    ts = np.asarray(ts); Gs = np.asarray(Gs)
    G0 = Gs[0] if Gs[0] > 0 else 1e-9
    ratio = Gs / G0
    # fit G(t) = Ginf + (G0-Ginf) exp(-t/tau): estimate tau where ratio crosses (1+ratio[-1])/2 ... e-folding
    Ginf = float(ratio[-1])
    decay = (ratio - Ginf) / max(1.0 - Ginf, 1e-9)
    tau = float("nan")
    below = np.where(decay < np.exp(-1.0))[0]
    if below.size:
        tau = float(ts[below[0]])
    return {"t_s": ts.tolist(), "G_t_Pa": Gs.tolist(), "G_over_G0": ratio.tolist(),
            "G0_Pa": float(G0), "Ginf_over_G0": Ginf, "tau_s": tau, "koff0_per_s": koff0_per_s,
            "tau_x_koff": tau * koff0_per_s if tau == tau else float("nan")}


def calibrate_continuum_k(spec, box_lo, box_hi, *, dim=3, node_spacing_um=2.0, probe="uniaxial",
                          n_steps=3000, device="cpu", rng=None) -> float:
    """One-shot calibration of the continuum-lattice bond stiffness: build with k=1, measure E-per-k via
    the affine probe, then return k = E_target / (E_per_k). This is the constitutive inverse (solving for
    the material property), NOT a fit to an outcome — the target E is the literature input."""
    from ffn_sim.ff.ecm_library import build_continuum_ecm
    ecm1 = build_continuum_ecm(spec, box_lo, box_hi, dim=dim, node_spacing_um=node_spacing_um,
                               k_bond_pN_um=1.0, rng=rng)
    if probe == "shear":
        E_per_k = shear_modulus(ecm1, n_steps=n_steps, device=device)["G_Pa"]
    else:
        E_per_k = uniaxial_modulus(ecm1, n_steps=n_steps, device=device)["E_Pa"]
    return spec.E_gel_Pa / max(E_per_k, 1e-12)


__all__ = ["shear_modulus", "uniaxial_modulus", "anisotropy", "indentation_modulus",
           "hertz_E_from_force", "calibrate_continuum_k", "shear_energy", "strain_stiffening",
           "stress_relaxation", "ecm_material_stress"]
