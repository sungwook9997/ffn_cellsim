"""FF motility — physical-η overdamped crawl on a substrate (piece-4 closed loop, real-time dynamics).

The going-forward answer to the PI (2026-07-06): "protrusion이 세포를 변형시키며 움직이는 것은 보지 못함 …
막이랑 필라멘트 다같이 모양 변형되는 것도 못 봤음." This module makes the whole cell CRAWL — deforming and
translocating in **real physical time**, membrane + cortex + filaments + nucleus co-moving — instead of the
prior quasi-static macro×relax split (equilibrate-between-ticks) that only showed a thin filopodial appendage.

Two changes vs the existing GPU whole-cell path:

  (B) REAL η-dynamics — replaces the numerical pseudo-time step ``x += (0.1/kmax)·F`` with the physical
      overdamped law ``x += (dt/γ_i)·F_i``, where ``γ_i`` is the NF2007 §5.2 per-node cytoplasm drag
      (:func:`ffn_sim.ff.units.fiber_point_drag`, η = 65.9 Pa·s MCF7 cytoplasm, Dessard 2024). The trajectory
      is now physical seconds, so the crawl SPEED is a measurable observable comparable to literature.

  (C) piece-4 closed loop — a polarized leading-edge protrusive stress (lamellipodial Brownian-ratchet push,
      Mogilner-Oster; KU-3.6) drives the FRONT cortex/membrane forward; the emergent opposing load (membrane
      tension + elastic resistance the front pushes against) throttles the ratchet ``v(f)=v0·e^{-fδ/kT}``;
      basal FA integrin catch-slip clutches (KU-2.5) grip the substrate → the protrusion is converted to net
      forward COM translocation by TRACTION. Newton's 3rd law: polymerization alone is an internal force, so
      the decisive control is that **without clutches the net COM drift is ≈ 0** (protrusion deforms the front
      in place); only with substrate traction does the cell body move.

Single physical-time loop: bending + crosslink + myosin + turgor + membrane(inward γ_mem) + nucleus + substrate
+ FA clutch + leading-edge protrusion ALL tick at the same ``dt``. The slow incompressible turgor/volume mode
is refreshed periodically (Guo closure) and held between refreshes — a legitimate multiscale treatment of a
genuinely fast mode, NOT the mechanical-deformation-hiding quasi-statics the PI objected to.

KU refs: KU-2.4 (Chan-Odde motor-clutch), KU-2.5 (integrin catch-slip), KU-3.6 (Mogilner-Oster ratchet),
KU-3.12 (retrograde flow). Reuses the validated kernels in ``network_warp`` / ``fa_clutch_warp``.
"""
from __future__ import annotations

import numpy as np
import warp as wp

from ffn_sim.ff import units as U
from ffn_sim.ff.units import ETA_CYTOPLASM, fiber_point_drag


@wp.kernel
def axpy_physical_kernel(pos: wp.array(dtype=wp.vec3d), dt: wp.float64,
                         gamma: wp.array(dtype=wp.float64), F: wp.array(dtype=wp.vec3d)):
    """Physical overdamped step ``x_i += (dt/γ_i)·F_i`` — γ_i is the per-node cytoplasm drag [pN·s/µm]
    (NF2007 §5.2, η = 65.9 Pa·s). Unlike ``axpy_kernel`` (whose scalar ``step = 0.1/kmax`` is a numerical
    gradient-descent step with γ≡1), the step here is PHYSICAL: dt is real seconds and the resulting velocity
    v_i = F_i/γ_i is the real overdamped velocity, so trajectory time and crawl speed are physical."""
    i = wp.tid()
    pos[i] = pos[i] + (dt / gamma[i]) * F[i]


@wp.kernel
def leading_edge_push_kernel(pos: wp.array(dtype=wp.vec3d), com: wp.vec3d, phat: wp.vec3d,
                             front_cos_R: wp.float64, f_pro: wp.float64, delta: wp.float64,
                             kT: wp.float64, force: wp.array(dtype=wp.vec3d),
                             total: wp.array(dtype=wp.float64)):
    """Lamellipodial leading-edge protrusion (front push half of the INTERNAL force pair). Cortex nodes on the
    FRONT cap (projection onto the polarity axis ``phat`` exceeds ``front_cos_R``) get a protrusive force
    ``f_pro·v(load)·phat`` along the crawl axis, and the total pushed magnitude is accumulated into ``total[0]``
    (consumed by :func:`protrusion_reaction_kernel` as the equal-and-opposite retrograde reaction, so the net
    protrusion force on the cell is ZERO — polymerization is INTERNAL).

    The Brownian-ratchet factor ``v = e^{-f_comp·δ/kT}`` (Mogilner-Oster, KU-3.6) throttles the push by the
    EMERGENT opposing load ``f_comp = -F·phat`` (membrane tension + elastic resistance already accumulated in
    ``force`` this step) — so this MUST be launched last, after all structural forces are assembled. High
    counter-load → stalled push (piece-4 emergent closed loop); zero load → full ``f_pro``."""
    i = wp.tid()
    d = pos[i] - com
    proj = wp.dot(d, phat)
    if proj > front_cos_R:
        f_comp = wp.max(-wp.dot(force[i], phat), wp.float64(0.0))   # emergent opposing (compressive) load along +phat
        v = wp.exp(-f_comp * delta / kT)                 # Mogilner-Oster ratchet throttle ∈ (0,1]
        fmag = f_pro * v
        wp.atomic_add(force, i, fmag * phat)
        wp.atomic_add(total, 0, fmag)


@wp.kernel
def protrusion_reaction_kernel(phat: wp.vec3d, total: wp.array(dtype=wp.float64), n_cortex: wp.float64,
                               force: wp.array(dtype=wp.vec3d)):
    """Retrograde reaction (the other half of the pair): the total leading-edge push ``total[0]`` is subtracted
    uniformly over the ``n_cortex`` cortex nodes, ``F_i −= (total[0]/n_cortex)·phat``. This is Newton's 3rd law —
    the actin polymerizing against the front membrane pushes the network BACKWARD (retrograde flow). Net
    protrusion force on the cell = 0, so protrusion alone CANNOT translocate the COM; only when the retrograde
    reaction on the basal network is resisted by substrate clutches (traction) does the cell move forward. Launch
    over cortex nodes only, AFTER ``leading_edge_push_kernel``."""
    i = wp.tid()
    wp.atomic_add(force, i, -(total[0] / n_cortex) * phat)


@wp.kernel
def gravity_kernel(fz_node: wp.float64, force: wp.array(dtype=wp.vec3d)):
    """Net sedimentation body force (gravity − buoyancy): every cortex node gets a downward z-force
    ``fz_node = −Δρ·g·v_node`` (Δρ = ρ_cell − ρ_medium > 0 ⇒ the cell sinks toward the dish). This is REAL
    physics (a derived body force ≈ 1 pN/cell at MCF7 size), NOT a tuned proxy — it rests the cell on the
    substrate. Ported from the DCM ``gravity_body_force_kernel`` (dcm_neighbor_warp); engine-agnostic."""
    i = wp.tid()
    f = force[i]
    force[i] = wp.vec3d(f[0], f[1], f[2] + fz_node)


@wp.kernel
def cortex_volume_kernel(pos: wp.array(dtype=wp.vec3d), faces: wp.array(dtype=wp.int32, ndim=2),
                         centre: wp.vec3d, vol_out: wp.array(dtype=wp.float64)):
    """Enclosed cortex volume via the divergence theorem: sum the signed tetra volumes {centre, v0, v1, v2}
    over a FIXED surface triangulation (the frame-0 hull faces; nodes move, faces stay). ``vol_out[0]`` must be
    zeroed before launch; take ``abs`` on read. Lets the osmotic volume feedback be applied EVERY step (no
    refresh lag) so the stiff incompressible-cytoplasm constraint is stable — this is what stops both the
    balloon (one-sided, lagged) and the collapse (bidirectional, lagged)."""
    m = wp.tid()
    a = pos[faces[m, 0]] - centre
    b = pos[faces[m, 1]] - centre
    c = pos[faces[m, 2]] - centre
    wp.atomic_add(vol_out, 0, wp.dot(a, wp.cross(b, c)) / wp.float64(6.0))


@wp.kernel
def xl_turnover_kernel(pos: wp.array(dtype=wp.vec3d), xl_ij: wp.array(dtype=wp.int32, ndim=2),
                       xl_rest: wp.array(dtype=wp.float64), frac: wp.float64):
    """Cortical crosslink turnover (α-actinin unbinds/rebinds, k_off): a rebinding crosslink is force-free at
    its CURRENT length, so its rest length Maxwell-relaxes toward the current bond length,
    ``rest ← rest + frac·(L − rest)`` with ``frac = 1 − e^{−k_off·Δt}``. This fluidizes the cortex (a crawling
    cell's cortex FLOWS — it is not a rigid Ferrer-stiff shell), letting the leading-edge protrusion deform the
    cell body over time. On-device (no host round-trip). Mechanistic: crosslinks are dynamic bonds, not permanent
    springs."""
    k = wp.tid()
    i = xl_ij[k, 0]; j = xl_ij[k, 1]
    L = wp.length(pos[i] - pos[j])
    xl_rest[k] = xl_rest[k] + frac * (L - xl_rest[k])


@wp.kernel
def spreading_push_kernel(pos: wp.array(dtype=wp.vec3d), cx: wp.float64, cy: wp.float64, z_sub: wp.float64,
                          h_basal: wp.float64, r_min: wp.float64, f_pro: wp.float64, delta: wp.float64,
                          kT: wp.float64, force: wp.array(dtype=wp.vec3d), total: wp.array(dtype=wp.float64)):
    """EMERGENT spreading: basal peripheral (rim) cortex nodes get a RADIALLY-OUTWARD (in the substrate plane)
    Brownian-ratchet push — the lamellipodial actin polymerizing against the membrane all around the contact
    edge, widening the footprint. Same Mogilner-Oster mechanism as the leading-edge crawl push, isotropic
    instead of polarized, so the cell SPREADS (contact area grows; at constant volume the cell flattens). NO
    wetting energy / adhesion-energy proxy — spreading is emergent from actin polymerization + (basal clutch)
    adhesion. ``total`` (2 floats: x,y) accumulates the pushed vector for the retrograde reaction (Newton)."""
    i = wp.tid()
    p = pos[i]
    if (p[2] - z_sub) < h_basal:
        rx = p[0] - cx
        ry = p[1] - cy
        r = wp.sqrt(rx * rx + ry * ry)
        if r > r_min:
            ux = rx / r
            uy = ry / r
            f_comp = wp.max(-(force[i][0] * ux + force[i][1] * uy), wp.float64(0.0))   # opposing radial load
            fmag = f_pro * wp.exp(-f_comp * delta / kT)                                 # ratchet throttle
            wp.atomic_add(force, i, wp.vec3d(fmag * ux, fmag * uy, wp.float64(0.0)))
            wp.atomic_add(total, 0, fmag * ux)
            wp.atomic_add(total, 1, fmag * uy)


@wp.kernel
def spreading_reaction_kernel(total: wp.array(dtype=wp.float64), n_cortex: wp.float64,
                              force: wp.array(dtype=wp.vec3d)):
    """Retrograde reaction for spreading (Newton): the net pushed vector ``(total[0], total[1])`` is subtracted
    uniformly over the cortex nodes, so the spreading push exerts ZERO net force on the cell (it widens the
    footprint, it does not translocate the COM). For a symmetric spread total≈0; this cancels any asymmetry."""
    i = wp.tid()
    wp.atomic_add(force, i, wp.vec3d(-total[0] / n_cortex, -total[1] / n_cortex, wp.float64(0.0)))


def physical_node_gammas(net, n_cortex: int, n_nuc: int, *, eta: float = ETA_CYTOPLASM,
                         bead_radius_um: float = 0.1) -> np.ndarray:
    """Per-node overdamped drag γ [pN·s/µm] for ``pos_all = [cortex nodes ; nucleus beads]``.

    Cortex/filament nodes: NF2007 §5.2 fiber log-drag :func:`units.fiber_point_drag` (γ = 1/((p+1)·μ),
    μ = log(L_h/δ)/(3πηL)), one value per fiber broadcast to its nodes — physically, a longer/more-resolved
    filament couples to more cytoplasm. Nucleus beads: Stokes ``6πη·a`` (a = ``bead_radius_um``). All from
    the single physiological viscosity η — no free drag parameter."""
    off = np.asarray(net.fiber_offsets, dtype=np.int64)
    nfib = len(off) - 1
    seg_per = np.diff(off) - 1
    seg_off = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int64)
    sr = np.asarray(net.seg_rest, dtype=np.float64)
    seg_mean = float(sr.mean()) if sr.size else 1.0
    g = np.empty(n_cortex + n_nuc, dtype=np.float64)
    for f in range(nfib):
        p = int(seg_per[f])
        if p > 0:
            L = float(sr[seg_off[f]:seg_off[f] + p].sum())
        else:
            L = seg_mean
        g[off[f]:off[f + 1]] = fiber_point_drag(max(L, 1e-3), max(p, 1), eta=eta)
    if n_nuc:
        g[n_cortex:] = 6.0 * np.pi * eta * bead_radius_um            # Stokes bead drag (nucleus cloud)
    return g


def volume_gradient(pcx: np.ndarray, faces: np.ndarray, centre: np.ndarray) -> np.ndarray:
    """∂V/∂x for the enclosed volume V = (1/6)·Σ_faces (a−c₀)·((b−c₀)×(c−c₀)) over the oriented surface
    triangulation ``faces`` (a,b,c = the 3 vertex ids, c₀ = ``centre``). Returns (Nc,3). This is the EXACT
    turgor-force direction (``f_i = ΔP·g_i``): net-zero over rigid translation (internal), outward for ΔP>0.
    Used so the osmotic volume constraint can be treated IMPLICITLY as the rank-1 stiffness ``k_vol·g·gᵀ``
    (g = flattened this) — the piece that lets the implicit crawl run at large dt without the volume blowing up."""
    g = np.zeros_like(pcx)
    a = pcx[faces[:, 0]] - centre
    b = pcx[faces[:, 1]] - centre
    c = pcx[faces[:, 2]] - centre
    np.add.at(g, faces[:, 0], np.cross(b, c) / 6.0)
    np.add.at(g, faces[:, 1], np.cross(c, a) / 6.0)
    np.add.at(g, faces[:, 2], np.cross(a, b) / 6.0)
    return g


def crawl_cfl_dt(gammas: np.ndarray, kmax: float, *, safety: float = 0.1) -> float:
    """Stable explicit-overdamped timestep: ``dt = safety·min_i(γ_i)/kmax`` so that ``(dt/γ)·k < safety`` for
    the stiffest force constant ``kmax`` (bending κ/seg³, crosslink, nucleus, clutch, substrate). The slow
    turgor/volume mode is EXCLUDED from ``kmax`` (held constant between periodic Guo refreshes — it is not a
    stiff per-step spring here). Returns physical seconds."""
    return float(safety) * float(np.min(gammas)) / float(kmax)


__all__ = ["axpy_physical_kernel", "leading_edge_push_kernel", "protrusion_reaction_kernel",
           "spreading_push_kernel", "spreading_reaction_kernel", "gravity_kernel", "cortex_volume_kernel",
           "xl_turnover_kernel", "volume_gradient", "physical_node_gammas", "crawl_cfl_dt"]
