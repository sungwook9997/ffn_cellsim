"""Warp port of the lamellipodial traction-tether force (Phase B, lamellipodium).

Ports ``cell.dcm_lamellipodium.LamellipodialTractionTether`` — the traction that
pulls each rim cell's LEADING basal membrane nodes toward the nearest OUTWARD
clutch-anchored actin bead (the emergent crawl traction that replaces a prescribed
body force). Per leading node: find the nearest basal ``actin_lamel`` bead that is
radially farther out (within ``tether_radius``), apply a force-capped harmonic pull
``F = k_tether·(r_actin − r_node)`` on the node (force only on the node; the
clutch-gripped actin is the fixed anchor).

Port split (mirrors B2 / DCM-turgor): the per-cell geometry (spheroid centroid,
per-rim-cell centroid + outward direction, leading+basal node selection) is a set
of cheap reductions computed in numpy EXACTLY as the reference; the parallel
per-leading-node nearest-actin SEARCH + capped force runs in a Warp kernel (one
thread per leading node). Nearest = first strict-minimum distance (matches
``np.argmin``). Each leading node belongs to one cell → own-row write, no atomics.
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()


@wp.kernel
def lamellipodium_tether_kernel(
    lead_rp: wp.array(dtype=wp.vec3d),       # (L,) leading-node positions
    lead_ccx: wp.array(dtype=wp.float64),    # (L,) cell centroid x
    lead_ccy: wp.array(dtype=wp.float64),    # (L,) cell centroid y
    lead_ox: wp.array(dtype=wp.float64),     # (L,) outward dir x
    lead_oy: wp.array(dtype=wp.float64),     # (L,) outward dir y
    lead_proj: wp.array(dtype=wp.float64),   # (L,) node outward projection
    lead_idx: wp.array(dtype=wp.int32),      # (L,) global node row
    actin: wp.array(dtype=wp.vec3d),         # (A,) basal actin positions
    n_actin: wp.int32,
    k_tether: wp.float64,
    force_cap: wp.float64,
    tether_radius: wp.float64,
    force: wp.array(dtype=wp.vec3d),         # (N,) out (own-row write)
):
    t = wp.tid()
    rp = lead_rp[t]
    ccx = lead_ccx[t]
    ccy = lead_ccy[t]
    ox = lead_ox[t]
    oy = lead_oy[t]
    npj = lead_proj[t]

    min_dist = wp.float64(1.0e300)
    jmin = wp.int32(-1)
    for a in range(n_actin):
        ap = actin[a]
        dx = ap[0] - rp[0]
        dy = ap[1] - rp[1]
        dz = ap[2] - rp[2]
        dist = wp.sqrt(dx * dx + dy * dy + dz * dz)
        act_proj = (ap[0] - ccx) * ox + (ap[1] - ccy) * oy
        if dist <= tether_radius and act_proj > npj:
            if dist < min_dist:
                min_dist = dist
                jmin = wp.int32(a)

    if jmin >= wp.int32(0):
        ap = actin[jmin]
        fvx = k_tether * (ap[0] - rp[0])
        fvy = k_tether * (ap[1] - rp[1])
        fvz = k_tether * (ap[2] - rp[2])
        fmag = wp.sqrt(fvx * fvx + fvy * fvy + fvz * fvz)
        if fmag > force_cap:
            s = force_cap / fmag
            fvx = fvx * s
            fvy = fvy * s
            fvz = fvz * s
        force[lead_idx[t]] = wp.vec3d(fvx, fvy, fvz)


def run_lamellipodium_tether_warp(
    *,
    pos: np.ndarray,            # (N, 3) tag-ordered
    typeid: np.ndarray,         # (N,)
    cell_of_node: np.ndarray,   # (n_mem,) cell of each membrane node
    rim_cells: np.ndarray,      # (R,) rim cell ids
    actin_typeid: int,
    mem_typeid: int,
    z_basal: float,
    basal_band: float,
    k_tether: float,
    force_cap: float,
    tether_radius: float,
    lead_frac: float = 0.0,
    device: str = "cpu",
) -> dict:
    """Lamellipodial traction-tether force (N,3) in Warp, matching
    dcm_lamellipodium.LamellipodialTractionTether (tag-ordered I/O)."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    tid = np.ascontiguousarray(typeid)
    cof = np.ascontiguousarray(cell_of_node, dtype=np.int64)
    N = pos.shape[0]
    F = np.zeros((N, 3), dtype=np.float64)

    # ---- host: per-cell geometry + leading-node selection (== the reference) ----
    actin_all = np.where(tid == actin_typeid)[0]
    if actin_all.size == 0:
        return {"force": F, "n_tethered": 0}
    actin_pos_all = pos[actin_all]
    actin_basal = np.abs(actin_pos_all[:, 2] - z_basal) <= basal_band
    actin_pos = actin_pos_all[actin_basal]
    if actin_pos.shape[0] == 0:
        return {"force": F, "n_tethered": 0}

    n_mem_block = cof.shape[0]
    mem_idx = np.where(tid == mem_typeid)[0]
    mem_idx = mem_idx[mem_idx < n_mem_block]
    if mem_idx.size == 0:
        return {"force": F, "n_tethered": 0}
    sph_centroid = pos[mem_idx].mean(axis=0)
    sph_centroid[2] = z_basal

    rp_l, ccx_l, ccy_l, ox_l, oy_l, proj_l, idx_l = [], [], [], [], [], [], []
    for c in rim_cells:
        c = int(c)
        cnodes = np.where(cof == c)[0]
        if cnodes.size == 0:
            continue
        cell_centroid = pos[cnodes].mean(axis=0)
        out = cell_centroid - sph_centroid
        out[2] = 0.0
        on = np.linalg.norm(out)
        if on < 1e-12:
            continue
        out = out / on
        npos = pos[cnodes]
        rel = npos - cell_centroid
        proj = rel[:, 0] * out[0] + rel[:, 1] * out[1]
        basal = np.abs(npos[:, 2] - z_basal) <= basal_band
        lead = (proj >= lead_frac * np.maximum(
            np.linalg.norm(rel[:, :2], axis=1), 1e-18)) & basal
        lead_nodes = cnodes[lead]
        lead_projs = proj[lead]
        for ln, npj in zip(lead_nodes, lead_projs):
            rp_l.append(pos[ln]); ccx_l.append(cell_centroid[0])
            ccy_l.append(cell_centroid[1]); ox_l.append(out[0]); oy_l.append(out[1])
            proj_l.append(float(npj)); idx_l.append(int(ln))

    L = len(idx_l)
    if L == 0:
        return {"force": F, "n_tethered": 0}

    # ---- Warp kernel: per-leading-node nearest-actin search + capped force ----
    lead_rp = wp.array(np.asarray(rp_l, dtype=np.float64), dtype=wp.vec3d, device=device)
    lead_ccx = wp.array(np.asarray(ccx_l, dtype=np.float64), dtype=wp.float64, device=device)
    lead_ccy = wp.array(np.asarray(ccy_l, dtype=np.float64), dtype=wp.float64, device=device)
    lead_ox = wp.array(np.asarray(ox_l, dtype=np.float64), dtype=wp.float64, device=device)
    lead_oy = wp.array(np.asarray(oy_l, dtype=np.float64), dtype=wp.float64, device=device)
    lead_proj = wp.array(np.asarray(proj_l, dtype=np.float64), dtype=wp.float64, device=device)
    lead_idx = wp.array(np.asarray(idx_l, dtype=np.int32), dtype=wp.int32, device=device)
    actin_d = wp.array(np.ascontiguousarray(actin_pos, dtype=np.float64),
                       dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)

    wp.launch(
        lamellipodium_tether_kernel, dim=L,
        inputs=[lead_rp, lead_ccx, lead_ccy, lead_ox, lead_oy, lead_proj, lead_idx,
                actin_d, wp.int32(actin_pos.shape[0]), wp.float64(k_tether),
                wp.float64(force_cap), wp.float64(tether_radius), force_d],
        device=device,
    )
    wp.synchronize_device(device)
    return {"force": force_d.numpy().astype(np.float64), "n_tethered": L}
