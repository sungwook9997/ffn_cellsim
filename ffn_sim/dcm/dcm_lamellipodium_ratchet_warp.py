"""GPU geometry for the lamellipodium ratchet — the per-batch leading-node build.

Phase C migration · M2 — GPU-port of the DETERMINISTIC geometry inside
:meth:`dcm_lamellipodium_host.LamellipodiumHost.update`, so the active-junction
stack becomes 100% GPU-resident. The filopodia PROBE was ported the same way
(:mod:`dcm_filopodia_probe_warp`); this module follows that STRUCTURE.

What `LamellipodiumHost.update(P)` does each batch tick splits cleanly into:

  DETERMINISTIC (pos-only, pure geometry — ported here)
    * per-cell centroid ``cc`` (mean of a cell's node positions, by ``cof``);
    * rim test ``(cc_z - z0) <= contact_band·R`` (host fallback to the lowest-third
      when the band is empty is NON-geometric and stays on the host — see below);
    * spheroid in-plane centroid ``sph`` (mean of all rim-cell node positions);
    * per-leading-node geometry: outward dir ``out = normalize(cc - sph)`` with z
      dropped, node projection ``proj = (r-cc)·out``, the basal+outward leading
      predicate, and the leading node's own position ``rp = P[ln]``.

  PROBABILISTIC + STATEFUL (RNG-driven actin SEED/ADVANCE over the append-only
  anchor pool — stays on the host, see :meth:`LamellipodiumHost._ratchet_from_geo`).
    The ratchet is sequential per cell (each seed mutates the actin view the next
    node reads) and append-only across ticks; it is cheap and not parallelism-shaped.

PATH CHOSEN: (a) **GPU-compute the geometry, feed it to the host ratchet.** This
removes the per-batch full GPU→CPU copy of all N node positions (``pos_d.numpy()``):
the host ratchet then runs only over the small per-leading-node geometry the GPU
returns (O(rim leading nodes), already on host as the kernel result). The RNG path
is kept bit-identical to the host (same ``np.random.default_rng`` stream), so the
opt-in GPU ratchet is byte-identical to the host path on its outputs — the only
thing that changed is WHERE the deterministic geometry was computed.

Why not (b) GPU-rng advance: it would need a parallel append into a shared pool
(atomic high-water mark) and a per-cell RNG substream, and could at best be
STATISTICALLY equivalent (not bit-exact) to the numpy stream — for a cheap,
sequential, append-only op that buys nothing, so (a) is strictly better here.

BIT-PARITY: the kernels use the SAME f64 formulas as the host numpy
(centroid = sum/count, ``out`` normalized after zeroing z, ``proj`` = dot in-plane,
basal ``|z - z_basal| <= basal_band``, lead ``proj >= lead_frac·max(rel_xy,1e-18)``).
The leading-node ITERATION ORDER is preserved (cells in ``rim_cells`` order, nodes
in ascending node-id order) so the host ratchet — which is order-sensitive — sees
an identical leading-node stream.

Positions are f64 ``vec3d``; the rim/spheroid reductions run in f64 (no f32 grid is
needed — there is no neighbour search here, only per-cell reductions over ``cof``).

SI units throughout.
"""

from __future__ import annotations

import numpy as np
import warp as wp

wp.init()


# ---------------------------------------------------------------------------
# Kernels — per-cell reductions + per-leading-node geometry, all f64
# ---------------------------------------------------------------------------
@wp.kernel
def cell_centroid_accum_kernel(
    pos: wp.array(dtype=wp.vec3d),          # (N,) f64 node positions
    cof: wp.array(dtype=wp.int32),          # (N,) owner cell of each node (−1 dormant)
    csum: wp.array(dtype=wp.vec3d),         # (n_cells,) OUT accumulator (zeroed by caller)
    ccnt: wp.array(dtype=wp.float64),       # (n_cells,) OUT count accumulator (zeroed)
):
    """Per-cell centroid via atomic reduction: ``csum[c] += pos[i]``, ``ccnt[c] += 1``
    for every live node ``i`` (``cof[i] >= 0``). Caller divides ``csum/ccnt`` to get the
    centroid — identical to the host ``np.add.at`` accumulation (sum then divide)."""
    i = wp.tid()
    c = cof[i]
    if c >= wp.int32(0):
        wp.atomic_add(csum, c, pos[i])
        wp.atomic_add(ccnt, c, wp.float64(1.0))


@wp.kernel
def cell_centroid_finalize_kernel(
    csum: wp.array(dtype=wp.vec3d),
    ccnt: wp.array(dtype=wp.float64),
    cc: wp.array(dtype=wp.vec3d),           # (n_cells,) OUT centroid
):
    """``cc[c] = csum[c] / max(ccnt[c], 1)`` — matches the host ``np.where(cnt>0, cnt, 1)``
    guard so an empty cell yields a zero centroid (count clamped to 1, sum is 0)."""
    c = wp.tid()
    n = ccnt[c]
    if n <= wp.float64(0.0):
        n = wp.float64(1.0)
    s = csum[c]
    cc[c] = wp.vec3d(s[0] / n, s[1] / n, s[2] / n)


@wp.kernel
def rim_test_kernel(
    cc: wp.array(dtype=wp.vec3d),           # (n_cells,) cell centroids
    ccnt: wp.array(dtype=wp.float64),       # (n_cells,) node counts (0 → cell absent)
    z0: wp.float64,
    band: wp.float64,                       # = contact_band · R
    is_rim: wp.array(dtype=wp.int32),       # (n_cells,) OUT 1 if rim else 0
):
    """Rim cell ⇔ ``(cc_z − z0) <= contact_band·R`` (host ``detect_rim_cells`` primary
    branch). Cells with no nodes (``ccnt==0``) are never rim. The host FALLBACK (lowest
    third by z when NO cell is within the band) is non-geometric and is applied on the
    host after this kernel, so the band-empty case stays bit-identical."""
    c = wp.tid()
    if ccnt[c] > wp.float64(0.0) and (cc[c][2] - z0) <= band:
        is_rim[c] = wp.int32(1)
    else:
        is_rim[c] = wp.int32(0)


@wp.kernel
def sph_centroid_accum_kernel(
    pos: wp.array(dtype=wp.vec3d),          # (N,) node positions
    cof: wp.array(dtype=wp.int32),          # (N,) owner cell
    cell_is_rim: wp.array(dtype=wp.int32),  # (n_cells,) 1 if cell is in the rim set
    ssum: wp.array(dtype=wp.vec3d),         # (1,) OUT sum of rim-cell node positions
    scnt: wp.array(dtype=wp.float64),       # (1,) OUT count of rim-cell nodes
):
    """Spheroid in-plane centroid accumulator: sum over every node whose owner cell is
    in the rim set (``rim_node_mask = np.isin(cof, rim_cells)``), then ``mean`` on the
    host. Matches ``P[rim_node_mask].mean(axis=0)`` (full 3D sum; the host only uses
    x,y but computes the 3D mean — kept identical)."""
    i = wp.tid()
    c = cof[i]
    if c >= wp.int32(0) and cell_is_rim[c] == wp.int32(1):
        wp.atomic_add(ssum, 0, pos[i])
        wp.atomic_add(scnt, 0, wp.float64(1.0))


@wp.kernel
def gather_f64_kernel(
    src: wp.array(dtype=wp.float64),        # (N,) per-node f64 value
    idx: wp.array(dtype=wp.int32),          # (L,) leading node ids
    out: wp.array(dtype=wp.float64),        # (L,) OUT gathered values
):
    """Device gather ``out[t] = src[idx[t]]`` for a per-node f64 geometry channel."""
    t = wp.tid()
    out[t] = src[idx[t]]


@wp.kernel
def gather_rp_kernel(
    pos: wp.array(dtype=wp.vec3d),          # (N,) node positions
    idx: wp.array(dtype=wp.int32),          # (L,) leading node ids
    out: wp.array(dtype=wp.vec3d),          # (L,) OUT leading node positions
):
    """Device gather ``out[t] = pos[idx[t]]`` so only the L leading-node positions cross
    the bus (not all N), keeping the host ratchet input O(L)."""
    t = wp.tid()
    out[t] = pos[idx[t]]


@wp.kernel
def leading_node_geo_kernel(
    pos: wp.array(dtype=wp.vec3d),          # (N,) node positions
    cof: wp.array(dtype=wp.int32),          # (N,) owner cell
    cc: wp.array(dtype=wp.vec3d),           # (n_cells,) per-cell centroid
    cell_is_rim: wp.array(dtype=wp.int32),  # (n_cells,) 1 if rim
    sx: wp.float64, sy: wp.float64,         # spheroid in-plane centroid
    z_basal: wp.float64, basal_band: wp.float64,
    lead_frac: wp.float64,
    # OUT, one slot per node i (compacted on the host in node-id order):
    is_lead: wp.array(dtype=wp.int32),      # 1 if node i is a leading node
    out_ox: wp.array(dtype=wp.float64),     # outward dir x of node i's cell
    out_oy: wp.array(dtype=wp.float64),     # outward dir y
    out_proj: wp.array(dtype=wp.float64),   # node i's outward projection from cc
    out_ccx: wp.array(dtype=wp.float64),    # node i's cell centroid x
    out_ccy: wp.array(dtype=wp.float64),    # node i's cell centroid y
):
    """Per-node leading-node predicate + geometry — exactly the host inner test.

    For each live node ``i`` of a rim cell: outward dir ``out = normalize((cc-sph) with
    z=0)`` (skip if ``|cc-sph|_xy < 1e-18``), projection ``proj = (r-cc)·out``,
    ``rel_xy = hypot(r-cc)``, basal ``|r_z - z_basal| <= basal_band``, and the leading
    predicate ``basal AND proj >= lead_frac·max(rel_xy, 1e-18)``. Writes per-node so the
    host can compact in node-id order — preserving the host loop's leading-node order
    (cells in rim order is enforced by the host's compaction key, nodes ascending)."""
    i = wp.tid()
    is_lead[i] = wp.int32(0)
    out_ox[i] = wp.float64(0.0)
    out_oy[i] = wp.float64(0.0)
    out_proj[i] = wp.float64(0.0)
    out_ccx[i] = wp.float64(0.0)
    out_ccy[i] = wp.float64(0.0)
    c = cof[i]
    if c < wp.int32(0):
        return
    if cell_is_rim[c] != wp.int32(1):
        return
    ccc = cc[c]
    ox = ccc[0] - sx
    oy = ccc[1] - sy
    on = wp.sqrt(ox * ox + oy * oy)
    if on < wp.float64(1.0e-18):
        return
    ox = ox / on
    oy = oy / on
    r = pos[i]
    relx = r[0] - ccc[0]
    rely = r[1] - ccc[1]
    proj = relx * ox + rely * oy
    rel_xy = wp.sqrt(relx * relx + rely * rely)
    basal = wp.abs(r[2] - z_basal) <= basal_band
    thresh = lead_frac * wp.max(rel_xy, wp.float64(1.0e-18))
    if basal and proj >= thresh:
        is_lead[i] = wp.int32(1)
        out_ox[i] = ox
        out_oy[i] = oy
        out_proj[i] = proj
        out_ccx[i] = ccc[0]
        out_ccy[i] = ccc[1]


# ---------------------------------------------------------------------------
# Persistent device-array cache (re-used / grown across batch ticks)
# ---------------------------------------------------------------------------
_GPU: dict = {}


def _gather_f64(src, idx_d, L: int, device):
    """Gather ``src[idx_d]`` into a fresh (L,) host array via the device gather kernel."""
    out_d = wp.zeros(L, dtype=wp.float64, device=device)
    wp.launch(gather_f64_kernel, dim=L, inputs=[src, idx_d, out_d], device=device)
    return out_d.numpy()


def _dev_cache(device, N: int, n_cells: int):
    """Per-device scratch bundle keyed by node/cell counts (geometric growth)."""
    key = str(device)
    d = _GPU.get(key)
    if d is None or d["N"] < N or d["n_cells"] < n_cells:
        Nn = max(N, d["N"] if d else 0)
        Cc = max(n_cells, d["n_cells"] if d else 0)
        d = {
            "N": Nn, "n_cells": Cc,
            "pos": wp.zeros(Nn, dtype=wp.vec3d, device=device),
            "cof": wp.zeros(Nn, dtype=wp.int32, device=device),
            "csum": wp.zeros(Cc, dtype=wp.vec3d, device=device),
            "ccnt": wp.zeros(Cc, dtype=wp.float64, device=device),
            "cc": wp.zeros(Cc, dtype=wp.vec3d, device=device),
            "is_rim": wp.zeros(Cc, dtype=wp.int32, device=device),
            "ssum": wp.zeros(1, dtype=wp.vec3d, device=device),
            "scnt": wp.zeros(1, dtype=wp.float64, device=device),
            "is_lead": wp.zeros(Nn, dtype=wp.int32, device=device),
            "out_ox": wp.zeros(Nn, dtype=wp.float64, device=device),
            "out_oy": wp.zeros(Nn, dtype=wp.float64, device=device),
            "out_proj": wp.zeros(Nn, dtype=wp.float64, device=device),
            "out_ccx": wp.zeros(Nn, dtype=wp.float64, device=device),
            "out_ccy": wp.zeros(Nn, dtype=wp.float64, device=device),
        }
        _GPU[key] = d
    return d


def compute_lamellipodium_geometry_gpu(
    *, pos: np.ndarray, cof: np.ndarray, n_cells: int,
    z0: float, R: float, contact_band: float,
    z_basal: float, basal_band: float, lead_frac: float,
    device: str = "cpu", pos_d=None,
):
    """GPU geometry for one lamellipodium ratchet tick — the deterministic half of
    :meth:`LamellipodiumHost.update`.

    Computes, fully on device (no per-node ``.numpy()`` of positions): per-cell
    centroids, the rim-cell set, the spheroid in-plane centroid, and the per-leading-node
    geometry (outward dir, projection, cell centroid). Returns small host arrays the host
    ratchet then consumes — the only data crossing the bus is O(n_cells + leading nodes),
    not O(N).

    Args:
        pos: ``(N,3)`` f64 node positions. Only used to upload when ``pos_d`` is None;
            when ``pos_d`` is given (a live device array) the positions never leave the GPU.
        cof: ``(N,)`` int owner cell of each node (−1 dormant).
        n_cells: number of cells.
        z0, R, contact_band: rim test ``(cc_z − z0) <= contact_band·R``.
        z_basal, basal_band, lead_frac: leading-node basal + outward predicate params.
        device: warp device ("cpu" or cuda).
        pos_d: optional pre-existing ``wp.array(vec3d)`` of positions ON device — when
            supplied, NO position copy happens (the GPU-resident fast path). ``cof`` is
            still small and is uploaded each tick.

    Returns:
        dict with keys:
          ``rim_cells`` (int64 ascending) — rim-cell ids (host fallback applied here),
          ``cell_cc`` ``(n_cells,3)`` f64 — per-cell centroids,
          ``sph`` ``(3,)`` f64 — spheroid centroid (or None if no rim node),
          ``lead_node_id`` (int64, node-id ascending) — leading node ids,
          ``lead_rp`` ``(L,3)`` f64 — leading nodes' own positions,
          ``lead_ccx/ccy`` ``(L,)`` f64, ``lead_ox/oy`` ``(L,)`` f64, ``lead_proj`` ``(L,)``.
        ``lead_*`` are ordered first by cell (rim order) then by node id, matching the
        host loop's leading-node emission order.
    """
    cof = np.ascontiguousarray(cof, dtype=np.int64).reshape(-1)
    N = cof.shape[0]
    band = float(contact_band) * float(R)

    d = _dev_cache(device, N, n_cells)

    # positions: reuse the live device array if given, else upload once
    if pos_d is not None:
        posarr = pos_d
    else:
        pos = np.ascontiguousarray(pos, dtype=np.float64).reshape(-1, 3)
        d["pos"].assign(pos)
        posarr = d["pos"]
    d["cof"].assign(cof.astype(np.int32))

    # --- per-cell centroid -------------------------------------------------
    d["csum"].zero_()
    d["ccnt"].zero_()
    wp.launch(cell_centroid_accum_kernel, dim=N,
              inputs=[posarr, d["cof"], d["csum"], d["ccnt"]], device=device)
    wp.launch(cell_centroid_finalize_kernel, dim=n_cells,
              inputs=[d["csum"], d["ccnt"], d["cc"]], device=device)

    # --- rim test (+ host fallback if band-empty) --------------------------
    wp.launch(rim_test_kernel, dim=n_cells,
              inputs=[d["cc"], d["ccnt"], wp.float64(z0), wp.float64(band), d["is_rim"]],
              device=device)
    cell_cc = d["cc"].numpy()[:n_cells].copy()
    is_rim_h = d["is_rim"].numpy()[:n_cells]
    rim_cells = np.flatnonzero(is_rim_h == 1).astype(np.int64)
    if rim_cells.size == 0:
        # host fallback: lowest third by z over PRESENT cells (ccnt>0), exactly
        # detect_rim_cells' ``np.argsort(centers[:,2])[: max(1, n//3)]``. The host
        # builds centers for ALL n_cells (absent cells have centroid 0 from the guard),
        # so we argsort the same n_cells-length z just like the host.
        order = np.argsort(cell_cc[:, 2])
        rim_cells = order[: max(1, n_cells // 3)].astype(np.int64)
        # reflect the fallback choice on device so the spheroid + leading kernels use it
        is_rim_h = np.zeros(n_cells, dtype=np.int32)
        is_rim_h[rim_cells] = 1
        d["is_rim"].assign(is_rim_h)

    # --- spheroid in-plane centroid over rim-cell nodes --------------------
    d["ssum"].zero_()
    d["scnt"].zero_()
    wp.launch(sph_centroid_accum_kernel, dim=N,
              inputs=[posarr, d["cof"], d["is_rim"], d["ssum"], d["scnt"]], device=device)
    scnt = float(d["scnt"].numpy()[0])
    if scnt <= 0.0:
        return {"rim_cells": rim_cells, "cell_cc": cell_cc, "sph": None,
                "lead_node_id": np.empty(0, np.int64), "lead_cell": np.empty(0, np.int64),
                "lead_rp": np.empty((0, 3)),
                "lead_ccx": np.empty(0), "lead_ccy": np.empty(0),
                "lead_ox": np.empty(0), "lead_oy": np.empty(0), "lead_proj": np.empty(0)}
    ssum = d["ssum"].numpy()[0]
    sph = np.asarray(ssum, dtype=np.float64) / scnt
    sx, sy = float(sph[0]), float(sph[1])

    # --- per-leading-node geometry -----------------------------------------
    wp.launch(leading_node_geo_kernel, dim=N,
              inputs=[posarr, d["cof"], d["cc"], d["is_rim"],
                      wp.float64(sx), wp.float64(sy), wp.float64(z_basal),
                      wp.float64(basal_band), wp.float64(lead_frac),
                      d["is_lead"], d["out_ox"], d["out_oy"], d["out_proj"],
                      d["out_ccx"], d["out_ccy"]],
              device=device)

    is_lead = d["is_lead"].numpy()[:N]
    lead_mask = is_lead == 1
    lead_node_id_natural = np.flatnonzero(lead_mask).astype(np.int64)  # ascending node id
    if lead_node_id_natural.size == 0:
        return {"rim_cells": rim_cells, "cell_cc": cell_cc, "sph": sph,
                "lead_node_id": np.empty(0, np.int64), "lead_cell": np.empty(0, np.int64),
                "lead_rp": np.empty((0, 3)),
                "lead_ccx": np.empty(0), "lead_ccy": np.empty(0),
                "lead_ox": np.empty(0), "lead_oy": np.empty(0), "lead_proj": np.empty(0)}

    # Re-order leading nodes by (cell rank in rim_cells, node id) to match the host's
    # emission order (host loops cells in rim_cells order, nodes ascending). Within a
    # cell, np.where already yields ascending node ids, so a stable sort on the cell rank
    # is sufficient.
    cof_lead = cof[lead_node_id_natural]
    rank = np.full(n_cells, np.iinfo(np.int64).max, dtype=np.int64)
    rank[rim_cells] = np.arange(rim_cells.size, dtype=np.int64)
    cell_rank = rank[cof_lead]
    order = np.argsort(cell_rank, kind="stable")
    lead_node_id = lead_node_id_natural[order]
    L = lead_node_id.size

    # Gather the L leading-node geometry slots + positions ON DEVICE, so only L-sized
    # arrays cross the bus (the per-node geometry stays GPU-resident). ``lead_node_id``
    # is the host compaction order; we upload it and gather each per-node device array.
    lead_idx_d = wp.array(lead_node_id.astype(np.int32), dtype=wp.int32, device=device)
    rp_d = wp.zeros(L, dtype=wp.vec3d, device=device)
    wp.launch(gather_rp_kernel, dim=L, inputs=[posarr, lead_idx_d, rp_d], device=device)
    ox = _gather_f64(d["out_ox"], lead_idx_d, L, device)
    oy = _gather_f64(d["out_oy"], lead_idx_d, L, device)
    proj = _gather_f64(d["out_proj"], lead_idx_d, L, device)
    ccx = _gather_f64(d["out_ccx"], lead_idx_d, L, device)
    ccy = _gather_f64(d["out_ccy"], lead_idx_d, L, device)
    lead_rp = np.ascontiguousarray(rp_d.numpy(), dtype=np.float64)

    return {
        "rim_cells": rim_cells,
        "cell_cc": cell_cc,
        "sph": sph,
        "lead_node_id": lead_node_id,
        "lead_cell": np.ascontiguousarray(cof[lead_node_id], dtype=np.int64),
        "lead_rp": lead_rp,
        "lead_ccx": np.ascontiguousarray(ccx, dtype=np.float64),
        "lead_ccy": np.ascontiguousarray(ccy, dtype=np.float64),
        "lead_ox": np.ascontiguousarray(ox, dtype=np.float64),
        "lead_oy": np.ascontiguousarray(oy, dtype=np.float64),
        "lead_proj": np.ascontiguousarray(proj, dtype=np.float64),
    }
