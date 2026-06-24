"""Parity test: GPU lamellipodium ratchet geometry vs the host numpy reference.

Builds a small multi-cell icosphere cluster sitting on a dish, makes two identical
:class:`LamellipodiumHost` instances (host path vs ``use_gpu_ratchet``), runs ONE
``update()`` from identical fresh state, and asserts:

  * rim-cell SET identical,
  * per-leading-node geometry (node ids / outward dir / projection / cc / rp) matches
    within 1e-9 (the GPU per-cell centroid is an atomic f64 reduction, so it matches the
    numpy ``mean`` to ULP-level ~1e-12, not bit-exact — well inside 1e-9),
  * the resulting DEVICE anchor arrays (the ones the tether kernel reads) match: the
    deterministic geometry exactly (to 1e-9) and, because the seed/advance RNG path is
    SHARED (same ``np.random.default_rng`` stream, fed pre-computed geometry), the actin
    pool is byte-identical too — so the advance COUNT/rate matches exactly, not just
    statistically. (We chose host-fed advance, path (a) — see the module docstring.)

Run:  conda activate ffn_sim && python -m ffn_sim.warp_port._parity_lamellipodium_ratchet
"""

from __future__ import annotations

import numpy as np

from ffn_sim.cell.dcm import icosphere_mesh
from ffn_sim.warp_port.dcm_lamellipodium_host import LamellipodiumHost, LamelParams


def _build_cluster(n_cells: int = 5, R: float = 7.5e-6, subdiv: int = 1, seed: int = 3):
    """A small cluster of icospheres placed near a dish plane z0, with light jitter so
    the rim test + leading-node selection are non-degenerate."""
    verts0, _edges, _tris = icosphere_mesh(R, subdivisions=subdiv)
    nv = verts0.shape[0]
    rng = np.random.default_rng(seed)
    # cluster centres on a rough hex ring + one centre cell, all resting on the dish
    z0 = 0.0
    centres = []
    centres.append(np.array([0.0, 0.0, R]))                      # central contact cell
    for k in range(n_cells - 1):
        th = 2.0 * np.pi * k / max(1, n_cells - 1)
        r = 1.9 * R
        centres.append(np.array([r * np.cos(th), r * np.sin(th), R]))
    centres = np.array(centres[:n_cells])

    pos = np.zeros((n_cells * nv, 3), dtype=np.float64)
    cof = np.zeros((n_cells * nv,), dtype=np.int64)
    for c in range(n_cells):
        v = verts0 + centres[c]
        v += rng.normal(0.0, 0.04e-6, size=v.shape)              # light jitter
        pos[c * nv:(c + 1) * nv] = v
        cof[c * nv:(c + 1) * nv] = c
    return pos, cof, z0, R, nv


def _max_abs(a, b):
    a = np.asarray(a, dtype=np.float64); b = np.asarray(b, dtype=np.float64)
    if a.size == 0 and b.size == 0:
        return 0.0
    if a.shape != b.shape:
        return np.inf
    return float(np.max(np.abs(a - b))) if a.size else 0.0


def run(device: str = "cpu") -> bool:
    import warp as wp
    wp.init()

    pos, cof, z0, R, nv = _build_cluster()
    n_cells = int(cof.max()) + 1
    # dt chosen so p_advance ~ 0.5 (v_front·batch·dt·S/l0) — this DRIVES the probabilistic
    # ADVANCE branch so we actually test count-parity of the RNG path (at the physiological
    # dt=8e-6, p_advance~8e-5 and advances ~never fire in a short test). It is a test knob,
    # not a physics change; the RNG branch is identical regardless of p_advance's value.
    dt = 5.0e-2
    params = LamelParams(substrate_clutch=False, batch_steps=50)

    # Multi-tick so the RNG-driven ADVANCE branch (not just first-tick SEED) is exercised:
    # the front exists after tick 1, and a small outward drift of the rim cells makes the
    # leading nodes catch up to their actin -> advance attempts (rng.uniform < p_advance).
    n_ticks = 6
    drift = np.zeros_like(pos)
    rim_for_drift = [1, 2, 3, 4]
    for c in rim_for_drift:                       # push the ring cells slightly outward each tick
        m = cof == c
        ctr = pos[m].mean(0)
        rad = ctr.copy(); rad[2] = 0.0
        rn = np.linalg.norm(rad)
        if rn > 0:
            drift[m] = (rad / rn) * 0.20e-6        # 0.2 um/tick outward

    # --- host reference -----------------------------------------------------
    lam_h = LamellipodiumHost(pos0=pos, cof=cof, n_cells=n_cells, z0=z0, R=R, dt=dt,
                              params=LamelParams(**params.__dict__), use_gpu_ratchet=False)
    # --- GPU geometry path (fed device positions, no full pos.numpy()) ------
    lam_g = LamellipodiumHost(pos0=pos, cof=cof, n_cells=n_cells, z0=z0, R=R, dt=dt,
                              params=LamelParams(**params.__dict__),
                              use_gpu_ratchet=True, device=device)

    Pc = pos.copy()
    for t in range(n_ticks):
        lam_h.update(Pc)
        pos_d = wp.array(np.ascontiguousarray(Pc), dtype=wp.vec3d, device=device)
        lam_g.update(pos_d=pos_d)
        Pc = Pc + drift
    dh = lam_h.upload("cpu")
    dg = lam_g.upload(device)

    print(f"\n=== Lamellipodium ratchet parity (device={device}) ===")
    print(f"  cluster: {n_cells} cells x {nv} nodes = {pos.shape[0]} nodes, R={R*1e6:.1f}um")

    ok = True

    # (1) rim-cell set identical
    rim_h = np.sort(lam_h.rim_cells); rim_g = np.sort(lam_g.rim_cells)
    rim_ok = (rim_h.shape == rim_g.shape) and bool(np.all(rim_h == rim_g))
    ok &= rim_ok
    print(f"  [1] rim set        host={rim_h.tolist()}  gpu={rim_g.tolist()}  "
          f"-> {'IDENTICAL' if rim_ok else 'MISMATCH'}")

    # (2) leading-node geometry within 1e-9
    lh, lg = lam_h._lead, lam_g._lead
    idx_ok = (lh["idx"].shape == lg["idx"].shape) and bool(np.all(lh["idx"] == lg["idx"]))
    cell_ok = (lh["cell"].shape == lg["cell"].shape) and bool(np.all(lh["cell"] == lg["cell"]))
    d_ox = _max_abs(lh["ox"], lg["ox"]); d_oy = _max_abs(lh["oy"], lg["oy"])
    d_pr = _max_abs(lh["proj"], lg["proj"])
    d_cx = _max_abs(lh["ccx"], lg["ccx"]); d_cy = _max_abs(lh["ccy"], lg["ccy"])
    geo_tol = 1e-9
    geo_ok = idx_ok and cell_ok and max(d_ox, d_oy, d_pr) < geo_tol and max(d_cx, d_cy) < (geo_tol * R)
    ok &= geo_ok
    print(f"  [2] leading nodes  host L={lh['idx'].size}  gpu L={lg['idx'].size}  "
          f"idx{'=' if idx_ok else '!'}=  cell{'=' if cell_ok else '!'}=")
    print(f"      max|Δout|={max(d_ox,d_oy):.2e}  max|Δproj|={d_pr:.2e} m  "
          f"max|Δcc|={max(d_cx,d_cy):.2e} m  -> {'OK' if geo_ok else 'FAIL'} (tol 1e-9)")

    # (3) device anchor arrays (what the tether reads) — deterministic geometry exact,
    #     and (shared RNG) the actin pool byte-identical => advance COUNT exact.
    nu_h, nu_g = dh["n_used"], dg["n_used"]
    nl_h, nl_g = dh["n_lead"], dg["n_lead"]
    d_actin = _max_abs(dh["actin"].numpy()[:nu_h], dg["actin"].numpy()[:nu_g]) if nu_h == nu_g else np.inf
    actin_cell_ok = (nu_h == nu_g) and bool(np.all(
        dh["actin_cell"].numpy()[:nu_h] == dg["actin_cell"].numpy()[:nu_g]))
    d_lrp_cc = max(_max_abs(dh["lead_ccx"].numpy()[:nl_h], dg["lead_ccx"].numpy()[:nl_g]),
                   _max_abs(dh["lead_ccy"].numpy()[:nl_h], dg["lead_ccy"].numpy()[:nl_g])) \
        if nl_h == nl_g else np.inf
    seed_ok = (lam_h.n_seeded == lam_g.n_seeded)
    adv_ok = (lam_h.n_advanced == lam_g.n_advanced)
    dev_ok = (nu_h == nu_g) and (nl_h == nl_g) and actin_cell_ok and d_actin < (1e-9) \
        and seed_ok and adv_ok and d_lrp_cc < (1e-9 * R)
    ok &= dev_ok
    print(f"  [3] device anchors n_used host={nu_h} gpu={nu_g}  n_lead host={nl_h} gpu={nl_g}")
    print(f"      seeded host={lam_h.n_seeded} gpu={lam_g.n_seeded}  "
          f"advanced host={lam_h.n_advanced} gpu={lam_g.n_advanced}")
    print(f"      max|Δactin|={d_actin:.2e} m  actin_cell {'=' if actin_cell_ok else '!'}=  "
          f"max|Δlead_cc|={d_lrp_cc:.2e} m  -> {'OK' if dev_ok else 'FAIL'}")
    print(f"      advance path = HOST-FED (path a): shared RNG stream -> COUNT exact, "
          f"not merely statistical.")

    print(f"\n  PARITY {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    import sys
    dev = sys.argv[1] if len(sys.argv) > 1 else "cpu"
    raise SystemExit(0 if run(dev) else 1)
