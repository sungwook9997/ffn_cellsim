"""Validate the GPU-graft mechanistic lamellipodium (CPU/numpy backend, small N).

Confirms the graft on ``build_gpu_dcm_simulation(lamellipodium=...)`` ENGAGES:
the ratchet activates actin beads (n_used > 0), the traction tether fires
(n_teth > 0), and the spheroid SPREADS (top-down A/A₀ rises, maxZ drops) — the
behaviour the body-force proxy could not produce. Runs on CPU (numpy path of the
device-dispatched forces) so it is testable here before the gbook-GPU production.

Starts from a substrate-TOUCHING small cluster (lowest node at z0), exactly like
production (``dcm_two_stage_production`` translates the aggregate down to touch).

Usage:
    python -m ffn_sim.scripts.lamel_gpu_validate --n 19 --dt 1e-4 \
        --steps 40000 --frames 20 --adh 1e7
"""

from __future__ import annotations

import argparse
import time

import numpy as np
from scipy.spatial import ConvexHull

from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_build import (
    ResolvedGpuDCM, build_gpu_dcm_snapshot, build_gpu_dcm_simulation)
from ffn_sim.archive.hoomd_legacy.cell.dcm_lamellipodium_gpu import ResolvedGpuLamellipodium


def _topdown(pos_mem: np.ndarray) -> float:
    pts = pos_mem[:, :2]
    if pts.shape[0] < 3:
        return 0.0
    try:
        return float(ConvexHull(pts).volume)
    except Exception:
        return 0.0


def _cell_volumes(pos: np.ndarray, faces: np.ndarray, face_cell: np.ndarray,
                  n_cells: int) -> np.ndarray:
    """Per-cell enclosed volume via the divergence (signed-tetra) sum over faces."""
    v0 = pos[faces[:, 0]]; v1 = pos[faces[:, 1]]; v2 = pos[faces[:, 2]]
    vol6 = np.einsum("ij,ij->i", v0, np.cross(v1, v2))
    out = np.zeros(n_cells)
    np.add.at(out, face_cell, vol6)
    return np.abs(out) / 6.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=19)
    ap.add_argument("--dt", type=float, default=1.0e-4)
    ap.add_argument("--steps", type=int, default=40000)
    ap.add_argument("--frames", type=int, default=20)
    ap.add_argument("--gate", type=int, default=2000)
    ap.add_argument("--adh", type=float, default=1.0e7, help="cohesion ω (weak)")
    ap.add_argument("--rep", type=float, default=4.0e7)
    ap.add_argument("--pool", type=int, default=40)
    ap.add_argument("--batch", type=int, default=50)
    ap.add_argument("--vfront-umin", type=float, default=6.0,
                    help="front velocity [µm/min]; raise to see spreading in fewer "
                         "steps during validation (production uses lit 3-12).")
    ap.add_argument("--contact-band", type=float, default=1.5)
    args = ap.parse_args()

    p = ResolvedGpuDCM(adh_strength=args.adh, rep_strength=args.rep, dt=args.dt)
    plam = ResolvedGpuLamellipodium(
        pool_per_cell=args.pool, batch_steps=args.batch,
        v_front=args.vfront_umin * 1e-6 / 60.0, rim_contact_band=args.contact_band)

    # substrate-touch start: build the cluster, translate lowest node to z0.
    b = build_gpu_dcm_snapshot(p, args.n)
    pos0 = np.asarray(b["snap"].particles.position, dtype=np.float64)
    pos0[:, 2] -= (pos0[:, 2].min() - p.z_substrate)

    print(f"[build] N={args.n} dt={args.dt:.1e} adh(ω)={args.adh:.1e} "
          f"pool/cell={args.pool} batch={args.batch}", flush=True)
    t0 = time.time()
    H = build_gpu_dcm_simulation(
        p, args.n, device=None, node_face_contact=True, settle_force=4.0e-10,
        init_pos=pos0, lamellipodium=plam)
    sim = H["sim"]
    st = H["lamel_state"]
    adv = H["lamel_advance"]
    tether = [f for f in sim.operations.integrator.forces
              if type(f).__name__ == "LamellipodialTractionTetherGPU"][0]
    ranges = H["ranges"]
    faces, face_cell = H["faces"], H["face_cell"]
    n_mem = int(ranges[-1][1])
    mem_tid = 0
    V0 = H["V0"]
    p_adv = min(1.0, plam.v_front * (plam.batch_steps * p.dt) * plam.S_kinetic
                / plam.actin_rest_length)
    print(f"[build] rim={H['lamel_rim'].size} actin_pool={st.n_pool} "
          f"p_advance={p_adv:.3e}/tick (build {time.time()-t0:.1f}s)", flush=True)

    def diag():
        snap = sim.state.get_snapshot()
        pos = np.asarray(snap.particles.position, dtype=np.float64)
        tid = np.asarray(snap.particles.typeid)
        mem = pos[:n_mem][tid[:n_mem] == mem_tid]
        vols = _cell_volumes(pos, faces, face_cell, H["n_cells"])
        return dict(area=_topdown(mem), maxZ=float(mem[:, 2].max()),
                    minZ=float(mem[:, 2].min()), vv0=float(vols.mean() / V0))

    sim.run(args.gate)
    d0 = diag(); A0 = d0["area"]
    print(f"[gate] step {args.gate}: A0={A0*1e12:.0f}µm² maxZ={d0['maxZ']*1e6:.1f}µm "
          f"V/V0={d0['vv0']:.3f} n_used={st.n_used} n_teth={tether.n_tethered}",
          flush=True)

    per = max(1, args.steps // args.frames)
    for f in range(1, args.frames + 1):
        sim.run(per)
        d = diag()
        aa = d["area"] / A0 if A0 > 0 else 0.0
        print(f"  [f{f:02d}] step {args.gate+f*per}: A/A0={aa:.3f} "
              f"maxZ={d['maxZ']*1e6:.1f}µm minZ={d['minZ']*1e6:.2f}µm "
              f"V/V0={d['vv0']:.3f} n_used={st.n_used}/{st.n_pool} "
              f"n_teth={tether.n_tethered} seed={adv.n_seeded} adv={adv.n_advanced}",
              flush=True)

    print(f"[done] {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
