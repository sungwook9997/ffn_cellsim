"""Stage 5.1 — CPU small-N validation of the mechanistic lamellipodium engine
on the NEW physiological foundation (γ=2.22e-4, large dt), BEFORE the GPU graft.

Question this answers: does the per-rim-cell lamellipodium engine
(``build_lamellipodium_spheroid``: protrusion + FA clutch + traction tether)
actually FLATTEN cells — top-down A/A₀ rises, maxZ drops, no mesh inversion — at
the new-foundation drag/dt? The engine was validated to A/A₀=2.03 on the OLD
foundation (γ=3.9e-10, S=6e5 accelerated clock); Stage 1 set S=1 + γ=2.22e-4, so
the kinetic budget only closes if dt rises ∝ γ (≈1e-3). This run tests that.

If cells flatten here → proceed to the GPU graft (Stages 2-4). If not → diagnose
the engine before the expensive graft.

Usage:
    python -m ffn_sim.scripts.lamel_validate_cpu --n 7 --dt 1e-3 \
        --steps 80000 --frames 20
"""

from __future__ import annotations

import argparse
import time

import numpy as np
from scipy.spatial import ConvexHull

from ffn_sim.cell.dcm_native_shell import ResolvedNativeDCM
from ffn_sim.cell.dcm_lamellipodium import (
    ResolvedLamellipodiumSpheroid, build_lamellipodium_spheroid)


def _topdown_area(pos: np.ndarray) -> float:
    """xy-silhouette convex-hull area of the given node positions (the PI assay)."""
    pts = pos[:, :2]
    if pts.shape[0] < 3:
        return 0.0
    try:
        return float(ConvexHull(pts).volume)  # 2D hull "volume" == area
    except Exception:
        return 0.0


def _diag(sim, n_mem: int, mem_tid: int, mean_edge0: float):
    """Top-down A (mem nodes), maxZ, minZ, mean-edge ratio (divergence guard)."""
    snap = sim.state.get_snapshot()
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    tid = np.asarray(snap.particles.typeid)
    mem = pos[:n_mem][tid[:n_mem] == mem_tid]
    area = _topdown_area(mem)
    return dict(area=area, maxZ=float(mem[:, 2].max()), minZ=float(mem[:, 2].min()),
                n_part=int(snap.particles.N))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=7)
    ap.add_argument("--dt", type=float, default=1.0e-3)
    ap.add_argument("--gamma", type=float, default=2.22e-4)
    ap.add_argument("--steps", type=int, default=80000)
    ap.add_argument("--frames", type=int, default=20)
    ap.add_argument("--gate", type=int, default=2000)
    ap.add_argument("--contact", action="store_true",
                    help="cells in contact (cohesion) vs gapped (engine-only)")
    args = ap.parse_args()

    base = ResolvedNativeDCM(
        R_cell=7.5e-6, subdivisions=1,
        gamma_node=args.gamma,    # physiological (new foundation)
        dt=args.dt,               # large dt enabled by physiological gamma
        spacing_factor=2.0 if args.contact else 2.6,
    )
    p = ResolvedLamellipodiumSpheroid(base=base)  # gamma_actin/ligand already 2.22e-4

    print(f"[build] N={args.n} dt={args.dt:.1e} gamma={args.gamma:.2e} "
          f"S_kinetic={p.S_kinetic} v_front={p.v_front:.2e} m/s "
          f"contact={args.contact}", flush=True)
    t0 = time.time()
    H = build_lamellipodium_spheroid(p, n_cells=args.n, device=None,
                                     contact=args.contact)
    sim = H["sim"]
    n_mem = H["n_mem"]
    mem_tid = H["mem_typeid"]
    mean_edge0 = H["mean_edge"]
    nuc = H["nucleator"]
    clutch = H["clutch"]
    tether = H["tether"]
    batch_dt = p.batch_steps * base.dt
    p_adv = min(1.0, p.v_front * batch_dt * p.S_kinetic / p.actin_rest_length)
    print(f"[build] rim_cells={H['rim_cells'].size} n_seed={H['n_seed']} "
          f"n_lig={H['n_lig']} p_advance={p_adv:.3e}/tick "
          f"(build {time.time()-t0:.1f}s)", flush=True)

    # finite gate: force-free seeding must not blow up before any long run
    sim.run(args.gate)
    d0 = _diag(sim, n_mem, mem_tid, mean_edge0)
    A0 = d0["area"]
    print(f"[gate] step {args.gate}: A0_topdown={A0*1e12:.0f}µm² "
          f"maxZ={d0['maxZ']*1e6:.1f}µm minZ={d0['minZ']*1e6:.2f}µm", flush=True)

    per = max(1, args.steps // args.frames)
    for f in range(1, args.frames + 1):
        sim.run(per)
        d = _diag(sim, n_mem, mem_tid, mean_edge0)
        aa = d["area"] / A0 if A0 > 0 else 0.0
        print(f"  [f{f:02d}] step {args.gate + f*per}: A/A0={aa:.3f} "
              f"maxZ={d['maxZ']*1e6:.1f}µm minZ={d['minZ']*1e6:.2f}µm "
              f"n_part={d['n_part']} n_grip={clutch.n_grips} "
              f"n_teth={tether.n_tethered} n_actin_add={nuc.n_promoted}",
              flush=True)

    print(f"[done] {time.time()-t0:.0f}s total", flush=True)


if __name__ == "__main__":
    main()
