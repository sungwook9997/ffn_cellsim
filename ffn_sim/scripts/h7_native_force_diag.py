"""Diagnose the ~620 us cupy compartment-force step: snapshot vs ops.

Times, at full-cell N, four set_forces variants on the same cloud:
  (a) snapshot-only      : enter gpu_local_snapshot, read pos, write zeros (no compute)
  (b) snapshot+centroid  : + the masked centroid reduction only
  (c) naive cupy nucleus : the current ~12-op NucleusConfinementGPU
  (d) fused RawKernel     : centroid reduction + ONE CUDA apply kernel
If (a) ~ (c) the gpu_local_snapshot context dominates -> need native C++ (ArrayHandle);
if (a) << (c) the elementwise launches dominate -> the fused single kernel (d) is the win.
"""

from __future__ import annotations

import argparse
import time

import hoomd
import hoomd.md as md
import numpy as np

from ffn_sim.cell.nucleus import resolve_nucleus
from ffn_sim.cell.nucleus_confinement_gpu import NucleusConfinementGPU

_APPLY_SRC = r'''
extern "C" __global__ void nuc_apply(
    const double* pos, const double* mask, const double* centroid,
    double R_nuc, double k_chrom, double k_lamin, double d_knee, double F_knee,
    double* force, double* U, unsigned int N) {
  unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i >= N) return;
  double m = mask[i];
  double dx = pos[3*i]-centroid[0], dy = pos[3*i+1]-centroid[1], dz = pos[3*i+2]-centroid[2];
  double r = sqrt(dx*dx+dy*dy+dz*dz);
  double rs = r > 0.0 ? r : 1.0;
  double nx = dx/rs, ny = dy/rs, nz = dz/rs;
  if (r <= 0.0) { nx = 0.0; ny = 0.0; nz = 0.0; }
  double d = r - R_nuc, ad = fabs(d), sgn = (double)((d>0)-(d<0));
  double e = ad - d_knee, uknee = 0.5*k_chrom*d_knee*d_knee, Fmag, Uu;
  if (ad <= d_knee) { Fmag = -k_chrom*d; Uu = 0.5*k_chrom*d*d; }
  else { Fmag = -sgn*(F_knee+(k_chrom+k_lamin)*e); Uu = uknee+F_knee*e+0.5*(k_chrom+k_lamin)*e*e; }
  force[3*i] = Fmag*nx*m; force[3*i+1] = Fmag*ny*m; force[3*i+2] = Fmag*nz*m;
  U[i] = Uu*m;
}
'''


class _Diag(md.force.Custom):
    def __init__(self, p, rng, mode):
        super().__init__()
        self.p = p; self.t0, self.t1 = rng; self.mode = mode
        self._cp = None; self._k = None

    def set_forces(self, timestep):
        import cupy as cp
        if self._cp is None:
            self._cp = cp
            self._k = cp.RawKernel(_APPLY_SRC, "nuc_apply")
        with self._state.gpu_local_snapshot as snap:
            pos = cp.asarray(snap.particles.position, dtype=cp.float64)
            tag = cp.asarray(snap.particles.tag)
            N = pos.shape[0]
            if self.mode == "snap":
                F = cp.zeros_like(pos); U = cp.zeros(N)
            else:
                mask = ((tag >= self.t0) & (tag < self.t1)).astype(cp.float64)
                mcount = cp.maximum(mask.sum(), 1.0)
                centroid = (pos * mask[:, None]).sum(axis=0) / mcount
                if self.mode == "centroid":
                    F = cp.zeros_like(pos); U = cp.zeros(N)
                else:  # fused
                    F = cp.empty_like(pos); U = cp.empty(N)
                    cen = cp.ascontiguousarray(centroid)
                    posc = cp.ascontiguousarray(pos)
                    blk = 256; grd = (N + blk - 1) // blk
                    self._k((grd,), (blk,), (posc, mask, cen,
                            cp.float64(self.p.R_nuc), cp.float64(self.p.k_chrom),
                            cp.float64(self.p.k_lamin), cp.float64(self.p.d_knee),
                            cp.float64(self.p.F_knee), F, U, cp.uint32(N)))
        with self.gpu_local_force_arrays as arrays:
            arrays.force[:] = F
            arrays.potential_energy[:] = U


def _sim(pos):
    n = pos.shape[0]
    snap = hoomd.Snapshot()
    if snap.communicator.rank == 0:
        L = 60e-6
        snap.configuration.box = [L, L, L, 0, 0, 0]
        snap.particles.N = n; snap.particles.types = ["shell"]
        snap.particles.position[:] = pos
        snap.particles.typeid[:] = np.zeros(n, np.int32)
    sim = hoomd.Simulation(device=hoomd.device.GPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.integrator = md.Integrator(dt=1e-6)
    return sim


def _time(sim, f, steps):
    sim.operations.integrator.forces = [f]
    sim.run(40)
    t0 = time.perf_counter(); sim.run(steps)
    return 1e6 * (time.perf_counter() - t0) / steps


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12840)
    ap.add_argument("--steps", type=int, default=600)
    args = ap.parse_args()
    p = resolve_nucleus({"nucleus": {"E_nuc": 5e3, "ratio_lamin": 3.0, "knee_strain": 0.10,
                                     "critical_pore_area_um2": 7.0}}, R_nuc=3e-6, n_beads=args.n)
    rng = np.random.default_rng(3)
    r = 3e-6 + (rng.random(args.n) - 0.5) * 6 * p.d_knee
    u = rng.normal(size=(args.n, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    pos = (r[:, None] * u).astype(np.float64)
    sim = _sim(pos)
    rg = (0, args.n)
    res = {}
    for mode in ("snap", "centroid", "fused"):
        res[mode] = _time(sim, _Diag(p, rg, mode), args.steps)
    res["naive_cupy"] = _time(sim, NucleusConfinementGPU(p, rg), args.steps)
    print(f"force diag: n={args.n}", flush=True)
    print(f"  (a) snapshot-only     {res['snap']:7.1f} us  (gpu_local_snapshot ctx + read + write)", flush=True)
    print(f"  (b) +centroid reduce  {res['centroid']:7.1f} us", flush=True)
    print(f"  (c) naive cupy (~12op){res['naive_cupy']:7.1f} us", flush=True)
    print(f"  (d) fused RawKernel   {res['fused']:7.1f} us", flush=True)
    print(f"=> snapshot floor {res['snap']:.0f}us; ops above floor: naive "
          f"{res['naive_cupy']-res['snap']:.0f}us -> fused {res['fused']-res['snap']:.0f}us", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
