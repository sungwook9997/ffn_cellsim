"""Build a single FF cell as a full compartment in the VALIDATED physiological RESTING state, and save the
geometry + fields as an npz for the resting viewer. This is the persistent visual record of the 2026-07-08
validated cortical-mechanics state: native cortex (woven actin) + 3000-bead nucleus (R_nuc=0.70R) + plasma
membrane (reservoir) + osmotic turgor at the interphase set-point (ΔP=40 Pa → γ=0.15 mN/m, pressure-borne) +
biphasic drained cytoplasm (K_drained=300 Pa). No compression (strain=0) — the resting baseline.

Run on the gbook A5000: python ff_resting_full_compartment.py <NF> <out.npz>
"""
import sys, time
import numpy as np
from scipy.spatial import ConvexHull
from ffn_sim.ff.gamma_floor import (CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex, TURGOR_DP0)
from ffn_sim.ff.network_warp import simulate_whole_cell_compression_on_device
from ffn_sim.ff.ff_virial_stress import cortex_node_stress
from ffn_sim.common.compartments import resolve_nucleus, resolve_membrane

NF = int(sys.argv[1]) if len(sys.argv) > 1 else 38000
OUT = sys.argv[2] if len(sys.argv) > 2 else "resting_cell.npz"
DEV = sys.argv[3] if len(sys.argv) > 3 else "cuda:0"
N_RELAX = int(sys.argv[4]) if len(sys.argv) > 4 else 6000

t0 = time.time()
cx = build_crosslinked_cortex(CortexParams(), n_filaments=NF, n_xl=NF, n_myo=NF // 10,
                              rng=np.random.default_rng(1))
cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
R0 = cx.R0_mean
nuc = resolve_nucleus(R_nuc_um=0.70 * R0, n_beads=3000)         # R_nuc = 0.70 R (physiological)
mem = resolve_membrane(f_excess=0.25)                          # reservoir-buffered plasma membrane
print(f"# build NF={NF} Nc={cx.net.n_nodes} R0={R0:.3f}µm build={time.time()-t0:.1f}s", flush=True)

t = time.time()
pos_all, m = simulate_whole_cell_compression_on_device(
    cx, NMIIA_MINIFIL_STALL_PN, strain=0.0, nucleus=nuc, membrane=mem,
    pressure_setpoint=float(TURGOR_DP0), K_drained_Pa=300.0, n_steps=N_RELAX,
    turgor_every=50, device=DEV)                               # resting, interphase turgor set-point
Nc, Ne, n_nuc = int(m["Nc"]), int(m["Ne"]), int(m["n_nuc"])
foff = np.asarray(cx.net.fiber_offsets, np.int64)             # cortex fiber offsets (MT off → Ne==Nc)
pcx = pos_all[:Nc]
faces = ConvexHull(pcx).simplices.astype(np.int64)
print(f"# relaxed {N_RELAX} steps {time.time()-t:.0f}s | ΔP={m['dP_turgor_Pa']:.1f} Pa γ_app={m['gamma_apparent_mN_m']:.3f} mN/m "
      f"V/V0={m['V_over_V0']:.3f} R_eq={m['R_eq_um']:.2f}µm", flush=True)

# resting cortex von-Mises field (deviatoric ~ γ-floor low + turgor isotropic) for the FEM-style colouring
V0 = (4.0 / 3.0) * np.pi * R0 ** 3
svm, _p = cortex_node_stress(pcx, np.full(Nc, V0 / Nc),
                             xl_ij=np.stack([cx.xl_i, cx.xl_j], 1).astype(np.int64), k_xl=cx.xl_k, r0_xl=cx.xl_rest,
                             myo_ij=np.stack([cx.myo_i, cx.myo_j], 1).astype(np.int64),
                             f_myo=float(NMIIA_MINIFIL_STALL_PN), dP=float(m["dP_turgor_Pa"]))

np.savez_compressed(OUT, frames=pos_all[None].astype(np.float32), Nc=Nc, Ne=Ne, n_nuc=n_nuc,
                    foff=foff, faces=faces, svm=svm[None].astype(np.float32), R=R0, R_nuc=nuc.R_nuc_um,
                    dP=float(m["dP_turgor_Pa"]), gamma_mN_m=float(m["gamma_apparent_mN_m"]))
print(f"# wrote {OUT}  (Nc={Nc}, n_nuc={n_nuc}, svm p95={np.percentile(svm,95):.1f} Pa)", flush=True)
print("# done", flush=True)
