"""DECISIVE aggregate-σ compaction test at PHYSIOLOGICAL TIMESCALE, using the finished #1 IPC.

Context (2026-07-08 reconciliation):
- The 2026-07-02 audit (DCM_AGGREGATE_COMPACTION §CORRECTION #2) fixed a 10⁶× unit bug
  (dP_agg = 2.0e6·σ/R → 2·σ/R; DCM is in METERS). After the fix, the N=100 aggregate-σ re-run
  showed NO compaction and CLEAN contact (pen ~0.3, not the buggy 123). BUT that null was over only
  ~0.06 s of physical time (8000 steps × 8e-6) — fidelity red-flag #2 (runs never reach min-hr
  mechanobiology). It never tested whether physical-magnitude Foty-Steinberg σ compacts over a
  physiological window.
- #1 (projected-Newton IPC + implicit large-dt, dt ceiling ~2 s) is exactly the tool that reaches
  min-hr physical time in feasible step counts. This run uses "more physical time" (a legit lever,
  NOT σ-amplification) to give the lit-anchored aggregate σ its actual timescale.

The controlled experiment: N=400 loose (gap 2.4) dispersed cluster, full compartment stack
(cadherin bundle-10, nucleus E_nuc=399, cortex γ), --ipc --ipc-newton implicit at accel_dt, run over
MINUTES of physical time. Two conditions differing ONLY in σ:
  SIGMA=0    → baseline (no aggregate tension)
  SIGMA>0    → + Foty-Steinberg aggregate σ (lit 1-20 mN/m; default 5)
The delta (baseline vs +σ) isolates the aggregate-tension effect and cancels any shared large-dt
soft-force lag. Question: does physical-magnitude aggregate σ drive CLEAN compaction (Rg↓, porosity↓,
pen~0) once given physiological time?

Env: NCELLS STEPS ACCEL_DT SIGMA GAP BUNDLE ENUC RNUC GAMMA NEWTON_MAX FRAMES TAG NUCLEUS.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.expanduser("~/ffn_cellsim"))
sys.path.insert(0, os.path.expanduser("~/ff_scratch"))
from scipy.spatial import cKDTree, ConvexHull

from ffn_sim.dcm.dcm_warp_decohesion import run_decohesion

OUT = os.path.expanduser("~/ff_scratch/_prod_out")
os.makedirs(OUT, exist_ok=True)

N        = int(os.environ.get("NCELLS", "400"))
STEPS    = int(os.environ.get("STEPS", "10000"))
ACCEL_DT = float(os.environ.get("ACCEL_DT", "8e-3"))    # implicit real-time dt [s]; 8e-3 = feasible ceiling
                                                        # for a loose start (Phase B 2026-07-09; 0.2 blows up)
SIGMA    = float(os.environ.get("SIGMA", "5.0e-3"))     # aggregate σ [N/m]; 5 mN/m = compaction sweet spot
                                                        # (Phase B non-monotonic; σ>5 = σ-dt over-shoot); 0 → baseline
GAP      = float(os.environ.get("GAP", "2.4"))          # loose start: centre spacing gap·R (>2 ⇒ void to compact)
BUNDLE   = float(os.environ.get("BUNDLE", "10"))        # cadherin N_cad bundle (~1.7 nN, lit KB-4.11; 40 diverges)
ENUC     = float(os.environ.get("ENUC", "399"))         # MCF7 nucleus in-situ [Pa] (audit#19)
RNUC     = float(os.environ.get("RNUC", "0.7"))         # nucleus radius factor (lit MCF7)
GAMMA    = float(os.environ.get("GAMMA", "5.0e-4"))     # cortical/surface tension [N/m]
NEWTONMAX= int(os.environ.get("NEWTON_MAX", "8"))
WARMUP   = int(os.environ.get("WARMUP", "2000"))        # soft-start steps at 0.1× dt
FRAMES   = int(os.environ.get("FRAMES", "40"))
NUCLEUS  = os.environ.get("NUCLEUS", "1") == "1"
MATURE   = os.environ.get("MATURE", "0") == "1"         # cadherin junction maturation (slow rearrangement viscosity)
TAU_MAT  = float(os.environ.get("TAU_MATURE", "600"))   # maturation timescale [s] (KB-4.11 5-30min)
MAT_LIFE = float(os.environ.get("MATURE_LIFETIME", "600"))  # mature junction lifetime [s] (KB-4.11)
CLUSTER  = os.environ.get("CLUSTER", "0") == "1"        # S3 load-sharing cluster break (junction=bundle_n parallel
                                                        # trans-dimers, m→0 death; +MATURE → contact-age maturation ENGAGES)
NNASCENT = int(os.environ.get("N_NASCENT", "4"))        # nascent cluster size (KB-4.3 controlled var)
BUILDER  = os.environ.get("BUILDER", "fcc")             # 'fcc' loose | 'confluent' Voronoi space-filling (feasible)
INSET    = float(os.environ.get("INSET", "0.0"))        # >0 → shrink cells at build so they START non-overlapping (G2 fix)
INIT_NPZ = os.environ.get("INIT_NPZ", "") or None       # restart from a saved aggregate npz (frames/faces/cof)
TAG      = os.environ.get("TAG", "")

npz = f"{OUT}/agg_compaction_n{N}_sig{SIGMA:.0e}_adt{ACCEL_DT:g}{TAG}.npz"

print(f"[AGG-COMPACT] N={N} steps={STEPS} accel_dt={ACCEL_DT}s (phys time≈{STEPS*ACCEL_DT:.0f}s="
      f"{STEPS*ACCEL_DT/60:.1f}min) sigma={SIGMA*1e3:.1f}mN/m gap={GAP} bundle={BUNDLE} "
      f"nucleus={NUCLEUS} -> {npz}", flush=True)

t0 = time.time()
r = run_decohesion(
    n_cells=N, subdiv=2, steps=STEPS, frames=FRAMES, device="cuda:0",
    builder=BUILDER, gap=GAP, init_npz=INIT_NPZ, inset=INSET,
    v0_from_init=(BUILDER == "confluent" or INIT_NPZ is not None), lloyd_iters=6,
    # --- contact: finished #1 projected-Newton IPC (log-barrier, guaranteed non-penetration) ---
    conservative_contact=False, ipc=True,
    ipc_newton=True, ipc_newton_max=NEWTONMAX, ipc_newton_tol=1e-4,
    integrator="implicit", accel_dt=ACCEL_DT,
    # --- cell-cell adhesion: explicit cadherin catch-bonds, bundle-10 (de-cohesion emergent) ---
    cadherin=True, cad_bundle=BUNDLE,
    cad_mature=MATURE, cad_tau_mature=TAU_MAT, cad_mature_lifetime=MAT_LIFE,
    cad_cluster=CLUSTER, cad_n_nascent=NNASCENT,
    # --- full compartment stack at physiological setpoints ---
    nucleus=NUCLEUS, E_nuc=ENUC, R_nuc_factor=RNUC, ratio_lamin=1.4,
    surface_tension=True, gamma_surf=GAMMA, diff_tension=True,
    # --- the driver under test: aggregate-level Foty-Steinberg σ (SIGMA=0 ⇒ off = baseline) ---
    aggregate_tension=(SIGMA > 0.0), sigma_agg=SIGMA, agg_every=25,
    warmup=WARMUP, save_frames=npz)

# ---- compaction metrics from the saved trajectory (DCM pos in METERS ×1e6 → µm) ----
d = np.load(npz)
fr = d["frames"]
faces = d["faces"].astype(int)
cof = d["cof"].astype(int)
nc = int(cof.max()) + 1
steps_arr = d["step"] if "step" in d.files else np.arange(len(fr))


def cell_volumes(pos):
    v0 = pos[faces[:, 0]]; v1 = pos[faces[:, 1]]; v2 = pos[faces[:, 2]]
    tv = np.einsum("ij,ij->i", v0, np.cross(v1, v2)) / 6.0
    V = np.zeros(nc)
    np.add.at(V, cof[faces[:, 0]], tv)
    return np.abs(V)


def com_rg(p):
    cen = np.array([p[cof == c].mean(0) for c in range(nc)])
    com = cen.mean(0)
    return float(np.sqrt(((cen - com) ** 2).sum(1).mean())) * 1e6


def porosity(p):
    try:
        Vh = float(ConvexHull(p).volume)
        Vc = float(cell_volumes(p).sum())
        return 1.0 - Vc / Vh
    except Exception:
        return float("nan")


def min_gap(p):
    # EFFICIENT (2026-07-09): ONE global KDTree + k-NN; per cell, the nearest node on ANOTHER cell = its
    # inter-cell gap. Was per-cell cKDTree over ~all other nodes → O(cells×nodes) (minutes at N=2000, it
    # stalled the overnight Phase A inline TREND). k-NN + node subsample is O(nodes·log) — seconds at N=2000.
    tree = cKDTree(p)
    K = min(32, len(p))
    gaps = []
    for c in range(nc):
        my = np.where(cof == c)[0]
        if my.size < 4:
            continue
        if my.size > 24:                                 # subsample nodes (cell min-gap robust to it)
            my = my[:: max(1, my.size // 24)]
        dd, ii = tree.query(p[my], k=K)
        best = np.inf
        for row_d, row_i in zip(dd, ii):
            mask = cof[row_i] != c                        # neighbours on another cell
            if mask.any():
                best = min(best, float(row_d[mask][0]))   # nearest other-cell node for this node
        if np.isfinite(best):
            gaps.append(best)
    return np.median(np.array(gaps)) * 1e6 if gaps else float("nan")


print(f"\n[TREND] N={N} cells={nc} sigma={SIGMA*1e3:.1f}mN/m  "
      f"({time.time()-t0:.0f}s wall, {STEPS*ACCEL_DT/60:.1f}min phys)", flush=True)
print(f"{'frame':>5} {'step':>8} {'Rg[um]':>8} {'poros':>7} {'gap[um]':>8}", flush=True)
for t in range(len(fr)):
    p = fr[t].astype(float)
    print(f"{t:>5} {int(steps_arr[t]):>8} {com_rg(p):>8.2f} {porosity(p):>7.3f} {min_gap(p):>8.3f}",
          flush=True)

rg0, rg1 = com_rg(fr[0].astype(float)), com_rg(fr[-1].astype(float))
p0, p1 = porosity(fr[0].astype(float)), porosity(fr[-1].astype(float))
print(f"\n[SUMMARY sigma={SIGMA*1e3:.1f}mN/m] Rg {rg0:.2f}->{rg1:.2f}um "
      f"({100*(rg1/max(rg0,1e-9)-1):+.2f}%)  porosity {p0:.3f}->{p1:.3f}  "
      f"vv0={r.get('vv0_final','?')}  pen_peak={r.get('pen_frac_peak','?')}  "
      f"pen_final={r.get('pen_frac_final','?')}", flush=True)
print("AGG-COMPACT DONE", flush=True)
