"""k_T1(s) calibration figure — grounds the T1 rate law in the small-dt fine-grained mechanics.

Measures the T1 (neighbour-exchange) rate vs shape index from small-dt references at several drives (v0), fits the
barrier stiffness B, and plots k_T1(s) against the grounded rate law + literature anchors (k_endo band KB-4.13,
s0*=5.41). Per PI: figure → ffn_sim/outputs/h_dcm_two_stage/figs/. See DCM_T1_RATE_COARSEGRAIN_DESIGN_2026-07-12.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/Users/sw1/ffn_cellsim")
from ffn_sim.dcm.dcm_t1_rate import (                                    # noqa: E402
    cell_centroids, auto_cutoff, cell_neighbor_set, measure_t1_rate,
    fit_barrier_stiffness, k_t1_arrhenius, K0_DEFAULT, S0_STAR_3D, K_ENDO_LO, K_ENDO_HI)
from ffn_sim.dcm.dcm_jamming_metrics import cell_shape_index_3d          # noqa: E402

V = sys.argv[1] if len(sys.argv) > 1 else "."
DT_FRAME = float(sys.argv[2]) if len(sys.argv) > 2 else 0.6              # 40 frames over 24s
# (label, npz path) — small-dt references at different drive v0 + the frozen/split points for context
RUNS = [
    ("v0=10", f"{V}/s17_kt_v010.npz"),
    ("v0=20", f"{V}/s16_bdf2_dt8e-4.npz"),      # the existing v0=20 small-dt reference (frame every 1.2s)
    ("v0=40", f"{V}/s17_kt_v040.npz"),
]
DT_FRAME_V20 = 1.2                                                        # the v0=20 ref has 20 frames over 24s


def per_cell_rate_and_s(frames, cof, faces, dt_frame):
    """Per-cell: (mean local shape index, T1-participation rate [1/s]). A cell 'participates' in a neighbour
    change if it is in the changed pair."""
    cof = np.asarray(cof); nc = int(cof.max()) + 1
    cens = [cell_centroids(frames[k], cof, nc) for k in range(len(frames))]
    cutoff = auto_cutoff(cens[0])
    sets = [cell_neighbor_set(c, cutoff) for c in cens]
    part = np.zeros(nc)
    for k in range(1, len(sets)):
        for (i, j) in sets[k] ^ sets[k - 1]:
            part[i] += 1; part[j] += 1
    T = dt_frame * (len(frames) - 1)
    # each T1 = loss+gain, and per-cell we counted both endpoints of each change → divide by 2 for events/cell
    rate = part / 2.0 / T
    s_mean = np.nanmean([cell_shape_index_3d(frames[k], faces, cof) for k in range(len(frames))], axis=0)
    return s_mean, rate


agg_s, agg_rate, labels = [], [], []
pc_s, pc_rate = [], []
for lbl, path in RUNS:
    try:
        d = np.load(path)
    except FileNotFoundError:
        print(f"  (skip {lbl}: {path} not found)"); continue
    fr = d["frames"]; cof = d["cof"].astype(int); faces = d["faces"].astype(int)
    dtf = DT_FRAME_V20 if lbl == "v0=20" else DT_FRAME
    r = measure_t1_rate(fr, cof, dtf)
    sf = np.nanmean(cell_shape_index_3d(fr[-1], faces, cof))
    agg_s.append(sf); agg_rate.append(r["rate_per_cell_s"]); labels.append(lbl)
    s_c, rate_c = per_cell_rate_and_s(fr, cof, faces, dtf)
    m = np.isfinite(s_c) & np.isfinite(rate_c)
    pc_s.append(s_c[m]); pc_rate.append(rate_c[m])
    print(f"  {lbl}: aggregate s={sf:.3f} rate={r['rate_per_cell_s']*1e3:.2f}e-3/cell/s")

agg_s = np.array(agg_s); agg_rate = np.array(agg_rate)
B = fit_barrier_stiffness(agg_s, agg_rate, k0=K0_DEFAULT) if (agg_rate > 0).any() else float("nan")
print(f"  fitted barrier stiffness B = {B:.2f}  (k0={K0_DEFAULT} /s, s0*={S0_STAR_3D})")

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
# left: per-cell scatter + aggregate + fitted law
cols = ["#4477aa", "#228833", "#ee6677"]
for i, lbl in enumerate(labels):
    ax[0].scatter(pc_s[i], np.array(pc_rate[i]) * 1e3, s=8, alpha=0.25, color=cols[i % 3], label=f"{lbl} per-cell")
ax[0].scatter(agg_s, agg_rate * 1e3, s=90, color="k", zorder=5, marker="D", label="aggregate")
ss = np.linspace(4.4, 5.6, 100)
if np.isfinite(B):
    ax[0].plot(ss, k_t1_arrhenius(ss, K0_DEFAULT, S0_STAR_3D, B) * 1e3, "k-", lw=2,
               label=f"k_T1=k0·exp(-{B:.1f}·(s0*-s)+)")
ax[0].axvspan(S0_STAR_3D - 0.005, S0_STAR_3D + 0.005, color="crimson", alpha=0.6)
ax[0].axvline(S0_STAR_3D, color="crimson", ls="--", lw=1, label="s0*=5.41 (unjam)")
ax[0].axhspan(K_ENDO_LO * 1e3, K_ENDO_HI * 1e3, color="gray", alpha=0.15, label="k_endo 0.01-0.1/s (KB-4.13)")
ax[0].set_xlabel("3D shape index s"); ax[0].set_ylabel("T1 rate [1e-3 /cell/s]")
ax[0].set_title("k_T1(s) grounded in the small-dt mechanics"); ax[0].legend(fontsize=7, loc="upper left")
ax[0].set_ylim(0, max(agg_rate.max() * 1e3 * 1.3, 60))
# right: log-scale barrier plot
ax[1].scatter(agg_s, agg_rate * 1e3, s=90, color="k", marker="D", zorder=5)
if np.isfinite(B):
    ax[1].plot(ss, k_t1_arrhenius(ss, K0_DEFAULT, S0_STAR_3D, B) * 1e3, "k-", lw=2)
ax[1].axvline(S0_STAR_3D, color="crimson", ls="--", lw=1)
ax[1].axhline(K0_DEFAULT * 1e3, color="gray", ls=":", label="k0 (gate)")
ax[1].set_yscale("log"); ax[1].set_xlabel("3D shape index s"); ax[1].set_ylabel("T1 rate [1e-3 /cell/s] (log)")
ax[1].set_title("barrier: rate → k0 at s0*, suppressed when jammed"); ax[1].legend(fontsize=8)
fig.suptitle("T1 rate calibration — physical rate MEASURED from small-dt mechanics (grounds the KMC rate law)",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = "/Users/sw1/ffn_cellsim/ffn_sim/outputs/h_dcm_two_stage/figs/dcm_t1_rate_calibration.png"
fig.savefig(out, dpi=130)
print("wrote", out, "| B =", round(B, 3))
