"""PI montage: 2-cell contact DOES flatten into a clean junction at every rep.

Corrects commit 9ec3ee1 ("cells overlap as rigid spheres, do NOT flatten"), which was a
SHORT-settle + 3D-trisurf-render artifact (the smooth-sphere render HIDES the flat junction).
At equilibrium (10k settle) the cross-section shows two truncated spheres meeting on a flat
contact face (oblate axial/lateral 0.90-0.94 < 1), with only sub-edge-scale (~0.5µm, 20% of
mean_edge) interface interleaving — NOT macroscopic rigid-sphere interpenetration.
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ffn_sim.scripts.twocell_overlap_diag import diagnose

REPS = ["4e+07", "2e+08", "1e+09", "5e+09"]
fig, axs = plt.subplots(1, 4, figsize=(20, 5.4), sharey=True)
for ax, rep in zip(axs, REPS):
    r = diagnose(f"/tmp/_2c_sweep_{rep}.npz")
    P = r["frame_nodes"]; cof = r["cof"]; axv = r["ax"]; mid = r["mid"]; R = r["R"]
    lat1 = np.cross(axv, [0, 0, 1.0]); lat1 /= np.linalg.norm(lat1)
    lat2 = np.cross(axv, lat1)
    x = (P - mid) @ axv / R; y = (P - mid) @ lat1 / R; z = (P - mid) @ lat2 / R
    slab = np.abs(z) < 0.28
    for cell, col in zip(np.unique(cof[cof >= 0]), ["#1f77b4", "#d62728"]):
        m = (cof == cell) & slab
        ax.scatter(x[m], y[m], s=22, color=col)
    ax.axvline(0, color="k", ls="--", lw=1)
    ax.set_aspect("equal"); ax.set_xlim(-2, 2)
    ax.set_title(f"rep={rep}\noblate(ax/lat)={r['oblate_axial_over_lateral']:.3f}  "
                 f"NN={r['NN_over_R']:.2f}R\ncross={r['cross_midplane_nodes']}/324  "
                 f"deep={r['deep_cross_R']:.3f}R", fontsize=10)
    ax.set_xlabel("contact axis (R)")
axs[0].set_ylabel("lateral (R)")
fig.suptitle("2-cell contact at equilibrium (10k settle): cells FLATTEN into a clean junction at "
             "EVERY rep (flat face at midplane, oblate<1).\n"
             "Crossing is sub-edge mesh interleaving (~0.5µm), NOT rigid-sphere overlap. "
             "Corrects 9ec3ee1's 3D-render 'rigid sphere' misread.", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.92])
out = "ffn_sim/outputs/warp_decohesion/figs/twocell_flattening_VERDICT.png"
fig.savefig(out, dpi=120); print("saved", out)
