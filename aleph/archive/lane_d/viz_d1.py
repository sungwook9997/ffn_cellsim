"""D1 visualization — force-path + timeseries figures (+ MP4) from d1_trajectory.npz.

    python -m aleph.archive.lane_d.viz_d1 --npz aleph/outputs/lane_d/d1_trajectory.npz \
        --figdir aleph/outputs/lane_d/figs

Matplotlib only. Projection: x (push axis) horizontal, z vertical.
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection


def _seg(pos, ia, ib):
    return np.stack([pos[ia][:, [0, 2]], pos[ib][:, [0, 2]]], axis=1)


def _mt_struts(pos, fiber_off, off_mt):
    segs = []
    fo = np.asarray(fiber_off)
    for i in range(len(fo) - 1):
        a = off_mt + int(fo[i]); b = off_mt + int(fo[i + 1])
        nodes = pos[a:b][:, [0, 2]]
        for k in range(len(nodes) - 1):
            segs.append([nodes[k], nodes[k + 1]])
    return np.asarray(segs)


def make_force_path(npz, figpath):
    d = np.load(npz, allow_pickle=True)
    posA = d["framesA"][-1]; posB = d["framesB"][-1]
    ns, ne = d["nuc_slice"]; as_, ae = d["anc_slice"]
    mem_n = int(d["mem_n"]); off_mt = int(d["off_mt"])
    linc_n = d["linc_n"]; linc_a = d["linc_a"]
    cap_ia = d["cap_ia"]; cap_ib = d["cap_ib"]

    fig, axes = plt.subplots(1, 2, figsize=(15, 7.2))
    for ax, pos, tag, chl, chi, react_label in (
        (axes[0], posA, "A: membrane push → nucleus", d["chA_linc"], d["chA_if"], "nucleus reaction"),
        (axes[1], posB, "B: nuclear pull → membrane", d["chB_linc"], d["chB_if"], "membrane reaction"),
    ):
        # membrane ring (nodes near z-mid band, ordered by angle)
        mem = pos[:mem_n]
        band = np.abs(mem[:, 1]) < 1.2
        if band.sum() < 20:
            band = np.argsort(np.abs(mem[:, 1]))[:40]
            band = np.isin(np.arange(mem_n), band)
        ring = mem[band][:, [0, 2]]
        ang = np.arctan2(ring[:, 1], ring[:, 0]); ring = ring[np.argsort(ang)]
        ax.plot(np.append(ring[:, 0], ring[0, 0]), np.append(ring[:, 1], ring[0, 1]),
                "-", color="0.5", lw=1.3, label="membrane frame")
        # nucleus node cloud
        nuc = pos[ns:ne][:, [0, 2]]
        ax.scatter(nuc[:, 0], nuc[:, 1], s=9, color="#8888dd", alpha=0.6, label="nucleus", zorder=2)
        # MT struts colored by |bending force| (compression/buckling resistance)
        mtseg = _mt_struts(pos, d["mt_fiber_off"], off_mt)
        if len(mtseg):
            ax.add_collection(LineCollection(mtseg, colors="#2ca02c", linewidths=1.0, alpha=0.55,
                                             label="MT struts"))
        # capture tethers (thin)
        capseg = _seg(pos, cap_ia, cap_ib)
        ax.add_collection(LineCollection(capseg, colors="0.75", linewidths=0.6, alpha=0.6))
        # LINC joints colored by tension (equal-opposite line nucleus↔anchor)
        lincseg = _seg(pos, linc_a, linc_n)
        tens = np.linalg.norm(chl[linc_a], axis=1)
        lc = LineCollection(lincseg, cmap="Reds", linewidths=2.2, zorder=3)
        lc.set_array(tens); lc.set_clim(0, max(1e-6, tens.max()))
        ax.add_collection(lc)
        cb = plt.colorbar(lc, ax=ax, fraction=0.045, pad=0.03)
        cb.set_label("LINC tension [pN]")
        bound = np.linalg.norm(chl[linc_a], axis=1) > 1e-9
        anc = pos[as_:ae][:, [0, 2]]
        ax.scatter(anc[:, 0], anc[:, 1], s=22, marker="s", facecolor="none",
                   edgecolor="k", linewidths=0.8, label="capture anchor", zorder=4)
        # reaction arrow on nucleus (or membrane) — net force direction
        f_react = chl + chi
        if tag.startswith("A"):
            tgt = pos[ns:ne][:, [0, 2]]; fr = np.sum(f_react[ns:ne][:, [0, 2]], axis=0)
            cen = tgt.mean(0)
        else:
            # pulse B reaction travels LINC/IF → anchors → (capture) → membrane frame; report the net
            # reaction on the anchor shell (what the pinned membrane must hold).
            fr = np.sum(f_react[as_:ae][:, [0, 2]], axis=0)
            cen = pos[as_:ae][:, [0, 2]].mean(0)
        frn = fr / (np.linalg.norm(fr) + 1e-12) * 2.5
        ax.annotate("", xy=(cen[0] + frn[0], cen[1] + frn[1]), xytext=(cen[0], cen[1]),
                    arrowprops=dict(arrowstyle="-|>", color="crimson", lw=2.5))
        ax.text(cen[0], cen[1] - 1.0, f"{react_label}\n|Σf|={np.linalg.norm(fr):.1f} pN",
                color="crimson", fontsize=9, ha="center")
        ax.set_title(tag, fontsize=11); ax.set_xlabel("x [µm]"); ax.set_ylabel("z [µm]")
        ax.set_aspect("equal"); ax.set_xlim(-9, 9); ax.set_ylim(-9, 9)
        ax.legend(loc="upper left", fontsize=7)
    fig.suptitle("Lane D1 — MT / IF / LINC internal force path  (NON-AUTHORITATIVE DIAGNOSTIC)\n"
                 "MT = bending/buckling (green), LINC = nonlinear tension-only cabling (red), "
                 "capture = grey", fontsize=12, y=1.03)
    fig.tight_layout()
    fig.savefig(figpath, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("wrote", figpath)


def make_timeseries(npz, figpath):
    d = np.load(npz, allow_pickle=True)
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.2))
    # pulse A
    ax = axes[0]
    x = d["A_applied_um"]
    ax.plot(x, d["A_load_at_patch_pn"], "-o", ms=3, color="C0", label="applied load @ patch")
    ax.plot(x, d["A_nuc_reaction_pn"], "-s", ms=3, color="crimson", label="nucleus reaction |Σf|")
    ax.plot(x, d["A_linc_pn"], "-^", ms=3, color="C3", label="max LINC tension")
    ax.plot(x, d["A_if_pn"], "-v", ms=3, color="C1", label="max IF tension")
    ax.plot(x, d["A_mt_pn"], "-d", ms=3, color="C2", label="max MT bending")
    ax.set_yscale("log"); ax.set_xlabel("prescribed patch push [µm]"); ax.set_ylabel("force [pN]")
    ax.set_title("Pulse A — membrane push → nucleus reaction"); ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
    # pulse B
    ax = axes[1]
    x = d["B_prescribed_um"]
    ax.plot(x, d["B_membrane_reaction_pn"], "-s", ms=3, color="crimson", label="membrane-anchor reaction |Σf|")
    ax.plot(x, d["B_anchor_reaction_pn"], "-o", ms=3, color="C0", label="anchor reaction |Σf|")
    ax.plot(x, d["B_linc_pn"], "-^", ms=3, color="C3", label="max LINC tension")
    ax.plot(x, d["B_if_pn"], "-v", ms=3, color="C1", label="max IF tension")
    ax.set_xlabel("prescribed nuclear displacement [µm]  (EXTERNAL input)")
    ax.set_ylabel("force [pN]")
    ax.set_title("Pulse B — nuclear pull → membrane reaction"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.suptitle("Lane D1 — force-path ledger (reciprocal load transfer)  NON-AUTHORITATIVE DIAGNOSTIC",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(figpath, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("wrote", figpath)


def make_movie(npz, moviepath):
    d = np.load(npz, allow_pickle=True)
    fA = d["framesA"]; fB = d["framesB"]
    ns, ne = d["nuc_slice"]; as_, ae = d["anc_slice"]; mem_n = int(d["mem_n"])
    off_mt = int(d["off_mt"]); linc_n = d["linc_n"]; linc_a = d["linc_a"]
    frames = [("A push", f) for f in fA] + [("B pull", f) for f in fB]

    fig, ax = plt.subplots(figsize=(6.6, 6.6))
    ax.set_xlim(-9, 9); ax.set_ylim(-9, 9); ax.set_aspect("equal")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("z [µm]")
    nucsc = ax.scatter([], [], s=9, color="#8888dd", alpha=0.6)
    mtlc = LineCollection([], colors="#2ca02c", linewidths=1.0, alpha=0.5); ax.add_collection(mtlc)
    linclc = LineCollection([], cmap="Reds", linewidths=2.2); ax.add_collection(linclc)
    ancsc = ax.scatter([], [], s=20, marker="s", facecolor="none", edgecolor="k", linewidths=0.7)
    ttl = ax.set_title("")

    def update(i):
        tag, pos = frames[i]
        nucsc.set_offsets(pos[ns:ne][:, [0, 2]])
        ancsc.set_offsets(pos[as_:ae][:, [0, 2]])
        mtlc.set_segments(_mt_struts(pos, d["mt_fiber_off"], off_mt))
        seg = _seg(pos, linc_a, linc_n); linclc.set_segments(seg)
        # LINC extension as tension proxy
        ext = np.linalg.norm(pos[linc_a] - pos[linc_n], axis=1)
        linclc.set_array(ext); linclc.set_clim(ext.min(), ext.max() + 1e-6)
        ttl.set_text(f"Lane D1 — pulse {tag}  (frame {i})")
        return nucsc, ancsc, mtlc, linclc, ttl

    anim = animation.FuncAnimation(fig, update, frames=len(frames), interval=140, blit=False)
    try:
        anim.save(moviepath, writer=animation.FFMpegWriter(fps=8, bitrate=1800)); print("wrote", moviepath)
    except Exception as e:
        gif = os.path.splitext(moviepath)[0] + ".gif"
        anim.save(gif, writer=animation.PillowWriter(fps=8)); print(f"ffmpeg unavailable ({e}); wrote", gif)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="aleph/outputs/lane_d/d1_trajectory.npz")
    ap.add_argument("--figdir", default="aleph/outputs/lane_d/figs")
    a = ap.parse_args()
    os.makedirs(a.figdir, exist_ok=True)
    make_force_path(a.npz, os.path.join(a.figdir, "mt_if_linc_force_path.png"))
    make_timeseries(a.npz, os.path.join(a.figdir, "mt_if_linc_timeseries.png"))
    make_movie(a.npz, os.path.join(a.figdir, "mt_if_linc_force_path.mp4"))
