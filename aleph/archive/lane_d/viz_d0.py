"""D0 visualization — membrane_cytosol_traction.png + MP4/GIF from d0_trajectory.npz.

Matplotlib only (no CUDA needed): run after pulling the trajectory back from the A5000.

    python -m aleph.archive.lane_d.viz_d0 --npz aleph/outputs/lane_d/d0_trajectory.npz \
        --figdir aleph/outputs/lane_d/figs
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation


def _midplane_extent(origin, dx, n):
    lo = origin[0]
    hi = origin[0] + dx * n
    return [lo, hi, lo, hi]


def _membrane_ring(pos, tol_z):
    """Nodes near the z=0 plane, ordered by angle — the membrane outline in the x-y cut.

    Adaptive: widen the band until at least ~24 nodes are captured (icosphere nodes rarely sit on z=0).
    """
    cz = pos[:, 2].mean()
    for scale in (1.0, 1.6, 2.4, 3.5, 5.0):
        m = np.abs(pos[:, 2] - cz) < tol_z * scale
        if int(m.sum()) >= 24:
            break
    xy = pos[m][:, :2]
    ang = np.arctan2(xy[:, 1] - pos[:, 1].mean(), xy[:, 0] - pos[:, 0].mean())
    order = np.argsort(ang)
    return xy[order], m, order


def make_static(npz, figpath):
    d = np.load(npz, allow_pickle=True)
    pos = d["pos"]; disp = d["disp"]; trac = d["traction"]; pmid = d["pmid"]
    origin = d["origin"]; dx = float(d["dx"]); patch = d["patch"]
    t = d["t"]; nsteps = pos.shape[0]
    pulse_steps = int(d["pulse_steps"])
    kf = min(nsteps - 1, int(pulse_steps))  # frame near end of pulse — strongest signal
    n = pmid.shape[1]
    extent = _midplane_extent(origin, dx, n)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4))
    # --- panel 1: midplane cytosol pressure heatmap + membrane traction vectors ---
    ax = axes[0]
    pm = pmid[kf].T  # (x,y) -> imshow expects [row=y, col=x] with origin lower
    vex = float(np.percentile(np.abs(pmid), 99.5)) or 1.0
    im = ax.imshow(pm, extent=extent, origin="lower", cmap="RdBu_r", aspect="equal",
                   vmin=-vex, vmax=vex)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="cytosol pressure excess  p − p_bar [Pa]")
    # membrane ring + traction-DEVIATION vectors (z≈0 slice). The uniform ~35 pN turgor traction is
    # subtracted so the arrows show how the perturbation REDISTRIBUTES traction (outward bump at the +x
    # patch, faint inward compensation elsewhere) rather than the huge uniform baseline.
    tol_z = 1.2 * dx
    xy, mmask, order = _membrane_ring(pos[kf], tol_z)
    ax.plot(np.append(xy[:, 0], xy[0, 0]), np.append(xy[:, 1], xy[0, 1]), "k-", lw=1.4, label="membrane")
    ring_pos = pos[kf][mmask][order]
    ctr = pos[kf].mean(axis=0)
    nrm = ring_pos - ctr
    nrm /= (np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-12)
    tr = trac[kf][mmask][order]
    outward = np.sum(tr * nrm, axis=1)
    dev = outward - np.mean(outward)            # traction deviation from the ring mean [pN]
    dvec = (dev[:, None] * nrm)[:, :2]
    amax = np.max(np.abs(dev)) + 1e-9
    alen = 2.2 / amax                            # scale so the largest deviation arrow is ~2.2 µm
    ax.quiver(xy[:, 0], xy[:, 1], dvec[:, 0], dvec[:, 1], color="k", scale_units="xy",
              angles="xy", scale=1.0 / alen, width=0.005,
              label=f"traction deviation (max {amax:.2f} pN)")
    ax.plot(patch[0], patch[1], "*", color="lime", ms=18, mec="k", label="perturbation patch")
    ax.set_title(f"Midplane cytosol pressure + membrane traction\n(t={t[kf]:.2f}s, end of pulse)")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("y [µm]"); ax.legend(loc="upper left", fontsize=8)

    # --- panel 2: membrane radial displacement (colored nodes) ---
    ax = axes[1]
    c = pos[0].mean(axis=0)
    normals = disp * 0.0  # placeholder; use radial direction
    rad_dir = pos[kf] - c
    rad_dir /= (np.linalg.norm(rad_dir, axis=1, keepdims=True) + 1e-12)
    radial_disp = np.sum(disp[kf] * rad_dir, axis=1)
    sctr = ax.scatter(pos[kf][:, 0], pos[kf][:, 1], c=radial_disp, cmap="viridis",
                      s=14, vmin=radial_disp.min(), vmax=radial_disp.max())
    plt.colorbar(sctr, ax=ax, fraction=0.046, pad=0.04, label="radial displacement [µm]")
    ax.plot(patch[0], patch[1], "*", color="red", ms=16, mec="k")
    ax.set_title(f"Membrane radial displacement (t={t[kf]:.2f}s)")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("y [µm]"); ax.set_aspect("equal")

    # --- panel 3: scalar ledger vs time ---
    ax = axes[2]
    ax.plot(t, d["p_excess_max"], "-o", ms=3, color="C3", label="max p_excess [Pa]")
    ax.plot(t, 1e3 * d["disp_max"], "-s", ms=3, color="C0", label="max |disp| [nm]")
    ax.plot(t, 1e3 * np.asarray(d["bulge_x"]), "-^", ms=3, color="C2", label="mean +x cap bulge [nm]")
    dtrac = np.asarray(d["trac_plus"]) - np.asarray(d["trac_minus"])
    ax.plot(t, dtrac, "-D", ms=3, color="C4", label="+x−(−x) traction asym [pN]")
    ax.axvspan(0, t[min(pulse_steps, len(t) - 1)], color="orange", alpha=0.12, label="pulse ON")
    ax.set_xlabel("time [s]"); ax.set_ylabel("magnitude")
    ax.set_title("D0 response ledger"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

    fig.suptitle("Lane D0 — membrane ⇄ Biot-Darcy cytosol traction rig  "
                 "(NON-AUTHORITATIVE DIAGNOSTIC)", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(figpath, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("wrote", figpath)


def make_movie(npz, moviepath):
    d = np.load(npz, allow_pickle=True)
    pos = d["pos"]; trac = d["traction"]; pmid = d["pmid"]
    origin = d["origin"]; dx = float(d["dx"]); patch = d["patch"]; t = d["t"]
    pulse_steps = int(d["pulse_steps"])
    n = pmid.shape[1]; extent = _midplane_extent(origin, dx, n)
    vmax = float(np.percentile(np.abs(pmid), 99.5)) or 1.0
    tol_z = 1.2 * dx

    fig, ax = plt.subplots(figsize=(6.5, 6))
    im = ax.imshow(pmid[0].T, extent=extent, origin="lower", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, aspect="equal")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="cytosol pressure excess  p − p_bar [Pa]")
    xy, mmask, order = _membrane_ring(pos[0], tol_z)
    ring, = ax.plot([], [], "k-", lw=1.4)
    q = ax.quiver(xy[:, 0], xy[:, 1], np.zeros(len(xy)), np.zeros(len(xy)), color="k",
                  scale_units="xy", angles="xy", scale=0.5, width=0.005)
    star, = ax.plot([patch[0]], [patch[1]], "*", color="lime", ms=16, mec="k")
    ttl = ax.set_title("")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("y [µm]")

    def _dev(f):
        xy, mmask, order = _membrane_ring(pos[f], tol_z)
        rp = pos[f][mmask][order]
        nrm = rp - pos[f].mean(axis=0)
        nrm /= (np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-12)
        outward = np.sum(trac[f][mmask][order] * nrm, axis=1)
        dev = outward - np.mean(outward)
        return xy, (dev[:, None] * nrm)[:, :2]

    def update(f):
        im.set_data(pmid[f].T)
        xy, mmask, order = _membrane_ring(pos[f], tol_z)
        ring.set_data(np.append(xy[:, 0], xy[0, 0]), np.append(xy[:, 1], xy[0, 1]))
        xy2, dvec = _dev(f)
        q.set_offsets(xy2)
        q.set_UVC(dvec[:, 0], dvec[:, 1])
        state = "pulse ON" if f < pulse_steps else "relaxing"
        ttl.set_text(f"Lane D0  t={t[f]:.2f}s  ({state})")
        return im, ring, q, star, ttl

    anim = animation.FuncAnimation(fig, update, frames=pos.shape[0], interval=120, blit=False)
    try:
        anim.save(moviepath, writer=animation.FFMpegWriter(fps=8, bitrate=1800))
        print("wrote", moviepath)
    except Exception as e:
        gif = os.path.splitext(moviepath)[0] + ".gif"
        anim.save(gif, writer=animation.PillowWriter(fps=8))
        print(f"ffmpeg unavailable ({e}); wrote", gif)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="aleph/outputs/lane_d/d0_trajectory.npz")
    ap.add_argument("--figdir", default="aleph/outputs/lane_d/figs")
    a = ap.parse_args()
    os.makedirs(a.figdir, exist_ok=True)
    make_static(a.npz, os.path.join(a.figdir, "membrane_cytosol_traction.png"))
    make_movie(a.npz, os.path.join(a.figdir, "membrane_cytosol_traction.mp4"))
