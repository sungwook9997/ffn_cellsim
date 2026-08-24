"""Developmental-trajectory curves for an FF crawl run — the QUANTITATIVE companion to the interactive
morphology viewer (`ff_crawl_viewer`). Reads the crawl `*_on.npz` (and, if present, the clutches-OFF audit
`*_off.npz`) and plots how the cell DEVELOPS over real time:

  1. COM displacement along the polarity axis phat vs time — clutches ON vs OFF (the audit: protrusion alone is
     internal ⇒ OFF ≈ 0; only traction translocates). Retrograde/crawl reality band 10-100 nm/s (KU-3.12).
  2. Shape POLARIZATION — cortex extent along phat / extent perpendicular; grows as the cell elongates into motion.
  3. Bound focal-adhesion fraction vs time (the clutch treadmill's engaged population).
  4. Cortical σ_vm p95 vs time — where the developing cell is loaded (front-protrusion + rear-traction stress).

Usage: python -m aleph.scripts.ff_dev_curves --npz outputs/ff/figs/<tag>_on.npz [--title ...]
The morphology (3D cell shape) is the primary record in the HTML viewer; these are the observable curves.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _perp_axis(phat: np.ndarray) -> np.ndarray:
    """A unit vector perpendicular to phat (for the lateral-extent / polarization measure)."""
    a = np.array([0.0, 0.0, 1.0]) if abs(phat[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e = a - (a @ phat) * phat
    return e / (np.linalg.norm(e) + 1e-12)


def _series(npz_path: str):
    d = np.load(npz_path)
    frames = np.asarray(d["frames"], np.float64)              # (T, N, 3)
    Nc = int(d["Nc"]); phat = np.asarray(d["phat"], np.float64); phat = phat / (np.linalg.norm(phat) + 1e-12)
    times = np.asarray(d["times"], np.float64) if "times" in d.files else np.arange(frames.shape[0], dtype=float)
    com = np.asarray(d["com"], np.float64) if "com" in d.files else frames[:, :Nc].mean(1)
    if com.shape[0] != frames.shape[0]:                        # com may be per-step; fall back to per-frame centroid
        com = frames[:, :Nc].mean(1)
    e2 = _perp_axis(phat)
    disp_along = (com - com[0]) @ phat                         # µm
    disp_perp = np.linalg.norm((com - com[0]) - np.outer(disp_along, phat), axis=1)
    pol, svm95 = [], []
    cortex = frames[:, :Nc]
    for t in range(frames.shape[0]):
        c = cortex[t] - cortex[t].mean(0)
        ext1 = np.ptp(c @ phat); ext2 = np.ptp(c @ e2)
        pol.append(ext1 / (ext2 + 1e-9))
    svm = np.asarray(d["svm"], np.float64) if "svm" in d.files else None
    svm95 = np.percentile(svm, 95, axis=1) if svm is not None else None
    bfrac = np.asarray(d["bound_frames"], np.float64).mean(1) if "bound_frames" in d.files else None
    return dict(t=times, disp_along=disp_along, disp_perp=disp_perp, pol=np.array(pol), svm95=svm95, bfrac=bfrac)


def build(npz_on: str, out: str, title: str | None = None):
    on = _series(npz_on)
    off_path = npz_on.replace("_on.npz", "_off.npz")
    off = _series(off_path) if Path(off_path).exists() else None
    meta = {}
    jp = npz_on.replace("_on.npz", ".json")
    if Path(jp).exists():
        meta = json.load(open(jp)).get("clutch_on", {})

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(title or "FF cell — developmental trajectory (how it develops over real time)", fontsize=13)

    # (1) COM translocation along phat — ON vs OFF audit
    a = ax[0, 0]
    a.plot(on["t"], on["disp_along"] * 1e3, "-o", color="#2166ac", lw=2, ms=4, label="clutches ON (traction)")
    if off is not None:
        a.plot(off["t"], off["disp_along"] * 1e3, "--s", color="#b2182b", lw=1.5, ms=3, label="clutches OFF (audit ≈ 0)")
    a.axhline(0, color="k", lw=0.6)
    a.set_xlabel("time [s]"); a.set_ylabel("COM displacement ∥ phat [nm]")
    a.set_title("translocation (the cell MOVES)"); a.legend(fontsize=8); a.grid(alpha=0.3)

    # (2) shape polarization (elongation into motion)
    a = ax[0, 1]
    a.plot(on["t"], on["pol"], "-o", color="#1a9850", lw=2, ms=4)
    a.axhline(1.0, color="k", lw=0.6, ls=":")
    a.set_xlabel("time [s]"); a.set_ylabel("extent∥ / extent⊥")
    a.set_title("shape polarization (elongates into motion)"); a.grid(alpha=0.3)

    # (3) bound FA fraction
    a = ax[1, 0]
    if on["bfrac"] is not None:
        a.plot(on["t"], on["bfrac"], "-o", color="#762a83", lw=2, ms=4, label="ON")
    if off is not None and off["bfrac"] is not None:
        a.plot(off["t"], off["bfrac"], "--s", color="#b2182b", lw=1.2, ms=3, label="OFF")
    a.set_ylim(-0.02, 1.05); a.set_xlabel("time [s]"); a.set_ylabel("bound FA fraction")
    a.set_title("clutch engagement"); a.legend(fontsize=8); a.grid(alpha=0.3)

    # (4) cortical stress
    a = ax[1, 1]
    if on["svm95"] is not None:
        a.plot(on["t"], on["svm95"], "-o", color="#e08214", lw=2, ms=4)
    a.set_xlabel("time [s]"); a.set_ylabel("cortex σ_vm p95 [Pa]")
    a.set_title("where the cell is loaded"); a.grid(alpha=0.3)

    v = float(meta.get("v_crawl_nm_s", np.nan)); tr = float(meta.get("traction_nN", np.nan))
    fig.text(0.5, 0.005, f"v_crawl≈{v:.1f} nm/s (reality band 10-100, KU-3.12) · traction≈{tr:.1f} nN · "
             f"net ∥-disp {on['disp_along'][-1]*1e3:+.1f} nm over {on['t'][-1]:.0f} s", ha="center", fontsize=9)
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])
    fig.savefig(out, dpi=130)
    return dict(out=out, disp_along_nm=float(on["disp_along"][-1] * 1e3), has_off=off is not None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True, help="crawl *_on.npz")
    ap.add_argument("--out", default=None)
    ap.add_argument("--title", default=None)
    args = ap.parse_args()
    out = args.out or args.npz.replace("_on.npz", "_devcurves.png").replace(".npz", "_devcurves.png")
    info = build(args.npz, out, title=args.title)
    print(f"wrote {out}  (net ∥-disp {info['disp_along_nm']:+.1f} nm; OFF-audit overlay: {info['has_off']})")


if __name__ == "__main__":
    main()
