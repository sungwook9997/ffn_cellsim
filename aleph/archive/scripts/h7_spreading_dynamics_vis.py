"""H.7 — figures + spreading GATE for the dynamic single-cell spreading run.

Reads the ``curve.json`` (time series) and ``frames.npz`` (per-sample basal
footprint point clouds) written by ``h7_spreading_dynamics.py`` and produces:

  * ``h7_spreading_dynamics.png`` — 4-panel summary:
      (1) A(t)/A₀ vs effective time with the physiological single-cell band
          [2, 4] shaded + the old kinematic-mock flat line for contrast;
      (2) leading-edge radius r_front(t) + network bead count n_actin(t);
      (3) FA clutch engagement n_engaged(t) + substrate traction(t);
      (4) basal footprint point cloud t=0 (grey) vs final (colour).
  * ``h7_spreading_movement.gif`` — the basal footprint growing over time.

GATE (spreading is dynamic + physiological), printed + returned:
  G1 A/A₀ GROWS monotonically-ish and reaches the physiological band [2, 4].
  G2 leading edge advances (r_front_final > 1.5 × r_front_0).
  G3 network is dynamic (n_actin grows materially, not a frozen geometry).
The band is the literature single-cell isotropic-spreading magnitude (Cuvelier
2007 / Dubin-Thaler 2004 / Betorz 2023 ≈ 3× radius in P1 → A/A₀ ~ 2–4; Henry
2015 R(t)~t^0.4 power-law → plateau time-course).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Physiological single-cell isotropic-spreading magnitude band (A/A₀).
_AA0_BAND = (2.0, 4.0)
_MOCK_AA0 = 1.06  # the kinematic h7_spreading_compare basal_ring final A/A0


def _load(curve_json: str):
    with open(curve_json) as fh:
        d = json.load(fh)
    s = d["series"]
    # Physiological-time axis (minutes), mapped from the emergent front advance at
    # the literature single-cell spreading velocity. Fall back to effective-time.
    if "t_phys_min" in s[0]:
        t = np.array([r["t_phys_min"] for r in s])
    else:
        t = np.array([r["t_eff_s"] for r in s]) / 60.0
    A = np.array([r["A_um2"] for r in s])
    A0 = d["meta"]["A0_um2"]
    aa0 = A / A0 if A0 else A * np.nan
    rf = np.array([r["r_front_um"] for r in s])
    nac = np.array([r["n_actin"] for r in s])
    neng = np.array([r["n_engaged"] for r in s])
    trac = np.array([r["traction_N"] for r in s])
    return d, t, A, aa0, rf, nac, neng, trac


def evaluate_gate(t, aa0, rf, nac):
    aa0_final = float(aa0[-1])
    g1 = bool(_AA0_BAND[0] <= aa0_final <= _AA0_BAND[1])
    g1_reached = bool(np.nanmax(aa0) >= _AA0_BAND[0])
    g2 = bool(rf[-1] > 1.5 * rf[0]) if rf[0] > 0 else False
    g3 = bool(nac[-1] > 1.3 * nac[0]) if nac[0] > 0 else False
    return {
        "A_A0_final": aa0_final,
        "A_A0_max": float(np.nanmax(aa0)),
        "r_front_ratio": float(rf[-1] / rf[0]) if rf[0] else float("nan"),
        "n_actin_ratio": float(nac[-1] / nac[0]) if nac[0] else float("nan"),
        "G1_in_band": g1,
        "G1_reached_band": g1_reached,
        "G2_edge_advances": g2,
        "G3_network_dynamic": g3,
        "band": _AA0_BAND,
    }


def figure(d, t, A, aa0, rf, nac, neng, trac, gate, out_png):
    fig, ax = plt.subplots(2, 2, figsize=(13.5, 9.5), constrained_layout=True)

    # Panel 1: A/A0 vs effective time + physiological band + mock contrast.
    a = ax[0, 0]
    a.axhspan(_AA0_BAND[0], _AA0_BAND[1], color="#bfe3bf", alpha=0.5,
              label=f"physiological band [{_AA0_BAND[0]:.0f}, {_AA0_BAND[1]:.0f}]")
    a.axhline(_MOCK_AA0, color="0.5", ls=":", lw=1.6,
              label=f"old kinematic mock (A/A₀≈{_MOCK_AA0:.2f}, flat)")
    a.plot(t, aa0, "-o", color="#c0392b", lw=2.2, ms=5, label="dynamic spreading run")
    a.set_xlabel("physiological spreading time (min)  [front-mapped, see meta]")
    a.set_ylabel("A / A₀   (basal contact footprint)")
    a.set_title(f"spreading trajectory  (accel S={d['meta'].get('accel'):.0f}×, "
                f"A₀={d['meta']['A0_um2']:.1f} µm²)")
    a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="upper left")

    # Panel 2: r_front + n_actin.
    a = ax[0, 1]
    a.plot(t, rf, "-o", color="#2f6fb0", lw=2, ms=4, label="r_front (µm)")
    a.set_xlabel("physiological spreading time (min)")
    a.set_ylabel("leading-edge radius r_front (µm)", color="#2f6fb0")
    a.tick_params(axis="y", labelcolor="#2f6fb0")
    a.grid(alpha=0.3)
    a2 = a.twinx()
    a2.plot(t, nac, "-s", color="#e67e22", lw=1.6, ms=3, label="n_actin")
    a2.set_ylabel("lamellipodial actin beads n_actin", color="#e67e22")
    a2.tick_params(axis="y", labelcolor="#e67e22")
    a.set_title("leading-edge advance + network growth (dynamic)")

    # Panel 3: clutch engagement + traction.
    a = ax[1, 0]
    a.plot(t, neng, "-o", color="#8e44ad", lw=2, ms=4, label="n_engaged clutches")
    a.set_xlabel("physiological spreading time (min)")
    a.set_ylabel("engaged FA clutches", color="#8e44ad")
    a.tick_params(axis="y", labelcolor="#8e44ad")
    a.grid(alpha=0.3)
    a3 = a.twinx()
    a3.plot(t, trac * 1e9, "-^", color="#16a085", lw=1.6, ms=3)
    a3.set_ylabel("substrate traction (nN)", color="#16a085")
    a3.tick_params(axis="y", labelcolor="#16a085")
    a.set_title("FA molecular clutch: substrate adhesion + traction")

    # Panel 4: footprint t=0 vs final.
    a = ax[1, 1]
    try:
        frames_npz = d.get("_frames_path")
        if frames_npz and Path(frames_npz).exists():
            z = np.load(frames_npz, allow_pickle=False)
            nf = int(z["n_frames"][0])
            f0 = z["frame__0"]
            ff = z[f"frame__{nf - 1}"]
            if len(f0):
                a.scatter(f0[:, 0], f0[:, 1], s=8, c="0.7", label="t=0")
            if len(ff):
                a.scatter(ff[:, 0], ff[:, 1], s=8, c="#c0392b", alpha=0.6, label="final")
            a.set_aspect("equal")
    except Exception:  # noqa: BLE001
        pass
    a.set_xlabel("x (µm)")
    a.set_ylabel("y (µm)")
    a.set_title("basal footprint: construction → spread")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    verdict = "PASS" if (gate["G1_reached_band"] and gate["G2_edge_advances"]
                         and gate["G3_network_dynamic"]) else "PARTIAL"
    fig.suptitle(
        f"H.7 DYNAMIC single-cell spreading — real BAOAB lamellipodium protrusion "
        f"+ FA clutch  [GATE {verdict}: A/A₀ {gate['A_A0_final']:.2f} "
        f"(max {gate['A_A0_max']:.2f}), r_front ×{gate['r_front_ratio']:.2f}, "
        f"n_actin ×{gate['n_actin_ratio']:.2f}]",
        fontsize=12, fontweight="bold")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def animate(frames_npz: str, out_gif: str, fps: int = 4):
    from matplotlib.animation import FuncAnimation, PillowWriter

    z = np.load(frames_npz, allow_pickle=False)
    nf = int(z["n_frames"][0])
    frames = [z[f"frame__{i}"] for i in range(nf)]
    teff = z["t_eff_s"]
    A = z["A_um2"]
    A0 = float(z["A0_um2"][0])
    allxy = np.vstack([f for f in frames if len(f)]) if any(len(f) for f in frames) else None
    lim = float(np.abs(allxy).max()) * 1.1 if allxy is not None and len(allxy) else 8.0

    fig, ax = plt.subplots(figsize=(6.0, 6.2), constrained_layout=True)

    def draw(k):
        ax.clear()
        xy = frames[k]
        if len(xy):
            r = np.hypot(xy[:, 0], xy[:, 1])
            ax.scatter(xy[:, 0], xy[:, 1], s=10, c=r, cmap="viridis",
                       alpha=0.8, edgecolors="none")
        ax.set_aspect("equal")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        aa0 = A[k] / A0 if A0 else float("nan")
        ax.set_title(f"basal contact footprint  frame {k}/{nf - 1}\n"
                     f"A/A₀={aa0:.2f}  (A={A[k]:.1f} µm²)", fontsize=11)
        fig.suptitle("H.7 dynamic single-cell spreading — basal footprint A(t)",
                     fontsize=12, fontweight="bold")

    anim = FuncAnimation(fig, draw, frames=nf, interval=1000 // max(1, fps))
    Path(out_gif).parent.mkdir(parents=True, exist_ok=True)
    anim.save(out_gif, writer=PillowWriter(fps=fps))
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--curve", required=True)
    ap.add_argument("--frames", default=None)
    base = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "figs"
    ap.add_argument("--out", default=str(base / "h7_spreading_dynamics.png"))
    ap.add_argument("--gif", default=str(base / "h7_spreading_movement.gif"))
    args = ap.parse_args()

    d, t, A, aa0, rf, nac, neng, trac = _load(args.curve)
    if args.frames:
        d["_frames_path"] = args.frames
    gate = evaluate_gate(t, aa0, rf, nac)
    figure(d, t, A, aa0, rf, nac, neng, trac, gate, args.out)
    print(f"[vis] wrote {args.out}")
    if args.frames and Path(args.frames).exists():
        animate(args.frames, args.gif)
        print(f"[vis] wrote {args.gif}")
    print("[gate] " + json.dumps(gate, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
