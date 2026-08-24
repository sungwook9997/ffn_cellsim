#!/usr/bin/env python
"""H.3 motor-engagement / xlink-bond time series from cell-anim capture log.

Parses the ``CAPTURE_FRAME k/K step=S myo_engaged=N xl_bonds=B`` lines
emitted by h3_cell_anim.py and renders a stats panel:

  fig_h3_cell_engagement_trace.png   — myo_engaged vs sim time + xlink
                                       bond count vs sim time + cumulative
                                       binding rate. Two-panel.

Use as a quick-look diagnostic during/after viz runs.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PKG = Path(__file__).resolve().parents[1]
OUT = PKG / "outputs" / "h3" / "figs"

LINE_RE = re.compile(
    r"CAPTURE_FRAME\s+(\d+)/\d+\s+step=(\d+)"
    r"(?:\s+myo_engaged=(\d+)\s+xl_bonds=(\d+))?"
)


def parse_log(log_path: Path, dt_s: float) -> dict:
    rows = []
    for line in log_path.read_text().splitlines():
        m = LINE_RE.search(line)
        if not m:
            continue
        k = int(m.group(1)); step = int(m.group(2))
        if m.group(3) is None:
            continue  # skip the t=0 line (no engagement data)
        rows.append(dict(
            k=k, step=step, t_ms=step * dt_s * 1e3,
            myo_engaged=int(m.group(3)),
            xl_bonds=int(m.group(4)),
        ))
    return dict(rows=rows)


def render(rows: list[dict], out_path: Path, *, n_fil: int, n_motors: int,
           dt_s: float) -> None:
    if not rows:
        raise SystemExit("no CAPTURE_FRAME rows parsed")
    t_ms = np.array([r["t_ms"] for r in rows])
    eng = np.array([r["myo_engaged"] for r in rows])
    xlb = np.array([r["xl_bonds"] for r in rows])
    n_heads = 2 * 10 * n_motors  # 2 sides × n_heads_per_side=10 × n_motors

    fig, axs = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    ax0, ax1 = axs

    # Top: motor engagement
    ax0.plot(t_ms, eng, marker="o", color="#e02020", lw=1.4, ms=3,
             label=f"engaged motor heads")
    ax0.axhline(n_heads, color="gray", ls="--", lw=0.7, alpha=0.5,
                label=f"total heads = {n_heads}")
    ax0.set_ylabel("# engaged heads")
    ax0.set_title(
        f"H.3 full assembly — motor binding kinetics over time "
        f"({n_fil} filaments, {n_motors} minifilaments)")
    ax0.grid(alpha=0.25); ax0.legend(fontsize=9, loc="upper left")
    # Right-axis: % engagement
    ax0r = ax0.twinx()
    ax0r.plot(t_ms, 100 * eng / n_heads, color="#e02020", alpha=0)  # invisible
    ax0r.set_ylim(0, 100 * max(eng.max(), 1) / n_heads * 1.05)
    ax0r.set_ylabel("% of total heads", color="#e02020")

    # Bottom: xlink bond count (total system bonds in current snapshot — gives
    # a stable structural baseline; large vs small swings reveal binding
    # events vs steady backbone).
    ax1.plot(t_ms, xlb, marker="s", color="#19c0ff", lw=1.4, ms=3,
             label=f"total bonds in snap (cortex+xlink+myosin+attach)")
    ax1.set_xlabel(r"sim time $t$ [ms]")
    ax1.set_ylabel("# bonds in snapshot")
    ax1.grid(alpha=0.25); ax1.legend(fontsize=9, loc="lower right")

    # Annotate inferred binding rate (mid-right to avoid legend overlap).
    if len(t_ms) > 1:
        dt_ms = t_ms[-1] - t_ms[0]
        d_eng = eng[-1] - eng[0]
        rate = d_eng / dt_ms if dt_ms > 0 else float("nan")
        # Bell-Evans theoretical equilibrium fraction (no load): k_on /
        # (k_on + k_off0) = 50 / 60 = 83%. Geometric availability is
        # much smaller (only motor heads within head_actin_capture_perp
        # of an actin segment can bind) — this is the simulation's
        # actual equilibrium, not 83%.
        bell_eq_pct = 50.0 / (50.0 + 10.0) * 100
        ax0.text(0.98, 0.40,
                 f"binding gain: {d_eng} heads in {dt_ms:.1f} ms\n"
                 f"net rate ≈ {rate:.1f} heads/ms\n"
                 f"engagement at end: {100*eng[-1]/n_heads:.1f}% of {n_heads} heads\n"
                 f"Bell-Evans no-load max: {bell_eq_pct:.0f}% of geometrically eligible heads",
                 transform=ax0.transAxes, va="center", ha="right", fontsize=9,
                 bbox=dict(boxstyle="round", fc="#fff7e0", ec="#888", alpha=0.92))
    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE_TRACE {out_path}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--log", type=Path,
                    default=Path("/tmp/mac_cell_anim_v4.log"))
    ap.add_argument("--dt", type=float, default=2.094e-7,
                    help="Simulation dt in seconds (constrained-BD dtc).")
    ap.add_argument("--n-fil", type=int, default=120)
    ap.add_argument("--n-motors", type=int, default=40)
    ap.add_argument("--out", type=Path,
                    default=OUT / "fig_h3_cell_engagement_trace.png")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    data = parse_log(args.log, args.dt)
    render(data["rows"], args.out, n_fil=args.n_fil, n_motors=args.n_motors,
           dt_s=args.dt)


if __name__ == "__main__":
    main()
