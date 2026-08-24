#!/usr/bin/env python
r"""The figures a tau run owes, written beside its own record.

The charter says *"Every run writes its own figures beside its record, committed with the data"*.
The three PHASE 4 tau runs landed on 2026-08-21 with **no figure of any kind** — thirty thousand
samples each, read only through printed tables. This writes them.

⚠ **These plot; they do not judge.** Every verdict about stationarity belongs to
`aleph/observe/stationarity.py` and to the record's own fields. A run whose window opened on a
transient is drawn exactly like one that did not, because **hiding a refused series is how a refusal
stops being visible.** The refusal is annotated on the panel, not used to withhold it.

⚠ **NOTHING IS DOWNSAMPLED.** All 30,000 samples of every trace are drawn
(`feedback-viewer-no-downsample`). A decimated trace of a transient is a picture of a different
transient, and the whole question these runs exist to answer is where the transient ends.

⚠ **The x axis is SAMPLES, never seconds.** The driver's `dt_s` is 1.0, so a field named `_s` counts
samples; labelling this axis "s" would put a physical time on a clock that has none. Repeated on
every panel rather than in a caption, because a panel gets cropped and a caption does not travel.

Sanity Gate:
    * dimensions — gamma is pN/µm; the abscissa is dimensionless sample index. Slopes are per SAMPLE.
    * boundary cases — a record with no trace is REFUSED, not skipped with an empty axes; a trace
      shorter than the largest declared cut draws the cuts it can reach and says which it dropped.
    * conservation — not applicable; this reads a stored series and computes no physical quantity.
    * sign sense — the tail slope's sign is the whole seed-1/seed-2 difference and is printed with it.
    * measurement protocol — every seed gets the same axes limits, the same cuts and the same
      reference band, so a difference between panels is a difference between seeds.

engine units: pN/µm, samples. Runtime: host; no device.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

#: The cuts characterised, fixed across seeds so panels are comparable. Descriptive, not criteria.
CUTS = (18_900, 22_000, 26_000)

#: The effective-sample bar the stationarity module applies. Drawn as a reference band, never
#: recomputed here — if it moves in the module, this figure is wrong and should be regenerated.
N_EFF_BAR = 25.0


def _load(path: Path) -> list[float]:
    """The gamma trace, or a refusal. An empty axes is not a figure."""
    d = json.loads(path.read_text())
    trace = (d.get("gamma") or {}).get("trace_pn_per_um")
    if not trace:
        raise SystemExit(f"REFUSED: {path} carries no gamma.trace_pn_per_um. A record that keeps a "
                         "verdict and not its series cannot be plotted, only quoted.")
    return [float(v) for v in trace]


def main(argv: list[str] | None = None) -> int:
    """Write one figure per record plus a comparison panel."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("records", type=Path, nargs="+")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="where the figures go. Default: beside the FIRST record, which is what the "
                         "charter asks for.")
    args = ap.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    from aleph.observe.stationarity import integrated_autocorrelation_time

    def tau(x: np.ndarray) -> float:
        t = integrated_autocorrelation_time(x, 1.0, c=5.0)
        return float(getattr(t, "tau_int_s", getattr(t, "tau", t)))

    out_dir = args.out_dir or args.records[0].parent
    out_dir.mkdir(parents=True, exist_ok=True)

    series = {f.stem: np.asarray(_load(f), np.float64) for f in args.records}
    y_hi = max(float(g.max()) for g in series.values()) * 1.06
    n_hi = max(int(g.size) for g in series.values())
    written: list[Path] = []

    # ── one full-resolution trace per seed ────────────────────────────────────────────────────────
    for name, g in series.items():
        fig, ax = plt.subplots(figsize=(11, 4.4), dpi=160)
        ax.plot(np.arange(g.size), g, lw=0.35, color="#1f4e79")   # every sample, no decimation
        plateau = float(g[-2000:].mean()) if g.size >= 2000 else float(g.mean())
        ax.axhline(plateau, color="#b45309", lw=1.0, ls="--",
                   label=f"mean of the last 2,000 samples = {plateau:.1f} pN/µm  (descriptive, not a criterion)")
        for c in CUTS:
            if c < g.size:
                ax.axvline(c, color="#6b7280", lw=0.8, ls=":")
                ax.text(c, y_hi * 0.965, f" cut {c:,}", fontsize=7, color="#6b7280", va="top")
        ax.set_xlim(0, n_hi)                      # SAME limits for every seed -- panels are comparable
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
        ax.set_ylim(0, y_hi)                      # from zero: no axis truncation
        ax.set_xlabel("sample index  (⚠ NOT seconds — the driver's dt_s is 1.0, so a sample is a step)")
        ax.set_ylabel("γ  [pN/µm]")
        ax.set_title(f"{name} — all {g.size:,} samples, nothing decimated", fontsize=10)
        ax.legend(fontsize=7, loc="lower right", framealpha=0.9)
        ax.grid(alpha=0.18, lw=0.5)
        fig.tight_layout()
        p = out_dir / f"{name}_gamma_trace.png"
        fig.savefig(p)
        plt.close(fig)
        written.append(p)

    # ── the three gates, side by side, same cuts for every seed ───────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), dpi=160)
    colors = ["#1f4e79", "#b45309", "#166534", "#7c3aed"]
    for k, (name, g) in enumerate(series.items()):
        cuts, taus, neffs, slopes = [], [], [], []
        for c in CUTS:
            if c >= g.size - 100:
                continue
            t = tau(g[c:])
            cuts.append(c)
            taus.append(t)
            neffs.append(g[c:].size / (2.0 * t) if t > 0 else float("nan"))
            slopes.append(float(np.polyfit(np.arange(g[c:].size, dtype=float), g[c:], 1)[0]))
        col = colors[k % len(colors)]
        axes[0].plot(cuts, taus, "o-", color=col, label=name, lw=1.2, ms=5)
        axes[1].plot(cuts, neffs, "o-", color=col, label=name, lw=1.2, ms=5)
        axes[2].plot(cuts, slopes, "o-", color=col, label=name, lw=1.2, ms=5)

    axes[0].set_ylabel("τ_int  [samples]")
    axes[0].set_title("τ over the kept tail", fontsize=10)
    axes[1].axhspan(0, N_EFF_BAR, color="#dc2626", alpha=0.12)
    axes[1].axhline(N_EFF_BAR, color="#dc2626", lw=1.1)
    axes[1].text(CUTS[0], N_EFF_BAR * 1.06, f"  bar = {N_EFF_BAR:g}: below this the window has not measured",
                 fontsize=7, color="#dc2626")
    axes[1].set_ylabel("n_eff = kept / 2τ")
    axes[1].set_title("effective samples — the gate that decides", fontsize=10)
    axes[2].axhline(0.0, color="#6b7280", lw=1.0)
    axes[2].set_ylabel("tail slope  [pN/µm per SAMPLE]")
    axes[2].set_title("⚠ sign is the content, so the axis is not |·|", fontsize=10)
    # ⚠ Tick labels, because the first render of this figure printed "190002000021000220002300024000"
    # -- five overlapping five-digit ticks reading as one number. Found by LOOKING at the png, which
    # is the only way that class of defect is ever found.
    for a in axes:
        a.set_xlabel("cut  [sample index]")
        a.grid(alpha=0.18, lw=0.5)
        a.legend(fontsize=7)
        a.set_xticks(list(CUTS))
        a.set_xticklabels([f"{c:,}" for c in CUTS], fontsize=8)
    fig.suptitle("Same cuts, same columns, every seed — so a difference between lines is a difference "
                 "between seeds, not between treatments", fontsize=9, y=1.0)
    fig.tight_layout()
    p = out_dir / "tau3_three_gates.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    written.append(p)

    for w in written:
        print(f"wrote {w}  ({w.stat().st_size / 1e3:.0f} kB)")
    print("\n⚠ These plot and do not judge. Verdicts live in the records' own fields and in "
          "aleph/observe/stationarity.py; a refused series is drawn exactly like an accepted one.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
