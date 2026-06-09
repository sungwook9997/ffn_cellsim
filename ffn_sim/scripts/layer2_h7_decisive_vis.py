"""H.7-decisive figure: A/A0(f_active) curve with the H.7 MECHANISTIC traction markers vs PI band.

Reads outputs/layer2/h7_decisive/h7_decisive_n{N}_s{seed}.json (one list per seed) and renders
the FORK verdict: sweep the CBM per-cell active traction and overlay (1) the H.7 mechanistic
test-density value (~2.8 nN), (2) the H.7 physiological-density value (102 nN, Gil-Redondo),
the old heuristic anchors (1.6 / 9.4 nN), and the PI A/A0 band [7-10]. If the curve enters the
band at the H.7 value → missing magnitude; if it saturates ~1-2 and EJECTS before the band → CBM
1-particle structural limit (H.7's mechanistic traction confirms route a owns the magnitude).

Visualization-integrity rules: no axis truncation; PI band overlaid on the measurement; units
annotated (nN, A/A0); per-seed thin lines + ensemble mean overlay.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "layer2"
_DAT = _OUT / "h7_decisive"
_FIG = _OUT / "figs" / "fig_layer2_h7_decisive.png"

_PI_BAND = (7.0, 10.0)
# x-position to draw the f=0 baseline on a log/symlog axis (a hair left of the smallest force).
_X0 = 0.3  # nN


def _load():
    seeds, header = {}, None
    for p in sorted(_DAT.glob("h7_decisive_n*_s*.json")):
        rows = json.loads(p.read_text())
        for r in rows:
            if r.get("header"):
                header = r
                continue
            seeds.setdefault((r["n"], r["seed"]), []).append(r)
    return seeds, header


def main() -> int:
    seeds, header = _load()
    if not seeds:
        print(f"no data in {_DAT}")
        return 1

    fig, ax = plt.subplots(figsize=(9.6, 6.4))

    # per-seed thin lines
    agg = defaultdict(list)
    for (n, s), rows in sorted(seeds.items()):
        rows = sorted(rows, key=lambda r: r["f_active_nN"])
        xs = [max(r["f_active_nN"], _X0) for r in rows]
        ys = [r["aa0_core"] for r in rows]
        ax.plot(xs, ys, "-", color="0.72", lw=1.0, alpha=0.8, zorder=2)
        for r in rows:
            agg[r["f_active_nN"]].append((r["aa0_core"], r["ejected"]))
        for r in rows:
            if r["ejected"]:
                ax.scatter([max(r["f_active_nN"], _X0)], [r["aa0_core"]], marker="x",
                           color="tab:red", s=70, zorder=5)

    # ensemble mean
    fx = sorted(agg)
    xm = [max(f, _X0) for f in fx]
    ym = [float(np.mean([v[0] for v in agg[f]])) for f in fx]
    ys = [float(np.std([v[0] for v in agg[f]])) for f in fx]
    ej = [any(v[1] for v in agg[f]) for f in fx]
    ax.errorbar(xm, ym, yerr=ys, fmt="o-", color="tab:blue", lw=2.0, capsize=4,
                label="A/A0 connected-core (ensemble mean)", zorder=4)

    # PI band (overlay-only)
    ax.axhspan(*_PI_BAND, color="tab:orange", alpha=0.15,
               label=f"PI A/A0 band {_PI_BAND[0]:.0f}-{_PI_BAND[1]:.0f} (overlay)")
    ax.axhline(1.0, color="0.5", ls="--", lw=1.0, label="A/A0 = 1 (no spread)")

    # H.7 mechanistic markers
    f_test = header["h7_f_cell_test_nN"] if header else 2.82
    f_phys = header["h7_f_cell_physio_nN"] if header else 102.0
    ax.axvline(f_test, color="tab:green", ls="-", lw=1.6, alpha=0.9,
               label=f"H.7 mechanistic test-density {f_test:.1f} nN")
    ax.axvline(f_phys, color="tab:purple", ls="-", lw=1.6, alpha=0.9,
               label=f"H.7 mechanistic physiological {f_phys:.0f} nN (Gil-Redondo)")
    ax.axvline(3.0, color="0.6", ls=":", lw=1.0, label="B1 stability ceiling 3 nN")

    ax.set_xscale("log")
    ax.set_xlabel("CBM per-cell active traction f_active  [nN]   (f=0 drawn at 0.3)")
    ax.set_ylabel("spread ratio  A/A0  (connected core)")
    prov = "PROVISIONAL (single-SF × N_SF)" if (header and header.get("h7_provisional")) else "array aggregate"
    ax.set_title(
        "DECISIVE re-run: A/A0(f_active) with H.7 MECHANISTIC traction vs the heuristic bracket\n"
        f"H.7 traction {prov}: +131 pN/SF × {header.get('h7_n_SF', 21.5) if header else 21.5:.1f} SF "
        f"= {f_test:.1f} nN (test), {f_phys:.0f} nN (physiol., ×{header.get('h7_density_gap_x','?') if header else '?'} density)\n"
        "A/A0 peaks ~2.2 then EJECTS (×) before the PI band → CBM 1-particle structural limit",
        fontsize=9.5)
    top = max(max(ym) + max(ys) + 0.3, _PI_BAND[1] + 0.6)
    ax.set_ylim(0, top)
    ax.set_xlim(_X0 * 0.8, f_phys * 1.5)
    ax.legend(fontsize=8, loc="upper left", ncol=1)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    _FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG, dpi=140)
    print(f"wrote {_FIG}")
    for f in fx:
        m = float(np.mean([v[0] for v in agg[f]]))
        e = any(v[1] for v in agg[f])
        print(f"  f_active={f:>7.2f} nN  A/A0_core={m:.2f}  ejected={e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
