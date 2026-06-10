"""H.7 — FA molecular-clutch ON vs OFF contrast for the dynamic spreading run.

Overlays two ``h7_spreading_dynamics.py`` curve.json runs (identical protrusion
drive, clutch tether ON vs OFF) to isolate the FA clutch's mechanical role in
spreading. The clutch (BasalAdhesionTether = the engaged integrin FA ensemble)
holds the protruded lamellipodial sheet on the substrate; with it OFF the cortex
/ membrane line tension lifts and retracts the front out of the basal plane and
the footprint stalls / collapses — the Chan-Odde / Elosegui-Artola motor-clutch
signature that adhesion is required to convert protrusion into spread area.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_AA0_BAND = (2.0, 4.0)
_MOCK_AA0 = 1.06


def _load(path):
    with open(path) as fh:
        d = json.load(fh)
    s = d["series"]
    t = np.array([r.get("t_phys_min", r["t_eff_s"] / 60.0) for r in s])
    A0 = d["meta"]["A0_um2"]
    aa0 = np.array([r["A_um2"] for r in s]) / A0
    rf = np.array([r["r_front_um"] for r in s])
    return t, aa0, rf, d["meta"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--on", required=True, help="clutch-ON curve.json")
    ap.add_argument("--off", required=True, help="clutch-OFF curve.json")
    base = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "figs"
    ap.add_argument("--out", default=str(base / "h7_spreading_clutch_compare.png"))
    args = ap.parse_args()

    ton, aon, rfon, mon = _load(args.on)
    toff, aoff, rfoff, moff = _load(args.off)

    fig, ax = plt.subplots(1, 2, figsize=(13.0, 5.2), constrained_layout=True)

    a = ax[0]
    a.axhspan(*_AA0_BAND, color="#bfe3bf", alpha=0.5,
              label=f"physiological band [{_AA0_BAND[0]:.0f}, {_AA0_BAND[1]:.0f}]")
    a.axhline(_MOCK_AA0, color="0.5", ls=":", lw=1.5,
              label=f"old kinematic mock ≈{_MOCK_AA0:.2f} (flat)")
    a.plot(ton, aon, "-o", color="#c0392b", lw=2.4, ms=5,
           label="FA clutch ON (adhered)")
    a.plot(toff, aoff, "-s", color="#2f6fb0", lw=2.0, ms=4,
           label="FA clutch OFF (front lifts)")
    a.set_xlabel("physiological spreading time (min)")
    a.set_ylabel("A / A₀  (basal contact footprint)")
    a.set_title("FA molecular clutch sets the spread area")
    a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="upper left")

    a = ax[1]
    a.plot(ton, rfon, "-o", color="#c0392b", lw=2.4, ms=5, label="clutch ON")
    a.plot(toff, rfoff, "-s", color="#2f6fb0", lw=2.0, ms=4, label="clutch OFF")
    a.set_xlabel("physiological spreading time (min)")
    a.set_ylabel("leading-edge radius r_front (µm)")
    a.set_title("leading-edge advance")
    a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="upper left")

    fig.suptitle(
        f"H.7 dynamic spreading — FA molecular-clutch ON vs OFF  "
        f"(A/A₀ final: ON {aon[-1]:.2f} vs OFF {aoff[-1]:.2f}; "
        f"identical protrusion drive)",
        fontsize=12, fontweight="bold")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    plt.close(fig)
    print(f"[compare] wrote {args.out}")
    print(f"[compare] clutch ON  A/A0_final={aon[-1]:.3f}  r_front={rfon[-1]:.2f}µm")
    print(f"[compare] clutch OFF A/A0_final={aoff[-1]:.3f}  r_front={rfoff[-1]:.2f}µm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
