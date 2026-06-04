"""Motility-bridge diagnostic figure: localise the magnitude gap to the active side, and show
the protrusion-anchored traction EXCEEDS cohesion (the CBM 1-particle structural limit).

Pure-arithmetic (no sim): reads the resolved Layer-2 params + the lamellipodium→CBM bridge.
One entry point: ``python -m ffn_sim.scripts.layer2_motility_bridge_vis``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml

from ffn_sim.spheroid.motility_bridge import resolve_active_traction
from ffn_sim.spheroid.params import resolve_layer2

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "layer2"
_FIG = _OUT / "figs" / "fig_layer2_motility_bridge.png"


def main() -> int:
    resolved = resolve_layer2(yaml.safe_load((_OUT.parents[1] / "configs" / "layer2_cbm.yaml").read_text()))
    m = resolve_active_traction(resolved)
    cohesion = float(resolved.deadhesion_force_mature)  # N (Iturri-2020 de-adhesion)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.0, 5.4))

    # ---- Panel A: the force ladder ----
    levels = [
        ("whole-cell anchor\n(MCF7 19 µm/h)", m.f_wholecell * 1e9, "tab:green"),
        ("B1 stable ceiling", m.ceiling * 1e9, "0.5"),
        ("cohesion F_detach\n(Iturri 6.5 nN)", cohesion * 1e9, "tab:red"),
        ("protrusion anchor\n(Bieling 1.88 µm/min)", m.f_protrusion * 1e9, "tab:blue"),
    ]
    for i, (lab, val, c) in enumerate(levels):
        axA.barh(i, val, color=c, alpha=0.85)
        axA.text(val + 0.15, i, f"{val:.1f} nN", va="center", fontsize=9)
    axA.set_yticks(range(len(levels)))
    axA.set_yticklabels([l for l, _, _ in levels], fontsize=8)
    axA.axvspan(0, m.ceiling * 1e9, color="tab:green", alpha=0.06)
    axA.axvspan(cohesion * 1e9, m.f_protrusion * 1e9 + 1, color="tab:red", alpha=0.06)
    axA.text(cohesion * 1e9 + 0.1, 3.35,
             "detachment regime\n(active > cohesion → cells tear apart,\nnot spread — CBM 1-particle limit)",
             color="tab:red", fontsize=7.5, va="top")
    axA.set_xlabel("per-cell force (nN)")
    axA.set_title("A — active-traction force ladder\nprotrusion anchor (9.4 nN) > cohesion (6.5 nN) > ceiling (3 nN) > whole-cell (1.6 nN)")
    axA.set_xlim(0, m.f_protrusion * 1e9 + 1.5)
    axA.grid(axis="x", alpha=0.25)

    # ---- Panel B: the ratio = the magnitude gap ----
    axB.bar(["whole-cell\n1.6 nN", "protrusion\n9.4 nN"],
            [m.f_wholecell * 1e9, m.f_protrusion * 1e9],
            color=["tab:green", "tab:blue"], alpha=0.85)
    axB.axhline(m.ceiling * 1e9, color="0.5", ls=":", label="B1 ceiling 3 nN")
    axB.axhline(cohesion * 1e9, color="tab:red", ls="--", label="cohesion 6.5 nN")
    axB.annotate("", xy=(1, m.f_protrusion * 1e9), xytext=(1, m.f_wholecell * 1e9),
                 arrowprops=dict(arrowstyle="<->", color="k"))
    axB.text(1.08, (m.f_protrusion + m.f_wholecell) * 1e9 / 2,
             f"×{m.active_ratio:.1f}\n≈ the observed\nA/A0 magnitude\ngap (~5–9×)", fontsize=8.5, va="center")
    axB.set_ylabel("anchored f_active (nN)")
    axB.set_title("B — the two single-cell anchors differ by ~6×\n= the magnitude gap, localised to the ACTIVE side\n(cohesion is anchored 2 ways → NOT the cause)")
    axB.legend(fontsize=8, loc="upper left")
    axB.grid(axis="y", alpha=0.25)

    fig.tight_layout()
    _FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG, dpi=140)
    print(f"wrote {_FIG}")
    print(f"f_wholecell={m.f_wholecell*1e9:.2f} nN  ceiling={m.ceiling*1e9:.1f}  "
          f"cohesion={cohesion*1e9:.1f}  f_protrusion={m.f_protrusion*1e9:.2f} nN  ratio={m.active_ratio:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
