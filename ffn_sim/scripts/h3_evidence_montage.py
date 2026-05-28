#!/usr/bin/env python
"""H.3 evidence montage — one-page overview of cortex unit's deliverables.

Assembles a 3×3 grid PNG combining existing key figures into a single
review-ready overview for PI ✅ DONE evidence:

  Row 1: structural validation
    [1.1] L_p FULL bending equipartition + KU-3.20 angle PDF + topology
  Row 2: dynamic subsystems
    [2.1] crosslinker network + myosin minifilament structure + force
          constants
  Row 3: full assembly + closure
    [3.1] full-cell 3D structure (this loop) + motor engagement trace +
          (placeholder) KU-3.5 tension sweep result

Designed for printing as a single landscape-A4 / 16:9 slide.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

PKG = Path(__file__).resolve().parents[1]
FIGS = PKG / "outputs" / "h3" / "figs"
OUT_PATH = FIGS / "fig_h3_evidence_montage.png"


# Source figures (path, caption). If a path doesn't exist, render a
# placeholder panel with the caption text.
SOURCES = [
    ("fig_h3_lp_bending_equipartition.png",
     "L_p MEDIUM ⟨E_bend⟩=0.99 kT (KU-1.1)"),
    ("fig_h3_lp_angle_pdf.png",
     "Angle PDF vs 3D Boltzmann (KS gate)"),
    ("fig_h3_topology_3d.png",
     "Cortex topology — spherical-shell ×40 mesoscopic"),
    ("fig_h3_full_bond_network.png",
     "Bond-network: cortex + xlink + myosin + ERM"),
    ("fig_h3_force_constants.png",
     "Force constants — Cortex spring/angle + xlink/myosin"),
    ("fig_h3_lp_filament_configs.png",
     "L_p FULL filament configurations"),
    ("fig_h3_cell_structure.png",
     "Full assembly cell structure (this loop)"),
    ("fig_h3_cell_engagement_trace.png",
     "Motor binding kinetics (this loop)"),
    ("fig_h3_ku35_tension_sweep.png",
     "KU-3.5 tension sweep (gbook canonical — landing later)"),
]


def main() -> None:
    fig = plt.figure(figsize=(16, 14))
    for i, (fname, caption) in enumerate(SOURCES):
        ax = fig.add_subplot(3, 3, i + 1)
        ax.set_xticks([]); ax.set_yticks([])
        path = FIGS / fname
        if path.exists():
            img = mpimg.imread(path)
            ax.imshow(img)
            ax.set_title(caption, fontsize=10)
        else:
            ax.set_facecolor("#f4f4f4")
            ax.text(0.5, 0.5, f"(pending)\n{caption}",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=11, color="#888")
            ax.set_title(caption, fontsize=10, color="#888")
    fig.suptitle(
        "H.3 cortex evidence montage — one-page review (Phase 1 unit deliverables)",
        fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(OUT_PATH, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE_EVIDENCE {OUT_PATH}", flush=True)


if __name__ == "__main__":
    main()
