"""Figure: active-motility UNJAMMING of the DCM aggregate — the jamming→unjamming transition.

The DCM spheroid is JAMMED (cells don't rearrange → confined, A/A0≈1). Per-cell active self-propulsion
(the Self-Propelled-Voronoi active-matter lever, Bi-Yang-Marchetti-Manning) FLUIDISES it: as F_active rises,
the aggregate spreads (A/A0 up) AND the per-cell 3D shape index s = S/V^(2/3) rises toward the Merkel-Manning
3D rigidity threshold s0*≈5.41 — the jamming→unjamming order parameter. Confined s=4.93 (jammed, < s0*);
motility drives it to 5.35 (approaching the fluid transition) = collective spreading. All GENUINE COLLECTIVE
SPREAD (peel-aware spread_eval, V/V0=1), physiological F_active (10-100 nN single-cell traction, KB-2.12).

Data measured 2026-07-10 (native N=100 on gbook A5000, visual-verified). Numbers are hard-coded from the
committed sweep runs (s3_motility_fa*) — see _historical/DCM_UNJAMMING_ACTIVE_MOTILITY_2026-07-10.md.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = "aleph/outputs/h_dcm_two_stage/figs"
os.makedirs(OUT, exist_ok=True)

# measured sweep (F_active nN, shape_index s, A/A0 topdown, drift µm, peel_idx)
fa = np.array([0.0, 10.0, 50.0, 100.0])
s = np.array([4.934, 4.988, 5.200, 5.351])
aa0 = np.array([1.006, 1.118, 1.556, 1.924])
drift = np.array([0.08, 1.02, 2.60, 3.80])
peel = np.array([np.nan, 3.4, 1.4, 1.1])
S0_STAR = 5.41   # Merkel-Manning 3D rigidity (unjamming) threshold
S_SPHERE = 4.836


def main():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))

    # (L) shape index vs F_active — the jamming order parameter
    a1.axhspan(S_SPHERE, S0_STAR, color="tab:blue", alpha=0.08)
    a1.axhline(S0_STAR, color="tab:red", ls="--", lw=1.5)
    a1.axhline(S_SPHERE, color="gray", ls=":", lw=1.0)
    a1.text(2, S0_STAR + 0.008, "s0* ≈ 5.41  (Merkel–Manning 3D unjamming)", color="tab:red", fontsize=9)
    a1.text(2, S_SPHERE - 0.02, "sphere 4.84", color="gray", fontsize=8)
    a1.text(2, 5.0, "JAMMED\n(solid)", color="tab:blue", fontsize=10, va="center")
    a1.plot(fa, s, "-o", color="tab:blue", lw=2.2, ms=7)
    a1.set_xlabel("active self-propulsion  F_active  [nN/cell]")
    a1.set_ylabel("per-cell 3D shape index  s = S/V^(2/3)", color="tab:blue")
    a1.set_title("Active motility raises the shape index toward s0*\n= the jamming→unjamming transition")
    a1.set_ylim(4.80, 5.46)
    a1.grid(alpha=0.25)

    # (R) the spread (A/A0) + drift, colored by GENUINE
    a1b = a2
    l1 = a1b.plot(fa, aa0, "-o", color="tab:green", lw=2.2, ms=7, label="A/A0 (top-down spread)")
    a1b.axhline(1.0, color="gray", ls=":", lw=1.0)
    a1b.set_xlabel("active self-propulsion  F_active  [nN/cell]")
    a1b.set_ylabel("A/A0  (top-down silhouette spread)", color="tab:green")
    a1b.set_title("Collective SPREAD emerges (all GENUINE, peel-aware)\nA/A0 1.0→1.92, drift 0.08→3.8 µm")
    a2b = a1b.twinx()
    l2 = a2b.plot(fa, drift, "-s", color="tab:orange", lw=1.8, ms=6, label="COM-relative drift [µm]")
    a2b.set_ylabel("cell drift [µm]", color="tab:orange")
    for x, y, p in zip(fa[1:], aa0[1:], peel[1:]):
        a1b.annotate(f"peel {p:.1f}\n(GENUINE)", (x, y), textcoords="offset points", xytext=(-8, 8),
                     fontsize=7, color="darkgreen")
    a1b.legend(l1 + l2, [h.get_label() for h in l1 + l2], fontsize=8, loc="upper left")
    a1b.grid(alpha=0.25)

    fig.tight_layout()
    path = f"{OUT}/dcm_unjamming_active_motility_doseresponse.png"
    fig.savefig(path, dpi=130)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
