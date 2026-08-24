"""Figure: the SPV persistence (τ_p) axis of the DCM unjamming transition.

The Self-Propelled-Voronoi (Bi-Yang-Marchetti-Manning) jamming transition is set by BOTH the self-propulsion
force/speed AND the directional persistence τ_p (=1/Dr). At fixed F_active=100 nN, sweeping τ_p: below a
critical persistence (~100 s) the motility is effectively DIFFUSIVE (fast reorientation → net displacement
cancels) and the tissue stays JAMMED (shape index ≈ 4.92, the no-motility confined value); above it (τ_p ≥
150 s, incl. the physiological ~10 min = 600 s) the motion is ballistic/directed → the tissue UNJAMS and the
shape index reaches s0*≈5.41. Together with the F_active dose-response, this maps both SPV axes on the DCM.

Data measured 2026-07-10 (native N=100, gbook A5000, seed 11; τ_p=30 s visually confirmed jammed).
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = "aleph/outputs/h_dcm_two_stage/figs"
os.makedirs(OUT, exist_ok=True)

tau = np.array([30.0, 150.0, 600.0, 3000.0])
aa0 = np.array([1.04, 2.58, 2.19, 2.39])
s = np.array([4.922, 5.397, 5.408, 5.411])
frac = np.array([0.02, 0.42, 0.45, 0.39])
S0 = 5.41
S_SPH = 4.836


def main():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))

    a1.axhline(S0, color="tab:red", ls="--", lw=1.5)
    a1.text(35, S0 + 0.008, "s0*≈5.41 (unjamming)", color="tab:red", fontsize=9)
    a1.axhline(S_SPH, color="gray", ls=":", lw=1.0)
    a1.text(35, S_SPH - 0.02, "sphere / confined 4.84–4.93 (JAMMED)", color="gray", fontsize=8)
    a1.axhspan(S_SPH, S0, color="tab:blue", alpha=0.07)
    a1.semilogx(tau, s, "-o", color="tab:blue", lw=2.2, ms=8)
    a1.axvline(600, color="green", ls=":", alpha=0.6)
    a1.text(620, 5.0, "physiological\n~10 min", color="green", fontsize=8)
    a1.set_xlabel("directional persistence  τ_p  [s]  (log)")
    a1.set_ylabel("per-cell 3D shape index  s")
    a1.set_title("SPV persistence axis (F_active=100 nN)\nlow τ_p = diffusive = JAMMED; τ_p≥150 s → s0*")
    a1.set_ylim(4.80, 5.46)
    a1.grid(True, which="both", alpha=0.25)

    a2.semilogx(tau, aa0, "-s", color="tab:green", lw=2.2, ms=8)
    a2.axhline(1.0, color="gray", ls=":")
    a2.axvline(600, color="green", ls=":", alpha=0.6)
    a2.set_xlabel("directional persistence  τ_p  [s]  (log)")
    a2.set_ylabel("A/A0  (aggregate spread)")
    a2.set_title("Spread requires PERSISTENT motility\nτ_p=30 s → A/A0=1.04 (jammed), ≥150 s → 2.2–2.6")
    a2.grid(True, which="both", alpha=0.25)

    fig.suptitle("DCM unjamming — the SPV persistence axis (a jammed tissue needs DIRECTED, not diffusive, motility)",
                 fontsize=11)
    fig.tight_layout()
    path = f"{OUT}/dcm_unjamming_spv_persistence.png"
    fig.savefig(path, dpi=130)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
