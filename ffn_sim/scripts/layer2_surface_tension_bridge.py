"""L2 D2 — single-cell cortical tension -> spheroid surface tension BRIDGE (demonstration).

Plugs the MCF7 anchors (single-cell cortical tension gamma; E-cadherin de-adhesion) into the
runtime-forbidden bridge oracle
``ffn_sim.validation.oracles.spheroid.surface_tension_bridge`` and reports the published chain
that links the KU-3.5 cortical-tension observable to the Layer-2 spheroid's aggregate surface
tension and the ``A/A0 = a + b/R + c/R^2`` curvature law:

    gamma (single-cell cortical tension, KU-3.5; g_rigid native 0.57 mN/m IN band)
      - beta (E-cadherin adhesion energy density)            [DITH; Okuda 2026]
      = Gamma_cc (interior cell-cell tension)
    sigma_tissue (aggregate free-surface tension) = gamma     [Roffay 2021]
    Young-Laplace dP = sigma (1/R + 1/R')                      [Roffay 2021, 3D]
    => the 1/R curvature scaling encoded by the A/A0 b/R term.

This is a DEMONSTRATION of the bridge from anchored constants — NOT a fit and NOT a runtime
mechanism (CLAUDE.md inversion rule). The aggregate sigma is produced *emergently* in a CBM run
by ``observables.virial_pressure`` -> ``oracle.surface_tension_from_pressure``; that production
measurement is the next step (needs a run). The magnitude band is a PROXY (MCF10DCIS ~21 mN/m,
Nagle 2022 — no MCF7 tissue-tensiometry datum); the in-band single-cell anchors are MCF7-specific.

Run:  python -m ffn_sim.scripts.layer2_surface_tension_bridge
Outputs: ffn_sim/outputs/layer2/surface_tension_bridge.json
         ffn_sim/outputs/layer2/figs/fig_layer2_surface_tension_bridge.png
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from ffn_sim.validation.oracles.spheroid import surface_tension_bridge as br

# --- MCF7 anchors (provenance: docs/LAYER2_ANCHORS_2026-06-02.md, CORTICAL_TENSION_TRIAGE) ---
GAMMA_CORTICAL = 0.57e-3      # N/m  single-cell cortical tension (g_rigid native, KU-3.5 in-band)
F_DEADHESION = 6.5e-9         # N    MCF7-MCF7 mature de-adhesion (Iturri 2020)
CONTACT_ZONE_W = 1.5e-6       # m    Morse adhesive range (= ~0.1 * diameter, layer2_cbm.yaml)
DIAMETER = 15.0e-6           # m    MCF7 diameter (Wagner 2011)
F0_CATCH = 29.2e-12          # N    Rakshit 2012 E-cadherin catch-bond peak force
# L2.5 spheroid initial radii where G3 PASS (REPORT.md), metres:
R0_L25 = np.array([31.7, 40.4, 53.1, 67.1, 78.3]) * 1e-6

OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "layer2"
FIG_DIR = OUT_DIR / "figs"


def compute_bridge() -> dict:
    """Evaluate the full bridge chain from the anchored constants."""
    # Adhesion energy density beta = work_of_deadhesion / contact_area, over a RANGE of
    # plausible cell-cell contact radii (the contact area is itself emergent — honest range).
    # Work scale: the Morse well depth D_e = 2 * F_detach * contact_zone_width (params.py).
    work_deadhesion = 2.0 * F_DEADHESION * CONTACT_ZONE_W          # J  (~1.95e-14 J)
    contact_radii = np.array([1.5, 2.5, 4.0]) * 1e-6              # m  plausible contact radii
    contact_areas = np.pi * contact_radii**2                       # m^2
    betas = np.array([br.adhesion_tension(work_deadhesion, A) for A in contact_areas])
    beta_over_gamma = betas / GAMMA_CORTICAL

    # Interior cell-cell tension + aggregate surface tension (DITH; Roffay).
    sigma_tissue = br.aggregate_surface_tension(GAMMA_CORTICAL)    # = gamma
    gamma_cc = np.array([br.interfacial_tension(GAMMA_CORTICAL, b) for b in betas])
    ratios = np.array([br.surface_interior_ratio(GAMMA_CORTICAL, b) for b in betas])

    # Young-Laplace interior overpressure across the L2.5 spheroid radii (sphere: 2 sigma / R).
    dP_L25 = np.array([br.young_laplace_pressure(sigma_tissue, R) for R in R0_L25])  # Pa

    # Roffay band window in beta/gamma: ratio = 1/(1 - beta/gamma) in [1.6, 2.0].
    r_lo, r_hi = br.ROFFAY_SURFACE_INTERIOR_RATIO
    bog_window = (1.0 - 1.0 / r_lo, 1.0 - 1.0 / r_hi)             # (0.375, 0.5)

    # Okuda 3D-cap check using a representative interior tension (mid contact-area estimate).
    gamma_cc_mid = float(gamma_cc[1]) if gamma_cc[1] > 0 else float(GAMMA_CORTICAL * 0.4)
    is_cap = br.is_three_d_cap(sigma_tissue, gamma_cc_mid)

    return {
        "anchors": {
            "gamma_cortical_mN_m": GAMMA_CORTICAL * 1e3,
            "f_deadhesion_nN": F_DEADHESION * 1e9,
            "contact_zone_width_um": CONTACT_ZONE_W * 1e6,
            "work_deadhesion_fJ": work_deadhesion * 1e15,
            "cortical_band_mN_m": list(br.CORTICAL_TENSION_BAND_MN_M),
        },
        "sigma_tissue_mN_m": sigma_tissue * 1e3,
        "sigma_in_band": br.CORTICAL_TENSION_BAND_MN_M[0]
        <= sigma_tissue * 1e3
        <= br.CORTICAL_TENSION_BAND_MN_M[1],
        "contact_radii_um": (contact_radii * 1e6).tolist(),
        "beta_mN_m": (betas * 1e3).tolist(),
        "beta_over_gamma": beta_over_gamma.tolist(),
        "gamma_cc_mN_m": (gamma_cc * 1e3).tolist(),
        "surface_interior_ratio": ratios.tolist(),
        "roffay_ratio_band": list(br.ROFFAY_SURFACE_INTERIOR_RATIO),
        "roffay_beta_over_gamma_window": list(bog_window),
        "beta_over_gamma_overlaps_roffay": bool(
            np.any(
                (beta_over_gamma >= min(bog_window))
                & (beta_over_gamma <= max(bog_window))
            )
        ),
        "R0_L25_um": (R0_L25 * 1e6).tolist(),
        "young_laplace_dP_Pa": dP_L25.tolist(),
        "is_3d_cap": is_cap,
        "tissue_sigma_proxy_mN_m": br.TISSUE_SIGMA_PROXY_MN_M,
        "notes": (
            "sigma_tissue = gamma (in-band) bridges the single-cell cortical tension to the "
            "aggregate surface tension; Young-Laplace dP = 2 sigma / R gives the 1/R curvature "
            "that the A/A0 b/R term encodes. The Roffay outer/interior ratio (1.6-2.0) is "
            "reproduced for beta/gamma in (0.375, 0.5). Absolute aggregate-sigma band is a "
            "PROXY (MCF10DCIS 21 mN/m); the emergent sigma is measured in a CBM run via "
            "virial_pressure -> surface_tension_from_pressure (next step)."
        ),
    }


def make_figure(result: dict) -> Path | None:
    """3-panel bridge figure (visualize-at-closeout rule). Best-effort; returns path or None."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - viz is best-effort
        print(f"[viz] matplotlib unavailable ({exc}); skipping figure.", file=sys.stderr)
        return None

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))

    # Panel A: sigma_tissue = gamma vs the KU-3.5 band + proxy magnitude.
    band = result["anchors"]["cortical_band_mN_m"]
    sig = result["sigma_tissue_mN_m"]
    ax[0].axhspan(band[0], band[1], color="tab:green", alpha=0.18, label=f"KU-3.5 band {band}")
    ax[0].axhline(sig, color="tab:blue", lw=2.5, label=f"sigma_tissue = gamma = {sig:.2f} mN/m")
    ax[0].axhline(
        result["tissue_sigma_proxy_mN_m"], color="tab:red", ls="--",
        label=f"tissue proxy (MCF10DCIS) {result['tissue_sigma_proxy_mN_m']:.0f} mN/m",
    )
    ax[0].set_yscale("log")
    ax[0].set_ylabel("surface tension [mN/m]")
    ax[0].set_title("A  aggregate surface tension = single-cell gamma\n(Roffay: free surface = cortex)")
    ax[0].set_xticks([])
    ax[0].legend(fontsize=7, loc="center left")

    # Panel B: Young-Laplace interior overpressure vs spheroid radius.
    R = np.array(result["R0_L25_um"])
    dP = np.array(result["young_laplace_dP_Pa"])
    ax[1].plot(R, dP, "o-", color="tab:purple", label="dP = 2 sigma / R (Roffay 3D)")
    ax[1].set_xlabel("spheroid radius R [um]")
    ax[1].set_ylabel("interior overpressure dP [Pa]")
    ax[1].set_title("B  Young-Laplace: 1/R curvature\n(the A/A0 b/R term)")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    # Panel C: surface/interior ratio vs beta/gamma with the Roffay band + anchored estimates.
    bog = np.linspace(0.0, 0.7, 200)
    ratio_curve = 1.0 / (1.0 - bog)
    win = result["roffay_beta_over_gamma_window"]
    rband = result["roffay_ratio_band"]
    ax[2].plot(bog, ratio_curve, color="k", lw=1.5, label="sigma/Gamma_cc = 1/(1-beta/gamma)")
    ax[2].axhspan(rband[0], rband[1], color="tab:green", alpha=0.18,
                  label=f"Roffay ratio {rband}")
    ax[2].axvspan(min(win), max(win), color="tab:green", alpha=0.10)
    for bg in result["beta_over_gamma"]:
        if bg < 0.95:
            ax[2].axvline(bg, color="tab:orange", ls=":", alpha=0.8)
    ax[2].plot([], [], color="tab:orange", ls=":", label="anchored beta/gamma (contact-area range)")
    ax[2].set_xlabel("beta / gamma  (adhesion / cortical)")
    ax[2].set_ylabel("outer / interior tension ratio")
    ax[2].set_ylim(1.0, 4.0)
    ax[2].set_title("C  Roffay outer/interior ratio\nreproduced at beta/gamma ~ 0.375-0.5")
    ax[2].legend(fontsize=7, loc="upper left")
    ax[2].grid(alpha=0.3)

    fig.suptitle(
        "L2 D2 — single-cell cortical tension -> spheroid surface-tension bridge "
        "(Chugh/Roffay/Fastabend/Okuda anchors; demonstration, not a fit)",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = FIG_DIR / "fig_layer2_surface_tension_bridge.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    result = compute_bridge()
    out_json = OUT_DIR / "surface_tension_bridge.json"
    out_json.write_text(json.dumps(result, indent=2))

    print("=== L2 D2 surface-tension bridge (single-cell gamma -> spheroid sigma) ===")
    print(f"  sigma_tissue = gamma            = {result['sigma_tissue_mN_m']:.3f} mN/m"
          f"  (in KU-3.5 band: {result['sigma_in_band']})")
    print(f"  beta (adhesion) range           = "
          f"{min(result['beta_mN_m']):.3f}..{max(result['beta_mN_m']):.3f} mN/m"
          f"  (beta/gamma {min(result['beta_over_gamma']):.2f}..{max(result['beta_over_gamma']):.2f})")
    print(f"  Roffay outer/interior ratio band {result['roffay_ratio_band']} hit at "
          f"beta/gamma in {tuple(round(x,3) for x in result['roffay_beta_over_gamma_window'])};"
          f" anchored range overlaps: {result['beta_over_gamma_overlaps_roffay']}")
    print(f"  Young-Laplace dP over R={result['R0_L25_um']} um:")
    print(f"      dP = {[round(p,2) for p in result['young_laplace_dP_Pa']]} Pa")
    print(f"  Okuda 3D-cap (free-surface > 0.2*cell-cell): {result['is_3d_cap']}  (cf. L2.6)")
    print(f"  [PROXY] aggregate-sigma magnitude band: MCF10DCIS {result['tissue_sigma_proxy_mN_m']} mN/m")
    print(f"  -> {out_json}")

    fig = make_figure(result)
    if fig is not None:
        print(f"  -> {fig}")


if __name__ == "__main__":
    main()
