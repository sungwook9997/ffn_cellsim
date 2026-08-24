#!/usr/bin/env python3
r"""Figures for the Oracle-B bleb-growth driver (:mod:`aleph.components.incumbent.bleb_growth`).

Two modes, both CPU-only (matplotlib, no CUDA):

* ``--npz PATH`` — the REAL growth figure from a completed gbook A5000 run: patch vs. control outward
  displacement over time, and the cap ΔP with the Laplace critical radius ``a_crit(t)=2γ/ΔP`` overlaid. Run
  this after the native growth run lands its ``growth.npz`` (once (b) resting convergence is in place).
* no ``--npz`` — the ANALYSIS-VALIDATION figure: a synthetic area-uniform membrane sphere with a known cap
  bulge, showing that the host analysis pipeline (:func:`cap_membership` / :func:`patch_outward_um`) recovers
  the imposed displacement. This validates the ANALYSIS only; it makes no claim about the physical bleb.

Visualization-integrity: SI/engine units annotated, the resting baseline shown as the zero reference, no axis
truncation. Writes into ``aleph/outputs/ac/bleb/figs/``.
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

matplotlib_backend = "Agg"


def _synthetic_analysis_figure(out_png: str) -> None:
    """Render the CPU analysis-validation figure from synthetic geometry (no CUDA, no native run)."""
    import matplotlib

    matplotlib.use(matplotlib_backend)
    import matplotlib.pyplot as plt

    from aleph.components.incumbent.bleb_growth import cap_membership, patch_outward_um
    from aleph.components.incumbent.bleb_perturbation import cap_cos_half_angle

    radius, centroid = 7.5, np.array([1.0, -2.0, 0.5])
    cos_theta = np.linspace(-1.0 + 1e-3, 1.0 - 1e-3, 200)
    phi = np.linspace(0.0, 2.0 * np.pi, 200, endpoint=False)
    cc, pp = np.meshgrid(cos_theta, phi, indexing="ij")
    ss = np.sqrt(1.0 - cc**2)
    unit = np.stack([ss * np.cos(pp), ss * np.sin(pp), cc], axis=-1).reshape(-1, 3)
    pos = unit * radius + centroid
    axis = np.array([0.0, 0.0, 1.0])
    cos_half = cap_cos_half_angle(30.0)
    in_cap, base_r = cap_membership(pos, centroid, axis, cos_half)
    base_patch_r, base_ctrl_r = float(base_r[in_cap].mean()), float(base_r[~in_cap].mean())

    imposed = np.linspace(0.0, 1.0, 11)
    recovered_patch = [patch_outward_um(np.linalg.norm(
        (pos + unit * (in_cap[:, None] * b)) - centroid, axis=1), base_patch_r, in_cap) for b in imposed]
    recovered_ctrl = [patch_outward_um(np.linalg.norm(
        (pos + unit * (in_cap[:, None] * b)) - centroid, axis=1), base_ctrl_r, ~in_cap) for b in imposed]

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(11, 4.4))
    axl.scatter(pos[~in_cap, 0], pos[~in_cap, 2], s=2, c="#4477aa", label="control (tethered)")
    axl.scatter(pos[in_cap, 0], pos[in_cap, 2], s=2, c="#cc3311", label=f"cap θ½=30° ({in_cap.mean()*100:.1f}%)")
    axl.set_aspect("equal")
    axl.set_xlabel("x [µm]")
    axl.set_ylabel("z [µm]")
    axl.set_title("area-uniform membrane sphere\n(+z ablation cap)")
    axl.legend(loc="lower right", fontsize=8)

    axr.plot(imposed, recovered_patch, "o-", c="#cc3311", label="patch (recovered)")
    axr.plot(imposed, recovered_ctrl, "s-", c="#4477aa", label="control (recovered)")
    axr.plot(imposed, imposed, "k--", lw=1, label="imposed (y=x)")
    axr.set_xlabel("imposed cap bulge [µm]")
    axr.set_ylabel("recovered outward displacement [µm]")
    axr.set_title("analysis pipeline recovers the imposed bulge\n(cap only; control ≈ 0)")
    axr.legend(fontsize=8)
    axr.grid(alpha=0.3)
    fig.suptitle("Oracle-B bleb-growth ANALYSIS validation (synthetic; CPU) — not a physical bleb claim",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    print(f"[bleb-vis] wrote analysis-validation figure {out_png}", flush=True)


def _growth_figure(npz_path: str, out_png: str) -> None:
    """Render the real growth figure from a completed native run's ``growth.npz``."""
    import matplotlib

    matplotlib.use(matplotlib_backend)
    import matplotlib.pyplot as plt

    d = np.load(npz_path, allow_pickle=True)
    report = json.loads(str(d["report_json"])) if "report_json" in d else {}
    t = d["t"]
    patch = d["patch_disp"]
    ctrl = d["ctrl_disp"]
    dp = d["dp_cap"]
    a_crit = d["a_crit"]
    a_patch = d["a_patch"]
    deadhered_step = int(d["deadhered_step"]) if "deadhered_step" in d else -1

    fig, (axu, axl) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    axu.axhline(0.0, c="k", lw=0.8, label="resting baseline")
    axu.plot(t, patch, "o-", c="#cc3311", label="patch (de-adhered cap)")
    axu.plot(t, ctrl, "s-", c="#4477aa", label="control (intact membrane)")
    if 0 <= deadhered_step < len(t):
        axu.axvline(t[deadhered_step], c="#228833", ls=":", label="de-adhesion (accepted)")
    axu.set_ylabel("outward displacement [µm]")
    axu.set_title(f"Bleb growth after commit-safe de-adhesion  "
                  f"(grows={report.get('bleb_grows')}, bulge={report.get('final_bulge_patch_minus_control_um', 0):.2e} µm)")
    axu.legend(fontsize=8)
    axu.grid(alpha=0.3)

    axl.plot(t, dp, "o-", c="#ee7733", label="cap ΔP [Pa]")
    axl.set_ylabel("cap ΔP [Pa]", color="#ee7733")
    axl.tick_params(axis="y", labelcolor="#ee7733")
    axl.set_xlabel("physical time [s]")
    axl.grid(alpha=0.3)
    axr = axl.twinx()
    axr.plot(t, a_crit, "^-", c="#aa3377", label="a_crit=2γ/ΔP [µm]")
    axr.plot(t, a_patch, "--", c="#000000", label="a_patch [µm]")
    axr.set_ylabel("radius [µm]")
    lines = axl.get_lines() + axr.get_lines()
    axl.legend(lines, [ln.get_label() for ln in lines], fontsize=8, loc="best")
    axl.set_title("Local ΔP non-equilibration + Laplace balance (bulges while a_patch > a_crit)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    print(f"[bleb-vis] wrote growth figure {out_png}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description="Oracle-B bleb-growth figures.")
    p.add_argument("--npz", type=str, default="", help="completed growth.npz (real figure); omit for analysis fig")
    p.add_argument("--out", type=str, default="")
    a = p.parse_args()
    figs_dir = os.path.join("aleph", "outputs", "ac", "bleb", "figs")
    os.makedirs(figs_dir, exist_ok=True)
    if a.npz:
        _growth_figure(a.npz, a.out or os.path.join(figs_dir, "bleb_growth.png"))
    else:
        _synthetic_analysis_figure(a.out or os.path.join(figs_dir, "bleb_growth_analysis_validation.png"))


if __name__ == "__main__":
    main()
