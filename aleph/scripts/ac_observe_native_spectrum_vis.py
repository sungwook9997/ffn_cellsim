#!/usr/bin/env python
r"""Figures for the native cortex spectrum record — this lane's own entry point, no GPU, no re-run.

Reads the committed JSON written by ``ac_observe_native_spectrum.py`` and writes three PNGs beside it.
Each carries the reference it is judged against on the SAME axes, because a measurement without its
reference overlaid is not a check (project visualization rule, PI 2026-05-21 / revised 2026-07-28):

  1. ``family_certification.png`` — each stiffness family's share of the assembled tangent, and each
     family's worst entrywise asymmetry against ITS OWN derived round-off floor.  A family that
     contributed exactly nothing is drawn as such rather than omitted, because "absent" and "zero" are
     different findings and only one of them is a coverage hole.
  2. ``native_ritz.png``          — the Ritz values with their rigorous Parlett residual bounds as error
     bars, for ``P K P`` and for ``K``.  The bounds are the point: they show convergence at the
     extremes and its absence in the interior, which is the method's honest shape.
  3. ``explicit_step.png``        — the dimensionless stability product ``dt_mu·λ_max`` at the
     production step and at FIRE's maximum multiplier.  The forward-Euler limit 2 is drawn, and only
     the plain-descent bar is judged against it: FIRE is damped MD with its own differently-scaled
     bound, so its bar is hatched rather than coloured pass/fail.

Usage (dev Mac, no CUDA needed)::

    python aleph/scripts/ac_observe_native_spectrum_vis.py \
        --record aleph/outputs/ac/observe/native_cortex_spectrum/record.json

Visualization integrity: no axis truncation; every log axis says so in its label; units are engine
units (pN/µm, µm/pN) named on every axis; the reference is drawn, not described in a caption.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

#: One palette across the figures, so a colour means the same thing everywhere.
COLOR_PROJECTED = "#3b6ea5"
COLOR_RAW = "#c1502e"
COLOR_REFERENCE = "#6b6b6b"
COLOR_UNEXERCISED = "#b0b0b0"
COLOR_LIMIT = "#a01d1d"


def _family_figure(record: dict[str, Any], out: Path) -> Path:
    """Plot per-family share of the assembled tangent, and per-family symmetry against its own floor."""
    certification = record["measurements"]["family_certification"]
    families = certification["families"]
    names = list(families)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))

    shares = [families[name]["fraction_of_assembled_norm"] or 0.0 for name in names]
    colours = [COLOR_PROJECTED if families[name]["exercised"] else COLOR_UNEXERCISED
               for name in names]
    positions = np.arange(len(names))
    axes[0].barh(positions, shares, color=colours)
    axes[0].set_yticks(positions, names)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("‖K_j‖_F / ‖K‖_F  (dimensionless, log axis)")
    axes[0].set_title("share of the assembled tangent, per stiffness family\n"
                      "grey = isolated tangent exactly zero, i.e. NOT certified by this run")
    for index, name in enumerate(names):
        if not families[name]["exercised"]:
            axes[0].text(min([s for s in shares if s > 0.0] or [1e-12]), index,
                         "  not exercised", va="center", fontsize=8, color=COLOR_UNEXERCISED)

    exercised = [name for name in names if families[name]["exercised"]]
    asymmetry = [families[name]["symmetry"]["max_abs_asymmetry_pN_per_um"] for name in exercised]
    floors = [families[name]["symmetry"]["round_off_floor_pN_per_um"] for name in exercised]
    positions = np.arange(len(exercised))
    # An asymmetry of EXACTLY zero is the best outcome available and a log axis cannot draw it, which
    # would make the strongest result on the panel look like missing data. Park those bars below the
    # smallest drawable value and label them, so "exactly 0" reads as a measurement.
    drawable = [value for value in asymmetry + floors if value > 0.0]
    floor_of_axis = min(drawable) / 100.0 if drawable else 1e-20
    plotted = [value if value > 0.0 else floor_of_axis for value in asymmetry]
    axes[1].barh(positions - 0.19, plotted, height=0.36, color=COLOR_RAW,
                 label="max|K_j − K_jᵀ| measured [pN/µm]")
    axes[1].barh(positions + 0.19, floors, height=0.36, color=COLOR_REFERENCE,
                 label="derived float64 accumulation floor γ_n·max|K_j| [pN/µm]")
    for index, value in enumerate(asymmetry):
        if value == 0.0:
            axes[1].text(floor_of_axis, index - 0.19, "  = 0 exactly",
                         va="center", fontsize=7.5, color=COLOR_RAW)
    axes[1].set_xlim(left=floor_of_axis / 2.0)
    axes[1].set_yticks(positions, exercised)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("entrywise asymmetry [pN/µm] (log axis)")
    axes[1].set_title("does each family admit a potential?\n"
                      "measured bar shorter than its own floor bar = symmetric to round-off\n"
                      "(a bar at the axis floor labelled '= 0 exactly' is symmetric to the last bit)",
                      fontsize=9.5)
    axes[1].legend(fontsize=8, loc="lower right", framealpha=0.95)

    axes[0].set_title(axes[0].get_title(), fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def _ritz_axis(axis: Any, spectrum: dict[str, Any], colour: str, label: str) -> None:
    """Draw one pass's Ritz values with their rigorous residual bounds as error bars."""
    values = np.asarray(spectrum["direct_pass"]["ritz_values_pN_per_um"], float)
    bounds = np.asarray(spectrum["direct_pass"]["residual_bounds_pN_per_um"], float)
    index = np.arange(values.size)
    axis.errorbar(index, values, yerr=bounds, fmt=".", ms=4.0, color=colour, ecolor=colour,
                  elinewidth=0.9, capsize=2.0, label=label)
    axis.set_yscale("log")


def _ritz_figure(record: dict[str, Any], out: Path) -> Path:
    """Plot both operators' Ritz values with the bounds that say which of them are converged."""
    native = record["measurements"].get("native_spectrum")
    if native is None:
        raise SystemExit("[vis] this record has no native stage; nothing to plot here")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))
    for axis, key, colour, title in (
        (axes[0], "projected_PKP", COLOR_PROJECTED,
         "P K P restricted to range(P) — what the inner solve inverts"),
        (axes[1], "unprojected_K", COLOR_RAW,
         "K unprojected — what conservativity is a statement about"),
    ):
        block = native[key]
        _ritz_axis(axis, block["at_sigma_2theta"], colour, "Ritz value ± Parlett bound [pN/µm]")
        lambda_max = block["at_sigma_2theta"]["lambda_max_pN_per_um"]
        axis.axhline(lambda_max, ls="-", lw=1.0, color="black",
                     label=f"λ_max = {lambda_max:.4g} pN/µm "
                           f"({'resolved' if block['at_sigma_2theta']['lambda_max_resolved'] else 'NOT resolved'})")
        lambda_min = block["at_sigma_2theta"]["lambda_min_pN_per_um"]
        resolved = block["at_sigma_2theta"]["lambda_min_resolved"]
        verdict = ("resolved" if resolved
                   else "NOT resolved — quoting it would quote the start vector")
        if lambda_min > 0.0:
            axis.axhline(lambda_min, ls="--", lw=1.0, color=COLOR_REFERENCE,
                         label=f"λ_min (folded) = {lambda_min:.4g} pN/µm ({verdict})")
        axis.set_xlabel("Ritz index (ascending) — the interior is NOT the spectrum's interior")
        axis.set_ylabel("eigenvalue [pN/µm] (log axis)")
        axis.set_title(title, fontsize=10)
        axis.legend(fontsize=8, loc="lower right")

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def _explicit_step_figure(record: dict[str, Any], out: Path) -> Path:
    """Plot the stability product against the limit the explicit update is actually bounded by."""
    native = record["measurements"].get("native_spectrum")
    if native is None:
        raise SystemExit("[vis] this record has no native stage; nothing to plot here")
    verdict = native["explicit_step_verdict"]

    labels = ["production step\ndt_mu = 0.1/kmax\n(plain descent)",
              f"FIRE maximum\n{verdict['fire_max_multiplier']:g}× dt_mu\n(damped MD — other bound)"]
    products = [verdict["stability_product_base"], verdict["stability_product_at_fire_max"]]
    limit = verdict["forward_euler_limit"]
    # only the plain-descent bar is judged against this limit; FIRE's bar is drawn hatched because the
    # line does not apply to it, and colouring it as pass/fail would adjudicate what was not measured
    colours = [COLOR_PROJECTED if products[0] < limit else COLOR_LIMIT, COLOR_REFERENCE]

    fig, axis = plt.subplots(figsize=(7.4, 5.4))
    bars = axis.bar(labels, products, color=colours, width=0.55)
    bars[1].set_hatch("//")
    axis.axhline(limit, ls="-", lw=1.6, color=COLOR_LIMIT,
                 label=f"forward-Euler limit dt·λ_max = {limit:g} — applies to the PLAIN descent only")
    axis.set_yscale("log")
    axis.set_ylabel("dt · λ_max  (dimensionless, log axis)")
    axis.set_title(
        "is the production inner step stable at the measured λ_max?\n"
        f"λ_max = {verdict['measured_lambda_max_pN_per_um']:.4g} pN/µm  = "
        f"{verdict['lambda_max_over_kmax']:.3g}× the build's own kmax", fontsize=10)
    for index, value in enumerate(products):
        axis.text(index, value, f"  {value:.3g}", ha="center", va="bottom", fontsize=9)
    axis.legend(fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main() -> None:
    """Regenerate every figure for one native-spectrum record."""
    parser = argparse.ArgumentParser(description="Figures for the native cortex spectrum record.")
    parser.add_argument("--record", default="aleph/outputs/ac/observe/native_cortex_spectrum/record.json")
    parser.add_argument("--figs", default="", help="output directory (default <record dir>/figs)")
    args = parser.parse_args()

    path = Path(args.record)
    record = json.loads(path.read_text(encoding="utf-8"))
    figs = Path(args.figs) if args.figs else path.parent / "figs"
    figs.mkdir(parents=True, exist_ok=True)

    written = [_family_figure(record, figs / "family_certification.png")]
    if record["measurements"].get("native_spectrum") is not None:
        written.append(_ritz_figure(record, figs / "native_ritz.png"))
        written.append(_explicit_step_figure(record, figs / "explicit_step.png"))
    for item in written:
        print(f"[vis] wrote {item}")


if __name__ == "__main__":
    main()
