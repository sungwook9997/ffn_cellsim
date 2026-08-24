#!/usr/bin/env python
r"""Figures for the ``ac/engine/observe`` run record — one entry point, regenerates everything.

Reads the observer's committed JSON + NPZ (no GPU, no re-run) and writes four PNGs.  Each figure
carries the structural PREDICTION it is judged against drawn on the same axes, because the project's
visualization rule is that a measurement without its reference overlaid is not a check:

  1. ``spectrum.png``       — the full eigenvalue field at both configurations, with the eigensolve's
                              own resolvability floor, the true ``λ_max`` and the Gershgorin bound.
  2. ``localization.png``   — participation ratio and spatial extent per mode: the stress decay length.
  3. ``sensitivity.png``    — the Fisher eigenvalue spectrum (sloppiness) and the per-parameter
                              response of the observable.
  4. ``loop_work.png``      — the loop residual against BOTH predicted power laws, for the passive
                              field, the active field, and the heads-detached control.

Usage (dev Mac, no CUDA needed)::

    python aleph/scripts/ac_observe_vis.py \
        --record aleph/outputs/ac/observe/sf_operator_observables.json

Visualization integrity (project rule): no axis truncation; every log axis is labelled as such in the
axis title; units are SI-derived engine units (pN, µm, pN/µm) and named on every axis; the reference
prediction is overlaid on the measurement rather than described in a caption.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

#: One palette across all four figures so a colour means the same configuration everywhere.
COLOR_T0 = "#3b6ea5"
COLOR_BOUND = "#c1502e"
COLOR_CONTROL = "#4f8a55"
COLOR_REFERENCE = "#6b6b6b"


def _spectrum_figure(record: dict, arrays: dict, out: Path) -> Path:
    """Plot the full eigenvalue field at both configurations against the floors that bound it."""
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))
    for axis, key, colour, title in (
        (axes[0], "operator_t0", COLOR_T0, "t0 — no head bound (passive)"),
        (axes[1], "operator_bound", COLOR_BOUND, "after 10 accepted steps (motor engaged)"),
    ):
        spectrum = record["measurements"][key]["spectrum"]
        values = np.sort(arrays["eigenvalues_t0" if key == "operator_t0" else "eigenvalues_bound"])
        index = np.arange(values.size)
        positive = values > 0.0
        axis.semilogy(index[positive], values[positive], ".", ms=2.5, color=colour,
                      label=f"eigenvalues (n_dof = {spectrum['n_dof']})")
        axis.axhline(spectrum["zero_tolerance_pN_per_um"], ls="--", lw=1.2, color=COLOR_REFERENCE,
                     label=f"resolvability floor n·eps·‖K‖ = {spectrum['zero_tolerance_pN_per_um']:.2e}")
        axis.axhline(spectrum["lambda_max_pN_per_um"], ls="-", lw=1.0, color="black",
                     label=f"true λ_max = {spectrum['lambda_max_pN_per_um']:.4g} pN/µm")
        axis.axhline(spectrum["gershgorin_bound_pN_per_um"], ls=":", lw=1.4, color="black",
                     label=(f"Gershgorin bound = {spectrum['gershgorin_bound_pN_per_um']:.4g} "
                            f"({spectrum['gershgorin_over_lambda_max']:.3g}× λ_max)"))
        axis.axvline(spectrum["n_zero_modes"], ls="-.", lw=1.0, color=colour,
                     label=(f"{spectrum['n_zero_modes']} modes below the floor "
                            f"(rigid-body bound {spectrum['rigid_body_zero_modes']}, "
                            f"excess {spectrum['n_excess_zero_modes']})"))
        axis.set_title(title)
        axis.set_xlabel("mode index (ascending eigenvalue)")
        axis.set_ylabel("eigenvalue of K = ∇²U  [pN/µm]  (log axis)")
        axis.legend(fontsize=7.0, loc="lower right")
        axis.grid(alpha=0.25)
    fig.suptitle(
        "nmii_sf_motor slice — the stiffness spectrum, measured (positive eigenvalues shown on a log "
        "axis; the count below the floor is marked, none is discarded)", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def _localization_figure(record: dict, out: Path) -> Path:
    """Plot how far each softest mode reaches — the stress decay length, measured per mode."""
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))
    for key, colour, label in (
        ("operator_t0", COLOR_T0, "t0 (passive)"),
        ("operator_bound", COLOR_BOUND, "motor engaged"),
    ):
        modes = record["measurements"][key]["softest_modes"]
        if not modes:
            continue
        values = [m["eigenvalue_pN_per_um"] for m in modes]
        axes[0].semilogx(values, [m["participation_ratio_nodes"] for m in modes], "o-", ms=4,
                         color=colour, label=label)
        axes[1].semilogx(values, [m["gyration_radius_um"] for m in modes], "o-", ms=4,
                         color=colour, label=label)
    n_nodes = record["census"]["n_nodes_concatenated"]
    axes[0].axhline(n_nodes, ls="--", lw=1.2, color=COLOR_REFERENCE,
                    label=f"fully delocalised = all {n_nodes} nodes")
    axes[0].axhline(1.0, ls=":", lw=1.2, color=COLOR_REFERENCE, label="single-node mode")
    axes[0].set_ylabel("participation ratio  [nodes occupied]")
    axes[1].set_ylabel("energy-weighted gyration radius  [µm]")
    for axis in axes:
        axis.set_xlabel("mode eigenvalue  [pN/µm]  (log axis)")
        axis.legend(fontsize=8.0)
        axis.grid(alpha=0.25)
    fig.suptitle("softest modes: how far a mode reaches through the slice = the stress decay length",
                 fontsize=10.5)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def _sensitivity_figure(record: dict, arrays: dict, out: Path) -> Path:
    """Plot the sloppiness spectrum and the per-parameter response of the observable."""
    fisher = record["measurements"]["sensitivity"]["fisher"]
    response = record["measurements"]["sensitivity"]["sensitivity"][
        "log_response_norm_by_parameter"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))

    values = np.asarray(fisher["eigenvalues"], np.float64)
    positive = values > 0.0
    axes[0].semilogy(np.arange(values.size)[positive] + 1, values[positive], "o-", ms=5,
                     color=COLOR_BOUND, label="Fisher eigenvalues (descending)")
    axes[0].axhline(fisher["resolvable_floor"], ls="--", lw=1.2, color=COLOR_REFERENCE,
                    label=f"resolvability floor p·eps·λ_max = {fisher['resolvable_floor']:.2e}")
    axes[0].set_xlabel("direction index")
    axes[0].set_ylabel("Fisher eigenvalue  [dimensionless]  (log axis)")
    axes[0].set_title(
        f"sloppiness = {fisher['sloppiness_decades']:.3g} decades; effective dimension "
        f"{fisher['effective_dimension']:.3g} of {len(fisher['parameters'])}")
    axes[0].legend(fontsize=8.0)
    axes[0].grid(alpha=0.25)

    names = list(response)
    axes[1].barh(range(len(names)), [response[name] for name in names], color=COLOR_T0)
    axes[1].set_yticks(range(len(names)))
    axes[1].set_yticklabels(names, fontsize=8.0)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("‖∂ log λ / ∂ log θ‖₂ over the observed spectrum  [dimensionless]")
    axes[1].set_title("how much the spectrum responds to each stiffness family")
    axes[1].grid(alpha=0.25, axis="x")

    fig.suptitle(
        "Fisher information of the spectrum — which stiffness parameters this observable could infer",
        fontsize=10.5)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def _loop_work_figure(record: dict, out: Path) -> Path:
    """Plot the loop residual against BOTH predicted power laws for all three configurations."""
    measurements = record["measurements"]
    series = (
        ("loop_work_t0", COLOR_T0, "t0 — passive"),
        ("loop_work_bound", COLOR_BOUND, "motor engaged"),
        ("loop_work_bound_heads_detached_control", COLOR_CONTROL,
         "control: same configuration, heads detached"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    # The passive series and the heads-detached control COINCIDE — that coincidence is the result, so
    # the two are drawn with different marker sizes and line widths rather than one hiding the other.
    styles = {
        "loop_work_t0": dict(marker="o", ms=11, mfc="none", mew=1.8, lw=3.0, alpha=0.9),
        "loop_work_bound": dict(marker="s", ms=6, lw=1.8),
        "loop_work_bound_heads_detached_control": dict(marker="^", ms=5, lw=1.2, ls="--"),
    }
    for key, colour, label in series:
        block = measurements.get(key)
        if not block:
            continue
        base, half = block["base"], block["half_amplitude"]
        refined = block["refined_segments"]
        verdict = block["convergence"]["verdict"]
        style = dict(styles[key], color=colour, label=f"{label} → {verdict}")
        axes[0].loglog([half["amplitude_um"], base["amplitude_um"]],
                       [abs(half["loop_work_pN_um"]), abs(base["loop_work_pN_um"])], **style)
        axes[1].loglog([base["n_segments_per_edge"], refined["n_segments_per_edge"]],
                       [abs(base["loop_work_pN_um"]), abs(refined["loop_work_pN_um"])], **style)

    # Reference slopes, anchored on the passive series so they are comparable by eye.
    anchor = measurements["loop_work_t0"]["base"]
    a0, w0 = anchor["amplitude_um"], abs(anchor["loop_work_pN_um"])
    amplitudes = np.array([0.5 * a0, a0])
    axes[0].loglog(amplitudes, w0 * (amplitudes / a0) ** 3.0, "--", lw=1.2, color=COLOR_REFERENCE,
                   label="predicted quadrature slope  a³")
    axes[0].loglog(amplitudes, w0 * (amplitudes / a0) ** 2.0, ":", lw=1.4, color=COLOR_REFERENCE,
                   label="predicted circulation slope  a²")
    n0 = anchor["n_segments_per_edge"]
    segments = np.array([n0, 2 * n0], np.float64)
    axes[1].loglog(segments, w0 * (segments / n0) ** -2.0, "--", lw=1.2, color=COLOR_REFERENCE,
                   label="predicted quadrature slope  N⁻²")
    axes[1].loglog(segments, np.full_like(segments, w0), ":", lw=1.4, color=COLOR_REFERENCE,
                   label="predicted circulation: N-independent")

    axes[0].set_xlabel("circuit edge amplitude  [µm]  (log axis)")
    axes[1].set_xlabel("midpoint segments per edge  (log axis)")
    for axis in axes:
        axis.set_ylabel("|∮F·dx|  [pN·µm]  (log axis)")
        axis.legend(fontsize=7.0)
        axis.grid(alpha=0.25, which="both")
    fig.suptitle(
        "does the assembled force field admit an energy?  The residual identifies itself by its own "
        "scaling — no tolerance is applied", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main() -> None:
    """Regenerate every observer figure from the committed record."""
    ap = argparse.ArgumentParser(description="Figures for the ac/engine/observe run record.")
    ap.add_argument("--record", type=str,
                    default="aleph/outputs/ac/observe/sf_operator_observables.json")
    ap.add_argument("--figs", type=str, default="")
    args = ap.parse_args()

    record_path = Path(args.record)
    record = json.loads(record_path.read_text())
    arrays = np.load(record_path.with_suffix(".npz"))
    figs = Path(args.figs) if args.figs else record_path.parent / "figs"
    figs.mkdir(parents=True, exist_ok=True)

    written = [
        _spectrum_figure(record, arrays, figs / "spectrum.png"),
        _localization_figure(record, figs / "localization.png"),
        _sensitivity_figure(record, arrays, figs / "sensitivity.png"),
        _loop_work_figure(record, figs / "loop_work.png"),
    ]
    for path in written:
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
