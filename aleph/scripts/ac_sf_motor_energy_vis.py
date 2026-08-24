#!/usr/bin/env python
r"""Figures for the ``nmii_sf_motor`` lane's acceptance predicate and its energy ledger.

WHY A SEPARATE ENTRY POINT FROM ``ac_sf_arc_vis.py``.  That script renders the lane's GEOMETRY (where the
minifilaments sit, which sarcomeres carry a motor).  This one renders its ACCOUNTING — the balance gate's
margin and the energy ledger's closure — which share no data, no axes and no scene with the geometry.  The
project's rule is one entry point per unit so that nothing is left un-regenerated; both are listed in the
lane's ``REPORT.md`` §Figures, and both read only the committed run record.

The visual check earns its place here rather than being a formality.  Two of the things this plot makes
immediately visible are things a table hides: whether the balance margin is a comfortable one or is sitting
just under its tolerance, and whether the energy terms actually TRACK each other step to step rather than
happening to net out once.

Reads: ``outputs/ac/gate_b_sf_motor/sf_motor_native_gate.json`` (committed).
Writes: ``outputs/ac/gate_b_sf_motor/figs/*.png``.

Visualization integrity (project rule): no axis truncation, no silent log/linear flip (both log axes are
labelled as such), units annotated on every axis, and the reference the measurement is judged against —
the derived tolerance, and the two predicted scaling exponents — drawn ON the measurement rather than
described in a caption.

Runs on the dev Mac: matplotlib only, no CUDA, no Warp.
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


def _load(path: Path) -> dict[str, Any]:
    """Return the committed run record."""
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def figure_balance_margin(record: dict[str, Any], destination: Path) -> Path:
    """Plot the acceptance gate's residual against the tolerance it is judged by, per step.

    Both series are forces [pN] on one log axis — the quantity and its threshold, never on separate
    plots, because the only thing worth seeing is the MARGIN between them.  The positive control is drawn
    on the same axis, which is what shows the gate has a working range at all rather than a floor.
    """
    trace = record["measurements"]["trace"]
    steps = np.array([row["step"] for row in trace], np.float64)
    residual = np.array([row["balance_residual_pN"] for row in trace], np.float64)
    tolerance = np.array([row["balance_tolerance_pN"] for row in trace], np.float64)
    control = record["measurements"]["acceptance"]["positive_control"]

    figure, axis = plt.subplots(figsize=(8.0, 5.0))
    axis.semilogy(steps, np.maximum(residual, 1e-32), "o-", color="#1f77b4",
                  label=r"measured $|\mathrm{reaction}+\mathrm{traction}|$")
    axis.semilogy(steps, tolerance, "s--", color="#d62728",
                  label=r"derived tolerance $\gamma_n \cdot \Sigma|f_i|$  (PI D8)")
    if control.get("ran"):
        axis.axhline(control["balance_residual_pN"], color="#2ca02c", linestyle=":", linewidth=2)
        axis.text(steps[0], control["balance_residual_pN"] * 1.6,
                  f"positive control: {control['injected_unbalanced_force_pN']:g} pN injected one-sided "
                  f"→ REJECTED", color="#2ca02c", fontsize=9)
    axis.set_xlabel("accepted physical step")
    axis.set_ylabel("force [pN]  (log axis)")
    axis.set_title("Acceptance is decided by the device force-balance gate\n"
                   "(both channels are the two never-merged force arrays; the sum must vanish)")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend(loc="center right", fontsize=9)
    figure.tight_layout()
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def figure_energy_ledger(record: dict[str, Any], destination: Path) -> Path:
    """Plot the four energy terms per measured step, and the residual's mobility scaling beside them.

    Left panel: the terms on a linear axis, so a reader can see that ``ΔU + dissipated`` really does track
    ``W_active`` step to step instead of netting out once.  Right panel: the residual at ``μ`` and ``μ/2``
    against BOTH predicted power laws, which is what identifies the residual as the integrator's own
    first-order consistency error rather than a leak — the same self-identifying argument the loop-work
    figure uses one level down.
    """
    energy = record["measurements"]["energy_ledger"]
    rows = energy["per_step"]
    steps = np.array([row["step"] for row in rows], np.float64)
    delta_u = np.array([row["delta_potential_pN_um"] for row in rows], np.float64)
    dissipated = np.array([row["dissipated_pN_um"] for row in rows], np.float64)
    active = np.array([row["active_input_pN_um"] for row in rows], np.float64)
    residual = np.array([row["residual_pN_um"] for row in rows], np.float64)

    figure, (left, right) = plt.subplots(1, 2, figsize=(13.0, 5.0))

    width = 0.26
    left.bar(steps - width, delta_u, width, label=r"$\Delta U$  ($-\int F_{passive}\cdot dx$)",
             color="#1f77b4")
    left.bar(steps, dissipated, width, label=r"dissipated  ($\Sigma|\Delta x|^2/\mu$, on device)",
             color="#ff7f0e")
    left.bar(steps + width, -active, width, label=r"$-W_{active}$  (bound $-$ detached)",
             color="#2ca02c")
    left.plot(steps, residual, "kD-", markersize=6, label="residual (sum of the three)")
    left.axhline(0.0, color="k", linewidth=0.8)
    left.set_xlabel("accepted physical step")
    left.set_ylabel(r"energy [pN$\cdot\mu$m]  (linear axis)")
    left.set_title("The balance, term by term — no energy function is implemented\n"
                   r"$\Delta U + \mathrm{dissipated} + \mathrm{event\ jump} = W_{active}$")
    left.set_xticks(steps)
    left.grid(True, axis="y", alpha=0.3)
    # Headroom for the legend rather than a smaller axis range: the bars must stay untruncated (project
    # visualization rule), so the room is made above them, never by clipping them.
    ceiling = float(np.max(np.abs(np.concatenate((delta_u, dissipated, active)))))
    left.set_ylim(-1.15 * ceiling, 1.75 * ceiling)
    left.legend(fontsize=8, loc="upper left")

    scaling = energy.get("mobility_scaling", {})
    if scaling.get("ran"):
        base = scaling["base"]
        fine = scaling["halved_mobility"]
        convergence = scaling["convergence"]
        mobility = np.array([base["mobility_um_per_pN"], fine["mobility_um_per_pN"]], np.float64)
        measured = np.abs(np.array([base["residual_pN_um"], fine["residual_pN_um"]], np.float64))
        right.loglog(mobility, measured, "ko-", markersize=8,
                     label=f"measured  (exponent {convergence['residual_exponent']:.4f})")
        span = np.array([mobility.min() * 0.7, mobility.max() * 1.4])
        anchor = measured[0] / mobility[0] ** 1.0
        right.loglog(span, anchor * span ** 1.0, "--", color="#1f77b4",
                     label=r"integrator consistency:  residual $\propto \mu^1$")
        right.loglog(span, np.full_like(span, measured[0]), "--", color="#d62728",
                     label=r"genuine leak:  residual independent of $\mu$")
        right.set_xlabel(r"mobility $\mu = dt/\gamma$ [$\mu$m/pN]  (log axis)")
        right.set_ylabel(r"$|$residual$|$ [pN$\cdot\mu$m]  (log axis)")
        right.set_title("What the residual IS, from its own scaling\n"
                        f"verdict: {convergence['verdict']}")
        right.grid(True, which="both", alpha=0.3)
        right.legend(fontsize=9, loc="lower right")
    else:
        right.text(0.5, 0.5, "mobility scaling not run", ha="center", va="center")
        right.set_axis_off()

    figure.tight_layout()
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Regenerate the nmii_sf_motor acceptance + energy-ledger figures.")
    parser.add_argument(
        "--record", type=Path,
        default=Path("aleph/outputs/ac/gate_b_sf_motor/sf_motor_native_gate.json"),
        help="the committed run record to plot")
    parser.add_argument(
        "--figs", type=Path, default=Path("aleph/outputs/ac/gate_b_sf_motor/figs"),
        help="destination directory")
    args = parser.parse_args()

    record = _load(args.record)
    for path in (
        figure_balance_margin(record, args.figs / "balance_gate_margin.png"),
        figure_energy_ledger(record, args.figs / "energy_ledger.png"),
    ):
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
