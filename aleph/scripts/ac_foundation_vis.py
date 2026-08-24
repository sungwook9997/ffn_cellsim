#!/usr/bin/env python3
"""Regenerate the AC foundation-hardening gate summary figure.

Reads only committed gate artifacts; it never launches simulation kernels.  All residuals are shown relative
to their predeclared structural limits, axes begin at zero unless explicitly labelled logarithmic, and the
native convergence panel overlays the fixed acceptance threshold.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from aleph.components.incumbent.compartments import MEMBRANE_DEFAULTS
from aleph.components.incumbent.preload_contract import evaluate_preload_capacity


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path,
        default=Path("aleph/outputs/ac/foundation-hardening"),
    )
    args = parser.parse_args()
    root = args.root
    gates = json.loads((root / "native_gates.json").read_text(encoding="utf-8"))
    native = json.loads((root / "native_resting_70686.json").read_text(encoding="utf-8"))
    arrays = np.load(root / "native_gates_arrays.npz")
    report = native["report"]

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.8), constrained_layout=True)
    fig.suptitle(
        "FF-AC foundation — conservation, retry, and NG-9 instrument pass; physiological t0 open",
        fontsize=15,
    )

    # Live mesh: report failure/mismatch channels explicitly, including closed-set surface ties.
    ax = axes[0, 0]
    live = gates["live_domain"]
    labels = ["rest\nquery fail", "rest\nvalid mismatch", "deformed\nquery fail", "deformed\nvalid mismatch"]
    values = [live["rest"]["query_failures"], live["rest"]["valid_cell_mismatches"],
              live["deformed"]["query_failures"], live["deformed"]["valid_cell_mismatches"]]
    ax.bar(np.arange(4), values, color="#3b82f6")
    ax.set_xticks(np.arange(4), labels)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("cells")
    ax.set_title("Live membrane−nucleus parity")
    for i, value in enumerate(values):
        ax.text(i, value, str(value), ha="center", va="bottom")
    rest_ties = live["rest"]["membrane_surface_ties"] + live["rest"]["nucleus_surface_ties"]
    deformed_ties = live["deformed"]["membrane_surface_ties"] + live["deformed"]["nucleus_surface_ties"]
    ax.text(
        0.02, 0.94, f"surface ties reported: rest {rest_ties}; deformed {deformed_ties}",
        transform=ax.transAxes, va="top",
    )

    # NG-2: normalize each result by its immutable tolerance so the acceptance band is common and linear.
    ax = axes[0, 1]
    ng2 = gates["ng2"]
    coupled = ng2["canonical_coupled_outer_step"]
    fsi = gates["fsi_conservation"]
    ng5 = gates["ng5_closed_system_force"]
    ng5_force = ng5["spatial_affine_pressure_bulk_plus_surface_net_force"]
    ng5_work = ng5["fsi_transfer_adjoint_work_identity"]
    normalized = np.array([
        ng2["impermeable"]["relative_content_drift"] / 1.0e-11,
        ng2["source_and_dilatation"]["relative_error"] / 1.0e-10,
        ng2["membrane_content_equals_flux"]["content_balance_relative_error"] / 1.0e-10,
        ng2["moving_domain"]["relative_error"] / 1.0e-10,
        ng2["device_terzaghi_dirichlet"]["degree_absolute_error"] / 3.0e-3,
        ng2["device_green_impulse"]["relative_error"] / 5.0e-2,
        coupled["physical_conservation_relative_error"] / 1.0e-10,
        fsi["spread_partition_conservation"]["momentum_relative_error"] / 1.0e-10,
        ng5_force["net_force_relative_error"] / 1.0e-10,
        ng5_work["relative_error"] / 1.0e-10,
    ])
    ax.bar(np.arange(normalized.size), normalized, color="#10b981")
    ax.axhline(1.0, color="#dc2626", linewidth=1.5, label="fixed limit")
    ax.set_xticks(
        np.arange(normalized.size),
        ["impermeable", "source/div", "membrane", "remap", "Terzaghi", "Green", "coupled", "FSI mom.",
         "NG-5 force", "NG-5 work"],
        rotation=28,
        ha="right",
        fontsize=9,
    )
    ax.set_ylabel("measured residual / fixed limit")
    ax.set_ylim(0, max(1.15, float(normalized.max()) * 1.2))
    ax.set_title("NG-2/NG-5 conservation + FSI transfer")
    ax.legend(frameon=False)
    coupled_status = coupled["status"]
    ax.text(
        0.02, 0.94,
        f"canonical ledger: {coupled_status}; FSI: {fsi['status']}; NG-5: {ng5['status']}",
        transform=ax.transAxes, va="top", color="#047857",
    )

    # Pressure device/reference identity; all points and the y=x reference share untruncated axes.
    ax = axes[0, 2]
    x = arrays["pressure_force_oracle_norm"]
    y = arrays["pressure_force_device_norm"]
    xl = arrays["pressure_linear_oracle_norm"]
    yl = arrays["pressure_linear_device_norm"]
    lim = float(max(x.max(), y.max(), xl.max(), yl.max())) * 1.03
    ax.scatter(x, y, s=10, alpha=0.4, color="#8b5cf6", edgecolors="none", label="uniform 40 Pa")
    ax.scatter(xl, yl, s=10, alpha=0.4, color="#0ea5e9", edgecolors="none", label="linear p(x)")
    ax.plot([0.0, lim], [0.0, lim], color="#111827", linewidth=1.2, label="device = oracle")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("oracle nodal force [pN]")
    ax.set_ylabel("CUDA nodal force [pN]")
    ax.set_title("NG-3 surface pressure trace")
    ax.legend(frameon=False)

    # Strict native convergence: log is stated; threshold remains visible. Values span >12 orders of magnitude.
    ax = axes[1, 0]
    tolerance = report["inner_tolerance_um"]
    ratios = np.array([
        report["inner_max_displacement_um"] / tolerance,
        report["inner_constraint_residual_um"] / tolerance,
    ])
    ax.bar([0, 1], ratios, color=["#dc2626", "#10b981"])
    ax.axhline(1.0, color="#111827", linewidth=1.5, label="convergence limit")
    ax.set_yscale("log")
    ax.set_xticks([0, 1], ["projected displacement", "segment constraint"])
    ax.set_ylabel("residual / tolerance (log scale)")
    ax.set_title("Full-native mechanics convergence")
    ax.legend(frameon=False)

    # Directional memcpy profiler with a non-vacuous positive control.
    ax = axes[1, 1]
    copies = [report["hot_loop_dtoh_copies"], report["profiler_positive_control_dtoh_copies"]]
    ax.bar([0, 1], copies, color=["#10b981", "#3b82f6"])
    ax.set_xticks([0, 1], ["physical-time hot loop", "post-loop .numpy control"])
    ax.set_ylabel("DtoH memcpy records")
    ax.set_ylim(0, max(1.2, max(copies) * 1.2))
    ax.set_title("NG-6 directional memcpy trace")
    for i, value in enumerate(copies):
        ax.text(i, value, str(value), ha="center", va="bottom")
    rollback = gates["outer_rejection_transaction"]
    retry = gates["device_retry_multistep"]
    memory = gates["ng9_native_measurement"]
    ax.text(
        0.02, 0.94,
        f"rejected-state rollback: {'PASS' if rollback['pass'] else 'FAIL'}\n"
        f"ERM candidate/commit: {'PASS' if rollback['erm_commit_only_rupture']['pass'] else 'FAIL'}\n"
        f"device retry / fail-stop: {retry['status']}\n"
        f"exact peak: {memory['status']}\n"
        f"NVML probe: {memory['memory_accounting_probe_status']}",
        transform=ax.transAxes, va="top",
    )

    # Physiological preload gap: the required continuum balance is a reference, not a fitted target.
    ax = axes[1, 2]
    required = gates["ng3"]["young_laplace_required_tension_pN_per_um"]
    plateau = float(MEMBRANE_DEFAULTS["gamma_mem"])
    ax.bar([0, 1], [required, plateau], color=["#dc2626", "#f59e0b"])
    ax.set_xticks([0, 1], ["40 Pa × 7.5 µm / 2", "membrane plateau"])
    ax.set_ylabel("surface tension [pN/µm]")
    ax.set_ylim(0, required * 1.18)
    ax.set_title("Why physiological t0 is not yet equilibrated")
    for i, value in enumerate([required, plateau]):
        ax.text(i, value, f"{value:g}", ha="center", va="bottom")

    out = root / "figs" / "foundation_ng2_ng3_ng6.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(out)

    preload_path = root / "preload_force_probe_proxy600_tension_only.json"
    if not preload_path.exists():
        preload_path = root / "preload_force_probe_proxy600.json"
    if preload_path.exists():
        proxy = json.loads(preload_path.read_text(encoding="utf-8"))
        long_path = root / "preload_force_probe_proxy600_tension_only_30k.json"
        proxy_long = json.loads(long_path.read_text(encoding="utf-8")) if long_path.exists() else None
        native_ledger = native.get("ledger", report.get("ledger", {}))
        proxy_ledger = proxy["ledger"]
        default_capacity = evaluate_preload_capacity(
            pressure_pa=float(native_ledger["resting_pressure_Pa"]),
            radius_um=7.5,
            membrane_area_um2=float(proxy_ledger["preload_membrane_area_um2"]),
            membrane_tension_pn_per_um=float(MEMBRANE_DEFAULTS["gamma_mem"]),
            erm_tether_count=int(native_ledger["membrane_erm_tethers"]),
            erm_rupture_force_pn=float(proxy_ledger["preload_erm_rupture_force_pn"]),
        )

        fig, axes_grid = plt.subplots(2, 2, figsize=(14.5, 9.2), constrained_layout=True)
        axes = axes_grid.ravel()
        fig.suptitle(
            "Physiological preload audit — capacity, unilateral ERM, and explicit-solver convergence",
            fontsize=14,
        )

        ax = axes[0]
        capacity_ratios = [
            default_capacity.pressure_capacity_ratio,
            float(proxy_ledger["preload_pressure_capacity_ratio"]),
        ]
        ax.scatter([0, 1], capacity_ratios, s=90, color=["#dc2626", "#3b82f6"], zorder=3)
        ax.vlines([0, 1], 0.1, capacity_ratios, colors=["#dc2626", "#3b82f6"], linewidth=5, alpha=0.55)
        ax.axhline(1.0, color="#111827", linewidth=1.5, label="40 Pa capacity requirement")
        ax.set_yscale("log")
        ax.set_ylim(0.1, max(capacity_ratios) * 1.8)
        ax.set_xticks([0, 1], ["642 ERM\nunsourced diagnostic", "422,094 ERM\n600/µm² proxy"])
        ax.set_ylabel("pressure-capacity ratio (log scale)")
        ax.set_title("Rupture-force upper bound")
        ax.legend(frameon=False, loc="upper left")
        for i, value in enumerate(capacity_ratios):
            ax.annotate(f"{value:.3g}×", (i, value), xytext=(0, 8), textcoords="offset points", ha="center")

        ax = axes[1]
        native_ratio = float(report["inner_max_displacement_um"]) / float(report["inner_tolerance_um"])
        proxy_ratio = (
            float(proxy["solver"]["max_displacement_um"]) / float(report["inner_tolerance_um"])
        )
        ratios = [native_ratio, proxy_ratio]
        labels = ["native default\n60 iters + WCA", "600/µm² proxy\n6k, no WCA"]
        colors = ["#dc2626", "#f59e0b"]
        if proxy_long is not None:
            ratios.append(float(proxy_long["solver"]["displacement_to_tolerance_ratio"]))
            labels.append("600/µm² proxy\n30k explicit")
            colors.append("#7c3aed")
        x_values = np.arange(len(ratios))
        ax.scatter(x_values, ratios, s=90, color=colors, zorder=3)
        ax.vlines(x_values, 1.0e-1, ratios, colors=colors, linewidth=5, alpha=0.55)
        ax.axhline(1.0, color="#111827", linewidth=1.5, label="fixed convergence limit")
        ax.set_yscale("log")
        ax.set_ylim(1.0e-1, max(ratios) * 2.0)
        ax.set_xticks(x_values, labels)
        ax.set_ylabel("maximum displacement / tolerance (log scale)")
        ax.set_title("Strict fixed-point residual")
        ax.legend(frameon=False)
        for i, value in enumerate(ratios):
            ax.annotate(f"{value:.3g}×", (i, value), xytext=(0, 8), textcoords="offset points", ha="center")

        ax = axes[2]
        history = proxy["solver"]["convergence_history"]
        ax.plot(
            history["iteration"],
            np.asarray(history["max_displacement_um"], dtype=float) / float(proxy["solver"]["tolerance_um"]),
            color="#f59e0b", linewidth=2, label="6k unilateral-ERM run",
        )
        if proxy_long is not None:
            history_long = proxy_long["solver"]["convergence_history"]
            iteration_long = np.asarray(history_long["iteration"], dtype=int)
            valid_long = iteration_long <= int(proxy_long["solver"]["iterations"])
            ax.plot(
                iteration_long[valid_long],
                np.asarray(history_long["max_displacement_um"], dtype=float)[valid_long]
                / float(proxy_long["solver"]["tolerance_um"]),
                color="#7c3aed", linewidth=1.5, alpha=0.85, label="30k explicit budget",
            )
        ax.axhline(1.0, color="#111827", linewidth=1.5, label="fixed convergence limit")
        ax.set_yscale("log")
        ax.set_xlabel("inner iteration")
        ax.set_ylabel("maximum displacement / tolerance (log scale)")
        ax.set_title("Late explicit mode stalls above the gate")
        ax.legend(frameon=False)

        ax = axes[3]
        quantile_labels = ["min", "median", "p90", "p99", "max"]
        tension_source = proxy_long if proxy_long is not None else proxy
        tension = np.asarray(tension_source["erm_candidate"]["tension_pN_quantiles"], dtype=float)
        rupture = float(proxy_ledger["preload_erm_rupture_force_pn"])
        ax.plot(np.arange(tension.size), tension, marker="o", linewidth=2, color="#3b82f6", label="ERM tension")
        ax.axhline(rupture, color="#dc2626", linewidth=1.5, label="derived rupture force")
        ax.set_xticks(np.arange(tension.size), quantile_labels)
        ax.set_ylim(0.0, rupture * 1.08)
        ax.set_ylabel("tether tension [pN]")
        ax.set_title("Proxy load remains sub-rupture")
        ax.legend(frameon=False)
        ax.text(
            0.02,
            0.92,
            f"0 ruptures; {tension_source['erm_candidate']['compressed_count']:,} compressed links carry 0 force\n"
            "proxy is zebrafish mesendoderm, not MCF7",
            transform=ax.transAxes,
            va="top",
        )

        preload_out = root / "figs" / "preload_capacity_and_relaxation.png"
        fig.savefig(preload_out, dpi=180)
        plt.close(fig)
        print(preload_out)


if __name__ == "__main__":
    main()
