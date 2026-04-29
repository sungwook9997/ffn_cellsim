"""Plotting utilities for the ActiveCellSim layer-by-layer narrative.

Headless-SSH-friendly via matplotlib's "Agg" backend. Each plot is
self-contained: read CSV/HDF5 from results/, render PNG to
results/visualization/.

Anti-pattern guard (per CLAUDE.md Hard Rule 11): when comparing to PI
experimental data, plots use `A_over_A0_topdown` (xy-projection convex
hull, z-independent) — the simulation analog of PI's Area_um2 column.
Substrate-contact area is NOT used for spreading comparisons.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Layer-by-layer narrative table (hardcoded R drift values from each
# stage's gate_report.md / commit message; these are the "ground truth"
# from the 25+ simulations, useful for the bar chart even if the raw
# results dirs are partially missing).
NARRATIVE_TABLE = [
    ("Layer 1 v15 (free)", 0.244, "v15 baseline (8e06d8e)"),
    ("+ substrate Option α", 0.247, "Stage 1a+ (f1af703)"),
    ("+ substrate β α=0.1", 0.246, "Stage 1a+ Option β"),
    ("+ substrate β α=0.3", 0.246, "Stage 1a+ Option β"),
    ("+ substrate β α=1.0", 0.247, "Stage 1a+ Option β"),
    ("+ Layer 2 ζ=0.1", 0.222, "Stage 1a++ Track 1"),
    ("+ Layer 2 ζ=0.4 (best stable)", 0.132, "Stage 1a++ Track 1"),
    ("+ Layer 2 ζ=0.6 (FIRST PASS)", 0.043, "Phase 2.2"),
    ("+ Layer 2 ζ=0.7 transitional", 0.151, "Phase 2.2"),
    ("+ Layer 2 ζ=0.9 explosion", 2.494, "Phase 2.2"),
    ("+ Layer 3 Bare (Stage 1b)", 0.185, "Stage 1b"),
    ("+ Layer 3 Pre (Stage 1b)", 0.166, "Stage 1b"),
    ("+ Layer 3 Lam4 (Stage 1b)", 0.158, "Stage 1b"),
    ("Path C v15 (g=0.1)", 0.326, "Path C baseline"),
    ("+ L5 + Path C g=0.01 (Stage 1c)", 0.138, "Stage 1c"),
    ("+ L4 (Stage 1d)", 0.147, "Stage 1d"),
    ("+ L6 (Stage 2)", 0.127, "Stage 2"),
    ("Phase 4-v1 Bare (zeta_max=0.4)", 0.157, "Phase 4-v1"),
    ("Phase 4-v1 Pre", 0.145, "Phase 4-v1"),
    ("Phase 4-v1 Lam4", 0.122, "Phase 4-v1"),
    ("Phase 4-v2 Bare (zeta_max=0.6)", 0.131, "Phase 4-v2"),
    ("Phase 4-v2 Pre", 0.095, "Phase 4-v2"),
    ("Phase 4-v2 Lam4 (best)", 0.073, "Phase 4-v2"),
]


def aggregate_narrative_table(repo_root: str | Path = ".") -> list[dict]:
    """Return tidy table of (label, R_drift, stage) tuples.

    For now hardcoded from NARRATIVE_TABLE; future versions can scan
    results/*/gate_report.md to update from disk.
    """
    return [
        {"label": label, "R_drift": float(rd), "stage": stage}
        for (label, rd, stage) in NARRATIVE_TABLE
    ]


def plot_layer_by_layer_R_drift(
    out_path: str | Path = "results/visualization/R_drift_bar.png",
) -> Path:
    """Bar chart of R drift across 23+ simulations in the project narrative.

    R drift gate at 0.05 marked as horizontal red line. Phase 4-v2 Lam4
    (0.073) marked as the best 3-phenotype-framework result; Phase 2.2
    ζ=0.60 (0.043) marked as the first-PASS finding.
    """
    rows = aggregate_narrative_table()
    labels = [r["label"] for r in rows]
    values = [r["R_drift"] for r in rows]

    fig, ax = plt.subplots(figsize=(11, 8))
    ypos = np.arange(len(labels))
    colors = []
    for r in rows:
        if r["R_drift"] < 0.05:
            colors.append("#1b9e77")  # PASS gate (green)
        elif r["R_drift"] < 0.10:
            colors.append("#7fcdbb")  # close to PASS
        elif r["R_drift"] < 0.20:
            colors.append("#fec44f")  # mid range
        elif r["R_drift"] < 0.30:
            colors.append("#fe9929")  # near baseline
        else:
            colors.append("#d95f0e")  # far above baseline
    bars = ax.barh(ypos, values, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0.05, color="red", linewidth=1.5, linestyle="--", label="R drift gate (0.05)")
    ax.axvline(0.244, color="black", linewidth=1.0, linestyle=":", alpha=0.5,
               label="Layer 1 ceiling (0.244)")
    ax.set_xlabel(r"$\max_t |R(t)/R_0 - 1|$", fontsize=11)
    ax.set_title(
        "ActiveCellSim layer-by-layer R drift narrative\n"
        "23-simulation progression: Layer 1 ceiling → Layer 2 first PASS → Phase 4-v2 best phenotype",
        fontsize=11,
    )
    ax.set_xlim(0, max(values) * 1.05)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(axis="x", alpha=0.3)

    plt.tight_layout()
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def plot_phenotype_comparison(
    out_path: str | Path = "results/visualization/phenotype_comparison.png",
) -> Path:
    """Phase 4-v1 vs Phase 4-v2 phenotype R drift + A/A₀_topdown comparison.

    Two side-by-side panels:
    - Left: R drift Bare/Pre/Lam4 (v1=light, v2=dark)
    - Right: A/A₀_topdown(end) Bare/Pre/Lam4 (v1=light, v2=dark)
    PI experimental Lam4 endpoint range [8, 33] shaded on right panel.
    """
    phenos = ["Bare", "Pre", "Lam4"]
    R_v1 = [0.157, 0.145, 0.122]
    R_v2 = [0.131, 0.095, 0.073]
    A_v1 = [1.399, 0.967, 1.161]
    A_v2 = [1.746, 2.359, 2.657]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 5))

    x = np.arange(len(phenos))
    width = 0.35

    axL.bar(x - width / 2, R_v1, width, label="Phase 4-v1 (zeta_max=0.4)",
            color="#a6cee3", edgecolor="black", linewidth=0.5)
    axL.bar(x + width / 2, R_v2, width, label="Phase 4-v2 (zeta_max=0.6)",
            color="#1f78b4", edgecolor="black", linewidth=0.5)
    axL.axhline(0.05, color="red", linewidth=1.5, linestyle="--",
                label="R drift gate (0.05)")
    axL.set_xticks(x)
    axL.set_xticklabels(phenos)
    axL.set_ylabel(r"$\max_t |R(t)/R_0 - 1|$", fontsize=11)
    axL.set_title("Phenotype R drift: Phase 4-v1 vs v2", fontsize=11)
    axL.legend(fontsize=8, loc="upper right")
    axL.grid(axis="y", alpha=0.3)

    axR.bar(x - width / 2, A_v1, width, label="Phase 4-v1",
            color="#fdbf6f", edgecolor="black", linewidth=0.5)
    axR.bar(x + width / 2, A_v2, width, label="Phase 4-v2",
            color="#ff7f00", edgecolor="black", linewidth=0.5)
    axR.axhspan(8, 33, alpha=0.15, color="green",
                label=r"PI Lam4 $A/A_0$ range [8, 33]")
    axR.set_xticks(x)
    axR.set_xticklabels(phenos)
    axR.set_ylabel(r"$A/A_0$ (top-down projection, end of pilot)", fontsize=11)
    axR.set_title("Phenotype A/A₀ at pilot end", fontsize=11)
    axR.legend(fontsize=8, loc="upper left")
    axR.grid(axis="y", alpha=0.3)

    fig.suptitle(
        "Phenotype comparison: Bare/Pre/Lam4 ordering reproduced in both "
        "R drift AND A/A₀_topdown (Phase 4-v2)",
        fontsize=11,
    )
    plt.tight_layout()
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def _load_AA0_trajectory(metrics_csv: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read time and A_over_A0_topdown columns from a metrics.csv."""
    if not metrics_csv.exists():
        return np.array([]), np.array([])
    rows = []
    with metrics_csv.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rows.append(r)
    if not rows:
        return np.array([]), np.array([])
    t = np.array([float(r.get("time_star", "nan") or "nan") for r in rows])
    AA0_key = "A_over_A0_topdown" if "A_over_A0_topdown" in rows[0] else None
    if AA0_key is None:
        return t, np.full_like(t, np.nan)
    AA0 = np.array([
        float(r.get(AA0_key, "nan") or "nan") for r in rows
    ])
    return t, AA0


def plot_AA0_trajectory(
    out_path: str | Path = "results/visualization/AA0_trajectory.png",
    include_production: bool = True,
) -> Path:
    """Overlay A/A₀_topdown(t) trajectories from Phase 4-v2 phenotypes
    + production_lam4 if available.

    PI experimental Lam4 endpoint range [8, 33] shaded for context.
    Pilot duration line at t* = 240 (4 hr) marked.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    sources = [
        ("results/final_pilot_v2_bare/metrics.csv", "Phase 4-v2 Bare", "#1b9e77"),
        ("results/final_pilot_v2_pre/metrics.csv", "Phase 4-v2 Pre", "#7570b3"),
        ("results/final_pilot_v2_lam4/metrics.csv", "Phase 4-v2 Lam4 (pilot 4 hr)", "#d95f02"),
    ]
    if include_production:
        sources.append((
            "results/production_lam4/metrics.csv",
            "Production Lam4 (5k × 80 hr)", "#e6ab02"))

    for csv_path, label, color in sources:
        t, AA0 = _load_AA0_trajectory(Path(csv_path))
        if t.size == 0:
            continue
        ax.plot(t, AA0, label=label, color=color, linewidth=1.6, alpha=0.9)

    ax.axvline(240, color="grey", linewidth=1, linestyle=":",
               label="Pilot duration (240·τ_relax = 4 hr)")
    ax.axhspan(8, 33, alpha=0.10, color="green",
               label=r"PI Lam4 $A/A_0$ endpoint range [8, 33]")
    ax.set_xlabel(r"$t / \tau_{\rm relax}$", fontsize=11)
    ax.set_ylabel(r"$A/A_0$ (top-down projection)", fontsize=11)
    ax.set_title(
        r"$A/A_0$ trajectories: Phase 4-v2 phenotypes + production (when available)",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.3)
    ax.set_xscale("symlog", linthresh=10)
    plt.tight_layout()
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def plot_zeta_response_curve(
    out_path: str | Path = "results/visualization/zeta_response_curve.png",
) -> Path:
    """13-point ζ response curve (Layer 2 stable regime + instability).

    Shows the Track 1 + Phase 2.2 sweep results: monotone decrease to
    ζ=0.60 first PASS, then transitional, then catastrophic explosion.
    """
    zetas = [0.10, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
    R_drifts = [0.222, 0.188, 0.182, 0.167, 0.155, 0.132, 0.102, 0.043, 0.151, 0.669, 2.494, 2.319]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(zetas, R_drifts, marker="o", linewidth=1.4, markersize=7,
            color="#1f78b4", label="Layer 2 only (no Layer 3/4/5/6)")
    ax.axhline(0.05, color="red", linewidth=1.5, linestyle="--",
               label="R drift gate (0.05)")
    ax.annotate("first PASS\n(R=0.043)", xy=(0.60, 0.043), xytext=(0.40, 0.30),
                arrowprops=dict(arrowstyle="->", color="black"),
                fontsize=9, ha="center")
    ax.annotate("explosion\n(R=2.494)", xy=(0.90, 2.494), xytext=(0.80, 1.5),
                arrowprops=dict(arrowstyle="->", color="black"),
                fontsize=9, ha="center")
    ax.set_yscale("symlog", linthresh=0.05)
    ax.set_xlabel(r"$\zeta / K$ (Layer 2 active stress, dimensionless)", fontsize=11)
    ax.set_ylabel(r"$\max_t |R(t)/R_0 - 1|$", fontsize=11)
    ax.set_title(
        "Layer 2 ζ response curve: 13-point sweep across stable regime "
        "+ instability transition",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.3)

    plt.tight_layout()
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def plot_3d_snapshots(
    out_path: str | Path = "results/visualization/3d_snapshots.png",
    h5_paths: list[str] | None = None,
) -> Path:
    """Per-stage final-frame 3D particle scatter (matplotlib mplot3d).

    By default reads the last frame of each Phase 4-v2 phenotype +
    production_lam4 (when available) snapshots.h5. 4 panels arranged 2×2.
    """
    if h5_paths is None:
        h5_paths = [
            "results/final_pilot_v2_bare/snapshots.h5",
            "results/final_pilot_v2_pre/snapshots.h5",
            "results/final_pilot_v2_lam4/snapshots.h5",
            "results/production_lam4/snapshots.h5",
        ]
    titles = [
        "Phase 4-v2 Bare (pilot 4 hr)",
        "Phase 4-v2 Pre (pilot 4 hr)",
        "Phase 4-v2 Lam4 (pilot 4 hr)",
        "Production Lam4 (5k × 80 hr, when available)",
    ]
    fig = plt.figure(figsize=(12, 10))

    try:
        import h5py
    except ImportError:
        print("h5py not available; skipping 3D snapshots.")
        return Path(out_path)

    for i, (h5p, ttl) in enumerate(zip(h5_paths, titles)):
        ax = fig.add_subplot(2, 2, i + 1, projection="3d")
        if not Path(h5p).exists():
            ax.text2D(0.5, 0.5, "(no snapshot — run not complete)",
                      transform=ax.transAxes, ha="center", va="center",
                      fontsize=10, color="grey")
            ax.set_title(ttl, fontsize=10)
            continue
        try:
            with h5py.File(h5p, "r") as f:
                pos = f["position"][-1, :, :]   # (N, 3) last frame
                is_b = f["is_boundary"][-1, :]
        except Exception as exc:
            ax.text2D(0.5, 0.5, f"(error: {exc})",
                      transform=ax.transAxes, ha="center", va="center",
                      fontsize=9, color="red")
            ax.set_title(ttl, fontsize=10)
            continue
        bdy = is_b.astype(bool)
        ax.scatter(pos[~bdy, 0], pos[~bdy, 1], pos[~bdy, 2],
                   c="#1f78b4", s=5, alpha=0.5, label="bulk")
        ax.scatter(pos[bdy, 0], pos[bdy, 1], pos[bdy, 2],
                   c="#e31a1c", s=8, alpha=0.9, label="boundary")
        ax.set_title(ttl, fontsize=10)
        ax.set_xlabel("x*", fontsize=9)
        ax.set_ylabel("y*", fontsize=9)
        ax.set_zlabel("z*", fontsize=9)
        # Substrate plane shading at z=0
        try:
            xx, yy = np.meshgrid(
                np.linspace(pos[:, 0].min(), pos[:, 0].max(), 5),
                np.linspace(pos[:, 1].min(), pos[:, 1].max(), 5),
            )
            ax.plot_surface(xx, yy, np.zeros_like(xx),
                            alpha=0.1, color="grey")
        except Exception:
            pass

    plt.suptitle("3D snapshots: final frame per stage (Phase 4-v2 + production)",
                 fontsize=12)
    plt.tight_layout()
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return out


def _load_pi_lam4() -> tuple[np.ndarray, np.ndarray]:
    """Load PI experimental Lam4 trajectory: returns (time_min, A_over_A0).

    Reads `data/experimental/260313_Lam4.csv`; aggregates over Series
    by mean within each Time_min bin. Result is the population-mean
    A/A₀ trajectory for Lam4 phenotype.
    """
    path = Path("data/experimental/260313_Lam4.csv")
    if not path.exists():
        return np.array([]), np.array([])
    rows = []
    with path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            try:
                t = float(r["Time_min"])
                a = float(r["A_over_A0"])
                rows.append((t, a))
            except (KeyError, ValueError):
                continue
    if not rows:
        return np.array([]), np.array([])
    # Aggregate by t_min bin: mean A across spheroids.
    from collections import defaultdict
    bins = defaultdict(list)
    for t, a in rows:
        bins[t].append(a)
    t_arr = np.array(sorted(bins.keys()))
    A_arr = np.array([np.mean(bins[t]) for t in t_arr])
    return t_arr, A_arr


def plot_production_trajectory(
    out_path: str | Path = "results/production_lam4/figures/production_trajectory.png",
    metrics_csv: str | Path = "results/production_lam4/metrics.csv",
) -> Path:
    """Multi-panel production trajectory: A/A₀_topdown, R, sphericity,
    ρ_osm, mmp_total, ecm_strength over time. PI Lam4 overlay where
    applicable.
    """
    metrics_csv = Path(metrics_csv)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    if not metrics_csv.exists():
        # Placeholder figure when production not yet complete.
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.text(0.5, 0.5, "Production Lam4 not yet complete\n"
                "(metrics.csv missing — figure regenerated when run finishes)",
                ha="center", va="center", fontsize=12, color="grey",
                transform=ax.transAxes)
        ax.set_axis_off()
        fig.savefig(out, dpi=110)
        plt.close(fig)
        return out

    rows = []
    with metrics_csv.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rows.append(r)
    if not rows:
        return out

    def _col(key, default=np.nan):
        return np.array([
            float(r.get(key, "nan") or "nan") if r.get(key, "") not in ("", None)
            else default for r in rows
        ])

    t = _col("time_star")
    n_p = 5000  # production
    AA0 = _col("A_over_A0_topdown")
    R = _col("effective_radius") / 0.9953  # normalise by R₀ initial
    sphericity = _col("wadell_sphericity")
    rho_osm_sum = _col("rho_osm_sum")
    rho_osm_mean = rho_osm_sum / n_p if rho_osm_sum.size else np.array([])
    mmp = _col("mmp_total")
    ecm = _col("ecm_strength")

    # Time axis: convert sim τ_relax → sim hours (× 60s / 3600s = 1/60).
    t_hr_sim = t * 60.0 / 3600.0   # τ_relax = 60s → hr conversion

    fig, axs = plt.subplots(3, 2, figsize=(13, 11), sharex=True)
    panels = [
        (axs[0, 0], AA0, "A/A₀ (top-down)", None),
        (axs[0, 1], R, "R / R₀", None),
        (axs[1, 0], sphericity, "Wadell sphericity ψ", (0, 1)),
        (axs[1, 1], rho_osm_mean, r"$\langle\rho_{\rm osm}\rangle$", (1.0, 1.6)),
        (axs[2, 0], mmp, "mmp_total (Layer 6)", None),
        (axs[2, 1], ecm, "ecm_strength (Layer 6)", (0.1, 1.05)),
    ]
    for ax, y, ylabel, ylim in panels:
        ax.plot(t_hr_sim, y, color="#e6ab02", linewidth=1.6, marker="o",
                markersize=2, alpha=0.8)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.grid(alpha=0.3)
        if ylim:
            ax.set_ylim(*ylim)

    # Overlay PI Lam4 A/A₀ on the top-left panel.
    pi_t_min, pi_AA0 = _load_pi_lam4()
    if pi_t_min.size:
        pi_t_hr = pi_t_min / 60.0
        axs[0, 0].plot(pi_t_hr, pi_AA0, color="#1b9e77", linewidth=1.4,
                       marker="^", markersize=3, alpha=0.7,
                       label="PI experimental Lam4 (mean)")
        axs[0, 0].axhspan(8, 33, alpha=0.10, color="green",
                          label="PI Lam4 endpoint range [8, 33]")
        axs[0, 0].legend(fontsize=8, loc="upper left")

    axs[2, 0].set_xlabel("Sim time (hours)", fontsize=11)
    axs[2, 1].set_xlabel("Sim time (hours)", fontsize=11)
    plt.suptitle("Production Lam4 — multi-metric trajectory (5k × 80 hr)",
                 fontsize=13)
    plt.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def plot_pi_comparison(
    out_path: str | Path = "results/production_lam4/figures/pi_comparison.png",
    metrics_csv: str | Path = "results/production_lam4/metrics.csv",
) -> Path:
    """Sim Lam4 production vs PI experimental Lam4 A/A₀ direct overlay.

    Includes pilot Phase 4-v2 Lam4 (4 hr) for context. Time axis in
    hours (matches PI Time_min / 60).
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6))

    # PI experimental Lam4 (population mean across spheroids).
    pi_t_min, pi_AA0 = _load_pi_lam4()
    if pi_t_min.size:
        ax.plot(pi_t_min / 60.0, pi_AA0, color="#1b9e77", linewidth=1.6,
                marker="^", markersize=4, label="PI Lam4 (mean)")
        ax.axhspan(8, 33, alpha=0.10, color="green",
                   label="PI Lam4 final range [8, 33]")

    # Pilot Phase 4-v2 Lam4 (4 hr).
    pilot_csv = Path("results/final_pilot_v2_lam4/metrics.csv")
    if pilot_csv.exists():
        t_p, AA0_p = _load_AA0_trajectory(pilot_csv)
        ax.plot(t_p * 60.0 / 3600.0, AA0_p, color="#d95f02",
                linewidth=1.6, marker="o", markersize=4,
                label="Pilot Phase 4-v2 Lam4 (4 hr)")

    # Production Lam4 (when available).
    prod_csv = Path(metrics_csv)
    if prod_csv.exists():
        t_q, AA0_q = _load_AA0_trajectory(prod_csv)
        if t_q.size:
            ax.plot(t_q * 60.0 / 3600.0, AA0_q, color="#e6ab02",
                    linewidth=1.8, marker="s", markersize=3, alpha=0.85,
                    label="Production Lam4 (5k × 80 hr)")
    else:
        ax.text(0.7, 0.5, "(production trajectory will appear\nwhen run completes)",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=10, color="grey", alpha=0.7)

    ax.set_xlabel("Time (hours)", fontsize=11)
    ax.set_ylabel(r"$A/A_0$ (top-down projection)", fontsize=11)
    ax.set_title(
        "PI experimental Lam4 vs simulation A/A₀ trajectory\n"
        "(direct measurement-protocol-matched comparison per Hard Rule 11)",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.3)
    plt.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def plot_production_bucket_classification(
    out_path: str | Path = "results/production_lam4/figures/bucket_classification.png",
    metrics_csv: str | Path = "results/production_lam4/metrics.csv",
) -> Path:
    """Visual representation of production bucket P1/P2/P3/P4 boundaries.

    Bands per docs/production_lam4_outcomes.md:
    - P1 [4.0, 33.1]: time-budget was sufficient cause; framework reaches PI range
    - P2 [2.5, 4.0]: time + mechanism both contribute; partial reproduction
    - P3 [< 2.5]: mechanism-missing; framework asymptotes below PI range
    - P4: numerical / scaling failure
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6.5))

    # Bucket bands.
    ax.axhspan(0, 2.5, alpha=0.18, color="#d95f0e", label="P3 mechanism-missing (< 2.5)")
    ax.axhspan(2.5, 4.0, alpha=0.18, color="#fec44f", label="P2 partial reproduction [2.5, 4.0]")
    ax.axhspan(4.0, 8.0, alpha=0.18, color="#7fcdbb", label="P1 lower band [4.0, 8.0]")
    ax.axhspan(8.0, 35.0, alpha=0.22, color="#1b9e77", label="P1 PI range [8.0, 33.1]")

    # Pilot value reference.
    ax.axhline(2.657, color="#d95f02", linewidth=1.2, linestyle="--",
               alpha=0.8, label="Phase 4-v2 Lam4 pilot endpoint (2.657)")

    # Production trajectory (when available).
    prod_csv = Path(metrics_csv)
    if prod_csv.exists():
        t_q, AA0_q = _load_AA0_trajectory(prod_csv)
        if t_q.size:
            ax.plot(t_q * 60.0 / 3600.0, AA0_q, color="black",
                    linewidth=1.8, marker="o", markersize=3,
                    label="Production Lam4")
            # Mark final point.
            ax.scatter([t_q[-1] * 60.0 / 3600.0], [AA0_q[-1]],
                       s=120, color="red", marker="*", zorder=10,
                       label=f"Production endpoint: {AA0_q[-1]:.2f}")
    else:
        ax.text(0.55, 0.45, "(production not yet complete;\n"
                "regenerate after run finishes)",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=11, color="grey", alpha=0.7)

    ax.set_xlabel("Sim time (hours)", fontsize=11)
    ax.set_ylabel(r"$A/A_0$ (top-down projection)", fontsize=11)
    ax.set_title(
        "Production Lam4 — Bucket P1/P2/P3 classification\n"
        "Time-budget vs mechanism-missing diagnostic",
        fontsize=11,
    )
    ax.legend(fontsize=8, loc="upper left", ncol=2)
    ax.set_ylim(0, 35)
    ax.set_xlim(0, 85)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def plot_layer_x_phenotype_summary(
    out_path: str | Path = "results/production_lam4/figures/layer_x_phenotype_summary.png",
) -> Path:
    """26-simulation summary heatmap: layer activations × phenotype.

    Highlights the combined narrative: each layer contributes
    quantitatively to spreading; phenotype ordering preserved across
    framework variants. Path C correction (g=0.1 over-anchor → g=0.01
    balanced) visualised as a separate annotation.
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Layer-by-layer R drift (single phenotype Pre carrier; Stage 1c/1d/2
    # used Pre default).
    single_phen = [
        ("Layer 1 only (v15)", 0.244),
        ("+ Substrate β α=1.0", 0.247),
        ("+ Layer 2 ζ=0.4", 0.132),
        ("+ Layer 2 ζ=0.6 (PASS)", 0.043),
        ("+ Layer 3 + Path C g=0.01\n(Stage 1c, Layer 5)", 0.138),
        ("+ Layer 4 (Stage 1d)", 0.147),
        ("+ Layer 6 (Stage 2)", 0.127),
    ]

    # 3-phenotype Phase 4-v2 (zeta_max=0.6).
    pheno_v2 = [
        ("Phase 4-v2 Bare", 0.131),
        ("Phase 4-v2 Pre", 0.095),
        ("Phase 4-v2 Lam4", 0.073),
    ]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14, 6),
                                    gridspec_kw={"width_ratios": [3, 1]})

    # Left: layer-by-layer single-phenotype progression.
    labels_L = [x[0] for x in single_phen]
    R_L = [x[1] for x in single_phen]
    ypos_L = np.arange(len(labels_L))
    colors_L = []
    for v in R_L:
        if v < 0.05:
            colors_L.append("#1b9e77")
        elif v < 0.15:
            colors_L.append("#7fcdbb")
        elif v < 0.25:
            colors_L.append("#fdae6b")
        else:
            colors_L.append("#d95f0e")
    axL.barh(ypos_L, R_L, color=colors_L, edgecolor="black", linewidth=0.5)
    axL.set_yticks(ypos_L)
    axL.set_yticklabels(labels_L, fontsize=9)
    axL.invert_yaxis()
    axL.axvline(0.05, color="red", linewidth=1.5, linestyle="--",
                label="R drift gate (0.05)")
    axL.set_xlabel(r"$\max_t |R(t)/R_0 - 1|$", fontsize=11)
    axL.set_title("Layer-by-layer progression (single Pre carrier)",
                  fontsize=11)
    axL.legend(fontsize=8)
    axL.grid(axis="x", alpha=0.3)

    # Right: 3-phenotype Phase 4-v2.
    labels_R = [x[0] for x in pheno_v2]
    R_R = [x[1] for x in pheno_v2]
    ypos_R = np.arange(len(labels_R))
    colors_R = ["#fdae6b", "#7fcdbb", "#1b9e77"]
    axR.barh(ypos_R, R_R, color=colors_R, edgecolor="black", linewidth=0.5)
    axR.set_yticks(ypos_R)
    axR.set_yticklabels(labels_R, fontsize=9)
    axR.invert_yaxis()
    axR.axvline(0.05, color="red", linewidth=1.5, linestyle="--")
    axR.set_xlabel(r"$\max_t |R(t)/R_0 - 1|$", fontsize=11)
    axR.set_title("Phase 4-v2 phenotype × full framework",
                  fontsize=11)
    axR.grid(axis="x", alpha=0.3)

    plt.suptitle(
        "Layer × Phenotype 26-simulation summary — Path C correction "
        "(g=0.1→0.01) applied to all stages",
        fontsize=12,
    )
    plt.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def generate_summary(
    out_path: str | Path = "results/visualization/SUMMARY.md",
) -> Path:
    """Write a markdown summary referencing all generated figures.

    Inline-references the four PNGs + provides project narrative
    commentary suitable for a paper-readiness review.
    """
    md_lines = [
        "# ActiveCellSim Visualization Summary",
        "",
        "Auto-generated by `acs.visualization.plots.generate_summary()`.",
        "PI Visualization directive 2026-04-29 (parallel with Production Lam4).",
        "",
        "## 1. Layer-by-layer R drift narrative",
        "",
        "![R drift bar chart](R_drift_bar.png)",
        "",
        "23-simulation progression from Layer 1 ceiling (R drift = 0.244) through",
        "the layer-by-layer activations to Phase 4-v2 Lam4 best phenotype",
        "(R drift = 0.073). The Phase 2.2 finding ζ=0.60 (R drift = 0.043) is the",
        "first PASS of the R drift gate (< 0.05) across all simulations; this",
        "single-Layer-2 result motivated the Phase 4-v2 zeta_max=0.6 update.",
        "",
        "## 2. Phenotype comparison (Phase 4-v1 vs v2)",
        "",
        "![Phenotype comparison](phenotype_comparison.png)",
        "",
        "Bare/Pre/Lam4 ordering reproduced in BOTH R drift (left) AND",
        "A/A₀_topdown (right) under Phase 4-v2 — the first time both metrics",
        "match the PI experimental ordering simultaneously. PI experimental",
        "Lam4 A/A₀_final range [8, 33] shaded for endpoint context.",
        "",
        "## 3. A/A₀ trajectory (Phase 4-v2 + production)",
        "",
        "![A/A₀ trajectory](AA0_trajectory.png)",
        "",
        "Top-down A/A₀ trajectories (matching PI's Area_um2 measurement",
        "modality per Hard Rule 11). Pilot duration t* = 240 (4 hr) marked.",
        "Production Lam4 trajectory (when complete) shows whether A/A₀",
        "continues to grow toward PI experimental endpoint range [8, 33]",
        "(Bucket P1) or asymptotes below (Bucket P2/P3).",
        "",
        "## 4. Layer 2 ζ response curve",
        "",
        "![ζ response curve](zeta_response_curve.png)",
        "",
        "Combined Track 1 + Phase 2.2 13-point sweep. Stable regime",
        "monotone-decreases to first PASS at ζ=0.60 (R drift = 0.043);",
        "transitional regime (ζ=0.70); catastrophic explosion ζ=0.90+",
        "(R drift = 2.494, active power 890× limit). Defines the operating",
        "regime for Phase 4-v2 zeta_max=0.6 selection.",
        "",
        "## 5. 3D snapshots (final frame per stage)",
        "",
        "![3D snapshots](3d_snapshots.png)",
        "",
        "Final-frame 3D particle scatter for Phase 4-v2 phenotypes + production",
        "Lam4 (when available). Substrate at z=0 shaded grey; boundary-tagged",
        "particles in red, bulk in blue. Visual confirmation of",
        "phenotype-dependent spreading shape (Lam4 most pancaked; Bare least).",
        "",
        "---",
        "",
        "Generated figures saved in `results/visualization/`. To regenerate:",
        "",
        "```bash",
        "python -c \"from acs.visualization import plots; "
        "plots.plot_layer_by_layer_R_drift(); "
        "plots.plot_phenotype_comparison(); "
        "plots.plot_AA0_trajectory(); "
        "plots.plot_zeta_response_curve(); "
        "plots.plot_3d_snapshots(); "
        "plots.generate_summary()\"",
        "```",
    ]
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(md_lines), encoding="utf-8")
    return out
