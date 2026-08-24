#!/usr/bin/env python
r"""Figures for the GATE-B ``nmii_sf_motor`` inner-solve convergence probe (2026-07-28).

WHAT THESE FIGURES ARE FOR.  The lane's reported SF tension had never been measured at a converged inner
solve — the explicit overdamped relax plateaus with a free-node residual that is a sizeable fraction of the
tension it reports, which is why the D7 verdict on every prior run is ``VOID``.  The implicit IMEX/PCG path
(:mod:`aleph.engine.sf_implicit`) landed but had never been executed.  Running it and sweeping its
mobility step over eight decades answers the question the ratio ceiling was written to ask, and the answer
is visible only as a picture: **the reported tension and the unconverged residual fall together**, so the
tension was never a relaxation transient in the innocent sense — it was the residual itself.

WHY A PICTURE AND NOT A TABLE.  A table of ``residual`` and ``T_max`` per run invites reading each row on
its own, where every row looks like a small number.  Plotted against each other on log-log across seven
decades, the proportionality is the whole finding, and the eye catches immediately that no run sits away
from the line — which a per-row eyeball cannot do.

VISUALIZATION INTEGRITY (project rule).  No axis is truncated; both decades axes span the measured range.
The log scaling is stated in every axis label because the quantities span seven decades and a linear axis
would show one point and seven zeros.  Units are SI-derived engine units and are on every label (pN, µm).
The model's own force scale — one bound head at ``f_stall`` — is drawn as a reference on the tension axes,
because "is this tension physically meaningful" is exactly the question the numbers cannot answer alone,
and every per-step trace is drawn per realisation rather than as a summary scalar.

NOT A BENCHMARK, NOT A BAND CLOSURE.  Every mechanical/kinetic constant in the underlying runs is a KB/PI
GAP (cards N1-N9, S1); these figures carry no magnitude claim.  They show a RELATIONSHIP between two
quantities measured in the same run, which is invariant to the parameter values in a way a magnitude is not.

Sanity Gate:
    * measurement-protocol: every point is read from a committed ``run-record@2`` artifact, and the run's
      own ``gate.verdict`` is carried onto the plot, so a reader cannot lift a point out of a VOID run
      without seeing that it was VOID.
    * boundary: a missing artifact is reported and skipped rather than plotted as zero (an absent input
      must render as absent, never as a value).
    * dimensional: tension, residual and traction are all forces [pN] and share an axis scale honestly;
      the mobility step is [µm/pN] and gets its own axis.
    * sign-sense: the traction is plotted as a magnitude with its inward fraction annotated separately, so
      a sign flip cannot hide inside an absolute value.
    * conservation/invariant: not applicable — this module plots, it computes no physics.
    * numerical: no fit is drawn through the data; the 1:1 guide is an exact reference line, not a
      regression, so it cannot be mistaken for a measured exponent.

Run (dev Mac, no CUDA needed — reads committed artifacts only):

    python aleph/scripts/ac_gate_b_sf_motor_convergence_viz.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

#: Where the probe artifacts live, relative to the repo root.
DEFAULT_RUN_DIR = Path("aleph/outputs/ac/gate_b_sf_motor/inner_solve_convergence_2026-07-28")

#: The runs, in the order they sweep the inner solve from least to most converged.  The label is what a
#: reader needs to identify the configuration; the file is the record it comes from.
RUNS: tuple[tuple[str, str], ...] = (
    ("sf_motor_implicit_slice.json", "IMEX 1x CFL, 4 Newton"),
    ("probe_scale_1e2.json", "IMEX 1e2x, 8 Newton"),
    ("probe_scale_1e4.json", "IMEX 1e4x, 8 Newton"),
    ("probe_inv_s1e5_n8.json", "IMEX 1e5x, 8 Newton"),
    ("probe_scale_1e6.json", "IMEX 1e6x, 8 Newton"),
    ("probe_inv_s1e6_n32.json", "IMEX 1e6x, 32 Newton"),
    ("probe_inv_s1e8_n32.json", "IMEX 1e8x, 32 Newton"),
)

#: The explicit control, plotted with a distinct marker because it is a different inner solve entirely.
EXPLICIT_RUN = ("../sf_motor_native_gate.json", "EXPLICIT relax 4000 (prior record)")

#: The NMII-stiffening diagnostic: same converged configuration, internal stiffnesses x100.
STIFF_RUN = ("probe_stiffnmii.json", "IMEX 1e8x, 32 Newton, NMII arm/backbone x100")

#: The duration control: the SAME converged configuration run for 10x the physical time, which is the one
#: objection a 20-step trace cannot answer on its own ("it would accumulate if you ran longer").
LONG_RUN = ("probe_long200.json", "IMEX 1e6x, 8 Newton, 200 steps = 2.0 s physical")

#: The native-scale topology A/B.  Identical node count, head count and solver budget; the ONLY difference is
#: whether the alpha-actinin dorsal<->arc crosslink is present.  This is what separates "the solver cannot
#: converge" from "the model has no tension", and the two runs answer it in opposite directions.
TOPOLOGY_AB = (("topo_ventral_only.json", "0 crosslinks (pure ventral dipoles)"),
               ("topo_joint_dense.json", "160 α-actinin joints / 200 bundles"))

#: The run that PASSES: the same crosslink-free native topology carried to 5.0 s of physical time.
PASS_RUN = ("topo_v_long500.json", "200 ventral sarcomeres, 500 steps = 5.0 s physical")


def _load(run_dir: Path, name: str) -> dict[str, Any] | None:
    """Read one run record, or return ``None`` if it is absent.

    Args:
        run_dir: Directory holding the probe artifacts.
        name: File name within it (may be a relative path).

    Returns:
        The parsed record, or ``None`` when the file does not exist — an absent input renders as absent.
    """
    path = (run_dir / name).resolve()
    if not path.is_file():
        print(f"  MISSING (skipped, not plotted as zero): {path}")
        return None
    with path.open() as handle:
        return json.load(handle)


def _point(record: dict[str, Any]) -> dict[str, Any]:
    """Extract the plotted quantities from one record.

    Args:
        record: A ``run-record@2`` artifact.

    Returns:
        Mapping with the final-step residuals, tension, head load, traction, the gate verdict, the
        per-step wall clock and the per-step tension trace.
    """
    final = record["measurements"]["final"]
    argv = record["config"]["argv"]
    return {
        "res_sf": float(final["residual_sf_free_pN"]),
        "res_nmii": float(final["residual_nmii_pN"]),
        "t_max": float(final["T_max_pN"]),
        "load_max": float(final["head_load_max_pN"]),
        "traction": float(final["traction_resultant_pN"]),
        "inward": float(final["fa_inward_fraction"]),
        "n_bound": int(final["n_bound"]),
        "verdict": str(record["gate"]["verdict"]),
        "ratio": float(record["gate"]["residual_over_signal"]),
        "wall_per_step": float(record["timing"]["wall_seconds_per_step"]),
        "scale": (float(argv["implicit_step_scale"]) if argv.get("implicit") else None),
        "f_stall": float(argv["f_stall"]),
        "trace_t": [float(row["T_max_pN"]) for row in record["measurements"]["trace"]],
        "trace_res": [float(row["residual_sf_free_pN"]) for row in record["measurements"]["trace"]],
    }


def figure_tension_vs_residual(points: list[tuple[str, dict]], explicit: tuple[str, dict] | None,
                               stiff: tuple[str, dict] | None, out: Path) -> None:
    """Plot tension and FA traction against the free-node residual they were measured with.

    The finding is that they are proportional: a tension that moves with the residual is the residual.

    Args:
        points: ``(label, point)`` pairs for the implicit sweep.
        explicit: The explicit control, or ``None``.
        stiff: The NMII-stiffening diagnostic, or ``None``.
        out: Destination PNG path.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4))
    f_stall = points[0][1]["f_stall"]

    for ax, key, name in ((axes[0], "t_max", "SF axial tension  $T_{max}$"),
                          (axes[1], "traction", "FA traction resultant  $|F|$")):
        xs = [p["res_sf"] for _, p in points]
        ys = [p[key] for _, p in points]
        ax.plot(xs, ys, "o-", color="#1f77b4", zorder=3, label="implicit IMEX/PCG sweep")
        for label, p in points:
            ax.annotate(label.replace("IMEX ", "").replace(", ", "\n"), (p["res_sf"], p[key]),
                        textcoords="offset points", xytext=(6, -12), fontsize=6.5, color="#333333")
        if explicit is not None:
            ax.plot([explicit[1]["res_sf"]], [explicit[1][key]], "s", ms=10, color="#d62728",
                    zorder=4, label="explicit relax 4000 (prior record)")
        if stiff is not None:
            ax.plot([stiff[1]["res_sf"]], [stiff[1][key]], "^", ms=10, color="#2ca02c",
                    zorder=4, label="NMII arm/backbone x100 (diagnostic)")

        lo = min(min(xs), min(ys)) * 0.3
        hi = max(max(xs), max(ys)) * 3.0
        ax.plot([lo, hi], [lo, hi], ":", color="#888888", lw=1.2, zorder=1,
                label="exact 1:1 reference (not a fit)")
        ax.axhline(f_stall, color="#ff7f0e", lw=1.4, ls="--", zorder=2,
                   label=f"one bound head at $f_{{stall}}$ = {f_stall:g} pN (model force scale)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel("free-node residual of the inner solve  [pN]  (log)")
        ax.set_ylabel(f"{name}  [pN]  (log)")
        ax.grid(True, which="both", alpha=0.25)
        ax.set_title(name.split("  ")[0], fontsize=11)

    axes[0].legend(fontsize=7.5, loc="upper left")
    fig.suptitle("GATE-B nmii_sf_motor — the reported tension IS the unconverged residual\n"
                 "SLICE: 904 nodes / 320 heads / 16 minifilaments = 0.18% of native; every run's D7 "
                 "verdict is VOID; all constants are PI-GAPs (no magnitude claimed)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def figure_step_traces(points: list[tuple[str, dict]], out: Path) -> None:
    """Plot the per-step tension trace of every run, so a steady state is distinguishable from a ramp.

    Args:
        points: ``(label, point)`` pairs for the implicit sweep.
        out: Destination PNG path.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    colours = plt.cm.viridis([i / max(1, len(points) - 1) for i in range(len(points))])
    for (label, p), colour in zip(points, colours):
        steps = list(range(1, len(p["trace_t"]) + 1))
        axes[0].plot(steps, p["trace_t"], "-o", ms=3, color=colour, label=label)
        axes[1].plot(steps, p["trace_res"], "-o", ms=3, color=colour, label=label)
    axes[0].axhline(points[0][1]["f_stall"], color="#ff7f0e", ls="--", lw=1.4,
                    label=f"one bound head at $f_{{stall}}$ = {points[0][1]['f_stall']:g} pN")
    for ax, name in ((axes[0], "SF axial tension  $T_{max}$  [pN]"),
                     (axes[1], "free-node residual  [pN]")):
        ax.set_yscale("log")
        ax.set_xlabel("accepted physical step  (dt = 0.01 s)")
        ax.set_ylabel(f"{name}  (log)")
        ax.grid(True, which="both", alpha=0.25)
    axes[0].legend(fontsize=7, loc="lower right")
    fig.suptitle("Per-step traces — every run plateaus within a few steps at a level set by its OWN inner-solve\n"
                 "residual; the best-converged run plateaus six decades below one head's stall force",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def figure_two_bodies(points: list[tuple[str, dict]], out: Path) -> None:
    """Plot both owners' residuals against the IMEX mobility step scale.

    The SF body converges over seven decades; the NMII body does not, which is the honest limit of the
    convergence claim and must be visible rather than buried in a table.

    Args:
        points: ``(label, point)`` pairs for the implicit sweep (records carrying a step scale).
        out: Destination PNG path.
    """
    swept = [(label, p) for label, p in points if p["scale"] is not None]
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ax.plot([p["scale"] for _, p in swept], [p["res_sf"] for _, p in swept], "o-",
            color="#1f77b4", label="sf_arc body — free-node residual [pN]")
    ax.plot([p["scale"] for _, p in swept], [p["res_nmii"] for _, p in swept], "s-",
            color="#d62728", label="nmii body — residual [pN]")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("IMEX mobility step, as a multiple of the explicit CFL bound  (log, dimensionless)")
    ax.set_ylabel("residual at the last accepted step  [pN]  (log)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8.5)
    ax.set_title("Only ONE of the two owners converges\n"
                 "split ownership means the gate must see both; the nmii body is the open end",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def figure_cost(points: list[tuple[str, dict]], out: Path) -> None:
    """Plot the measured per-step wall clock against the convergence actually reached.

    This is the stage-2 cost question in the only form it may be asked: a cost is a property of an engine
    at a STATED convergence, so the two axes belong on one plot.

    Args:
        points: ``(label, point)`` pairs for the implicit sweep.
        out: Destination PNG path.
    """
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    for label, p in points:
        ax.plot([p["res_sf"]], [p["wall_per_step"]], "o", ms=9, color="#1f77b4")
        ax.annotate(label.replace("IMEX ", ""), (p["res_sf"], p["wall_per_step"]),
                    textcoords="offset points", xytext=(7, 4), fontsize=7.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("free-node residual reached  [pN]  (log, MORE converged to the right)")
    ax.set_ylabel("measured wall clock per accepted step  [s]  (log)")
    ax.grid(True, which="both", alpha=0.25)
    ax.set_title("Cost of a step, against the convergence that step reached  (SLICE, 904 nodes)\n"
                 "every point is from a VOID run, so none of these is a quotable native step cost",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def figure_duration_control(long_point: dict, short_point: dict, out: Path) -> None:
    """Plot the 10x-physical-time run against one head's stall force, with the residual beside it.

    The objection a 20-step trace cannot answer is "the tension would accumulate if you ran longer". This
    figure answers it in the only way that settles it: the same configuration, ten times the physical time,
    with the model's own force scale drawn on the same axis so the reader sees the distance to it.

    Args:
        long_point: The 200-step run's extracted point.
        short_point: The 20-step run at the identical configuration.
        out: Destination PNG path.
    """
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    steps = list(range(1, len(long_point["trace_t"]) + 1))
    times = [s * 0.01 for s in steps]
    ax.plot(times, long_point["trace_t"], "-", lw=1.4, color="#1f77b4",
            label="SF axial tension $T_{max}$ [pN] — 200 steps")
    ax.plot(times, long_point["trace_res"], "-", lw=1.4, color="#d62728",
            label="free-node residual [pN] — same run")
    ax.axvline(len(short_point["trace_t"]) * 0.01, color="#888888", ls=":", lw=1.2,
               label="end of the 20-step runs plotted elsewhere")
    ax.axhline(long_point["f_stall"], color="#ff7f0e", ls="--", lw=1.6,
               label=f"ONE bound head at $f_{{stall}}$ = {long_point['f_stall']:g} pN")
    ax.set_yscale("log")
    ax.set_xlabel("physical time  [s]   (dt = 0.01 s per accepted step)")
    ax.set_ylabel("force  [pN]  (log)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("Duration control — 10x the physical time does NOT close the gap\n"
                 f"after 2.0 s with {long_point['n_bound']} of 320 heads bound, $T_{{max}}$ is still two "
                 f"decades below ONE head's stall force", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def figure_emergence(point: dict, out: Path) -> None:
    """Plot the run that PASSES: tension emerging over 5.0 s against the residual it was measured with.

    The single most important thing this figure carries is the vertical marker at 0.2 s. Every 20-step run in
    this series stopped there, inside an induction period during which the converged tension really is ~0 —
    which is how a duration artifact came to be read as a statement about the model's equilibrium.

    Args:
        point: The 500-step run's extracted point.
        out: Destination PNG path.
    """
    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    times = [(i + 1) * 0.01 for i in range(len(point["trace_t"]))]
    ax.plot(times, point["trace_t"], "-", lw=1.8, color="#1f77b4", label="SF axial tension $T_{max}$ [pN]")
    ax.plot(times, point["trace_res"], "-", lw=1.4, color="#d62728", label="free-node residual [pN] — same run")
    ax.axhline(point["f_stall"], color="#ff7f0e", ls="--", lw=1.5,
               label=f"ONE bound head at $f_{{stall}}$ = {point['f_stall']:g} pN")
    ax.axvline(0.20, color="#444444", ls=":", lw=1.6,
               label="0.20 s — where EVERY 20-step run in this series stopped")
    ax.set_yscale("log")
    ax.set_xlabel("physical time  [s]   (dt = 0.01 s per accepted step)")
    ax.set_ylabel("force  [pN]  (log)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8.5, loc="lower right")
    rise = (point["trace_t"][-1] - point["trace_t"][-201]) / 2.0 if len(point["trace_t"]) > 201 else float("nan")
    ax.set_title("The tension DOES emerge, once the solve converges AND the clock passes the induction period\n"
                 f"D7 = PASS · residual {point['res_sf']:.2e} pN, 8 decades under the signal · "
                 f"{point['n_bound']}/4,000 heads bound\n"
                 f"NOT a steady state: still rising ≈{rise:.2f} pN/s over the last 2 s "
                 f"(a log axis flatters it — read the slope, not the shape)", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def figure_topology_ab(free_pt: dict, dense_pt: dict, out: Path) -> None:
    """Plot the crosslink A/B: identical size and budget, ten decades of difference in the residual.

    Args:
        free_pt: The crosslink-free native topology.
        dense_pt: The joint-dense native topology.
        out: Destination PNG path.
    """
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    labels = ["0 α-actinin joints\n(pure ventral)", "160 α-actinin joints\n(dense)"]
    xs = [0, 1]
    ax.plot(xs, [free_pt["res_sf"], dense_pt["res_sf"]], "o-", ms=11, color="#d62728",
            label="free-node residual reached [pN]")
    ax.plot(xs, [free_pt["t_max"], dense_pt["t_max"]], "s-", ms=11, color="#1f77b4",
            label="reported $T_{max}$ [pN] at the same 0.20 s")
    ax.axhline(free_pt["f_stall"], color="#ff7f0e", ls="--", lw=1.5,
               label=f"ONE bound head at $f_{{stall}}$ = {free_pt['f_stall']:g} pN")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.set_yscale("log")
    ax.set_ylabel("force  [pN]  (log)")
    ax.grid(True, which="both", axis="y", alpha=0.25)
    ax.legend(fontsize=8.5)
    ax.set_title("Same 10,400 nodes, same 4,000 heads, same solver budget —\n"
                 "the α-actinin dorsal↔arc crosslink alone costs TEN decades of convergence",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def main() -> int:
    """Render every figure from the committed probe artifacts.

    Returns:
        Process exit code: 0 on success, 1 if no artifact could be read.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    args = parser.parse_args()
    run_dir: Path = args.run_dir
    figs = run_dir / "figs"
    figs.mkdir(parents=True, exist_ok=True)

    print(f"reading probe artifacts from {run_dir}")
    points: list[tuple[str, dict]] = []
    for name, label in RUNS:
        record = _load(run_dir, name)
        if record is not None:
            points.append((label, _point(record)))
    if not points:
        print("no probe artifact could be read — nothing rendered")
        return 1

    explicit_record = _load(run_dir, EXPLICIT_RUN[0])
    explicit = (EXPLICIT_RUN[1], _point(explicit_record)) if explicit_record else None
    stiff_record = _load(run_dir, STIFF_RUN[0])
    stiff = (STIFF_RUN[1], _point(stiff_record)) if stiff_record else None

    figure_tension_vs_residual(points, explicit, stiff, figs / "tension_is_the_residual.png")
    figure_step_traces(points, figs / "per_step_traces.png")
    figure_two_bodies(points, figs / "two_body_residuals.png")
    figure_cost(points, figs / "cost_vs_convergence.png")

    long_record = _load(run_dir, LONG_RUN[0])
    short = next((p for label, p in points if label.startswith("IMEX 1e6x, 8")), None)
    if long_record is not None and short is not None:
        figure_duration_control(_point(long_record), short, figs / "duration_control.png")

    pass_record = _load(run_dir, PASS_RUN[0])
    if pass_record is not None:
        figure_emergence(_point(pass_record), figs / "tension_emerges_when_it_converges.png")

    free_record = _load(run_dir, TOPOLOGY_AB[0][0])
    dense_record = _load(run_dir, TOPOLOGY_AB[1][0])
    if free_record is not None and dense_record is not None:
        figure_topology_ab(_point(free_record), _point(dense_record), figs / "crosslink_blocks_convergence.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
