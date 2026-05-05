"""B-tier visual evidence dashboard for the Phase E v2 pilot runner.

This module is **read-only plumbing**. It consumes the artifacts emitted
by ``acs.v2.phase_e_v2_pilot_runner`` (commit ``38b679f``) and produces
a multi-panel PNG dashboard plus an optional self-contained HTML rollup.

Per Codex ``id=1694`` 5 guardrails (B-tier scope discipline):

1. Consume only runner artifacts: ``metadata.json``, ``diagnostics.csv``,
   and (optionally) ``frame_*.h5`` via the existing
   :func:`acs.v2.output.frame_dump.read_frame` path.
2. No new biological / physics claims. Plots are labeled as
   *existing diagnostic summaries*, NOT validation metrics.
3. Tests assert artifact-schema correctness (column names / count / file
   non-emptiness), NOT pixel-perfect image equality.
4. CLI smoke generates non-empty PNG / HTML outputs from a tiny temp run.
5. No new measurement extraction from frames — frame replay is read-only
   visualization only.

The dashboard does NOT introduce new physics, NOT establish satisfaction
claims, and NOT produce content that should be cited as scientific
evidence beyond the ``diagnostics.csv`` schema already documented in the
runner module. It is a viewer over runner outputs, not an analysis layer.
"""

from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass
from typing import Optional

import matplotlib

matplotlib.use("Agg")  # headless / SSH-friendly
import matplotlib.pyplot as plt  # noqa: E402


_DIAGNOSTIC_COLUMNS = (
    "step_index",
    "time_s",
    "v_active_um2",
    "multiplier_min",
    "multiplier_max",
    "multiplier_mean",
    "max_traction_norm_nN_per_um2",
    "orientation_tensor_abs_max",
)


_DASHBOARD_PNG_NAME = "dashboard.png"
_DASHBOARD_HTML_NAME = "dashboard.html"


@dataclass(frozen=True, slots=True)
class PilotDashboardArtifacts:
    """File paths emitted by :func:`render_phase_e_v2_pilot_dashboard`."""

    run_dir: str
    png_path: str
    html_path: Optional[str]
    n_diagnostic_rows: int
    n_frames: int


def _read_diagnostics_csv(csv_path: str) -> list[dict[str, float]]:
    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"diagnostics.csv missing at {csv_path!r}; expected runner artifact"
        )
    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError(f"{csv_path!r} has no header row")
        missing = set(_DIAGNOSTIC_COLUMNS) - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"{csv_path!r} missing expected diagnostic columns: "
                f"{sorted(missing)}; runner artifact schema mismatch"
            )
        rows: list[dict[str, float]] = []
        for raw in reader:
            row: dict[str, float] = {}
            for col in _DIAGNOSTIC_COLUMNS:
                if col == "step_index":
                    row[col] = float(int(raw[col]))
                else:
                    row[col] = float(raw[col])
            rows.append(row)
    return rows


def _read_metadata_json(metadata_path: str) -> dict[str, object]:
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(
            f"metadata.json missing at {metadata_path!r}; expected runner artifact"
        )
    with open(metadata_path, encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"{metadata_path!r} top-level is not a JSON object")
    runner = data.get("runner")
    accepted_runners = ("phase_e_v2_pilot", "phase_f_minimal_motility_pilot")
    if runner not in accepted_runners:
        raise ValueError(
            f"{metadata_path!r} runner field {runner!r} is not in "
            f"{accepted_runners}; refusing to render unrelated artifacts"
        )
    return data


def _series(rows: list[dict[str, float]], col: str) -> list[float]:
    return [r[col] for r in rows]


def _render_png(
    rows: list[dict[str, float]], png_path: str, *, run_label: str
) -> None:
    if not rows:
        # Render an explicit "no data" placeholder rather than raising,
        # so a zero-step pilot run still produces a valid PNG artifact.
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        ax.text(
            0.5,
            0.5,
            "no diagnostic rows recorded\n(zero-step pilot run)",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_axis_off()
        fig.suptitle(
            f"Phase E v2 pilot · {run_label} (display-only summary)",
            fontsize=10,
        )
        fig.tight_layout()
        fig.savefig(png_path, dpi=120)
        plt.close(fig)
        return

    steps = _series(rows, "step_index")
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), sharex=True)

    ax = axes[0, 0]
    ax.plot(steps, _series(rows, "v_active_um2"), marker="o")
    ax.set_ylabel("v_active_um2")
    ax.set_title("HB#5 Lyapunov-like (display-only diagnostic)", fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(steps, _series(rows, "multiplier_min"), marker="o", label="min")
    ax.plot(steps, _series(rows, "multiplier_mean"), marker="s", label="mean")
    ax.plot(steps, _series(rows, "multiplier_max"), marker="^", label="max")
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=0.8, label="neutral=1.0")
    ax.set_ylabel("ECM->FA multiplier")
    ax.set_title("HB#4-active multiplier summary (display-only)", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc="best")

    ax = axes[1, 0]
    ax.plot(steps, _series(rows, "max_traction_norm_nN_per_um2"), marker="o")
    ax.set_ylabel("max traction (nN/um^2)")
    ax.set_xlabel("step_index")
    ax.set_title("HB#3 scatter peak (display-only diagnostic)", fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(steps, _series(rows, "orientation_tensor_abs_max"), marker="o")
    ax.set_ylabel("|T|_inf")
    ax.set_xlabel("step_index")
    ax.set_title("ECM orientation L-inf (display-only diagnostic)", fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"Phase E v2 pilot · {run_label} (display-only diagnostic summary, "
        f"NOT a validation metric)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(png_path, dpi=120)
    plt.close(fig)


def _render_html(
    rows: list[dict[str, float]],
    metadata: dict[str, object],
    html_path: str,
    *,
    png_basename: str,
    run_label: str,
) -> None:
    config = metadata.get("config", {})
    notes = metadata.get("notes", [])
    n_frames = metadata.get("n_frames", 0)
    git_commit_hash = metadata.get("git_commit_hash", "unrecorded")
    final = metadata.get("final")

    config_rows = "".join(
        f"<tr><th>{k}</th><td>{v}</td></tr>"
        for k, v in (config.items() if isinstance(config, dict) else [])
    )
    notes_items = "".join(
        f"<li>{n}</li>" for n in (notes if isinstance(notes, list) else [])
    )
    final_block = ""
    if isinstance(final, dict):
        final_rows = "".join(
            f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in final.items()
        )
        final_block = (
            "<h3>Final-step diagnostic snapshot</h3>"
            f"<table>{final_rows}</table>"
        )
    diagnostic_table_rows = "".join(
        "<tr>"
        + "".join(f"<td>{row[col]:g}</td>" for col in _DIAGNOSTIC_COLUMNS)
        + "</tr>"
        for row in rows
    )
    diagnostic_header = "".join(f"<th>{c}</th>" for c in _DIAGNOSTIC_COLUMNS)
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Phase E v2 pilot dashboard - {run_label}</title>
<style>
body {{ font-family: -apple-system, sans-serif; margin: 24px; color: #222; }}
h1 {{ font-size: 18px; }}
h2, h3 {{ font-size: 14px; margin-top: 18px; }}
table {{ border-collapse: collapse; font-size: 12px; }}
th, td {{ border: 1px solid #bbb; padding: 4px 8px; text-align: left; }}
.note {{ color: #777; font-size: 11px; }}
img {{ max-width: 100%; border: 1px solid #ddd; }}
</style>
</head>
<body>
<h1>Phase E v2 pilot dashboard - {run_label}</h1>
<p class="note">Display-only diagnostic summary. NOT a validation metric.
NOT scientific evidence beyond the runner's CSV schema. B-tier visualization
over locked composition wrapper output.</p>

<h2>Run metadata</h2>
<table><tr><th>git_commit_hash</th><td>{git_commit_hash}</td></tr>
<tr><th>n_diagnostic_rows</th><td>{len(rows)}</td></tr>
<tr><th>n_frames</th><td>{n_frames}</td></tr></table>

<h2>Configuration</h2>
<table>{config_rows}</table>

<h2>Trajectory panel</h2>
<img src="{png_basename}" alt="Phase E v2 pilot trajectory panel">

<h2>Per-step diagnostics</h2>
<table><tr>{diagnostic_header}</tr>{diagnostic_table_rows}</table>

{final_block}

<h2>Scope notes</h2>
<ul>{notes_items}</ul>
</body>
</html>
"""
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(html)


def render_phase_e_v2_pilot_dashboard(
    run_dir: str,
    *,
    write_html: bool = True,
) -> PilotDashboardArtifacts:
    """Render dashboard PNG (and optional HTML) from a pilot runner's run dir.

    Reads ``<run_dir>/diagnostics.csv`` and ``<run_dir>/metadata.json`` and
    writes ``<run_dir>/dashboard.png`` (and optionally
    ``<run_dir>/dashboard.html``). The PNG and HTML are display-only
    artifacts; they do not modify any runner output.

    Args:
        run_dir: A directory previously produced by
            :func:`acs.v2.phase_e_v2_pilot_runner.run_phase_e_v2_pilot`.
        write_html: If True, also emit a self-contained HTML rollup that
            embeds the trajectory PNG by relative path and shows the run
            metadata + per-step diagnostic table.

    Returns:
        :class:`PilotDashboardArtifacts` with the absolute paths of the
        generated dashboard files.
    """

    if not os.path.isdir(run_dir):
        raise FileNotFoundError(
            f"run_dir {run_dir!r} is not an existing directory"
        )
    abs_run_dir = os.path.abspath(run_dir)
    csv_path = os.path.join(abs_run_dir, "diagnostics.csv")
    metadata_path = os.path.join(abs_run_dir, "metadata.json")

    rows = _read_diagnostics_csv(csv_path)
    metadata = _read_metadata_json(metadata_path)

    run_label = os.path.basename(abs_run_dir.rstrip(os.sep)) or abs_run_dir
    png_path = os.path.join(abs_run_dir, _DASHBOARD_PNG_NAME)
    _render_png(rows, png_path, run_label=run_label)

    html_path: Optional[str] = None
    if write_html:
        html_path = os.path.join(abs_run_dir, _DASHBOARD_HTML_NAME)
        _render_html(
            rows,
            metadata,
            html_path,
            png_basename=_DASHBOARD_PNG_NAME,
            run_label=run_label,
        )

    n_frames = int(metadata.get("n_frames", 0)) if isinstance(metadata, dict) else 0
    return PilotDashboardArtifacts(
        run_dir=abs_run_dir,
        png_path=png_path,
        html_path=html_path,
        n_diagnostic_rows=len(rows),
        n_frames=n_frames,
    )
