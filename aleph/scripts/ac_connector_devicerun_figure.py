#!/usr/bin/env python
r"""Render the per-connector device-run census beside its record.

One figure, two panels, and both are counts — there is no magnitude in this census to plot and drawing
one would invent evidence the run does not carry.

  LEFT   every declared connector as one cell, coloured by verdict, in architecture order. The whole
         population is shown; nothing is ranked, truncated or sampled, so "37 of 38 are one colour" is
         read off the picture rather than taken on trust.
  RIGHT  the verdict counts as a bar chart on a LINEAR axis starting at zero. No log axis: the counts
         span 1..38 and a log axis would make 1 and 37 look comparable, which is the entire finding.

The subtitle carries the population, the device and the control outcome, because a census whose controls
did not pass is not a census and the figure must not be separable from that fact.

Runs on the dev machine from the committed record — it reads JSON and touches no device.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

#: Verdict → colour. Ordered weakest-to-strongest so the legend reads as a ladder.
_ORDER = ("NOT_BOUND", "NOT_INSTRUMENTABLE", "NOT_CALLED", "STUB", "DEVICE_RUN")
_COLOUR = {
    "NOT_BOUND": "#d9d9d9",
    "NOT_INSTRUMENTABLE": "#bdbdbd",
    "NOT_CALLED": "#fdd0a2",
    "STUB": "#fdae6b",
    "DEVICE_RUN": "#2c7fb8",
}


def render(record_path: Path, out_path: Path, negative_path: Path | None = None) -> None:
    """Write the two-panel census figure for one ``run-record@2`` census record."""
    rec = json.loads(record_path.read_text(encoding="utf-8"))
    rows = rec["connectors"]
    counts = {v: sum(1 for r in rows if r["verdict"] == v) for v in _ORDER}

    controls = rec.get("controls", {})
    pos = controls.get("positive", {})
    neg_note = "negative: not in this pass"
    if negative_path and negative_path.is_file():
        nrec = json.loads(negative_path.read_text(encoding="utf-8"))
        nneg = nrec.get("controls", {}).get("negative", {})
        neg_note = (f"negative (self-null {nneg.get('connector')}): "
                    f"{nneg.get('verdict')} — {'PASS' if nneg.get('ok') else 'FAIL'}")

    fig, (ax_map, ax_bar) = plt.subplots(
        1, 2, figsize=(15.5, 7.2), gridspec_kw={"width_ratios": [2.6, 1.0]},
    )

    # LEFT — every declared connector, architecture order, no omission.
    n_cols = 4
    for i, row in enumerate(rows):
        col, line = i % n_cols, i // n_cols
        ax_map.add_patch(mpatches.Rectangle(
            (col, -line), 0.96, 0.92, facecolor=_COLOUR.get(row["verdict"], "#ffffff"),
            edgecolor="#4d4d4d", linewidth=0.6,
        ))
        label = row["connector"]
        ax_map.text(col + 0.48, -line + 0.46, label if len(label) <= 26 else label[:25] + "…",
                    ha="center", va="center", fontsize=7.1,
                    color="white" if row["verdict"] == "DEVICE_RUN" else "#1a1a1a")
    n_lines = (len(rows) + n_cols - 1) // n_cols
    ax_map.set_xlim(-0.1, n_cols + 0.05)
    ax_map.set_ylim(-n_lines + 0.9, 1.05)
    ax_map.axis("off")
    ax_map.set_title(f"All {len(rows)} declared connectors — none omitted, none ranked", fontsize=11)
    ax_map.legend(
        handles=[mpatches.Patch(facecolor=_COLOUR[v], edgecolor="#4d4d4d", label=f"{v} ({counts[v]})")
                 for v in _ORDER],
        loc="lower center", bbox_to_anchor=(0.5, -0.10), ncol=3, fontsize=8.5, frameon=False,
    )

    # RIGHT — linear counts from zero. Units: connectors (a count, not a magnitude).
    present = [v for v in _ORDER if counts[v]]
    ax_bar.barh(range(len(present)), [counts[v] for v in present],
                color=[_COLOUR[v] for v in present], edgecolor="#4d4d4d", linewidth=0.6)
    ax_bar.set_yticks(range(len(present)), present, fontsize=9)
    ax_bar.invert_yaxis()
    ax_bar.set_xlim(0, len(rows))
    ax_bar.set_xlabel("declared connectors [count]", fontsize=9.5)
    ax_bar.set_title("Verdict counts (linear axis from 0)", fontsize=11)
    for i, v in enumerate(present):
        ax_bar.text(counts[v] + 0.5, i, str(counts[v]), va="center", fontsize=9.5)
    ax_bar.spines[["top", "right"]].set_visible(False)

    census = rec["census"]
    fig.suptitle(
        "Per-connector DEVICE-RUN census — composed native cell\n"
        f"{census['cortex_filaments']:,} cortical filaments / {census['n_actin_nodes']:,} actin nodes / "
        f"{census['n_total_nodes']:,} total nodes ({census['fraction_of_native']:.0%} of native) · "
        f"{rec['device']} · build {rec['build']['commit']} ({rec['build']['source']})\n"
        f"components bound: {len(rec['registered_components'])}/14 "
        f"{rec['registered_components']} · "
        f"positive control ({pos.get('connector')}): {pos.get('verdict')} — "
        f"{'PASS' if pos.get('ok') else 'FAIL'} · {neg_note}",
        fontsize=10, y=0.99,
    )
    fig.text(0.5, 0.015,
             "Counts only — this census carries no force magnitude, and the verdict is the kernel-LAUNCH "
             "observation alone. Verdicts describe THIS composed configuration on THIS build.",
             ha="center", fontsize=8.5, style="italic", color="#4d4d4d")
    fig.tight_layout(rect=(0, 0.035, 1, 0.88))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    print(f"wrote {out_path}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Render the connector device-run census figure.")
    ap.add_argument("--record", required=True)
    ap.add_argument("--negative-record", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    render(Path(args.record), Path(args.out),
           Path(args.negative_record) if args.negative_record else None)


if __name__ == "__main__":
    main()
