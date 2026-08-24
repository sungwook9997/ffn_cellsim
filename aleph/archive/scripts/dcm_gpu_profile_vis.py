"""Figures for the GPU DCM per-step profile + active-traction vectorization win.

Reads the JSON artifacts written by ``dcm_gpu_profile.py`` (per-step breakdown)
and ``dcm_gpu_ab_wall.py`` (A/B wall: original loop vs vectorized vs fused) from
``outputs/h_dcm_gpu_lod/profile/`` and renders:

  (1) per-step time breakdown (stacked: each custom force + baoab + HOOMD internals)
      for the active (production) profile — shows the active traction dominating;
  (2) A/B wall bars (ms/step) for original-loop / vectorized / fused with speedups.

Run on gbook (after the profile + ab_wall runs) or wherever the JSONs landed:
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python scripts/dcm_gpu_profile_vis.py
"""

from __future__ import annotations

import json
from pathlib import Path

_DIR = Path("aleph/outputs/h_dcm_gpu_lod/profile")
_FIGS = _DIR / "figs"


def _load(name: str):
    fp = _DIR / name
    return json.loads(fp.read_text()) if fp.exists() else None


def main() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _FIGS.mkdir(parents=True, exist_ok=True)
    # prefer the tagged active profile; fall back to the legacy untagged name.
    prof = (_load("profile_n200_active_GPU.json")
            or _load("profile_n200_GPU.json"))
    ab = _load("ab_wall_n200_GPU.json")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (1) per-step breakdown
    ax = axes[0]
    if prof is not None:
        ps = prof["per_step_ms"]
        order = ["set_forces_active", "set_forces_tent", "set_forces_turgor",
                 "set_forces_substrate", "baoab_act"]
        labels, vals, colors = [], [], []
        palette = {"set_forces_active": "#d7301f", "set_forces_tent": "#fc8d59",
                   "set_forces_turgor": "#fdcc8a", "set_forces_substrate": "#fef0d9",
                   "baoab_act": "#74a9cf"}
        for k in order:
            if k in ps:
                labels.append(k.replace("set_forces_", ""))
                vals.append(ps[k])
                colors.append(palette[k])
        labels.append("HOOMD internals")
        vals.append(prof["hoomd_internal_ms"])
        colors.append("#bdbdbd")
        ax.bar(labels, vals, color=colors)
        ax.set_ylabel("ms / step")
        ax.set_title(f"GPU per-step breakdown (A5000, N={prof['n_cells']}, "
                     f"{prof['nodes']} nodes)\ntotal "
                     f"{prof['total_ms_per_step']:.1f} ms/step  |  empty gpu_local "
                     f"ctx floor {prof['empty_ctx_floor_ms']:.2f} ms/force")
        ax.tick_params(axis="x", rotation=30)
        for i, v in enumerate(vals):
            ax.text(i, v, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
    else:
        ax.text(0.5, 0.5, "profile_n200_GPU.json not found", ha="center")

    # (2) A/B wall
    ax = axes[1]
    if ab is not None:
        names = ["A: original\nloop active", "B: vectorized\nactive",
                 "C: fused\n(vec)"]
        vals = [ab["A_original_loop"]["ms_per_step"],
                ab["B_vectorized"]["ms_per_step"],
                ab["C_fused"]["ms_per_step"]]
        colors = ["#bdbdbd", "#2c7fb8", "#238443"]
        ax.bar(names, vals, color=colors)
        ax.set_ylabel("ms / step")
        spB = ab["B_vectorized"]["speedup_vs_A"]
        spC = ab["C_fused"]["speedup_vs_A"]
        ax.set_title(f"A/B wall (A5000, N={ab['n_cells']}, {ab['steps']} steps)\n"
                     f"vectorize {spB:.2f}x  |  +fuse total {spC:.2f}x")
        for i, v in enumerate(vals):
            ax.text(i, v, f"{v:.1f}", ha="center", va="bottom", fontsize=9)
    else:
        ax.text(0.5, 0.5, "ab_wall_n200_GPU.json not found", ha="center")

    fig.suptitle("GPU DCM per-step profile + active-traction vectorization "
                 "(gbook RTX A5000, 2026-06-11)")
    fig.tight_layout()
    out = _FIGS / "dcm_gpu_profile_ab.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"[fig] {out}")


if __name__ == "__main__":
    main()
