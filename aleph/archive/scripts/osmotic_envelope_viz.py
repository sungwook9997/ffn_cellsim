"""Render the osmotic / volume-response prediction-envelope oracle to a figure.

Pure-Python (NumPy + matplotlib), GPU-agnostic — this is the ORACLE layer, not the runtime.
Emits `aleph/outputs/ac/osmotic/osmotic_envelope.png` with four panels, each following the
CLAUDE.md viz-integrity rules: per-realisation thin lines + band overlay, units annotated, no
axis truncation, literature-range envelopes (never a single deterministic curve).

Usage:  python -m aleph.scripts.osmotic_envelope_viz
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from aleph.validation.oracles.osmotic_response import (  # noqa: E402
    B_GRID,
    CONDITIONS,
    permeant_envelope,
    pbvh_v_eq,
    shape_envelope,
    stiffness_envelope,
    timecourse_envelope,
)

OUT = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "osmotic"


def _band(ax, env, color, label, alpha_fill=0.18):
    for row in env.curves:
        ax.plot(env.x, row, color=color, lw=0.5, alpha=0.35)
    ax.fill_between(env.x, env.lo, env.hi, color=color, alpha=alpha_fill,
                    label=f"{label} band")
    ax.plot(env.x, env.med, color=color, lw=1.8)


def main() -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # (1) PBvH passive equilibrium vs r, swept over b -------------------------
    ax = axes[0, 0]
    r = np.linspace(0.6, 1.5, 100)
    for b in B_GRID:
        ax.plot(r, pbvh_v_eq(r, b), lw=1.4, label=f"b={b}")
    ax.fill_between(r, pbvh_v_eq(r, min(B_GRID)), pbvh_v_eq(r, max(B_GRID)),
                    color="C0", alpha=0.15)
    ax.axhline(1.0, color="k", lw=0.6, ls=":")
    ax.axvline(1.0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("medium ratio  r = O_out/O₀  [–]")
    ax.set_ylabel("v_eq = V_eq/V₀  [–]")
    ax.set_title("Ponder–Boyle–van 't Hoff passive equilibrium\n(band = b ∈ {0.1,0.2,0.3})")
    ax.legend(fontsize=8)

    # (2) impermeant time course for strong hyper, RVD/RVI fan ----------------
    ax = axes[0, 1]
    t = np.linspace(0.0, 8.0, 200)
    env = timecourse_envelope(t, CONDITIONS["strong_hyper"])
    _band(ax, env, "C3", "r=1.30 (b,R,τ_reg/τ_w)")
    ax.axhline(1.0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("time  t/τ_w  [–]")
    ax.set_ylabel("v(t) = V/V₀  [–]")
    ax.set_title("Two-timescale RVD/RVI time course, strong hyper\n(fan spans R=0 → 1 recovery)")
    ax.legend(fontsize=8)

    # (3) permeant solute shrink-then-re-swell --------------------------------
    ax = axes[1, 0]
    tp = np.linspace(0.0, 120.0, 300)
    envp = permeant_envelope(tp, CONDITIONS["strong_hyper"])
    _band(ax, envp, "C2", "permeant r=1.30 (b,τ_s/τ_w)")
    ax.axhline(1.0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("time  t/τ_w  [–]")
    ax.set_ylabel("v(t) = V/V₀  [–]")
    ax.set_title("Permeant solute σ(t): shrink-then-re-swell\n(glycerol) vs sustained (mannitol)")
    ax.legend(fontsize=8)

    # (4) stiffness sensitivity band vs volume ratio --------------------------
    ax = axes[1, 1]
    v = np.linspace(0.6, 1.4, 100)
    envs = stiffness_envelope(v)
    for row, lab in zip(envs.curves, envs.labels):
        ax.plot(v, row, lw=1.2, label=lab)
    ax.fill_between(v, envs.lo, envs.hi, color="C4", alpha=0.15)
    ax.axhline(1.0, color="k", lw=0.6, ls=":")
    ax.axvline(1.0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("volume ratio  v = V/V₀  [–]")
    ax.set_ylabel("E/E₀  [–]")
    ax.set_title("Apparent-stiffness sensitivity E/E₀ = v^(−m)\n(m ∈ {0,1,2}; Guo 2017 ≈ v^−2)")
    ax.legend(fontsize=8)

    fig.suptitle("Osmotic / volume-response prediction-envelope ORACLE "
                 "(dimensionless, literature-range) — closed-form, not runtime",
                 fontsize=12, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    path = OUT / "osmotic_envelope.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"wrote {path}")
    return path


if __name__ == "__main__":
    main()
