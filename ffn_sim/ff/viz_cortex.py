"""Visualize the FF cortex fiber network on the sphere (Stage 6c-b sanity figure).

Renders the assembled :class:`~ff.cortex_assembly.FiberNetwork` as 3D fiber polylines on the shell
+ a segment-length histogram (the inextensibility quantity) so the geometry can be eyeballed.
Writes to ``ffn_sim/outputs/ff/figs/``. Run: ``python -m ffn_sim.ff.viz_cortex``.
"""

from __future__ import annotations

import os

import numpy as np

from ffn_sim.ff.constraints import segment_lengths
from ffn_sim.ff.cortex_assembly import CortexParams, build_cortex_network

OUTDIR = os.path.join(os.path.dirname(__file__), "..", "outputs", "ff", "figs")


def render(n_filaments: int = 200, seed: int = 0, outdir: str = OUTDIR) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    params = CortexParams()
    net, meta = build_cortex_network(params, rng=np.random.default_rng(seed),
                                     n_filaments=n_filaments)
    os.makedirs(outdir, exist_ok=True)

    fig = plt.figure(figsize=(13, 6))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    off = net.fiber_offsets
    for f in range(net.n_fibers):
        sl = slice(int(off[f]), int(off[f + 1]))
        seg = net.pos[sl]
        ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], lw=0.6, color="tab:blue", alpha=0.6)
    ax.set_title(f"FF cortex: {n_filaments} actin fibers on R={params.R_um:.0f}µm shell\n"
                 f"(L={params.L_filament_um:.1f}µm, {params.beads_per_filament} beads, "
                 f"ℓ₀={params.seg_um:.1f}µm, κ={params.kappa:.3f} pN·µm²)")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("y [µm]"); ax.set_zlabel("z [µm]")
    R = params.R_um
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_zlim(-R, R)

    ax2 = fig.add_subplot(1, 2, 2)
    L = segment_lengths(net) * 1000.0   # nm
    ax2.hist(L, bins=40, color="tab:green", alpha=0.8)
    ax2.axvline(params.seg_um * 1000.0, color="k", ls="--",
                label=f"ℓ₀ = {params.seg_um * 1000:.0f} nm")
    ax2.set_xlabel("segment length [nm]"); ax2.set_ylabel("count")
    full_density = params.areal_density_um2
    ax2.set_title(f"segment lengths ({n_filaments} rendered; full cortex "
                  f"{params.n_filaments} fib → {full_density:.2f} µm⁻²)\n"
                  f"on-shell residual {meta['on_shell_residual_um'] * 1000:.2e} nm")
    ax2.legend()

    fig.tight_layout()
    path = os.path.join(outdir, "cortex_assembly.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"wrote {path}")
    print(f"meta: {meta}")
    return path


if __name__ == "__main__":
    render()
