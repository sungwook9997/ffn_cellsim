"""Interactive 3D viewer — the FF cortex filament network at varying filament count (NF).

Shows what "increasing / decreasing the cortex filaments" actually looks like: the woven
actin cortex rendered as full-resolution line segments (NO downsampling, PI rule), one SCENE
per NF selectable from the dropdown. Self-contained double-clickable HTML (three.js viewer;
opens in a browser, not a CSP-locked artifact).

Usage: python -m aleph.scripts.ff_cortex_density_viewer [--nfs 500,2000,8000,38000,70686]
"""

from __future__ import annotations

import argparse

import numpy as np

from aleph.laws.gamma_floor import CortexParams, build_crosslinked_cortex
from aleph.scripts.ff_viewer_html import build_viewer


def cortex_segments(n_filaments: int, seed: int = 1):
    """Build a fresh cortex at ``n_filaments`` and return (segment-endpoint verts, n_segments)."""
    cx = build_crosslinked_cortex(
        CortexParams(), n_filaments=n_filaments, n_xl=n_filaments,
        n_myo=max(1, n_filaments // 10), rng=np.random.default_rng(seed),
    )
    pos = cx.net.pos
    seg = cx.net.segments
    return pos[seg].reshape(-1, 3), int(seg.shape[0])


def main() -> None:
    ap = argparse.ArgumentParser(description="FF cortex filament-density 3D viewer (per NF).")
    ap.add_argument("--nfs", default="500,2000,8000,38000,70686",
                    help="comma-separated filament counts (sparse → native)")
    ap.add_argument("--out", default="aleph/outputs/mech_hier/figs/cortex_density_3d.html")
    a = ap.parse_args()

    nfs = [int(x) for x in a.nfs.split(",")]
    scenes = {}
    for nf in nfs:
        verts, nseg = cortex_segments(nf)
        tag = "native" if nf >= 38000 else "coarse"
        scenes[f"NF={nf} · {nseg} segments · {tag}"] = [{
            "name": f"cortex filaments (NF={nf}, FULL-RES)",
            "kind": "lines", "verts": verts, "color": "#8fbff0", "size": 1.2, "opacity": 0.55,
        }]
        print(f"NF={nf:6d}: {nseg} segments", flush=True)

    build_viewer(scenes, out=a.out, title="FF cortex — filament density (NF sweep)")
    print(f"wrote {a.out}", flush=True)


if __name__ == "__main__":
    main()
