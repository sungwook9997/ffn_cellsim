"""Interactive 3-D render of the staged Active-Cell (ac/engine) build, per COMPONENT + per CONNECTOR.

The staged ``ac/engine`` build is a component/connector graph: components own a DISJOINT node population;
connectors are the ONLY mechanical joints (co-location in an array is NEVER a connection). This viewer makes
that build visually readable — one isolation scene per component, one coloured load-path scene per connector
FAMILY (the explicit joints, drawn so the force PATH is visible, not mere proximity), a CUMULATIVE composition
family that grows one component at a time as slices bind, and the composed whole cell coloured by the composed
per-node force |F|.

It is fully **data-driven**: the scene set is assembled from a state **dump** with NO per-component hardcoding
(see :mod:`ac_viz_common`). Two dump schemas are accepted and normalised to one model:

* **v2 generic** — written by ``scripts/dump_state.py`` from the ``ac/engine`` graph (or the incumbent
  ``build_cell`` re-expressed as components/connectors) on the gbook A5000. The moment a new component (SF, MT,
  IF, ECM, lamellipodium, filopodium) or connector family (ERM, FA-clutch, LINC, MOTOR, plectin, spectraplakin,
  transient-actin, immersed-transfer) appears in the dump, it gets its scenes automatically.
* **legacy** — the incumbent ``assembled_state.npz`` (cortex F-actin + NMII + deformable nucleus + membrane +
  LINC), re-expressed as the same component/connector model so the generic viewer works on the real dump that
  already exists on disk.

Rendering is numpy/matplotlib only (dev Mac — NO Warp/CUDA). Colour = composed per-node |F| on a perceptually
ordered **turbo** ramp, **log₁₀** axis (the field spans >4 decades: resting ~0.02 pN vs excluded-volume
hotspots 10²–10³ pN — a linear scale would show only the hotspots; log is annotated on every colour bar).
FULL native resolution — NO downsampling.

    PYTHONPATH=/Users/sw1/ffn_cellsim python aleph/scripts/ac_cell_assembled_viz.py [--npz PATH] [--out PATH]
        -> aleph/outputs/ac/cell_assembled/ac_cell_assembled.html  (+ browser-verified screenshots)
"""
from __future__ import annotations

import argparse
from pathlib import Path

from aleph.scripts.ac_viz_common import build_scenes, compose_title, load_dump
from aleph.scripts.ff_viewer_html import build_viewer

_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_NPZ = _ROOT / "outputs" / "ac" / "cell_assembled" / "assembled_state.npz"
DEFAULT_OUT = _ROOT / "outputs" / "ac" / "cell_assembled" / "ac_cell_assembled.html"


def build(npz: Path, out_html: Path) -> tuple[str, dict]:
    """Read a state dump (v2 or legacy) and emit the data-driven per-component/connector interactive HTML."""
    dump = load_dump(npz)
    scenes, cbars = build_scenes(dump)
    title = compose_title(dump)
    build_viewer(scenes, out=str(out_html), title=title, cbars=cbars)
    counts = {
        "n_components": len(dump.components),
        "n_connectors": len(dump.connectors),
        "n_scenes": len(scenes),
        **{f"comp:{c.key}_nodes": c.n_nodes for c in dump.components},
        **{f"comp:{c.key}_actors": c.n_actors for c in dump.components},
        **{f"conn:{cn.family}_{cn.a_component}-{cn.b_component}": cn.n for cn in dump.connectors},
    }
    return str(out_html), counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--npz", type=Path, default=DEFAULT_NPZ, help="state dump (.npz): v2 generic or legacy")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output interactive HTML")
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    path, counts = build(a.npz, a.out)
    print(f"wrote {path}  (from {a.npz})")
    for k, v in counts.items():
        print(f"  {k:34s} {v:,}" if isinstance(v, int) else f"  {k:34s} {v}")


if __name__ == "__main__":
    main()
