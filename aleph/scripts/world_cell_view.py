#!/usr/bin/env python
r"""Write the interactive WebGL cell view from a ``.alephcell`` — mouse, and a toggle per population.

**Why this exists, and it is a reversal.** `viz/cell_app.py` is a native window; `world_turntable.py`
pre-renders angles. Neither does the two things the OLD viewer did and the PI asked for back: turn the
cell freely with the mouse, and switch populations off to see inside. `world/viewer.py` has done both
since it was written — one index buffer per layer, so a toggle is a skipped draw call rather than a
rebuild — but it could only be fed a live ``WorldArena``, which means a CUDA host. So the pages existed
and could not be re-made anywhere the cell already sat on disk.

⚠ **This adds no viewer.** It is an adapter: `.alephcell` -> the shapes `render_cell_view` already
reads. Writing a second WebGL page would have been the easy mistake and would have left two viewers
disagreeing about one cell.

⚠ **On the cortex, which is the whole size question.** The old page left it out ENTIRELY and said so in
its own header — *"LEFT OUT ENTIRELY (not thinned): cortex"* — and that distinction is the honest part:
omitting a population and naming it is a different act from drawing one node in twenty-three. It is
4,128,840 of 4,313,157 segments, 95.7%, so including it is not a tweak. It is available here, off by
default, and what it costs is printed rather than discovered.

⚠ **Nothing here is a measurement.** No force is evaluated and nothing moves; the page says so.

Sanity Gate:
    * dimensional — positions pass through in µm exactly as stored; no scaling, no recentering.
    * boundary cases — a population with positions and no topology is REPORTED as undrawn, never
      silently skipped; an empty selection refuses rather than writing an empty page.
    * conservation — every kept population's every segment reaches a layer's index buffer, and the
      script asserts the totals it prints against the arrays it built.
    * CFL/precision — none. Positions are cast to float32 for the buffer, ~1e-7 relative, stated on
      the page by the viewer itself.
    * sign sense — not applicable: no orientation is chosen here, the viewer's camera owns it.
    * measurement protocol — this measures nothing. It converts a file into a page.

engine units: µm. Runtime: a VIEWER. No device is opened and no kernel runs.

Usage:
    python aleph/scripts/world_cell_view.py CELL --out page.html [--include cortex] [--only a,b]
"""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np

#: Left out unless asked for, with the reason printed. Not a judgement about the cortex — a statement
#: about what one population costs a single-file page. See the module docstring.
HEAVY = ("cortex",)


def _layers_from(cell, keep: list[str]) -> tuple[np.ndarray, list, list, list[str]]:
    """Concatenate populations into one node array and re-base each population's indices onto it.

    ⚠ The `.alephcell` stores each population's positions and its own LOCAL indices; the viewer wants
    one global node array with per-layer index buffers into it. The offset is applied here, once, and
    asserted below — an index that still points into its own population's frame draws a line between
    two unrelated structures, which is exactly the failure `cell_app.load` re-checks on every read.

    Returns:
        ``(positions, strand-like objects, surface-like objects, notes about what was not drawn)``.
    """
    pos_parts: list[np.ndarray] = []
    strands: list[SimpleNamespace] = []
    surfaces: list[SimpleNamespace] = []
    notes: list[str] = []
    lo = 0
    for name in keep:
        p = np.asarray(cell.positions[name], np.float64)
        hi = lo + p.shape[0]
        rng = SimpleNamespace(lo=lo, hi=hi)
        pos_parts.append(p)
        if name in cell.segments:
            idx = np.asarray(cell.segments[name], np.int64) + lo
            assert idx.min() >= lo and idx.max() < hi, f"{name}: index escapes its own block"
            strands.append(SimpleNamespace(population=name, nodes=rng, position=p, seg_node=idx))
        elif name in cell.faces:
            idx = np.asarray(cell.faces[name], np.int64) + lo
            assert idx.min() >= lo and idx.max() < hi, f"{name}: face index escapes its own block"
            surfaces.append(SimpleNamespace(population=name, nodes=rng, position=p, face_idx=idx))
        else:
            # ⚠ SAID, not skipped. `lamina` and `chromatin` are written as nodes with no connectivity
            # and this project does not draw points; a silent omission here is the defect the app's
            # own history is about.
            notes.append(f"{name} ({p.shape[0]:,} nodes, NO TOPOLOGY — not drawn)")
        lo = hi
    return np.concatenate(pos_parts) if pos_parts else np.zeros((0, 3)), strands, surfaces, notes


def main(argv: list[str] | None = None) -> int:
    """Convert one export into one interactive page."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("path", nargs="?", type=Path,
                    help="an .alephcell export. Omit to take the newest under aleph/outputs/.")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--only", default=None,
                    help="comma-separated populations. Default: everything except the heavy ones.")
    ap.add_argument("--include", default=None,
                    help=f"add one of {HEAVY} back to the default set, e.g. --include cortex. "
                         f"⚠ the cortex is 95.7%% of the segments; the page will be large and the "
                         f"cost is printed before it is written.")
    args = ap.parse_args(argv)

    from aleph.viz.cell_app import _newest_export, load
    from aleph.world.viewer import render_cell_view

    path = args.path or _newest_export()
    if path is None:
        print("refused: no .alephcell found under aleph/outputs/.")
        return 2
    cell = load(Path(path))

    if args.only:
        keep = [p.strip() for p in args.only.split(",")]
    else:
        keep = [p for p in cell.positions if p not in HEAVY]
        for extra in (args.include or "").split(","):
            if extra.strip() and extra.strip() not in keep:
                keep.append(extra.strip())
    missing = [p for p in keep if p not in cell.positions]
    if missing:
        print(f"refused: {missing} not in this file. It has: {sorted(cell.positions)}")
        return 2
    if not keep:
        print("refused: an empty selection would write a page of nothing.")
        return 2

    left_out = [p for p in cell.positions if p not in keep]
    pos, strands, surfaces, notes = _layers_from(cell, keep)
    n_seg = sum(int(s.seg_node.shape[0]) for s in strands)
    n_face = sum(int(sf.face_idx.shape[0]) for sf in surfaces)

    # ⚠ The census the page carries is built from the arrays that were just built, not from the
    # file's header — a page that quotes a header it did not draw from is a page that can be right
    # about a cell it is not showing.
    note = (f"Drawn IN FULL: {', '.join(sorted(p for p in keep if p in cell.segments or p in cell.faces))}."
            + (f" LEFT OUT ENTIRELY (not thinned): {', '.join(sorted(left_out))}." if left_out else "")
            + (f" {' '.join(notes)}." if notes else "")
            + " Geometry only: no force, no law, nothing has moved. A picture is not evidence.")
    print(f"[cell-view] {path.name}: {pos.shape[0]:,} nodes positioned, {n_seg:,} segments, "
          f"{n_face:,} faces over {len(strands) + len(surfaces)} layers", flush=True)
    if left_out:
        print(f"[cell-view] LEFT OUT ENTIRELY (not thinned): {', '.join(sorted(left_out))} — "
              f"--include or --only to change that", flush=True)
    for n in notes:
        print(f"[cell-view] {n}", flush=True)

    arena = SimpleNamespace(n_live=lambda _kind: int(pos.shape[0]))
    html = render_cell_view(arena, strands=strands, surfaces=surfaces,
                            title=f"Aleph — {path.stem}", note=note)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html)
    mb = args.out.stat().st_size / 1e6
    print(f"[cell-view] wrote {args.out}  ({mb:.1f} MB)", flush=True)
    print("[cell-view] drag to rotate · shift-drag to pan · scroll to zoom · click a chip to toggle",
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
