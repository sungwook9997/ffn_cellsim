"""One command from the moving-boundary kernels to a figure of what an interior gate cannot see.

`ALEPH-PORT-3617`.

    python -m aleph.viz.boundary_figure

No arguments a person has to invent. It launches the kernels on CPU-Warp, grades them against the
frozen host law, and writes an SVG. A refusal exits **2** and says what is missing; exit **1** is an
unexpected fault, so a caller can tell a refusal from a crash — the convention `gpu_preflight.py`
established and `evidence_figure.py` carried, and the reason carries over exactly: a tool that
refuses and exits 0 lets ``make-figure && publish`` treat "nothing was drawn" as success.

What the figure shows that reading the numbers does not
-------------------------------------------------------
**The blind spot has a shape, and the shape is the argument.** Panel D is `|true - centre_ghost|` on
the same slice as B and C. It is **exactly zero** across the interior and lights up in a one-cell
shell against the mask. A table saying "interior 0.0, boundary 11.4" is the same measurement; the
picture is what makes it obvious that an interior-only gate is not merely *unlucky* — it is
sampling the one region where the two operators agree by construction.

**The mask is drawn, so the reader can check the claim.** Every panel carries the fluid extent as a
stroked outline. The defect lives exactly one cell inside that outline, which is visible rather than
asserted.

**Node marks are sized by how hard the mask bites them.** `W_n = 1` is the unclipped case — the one
on which `WRONG_UNMASKED_NORMALISATION` is bit-identical to the truth and the gate is blind. A state
that admitted such a node would show a full-size ring, and there are none, which is the state
builder's rejection rule made visible instead of trusted.

**Exact zero is drawn as absence, not as a small colour.** Panel D uses white for `0.0` and starts
its ramp above it. A diverging map with a near-zero midpoint would render round-off and structural
zero identically, which is the distinction the whole lane turns on.

Units: `1/s` for divergence, dimensionless for `W_n` and for the stencil weights.
"""

from __future__ import annotations

import argparse
import html
import sys
from pathlib import Path

import numpy as np

__all__ = [
    "EXIT_OK",
    "EXIT_REFUSED",
    "BoundaryEvidence",
    "collect",
    "main",
    "render_boundary_svg",
]

EXIT_OK = 0
EXIT_FAULT = 1
EXIT_REFUSED = 2

#: The slice drawn. The fluid block is `[2, 6]` on each axis, so `4` is its middle plane — the one
#: where the boundary shell is a closed ring rather than a corner.
SLICE_K = 4

CELL_PX = 30.0
PANEL_GAP = 46.0
MARGIN = 28.0
HEADER = 92.0
#: Panel title, subtitle, and the gap above the cells. Reserved per panel ROW, not once.
PANEL_HEAD = 42.0
#: Per caption entry: the measurement line, its gloss, and the space to the next entry.
CAPTION_STEP = 38.0


class BoundaryFigureRefused(RuntimeError):
    """Raised when the figure cannot be a picture of a measurement. Exits 2, never 0."""


class BoundaryEvidence:
    """Everything the figure draws, measured once. No panel re-derives a number it renders."""

    def __init__(self, *, grid, nodes, node_weight, raw_weight, true_div, ghost_div, spacing):
        self.grid = grid
        self.nodes = nodes
        self.node_weight = node_weight
        self.raw_weight = raw_weight
        self.true_div = true_div
        self.ghost_div = ghost_div
        self.spacing = spacing

        self.fluid = np.asarray(grid.fluid_mask, dtype=bool)
        self.gap = np.abs(true_div - ghost_div)

        # "Interior" means every one of the six axis neighbours is fluid — the region an
        # interior-only gate restricts to, defined from the mask rather than from a hand-picked box.
        interior = self.fluid.copy()
        for axis in range(3):
            for shift in (1, -1):
                interior &= np.roll(self.fluid, shift, axis=axis)
        self.interior = interior
        self.shell = self.fluid & ~interior

        self.interior_gap = float(np.max(self.gap[interior])) if interior.any() else float("nan")
        self.shell_gap = float(np.max(self.gap[self.shell])) if self.shell.any() else float("nan")
        self.unclipped = int(np.sum(node_weight >= raw_weight - 1e-12))

    def faults(self) -> list[str]:
        """What would make this figure a picture of something other than the measurement."""
        out: list[str] = []
        if not self.interior.any():
            out.append("no cell has six fluid neighbours, so there is no interior to compare")
        if not self.shell.any():
            out.append("the fluid region has no boundary shell, so the defect has nowhere to live")
        if self.interior_gap != 0.0:
            out.append(
                f"the two divergence forms differ by {self.interior_gap:.6e} in the interior. "
                "The blind spot has moved and panel D no longer shows what it claims"
            )
        if self.shell_gap <= 1.0:
            out.append(
                f"the boundary gap is only {self.shell_gap:.6e}; a defect that small is not "
                "evidence that an interior-only gate missed anything"
            )
        if self.unclipped:
            out.append(
                f"{self.unclipped} node(s) have an unclipped stencil, on which the unmasked "
                "normalisation is bit-identical to the truth -- the state cannot grade it"
            )
        return out


def collect(*, shape=(9, 9, 9), spacing_um: float = 0.25, fluid=slice(2, 7), seed: int = 3181):
    """Launch the kernels and grade them against the frozen law. Refuses rather than drawing a
    figure of a run that did not happen."""
    try:
        from aleph.runtime import law_cases, law_kernels, moving_boundary as mb
        from aleph.runtime.backend import WarpBackend
    except ImportError as exc:  # pragma: no cover - exercised by the import-boundary test
        raise BoundaryFigureRefused(f"the runtime is not importable: {exc}") from exc

    try:
        kernels = law_kernels.load()
    except Exception as exc:
        raise BoundaryFigureRefused(
            f"NVIDIA Warp did not load, so there are no kernels to draw: {exc}"
        ) from exc

    wp = law_kernels.loaded_warp()
    # CPU device: unguarded by construction. This figure is about a stencil, not about hardware, and
    # the GPU grant lapsed at 2026-08-04T00:00. NOTHING HERE WRITES AN AUTHORIZATION RECORD.
    backend = WarpBackend(device="cpu", dtype=np.float64)
    device = backend.warp_device

    state = law_cases.moving_boundary_state(shape=shape, spacing_um=spacing_um, fluid=fluid)
    grid = law_cases._mb_grid(state)
    nodes = np.asarray(state["nodes_um"], dtype=np.float64)
    count = nodes.shape[0]
    cells = grid.cell_count
    mask_flat = np.asarray(grid.fluid_mask, dtype=bool).reshape(-1)

    def as_int32(a):
        return wp.array(np.ascontiguousarray(a, dtype=np.int32), dtype=wp.int32, device=device)

    d_cell = as_int32(np.zeros((count, 64)))
    d_weight = backend.zeros((count, 64), dtype=np.float64)
    wp.launch(
        kernels["mb_block"],
        dim=(count, 64),
        inputs=[
            backend.array(nodes, dtype=np.float64),
            backend.array(state["grid_origin_um"], dtype=np.float64),
            wp.float64(float(spacing_um)),
            as_int32(shape),
            d_cell,
            d_weight,
        ],
        device=device,
    )
    d_node_weight = backend.zeros(count, dtype=np.float64)
    wp.launch(
        kernels["mb_node_weight_true"],
        dim=count,
        inputs=[d_cell, d_weight, as_int32(mask_flat), d_node_weight],
        device=device,
    )

    rng = np.random.default_rng(seed)
    velocity = rng.standard_normal((*shape, 3))
    d_velocity = backend.array(velocity.reshape(cells, 3), dtype=np.float64)

    def divergence(name: str) -> np.ndarray:
        out = backend.zeros(cells, dtype=np.float64)
        wp.launch(
            kernels[name],
            dim=cells,
            inputs=[d_velocity, as_int32(mask_flat), as_int32(shape), wp.float64(spacing_um), out],
            device=device,
        )
        backend.synchronize()
        return backend.to_host(out).reshape(shape)

    backend.synchronize()
    true_div = divergence("mb_divergence_true")
    ghost_div = divergence("mb_divergence_wrong_centre_ghost")

    # The kernels are graded here rather than trusted. A figure drawn from a kernel that does not
    # reproduce the frozen law is a picture of a bug wearing the law's name.
    host_true = mb.boundary_aware_divergence(velocity, grid)
    drift = float(np.max(np.abs(true_div - host_true)))
    if drift != 0.0:
        raise BoundaryFigureRefused(
            f"the divergence kernel is {drift:.6e} from the frozen host law; it must be "
            "bit-identical, and a figure drawn from it would not be a picture of the law"
        )

    return BoundaryEvidence(
        grid=grid,
        nodes=nodes,
        node_weight=backend.to_host(d_node_weight),
        raw_weight=np.sum(backend.to_host(d_weight), axis=1),
        true_div=true_div,
        ghost_div=ghost_div,
        spacing=float(spacing_um),
    )


def _esc(text: object) -> str:
    return html.escape(str(text), quote=True)


def _diverging(value: float, scale: float) -> str:
    """Blue-white-red about zero. Used for the two divergence fields, which are signed."""
    if scale <= 0.0:
        return "#ffffff"
    t = max(-1.0, min(1.0, value / scale))
    if t >= 0.0:
        return f"rgb(255,{int(255 * (1 - t) + 0.5)},{int(255 * (1 - t) + 0.5)})"
    return f"rgb({int(255 * (1 + t) + 0.5)},{int(255 * (1 + t) + 0.5)},255)"


def _magnitude(value: float, scale: float) -> str:
    """White at EXACTLY zero, then a ramp that starts visibly above it.

    The gap between structural zero and round-off is the distinction this whole lane turns on, so
    the colour map refuses to render them alike: `0.0` is the paper, and anything nonzero starts at
    a colour a reader can see.
    """
    if value == 0.0 or scale <= 0.0:
        return "#ffffff"
    t = 0.18 + 0.82 * min(1.0, value / scale)
    return f"rgb({int(255 * t + 0.5)},{int(255 * (1 - 0.75 * t) + 0.5)},{int(255 * (1 - t) + 0.5)})"


def _panel(evidence, field, x0, y0, *, title, subtitle, colour, scale, nodes=False):
    n0, n1, _ = evidence.grid.shape
    out = [
        f'<text x="{x0:.1f}" y="{y0 - 26:.1f}" class="ptitle">{_esc(title)}</text>',
        f'<text x="{x0:.1f}" y="{y0 - 10:.1f}" class="psub">{_esc(subtitle)}</text>',
    ]
    for i in range(n0):
        for j in range(n1):
            fill = colour(float(field[i, j, SLICE_K]), scale)
            out.append(
                f'<rect x="{x0 + j * CELL_PX:.1f}" y="{y0 + i * CELL_PX:.1f}" '
                f'width="{CELL_PX:.1f}" height="{CELL_PX:.1f}" fill="{fill}" '
                'stroke="#eceff1" stroke-width="0.5"/>'
            )
    # The fluid extent, stroked. Every claim in the caption is about position relative to this.
    fluid_slice = evidence.fluid[:, :, SLICE_K]
    rows = np.flatnonzero(fluid_slice.any(axis=1))
    cols = np.flatnonzero(fluid_slice.any(axis=0))
    if rows.size and cols.size:
        out.append(
            f'<rect x="{x0 + cols[0] * CELL_PX:.1f}" y="{y0 + rows[0] * CELL_PX:.1f}" '
            f'width="{(cols[-1] - cols[0] + 1) * CELL_PX:.1f}" '
            f'height="{(rows[-1] - rows[0] + 1) * CELL_PX:.1f}" fill="none" '
            'stroke="#37474f" stroke-width="1.6" stroke-dasharray="4 2"/>'
        )
    if nodes:
        for n in range(evidence.nodes.shape[0]):
            pos = evidence.nodes[n] / evidence.spacing
            if abs(pos[2] - SLICE_K) > 1.0:
                continue
            ratio = float(evidence.node_weight[n] / max(evidence.raw_weight[n], 1e-30))
            out.append(
                f'<circle cx="{x0 + pos[1] * CELL_PX:.1f}" cy="{y0 + pos[0] * CELL_PX:.1f}" '
                f'r="{2.5 + 6.0 * ratio:.2f}" fill="none" stroke="#00695c" stroke-width="1.6"/>'
            )
    return out


def _caption_lines(evidence: BoundaryEvidence) -> list[tuple[str, str]]:
    """`(measurement, gloss)` per line. One source, read by both the width computation and the
    drawing — sizing the canvas from text that is not the text drawn is how a caption ends up
    off the right edge while every number in it is correct."""
    return [
        (
            f"interior   max |B - C| = {evidence.interior_gap:.1f}"
            f"   over {int(evidence.interior.sum())} cells (every axis neighbour fluid)",
            "the two operators agree BY CONSTRUCTION here - an interior-only gate is not unlucky,"
            " it samples the one region where the defect cannot appear",
        ),
        (
            f"boundary   max |B - C| = {evidence.shell_gap:.4f}"
            f"   over {int(evidence.shell.sum())} cells (the one-cell shell)",
            "the numerator is right and only the measure is wrong, so the defect is a clean factor"
            " of two rather than noise",
        ),
        (
            f"unclipped nodes   {evidence.unclipped} of {evidence.nodes.shape[0]}",
            "on an unclipped node the unmasked normalisation is bit-identical to the truth and"
            " cannot be graded at all",
        ),
    ]


def render_boundary_svg(evidence: BoundaryEvidence, *, title: str | None = None) -> str:
    """The four panels, the mask, and a caption that states what was measured and where."""
    n0, n1, _ = evidence.grid.shape
    panel_w = n1 * CELL_PX
    panel_h = n0 * CELL_PX
    # TWO BY TWO, not one by four. A single row of four made the canvas 2.4:1, and every renderer
    # that fits a wide figure to a square viewport dropped the fourth panel off the right edge --
    # which is the panel carrying the finding. A layout whose most important element is the one most
    # likely to be cropped is a bad layout however correct its numbers are.
    caption = 30.0 + CAPTION_STEP * len(_caption_lines(evidence))
    panel_row = MARGIN * 2 + panel_w * 2 + PANEL_GAP
    # The gloss is 11px system-ui; ~5.6px per character is that face's mean advance. The canvas is
    # sized from the WIDER of the two, because text drawn past the right edge is invisible in every
    # renderer and no test that checks the numbers would notice.
    longest = max(len(gloss) for _, gloss in _caption_lines(evidence))
    width = max(panel_row, MARGIN * 2 + 18 + longest * 5.6)
    height = HEADER + (panel_h + PANEL_HEAD) * 2 + caption

    div_scale = float(np.max(np.abs(evidence.true_div))) or 1.0
    gap_scale = float(np.max(evidence.gap)) or 1.0

    body: list[str] = [
        f'<rect width="{width:.1f}" height="{height:.1f}" fill="#ffffff"/>',
        f'<text x="{MARGIN:.1f}" y="34" class="title">'
        f"{_esc(title or 'The moving-boundary element: what an interior-only gate cannot see')}"
        "</text>",
        f'<text x="{MARGIN:.1f}" y="55" class="sub">ALEPH-PORT-3617 · '
        f"k = {SLICE_K} slice · kernels on CPU-Warp, bit-identical to the frozen host law"
        "</text>",
    ]

    panels = [
        (
            evidence.fluid.astype(float),
            "A · the mask and the nodes",
            "ring size = W_n / sum Phi",
            lambda v, s: "#e0f2f1" if v > 0.5 else "#fafafa",
            1.0,
            True,
        ),
        (
            evidence.true_div,
            "B · div(V), boundary-aware",
            "one-sided over h at the edge",
            _diverging,
            div_scale,
            False,
        ),
        (
            evidence.ghost_div,
            "C · div(V), inherited centre-ghost",
            "centre as ghost, still over 2h",
            _diverging,
            div_scale,
            False,
        ),
        (
            evidence.gap,
            "D · |B - C|",
            "white is EXACTLY zero",
            _magnitude,
            gap_scale,
            False,
        ),
    ]
    for index, (field, ptitle, psub, colour, scale, nodes) in enumerate(panels):
        x0 = MARGIN + (index % 2) * (panel_w + PANEL_GAP)
        y0 = HEADER + (index // 2) * (panel_h + PANEL_HEAD) + PANEL_HEAD
        body.extend(
            _panel(
                evidence,
                field,
                x0,
                y0,
                title=ptitle,
                subtitle=psub,
                colour=colour,
                scale=scale,
                nodes=nodes,
            )
        )

    y = HEADER + (panel_h + PANEL_HEAD) * 2 + 30
    for measurement, gloss in _caption_lines(evidence):
        body.append(f'<text x="{MARGIN:.1f}" y="{y:.1f}" class="cap">{_esc(measurement)}</text>')
        body.append(
            f'<text x="{MARGIN + 18:.1f}" y="{y + 15:.1f}" class="gloss">{_esc(gloss)}</text>'
        )
        y += 38

    style = (
        "<style>"
        ".title{font:600 17px system-ui,sans-serif;fill:#263238}"
        ".sub{font:12px system-ui,sans-serif;fill:#607d8b}"
        ".ptitle{font:600 12px system-ui,sans-serif;fill:#37474f}"
        ".psub{font:10px system-ui,sans-serif;fill:#78909c}"
        ".cap{font:11px ui-monospace,SFMono-Regular,monospace;fill:#263238}"
        ".gloss{font:11px system-ui,sans-serif;fill:#78909c}"
        "</style>"
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}">{style}' + "".join(body) + "</svg>"
    )


def main(argv: "list[str] | None" = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("docs/results/2026-08-04-moving-boundary/interior_blind_spot.svg"),
        help="where to write the SVG",
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    try:
        evidence = collect()
    except BoundaryFigureRefused as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return EXIT_REFUSED

    faults = evidence.faults()
    if faults:
        for fault in faults:
            print(f"REFUSED: {fault}", file=sys.stderr)
        return EXIT_REFUSED

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(render_boundary_svg(evidence), encoding="utf-8")

    if not args.quiet:
        print(f"wrote {args.out}")
        print(
            f"  interior  max |true - ghost| = {evidence.interior_gap:.1f}"
            f" over {int(evidence.interior.sum())} cells"
        )
        print(
            f"  boundary  max |true - ghost| = {evidence.shell_gap:.4f}"
            f" over {int(evidence.shell.sum())} cells"
        )
        print(
            f"  {evidence.unclipped} of {evidence.nodes.shape[0]} nodes unclipped"
            " — an unclipped node cannot grade the unmasked normalisation"
        )
    return EXIT_OK


if __name__ == "__main__":  # pragma: no cover - exercised through `main` in the tests
    raise SystemExit(main())
