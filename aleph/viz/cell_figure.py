"""A composed cell as the thing it is: surfaces, curves and bonds — not a point per node.

:mod:`aleph.viz.render` emits one ``<circle>`` per node and says so at its line 28: *"Individual
elements, not a shaded surface."* That decision throws the topology away, and the topology **is** the
physics. The owners publish far more than positions —

    membrane   162 positions, 320 TRIANGLES
    cortex     486 positions, 324 SEGMENTS, 162 TRIPLES, 480 CROSSLINKS

— the axial law lives on segments, the bending law on triples, and the crosslinks are what connect
the network at all. **A stretched segment and a slack one are the same two dots.** On 2026-08-05 it
took a Hessian eigendecomposition to find that every cortex crosslink attached at one material
coordinate; drawn as bonds, the crosslink graph *is* the visible shell and the finding is in the
picture.

Four rules this module enforces rather than offers, each because its absence hid something real.

**A ground plane.** The cell sat 1.34 µm above its substrate and no render showed it until one was
drawn in side elevation. A visible ground costs one quad.

**A scale bar in µm.** Taken from the reference viewer's own source comment: *"the viewer is in
physical µm, so a zoom-reactive ruler keeps every screenshot self-describing — visualization rule:
always annotate units."*

**A logarithmic colour scale.** Forces here span four decades — one adhesion at ~10³ pN against a
cortex node at ~10⁻² — so a linear map puts every filament in the first colour and the figure reports
one outlier.

**Orthographic projection.** Lengths on the page are proportional to lengths in the cell, so the
figure can be measured. Perspective would make the scale bar a lie away from its own depth.

Layers are toggled by **class, not by group**. The depth order is a global painter's-algorithm sort,
and putting each layer in its own ``<g>`` would reorder the document and destroy it; every element
instead carries ``class="L-<layer>"`` and :func:`layer_toggle_page` hides them with CSS.

**This is a projection, not evidence.** It shows what the arrays contain. Whether the arrays are
right is what the controls are for.
"""

from __future__ import annotations

import base64
import json
from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = [
    "FORCE_DECADES",
    "CellFigure",
    "CellScene",
    "interactive_page",
    "layer_toggle_page",
    "log_norm",
    "render_cell_scene",
    "render_cell_svg",
]

#: How many decades the force colour ramp spans below its maximum. Four covers the measured spread
#: (adhesions ~10³ pN, cortex nodes ~10⁻² pN) without giving round-off its own colour band.
FORCE_DECADES = 4.0

#: Filament owners: colour, stroke width, label. Drawn along the owner's own ``segments``.
FILAMENT_STYLE: dict[str, tuple[str, float, str]] = {
    "cortex": ("#8FB8D6", 1.1, "cortex F-actin"),
    "sf_arc": ("#E85D4E", 2.8, "stress fibres"),
    "microtubule": ("#3FA88F", 2.6, "microtubule"),
    "intermediate_filament": ("#8898AC", 2.3, "interm. filament"),
    "lamellipodium": ("#F0B23C", 2.8, "lamellipodium"),
    "filopodium": ("#E8873A", 2.8, "filopodium"),
    "ecm": ("#B08A3E", 1.6, "ECM fibres"),
    "nmii": ("#D2568C", 2.5, "NMII"),
}

#: Closed surfaces: base colour, opacity, label. The membrane is drawn **far side only** so the cell
#: is a body you can see into rather than an opaque ball with everything hidden behind it.
SURFACE_STYLE: dict[str, tuple[tuple[float, float, float], float, str]] = {
    "membrane": ((0.34, 0.56, 0.76), 0.26, "plasma membrane"),
    "nucleus": ((0.55, 0.38, 0.66), 0.80, "nucleus"),
}

_ADHESION_TAGS = ("nascent_fa", "collagen_series", "membrane_ecm")

_TURBO = np.array([
    (0.19, 0.07, 0.23), (0.28, 0.31, 0.72), (0.16, 0.62, 0.85), (0.20, 0.83, 0.58),
    (0.65, 0.90, 0.25), (0.96, 0.75, 0.16), (0.95, 0.42, 0.13), (0.72, 0.11, 0.06),
])


def log_norm(force: np.ndarray, ceiling: float, decades: float = FORCE_DECADES) -> np.ndarray:
    """Map ``force`` to ``[0, 1]`` logarithmically over ``decades`` below ``ceiling``.

    Anything at or below the floor maps to exactly 0 rather than to a negative number, so a slack
    element and an unloaded one are the same colour — which is true, and which a clipped linear map
    would also give but only after crushing four decades of real variation into it.
    """
    top = max(float(ceiling), 1.0e-300)
    floor = top * 10.0 ** (-float(decades))
    t = (np.log10(np.maximum(np.asarray(force, dtype=np.float64), floor)) - np.log10(floor))
    return np.clip(t / float(decades), 0.0, 1.0)


def _turbo(t: np.ndarray) -> list[str]:
    t = np.clip(np.asarray(t, dtype=np.float64), 0.0, 1.0) * (len(_TURBO) - 1)
    lo = np.floor(t).astype(int)
    hi = np.minimum(lo + 1, len(_TURBO) - 1)
    f = (t - lo)[:, None]
    rgb = _TURBO[lo] * (1 - f) + _TURBO[hi] * f
    return ["#%02x%02x%02x" % tuple(int(round(255 * c)) for c in row) for row in rgb]


@dataclass(frozen=True, slots=True)
class CellFigure:
    """A rendered figure plus what a viewer needs to drive it.

    Attributes:
        svg: the standalone document.
        layers: one entry per toggleable layer — ``cls``, ``label``, ``colour``.
        force_ceiling_pn: the top of the colour ramp, so a caller can state it.
        spans_um: adhesion rows measured on their own paired sites, ``(name, min, max)``.
    """

    svg: str
    layers: tuple[dict[str, str], ...]
    force_ceiling_pn: float
    spans_um: tuple[tuple[str, float, float], ...]


def render_cell_svg(
    cell: Any,
    *,
    azimuth_deg: float = 34.0,
    elevation_deg: float = 20.0,
    width_px: float = 660.0,
    subtitle: str = "",
) -> CellFigure:
    """Render a :class:`~aleph.scenarios.whole_cell.WholeCell` in side-oblique orthographic view.

    Takes a **built** cell rather than build arguments: the figure is a view of a world somebody else
    composed, and a renderer that builds its own subject can disagree with the run it is illustrating.
    """
    world, owners = cell.world, cell.owners
    forces = {
        name: np.linalg.norm(np.asarray(block, dtype=np.float64), axis=1)
        for name, block in world.evaluate_forces().items()
    }

    az, el = np.radians(azimuth_deg), np.radians(elevation_deg)
    fwd = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    right = np.cross(np.array([0.0, 0.0, 1.0]), fwd)
    right /= np.linalg.norm(right)
    up = np.cross(fwd, right)
    light = fwd * 0.30 + up * 0.80 + right * 0.50
    light /= np.linalg.norm(light)

    def view(p):
        p = np.atleast_2d(np.asarray(p, dtype=np.float64))
        return np.stack([p @ right, p @ up], axis=1), p @ fwd

    blocks = [
        np.asarray(o.positions, dtype=np.float64)
        for o in owners.values()
        if getattr(o, "positions", None) is not None
        and np.asarray(o.positions).ndim == 2
        and np.asarray(o.positions).shape[1] == 3
    ]
    if not blocks:
        raise ValueError("no owner in this cell publishes an (N, 3) position block to draw")
    flat, _ = view(np.concatenate(blocks))
    lo, hi = flat.min(axis=0), flat.max(axis=0)
    scale = float(width_px) / float(max(hi - lo))
    pad_l, pad_t, pad_r, pad_b = 44, 104, 276, 78
    w = int((hi[0] - lo[0]) * scale) + pad_l + pad_r
    h = int((hi[1] - lo[1]) * scale) + pad_t + pad_b

    def screen(p):
        xy, d = view(p)
        out = np.empty_like(xy)
        out[:, 0] = pad_l + (xy[:, 0] - lo[0]) * scale
        out[:, 1] = pad_t + (hi[1] - xy[:, 1]) * scale
        return out, d

    def shade(base, normal, ambient=0.32):
        k = ambient + (1.0 - ambient) * max(0.0, float(normal @ light))
        return "#%02x%02x%02x" % tuple(int(round(255 * min(1.0, c * k))) for c in base)

    draw: list[tuple[float, str]] = []
    layers: list[dict[str, str]] = []

    # -- the ground, so "on a dish" needs no explaining -------------------------------------------
    matrix = owners.get("ecm")
    membrane = owners.get("membrane")
    if matrix is not None and membrane is not None:
        top = float(np.asarray(matrix.positions)[:, 2].max())
        reach = 1.45 * float(np.abs(np.asarray(membrane.positions)[:, :2]).max())
        quad, _ = screen(np.array([[-reach, -reach, top], [reach, -reach, top],
                                   [reach, reach, top], [-reach, reach, top]]))
        draw.append((-1.0e6, '<path d="M' + "L".join(f"{p[0]:.1f},{p[1]:.1f}" for p in quad)
                     + 'Z" fill="#152029" stroke="#2B4053" stroke-width="1.2" class="L-ground"/>'))
        layers.append({"cls": "L-ground", "label": "substrate plane", "colour": "#2B4053"})

    # -- closed surfaces -------------------------------------------------------------------------
    for name, (base, alpha, label) in SURFACE_STYLE.items():
        owner = owners.get(name)
        tri = getattr(owner, "triangles", None) if owner is not None else None
        if owner is None or tri is None:
            continue
        pos = np.asarray(owner.positions, dtype=np.float64)
        tri = np.asarray(tri, dtype=np.int64)
        a, b, c = pos[tri[:, 0]], pos[tri[:, 1]], pos[tri[:, 2]]
        n = np.cross(b - a, c - a)
        length = np.linalg.norm(n, axis=1, keepdims=True)
        n = n / np.where(length == 0.0, 1.0, length)
        mid = (a + b + c) / 3.0
        n[np.einsum("ij,ij->i", n, mid - pos.mean(axis=0)) < 0.0] *= -1.0
        # Far side only for the membrane: an opaque near wall hides the twelve compartments the
        # figure exists to show, and the far wall alone still reads as a closed body.
        keep = (n @ fwd) > 0.0 if name == "membrane" else np.ones(tri.shape[0], dtype=bool)
        pa, _ = screen(a)
        pb, _ = screen(b)
        pc, _ = screen(c)
        _, dm = screen(mid)
        for i in np.flatnonzero(keep):
            draw.append((float(dm[i]),
                         f'<path d="M{pa[i,0]:.2f},{pa[i,1]:.2f}L{pb[i,0]:.2f},{pb[i,1]:.2f}'
                         f'L{pc[i,0]:.2f},{pc[i,1]:.2f}Z" fill="{shade(base, n[i])}" '
                         f'fill-opacity="{alpha}" stroke="none" class="L-{name}"/>'))
        layers.append({"cls": f"L-{name}", "label": f"{label} · {tri.shape[0]} faces",
                       "colour": "#%02x%02x%02x" % tuple(int(255 * x) for x in base)})

    ceiling = max(
        (float(forces[k].max()) for k in FILAMENT_STYLE if k in forces and forces[k].size),
        default=1.0,
    )
    ceiling = max(ceiling, 1.0e-300)

    # -- filaments, along their own segments, coloured by the force they carry --------------------
    for name, (colour, stroke, label) in FILAMENT_STYLE.items():
        owner = owners.get(name)
        if owner is None:
            continue
        pos, seg = getattr(owner, "positions", None), getattr(owner, "segments", None)
        if pos is None or seg is None:
            continue
        pos = np.asarray(pos, dtype=np.float64)
        seg = np.asarray(seg, dtype=np.int64).reshape(-1, 2)
        if not seg.size:
            continue
        f = forces.get(name)
        raw = (0.5 * (f[seg[:, 0]] + f[seg[:, 1]])
               if f is not None and f.size == pos.shape[0] else np.zeros(seg.shape[0]))
        tint = _turbo(log_norm(raw, ceiling))
        p0, d0 = screen(pos[seg[:, 0]])
        p1, d1 = screen(pos[seg[:, 1]])
        for i in range(seg.shape[0]):
            draw.append((float(0.5 * (d0[i] + d1[i])),
                         f'<line x1="{p0[i,0]:.2f}" y1="{p0[i,1]:.2f}" x2="{p1[i,0]:.2f}" '
                         f'y2="{p1[i,1]:.2f}" stroke="{tint[i] if raw[i] > 0.0 else colour}" '
                         f'stroke-width="{stroke}" stroke-linecap="round" opacity="0.94" '
                         f'class="L-{name}"/>'))
        layers.append({"cls": f"L-{name}", "label": f"{label} · {seg.shape[0]} segments",
                       "colour": colour})

    # -- crosslinks: where a bond attaches is a physical fact --------------------------------------
    for name in ("cortex", "ecm"):
        owner = owners.get(name)
        links = list(getattr(owner, "crosslinks", ()) or ()) if owner is not None else []
        drawn = 0
        for link in links:
            pair = None
            for fa, ca, fb, cb in (
                ("filament_a", "material_coordinate_a_um", "filament_b", "material_coordinate_b_um"),
                ("fibre_a", "coordinate_a_um", "fibre_b", "coordinate_b_um"),
            ):
                if not hasattr(link, fa):
                    continue
                try:
                    pair = (
                        owner.resolve_material_point(getattr(link, fa), getattr(link, ca)).position,
                        owner.resolve_material_point(getattr(link, fb), getattr(link, cb)).position,
                    )
                except Exception:  # noqa: BLE001 — a bond that cannot resolve is not drawn
                    pair = None
                break
            if pair is None:
                continue
            q, d = screen(np.stack(pair))
            draw.append((float(d.mean()),
                         f'<line x1="{q[0,0]:.2f}" y1="{q[0,1]:.2f}" x2="{q[1,0]:.2f}" '
                         f'y2="{q[1,1]:.2f}" stroke="#6FE0C4" stroke-width="0.85" opacity="0.6" '
                         f'class="L-{name}xl"/>'))
            drawn += 1
        if drawn:
            layers.append({"cls": f"L-{name}xl", "label": f"{name} crosslinks · {drawn}",
                           "colour": "#6FE0C4"})

    # -- adhesion at true length, always on top ----------------------------------------------------
    spans: list[tuple[str, float, float]] = []
    for slot in world.connectors:
        if not any(tag in slot.name for tag in _ADHESION_TAGS):
            continue
        paired = slot.owners
        if paired is None or len(paired) != 2:
            continue
        try:
            qa = np.asarray(paired[0].positions, dtype=np.float64)
            qb = np.asarray(paired[1].positions, dtype=np.float64)
        except Exception:  # noqa: BLE001
            continue
        if qa.shape != qb.shape or not qa.size:
            continue
        d = np.linalg.norm(qb - qa, axis=1)
        spans.append((slot.name, float(d.min()), float(d.max())))
        ra, _ = screen(qa)
        rb, _ = screen(qb)
        for j in range(qa.shape[0]):
            draw.append((1.0e6,
                         f'<line x1="{ra[j,0]:.2f}" y1="{ra[j,1]:.2f}" x2="{rb[j,0]:.2f}" '
                         f'y2="{rb[j,1]:.2f}" stroke="#FF5C7A" stroke-width="2.4" '
                         f'stroke-linecap="round" opacity="0.96" class="L-adhesion"/>'))
    if spans:
        layers.append({"cls": "L-adhesion", "label": "adhesion, true length", "colour": "#FF5C7A"})

    # -- emit ---------------------------------------------------------------------------------------
    nodes = sum(int(b.shape[0]) for b in blocks)
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="100%" '
        f'style="max-width:{w}px;height:auto" '
        f'font-family="ui-monospace,SFMono-Regular,Menlo,monospace">',
        f'<rect width="{w}" height="{h}" fill="#0C1116"/>',
    ]
    svg += [s for _, s in sorted(draw, key=lambda kv: kv[0])]

    svg.append(f'<text x="{pad_l}" y="28" font-size="15" font-weight="700" fill="#E6EDF3">'
               f'{subtitle or "Whole cell"}</text>')
    svg.append(f'<text x="{pad_l}" y="50" font-size="11.5" fill="#8FA3B0">'
               f'{len(cell.present)} compartments · {len(cell.couplings)} connectors · {nodes} '
               f'nodes · surfaces as faces, filaments as segments, crosslinks as bonds · '
               f'orthographic, az {azimuth_deg:.0f}° el {elevation_deg:.0f}°</text>')
    svg.append(f'<text x="{pad_l}" y="70" font-size="11.5" fill="#8FA3B0">'
               f'filament colour = per-node |F|, log over {FORCE_DECADES:.0f} decades to '
               f'{ceiling:.4g} pN · FULL-RES, no downsample</text>')

    y = pad_t + 4
    svg.append(f'<text x="{w - pad_r + 14}" y="{y}" font-size="11.5" font-weight="700" '
               f'fill="#E6EDF3">layers</text>')
    for layer in layers:
        y += 19
        svg.append(f'<rect x="{w - pad_r + 14}" y="{y - 9}" width="11" height="11" rx="2" '
                   f'fill="{layer["colour"]}"/>')
        svg.append(f'<text x="{w - pad_r + 32}" y="{y}" font-size="10.5" fill="#C9D6DF">'
                   f'{layer["label"]}</text>')

    if spans:
        y += 26
        svg.append(f'<text x="{w - pad_r + 14}" y="{y}" font-size="11.5" font-weight="700" '
                   f'fill="#E6EDF3">adhesion span, µm</text>')
        for name, a, b in sorted(spans):
            y += 16
            svg.append(f'<text x="{w - pad_r + 14}" y="{y}" font-size="9.5" fill="#B4707F">'
                       f'{name[:26]}</text>')
            y += 13
            svg.append(f'<text x="{w - pad_r + 24}" y="{y}" font-size="11" font-weight="700" '
                       f'fill="#FF5C7A">{a:.3f} – {b:.3f}</text>')
        y += 22
        svg.append(f'<text x="{w - pad_r + 14}" y="{y}" font-size="10" fill="#8FA3B0">'
                   f'integrin, real: 0.020 – 0.030</text>')

    bx, by, bh = w - pad_r + 16, h - pad_b - 150, 130
    for k in range(60):
        svg.append(f'<rect x="{bx}" y="{by + bh - (k + 1) * bh / 60:.2f}" width="13" '
                   f'height="{bh / 60 + 0.6:.2f}" fill="{_turbo(np.array([k / 59.0]))[0]}"/>')
    svg.append(f'<text x="{bx}" y="{by - 8}" font-size="10" fill="#C9D6DF">'
               f'per-node |F| [pN, log]</text>')
    for k in range(int(FORCE_DECADES) + 1):
        t = k / FORCE_DECADES
        svg.append(f'<text x="{bx + 18}" y="{by + bh - t * bh + 3.5:.1f}" font-size="9.5" '
                   f'fill="#8FA3B0">{ceiling * 10.0 ** (k - FORCE_DECADES):.3g}</text>')

    # The scale bar is not optional; see the module docstring.
    bar_um = 10.0 ** np.floor(np.log10(0.3 * float(max(hi - lo))))
    if bar_um * scale < 60.0:
        bar_um *= 2.0
    px = bar_um * scale
    sx, sy = pad_l + 8, h - 34
    svg.append(f'<rect x="{sx}" y="{sy}" width="{px:.1f}" height="3" rx="1.5" fill="#C9D6DF"/>')
    svg.append(f'<text x="{sx + px / 2:.1f}" y="{sy + 18}" font-size="11" fill="#C9D6DF" '
               f'text-anchor="middle">{bar_um:g} µm</text>')
    svg.append("</svg>")

    return CellFigure(
        svg="\n".join(svg),
        layers=tuple(layers),
        force_ceiling_pn=float(ceiling),
        spans_um=tuple(spans),
    )


def layer_toggle_page(figure: CellFigure, *, title: str = "Aleph — whole cell") -> str:
    """A self-contained page with one checkbox per layer. No external asset of any kind.

    Toggling is by CSS class rather than by regrouping the SVG: the depth order is a **global**
    painter's-algorithm sort, and moving each layer into its own ``<g>`` would reorder the document
    and put near geometry behind far geometry.
    """
    rows = "\n".join(
        f'<label><input type="checkbox" checked data-cls="{layer["cls"]}">'
        f'<span class="sw" style="background:{layer["colour"]}"></span>{layer["label"]}</label>'
        for layer in figure.layers
    )
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>{title}</title>
<style>
 body{{margin:0;background:#0C1116;color:#C9D6DF;
      font:13px ui-monospace,SFMono-Regular,Menlo,monospace}}
 #panel{{position:fixed;top:10px;left:10px;background:#141B21ee;border:1px solid #26313A;
        border-radius:8px;padding:10px 12px;max-height:92vh;overflow:auto;z-index:9}}
 #panel b{{display:block;margin-bottom:6px}}
 label{{display:flex;align-items:center;gap:6px;margin:3px 0;cursor:pointer;white-space:nowrap}}
 .sw{{width:11px;height:11px;border-radius:2px;flex:none}}
 svg{{margin-left:264px}}
 a{{color:#58a6ff;cursor:pointer}}
</style></head><body>
<div id="panel"><b>layers &nbsp;<a id="all">all</a> / <a id="none">none</a></b>{rows}</div>
{figure.svg}
<script>
const set=(c,on)=>document.querySelectorAll('.'+c).forEach(e=>e.style.display=on?'':'none');
const boxes=document.querySelectorAll('#panel input');
boxes.forEach(i=>i.onchange=()=>set(i.dataset.cls,i.checked));
document.getElementById('all').onclick=()=>boxes.forEach(i=>{{i.checked=true;set(i.dataset.cls,true)}});
document.getElementById('none').onclick=()=>boxes.forEach(i=>{{i.checked=false;set(i.dataset.cls,false)}});
</script></body></html>"""


def layers_json(figure: CellFigure) -> str:
    """The layer list, for a caller building its own viewer."""
    return json.dumps(list(figure.layers), indent=1)


# ==================================================================================================
# The interactive scene
#
# `render_cell_svg` emits one element per segment. At 324 cortex segments that is a document; at the
# sourced density it is **37 MB and 31,165 drawable elements**, measured 2026-08-07, and no browser
# opens it. The bytes were never the problem — 400 MB artifacts have been opened on this project.
# **The element count is**, because each one becomes a laid-out, hit-tested object.
#
# Decimation would also make the page open, and would be the wrong repair: the density is the thing
# being modelled, so thinning it produces a picture of a cell nobody ran. What changes here is the
# **representation** — geometry travels as base64 `Float32Array`/`Int32Array` and **one** canvas
# draws it, so the DOM cost is constant in the geometry and a 400 MB payload is merely a large file.
#
# Carrying 3-D rather than a projection is not a bonus, it is the point of doing it at all: the SVG
# has one viewpoint baked in and cannot be interrogated, and *"a stretched segment and a slack one
# are the same two dots"* is only half the problem — a segment hidden behind the nucleus is no dots
# at all until the camera moves.
# ==================================================================================================


@dataclass(frozen=True, slots=True)
class CellScene:
    """Geometry as data, in µm, with no camera applied.

    Attributes:
        payload: the JSON a viewer consumes. Embedded verbatim by :func:`interactive_page`.
        force_ceiling_pn: the top of the colour ramp — **the same value the SVG uses**, so the two
            renderings of one cell cannot disagree about what red means.
        segment_count: total segments carried, before any drawing.
        node_count: total published positions.
        spans_um: adhesion rows, ``(name, min, max)``, as :class:`CellFigure` reports them.
    """

    payload: str
    force_ceiling_pn: float
    segment_count: int
    node_count: int
    spans_um: tuple[tuple[str, float, float], ...]


def _b64(array: np.ndarray, dtype) -> str:
    return base64.b64encode(np.ascontiguousarray(array, dtype=dtype).tobytes()).decode("ascii")


def render_cell_scene(cell: Any) -> CellScene:
    """Extract a :class:`~aleph.scenarios.whole_cell.WholeCell` as drawable data.

    Takes a **built** cell for the reason :func:`render_cell_svg` does: a renderer that builds its
    own subject can disagree with the run it is illustrating.

    Nothing is thinned, sorted or projected. Depth sorting is a property of a viewpoint and there is
    no viewpoint here; the viewer owns the camera and does its own ordering.
    """
    world, owners = cell.world, cell.owners
    forces = {
        name: np.linalg.norm(np.asarray(block, dtype=np.float64), axis=1)
        for name, block in world.evaluate_forces().items()
    }
    ceiling = max(
        (float(forces[k].max()) for k in FILAMENT_STYLE if k in forces and forces[k].size),
        default=1.0,
    )
    ceiling = max(ceiling, 1.0e-300)

    layers: list[dict[str, Any]] = []
    segments = 0
    nodes = 0

    # -- the ground, as a quad. Same reason as the SVG: "on a dish" should need no explaining. -----
    matrix, membrane = owners.get("ecm"), owners.get("membrane")
    if matrix is not None and membrane is not None:
        top = float(np.asarray(matrix.positions)[:, 2].max())
        reach = 1.45 * float(np.abs(np.asarray(membrane.positions)[:, :2]).max())
        quad = np.array([[-reach, -reach, top], [reach, -reach, top],
                         [reach, reach, top], [-reach, reach, top]])
        layers.append({
            "cls": "L-ground", "label": "substrate plane", "colour": "#2B4053", "kind": "quad",
            "positions": _b64(quad, np.float32),
        })

    # -- closed surfaces, as triangles. Normals are computed in the viewer, which has the camera. --
    for name, (base, alpha, label) in SURFACE_STYLE.items():
        owner = owners.get(name)
        tri = getattr(owner, "triangles", None) if owner is not None else None
        if owner is None or tri is None:
            continue
        pos = np.asarray(owner.positions, dtype=np.float64)
        tri = np.asarray(tri, dtype=np.int64).reshape(-1, 3)
        nodes += int(pos.shape[0])
        layers.append({
            "cls": f"L-{name}", "label": f"{label} · {tri.shape[0]} faces", "kind": "surface",
            "colour": "#%02x%02x%02x" % tuple(int(255 * x) for x in base),
            "rgb": [float(x) for x in base], "alpha": float(alpha),
            # The membrane is drawn far-side-only, but *which* side is far depends on the camera, so
            # the flag travels and the viewer decides per frame. The SVG could pre-compute it; a
            # scene that can be rotated cannot.
            "cull": bool(name == "membrane"),
            "positions": _b64(pos, np.float32), "triangles": _b64(tri, np.int32),
        })

    # -- filaments, with the per-segment force so the viewer can colour without recomputing --------
    for name, (colour, stroke, label) in FILAMENT_STYLE.items():
        owner = owners.get(name)
        if owner is None:
            continue
        pos, seg = getattr(owner, "positions", None), getattr(owner, "segments", None)
        if pos is None or seg is None:
            continue
        pos = np.asarray(pos, dtype=np.float64)
        seg = np.asarray(seg, dtype=np.int64).reshape(-1, 2)
        if not seg.size:
            continue
        f = forces.get(name)
        raw = (0.5 * (f[seg[:, 0]] + f[seg[:, 1]])
               if f is not None and f.size == pos.shape[0] else np.zeros(seg.shape[0]))
        segments += int(seg.shape[0])
        nodes += int(pos.shape[0])
        layers.append({
            "cls": f"L-{name}", "label": f"{label} · {seg.shape[0]} segments", "kind": "filament",
            "colour": colour, "width": float(stroke),
            "positions": _b64(pos, np.float32), "segments": _b64(seg, np.int32),
            "force_pn": _b64(raw, np.float32),
        })

    # -- crosslinks: where a bond attaches is a physical fact ---------------------------------------
    for name in ("cortex", "ecm"):
        owner = owners.get(name)
        links = list(getattr(owner, "crosslinks", ()) or ()) if owner is not None else []
        ends: list[np.ndarray] = []
        for link in links:
            for fa, ca, fb, cb in (
                ("filament_a", "material_coordinate_a_um", "filament_b", "material_coordinate_b_um"),
                ("fibre_a", "coordinate_a_um", "fibre_b", "coordinate_b_um"),
            ):
                if not hasattr(link, fa):
                    continue
                try:
                    ends.append(np.stack([
                        owner.resolve_material_point(getattr(link, fa), getattr(link, ca)).position,
                        owner.resolve_material_point(getattr(link, fb), getattr(link, cb)).position,
                    ]))
                except Exception:  # noqa: BLE001 — a bond that cannot resolve is not drawn
                    pass
                break
        if ends:
            stacked = np.concatenate(ends).reshape(-1, 3)
            pairs = np.arange(stacked.shape[0], dtype=np.int64).reshape(-1, 2)
            segments += int(pairs.shape[0])
            layers.append({
                "cls": f"L-{name}xl", "label": f"{name} crosslinks · {pairs.shape[0]}",
                "kind": "filament", "colour": "#6FE0C4", "width": 0.85, "opacity": 0.6,
                "positions": _b64(stacked, np.float32), "segments": _b64(pairs, np.int32),
                "force_pn": _b64(np.zeros(pairs.shape[0]), np.float32),
            })

    # -- adhesion at true length ---------------------------------------------------------------------
    spans: list[tuple[str, float, float]] = []
    ends = []
    for slot in world.connectors:
        if not any(tag in slot.name for tag in _ADHESION_TAGS):
            continue
        paired = slot.owners
        if paired is None or len(paired) != 2:
            continue
        try:
            qa = np.asarray(paired[0].positions, dtype=np.float64)
            qb = np.asarray(paired[1].positions, dtype=np.float64)
        except Exception:  # noqa: BLE001
            continue
        if qa.shape != qb.shape or not qa.size:
            continue
        spans.append((slot.name, float(np.linalg.norm(qb - qa, axis=1).min()),
                      float(np.linalg.norm(qb - qa, axis=1).max())))
        ends.append(np.stack([qa, qb], axis=1).reshape(-1, 3))
    if ends:
        stacked = np.concatenate(ends)
        pairs = np.arange(stacked.shape[0], dtype=np.int64).reshape(-1, 2)
        segments += int(pairs.shape[0])
        layers.append({
            "cls": "L-adhesion", "label": "adhesion, true length", "kind": "filament",
            "colour": "#FF5C7A", "width": 2.4, "always_on_top": True,
            "positions": _b64(stacked, np.float32), "segments": _b64(pairs, np.int32),
            "force_pn": _b64(np.zeros(pairs.shape[0]), np.float32),
        })

    payload = json.dumps({
        # No `azimuth_deg`, no `width_px`. A projected scene is a picture; the viewer owns the camera.
        "units": "um",
        "force_decades": FORCE_DECADES,
        "force_ceiling_pn": float(ceiling),
        "turbo": [[float(c) for c in row] for row in _TURBO],
        "compartments": len(cell.present),
        "connectors": len(cell.couplings),
        "nodes": int(nodes),
        "segments": int(segments),
        "spans_um": [[n, a, b] for n, a, b in spans],
        "layers": layers,
    }, separators=(",", ":"))

    return CellScene(
        payload=payload,
        force_ceiling_pn=float(ceiling),
        segment_count=int(segments),
        node_count=int(nodes),
        spans_um=tuple(spans),
    )


#: Everything drawable goes on the canvas. What is left is chrome, and it is a **fixed** list —
#: which is the property `test_a_page_costs_a_bounded_number_of_elements` pins.
_VIEWER_JS = r"""
const S = JSON.parse(document.getElementById('scene').textContent);
const cv = document.getElementById('cv'), cx = cv.getContext('2d');
const dec = (s, T) => { const b = atob(s), u = new Uint8Array(b.length);
  for (let i = 0; i < b.length; i++) u[i] = b.charCodeAt(i); return new T(u.buffer); };
for (const L of S.layers) {
  if (L.positions) L.P = dec(L.positions, Float32Array);
  if (L.segments)  L.S = dec(L.segments, Int32Array);
  if (L.triangles) L.T = dec(L.triangles, Int32Array);
  if (L.force_pn)  L.F = dec(L.force_pn, Float32Array);
  L.on = true;
}
// Bounds from every published position, so the fit does not depend on which layers are visible.
let lo = [1e30, 1e30, 1e30], hi = [-1e30, -1e30, -1e30];
for (const L of S.layers) for (let i = 0; L.P && i < L.P.length; i += 3)
  for (let k = 0; k < 3; k++) { const v = L.P[i + k];
    if (v < lo[k]) lo[k] = v; if (v > hi[k]) hi[k] = v; }
const mid = [0, 1, 2].map(k => 0.5 * (lo[k] + hi[k]));
const extent = Math.max(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1;

let az = 34 * Math.PI / 180, el = 20 * Math.PI / 180, zoom = 0.86, panx = 0, pany = 0;
const TB = S.turbo, DEC = S.force_decades, CEIL = S.force_ceiling_pn;
// The same log map as `log_norm`, to the same ceiling over the same decades. Two renderings of one
// cell that disagree about what red means are worse than one rendering.
function tint(f) {
  if (!(f > 0)) return null;
  let t = 1 + Math.log10(f / CEIL) / DEC;
  t = t < 0 ? 0 : t > 1 ? 1 : t;
  const x = t * (TB.length - 1), i = Math.min(TB.length - 2, Math.floor(x)), u = x - i;
  const c = k => Math.round(255 * (TB[i][k] * (1 - u) + TB[i + 1][k] * u));
  return 'rgb(' + c(0) + ',' + c(1) + ',' + c(2) + ')';
}

function draw() {
  const W = cv.width = cv.clientWidth * devicePixelRatio;
  const H = cv.height = cv.clientHeight * devicePixelRatio;
  const ca = Math.cos(az), sa = Math.sin(az), ce = Math.cos(el), se = Math.sin(el);
  const fwd = [ce * ca, ce * sa, se];
  let rt = [-fwd[1], fwd[0], 0]; const rn = Math.hypot(rt[0], rt[1]) || 1;
  rt = [rt[0] / rn, rt[1] / rn, 0];
  const up = [fwd[1] * rt[2] - fwd[2] * rt[1], fwd[2] * rt[0] - fwd[0] * rt[2],
              fwd[0] * rt[1] - fwd[1] * rt[0]];
  const s = zoom * Math.min(W, H) / extent;
  const px = (P, i) => { const x = P[i] - mid[0], y = P[i+1] - mid[1], z = P[i+2] - mid[2];
    return [W / 2 + panx + s * (x * rt[0] + y * rt[1] + z * rt[2]),
            H / 2 + pany - s * (x * up[0] + y * up[1] + z * up[2]),
            x * fwd[0] + y * fwd[1] + z * fwd[2]]; };
  const light = (() => { const l = [fwd[0]*.3 + up[0]*.8 + rt[0]*.5, fwd[1]*.3 + up[1]*.8 + rt[1]*.5,
    fwd[2]*.3 + up[2]*.8 + rt[2]*.5]; const n = Math.hypot(l[0], l[1], l[2]);
    return [l[0]/n, l[1]/n, l[2]/n]; })();

  cx.fillStyle = '#0C1116'; cx.fillRect(0, 0, W, H);
  cx.lineCap = 'round';
  const top = [];
  for (const L of S.layers) {
    if (!L.on) continue;
    if (L.always_on_top) { top.push(L); continue; }
    paint(L);
  }
  for (const L of top) paint(L);

  function paint(L) {
    if (L.kind === 'quad') {
      cx.beginPath();
      for (let i = 0; i < L.P.length; i += 3) { const q = px(L.P, i);
        i ? cx.lineTo(q[0], q[1]) : cx.moveTo(q[0], q[1]); }
      cx.closePath(); cx.fillStyle = '#152029'; cx.fill();
      cx.strokeStyle = '#2B4053'; cx.lineWidth = 1.2 * devicePixelRatio; cx.stroke(); return;
    }
    if (L.kind === 'surface') {
      // Depth-sorted per frame. The painter's order is a property of the camera, so it cannot be
      // baked into the payload the way the SVG baked it into the document.
      const n = L.T.length / 3, ord = new Array(n), dep = new Float32Array(n);
      for (let t = 0; t < n; t++) {
        const a = 3*L.T[3*t], b = 3*L.T[3*t+1], c = 3*L.T[3*t+2];
        dep[t] = (px(L.P,a)[2] + px(L.P,b)[2] + px(L.P,c)[2]) / 3; ord[t] = t;
      }
      ord.sort((p, q) => dep[p] - dep[q]);
      cx.globalAlpha = L.alpha;
      for (const t of ord) {
        const a = 3*L.T[3*t], b = 3*L.T[3*t+1], c = 3*L.T[3*t+2];
        const A = px(L.P,a), B = px(L.P,b), C = px(L.P,c);
        const ux = L.P[b]-L.P[a], uy = L.P[b+1]-L.P[a+1], uz = L.P[b+2]-L.P[a+2];
        const vx = L.P[c]-L.P[a], vy = L.P[c+1]-L.P[a+1], vz = L.P[c+2]-L.P[a+2];
        let nx = uy*vz - uz*vy, ny = uz*vx - ux*vz, nz = ux*vy - uy*vx;
        const nl = Math.hypot(nx, ny, nz) || 1; nx/=nl; ny/=nl; nz/=nl;
        const cxm = (L.P[a]+L.P[b]+L.P[c])/3 - mid[0], cym = (L.P[a+1]+L.P[b+1]+L.P[c+1])/3 - mid[1],
              czm = (L.P[a+2]+L.P[b+2]+L.P[c+2])/3 - mid[2];
        if (nx*cxm + ny*cym + nz*czm < 0) { nx=-nx; ny=-ny; nz=-nz; }
        if (L.cull && (nx*fwd[0] + ny*fwd[1] + nz*fwd[2]) <= 0) continue;
        const k = 0.32 + 0.68 * Math.max(0, nx*light[0] + ny*light[1] + nz*light[2]);
        cx.fillStyle = 'rgb(' + L.rgb.map(v => Math.round(255 * Math.min(1, v * k))).join(',') + ')';
        cx.beginPath(); cx.moveTo(A[0],A[1]); cx.lineTo(B[0],B[1]); cx.lineTo(C[0],C[1]);
        cx.closePath(); cx.fill();
      }
      cx.globalAlpha = 1; return;
    }
    // filaments — one path per colour bucket, not one per segment
    cx.globalAlpha = L.opacity || 0.94;
    cx.lineWidth = L.width * devicePixelRatio * (zoom / 0.86 > 2 ? 1.6 : 1);
    const buckets = new Map();
    for (let e = 0; e < L.S.length; e += 2) {
      const col = tint(L.F ? L.F[e >> 1] : 0) || L.colour;
      let arr = buckets.get(col); if (!arr) buckets.set(col, arr = []);
      arr.push(e);
    }
    for (const [col, list] of buckets) {
      cx.strokeStyle = col; cx.beginPath();
      for (const e of list) {
        const A = px(L.P, 3 * L.S[e]), B = px(L.P, 3 * L.S[e + 1]);
        cx.moveTo(A[0], A[1]); cx.lineTo(B[0], B[1]);
      }
      cx.stroke();
    }
    cx.globalAlpha = 1;
  }

  // The scale bar is not optional; see the module docstring. It is redrawn per frame because zoom
  // changes what a micron is worth in pixels, and a ruler that does not track the zoom is a lie.
  let bar = Math.pow(10, Math.floor(Math.log10(0.3 * extent)));
  while (bar * s < 60 * devicePixelRatio) bar *= 2;
  const bx = 22 * devicePixelRatio, by = H - 34 * devicePixelRatio;
  cx.fillStyle = '#C9D6DF'; cx.fillRect(bx, by, bar * s, 3 * devicePixelRatio);
  cx.font = (12 * devicePixelRatio) + 'px ui-monospace,Menlo,monospace';
  cx.textAlign = 'center';
  cx.fillText(bar + ' µm', bx + bar * s / 2, by + 18 * devicePixelRatio);
  cx.textAlign = 'left';
  document.getElementById('cam').textContent =
    'az ' + (az * 180 / Math.PI).toFixed(0) + '°  el ' + (el * 180 / Math.PI).toFixed(0) +
    '°  zoom ' + zoom.toFixed(2) + '×';
}

let drag = null;
cv.addEventListener('pointerdown', e => { drag = [e.clientX, e.clientY, e.shiftKey];
  cv.setPointerCapture(e.pointerId); });
cv.addEventListener('pointerup', () => drag = null);
cv.addEventListener('pointermove', e => {
  if (!drag) return;
  const dx = e.clientX - drag[0], dy = e.clientY - drag[1];
  if (drag[2]) { panx += dx * devicePixelRatio; pany += dy * devicePixelRatio; }
  else { az -= dx * 0.008; el = Math.max(-1.53, Math.min(1.53, el + dy * 0.008)); }
  drag = [e.clientX, e.clientY, drag[2]]; draw();
});
cv.addEventListener('wheel', e => { e.preventDefault();
  zoom *= Math.exp(-e.deltaY * 0.0014); zoom = Math.max(0.05, Math.min(400, zoom)); draw(); },
  {passive: false});
document.querySelectorAll('#panel input').forEach(i => i.onchange = () => {
  S.layers.find(L => L.cls === i.dataset.cls).on = i.checked; draw(); });
const flip = v => document.querySelectorAll('#panel input').forEach(i => {
  i.checked = v; S.layers.find(L => L.cls === i.dataset.cls).on = v; draw(); });
document.getElementById('all').onclick = () => flip(true);
document.getElementById('none').onclick = () => flip(false);
document.getElementById('reset').onclick = () => {
  az = 34 * Math.PI / 180; el = 20 * Math.PI / 180; zoom = 0.86; panx = pany = 0; draw(); };
addEventListener('resize', draw);
draw();
"""


def interactive_page(scene: CellScene, *, title: str = "Aleph — whole cell") -> str:
    """A standalone, rotatable page. One canvas, no fetch, no external asset of any kind.

    The DOM cost is **constant in the geometry**: a cortex of 324 segments and one of 31,165 produce
    the same number of elements, which is why this opens at a density the SVG cannot.
    """
    meta = json.loads(scene.payload)
    rows = "\n".join(
        f'<label><input type="checkbox" checked data-cls="{lay["cls"]}">'
        f'<span class="sw" style="background:{lay["colour"]}"></span>{lay["label"]}</label>'
        for lay in meta["layers"]
    )
    spans = "".join(
        f'<div class="sp"><i>{n[:30]}</i><b>{a:.3f} – {b:.3f}</b></div>'
        for n, a, b in sorted(meta["spans_um"])
    )
    ramp = "".join(
        f'<i style="background:{_turbo(np.array([k / 8.0]))[0]}"></i>' for k in range(9)
    )
    ceil = meta["force_ceiling_pn"]
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>{title}</title>
<style>
 html,body{{margin:0;height:100%;background:#0C1116;color:#C9D6DF;overflow:hidden;
      font:13px ui-monospace,SFMono-Regular,Menlo,monospace}}
 #cv{{position:fixed;inset:0;width:100%;height:100%;cursor:grab;touch-action:none}}
 #cv:active{{cursor:grabbing}}
 #panel{{position:fixed;top:10px;left:10px;background:#141B21ee;border:1px solid #26313A;
        border-radius:8px;padding:10px 12px;max-height:94vh;overflow:auto;z-index:9;
        backdrop-filter:blur(3px)}}
 #panel b.h{{display:block;margin:8px 0 5px}}
 label{{display:flex;align-items:center;gap:6px;margin:3px 0;cursor:pointer;white-space:nowrap}}
 .sw{{width:11px;height:11px;border-radius:2px;flex:none}}
 a{{color:#58a6ff;cursor:pointer}}
 .sp{{margin:2px 0;font-size:11px}} .sp i{{color:#B4707F;font-style:normal}}
 .sp b{{display:block;color:#FF5C7A;margin-left:8px}}
 #ramp{{display:flex;height:10px;border-radius:2px;overflow:hidden;margin:4px 0 2px}}
 #ramp i{{flex:1}}
 .dim{{color:#8FA3B0;font-size:10.5px;line-height:1.45}}
 #hud{{position:fixed;right:12px;bottom:12px;color:#8FA3B0;font-size:11px;text-align:right;z-index:9}}
</style></head><body>
<canvas id="cv"></canvas>
<div id="panel">
 <b class="h">layers &nbsp;<a id="all">all</a> / <a id="none">none</a> / <a id="reset">reset view</a></b>
 {rows}
 <b class="h">per-node |F| [pN, log]</b>
 <div id="ramp">{ramp}</div>
 <div class="dim">{ceil * 10.0 ** -FORCE_DECADES:.3g} &rarr; {ceil:.4g} pN
  &middot; {FORCE_DECADES:.0f} decades</div>
 <b class="h">adhesion span, &micro;m</b>{spans}
 <div class="dim">integrin, real: 0.020 – 0.030</div>
 <b class="h">scene</b>
 <div class="dim">{meta["compartments"]} compartments &middot; {meta["connectors"]} connectors<br>
  {meta["nodes"]} nodes &middot; {meta["segments"]} segments<br>
  FULL-RES, no downsample</div>
 <div class="dim" style="margin-top:7px">drag rotate &middot; shift-drag pan &middot; wheel zoom</div>
</div>
<div id="hud"><span id="cam"></span></div>
<script type="application/json" id="scene">{scene.payload}</script>
<script>{_VIEWER_JS}</script>
</body></html>"""
