"""Draw the load path: every owner from its own arrays, and every isolation state distinctly.

`ALEPH-PORT-3602`.

Two failures this module is shaped to make impossible
-----------------------------------------------------
**1. An owner that is silently omitted.** A renderer that cannot read an owner and quietly draws the
rest produces a picture that looks complete. The reader has no way to tell the difference between "this
owner holds nothing" and "the renderer could not read it", and the second is the interesting case. So
:func:`draw_scene` **partitions** the scene's components: every one is either drawn from its own arrays
or named in :attr:`DrawnScene.undrawn` with a reason. Neither set may be missing a component and the two
may not overlap -- asserted by
``tests/viz/test_render.py::test_every_component_is_either_drawn_or_named_undrawn``.

**2. A defect drawn to look like a decision.** :class:`~aleph.viz.scene.IsolationState` distinguishes
five situations, and the reason it does is that ``NOT_SCHEDULED`` -- registered and never evaluated --
is a *hole in the schedule* while ``ISOLATED`` is somebody's deliberate choice. `scene.py` calls the
first "the failure this project exists to make visible". A viewer that greys both out has hidden the
only one of the five that is a defect. So the five styles differ pairwise on every visual channel this
module has, and ``NOT_SCHEDULED`` versus ``ISOLATED`` is checked separately and harder.

``ENGAGED_ZERO`` is the subtle one. It means a unilateral element sitting on its inactive branch --
a tensile-only tether under compression -- carrying **exactly** zero rather than a small number. The
exactness is the evidence, so it is drawn with its own glyph and labelled ``0`` rather than being left
to look like a very thin ``COUPLED``. ``scene.ConnectorNode`` refuses the state outright if either
force is non-zero, so by the time anything reaches this module the claim has already been enforced.

Individual elements, not a shaded surface
-----------------------------------------
Filaments, fibres, crosslinks and clutches are drawn one by one from the arrays their owner holds.
A shaded surface is a picture of one owner; this draws all of them, and the legend carries the element
counts so the quantitative record is a number and not a pixel.

Output is **SVG**: it is text, so it diffs and reviews; it needs no rasteriser, and this package has
none; and matplotlib is neither installed nor declared in ``pyproject.toml``. Lines are what SVG is
best at, and a filament is a line.

Honest limits
-------------
The projection is a fixed axis-aligned drop of one coordinate, with depth mapped to opacity. **This is
a legibility aid and not a hidden-line computation** -- a dense population overdraws, and when it does,
the legend's element counts are the record and the picture is not. Units: um throughout.
"""

from __future__ import annotations

import html
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from aleph.viz.scene import IsolationState, SceneTree

__all__ = [
    "AXIS_CHOICES",
    "ISOLATION_STYLE",
    "DrawnElement",
    "DrawnScene",
    "IsolationStyle",
    "OwnerElements",
    "draw_scene",
    "owner_elements_from",
    "render_scene_svg",
]


@dataclass(frozen=True, slots=True)
class IsolationStyle:
    """How one isolation state is drawn.

    Every field is a separate visual channel, and the five states are required to differ on more than
    one of them. Colour alone would be a single point of failure -- a greyscale print, a colour-blind
    reader, or a viewer that themes the palette would collapse the distinction that matters most.

    Attributes:
        stroke: Line colour.
        dash: SVG ``stroke-dasharray``; the empty string is a solid line.
        width: Stroke width in px.
        marker: A short glyph drawn at the element's midpoint; the empty string draws none.
        opacity: Base opacity before any depth cue.
        label: Legend text.
        is_fault: Whether this state is a defect rather than a decision.
    """

    stroke: str
    dash: str
    width: float
    marker: str
    opacity: float
    label: str
    is_fault: bool = False

    def channels(self) -> tuple[str, str, float, str]:
        """The four channels distinctness is judged on."""
        return (self.stroke, self.dash, self.width, self.marker)


#: One style per isolation state. See the module docstring for why the differences are redundant.
ISOLATION_STYLE: dict[IsolationState, IsolationStyle] = {
    IsolationState.COUPLED: IsolationStyle(
        stroke="#2e7d32",
        dash="",
        width=2.2,
        marker="",
        opacity=1.0,
        label="coupled — evaluated, carrying load",
    ),
    IsolationState.ENGAGED_ZERO: IsolationStyle(
        # Solid, because it WAS evaluated -- the schedule did its job. Distinguished from `coupled`
        # by hue, weight and an explicit `0`, because "carried exactly nothing" is a measurement and
        # deserves to be shown as one rather than as a faint version of carrying something.
        stroke="#0277bd",
        dash="",
        width=1.2,
        marker="0",
        opacity=0.95,
        label="engaged_zero — evaluated, exactly zero (inactive branch)",
    ),
    IsolationState.ISOLATED: IsolationStyle(
        stroke="#8d6e63",
        dash="7 4",
        width=1.6,
        marker="/",
        opacity=0.75,
        label="isolated — deliberately decoupled, with a reason",
    ),
    IsolationState.NOT_SCHEDULED: IsolationStyle(
        # The one state that is a defect. Loud on purpose, and different from `isolated` on stroke,
        # dash, width AND marker -- see `test_not_scheduled_never_renders_like_isolated`.
        stroke="#c62828",
        dash="2 3",
        width=3.4,
        marker="!",
        opacity=1.0,
        label="NOT_SCHEDULED — registered and never evaluated (a defect)",
        is_fault=True,
    ),
    IsolationState.EXCLUDED: IsolationStyle(
        stroke="#90a4ae",
        dash="1 5",
        width=1.0,
        marker="x",
        opacity=0.55,
        label="excluded — declared out of scope, not a fault",
    ),
}

#: Which pair of world axes the projection keeps. Index into (x, y, z); the third becomes depth.
AXIS_CHOICES: dict[str, tuple[int, int, int]] = {
    "xy": (0, 1, 2),
    "xz": (0, 2, 1),
    "yz": (1, 2, 0),
}


@dataclass(frozen=True, slots=True)
class OwnerElements:
    """One owner's own arrays, as the renderer reads them.

    Attributes:
        name: The owner's registration name.
        positions: ``(N, 3)`` node positions [um], or ``None`` when the owner is not spatially
            resolved. A field owner legitimately has no nodes; that is a property of the component and
            is reported as an undrawn reason rather than treated as a gap in the renderer.
        segments: ``(M, 2)`` node-index pairs -- filaments, fibres, struts.
        crosslinks: ``(K, 2)`` node-index pairs drawn as crosslinks.
        clutches: ``(J,)`` node indices drawn as point elements.
        colour: Stroke colour for this owner's elements.
    """

    name: str
    positions: np.ndarray | None = None
    segments: np.ndarray | None = None
    crosslinks: np.ndarray | None = None
    clutches: np.ndarray | None = None
    colour: str = "#546e7a"


@dataclass(frozen=True, slots=True)
class DrawnElement:
    """One drawable primitive in world coordinates [um]."""

    owner: str
    kind: str
    start: np.ndarray
    end: np.ndarray | None
    colour: str


@dataclass(slots=True)
class DrawnScene:
    """What could be drawn, and -- named rather than omitted -- what could not.

    Attributes:
        tree: The scene tree this was drawn against.
        elements: Every drawable primitive.
        undrawn: ``(owner, reason)`` for every component that could not be drawn from its arrays.
        emergence: Optional per-owner emergence observations, as JSON-able objects.
    """

    tree: SceneTree
    elements: list[DrawnElement] = field(default_factory=list)
    undrawn: list[tuple[str, str]] = field(default_factory=list)
    emergence: dict[str, Any] = field(default_factory=dict)

    @property
    def drawn_owners(self) -> set[str]:
        """Owners that contributed at least one element."""
        return {element.owner for element in self.elements}

    @property
    def undrawn_owners(self) -> set[str]:
        """Owners named as undrawn."""
        return {name for name, _reason in self.undrawn}

    def counts(self) -> dict[str, int]:
        """Element count per ``owner/kind`` -- the quantitative record the picture is not."""
        out: dict[str, int] = {}
        for element in self.elements:
            key = f"{element.owner}/{element.kind}"
            out[key] = out.get(key, 0) + 1
        return out

    def coverage_problems(self) -> list[str]:
        """Every way the drawn/undrawn split fails to partition the tree's components.

        Returned rather than raised so a caller can render the failure into the figure instead of
        crashing on it: a picture that says "this owner was neither drawn nor explained" is more use
        to the reader than a traceback.
        """
        components = {node.name for node in self.tree.components}
        problems: list[str] = []
        missing = sorted(components - self.drawn_owners - self.undrawn_owners)
        for name in missing:
            problems.append(f"{name}: neither drawn nor listed undrawn — silently omitted")
        for name in sorted(self.drawn_owners & self.undrawn_owners):
            problems.append(f"{name}: both drawn and listed undrawn — the split is not a partition")
        for name in sorted(self.drawn_owners - components):
            problems.append(f"{name}: drawn but is not a component of this tree")
        return problems


def _as_pairs(topology: Any) -> np.ndarray | None:
    """Read an ``(M, 2)``-ish index array, or return ``None`` if it is not one.

    Tolerant about dtype and strict about meaning: a segment array is a list of node-index pairs
    however it arrived, and anything that is not a pair of indices is not a segment and is reported.
    """
    if topology is None:
        return None
    try:
        array = np.asarray(topology)
    except Exception:  # noqa: BLE001 - any conversion failure means "not readable", which we report
        return None
    if array.ndim != 2 or array.shape[1] < 2 or array.shape[0] == 0:
        return None
    if not np.issubdtype(array.dtype, np.number):
        return None
    return array[:, :2].astype(np.int64)


def owner_elements_from(name: str, source: Any, colour: str = "#546e7a") -> OwnerElements:
    """Read an owner object's own arrays by duck-typing the names this package already uses."""
    return OwnerElements(
        name=name,
        positions=getattr(source, "positions", None),
        segments=getattr(source, "segments", None),
        crosslinks=getattr(source, "crosslinks", None),
        clutches=getattr(source, "clutches", None),
        colour=colour,
    )


def draw_scene(
    tree: SceneTree,
    owners: Mapping[str, OwnerElements] | Iterable[OwnerElements] = (),
    emergence: Mapping[str, Any] | None = None,
) -> DrawnScene:
    """Collect every drawable element, naming every component that could not be drawn.

    **Every component of ``tree`` ends up in exactly one of two places.** An owner with no entry in
    ``owners``, with no positions, with malformed positions, or with no readable topology is appended
    to :attr:`DrawnScene.undrawn` with a reason saying which. None of those is a silent skip, and none
    of them is an error either -- a field owner genuinely has no nodes.

    Args:
        tree: The scene tree, built from the runtime's own census and adjoint pairs.
        owners: Per-owner arrays, keyed by owner name or as an iterable of :class:`OwnerElements`.
        emergence: Optional per-owner observation objects to annotate the figure with.

    Returns:
        The :class:`DrawnScene`.
    """
    if isinstance(owners, Mapping):
        supplied = {str(key): value for key, value in owners.items()}
    else:
        supplied = {element.name: element for element in owners}

    scene = DrawnScene(tree=tree, emergence=dict(emergence or {}))

    for node in tree.components:
        source = supplied.get(node.name)
        if source is None:
            # If the tree already carries a reason for this component, use it. The tree was built by
            # whoever read the artefact and knows *why* -- "array of shape (5,) is not (N, 3)" is
            # worth more to the reader than this function's generic "nothing was supplied".
            scene.undrawn.append(
                (
                    node.name,
                    node.isolation_reason.strip()
                    or "no arrays were supplied for this owner, so nothing could be read",
                )
            )
            continue
        if source.positions is None:
            scene.undrawn.append(
                (node.name, "no `positions`: not a spatially resolved owner (a field or a ledger)")
            )
            continue
        try:
            nodes = np.asarray(source.positions, dtype=np.float64)
        except (TypeError, ValueError):
            scene.undrawn.append((node.name, "`positions` could not be read as a numeric array"))
            continue
        if nodes.ndim != 2 or nodes.shape[1] != 3:
            scene.undrawn.append(
                (node.name, f"`positions` has shape {nodes.shape}, expected (N, 3) in um")
            )
            continue
        if nodes.shape[0] == 0:
            scene.undrawn.append((node.name, "`positions` is empty: the owner holds no nodes"))
            continue
        if not np.isfinite(nodes).all():
            scene.undrawn.append((node.name, "`positions` carries a non-finite value"))
            continue

        before = len(scene.elements)
        for kind, topology in (
            ("segment", source.segments),
            ("crosslink", source.crosslinks),
        ):
            pairs = _as_pairs(topology)
            if pairs is None:
                continue
            valid = (pairs >= 0).all(axis=1) & (pairs < nodes.shape[0]).all(axis=1)
            for a, b in pairs[valid]:
                scene.elements.append(
                    DrawnElement(node.name, kind, nodes[a], nodes[b], source.colour)
                )
            if not valid.all():
                scene.undrawn.append(
                    (
                        node.name,
                        f"{int((~valid).sum())} {kind}(s) index a node this owner does not have, "
                        "so they were not drawn",
                    )
                )

        clutches = source.clutches
        if clutches is not None:
            indices = np.asarray(clutches).astype(np.int64).ravel()
            for index in indices[(indices >= 0) & (indices < nodes.shape[0])]:
                scene.elements.append(
                    DrawnElement(node.name, "clutch", nodes[index], None, source.colour)
                )

        if len(scene.elements) == before:
            # Positions were readable but no topology was: draw the nodes themselves rather than
            # nothing. The owner is still drawn from its own arrays, which is the requirement.
            for point in nodes:
                scene.elements.append(DrawnElement(node.name, "node", point, None, source.colour))

    return scene


def _project(
    points: np.ndarray, axes: tuple[int, int, int]
) -> tuple[np.ndarray, np.ndarray]:
    """Drop one world axis; return the kept pair and the depth coordinate."""
    horizontal, vertical, depth = axes
    return points[:, [horizontal, vertical]], points[:, depth]


def render_scene_svg(
    scene: DrawnScene,
    *,
    width: int = 1180,
    height: int = 760,
    axes: str = "xy",
    title: str = "",
) -> str:
    """Render a :class:`DrawnScene` to standalone SVG text.

    The figure carries four blocks, and the last two are the ones that make it evidence rather than
    decoration: the elements, the connector load path in its isolation styles, the **undrawn list**,
    and the element counts.

    Args:
        scene: What to draw.
        width: Canvas width in px.
        height: Canvas height in px.
        axes: Which world plane to keep -- one of :data:`AXIS_CHOICES`.
        title: Figure title.

    Returns:
        SVG text, newline-terminated.
    """
    if axes not in AXIS_CHOICES:
        raise ValueError(f"axes must be one of {sorted(AXIS_CHOICES)}; got {axes!r}")
    chosen = AXIS_CHOICES[axes]

    margin = 46.0
    panel = width * 0.66
    out: list[str] = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
        f"viewBox='0 0 {width} {height}' font-family='ui-monospace,monospace' font-size='11'>",
        f"<rect width='{width}' height='{height}' fill='#fbfbfd'/>",
    ]
    heading = title or f"Aleph scene — accepted step {scene.tree.accepted_step}"
    out.append(
        f"<text x='{margin}' y='26' font-size='15' fill='#111'>{html.escape(heading)}</text>"
    )
    out.append(
        f"<text x='{margin}' y='42' fill='#555'>topology epoch {scene.tree.topology_epoch} · "
        f"{len(scene.tree.components)} components · {len(scene.tree.connectors)} connectors · "
        f"{len(scene.elements)} elements drawn</text>"
    )

    # ---- the elements -------------------------------------------------------------------------
    owner_centre: dict[str, tuple[float, float]] = {}
    if scene.elements:
        starts = np.array([element.start for element in scene.elements], dtype=np.float64)
        ends = np.array(
            [
                element.end if element.end is not None else element.start
                for element in scene.elements
            ],
            dtype=np.float64,
        )
        all_points = np.concatenate([starts, ends], axis=0)
        flat, depth = _project(all_points, chosen)
        low = flat.min(axis=0)
        high = flat.max(axis=0)
        span = np.maximum(high - low, 1.0e-9)
        scale = min((panel - 2 * margin) / span[0], (height - 2 * margin - 40) / span[1])

        def place(point: np.ndarray) -> tuple[float, float]:
            horizontal = float(point[chosen[0]])
            vertical = float(point[chosen[1]])
            x = margin + (horizontal - low[0]) * scale
            # SVG's y grows downward; flip so +y is up, which is what a reader expects of um.
            y = height - margin - (vertical - low[1]) * scale
            return float(x), float(y)

        depth_low = float(depth.min())
        depth_span = max(float(depth.max()) - depth_low, 1.0e-9)

        accumulator: dict[str, list[tuple[float, float]]] = {}
        # Points last, so a clutch is never hidden under a filament it sits on.
        ordered = sorted(scene.elements, key=lambda e: (e.end is None, e.owner))
        for element in ordered:
            x1, y1 = place(element.start)
            accumulator.setdefault(element.owner, []).append((x1, y1))
            relative = (float(element.start[chosen[2]]) - depth_low) / depth_span
            opacity = 0.42 + 0.5 * relative
            colour = html.escape(element.colour)
            if element.end is None:
                radius = 3.0 if element.kind == "clutch" else 1.5
                out.append(
                    f"<circle cx='{x1:.2f}' cy='{y1:.2f}' r='{radius}' fill='{colour}' "
                    f"fill-opacity='{opacity:.2f}'/>"
                )
            else:
                x2, y2 = place(element.end)
                accumulator[element.owner].append((x2, y2))
                stroke_width = 1.7 if element.kind == "segment" else 0.9
                dash = " stroke-dasharray='3 2'" if element.kind == "crosslink" else ""
                out.append(
                    f"<line x1='{x1:.2f}' y1='{y1:.2f}' x2='{x2:.2f}' y2='{y2:.2f}' "
                    f"stroke='{colour}' stroke-width='{stroke_width}' "
                    f"stroke-opacity='{opacity:.2f}'{dash}/>"
                )
        for name, coordinates in accumulator.items():
            array = np.array(coordinates, dtype=np.float64)
            owner_centre[name] = (float(array[:, 0].mean()), float(array[:, 1].mean()))
    else:
        out.append(
            f"<text x='{margin}' y='{height / 2:.0f}' fill='#c62828'>no elements drawn — "
            "every owner is listed undrawn below</text>"
        )

    # ---- the connector load path, in isolation styles ------------------------------------------
    fallback_y = 70.0
    for index, node in enumerate(scene.tree.connectors):
        style = ISOLATION_STYLE[node.isolation]
        a = owner_centre.get(node.endpoint_a)
        b = owner_centre.get(node.endpoint_b)
        if a is None or b is None:
            # An endpoint whose owner was not drawn still gets its connector shown, on a stub rail.
            # Dropping it would hide exactly the NOT_SCHEDULED case, which is the point of the module.
            a = (panel + 30.0, fallback_y + 18.0 * index)
            b = (panel + 110.0, fallback_y + 18.0 * index)
        mid = ((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)
        dash = f" stroke-dasharray='{style.dash}'" if style.dash else ""
        out.append(
            f"<line x1='{a[0]:.2f}' y1='{a[1]:.2f}' x2='{b[0]:.2f}' y2='{b[1]:.2f}' "
            f"stroke='{style.stroke}' stroke-width='{style.width}' "
            f"stroke-opacity='{style.opacity:.2f}'{dash}/>"
        )
        if style.marker:
            out.append(
                f"<text x='{mid[0]:.2f}' y='{mid[1]:.2f}' fill='{style.stroke}' font-size='12' "
                f"text-anchor='middle'>{html.escape(style.marker)}</text>"
            )
        out.append(
            f"<title>{html.escape(node.name)} — {html.escape(node.isolation.value)} — "
            f"{html.escape(node.isolation_reason)}</title>"
        )

    # ---- the side panel: legend, undrawn list, counts, emergence -------------------------------
    x = panel + 24.0
    y = 70.0
    out.append(f"<text x='{x}' y='{y}' font-size='12' fill='#111'>isolation states</text>")
    y += 16.0
    for state in IsolationState:
        style = ISOLATION_STYLE[state]
        dash = f" stroke-dasharray='{style.dash}'" if style.dash else ""
        out.append(
            f"<line x1='{x}' y1='{y - 4:.0f}' x2='{x + 26}' y2='{y - 4:.0f}' "
            f"stroke='{style.stroke}' stroke-width='{style.width}'{dash}/>"
        )
        if style.marker:
            out.append(
                f"<text x='{x + 13}' y='{y - 6:.0f}' fill='{style.stroke}' font-size='10' "
                f"text-anchor='middle'>{html.escape(style.marker)}</text>"
            )
        weight = " font-weight='bold'" if style.is_fault else ""
        colour = "#c62828" if style.is_fault else "#333"
        out.append(
            f"<text x='{x + 34}' y='{y:.0f}' fill='{colour}'{weight}>"
            f"{html.escape(style.label)}</text>"
        )
        y += 15.0

    faults = scene.tree.faults
    y += 8.0
    out.append(
        f"<text x='{x}' y='{y:.0f}' font-size='12' fill="
        f"'{'#c62828' if faults else '#111'}'>faults: {len(faults)}</text>"
    )
    y += 15.0
    for entry in faults[:6]:
        out.append(f"<text x='{x + 8}' y='{y:.0f}' fill='#c62828'>{html.escape(entry[:64])}</text>")
        y += 13.0

    y += 8.0
    problems = scene.coverage_problems()
    out.append(
        f"<text x='{x}' y='{y:.0f}' font-size='12' fill='#111'>undrawn owners: "
        f"{len(scene.undrawn)}</text>"
    )
    y += 15.0
    for name, reason in scene.undrawn[:8]:
        out.append(
            f"<text x='{x + 8}' y='{y:.0f}' fill='#6d4c41'>{html.escape(name)}: "
            f"{html.escape(reason[:58])}</text>"
        )
        y += 13.0
    if problems:
        y += 4.0
        out.append(
            f"<text x='{x}' y='{y:.0f}' font-size='12' fill='#c62828' font-weight='bold'>"
            f"COVERAGE DEFECT: {len(problems)}</text>"
        )
        y += 14.0
        for problem in problems[:4]:
            out.append(
                f"<text x='{x + 8}' y='{y:.0f}' fill='#c62828'>{html.escape(problem[:60])}</text>"
            )
            y += 13.0

    counts = scene.counts()
    if counts:
        y += 10.0
        out.append(f"<text x='{x}' y='{y:.0f}' font-size='12' fill='#111'>elements</text>")
        y += 15.0
        for key, value in sorted(counts.items())[:10]:
            out.append(
                f"<text x='{x + 8}' y='{y:.0f}' fill='#333'>{html.escape(key)}: {value}</text>"
            )
            y += 13.0

    if scene.emergence:
        y += 10.0
        out.append(f"<text x='{x}' y='{y:.0f}' font-size='12' fill='#111'>emergence</text>")
        y += 15.0
        for name, observation in sorted(scene.emergence.items()):
            block = (
                observation.as_json_obj()
                if hasattr(observation, "as_json_obj")
                else dict(observation)
            )
            verdict = "ORDERED" if block.get("is_ordered") else "not ordered"
            out.append(
                f"<text x='{x + 8}' y='{y:.0f}' fill='#333'>{html.escape(name)}: S="
                f"{block.get('order', float('nan')):.3f} z={block.get('z', float('nan')):.1f} "
                f"{verdict}, {block.get('n_bundles', 0)} bundle(s)</text>"
            )
            y += 13.0

    out.append("</svg>")
    return "\n".join(out) + "\n"
