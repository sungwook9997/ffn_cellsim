#!/usr/bin/env python3
"""Native (non-browser) debug viewer for an `ac` state dump — the fast inner loop.

The WebGL/HTML viewer (`ff_viewer_html.build_viewer`) is the SHARE artifact: self-contained,
sendable, browser-verified, and the thing PI opens. It is not a good *debugging* loop — a
native-population scene is a ~100 MB inline payload that must be written, opened, and
browser-checked before you can rotate it once.

This module is the other half: the same `CellDump` data layer (`ac_viz_common.load_dump`)
rendered straight into a native OpenGL window, so "is the cortex interpenetrating? is that
force field the right sign?" is answered in seconds. It renders the SAME arrays the HTML
viewer does — full native resolution, no downsampling (PI 2026-07-07) — so what you inspect
here is what ships there.

Scope: a debug/inspection tool. It is NOT a milestone artifact — the per-stage visualization
gate (PI 2026-07-22) is still satisfied by the browser-verified HTML, never by a screenshot
from this script.

Usage:
    python -m aleph.scripts.ac_polyscope_view --dump outputs/ac/cell_assembled/assembled_state.npz
    python -m aleph.scripts.ac_polyscope_view --dump <npz> --screenshot /tmp/scene.png   # headless-ish
    python -m aleph.scripts.ac_polyscope_view --dump <npz> --no-connectors --force-field
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from aleph.scripts.ac_viz_common import (
    CellDump,
    Component,
    Connector,
    component_color,
    connector_color,
    load_dump,
    seg_pairs_from_offsets,
)


def hex_to_rgb(color: str) -> tuple[float, float, float]:
    """Convert a ``#rrggbb`` string to a polyscope ``(r, g, b)`` float triple in [0, 1].

    Args:
        color: Hex colour string, with or without the leading ``#``.

    Returns:
        The colour as three floats in [0, 1]. Unparseable input falls back to mid grey.
    """
    s = color.lstrip("#")
    if len(s) != 6:
        return (0.6, 0.65, 0.75)
    try:
        return tuple(int(s[i:i + 2], 16) / 255.0 for i in (0, 2, 4))  # type: ignore[return-value]
    except ValueError:
        return (0.6, 0.65, 0.75)


def _component_edges(c: Component) -> np.ndarray | None:
    """Return the (M, 2) node-index pairs that draw a filament component, or None.

    Prefers the component's explicit ``seg``; falls back to per-filament offsets in ``meta``
    (the same convention `ac_viz_common` uses to rebuild segments for the HTML layers).
    """
    if c.seg is not None and len(c.seg):
        return np.asarray(c.seg, dtype=np.int64).reshape(-1, 2)
    offsets = c.meta.get("offsets")
    if offsets is not None:
        seg = seg_pairs_from_offsets(np.asarray(offsets), c.n_nodes)
        return np.asarray(seg, dtype=np.int64).reshape(-1, 2) if len(seg) else None
    return None


def register_component(ps, c: Component, *, force_field: bool) -> None:
    """Register one component with polyscope as a curve network, mesh, or point cloud.

    Args:
        ps: The imported ``polyscope`` module.
        c: The component to draw, at full native resolution.
        force_field: When True and the component carries per-node ``fmag``, attach it as a
            scalar quantity (log10 |F| in pN) and enable it, matching the HTML viewer's
            force colouring.
    """
    pos = np.ascontiguousarray(c.pos, dtype=np.float64).reshape(-1, 3)
    if pos.size == 0:
        return
    rgb = hex_to_rgb(c.color or component_color(c.key))
    name = f"{c.key} · {c.label}"[:120]

    struct = None
    if c.kind == "mesh" and c.faces is not None and len(c.faces):
        struct = ps.register_surface_mesh(name, pos, np.asarray(c.faces, dtype=np.int64), smooth_shade=True)
    else:
        edges = _component_edges(c)
        if edges is not None and len(edges):
            struct = ps.register_curve_network(name, pos, edges, radius=0.0012)
        else:
            struct = ps.register_point_cloud(name, pos, radius=0.0016)
    struct.set_color(rgb)

    if force_field and c.fmag is not None and len(c.fmag) == pos.shape[0]:
        f = np.asarray(c.fmag, dtype=np.float64)
        struct.add_scalar_quantity("log10 |F| [pN]", np.log10(np.maximum(f, 1e-6)),
                                   enabled=True, cmap="turbo")


def register_connector(ps, cn: Connector) -> None:
    """Register one connector family as a curve network of its explicit joints.

    Connectors are the ONLY mechanical connection in the engine graph (co-location is never a
    connection), so they are drawn as their own structures rather than folded into a component.
    """
    ep = np.ascontiguousarray(cn.endpoints, dtype=np.float64).reshape(-1, 3)
    if ep.shape[0] < 2:
        return
    n_pairs = ep.shape[0] // 2
    edges = np.arange(n_pairs * 2, dtype=np.int64).reshape(-1, 2)
    net = ps.register_curve_network(f"⇄ {cn.family} · {cn.label}"[:120], ep[:n_pairs * 2], edges,
                                    radius=0.0018)
    net.set_color(hex_to_rgb(cn.color or connector_color(cn.family)))
    if cn.load is not None and len(cn.load) == n_pairs:
        net.add_scalar_quantity("load [pN]", np.asarray(cn.load, dtype=np.float64),
                                defined_on="edges", enabled=True, cmap="turbo")


def build_scene(dump: CellDump, *, connectors: bool = True, force_field: bool = False) -> tuple[int, int]:
    """Register every component (and optionally every connector) of a dump with polyscope.

    Args:
        dump: The normalized cell dump to draw.
        connectors: Whether to draw connector families as separate structures.
        force_field: Whether to attach per-node |F| scalar quantities.

    Returns:
        ``(n_components, n_connectors)`` actually registered.
    """
    import polyscope as ps

    n_c = 0
    for c in dump.components:
        register_component(ps, c, force_field=force_field)
        n_c += 1
    n_k = 0
    if connectors:
        for cn in dump.connectors:
            register_connector(ps, cn)
            n_k += 1
    return n_c, n_k


def describe(dump: CellDump) -> str:
    """Return a one-line population summary (nodes and actors are reported separately)."""
    nodes = sum(c.n_nodes for c in dump.components)
    actors = sum(c.n_actors for c in dump.components)
    joints = sum(cn.n for cn in dump.connectors)
    return (f"{len(dump.components)} components · {nodes:,} nodes · {actors:,} actors | "
            f"{len(dump.connectors)} connector families · {joints:,} joints")


def main() -> None:
    """CLI entry point."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--dump", required=True, type=Path, help="state dump .npz (v2 or legacy)")
    ap.add_argument("--screenshot", type=Path, default=None,
                    help="render one frame to this PNG instead of opening the interactive window")
    ap.add_argument("--no-connectors", action="store_true", help="components only")
    ap.add_argument("--force-field", action="store_true",
                    help="attach per-node log10|F| scalar quantities and enable them")
    args = ap.parse_args()

    dump = load_dump(args.dump)
    print(f"[polyscope] {args.dump}")
    print(f"[polyscope] {describe(dump)}")

    import polyscope as ps

    ps.set_program_name("ffn_cellsim — ac state")
    ps.set_up_dir("z_up")
    ps.set_ground_plane_mode("shadow_only")
    ps.init()

    n_c, n_k = build_scene(dump, connectors=not args.no_connectors, force_field=args.force_field)
    print(f"[polyscope] registered {n_c} components, {n_k} connectors (full native resolution)")

    if args.screenshot is not None:
        args.screenshot.parent.mkdir(parents=True, exist_ok=True)
        ps.set_screenshot_extension(".png")
        ps.screenshot(str(args.screenshot))
        print(f"[polyscope] wrote {args.screenshot}")
        return
    ps.show()


if __name__ == "__main__":
    main()
