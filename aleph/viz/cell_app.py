r"""A native viewer for the whole cell — every node, every segment, no page and no thinning.

**Why this is not a web page.** The HTML viewer works and is the right tool at a slice, but a page of
eight populations is 78–124 MB of markup and the browser parses every line before it draws anything.
The cell is not a large amount of DATA — 4.6 M nodes and 4.3 M segments is 94 MB as raw float32/int32
— it is only large as text. Above about a million lines a page stops being slow and starts being
unopenable, and the standing instruction (`feedback-viewer-no-downsample`) forbids the obvious escape.

So this reads the binary :mod:`aleph.scripts.world_export_cell` writes and hands it to the GPU. The
segments go up once as a vertex buffer and stay there; orbiting re-uses them. **That is the whole
trick** — the geometry never crosses back to the host, which is the same principle the arena itself is
built on.

**Fat lines, never dots.** `feedback-viewer-fat-lines-no-dots` is a standing instruction and it is
physical, not aesthetic: a filament has a real thickness and a point cloud silently discards both the
connectivity and the width. Warp's ``render_line_list`` draws capsules of a stated radius, so a
cortex filament is drawn at ITS radius and a stress fibre at ITS OWN, and the two do not look alike
because they are not.

**Surfaces are meshes, wire cages are not surfaces.** The membrane and the nuclear envelope arrive as
faces and are drawn as faces, so an interior structure is genuinely occluded rather than showing
through a lattice. The exporter writes them as ``faces`` blocks for exactly this reason.

⚠ **This viewer draws geometry and claims nothing.** No force is evaluated, nothing moves, and a
picture is not evidence. Its own status line says so, because a screenshot of a cell is the single
most quotable artifact this project produces and the caption has to travel with it.

Usage::

    python -m aleph.viz.cell_app aleph/outputs/ac/world_phase1/cell_full.alephcell
    python -m aleph.viz.cell_app cell_full.alephcell --only membrane,microtubule --headless shot.png

Controls are Warp's: drag to orbit, scroll to zoom, ``W``/``A``/``S``/``D`` to fly, ``X`` to
screenshot, ``ESC`` to quit.

Sanity Gate:
    * dimensions — µm throughout; radii are µm and are stated per population below.
    * boundary cases — a missing file, a wrong magic, a truncated payload and an empty ``--only``
      selection each fail with what was wrong, not with a traceback into a parser.
    * conservation — the loader re-checks every segment index against its own block, so a file that
      would draw a line between unrelated structures is refused HERE as well as at export.
    * precision — float32 positions; the header carries the error bound and it is printed.
    * measurement protocol — the header's population counts are printed next to what was drawn, so a
      selection that silently dropped something is visible in the same breath.
"""

from __future__ import annotations

import argparse
import json
import struct
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

MAGIC = b"ALEPHCEL"

#: Draw radius per population [µm], and a word on where each comes from. **These are DRAW radii, not
#: claims about the structure.** Where a real thickness is known it is used; where it is not, the
#: value is chosen so the population is legible beside its neighbours and is marked ``legibility``.
#: A viewer that pretended these were measurements would be doing what this project spends its time
#: catching, so the table says which is which and :func:`describe_radii` prints it on request.
DRAW_RADIUS_UM: dict[str, tuple[float, str]] = {
    "cortex": (0.0035, "actin filament radius ~3.5 nm — a real thickness"),
    "microtubule": (0.0125, "MT outer radius ~12.5 nm — a real thickness"),
    "intermediate_filament": (0.005, "IF radius ~5 nm — a real thickness"),
    "filopodium": (0.004, "legibility: a bundle drawn at single-filament scale"),
    "lamellipodium": (0.004, "legibility"),
    "stress_fiber": (0.006, "legibility: an SF bundle is far thicker than one filament"),
    "nmii": (0.008, "legibility: a minifilament backbone drawn thick enough to find"),
    "sf_arc": (0.006, "legibility: a transverse-arc bundle, as for stress_fiber"),
    "lamina": (0.004, "legibility"),
    "chromatin": (0.010, "legibility"),
}

#: Colour per population. Chosen for separability against a light background, not for meaning.
COLOUR: dict[str, tuple[float, float, float]] = {
    "cortex": (0.94, 0.45, 0.33),
    "membrane": (0.98, 0.62, 0.48),
    "nuclear_envelope": (0.85, 0.55, 0.80),
    "microtubule": (0.35, 0.55, 0.85),
    "intermediate_filament": (0.45, 0.75, 0.55),
    "filopodium": (0.20, 0.70, 0.68),
    "lamellipodium": (0.55, 0.75, 0.30),
    "stress_fiber": (0.60, 0.80, 0.25),
    "nmii": (0.15, 0.75, 0.85),
    "sf_arc": (0.90, 0.70, 0.20),
    "lamina": (0.70, 0.45, 0.85),
    "chromatin": (0.55, 0.35, 0.70),
}


#: Named camera positions, because the default one hides the half of the cell that matters most.
#: The basal structures — ventral stress fibres, transverse arcs, the lamellipodium — all sit at
#: ``z = -7 µm``, so an axis-on view puts every one of them behind the nucleus. That is not a
#: rendering problem; it is a fact about where they are, and a viewer whose only camera hides them
#: teaches the wrong shape. ``basal`` looks up at the contact disc from below, which is the view an
#: adherent cell is normally imaged in.
VIEWS: dict[str, tuple[tuple[float, float, float], tuple[float, float, float]]] = {
    "front":  ((0.0, 4.0, 34.0), (0.0, -0.1, -1.0)),
    "side":   ((34.0, 4.0, 0.0), (-1.0, -0.1, 0.0)),
    "top":    ((0.0, 34.0, 0.1), (0.0, -1.0, 0.0)),
    "basal":  ((0.0, -26.0, 20.0), (0.0, 0.72, -0.70)),
    "oblique": ((22.0, -12.0, 22.0), (-0.62, 0.34, -0.62)),
}


def camera_is_on_the_kept_side(view: str, cut: tuple[int, float, int] | None) -> float | None:
    """How far the camera sits INTO the half a cut keeps, in µm. ``None`` when there is no cut.

    ⚠ **A cut only opens the cell if the camera is on the side the cut REMOVED.** Otherwise the
    retained half presents its closed outer surface and the cut is invisible — the render looks like
    an uncut cell and the element counts in the status line are the only sign anything was removed.

    This was found by predicting the opposite and being wrong. The first diagnosis was that `--cut`
    and the camera are unrelated, so `z-` hid the interior; it does not — `oblique` sits at
    ``z = +22`` and `--cut z-` keeps ``z <= 0``, so the camera is already on the removed side and the
    cut face IS toward it. What hides the interior there is **the cortex's own inner surface**: 4.2 M
    filaments on a shell are opaque from inside as well as outside, which this module's own header
    already said. `--cut z+` with the same camera is the case that is actually wrong, and it is the
    one this function computes.

    Returns:
        Signed distance [µm] from the camera to the cut plane, positive when the camera is inside the
        KEPT half — i.e. when the closed side faces it. ``None`` if ``cut`` is ``None``.
    """
    if cut is None:
        return None
    axis, offset, sign = cut
    return (VIEWS[view][0][axis] - offset) * float(sign)


@dataclass(frozen=True, slots=True)
class CellFile:
    """A loaded export: the header, and each population's arrays."""

    header: dict
    positions: dict[str, np.ndarray]
    segments: dict[str, np.ndarray]
    faces: dict[str, np.ndarray]
    #: frame index -> population -> positions. A static export has exactly one frame, ``0``, and it
    #: is the same array as :attr:`positions` — so nothing about a single-frame file changed when
    #: sequences arrived.
    frames: dict[int, dict[str, np.ndarray]] = field(default_factory=dict)
    #: Where it was read from. Carried so the app can find the sibling run record; ``None`` for a
    #: CellFile built in a test. ⚠ It is NOT optional in practice — `getattr(cell, "path", None)`
    #: was the first attempt and silently returned None forever, which is the app doing nothing
    #: while looking like it works. The field exists so the absence has to be constructed.
    path: Path | None = None

    @property
    def populations(self) -> list[str]:
        return list(self.positions)

    def summary(self) -> str:
        n = sum(int(p.shape[0]) for p in self.positions.values())
        s = sum(int(x.shape[0]) for x in self.segments.values())
        f = sum(int(x.shape[0]) for x in self.faces.values())
        base = (f"{n:,} nodes · {s:,} segments · {f:,} faces · "
                f"{len(self.positions)} populations")
        return base if len(self.frames) <= 1 else f"{base} · {len(self.frames)} FRAMES"


def load(path: Path) -> CellFile:
    """Read an ``.alephcell`` file, re-checking what the exporter checked.

    Raises:
        FileNotFoundError: if the path does not exist.
        ValueError: on a wrong magic, a truncated payload, or a segment index that leaves its own
            population's block. The last one is re-checked here on purpose: a file is a thing that
            gets copied, and a viewer that trusts its input to have been validated once is a viewer
            that will one day draw a line between two unrelated structures.
    """
    raw = path.read_bytes()
    if raw[:8] != MAGIC:
        raise ValueError(f"{path} is not an aleph cell export (magic {raw[:8]!r}, expected {MAGIC!r})")
    (hlen,) = struct.unpack("<I", raw[8:12])
    header = json.loads(raw[12:12 + hlen].decode("utf-8"))
    body = raw[12 + hlen:]

    positions: dict[str, np.ndarray] = {}
    segments: dict[str, np.ndarray] = {}
    faces: dict[str, np.ndarray] = {}
    # frame -> population -> positions. Frame 0 IS `positions`, so a single-frame file behaves exactly
    # as it did before this existed and an old file still loads.
    frames: dict[int, dict[str, np.ndarray]] = {}
    for blk in header["blocks"]:
        off, nb = int(blk["offset"]), int(blk["nbytes"])
        if off + nb > len(body):
            raise ValueError(
                f"{path} is truncated: block {blk['population']}/{blk['kind']} wants bytes "
                f"[{off}, {off + nb}) and the payload is {len(body)}."
            )
        arr = np.frombuffer(body, dtype=np.dtype(blk["dtype"]), count=nb // np.dtype(blk["dtype"]).itemsize,
                            offset=off).reshape(blk["shape"])
        kind = blk["kind"]
        if kind == "positions_um":
            fr = int(blk.get("frame", 0))
            frames.setdefault(fr, {})[blk["population"]] = arr
            if fr == 0:
                positions[blk["population"]] = arr
        else:
            {"segments": segments, "faces": faces}[kind][blk["population"]] = arr

    for name, idx in list(segments.items()) + list(faces.items()):
        n = int(positions[name].shape[0])
        if idx.size and (int(idx.min()) < 0 or int(idx.max()) >= n):
            raise ValueError(
                f"{path}: {name} has an index outside its own block "
                f"([{idx.min()}, {idx.max()}] against {n} nodes). Refusing rather than drawing a line "
                "to an unrelated structure."
            )
    return CellFile(header=header, positions=positions, segments=segments, faces=faces,
                    frames=frames, path=path)


def describe_radii() -> str:
    """The draw-radius table, with which entries are real thicknesses and which are legibility."""
    rows = [f"  {k:<24} {v:.4f} µm   {why}" for k, (v, why) in DRAW_RADIUS_UM.items()]
    return ("draw radii — NOT claims about the structure:\n" + "\n".join(rows))


#: Axis letter -> column of the position array. Y is up in this renderer.
_CUT_AXIS = {"x": 0, "y": 1, "z": 2}


def parse_cut(spec: str) -> tuple[int, float, int]:
    """Parse ``AXIS[+-][:OFFSET]`` into ``(axis index, offset um, sign)``.

    ⚠ A malformed cut RAISES, never silently yields None. A viewer that quietly drops a cut draws the
    whole cell and prints a summary of a picture the caller did not ask for — the same shape as the
    empty index list this file's history is about.

    ⚠ **Shared rather than copied.** It was inline in :func:`main` until `world_turntable.py` needed
    the same grammar; two parsers for one spelling is how the two placement copies in `world/` drifted
    on 2026-08-20, and a cut that means different halves in two tools is worse than either.

    Raises:
        ValueError: on a bad axis, a missing or bad sign, or a non-numeric offset.
    """
    text = spec.strip().lower()
    axis_c, rest = (text[:1], text[1:])
    if axis_c not in _CUT_AXIS or not rest or rest[0] not in "+-":
        raise ValueError(f"--cut {spec!r} is not AXIS[+-][:OFFSET], e.g. 'y-' or 'z+:1.5'. axis must "
                         f"be one of {sorted(_CUT_AXIS)} and the sign says which half to KEEP.")
    tail = rest[1:].lstrip(":")
    try:
        offset = float(tail) if tail else 0.0
    except ValueError:
        raise ValueError(f"--cut offset {tail!r} is not a number.") from None
    return (_CUT_AXIS[axis_c], offset, 1 if rest[0] == "+" else -1)


def _apply_cut(pos: np.ndarray, idx: np.ndarray, cut: tuple[int, float, int] | None
               ) -> tuple[np.ndarray, int]:
    """Keep only the elements lying wholly on one side of a plane. Returns ``(kept, n_hidden)``.

    ⚠ **This is a SPATIAL SELECTION, never a thinning.** The distinction is the whole reason the flag
    exists: a thinned draw shows *some of* the structure everywhere and misrepresents density, which
    this project forbids; a cut shows *all of* the structure in half the volume and misrepresents
    nothing, provided the reader is told where the plane is. The caller prints the plane and the
    hidden count with every frame.

    An element is kept only when **every** vertex of it is on the kept side. Clipping an element
    partway would invent a filament end that the build does not have, and a viewer that invents ends
    is a viewer that can be quoted about lengths.

    Args:
        pos: ``(N, 3)`` positions [µm].
        idx: ``(M, k)`` element vertex indices — segment pairs or triangle/quad faces.
        cut: ``(axis, offset_um, sign)`` or ``None``. ``sign`` is +1 to keep the side above
            ``offset`` and -1 to keep the side below.

    Returns:
        The kept index rows, and how many rows were hidden.
    """
    if cut is None:
        return idx, 0
    axis, offset, sign = cut
    side = (pos[:, axis] - offset) * sign >= 0.0
    keep = side[idx].all(axis=1)
    return idx[keep], int((~keep).sum())


def _draw_population(renderer, name: str, pos: np.ndarray, cell: CellFile, thicken: float,
                     log: list[str], cut: tuple[int, float, int] | None = None) -> None:
    """Draw one population at the given positions. One implementation, used by both draw paths.

    Two paths drawing the same population differently is the shape of the defect this project spent
    2026-08-20 removing, so playback and the first frame share this rather than each having their own.

    ⚠ **That sentence was false when it was written.** This function existed and said it, and the
    static first-frame path had its own inline copy of the same fifteen lines regardless — found
    2026-08-21 while adding ``--cut``, which would otherwise have landed in one of the two paths and
    been absent from the other. The static path now calls this. **A docstring asserting a property
    the code does not have is the same defect as a gate reporting a verdict it did not reach.**
    """
    if name in cell.faces:
        faces = cell.faces[name].astype(np.int32)
        kept, hidden = _apply_cut(pos, faces, cut)
        renderer.render_mesh(name=name, points=pos.astype(np.float32),
                             indices=kept.reshape(-1),
                             colors=COLOUR.get(name, (0.6, 0.6, 0.6)))
        log.append(f"{name} (mesh, {kept.shape[0]:,} faces"
                   + (f", {hidden:,} cut away)" if hidden else ")"))
        return
    if name not in cell.segments:
        # ⚠ POSITIONS WITHOUT TOPOLOGY. `lamina` and `chromatin` are written as nodes with nothing
        # connecting them, because their builders record no n_strands/nodes_per_strand and the
        # exporter refuses to invent a connectivity the build does not have. A viewer could draw them
        # as points; this project does not draw points, because a dot cloud reads as a density map
        # and these are filaments and a chain.
        #
        # ⚠ SAID, not skipped. The first version of this viewer drew nothing for an empty index list
        # and printed a confident summary of what it had not drawn; silence here would be the same
        # defect wearing the word "unsupported".
        log.append(f"{name} ({pos.shape[0]:,} nodes, NO TOPOLOGY — not drawn)")
        return
    seg = cell.segments[name].astype(np.int32)
    kept, hidden = _apply_cut(pos, seg, cut)
    radius = DRAW_RADIUS_UM.get(name, (0.004, "legibility"))[0] * thicken
    colour = COLOUR.get(name, (0.6, 0.6, 0.6))
    # ⚠ NOT `render_line_list`, and the reason is measured. That function builds its capsule array
    # with a PYTHON LOOP over every segment — `for i in range(len(indices)//2): lines.append(...)` —
    # so issuing the cortex alone costs ~3.0 s, and this app re-issues EVERY population on EVERY
    # keystroke. A full 12-population render took 24.5 s before this line changed. `kept` is already
    # (n, 2), so `pos[kept]` is the same array in one fancy-index: 40x on the same data.
    # ponytail: calls warp's private `_render_lines`. The public wrapper does nothing else — it only
    # builds that array and forwards — so if warp renames it, the fallback below is byte-identical
    # in output and merely slow again.
    lines = pos.astype(np.float32)[kept]
    if hasattr(renderer, "_render_lines"):
        renderer._render_lines(name, lines, colour, radius)
    else:
        renderer.render_line_list(name=name, vertices=pos.astype(np.float32),
                                  indices=kept.reshape(-1), color=colour, radius=radius)
    log.append(f"{name} ({kept.shape[0]:,} segments @ {radius:.4f} um"
               + (f", {hidden:,} cut away)" if hidden else ")"))


def _write_png(path: Path, rgb: np.ndarray) -> None:
    """Write an RGB uint8 array as a PNG, with no image library.

    Deliberately dependency-free: this project already installs one package for the window and adding
    a second for a screenshot would be paying twice for something the standard library does. PNG's
    scanline format is a filter byte then the row, zlib-compressed — about ten lines.
    """
    import binascii
    import zlib

    h, w, _ = rgb.shape
    raw = b"".join(b"\x00" + rgb[y].tobytes() for y in range(h))

    def chunk(tag: bytes, data: bytes) -> bytes:
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", binascii.crc32(tag + data) & 0xFFFFFFFF))

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(raw, 6))
        + chunk(b"IEND", b""))


#: Lines the HUD paints, and the ones it paints in warning red. A ⚠ line is not decoration: it is
#: the sentence that stops the picture being quoted for more than it shows.
_HUD_WARN = "\u26a0"


def _caption(cell: CellFile, view: str, cut, into) -> list[str]:
    """The lines the PICTURE carries, top first.

    ⚠ Built from the SAME values the terminal block prints — `cell.summary()`, `view`, `cut`, `into`
    and the header's own `thinning` — so no number is retyped and the two surfaces cannot drift on a
    fact. They differ only in LENGTH, deliberately: the terminal carries the full argument and the
    picture carries the sentence that has to travel with it, because the picture is what leaves.
    """
    lines = [
        f"{cell.path.name if cell.path is not None else 'cell'} — {cell.summary()}",
        f"view {view}"
        + (f"  ·  cut: keeping {'xyz'[cut[0]]} {'>=' if cut[2] > 0 else '<='} {cut[1]:g} um"
           if cut is not None else "  ·  no cut"),
        f"thinning: {cell.header.get('thinning')} — a cut is a SPATIAL SELECTION, never a thinning",
        "⚠ geometry only — no force is evaluated, nothing has moved, and a picture is not evidence",
    ]
    if into is not None and into > 0.0:
        lines.insert(2, f"⚠ THE CUT FACES AWAY FROM THIS CAMERA ({into:.1f} um inside the kept "
                        f"half) — this looks like an UNCUT cell. Not auto-corrected.")
    return lines


def _attach_hud(renderer, get_lines, *, on_screen: bool = True) -> None:
    """Paint the app's own state ON the picture, not only into the terminal.

    ⚠ **WHY THIS EXISTS, and it is not ergonomics.** Until this, every word the app said about itself
    went to `print()` — which populations are on, how many elements the cut hid, `thinning: NONE`, the
    cut-facing-away warning, and the standing caveat that a picture is not evidence. The window showed
    geometry and nothing else, so **the screenshot `p` writes carried no caveat at all.** The PNG is
    the thing that leaves this machine and gets pasted into a report; the terminal scrollback is not.
    A render whose caveats live somewhere the render does not go is a render that will one day be
    shown without them.

    ⚠ **It also means the window can be read alone.** The comment this replaces said "the terminal IS
    this app's status bar", which is true of a tool you drive from a terminal and false of an app —
    an app tells you its state in itself, and the PI asked for an app.

    The hook is ``renderer._draw`` rather than a pyglet ``on_draw`` handler because pushed handlers
    dispatch BEFORE the renderer's own, and the HUD has to land on top of the scene. It is also the
    ONE seam both paths share: the interactive window reaches it through pyglet's redraw, and the
    headless PNG calls it directly — so the caption is on the screenshot by construction rather than
    by a second code path that could drift from this one.

    Args:
        renderer: A ``warp.render.OpenGLRenderer``. Its private ``_draw`` is wrapped in place.
        get_lines: Called every frame; returns the lines to paint, top line first. Called late so a
            live state line is current rather than a snapshot taken at attach time.
    """
    # ponytail: wraps a private `_draw`. warp exposes no public post-scene hook and its pyglet
    # `on_draw` dispatches handlers before the renderer's own. If warp renames it, the HUD vanishes
    # and the geometry still draws — so the failure is a missing caption, never a black window.
    import pyglet
    from pyglet import gl
    from pyglet.math import Mat4

    window = getattr(renderer, "window", None)
    original = renderer._draw
    if window is None or not isinstance(window, pyglet.window.Window):
        return                                    # no window to paint on; the prints still happen

    # ⚠ Labels are CREATED ONCE and reused. pyglet 2.1 tears a Label's glyph boxes down on garbage
    # collection, so building them per frame raised `'Label' object has no attribute '_boxes'` on the
    # second frame — found by rendering, not by reading. Reuse also means the caption costs nothing
    # per frame, which matters because this loop is already re-issuing 4.3 M segments.
    cache: list = []
    batch = pyglet.graphics.Batch()

    def _paint(lines: list[str], flip: bool, fbo) -> None:
        """Paint the caption once, into one target, in that target's own orientation."""
        w, h = window.width, window.height
        size = max(11, int(h / 78))
        pad = size
        while len(cache) < len(lines):
            cache.append(pyglet.text.Label("", x=pad, y=0, font_size=size, batch=batch))
        for k, label in enumerate(cache):
            if k < len(lines):
                label.text = lines[k]
                label.x, label.y = pad, int(pad + size * 1.55 * (len(lines) - 1 - k))
                label.color = (176, 32, 32, 255) if _HUD_WARN in lines[k] else (34, 34, 34, 235)
            else:
                label.text = ""          # emptied, never destroyed — see the note above
        top, bottom = (h, 0) if flip else (h, 0)   # see _attach_hud: both targets now agree
        window.projection = Mat4.orthogonal_projection(0, w, bottom, top, -255, 255)
        window.view = Mat4()
        # ⚠ The renderer leaves depth testing and face culling ON and blending OFF — its state for a
        # 3D scene. Under those a 2D caption loses the depth test against the scene it is meant to sit
        # on, loses its alpha, and — when the projection is flipped — has its own glyph quads culled
        # as back faces. Each of those produced a render that succeeded and showed nothing.
        gl.glDisable(gl.GL_DEPTH_TEST)
        gl.glDisable(gl.GL_CULL_FACE)
        gl.glEnable(gl.GL_BLEND)
        gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        gl.glViewport(0, 0, w, h)
        if fbo is not None:
            gl.glBindFramebuffer(gl.GL_FRAMEBUFFER, fbo)
        batch.draw()
        if fbo is not None:
            gl.glBindFramebuffer(gl.GL_FRAMEBUFFER, 0)
        # ⚠⚠ RESTORED, AND THIS IS NOT HYGIENE. GL state is global, and the caption was leaving depth
        # testing and face culling OFF for whatever drew next. In a ONE-frame render nothing draws
        # next, so `--headless` and every figure taken with it were correct and the bug was invisible.
        # In a MULTI-frame render every frame after the first drew with no depth test: the turntable's
        # second frame showed the nucleus THROUGH the membrane, and the mean colour went orange to
        # pink while every count in the caption stayed right. It looks like a plausible picture, which
        # is why it was only caught by rendering three frames of a ROTATIONALLY SYMMETRIC view and
        # comparing them numerically — three pictures that must be identical, and were not.
        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glEnable(gl.GL_CULL_FACE)
        gl.glDisable(gl.GL_BLEND)

    def draw_with_hud() -> None:
        original()
        lines = list(get_lines())
        if not lines:
            return
        # ⚠ BOTH TARGETS, ALWAYS, and they need OPPOSITE orientations.
        #   `_frame_fbo` -> `_frame_texture` is what `get_pixels` reads, and BOTH file paths go
        #     through it: `--headless` and the interactive `p` screenshot. It is read raw and then
        #     reversed row-wise, so a caption drawn bottom-up lands upside down — it did, legibly.
        #   framebuffer 0 is what the window presents, in GL's own bottom-left origin.
        # Painting only the one you happen to be testing is how `p` would have gone on writing
        # caption-free PNGs from a window that showed the caption the whole time.
        fbo = getattr(renderer, "_frame_fbo", None)
        if fbo is not None:
            _paint(lines, flip=True, fbo=fbo)
        if on_screen:
            _paint(lines, flip=False, fbo=None)

    renderer._draw = draw_with_hud


def _render(cell: CellFile, only: list[str], headless: Path | None, thicken: float,
            view: str, frame: int = 0, cut: tuple[int, float, int] | None = None) -> int:
    import warp as wp
    import warp.render

    wp.init()
    title = "Aleph — native cell"
    pos, front = VIEWS[view]
    kw = dict(title=f"{title} — {view}", screen_width=1600, screen_height=1000,
              camera_pos=pos, camera_front=front, up_axis="Y",
              background_color=(0.97, 0.96, 0.94), near_plane=0.05, far_plane=400.0)
    # Feature-detect rather than pin a warp version: this project runs one warp on the dev machine and
    # another on the GPU host, and a viewer that hard-fails on the older one is a viewer nobody opens.
    for opt, val in (("draw_grid", False), ("draw_axis", False),
                     # ⚠ warp defaults `vsync=False`, so the loop below presents as fast as the
                     # machine allows. With this scene that is not a frame-rate question — see the
                     # note on the idle branch.
                     ("vsync", True),
                     ("headless", headless is not None)):
        if _accepts(warp.render.OpenGLRenderer, opt):
            kw[opt] = val
    renderer = warp.render.OpenGLRenderer(**kw)

    def _draw(frame: int, t: float) -> list[str]:
        """Draw one frame. Called once for a static export and per frame for a sequence."""
        pos_by_pop = cell.frames.get(frame, cell.positions)
        out: list[str] = []
        renderer.begin_frame(t)
        for nm in only:
            _draw_population(renderer, nm, pos_by_pop[nm], cell, thicken, out, cut)
        renderer.end_frame()
        return out

    drawn: list[str] = []
    pos_by_pop = cell.frames.get(frame, cell.positions)
    into = camera_is_on_the_kept_side(view, cut)
    caption = _caption(cell, view, cut, into)
    hud_lines = list(caption)          # ONE holder, attached once; the live loop refreshes it
    # ⚠ ATTACHED BEFORE THE ONLY DRAW. The first version attached it afterwards and then drew every
    # population a SECOND time so the wrapped `_draw` would run — which doubled a 4.5 s render and
    # made the app markedly slower than the thing it was fixing. Measured, not noticed: 24.5 s for a
    # full render, of which ~9 s was the duplicate pass. The caption needs nothing the draw produces.
    _attach_hud(renderer, lambda: hud_lines, on_screen=headless is None)
    renderer.begin_frame(0.0)
    # Vertices are the population's nodes and indices are its segment pairs, flattened. Passing an
    # empty index list draws nothing and does not complain — the first version of this viewer did
    # exactly that and printed a confident summary of what it had not drawn. That lesson lives in
    # `_draw_population` now, and this path calls it rather than restating it.
    for name in only:
        _draw_population(renderer, name, pos_by_pop[name], cell, thicken, drawn, cut)
    renderer.end_frame()

    print(f"[app] {cell.summary()}" + (f"  ·  frame {frame}" if len(cell.frames) > 1 else ""), flush=True)
    print(f"[app] drawn: {', '.join(drawn)}", flush=True)
    if into is not None and into > 0.0:
        axis, offset, sign = cut
        print(f"[app] ⚠ THE CUT IS FACING AWAY FROM YOU. The `{view}` camera sits {into:.1f} um INSIDE "
              f"the half `--cut {'xyz'[axis]}{'+' if sign > 0 else '-'}` keeps, so you are looking at "
              f"that half's closed outer surface and the cut face is behind it. The picture will look "
              f"like an UNCUT cell; only the element counts above say anything was removed. "
              f"Use `--cut {'xyz'[axis]}{'-' if sign > 0 else '+'}` with this camera, or a camera on "
              f"the other side. NOT auto-corrected: a viewer that quietly changes your cut is a "
              f"viewer that can be quoted about a picture you did not ask for.", flush=True)
    print(f"[app] thinning: {cell.header.get('thinning')}", flush=True)
    if cut is not None:
        ax, off, sign = cut
        print(f"[app] CUT: keeping {'xyz'[ax]} {'>=' if sign > 0 else '<='} {off:g} um. "
              "A spatial selection, NOT a thinning: everything on the kept side is drawn in full, "
              "and the counts above say how many elements the plane hid.", flush=True)
    # ⚠ flush=True on every line the app prints about itself. The terminal is this app's LOG — the
    # state line, the cut plane, the hidden counts, and the full text of every warning the HUD
    # abbreviates onto the picture — and Python buffers stdout when it is not a
    # tty, so launching it with output piped anywhere showed NOTHING for as long as it ran. Found by
    # actually starting the window rather than by reading the code: 26 s alive, 0 bytes out.
    print("[app] ⚠ geometry only — no force is evaluated, nothing has moved, and a picture is not "
          "evidence.", flush=True)

    if headless is not None:
        # ⚠ `OpenGLRenderer.save()` takes no path — it is the interactive screenshot, not a file
        # writer, and calling it with one raises. The pixels come back through `get_pixels` into a
        # warp array, which is also the only route that works when there is no window at all.
        buf = wp.zeros((renderer.screen_height, renderer.screen_width, 3), dtype=wp.uint8)
        renderer.get_pixels(buf, split_up_tiles=False, mode="rgb", use_uint8=True)
        # ⚠ NO ROW REVERSAL. `get_pixels` already hands back a TOP-LEFT-origin image; the
        # `[::-1]` that stood here for the life of this file flipped every PNG it ever wrote
        # upside down — headless renders, `make cell-png`, and the `p` screenshot alike. Measured
        # with a thick marker bar lying entirely at +Y and a camera on the axis: raw centroid row
        # 215/800 (TOP, correct), reversed 584/800 (BOTTOM). "GL origin is bottom-left" is true of
        # the default framebuffer and not of what this call returns.
        img = buf.numpy()
        _write_png(headless, img)
        print(f"[app] wrote {headless}  ({img.shape[1]}x{img.shape[0]})", flush=True)
        return 0

    # ── the interactive app ──────────────────────────────────────────────────────────────────────
    # ⚠ Before this, every option was a command-line flag, so seeing the cell a second way meant
    # quitting and reloading — 16 s to re-read 94 MB and re-upload 4.3 M segments, for a keystroke's
    # worth of question. That is a renderer with arguments, not an app. The state below is live and
    # the keys mutate it; nothing is re-read from disk and no geometry is re-parsed.
    series = load_run_series(cell.path, cell.header) if cell.path is not None else None
    state = _AppState(only=list(only), cut=cut,
                      frames=sorted(cell.frames) if len(cell.frames) > 1 else [frame],
                      all_pops=[n for n in cell.positions], series=series)
    if series is None:
        # ⚠ SAID, not skipped. An absent panel and an absent series look identical on screen, and
        # this app has already been wrong once by drawing nothing for an empty index list in silence.
        print("[app] no gamma series: the sibling .json is missing, carries no "
              "`gamma.trace_pn_per_um`, or records no --snapshot-every to index it by. Nothing is "
              "plotted rather than something being plotted against a guessed stride.", flush=True)
    _install_keys(renderer, state, cell, headless_dir=Path("."))
    # The live half. `_state_line` already produces the status bar; it is now painted as well as
    # printed, so the window answers "what am I looking at" without a second surface. The HUD was
    # attached above and reads `hud_lines`, so this refreshes the holder rather than hooking twice.
    def _refresh_hud() -> None:
        hud_lines[:] = [ln.removeprefix("[app] ")
                        for ln in _state_line(state).splitlines()] + caption

    _refresh_hud()

    if len(state.frames) > 1:
        # ⚠ Playback re-sends POSITIONS only. The topology went up once and does not change, which is
        # what makes a sequence affordable at all: a frame of 4.59 M positions is 55 MB against 39 MB
        # of topology written once. Re-uploading the topology per frame would triple the file and
        # teach the viewer that the connectivity is time-dependent, which it is not.
        print(f"[app] sequence: {len(state.frames)} frames. The cell MOVES; the topology does not.", flush=True)
    print(_help_text(state, cell), flush=True)

    # What each population was last issued with. See the loop's note on why this is
    # sound rather than a stale-state risk.
    issued: dict[str, tuple] = {}
    last_cam = (renderer._camera_pos, renderer._camera_front, renderer.camera_fov)
    i = 0
    while renderer.is_running():
        if state.dirty or (state.playing and len(state.frames) > 1):
            if state.playing and len(state.frames) > 1:
                state.frame_i = (state.frame_i + 1) % len(state.frames)
            out: list[str] = []
            renderer.begin_frame(i * 0.05)
            pos_by_pop = cell.frames.get(state.frames[state.frame_i], cell.positions)
            # ⚠ HIDING IS AN EMPTY INDEX LIST, NEVER A SKIP. A name that is simply not issued keeps
            # last frame's geometry on screen under that name, so the picture would show a population
            # the state says is off. `_CUT_NOTHING` is what actually removes it.
            #
            # ⚠ AND ONLY WHAT CHANGED IS RE-ISSUED. Re-issuing all of them cost 4.5 s PER KEYSTROKE —
            # measured — because warp rebuilds a capsule instancer per call and the cortex alone is
            # 4.1 M segments. The instancer persists across frames, which is the same property the
            # paragraph above warns about; used deliberately it makes a population toggle cost one
            # population instead of twelve. The key is everything `_draw_population` reads, so a
            # change it would have rendered cannot be skipped.
            for nm in state.all_pops:
                shown = nm in state.only
                this_cut = state.cut if shown else _CUT_NOTHING
                key = (shown, this_cut, state.frames[state.frame_i], thicken)
                if issued.get(nm) == key:
                    continue
                issued[nm] = key
                _draw_population(renderer, nm, pos_by_pop[nm], cell, thicken, out, this_cut)
            renderer.end_frame()
            _refresh_hud()
            if state.dirty:
                print(f"[app] {_state_line(state)}", flush=True)
            state.dirty = False
        else:
            # ⚠⚠ NOTHING CHANGED, SO NOTHING IS DRAWN. THIS IS THE FIX FOR "the whole screen
            # freezes", and it was not a frame-rate problem. Measured on the native cell: ONE idle
            # redraw of 4.1 M capsule instances takes **7.4 s**, and this branch used to run it in a
            # tight loop forever, with `vsync=False` and no sleep. The window was therefore holding
            # the GPU at 100% for as long as it was open, whether or not anyone touched it, and the
            # compositor starved with it — which is why the desktop stopped, not just the app.
            #
            # A static picture does not need redrawing. So the branch pumps events, lets warp move
            # the camera, and re-renders ONLY if the camera actually moved. Idle cost goes from
            # 7.4 s of GPU per frame, forever, to zero.
            #
            # ⚠ The event pump and `_process_inputs` are BOTH required and neither is optional:
            # warp calls `_process_inputs` from `_draw`, so skipping the draw would silently kill
            # w/a/s/d while leaving the mouse working — a half-dead camera is worse than a slow one.
            renderer.app.platform_event_loop.step(0.03)
            if hasattr(renderer, "_process_inputs"):
                renderer._process_inputs()
            cam = (renderer._camera_pos, renderer._camera_front, renderer.camera_fov)
            if cam != last_cam:
                last_cam = cam
                renderer.begin_frame(i * 0.05)
                renderer.end_frame()
        i += 1
    return 0


#: A cut that keeps nothing. Hiding a population re-issues it with this rather than skipping it, so
#: the previous frame's geometry cannot survive under the same name.
_CUT_NOTHING = (0, float("inf"), 1)


def load_run_series(path: Path, header: dict | None = None) -> tuple[list[float], int] | None:
    r"""The gamma trace of the run that produced ``path``, and its step stride, or ``None``.

    A sequence is scrubbed frame by frame and the frames mean nothing on their own: they are steps of
    a run whose measured quantity lives in the sibling record. **Scrubbing geometry with no idea where
    you are on the curve is half an app**, so if ``<stem>.json`` is beside the ``.alephcell`` and
    carries ``gamma.trace_pn_per_um``, it is read and shown under the frame counter.

    ⚠ **This does NOT make the viewer a measuring instrument.** The series is read from the record and
    displayed; nothing here computes it, judges it, or decides whether it may be quoted. A record whose
    verdict says its window opened on a transient still plots — the plot is the data, and the record's
    own fields are where the verdict lives. The status line says so.

    ⚠ And the stride is READ, never assumed. A sequence written with ``--snapshot-every 6000`` has
    frames 6,000 steps apart, and a viewer that assumed one sample per frame would index the trace by
    the frame number and point at the wrong place with total confidence.

    Returns:
        ``(trace, steps_per_frame)``, or ``None`` when there is no sibling record, no trace in it, or
        no recorded stride to index it by. **Every one of those is a silent absence, so the caller
        prints what was missing rather than simply omitting the panel.**
    """
    rec = path.with_suffix(".json")
    if not rec.exists():
        return None
    try:
        d = json.loads(rec.read_text())
    except (OSError, ValueError):
        return None
    trace = (d.get("gamma") or {}).get("trace_pn_per_um")
    # ⚠ The .alephcell HEADER is the first place to look, not the sibling record. The stride indexes
    # the frames, and the frames are in this file: a header value cannot disagree with the frames it
    # ships beside, whereas a sibling record can be from a different run with the same stem. The
    # record is the fallback, not the source.
    stride = (header or {}).get("snapshot_every")
    stride = stride or (d.get("argv_parsed") or {}).get("snapshot_every") or d.get("snapshot_every")
    if not stride:
        # The record may only carry the raw argv. Reading it is better than assuming 1.
        #
        # ⚠ `argv` is stored as a STRING, not a list, and the first version of this did
        # `if "--snapshot-every" in argv: int(argv[argv.index("--snapshot-every") + 1])`. On a string
        # both of those still work and both mean something else: `in` is a substring test and
        # `.index` is a CHARACTER offset, so `argv[i + 1]` is one letter. It raised ValueError, was
        # caught, and the panel silently never appeared. **Split first, then the list operations mean
        # what they read as.**
        argv = (d.get("provenance") or {}).get("argv") or []
        tokens = argv.split() if isinstance(argv, str) else list(argv)
        if "--snapshot-every" in tokens:
            try:
                stride = int(tokens[tokens.index("--snapshot-every") + 1])
            except (IndexError, ValueError):
                stride = None
    if not trace or not stride:
        return None
    return [float(v) for v in trace], int(stride)


#: How wide the terminal sparkline is. Purely cosmetic; nothing is decided by it.
_SPARK_W = 56
_SPARK = "▁▂▃▄▅▆▇█"


def sparkline(trace: list[float], mark: int | None = None, width: int = _SPARK_W) -> str:
    r"""A fixed-width sparkline over ``trace``, with ``mark`` (an index into it) shown as ``|``.

    ⚠ **A sparkline is not a figure and may not be read as one.** It has no axis, its vertical range is
    the window's own min and max so a flat series and a doubling series look identical, and each column
    is a MEAN over ``len(trace)/width`` samples, which hides everything faster than that. The run's real
    figures are written beside its record. This exists to answer one question — *where in the run am I*
    — and the caller prints the actual endpoint values next to it so the shape is never the only thing.
    """
    if not trace:
        return ""
    n = len(trace)
    cols: list[float] = []
    for i in range(width):
        lo, hi = i * n // width, max((i + 1) * n // width, i * n // width + 1)
        chunk = trace[lo:hi]
        cols.append(sum(chunk) / len(chunk))
    a, b = min(cols), max(cols)
    span = (b - a) or 1.0
    out = [_SPARK[min(int((c - a) / span * len(_SPARK)), len(_SPARK) - 1)] for c in cols]
    if mark is not None and n:
        out[min(max(mark * width // n, 0), width - 1)] = "|"
    return "".join(out)


@dataclass
class _AppState:
    """What the app is showing right now. Mutated by keys, read by the draw loop."""

    only: list[str]
    cut: tuple[int, float, int] | None
    frames: list[int]
    all_pops: list[str]
    frame_i: int = 0
    series: tuple[list[float], int] | None = None   # (gamma trace, steps per frame) from the record
    playing: bool = True
    dirty: bool = True
    shot: int = 0


def _state_line(st: _AppState) -> str:
    """One line saying exactly what is on screen. Painted ON the window by the HUD and printed too.

    ⚠ It used to say "the app's status bar is the terminal", and that was the defect: an app tells
    you its state in itself. `_attach_hud` paints this line onto the picture; the print stays because
    a terminal keeps scrollback and a window does not.
    """
    cut = "none" if st.cut is None else (
        f"{'xyz'[st.cut[0]]} {'>=' if st.cut[2] > 0 else '<='} {st.cut[1]:g} um")
    hidden = [n for n in st.all_pops if n not in st.only]
    frame = (f" · frame {st.frames[st.frame_i]} of {len(st.frames)}"
             f"{' (playing)' if st.playing else ' (paused)'}" if len(st.frames) > 1 else "")
    line = (f"showing {len(st.only)}/{len(st.all_pops)} populations · cut {cut}{frame}"
            + (f" · hidden: {', '.join(hidden)}" if hidden else ""))
    if st.series is None:
        return line
    trace, stride = st.series
    i = min(st.frames[st.frame_i] * stride if len(st.frames) > 1 else 0, len(trace) - 1)
    return (f"{line}\n[app] γ {sparkline(trace, i)}  step {i:,}/{len(trace):,}  "
            f"γ={trace[i]:.2f} pN/µm  (run {min(trace):.1f}…{max(trace):.1f}) "
            f"— READ from the record, judged nowhere here")


def _help_text(st: _AppState, cell: CellFile) -> str:
    """The key map, printed at startup and on `h`."""
    pops = "  ".join(f"{i + 1 if i < 9 else 0}:{n}" for i, n in enumerate(st.all_pops[:10]))
    return (
        "\n[app] ── keys ─────────────────────────────────────────────────────────────\n"
        f"[app]  1-9,0  toggle a population   {pops}\n"
        "[app]  f / n   show all / show none\n"
        "[app]  k       cycle the cut: off -> x- -> x+ -> y- -> y+ -> z- -> z+ -> off\n"
        "[app]  [ / ]   move the cut plane by 0.5 um\n"
        "[app]  r       play / pause a sequence      , / .   step one frame\n"
        "[app]  p       write a PNG of what is on screen\n"
        "[app]  h       this help, and the current state\n"
        "[app]  o       what this cell HAS -- and what it does NOT, which no census says\n"
        "[app] ── the CAMERA is the renderer's, not this app's ─────────────────────\n"
        "[app]  mouse   drag to look, scroll to zoom\n"
        "[app]  w/a/s/d or arrows   move the camera\n"
        "[app]  ⚠ warp also owns c g i x t b TAB space esc. This app's keys were moved\n"
        "[app]     off them (a->f, c->k, i->o, space->r); space is warp's PAUSE and it\n"
        "[app]     spins inside end_frame, which looks exactly like a hang.\n"
        "[app] ─────────────────────────────────────────────────────────────────────\n"
        f"[app] {_state_line(st)}")


#: The cut ring: OFF, then each axis in each direction. Seven states, not a boolean — "which half"
#: is the question, and which half is interesting is not known in advance.
_CUT_RING: tuple[tuple[int, float, int] | None, ...] = (
    None, (0, 0.0, -1), (0, 0.0, 1), (1, 0.0, -1), (1, 0.0, 1), (2, 0.0, -1), (2, 0.0, 1))


def _apply_key(st: _AppState, ch: str) -> bool:
    """Apply one keystroke to the state. Returns whether anything changed.

    ⚠ **NONE OF THESE MAY COLLIDE WITH THE RENDERER'S OWN KEYS**, and four of them did. warp binds
    `W A S D` and the arrows to camera movement, `C` to the axis, `G` to the grid, `I` to its own info
    overlay, `X` wireframe, `T` depth, `B` culling, `TAB` skip-rendering, `SPACE` PAUSE and `ESC`
    close — all through its own `on_key_press`, which this app's callback is registered ALONGSIDE
    rather than instead of. So `a` both showed every population and strafed the camera left, `c` cut
    the cell and toggled the axis, `i` printed the absent list and turned on warp's overlay, and
    `space` did the app's play/pause AND warp's pause, which spins inside `end_frame` and looks
    exactly like a hang. Moved to `f` / `k` / `o` / `r`, which warp does not claim.

    ⚠ **Separated from the callback so it can be checked on a machine that cannot open a window.**
    The dev machine has no display and the GPU host has no runner; a state machine reachable only
    through a live GL window is a state machine nobody verifies — which is the same hole as a
    CUDA-gated test that skips everywhere, found on 2026-08-21 one layer down.
    """
    if ch and ch in "1234567890":
        i = (int(ch) - 1) if ch != "0" else 9
        if i >= len(st.all_pops):
            return False
        nm = st.all_pops[i]
        st.only = [n for n in st.only if n != nm] if nm in st.only else st.only + [nm]
    elif ch == "f":
        st.only = list(st.all_pops)
    elif ch == "n":
        st.only = []
    elif ch == "k":
        here = next((k for k, v in enumerate(_CUT_RING) if v == st.cut), 0)
        st.cut = _CUT_RING[(here + 1) % len(_CUT_RING)]
    elif ch in "[]" and st.cut is not None:
        ax, off, sg = st.cut
        st.cut = (ax, off + (0.5 if ch == "]" else -0.5), sg)
    elif ch == "r":
        st.playing = not st.playing
    elif ch in ",." and len(st.frames) > 1:
        st.playing = False
        st.frame_i = (st.frame_i + (1 if ch == "." else -1)) % len(st.frames)
    else:
        return False
    return True


def describe_cell(cell: CellFile) -> str:
    r"""What this file HAS, and — the point of the function — **what it does not have.**

    ⚠ **A census of eleven names does not tell a reader it is missing a twelfth.** You have to already
    know the twelve. That is how the PHASE 4 tau runs were stepped and read for a full night before
    anyone noticed the cell has **no myosin in it**: the records list what they contain, correctly, and
    nothing lists what they lack.

    So the absent list is computed, not the present one. The reference is **derived, never typed**: the
    builder modules that exist in `aleph/world/build/`, read from the filesystem. If a builder is added
    the reference grows with no edit here, and if this list is wrong it is wrong about the tree rather
    than about a constant somebody forgot.

    ⚠ **A name in "not in this file" is not an accusation.** `--core-only` builds three populations on
    purpose and every other builder is legitimately absent; a driver may have no reason to stand one.
    The function reports a difference and says nothing about whether the difference is a defect. What it
    refuses to do is let the difference stay invisible.
    """
    lines = [f"[app] {cell.path.name if cell.path else '<in memory>'} — {cell.summary()}"]
    have = sorted(cell.positions)
    lines.append(f"[app] HAS ({len(have)}): {', '.join(have)}")

    build_dir = Path(__file__).resolve().parents[1] / "world" / "build"
    known = sorted(f.stem for f in build_dir.glob("*.py") if not f.stem.startswith("_"))
    if not known:
        lines.append(f"[app] ⚠ cannot report absences: no builder modules found under {build_dir}. "
                     "An empty reference is not an empty difference.")
        return "\n".join(lines)

    # `envelope.py` builds `nuclear_envelope`; match on containment so a builder stem is credited when
    # the population it stands is named after it. ⚠ Same substring compromise the constraint scanner
    # declares, and it can over-credit -- reported as a difference, never as a verdict.
    missing = [k for k in known if not any(k in pop for pop in have)]
    if missing:
        lines.append(f"[app] ⚠ NOT IN THIS FILE ({len(missing)} of {len(known)} builders in the tree): "
                     f"{', '.join(missing)}")
        lines.append("[app]   Not an accusation — a driver may have had no reason to stand one, and "
                     "`--core-only` omits nine deliberately. It is here so the omission is not invisible.")
    else:
        lines.append(f"[app] every one of the {len(known)} builders in the tree has a population here.")

    rec = cell.path.with_suffix(".json") if cell.path else None
    if rec and rec.exists():
        try:
            d = json.loads(rec.read_text())
        except (OSError, ValueError):
            d = {}
        prov = d.get("provenance") or {}
        for k in ("build_commit", "closure_digest", "warp", "env_prefix"):
            if prov.get(k):
                lines.append(f"[app]   {k}: {prov[k]}")
        if prov.get("argv"):
            lines.append(f"[app]   argv: {prov['argv']}")
        for k in ("quantitative_claim_status", "not_a_claim"):
            if d.get(k):
                lines.append(f"[app]   ⚠ {k}: {d[k]}")
    elif rec:
        lines.append(f"[app]   no sibling record at {rec.name} — provenance unknown, which is not the "
                     "same as provenance clean.")
    return "\n".join(lines)


def _install_keys(renderer, st: _AppState, cell: CellFile, headless_dir: Path) -> None:
    """Bind the keys. Every binding mutates ``st`` and sets ``dirty``; none touches the disk.

    ⚠ The cut cycle is a ring of SEVEN, not a boolean. A two-state toggle would make "which half"
    another flag to remember, and the whole reason the cut exists is that the interesting half is not
    known in advance — it depends on which population you are looking for.
    """

    def on_key(symbol: int, modifiers: int) -> None:
        ch = chr(symbol) if 32 <= symbol < 127 else ""
        if ch == "p":
            st.shot += 1
            _snapshot(renderer, headless_dir / f"cell_shot{st.shot:02d}.png")
            return
        if ch == "h":
            print(_help_text(st, cell), flush=True)
            return
        if ch == "o":
            print(describe_cell(cell), flush=True)
            return
        if _apply_key(st, ch):
            st.dirty = True

    renderer.register_key_press_callback(on_key)


def _snapshot(renderer, path: Path) -> None:
    """Write what is on screen. Same read-back route as ``--headless``, so they cannot diverge."""
    import warp as wp

    buf = wp.zeros((renderer.screen_height, renderer.screen_width, 3), dtype=wp.uint8)
    renderer.get_pixels(buf, split_up_tiles=False, mode="rgb", use_uint8=True)
    _write_png(path, buf.numpy())
    print(f"[app] wrote {path}", flush=True)


def _accepts(fn, name: str) -> bool:
    import inspect

    try:
        return name in inspect.signature(fn.__init__).parameters
    except (TypeError, ValueError):                      # pragma: no cover
        return False


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("path", type=Path, nargs="?", default=None,
                    help="an .alephcell export. OMIT IT to open the most recent one under "
                         "aleph/outputs/ — the tool prints which file it chose and why, because a "
                         "viewer that silently picks a file is a viewer that can show you the wrong "
                         "cell and let you believe it is the current one.")
    ap.add_argument("--only", default="",
                    help="comma-separated populations to draw. Default: all of them. Named "
                         "populations are drawn IN FULL; there is no thinning anywhere in this tool.")
    ap.add_argument("--view", default="front", choices=sorted(VIEWS),
                    help="named camera. `basal` looks up at the contact disc, which is where the "
                         "ventral fibres, the transverse arcs and the lamellipodium all are — an "
                         "axis-on view puts every one of them behind the nucleus.")
    ap.add_argument("--thicken", type=float, default=1.0,
                    help="multiply every draw radius. For legibility at a distance; the radii "
                         "themselves are stated by --radii and are not claims.")
    ap.add_argument("--radii", action="store_true", help="print the draw-radius table and exit")
    ap.add_argument("--frame", type=int, default=0,
                    help="which frame of a sequence to draw. Only meaningful with --headless; "
                         "interactive playback walks them all. A frame this file does not have is "
                         "REFUSED rather than clamped to the nearest — silently drawing a different "
                         "frame than the one asked for is how a picture stops being of what it says.")
    ap.add_argument("--cut", default="", metavar="AXIS[:OFFSET]",
                    help="cut the cell open on a plane and draw the half that remains — e.g. "
                         "`y-`, `x+`, `z-:1.5`. AXIS is x|y|z, the sign says which half to KEEP "
                         "(`-` keeps below the plane, `+` above), OFFSET defaults to 0. "
                         "A SPATIAL SELECTION, never a thinning: every element on the kept side is "
                         "drawn in full and the hidden count is printed. At native density the "
                         "cortex is 4.1 M filaments in a 0.2 um shell and hides everything inside "
                         "it, including itself — a cut is how the interior is seen without "
                         "misrepresenting the density, which downsampling would.")
    ap.add_argument("--headless", type=Path, default=None, help="render one frame to a PNG and exit")
    ap.add_argument("--info", action="store_true", help="print the header and exit, drawing nothing")
    ap.add_argument("--list", action="store_true", dest="list_exports",
                    help="what can be opened, and what each one is. Needs no path.")
    args = ap.parse_args(argv)

    if args.radii:
        print(describe_radii())
        return 0

    if args.list_exports:
        print(list_exports(), flush=True)
        return 0

    if args.path is None:
        args.path = _newest_export()
        if args.path is None:
            print("refused: no .alephcell found under aleph/outputs/. Make one with "
                  "`world_export_cell.py`, or name a file.")
            return 2

    cut: tuple[int, float, int] | None = None
    if args.cut:
        try:
            cut = parse_cut(args.cut)
        except ValueError as exc:
            print(f"refused: {exc}")
            return 2
    if not args.path.exists():
        raise SystemExit(f"no such file: {args.path}")

    cell = load(args.path)
    if args.info:
        h = dict(cell.header)
        h.pop("blocks", None)
        print(json.dumps(h, indent=1))
        print()
        # ⚠ The header is what the writer chose to say. `describe_cell` is what the file does NOT
        # contain, which no header has ever said, and this is the only path to it on a machine with
        # no display -- the GPU host has none, and it is where these files are written.
        print(describe_cell(cell))
        return 0

    only = [w.strip() for w in args.only.split(",") if w.strip()] or cell.populations
    unknown = [w for w in only if w not in cell.positions]
    if unknown:
        raise SystemExit(f"refused: {unknown} are not in this file. Present: {cell.populations}")
    if not only:
        raise SystemExit("refused: an empty selection would open a black window")
    if args.frame and args.frame not in cell.frames:
        raise SystemExit(
            f"refused: --frame {args.frame} is not in this file. It has "
            f"{sorted(cell.frames) if len(cell.frames) > 1 else 'one frame (0)'}.")
    return _render(cell, only, args.headless, args.thicken, args.view, args.frame, cut)


def list_exports(root: Path = Path("aleph/outputs")) -> str:
    r"""Every openable cell under ``root``, with what each one is — so the app can be opened blind.

    ⚠ **Before this, using the app meant already knowing a path.** That is a renderer with an
    argument, not an application: the one thing a person opening it does not have is the name of the
    file they want. `--list` is the answer to *"what is there"*, and it answers with the facts that
    decide which one you want — how many populations, how many frames, whether a gamma series is
    beside it — rather than with a directory listing they could have got from `ls`.

    ⚠ **Every row is READ, not inferred.** The frame count and the population names come from the
    file's own header, and the gamma column says whether `load_run_series` would actually find a
    series, not whether a `.json` happens to sit next to it. A listing that guessed would send
    somebody to open a 367 MB file to find out it is not the one.
    """
    import time as _t

    found = sorted(root.rglob("*.alephcell"), key=lambda f: f.stat().st_mtime, reverse=True)
    if not found:
        return (f"no .alephcell under {root}. Make one with `world_export_cell.py`, or with a "
                "driver's --snapshot-every. ⚠ An empty listing is a statement about this directory, "
                "not about whether any cell has ever been exported.")
    rows = [f"{len(found)} openable under {root}, newest first:",
            f"  {'file':<44} {'MB':>7} {'pops':>5} {'frames':>7} {'γ':>3}  age"]
    for f in found:
        try:
            cell = load(f)
            pops, frames = len(cell.positions), max(len(cell.frames), 1)
            gamma = "yes" if load_run_series(f, cell.header) else " no"
        except Exception as exc:                       # a broken file is LISTED, not hidden
            rows.append(f"  {str(f.relative_to(root)):<44}    ⚠ unreadable: {type(exc).__name__}")
            continue
        age_h = (_t.time() - f.stat().st_mtime) / 3600.0
        age = f"{age_h * 60:.0f} min" if age_h < 2 else f"{age_h:.1f} h"
        rows.append(f"  {str(f.relative_to(root)):<44} {f.stat().st_size / 1e6:>7.0f} "
                    f"{pops:>5} {frames:>7} {gamma:>3}  {age}")
    rows.append("  open one with:  make cell ARGS=<path>      inspect without a window:  --info")
    return "\n".join(rows)


def _newest_export(root: Path = Path("aleph/outputs")) -> Path | None:
    """The most recently modified ``.alephcell`` under ``root``, or None.

    ⚠ **It says which file it picked and how old it is.** A viewer that silently opens "the latest"
    is a viewer that shows a stale cell to someone who believes it is current — the same failure as
    the run host's hand-written commit stamp, which was 101 commits old and said nothing. Choosing
    for the user is fine; choosing silently is not.
    """
    import time as _t

    found = sorted(root.rglob("*.alephcell"), key=lambda f: f.stat().st_mtime, reverse=True)
    if not found:
        return None
    best = found[0]
    age_min = (_t.time() - best.stat().st_mtime) / 60.0
    print(f"[app] no file named; opening the newest of {len(found)} exports under {root}:",
          flush=True)
    print(f"[app]   {best}  ({best.stat().st_size / 1e6:.0f} MB, modified "
          f"{age_min:.0f} min ago)", flush=True)
    if age_min > 120:
        print(f"[app]   ⚠ that is {age_min / 60:.1f} HOURS old — it may not be the cell you just "
              "built. Name the file if you meant a different one.", flush=True)
    return best


def _demo() -> None:
    """Sanity Gate for the cut — the properties that separate a cutaway from a thinning.

    ⚠ Placed ABOVE ``__main__`` deliberately. `world/observe_gamma.py` shipped its device section
    *below* its own entry point on 2026-08-21, so the documented invocation bound none of it and only
    a static check could see that. Nothing is defined after the block below.

    What is checked without a renderer: that the cut keeps whole elements, that it is exhaustive and
    disjoint, and that "no cut" is not quietly a cut.
    """
    pos = np.array([[0.0, 0.0, -2.0], [0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 2.0]])
    seg = np.array([[0, 1], [1, 2], [2, 3]], np.int32)      # below, STRADDLING, above

    # 1. No cut is not a cut. Every element survives and the hidden count is zero.
    kept, hidden = _apply_cut(pos, seg, None)
    assert kept.shape == seg.shape and hidden == 0

    # 2. An element with one vertex on each side is DROPPED, never clipped. Clipping would invent a
    #    filament end the build does not have, and a viewer that invents ends can be quoted about
    #    lengths. Below-keeping gets [0,1]; above-keeping gets [2,3]; the straddler goes to neither.
    below, hid_b = _apply_cut(pos, seg, (2, 0.0, -1))
    above, hid_a = _apply_cut(pos, seg, (2, 0.0, +1))
    assert below.tolist() == [[0, 1]], below.tolist()
    assert above.tolist() == [[2, 3]], above.tolist()
    assert hid_b == 2 and hid_a == 2

    # 3. Exhaustive and disjoint over the elements that do NOT straddle: nothing is drawn twice by
    #    the two halves, and nothing whole is lost between them.
    pair = {tuple(r) for r in below.tolist()} , {tuple(r) for r in above.tolist()}
    assert not (pair[0] & pair[1]), "an element appears in both halves"
    whole = {tuple(r) for r in seg.tolist() if abs(pos[r[0], 2]) and pos[r[0], 2] * pos[r[1], 2] > 0}
    assert (pair[0] | pair[1]) == whole, "a non-straddling element is in neither half"

    # 4. The offset moves the plane, and the sign says which half is KEPT — independently.
    #    ⚠ Written without an `or True` tail. The first draft had one, which makes an assert that
    #    cannot fail, which is the defect this whole file's sibling documents are about.
    #    Node 1 sits at z = -1.0, so `z <= -1.5` holds for node 0 alone and segment [0,1] — which
    #    needs BOTH ends — is dropped. The first draft asserted [[0, 1]] here and hid the wrong
    #    expectation behind `or True`; removing the tail is what surfaced it.
    assert _apply_cut(pos, seg, (2, -1.5, -1))[0].tolist() == [], "z <= -1.5 keeps no whole segment"
    assert _apply_cut(pos, seg, (2, 99.0, +1))[0].shape[0] == 0, "nothing is above z=99"
    assert _apply_cut(pos, seg, (2, 99.0, -1))[0].shape[0] == 3, "everything is below z=99"

    # 5. Meshes go through the same predicate: a face is kept only if every vertex survives.
    tri = np.array([[0, 1, 2]], np.int32)                    # straddles z=0
    assert _apply_cut(pos, tri, (2, 0.0, -1))[0].shape[0] == 0

    # ── the app's state machine, checked without a window ────────────────────────────────────────
    pops = ["cortex", "membrane", "microtubule", "nmii"]
    st = _AppState(only=list(pops), cut=None, frames=[0, 1, 2], all_pops=list(pops))

    # 6. A population toggles off and back on, and the toggle is by INDEX into all_pops — so the key
    #    map cannot drift from what is on screen just because `only` was reordered.
    assert _apply_key(st, "3") and "microtubule" not in st.only
    assert _apply_key(st, "3") and "microtubule" in st.only
    assert not _apply_key(st, "9"), "a key past the population count must change nothing"

    # 7. `f` and `n` are total, and `n` leaves an EMPTY list rather than a sentinel — the draw loop
    #    re-issues hidden populations with a keep-nothing cut, so "none" must really mean none.
    assert _apply_key(st, "n") and st.only == []
    assert _apply_key(st, "f") and st.only == pops

    # 7b. ⚠ NO APP KEY MAY BE ONE THE RENDERER ALREADY OWNS. Four of them were, and `space` was the
    #     worst: it is warp's PAUSE, which spins inside `end_frame` and reads as a hang. Pinned by
    #     value rather than by comment, so re-introducing one fails here instead of in someone's
    #     window. The list is warp's `on_key_press` plus its WASD/arrow movement.
    warp_owns = set("wasdcgixtb") | {" ", "\t"}
    mine = {"f", "n", "k", "[", "]", "r", ",", ".", "p", "h", "o", *"1234567890"}
    assert not (mine & warp_owns), f"app keys collide with the renderer's: {sorted(mine & warp_owns)}"

    # 8. The cut ring is a RING: seven states, and cycling all seven returns to where it started.
    seen = [st.cut]
    for _ in range(len(_CUT_RING) - 1):
        assert _apply_key(st, "k")
        seen.append(st.cut)
    assert len(set(map(str, seen))) == len(_CUT_RING), seen
    assert _apply_key(st, "k") and st.cut is None, "the ring must close"

    # 9. The plane only moves when there IS a plane, and `[`/`]` are inverses.
    assert not _apply_key(st, "]"), "moving a plane that is off must change nothing"
    _apply_key(st, "k")
    before = st.cut
    assert _apply_key(st, "]") and _apply_key(st, "[") and st.cut == before

    # 10. Stepping a frame pauses playback — a viewer that keeps playing while you step is a viewer
    #     that shows a different frame from the one you asked for, which `--frame` already refuses.
    st.playing = True
    assert _apply_key(st, ".") and st.playing is False and st.frame_i == 1
    assert _apply_key(st, ",") and st.frame_i == 0

    # 11. Hiding is done by a cut that keeps NOTHING, not by skipping the population. A skipped name
    #     leaves the previous frame's geometry on screen under that name.
    tri = np.array([[0, 1, 2]], np.int32)
    assert _apply_cut(pos, tri, _CUT_NOTHING)[0].shape[0] == 0

    # 12. The status line names what is hidden, so the terminal and the window cannot disagree.
    st.only = [n for n in pops if n != "nmii"]
    assert "hidden: nmii" in _state_line(st), _state_line(st)

    print("cell_app self-check OK — the cut keeps whole elements, halves are disjoint and "
          "exhaustive, no-cut is not a cut, and the app's state machine holds without a window")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--self-check":
        _demo()
        raise SystemExit(0)
    raise SystemExit(main(sys.argv[1:]))
