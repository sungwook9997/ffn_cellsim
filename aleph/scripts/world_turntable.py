#!/usr/bin/env python
r"""Pre-render the cell once, then look at it in a browser with no renderer running at all.

**Why this exists.** `aleph/viz/cell_app.py` opens the native cell at full density and lets you turn
it, and that is the right tool for asking a question of the geometry. It is the wrong tool for
*looking*: one redraw of 4.1 M capsule instances costs about 7.4 s, so a live window spends real GPU
on every frame it shows, and the PI's machine stopped while it was open. The fix in the app was to
stop redrawing a picture that had not changed. This goes further: **render every angle ONCE, then
view them as images.** Turning the cell in the browser costs nothing, because nothing is rendering.

⚠ **THIS IS NOT DOWNSAMPLING AND MUST NOT BECOME IT.** Every frame is a full-native render — every
node and every segment of every requested population, `thinning: NONE`, exactly what the app draws.
What is reduced is the number of ANGLES you can look from, which is a property of the viewer and not
of the cell. `feedback-viewer-no-downsample` forbids showing one node in twenty-three; it does not
require that a picture be recomputed while you look at it.

⚠ **What you give up, stated rather than discovered.** The angles are the ones rendered and no
others; the cut is fixed for a sweep; and populations cannot be toggled, because compositing them as
transparent layers would put a near filament behind a far one — depth does not survive being flattened
into a PNG. Those questions belong in the app. This answers "show me the cell", not "let me
interrogate it".

⚠ **Each frame carries its own caption**, from `cell_app._caption`, including the cut-facing-away
warning — so the half of a turn where the cut is invisible says so on the picture rather than in a
scrollback the image will be separated from.

Sanity Gate:
    * dimensional — the camera orbits at a radius in µm about the cell centre; elevation is degrees.
      The cell is not moved and no geometry is recomputed between frames: only the view matrix.
    * boundary cases — one frame is a legal sweep and produces a viewer that simply does not turn; a
      cut whose kept half faces the camera for part of the orbit is REPORTED per frame, never fixed.
    * conservation — the geometry is issued ONCE and every frame draws the same instancers, so all
      frames are of the same cell by construction rather than by assertion.
    * CFL/precision — not applicable: no time integration, no force, no solver.
    * sign sense — azimuth increases counter-clockwise seen from +Y, and the viewer's drag follows the
      same sign, so dragging right turns the cell the way it turns in the app.
    * measurement protocol — this measures nothing. It renders. `quantitative_claim_status` is BLOCKED
      in the manifest it writes, and the caption on every frame says a picture is not evidence.

engine units: µm. Runtime: this is a VIEWER. It evaluates no force and advances no time.

Usage:
    python aleph/scripts/world_turntable.py [CELL] --out DIR [--frames N] [--cut y+] [--only a,b]
    open DIR/index.html
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

#: Frames in a full turn. 120 is one every 3 degrees. The cost of an angle is a file on disk, not
#: memory: the viewer loads a window around the frame you are on, so the number here is bounded by
#: patience at render time (about 3.8 s/frame at full native) rather than by anything the browser has
#: to hold. Declared here rather than defaulted in the parser so it has a stated reason beside it.
DEFAULT_FRAMES = 120

#: Render size. Smaller than the app's 1600x1000 on purpose: this writes N of them, and the whole
#: point is a directory a browser can hold at once.
DEFAULT_WIDTH, DEFAULT_HEIGHT = 1200, 800


def _camera(azimuth_deg: float, elevation_deg: float, radius_um: float
            ) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Camera position and unit front vector for one point on the orbit, looking at the origin."""
    a, e = math.radians(azimuth_deg), math.radians(elevation_deg)
    pos = (radius_um * math.cos(e) * math.sin(a),
           radius_um * math.sin(e),
           radius_um * math.cos(e) * math.cos(a))
    n = math.sqrt(sum(c * c for c in pos)) or 1.0
    return pos, (-pos[0] / n, -pos[1] / n, -pos[2] / n)


def main(argv: list[str] | None = None) -> int:
    """Render one full turn and write a viewer beside it."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("path", nargs="?", type=Path,
                    help="an .alephcell export. Omit to take the newest under aleph/outputs/.")
    ap.add_argument("--out", type=Path, required=True, help="directory for the frames + index.html")
    ap.add_argument("--frames", type=int, default=DEFAULT_FRAMES,
                    help=f"frames in a full 360 turn (default {DEFAULT_FRAMES}, one every 10 deg)")
    ap.add_argument("--elevation", type=float, default=18.0, help="camera elevation [deg]")
    ap.add_argument("--radius", type=float, default=32.0, help="orbit radius [um]")
    ap.add_argument("--width", type=int, default=DEFAULT_WIDTH)
    ap.add_argument("--height", type=int, default=DEFAULT_HEIGHT)
    ap.add_argument("--only", default=None, help="comma-separated populations; default all")
    ap.add_argument("--cut", default=None,
                    help="AXIS[+-][:OFFSET] as in cell_app, e.g. y+ — FIXED for the whole sweep")
    ap.add_argument("--thicken", type=float, default=1.0)
    ap.add_argument("--format", choices=("webp", "png"), default="webp",
                    help="webp (default) is ~3x smaller than png for the same picture. ⚠ It is LOSSY: "
                         "measured on a native frame at the default quality, mean |delta| is 1.3 of "
                         "255 (0.5%%) and the caption stays legible. That is a viewer's trade and it "
                         "is recorded in turntable.json; pass png when the frames are the artifact.")
    ap.add_argument("--quality", type=int, default=82, help="webp quality; ignored for png")
    args = ap.parse_args(argv)

    import warp as wp
    import warp.render

    from aleph.viz.cell_app import (
        _accepts,
        _attach_hud,
        _caption,
        _draw_population,
        _newest_export,
        load,
        parse_cut,
    )

    path = args.path or _newest_export()
    if path is None:
        print("refused: no .alephcell found under aleph/outputs/.")
        return 2
    cell = load(Path(path))
    only = [p.strip() for p in args.only.split(",")] if args.only else list(cell.positions)
    missing = [p for p in only if p not in cell.positions]
    if missing:
        print(f"refused: {missing} not in this file. It has: {sorted(cell.positions)}")
        return 2
    try:
        cut = parse_cut(args.cut) if args.cut else None
    except ValueError as exc:
        print(f"refused: {exc}")
        return 2
    args.out.mkdir(parents=True, exist_ok=True)

    wp.init()
    kw = dict(title="Aleph — turntable", screen_width=args.width, screen_height=args.height,
              up_axis="Y", background_color=(0.97, 0.96, 0.94), near_plane=0.05, far_plane=400.0)
    for opt, val in (("draw_grid", False), ("draw_axis", False), ("headless", True)):
        if _accepts(warp.render.OpenGLRenderer, opt):
            kw[opt] = val
    renderer = warp.render.OpenGLRenderer(**kw)

    # ⚠ ISSUED ONCE, OUTSIDE THE LOOP. That is the whole economy of this script: the instancers
    # persist across frames, so every frame after the first costs a view matrix and a draw rather than a
    # rebuild of 4.1 M capsules. It is also what makes every frame provably the same cell.
    drawn: list[str] = []
    pos_by_pop = cell.positions
    renderer.begin_frame(0.0)
    for name in only:
        _draw_population(renderer, name, pos_by_pop[name], cell, args.thicken, drawn, cut)
    renderer.end_frame()
    print(f"[turntable] {cell.summary()}", flush=True)
    print(f"[turntable] drawn: {', '.join(drawn)}", flush=True)
    print(f"[turntable] thinning: {cell.header.get('thinning')} — every frame is a FULL render; "
          f"what is reduced is the number of angles, not the cell", flush=True)

    caption: list[str] = []
    _attach_hud(renderer, lambda: caption, on_screen=False)

    buf = wp.zeros((renderer.screen_height, renderer.screen_width, 3), dtype=wp.uint8)
    frames: list[dict] = []
    t0 = time.perf_counter()
    for k in range(args.frames):
        azimuth = 360.0 * k / args.frames
        pos, front = _camera(azimuth, args.elevation, args.radius)
        renderer.camera_pos, renderer.camera_front = pos, front
        renderer.update_view_matrix()
        # The caption is rebuilt per frame because the cut-facing-away test depends on the camera,
        # and that is exactly the fact a still image must not lose.
        into = _into_for(pos, cut)
        caption[:] = _caption(cell, f"turntable {azimuth:.0f}°", cut, into)
        renderer.begin_frame(0.0)
        renderer.end_frame()
        renderer.get_pixels(buf, split_up_tiles=False, mode="rgb", use_uint8=True)
        name = f"frame_{k:03d}.{args.format}"
        _save(args.out / name, buf.numpy(), args.format, args.quality)
        frames.append({"file": name, "azimuth_deg": round(azimuth, 2),
                       "cut_faces_away": bool(into is not None and into > 0.0)})
        print(f"[turntable] {k + 1}/{args.frames}  azimuth {azimuth:6.1f}°  {name}", flush=True)
    wall = time.perf_counter() - t0

    total_mb = sum((args.out / f["file"]).stat().st_size for f in frames) / 1e6
    manifest = {
        "schema": "ffn-world-turntable@1",
        "kind": "viewer",
        "quantitative_claim_status": "BLOCKED",
        "not_a_claim": ("A pre-rendered turntable. No force is evaluated, nothing moves, and a "
                        "picture is not evidence. Every frame is a FULL-density render of the same "
                        "issued geometry; the reduction is in ANGLES, never in nodes."),
        "source": str(path), "populations": only, "n_nodes": int(sum(v.shape[0] for v in cell.positions.values())),
        "cut": args.cut, "elevation_deg": args.elevation, "radius_um": args.radius,
        "size_px": [args.width, args.height], "frames": frames,
        "encoding": {"format": args.format, "quality": args.quality if args.format != "png" else None,
                     "lossy": args.format != "png",
                     "measured_error": ("mean |delta| 1.3/255 (0.5%) at quality 82 on a native frame; "
                                        "max 163 at single-pixel filament edges; caption legible"
                                        if args.format != "png" else "none — png is lossless")},
        "render_s": round(wall, 1), "total_mb": round(total_mb, 1),
    }
    (args.out / "turntable.json").write_text(json.dumps(manifest, indent=1))
    (args.out / "index.html").write_text(_viewer_html(manifest))
    print(f"[turntable] {args.frames} frames, {total_mb:.1f} MB, {wall:.0f} s "
          f"({wall / args.frames:.1f} s/frame)", flush=True)
    print(f"[turntable] open {args.out / 'index.html'}", flush=True)
    return 0


def _save(path: Path, rgb, fmt: str, quality: int) -> None:
    """Write one frame. PNG goes through the viewer's own writer; webp through Pillow.

    ⚠ The PNG path is `cell_app._write_png` rather than Pillow, so a `--format png` turntable is
    byte-comparable with what `--headless` writes. Two encoders for one format is how two pictures of
    one thing come to differ for a reason nobody can name.
    """
    if fmt == "png":
        from aleph.viz.cell_app import _write_png
        _write_png(path, rgb)
        return
    from PIL import Image
    Image.fromarray(rgb, "RGB").save(path, "WEBP", quality=quality, method=4)


def _into_for(pos: tuple[float, float, float], cut) -> float | None:
    """How far the camera at ``pos`` sits into the half ``cut`` keeps [µm], or None with no cut.

    ⚠ The app's :func:`camera_is_on_the_kept_side` takes a NAMED view, and a turntable has none — its
    camera is a continuum. Same arithmetic, taken from the position directly, so the two cannot give
    different answers about the same geometry.
    """
    if cut is None:
        return None
    axis, offset, sign = cut
    return (pos[axis] - offset) * sign


def _viewer_html(manifest: dict) -> str:
    """The page. One file, no framework, no network, and no renderer.

    Drag or use the arrow keys to turn; every frame is a static PNG that was rendered once. The page
    states what it is and what it is not, because the picture will outlive the terminal that made it.
    """
    data = json.dumps(manifest["frames"])
    cut = manifest["cut"] or "none"
    e = manifest["encoding"]
    enc = (f'<span title="{e["measured_error"]}">({e["format"]} q{e["quality"]}, lossy)</span>'
           if e["lossy"] else "(png, lossless)")
    return f"""<!doctype html>
<meta charset="utf-8">
<title>Aleph — cell turntable</title>
<style>
 :root {{ color-scheme: light dark; }}
 body {{ margin:0; font:14px/1.5 ui-sans-serif,system-ui,sans-serif;
        background:#f7f6f4; color:#222; display:flex; flex-direction:column; height:100vh; }}
 @media (prefers-color-scheme: dark) {{ body {{ background:#15161a; color:#e7e7e7; }} }}
 header {{ padding:.6rem 1rem; border-bottom:1px solid #8883; }}
 h1 {{ font-size:1rem; margin:0 0 .15rem; }}
 .meta {{ font-size:.78rem; opacity:.75; }}
 .warn {{ color:#b02020; font-weight:600; }}
 #stage {{ flex:1; display:flex; align-items:center; justify-content:center; overflow:hidden;
           cursor:grab; touch-action:none; }}
 #stage.drag {{ cursor:grabbing; }}
 img {{ max-width:100%; max-height:100%; object-fit:contain; user-select:none; -webkit-user-drag:none; }}
 footer {{ padding:.5rem 1rem; border-top:1px solid #8883; font-size:.78rem; display:flex;
           gap:1rem; align-items:center; flex-wrap:wrap; }}
 input[type=range] {{ flex:1; min-width:180px; }}
</style>
<header>
  <h1>Aleph — cell turntable</h1>
  <div class="meta">
    {manifest["n_nodes"]:,} nodes · {len(manifest["frames"])} pre-rendered angles ·
    cut {cut} · {manifest["total_mb"]} MB {enc} ·
    <span class="warn">geometry only — no force is evaluated, nothing has moved,
    and a picture is not evidence</span>
  </div>
</header>
<div id="stage"><img id="f" alt="cell"></div>
<footer>
  <span id="ang"></span>
  <input type="range" id="sl" min="0" max="{len(manifest["frames"]) - 1}" value="0">
  <span id="note"></span>
</footer>
<script>
const F = {data};
const img = document.getElementById('f'), sl = document.getElementById('sl');
const ang = document.getElementById('ang'), note = document.getElementById('note');
let i = 0;
// ⚠ A WINDOW, NOT THE WHOLE SET. The first version created an Image for every frame up front, which
// made the angle count a MEMORY budget: 360 angles would have been ~360 decoded bitmaps resident at
// once, hundreds of MB, for a viewer showing one. Only the neighbourhood you can reach in the next
// moment is prefetched, so the number of angles costs disk and render time and nothing else.
const WINDOW = 8, KEEP = 40;
const cache = new Map();
function warm(k) {{
  for (let d = -WINDOW; d <= WINDOW; d++) {{
    const j = ((k + d) % F.length + F.length) % F.length;
    if (!cache.has(j)) {{ const im = new Image(); im.src = F[j].file; cache.set(j, im); }}
  }}
  // Bounded, and evicting the OLDEST insertion rather than the farthest: a drag walks forward, so
  // insertion order already tracks distance and the cheap rule is the right one.
  while (cache.size > KEEP) cache.delete(cache.keys().next().value);
}}
function show(k) {{
  i = ((k % F.length) + F.length) % F.length;
  img.src = F[i].file; sl.value = i; warm(i);
  ang.textContent = F[i].azimuth_deg.toFixed(0) + '\\u00b0';
  note.innerHTML = F[i].cut_faces_away
    ? '<span class="warn">\\u26a0 the cut faces away from this angle \\u2014 this looks like an UNCUT cell</span>'
    : '';
}}
show(0);
sl.addEventListener('input', e => show(+e.target.value));
const stage = document.getElementById('stage');
let dragging = false, x0 = 0, i0 = 0;
stage.addEventListener('pointerdown', e => {{
  dragging = true; x0 = e.clientX; i0 = i; stage.classList.add('drag');
  stage.setPointerCapture(e.pointerId);
}});
stage.addEventListener('pointermove', e => {{
  if (!dragging) return;
  // One full turn per stage width, so a drag across the picture is a lap.
  show(i0 + Math.round((e.clientX - x0) / stage.clientWidth * F.length));
}});
for (const ev of ['pointerup', 'pointercancel'])
  stage.addEventListener(ev, () => {{ dragging = false; stage.classList.remove('drag'); }});
addEventListener('keydown', e => {{
  if (e.key === 'ArrowLeft') show(i - 1);
  if (e.key === 'ArrowRight') show(i + 1);
}});
</script>
"""


if __name__ == "__main__":
    raise SystemExit(main())
