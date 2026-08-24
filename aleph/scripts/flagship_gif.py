#!/usr/bin/env python3
"""Extract animated GIFs (per view: shape / stress / strain / …) from a flagship viewer HTML.

The native-resolution FF/DCM morph viewers (``*_morph.html``) are 100–500 MB single files —
far too large for GitHub Pages. This driver renders them in a REAL headless Chrome, drives the
viewer's own JS (``showScene`` / ``setFrame`` / OrbitControls camera) over the Chrome DevTools
Protocol, screenshots every frame, and assembles per-scene GIFs with ffmpeg (lanczos + palette).

Why CDP and not ``browser_check.py``'s one-shot ``--screenshot``: a GIF needs many frames of the
SAME loaded document. Re-launching Chrome per frame would re-parse the 500 MB payload every time.
So we load once, keep a persistent CDP websocket, and step the viewer in place.

Leak-proofing (same discipline as ``browser_check.py``, the 2026-07-08 CPU-leak fix): Chrome runs
in its own process group with a scratch ``--user-data-dir``; the whole group is SIGKILLed on exit,
on a hard deadline, and via ``atexit`` — Chrome-149 ``--headless=new`` lingers after work, so we
never wait for it to self-exit. stderr goes to a FILE, never a pipe (a full 64 KB pipe deadlocks
Chrome). The scratch profile is always removed.

Style contract (matches the ESTABLISHED flagship GIFs in the gh-pages gallery, 2026-07-10): the
on-canvas control panel (``#ui``) and object legend (``#leg``) are HIDDEN, the turbo colorbar
(``#cbar``) is KEPT, the camera is dollied in so the cell fills the frame, and the output is
760 px wide at 150 ms/frame — the same look as ``ff_afm_*.gif`` / ``ff_mt_*.gif``.

Viewer contract (both ``ff_viewer_html.py`` and ``dcm_mesh_viewer_html.py`` satisfy it): global
``P`` (payload; ``P.scenes`` object, ``P.n_frames`` int), ``renderer`` / ``scene`` / ``cam`` /
``ctr`` (Three.js r128 + OrbitControls), and functions ``showScene(name)`` / ``setFrame(i)``.
Missing pieces degrade gracefully (no ``showScene`` → single view; ``n_frames<=1`` → turntable).

Usage
-----
    # one GIF per scene, animated over the morph frames (front view):
    python aleph/scripts/flagship_gif.py outputs/ff/figs/afm_fullcompartment_morph.html \
        --out-dir /tmp/gifs --prefix ff_afm

    # turntable (rotate camera) instead of frame-stepping, e.g. for a static confluent spheroid:
    python aleph/scripts/flagship_gif.py outputs/.../N2000_confluent...FULL.html \
        --out-dir /tmp/gifs --prefix dcm_n2000 --mode turntable --scenes stress,junction

Exit 0 = at least one GIF written. Non-zero (message on stderr) otherwise.
"""
from __future__ import annotations

import argparse
import atexit
import base64
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path
from typing import Any

_CHROME_CANDIDATES = (
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Chromium.app/Contents/MacOS/Chromium",
    "/Applications/Google Chrome Canary.app/Contents/MacOS/Google Chrome Canary",
    "google-chrome",
    "chromium",
    "chromium-browser",
)


def find_chrome() -> str:
    """Return a usable Chrome/Chromium binary path or raise (``$CHROME_BIN`` wins)."""
    env = os.environ.get("CHROME_BIN")
    if env:
        if Path(env).exists() or shutil.which(env):
            return env
        raise FileNotFoundError(f"$CHROME_BIN={env!r} does not exist")
    for cand in _CHROME_CANDIDATES:
        hit = cand if Path(cand).exists() else shutil.which(cand)
        if hit:
            return hit
    raise FileNotFoundError("No Chrome/Chromium found; set $CHROME_BIN.")


# ── injected once into the page: set scene + frame + camera azimuth, then force a
#    SYNCHRONOUS render so the very next CDP screenshot captures this exact state
#    (the viewer's rAF loop is event-driven/idle, so we cannot rely on it firing). ──
_CAP_JS = r"""
window.__cap = function(sc, frac, azi){
  try {
    if (sc && typeof showScene === 'function') showScene(sc);
    if (typeof P !== 'undefined' && P.n_frames > 1 && typeof setFrame === 'function')
      setFrame(Math.round(frac * (P.n_frames - 1)));
    if (azi !== null) {
      var t = ctr.target, cp = cam.position;
      var dx = cp.x - t.x, dy = cp.y - t.y, dz = cp.z - t.z;
      var r = Math.sqrt(dx*dx + dy*dy + dz*dz);
      var el = Math.asin(Math.max(-1, Math.min(1, dy / r)));   // keep the viewer's own elevation
      var rc = r * Math.cos(el);
      cam.position.set(t.x + rc*Math.sin(azi), t.y + r*Math.sin(el), t.z + rc*Math.cos(azi));
      cam.up.set(0,1,0); cam.lookAt(t); cam.updateProjectionMatrix();
    }
    renderer.render(scene, cam);
    return { ok: true, n: (typeof P !== 'undefined' ? P.n_frames : 1) };
  } catch (e) { return { ok: false, err: String(e) }; }
};
// hide the on-canvas control panel + object legend, KEEP the colorbar — the flagship look.
window.__style = function(hideIds){
  hideIds.forEach(function(id){ var e = document.getElementById(id); if (e) e.style.display='none'; });
  return true;
};
// dolly the camera toward its target so the cell fills the frame (factor<1 = closer). One-shot.
window.__zoom = function(f){
  var t = ctr.target, cp = cam.position;
  cam.position.set(t.x + (cp.x-t.x)*f, t.y + (cp.y-t.y)*f, t.z + (cp.z-t.z)*f);
  cam.lookAt(t); cam.updateProjectionMatrix(); return true;
};
true
"""


class Chrome:
    """A persistent headless-Chrome + CDP websocket session with guaranteed teardown."""

    def __init__(self, url: str, window: tuple[int, int] = (1100, 900)) -> None:
        import websocket  # websocket-client; imported here so --help works without it

        self._ws_mod = websocket
        self.profile = Path(tempfile.mkdtemp(prefix="flgif_"))
        self.err_path = Path(tempfile.mkstemp(prefix="flgif_", suffix=".log")[1])
        self.proc: subprocess.Popen | None = None
        self.ws: Any = None
        self._id = 0
        atexit.register(self.close)
        self._launch(url, window)

    def _launch(self, url: str, window: tuple[int, int]) -> None:
        args = [
            find_chrome(), "--headless=new", "--no-sandbox",
            f"--user-data-dir={self.profile}", "--allow-file-access-from-files",
            "--use-gl=angle", "--use-angle=swiftshader", "--hide-scrollbars",
            "--force-color-profile=srgb", "--disable-gpu-vsync",
            f"--window-size={window[0]},{window[1]}",
            "--remote-debugging-port=0", "--remote-allow-origins=*", url,
        ]
        err_fh = open(self.err_path, "wb")
        try:
            self.proc = subprocess.Popen(
                args, stdout=subprocess.DEVNULL, stderr=err_fh, start_new_session=True,
            )
        finally:
            err_fh.close()
        port = self._await_port()
        ws_url = self._page_ws(port)
        self.ws = self._ws_mod.create_connection(ws_url, max_size=None, timeout=120)
        self.call("Page.enable")
        self.call("Runtime.enable")

    def _await_port(self, timeout_s: float = 30.0) -> int:
        f = self.profile / "DevToolsActivePort"
        deadline = time.monotonic() + timeout_s
        while time.monotonic() < deadline:
            if self.proc and self.proc.poll() is not None:
                raise RuntimeError(f"Chrome exited early; see {self.err_path}")
            if f.exists():
                txt = f.read_text().splitlines()
                if txt and txt[0].strip().isdigit():
                    return int(txt[0].strip())
            time.sleep(0.1)
        raise TimeoutError("Chrome never wrote DevToolsActivePort")

    def _page_ws(self, port: int, timeout_s: float = 20.0) -> str:
        deadline = time.monotonic() + timeout_s
        while time.monotonic() < deadline:
            try:
                data = json.loads(
                    urllib.request.urlopen(f"http://127.0.0.1:{port}/json", timeout=5).read()
                )
                for t in data:
                    if t.get("type") == "page" and t.get("webSocketDebuggerUrl"):
                        return t["webSocketDebuggerUrl"]
            except Exception:
                pass
            time.sleep(0.2)
        raise TimeoutError("no CDP page target appeared")

    def call(self, method: str, params: dict | None = None, timeout: float = 120.0) -> dict:
        self._id += 1
        mid = self._id
        self.ws.settimeout(timeout)
        self.ws.send(json.dumps({"id": mid, "method": method, "params": params or {}}))
        while True:
            msg = json.loads(self.ws.recv())
            if msg.get("id") == mid:
                if "error" in msg:
                    raise RuntimeError(f"{method}: {msg['error']}")
                return msg.get("result", {})
            # otherwise it's an event; ignore

    def evaluate(self, expr: str, timeout: float = 120.0) -> Any:
        r = self.call(
            "Runtime.evaluate",
            {"expression": expr, "returnByValue": True, "awaitPromise": True},
            timeout=timeout,
        )
        if r.get("exceptionDetails"):
            raise RuntimeError(f"JS: {r['exceptionDetails'].get('text')}")
        return r.get("result", {}).get("value")

    def screenshot(self) -> bytes:
        r = self.call("Page.captureScreenshot", {"format": "png", "captureBeyondViewport": False})
        return base64.b64decode(r["data"])

    def close(self) -> None:
        try:
            if self.ws is not None:
                self.ws.close()
        except Exception:
            pass
        self.ws = None
        if self.proc is not None and self.proc.poll() is None:
            try:
                os.killpg(os.getpgid(self.proc.pid), signal.SIGKILL)
            except Exception:
                self.proc.kill()
        self.proc = None
        shutil.rmtree(self.profile, ignore_errors=True)
        self.err_path.unlink(missing_ok=True)


def _slug(s: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", s).strip("_").lower() or "view"


def _classify(name: str) -> str:
    """Map a viewer scene name to a friendly gallery key (shape / stress / strain / …)."""
    n = name.lower()
    if "strain" in n:
        return "strain"
    if "stress" in n or "vm" in n or "σ" in n or "sigma" in n or "von" in n:
        return "stress"
    if "junction" in n or "cadher" in n:
        return "junction"
    if any(k in n for k in ("shape", "morph", "geom", "cell", "native", "cortex")):
        return "shape"
    return _slug(name)


def _assemble_gif(frame_dir: Path, out_gif: Path, width: int, duration_ms: int) -> bool:
    """Assemble ``frame_%04d.png`` in ``frame_dir`` into an optimized GIF; ffmpeg then PIL."""
    framerate = 1000.0 / max(1, duration_ms)
    if shutil.which("ffmpeg"):
        pal = frame_dir / "_palette.png"
        vf = f"scale={width}:-1:flags=lanczos"
        try:
            subprocess.run(
                ["ffmpeg", "-y", "-v", "error", "-framerate", f"{framerate:.4f}",
                 "-i", str(frame_dir / "frame_%04d.png"),
                 "-vf", f"{vf},palettegen=stats_mode=diff", str(pal)],
                check=True,
            )
            subprocess.run(
                ["ffmpeg", "-y", "-v", "error", "-framerate", f"{framerate:.4f}",
                 "-i", str(frame_dir / "frame_%04d.png"), "-i", str(pal),
                 "-filter_complex", f"{vf}[x];[x][1:v]paletteuse=dither=bayer:bayer_scale=3",
                 "-loop", "0", str(out_gif)],
                check=True,
            )
            return out_gif.exists()
        except subprocess.CalledProcessError:
            pass
    # PIL fallback
    try:
        from PIL import Image
    except Exception:
        return False
    pngs = sorted(frame_dir.glob("frame_*.png"))
    if not pngs:
        return False
    imgs = []
    for p in pngs:
        im = Image.open(p).convert("RGB")
        h = int(im.height * width / im.width)
        imgs.append(im.resize((width, h), Image.LANCZOS).quantize(colors=128, method=Image.MEDIANCUT))
    imgs[0].save(out_gif, save_all=True, append_images=imgs[1:],
                 duration=duration_ms, loop=0, optimize=True, disposal=2)
    return out_gif.exists()


def extract(
    html: Path, out_dir: Path, prefix: str, mode: str, scenes_filter: list[str] | None,
    n_frames: int, width: int, duration_ms: int, ready_timeout: float, window: tuple[int, int],
    still: bool, zoom: float, hide_ui: bool,
) -> list[Path]:
    """Render ``html`` and write one GIF per selected scene into ``out_dir``. Returns the GIFs."""
    out_dir.mkdir(parents=True, exist_ok=True)
    ch = Chrome(html.resolve().as_uri(), window=window)
    written: list[Path] = []
    try:
        # wait for the viewer to finish parsing the (possibly 500 MB) payload + first build
        deadline = time.monotonic() + ready_timeout
        while time.monotonic() < deadline:
            ready = ch.evaluate(
                "typeof P!=='undefined' && typeof renderer!=='undefined' && "
                "typeof scene!=='undefined' && typeof cam!=='undefined'"
            )
            if ready:
                break
            time.sleep(0.5)
        else:
            raise TimeoutError(f"viewer never became ready within {ready_timeout:.0f}s")
        ch.evaluate(_CAP_JS)
        if hide_ui:  # flagship look: drop the control panel + legend, keep the colorbar
            ch.evaluate("__style(['ui','leg','title','hud','panel','legend','help','info'])")
        if zoom and zoom != 1.0:  # dolly in so the cell fills the frame
            ch.evaluate(f"__zoom({zoom})")

        scene_names: list[str] = ch.evaluate(
            "(typeof P!=='undefined' && P.scenes) ? Object.keys(P.scenes) : []"
        ) or [""]
        has_show = bool(ch.evaluate("typeof showScene==='function'"))
        n_data = int(ch.evaluate("(typeof P!=='undefined' && P.n_frames) ? P.n_frames : 1") or 1)
        if not has_show:
            scene_names = [""]  # single view; camera/frames only

        # choose which scenes to render
        chosen: list[tuple[str, str]] = []  # (scene_name, gallery_key)
        for sn in scene_names:
            key = _classify(sn) if sn else _slug(prefix)
            if scenes_filter and not any(f.lower() in (sn.lower() + " " + key) for f in scenes_filter):
                continue
            chosen.append((sn, key))
        if not chosen:
            chosen = [(scene_names[0], _classify(scene_names[0]) if scene_names[0] else _slug(prefix))]

        use_turntable = mode == "turntable" or (mode == "auto" and n_data <= 1)
        print(f"[flagship_gif] {html.name}: scenes={scene_names} n_frames={n_data} "
              f"mode={'turntable' if use_turntable else 'frames'} → {[k for _,k in chosen]}")

        seen_keys: dict[str, int] = {}
        for sname, key in chosen:
            if key in seen_keys:
                seen_keys[key] += 1
                key = f"{key}{seen_keys[key]}"
            else:
                seen_keys[key] = 0
            frame_dir = Path(tempfile.mkdtemp(prefix="flgif_frames_"))
            try:
                for i in range(n_frames):
                    frac = i / max(1, n_frames - 1)
                    if use_turntable:
                        azi = 2 * 3.141592653589793 * i / n_frames  # full 360° orbit
                        r = ch.evaluate(f"__cap({json.dumps(sname)}, 0.0, {azi})")
                    else:
                        r = ch.evaluate(f"__cap({json.dumps(sname)}, {frac}, null)")
                    if not (isinstance(r, dict) and r.get("ok")):
                        raise RuntimeError(f"__cap failed on {sname!r}: {r}")
                    (frame_dir / f"frame_{i:04d}.png").write_bytes(ch.screenshot())
                out_gif = out_dir / f"{prefix}_{key}.gif"
                if _assemble_gif(frame_dir, out_gif, width, duration_ms):
                    written.append(out_gif)
                    kb = out_gif.stat().st_size / 1024
                    print(f"  ✓ {out_gif.name}  ({kb:.0f} KB, {n_frames}f)")
                    if still:  # representative mid still
                        mid = frame_dir / f"frame_{n_frames//2:04d}.png"
                        if mid.exists():
                            shutil.copy(mid, out_dir / f"{prefix}_{key}.png")
                else:
                    print(f"  ✗ GIF assembly failed for {key}", file=sys.stderr)
            finally:
                shutil.rmtree(frame_dir, ignore_errors=True)
    finally:
        ch.close()
    return written


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("html", type=Path, help="flagship viewer HTML")
    ap.add_argument("--out-dir", type=Path, required=True, help="where GIFs are written")
    ap.add_argument("--prefix", required=True, help="output basename prefix (e.g. ff_afm)")
    ap.add_argument("--mode", choices=["auto", "frames", "turntable"], default="auto",
                    help="auto: frame-step if animated else turntable")
    ap.add_argument("--scenes", default=None,
                    help="comma-list substring filter of scenes/keys to render (default: all)")
    ap.add_argument("--frames", type=int, default=18, help="frames in each GIF (flagship default 18)")
    ap.add_argument("--width", type=int, default=760, help="GIF width px (flagship default 760)")
    ap.add_argument("--duration-ms", type=int, default=150, help="ms per frame (flagship default 150)")
    ap.add_argument("--zoom", type=float, default=0.72, help="camera dolly factor (<1 = closer/fills frame)")
    ap.add_argument("--keep-ui", action="store_true", help="keep the control panel + legend (default hides them)")
    ap.add_argument("--ready-timeout", type=float, default=300.0, help="max s to load the payload")
    ap.add_argument("--window", default="960x960", help="framebuffer WxH (square fills like the flagships)")
    ap.add_argument("--no-still", action="store_true", help="do not also write a mid-frame PNG")
    args = ap.parse_args()

    if not args.html.exists():
        print(f"[FAIL] no such file: {args.html}", file=sys.stderr)
        return 6
    try:
        w, h = (int(x) for x in args.window.lower().split("x"))
    except Exception:
        print(f"[FAIL] bad --window {args.window!r} (want WxH)", file=sys.stderr)
        return 6
    scenes = [s.strip() for s in args.scenes.split(",")] if args.scenes else None
    try:
        gifs = extract(args.html, args.out_dir, args.prefix, args.mode, scenes,
                       args.frames, args.width, args.duration_ms, args.ready_timeout, (w, h),
                       still=not args.no_still, zoom=args.zoom, hide_ui=not args.keep_ui)
    except Exception as e:
        print(f"[FAIL] {type(e).__name__}: {e}", file=sys.stderr)
        return 1
    if not gifs:
        print("[FAIL] no GIF written", file=sys.stderr)
        return 2
    print(f"[OK] {len(gifs)} GIF(s) → {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
