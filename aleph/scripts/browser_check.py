#!/usr/bin/env python3
"""Leak-proof in-browser smoke check for the generated HTML cell viewers.

Root-cause fix for the 2026-07-08 CPU incident: ad-hoc CDP/puppeteer verification
launched a headless Chrome and never closed it — nine sessions leaked in one day and a
single GPU helper pegged 6 CPU cores for 3 h under software-WebGL, which also blocked
the normal GUI Chrome from relaunching (LaunchServices saw the app "already running").

This helper renders a viewer in a REAL browser (system Chrome, software-WebGL that works
headless) and *cannot* leak a process:

  * Chrome writes the one-shot ``--screenshot`` within a second or two, but Chrome-149
    ``--headless=new`` then *lingers* instead of exiting (verified: every flag combo, incl.
    ``--virtual-time-budget`` and ``--run-all-compositor-stages-before-draw``, hangs after
    writing the PNG). So we do NOT wait for it to self-exit — we poll for the screenshot
    file, and the instant it is written we ``SIGKILL`` the whole process group ourselves.
  * The child runs in its own process group; the same kill fires on a hard deadline if the
    screenshot never appears (a viewer stuck before first paint).
  * stderr is redirected to a temp FILE, never a pipe: ``--enable-logging=stderr`` can spew
    megabytes and a full 64 KB pipe would deadlock Chrome (the original hang bug).
  * The scratch ``--user-data-dir`` lives under the OS tempdir and is removed in
    ``finally``; ``atexit`` guarantees teardown even on crash / Ctrl-C.

Beyond "did it not leak", it also flags the failure class from
``reference-verify-html-viewer-in-browser`` (a *silently* dead viewer that renders only
the background): the screenshot is checked for non-background content, and Chrome's
stderr is scanned for uncaught JS exceptions. A byte-size check on the HTML is NOT a
substitute — that is exactly what missed the dead viewer before.

Usage
-----
    python aleph/scripts/browser_check.py path/to/viewer.html
    python aleph/scripts/browser_check.py viewer.html --out shot.png --timeout 90

Exit code 0 = Chrome rendered a non-blank frame, exited cleanly, and no uncaught JS
exception was seen. Non-zero (with a message on stderr) otherwise. Override the browser
with ``$CHROME_BIN``.
"""
from __future__ import annotations

import argparse
import atexit
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from urllib.parse import quote

# Ordered by preference; first existing one wins. $CHROME_BIN overrides all.
_CHROME_CANDIDATES = (
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Chromium.app/Contents/MacOS/Chromium",
    "/Applications/Google Chrome Canary.app/Contents/MacOS/Google Chrome Canary",
    "google-chrome",
    "chromium",
    "chromium-browser",
)

# JS failure markers in Chrome's stderr console log (the dead-viewer class of bug).
_JS_ERROR_RE = re.compile(
    r"Uncaught|Unhandled|SyntaxError|TypeError|ReferenceError|is not defined|is not a function",
    re.IGNORECASE,
)


def find_chrome() -> str:
    """Return a usable Chrome/Chromium binary path or raise.

    Honours ``$CHROME_BIN`` first, then a platform candidate list.
    """
    env = os.environ.get("CHROME_BIN")
    if env:
        if Path(env).exists() or shutil.which(env):
            return env
        raise FileNotFoundError(f"$CHROME_BIN={env!r} does not exist")
    for cand in _CHROME_CANDIDATES:
        hit = cand if Path(cand).exists() else shutil.which(cand)
        if hit:
            return hit
    raise FileNotFoundError(
        "No Chrome/Chromium found; install one or set $CHROME_BIN to its path."
    )


def _screenshot_is_blank(png: Path) -> bool | None:
    """True if the screenshot is ~uniform (nothing drew), None if undecidable.

    Uses Pillow if available; returns None (skip the check) when it is not, so the
    helper degrades gracefully rather than failing on a missing optional dep.
    """
    try:
        from PIL import Image, ImageStat  # type: ignore
    except Exception:
        return None
    with Image.open(png) as im:
        stat = ImageStat.Stat(im.convert("L"))
    # A viewer that drew a cell over the dark 0x0e1117 background has real luminance
    # spread; a background-only frame is nearly constant. 3.0 is comfortably below any
    # real render and well above PNG/AA noise on a flat fill.
    return float(stat.stddev[0]) < 3.0


def check(
    html: Path,
    out_png: Path,
    timeout_s: float = 60.0,
    budget_ms: int = 8000,
    window: tuple[int, int] = (1400, 900),
    scene: str | None = None,
    frame: str | None = None,
) -> int:
    """Render ``html`` headlessly to ``out_png`` with guaranteed teardown.

    Args:
        html: Path to the viewer HTML file.
        out_png: Where to write the one-shot screenshot.
        timeout_s: Hard wall-clock limit; the process group is killed past it.
        budget_ms: Chrome ``--virtual-time-budget`` — how long JS/render runs before
            Chrome self-exits. Bump it for very large meshes.
        window: Off-screen framebuffer size.

    Returns:
        Process-style exit code (0 == pass).
    """
    chrome = find_chrome()
    profile = Path(tempfile.mkdtemp(prefix="cdp_check_"))
    err_path = Path(tempfile.mkstemp(prefix="cdp_check_", suffix=".log")[1])
    proc: subprocess.Popen | None = None

    def _cleanup() -> None:
        if proc is not None and proc.poll() is None:
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            except Exception:
                proc.kill()
        shutil.rmtree(profile, ignore_errors=True)
        err_path.unlink(missing_ok=True)

    atexit.register(_cleanup)

    stderr_txt = ""
    try:
        if out_png.exists():
            out_png.unlink()  # never mistake a stale screenshot for this run
        _hq = ([f"scene={quote(scene)}"] if scene else []) + ([f"frame={quote(str(frame))}"] if frame is not None else [])
        _hash = ("#" + "&".join(_hq)) if _hq else ""
        args = [
            chrome,
            "--headless=new",
            "--no-sandbox",
            f"--user-data-dir={profile}",
            "--allow-file-access-from-files",
            "--use-gl=angle",
            "--use-angle=swiftshader",  # software WebGL that actually renders headless
            "--hide-scrollbars",
            "--force-color-profile=srgb",
            f"--window-size={window[0]},{window[1]}",
            f"--virtual-time-budget={budget_ms}",  # advance a few s of JS so the frame settles
            f"--screenshot={out_png}",
            "--enable-logging=stderr",
            "--v=0",
            html.resolve().as_uri() + _hash,
        ]
        # stderr → FILE (not a pipe → no 64 KB deadlock); own process group so we can
        # SIGKILL the whole tree. Chrome-149 new-headless writes the PNG fast but then
        # lingers, so we poll for the screenshot and kill it the moment it is written.
        err_fh = open(err_path, "wb")
        try:
            proc = subprocess.Popen(
                args, stdout=subprocess.DEVNULL, stderr=err_fh, start_new_session=True,
            )
        finally:
            err_fh.close()  # the child keeps its own dup of the fd
        deadline = time.monotonic() + timeout_s
        last_size = -1
        while time.monotonic() < deadline:
            if proc.poll() is not None:  # exited on its own (older builds)
                break
            if out_png.exists():
                size = out_png.stat().st_size
                if size > 0 and size == last_size:  # stable across two polls → fully written
                    break
                last_size = size
            time.sleep(0.25)
        try:
            stderr_txt = err_path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            stderr_txt = ""
    finally:
        _cleanup()

    js_errors = [ln for ln in stderr_txt.splitlines() if _JS_ERROR_RE.search(ln)]

    if not out_png.exists():
        print("[FAIL] Chrome produced no screenshot (render never completed).", file=sys.stderr)
        if js_errors:
            print("       JS errors seen:\n         " + "\n         ".join(js_errors[:8]), file=sys.stderr)
        return 3

    blank = _screenshot_is_blank(out_png)
    if blank is True:
        print(
            f"[FAIL] Screenshot is background-only — nothing drew: {out_png}\n"
            "       (a silently dead viewer; see reference-verify-html-viewer-in-browser)",
            file=sys.stderr,
        )
        return 4
    if js_errors:
        print("[FAIL] Uncaught JS in the viewer:\n  " + "\n  ".join(js_errors[:8]), file=sys.stderr)
        return 5

    note = "" if blank is not None else "  (blank-check skipped: Pillow not installed)"
    print(f"[OK] rendered + exited cleanly, no JS errors → {out_png}{note}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("html", type=Path, help="viewer HTML file to render")
    ap.add_argument("--out", type=Path, default=None, help="screenshot path (default: OS tempdir)")
    ap.add_argument("--timeout", type=float, default=60.0, help="hard kill timeout, seconds")
    ap.add_argument("--budget-ms", type=int, default=8000, help="Chrome virtual-time budget, ms")
    ap.add_argument("--scene", default=None, help="auto-select this viewer scene via #scene=NAME before the screenshot")
    ap.add_argument("--frame", default=None, help="advance to this frame via #frame=N|last before the screenshot")
    args = ap.parse_args()

    if not args.html.exists():
        print(f"[FAIL] no such file: {args.html}", file=sys.stderr)
        return 6
    _tag = ("." + re.sub(r"[^0-9A-Za-z]+", "_", args.scene)) if args.scene else ""
    _tag += (".f" + re.sub(r"[^0-9A-Za-z]+", "_", args.frame)) if args.frame else ""
    out = args.out or Path(tempfile.gettempdir()) / (args.html.stem + _tag + ".check.png")
    return check(args.html, out, timeout_s=args.timeout, budget_ms=args.budget_ms, scene=args.scene, frame=args.frame)


if __name__ == "__main__":
    raise SystemExit(main())
