"""The native viewer's gates, run with everything else rather than when someone remembers.

`aleph/viz/cell_app.py` carries its own `_demo()` because it must be runnable on a machine with no
pytest — that is the shape `world/observe_gamma.py` adopted on 2026-08-21 after the production
environment on the GPU host turned out to have no test runner. **A self-check that only runs when
invoked by hand is a self-check that stops running**, so this delegates to it.

⚠ **What is NOT covered here, stated rather than implied.** These exercise the cut arithmetic and the
app's state machine, both of which are pure. The **GL binding is not exercised**: no window is opened,
`register_key_press_callback` is not called, and nothing here would notice if warp changed the
callback signature. That gap is why `_apply_key` was split out of the callback in the first place —
the logic is checkable everywhere and the binding is a two-line adapter — and it is a gap, not a
solved problem.
"""

from __future__ import annotations

import inspect
from pathlib import Path

from aleph.viz import cell_app


def test_the_module_self_check_passes() -> None:
    """Delegate to the module's own gate, so the two cannot drift apart."""
    cell_app._demo()


def test_the_key_callback_matches_what_warp_will_call() -> None:
    """The one thing about the GL binding that IS checkable without a display: its arity.

    Warp invokes the registered callback as ``cb(symbol, modifiers)``. A signature mismatch is a
    `TypeError` at the first keystroke — in a window, on a machine that may have no terminal in view.
    """
    src = inspect.getsource(cell_app._install_keys)
    assert "def on_key(symbol: int, modifiers: int) -> None:" in src, (
        "the callback warp calls takes (symbol, modifiers); a mismatch fails at the first keypress")
    assert "register_key_press_callback(on_key)" in src


def test_nothing_is_defined_below_the_main_block() -> None:
    """The defect that shipped in `observe_gamma.py` this morning, ratcheted here.

    Its device section sat below `if __name__ == "__main__"`, so `python -m` bound none of it and the
    documented invocation failed while both an import-based test and a subprocess test on a machine
    without CUDA reported success. Only a static check can see it.
    """
    import ast
    import pathlib

    tree = ast.parse(pathlib.Path(cell_app.__file__).read_text())
    main_line = next(n.lineno for n in ast.walk(tree)
                     if isinstance(n, ast.If) and "__main__" in ast.dump(n.test))
    after = [n.name for n in tree.body
             if isinstance(n, (ast.FunctionDef, ast.ClassDef)) and n.lineno > main_line]
    assert not after, f"defined below __main__ and unreachable from `python -m`: {after}"


def test_a_cut_that_faces_away_from_the_camera_is_detected() -> None:
    """⚠ The sign is the whole content, so both directions are pinned.

    A cut only opens the cell if the camera is on the side the cut REMOVED. `oblique` sits at
    z = +22, so `--cut z-` (keep z <= 0) opens toward it and `--cut z+` presents the closed dome —
    which is what `--cut z+ --view oblique` actually rendered on 2026-08-21: an opaque membrane cap
    and nothing else.
    """
    from aleph.viz.cell_app import camera_is_on_the_kept_side as facing

    assert facing("oblique", (2, 0.0, -1)) == -22.0     # keep z<=0: camera on the removed side, opens
    assert facing("oblique", (2, 0.0, +1)) == +22.0     # keep z>=0: camera inside the kept half, closed
    assert facing("front", None) is None                # no cut, no verdict -- not 0.0, which sorts

    # An offset plane moves the answer, or the check is reading the axis and not the plane.
    assert facing("oblique", (2, 10.0, +1)) == +12.0
    assert facing("oblique", (2, 30.0, +1)) == -8.0     # camera now OUTSIDE the kept half


def test_the_sparkline_says_where_you_are_and_nothing_else() -> None:
    """It is a position indicator, not a figure — and it must not lie about position."""
    from aleph.viz.cell_app import sparkline

    rising = [float(i) for i in range(400)]
    assert sparkline(rising).startswith("▁") and sparkline(rising).endswith("█")
    assert sparkline([]) == ""                          # no series, no picture
    assert len(sparkline(rising)) == 56

    # The mark tracks the index, at both ends and in between.
    assert sparkline(rising, 0).index("|") == 0
    assert sparkline(rising, 399).index("|") == 55
    assert 26 <= sparkline(rising, 200).index("|") <= 29

    # A flat series must not render as noise: every column equal -> every bar the same.
    assert len(set(sparkline([7.0] * 400))) == 1


def test_a_run_series_needs_a_stride_and_argv_is_a_STRING(tmp_path: Path) -> None:
    """⚠ The record stores `argv` as a string, not a list, and both read the same in Python.

    `"--snapshot-every" in argv` is a substring test on a string and `.index` is a character offset,
    so the first version pulled one letter, raised ValueError, swallowed it, and the panel silently
    never appeared. Both shapes are pinned so the split cannot be dropped.
    """
    import json

    from aleph.viz.cell_app import load_run_series

    cellfile = tmp_path / "run.alephcell"
    cellfile.write_bytes(b"")
    rec = tmp_path / "run.json"

    def write(argv):
        rec.write_text(json.dumps({"gamma": {"trace_pn_per_um": [1.0, 2.0, 3.0]},
                                   "provenance": {"argv": argv}}))

    write("driver.py --steps 30000 --snapshot-every 6000 --out x.json")
    assert load_run_series(cellfile) == ([1.0, 2.0, 3.0], 6000)

    write(["driver.py", "--steps", "30000", "--snapshot-every", "6000"])
    assert load_run_series(cellfile) == ([1.0, 2.0, 3.0], 6000)

    # No stride recorded -> None, NOT a guessed 1. Indexing a trace by frame number when the frames
    # are 6,000 steps apart points at the wrong place with total confidence.
    write("driver.py --steps 30000")
    assert load_run_series(cellfile) is None

    # No sibling record at all, and a record with no trace.
    assert load_run_series(tmp_path / "absent.alephcell") is None
    rec.write_text(json.dumps({"provenance": {"argv": "driver.py --snapshot-every 6000"}}))
    assert load_run_series(cellfile) is None


def test_describe_cell_reports_what_is_ABSENT(tmp_path: Path) -> None:
    """⚠ The present list is the easy half and is not what this function is for.

    A census of eleven names does not tell a reader it is missing a twelfth — that is how the PHASE 4
    tau runs were stepped and read for a night before anyone noticed the cell has no myosin. So the
    test pins the ABSENT list, and pins that the reference is derived from the builder modules rather
    than typed.
    """
    import numpy as np

    from aleph.viz.cell_app import CellFile, describe_cell

    cell = CellFile(header={}, positions={"cortex": np.zeros((3, 3)), "membrane": np.zeros((2, 3))},
                    segments={}, faces={}, frames={}, path=tmp_path / "x.alephcell")
    out = describe_cell(cell)
    assert "HAS (2): cortex, membrane" in out
    assert "NOT IN THIS FILE" in out
    assert "nmii" in out.split("NOT IN THIS FILE")[1]
    # An absence is reported, never accused.
    assert "Not an accusation" in out
    # No sibling record must say so rather than pass silently.
    assert "provenance unknown, which is not the same as provenance clean" in out


def test_describe_cell_says_so_when_nothing_is_missing(tmp_path: Path) -> None:
    """The all-present branch must exist and be reachable, or the absent branch is untested by half."""
    import numpy as np

    from aleph.viz.cell_app import CellFile, describe_cell

    build_dir = Path(__file__).resolve().parents[2] / "world" / "build"
    stems = [f.stem for f in build_dir.glob("*.py") if not f.stem.startswith("_")]
    assert stems, "no builder modules found — the reference this function derives is empty"
    cell = CellFile(header={}, positions={s: np.zeros((1, 3)) for s in stems},
                    segments={}, faces={}, frames={}, path=None)
    out = describe_cell(cell)
    assert "every one of the" in out and "NOT IN THIS FILE" not in out


def test_list_exports_answers_what_is_there_without_a_path(tmp_path: Path) -> None:
    """⚠ An empty listing must be a statement about the directory, not about the world."""
    from aleph.viz.cell_app import list_exports

    out = list_exports(tmp_path)
    assert "no .alephcell under" in out
    assert "not about whether any cell has ever been exported" in out

    # A file that is not a valid export is LISTED as unreadable, never hidden. A listing that
    # silently drops a broken file tells you it does not exist.
    (tmp_path / "broken.alephcell").write_bytes(b"not an export")
    out = list_exports(tmp_path)
    assert "broken.alephcell" in out and "unreadable" in out


def test_the_caption_reaches_the_SAVED_PIXELS_and_not_only_the_window(tmp_path: Path) -> None:
    """⚠ The one gate that would have caught every way the HUD failed, and it failed three ways.

    This file's own header says the GL binding is not exercised. That gap is exactly where the HUD
    lives, and while it was being written it produced a render that succeeded and showed nothing
    **three times**, for three unrelated reasons: the caption went to framebuffer 0 while
    `get_pixels` reads `_frame_texture`; the renderer's depth test discarded it; and a flipped
    projection turned its own glyph quads into culled back faces. Every one of those printed
    `wrote ... 3200x2000` and wrote a clean, caption-free file. **Only opening the PNG showed it**,
    so this opens the PNG.

    It asserts INK IN THE CAPTION BAND, not the strings — the strings were right every time. A test
    on `caption` would have passed on all three broken renders.

    Skipped, with a reason, where no GL context can be made: the GPU host runs inside a Slurm
    allocation with no display. A skip that says why is not the same as a pass.
    """
    import numpy as np
    import pytest
    from PIL import Image

    from aleph.viz.cell_app import CellFile, _render

    cell = CellFile(
        header={"thinning": "NONE"},
        positions={"membrane": np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])},
        segments={"membrane": np.array([[0, 1], [1, 2]], dtype=np.int64)},
        faces={}, frames={}, path=tmp_path / "tiny.alephcell")
    out = tmp_path / "shot.png"
    try:
        _render(cell, ["membrane"], out, thicken=1.0, view="oblique")
    except Exception as exc:                       # noqa: BLE001 — no display is a SKIP, not a fail
        pytest.skip(f"no GL context here: {type(exc).__name__}: {exc}")

    img = np.asarray(Image.open(out).convert("L"))
    band = img[int(img.shape[0] * 0.88):, :]       # the caption sits along the bottom edge
    # The background is near-white (0.97, 0.96, 0.94). Text is dark; anything below this is ink and
    # nothing else in the scene reaches down here on a three-node cell.
    assert (band < 128).sum() > 200, (
        f"the caption band carries {(band < 128).sum()} dark pixels — the HUD drew nothing into the "
        f"file that leaves this machine, which is the failure it exists to prevent")


def test_the_caption_leaves_no_GL_STATE_BEHIND_for_the_next_frame(tmp_path: Path) -> None:
    """⚠ The caption is drawn with depth testing and culling OFF. It must put them back.

    GL state is global. While the HUD was being written it disabled both and never restored them, so
    everything drawn AFTER a caption drew with no depth test. In a one-frame render nothing is drawn
    after, which is why `--headless` and every figure taken with it were correct and this was
    invisible: the bug needs a SECOND frame to show. On the second frame the nucleus rendered THROUGH
    the membrane and the mean colour went orange to pink, while every count in the caption stayed
    exactly right — a plausible picture, which is the worst kind of wrong.

    So the check renders the SAME camera three times and demands the frames be bit-identical. Not
    "similar": identical. An earlier version of this compared three angles of a rotationally
    symmetric view and had to guess a tolerance, which cannot distinguish a state leak from the
    lighting moving with the camera — and did not.
    """
    import numpy as np
    import pytest
    import warp as wp
    import warp.render

    from aleph.viz.cell_app import CellFile, _attach_hud, _caption, _draw_population

    cell = CellFile(
        header={"thinning": "NONE"},
        positions={"membrane": np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])},
        segments={"membrane": np.array([[0, 1], [1, 2]], dtype=np.int64)},
        faces={}, frames={}, path=tmp_path / "tiny.alephcell")
    try:
        wp.init()
        r = warp.render.OpenGLRenderer(screen_width=320, screen_height=240, up_axis="Y",
                                       headless=True, camera_pos=(0.0, 1.0, 6.0),
                                       camera_front=(0.0, -0.1, -1.0))
    except Exception as exc:                       # noqa: BLE001 — no display is a SKIP, not a fail
        pytest.skip(f"no GL context here: {type(exc).__name__}: {exc}")

    r.begin_frame(0.0)
    _draw_population(r, "membrane", cell.positions["membrane"], cell, 1.0, [], None)
    r.end_frame()
    _attach_hud(r, lambda: _caption(cell, "leak-test", None, None), on_screen=False)

    buf = wp.zeros((r.screen_height, r.screen_width, 3), dtype=wp.uint8)
    shots = []
    for _ in range(3):
        r.begin_frame(0.0)
        r.end_frame()
        r.get_pixels(buf, split_up_tiles=False, mode="rgb", use_uint8=True)
        shots.append(buf.numpy().copy())
    for k in (1, 2):
        assert np.array_equal(shots[0], shots[k]), (
            f"frame {k} differs from frame 0 at the SAME camera by "
            f"{np.abs(shots[0].astype(int) - shots[k].astype(int)).max()} — the caption left GL "
            f"state behind and everything after it is drawing wrong")
def test_UP_IS_UP_in_the_written_png(tmp_path: Path) -> None:
    """⚠ Every PNG this viewer wrote was upside down, for the life of the file.

    `_write_png` was fed `buf.numpy()[::-1]` under the comment *"GL origin is bottom-left; images are
    top-left"*. That is true of the default framebuffer and NOT of what `get_pixels` returns, which is
    already top-left. So headless renders, `make cell-png` and the `p` screenshot were all flipped —
    and nothing noticed, because a cell is nearly symmetric about the horizontal and a flipped dome
    reads as a plausible bowl.

    ⚠ **It was seen once and EXPLAINED AWAY.** When the caption came out inverted, the fix was to flip
    the caption's projection to match — compensating for the bug and hiding it for another hour. A
    correction that makes a symptom go away without naming a cause is a second defect on the first.

    ⚠ **This goes through `_render` and opens the FILE**, which the first version of it did not: that
    one built the buffer itself, so it pinned the convention and could not see the write path. Its
    positive control PASSED with the flip put back, which is the only reason the hole was found.
    Orientation is pinned by a MARKER, never by looking: a thick bar lying entirely at +Y.
    """
    import numpy as np
    import pytest
    from PIL import Image

    from aleph.viz.cell_app import CellFile, _render

    pos = np.array([[0.0, 3.0, 0.0], [0.0, 7.0, 0.0]])
    cell = CellFile(header={"thinning": "NONE"}, positions={"cortex": pos},
                    segments={"cortex": np.array([[0, 1]], dtype=np.int64)},
                    faces={}, frames={}, path=tmp_path / "marker.alephcell")
    out = tmp_path / "marker.png"
    try:
        _render(cell, ["cortex"], out, thicken=200.0, view="front")
    except Exception as exc:                       # noqa: BLE001 — no display is a SKIP, not a fail
        pytest.skip(f"no GL context here: {type(exc).__name__}: {exc}")

    img = np.asarray(Image.open(out).convert("L"))
    body = img[: int(img.shape[0] * 0.85)]         # the caption band is not the subject
    rows, _ = np.nonzero(body < 200)
    assert rows.size > 100, "the marker did not render; this test cannot say anything"
    centroid, middle = rows.mean(), body.shape[0] / 2
    assert centroid < middle, (
        f"a bar lying entirely at +Y has its ink centroid at row {centroid:.0f} of {body.shape[0]} — "
        f"below the middle, so the WRITTEN image is upside down")
