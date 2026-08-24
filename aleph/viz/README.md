# `cell_app` — the native cell application

One command, no arguments:

```
make cell
```

It finds the newest export, says which one it chose and how old it is, opens it at **full native
resolution**, and prints its keys.

## ⚠ Why this is not a web page — and the version of this section that stood until 2026-08-22 was WRONG

The old text argued from the **16 MB artifact cap**: 367 MB of sequence is 23× over it, so fitting
would mean one node in twenty-three, which `feedback-viewer-no-downsample` forbids. **That argument
does not apply to a local HTML file, which has no such cap**, and the PI was right to ask why the old
viewers worked. They did work, and one of them is still on disk and still opens:
`outputs/ac/world_phase1/cell_inside.html`, 78 MB, hand-written canvas, drag to turn.

**The real reason is one population.** That viewer's own header says it:

```
4,586,067 nodes · 388,413 positioned · 795,666 lines · every one drawn, none sampled
Drawn IN FULL: membrane, nuclear_envelope, microtubule, ... LEFT OUT ENTIRELY (not thinned): cortex
```

The cortex is **4,128,840 of 4,313,157 segments — 95.7%**. 795,666 lines cost 78 MB, so including it
is ~6.2× that (~480 MB in one file) and asks the browser to project **4,924,506** segments per frame
in JS. The old route is not beaten by the cell; it is beaten by the cortex, and it remains the right
tool for the other eleven populations.

⚠ The distinction the old header draws is the load-bearing one and it is why that viewer is honest:
the cortex was **omitted, not thinned**. Leaving a population out and saying so is a different act
from drawing one node in twenty-three.

## What it does

```
--list          what can be opened, and what each one is. Needs no path.
--info          the header AND what the file does NOT contain. Draws nothing, opens no window.
--only A,B      show these populations
--cut z-[:1.5]  keep one side of a plane. A SELECTION, never a thinning.
--view NAME     front | side | top | basal | oblique
--frame N       a frame of a sequence
--headless P    render one frame to a PNG and exit — for a machine with no display
--thicken F     scale line width; the widths are real radii, not a style
```

Keys, live, with nothing re-read from disk:

```
1-9,0   toggle a population        f / n   show all / none
k       cycle the cut (a ring of seven, because which half matters is not known in advance)
[ ]     slide the cut plane        r       play / pause a sequence      , .   step frames
o       what this cell HAS — and what it does NOT
h       help + current state       p       screenshot
mouse   drag to look, scroll to zoom       w/a/s/d or arrows   move the camera
```

⚠ **The camera belongs to the renderer, not to this app**, and four of this app's keys used to collide
with warp's own: `a` both showed every population and strafed left, `c` cut the cell and toggled the
axis, `i` printed the absent list and turned on warp's overlay, and **`space` was warp's PAUSE, which
spins inside `end_frame` and reads exactly like a hang.** They are now `f` / `k` / `o` / `r`, and the
module's self-check pins the whole app key set as disjoint from warp's so it cannot regress.

## The picture carries its own caption

⚠ **Until 2026-08-22 every word the app said about itself went to the terminal**, and the window
showed geometry and nothing else. That is a tool you drive from a terminal, not an app — and worse,
**the screenshot `p` writes carried no caveat at all.** The PNG is what leaves this machine and gets
pasted into a report; the scrollback is not.

The caption is now painted onto the render itself — the file, the view, the cut, `thinning: NONE`,
the cut-facing-away warning in red, and *"geometry only — a picture is not evidence"*. It is built
from the same variables the terminal block prints, so the two cannot drift on a fact; they differ
only in length, because the terminal carries the full argument and the picture carries the sentence
that has to travel with it.

It is painted into **both** the window and `_frame_texture`, in opposite orientations, because
`--headless` and the interactive `p` both read the texture and reverse it row-wise. Painting only the
one you happen to be testing is how `p` would have gone on writing caption-free PNGs from a window
that showed the caption the whole time.

⚠ **It failed silently three times while being written** — wrong framebuffer, depth test, and culled
back faces under a flipped projection — and every one printed `wrote … 3200x2000` and produced a
clean file. `test_the_caption_reaches_the_SAVED_PIXELS_and_not_only_the_window` opens the PNG and
counts ink in the caption band, because a test on the strings passed on all three broken renders.

## The three things it refuses to do

* **It never thins.** Every node of every population is present in every frame, and it says so in the
  file header and on every run. A cut removes elements by position and reports how many; that is a
  spatial selection and the status line insists on the distinction.
* **It never claims.** No force is evaluated and nothing moves. Every render prints *"geometry only —
  no force is evaluated, nothing has moved, and a picture is not evidence."* The γ panel plots the
  series **read from the record** and judges nothing: a run whose window opened on a transient draws
  exactly like one that did not.
* **It never silently omits.** A population with no topology is named as undrawn rather than skipped.
  An absent γ series prints why. A cut facing away from the camera says so and **does not
  auto-correct** — a viewer that quietly changes your cut is a viewer that can be quoted about a
  picture you did not ask for.

## The one it insists on doing

`i` (and `--info`) prints **what the file does not contain**:

```
HAS (11): chromatin, cortex, filopodium, intermediate_filament, lamellipodium,
          lamina, membrane, microtubule, nuclear_envelope, sf_arc, stress_fiber
⚠ NOT IN THIS FILE (2 of 13 builders in the tree): cytosol, nmii
```

⚠ **A census of eleven names does not tell a reader it is missing a twelfth** — you have to already
know the twelve. That is how the PHASE 4 τ runs were stepped and read for a full night before anyone
noticed the cell has **no myosin in it**. The records list what they contain, correctly, and nothing
listed what they lack. So the absent list is the computed one, and its reference is **derived from
the builder modules on disk, never typed**.

It is a difference, not an accusation: `--core-only` omits nine on purpose and a driver may have had
no reason to stand one. The line says so.

## ⚠ Why it used to stop the whole desktop

Not a frame-rate problem. **One idle redraw of the native cell costs 7.4 s** — 4.1 M capsule
instances — and the loop ran that continuously, with warp's default `vsync=False` and no sleep, for
as long as the window was open. The GPU was pinned at 100% whether or not anyone touched the app, and
the compositor starved with it.

A static picture does not need redrawing. The idle branch now pumps events, lets the renderer move the
camera, and re-renders **only when the camera actually moved**. Measured with the window open and
untouched: **0.7% CPU**. `vsync=True` is passed as well.

Two things that are NOT optional in that branch, both learned by reading warp rather than guessing:
the event pump, and `_process_inputs()` — warp calls the latter from `_draw`, so skipping the draw
would silently kill `w/a/s/d` while leaving the mouse working, and a half-dead camera is worse than a
slow one.

## The turntable — the other way past the size limit

`aleph/scripts/world_turntable.py` pre-renders a full turn and writes a viewer beside it. Turning the
cell in the browser is then an array index: **nothing renders while you look**, which is the point,
because a live window cost 7.4 s of GPU per idle frame.

Two things make the angle count cheap rather than the binding constraint:

* **webp, not png** — 3.0x smaller measured on a native frame. ⚠ Lossy: mean |delta| 1.3 of 255
  (0.5%), max 163 at single-pixel filament edges, caption legible. That trade is a viewer's to make
  and it is recorded in `turntable.json` and printed on the page. `--format png` when the frames are
  the artifact rather than the view.
* **a loading window, not a preload** — the first viewer created an Image for every frame up front,
  which made the angle count a MEMORY budget. Only the neighbourhood you can reach next is fetched,
  so angles cost disk and render time and nothing else.

What it cannot do, stated rather than discovered: arbitrary angles, a moving cut, or population
toggles — depth does not survive being flattened into a PNG, so layering populations would put a near
filament behind a far one. Those questions belong in the app.

## Where the figures go

`make cell-png` writes one PNG. The τ runs' own figures are `make tau3`, which reads the seeds in the
contract's declared order **first** and then draws — because contract III says nothing later may be
quoted if something earlier refuses, and the figures are drawn either way.

## Testing

`aleph/tests/world/test_cell_app.py` delegates to the module's own `_demo()` and pins the parts that
need no window: the cut arithmetic, the key state machine, the cut-facing sign, the sparkline, the
absence report and the listing. The module carries `_demo()` because the GPU host — where these files
are written — has no test runner, and **a self-check that only runs when invoked by hand is a
self-check that stops running.**
