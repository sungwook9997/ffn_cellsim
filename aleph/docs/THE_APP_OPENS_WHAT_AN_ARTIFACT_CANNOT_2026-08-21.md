# The app opens what an artifact cannot — measured, not asserted

**2026-08-21, session `d380041d`.** The PI's instruction was an app rather than a web page, and gave
the reason: *"웹 html 방식 더 싫은 이유는 용량 제한 때문에 나는 풀로 보고 싶은 것인데"* — the size cap,
because the whole cell is what is wanted. This records the measurement that the choice was right,
rather than repeating the reason.

## What was opened

`aleph/viz/cell_app.py`, `load()` path, on the two real files that exist — not a fixture.

```
PHASE 1 single cell            94.5 MB    0.02 s    1 frame    4,591,262 nodes   12 populations
PHASE 4 tau seed 2 sequence   367.4 MB    0.21 s    6 frames   4,558,554 nodes   11 populations
                                                    peak RSS 0.75 GB
```

Every node is present. Nothing is decimated on the way in:

```
cortex                 4,197,654
membrane                 163,842
microtubule               74,001
filopodium                61,000
nuclear_envelope          40,962
nmii                      32,708    (PHASE 1 only -- see below)
stress_fiber               9,240
sf_arc                     2,600
intermediate_filament      2,460
lamina                     2,043
chromatin                    552
```

## Against the cap

An Artifact renders at **16 MB or smaller**, embedded assets included.

```
367.4 MB  /  16 MB  =  23.0x over
```

To fit, 4,558,554 nodes would have to become about **198,000** — one node in twenty-three. The cortex
alone is 4.2 M of them. **That is the downsampling the PI has forbidden in viewers since the
`feedback-viewer-no-downsample` ruling**, so the web route does not have a version that satisfies the
instruction; it has a version that violates a different one.

⚠ **This is not an argument that HTML is bad.** It is a measurement that this particular file does not
fit, by a factor with two digits in it. A 4 MB cell would fit fine, and none of this would apply.

## What is verified and what is not

* **Verified**: the load path, at full scale, on both real files, including the six-frame sequence —
  topology written once and positions per frame, which is why 6 frames of 4.56 M nodes is 367 MB and
  not 6 × 94 MB.
* **Verified**: the cell moves between frames, so the sequence is a sequence and not six copies:

  ```
  frame 1   max |dx| vs frame 0 = 0.007102 um    mean 0.001624
  frame 2                         0.010596              0.002068
  frame 3                         0.012830              0.002386
  frame 4                         0.013467              0.002640
  frame 5                         0.014562              0.002859
  ```

  ⚠ **Descriptive only.** Six points at 6,000-step spacing is not a relaxation measurement, and this
  is seed 2 — the replicate that does not reach stationarity inside its own run. The displacement is
  still growing at the last frame, which is consistent with that and is not evidence for it.
* **NOT verified**: the render path at 367 MB with a display. It needs the file and a display in the
  same place; the file is produced on the GPU host, which has no display, and the dev machine has no
  CUDA. The render path was exercised earlier tonight on the single-frame file. **The multi-frame
  render at full size is untested and is recorded here as untested.**

## The thing the census does not say out loud

The PHASE 4 sequence carries **11** populations. The PHASE 1 cell carries **12**. The one that is
missing is `nmii` — 32,708 nodes, the cell's only active element. A reader who opens the sequence
sees eleven population names and no gap where a twelfth would be.

See `THE_MYOSIN_COUNT_RESOLVES_AGAINST_THE_WRONG_SPHERE_2026-08-21.md` and PI queue item 19.
