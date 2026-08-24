# Figures — 2026-08-21, session `d380041d`

Both are renders of `aleph/outputs/ac/world_phase4/tau3_seed2.alephcell`, **frame 5 of 6**
(step 30,000), produced by `aleph/viz/cell_app.py` at full resolution.

⚠ **A picture is not evidence, and the app says so itself on every run.** These accompany
`THE_CENTROSOME_IS_INSIDE_THE_NUCLEUS_2026-08-21.md`, where the claim rests on the measured radii,
not on either image. The images are here because the first one is *how the question got asked*.

⚠ **Rendered on a CPU Warp build** (`Warp 1.14.0, CUDA not enabled, device "cpu"`) on the dev
machine, which has no CUDA. **No physics was computed** — `load()` reads stored positions and the
renderer draws them, nothing is stepped and no force is evaluated. A CPU number would not be a
result; a CPU *drawing of stored positions* is not a number.

## `cell_cut_zminus_seed2_f5.png` — 3200×2000

`--cut z- --view oblique`. Everything, cut at z = 0, keeping z ≤ 0.

⚠ **This picture is the reason the next one exists, and it is nearly uninformative — but the first
explanation written here for WHY was wrong, and is corrected in place.**

It said the cut and the camera are unrelated, so the retained hemisphere's convex surface faced the
viewer. **It did not.** `oblique` sits at z = +22 and `--cut z-` keeps z ≤ 0, so the camera is
already on the side the cut removed and the cut face IS toward it. What hides the interior here is
**the cortex's own inner surface**: 4.2 M filaments on a shell are opaque from inside as well as
outside, which `cell_app.py`'s header already said when `--cut` was added.

⚠ **`--cut z+` is the case that really is wrong**, and it is computable: the camera sits 22 µm inside
the half that cut keeps, so the closed dome is in front and the picture looks like an uncut cell.
`camera_is_on_the_kept_side()` now says so at render time, in the app, with the cut to use instead —
**and does not auto-correct**, because a viewer that quietly changes your cut is a viewer that can be
quoted about a picture you did not ask for.

## `centrosome_inside_nucleus_seed2_f5.png` — 3200×2000

Same cut and camera, `--only microtubule,nuclear_envelope,stress_fiber,sf_arc,lamellipodium,`
`intermediate_filament,chromatin` — cortex and membrane dropped so the interior is visible.

The microtubule aster (blue) converges to a single point. The nuclear envelope (purple) reads as a
**crescent** rather than the hemispherical bowl a centred sphere cut at z = 0 should give. That is
what prompted the measurement:

```
microtubule        74,001 nodes   r 0.000 .. 7.400 um   about the origin
nuclear_envelope   40,962 nodes   r 5.100 .. 5.100 um   exact sphere, on the origin

50,735 / 74,001  =  68.6%  of MT nodes lie inside the nuclear envelope
```

Units: µm throughout. Nothing is downsampled — the app prints
*"thinning: NONE. Every node of every population is present in every frame."*
and the per-population element counts it drew are in the run logs.
