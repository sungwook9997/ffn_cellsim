# Overnight autonomous run — summary (2026-06-11, PI asleep ~4 h)

Branch `h7/compartment-platform`. All results committed (NOT pushed). gbook `~/ffn_cellsim`
(dirty) untouched; all work in the freshly-synced `~/ffn_cellsim-platform`. Runs in detached
`screen` (survive wifi drops); 10-min night-loop `ad9c0c16` as backup.

## What landed (committed)
1. **GPU UNBLOCKED** (`3a7e527`) — K1-turgor kernel replaces native md.mesh per-cell-types
   (which stalled the A5000 at 0% util); N=60 spheroid runs finite in 97 s. + activity-LOD
   (freeze inert core).
2. **a+b/R+c/R² law FORM reproduced on the A5000** (`55f23c0`) — active size sweep R=31-78 µm:
   fit a=−0.27, b=154, c=−3152 (signs + order MATCH the PI law −0.33/188.7/−2655), A/A₀
   DECREASES with R (corr −0.67). The robust scientific result.
3. **REAL-TIME** (`f6916df`) — `common/sim_realtime.py`: every production run + viz now shows
   t_sim AND a real-time-equivalent via accel S=6e5 (lamellipodium-velocity basis). Honest
   ceiling: ~hours wall = ~seconds real (BAOAB-sync step rate), so accelerated sim covers the
   EARLY spreading episode.
4. **Large necrosis-ON 3-zone** (`...`) — N=400 coarse-patch (R_cell=22 µm) → R>150 µm →
   necrosis ON: rim 225 / quiescent 175 / **necrotic core 24**, A/A₀=3.15.
5. **Magnitude gap RESOLVED as a TIMESCALE truth** (`113aadc`) — fixed the prolif build's
   active classification (now = plain build) + ANCHORED the division rate to the cell cycle
   (`--auto-pdiv`, p_div=S·dt·div_every/T_cycle). A single 7.2 s-real episode = 1e-4 of a 20 h
   cycle → ~0 divisions (correct physics, not a bug). The PI law's b/R is a DAYS-LONG-assay
   proliferation integral the accelerated single-episode sim cannot reach in feasible
   wall-time; the FORM is reproduced per-episode; the magnitude must NOT be faked by tuning.
6. **Size-dependent necrosis** (`...`) — sweep N=200-800 (R=180-275 µm): **necrotic core grows
   3→24→64→101** (1.5→11.2 %), proliferating rim fraction falls 0.64→0.41, quiescent rises
   0.34→0.48 — Greenspan surface/volume; the falling rim fraction IS the b/R mechanism.
   `figs/necrosis_size_scaling.png`.

## Known issues / honest limits
- **A/A₀ unstable at large N** (N=600 gave 14.5, a basal convex-hull footprint artifact when
  peripheral cells disperse); the zone/necrotic COUNTS are clean + monotonic.
- **N=1500 run errored** (`IndexError: index 1372 out of bounds size 1372`) — a large-N
  indexing bug (some cells dropped/merged vs the assumed count); needs a fix before N>~1000.
  Realistic spheroid scale (R 180-275 µm) is already covered by N=200-800, so not critical.
- **Fine-grained R>150 µm** (thousands of 7.5 µm cells) is hours-to-days of wall-time — needs
  the BAOAB GPU port (frozen, PI-gated). Coarse patches are the feasible route used here.
- **Notion MCP disconnected** overnight — git is the authoritative record; Notion to be updated
  on reconnect.

## Bottom line for the PI
The spheroid platform now runs GPU-resident at scale (hundreds of cells), reproduces the
spreading-law FORM + size-dependence, shows the full size-dependent 3-zone necrosis, carries
real-time labels, and the proliferation magnitude question is resolved honestly (a multi-day
timescale integral, not a tunable knob). Next decisions for you: (a) fix the large-N index bug
+ push bigger; (b) the BAOAB GPU port (frozen) for fine-grained R>150 µm + longer real-time;
(c) the remaining refinements (junction-switch dynamics, explicit ECM, fully integrated
spheroid).
