# Handoff — next Lead session, 2026-08-21 14:05 KST

Paste the block below as your opening prompt. Everything in it was measured, not recalled; where a
claim was later retracted, the retraction is included rather than the claim being deleted.

---

## Where you are starting

Branch `engine/main`, HEAD `0df18bc5`. **114 commits** since `c41d4e95` (05:15 today), all with a
pathspec. `aleph/tests/world` green, ruff clean, `make kb-check` **no DRIFT**, boot 47,943 / 48,000.

**GPU grant in force**: job 73, 4090-1, **2:54 remaining** (ends 17:00:46). The τ runs are finished —
**the card is free for the next thing, inside that window.** ⚠ Beyond 17:00 you must ask the PI. You
may transcribe a grant with transcript path and timestamp; **you may never write one.**

## What landed

**PHASE 4 τ run III is complete.** Three seeds × 30,000 samples, all rc 0, all carrying the **same
`closure_digest`** — the same code ran three times.

```
status                     REFUSED
min_tau_reliability        10.95     (bar 50)   <- refuses at the FIRST field
windows_open_on_transient  [1, 2, 3]            <- all three seeds
equilibration_samples      14,625    window_samples 15,375
```

**The run is good and the analysis cannot read it.** The module's equilibration search cuts all three
seeds at the same 14,625 because `stationarity.py:607` is `limit = max(1, x.size // 2)`, so every τ is
measured across a rise. **That is PI queue item 18, now confirmed at three seeds.**

⚠ **Do not patch `stationarity.py:607` if you are holding the result it would rescue.** That was this
session's recorded reason for leaving it: *a shared threshold must not be moved by the person holding
the result.* If the PI rules, patch it.

**What the three runs establish**: settling behaviour differs between seeds at fixed configuration.
That is an input to **D-3 (initial-condition independence)**, not D-2. **One seed of three produced a
τ. There is no τ scatter.** No magnitude is quotable — `(c) 3` and `(c) 17` stand.

⚠ **Two claims this session made about seed 3 are RETRACTED** (`0df18bc5`): *"two of three seeds are
measurements"* and *"the first seed-to-seed τ comparison, ratio 1.30"*. τ at cuts 21k–26k is
`115.8/97.7…` for seed 1 (spread 16.7%, converged) and `89.2 → 49.0` for seed 3 (spread 61.2%,
monotone falling); the ratio moves **1.11 → 2.28** with the cut, so it is a property of the cut. And
**`n_eff = N/(2τ)`** — a falling τ mechanically raises n_eff, so seed 3 passing the length gate is not
independent evidence. **Read that document before quoting anything about seed 3.**

## The queue is 23 items and that is the work

`aleph/docs/v2_audit/PI_QUEUE_2026-08-21_MORNING.md`, index at the head.

**BLOCKING, and 18 is upstream of most of the rest:**

| | |
|---|---|
| **18** | τ data is good, the analysis cannot read it. **Everything chemical is downstream** — a kinetic connector may only commit on an accepted step and **zero steps have ever been accepted.** |
| **14 → 11** | Motors sit **0.250 µm** from actin against a 0.210 µm capture radius. **No value of the eight crossbridge constants produces a single bond.** Read 14 before 11. |
| **3 ← 20** | Archetype: sphere-in-contact or spread cell. Item 20 hangs off it — **the MTOC is at (0,0,0), the centre of the nucleus, with 68.6% of MT nodes inside the envelope.** |
| **19** | Myosin count resolves against the membrane sphere (706.86 = 4π·7.5²) while `plan_nmii`'s own docstring asks for the cortical shell. **442 vs 379, +16.5%.** |
| **21** | `STATE.md` §(c) says it *only grows* and the header caps the file at 22,000 B. **39 B free, next row needs 141.** `c19` is the first entry with no row. |
| **22** | **The cell is attached to nothing** — `world/` binds 11 of 46 laws and **not one is adhesion or ECM.** |
| **23** | **0 of 77 corpus dossiers are citable.** 72 unregistered, 5 registered all `CHECK`, none OK. 7 preprints (p13–p19, consecutive). |

## The adversarial audit answers the representation question — read it before adding any protein

`aleph/docs/ADVERSARIAL_AUDIT_HOW_TO_REPRESENT_A_PROTEIN_2026-08-21.md`

**`arena.py:76-78` already decided it**, in the argument that closed the primitive set: a molecular
attachment is a **`BOND`** with its far address on a `SEGMENT`; a binding site is **not an allocation**;
a species pool is a **`GRID_CELL` channel**. **Both of those primitives are at zero** — bond 0/1,000,000
and grid_cell 0/4,000,000.

**A node is 50 nm** (measured), ~9 actin monomers. Integrin (10–20 nm) and α-actinin (35 nm) are
**sub-node**, so giving them nodes resolves a receptor more finely than the actin it binds. **Cost is
not the obstacle**: a whole-cell basement membrane is 565,520 nodes against 7,441,446 free.

⚠ **The sharpest gap found, and it is live rather than absent**: `G_actin_uM` is a **declared scalar**
passed to a law resolver while `grid_cell` live is 0 — a spatially uniform constant standing in for a
diffusing pool. **That is the one place a "probability cloud" is not an option but the correct
answer**, and it is a lumping the charter forbids, today.

## Loose ends someone must decide

* ⚠ **`aleph/outputs/tag_kb/snapshots/notion_snapshot.json` is modified and uncommitted.** Written by
  `refresh.sh` stage 1b when session 75 built `kb.duckdb`. It is the **`outer` lane's** — neither
  session 75 nor this one committed it. Leaving it makes the next session read it as their own change.
* `STATE.md` carries no τ run III row and **cannot** until item 21 is ruled. Session 21 owns the
  stationarity reporting; coordinate rather than duplicating.
* The PI ratified *"use a physiological range where no source exists"*. **It was not acted on.** The
  highest-value target is **`mobility`**: it is what turns step count into physical time, and until it
  is sourced the file's own header is right that `dt` is *"a step index wearing a unit"* — **there is
  no real-time ratio to quote.**

## The app

`aleph/viz/cell_app.py` + `aleph/viz/README.md`. `make cell` opens the newest export; `--list` says what
can be opened without a path; `i` / `--info` prints **what the file does NOT contain**; a γ sparkline
puts the current frame on the run's curve; a warning fires when the cut faces away from the camera.
Opens **367 MB / 6 frames / 4.56 M nodes with nothing downsampled** — 23× an artifact's 16 MB cap,
which is the measured reason it is not a web page. Tests in `aleph/tests/world/test_cell_app.py`.

## What worked, and it is worth keeping

* **Reproduce before claiming, not after.** Every claim caught before sending was calculated first;
  every retraction was written first and calculated after.
* **Hand over more columns than your conclusion needs.** The 1.30 refutation was only possible because
  three cuts were sent without interpretation. One cut would have entered the record as a seed property.
* **Execute the thing.** The app's status bar printed nothing for as long as it ran; a figure's five
  tick labels overlapped into one number; a positive control passed until the matcher was deliberately
  broken. None of those is findable by reading.
* ⚠ **And the standing one, which this session then violated**: a check must not be entangled with its
  own subject. `n_eff` is a function of `τ`, and it was used to certify `τ`.

---

**Read `STATE.md` and `CLAUDE.md` first, then the queue index. Do not quote a magnitude.**
