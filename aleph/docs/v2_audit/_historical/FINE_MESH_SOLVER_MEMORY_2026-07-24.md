---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# Fine physiological cortex mesh — solver + memory analysis (2026-07-24)

**Read-only analysis. No sim code changed.** Follows `CORTEX_MESH_FIDELITY_2026-07-24.md` §9, which found the
FINE cortex (`cortex_seg_um=0.075`, `cortex_density_per_fil=40`, **2,898,126 actin nodes / 3,569,158 total**)
does not converge: pathA (matrix-free, weak) plateaus at `max|PF|=3.46`, pathB (explicit BSR `A_c`, strong)
**OOMs** (an 85 MB alloc fails ⇒ device already full). This memo estimates the peak GPU bytes from the code,
ranks the cheapest ways to fit a strong solve in the A5000 16 GB, and gives the honest single-vs-multi-GPU verdict.

## Headline verdict

**Single-A5000-FEASIBLE — the 75 nm / 2.9 M-node mesh does NOT genuinely need multi-GPU to FIT.** The OOM is
**not a fundamental capacity wall**; it is two implementation inefficiencies in pathB, both avoidable:

1. **Materialising 11.4 M explicit `MAT66` triplets and sorting them** (`vals_d` 3.28 GB + the
   `bsr_set_from_triplets` sort transient ~2–4 GB ≈ **5–7 GB**), and
2. **A dense per-fiber Cholesky that stores ~94 % zeros** (`fiber_block_factor` **2.85 GB**, scaling as
   `max_fiber_nodes²` — a 34× blow-up from the coarse mesh).

The *intrinsic* fine-mesh footprint — state (~0.65 GB) + shared solver work vectors (~1.8 GB) + an
**efficiently-assembled** explicit `A_c` (~1.4 GB) — is **~4–5 GB, a third of the card.** Multi-GPU / a bigger
card is only genuinely required at the **30–50 nm NATIVE stretch** (seg 0.03–0.05 ⇒ 4.3–8.5 M nodes), where
state alone approaches ~2 GB, an efficient `A_c` grows past ~5 GB, and conditioning worsens further.

**Separate, still-OPEN question — convergence.** The memory fixes let pathB *run* at the fine mesh; whether the
strong deep-block-Jacobi CG actually drives `max|PF|` from 3.46 down to the ~0.2 gate at 2.9 M nodes is
**unproven natively** (the CPU prototype `_pathb_selfcheck` passes on a *downscaled* curved sub-isostatic
operator — supporting, not native, evidence). So: **memory = single-A5000-solvable; native fine convergence =
plausible-but-must-be-measured once it fits.**

## 1. Memory breakdown at 2.9 M nodes (fine, membrane subdiv 8)

Sizes: `n_total = 3,569,158`, `n_actin = 2,898,126`, `n_xl = 2,827,440`, `n_fibers = 70,686`,
`max_fiber_nodes = 41` (fine; 7 coarse), `n_tri = 39×70,686 = 2,756,754`, `n_seg = 40×70,686 = 2,827,440`.
`vec3d`(f64) = 24 B; `MAT66`(f64) = 288 B; `VEC6`(f64) = 48 B. MB = ÷10⁶.

### (a) Fine cortex STATE (persistent) — NOT the bottleneck

| array | dtype | count | MB |
|---|---|---|---|
| `pos_d`, `f_d` (`assemble.py:834`) | vec3d | 2×n_total | 171 |
| `node_volume_d` f64 / `solid_active_d` i32 | — | n_total | 43 |
| `node_fiber_d` i32 | i32 | n_actin | 12 |
| crosslinks `xl_d`(2i32)+`kxl_d`+`r0xl_d`(2×f64) (`:841`) | — | n_xl | 68 |
| bending `tri_d`(3i32)+`alpha_d`(f64) | — | n_tri | 55 |
| `srest_d`(f64,n_seg) + foff/toff/soff/barbed | — | — | 24 |
| steric hash grid (256³ cells + qpts/fiber_id/active) | — | — | ~140 |
| membrane subdiv8 (655,362 v, ~1.31 M f) + ERM | — | — | ~80 |
| nucleus (~15.7 k v) + Biot fluid grid (43³) | — | — | ~10 |
| **STATE subtotal** | | | **~600–650 MB** |

State grows ~linearly with nodes/xl/tri but is small; it is **not** the OOM. (Matches the audit's subdiv-6
≈ 592 MB anchor.)

### (b) pathA solver (`ProjectedAnalyticCG` + probe) — fits ~5.8 GB

| buffer | where | size | MB |
|---|---|---|---|
| 13 × n_total vec3d work vectors (`r,z,projected_z,p,ap,dx,k_input,k_output,projected_k,preconditioner,coarse_external_diagonal,coarse_fine,coarse_fine_action`) | `implicit_mechanics.py:1504–1528` | 13×n_total | 1,114 |
| `fiber_block_work` (3·n_total f64) | `:1520` | n_total-equiv | 86 |
| `FiberQuotientCoarse.fine`+`fine_action` (pathA only) | `fiber_quotient_coarse.py:442–443` | 2×n_total | 171 |
| FQ block_factor (n_fibers·max_rank² f64) + coarse vecs | `:439` | — | ~20 |
| **probe scratch** `f,projected,trial_pf,pos_prev,pos_trial,pos_best`(+`pos0`) | `ac_gate_a_fq_coarse_test.py:93–98,194` | 7×n_total | 600 |
| **`fiber_block_factor`** (n_fibers·3·max_fiber_nodes² f64) | `:1515–1519` | 41² dense | **2,852** |
| **pathA solver subtotal** | | | **~4,840** |
| **pathA TOTAL (+state)** | | | **~5.5 GB → FITS** (plateaus, no OOM) ✓ |

`fiber_block_factor` alone is 2.85 GB and is the single largest solver array — at coarse (7 nodes) it is only
83 MB; the 41²/7² = **34×** blow-up is pure `max_fiber_nodes²` dense storage.

### (c) pathB extra (`FiberQuotientCoarsePathB` replaces the pathA inner) — OOMs

pathB still allocates the shared block-b buffers **minus** the two pathA FQ fines, i.e. 13 work vecs
(1,114) + `fiber_block_work` (86) + probe (600) + **`fiber_block_factor` 2,852** ≈ 4,650 MB, **plus**:

| pathB-specific | where | size | MB |
|---|---|---|---|
| **`vals_d`** `MAT66` triplets (`4·n_xl+n_fibers = 11,380,446`) | `fiber_quotient_coarse.py:814` | 288 B ea | **3,277** |
| `rows_d`+`cols_d` i32 triplets | `:812–813` | 2×11.4 M | 91 |
| **`A_c` BSR values** (deduped ~4–5 M blocks) + col idx | `:815, :843` | 288 B ea | ~1,400 |
| `diag_blocks`+`bj_factor` (`MAT66`×n_fibers) + vec6 workspaces | `:816–824` | — | 60 |
| **`bsr_set_from_triplets` sort/scan TRANSIENT** (11.4 M × 288 B) | `:843` | — | ~2,000–4,000 |
| **pathB extra subtotal** | | | **~6,800–8,800** |
| **pathB PEAK (+state+shared)** | | | **~12–14 GB steady + GB-scale transient spike → 16 GB → OOM** ✓ |

**Which dominates the OOM?** Not the state (~0.65 GB), not the work vectors (~1.1 GB). It is the **explicit-`A_c`
machinery** (`vals_d` 3.28 GB + BSR ~1.4 GB + sort transient 2–4 GB ≈ **7–9 GB**) stacked on the **dense
`fiber_block_factor` 2.85 GB**. pathA fit at ~5.5 GB; pathB adds ~7–9 GB on top and the triplet-sort transient
tips it over — exactly why an 85 MB alloc (`fiber_quotient_coarse.py:843`) is what fails.

## 2. Cheapest wins to fit pathB (or a strong matrix-free solve) in 16 GB

| # | win | code? | frees | note |
|---|---|---|---|---|
| **1** | **`disable_fiber_block=True` for the fine run** — drops `fiber_block_factor` entirely | 1-line solver kwarg (`implicit_mechanics.py:1487–1489,1515–1519`; `driver.py:462`) | **2.85 GB** | The 2026-07-23 diagnosis already found the per-fiber block OVER-STEPS; the pathB coarse supplies the inter-fiber correction, so a pure node-Jacobi smoother is the intended pairing. `use_fiber_block=False ⇒ factor_size=0`. |
| **2** | **membrane subdiv 8→6** (655 k→41 k verts) | CLI `--membrane-subdiv 6` (`ac_gate_a_fq_coarse_test.py:163`) | ~0.45 GB | Audit already blesses this for isolating the CORTEX effect. n_total −17 % shrinks every n_total vec3d. Small. |
| **3** | **alias the probe's 7 n_total vec3d scratch → ~4** | trivial script edit (`:93–98`) | ~0.25 GB | `f/projected/trial_pf` and `pos_trial/pos_best` can reuse. |
| **4** | **assemble `A_c` on a PRE-BUILT static fiber-adjacency BSR pattern via atomic scatter-add** — eliminate `rows_d/cols_d/vals_d` + the `bsr_set_from_triplets` sort | MODERATE (build the fiber-adjacency CSR once at setup from the static crosslink→(f,g) map; per build zero `A_c.values` and atomic-add the 4 projected 6×6 blocks/crosslink into their CSR slots) | **~5–7 GB** | **THE decisive structural fix.** Net pathB coarse memory → `A_c` values (~1.4 GB) + `bj_factor` (40 MB) only. Topology is static (crosslinks don't rewire in the resting solve), so the pattern is built once. |
| **5** | **banded per-fiber Cholesky** — store `max_fiber_nodes × bandwidth` not `max_fiber_nodes²` | MODERATE (`build_fiber_block_cholesky_kernel`) | ~2.7 GB | Bending couples i,i±1,i±2 ⇒ half-bandwidth 2, so a dense 41² block is ~94 % zeros → ~160 MB banded. Alternative to #1 if the fiber block turns out needed for convergence. |
| **6** | **single-precision (f32) solver scratch** | LARGER | ~1–3 GB | Halves work vecs / `vals` / `A_c`. **Deferred** — risks the SPD Cholesky conditioning and the `sqrt(eps64)` tolerance contract; not first-line. |

**Combination that fits 2.9 M fine + a strong solve in 16 GB:**

- **No-code (borderline):** #1 + #2 + #3 free ~3.5 GB → pathB steady ~9–11 GB, but the `bsr_set_from_triplets`
  transient (2–4 GB) still spikes on top. We were only ~85 MB short, so this **may** clear the OOM — but the
  transient makes it risky.
- **Robust (the real fix):** **#4 alone** frees ~5–7 GB → pathB peak **~6–8 GB even at subdiv 8** → fits with
  wide headroom. Add #1/#5 and it drops another ~2.7 GB.

On the "strong matrix-free without explicit `A_c`" option (prompt 2c): pathA is *already* matrix-free but weak
**for a time/iteration reason, not memory** — each inner iter is a full fine-operator apply
(`fiber_quotient_coarse.py:461–472`), so only ~40–200 inner iters are affordable and per-fiber block-Jacobi
inner-CG cannot converge the ill-conditioned 424 k inter-fiber `A_c` in that budget (`:522–533`). Raising
`fq_coarse_iterations` on pathA is cheap in memory but each iter is expensive in *time* and the preconditioner is
too weak — so more pathA iters alone will not close 3.46 → 0.2 in a practical budget. The explicit `A_c` (cheap
SpMV ⇒ 600–1500 iters) is what makes the strong solve tractable; win #4 keeps that strong SpMV while removing its
memory cost. That is strictly better than either raw pathA or triplet-pathB.

## 3. Verdict + finest feasible rung

- **75 nm / 2.9 M nodes / density-40 (this target): single-A5000-FEASIBLE-with-win-#4** (pattern-atomic `A_c`),
  or borderline-feasible no-code via #1+#2+#3. Convergence to be MEASURED once it fits.
- **100 nm rung (seg 0.1 ⇒ ~2.1 M nodes):** an even safer fallback if 75 nm's conditioning proves too stiff for
  the strong solver's budget — same memory story, smaller.
- **30–50 nm NATIVE (seg 0.03–0.05 ⇒ 4.3–8.5 M nodes):** genuinely pressures 16 GB (state ~2 GB, efficient
  `A_c` >5 GB, conditioning worse) — **this is the multi-GPU / bigger-card case** (per the HARD rule: never lower
  density to fit memory), the stretch/validation resolution, not the first fine build.

**Resolution ladder:** 500 nm (validated) → **75 nm 2.9 M (A5000-feasible with #4; convergence TBD)** → 30–50 nm
native (multi-GPU / bigger card).

## Recommended single next native experiment

Cheapest run that yields a fine-mesh **convergence number** while directly testing the memory fix — free
~3.3–3.5 GB (the OOM shortfall) with a one-line solver-config addition (not sim physics):

```
# ac_gate_a_fq_coarse_test.py: pass disable_fiber_block=True into ProjectedAnalyticCG(...) (~line 210)
python aleph/scripts/ac_gate_a_fq_coarse_test.py --native --outer 40 \
    --cortex-seg-um 0.075 --cortex-density 40 \
    --membrane-subdiv 6 --mode pathB --fq-iters 1500 --cg-iters 60
```

Log peak bytes from the assemble ledger (`gpu_bytes_used`, `warp_mempool_used_high_bytes`,
`assemble.py:1011–1016`) and the fine `max|PF|` trajectory.

- **If it fits AND converges** → the 75 nm rung is closed on one A5000; proceed to GATE-B dynamic slice at fine.
- **If it fits but plateaus** → conditioning, not memory, is the wall → raise `--fq-iters` / add a stronger inner
  preconditioner (the convergence question, separate from this memo).
- **If it still OOMs on the `bsr_set_from_triplets` transient** → win #4 (pre-built-pattern atomic `A_c`
  assembly) is REQUIRED, after which subdiv 8 fits with headroom.
