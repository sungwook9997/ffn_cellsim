# STATE_NONQUOTABLE.md — the full reasons a result may NOT be cited

Split out of `STATE.md` §(c) on **2026-07-28 (PI decision)**. The reason is structural, not cosmetic:
this list only ever GROWS. A retraction is permanent — nothing is ever removed once a number has been
shown not to transfer — so a fixed line budget on `STATE.md` was guaranteed to collide with it every
few sessions, and the collision was resolving the wrong way: the budget was editing the content, and
in the session that prompted the split it had already forced an already-CLOSED entry to be compressed
to make room for a fresh retraction.

`STATE.md` §(c) keeps one line per entry — the claim and the one-line reason — and links here for the
evidence. **The short form is the contract; this file is the working.** If the two ever disagree,
`STATE.md` is the one a reader is assumed to have read, so fix it first.

Anchors are `#cN`. `check_state_md.py` asserts every `STATE.md` §(c) pointer resolves to one that
exists here, because a link that does not open is worse than the inline text it replaced.


<a id="c1"></a>

**1.** **GATE A "CLOSED" / `max|PF|` 0.1398** — status is **PROVISIONAL, PI ratification pending** (Card G-A): the observable (raw `|F|` → projected `|PF|`) *and* the configuration (resting myosin removed from the static baseline) both changed between the FAIL and the PASS; no data artifact exists (one PNG, `df8404cb`) and the closure config is not reconstructible from disk. The `0.21 pN` threshold is also **not a physics threshold**: `GATE_A_THRESHOLD_DERIVATION_2026-07-25.md` identifies it as `tolerance_um/dt_mu = 10·√eps64·ℓ·k_max = 0.2107 pN`, the runtime's own float64 round-off *displacement* tolerance evaluated at one mesh and then frozen by hand into 54 sites — so it scales ∝ ℓ, and scoring the 200 nm / 75 nm rungs against the frozen literal loosens the gate 2.5× / 6.6×. Neither that derivation nor any replacement metric is PI-ratified.
  **REFRAMED 2026-07-28 (night), and this is the substantive half:** the configuration change was not a convenience — `GATE_A_RESTING_CONVERGENCE_2026-07-23.md` §10 (2026-07-24) reports the decisive control. The identical cell with resting myosin **OFF** converges seed-robustly (`max|PF|` 0.1398 / 0.1416 / 0.1586 across three cortex-construction seeds, ≥24% under the floor), while seeding 4,420 discrete isometric heads injects a ~1.5 residual that **no solver removes**: `PF ≈ f_head`, ∝ f_head, concentrated on the ~1% myosin sites, **inherent to the discrete point loads rather than a conditioning artefact** — closing it would need ~7× more heads, which is unphysical. So the static solve is not failing; it is being asked for an object that does not exist. A physiological resting cortex's tension is maintained **dynamically** by myosin turnover — a dynamic steady state, not a static equilibrium — and forcing it into a fixed isometric seed is the **category error** the residual reports. Consequences already visible elsewhere in this file: the 6.30×-above-floor native stall in (b) is that same irreducible active force and is **NOT a preconditioner target**; and the whole gate vocabulary is static-shaped — no stationarity criterion exists anywhere in the repo, and GATE B's tension "emergence" was observed over 0.3 s of physical time against a 2.5 s bound-head lifetime (`k_off0` = 0.4/s, Nagy 2013), i.e. 12% of ONE turnover, so it measured a transient. Lead recorded this after quoting the superseded 2026-07-22e conditioning diagnosis for several turns.
<a id="c2"></a>

**2.** **γ-floor "~530× / 100–2300× under band"**, and every full-native `ff/` headline with it (`engine_reval` P6.1 γ = 0.150 mN/m, `mech_hier` S1/S2/S3, S1 Hertz `E_fit` 3624 Pa, Biot FSI τ_p 1.15×, SC0) — measured on the pre-2026-07-23 cortex build, **no artifact records which build**, not re-run, and three of that audit's five defects are still live in the γ-floor production path. Banners: `ff/ENGINE.md`, `outputs/mech_hier/REPORT.md`, `outputs/engine_reval/RETIRED_BUILD_BANNER.md`.
<a id="c3"></a>

**3.** **γ = 3.72 pN/µm as "real power-stroke tension"**, and with it **every SF/cortex motor tension this project has reported** — the crossbridge term is 0.0191 (0.51%); the rest is turgor-loaded passive network; no t0 row is measured, every step is force-accepted, `max|PF|` plateaus ≈0.5 against the 0.21 gate, and every NMII magnitude is a PI-GAP. **Sharpened, then SCOPED the same night (2026-07-28)**: swept against a converging inner solve, the SF lane's reported tension tracks the residual 1:1 over seven decades at a fixed bound-head count — so at that sampling instant the magnitude was the un-relaxed residual, not a transient on its way anywhere. **The scope is the point and it was missed once**: every run in that sweep was 20 steps = **0.20 s**, which is inside an induction period; the converged tension there is 3.7e-06 pN and by 0.50 s it is 0.947. Past the induction period the tension is real and the gate PASSES (see (b)). So this entry retires the MAGNITUDES, not the mechanism — and the reason they are non-quotable is that they were read at a time and a convergence at which nothing had happened yet.
<a id="c4"></a>

**4.** **Any `CONNECTED` rung, and the "first composed native cell"** — but **NOT for the reason this entry recorded until 2026-07-28 (night)**, and the old reason was wrong in the direction that matters. It said "2 of 32 connectors are ever dispatched" and "`ac_composed_world_dump.py` currently raises". Checked against the code: **8 canonical facades claim 36 of 36 connectors** and `validate_canonical_facade_claims(reference_cell_architecture())` **PASSES**, so no connector lacks an owner; and the dump does **not** raise — it writes `composed_world_v2.npz` with 14 components / 36 connectors / 71,624 filaments. The wiring is not the gap. **The gap is that nothing solves.** `composed_native.py::solve_candidate` is one force-assembly pass — zero the force arrays, run the exact-once pipeline, compute head loads — with no iteration, no Newton, no relaxation; and `CandidatePhase.SOLVE_COUPLED` has **zero production tasks** assigned anywhere in the repo (one test file only). So a composed candidate transfers forces correctly ONCE and stops at whatever residual that leaves, which is why no magnitude from it is a physical one — not because the pieces are unconnected but because the assembled system is never brought to force balance. This splits the two engines cleanly: `ac/cell` HAS an iterative solve and stalls at 6.30x its own derived floor (see (b)); `ac/engine` does not stall because it does not try. Still true and still blocking: the composed path has never been stepped on CUDA (no composed run in the run index). **Consequence for the roadmap: the remaining work is not 34 connectors, it is ONE coupled solver in the engine** — and it has to be one that survives the conditioning that already defeats the incumbent, i.e. the backbone-aware preconditioner belongs in `SOLVE_COUPLED`. Corrected by Lead after quoting the stale version of this entry twice in one session.
<a id="c5"></a>

**5.** **"22 CUDA tests pass on the A5000"** (= 6 CUDA-gated + 10 fake-device + 6 host), **"cortex/membrane seam CUDA_UNIT"** (`test_surface_body.py` has no CUDA-gated test), **"interior column native-runs 40 steps"**, GATE A's 3-seed and GATE B's 5-seed ensembles — no artifact for any of them.
<a id="c6"></a>

**6.** **`engine_reval` "19/21 native PASS" / "28/28"** and its **P3.1 percolation PASS** — the ledger is unreconstructable (20/28 rows record no device; P7.F.5/6/7 are `device: cpu`, one at 4% of native), and P3.1 tests giant ≫ 2nd rather than absolute connectivity on a cortex silently built with 300 crosslinks / 100 myosin, i.e. a gate a broken build cannot fail.
<a id="c7"></a>

**7.** **"6/6 ECM materials real-Pa"** — four of six recover a material *input* to ~27% against 100–400×-wide bands; the fibrillar law has concentration exponent ≈1.05 vs literature 2.05 and a **wrong-sign** N₁, both inherited by `ac/engine/ecm_mechanics.py`'s default collagen card.
<a id="c8"></a>

**8.** **"Physiological crawl 28–45 nm/s"** — coarse only; native is 0.12 nm/s and the two-way S6 artifact gives `v_crawl` 0.0038 nm/s. There is no crawl and no motility gate in `ac/` at all.
<a id="c9"></a>

**9.** **The dcm July claim set** — aggregate-σ compaction −8.1%, N=2000 large-dt stability, the T1-KMC heroes, "biology-time gap bridged", spreading A/A₀ 1.62–1.68: no committed numeric artifact, or disowned by their own docs, or a superseded pre-fix FAIL (use 3.38/3.62/3.44).
<a id="c10"></a>

**10.** **"266k particles in the native MCF7 cell"** and the "~530× short" line — both retired counts/results, and both still on the public gh-pages gallery.
<a id="c11"></a>

**11.** **The two KB self-reports.** *Benchmark "C4 = 100%"* — the gold answers and the C4 context are the same SQL over the same DuckDB, and the run is no longer reproducible. *CLAUDE.md's scale + citation figures* (329 rows / OK 177; vault 640; edges 808; PDFs 95; chunks 4054) — actual is 565 / 327 and 1648 / 1755 / 278 / 12810, and `CHECK` is ~95% a year-parsing regex artifact, not a citation-quality statement.
<a id="c12"></a>

**12.** **The `lane_d` D0/D1 gate results** (`ac/mt-if-linc-vertical@d2f5ba95`'s `outputs/lane_d/`, "6 CUDA-green gates") — the RIG and its gates were ported in on 2026-07-28 (PI D6 = port, not merge) but **its result artifacts deliberately were not**; they remain on the branch. The scenes are coarse diagnostic rigs by construction, so no magnitude from them is quotable however green the gates run — rung `CUDA_UNIT` at best, `QuantitativeClaim.BLOCKED` unconditionally. Porting rather than merging also caught a hard-coded CUDA device ordinal that a merge would have carried in silently.
<a id="c13"></a>

**13.** **"HOOMD is never imported or executed" — CLOSED 2026-07-28** (PI: delete, not archive; `archive/` leaves a runner runnable). 320 files deleted across `2b666f43` / `00858821` / `d06f7ea2`, including `h1_h2_vis.py` once the viz-entry-point decision retired the charter clause naming it. Outside `archive/hoomd_legacy/` the importer count is **zero**, verified by AST, while `dcm/fixtures/warp_parity_results*.json` + the 15 `tests/dcm` files survive their producer and still pass. Residual caveat only: HOOMD 7.0.1 is still installed in the dev env, so the last word is "not installed", not "cannot run". **FINISHED 2026-07-29 (PI: "HOOMD 관련 내용은 전부 삭제")** — the 158-file / 4.3 MB `archive/hoomd_legacy/` tree that the 07-28 pass had kept as a "read-only port specification" is now DELETED, together with `dcm/fixtures/warp_parity_results{,_cuda0}.json` and the 15 `tests/dcm` files that read them. **What that cost, stated because it is the real price of the instruction: 25 declared parity results are now permanently uncheckable** (`warp-B1-baoab-*`, `warp-B2-radial-*`, `warp-MSHAKE-*`, `warp-Fixman-*`, `warp-compartment-*`, `warp-DCM-*`, `warp-lamellipodium-tether`, `warp-B4-differentiability`, and their five `warp-gpu-*` counterparts). They were removed from `results_manifest.yaml` rather than relabelled: `needs-regen` would read "regenerate + commit", which is impossible because regenerating requires the runtime this repo forbids, and `retract` would read "disk contradicts the claim", which is not what happened — the evidence was deliberately removed. A manifest is a list of CHECKABLE claims; these are not, and their history belongs here instead. **What was deliberately NOT deleted, because deleting it would UN-forbid HOOMD:** `tests/ac/test_warp_only_contract.py` (the AST guard), the `conftest.py` import hook, and the "HOOMD is never imported" assertions in `ac/` docstrings — these name HOOMD precisely because they prohibit it. Nor the `references/analysis/` paper summaries, where HOOMD appears as a fact about a published study; editing those would falsify a literature record. Ported-topology provenance in `ac/motor/minifilament_{topology,warp}.py` and `params_i0b3.yaml` was re-anchored from a dead path to `git show 1a9ded66^:...` so the derivation stays auditable.
<a id="c14"></a>

**14.** **Every population-DEPENDENT number from the 2026-07-28 `observe` run** — the Gershgorin-over-true-`λ_max` factor **1.871×**, the condition number **2.03e11 / 11.31 decades**, the floppy-mode **fraction 439 of 2,712 DOF**, and the whole Fisher block (**effective dimension 5.42 / 8**, sloppiness 2.72 decades, α-actinin ranked least identifiable). Retired from (b) on **2026-07-28 by PI pushback**, one session after landing. All four are functions of the connectivity graph, of system size, or of the slice's own composition, and were measured on **904 nodes = 0.18% of the 511,114-node native cell** with 4 dorsal-arc joints against native's 1,413,720 crosslinks. `λ_min` in particular is a GLOBAL mode, so a small system UNDERSTATES the separation — the implicit-solver argument survives, its number does not. Precedent that slice results do not transfer here: (b)'s fiber-quotient negative result, which passed synthetically and failed at native. The instrument is sound; only the denominator was wrong. Re-measure at native (needs Lanczos — dense column probing refuses above 8,192 DOF by design), and restate the floppiness as a PER-MINIFILAMENT count, which is a single-object property and does transfer. **Both unblocked 2026-07-28 (`49d2215a`):** the per-minifilament restatement is measured and now sits in (b), and `observe/spectrum.py::lanczos_extremal_spectrum` gives the matrix-free path to the remaining three — so they are re-measurable rather than merely withdrawn. Honest limit of that path, recorded in the record itself: `λ_max` comes out exact, while the solve-free folded small end often does NOT converge, and it reports resolution per end via Parlett's bound instead of quoting the iteration count as physics.

<a id="c15"></a>

**15.** **The `nmii_sf_motor` cost multiplier "1.24 → 24.9 s/step, ≈20× for six decades"** — withdrawn the day it landed — **the lane caught it first** (`ac62041b`, 16:57) with Lead review converging independently four minutes later, which is the coordination working rather than a catch — and unlike (c) 14 the reason is not population but **confounding inside the slice**: the cheap and expensive runs differ in the solver budget the operator set, not only in the residual reached. Across all nine `config.argv`, wall-clock tracks `cg_iterations × newton` and is nearly independent of residual — budget 800 → 0.313 s/step, 3,200 → 1.170–1.267, 64,000 → 24.4–25.1. The quoted ~20× IS the 20× budget rise (`cg` 400→2000 × `newton` 8→32) and bought only 2.2× residual (1.311e-04 → 5.891e-05 at the same `implicit_step_scale`). **The same data read correctly points the other way and is the better result:** at fixed budget 3,200, `implicit_step_scale` 1e2 → 1e6 moved the residual 1.035e+00 → 1.311e-04 — **four decades for 8% wall-clock**. The cheap lever is the IMEX step scale, not the iteration budget. Two further blocks: all nine gates are `VOID`, so the records had already stamped all nine timings `comparable: false` — the row quoted a ratio the artifacts refuse to compare — and 904 nodes with 4 α-actinin joints cannot carry a cost law when conditioning worsens with density (2026-07-22e: 8k 2× → native 8×). The corrected statement is a (c) entry too. Re-measure at native; `ac_resting_cost_curve.py` does this on the resting lane.


<a id="c16"></a>

**16.** **The 2026-07-28 blebbistatin dose response, ALL SIX points** — "myosin is not a passenger: `k_xb` 1000 → 100 drops γ 4.32 → 3.47 at fixed 98.3% bound", and the non-monotonic `k_xb`=3 point offered alongside it as a protocol artifact. The commit that landed them recorded that **one** point (`k_xb`=1000, "tail drift 0.6 std") had settled and four had not. Re-judged 2026-07-28 (night, 4) through `observe/stationarity.py`, which is now wired into the driver: **none of the six settled, and none of them could have.** Two independent reasons, either sufficient. (i) The driving time was wrong: the driver stamped `bound_head_lifetime_s = 1/0.4 = 2.5 s`, but this repo's `NMII_KOFF0` is **0.35** and the sweep ran `--catch-slip`, whose zero-load off-rate is the SUM of the two Pereverzev pathways, `k_catch0 + k_slip0` = 0.70/s → **1.43 s** — so the required analysable span is 7.14 s and every point ran **0.60 s**. No `min_windows` ≥ 1 admits any of them; the verdict is `TOO_SHORT`, not `DRIFTING`, i.e. these series cannot even be said to have been going somewhere. (ii) The k_xb=1000 tail drift is **1.46σ** over the second half, not 0.6σ. And the reading itself does not survive: the three softest points, `k_xb` = 3 / 1 / 0.1, give γ = 4.365 / 4.326 / 4.291 against `k_xb`=1000's **4.321** — i.e. in the limit where the crossbridge is four decades softer the tension is the SAME, which read literally says myosin IS a passenger at 0.6 s and the 100/10 dip is transient. **Both readings are void**; what is quotable is the protocol finding, which is (b)-eligible and independent of any magnitude: a fixed step count is the wrong protocol for a sweep whose swept parameter sets the relaxation time. The correct ratio `R = γ(k_xb→0)/γ(1000)` against the experimental 0.1–0.5 band is NOT computed here and must not be inferred from these numbers.


<a id="c17"></a>

**17.** **Every cortex-motor γ measured at `dt` = 0.01 — the magnitudes AND the ratios built from them.** A `dt`-independence test at the one physiological dose (`k_xb`=100, whose 5.0 nm working stroke is the only one inside the model's 5–20 nm window) was queued with its prediction committed first: γ(dt=0.005) must equal 3.4807 within 3 combined sems. **Refuted.** γ(dt=0.005) = **3.1998 ± 0.0011** against **3.4807 ± 0.0017** — **137 σ, 8.07%** — with both points STATIONARY over comparable physical time (8.8 s vs 8.0 s), so it is not a duration artifact. `dt` = 0.01 is **not time-converged**, and every γ this lane has reported was measured there. **The shift is not in the kinetics**: decomposed, γ_network moves **8.13%** (135.7 σ) while `bound_fraction` moves **0.25%** and γ_source **1.16%** — the entire `dt` dependence sits in the mechanical relaxation term, and `max_f_cortex` moves 5.61% with it. **Why the ratios die too, which is the part that is easy to get wrong**: the standing argument for quoting `R` was that a ratio of one observable to itself survives the per-parameter PI-GAPs. That is true of an error that CANCELS. Nothing shows the `dt` error is equal at two different `k_xb` — it has been measured at exactly one — so the cancellation is an assumption, not a result, and `R` is blocked by this entry independently of (e) 2's separate finding that `k_xb` is the wrong knob. **What is NOT blocked**: the engine-cortex A/B (2.06e-16 over a settled trajectory) and the cross-build parity (bit-identical to 15 digits) are comparisons at the *same* `dt` between two code paths, so a common `dt` error cancels there by construction — this entry constrains absolute values and the ratios across `k_xb`, not whether two implementations agree. **Consequence for the roadmap**: `ROADMAP.md` stage 2 asks for `dt` to be **RAISED** and the answer shown `dt`-independent; `dt` cannot be raised from a point that is not converged, so stage 2 is further away than the plan assumed, not closer. ⚠⚠ **RE-JUDGED 2026-07-29 AGAINST A MEASURED SEED SCATTER, AND MOST OF THE SIGMA FIGURES BELOW WERE INFLATED.** Every run in this sweep used `seed = 0`, so the run-to-run scatter of the settled γ had never been measured and every σ quoted here used the WITHIN-RUN `sem`. Three replicates at `k_xb`=100, `dt`=0.01 (identical `t ≥ 2.0 s` window) give **3.48131 / 3.46411 / 3.35797** → scatter **S = 0.0668 pN/µm = 1.95% of γ, which is 67× the within-run sem (0.0010)**. The uncertainty on the difference of two single-seed values is `√2·S = 0.094`. Re-expressed: `k_xb` 10 vs 3 = **99 σ**, 1000 vs 100 = **8.9 σ**, 100 vs 10 = **6.3 σ** (was "168 σ"), `dt` 0.01 vs 0.005 = **3.0 σ** (was "137 σ"), `dt` 0.005 vs 0.0025 = **2.3 σ**. **CONSEQUENCES: (i)** the `k_xb` findings are unaffected and the runaway is overwhelming; **(ii)** the `dt` non-convergence is REAL but marginal, not overwhelming; **(iii)** ⚠ **the non-monotonicity claimed below is RETRACTED — at 2.3 σ it is not distinguishable from seed noise**, and with one seed per `dt` point "non-monotone `dt` dependence" and "monotone `dt` dependence plus scatter" cannot be told apart by these data, so "no order fits" and "Richardson does not apply" are withdrawn as unsupported; **(iv)** single-seed runs cannot resolve γ differences below ~5%, which is a standing limit on every sweep this project runs, stage 4 included. **NOT affected**: the engine-cortex A/B and the cross-build parity compare two code paths at the SAME seed and config, so the scatter cancels by construction.

**SEPARATED — it is (B), the outer integrator.** Two readings were possible and both were declared before the separating run: (A) relaxation depth — halving `dt` doubles inner iterations per physical second, so the reported tension would contain unrelaxed residual, the same failure mode as (c) 3; (B) outer time-integration order. A run at `dt`=0.01 with `--outer` **40 → 80** matches (A)'s relaxation budget without touching `dt`. **Doubling the budget moved γ by 6.0e-08 relative** (3.4807438591 → 3.4807439934) against the 8.07% the `dt` halving produced. The flag is plumbed and the work was done — inner iterations 32,000 → 64,000, wall ×1.59 — and `max|PF|` **did** improve, 0.5291 → 0.4998 (5.5%); `bound_fraction` and `n_bound` are bit-identical. **So the residual falls with more iterations while γ does not move.** Two consequences: γ is **not** reading unrelaxed residual — direct evidence for a conclusion that had rested only on `max|PF|` staying flat, a weak argument that happened to reach the right answer (`max|PF|` is a max over nodes, γ_network an integral, and the max can hold while the integral moves) — and **the `dt` dependence lives in the OUTER physical-time integration, not the inner solve.** That last point redirects effort: this lane's convergence work has been inner-solver conditioning (fiber-arclength multigrid, ERM-Schwarz preconditioner, IMEX step scale) and **none of it addresses this**. **The order was tested and BOTH predictions failed, so it is not merely unknown — it is not defined by these data.** `dt`=0.0025 was run with both predictions declared first (first order ⇒ 3.0594; second order ⇒ 3.1296). Measured, after correcting a detector defect described below: **3.4148**. The sequence 3.4813 → 3.2003 → 3.4148 (all re-measured on an identical `t ≥ 2.0 s` window, since each run's auto-selected transient differed) **is not monotone** — it falls then rises — so no convergence order fits and Richardson extrapolation does not apply. ⚠ **RETRACTED: the "either way the `dt`=0.01 values are 10–19% high" quantification above.** It assumed first or second order; both are refuted, so the true `dt`→0 value is simply **unknown**. That does not weaken this entry, it strengthens it — an error of unknown size and unknown sign cannot be corrected for. **A detector defect surfaced in the same run and is fixed**: `equilibration_point` maximises `n_eff = T/(2·τ_int)` and chose `start = 0` on a series that rose 0.8 → 3.4 in 0.5 s, so the "steady-state mean" was taken ACROSS the rise (3.373 rather than the plateau's 3.415). It self-reinforced — keeping the transient inflated `std` twentyfold, and the drift verdict is `|drift| < σ·std`, so a bigger transient makes the test EASIER to pass; that run was certified STATIONARY at 0.94 σ. `window_opens_on_a_transient` now compares the retained window's leading tenth against its trailing tenth, normalised by the TRAILING tenth's scatter (transient-free by construction, so it cannot inherit the same loop), and returns DRIFTING with no mean when they differ by more than one standard deviation. It reads 0.0 / 0.8 / 33.5 on the three runs. The first draft tested the first SAMPLE and rejected all three — a transient is a property of a segment, not a point.

<a id="e4"></a>

**(e) 4.** **`STRUCTURE.md` names directories that no longer exist.** Measured 2026-07-29 by resolving all 25
backticked directory entries against four bases (repo root, `aleph/`, `aleph/ac/`, `aleph/docs/`):
16 resolve, **9 do not** — `cell_engine/`, `warp_port/`, and seven under the heading "1. `aleph/` —
historical 2026-05-20 HOOMD layout (non-executable)": `aleph/{ecm,cell,cortex,bridge,junction,integrator,warp_port}/`.
**Why this is a PI question and not a defect to fix.** Seven of the nine sit under a heading that says
*historical*, so listing directories that were later deleted may be deliberate record-keeping rather than
drift. Two things cut against reading it that way: `aleph/warp_port/` carries "adopted 2026-06-21 (G1 GPU
parity, G2 substrate)", which reads as live status, and `aleph/tests/` is described as "HOOMD runtime
테스트 (현재 `__init__.py`만)" while the real `aleph/tests/` holds `ac/`, `ff/`, `mechanics/`,
`coordination/` and ~20 top-level modules. **Not edited**: `STRUCTURE.md` is read-only during normal work
(`CLAUDE.md` §File ownership), and the 2026-07-29 PI instruction "HOOMD 관련 내용은 전부 삭제" is about
HOOMD content, which is not the same question as whether a historical layout section should be kept.
**Re-measure, do not quote this count** — it changes whenever the tree does. The checker:

```python
import re, pathlib
root = pathlib.Path(".")
txt = (root / "STRUCTURE.md").read_text(encoding="utf-8")
bases = ("", "aleph/", "aleph/ac/", "aleph/docs/")
for p in sorted(set(re.findall(r"`([A-Za-z0-9_./-]+/)`", txt))):
    if not any((root / (b + p)).exists() for b in bases):
        print("MISSING:", p)
```

<a id="closure2026-07-29"></a>

**Import closure of the production driver — re-measured 2026-07-29, and the parked trees stay unreachable.**
Recorded here rather than in `STATE.md` because it is a maintenance measurement, not a claim, and because
the boot budget has eight bytes ((e) 3). AST closure from `scripts/ac_gate_b_cortex_motor_native.py`:
**120 modules / 2,576,592 B**, split `ac/` 92, `ff/` 20, `scripts/` 5, `common/` 3. **Reachable modules under
`dcm/` or `archive/`: ZERO** — which is the number that matters, because the 2026-07-29 overnight cleanup
removed a seven-hop chain that dragged parked `dcm/` into the closure through two constants, and this
confirms it has not grown back. Kernel share: **51 of 120 (42%)** define `@wp.kernel` or call `wp.launch`.

**The 69 that fire no kernel are not dead weight, and this is the part worth writing down**, because "modules
that are imported but use no kernel" reads like a cleanup target and mostly is not. Audited by importer:
`common/{compartments,surface_manifold,turgor_pi0}.py` are physical constants pulled in by `ff/membrane_surface.py`
and `ff/turgor_constants.py`; the four `scripts/` entries are the **mandated** auto-viz path
(`ac_viz_common → ff_cell_morphology → ff_viewer_html`) plus `dump_state.py`, which exist because
`CLAUDE.md` requires every run to write its own figures beside its record. Removing them would break a hard
rule, not tidy anything. The rest are contracts, dataclasses and observability — kernel-free by design.

So the eager-import question is **CLOSED**: what remains has a justification chain per module. Re-measure
rather than quoting these counts; the checker walks `ast.Import`/`ast.ImportFrom` from the driver and
resolves names against the repo root.

<a id="gpu2026-07-29"></a>

**A second GPU host exists, and all three of its cards are verified — measured 2026-07-29.**
Recorded here, not in `STATE.md`, for the same reason as the closure above: it is a device inventory,
not a physics claim, and the boot budget has 30 bytes ((e) 3).

**The host.** `DESKTOP-EEUU4RQ` (Tailscale `100.104.168.13`, KAIST direct), owned by **함석훈**
(`tjrgns4576@gmail.com`) and **shared with five other users** — `ham`, `jimin`, `yerim`, `yujin`,
`nari` besides us. **WSL2** (`5.15.167.4-microsoft-standard-WSL2`) on Ubuntu 24.04.1, driver
`610.43.02`, CUDA Toolkit 12.9 / Driver API 13.3 as Warp reports them, 64 CPU threads, 62 GiB RAM
(Slurm `RealMemory` 54,579 MB), `/` 84% full with 40 GB free.

**The cards** (PCI instance IDs from the Windows device log, three distinct bus locations):
**2× RTX 4090** (`DEV_2684`, GIGABYTE, **sm_89**) + **1× RTX 3090** (`DEV_2204`, ASUS, **sm_86**),
24 GiB each. Contrast gbook, which is a single **RTX A5000 *Laptop* GPU, 16 GiB, sm_86** — the
charter and `STATE.md` both say only "A5000", which is narrower than the truth in both directions.

**Access is NOT ours to take, and is already governed.** `/dev/dxg` is ACL'd to `ham` alone, so no
CUDA context can be created directly (`cuInit -> 100`, `NO_DEVICE`, measured). The admin runs a
token system instead: Slurm **licenses** (`rtx3090`, `rtx4090_0`, `rtx4090_1`, Total=1 each) rather
than GRES — `Gres=(null)`, which reads as "unmanaged" and is not. `gpu-check` / `gpu-shell` /
`gpu-submit` allocate a token, `gpu-dxg-gate` grants `/dev/dxg` for the job's lifetime and revokes it
at exit, and `CUDA_VISIBLE_DEVICES` is set by the gate so the job sees exactly one card as `cuda:0`.
**A `gpu-bypass-watch` daemon exists**, so building our own lease for this host would be, functionally,
a bypass. `ffn_gpu.py`'s lease stays correct for gbook and must WRAP `gpu-submit` here, not replace it.

**Verified by three submitted jobs** (5 / 6 / 7; logs under `~sungwook/slurm_logs/` on that host),
each compiling and running an atomic float32/float64 scatter of 4.19 M pairs against a host reference:

| card | arch | f32 | f64 | f64/f32 | float64 rel. err |
|---|---|---|---|---|---|
| RTX 4090 | sm_89 | 2.247 ms | 3.851 ms | 1.71x | 1.64e-16 PASS |
| RTX 3090 | sm_86 | 5.393 ms | 8.415 ms | 1.56x | 3.27e-16 PASS |

**These are microbenchmarks of one synthetic kernel and may not be quoted as engine performance** —
no population, no build, no gate. What they license is narrower and is the point: Warp 1.14.0 resolves
CUDA under WSL2, **sm_89 is a target this project had never compiled for** and it compiles and computes
correctly, and float64 atomics cost only ~1.6-1.7x float32 *for a scatter-bound kernel* — nothing like
the 1/64 FP64 ALU ratio a consumer board's spec sheet implies. The force assembly is scatter-bound too,
which is why this was the shape worth probing first.

**The mixed-architecture cache trap is CLOSED, by experiment rather than by argument.** Warp keeps
per-arch PTX *inside* one module directory — `wp___main___2df8541.sm86.ptx` beside `...sm89.ptx` — so
the arch is a FILE component, not a directory one, and a check for an arch-tagged directory (which is
what I wrote first) asks the wrong question and answers False. Run in the order 89 -> 86 -> 89, the
third job hit the cache (`4.11 ms (cached)` against `1317.85 ms (compiled)`) **and stayed correct at
1.64e-16**. So one cache directory across these cards is safe, and no cache separation is needed.

**One retraction from this session, since it is exactly the failure this file exists for.** The first
probe reported `float64 verdict: FAIL` at relative error `5.000e+00`. That was the PROBE, not the card:
the timing loop accumulated six launches into the array the correctness check then compared against a
one-launch reference, so the exact 5.0 was the kernel being right six times. Fixed by zeroing and
relaunching once after timing. A verification harness that cannot fail cleanly is worth no more than
the gate it feeds.

<a id="gpuspeed2026-07-29"></a>

**RTX 4090 vs the A5000, measured on the REAL engine at full native — 2026-07-29.** Four runs, identical
code (`384dc1f3`), identical command, `seed=0`, **70,686 cortex filaments / 494,802 actin / 1,165,834 total
nodes / 8,840 heads, `fraction_of_native = 1.0`**, 50 accepted steps, `k_xb`=100 (the one physiological dose),
`dt`=0.01, `n_inner`=40. Jobs 8/9 on the workstation's 4090 (sm_89), two leased gbook runs on the A5000
Laptop (sm_86).

| telemetry | 4090 s/step | A5000 s/step | ratio |
|---|---|---|---|
| every step (`--gamma-every 1`, the default) | 6.847 | 9.776 | **1.43x** |
| every 25 steps (`--gamma-every 25`) | **1.124** | **4.655** | **4.14x** |

**The two rows are the finding, not the second row.** The synthetic atomic-scatter probe ([gpu2026-07-29])
put the 4090 at 4.25-4.52x. The first real run said 1.43x and looked like the probe overselling the
hardware. It was not: at the default stride the step is **telemetry-bound**, and telemetry is host work
that both machines do at the same speed — **5.72 s/step on the 4090, 5.12 s/step on the A5000**, i.e.
machine-independent to within the hosts' own CPU difference. Stride it out and the device ratio returns to
**4.14x**, which agrees with the probe. The `timing.device_note` in every record already said "at native it
dominated the per-step cost"; this measures how much: **84% of the 4090's default step and 52% of the
A5000's.**

**The free 2.1x.** Raising the telemetry stride costs no hardware and no physics: on the A5000 alone it is
9.776 -> 4.655 s/step. `ROADMAP.md`'s sweep budget is denominated in s/step, so this moves that arithmetic
before any new card is used. What the stride trades away is sampling density of `gamma_*` and
`bound_fraction`, which a stationarity judgement needs — so the right stride is a function of the
correlation time, not a constant, and picking one is a contract question, not a tuning knob.

**NONE of these four timings is quotable as the engine's speed**, and the records say so themselves:
`timing.comparable = False`, `not_comparable_reason` = "the gate verdict is absent, so this run did not
establish the convergence its cost would have to be quoted at." All four are `TOO_SHORT` on every
stationarity observable (0.5 s physical = **0.175 bound-head lifetimes**). What they license is a RATIO
between two devices running provably the same work, and the free-stride result — not a cost-per-converged-step,
which is what stage 2 still asks for and still does not have.

**Cross-architecture agreement, measured rather than assumed.** Comparing the 50-step trajectories field by
field, sm_89 vs sm_86: `n_bound` is identical at every step, the census is identical, and every float64
observable agrees to **max relative 5.10e-09** (248 of 450 scalars differ at all). That is far above float64
round-off and is the expected amplification of atomic accumulation order through 2,000 inner iterations —
the same "float64 atomic arrival order" effect `engine_cortex_ab` names. **The two cards do the same physics
to nine significant figures over a full-native trajectory.** ⚠ First read of the console table suggested
bit-identity; that was the table printing four significant figures, and the field-by-field diff is what
corrected it. A comparison whose resolution is coarser than the effect it reports proves nothing.

<a id="native2026-07-29"></a>

**`fraction_of_native: 1.0` counts the CORTEX ONLY, and the membrane behind `ROADMAP.md`'s cost
arithmetic is at a quarter of the driver's own default — found 2026-07-29 while reconciling a timing.**

**The measurement.** One variable changed, everything else identical (4090, 50 steps, `k_xb`=100, seed 0,
`--gamma-every 25`, SLIP, 70,686 filaments):

| `--membrane-subdiv` | membrane verts | `n_total_nodes` | s/step | `fraction_of_native` stamped |
|---|---|---|---|---|
| 7 | 163,842 | 674,314 | 0.5779 | **1.0** |
| 8 (the driver default) | 655,362 | 1,165,834 | 1.1244 | **1.0** |

4x the membrane mesh -> 1.73x the nodes -> **1.95x the step cost**, and **both runs stamp
`fraction_of_native: 1.0`** because line 837 computes it as `args.filaments / 70686.0` — cortex filaments
and nothing else. A run at subdiv 6 (40,962 verts) would stamp 1.0 too. The census field that exists to
stop a reduced run being read as native does not see the compartment that CLAUDE.md's own architectural
principle names as the cost: "the cost is the fluid grid and the moving membrane boundary."

**What it reconciles.** The five settled runs `ROADMAP.md` quotes at **2.87-2.90 s/step** ran
`--membrane-subdiv 7`; the default was 8 at their build (`efe74e52`) too, so it was passed deliberately.
Today's A5000 run at the default gives 4.655 s/step; divided by the measured 1.95x membrane factor that is
2.393, against the settled 2.900 — a 1.21x residual in the right direction, since the settled runs also do
5x more telemetry readbacks (`gamma_every` 5 vs 25) and use CATCH_SLIP rather than SLIP. **The arithmetic
closes.** ⚠ The 1.95x factor is measured on the 4090 only; it is not yet confirmed on the A5000.

**⚠ CORRECTED WITHIN THE HOUR — the framing above about WHICH run is reduced was backwards, and the
correction is the more useful finding.** `CellConfig.membrane_subdivisions` (`ac/cell/assemble.py:267`)
carries a **PI ratification dated 2026-07-22**: level 3 = 642 verts is the smooth resting sphere kept for
backward-compat and gamma-validation; **6 = 40,962 verts (~131 nm) is the PRODUCTION BASELINE** for dynamic
morphology (blebs ~0.5 um give ~4 nodes per neck; ERM 58/um^2); **7 = 163,842 (~66 nm) is the VALIDATION
resolution** (~8 nodes/neck, ERM 235/um^2, physiological). So the settled runs' `--membrane-subdiv 7` is
**the ratified validation resolution, not a memory reduction** — they were right. The anomaly is the other
way round: **`ac_gate_b_cortex_motor_native.py` defaults to 8**, a level that appears in no ratification,
and its help text calls anything lower "an A5000 memory safety valve" while the ratified comment says in
terms that **"Memory is a non-issue (subdiv 6 full-native = 592 MB on the A5000 16 GB)"**. A production
driver whose default exceeds the ratified validation resolution, justified by a memory argument the
ratification explicitly rejects, is the defect here — and today's runs at that default are the ones that
need re-labelling, not the settled ones. `ROADMAP.md`'s 2.87-2.90 s/step therefore stands on the ratified
resolution and needs no (c) row on this ground.

**What survives the correction, unchanged, because it was measured rather than inferred:** the 1.95x
membrane cost factor; that `fraction_of_native` cannot see the membrane at all; and that two runs 1.73x
apart in node count both stamp 1.0. Those are facts about the census. What was wrong was my reading of
which configuration the census was hiding.

**Recommended to PI, not taken unilaterally, because declaring a `ROADMAP.md` budget figure non-quotable is
a substantive call:** (i) a (c) row against "2.87-2.90 s/step at full native"; (ii) `fraction_of_native`
either renamed to `fraction_of_native_cortex` or computed across compartments — a scalar that says 1.0 for
two runs 1.95x apart in cost is worse than no scalar; (iii) a PI answer to what the membrane's PHYSIOLOGICAL
subdivision is, since neither 7 nor 8 is sourced anywhere — 8 is a default, not a derivation.

**Retracted from earlier the same session:** the hypothesis that the gap was in-loop Warp kernel
recompilation on a cache my own rsync invalidated. Wrong, and it was reasoning where a config diff was
available for free. `wall_seconds` is `loop_s`, timed from after the cell is built
(`ac_gate_b_cortex_motor_native.py:723` vs `:580`), so setup was never in it either.

<a id="popkinetics2026-07-29"></a>

**Population kinetics readiness, per component — audited 2026-07-29, because "who signs the MT count?" was
the wrong question.** The counts (MT, IF, filopodia, FA) are OUTPUTS the sweep infers, not inputs a PI
ratifies: `CLAUDE.md` §What this is says "Parameters are OUTPUTS; literature supplies priors", and in a real
cell those numbers are not fixed anyway — MT count moves with dynamic instability, filopodia form and
retract, FA assemble and disassemble. `ac/engine/population.py` already refuses the wrong framing in code:
"**count changes only via the free-list, never by resampling to a target** … there is no 'draw N filaments
to hit a density' path". So a component's population is admissible only if EVENTS can move it, and what
needs a literature prior is the RATE, never the count.

| component | mechanism | candidate prior | what is actually missing |
|---|---|---|---|
| `ecm` | events wired — the only component calling `PopulationLedger.allocate/release` (`ac/ecm/mikado_topology.py`, `remodel_runtime.py`) | — | nothing structural |
| `microtubule` | **stochastic DI implemented** — `_mt_di_stochastic_proposal_kernel` + commit/rollback kernels + `propose_dynamic_instability_stochastic` | **`DynamicInstabilityRateCard.rusan_2001_llcpk_proxy()`**, `ratified=False`, provenance recorded | **load dependence**: the kernel takes the two rates as SCALARS and no force array, so rates cannot respond to tip load |
| `intermediate_filament` | mechanics delegated (same pattern as MT bending) via `LegacyIntermediateFilamentCageAdapter` | **NONE — and the corpus does not contain one** | a sourced turnover prior; see the negative result below |
| `protrusion` (lamellipodium / filopodium) | 1 kernel, 9 gap markers, 2 kinetics mentions | — | initiation/retraction events |
| `cortex`, `nmii`, `sf_arc` | static populations built once | — | no birth/death path |

**Negative literature result, recorded so it is not re-searched.** `IntermediateFilamentTurnoverCard` is an
UNSET slot whose own docstring says no IF-specific `k_off0` / Bell force scale is sourced. Searched the KB
content layer directly: chunks mentioning vimentin/keratin/intermediate-filament together with
exchange/turnover/half-life/FRAP/recovery return **9 candidates and ZERO extractable rates**. `KB-DRAFT-3.B-11`
carries "subunit exchange timescale = minutes (order)" but is status `draft` with an EMPTY evidence field, so
it cannot seed a candidate card. The IF prior is a literature-ACQUISITION task (Lorenz 2019/2023 are recorded
as paywalled), not a code task, and MT's `rusan_2001_llcpk_proxy` is the pattern to copy once a source exists.

**One stale comment corrected in the same pass.** `microtubule_rig.py`'s header said "catastrophe/rescue RNG
transitions **and** load dependence are deferred". Half of that was false — the stochastic kernel has existed
since the rig landed. Verified by reading the kernel signature rather than the note: the transitions are
there, the force argument is not. The note now says which half is which.

**Why this matters beyond bookkeeping.** The composed world assembles **14 components / 36 connectors across
13 connector families** and its disjointness invariants hold — the wiring is genuinely complete, exactly as
(c) 4 says. But `ac_composed_world_dump.py`'s own docstring marks node GEOMETRY and the non-cortex counts as
"**small synthetic render placeholders, NOT physiological claims**", and the dump bears that out: cortex
appears as 4,000 nodes against the 494,802 a native cortex build produces. The placeholders are not there
because someone failed to type a number — they are there because, for every component except ECM, **no event
can create or destroy a member**, so there is nothing for a population to emerge from.

---

<a id="c18"></a>

### c18 — the ECM row's "zero crosslink bonds" is about the code, not the configuration

**Do not quote: "zero crosslink bonds" as the ECM tier-(a) row's population.**

The row reads *"200–220 collagen fibers, **zero crosslink bonds**"*, and that zero has been read as a
scoping choice — a deliberately un-crosslinked bed, existence-and-sign only. Session `97fb3916`
found on 2026-08-09 that it was not a choice.

`ECMWorld` has carried a `crosslink_connector` field and validated its name, its internal endpoints
and its adjoint requirement **for as long as the row has existed — with nothing satisfying the
Protocol.** There were zero crosslink bonds because no implementation existed to make one, not
because the measurement excluded them.

**Why this is a (c) entry and not a correction.** The number is still true: that run had no crosslink
bonds. What may no longer be quoted is the *population descriptor*, because a reader takes "zero
crosslink bonds" as a statement about the configuration under test and it is a statement about the
code. A missing implementation leaking into a population is the same shape as
`require_complete=True` having only ever been satisfied by `_CensusRuntime` stubs — a contract
reporting satisfaction it never obtained.

Requoting requires a run on `ecm_crosslink` with a real runtime (`94775e5d`), with its own population.

**Source:** session `97fb3916`, lane `engine`, while authoring the `ecm_crosslink` connector — the
only one of five that was a genuine authoring job rather than a rewiring.

<a id="c19"></a>

### c19 — the MTOC sits at the centre of the nucleus, so no MT tensegrity result on this cell is about a cell

✅ **PUSHED 2026-08-21.** For part of one day this entry was pull-only, because `STATE.md` had 39 free
bytes and the row needed 141 — PI queue item 21, the collision between §(c)'s *"only GROWS"* and a fixed
byte cap. The PI ruled the cap to EXCLUDE §(c) rather than raising it, and the boot total was paid for by
retiring five superseded memories rather than by shortening this list. The row is in §(c).

**Do not quote: any MT-based tensegrity, MT↔nucleus, or MT↔cortex result measured on this cell.**

Measured on `aleph/outputs/ac/world_phase4/tau3_seed2.alephcell`, frame 0, radii about the origin [µm]:

```
microtubule            74,001    r 0.000 .. 7.400    centroid (0,0,0)
nuclear_envelope       40,962    r 5.100 .. 5.100    exact sphere, on the origin
intermediate_filament   2,460    r 5.250 .. 7.250    -- starts OUTSIDE, correctly

microtubule nodes inside the nuclear envelope:  50,735 / 74,001  =  68.6%
```

`build_microtubules`' Args block calls `centre` *"the MTOC position"*; `populations.py:104` passes
`(0.0, 0.0, 0.0)`, which is the geometric centre of the nucleus. A centrosome is a cytoplasmic
organelle. `aleph/laws/microtubule.py:40` states as physical fact that a microtubule *"cannot
interpenetrate the nucleus or other filaments"*, and line 10 that the aster *"positions the nucleus"* —
which a centrosome at the nucleus's own centre has no direction to do.

**Why this is a (c) entry and not a correction.** No number is wrong. γ is method-of-planes on the
cortex and microtubules carry no bound law in these runs, so tonight's γ and τ are unaffected — session
21 confirmed this on both the host and device paths rather than asserting it. What may not be quoted is
any conclusion *about the MT compartment as a mechanical element*: the strut's base is inside the object
it is supposed to brace, so a result on this geometry is a result about a configuration that does not
occur.

**Why nothing caught it.** `assert_inside_membrane` asks whether every node is inside R = 7.5 µm, and
7.400 is. `assert_partitioned` is a property of ID ranges. No test compares the two populations. No law
evaluates MT–nucleus sterics, so the interpenetration costs no energy and raises nothing. **No number in
the tree contradicted any other number** — it was found by rendering the cell and measuring what the
picture suggested. Full working: `aleph/docs/THE_CENTROSOME_IS_INSIDE_THE_NUCLEUS_2026-08-21.md`.

**When it may be quoted again.** When the archetype is ruled (item 3) and the MTOC is placed off-centre,
adjacent to the envelope, with `build_microtubules` taking `R_nuc_um` the way
`build_intermediate_filaments` already does. Not before: moving the MTOC changes the built cell.


<a id="c20"></a>

### c20 — the NMII minifilaments stand clear of the actin, so nothing built on this cell is contractile

✅ **CLOSED 2026-08-22 for builds from `ba0160f9` onward — and this closure is NARROWER than c13's.**
c13 closed because its claim turned out to be TRUE when checked. This one closes because the
CONFIGURATION WAS REPAIRED, which is a different act: the PI ruled queue 14 (a), the excursion guard
was made to promise something true, and `geometry.py`'s second derivation of the NMII shell was
deleted so `build.cortex_shell()` owns it. Re-measured with no override on the built default:

```
shell radius 7.40 = the cortical shell   nmii span [7.300, 7.500] = cortex span
nearest actin   min 0.0028   p95 0.0289 um       stations @ reach 0.20   442 / 442
nodes outside the membrane   0
```

⚠ **RESULTS MEASURED ON THE CELL AS IT WAS BUILT BEFORE `ba0160f9` REMAIN BLOCKED.** Repairing a
configuration does not make the numbers taken from the broken one quotable — it only stops new ones
joining them. Anything whose build commit predates `ba0160f9` is still covered by everything below.

**Do not quote: any NMII binding, crossbridge, cortical-motor tension or contractility result measured
on this built cell** — including any `gamma_source` read from it.

Measured at full native on 4090-1, build `13233392`, closure 22/22 hash-matched against the run host
before launch. 442 minifilaments / 32,708 NMII nodes; ten-point radius ladder declared before the run:

```
shell radius   nearest cortex node   stations @ reach 0.20   assert_inside_membrane
   6.95 (built)   min 0.2501 um            0 / 442            PASS
   7.20           min 0.0058               422 / 442          PASS
   7.35           min 0.0017               442 / 442          PASS
   7.40           min 0.0013               442 / 442          REFUSED, 414/32,708 nodes out
```

**At the built shell the CLOSEST minifilament in the population is farther from actin than the
`hand_kmc` capture proxy.** Not "most cannot bind" — none can, and no rate constant reaches across a
gap. So a contractility result on this cell is not weak evidence; it is evidence about a cell whose
motors are not touching what they pull.

**Two things this entry does NOT say.** It does not retire the mechanism — `world/` never bound a
crossbridge here, so there is no magnitude to withdraw, only a configuration to stop reporting against.
And it does not name a replacement radius: the two constraints leave a window, and choosing inside it is
PI queue item 14.

⚠ **Upstream of that choice, and the reason the window's edge is not where it reads.** `build/nmii.py`
raises when `excursion > thickness/2` — a FIT test — while `nmii.py:156` draws the per-filament radius
over the FULL thickness with no inset, so a filament at the outer edge sits that excursion past the
shell whatever the guard returns. The measured overshoot is CONSTANT across the ladder and SMALLER than
the guard's own prediction: the number was never wrong, the question was. Same family as `balance_ok`
and `descent_ratio`, which `geometry.py`'s own self-check already names.

⚠ **And PI queue item 14's arithmetic had the sign backwards** — it reasoned the heads land a few nm
INSIDE the membrane at 7.40. They land outside, and the gate refuses.

**Source:** Lead session 2026-08-21 PM, `615312cd`. Artifacts:
`aleph/outputs/ac/world_tier0/sweep_2026-08-21/` (10 records + `REPORT.md` + 1 fig),
`aleph/outputs/ac/world_phase1/phase1_inside_2026-08-21_r7.{20,35}.json`.

<a id="c21"></a>

### c21 — the PHASE 4 cell carried two of its eleven populations' laws, and its cortex was a tenth of actin

**Do not quote: any γ, τ, reliability, seed scatter or residual taken from a PHASE 4 run before
2026-08-24.** This is not a second reason to doubt the numbers already blocked by c3 and c17 — those
are about the STATISTIC (an induction period; time convergence). **This one is about the CELL**, and it
reaches every number the cell produced, including ones no earlier entry covers.

Measured 2026-08-24 on `world_phase4/tau3_seed1.json`, the record the three τ seeds were taken from:

```
populations standing                            11
populations carrying ANY law                     2   cortex (axial+bending), membrane (tension+Helfrich)
  microtubule   74,001 nodes / 73,500 segments   0 laws
  filopodium · lamellipodium · stress_fiber · sf_arc · intermediate_filament · chromatin   0 laws
  nuclear_envelope  81,920 faces                 0 laws
cortex k_axial                        88,000 pN/µm   = 1/10 of EA_actin / L_seg = 880,000
census.live.bond / grid_cell             0 / 0
residual at step 0                    2.3e-13 pN     built AT rest: nothing to relax
motion per step                          20.7 pm     thermostat only; the drift term is 414× smaller
motion over 30,000 steps                  3.58 nm    7.2% of ONE 50 nm bond
```

**Three separate things make the cell the wrong object, and each alone would be enough.**

1. **`build/nmii.py` never runs in production.** `build_all` and `build_remaining_populations` contain
   zero references to it; the only callers are two TIER driver scripts. **There is no motor in the cell
   the τ series was measured on.** The finding was already written down — in the docstring of
   `test_observe_gamma.py`'s `test_an_absent_family_reports_none_not_zero`, which exists so
   `gamma_source` reports `None` and not `0.0` because *"the motors contributed nothing"* and *"there
   are no motors"* are different findings. The distinction stood; what it distinguished reached no
   record. Closed by `21fae96f` (`builders_not_standing` + a ratchet).

2. **Nine of eleven populations were inert.** A microtubule with 73,500 segments and no axial or
   bending law is geometry, not a microtubule. `--bind` existed the whole time and every τ record shows
   `extra_bindings: []` — nobody passed the flag.

3. **The cortex was a tenth of actin.** 88,000 pN/µm was honestly recorded as a TEST POINT, so it is
   not a hidden defect; but the charter requires starting at the physiological operating point, and
   `EA_actin / L_seg` is 880,000. Raised by PI ruling 2026-08-24 (`f4349abe`).

⚠ **What this entry does NOT say.** It does not retire the machinery. The step ran, the partition
invariant held every step, the engine is bitwise deterministic, and γ is emitted on the resting path —
all of that is (b) and stands. It says the OBJECT those instruments were pointed at was not the cell
the conclusions were about. A run after the binding lands is a different cell and its numbers do not
inherit this entry; they also do not inherit c3 or c17's clearance, which were never granted.
