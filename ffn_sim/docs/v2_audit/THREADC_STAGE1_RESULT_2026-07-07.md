# Thread-C Stage 1 result — MT aster merged into the crawl driver's implicit K (2026-07-07)

Wires the verified Stage-1 building blocks (`concat_fiber_networks`, `merge_aster_into_cortex`, committed a9212ee)
into the production crawl driver `ffn_sim/scripts/ff_crawl_on_substrate.py`, per `THREADC_UNIFIED_PLAN_2026-07-07.md`
§2 Stage 1. The MT aster (n_mt tubes, κ=KAPPA_MT=20 pN·µm², radiating from one MTOC) now rides the **same**
implicit `(γ/dt·I + K)Δx = F` solve as the cortex at 300× the bending stiffness — no second solver.

## What changed

- **`build(microtubules=, n_mt=, L_mt_um=)`** — builds the aster + `merge_aster_into_cortex`, returns the merged
  elastic net + `Ne`/`mtoc_idx`/hub-crosslink arrays. `pos_all = [cortex(Nc) ; MT_arms(Nmt) ; MTOC(1) ; nucleus(n_nuc)]`.
  `microtubules=False` → no-op passthrough (Ne==Nc, empty hub) → OFF path bit-identical.
- **`run()` `Nc→Nc/Ne/N` split** — cortex surface physics (faces, turgor `F[:3Nc]`, substrate, clutch, protrusion,
  gravity, centroid/mean-R reductions, membrane proxy) stays on `Nc`; bending/crosslink/reshape/K/γ/pos span `Ne`;
  solve vector `N=Ne+n_nuc`. Hub crosslinks (MTOC↔arm-base, finite rest = seg_um) concatenated into `xl_ij/kxl/r0`
  (host + device + cupy). Cortex `ConvexHull`/V0 taken from `net.pos[:Nc]`. Crosslink turnover excludes the hub;
  actin assembly / barbed-end growth restricted to cortex fibers; reshape covers all fibers (MT arms inextensible).
- **MT-node drag (R-γ / R7)** — the GPU solver uses a scalar `gamma_rep`; MT arms + MTOC (nodes `[Nc,Ne)`) get the
  correct fiber-log/Stokes drag via `diag[3i(+1,+2)] += (γ_i−γ_rep)/dt`, folded into the existing `diag_extra`
  accumulator. Cortex + nucleus untouched (the pre-existing scalar-γ nucleus mismatch is **left for PI** — see flags).
- **Nucleus offset** — the single moved index: `nucleus_shell_kernel` base `wp.int32(Nc)` → `wp.int32(Ne)` (3 call sites).
- **Viewer** — `ff_crawl_viewer.py` splits cortex / MT-arms / MTOC / nucleus by the `Ne`/`n_mt`/`mtoc_idx` npz keys and
  draws the aster on-top (through the cortex cage); `ff_viewer_html.py` line layers gained `opacity`/`on_top`.

## Gates (all met)

**(a) Bit-identity — `--microtubules` OFF ⇒ byte-identical.**
- CPU exact: `np.array_equal(frames, golden)` PASS on **both** explicit and host-implicit, `max|Δ|=0.00e+00`,
  metrics equal (CPU verified deterministic run-to-run). Ne==Nc when off.
- GPU: new-OFF vs HEAD-baseline-OFF `max|Δ|=4.77e-7` = exactly the GPU run-to-run atomic-add noise floor (rel 6.4e-8)
  → the edit adds no logic to the OFF path.

**(b) Merged-K buckling + decoupling.** Single MT fiber's discrete buckling eigenvalue in the **merged** K:
**8.29 pN vs Euler π²EI/L² = 7.90 pN (5.0% at seg 0.25 µm), converging 5.0% → 2.5% → 1.3%** — matches the
standalone `microtubule.py` sanity numbers exactly. MT and cortex bending sub-blocks are **bit-identical**
merged-vs-standalone (true decoupling — the concat does not corrupt per-fiber bending).

**(c) V/V0 invariant + CG finite — native scale on the A5000.**
- V/V0 = **1.0004** at **native N=486,481** (Nc=483,000 + 40-tube aster + 3,000-bead nucleus), 40 implicit steps,
  no divergence, no OOM on the 16 GB A5000. Cortex trajectory unperturbed by the load-free aster (as designed —
  MT contributes bending harmlessly with no tip contact; tip↔cortex contact is a deferred post-Stage-1 step).
- CG iters **finite/bounded**: native OFF 113 → ON 178 (converged, well under the 400 maxiter). Merging the aster
  costs a bounded **~1.5×** — the plan's anticipated "only numerical cost of the merge".

### R-CG: Jacobi preconditioner TESTED and REJECTED

The plan (R-CG guardrail) suggested a Jacobi preconditioner "if κ_MT inflates CG iters". Measured on the A5000, a
diagonal (Jacobi) preconditioner **inflates** iterations instead (mid 115→279 = 4.4×; native 171→291 = 2.6×). Cause:
K is dominated by the rank-1 crosslink blocks `k·ûûᵀ` (k_xl≈4.6e5) which are **not diagonally dominant**, so `1/diag`
mis-scales them. Un-preconditioned CG is therefore retained (bounded ~1.5×). `ff_implicit_step_gpu` is unchanged
except a docstring recording this. A **block-Jacobi (3×3) / IC(0)** preconditioner is the flagged follow-up **iff**
native CG throughput becomes a production bottleneck — out of Stage-1 scope.

## Viewer (browser-verified)

`ffn_sim/outputs/ff/figs/mt_stage1_morph.html` (+ `mt_stage1_aster.png`). Verified in **real headless Chrome**
(puppeteer-core + system Chrome, WebGL 2.0 SwiftShader): no JS/page errors, GL context alive (not lost), scene
renders (28.9% non-background), and the **40 orange aster arms genuinely radiate** from the MTOC out to the woven
cortex cage (2449 scene-orange px spread 222×225, 1698 beyond the MTOC block) — not a grep/byte check.

## Surface-to-PI flags (no magic numbers were tuned to pass a gate)

1. **`k_hub`** (centrosome/PCM rigidity) — interim = cortex `xl_k.max()`; needs a KB anchor before it is "physiological".
2. **`n_mt`, `L_mt_um`** (aster size/reach) — MCF7 aster geometry must be KB-grounded, not chosen to reach the cortex.
3. **Nucleus scalar-γ drag mismatch** — pre-existing (nucleus beads run at `gamma_rep`, not their Stokes drag);
   correcting it perturbs the validated baseline, so it is deferred as a **separate** PI decision (not folded here).
4. **Native CG cost** — un-preconditioned merge is ~1.5× CG iters; block-Jacobi/IC(0) is the follow-up if it matters.

## Viewer + native-crawl note (PI feedback 2026-07-07)

PI reviewed the first (toy `cortex_fil=150`) viewer and flagged it as "not native, not crawling, looks like it's
splitting." Findings after investigation:
- **"Splitting" = sparsity artifact, NOT physics.** A 150-filament cortex is under-crosslinked/floppy and deforms
  irregularly under the crawl forces (end radius std 1.88, max 1.7×R₀), and the fixed-topology hull mesh distorts.
  At **native (Nc=483k)** the cell is a near-perfect coherent sphere (radius std **0.07**, max ~R₀, V/V0=1.0001).
- **Native full-resolution viewer** rendered (all 69k fibers, no downsampling — PI standing rule) and browser-verified.
- **Does NOT translocate.** Over 60 s (6000 steps, 3 full clutch-treadmill/KMC cycles) at native scale, disp∥=0,
  COM flat, `bound=1.00` throughout. The leading-edge protrusion is INTERNAL (net-zero reaction spread over all
  cortex nodes) and acts on the upper-front cap, while the adhesion clutches sit on the *basal* cap and stay
  uniformly bound to fixed substrate anchors → the adhesion footprint never advances → no net motion. This is the
  **pre-existing "traction piece" (piece-3) gap** — the crawl needs a polarized adhesion mechanism (front nascent
  adhesion ahead of the footprint + rear release) to translocate. It is **independent of this Stage-1 MT merge**
  (MT is load-free, tip↔cortex contact deferred) and is a physics decision surfaced to PI, not a parameter tweak.
- Viewer downsampling **removed by default** (`max_fibers=0` = show all); committed figure regenerated at a dense
  coherent scale (Nc=14000, all fibers). The full 483k-node HTML (452 MB) is regenerable, not committed (size).

## Next (plan §2)

Stage 2 (plasma membrane as a separate explicit array — area-tension + ERM + containment via `diag_extra`) and
Stage 3 (nuclear-envelope surface), promoting the verified separate-array pattern from `scripts/ff_membrane_cell.py`.
