# Thread-C unified-cell integration plan (2026-07-07)

Design workflow (3 agents + synthesis), code-grounded. Execute stage-by-stage; each stage bit-identical when its flag is OFF.

All three designs verified against the live code. Every load-bearing claim checks out, and the code resolves the three real disagreements between the designs. Here is the synthesized plan.

---

# Thread-C: Unified FF cell — one implementation plan

Node-adding compartments (MT aster, plasma membrane, nuclear envelope) into the production crawl driver `ffn_sim/scripts/ff_crawl_on_substrate.py`. Grounded in the live code: `implicit_ff.py` (scalar-γ solver + `diag_extra`), the crawl driver (`gamma_rep = median(gammas[:Nc])`, existing `diag` accumulation for clutch/substrate, degenerate-bond skip `L<1e-9`), and the proven separate-array pattern in `scripts/ff_membrane_cell.py`.

## Naming (fixed for the whole plan — resolves Design-1 `Ncx/Nfib` vs Design-3 `Nc/Ne`)

- **`Nc`** — cortex **surface** node count, block `[0, Nc)`. Owns everything geometric: ConvexHull `faces`, turgor/osmotic `F[:3Nc]`, membrane proxy, substrate, clutch, protrusion, gravity, centroid/mean-radius reductions.
- **`Ne = Nc + Nmt`** — **elastic** node count, block `[0, Ne)`. Owns bending, crosslinks, reshape, the implicit `K`, per-node γ.
- **`N = Ne + n_nuc`** — total, for `pos_all` and the solve vector.
- Node layout: `pos_all = [cortex(Nc) ; MT+MTOC(Nmt) ; nucleus_beads(n_nuc)]`.
- The single driver offset that moves: `nucleus_shell_kernel`'s base index `wp.int32(Nc)` → `wp.int32(Ne)`. Everything the code calls `Nc` today that means *cortex surface* **stays `Nc`**; only the nucleus offset and the K/γ/pos spans become `Ne`/`N`.

---

## (1) Final architecture — what MERGES, what steps SEPARATE, and why

| Compartment | Placement | Why |
|---|---|---|
| Cortex (bend + α-actinin xl + osmotic rank-1 `k_vol·g·gᵀ`) | **implicit K** (existing) | α-actinin k≈4.6e5 pN/µm sets the explicit CFL; hidden in K → dt is accuracy-bound |
| **MT aster** | **MERGE into the same implicit K** | hypothesis **CONFIRMED below** — MT arm is just a fiber with κ=20 |
| MTOC hub, MT node drag | **`K` via `diag_extra`** | reuse the existing diagonal-enrichment path (crosslink block for the hub; per-node γ correction for MT) |
| Clutch, substrate, **membrane containment** | **`K` via `diag_extra`** | stiff one-sided penalties on cortex DOFs; **Design-3 wins over Design-2 here** (see §3, R-CONT) |
| **Plasma membrane** (area-tension + ERM + reservoir) | **SEPARATE explicit array**, stepped AFTER the cortex solve | soft in the physiological plateau (γ_mem→dt_stab≈1.5 s ≫ cortex dt); proven in `ff_membrane_cell.py` |
| **Nuclear envelope** (bilinear lamina) | **SEPARATE explicit array**, augments then replaces the bead cloud | even softer (K_chrom≈6, K_lamin≈18 pN/µm) → explicit-stable at cortex dt |

**MT-merge hypothesis — CONFIRMED (from code).** `assemble_K_current_cupy(pos, bend_triples, alpha, xl_ij, k_xl, N, diag_extra)` builds `K_bend` purely from `bend_triples` + per-triple `alpha`, and `alpha = _per_triple_alpha(net)` reads **per-fiber κ**. So concatenating MT fibers into the cortex `FiberNetwork` with `kappa=KAPPA_MT` makes their triples land in the *same* `(γ/dt·I + K)Δx=F` solve at 300× the bending stiffness — unconditionally stable, **no second solver**. MT `κ/seg³ = 20/0.5³ = 160 pN/µm` ≪ the crosslink `k_max≈4.6e5` that already sets any explicit dt, and in K it is dt-irrelevant → **MT changes neither path's dt.** Axial inextensibility is the existing per-fiber `reshape_kernel` (dt-independent). MT requires the `--implicit` GPU path; the explicit path is dev-only.

**The invariant that shapes all coupling (bit-identity contract).** Every new compartment is purely additive and **default-OFF**, so all-flags-off ⇒ the kernel-launch sequence is byte-identical to today. This is why membrane/envelope are *separate arrays riding outside their host*, coupled only by soft ERM + one-sided containment: the cortex keeps its own validated osmotic envelope. We do **not** hand osmosis to the membrane the way the `ff_membrane_cell` *demo* does (that demo zeroed cortex turgor — it would perturb the validated path). The existing lumped-Laplace membrane proxy (`dP_mem_area` via `turgor_kernel`, driver L306/L326) is **retained when `--membrane` is off** (bit-identity) and **replaced** by the explicit surface when on.

---

## (2) Ordered build steps + incremental rollout + per-stage gate

### New functions

`ffn_sim/ff/fiber_network.py` (engine-agnostic, reusable):
```python
def concat_fiber_networks(nets: list[FiberNetwork]) -> tuple[FiberNetwork, np.ndarray]:
    """Stack pos/segments/bend_triples/fiber_offsets/seg_rest/kappa with index offsets.
    Returns (merged, node_off); sub-net k owns merged.pos[node_off[k]:node_off[k+1])."""
```

`ffn_sim/ff/microtubule.py` (MT-specific hub wiring):
```python
def merge_aster_into_cortex(cortex_net, aster: MicrotubuleAster | None, *, k_hub_pn_um: float) -> dict:
    """Returns {'net','Nc','Ne','mtoc_idx','hub_xl_i','hub_xl_j','hub_xl_k','hub_xl_rest','mt_gammas'}.
    aster is None → no-op passthrough: net=cortex_net, Ne==Nc, empty hub arrays  (guarantees MT-off bit-identity)."""
```
Also tweak `build_microtubule_aster`: run beads from `seg_um` out to `L` (**drop the coincident centre bead**) and return the MTOC as a separate node/index. Keep `L_mt_um` fixed so the verified buckling gate is unaffected.

### Rollout — 5 stages, each gated. Freeze a golden trajectory at Stage 0; every later OFF-path must match it.

**Stage 0 — baseline lock.** Current cortex + nucleus-bead crawl. Hash the golden `_on.npz` `frames` + `disp_along_um`/`traction_nN`/`v_crawl_nm_s`. This is the bit-identity oracle.

**Stage 1 — MT aster into implicit K** (lowest risk: no new array, no soft coupling).
- `build()`: after `build_crosslinked_cortex`, if `--microtubules`: `aster = build_microtubule_aster(centre=c, n_mt=…, L_mt_um=…)`; `m = merge_aster_into_cortex(cx.net, aster, k_hub_pn_um=cx.xl_k.max())`. Replace `cx.net` with `m['net']` for the elastic path; keep the cortex block `[:Nc]` for surface physics. Return `Nc, Ne, mtoc_idx` + hub-crosslink arrays.
- `run()`: `Nc=m['Nc']`, `Ne=m['Ne']`, `N=Ne+n_nuc`. Concatenate hub crosslinks into `xl_ij/kxl/r0` (host + `xl_d/kxl_d/r0_d` device + cupy `xlij_cp/kxl_cp`). `pos_all = concat([merged.pos, nuc_pos])`. `alpha/tri/foff/soff/gammas` all from the merged net (`_per_triple_alpha` auto-routes KAPPA_MT).
- **MTOC hub (Design-1 resolution — verified):** the code skips any bond with `L<1e-9` (`implicit_ff.py:65`) and `build_fiber_network` duplicates coincident bases into *separate* nodes — so a weld-at-centre hub is **silently dropped** and the aster becomes free rods. Instead add **one MTOC node** and bond it to each arm's innermost real bead (at `seg_um`) with a **finite-rest** stiff crosslink (`rest=seg_um`). Reuses `link_spring_kernel` (force) + the `k·ûûᵀ` block of K (stiffness) — no new kernel. `k_hub = cx.xl_k.max()` interim (flag to PI, §4).
- **MT node drag (Design-3 R7 resolution — verified real):** the solver uses a **scalar** `gamma_rep = median(gammas[:Nc])` (`driver:123`), so an added node with different drag gets the wrong timescale. Fix **only for MT nodes**, via the existing `diag` accumulator: `diag[3*mt] += (γ_mt − gamma_rep)/dt` (x/y/z). This gives MT correct fiber-log drag **and preserves cortex+nucleus bit-identity** (the pre-existing scalar-γ nucleus behavior is left untouched — flag that latent nucleus-drag mismatch to PI as a *separate* decision, do not fix it here).
- **Gate:** (a) `--microtubules` OFF ⇒ `np.array_equal(frames, golden)`. (b) merged-K single-fiber discrete buckling eigenvalue reproduces `F_crit=π²·EI/L²=euler_buckling_load(KAPPA_MT,L)` (7.90 pN at L=5 µm) within the same ~1–5% discretization error `microtubule.py` logs; cortex fibers' own bending eigenvalues unchanged (decoupling). (c) V/V0 unchanged; CG iters/step bounded (add Jacobi preconditioner to the cupy CG if κ_MT inflates them — the only numerical cost of the merge).

**Stage 2 — plasma membrane as separate explicit array + ERM + containment.**
- Arrays (build in `build()`, upload in `run()`, mirror `ff_membrane_cell.run` L86–100): `mem_mesh = build_membrane_mesh(R+gap, subdivisions=4)` (2562 nodes/5120 faces, `gap≈0.3 µm`); `mpos_d, mfaces_d, mf_d, marea_d[1]`; ERM map `erm_c_d = cKDTree(cortex_pos).query(mem_verts)[1]`, `erm_rest_d = ‖mvert−cortex‖` (force-free at formation), `bound_d = ones`. Scalars `k_erm=50`, `f_rupt=2π√(2κ(γ_mem+γ_MCA))≈5–40 pN`, `k_wall=5e3`, `m_gamma=median(gammas)`, `renv_d[1]`.
- Loop wiring per cortex step (`gpu_impl` branch):
  1. reduce mean membrane radius → `renv_d` (frozen for the whole cortex solve).
  2. **cortex implicit solve** — with ERM-reaction + containment folded in (next bullet). `mpos_d` frozen here.
  3. **membrane sub-step** (explicit, AFTER the cortex step), `n_sub` times (§3 CFL).
  4. **rupture eval once per step** (bleb emergence): evaluate ERM tension vs the real `f_rupt`, mutate `bound_d` — kept out of the RHS.
- **Coupling into the cortex force** (`gpu_force_fn`, under `if membrane:`, against the *frozen* `mpos_d`/`renv`):
  - **Containment → `diag_extra`, NOT explicit-in-F** (resolves Design-2 vs Design-3 in favor of **Design-3 R4**): verified `γ/k_wall = 45/5000 ≈ 9e-3 s`, and the implicit path runs at a much larger dt, so explicit containment is unstable. Add it to the same `diag` the driver already builds for clutch/substrate: for cortex nodes with `r>renv`, `diag[3*i (+1,+2)] += k_wall`. State-dependent (recompute the penetrating set each step) exactly like the bound-basal clutch set.
  - **ERM reaction** stays explicit-in-F (soft, k_erm=50, force-free at formation → keeping it out of K keeps M SPD/reusable): launch `erm_tether_kernel` writing the cortex half into `f_d[:Nc]`, membrane half into an ignored scratch, with `f_rupt=+inf` so no `bound` mutation inside the RHS.
- **Membrane sub-step** (its own overdamped Euler, per `ff_membrane_cell`): zero `mf_d` → area-reduce → `sigma = reservoir_tension(A,A0,mem)` → `membrane_area_kernel` (−σ∂A/∂x) → **self-Laplace** `turgor_kernel(mpos_d, centre, 2σ/R_mem·A/Nm)` (inflates the sheet to R_mem on its *own* nodes → no double-count of cortex turgor) → ERM membrane-half (`f_rupt=+inf`) → `_step_scalar_gamma(mpos_d, dt_sub, m_gamma)`.
- When `--membrane` on, **disable the lumped `dP_mem_area` proxy** (the explicit surface replaces it); when off, keep the proxy → bit-identity.
- **Gate:** OFF ⇒ bit-identical to Stage 1; Young-Laplace ΔP=2σ/R on the membrane; reservoir plateau→K_A upturn; ERM rupture f_t≈5–40 pN emergent; **no cortex node leaks past the envelope** (containment gate); no divergence when driven to lysis.

**Stage 3a — nuclear envelope AUGMENT (keep bead cloud, add lamina shell)** — Design-2 conservative staging.
- `ne = resolve_nuclear_envelope(R_nuc=0.70R, subdivisions=3)` (642 nodes/1280 faces). `K_chrom = E_nuc·h = 399·0.015 = 5.99 pN/µm`, `K_lamin = 3·K_chrom = 17.96`. Arrays `npos_d, nfaces_d, nf_d, narea_d[1], ruptured_d[1280]`, `n_gamma = 6πη·a` (bead-scale).
- Coupling (stepped after the cortex solve, alongside the membrane): envelope self-Laplace `turgor_kernel(npos_d, nuc_centre, 2σ_lam/R_nuc·A/Nne)` with `σ_lam = lamina_tension(A,A0,ne)`; **containment of chromatin beads** `membrane_containment_kernel(beads, nuc_centre, r_env(frozen), k_wall_nuc)`; add a **nucleus volume term + collapse guard** (Design-3 R5): reuse the cortex `k_vol·g·gᵀ` over envelope nodes and the existing `V ≤ 1.02·vmin` truncation as `V_nuc ≤ 1.02·vmin_nuc` (bilinear area tension alone has no volume preservation → would deflate under compression).
- **Gate:** OFF ⇒ bit-identical to Stage 2; Laplace ΔP=2σ/R_nuc; strain-stiffening knee ×ratio at ε=10%; rupture emerges >ε_rupture, absent below; V_nuc/V_nuc0≈1.

**Stage 3b — envelope REPLACE (optional, PI-gated).** Drop the bead cloud, resolve chromatin-only (`ratio_lamin=1.0`) so the lamina lives on the envelope, not doubled in the bead bilinear. Only after 3a is green.

**Stage 4 — unified full-cell, all additive units on.** Gate: whole-cell V/V0∈[0.8,1.25], basal-gap<0.4, bound_frac>0.5 (existing STABLE-ADHERED band); every compartment's own gate still holds in-situ; crawl speed in the retrograde-flow band (10–100 nm/s); the adversarial clutch-OFF audit still ≈0 net drift.

### Driver `Nc→Nc/Ne/N` mechanical edits (`ff_crawl_on_substrate.py`)
Replace `Nc` with `Ne` **only** in: the implicit solve vector `N`, `gammas`/`gamma_rep` span, `bt_cp/al_cp` dims, `reshape dim=merged_foff-1`, the `diag`/`vol_g` allocation `3*N`. Keep `Nc` in: `turgor_kernel`, `gravity_kernel`, `substrate_plane_kernel`, `spreading_*`, `leading_edge_*`, `sum_pos_kernel`, `sum_radius_kernel`, `cortex_volume_kernel`, `faces`, both osmotic `F[:3*Nc]`. `nucleus_shell_kernel` offset `wp.int32(Nc)`→`wp.int32(Ne)`. New argparse: `--microtubules --n-mt --membrane --envelope --gap`. Record `mpos_d/npos_d/bound_d/ruptured_d` into the frames npz for the HTML morphology viewer (cortex + translucent membrane sheet + nucleus envelope), per the viz rule.

---

## (3) Stability guardrails

| # | Risk | Number | Guardrail |
|---|---|---|---|
| R-MT | MT stiff bending destabilizes explicit coupling | α_MT=160 pN/µm ≪ k_xl=4.6e5 | Bending goes in K (unconditional); axial via reshape projection (dt-independent) |
| R-HUB | Aster arms fly apart (silent degenerate-bond skip) | `L<1e-9`→`continue`, verified `implicit_ff.py:65` | Dedicated MTOC node + **finite-rest** hub crosslinks; enters K_xl, no explicit stiff hub spring |
| R-γ (R7) | Scalar `gamma_rep` gives added nodes wrong drag | `gamma_rep=median(gammas[:Nc])`, verified `driver:123` | Per-node `diag[3i]+= (γ_i−γ_rep)/dt` for **MT only**; leaves cortex+nucleus bit-identity intact |
| R-CONT | Explicit containment unstable at implicit dt | `γ/k_wall=45/5000≈9e-3 s ≪ dt` | **k_wall → `diag_extra`** (the verified clutch/substrate pattern, `driver:363–366`), not explicit-in-F |
| R-MEMB | Membrane reservoir/lysis upturn (stiff) | plateau σ=10→dt_stab≈1.5 s; K_A=2.35e5→dt_stab≈6e-4 s | Adaptive `n_sub = ceil(dt/(0.5·15/σ))` from current σ, cap `n_sub_max=64`; reservoir **DEFAULT-OFF** for caveolae-deficient MCF7 ⇒ n_sub=1 in production; if a caveolae-competent type ever needs the upturn beyond the cap → give the membrane its own implicit area step (PI-gate), never sub-cycle thousands |
| R-NUC | Envelope deflates (area tension has no volume preservation) | — | Nucleus `k_vol·g·gᵀ` volume term + reuse `V_nuc≤1.02·vmin_nuc` collapse guard |
| R-RUPT | Envelope/ERM rupture transient | tension → 0 | Non-divergent by construction (rupture only *removes* force); evaluate at KMC cadence, out of the RHS |
| R-CG | κ_MT (275×) inflates CG iterations | — | Monitor cg_iters/step across the merge; add Jacobi preconditioner to cupy CG if it climbs |

**CFL / dt budget.** Implicit cortex+MT: no stiff term in the RHS → dt accuracy-bound (production). Membrane plateau (γ_mem): dt_stab≈15/σ≈1.5 s → n_sub=1. Nuclear envelope: K_lamin=18→dt_stab≈0.8 s, K_chrom→≈25 s → n_sub=1 unconditionally. **All soft surfaces are explicit-stable at production cortex dt with n_sub=1; the K_A reservoir upturn is the only stiff bound, and it is default-off.**

**Bit-identity gate (run before every commit; red ⇒ halt/surface).** (1) *Launch-count* (structural): with all new flags off, the set of launched kernels + dims equals the pre-change baseline exactly. (2) *Trajectory* (numerical): on the **CPU dev path** (deterministic reductions), `np.array_equal(frames_new, golden)` and `com` equal. GPU is not bit-reproducible run-to-run (atomic_add), so the GPU gate is launch-count + ≤1e-10 relative agreement; the exact bit gate lives on CPU.

---

## (4) Deferred, and why

- **MT-tip ↔ cortex contact — DEFERRED to a post-Stage-1 step (Design-1 resolution).** Not needed to bring MT into K: with no external load the stiff arms stay straight and stable, contributing bending harmlessly; the MT-off bit-identity + Euler gates pass without any tip force. When added (so the aster bears actomyosin load / positions the nucleus — the tensegrity point), it is a **soft explicit one-sided repulsion in `F_total`** on tip nodes penetrating the cortex surface (reaction transmitted through the arm to the MTOC and thence the nucleus), riding the `gpu_force_fn` pattern — **it does NOT go into K.**
- **Envelope REPLACE mode (Stage 3b) — DEFERRED behind AUGMENT (Stage 3a).** Augment (bead cloud + lamina shell) perturbs the validated nucleus least and keeps a clean bit-identity fallback; only collapse the bead cloud into the envelope after 3a is green and PI approves.
- **Fixing the pre-existing scalar-γ nucleus drag — DEFERRED (surface to PI).** The scalar `gamma_rep` already gives the nucleus beads cortex drag rather than their Stokes drag; correcting it changes the *validated* baseline trajectory, so it is a separate PI-gated decision, not folded into this merge.
- **MT dynamic instability (Mitchison-Kirschner)** — stays default-off; this is a **stable-aster baseline** pending KB DI rates (already noted in `microtubule.py`).

### Surface-to-PI flags (no magic numbers)
- `k_hub` (centrosome/PCM rigidity): interim = cortex `xl_k`; needs a KB anchor before it is called physiological.
- `n_mt`, `L_mt_um`: MCF7 aster size/reach must be KB-grounded, not chosen to reach the cortex.
- `RESERVOIR_STRAIN` flagged unsourced + MCF7 caveolae-deficient (`membrane_surface.py:34–39`) → reservoir off in production.

### Files carrying the change (absolute)
- `/Users/sw1/ffn_cellsim/ffn_sim/ff/fiber_network.py` — add `concat_fiber_networks`
- `/Users/sw1/ffn_cellsim/ffn_sim/ff/microtubule.py` — add `merge_aster_into_cortex`; drop coincident centre bead + expose MTOC node
- `/Users/sw1/ffn_cellsim/ffn_sim/scripts/ff_crawl_on_substrate.py` — `Nc→Nc/Ne/N` split, hub-crosslink concat, membrane/envelope arrays + loop wiring, `diag_extra` for MT-γ + containment, new flags
- `/Users/sw1/ffn_cellsim/ffn_sim/ff/{membrane_surface,nucleus_envelope}.py` — reused as-is (area machinery, tensions, rupture)
- **No change to `/Users/sw1/ffn_cellsim/ffn_sim/ff/implicit_ff.py`** — MT bending rides per-fiber κ; hub, MT-γ, and containment all ride the existing `diag_extra`.