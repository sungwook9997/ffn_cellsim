# Phase C migration · M2 — per-cell lamellipodium → Warp loop (BOOT + PORT SPEC)

**This is the boot artifact for a FRESH session to do M2** (the de-cohesion lamellipodium
port). Context got long in the session that produced this; M2 is a physics-subtle port best
done in a clean context. Read this top-to-bottom, then implement §M2b.

---

## BOOT STATE (2026-06-21) — where Phase C migration stands

**Goal:** run the de-cohesion verdict on the fast Warp DCM engine (not the old per-step
`cpu_local_snapshot` HOOMD driver). The de-cohesion PLATEAU verdict was RETRACTED
(`DCM_LANDMINE_REGISTER_2026-06-21.md`, branch `dcm/decohesion`): the committed A/A0=1.958/2.286
was driven by the deleted `SettlingForce` body-force proxy; the only mechanistic-isolating
control (ctrl_C, lamel-only) reached A/A0 **3.0 above the claimed ceiling** but diverged past 120k
→ verdict is **UNDECIDED**. We settle it cleanly on the migrated engine.

**Done (all committed on `h7/compartment-platform`):**
- Engine adopted (`ecf130d`), G1 GPU-backend parity 15/15 (`13435b6`), G2 substrate forces
  ported (`d8ab9c7`), G3/4/5 engine API+docs+CI (`f055108`), STRUCTURE/CLAUDE diffs (`e160b08`).
- **M1** (`5e20c54`): `ffn_sim/warp_port/dcm_warp_decohesion.py` — substrate (z-well + wetting)
  wired into the multicell Warp loop + cleanball-on-substrate build + per-frame metrics
  (top-down A/A0 + maxZ + V/V0 + COM drift) + divergence guard.
- **M1.5** (`5b61894`): hash-grid path (cohesion+contact O(N·k)). **N=100 wetting run on A5000:
  A/A0 1.0→1.45 (genuine, maxZ held 81µm, V/V0 1.04), ~3 min @ 323 steps/s.**
  Artifact `outputs/warp_decohesion/n100_wetting_a5000.json`.

**Established conventions (USE THESE — already validated in `dcm_warp_decohesion.py`):**
- **Physiological per-node drag** `γ = 6π·η·R/nv`, η=65.9 Pa·s (NOT legacy ResolvedDCM 3.9e-10 —
  that blows up). With it the Warp loop is stable at **dt=8e-6** (same as the HOOMD driver) ⇒ the
  30× per-step win is a true wall-clock win.
- cleanball-on-substrate build (`build_multicell` gap=2.2, rested so min-z = z0).
- warmup soft-start (1000 steps at 0.1·dt). hash-grid every step. divergence guard truncates.
- params: k_vol 7.73e5, rep 2e8, adh 1e7, w_cs_jm2 2.85e-3, force_cap 5e-8.

**M2 = add the per-cell lamellipodium to `dcm_warp_decohesion.py`'s loop.** Then M3 (junction-switch,
small) + the decisive N=100 full-stack run. The N=100 wetting A/A0 1.45 is the wetting-only
contribution; the lamellipodium crawl (ctrl_C ~3.0) is additive — that is what M2 adds.

**Boot checklist for the fresh session:**
1. `git -C /Users/sw1/ffn_cellsim-platform log --oneline -8` (confirm the commits above).
2. Read `ffn_sim/warp_port/dcm_warp_decohesion.py` (the driver to extend) + `ffn_sim/warp_port/ENGINE.md`.
3. Read this spec §A–§G + §M2b.
4. The HOOMD reference is on branch `dcm/decohesion` (`git show dcm/decohesion:ffn_sim/cell/dcm_lamellipodium_gpu.py`).
   The existing Warp tether kernel is `ffn_sim/warp_port/dcm_neighbor_warp.py::lamellipodium_tether_accum`
   + `gather_lead_pos` (single-cell; see §B/§E caveats for the multi-cell deltas).
5. gbook A5000: rsync `ffn_sim` → `gbook:/home/sungwook/ffn_phase_c/` then
   `cd ~/ffn_phase_c && PYTHONPATH=. python -m ffn_sim.warp_port.dcm_warp_decohesion --device cuda:0 ...`
   (Syncthing to this worktree is STALE — use rsync. `source ~/miniconda3/etc/profile.d/conda.sh;
   conda activate ffn_sim; export CUDA_HOME=$CONDA_PREFIX`.)

---

## M2b — implementation plan (the leanest faithful port)

Per §F: **actin = host-managed advancing anchor points (clutch ~rigid), NOT BAOAB DOFs.** Justified:
`k_clutch 5e-2 ≫ k_tether 4e-3`, tether capped 5e-9 ⇒ actin slip ≤ tether_cap/k_clutch = 0.1µm = 0.2·ℓ₀.

1. **Host, low cadence (every `batch_steps=50` steps):**
   - rim cells: one-shot at build, `centroid_z − z0 ≤ rim_contact_band·R` (§A).
   - per rim cell: cell centroid xy (segment-mean over cof) + spheroid-centroid xy (mean of rim-cell
     centroids) → outward unit dir `out = norm(cell_cen − sph_cen)` (xy, z=0); leading basal nodes
     (§B predicate: in_rim & basal-band & `proj ≥ lead_frac·|rel_xy|`).
   - ratchet SEED/ADVANCE (§D): per leading node, find cell's outermost front bead within
     tether_radius; SEED if none (one seed_offset·ℓ₀ outward); else ADVANCE one ℓ₀ outward when
     `(fb_proj − node_proj) < catch_factor·ℓ₀` AND `rng.uniform() < p_advance`, with
     `p_advance = min(1, v_front·(50·dt)·S_kinetic/ℓ₀)` (computed once). Seed RNG with 7 for parity.
   - rewrite the device `actin` anchor array + the per-leading-node geometry arrays
     (lead_idx, ccx, ccy, ox, oy, proj).
2. **Device, every step:** `gather_lead_pos` → `lamellipodium_tether_accum` into `force_d`.
   **Two fixes vs the existing kernel (§E):** (a) add a **same-cell mask** (pass per-actin cell id +
   per-leading-node cell id, or pass cell-local actin slices); (b) **decide Newton-3 reaction** — drop
   it (formalises the rigid-pin; recommended for the lean port) OR add `atomic_add(force, a0+aj, −f)`
   if integrating actin. Document the choice.
3. **Flags on `run_decohesion`:** `lamellipodium: bool`, plus the §G params (defaults from the table).
   Wire CLI `--lamellipodium`.
4. **Validate:** CPU smoke (small N, rim cells crawl: footprint grows at held maxZ, finite, V/V0~1);
   then gbook A5000 N=100. Co-track A/A0 + maxZ + V/V0 + COM drift (NEVER A/A0 alone — the landmine
   rule). Watch for divergence (ctrl_C diverged past 120k on the old engine; the guard truncates).
5. **M3 (after M2):** junction-switch = host cad_mult reduction under crowding at cadence 500
   (`attach_junction_switch` in `dcm_gpu_build.py` on `dcm/decohesion`; secondary effect — settle
   dominated, per the register). Then the **decisive run**: N=100 cleanball, lamellipodium +
   junction + wetting, NO proxy, headline A/A0 + maxZ + V/V0 + crawl_residual + a figure
   (visualize-at-closeout rule), commit a results_manifest claim, Notion wrap.

---

## PORTING SPEC (from a focused read of the `dcm/decohesion` reference)

> Quotes are verbatim from `dcm_lamellipodium_gpu.py` / `dcm_gpu_build.py` on `dcm/decohesion`,
> and `dcm_neighbor_warp.py` / `lamellipodium_warp.py` / `dcm_warp_crawl.py` in the working tree.

### A. Rim-cell detection (`detect_rim_cells`)
Criterion (`dcm_lamellipodium_gpu.py:272`): `basal = np.flatnonzero((centers[:, 2] - z0) <= contact_band * R)`
— a cell crawls iff **centroid z within `contact_band·R` of z0**. Default `rim_contact_band = 1.5`.
Fallback if empty (`:274`): lowest third by z — `np.argsort(centers[:,2])[: max(1, n//3)]`.
**One-shot at build** (`dcm_gpu_build.py:415`), frozen, passed to tether + ratchet; never recomputed.
`neighbour_factor`/`max_neighbours` are accepted but **unused** — omit.

### B. Per-cell leading-node selection
Predicate (`:427–442`): leading iff `in_rim & basal & (proj >= lead_frac·|rel_xy|)` where
`basal = |z_node − z_basal| ≤ basal_band`, `proj = (node_xy − cell_centroid_xy)·out_hat`.
**Two centroids:** per-cell centroid `(cen_x[c],cen_y[c])` (segment xy-mean over `cof`) for `rel`/`proj`;
**spheroid** centroid `(sx,sy)` = mean of **rim-cell** centroids (`:421–424`) for the outward dir
`out_hat = normalize(cell_centroid_xy − spheroid_centroid_xy)`, **z dropped** (2D radial; updater uses
`[cc-sph, 0.0]`, `:575`). `lead_frac=0.0` ⇒ outward hemisphere; `basal_band = basal_band_factor·R`,
`basal_band_factor=0.3`. Per-node arrays to device: lead_rp (gathered/step), lead_ccx/ccy (cell
centroid), lead_ox/oy (cell outward), lead_proj (node `npj`), lead_idx (global row).
⚠️ The single-cell helpers (`_leading_geometry`, `lamellipodium_warp`) use one global centroid +
radial-distance selection — **NOT** the multi-cell two-centroid predicate. Port §B fresh.

### C. Actin anchor pool + clutch
Pool (`:202–246`): type `actin_lamel` (typeid 2), parked far off-cluster (force-free), `cof=−1`,
contiguous trailing tag block `[a0, a0+n_pool)`, `n_actin = n_rim·pool_per_cell`, `pool_per_cell=60`.
Activated **basal** at `z_basal = z0 + seed_basal_offset` (`seed_basal_offset=0.3e-6`).
Clutch (`ActinClutchAnchorGPU`, `:281–329`): `F = −k_clutch·(r − anchor)` capped `clutch_cap`;
`k_clutch=5e-2`, `clutch_cap=1e-7`. **Actin can be host-pinned anchors** (clutch ≫ tether ⇒ slip
≤0.2·ℓ₀; the standalone Warp tether already does this — `lamellipodium_warp.py:8-9` "the
clutch-gripped actin is the fixed anchor").

### D. Ratchet advance (`GpuLamellipodiumAdvance`)
Cadence: every `batch_steps=50` steps. `p_advance = min(1.0, v_front·(batch_steps·dt)·S_kinetic/ℓ₀)`
(`dcm_gpu_build.py:594`), computed once. Per leading node (`:597–628`):
front bead = cell's active actin within tether_radius with **LARGEST outward projection** (NOT
strictly-outward-of-node — that re-seeds); **SEED** (`fb=None`) one `seed_offset·ℓ₀` outward of node
(`seed_offset=1.0`); **ADVANCE** else if `(fb_proj − node_proj) < catch_factor·ℓ₀` (`catch_factor=1.2`)
AND `rng.uniform() < p_advance` → new bead at `cell_cen_xy + (fb_proj + ℓ₀)·out_hat`, z=z_basal. Stop
on pool exhaustion. RNG `np.random.default_rng(7)`. Advance distance always ℓ₀; v_front·S only set the
per-tick probability.

### E. Traction tether
Force (`:479–488`): `F_node = +k_tether·(r_actin − r_node)` capped `tether_cap`, Newton-3 `−f` on the
actin bead. `k_tether=4e-3`, `tether_cap=5e-9`, `tether_radius=3e-6`. Actin selection in the FORCE:
same-cell, **strictly outward of node** (`act_proj>0`), **nearest** within radius.
Existing Warp `lamellipodium_tether_accum` (`dcm_neighbor_warp.py:162–204`): one thread/leading-node,
loops all actin, picks nearest with `dist≤radius and act_proj>npj` (uses cell centroid + lead_proj),
applies capped pull `atomic_add(force, lead_idx[t], +f)`. **Missing vs reference:** (1) no same-cell
mask (assumes the passed actin array is this cell's), (2) no Newton-3 `−f` on actin. `gather_lead_pos`
(`:153–159`) refreshes leading positions on-device each step.

### F. Faithful-vs-simplified
Module note (`:51–60`): explicit (not lumped) = actin DOFs + clutch spring + tether + ratchet;
simplified = "born gripped" (anchors at activation) + per-leading-node radial lanes (not free Arp2/3).
Recommendation: **host advancing anchors + device tether, actin pinned (no DOF integration)** — slip
≤0.2·ℓ₀; lost = bit-faithfulness w/ dcm/decohesion, clutch-as-measurable (slip/rupture), actin thermal
motion — all immaterial for a spreading verdict. Add the same-cell mask + decide Newton-3 (§E).

### G. Params (defaults `ResolvedGpuLamellipodium`)
v_front 6e-6/60 m/s (6 µm/min, lit 3–12) · ℓ₀ actin_rest_length 0.5e-6 m · S_kinetic 1.0 ·
catch_factor 1.2(×ℓ₀) · seed_offset 1.0(×ℓ₀) · rim_contact_band 1.5(×R) · lead_frac 0.0 ·
basal_band_factor 0.3(×R) · seed_basal_offset 0.3e-6 m · k_clutch 5e-2 N/m · clutch_cap 1e-7 N ·
k_tether 4e-3 N/m · tether_cap 5e-9 N (5·60Pa·area_node, Gil-Redondo 2023) · tether_radius 3e-6 m ·
pool_per_cell 60 · batch_steps 50 · seed 7. Derived: z_basal=z0+seed_basal_offset,
basal_band=0.3·R, p_advance=min(1,v_front·50·dt·S/ℓ₀), gamma[actin]=6πηR/nv (only if integrating actin).
CFL (`:134`): clutch dt<2γ/k → at γ=2.22e-4,k=5e-2 → dt<8.9e-3 (spread dt safe).
