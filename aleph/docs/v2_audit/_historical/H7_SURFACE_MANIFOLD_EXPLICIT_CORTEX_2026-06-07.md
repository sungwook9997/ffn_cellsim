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

# H.7 — Surface-Manifold-as-Spatial-Substrate for an EXPLICIT Cortex
**`aleph/docs/v2_audit/H7_SURFACE_MANIFOLD_EXPLICIT_CORTEX_2026-06-07.md`**
Status: DESIGN (pre-implementation). Branch `h7/full-cell-integration`. Synthesized from 4 parallel investigations (literature / model / systems / critique). Companion to `H7_CORTEX_AS_MESH_2026-06-07.md` (which rejected mesh-edges-as-cortex), `H7_DIRECTION_FOR_CODEX_2026-06-07.md` and `H7_SIGMA_A_NETWORK_INVESTIGATION_2026-06-07.md` (the γ-floor diagnosis).

> **One-line honest verdict.** A deformable finite-thickness cortical *manifold* used **only** as a geometric/locality/soft-confinement substrate — with actin, motors, crosslinkers, adhesions staying explicit — is **fine-grained-rule-compliant** and a genuine win for **spatiality, search efficiency, connectivity organization, and spreading coordination**. It is **ORTHOGONAL to the active-γ floor** (Gate A: heads don't walk, s_grip≈0; Gate B: M-SHAKE forbids buckling). It MUST NOT be sold, scoped, or tuned as a γ-floor fix.

---

## 1. Revised concept

The cortex is represented in **two decoupled layers**:

- **Layer 1 — explicit mechanistic objects (unchanged fidelity).** Actin filaments stay explicit bead-spring/WLC chains with stretch, bend, compression, **buckling**, and turnover. Myosin minifilaments (Stam-Hocky bipolar), crosslinkers (catch/slip), Arp2/3 branches, and FA clutches remain explicit particles/bonds. These generate **all** force. γ is measured **only** from their realized bond tensions.
- **Layer 2 — a deformable triangulated MANIFOLD as the cortical *spatial substrate*.** A finite-thickness (`h_cortex ≈ 200 nm`) deforming surface that supplies: (a) the 3D cell-surface shape and a shared coordinate system for membrane/turgor/FA/lamellipodium; (b) per-query local tangent/normal frames `(n̂, e1, e2)`; (c) signed normal distance `d_normal` and band center `d0` for soft confinement; (d) patch identity for neighbor-search acceleration. The manifold is **slaved to** the explicit bead cloud (beads → manifold geometry) and couples **back** to beads through **exactly one** force: a soft normal confinement `U_conf = ½ k_conf (d_normal − d0)²`.

This is explicitly **NOT** the rejected triangulated-surface-as-cortex idea. Mesh edges are **not** actin. The manifold carries **no** mechanical tension and **no** elastic DOF that enters the force budget. It is geometry + bookkeeping + one localization spring.

**Slaving direction (the fidelity-preserving invariant):** beads → manifold (re-fit/advect every K steps); manifold → beads (soft normal confinement only). The manifold never relaxes toward its own elastic minimum and pulls the cortex.

---

## 2. Motivation

Three concrete deficiencies of the current pure-3D `connected_mesh` build motivate a manifold:

1. **Search cost.** Every binder (`crosslinkers.py:1060-1064` `XlinkBondUpdater.act`, `myosin.py:1099` legacy + `myosin.py:1193` segment-projection, plus one-shot seeding `crosslinkers.py:613/647` and Arp2/3 `cortex.py:1332`) rebuilds a **global `scipy.spatial.cKDTree`** over all ~7000 cortex beads **every batch tick** and throws it away. Patch-local geodesic k-ring search replaces O(N log N) per query with O(k) — an analytic ~50–250× search-cost win (to be **measured** by the prototype, not asserted).
2. **No deformable surface / shared coordinate frame.** Today `cortex.py:1246 _project_to_shell_band` and `cortex.py:569 _tangent_plane_basis` give a **static, sphere-only** band + frames; `enclosed_volume.py` and `membrane_surface.py` assume a **centroid-sphere** for normals. A spreading / non-spherical cell needs a deformable surface that membrane, turgor, FA, and lamellipodium can all reference for "where is the surface here."
3. **Connectivity organization.** The manifold gives a principled locality structure for the cortical network on a curved surface, and local frames that make the curved-surface γ measurement (§9) shape-agnostic.

**What it does NOT motivate:** the active-γ floor. Per all four H7 docs the floor is **generation-limited** (Gate A s_grip walking + Gate B M-SHAKE buckling), upstream of both topology and measurement. A coordinate manifold adds zero myosin walking and zero bond extension. See §13 / honest crux.

---

## 3. FORBIDDEN vs ALLOWED mesh usage (the explicit guardrail table)

| # | ALLOWED (geometry / bookkeeping / soft confinement) | FORBIDDEN (mesh carries mechanics) |
|---|---|---|
| 1 | Manifold defines **rest shape & current deformed shape** of the cell surface | Mesh **edges represent actin filaments** (edge = filament) |
| 2 | Per-query **local tangent/normal frame** `(n̂, e1, e2)` | **Edge-length springs** on the manifold storing/transmitting force |
| 3 | **Signed normal distance** `d_normal`, band center `d0`, finite thickness `h_cortex` | **Area / dilation springs** (esp. tuned to produce γ) |
| 4 | **Patch identity** for neighbor-search acceleration (binding-candidate lookup) | **Bending/curvature elastic DOF** that enters the force budget |
| 5 | **Shared coordinate system** for membrane / turgor / FA / lamellipodium | **Motors/crosslinkers replaced by a surface active-stress field** |
| 6 | **ONE soft normal confinement** `U_conf = ½ k_conf (d_normal − d0)²` (keeps beads in shell, allows tangential motion + buckling + limited normal fluctuation) | **Hard projection** onto the surface (snaps beads to 2D, kills buckling) |
| 7 | Manifold **slaved** to explicit beads (re-fit/advected) | Manifold **independently relaxed** to its own minimum, pulling beads |
| 8 | γ **measured** from explicit filament/motor/xlink bond stress | γ **defined as a mesh property** (edge tension, area stress) |
| 9 | Manifold **supplies normals** that `enclosed_volume`/`membrane_surface` read | Manifold adds a **THIRD** Laplace/normal force on the same beads |
| 10 | Mesh resolution chosen so **patch_size ≥ physical_reach** (a correctness bound) | **Mesh resolution tuned to pass the γ gate** (or to move connectivity) |

**The operational line (single test):** *does any manifold DOF appear in the γ budget?* γ is `Σ` over **explicit** bond tensions `T = k(L−r0)` (soft MOP, `cortical_tension.py:289-393`), `λ·r0/Δt` (M-SHAKE rigid, 396-461), and the IK virial `Σ T·L·(1−(r̂·û)²)/(8πR²)` (`cortical_tension.py:380-386`). If a manifold edge/area/curvature term ever enters that sum → forbidden, halt.

---

## 4. Proposed DATA MODEL — **recommend (a) global-Cartesian + soft confinement** (reject barycentric)

**Recommendation: GLOBAL Cartesian beads + a SOFT manifold-confinement force.** Reject local triangle-id/barycentric/normal-offset state.

**The HOOMD reason.** HOOMD integrates **global Cartesian** positions. The constrained L-M BAOAB action (`constrained_baoab.py:827-955`) reads `pos + net_force`, adds the Fixman correction, predicts, and M-SHAKE-projects. A manifold confinement is then **just another `md.force.Custom` landing in `net_force`** — exactly like `erm.py` (radial pin) and `enclosed_volume.py` (Laplace normal pressure) — so it requires **zero integrator change** and is GPU-port compatible.

A barycentric/triangle-local representation would require a per-step **local→global before** and **global→local after** every HOOMD integration step — a CPU projection that (a) fights the Cartesian core, (b) **re-adds the per-step host GPU↔CPU sync the GPU-main port is removing** (`constrained_baoab.py:789-825`), and (c) complicates buckling (in-plane buckling is awkward in surface coordinates). Verdict: integrated state stays global Cartesian; the manifold is an **auxiliary field** queried for `d_normal`, the local frame, and patch id.

**Manifold-side structures (auxiliary, not integrated state):**
- `verts (n_vert, 3)`, `tris (n_tri, 3)` — triangulated icosphere/Fibonacci-sphere at resolution `n_tri ≈ n_filaments` (~1 filament-COM per patch).
- `tri_adj (n_tri, 3)` — edge-neighbor adjacency for geodesic BFS / k-ring.
- `bead_tri (n_cortex_actin,)` — bead→triangle map (nearest-centroid or barycentric), refreshed **lazily** every K batches (beads move ≪ patch size between rebuilds).
- per-face local frame `(n̂, e1, e2)` reusing `cortex.py:569-587 _tangent_plane_basis`.

**Tag-append invariant (hard migration constraint, `cell.py:914-993`):** the manifold must be **additive / appended last** and must **not** renumber cortex actin tags `[0, n_cortex_actin)`; every subsystem (xlink, myosin, lamellipodium, FA, nucleus) is appended last so absolute-tag Updater bookkeeping stays valid.

---

## 5. FORCE MODEL

- **Filament / motor / crosslinker / FA forces: UNCHANGED.** Global Cartesian, into `net_force`. No edits to AFINES filament mechanics, Stam-Hocky myosin, catch/slip crosslinkers, or FA clutches.
- **ADD exactly one manifold→bead coupling:** `F_conf = −k_conf (d_normal − d0) · n̂_local`, where `n̂_local` is the manifold normal at the bead's nearest surface point and `d_normal` is the signed distance. This is the **only** new force. It is **normal-only** (no tangential / edge-length spring on actin beads).
- **Manifold deformation is DRIVEN by explicit forces (two-way, Newton's 3rd law):** the confinement reaction `+F_conf` acts on the manifold node, and the manifold's shape is updated as a passive geometric tracker slaved to the bead cloud. The manifold has **no independent elastic dynamics** that feed the force budget. A fixed sphere is a **non-production dev fallback** only.
- **Laplace / surface ownership stays with the existing modules.** `enclosed_volume.py EnclosedVolumePressure` (turgor `F_i = ΔP·A_i·n̂_i`, lines 466-496) and `membrane_surface.py MembraneSurfaceTension` (membrane Laplace `F_i = −ΔP·(S/N)·n̂_i`, lines 34-47) remain the **sole** owners of the normal/Laplace pressure. They now read their **normals from the manifold** instead of a centroid-sphere assumption. The manifold adds **no** Laplace force itself (eliminates the triple-count, §9).

---

## 6. CONFINEMENT MODEL — soft `U_conf`, buckling-preserving

**Potential:** `U_conf = ½ k_conf (d_normal − d0)²`, `d_normal` = signed distance to the manifold's nearest point, `d0` = shell mid-plane, normal half-width from `h_cortex/2` (`h_cortex = 200 nm`, KU-3.17).

**`k_conf` anchor — `[MAGIC-NUMBER BLOCK REQUIRED]`.** Anchor to the **ERM tether** `k_ERM = 0.1 N/m` (KU-3.18, `erm.py:7-18,115`): the ERM/membrane tether **is physically** the cortex-to-shell confinement spring, and `U_conf` is literally `erm.py`'s `U = ½ k_ERM(|r|−R_cell)²` generalized to a **deformable, banded** target. So `k_conf` is the **same knob class** as `k_ERM` and must be anchored, not freely tuned. With `k_ERM = 0.1 N/m`, `σ_n = √(kT/k_conf) = 0.65 nm` (a thin pin = today's `erm.py`). For band-filling normal undulation, derive a **softer** `k_conf` as a *stated fraction of `h_cortex/2`* and **surface to PI** (sign-off). Cross-checks: `κ_m = 1e-19 J` (Rawicz 2000, KU-3.B1) and `γ_MCA = 1e-5 J/m²` (KU-3.B1).

**Why SOFT, not hard (buckling-preserving).** A normal-only penalty leaves **tangential motion free**, so filaments keep in-plane **buckling** and condensation (`r/r0 < 1`, the Gate-B ~10× amplification lever) plus limited normal undulation. A **hard projection** (the Janzen-Sknepnek 2026 mode) snaps beads to a 2D surface and **suppresses the very buckling** the team rejected the triangulated mesh for. The target is a **thin 3D shell, not a 2D sheet.**

**Honest caveat (Gate B).** Soft confinement makes buckling *possible* but does **not by itself unlock it**: the rigid backbone is currently held by **M-SHAKE** (`constrained_baoab.py`, `r/r0 = 1.000`). The actual lever is **relaxing the M-SHAKE backbone constraint** (a FROZEN-integrator change → **PI sign-off**); the existing independent filament lines could buckle the instant M-SHAKE is relaxed, **with no manifold at all**. The manifold is at most a convenient host for a buckling-capable cortex, not the cause of buckling.

**CFL / escape guard.** Attach gate `τ_conf = γ_b / k_conf`, `dt ≤ safety·τ_conf` (mirror `attach_erm` `erm.py:282-298` and `attach_enclosed_volume:559-575`); optional steeper half-harmonic/WCA wall **only** beyond `±h_cortex/2`; per-step inside-band/finite assertion (sibling of `constrained_baoab.py:849-856`).

---

## 7. NEIGHBOR-SEARCH MODEL (patch / geodesic)

**Current state (to replace):** there is **no** triangle adjacency or bead→triangle map anywhere. The only spatial metadata is the flat per-bead `filament_idx` (`cortex.py:1078`) and a **1D** per-bead→segment chain adjacency (`myosin.py:857-863 self._bead_to_segs`) — neither is a 2D patch map. Every binder rebuilds a global cKDTree each batch tick.

**Proposed:** candidate search = "beads in `bead_tri[t]` ∪ its **k-ring** neighbor triangles" via **BFS over `tri_adj`** → O(k) instead of O(N log N) per query, with no per-tick tree teardown/rebuild. Keep the index **topology-agnostic** so γ and myosin need zero edits.

**Correctness bound (NOT a tuning knob):** `patch_size ≥ physical_reach` must hold at construction so a single k-ring always covers the reach ball at every resolution. At fine `n_tri` with the `√(A/n)` compensation reach (~1 µm) a triangle can be smaller than the reach → **k must scale with resolution** (k-ring grows). Enforce this as a construction-time assertion; it is the guard that keeps candidate sets resolution-invariant (§10).

**Geodesic vs chord caveat:** a 2D UV/geodesic binding metric could marginally change which partners are "within reach" vs the 3D-chord metric on a curved shell (chord < arc). Judged second-order at `R_cell ≈ 7.5 µm` with 60 nm reach; quantify in the prototype, do not assume.

---

## 8. MOTOR/CROSSLINK BINDING MODEL — patch-local, and **honest about the `√(A/n)` hack**

Binding becomes **patch-local at the physical reach**: candidate partners are filtered by the **physical** crosslinker reach (`max_bind_dist = 60 nm`, `phase1_h3.yaml:124,404`), restricted to the bead's patch ∪ k-ring.

**Does this remove the `√(A/n)` bridge-reach hack? NO — and this is the load-bearing honest answer.** The manifold changes only the **search COST**, not the **physical bead spacing**. The ×40 sparsity is geometric: 1000 filaments on a `4πR²` shell sit `~√(A/n) ≈ 0.9–1.9 µm` apart **regardless of the spatial index**. The 2D test (`H7_SIGMA_A_NETWORK_INVESTIGATION:71-78`, `scripts/crosslink_2d_test.py`) is decisive: **60 nm dynamic reach → z = 0.20, giant 1.5% (FRAGMENTED); `√(A/n)` reach → z = 9.4, giant 99% (PERCOLATED).** A manifold/k-ring query at 60 nm returns the **same empty candidate sets** the current 60 nm cKDTree returns — just faster. Running binding in 2D manifold/UV coordinates only swaps the 3D-chord metric for the in-surface geodesic metric; the in-shell nearest-neighbor distance is **still `~√(A/n) ≫ 60 nm`** at `n_filaments = 1000`.

**Conclusion:** **keep `√(A/n)` as the documented, sanctioned coarse-graining-COMPENSATION device** (`crosslinkers.py:507-509`, `connected_mesh.py:57-62` — the crosslinker analogue of the mesoscale myosin force-scaling, the geometric dual of the ×40 areal coarse-graining, **not** a tuned `bind_scale`). The manifold **reorganizes and accelerates** the sparsity; it does not eliminate it. The only ways to remove `√(A/n)` are to raise `n_filaments` toward native (defeats the one sanctioned ×40 coarse-graining) — so we don't.

> For the **connectivity** sub-goal specifically, a **Mikado-on-sphere** weave (`cortex/mikado_sphere.py` reusing `ecm/mikado.py`, per `H7_CORTEX_AS_MESH:65-68`) is a better-scoped tool than a manifold: it crosslinks explicit fibers at geodesic intersections at the physical 60 nm scale, is buckling-capable with soft edges, and needs zero γ/myosin edits. A manifold and a Mikado-sphere solve **different** problems (spreading-coordination/efficiency vs woven connectivity). Don't conflate them.

---

## 9. CORTICAL-TENSION MEASUREMENT on the manifold — MOP vs IK, curved-surface hoop

- **Promote IK whole-shell virial `γ_ik` to PRIMARY** curved-surface estimator. It is **already implemented and validated** (`cortical_tension.py:380-386`): `γ_ik = Σ T·L·(1 − (r̂·û)²) / (8πR²)`, **<0.5%** vs corrected MOP, **~12× lower** sampling variance. It sums **every** bond's tangential projection over the whole shell (not only bonds crossing one flat cut-plane), so it degrades far more gracefully on a deformed surface. It is the **same realized bond-tension state** as MOP — it de-noises, it does not fix an undercount the MOP doesn't have.
- **Demote flat-plane MOP to cross-check.** `_method_of_planes_gamma` (`cortical_tension.py:213-258`) uses **flat through-center** cut-planes and `2πR_cell` circumference — both become ill-defined on a spread/deformed cell. (Keep the required `|û·n̂|` absolute value, docstring:52-58 — without it the per-plane sum cancels as `√N` and fabricates a ~97× suppressed floor.)
- **Curved-surface hoop, made shape-agnostic via the manifold:** at each bond midpoint use the **local manifold normal** `n̂` (instead of centroid-radial `r̂`) to build the tangent plane, project the bond into it, and normalize by the **realized manifold area** (`Σ` triangle areas) instead of `4πR²`. The hoop direction is then defined locally and the estimator is shape-agnostic.
- **Third independent cross-check:** HOOMD `pressure_tensor` / `σ_ab` virial (already planned, `H7_SIGMA_A:57`).
- **Guardrail:** the manifold confinement force `F_conf` is a **localization** spring and is **EXCLUDED** from every γ sum. γ stays a sum over explicit filament/motor/xlink bond tensions only.

---

## 10. VALIDATION GATES

1. **Manifold-OFF bit-identity.** With the manifold disabled (default-OFF), the build is **bit-identical** to current `connected_mesh` (additive force pattern, mirror ERM/enclosed_volume append).
2. **Zero-tension-injection invariant (decisive).** With **myosin OFF + turgor OFF**, `γ_soft + γ_rigid + γ_ik` measured with the manifold **ON** must equal the value with the manifold **OFF** to **<1%**. The manifold must inject **zero** cortical tension.
3. **No manifold bond types.** Sanity-Gate test asserting `bonds.types` contains **no** manifold/edge/area type, so `cortical_bond_typeid_mask` (`cortical_tension.py:150-180`) never sees a manifold edge.
4. **Mesh-resolution independence (grid-convergence, the falsification tripwire).** Hold the **physical** cell fixed (same `n_filaments`, `n_xl`, `n_motors`, seed); refine `n_tri ∈ {N/4, N/2, N, 2N, 4N}`; assert `γ_soft`, `γ_ik`, coordination `z`, giant-component fraction, and total bind-count **converge to a plateau** (successive relative change < a few %). **HARD invariance guard:** because candidates are filtered by the **physical** reach, these MUST be invariant under `n_tri`. **If any tracks `n_tri`, mesh resolution has leaked into the physics → invalid build, halt.** Enforce `patch_size ≥ physical_reach` so a k-ring always covers the reach ball at every resolution.
5. **Candidate-set identity (prototype gate).** Patch k-ring candidate **sets** == global cKDTree candidate sets at a fixed reach (proves the index is faithful, no physics change).
6. **Connectivity gates retained** (from `connected_mesh`): `z ∈ [3.0, 3.5]`, `giant ≥ 0.9`, `L/lc ≥ 5.9`.
7. **No double-count.** Decomposition test asserting `F_conf` does not duplicate the membrane `2γ_mem/R` Laplace (`membrane_surface.py:34-47`) or the enclosed-volume turgor (`enclosed_volume.py:484`) on the same tags (sibling of the existing composite-tension test).
8. **CFL / escape.** Per-step inside-band assertion; `dt ≤ safety·τ_conf`.

**No gate-loosening, no magic numbers without a block, no tuning of `k_conf`/`d0`/`n_tri` to pass γ or connectivity.**

---

## 11. RESOURCE ESTIMATE

- **Search cost (analytic; prototype must MEASURE).** Current: per batch tick, cKDTree build O(N log N) over ~7000 beads + `query_ball_point` of `n_heads` (xlink heads up to ~2000; myosin heads `2·28·100 = 5600` MCF7) — rebuilt every `batch_steps = 100` and discarded. Patch k-ring: O(k) per query, no rebuild. Analytic ~50–250× **search** speedup; **not** a whole-step speedup (BAOAB + Fixman + M-SHAKE dominate per-step cost).
- **Memory.** `tri_adj (n_tri,3)`, `bead_tri (n_cortex_actin,)`, per-face frames `(n_tri,3,3)` — all O(N), negligible vs the ~7000-bead state. **GPU caveat:** a naive CPU nearest-point query re-adds the host sync the GPU port removes; `bead_tri` and the patch index must be **GPU-resident** (cupy) for production. Lazy refresh (every K batches) keeps re-fit cost amortized.
- **Implementation cost.** ~1 file for the manifold + index; one `md.force.Custom` for `U_conf`; reuse `_tangent_plane_basis`. No integrator change, no AFINES/myosin/xlink mechanism edits.

---

## 12. MIGRATION PLAN from `connected_mesh`

Incremental, additive, each step proven before the next:

- **Step 0 — Prototype (no simulator change).** `scripts/manifold_search_prototype.py` (§ smallest prototype below). Validates data model + search win + the §8 honesty claim + a slice of the §10 resolution gate in one file.
- **Step 1 — Faithful index refactor (zero physics change).** Wrap `seed_connected_mesh_xlinks` cKDTree (`crosslinkers.py:613`) behind a `ManifoldIndex` returning **byte-identical** candidate sets. Prove with the connectivity gates (gate 6).
- **Step 2 — Per-tick search swap.** Replace the per-tick cKDTree in `XlinkBondUpdater.act` + `MyosinStepUpdater.act` with `ManifoldIndex.kring_query`; assert **identical bind/break counts** on a fixed seed (search-cost win, no physics change). Promote the full faithful `connected_mesh` build to production (`H7_CORTEX_AS_MESH` rec 1a) at the same time.
- **Step 3 — Soft normal confinement.** Add `U_conf` as an **additive default-OFF** `md.force.Custom`; mirror ERM/enclosed_volume append; preserve tag-append invariant (`cell.py:914`). Run gates 1, 2, 7, 8.
- **Step 4 — Deformable manifold + shape coupling.** Have `enclosed_volume`/`membrane_surface` read normals from the manifold; enable two-way slaving. Run gate 4 (resolution independence) as the production acceptance test.

`√(A/n)` compensation is **retained throughout** (it is not a migration target).

---

## 13. RISKS + OPEN PI DECISIONS

**Risks**
- **Marketed as a γ fix (highest-stakes).** It is not. The floor is Gate A + Gate B, orthogonal to a coordinate manifold. Selling it as a fix repeats the refuted topology/measurement-blaming pattern (`H7_SIGMA_A:7-16`).
- **Triple-counted surface tension.** A manifold normal force would be a **third** Laplace force on the same beads, atop `enclosed_volume` (turgor) and `membrane_surface` (membrane γ). Mitigation: manifold supplies normals only; existing modules stay sole force owners (§5, gate 7).
- **`k_conf` magic-number creep.** Same knob class as `k_ERM`; strong temptation to tune it to "help connectivity" or nudge γ. Needs the Magic-Number Block or PI surface (§6).
- **Resolution leaking into γ.** If γ moves with `n_tri`, the manifold is in the tension budget — gate 4 is the decisive falsification.
- **Decoupling / lag.** Stale `d_normal` if re-fit cadence is slower than bead dynamics → beads escape or snap back. Mitigation: two-way coupling, `τ_manifold ≳ τ_conf`, re-query nearest point each step (cheap with patch id).
- **Deforming-manifold vs remodeling-cortex.** Active spreading is polymerization growth (adds area); a fixed-topology manifold resists it → needs dynamic re-meshing/node insertion at the front (`H7_CORTEX_AS_MESH:41-49`).
- **GPU host-sync regression** from naive CPU nearest-point query (§11).
- **Unverified citations.** Ronceray-Broedersz-Lenz 2016 (10.1073/pnas.1514208113) and Chugh 2017 appear cited-from-memory in source docs; the project has a confirmed-hallucination history. Verify before any Notion SourceEvidence registration. Two literature DOIs not fully fetched (Membrane-MEDYAN 10.1021/acs.jpcb.1c02336 paywalled; Cytosim 2007 NJP via secondary) — re-verify.

**Open PI decisions**
1. **Build now or DEFER?** The critique investigation recommends **DEFER as a γ effort** and attack Gate A (construction-time binding pre-equilibration) + Gate B (M-SHAKE relax, integrator-freeze sign-off) on the existing explicit network first; build the manifold later strictly as a coordinate/search/spreading substrate. The systems/model investigations support building it now as a scoped efficiency/spatiality layer. **PI call needed.**
2. **`k_conf` value** — thin pin (`= k_ERM`, σ_n 0.65 nm) vs a derived softer band-filling value (fraction of `h_cortex/2`). Magic-Number Block + sign-off.
3. **M-SHAKE relaxation** for buckling — FROZEN integrator, requires PI sign-off; the manifold does not unlock buckling on its own.
4. **Manifold vs Mikado-on-sphere** for connectivity — different problems; pick scope.
5. **Single owner of the surface Laplace pressure** — confirm `enclosed_volume` + `membrane_surface` remain sole owners, manifold passive.

---

## Answers to ALL 10 PI questions

1. **Does the manifold-with-explicit-filaments design satisfy the fine-grained, mechanistic hard rule?** **YES — conditionally.** From literature precedent (Cytosim soft-wall; Bidone-Tang-Vavylonis 2014 soft normal cortical force, 10.1016/j.bpj.2014.10.034 — closest single precedent; Membrane-MEDYAN 10.1021/acs.jpcb.1c02336 deformable triangulated substrate + explicit filaments; FreeDTS 10.1038/s41467-024-44819-w; Janzen-Sknepnek 2026 10.1039/D5SM00884K explicit filaments on a curved manifold) every ingredient is precedented; the integrated whole is a **novel union**, not an off-the-shelf method. It is categorically distinct from the **forbidden** lumped continuum active-surface theory (Salbreux-Julicher 2017, PMID 29346890). **Compliant IF AND ONLY IF** the manifold supplies only geometry/locality/soft-confinement and **every** force including γ comes from explicit objects and is measured from explicit bond stress. It becomes a forbidden second coarse-graining the instant any manifold DOF carries mechanical load. (The one sanctioned coarse-graining stays the ×40 mesoscale.)

2. **Minimal architecture that avoids replacing actin mechanics with mesh mechanics?** The manifold gets **exactly three** responsibilities, all coordinate/bookkeeping, zero mechanical: (1) **geometric tracker** (deforming surface exposing `n̂`, `d_normal`, `d0`, `(e1,e2)` — generalizes `_project_to_shell_band` + `_tangent_plane_basis` to a deformable surface); (2) **locality/patch index** (accelerated binding-candidate search + shared "where is the surface" frame for membrane/turgor/FA/lamellipodium); (3) **one soft normal confinement** `U_conf`. It does NOT carry the Laplace pressure (stays in `enclosed_volume` + `membrane_surface`, now reading manifold normals), owns no edge/area elastic energy in the force budget, and does not define γ. Manifold shape is **slaved to the bead cloud**.

3. **Bead representation — (a) global Cartesian + soft confinement, or (b) local triangle-id/barycentric?** **(a) global Cartesian + soft confinement.** HOOMD integrates global Cartesian; the constrained BAOAB action reads `pos + net_force`, so confinement is just another `md.force.Custom` (like `erm.py`/`enclosed_volume.py`) — zero integrator change, GPU-port compatible. (b) needs per-step local↔global projection that fights the Cartesian core, re-adds the GPU host sync, and complicates buckling. The manifold is an auxiliary field; integrated state stays Cartesian. (§4)

4. **How to keep beads near the shell while permitting buckling?** Soft normal-only `U_conf = ½ k_conf (d_normal − d0)²`, `d0` at shell mid-plane, half-width `h_cortex/2 = 100 nm`. **`k_conf` anchored to ERM `k_ERM = 0.1 N/m` (KU-3.18)** — same knob class, **[MAGIC-NUMBER BLOCK REQUIRED]**; softer band-filling value derived as a stated fraction of `h_cortex/2` if needed (PI sign-off). SOFT (normal-only) leaves tangential motion free → in-plane buckling + condensation + normal undulation survive; hard projection would kill them. **Caveat:** buckling is only *permitted*, not *unlocked* — the lever is relaxing the M-SHAKE backbone (integrator-freeze, PI sign-off), independent of the manifold. (§6)

5. **Current neighbor-search, and what's needed for patch/geodesic search?** Today: **no** triangle adjacency / bead→triangle map; every binder rebuilds a **global cKDTree** over ~7000 beads each batch tick (`crosslinkers.py:1060-1064`, `myosin.py:1099/1193`, seeding `crosslinkers.py:613/647`, Arp2/3 `cortex.py:1332`). Only spatial metadata = flat `filament_idx` + a 1D bead→segment chain (`myosin.py:857-863`). Needed: `tri_adj (n_tri,3)`, `bead_tri` (lazy-refresh), per-face frames (reuse `_tangent_plane_basis`); search = patch ∪ k-ring via BFS = O(k) vs O(N log N). Construction bound `patch_size ≥ physical_reach`. (§7)

6. **Does patch-local search at physical 60 nm reach replace the `√(A/n)` bridge hack? (load-bearing)** **NO.** The manifold changes search **cost**, not physical **spacing**. ×40 sparsity is geometric (`~√(A/n) ≈ 0.9–1.9 µm` apart, index-independent). At 60 nm the k-ring returns the **same empty sets** as the current cKDTree — just faster. 2D test is decisive: 60 nm → z 0.20 / giant 1.5% (fragmented) vs `√(A/n)` → z 9.4 / giant 99% (percolated). UV/geodesic binding only swaps chord for arc; in-shell nearest-neighbor is still `~√(A/n) ≫ 60 nm`. **Keep `√(A/n)`** as the documented sanctioned coarse-graining compensation. (§8)

7. **Cortical-tension measurement on the manifold — MOP vs IK; curved hoop?** **Promote IK whole-shell virial `γ_ik` to primary** (already implemented `cortical_tension.py:380-386`, <0.5% vs MOP, ~12× lower variance); demote flat through-center MOP to cross-check (its "center"/`R_cell`/`2πR` degrade on a deformed cell). Make `γ_ik` shape-agnostic: at each bond midpoint use the **local manifold normal** for the tangent plane, normalize by **realized manifold area** (`Σ` triangle areas) not `4πR²`. Add HOOMD `σ_ab` virial as a third cross-check. Manifold confinement force is **excluded** from every γ sum. (§9)

8. **Mesh-resolution independence — how validated?** Grid-convergence study: fix the physical cell, refine `n_tri ∈ {N/4…4N}`, assert `γ_soft`, `γ_ik`, `z`, giant-fraction, bind-count converge to a plateau (successive change < a few %). **HARD guard:** candidates are filtered by physical reach, so these MUST be invariant under `n_tri`; if any tracks `n_tri`, resolution leaked into physics → halt. Enforce `patch_size ≥ physical_reach` (k scales with resolution). This is the decisive falsification that the manifold is coordinate-only. (§10 gate 4)

9. **New failure modes the manifold introduces?** (1) **Triple-counted surface tension** — manifold normal force atop `enclosed_volume` turgor (`F_i=ΔP·A_i·n̂_i`) + `membrane_surface` Laplace (`F_i=−ΔP·(S/N)·n̂_i`) on the same beads; fix = manifold supplies normals only. (2) **Manifold-bead decoupling/lag** (stale `d_normal`) → escape or snap-back; fix = two-way coupling + `τ_manifold ≳ τ_conf` + per-step re-query. (3) **Manifold as unphysical tension source** if it has independent elastic DOF; forbidden. (4) **`k_conf` hidden knob** (too stiff → re-creates M-SHAKE buckling kill; too soft → leak); forbidden to tune for γ/connectivity. (5) **Resolution affecting γ** (gate 4 tripwire). (6) **Deforming manifold fighting remodeling cortex** during spreading (needs re-meshing). Plus CFL/escape under strong active load + large `dt`. (§9-critique, §13)

10. **Smallest prototype + comparison metrics?** One standalone `scripts/manifold_search_prototype.py` (no Updater/integrator/γ-gate edits) — full spec below.

---

## Smallest prototype spec

**File:** `scripts/manifold_search_prototype.py` (single file; reuses `generate_bimodal_cortex_layout` + `connected_mesh` diagnostics; touches **no** Updater, **no** integrator, **no** γ gate).

**Steps:**
1. Generate the existing bimodal cortex bead cloud via `generate_bimodal_cortex_layout` (MCF7: `R_cell = 7.5 µm`, `n_filaments = 1000`).
2. Build the triangulated manifold (icosphere/Fibonacci, `n_tri ≈ n_filaments`) + `bead_tri` map + `tri_adj`.
3. Run **BOTH** the current global cKDTree candidate query **AND** the patch k-ring query for the crosslinker-seeding pass.
4. **Assert candidate SETS are identical** at a fixed reach (proves faithful index, no physics change); report wall-time speedup + candidate-count distribution.
5. Sweep reach ∈ {60 nm, `√(A/n)`} to reproduce z/giant fragment-vs-percolate, numerically demonstrating the §8 claim that the manifold does not change spacing.
6. Sweep `n_tri ∈ {N/4, N/2, N, 2N, 4N}` at fixed reach → show candidate sets / z / giant invariant (slice of gate 4).

**Comparison metrics vs current pure-3D `connected_mesh`:**

| Metric | Current `connected_mesh` (cKDTree) | Manifold prototype | Pass criterion |
|---|---|---|---|
| Connectivity `z` | ~3.3 | same | identical at fixed reach |
| Giant-component fraction | ~99% | same | identical at fixed reach |
| Candidate count distribution | baseline | report | sets **identical** at fixed reach |
| Search wall-time per query | O(N log N) build + query | O(k) | measured speedup (target ~50–250×) |
| `γ_soft` / `γ_ik` | baseline | (in full build, later) | gate 2: <1% with manifold ON vs OFF (myosin/turgor OFF) |
| Buckling `r/r0` | 1.000 (M-SHAKE) | 1.000 (unchanged) | **unchanged** — manifold does NOT alter it |
| Mesh-resolution sensitivity | n/a | sweep `n_tri` | candidate sets / z / giant **invariant** under `n_tri` |

---

## Honesty statement (acceptance criteria, enforced)

This design **improves spatiality** (a real 3D deforming surface + local frames), **efficiency** (O(k) patch search vs per-tick global cKDTree rebuild), **connectivity organization** (principled curved-surface locality), and **spreading coordination** (shared coordinate frame for membrane/turgor/FA/lamellipodium). It is **ORTHOGONAL to the active-γ floor**: the floor is Gate A (heads load to 0.2–0.44 F_stall but `s_grip ≈ 0.012 nm ≪ ℓ0 = 500 nm` → don't walk → bond extensions never grow) and Gate B (M-SHAKE `r/r0 = 1.000` forbids the ~10× buckling amplification). A coordinate manifold adds **zero** myosin walking and **zero** bond extension; `γ_soft` is a sum over explicit bond tensions and a coordinate manifold contributes none. **It must not be sold or tuned as a γ-floor fix.** Acceptance criteria enforced: explicit objects stay explicit (§5); buckling possible (soft confinement, §6); γ from explicit stress only (§9, gate 2/3); resolution-independence testable (gate 4); no magic numbers without a block (`k_conf` flagged `[MAGIC-NUMBER BLOCK REQUIRED]`, §6); no gate-loosening (§10). The active-γ floor is attacked elsewhere (Gate A construction-time binding pre-equilibration; Gate B M-SHAKE relaxation with PI integrator-freeze sign-off) — not here.

---

**Synthesis notes for the caller:** All four investigations are mutually consistent and well-grounded in the actual code (I verified `cortical_tension.py:380-386` IK virial, `connected_mesh.py:57-62` `√(A/n)`, `erm.py` `k_ERM=0.1 N/m`). The one substantive divergence the PI must resolve is **timing**: the critique recommends DEFER (attack Gate A/B first); systems/model support building the manifold now as a scoped efficiency/spatiality layer. Both agree it is orthogonal to γ and must be guardrailed. The single highest-value, lowest-risk next action regardless of that decision is the **Step-0 prototype** (`scripts/manifold_search_prototype.py`), which validates the data model, the search win, and the load-bearing §8 honesty claim in one file with no simulator changes. Two value flags carry a `[MAGIC-NUMBER BLOCK REQUIRED]` (`k_conf`) and a verify-before-SE warning (Ronceray-Broedersz-Lenz 2016, Chugh 2017, Membrane-MEDYAN/Cytosim DOIs) given the project's confirmed-hallucination history.