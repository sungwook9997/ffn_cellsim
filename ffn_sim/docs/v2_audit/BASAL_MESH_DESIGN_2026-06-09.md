# Connected Basal Contractile Mesh — design (PI directive (a), 2026-06-09)

> **2026-06-09 REFRAME (PI clarification).** The intent is a **TWO-LAYER**
> architecture, and it lands EXACTLY on the already-sanctioned manifold-guided
> design (`H1_H5_H7_MANIFOLD_CONTACT_ARCHITECTURE_2026-06-07.md`) — NOT on the
> force-bearing surface that `H7_CORTEX_AS_MESH_2026-06-07.md` (commit `7b19276`)
> rejected:
>
> - **S-layer — 2D SURFACE MESH (SimuCell3D-style, `cortex/surface_manifold.py`).**
>   Its job is to give the FA particles a real 2D surface to be POSITIONED on (FA is
>   otherwise hard to attach — the spatial-disjoint integrin problem), plus SOME
>   surface connectivity ("표면 연락 일부" — patch adjacency + an optional soft
>   **normal-only** confinement `U_conf`). **The surface carries NO in-plane /
>   edge force** — a force-bearing triangulation edge is the 2nd unsanctioned
>   coarse-graining and over-constrains single-filament buckling (the Gate-B lever),
>   rejected by `7b19276`. The surface is geometry / positioning / contact ONLY,
>   invisible to the γ estimator (not a bond).
> - **F-layer — explicit FILAMENT NETWORK ON the surface (the REAL mechanics).**
>   The actual FA forces/traction come from a fine-grained actin filament network
>   (this doc's B1/B2) anchored to the FA particles ON the surface, contracted by
>   `sf_myosin_` NMII. "실제 FA 메카닉스는 그 표면 위에서 필라멘트 네트워크를
>   구축해서 거기서 실제 힘을 본다" (PI). This is the sanctioned
>   "explicit fibers woven on the manifold" option; the manifold consumes its
>   normal/patch services unchanged.
>
> Bright line (from `7b19276`, restated): **a surface edge that PULLS two nodes is
> physics (FORBIDDEN); a surface that POSITIONS particles, says "search this patch",
> and softly keeps a bead in the band is geometry/broad-phase (ALLOWED).** The B1/B2
> filament work below is the F-layer (reusable); what changes is it is placed ON the
> S-layer surface and anchored to FA positioned on it.

PI chose **(a)**: lay down ONE connected actin mesh on the basal plane, weave the
ventral stress-fiber cables INTO it, anchor it to the substrate through FA, and let
NMII contract the whole mesh to generate traction. This resolves the single-chain /
bipolar-gate blocker (a connected mesh has many filaments → the Stam-Hocky bipolar
NMII gate works unchanged: +/− head-sets grip DIFFERENT basal filaments).

## The key reuse — basal mesh = cortex connected-mesh on a FLAT DISK

The cortex is ALREADY a connected spanning mesh built from:
* `generate_bimodal_cortex_layout` — BIMODAL filament length (short Arp2/3 infill +
  long formin backbone) + isotropic orientation, projected onto the **sphere shell**.
* `seed_connected_mesh_xlinks` — BRIDGE-DIFFERENT-FILAMENT crosslinkers, per-filament
  z≈3-4, bundling → giant ≥ 0.9, L/lc ≥ 5.9.

The basal mesh is the SAME construction projected onto a **flat basal disk** instead
of the sphere shell. The biology maps cleanly:
* **short infill filaments**  → the connecting basal actin meshwork (the "mesh").
* **long filaments**          → the ventral stress-fiber CABLES (FA→FA), "woven in".
* **bridge xlinks**           → α-actinin/filamin cross-bundle bridges (connectivity).
* **NMII minifilaments**      → `sf_myosin_*` (prefix-split, ①a) grip the connected
                                basal sf_actin → bipolar contraction → traction.

So (a) is largely a GEOMETRY swap (sphere-shell projection → basal-disk projection) +
re-anchoring to FA + the `sf_` γ-denylist, NOT a from-scratch mechanism. Every new
particle/bond stays `sf_`-prefixed → excluded from cortical hoop γ (KU-3.5) by the
registry denylist; the basal contraction is reported as its own observable.

## Geometry (the one genuinely new piece)

```
   side view (x–z):                         top view (x–y), basal disk z≈z_basal:
                                              FA●          ●FA
     cortex sphere shell                        \  ____  /
         ___________                             \/    \/      long cables = SF
        /           \                            /\    /\      (FA→FA, long-axis)
       |             |                           /  ‾‾‾‾  \     short infill = mesh
       |    nucleus  |                         FA●  ·····  ●FA  (connecting network)
        \___________/   ← basal cap
   ======●==●==●==●====  FA integrins (clutch → cortex south cap)
   ----- substrate z=0 (pinned ligands) -----
```

* The basal mesh is a **flat disk** at `z ≈ z_basal` (the south-cap contact zone),
  radius = the FA contact-footprint radius. Filaments lie IN the basal plane
  (in-plane isotropic orientation; long cables biased along the cell long axis =
  the existing `select_aligned_fa_pairs` PCA axis).
* **FA anchoring**: cable ENDS anchor to FA integrin clutches (`sf_anchor`, force-free,
  reusing the current SF→FA path). Interior mesh nodes are NOT FA-anchored (free to be
  pulled inward — that inward pull IS the traction the FA clutches resist).
* **NMII**: `sf_myosin_` minifilaments placed on the basal sf_actin via the actin-aware
  layout (①b actin_pool_tags = the basal sf_actin tags); grip_walk + the bipolar gate
  now satisfiable (many basal filaments).

## Observables / gates (stated before any run)

1. **Connectivity** (build-time): giant-component fraction ≥ 0.9, z ∈ [3, 3.5] on the
   basal sf_actin mesh (same gate as cortex connected-mesh).
2. **No cortical-γ contamination** (build-time): cortical γ_soft IDENTICAL basal-mesh
   OFF vs ON (all basal bonds `sf_`-denylisted).
3. **Force-free construction** (build-time): every basal bond born at its exact r0.
4. **Single-SF tension** (equilibrated, DEFERRED): Kumar 2006 10-30 nN per cable.
5. **Traction** (equilibrated, DEFERRED): basal stress through FA ≈ Balaban 5.5 nN/µm².
   (4) + (5) need the equilibration prelude + N_filaments/k_actin ratification.

## Implementation increments (each its own commit + build-time gate)

**F-layer (explicit filament network — the real mechanics):**
* **B1 ✅ DONE** (`6f67b99`) — basal-disk filament layout (bimodal cables + infill;
  `cell/basal_mesh.py:generate_basal_mesh_layout`). Gate PASS: planar, force-free,
  cables |cos|=1.0, infill |cos|=0.639. (Currently a free disk; B3 places it ON the
  S-layer surface.)
* **B2 ✅ DONE** — connect the filament network (`connect_basal_mesh`, reuses
  `seed_connected_mesh_xlinks` with the disk reach). Gate PASS: giant=0.998,
  z=3.087∈[3,3.5], L/lc=6.18≥5.9, n_xl=1850.

**S-layer (surface manifold — positioning + connectivity, NO force):**
* **S1 ✅ DONE (FLAT, PI 2026-06-09)** — flat ventral surface layer
  (`cell/basal_surface.py`): an adherent cell FLATTENS its ventral surface against the
  substrate, so the basal surface is a **FLAT 2D triangulated disk** (Delaunay over
  concentric-ring disk points at z=z_basal, +z normal), NOT a curved sphere cap (=
  suspended geometry, wrong here; also flat eases the lamellipodium). POSITIONS FA
  particles ON the flat surface (z=z_basal) at triangle centroids + exposes triangle
  adjacency = the "표면 연락 일부". NO force, not a bond → γ-invisible. Gate PASS: 629
  tris, 60 FA planar + within-disk + connected, area 99.8% of disk; VG-1 master gate —
  disk area invariant to ring resolution (rel spread 0.4%). This consistency-fixes the
  earlier curved-cap S1 (the surface and the planar F-layer filaments now share one
  flat geometry — no flat-vs-curved dichotomy). `disk_area_resolution_invariance` =
  VG-1 check. (Optional soft normal-only `U_conf` deferred to integration.)

**Integration + active:**
* **B3** anchor the F-layer ON the S-layer — filament network beads sit on the
  surface; cable ends → FA integrin clutches positioned on the surface (`sf_anchor`).
* **B4** `sf_myosin_` NMII placement on the surface-borne filament network (①a/①b).
  Build-time: assembled, force-free, no-contam (sf_ γ-denylisted).
* **B5** equilibrated active gate — equilibration prelude → Kumar single-SF tension +
  Balaban traction. (Needs PI N_filaments/k_actin + likely gbook GPU run.)

## Open sub-decisions for PI (small)
* **Q1.** Dedicated basal sf_actin mesh layer (this design) vs. NMII gripping the
  EXISTING basal cortex shell beads. Recommend the DEDICATED layer (clean γ separation:
  gripping cortex actin would tense cortex backbone bonds → leak into hoop γ; a separate
  `sf_` mesh keeps traction and cortical tension cleanly distinct). PI wording
  ("바닥면에 연결된 액틴 mesh를 하나 깔고") already implies the dedicated layer.
* **Q2.** Is this a NEW compartment (`basal_actin_mesh`) in the registry, or an
  EXTENSION of `ventral_stress_fibers` (the mesh subsumes the cables)? Recommend
  EXTENDING `ventral_stress_fibers` (the cables ARE the long filaments of this mesh;
  one concept, one module) — rename intent: "ventral_stress_fibers" = the basal
  contractile mesh incl. its cables.
