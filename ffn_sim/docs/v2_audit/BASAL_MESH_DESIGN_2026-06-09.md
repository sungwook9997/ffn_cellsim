# Connected Basal Contractile Mesh — design (PI directive (a), 2026-06-09)

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

* **B1** basal-disk layout generator — bimodal in-plane filaments on the basal disk
  (reuse `generate_bimodal_cortex_layout` math with a planar projection; long cables
  along the long axis). Build-time: filaments assembled, in-plane, force-free.
* **B2** connect the basal mesh — reuse `seed_connected_mesh_xlinks` on the basal
  filaments. Gate: giant ≥ 0.9, z ∈ [3, 3.5].
* **B3** FA anchoring — cable ends → FA integrin clutches (reuse `sf_anchor`).
* **B4** `sf_myosin_` NMII placement on the basal mesh (①a/①b). Build-time: assembled,
  force-free, no-contam.
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
