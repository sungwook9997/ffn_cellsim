# H.7 — Lamellipodium ↔ Spherical-Cortex Integration (active-spreading enhancement)

**Status**: DESIGN (not yet implemented). **Branch**: `h7/full-cell-integration`.
**Owner**: Lead. **Scope (PI 2026-06-06)**: ACTIVE-SPREADING track for the full physiological MCF7 cell.
**Reads / grounding**: `aleph/cell/lamellipodium.py`, `aleph/cell/membrane.py`,
`aleph/cell/cell.py`, `aleph/cortex/cortex.py`, `aleph/configs/{phase1_h5,mcf7_baseline}.yaml`.

> This document is a geometry-integration design problem statement + candidate
> approaches + recommendation. It writes NO code. Every claim is grounded at `file:line`.

---

## 0. TL;DR

The H.5 lamellipodium is a **flat leading-edge model**: a WAVE/NPF plane at the box top
`y = +Y_max` (default `0.45·L_box ≈ +13.5 µm`), spanning `wave_area = L_box²`, with mother
barbed ends seeded one rest-length below and growing toward `−ŷ` *into* the cytosol; the
membrane it pushes is a tension plane with normal `+ŷ`. The full cell is a **closed
spherical cortex** (`R_cell = 7.5 µm`, origin-centred) adhered to a substrate plane at
`z ≈ 0` through FA at its **south cap** (`z ≤ z_min + cap_depth`, `z_min ≈ −R_cell`). The
flat-patch-in-a-box geometry does not map onto a sphere spreading on a substrate: the WAVE
plane sits at `+y` (top of the box, several µm of vacuum away from any cortex bead and
orthogonal to the substrate), whereas a spreading cell's true leading edge is the **basal
contact ring / lamellar zone at the cell–substrate interface** (around the south cap, near
`z = 0`, growing **radially outward in the basal plane**). Recommendation: **Approach (a) —
a basal lamellar patch at the contact ring**: relocate the WAVE plane to the basal substrate
plane and reorient barbed-end growth to radially-outward-in-the-basal-plane, keeping the full
Arp2/3 mechanistic D2 machinery. This is the highest-fidelity option and the smallest
*conceptual* change to the existing mechanism (it re-parameterises geometry, not physics).

This integration is the **active-spreading enhancement track and does NOT block KU-3.5 /
GATE-B** (see §5).

---

## 1. Geometry mismatch (precise)

### 1.1 The flat WAVE-plane patch (current H.5)

- **WAVE plane position.** `resolve_h5_lamellipodium` sets `Y_max = cfg.get("Y_max",
  0.45 * L_box)` (`lamellipodium.py:227`) and validates `|Y_max| < 0.5·L_box`
  (`lamellipodium.py:279-283`) — i.e. the plane is a horizontal slice near the **top of the
  cubic box**, normal `+ŷ`. `phase1_h5.yaml:26` pins `Y_max = 13.5e-6 m = 0.45·L_box`
  (L_box = 30 µm).
- **Patch area = whole box face.** `wave_area = cfg.get("wave_area", L_box * L_box)`
  (`lamellipodium.py:228`; `phase1_h5.yaml:27` = `(30 µm)² = 9·10⁻¹⁰ m²`). WAVE count
  `n_WAVE = round(sigma_NPF · wave_area)` (`lamellipodium.py:230`).
- **WAVE placement = uniform on the (x,z) plane.** `generate_lamellipodium_layout` scatters
  WAVE beads at `x,z ∈ U(−0.475·L_box, +0.475·L_box)`, `y = Y_max`
  (`lamellipodium.py:352-354`). Each WAVE gets ONE mother seed one rest-length below at
  `y = Y_max − ℓ₀` (`lamellipodium.py:356-358`).
- **Barbed ends grow toward −ŷ (into cytosol, away from membrane).** Mother tangent is
  hard-coded `[0, −1, 0]` at construction (`lamellipodium.py:1016`,
  `lamellipodium.py:1297`). Elongation appends a bead at `r_be + ℓ₀·tangent`
  (`lamellipodium.py:539`).
- **WAVE pin is a 1-D harmonic to the y-plane.** `WaveMembranePin` applies
  `F_y = −k_wave_pin·(y − Y_max)`, force only on the y-component
  (`lamellipodium.py:400-408`) — i.e. it confines WAVE to the *plane* `y = Y_max`, free in
  x,z.
- **Membrane load lives on the same +ŷ plane.** `resolve_membrane` defaults `y_plane = Y_max`
  and `A_mem = wave_area` (`membrane.py:327-331`); the membrane normal default is `(0,1,0)`
  (`membrane.py:557`, `MembraneReactionForce` `:735`); the contact test is `y ≥ y_plane −
  contact_range` (`membrane.py:646-649`) and the reaction is `−n̂·F_per` along `−ŷ`
  (`membrane.py:774`). So the WHOLE H.5 force-velocity story is "barbed ends push a flat
  ceiling at the top of the box."

**Net:** the H.5 model is a Bieling/Funk *in-vitro* reconstitution geometry — a flat NPF-coated
coverslip with dendritic actin growing off it against a flat membrane — embedded in the +y
face of a cube. It has no notion of a cell body, a substrate, or a contact ring.

### 1.2 The closed spherical cortex on a substrate (full cell)

- **Shell is origin-centred, radius R_cell.** `_sample_sphere_surface` (Marsaglia)
  returns CoMs with `|p| = R_cell` centred on the origin (`cortex.py:542-566`); filaments are
  laid tangent on that shell (`cortex.py:659-662`). Box `L_box = box_factor·R_cell`
  (`cortex.py:438`; default `box_factor = 3`). So the cortex spans `z ∈ [−R_cell, +R_cell]`
  about the origin; the south pole is at `z ≈ −R_cell`.
- **Substrate is the plane z = 0.** FA substrate ligands are placed at `z = 0` and held
  immobile by `SubstrateLigandPin` (`cell.py:314`, docstring; "ligands placed at z=0 (the
  mechanical ground)").
- **Adhesion is at the SOUTH CAP.** `_extend_snapshot_with_fa` computes `z_min = min(cortex
  z)` (`cell.py:346`) and the contact footprint = cortex beads with `z ≤ z_min + cap_depth`
  (`cell.py:362-363`), seeds integrins under those cap beads, and forms clutch bonds to the
  nearest within-capture south-cap actin bead (`cell.py:452-465`). The module docstring
  states the honest geometry: "the cortex shell is origin-centred (south pole at z ~ −R_cell)
  while the substrate plane is at z ~ 0" (`cell.py:332-334`).

### 1.3 Where the leading edge physically sits when a sphere spreads on a substrate

A cell settling/spreading on a flat substrate is **not** pushing a membrane at its apex (+y or
+z). The protrusive machinery (lamellipodium) lives at the **basal periphery**: the thin
lamellar zone where the ventral membrane meets the substrate, i.e. the **contact ring** at the
edge of the cell–substrate footprint, near `z = 0`, around the south cap that the FA already
adheres (`cell.py:362-363`). There, dendritic actin grows **radially outward in the basal
plane** (in the x–y plane at `z ≈ 0`), advancing the contact ring and increasing the footprint
area `A/A₀` — the spreading observable. The membrane the basal barbed ends push is the
**advancing edge of the ventral/peripheral membrane**, whose local outward normal is
**radial-in-the-basal-plane** (pointing away from the cell centre, parallel to the substrate),
NOT `+ŷ` and NOT the substrate normal `+ẑ`.

### 1.4 The three concrete mismatches

| Axis | Flat H.5 patch | Sphere-on-substrate reality | File evidence |
|---|---|---|---|
| **Leading-edge location** | plane `y = +0.45·L_box` (box top, vacuum) | basal contact ring, `z ≈ 0`, around south cap | `lamellipodium.py:227,354` vs `cell.py:362,314` |
| **Growth direction** | all mothers `−ŷ` (one global direction) | radially-outward in basal plane (per-WAVE, az-symmetric ring) | `lamellipodium.py:1016,1297` vs §1.3 |
| **Membrane normal / load axis** | `+ŷ`, flat ceiling | radial-in-basal-plane (per edge segment) | `membrane.py:557,774` vs §1.3 |
| **Patch ↔ cell coupling** | WAVE plane disjoint from cortex (top of box); WAVE↔mother anchor only (`lamel_wave_anchor`) | basal lamella must be mechanically continuous with the south-cap cortex it protrudes from | `lamellipodium.py:1258` vs `cell.py:362` |

There is also a **spatial-disjointness hazard already documented for FA** that this design
inherits: at construction the cross-subsystem WCA is disabled precisely because the substrate
(z≈0) and cortex (r≈R_cell) layers are disjoint (`cell.py:1106-1109`), and the WAVE-plane
beads (at +y) are likewise WCA-decoupled from everything (`cell.py:1086-1098`). Any relocation
of the WAVE plane to the basal zone must preserve force-free / non-exploding construction
(the per-r0-bin / WCA-off discipline the FA clutch uses, `cell.py:429-435`).

---

## 2. Candidate integration approaches

All four keep the integrator frozen (BAOAB, per CLAUDE.md) and reuse the existing D2-batched
updaters where possible. "Code-change surface" lists the functions actually touched.

### (a) Basal lamellar patch at the contact ring (RECOMMENDED — see §3)

**Mechanism sketch.** Move the WAVE/NPF reservoir from the box-top plane to the **basal
substrate plane** (`y_plane → z = 0`, or a thin offset above it), and seed WAVE beads on a
**ring** at the cell-substrate contact radius rather than uniformly across a square. Each WAVE
mother grows **radially outward in the basal plane** (per-WAVE tangent = the outward radial
direction in x–y at that azimuth), so the dendritic array advances the contact ring outward —
literally the spreading mechanism. The membrane the tips push is the advancing peripheral
membrane, normal = local outward radial direction (per-WAVE), fed to the Bell-Evans
elongation/capping rate exactly as today via `set_network_load`. Arp2/3 branching, capping,
and the membrane-load law are **physically unchanged** — only the *frame* (plane → ring,
global −ŷ → per-WAVE radial) changes.

**Fidelity vs the real lamellipodium.** Highest. Keeps the full Bieling/Funk dendritic
mechanism (`BarbedEndElongationUpdater` / `ArpBranchingUpdater` / `CappingUpdater`,
`lamellipodium.py:466,574,788`) AND places it where the real lamellipodium is (basal periphery).
Reproduces the spreading observable (`A/A₀` = footprint growth) as an emergent consequence of
barbed-end advance, not a lumped force.

**Code-change surface.**
- `lamellipodium.py::generate_lamellipodium_layout` (`:329`) — replace the `y = Y_max` square
  scatter (`:352-354`) with a ring at the contact radius in the basal plane; mother seeds one
  ℓ₀ *inward* (`:356-358`).
- `lamellipodium.py::extend_cortex_snapshot_with_lamellipodium` (`:1141`) — set per-WAVE
  mother tangent to the outward radial direction instead of the hard-coded `[0,−1,0]`
  (`:1297`); store per-WAVE tangents in `LamellipodiumState.tangent_of` (already per-tag,
  `:436`).
- `lamellipodium.py::WaveMembranePin` (`:379`) — generalise the 1-D `y`-plane harmonic
  (`:400-408`) to pin WAVE to the basal plane (a z-plane harmonic) — or to the ring — without
  constraining the in-plane radial coordinate.
- `membrane.py` — the load law is geometry-light (`F = γ_mem·s_fil`, `:413-452`) so it largely
  carries over, BUT `_contact_set` (`membrane.py:633-653`) and `MembraneReactionForce`
  (`:746-790`) assume a single global normal and a planar `y ≥ y_plane − contact_range` test;
  these need a **per-WAVE / per-edge-segment radial** contact test + normal. New resolver
  fields (basal plane, contact radius) on `ResolvedMembrane` (`:272`).
- `cell.py` — the lamellipodium wiring at `:880-892`, the WAVE-pin force append at
  `:1202-1209`, and the membrane attach at `:1428-1440` pass through unchanged in shape; only
  the resolved-geometry inputs change. WCA pairing (`:1085-1098`) unchanged.
- New config block under `cortex.lamellipodium` for the contact-ring geometry
  (radius, basal-plane z, ring azimuthal count) — additive, replacing `Y_max`/`wave_area`
  usage in the full-cell config; the standalone KU-5.x reconstitution config (`phase1_h5.yaml`)
  keeps the flat-plane defaults.

**CFL / stability risk.** Low. The advance is still a Bell-Evans *rate* (no new stiff force);
the WAVE pin keeps its ERM-style CFL gate (`lamellipodium.py:99-101`). The radial-growth tips
live near z=0 where the substrate ligands and south-cap FA already are — the construction-time
spatial-disjointness that motivated WCA-off (`cell.py:1106`) is *reduced*, not worsened. Main
care item: per-WAVE radial tangents must be recomputed as tips advance around curvature (a tip
that grows straight will drift off the true radial as the ring expands) — but at lamellipodial
length scales over a single spreading episode this is a small correction (the same
straight-tangent approximation the cortex tangent-plane placement already accepts,
`cortex.py:24-26`).

**Composition with FA + membrane-load.** Natural. The basal lamella sits in the **same**
south-cap / z≈0 zone the FA clutch occupies (`cell.py:362`); spreading (lamella advancing the
ring) and adhesion (FA gripping the substrate) are co-located and mutually reinforcing — the
biological picture (protrude → adhere → advance). Membrane-load (KU-5.x, PI-approved
`membrane.py:217`) loads the radial barbed ends exactly as it loads the −ŷ tips today, once
its contact test is made radial.

### (b) Keep the flat-patch model, reposition/reorient to the basal zone (MINIMAL)

**Mechanism sketch.** Leave the planar WAVE model intact but **rigidly re-place** it: set the
WAVE plane to a vertical plane tangent to one side of the cell at the contact zone (or
translate `Y_max`'s role to a `z`-offset patch at the south flank) and rotate the global growth
direction so mothers grow toward that side's outward direction. Still ONE global growth axis
and ONE flat membrane normal — just rotated/translated from `+ŷ` to a basal-edge direction.

**Fidelity.** Medium. Models a **single migrating front** (one lamellipodium on one side of the
cell, as in directed migration) rather than the axisymmetric spreading ring. Faithful to a
polarised, motile cell; NOT faithful to symmetric spreading (which needs a ring, approach (a)).

**Code-change surface.** Smallest. `generate_lamellipodium_layout` (`:329-373`) repositions the
square patch; the mother tangent (`:1016,1297`) and the membrane normal
(`membrane.py:557`) become a config-supplied unit vector instead of `±ŷ`. `WaveMembranePin`
(`:379-408`) generalises to pin to an arbitrary-orientation plane (dot-product form). `cell.py`
wiring unchanged.

**CFL / stability risk.** Low — identical force structure to today, just rotated. The single
flat patch on the cell flank may overlap the cortex shell (the patch is a plane through the box;
the cell occupies its centre) → must keep WCA cross-pairs off (already the case,
`cell.py:1086-1098`) and choose the plane offset so mother seeds don't construct *inside* the
cortex.

**Composition with FA + membrane-load.** Works, but the single front sits on ONE flank; only the
FA under that flank co-locates. Membrane-load carries over directly (still one flat normal).

### (c) Curved / spherical-cap leading edge

**Mechanism sketch.** Replace the flat WAVE plane with a **spherical-cap** WAVE shell (a patch
of a sphere of radius R_edge at the basal periphery) so the membrane the tips push is curved,
matching the real cell edge. Per-WAVE normals point radially off the cap; barbed ends grow along
local outward cap normals.

**Fidelity.** High in principle (curvature is explicit), but it conflates two things the
Bieling/Funk model deliberately separates: the **dendritic-array mechanism** (local, flat at the
filament scale — the lamellipodium IS locally flat, ~200 nm thick) and the **µm-scale cell
curvature**. The real lamellipodial leading edge is locally near-flat; the curvature that matters
for spreading is the in-plane contact-ring curvature, which approach (a) captures more cheaply.

**Code-change surface.** Large. `generate_lamellipodium_layout` (`:329`) gets cap sampling;
`WaveMembranePin` (`:379`) becomes a radial-shell pin (like `cortex/erm.py`); `membrane.py`
contact test + normal become per-tip radial (`:633-790`); the elongation tangent bookkeeping
(`:537-547`) must track curved growth. Highest test burden (new Sanity Gate geometry cases).

**CFL / stability risk.** Medium — a radial-shell pin is stiffer-feeling than a 1-D plane pin;
the ERM CFL gate (`lamellipodium.py:99-101`) must be re-derived for the cap. More
snapshot-rebuild churn as curved tips are appended.

**Composition.** Same co-location benefit as (a) if the cap sits at the basal periphery, but at
much higher cost for marginal added realism over (a).

### (d) Decouple — basal protrusion force on the cortex (LUMPED, fallback)

**Mechanism sketch.** Do NOT instantiate the Arp2/3 patch in the full cell at all. Instead model
spreading as an **outward radial protrusion force** applied to the south-cap cortex beads (a
`md.force.Custom` pushing the contact-ring beads radially outward in the basal plane), with the
force magnitude taken from the *standalone* KU-5.x lamellipodium force-velocity result. The
full-fidelity dendritic patch stays a separate KU-5.x reconstitution sim; the full cell consumes
its *output* (a protrusion stress) as a boundary condition.

**Fidelity.** Lowest — this is a **lumped/proxy** mechanism, exactly what CLAUDE.md's
architectural principle says to avoid ("does this replace a mechanistic process with a lumped
one? If yes, prefer the mechanistic alternative"). It replaces explicit barbed-end ratcheting
with a prescribed force.

**Code-change surface.** A new `md.force.Custom` over the south-cap shell tags (pattern:
`enclosed_volume`/`membrane_surface`, attached like `cell.py:1222-1229`); NO change to
`lamellipodium.py`/`membrane.py`. The lamellipodium subsystem is simply left OFF in the
full-cell config.

**CFL / stability risk.** Low (soft applied force, like the enclosed-volume term). But it
forfeits the mechanistic dendritic density / abortive-branching observables (KU-5.1/5.3) in the
integrated cell.

**Composition.** Trivial — it's just another south-cap custom force alongside FA and
enclosed-volume. Useful as a **scaffold / sanity baseline** while (a) is built, or as a
scale-bridge analogue (cf. the Layer-2 "magnitude → fine-grained via scale-bridge" pattern),
but not the mechanistic endpoint.

---

## 3. Recommendation

**Adopt Approach (a): a basal lamellar patch at the cell–substrate contact ring**, with the
WAVE/NPF reservoir relocated to the basal plane (z ≈ 0) on a ring at the contact radius, mother
barbed ends growing radially outward in the basal plane, and the membrane-load contact test made
per-WAVE radial.

**Justification.**
- **Mechanistic fidelity (CLAUDE.md hard rule).** It preserves the full fine-grained Bieling/Funk
  dendritic machinery verbatim — every barbed end, branch, and cap stays an explicit event
  (`lamellipodium.py:466,574,788`) — and only re-frames the *geometry*. This is exactly the
  "pick the mechanistic option" call; approach (d) is the lumped alternative the principle
  rejects.
- **Physical correctness.** It puts the leading edge where it actually is in a spreading cell
  (basal contact ring, z≈0, around the south cap) rather than at a box-top ceiling, and makes
  the spreading observable (`A/A₀` footprint growth) emerge from radial barbed-end advance.
- **Composes with what's already built.** The basal lamella co-locates with the FA south cap
  (`cell.py:362`) and the substrate (`cell.py:314`); the membrane-load hook (PI-approved,
  `membrane.py:217`) and the cell.py wiring (`:880-892`, `:1428-1440`) carry over in shape.
- **Cost is honest but bounded.** The real work is in three functions (`generate_lamellipodium_
  layout`, `extend_cortex_snapshot_with_lamellipodium`, `WaveMembranePin`) plus a per-WAVE-radial
  generalisation of the membrane contact test (`membrane.py:633-790`). No integrator change, no
  new particle physics, modest new Sanity-Gate cases (radial sign/sense, ring placement,
  per-WAVE normal). Larger than (b), much smaller than (c).

**Cost call-out.** The membrane-load module (`membrane.py`) currently bakes in a *single global
planar normal* in two places — `MembraneLoad._contact_set` (`:633-653`) and
`MembraneReactionForce.set_forces` (`:746-790`). Making these per-WAVE radial is the main
non-trivial edit. The load *magnitude* law (`per_filament_load*`, `:413-452`) and the Helfrich
span cap (`:391-407`) are geometry-agnostic and need no change.

**Concrete first implementation step.** Before touching `membrane.py`, do the cheap,
high-information geometry step: **add a `geometry: "basal_ring"` option (default `"flat_plane"`)
to `generate_lamellipodium_layout` (`lamellipodium.py:329`)** that, given the FA contact radius
and basal-plane z, places WAVE beads on the contact ring with per-WAVE outward-radial mother
tangents (returned in `LamellipodiumState.tangent_of`, `:436`), keeping `geometry="flat_plane"`
bit-for-bit identical to today (so `phase1_h5.yaml` and all KU-5.x tests are untouched). Write
the Sanity-Gate sign/sense test for radial growth (tips advance outward, footprint area
increases) FIRST, on the standalone `build_lamellipodium_simulation` path (`:999`), with the
membrane OFF. Only once the dry (load-free) basal ring advances correctly do we generalise the
membrane contact test (step 2) and wire it into `build_cortex_full_simulation`. This staged
order keeps each change independently testable and never breaks the existing flat-plane gate.

---

## 4. Open questions for PI

1. **Migrating vs symmetric-spreading cell.** Does the platform need a *polarised migrating*
   cell (one lamellipodium on one flank → approach (b)'s single rotated patch is the natural
   fit) or a *symmetric spreading* cell (axisymmetric contact ring → approach (a))? The MCF7
   spreading observable (`A/A₀`, the PI experimental overlay) reads as symmetric spreading, so
   (a) is the default — please confirm. (If both are eventually needed, (a)'s `geometry` switch
   can carry a `"basal_arc"` variant for a polarised front.)
2. **Contact-ring radius / basal-plane offset.** What sets the contact radius and the lamella's
   z-offset above the substrate? Candidate: derive from the FA south-cap footprint (`z_min +
   cap_depth`, `cell.py:362`) so the lamella ring = the current footprint edge. Is that the
   intended coupling, or should the lamella define the footprint and FA follow?
3. **Substrate stiffness k_sub.** The substrate ligands are currently rigidly pinned at z=0
   (`SubstrateLigandPin`, `cell.py:314`). Does active spreading need a compliant substrate
   (finite `k_sub`, as in the molecular-clutch / traction literature) for the force-velocity
   balance to be physical, or is the rigid-substrate idealisation adequate for the operating
   point? (This affects whether the radial barbed-end load is realistic.)
4. **Does basal lamella suffice for the operating point?** GATE-B / KU-3.5 is cortex+adhesion
   (§5); is the lamellipodium required to *quantitatively* hit the spreading-area operating
   point, or is it the qualitative active-spreading mechanism with FA setting the magnitude?
   (Determines whether (a) must be production-tuned or can stay a mechanism demonstration.)
5. **Membrane-load fidelity for the radial geometry.** Approve generalising the membrane
   contact test to per-WAVE radial (the only `membrane.py` physics-shape change) — the load
   magnitude law and Helfrich cap are untouched, but the contact *geometry* edit touches the
   PI-approved KU-5.x module (`membrane.py:217`).

---

## 5. Non-goals / what this does NOT block

- **KU-3.5 / GATE-B does NOT depend on this** (PI 2026-06-06). Cortical tension is a **cortex
  property**; the operating point is set by **FA-adhesion**, not the lamellipodium. The
  `mcf7_baseline.yaml` config states this inline: the GATE-B operating point "needs adhesion,
  NOT the lamellipodium (cortex tension is a cortex property; lamellipodium is a parallel
  enhancement)" (`mcf7_baseline.yaml:81-82`), and the lamellipodium block is annotated "the
  active-spreading machine" (`:87`). This design is therefore the **active-spreading enhancement
  track**, sequenced *after* / *in parallel with* the cortex-tension + FA operating-point work —
  it must not gate, and is not gated by, GATE-B.
- **Not a rebuild of the Arp2/3 mechanism.** The Bieling/Funk D2 updaters, the γ-branching
  free-NPF coupling (`lamellipodium.py:574-615`), and the membrane-load law (`membrane.py`)
  are kept; only their geometric frame changes.
- **Not an integrator change.** BAOAB stays frozen; all advance remains a Bell-Evans rate +
  soft custom-force pins (no new stiff forces, no new CFL regime beyond the existing ERM-style
  WAVE-pin gate, `lamellipodium.py:99-101`).
- **Standalone KU-5.x reconstitution is preserved.** The flat-plane geometry
  (`phase1_h5.yaml`, `build_lamellipodium_simulation`, `lamellipodium.py:999`) remains the
  in-vitro Bieling/Funk acceptance harness; the basal-ring geometry is an additive,
  default-off-for-the-standalone-path option selected only in the full-cell config.

---

## References (in-repo, grounded)

- `aleph/cell/lamellipodium.py` — WAVE plane `:227,279-283,352-358`; mother growth `−ŷ`
  `:1016,1297`; updaters `:466,574,788`; WAVE pin `:379-408`; layout `:329`; cortex-extension
  `:1141`; standalone builder `:999`; attach `:1307`.
- `aleph/cell/membrane.py` — load law `:413-452`; Helfrich cap `:391-407`; defaults from
  `Y_max`/`wave_area` `:327-331`; planar contact test + `+ŷ` normal `:557,633-653,746-790`;
  KU-5.x PI-approval `:217`.
- `aleph/cell/cell.py` — FA south-cap `:332-363,452-465`; substrate z=0 `:314`; lamellipodium
  wiring `:880-892,1202-1209`; membrane attach `:1428-1440`; cross-subsystem WCA-off
  `:1085-1098,1106-1109`; full builder signature/docstring `:689-794`.
- `aleph/cortex/cortex.py` — origin-centred shell `:542-566,659-662`; `L_box =
  box_factor·R_cell` `:438`; straight-tangent placement note `:24-26`.
- `aleph/configs/phase1_h5.yaml` — flat-plane defaults `:26-29`. `aleph/configs/
  mcf7_baseline.yaml` — active-spreading framing + GATE-B non-dependence `:81-93`.
- `aleph/docs/briefs/H5_lamellipodium.md`, `H5_GAMMA_BRANCHING_DESIGN.md` — H.5 spec + γ
  branching design.
