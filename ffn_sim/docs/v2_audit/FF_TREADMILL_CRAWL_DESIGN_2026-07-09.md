# FF crawl — molecular-clutch TREADMILL design (fixed-node material treadmill)

**Date:** 2026-07-09 · **Owner:** Lead session · **Status:** DESIGN (PI asked to design before implementing)
**Prereq:** [FF_CRAWL_PROTRUSION_REBUILD_2026-07-09.md](FF_CRAWL_PROTRUSION_REBUILD_2026-07-09.md) — body-force
proxy retired (commit fa9042c); cell is COMPACT but the crawl is not yet protrusion-driven.

---

## 0. Why we are here (the finding that forced this design)

Retiring the body-force protrusion revealed that **no "front push" gives compact translocation in a FIXED
network**: smear reaction → tears (dipole, +100 µm); no reaction → net-force runaway (+251 µm); pure tip
growth → too weak (capped/sparse, ≈0). The reason is physical: a real crawl is a **molecular clutch** — actin
polymerizes at the leading edge, the network **flows backward (retrograde flow)**, depolymerizes at the rear,
and clutches **grip the flowing network** → traction → translocation. Our cortex is a *fixed* Lagrangian mesh
with *no flow*, so clutches grip a stationary material point and cannot convert protrusion into crawl traction.

## 1. Target physics — the Chan-Odde / Bangasser-Odde molecular clutch (KB-2.4, KB-PIV-4)

Three governing relations (KB-2.4, ChanOdde2008 *Science* 10.1126/science.1163595; Bangasser2013 *BiophysJ* 105:581):

1. **Retrograde actin flow** (Hill force-velocity): `v_actin = v_unloaded·(1 − F_total/(N_m·F_stall^m))`.
   Physiological **10–100 nm/s** (ChanOdde2008, Vallotton2003); ~30–80 nm/s at E_sub=5 kPa. Driven by
   leading-edge polymerization pushing the network back + myosin pulling it back.
2. **Clutch force** (linear spring): `F_i = k_int·(x_actin − x_sub,i)` — each engaged clutch is a Hookean bond
   from a flowing-actin point to a fixed substrate anchor.
3. **Clutch binding** (Bell/Pereverzev force-dependent off-rate): `dN_eng/dt = k_on·(N_tot−N_eng) − ⟨k_off(F)⟩·N_eng`.

Net translocation = protrusion (front adds) − retraction (rear removes); traction = clutches gripping the
retrograde flow. Parameter anchors already resolved in-code: `k_int=1000 pN/µm` (INTEGRIN_A5B1, Kong 2009,
KB-2.5), catch-slip `F*=7 pN` (Pereverzev; unreconciled vs Kong's ~30 pN, **kept AS-RECORDED per PI** — do NOT
retune x_c/x_s), per-clutch traction 5–20 pN (KB-2.12), barbed-end ratchet v0≈0.62 µm/s (Pollard 1986).

## 2. Representation decision — **fixed-node MATERIAL treadmill**, NOT topological remesh

The hard question: how to represent retrograde flow + front-add/rear-remove in the FF mesh. Audit verdict:

**FF has NO remesh, and a topology change is a re-build()-per-event.** All connectivity is flat integer index
arrays into a fixed-length `pos` (fiber_offsets, segments, seg_rest+seg_off, bend_triples, xl_ij(+k,+rest),
myo_ij) plus a frame-0 ConvexHull `faces` list, the frozen volume target `V0`, rest areas `_area0`, the
nucleus block offset `Ne`, and ~14 fixed-dim Warp/cupy device arrays + the implicit-CG operator. A single
front-insert/rear-remove would force: re-index every flat array, re-triangulate + **re-baseline V0/A0**,
reallocate every device buffer, and rebuild the CG operator — **breaking GPU-residency every event**.

**DCM's remesher does NOT port.** `dcm/dcm_remesh.py` (split/collapse/swap) operates on a closed *surface*
mesh where the triangle IS the mechanical element. FF's cortex mechanics live in a **fiber network** (segments
+ bend triples + crosslink node-pairs); its surface triangulation is a *separate bookkeeping hull* used only for
enclosed volume. Porting DCM split/collapse would remesh the hull, not insert/remove material in the fibers —
the cortex would not treadmill. DCM remesh is at most an optional **hull-quality safety net** (if crawl
deformation slivers the hull), never the treadmill mechanism.

**Decision (both agents concur): fixed-node MATERIAL treadmill.** Represent retrograde flow as a **material
field advected through the fixed node lattice** — nodes stay put; the actin *material* (rest-length profile +
crosslink rest + clutch grip label) flows rearward. This keeps node count, every index array, and the hull-face
volume topology FIXED — no re-index, no realloc, no CG rebuild, GPU-resident — and is *more* mechanistically
faithful (mass flux through a material, not discrete node teleportation). It also reuses the validated primitives.

## 3. The mechanism — concrete

A fiber is a chain `[a, b)`; barbed end (front-facing) = `b−1`, pointed end (rear-facing) = `a`. "Material
coordinate" = which segment/node a piece of actin currently occupies. The treadmill advects that coordinate
rearward (barbed→pointed) at `v_retro`, while the *nodes* move only under force.

| piece | what | status |
|---|---|---|
| **(a) Front polymerization** | grow leading-cap barbed-end `seg_rest` (Mogilner-Oster ratchet, capped) | ✅ HAVE `directed_front_growth_kernel` |
| **(b) Rear depolymerization** | shrink pointed-end `seg_rest` at the pointed end of rear-facing fibers (mirror of (a); pointed-end disassembly rate) | 🔨 NEW mirror kernel |
| **(c) Retrograde material advection** | shift the rest-length profile rearward along each fiber (segment k inherits k+1's rest length) at `v_retro`; advect crosslink rest + clutch grip labels the same way | 🔨 NEW advection kernel |
| **(d) Clutch↔flow coupling** | the clutch grips a **material label**, re-indexed rearward each tick; its spring pulls the *current* grip node toward the *fixed* anchor → as material flows back, the grip point falls behind the anchor → `L` grows from FLOW (not just deformation) → load builds | 🔨 NEW coupling (re-index `ac_d` rearward) |
| **(e) Front clutch nucleation** | redirect the nascent-rebind block to plant NEW clutches on freshly-polymerized FRONT material (the rejected front-bias failed because it re-formed OLD nodes that fall behind the COM) | 🔧 MODIFY existing rebind |
| **(f) Rear clutch release** | catch-slip at `F*=7 pN` — becomes rear-release automatically once `L` is flow-driven (load peaks as the grip reaches the pointed end) | ✅ HAVE `clutch_catchslip_kmc_kernel` |
| **(g) Mass conservation** | front polymerization rate ≈ rear depolymerization rate → total actin ~const → cell stays compact | 🔨 NEW balance (rate-match (a)↔(b)) |
| **(h) Fluidization** | crosslink Maxwell-relax so the cortex flows, not Ferrer-stiff | ✅ HAVE `xl_turnover_kernel` |

**Traction sign (the payoff).** As material flows rearward, the clutch's gripped node falls *behind* its fixed
anchor → `f = k_int·(anchor − actin) > 0` (forward) → the clutch pulls the network **forward** → net forward
traction on the cell body → **compact translocation** (front adds = rear removes, so the cell moves without
elongating). Clutches OFF ⇒ flow is free ⇒ no traction ⇒ ≈0 net drift (the built-in `--audit` gate).

**Volume/area re-baselining: NOT needed** for the material treadmill (unlike a topology remesh). `cortex_volume_kernel`
recomputes V every step over current node positions; because mass is conserved (a↔b balanced), total rest length
and the equilibrium shape stay ~constant, so the frame-0 `V0`/`_area0` targets remain valid. (A topology remesh
*would* need re-baselining — another reason to prefer the material treadmill.)

## 4. Validation gates (write before implementing; do not loosen)

- **G1 — retrograde flow speed:** the material-advection `v_retro` lands in **10–100 nm/s** (ChanOdde) at the
  physiological setpoint, emergent from polymerization + myosin, not tuned.
- **G2 — compact:** front-rear extent Δ ≤ ~1 µm over the run (no tear, no monotonic elongation).
- **G3 — whole-cell translocation:** front AND rear COM advance together (Δ_front ≈ Δ_rear), not front-only.
- **G4 — traction-driven:** `--audit` clutches-OFF ⇒ net COM drift ≈ 0 (the crawl comes from clutch-gripped flow).
- **G5 — mass conserved:** Σ front polymerization ≈ Σ rear depolymerization each window; total actin length ~const.
- **G6 — biphasic traction (cross-check):** traction vs substrate stiffness peaks at the Chan-Odde optimum
  (2–300 kPa band, Bangasser2013) — a *cross-check* against analytic/lit, not a tuning target.
- **G7 — seed-robust:** compact + forward across seeds 7/11/17/23 (ensemble, not single-seed).
- **G8 — NATIVE:** Nc≈266k on A5000 — compact + physiological speed + stable (no NaN/divergence).

## 5. Implementation staging

- **S1 — rear depolymerization** (`pointed_end_depoly_kernel`, mirror of `directed_front_growth`): shrink
  pointed-end `seg_rest` on rear-facing fibers; rate-matched to front growth (G5). Gate: mass conserved, compact.
- **S2 — retrograde material advection** (`retrograde_advect_kernel`): shift `seg_rest` (+ `xl_rest`) rearward
  along fibers at `v_retro`. Gate: G1 flow speed physiological, still compact (G2).
- **S3 — clutch↔flow coupling** (`ac_d` rearward re-index + front nucleation redirect): the clutch load builds
  from flow; catch-slip releases at the rear. Gate: G4 traction-driven, G3 whole-cell translocation.
- **S4 — coarse validation:** all gates G1–G7 on the coarse cell (host + small GPU). Tune ONLY physiological
  setpoints (G-actin, myosin, k_on) at their lit values — never to pass a gate.
- **S5 — NATIVE (A5000):** G8. Watch solver stability (advection + flow + catch-slip + implicit CG).

## 6. Reuse vs new (from the audit)

**REUSE unchanged:** `clutch_spring_kernel`, `clutch_catchslip_kmc_kernel` (material-agnostic — F=k_int·|L−rest|),
`fa_maturation` stack, `xl_turnover_kernel`, `reshape_kernel`, `directed_front_growth_kernel`, `resolve_polymerization`.
**NEW:** rear depolymerization kernel; retrograde material-advection kernel; clutch material-label re-index +
front-material nucleation. **MODIFY:** nascent-rebind (front-material target, not fixed basal node).

## 7. Risks + open questions (surface to PI)

1. **Bounded per-node advection.** Rest-length advection through fixed nodes is bounded (seg cap); a *sustained*
   treadmill over many cell-lengths may need periodic material **re-seeding** (reset the rest-length profile) —
   or, only if that proves insufficient, escalate to Option 2 (periodic DCM-style remesh at low cadence). Decide
   the re-seeding scheme in S2.
2. **Clutch re-index discreteness.** Walking `ac_d` rearward is a discrete jump per tick; sub-node advection
   needs either a fractional material coordinate or a fine enough tick. Watch for load-jump artifacts (G4).
3. **F* = 7 pN unreconciled** (vs Kong ~30 pN) — kept AS-RECORDED per PI; the treadmill inherits it. Flag if the
   rear-release timing looks wrong; do NOT retune x_c/x_s.
4. **Coupling stability** at large implicit dt (advection + flow + catch-slip KMC + CG). May cap dt or need
   sub-stepping the advection. Watch G8 native.
5. **Myosin's role in flow.** Retrograde flow is polymerization + myosin-driven; the crawl already carries myosin
   (Stam-Hocky). Confirm the flow (G1) is emergent from both, not from a tuned `v_retro`.

## 8. One-line summary
Represent the actin treadmill as a **material field advected rearward through fixed nodes** (front polymerize +
rear depolymerize + rest-length/clutch-label advection), and **couple the clutch to the flowing material** (grip
label re-indexed rearward) so retrograde flow builds clutch load → forward traction → **compact translocation** —
all with the node topology, index arrays, and hull-volume FIXED (GPU-resident), no remesh.
