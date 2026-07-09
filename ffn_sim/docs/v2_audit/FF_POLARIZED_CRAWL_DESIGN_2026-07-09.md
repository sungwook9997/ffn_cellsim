# FF crawl — POLARIZED SHAPE-CHANGE design (crawling ≠ rigid translation)

**Date:** 2026-07-09 · **Owner:** Lead session · **Status:** DESIGN (reset after PI rejected rigid-translation)
**Supersedes** the force-driven "crawl" line (body-force → anchor-drift → slip-traction), all of which the PI
correctly rejected.

---

## 0. Why every previous attempt was wrong (own it)

| attempt | what it did | why it's NOT a crawl |
|---|---|---|
| body-force (front+/rear− smear) | pushed the front out | UNBOUNDED → tore the cell into a 136 µm tube |
| anchor-drift | slid the substrate anchors forward | the FLOOR moved, cell rode along |
| slip-traction (fixed anchor) | uniform forward force on basal nodes | **RIGID TRANSLATION** — a sphere slid across the substrate, **no shape change, no polarity, MTOC static** |

All three **imposed a FORCE to move the COM**. That is not crawling. **Crawling is a SHAPE-CHANGE process** — the
front protrudes, the rear retracts, the cell is polarized, and the COM moves *as a consequence*. The PI's three
critiques (no deformation / MTOC unchanged / cut fixed at initial coords) are all one thing: **the cell was a
rigid sphere being dragged, not a crawling cell.**

## 1. What a real crawl IS (the physics we must reproduce)

A migrating cell is **polarized** and its shape **continuously changes**:
1. **Front protrusion** — actin polymerizes at the leading edge → the front cortex/membrane EXTENDS forward
   (lamellipodium/pseudopod). The front broadens/flattens. *Shape change #1.*
2. **Front adhesion** — nascent clutches form under the extended front and MATURE under load (strong, persistent).
3. **Contraction** — actomyosin contracts; with the front anchored (strong adhesions), the contraction pulls the
   cell BODY forward toward the front.
4. **Rear release + retraction** — rear adhesions release (force-gated catch-slip); the rear cortex RETRACTS
   forward (contraction pulls it) and depolymerizes → a tapered uropod. *Shape change #2.*
5. **Internal polarity** — the cell is asymmetric (broad front, tapered rear; aspect > 1 along phat); the **MTOC
   repositions toward the leading edge** and the nucleus sits rearward.
6. **Net translocation EMERGES** — front extends + rear retracts + body pulled forward. Mass conserved (front
   polymerizes = rear depolymerizes). **No force is imposed on the COM.**

**The needle to thread:** the protrusion must EXTEND the front (shape change) but stay BOUNDED (not tear) and NOT
be a uniform push (not rigid translation). Bounded by cortical tension + membrane + volume; localized to the front.

## 2. Mechanism design — orchestrate the FF pieces (no new force shortcuts)

We already have every physical piece; the failure was orchestration, not missing kernels.

**HAVE:** cortex mesh + crosslinks + reshape (structure) · myosin Stam-Hocky minifilaments (contraction) ·
clutch_spring + catch-slip + nascent-rebind + fa_maturation (adhesion) · directed_front_growth (front
polymerization) · pointed_end_depoly (rear depoly) · hard volume V=V0 · membrane inward tension · MT aster.

**The crawl cycle (all EMERGENT, no imposed COM force):**

1. **Protrude (front polymerization, bounded).** `directed_front_growth` grows the FRONT-cap forward barbed ends
   → the front cortex extends. Bounded by (a) growth cap seg_max, (b) membrane inward tension + crosslinks
   resisting, (c) volume V=V0. The front becomes a lamellipodium (extended, flatter). **Do NOT** advance nodes
   kinematically or add a body force — let the added rest length + reshape + the anchored front produce the
   extension.
2. **Adhere the front.** Nascent-rebind forms clutches **under the newly-extended leading edge** (front-biased
   RELATIVE to the current leading edge, not an absolute x); `fa_maturation` strengthens them under load →
   strong, persistent front adhesions that ANCHOR the extended front (so polymerization advances it, not retracts).
3. **Contract (the forward engine).** Myosin contracts the cortex (already on). With the FRONT anchored and the
   REAR released, the isotropic contraction pulls the cell body FORWARD toward the front adhesions — this is what
   converts front-back adhesion POLARITY into forward COM motion. **This replaces the slip-traction: the forward
   pull EMERGES from contraction + adhesion polarity, it is not injected.**
4. **Release + retract the rear.** Rear clutches, under the load of the forward pull, exceed F* and catch-slip
   release; the released rear cortex retracts forward (contraction) and `pointed_end_depoly` shrinks it → a
   tapered uropod. Mass conserved (front adds = rear removes).
5. **Reposition the MTOC/nucleus (polarity).** Couple the MTOC to the polarity: the MT aster is pulled toward the
   leading edge (MT-cortex dynein-like pull, or a polarity spring toward the front centroid), the nucleus dragged
   rearward. This gives the internal polarity the PI flagged (MTOC must NOT sit static).

**Polarity source:** the initial cue (phat) seeds front-cap polymerization + front-biased adhesion maturation;
the polarity then SELF-REINFORCES (front extends → more front adhesion → stronger forward pull → more front
extension). Emergent, not scripted per-step.

## 3. What to REMOVE / change

- **DELETE** `clutch_slip_traction` + `clutch_slip_accumulate` (the uniform forward force = rigid translation).
- **DELETE** `com_drag`/modal drag *as a crawl driver* — the COM must move because the shape changes, not because
  a soft COM mode lets an injected force slide it. (Physical per-node drag stays; the crawl speed emerges.)
- **KEEP + re-aim** directed_front_growth (front, bounded), pointed_end_depoly (rear), fa_maturation (front-mature),
  catch-slip (rear-release), myosin (contraction engine).
- **FIX the viewer**: the cut plane + camera must follow the cell (track the live COM), not sit at frame-0 coords.

## 4. Validation gates — VISUAL-FIRST (shape, not just COM)

The lesson: a COM number hid the rigid translation. Gate on SHAPE + POLARITY, on the NATIVE full cell.

- **G-shape:** the cell is POLARIZED — aspect ratio along phat > ~1.3 (broad front lamellipodium, tapered rear
  uropod); NOT a translating sphere. Measured on the native cortex, visualized.
- **G-protrude:** the front leading edge advances RELATIVE TO THE COM (the front extends, not the whole cell).
- **G-retract:** the rear advances relative to the COM faster than the front (the rear tapers/catches up).
- **G-MTOC:** the MTOC repositions toward the leading edge relative to the nucleus/COM (polarity index changes).
- **G-emergent:** NO imposed COM force — the COM displacement equals the integral of the (adhesion) external
  forces only; turning off myosin OR front polymerization kills the crawl (it is not a free slide).
- **G-traction:** clutch-OFF → no sustained translocation (adhesion-dependent).
- **G-compact/mass:** no tearing (extent bounded), mass conserved (Σ front growth ≈ Σ rear depoly).
- **G-native:** Nc≈266k, all compartments, stable — and the shape change is visible in the interactive HTML
  (camera tracking the cell, cut plane following it).

## 5. Implementation staging (native + full every step, per the HARD rule)

- **S1 — polarity scaffold + viewer fix.** Front-biased adhesion maturation (relative to live leading edge) + MTOC
  polarity coupling; fix the viewer camera/cut to track the COM. Gate: MTOC repositions, viewer follows.
- **S2 — protrusion + contraction engine.** Bounded front polymerization + myosin contraction toward the anchored
  front. Gate: G-protrude + the body is pulled forward by contraction (kill-myosin → no pull).
- **S3 — rear release + retraction.** Force-gated rear catch-slip + pointed-end depoly. Gate: G-retract, G-shape
  (polarized), mass conserved.
- **S4 — emergent translocation, native.** Remove all injected COM forces; confirm the crawl EMERGES (G-emergent,
  G-traction) with a POLARIZED shape (G-shape) at Nc≈266k, visualized with a tracking camera. This is the real
  deliverable.
- **S5 — ECM Mikado** (compliant + remodeling) once the shape-change crawl is real on the rigid substrate.

## 6. Honest challenges (where this is hard / may need iteration)

1. **Protrusion vs tearing vs rigid-translation** — the front must extend a BOUNDED lamellipodium. Too strong →
   tears (my first failure); a uniform push → rigid slide (my last failure). The balance (polymerization rate vs
   cortical tension vs volume) is the crux and will need tuning at PHYSIOLOGICAL values (not tuned-to-pass).
2. **Contraction-driven forward pull** needs a real front-back adhesion asymmetry. An earlier front-bias rebind
   was rejected because it biased an ABSOLUTE x and fell behind the advancing COM — here it must be RELATIVE to
   the live leading edge, and re-derived each turnover.
3. **MCF7 is a poor migrator** (epithelial). The physiological crawl mode may be slow/bleb-based, not a fast
   keratocyte lamellipodium — the shape change may be subtle. Set the mode/params to MCF7 physiology; surface to
   PI if a value is unknown.
4. **MTOC repositioning** requires coupling the currently-decoupled MT aster to the polarity — a real mechanism
   (dynein cortical pull) or a minimal polarity spring; pick the mechanistic one.
5. **Everything native + full** (HARD rule) — slower iteration, but coarse/stripped manufacture artifacts.

## 7. One-line summary
Stop imposing a force to move the COM. Build the **polarized shape-change cycle** — bounded front protrusion +
front adhesion maturation + myosin contraction toward the anchored front + rear release/retraction + MTOC
repositioning — so the cell **deforms into a polarized crawler** and the COM translocates *as a consequence*,
validated by SHAPE/polarity gates on the native full cell with a tracking camera.
