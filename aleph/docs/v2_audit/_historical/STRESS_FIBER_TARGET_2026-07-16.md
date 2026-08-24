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

# Stress-fiber architecture — the target (PI reference images, 2026-07-16)

PI-supplied reference (3 figures) fixing the SF architecture to build AFTER the membrane native gate.
This is the DESIGN TARGET, not yet built. It reframes the SF gap surfaced this session.

## 1. The SF taxonomy (must be distinct, not one lumped "bundle")
- **Ventral stress fibers** — contractile actomyosin bundle with a **focal adhesion at BOTH ends** (FA→SF→FA).
  The canonical traction generator (Kumar 2006 10–30 nN). α-actinin Z-bodies periodic + NMII bands.
- **Dorsal stress fibers** — **FA at one end**, rise dorsally, terminate in a transverse arc (no FA at the far
  end). Formin-nucleated, initially little myosin (assembles myosin as it matures).
- **Transverse arcs** — curved actomyosin bundles PARALLEL to the leading edge, **not FA-anchored** directly;
  connected to dorsal SFs, flow retrogradely. Myosin + α-actinin.
- **Perinuclear actin cap** — ventral-SF-like bundles arcing OVER the nucleus, connected to the nuclear
  envelope by **LINC (nesprin/SUN)**. ⭐ This is a DIRECT nucleus deformation load path — the same physics I
  chased with the IF cage. The actin cap (not just IF) is a dominant apical nucleus-shaping structure.

## 2. SF FORMATION is EMERGENT, not scripted (Nature-figure b/c/d)
Mechanistic self-organization (the CLAUDE.md fine-grained rule demands this, NOT hand-placed bundles):
  myosin turnover → LOCAL myosin-density increase → myosin motor activity ALIGNS filaments into a bundle →
  new filaments nucleated + motors recruited → a stress fiber. (Signalling can locally template it, panel c.)
So the correct build is: seed the actomyosin network (filaments + NMII + α-actinin + FA clutches) and let the
SF EMERGE from myosin-driven alignment + FA anchoring — do NOT paint a pre-made sarcomeric bundle and call it
a stress fiber (that would be the lumped/scripted anti-pattern).

## 3. What this reframes (the SF gap, honest)
- Current `architecture_spec.STRESS_FIBER` is a PASSIVE backbone (NMII deferred, h7_stress_fibers_activation_
  gate.py) + PI-gated magic numbers (sarcomere_um / N_fil / μ_SF). The ACTIVE contractility + Kumar 10–30 nN
  gate is DEFERRED; the H.7 SF-array traction '+131pN' did NOT reproduce across seeds (HALTED).
- FA/traction was studied via imposed retrograde flow (ff_clutch_traction) or the CORTEX as substrate
  (h7_manifold_traction) — NOT contractile ventral SF. So "FA without functional SF" is a fair description.
- The actin cap connection means the SF line and the nucleus-coupling line (IF cage / re-frame) are the SAME
  problem from two structures — the apical actin cap + LINC is likely the dominant nucleus-flatten load path
  that the (polar-radial, buckling) IF cage failed to be.

## 4. Build order (post-membrane-gate)
1. Ventral SF: FA→contractile-NMII-bundle→FA, EMERGENT from the actomyosin network (not scripted), on the
   basal footprint of the `--from-resting` adherent cell. Active gate = Kumar 2006 10–30 nN single-SF tension.
2. Perinuclear actin cap: ventral-SF-like bundles over the nucleus, LINC-tethered → the direct nucleus load
   path (re-test the strain>0.34 nucleus-flatten question with the cap, not the IF cage).
3. Dorsal SF + transverse arcs: the leading-edge organization (lower priority; needs the lamellipodium).
PI-gated magic numbers (sarcomere_um, N_fil, μ_SF, cap bundle count) surfaced, not tuned.

## Change log
- 2026-07-16: created from PI reference images while the membrane native gate is in flight. Design target only.
