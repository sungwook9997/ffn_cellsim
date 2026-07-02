# FF active movement — piece 1/5: barbed-end polymerization ratchet (2026-07-02)

**Date:** 2026-07-02  **Engine:** FF (Warp)  **Branch:** dcm/main
**PI directive:** "ff엔진 능동적 세포 움직임 구현가능한지 … full compartment." FF_ACTIVE_MOVEMENT_ASSESSMENT
ranked the missing propulsive pieces; #1 (decisive) = actin polymerization (protrusion engine). This is that
piece, built + analytic-validated. **Scope: 1 of 5** — NOT migration (pieces 2-5 below still missing).

## What this is

`ff/polymerization_warp.py` — a load-dependent barbed-end ELONGATION ratchet (Mogilner-Oster form): each step
the barbed-end segment grows by v(f)·dt, v(f)=v0·exp(−f·δ/kBT); the reshape then advances the barbed node →
protrusion. Load shares as f/N across barbed ends (`load_share`).

## Lit-anchored constants (DOIs verified; flagged for PI KB-registration)

- k_on=11.6 µM⁻¹s⁻¹, k_off=1.4 s⁻¹ (barbed, ATP-actin) — **Pollard 1986** JCB 10.1083/jcb.103.6.2747 (canonical).
- δ_full=2.7 nm; ratchet δ=1.35 nm — **Mogilner-Oster 1996** 10.1016/S0006-3495(96)79496-1.
- kBT=4.28 pN·nm @310 K; v0=(k_on·M−k_off)·δ_full → **0.62 µm/s at M=20 µM** G-actin.
- Exponential Brownian-ratchet form — Mogilner-Oster 2003 (in corpus). PI KB-registration: Pollard 1986,
  Footer 2007, Kovar 2004/2006 (SourceEvidence rows) — the corpus has only the M-O theory.

## Validation (analytic ground truth — done BEFORE any claim, per the unit-bug lessons)

`tests/ff/test_polymerization.py` (3/3): the kernel's v(f) equals the analytic v0·exp(−fδ/kBT) to rel_err
0.0e+00 over 0–15 pN; v0=0.6226 µm/s; e-fold force kBT/δ=3.17 pN (=v0/e at f=efold); load-sharing f=10/N=10 ≡
f=1; barbed segment elongates by v·dt. All the always-on analytic checks pass exactly.

## Honest scope + caveats

- **This is the protrusion FORCE-VELOCITY + elongation only.** Net cell MIGRATION still needs pieces 2-5
  (FF_ACTIVE_MOVEMENT_ASSESSMENT): (2) polarization/symmetry-breaking, (3) wired substrate FA-clutch traction,
  (4) retrograde-flow↔clutch↔protrusion closed loop, (5) LINC nucleus coupling + confinement. Not migration yet.
- **Single free filament caps ~1 pN** (Footer 2007, buckling) — real protrusion needs bundling/branching so the
  membrane load (shared f/N) stays sub-pN per filament; the kernel's load_share models this.
- **Elongation is kinematic rest-length growth** (no node-splitting topology yet) — valid for the force-velocity
  + moderate protrusion; discrete monomer-add + capping + finite G-actin pool are refinements.
- Analytic ground truth (kBT/δ, Pollard v0) is the runtime check; Footer/M-O stall values are cross-checks, not
  hard-coded (oracle-is-crosscheck rule).

## Piece-1 coupled into a cell (two demos — cortex vs the natural protrusion structure)

**(1) Cortex patch (protrusion_piece1.png).** Wired the kernel into a full-compartment cortex (turgor +
membrane, no plates), polymerizing a +x leading-edge patch (1231 barbed ends), physical-time-stepped (grow
v(f)·dt_phys, relax between). The leading edge advances **+15.5 nm with polymerization vs +0.1 nm control** over
2 s — piece-1 DOES couple. But small, because the **crosslinked cortex network + reshape (COG-conserving) absorb
most of the ~1.2 µm of segment growth** — physically correct: **the cortex is not a protrusion structure.**

**(2) Filopodium bundle (filopodium_protrusion.png + .html) — the VISIBLE protrusion.** A base-anchored parallel
actin bundle (NF=24, 30 nm hex spacing) is the natural filopodial architecture: the poly kernel grows each
barbed tip by v(f)·dt_phys and positions propagate from the FIXED pointed-end base (clean base-anchored
propagation — no COG-reshape artifact). Result: the bundle tip advances **free +1.84 µm (v=0.623=v0 exactly) vs
loaded +0.68 µm at the e-fold load 3.17 pN (v=0.229=v0/e exactly)**, both matching the analytic ∫v(f)dt to
rtol 0.02. This is a genuine µm-scale VISIBLE protrusion with the exact Mogilner-Oster force-velocity — the same
validated kernel, now in the architecture that actually protrudes. **The tip load is IMPOSED here** (0 vs the
e-fold 3.17 pN) to exercise the force-velocity cleanly; the *emergent* membrane counter-load (protrusion pushes
membrane, membrane resists, velocity settles self-consistently) is the piece-4 retrograde-flow↔clutch closed
loop — not claimed yet.

*(An earlier filopodium attempt re-anchored the base AFTER a COG-conserving reshape, injecting a strain load so
v_med≠v0; that version was correctly NOT committed. The fix is base-anchored propagation, above.)*

## Files
- `ff/polymerization_warp.py` — `resolve_polymerization`, `polymerization_kernel`, `ratchet_velocity_np`.
- `tests/ff/test_polymerization.py`.
- `outputs/ff/figs/protrusion_piece1.png` — leading-edge advance (poly vs control) in the full-compartment cortex.
- `outputs/ff/figs/filopodium_protrusion.png` — filopodium bundle: tip-advance curves + force-velocity (sim ON
  the analytic) + bundle side-view elongation (free vs loaded).
- `outputs/ff/figs/filopodium_protrusion.html` — **interactive animated** viewer (▶ play / frame slider,
  OrbitControls): the actin bundle visibly ELONGATING over the 60 polymerization steps, free vs loaded scenes.
- `scripts/ff_viewer_html.py` — extended with frame playback (▶ + slider) for animated FF line/point layers.

Related: FF_ACTIVE_MOVEMENT_ASSESSMENT_2026-07-02 (the 5-piece plan), FF_STAGE6L (unified architecture).
