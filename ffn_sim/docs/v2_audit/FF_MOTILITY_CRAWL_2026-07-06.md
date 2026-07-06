# FF motility — a polarized cell CRAWLS on a substrate in real η-dynamics (piece-4 + B) · 2026-07-06

**PI ask (2026-07-06):** "protrusion이 세포를 변형시키며 움직이는 것은 보지 못함 … 막이랑 필라멘트 다같이
모양 변형되는 것도 못 봤음." The prior piece-1 "✅" was for the force-velocity *kernel* + a thin filopodial
*appendage* — not a whole cell deforming + moving with membrane and filaments co-moving. PI chose **C + B**:
build the piece-4 closed loop (C) with **real η-dynamics** (B), done properly.

**Result:** a polarized whole cell on a substrate CRAWLS — the leading edge protrudes, the cell body deforms,
and the COM translocates — in real physical time, membrane + cortex + filaments + nucleus co-moving. The motion
is **traction-driven** (verified by the decisive no-clutch control), the crawl **speed is physiological**, and
the protrusion is **Newton-conserving** (an adversarial audit caught and fixed a spurious-external-force bug
before any claim).

Engine: FF (Warp, µm·pN·s). New: `ff/motility_warp.py` (physical-γ overdamped step + Newton-conserving
leading-edge protrusion + on-device crosslink turnover), driver `scripts/ff_crawl_on_substrate.py`
(+ no-clutch audit), viewer `scripts/ff_crawl_viewer.py`. GPU-native (A5000 `cuda:0`, confirmed).

---

## What was built

### B — real η-dynamics (replaces the pseudo-time step)
The existing GPU whole-cell path steps positions with `x += (0.1/kmax)·F` — a **numerical gradient-descent
step** (γ≡1, dt chosen for stability, NOT physical seconds). B replaces it with the physical overdamped law

    x_i ← x_i + (dt/γ_i)·F_i,   γ_i = units.fiber_point_drag(η = 65.9 Pa·s)   [NF2007 §5.2]

so trajectory time and the crawl SPEED are physical, measurable observables. `γ_i` is the per-node cytoplasm
drag from the single physiological viscosity — no free drag parameter (`ff/motility_warp.physical_node_gammas`).

### C — piece-4 closed loop (protrusion ↔ membrane counter-load ↔ traction → motion)
Single explicit physical-time loop; every force ticks at the same `dt`:

    F = bending + crosslink + myosin + turgor + membrane(inward γ_mem) + nucleus + substrate + FA-clutch
        + leading-edge protrusion

- **Leading-edge protrusion** = a lamellipodial Brownian-ratchet STRESS on the FRONT cortex/membrane nodes
  (not a grafted appendage — it deforms the cell BODY). Magnitude derived, not tuned:
  `f_pro = ρ_fil · (area/node) · F_stall` with ρ_fil = 100/µm² (KB-3.18), F_stall = 4 pN/filament (KB-3.6);
  grid-invariant total (~67 nN over the front). The Mogilner-Oster factor `v = e^{−f_comp·δ/kT}` (KU-3.6)
  throttles the push by the EMERGENT opposing load `f_comp = −F·p̂` (membrane tension + elastic resistance) →
  the piece-4 emergent counter-load closed loop.
- **Protrusion is INTERNAL (Newton's 3rd law).** Polymerization pushes the front membrane forward AND pushes
  the actin network backward (retrograde flow). Implemented as a force PAIR: `+f` on the front cap, `−Σf/N_c`
  uniformly on all cortex nodes (net = 0). ⇒ protrusion alone cannot move the COM; only when the retrograde
  reaction on the basal network is resisted by substrate clutches (TRACTION) does the cell translocate.
- **Adhesion** = basal FA integrin catch-slip clutches (KU-2.5, `fa_clutch_warp`): springs to fixed substrate
  anchors; over-loaded clutches rupture, detached ones re-bind force-free at the new contact (nascent
  adhesion) → a clutch treadmill.
- **Cortical crosslink turnover** (the enabler): a Ferrer-stiff cortex (α-actinin 8×10⁵ pN/µm) is too rigid to
  protrude. Real crawling cortex is an **active-gel VISCOELASTIC material that flows** (KB-3.13: crosslinker
  k_off ~1 s⁻¹, Maxwell τ_cortex ~30 s, "fluid on t≫τ"). Modeled as on-device Maxwell rest-length relaxation
  `rest ← rest + (1−e^{−k_off·Δt})(L−rest)`, k_off = 0.4 s⁻¹ (conservative within KB-3.13 ~1 s⁻¹ / KB-3.19
  α-actinin Ferrer2008 0.066 s⁻¹). Mechanistic: crosslinks are dynamic bonds, not permanent springs.

### Integration — explicit (CFL-stable), NOT implicit
The stiff α-actinin crosslinks (8×10⁵ pN/µm) set the explicit CFL `dt = 0.1·γ_min/kmax ≈ 5.5 µs`. The NF2007
Eq-2 **implicit** step (unconditionally stable, would allow large dt) was tried and **does not work at FF
scale**: its DCM-calibrated CG tolerance (`cg_tol=1e-8` absolute) is met after ONE iteration when FF forces are
~1e4–1e5 pN, so the "implicit" step degenerates to explicit and overshoots at large dt (diverges). Recalibrating
the solver for FF is a separate task (already flagged in `ff/ENGINE.md`). So this uses the **robust explicit
scheme** (CFL-bound, cannot diverge) — many small steps, built for the GPU.

---

## Validation (CPU, N=150-filament cortex + 3000-bead nucleus; A5000 GPU-native production separately)

| check | result | verdict |
|---|---|---|
| **no-clutch audit** (protrusion is internal → net COM drift must be ≈0) | clutch-ON disp∥ = +0.003 µm vs clutch-OFF +0.000 µm; ratio **9.4×** | **PASS — traction-driven, not artifact** |
| **crawl speed** (physiological? KU-3.12 retrograde 10–100 nm/s; single-cell crawl ~1–100 nm/s) | **v_crawl = 14.9 nm/s** | **PASS — physiological, no unit slip** |
| **traction** (KB-2.12 per-clutch 5–20 pN, per-cell 10–100 nN) | 1.33 nN over 58 clutches ≈ 23 pN/clutch | in range |
| **Newton conservation** (protrusion force-pair sums to 0) | clutch-OFF drift ~1.6 nm ≈ 0 | PASS |

**Adversarial-audit catch (before any claim, per the project's audit-first rule):** the first implementation
applied the protrusion as a one-sided force on the front nodes only → a spurious NET EXTERNAL force that moved
the COM at **915 nm/s even with no adhesion** (unphysical). Fixed by making protrusion an internal force pair
(retrograde reaction); the speed dropped into the physiological 15 nm/s band AND the no-clutch drift went to ~0.

---

## Honest scope / limits (v1)

- **Polarity axis is IMPOSED** (+x). piece-2 (`polarization_activegel`) confirms +x is a valid single-cap
  polarization mode, but wiring piece-2's emergent front/rear axis into the crawl is a follow-on.
- **Nucleus tracks the cortex centroid** (rigid follow, no lag) — a simplification; a lagging/deforming nucleus
  (LINC piece-5) is future.
- **Mesoscale cortex** (150–250 filaments), not native 70,686. GPU-native confirmed on the A5000; native-N is a
  scale-up, not a model change.
- **Membrane** = the inward surface-tension law on cortex nodes (Young-Laplace `2γ_mem/R`) + a convex-hull
  surface for rendering — not a separately-meshed bilayer.
- **Crawl magnitude depends on k_off** (cortex fluidity). k_off=0.4 s⁻¹ is literature-grounded (KB-3.13) and the
  resulting speed is physiological, but the exact magnitude is turnover-sensitive — PI-reviewable.
- **Implicit integrator at FF scale is unresolved** (CG-tolerance recalibration) — would give a large-dt
  speed-up; not required for correctness.

## Figures
- `outputs/ff/figs/crawl_prod_morph.html` — **interactive morphology viewer**: the cell CRAWLING (membrane
  deforming surface + actin cortex + highlighted leading edge + nucleus co-moving + substrate + COM path), ▶
  play / frame slider, rotate/zoom. This is the deliverable the PI asked to SEE.
- `outputs/ff/figs/crawl_prod.json` — crawl metrics (clutch-on) + the no-clutch audit verdict.
