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

# FF motility — adhesion/contact FOUNDATION for a substrate-migrating cell (2026-07-06)

**Goal (PI 2026-07-06):** a cell that behaves like a real one — **adheres, generates traction, and MOVES**.
Not "flattening" for its own sake; realism is the target. This doc is the **foundation** toward that:
mechanistic, **no magic numbers**, verified incrementally. It also records what the PI's skeptical review
caught and what was corrected — because several first-pass claims did NOT hold up.

Engine: FF (Warp, µm·pN·s). Code: `ff/motility_warp.py` (kernels), `scripts/ff_crawl_on_substrate.py`
(driver: `--static` / `--spread` / crawl, `--mature` FA, `--audit`), `scripts/ff_crawl_viewer.py` (clean
morphology viewer — nested surfaces, not point clouds).

---

## VERIFIED (mechanistic, no magic numbers)

| piece | mechanism | verification |
|---|---|---|
| **Volume conservation** | incompressible cytoplasm: **per-step, on-device** osmotic ΔP (Guo closure, physiological Π_in0) computed from the enclosed volume every step (`cortex_volume_kernel`, oriented-face divergence theorem — matches ConvexHull to 1e-4) | static **V/V0 = 1.003**; under protrusion **V/V0 = 1.026** (was **ballooning to 4.06** before the fix). No balloon, no collapse. |
| **Gravity − buoyancy** | real body force `−Δρ·g·v_node` (Δρ = 55 kg/m³, ≈ 1 pN/cell), ported from DCM `gravity_body_force_kernel` — DERIVED, not tuned | rests the cell on the dish |
| **Stable adhesion** | basal integrin FA clutches; mature (stable) FA hold, nascent catch-slip cycle | **bound = 1.0**, basal flat contact (gap 0.023 µm), STABLE ADHERED |
| **Emergent spreading** | peripheral (basal-rim) actin-polymerization ratchet (`spreading_push_kernel`, Newton-conserving) + basal clutches — the cell widens its footprint from actin + adhesion, **NO wetting energy** | contact radius **3.1 → 5.2 µm** at constant volume |
| **Traction** | clutch springs to substrate anchors, load from cortex/protrusion | **1.25 nN**, and the no-clutch audit confirms motion is **traction-driven** (Newton-conserving protrusion) |
| **Real η-dynamics (B)** | `x += (dt/γ)·F`, γ = `units.fiber_point_drag` (η = 65.9 Pa·s) — physical time, not pseudo-time | crawl speed is a measurable observable (~15–21 nm/s, physiological) |

## What the PI's review CAUGHT and what was corrected (honest record)

- **Ballooning (PI-caught):** the first "crawl works" claim ran on a foundation whose volume was **not
  conserved** (one-sided turgor) — over ~16 s it inflated to V/V0 = 4. Root cause: no inward restoring when
  V>V0, and the membrane K_A was deferred (a physiological-baseline violation). **Fixed** by the per-step
  bidirectional osmotic volume constraint above.
- **Wetting REJECTED (PI-caught):** reusing DCM's substrate-wetting to spread the cell = **magic-number
  tuning** (a dialed adhesion-energy `W_cs` proxy) — it violates the fine-grained-over-lumped HARD rule.
  Wetting is a DCM (coarse-mesh) mechanism; the FF (fine-grained) engine must get spreading **emergently**
  from actin polymerization + clutches, which is what `spreading_push_kernel` does.
- **Contact was not proper (PI-caught):** early runs showed a floating/inflating ball. Fixed by gravity
  (settling) + volume conservation (no inflation) + stable adhesion.
- **Adversarial-audit catch (self):** the protrusion was first a one-sided force on the front nodes → a
  spurious NET external force that moved the COM even with no adhesion (915 nm/s, unphysical). Fixed to an
  internal force pair (retrograde reaction); speed dropped into the physiological band and the no-clutch
  drift went to ~0.

## HONEST limits — what is NOT done yet

- **The cell does not fully FLATTEN into a spread pancake.** A flat disk needs ~2× the surface area of the
  equal-volume sphere; the FF cortex area is ~fixed (segment inextensibility + limited crosslink turnover), so
  it cannot flatten without **membrane/cortex AREA GROWTH** (reservoir unfold + actin assembly) — not yet
  modelled. The cell stays domed with a widened base.
- **No SUSTAINED migration yet.** Traction + initial movement are shown, but translocating µm-scale over a
  real migration (minutes) needs (a) a working clutch **treadmill** (front engage / rear release) and (b) a
  **timescale** unreachable by explicit stepping (dt = 5.5 µs, set by the stiff α-actinin crosslinks).

## NEXT — FF-specific implicit solver (the real unlock; verified diagnosis)

Sustained migration needs large physical dt. **DCM's implicit step does NOT transfer to FF scale** — verified:
its matrix-free finite-difference-Jacobian CG **blows up at large dt regardless of the FD step ε** (1e-9→1e-2
all identical: CG bails at 1 iteration, |x|→1.8e6). The FD JVP is unreliable at FF force scales (~1e6 pN). The
FF-implicit must use the **NF2007 approach — assemble the ANALYTIC elastic stiffness** (banded bending 4th-diff
`[−1,4,−6,4,−1]` + crosslink-spring + volume Jacobian) and solve — the "key contribution" / Stage-6b piece
ENGINE.md already flags. A substantial, careful numerical build; do it as a focused effort, not a calibration.

## Figures
- `outputs/ff/figs/found_spread_morph.html` — interactive morphology (membrane surface + nucleus + substrate,
  ▶ play): the cell adhered on the substrate with volume conserved, footprint spreading 3.1→5.2 µm. This is the
  verified foundation state (a domed adherent cell, not yet a flat pancake, not yet sustained migration).
