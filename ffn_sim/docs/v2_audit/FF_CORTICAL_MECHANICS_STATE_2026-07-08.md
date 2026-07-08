# FF cortical mechanics — consolidated state (2026-07-08)

**Purpose.** One reference reconciling the scattered pieces after the 2026-07-08 AFM investigation, and correcting
a stale framing. Supersedes the "γ-floor is the open gap / needs an engaged-NMII density datum" narrative.

## The one-paragraph state

A whole-cell mechanical measurement resolves into **two distinct quantities**, and the model now handles both:
**cortical TENSION** (a surface quantity, γ [mN/m]) and the **apparent MODULUS** (a bulk AFM quantity, E [Pa]).
The **tension was resolved 2026-07-01 (FF_STAGE6V): it is osmotic-pressure-borne** (γ=ΔP·R/2, ΔP set by
osmoregulation), reproducing both interphase (0.15≈0.17) and mitotic (1.5≈1.6) states. The **modulus was resolved
2026-07-08 (this session): the apparent "too stiff" was three compounding artifacts**, and with the biphasic
cytoplasm + physical press-speed + η-calibrated drag the model lands in the MCF7 band (~249 Pa) at a physical
sub-µm/s press. **Neither quantity is a genuine "too soft / too stiff" gap.** The real open lever is the FF network
**contractile-transmission (screening)** — a physics-vs-artifact question, not a missing datum.

## Quantity 1 — cortical TENSION γ (RESOLVED, FF_STAGE6V/6T/6Z, 2026-07-01/02)

- Cortical tension = Young-Laplace of the osmotically-pressurised cell: **γ = ΔP·R/2**, with ΔP set by
  **osmoregulation** (ion/water flux), NOT by myosin material contraction. Volume regulation (`pressure_setpoint`)
  reproduces both states + Fischer-Friedrich confinement-independence:

  | state | ΔP [Pa] | model γ [mN/m] | measured [mN/m] |
  |---|---|---|---|
  | interphase | 40 | 0.150 | ~0.17 |
  | mitotic | 400 | 1.503 | ~1.6 |

- **Myosin = container/modulator, not generator**: it bears ΔP via Laplace (blebbistatin → cortex fails →
  tension −75%), and cannot itself compress the ~7e5 Pa osmotic modulus to drive the ΔP ramp. γ_myo (direct
  material tension) ≈ 6e-5–4e-4 mN/m — the old "~530× γ-floor" is REFUTE **for γ_myo**, but γ_myo ≠ cortical tension.
- **Myosin protein is present** (FF_STAGE6T proteomics, Itzhak/Hein MYH9 copy# → 10–90 minifil/µm², spans the band;
  Nie's 0.6/µm² imaging undercounts). The "missing engaged-NMII density datum" framing is **superseded**.
- Force scale: FF_STAGE6Z — the ~10³× plate force was a soft-plate artifact; a **rigid plate gives realistic 3–39 nN**.

## Quantity 2 — apparent MODULUS E (this session, 2026-07-08)

The AFM "cell too stiff vs 249 Pa" (reported 857× → 16.5×) was **three compounding artifacts**, none real:

1. **Category error** — 857× divided a large-strain nominal σ=F/πR² by a small-strain fixed-E Hertz; matched Sneddon
   inversion → ~16.5×.
2. **Under-relaxation transient** — the solve was not converged; the small-strain modulus drops monotonically with
   relaxation (NATIVE Nc=266k: 500→64000 steps = 21900× → 0.73×, through the MCF7 band to the soft side). The
   reported "stiffness" read the force before the cell finished deforming (PI's "누르는 속도" insight).
3. **Press-speed / uncalibrated drag** — instant strain + a numerical (non-physical) dt. With physical η=65.9 Pa·s
   drag, E ∝ v_press = a VISCOUS transient. Bulk-η-calibrated drag γ_node=6πηR/Nc (canonical Stokes, no tuning) →
   η_eff validated ~2× of 65.9 (Newtonian); the model then produces a rate-dependent modulus that **passes through
   the MCF7 band at ~0.5–1 µm/s** (viscous-stiff above, tension-soft equilibrium below).

New mechanisms added (all additive, off→bit-identical, 12 sanity gates; `ff/network_warp.py` + `common/compartments.py`):
biphasic poroelastic cytoplasm (Kedem-Katchalsky drainage + Terzaghi drained solid), membrane reservoir K_A upturn,
Bell-slip crosslink turnover, physical press-speed ramp, bulk-η drag, Hertz harness. Figures:
`outputs/ff/figs/{ff_drag_calibrated,ff_native_convergence,ff_press_speed,ff_poroelastic_validation}.png`.

## Reconciliation of the two quantities

The "13 Pa too soft" I kept invoking this session **is** the pure tension-governed relaxed modulus (≈ the correct
0.15 mN/m interphase tension → ~13–20 Pa apparent), NOT a deficit. The full AFM modulus (with the cytoplasm drained
solid + physical press this session added) is what reaches the MCF7 band. **Tension (Quantity 1) and modulus
(Quantity 2) are separately correct; the "γ-floor / too soft" narrative conflated them.**

## The real open levers (NOT a datum)

1. **⭐ Network contractile-transmission (screening), ~55–74×** (FF_STAGE6H/6T) — is the FF connected-inextensible
   network's screening of contraction (σ/σ_dipole ≈ 0.2–0.6) PHYSICAL (real cortex screens too) or an FF transmission
   artifact (real cortex transmits better)? **Validate FF against a reconstituted actomyosin network of KNOWN density**
   (Murrell/Gardel/Koenderink; Linsmeier 2016 R≈0.085) — does a known-density FF network generate the measured
   reconstituted stress, or less? Non-circular. **This is the next task (PI-directed (a)).**
2. Cortex still ~6–10× over-stiff at high strain (FF_STAGE6Z residual) — not a clean liquid drop.
3. Osmoregulation as an explicit dynamic driver (finite Lp + ion-pump ΔP setpoint) — currently an imposed setpoint.
4. No direct MCF7 cortical-TENSION datum exists (only modulus, Zbiral 249 Pa) — the tension target is poorly grounded.

## Detailed source docs

FF_STAGE6V (tension pressure-borne) · FF_STAGE6T (proteomics density) · FF_STAGE6Z (rigid plate) · FF_STAGE6S
(engagement + density litsweep) · FF_STAGE6H (buckling network / screening) · FF_RESULTS_LOG (this session's AFM
modulus work, dated entries 2026-07-07/08). Memory: [[project-ff-afm-overshoot-was-artifact]] (with the correction),
[[project-gamma-floor-likely-deficit]] (⚠️ pre-resolution, superseded by FF_STAGE6V).
