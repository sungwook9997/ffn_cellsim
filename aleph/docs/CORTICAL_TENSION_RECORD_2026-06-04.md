---
kb_record:
  topic: KU-3.5-cortical-tension
  claim: KB-3.5
  gate: VG-H3-KU35-cortex-tension
  status: superseded
  primary: false
  authoritative_as_of: 2026-06-04
  superseded_by:
    - CORTICAL_TENSION_RECORD_2026-06-30
  aliases: [KU-3.5, KU3.5, cortical tension, cortical-tension floor, g_soft, gamma-floor]
  supersedes:
    - KU35_FLOOR_ROOT_CAUSE_2026-05-31
    - CORTICAL_TENSION_TRIAGE_2026-06-03
    - STAGE2-v0accel-RESOLVED@18ab034
  conclusion: >-
    The KU-3.5 active cortical-tension (g_soft) floor is a force GENERATION /
    AGGREGATION problem: per-head myosin force does not aggregate into a sustained
    shell tension (r/r0=1.0, confirmed by the 2026-06-04 compliant-backbone and
    FA-anchored ON/OFF controls). The dead --v0-accel flag WAS a real bug but does
    NOT close KU-3.5 — best case 0.030 mN/m (~11.6x under the 0.35 band) and only
    by over-driving per-head force to 2.7x stall (non-Hill, unphysical). KU-3.5
    splits into ACTIVE (g_soft, still OPEN) vs PASSIVE/composite (g_rigid native
    0.57 mN/m, in-band but myosin-INDEPENDENT) vs Layer-2 BRIDGE anchor. Cortical
    tension is measurement-protocol / timescale / adhesion-context dependent, not
    a single number. The earlier "--v0-accel resolved / KU-3.5 resolved" closure
    (commit 18ab034) is formally SUPERSEDED.
---

# Cortical-tension (KU-3.5) authoritative record + measurement reframe — 2026-06-04

> PI-ratified actions from the 2026-06-04 ultracode review (decisions 1A / 2C / 3A) +
> the H.4 re-check. Consolidates the honest state of the active cortical-tension (g_soft)
> floor after reconciling two parallel sessions' conclusions and three controlled
> experiments run this day. Companion to
> [`CORTICAL_TENSION_TRIAGE_2026-06-03.md`](CORTICAL_TENSION_TRIAGE_2026-06-03.md),
> [`STAGE2_TRANSMISSION_ELIMINATION_2026-06-03.md`](STAGE2_TRANSMISSION_ELIMINATION_2026-06-03.md),
> and [`v2_audit/_historical/MECHANISM_AUDIT_2026-05-30.md`](v2_audit/_historical/MECHANISM_AUDIT_2026-05-30.md).

## 1A — "RESOLVED" is SUPERSEDED; force-aggregation is authoritative

The STAGE-2 commit **`18ab034`** ("RESOLVED — dead --v0-accel flag was the floor") on
branch `test/mcf7-fullcell` is **formally superseded** by the force-aggregation finding
**`97f6b3e`** (Track-1 verification). The "RESOLVED" claim was over-stated and was already
internally contradicted on its own branch (the doc body says the climb "extrapolates to
~40 M steps to reach band — a real floor", and the later figure commit `0f62d47` walked it
back with "CH2 over-driving revealed").

**Authoritative state of the active-γ (g_soft) floor:**

- The dead `--v0-accel` flag **was a real bug** (myosin v0 never applied → motors acted as
  static crosslinks); wiring it lifts g_soft 17–150× off the floor. **But it does NOT reach
  the KU-3.5 band**: best case v0×3000 = **0.030 mN/m (≈11.6× under the 0.35 floor)**, and
  only by **over-driving per-head force to 2.7× stall** (the "CH2 over-driving" — non-Hill,
  unphysical). v0_accel is an un-derived tuning knob; chasing the band with it violates the
  no-magic-number / no-gate-loosening hard rules.
- The wall is **force GENERATION / AGGREGATION**: per-head myosin force does not sum into
  sustained shell tension. Confirmed across **every** lever — STAGE-2's six negatives
  (percolation, turnover, channel, coherence/meridional, buckling, soft-coupling ceiling),
  Track-1 (binding-count refuted, corr=−0.10), and the two controlled experiments below.
- The structural/passive channel **g_rigid is in-band (native 0.57 mN/m)** but is
  myosin-independent — NOT the KU-3.5 active target.

### Two controlled experiments run 2026-06-04 (both confirm force-generation, not boundary)

| Experiment | Commit | Result |
|---|---|---|
| **Compliant backbone** (`stage2_compliant_backbone.py`, motors ON/OFF control) | `48ea696` | r/r0=1.00000 **even compliant** (⇒ r/r0=1.0 is pressurized-sphere physics, NOT rigid-constraint suppression — refutes the 4a premise); **Δ_active = g_soft(ON)−g_soft(OFF) = 1.4e-4 mN/m = the same floor**. Compliant backbone does NOT lift the active floor. |
| **FA-anchored** (`h3_ku35_v4_fa.py`, n_fil=120) | run `ku35_v4_fa_ONOFF` | FA clutch **engages (86 bonds, z=0 substrate)**, yet **g_soft = 3.1e-5 mN/m (floored), r/r0=1.0**. FA grounds the *structural* channel, does NOT lift the active floor. ⚠️ Caveat: step_advances=0 / int_bound=0 (myosin not stepping at the FA-CFL dt, no v0-accel) → confirms FA doesn't help, but is not a clean active-myosin test. |

## 2C — physical-velocity acceptance filter (no super-stall g_soft)

`scripts/track1_gsoft_verify.py` now reports a **`physical_validity`** block: g_soft over
**only the Hill-bounded samples (mean F/F_stall ≤ 1)**, plus `overdriven_fraction`. A myosin
head cannot deliver more than its stall force, so any g_soft accumulated while F/F_stall > 1
is bought with non-physical (super-stall) force and is **not** a real tension. For a
literal-v0 run (F/F_stall ≈ 0.2) this is a no-op; it only bites over-driven (high-v0_accel)
runs — converting the v0-accel line from band-chasing into a falsifiable, physical readout.
(This is a measurement-protocol clarification, not a gate-contract change.)

## 3A — the KU-3.5 target is MEASUREMENT-PROTOCOL / TIMESCALE / ADHESION-CONTEXT resolved

Cortical tension is **not a single number**. The new literature (2026-06-04 batch) shows it
spans ~3 orders of magnitude by measurement protocol, force-application timescale, and
adhesion context. Report the **protocol-matched** value, not a single fixed band:

| Context / protocol | Cortical tension | Source |
|---|---|---|
| Micropipette aspiration / AFM, rounded cell (slow) | **0.03 – 3 mN/m** (KU-3.5 band [0.35,0.65] sits here) | classical; electrodeformation paper reconciliation (SE309) |
| Mechanistic-cortex model peak (coherent, intermediate filament length) | **~0.37 mN/m** (= KU-3.5 band floor) | Chugh 2017 (SE276), T₀=230 pN/µm |
| Optical-stretcher, single **epithelial** cell (MCF-10A), suspended | **~0.013 mN/m** (blebbistatin → ~0.0012) | Warmt 2021 (SE307) |
| Electrodeformation, fast pulse (0.01 s) | **~10⁻² N/m**; tp-dependent over 0.01–10 s | SE309 (all 4 platform lines) |
| Active-tension lift law | T = T₀ + linear(myosin recruitment), dissipated by turnover ~300 s | Bohec 2026 (SE312) |

**Implications:** (i) the platform's Hill-bounded plateau (~0.003–0.015 mN/m) is of the same
order as Warmt's single-epithelial γ_CST — the model may be physically correct while "under
band", because the band is a *slow-micropipette-rounded-cell* number. (ii) the measurement
context must be specified: **suspended** (FA-off, internal-pressure grounded — the canonical
band geometry, what the current model is) vs **adherent** (FA-on, substrate-coupled — the PI
experiment). Both are now runnable (FA wired via `p_fa`). The KU-3.5 band [0.35,0.65] is
retained as the **slow-timescale rounded-cell oracle**, with the timescale/context axis
flagged for any future contract revision (PI-gated).

## H.4 re-check (does FA belong in a "proper" cortical-tension implementation?)

- **Suspended-cell cortical-tension MAGNITUDE**: FA **not** required — the 2026-05-30 audit's
  "floating shell → no ground → floor" diagnosis was real but was resolved by
  **enclosed-volume internal pressure** (g_rigid → in-band), NOT FA. Canonical measurements
  are FA-free suspended cells.
- **Active g_soft floor**: FA does **not** fix it — empirically shown today (FA clutch
  engaged, g_soft floored; consistent with the STAGE-2 meridional falsification + compliant
  refutation). The floor is force-generation, upstream of any boundary.
- **Adherent-MCF7 / spreading / bridge fidelity**: FA **should** be ON — and **is already
  wired** into `Cell.build` (`p_fa`, α-restart S0/S1/S2; `h3_ku35_v4_fa.py` driver), just
  off in the suspended-cell measurement runs. This is the PI-experiment-matching axis, not
  the active-floor fix the α-restart originally bet on.

## Open / next (PI-gated for the floor itself)

The active floor remains genuinely open and is **force generation** (per-head force not
aggregating into sustained shell tension, r/r0=1.0). Candidate levers, all PI-gated core
myosin: (a) enforced bipolar antiparallel sidedness → net contractile dipole; (b) sustained
myosin recruitment > turnover (Bohec 2026), with v0 NOT the recruitment axis; (c) a clean
FA-on-vs-off test **with myosin actually stepping** (needs v0-accel wired into the v4 driver).
