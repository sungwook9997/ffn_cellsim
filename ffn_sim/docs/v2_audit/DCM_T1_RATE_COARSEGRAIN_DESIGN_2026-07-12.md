# DCM rate/event-driven T1 coarse-graining — DESIGN (PI 2026-07-12)

## Why (from the timescale attack)
`DCM_TIMESCALE_ATTACK_2026-07-11` proved the jammed-unjamming biology-time gap is FUNDAMENTAL: tissue flow/unjamming
comes only from fast, sub-timestep T1 rearrangement events, so ANY large-dt scheme on the fine-grained mechanics
correctly freezes (BDF2 confirmed — a better integrator freezes cleaner). PI's decision: model the slow rearrangement as
a **physical RATE process** (KMC of T1 events) at large dt, layered on the fine-grained DCM mechanics — the only rigorous
route to biology-time. Matches the Kim-corpus τ_p=2^N·τ compaction-rate paradigm.

## The physics to coarse-grain — the T1 event
A T1 is the elementary unit of tissue rearrangement: two neighbouring cells (A,B) lose contact while two previously
non-neighbouring cells (C,D) gain contact (a junction flip). It is the mechanistic origin of tissue fluidity; its RATE
sets the effective tissue viscosity. In a jammed packing the T1 requires crossing an energy barrier ΔE (the cells must
deform to swap); near the unjamming shape index s0*=5.41 the barrier → 0.

## Rate law (PHYSICAL, not Metropolis — consistent with the mechanistic hard rule)
This is a coarse-grained rate of a mechanical process (like Bell-Evans for bonds, Hill for motors — NOT a
detailed-balance MC acceptance). Arrhenius barrier-crossing under active drive:

  k_T1(pair) = k0 · exp( − ΔE_barrier(local) / E_active )

- **ΔE_barrier(local)** from the local junction state — the deformation energy to execute the swap, a function of the
  local shape index / junction tension / how far below s0*=5.41. → 0 as s → s0* (barrier vanishes at unjamming). This
  is READ from the fine-grained mechanics (the virial junction stress + cell shape already computed), NOT a free
  parameter.
- **E_active** = the active energy scale that crosses the barrier: from the motility (½·drag·v0²·τ_p), junction
  turnover, or cortical-tension fluctuation — again read from the physiological inputs, not tuned.
- **k0** = attempt frequency, set by the fast mechanical relaxation time (cortex/contact), NOT fitted.

Every quantity is grounded in the fine-grained mechanics or physiological inputs (Magic-Number Block clean). The rate
law's shape (barrier → 0 at s0*) must reproduce the SPV/Bi-Manning unjamming transition as a consistency gate.

**TAG literature validation (2026-07-12).** The design + G3 grounding are corroborated by the KB:
- **KB-5.12 (T1 topological transitions):** "a T1 fires when an edge L < L_min; the energy barrier depends on p0;
  **T1 frequency is high in the fluid phase and zero in the solid phase**." This is EXACTLY the rate-law shape
  (barrier → 0 at unjamming; rate → 0 when jammed) — and it matches G3 (frozen/jammed state = 0 T1; fluid state =
  6.9×10⁻³ /cell/s).
- **KB-4.13 (junction remodeling):** endocytic junction turnover **k_endo ≈ 0.01–0.1 s⁻¹**, "T1 = junction shrinks to
  zero, neighbour swap." This is the literature anchor for **k0** (the attempt/gating frequency) — NOT a fitted number.
  The measured T1 rate 6.9×10⁻³ s⁻¹ sits just below this gate (a T1 needs a junction to remodel, so k_T1 ≤ k_endo ✓).
- **KB-5.8 (RAB5A) + motility-driven-glass:** T1 rate scales with motile force v0 (higher activity → above threshold →
  more T1s) — grounds the E_active(v0) dependence; the k_T1(s) calibration sweep (v0=10/20/40) tests this scaling.
- **s0* caveat (from TAG):** KB carries the **2D** threshold p0*=3.81 (Bi/Manning vertex); our s0*=**5.41** is the 3D
  Voronoi analog (s=S/V^{2/3}). Don't conflate — the DCM is 3D, so 5.41 is correct here.

So k0 is literature-anchored (k_endo, KB-4.13), the rate-law shape is literature-validated (KB-5.12), and the v0
dependence is literature-supported (KB-5.8) — the rate law is fully grounded, no magic numbers.

## Algorithm — KMC of T1 layered on the BDF2 fine-grained mechanics
1. Run the fine-grained DCM at large dt (BDF2 — stable, resolves the fast elastic/contact relaxation; it FREEZES the
   slow creep, which is exactly the part we now supply by rate).
2. At a rearrangement cadence Δt_rearr, for each candidate neighbour pair (cells sharing a small/high-stress junction —
   a near-T1 configuration, detected from the contact set + junction stress):
   a. compute k_T1 from the local state (above);
   b. draw the event with probability p = 1 − exp(−k_T1·Δt_rearr);
   c. if fired, EXECUTE the T1: (topology) break the A–B cadherin junction + seed the C–D junction; (geometry) apply
      the associated node displacement toward the post-swap configuration (A,B separate; C,D touch) — a discrete move,
      not a sub-dt mechanical resolution;
   d. let the BDF2 mechanics relax the post-T1 state over the next step(s).
3. The tissue now FLOWS at the physical T1 rate → biology-time unjamming/rearrangement at large, feasible dt.

## Grounding + validation (Sanity Gates — before first production)
- **G1 dimensional:** k_T1 in s⁻¹; p ∈ [0,1]; Δx of the T1 move in metres.
- **G2 barrier limit:** ΔE_barrier → 0 as s → s0*=5.41 (recovers the SPV unjamming threshold); large when deep-jammed.
- **G3 rate grounding (KEY, no magic number):** MEASURE the actual T1 rate from the small-dt fine-grained reference
  (the dt=8e-4 runs that DO resolve the creep) by tracking the cell-cell neighbour set over time and counting exchanges;
  the coarse-grained k_T1 law must reproduce THAT measured rate. The rate comes FROM the mechanics.
- **G4 convergence:** a large-dt (BDF2 + KMC-T1) run must reproduce the small-dt reference's unjamming (s→s0*, per-cell
  flow, A/A0) at the correct RATE — the whole point. Validate against the dt=8e-4 creep.
- **G5 conservation:** V/V0 = 1 through T1 moves (the swap must not create/destroy volume); no interpenetration (IPC
  re-resolves post-move).
- **G6 detailed-balance-free check:** confirm the rate is a physical forward rate (barrier + active drive), not an
  equilibrium MC — reversing time must NOT give the reverse rate for free (it's driven).

## First increments (this is a multi-step module — honest scope)
1. **G3 grounding measurement (do first):** from the committed small-dt reference frames, track the cell neighbour set
   (from the cadherin contact set / node proximity) and measure the neighbour-exchange (T1) rate vs shape index. This
   grounds the rate law + gives the k_T1(s) curve empirically.
2. T1 detection on the 3D node mesh (candidate near-T1 pairs from contact + junction stress).
3. The T1 geometric move (node displacement for a swap) + topology update (junction break/form).
4. Wire the KMC cadence into the BDF2 step loop behind a flag (`--t1-rate`); Sanity Gates G1–G6.
5. Convergence run (G4) vs the small-dt reference; then a physiological biology-time hero run.

## G3 grounding measurement — DONE (the T1 rate is real, and it IS the freeze/flow order parameter)
Tracked the cell–cell neighbour set (auto-cutoff = 1.15× median nearest-neighbour centroid distance = 17.7 µm) frame to
frame and counted exchanges (symmetric difference / 2 = T1 events) over the committed references (v0=20 µm/min, 24 s):

| reference | per-cell disp | shape index s | neighbour changes | **T1 rate** |
|---|---|---|---|---|
| BDF2 dt=8e-4 (small-dt truth) | 3.1 µm | →5.00 | 35 | **6.9×10⁻³ /cell/s** |
| operator-split dt=2e-2 | 6.7 µm | →5.31 | 227 | 45×10⁻³ /cell/s (6.5× too many) |
| BDF2 dt=2e-2 (FROZEN) | 0.16 µm | →4.93 | **0** | **0** |

**Three things this nails down:** (1) the freeze is EXACTLY zero T1 events (BDF2 2e-2 → 0 exchanges) — the mechanism of
the gap, now quantitative; (2) the correct physical T1 rate at this state is **6.9×10⁻³ /cell/s** (the small-dt
reference) — the target the coarse-grained rate must reproduce, MEASURED from the mechanics (no magic number); (3) the
T1 rate rises with shape index (frozen s=4.93 → 0; reference s=5.00 → 6.9e-3; split s=5.31 → 45e-3), empirically
tracing the k_T1(s) curve the rate law needs (barrier → 0 as s → s0*). The operator-split's 6.5× excess T1 rate is the
quantitative signature of its ~2× over-drift.

## Status
Design + Sanity Gates + G3 grounding DONE (this doc). The rate model has a MEASURED target (6.9e-3 /cell/s) and a
grounded k_T1(s) trend. Next increments (multi-step module): (2) 3D near-T1 pair detection, (3) the T1 geometric move +
junction topology update, (4) KMC cadence in the BDF2 loop behind `--t1-rate`, (5) G4 convergence vs the small-dt
reference, then a physiological biology-time hero run.
