---
kb_record:
  topic: KU-3.5-cortical-tension
  claim: KB-3.5
  status: superseded
  authoritative_as_of: 2026-05-31
  superseded_by: CORTICAL_TENSION_RECORD_2026-06-04
  aliases: [KU-3.5, KU3.5, cortical tension, binned-r0 ratchet, grip-walk]
  conclusion: >-
    Diagnosed the floor as the binned-r0 Hill ratchet (a lumped proxy that
    relabels a rest length instead of transporting actin material, so it cannot
    sustain contractile force). The force-generation root cause is CONFIRMED by
    the 2026-06-04 authoritative record, but the single-fix-closure framing
    (grip-walk alone resolves KU-3.5) was SUPERSEDED: grip-walk + --v0-accel do
    not reach the band; the wall is force aggregation per the 6/4 record.
---

# KU-3.5 cortical-tension floor — ROOT CAUSE (2026-05-31, autonomous session)

> **The months-long KU-3.5 γ floor (~3–4×10⁻⁴ mN/m, ~1000× under the [0.35, 0.65] band)
> is now fully diagnosed.** It is NOT the FA clutch, NOT myosin failing to step, NOT the
> long-run crash, NOT run length, and NOT the measurement. **It is the myosin contractile
> mechanism: the binned-r0 Hill ratchet is a LUMPED PROXY that does not transport actin
> material, so it cannot sustain contractile force.** The fix is core contractile physics
> (AFINES grip-point walking) → **PI-gated**. This doc is the decision brief.
>
> Derived by two read-only investigation agents + a crash-free 2.5e6-step diagnostic run.
> Branch `phase1/h3-cortex`; see `AUTONOMOUS_LOG_2026-05-31.md` for the full session trail.

## 0. What was ruled OUT (each verified, not assumed)

| Hypothesis | Status | Evidence |
|---|---|---|
| FA clutch never loads (blocker a) | **FIXED (G1)** | contact-footprint seeding → clutch 9→45; commits 75f464e…3d70e6b |
| Myosin never steps (blocker b) | **FIXED / was artifact (G2)** | `step_advances` 203→1396 over 2.5e6 steps; not a dead gate |
| Long-run crash blocks the measurement | **FIXED (pre-existing)** | HOOMD 7-bond/particle exclusion cap; degree-cap fix `f2f4e75`; 2.5e6 steps crash-free |
| Needs a longer run for contraction | **NO** | 5 samples over 2.5e6 steps: γ flat ~4e-4, r/r0=1.000 throughout — no trend |
| Measurement misses active stress | **NO** | method-of-planes sums k·(L−r0) over ALL bonds incl. myosin attach + the SHAKE-rigid Lagrange stress; faithful. γ reads floor because there is genuinely no strain |

**A/B (clutch on vs off):** γ 4.84e-4 (G1-on, 45 clutch) vs 4.20e-4 (G1-off, 9 clutch) — the
loaded clutch raises γ ~15%, confirming the clutch works but is NOT the floor.

## 1. ROOT CAUSE — the binned-r0 ratchet loses the contractile force

The "step" (`cortex/myosin.py:985-1022`) advances a bound head's bin by decrementing `bond_bins`,
which **relabels the head-actin bond to the next-lower-r0 bin type** (`r0 ∈ {16.5 … 313.5} nm`,
`bin_width=33 nm`, registered at `myosin.py:646-660`). A `md.bond.Harmonic` with a reduced r0
*does* pull (`F = k·(r−r0)`, k=1e-6 N/m). The mechanism is wired — **but it cannot sustain a load:**

- **(A) It transports no actin material — it relabels a rest length that relaxes out.** A real
  walking motor (AFINES `pos_a_end`, `AFINES_ALGORITHM_NOTES.md:466-565`) advances the head's
  **grip point along the actin filament**, dragging actin material a net distance each step, so the
  head spring is continually re-stretched against the load → **sustained force**. Here, advancing a
  bin only shrinks an *idealized* rest length between the head and the **same** actin bead. Nothing
  holds that bead at its old position, so each tick the head spring + the mobile actin relax toward
  the shorter r0 and the tension decays to ~0. **No strain accumulates** → γ_soft ≈ 1.4e-5 mN/m
  (~25,000× under band); even γ_rigid is on the floor because the shell is never strained.
- **(B) The single-head force ceiling is below F_stall.** r0 is floored at bin 0 = 16.5 nm
  (`np.clip(...,0,n_bins-1)`); the bin ceiling `head_actin_max_bind_dist = 330 nm` caps a single
  bond at `k·330nm = 0.33 pN < F_stall = 0.5 pN` (which needs 500 nm of stretch at k=1e-6). In
  series with the head-backbone spring the effective pulling stiffness halves to 5e-7 N/m. Heads
  hit the r0 floor or are stripped by the Bell-Evans slip bond (k_off0=10/s) before stall-loading.
- **(C) Bipolar contraction geometry is absent in effect.** Binding (`myosin.py:859-983`) attaches
  each head to its nearest bead with **no polarity/sidedness constraint** — the two opposite-polarity
  head sets are not required to grip antiparallel filaments, so even a sustained per-bond force would
  not organize into net minifilament contraction (the Stam–Hocky mechanism the brief mandates).

For scale: 555 heads × 0.33 pN (all maximally stretched, perfectly aligned) = 183 pN, vs the
~31,000 pN great-circle tension the KU-3.5 band implies — **2 orders short before any cancellation.**
The motor never generates the force; resistance is not the issue.

## 2. This violates the project's core architectural principle

CLAUDE.md (PI 2026-05-19, hard rule): *"pick the full-fidelity, fine-grained, mechanistic option
over abstracted, lumped, or proxy mechanisms."* The binned-r0 ratchet is exactly a **lumped proxy**
for the Stam–Hocky/AFINES walking motor (`briefs/H3_cortex.md:56-65` specifies AFINES §4.3 stepping).
It discards the one part that makes myosin contractile (grip-point material transport). So the KU-3.5
floor is the predicted failure mode of a sanctioned-elsewhere lumping that slipped into the runtime
contractile mechanism.

## 3. Proposed fix — TRUE grip-point walking (PI-gated; NOT implemented)

Convert the proxy into a real motor (AFINES `pos_a_end` on HOOMD):

1. **Replace bin-r0-decrement with grip-point transport.** On each step advance the head's
   **attachment bead index** along its bound actin filament toward the minus end by
   `n = floor(accum)` beads (the fractional accumulator already exists), and **re-target the
   head-actin bond to the new downstream bead while keeping the head where it is** → the bond is
   stretched by ~one bead spacing and pulls. Keep r0 ≈ 0 (idealized cross-bridge — the
   `bridge/motor.py:209` `motor_head_actin r0=0` precedent), so the *stretch*, not a shrinking r0,
   carries the force. **This single change converts the proxy into a real motor.**
2. **Restore a force ceiling that reaches F_stall** (raise k_head_actin, or let grip-advance build
   multi-bead stretch; the 330 nm bin ceiling currently caps force at 0.33 pN < F_stall).
3. **Enforce bipolar sidedness at binding** so +/− polarity head sets grip antiparallel filaments
   (net contractile dipole about the rigid rod).
4. **Keep the dt coupling rate-correct** via the existing fractional accumulator (advance in integer
   actin beads), since at the FA-limited dt the per-tick advance is sub-picometre.

This is a redesign of the runtime myosin contractile mechanism — **PI sign-off required** (it changes
core physics and the KU-3.5 mechanism). It is the load-bearing next step for KU-3.5.

## 4. State at this finding

- **Committed (production-capable, regression-green):** G1 FA contact-footprint + fan-in/exclusion
  fixes (`75f464e`,`c960406`,`d33230b`,`3d70e6b`); EXTEND modules `db1c756`; degree-cap crash fix
  `f2f4e75`. The full cell now builds + runs long (2.5e6 steps) crash-free.
- **The remaining KU-3.5 work is the §3 myosin redesign (PI-gated).** Once that lands, re-run the
  (now crash-free) v4 production to confirm γ lifts into [0.35, 0.65].

## 5. Key files
- `cortex/myosin.py` — stepping kernel `:985-1022` (relabel-not-transport), bin-r0 `:646-660`, binding `:859-983`, Bell-Evans strip `:829-854`.
- `bridge/motor.py` — `hill_velocity_clamped :134-143`; `motor_head_actin r0=0` precedent `:209`.
- `scripts/h3_ku35_tension.py` — method-of-planes `:51-115` (soft) + `:192-259` (rigid Lagrange); both correct.
- `docs/AFINES_ALGORITHM_NOTES.md:466-565` — the true walking mechanism. `docs/briefs/H3_cortex.md:56-65` — the spec the proxy departs from.
