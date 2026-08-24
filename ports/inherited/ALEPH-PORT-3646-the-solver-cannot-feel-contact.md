# ALEPH-PORT-3646 — the solver walks filaments through each other, because nothing stops it

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3646` |
| Lane | `46143f30` Lane W2 (cortex) |
| Status | `PROPOSED` |
| Written | `2026-08-07` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **RE-DERIVATION onto the device.** `segment_steric_energy_and_forces` already exists in NumPy and is its own reference; this writes the same law as Warp kernels. Nothing is read from a provider. |
| Aleph target | `aleph/runtime/cortex_relax_kernels.py`, `aleph/runtime/resident_cortex.py` |
| Depends on | `ALEPH-PORT-3643` (the law), `-3637` (the resident state), `-3641` (the CG driver) |
| Exists because | A shell built with **zero** segment overlaps ends one relaxation with **~200** of them at ρ = 100. |

---

## 0. The measurement, on `4090-1`

Each shell counted at build and again on the coordinates the solver returned. Job 43, and the
nine-run convergence sweep of job 48:

| ρ | resolution | segment pairs @ build | **@ relaxed** | closest, relaxed | max \|F\| |
|---:|---|---:|---:|---:|---:|
| 8 | segment | **0** | 2 | 6.81e-3 µm | 0.24 pN |
| 24 | segment | **0** | 13 | 6.61e-3 µm | 0.40 pN |
| 100 | segment | **0** | **193–206** *(nine runs)* | 1.2e-3 – 5.9e-3 µm | — |
| 100 | node | 5,159 | 5,189 | 6.47e-6 µm | 7.4e+38 pN |

**The relaxation creates interpenetration.** Not as a numerical accident at one setting: across nine
independent runs at ρ = 100 varying the tolerance over three decades and the initial step over three,
the count lands between 193 and 206 every time. It is a property of the operator, not of a run.

## 1. Why: `assemble_forces` has no contact term

`ALEPH-PORT-3637` put the axial, bending and crosslink laws on the device and **left the steric law
in NumPy**. `ResidentCortexRelaxation.assemble_forces` sums three terms. The energy the line search
minimises therefore contains no excluded volume at all, so two filaments sliding through each other
is not merely permitted — it is **rewarded** whenever it lowers a crosslink or bending term, and the
CG line search will take exactly that step because nothing opposes it.

`-3637` guarded against this by *refusing a shell with steric pairs at construction* and never
looking again. That guard was written when the only law was node-pair, so it could not have caught
this; and it is a **precondition check**, not a force. A precondition cannot constrain a trajectory.

### This is why the whole `-3643`/`-3644`/`-3645` arc was a precondition and not a result

`-3644` G-6 measured the relaxed modulus changing by under 3% when 5,159 interpenetrations were
removed at construction — and `-3644`'s retraction records why that number could not have been
larger: **the solver could not tell the two shells apart on contact, because it evaluates no contact
term.** Fixing the construction while the solver is blind buys a clean `t = 0` and nothing else.

## 2. What goes on the device

Three pieces, and the second is the one with a real design decision in it.

### 2a. The kernel

`segment_steric_energy_and_forces` (`aleph/vertical/cortex_filaments.py`) is the reference: closest
approach between two segments, WCA at that point, force split to the four endpoints by the
interpolation weights, deduplicated by contact position. The Warp kernel computes the same thing per
segment pair and accumulates with atomics.

**`aleph/runtime/cortex_relax_kernels.py` must stay pure ASCII.** Warp parses kernel source and one
em-dash in a docstring is a `SyntaxError` in codegen pointing at the wrong line
(`ALEPH-PORT-3641`'s amendment records the hour that cost).

### 2b. The neighbour list, which is where the honesty is

The pair list is built on the host with `cKDTree` on segment midpoints. It cannot be rebuilt every
iteration — a 93,495-node shell runs 30,000–45,000 CG iterations, and a host-side rebuild in the hot
loop would both destroy the wall clock and **break the residency invariant** (`-3637`: *"between
`push` and `pull`, no array whose element count scales with the mesh crosses the host boundary"*).

So the list is built once and reused, and **a stale neighbour list is a wrong force, not a slow
one**. Two options and this entry does not pre-judge between them:

* **A skin.** Search at `r_c + δ` and rebuild when any node has moved more than `δ/2` since the last
  build. The rebuild trigger is a scalar reduction, which is already a supported readback.
* **A device-resident rebuild.** A uniform grid on the device, no host crossing at all.

**A is proposed first** because it is small and because the measurement below decides whether B is
needed. What may **not** happen is a fixed list with no staleness check: that is a force law that
silently stops being one, and this project has a proposal on file about defaults that are known
wrong.

### 2c. The line search

The energy readback already exists. Adding a term changes what is minimised, so
`ALEPH-PORT-3641`'s descent guard `⟨d, r⟩ > 0` and the Polak-Ribière reset are both re-tested rather
than assumed: a stiff short-ranged term is exactly what makes a line search fail, and `r⁻¹²` at
contact is the stiffest thing in this engine by many orders of magnitude.

## 3. The gates — written before the code

| # | Gate | Threshold |
|---|---|---|
| **I-1** | **The device law equals the NumPy law.** Energy and per-node force on the same configuration, CPU-Warp against `segment_steric_energy_and_forces`. | `≤ 1e-10` relative |
| **I-2** | **The shell that starts clean stays clean.** ρ = 8, 24, 100, segment-resolved, one relaxation. | **0** segment pairs at the end, where it is now 2 / 13 / ~200 |
| **I-3** | **Residency.** Mesh-sized host crossings between `push` and `pull`, including every neighbour-list rebuild. | exactly `0`, or the rebuild count and its cost reported and the invariant amended in the open |
| **I-4** | **The list is not stale.** Max node displacement since the last rebuild, against the skin `δ`. | never exceeds `δ/2`, asserted |
| **I-5** | **The solver still converges.** Nine-run tolerance × initial-step sweep at ρ = 100, as `-3645`'s W-3 run. | spread comparable to the **0.435%** measured without the term |
| **I-6** | **What excluded volume costs the modulus.** Relaxed `K_A` at ρ = 8, 24, 100 with and without the device term. | **measured, no threshold** |
| **I-7** | **The retained control.** Every committed relaxed number is still reachable with the term off. | bitwise |

**I-2 is the port's reason and can fail honestly.** If contacts still form with the force present, the
term is too soft, the list is stale, or the line search is stepping past the barrier — three different
failures with three different fixes, and which one it is would be the next measurement rather than a
guess.

**I-5 can fail honestly and is the risk that would sink this.** `r⁻¹²` is not a gentle addition to a
minimisation. If adding it makes CG stop at a few hundred iterations the way the ρ = 100
segment-resolved shell already does at one setting, then a static minimisation of a hard-core
potential needs a different method — a capped or shifted potential during the solve, or an
augmented-Lagrangian treatment of the contacts — and that is a different entry.

**I-6 carries no threshold deliberately.** This lane has now been wrong about the size of an effect
**seven** times on this project, most recently in `-3644`'s retracted amendment, where the arithmetic
was right and the attribution was not.

## 4. What is NOT claimed

- **Not that this makes the cortex correct.** It makes the solver unable to do one specific wrong
  thing. The mesh is still 3 nodes, the crosslinks still do not unbind, and there is still no thermal
  motion.
- **Not that the neighbour-list choice is settled.** §2b names two and proposes one; if the skin's
  rebuild rate at ρ = 100 is high enough to dominate the wall clock, B is the answer and this entry
  will say so with the rate in it.
- **Not that ε = 50 pN/µm is right.** `steric_contact_stiffness_pn_per_um` is row 23 of the axis
  registry, class `axis`, marked **UNSOURCED**, and porting the law to the device does not source it.
- **Not that the construction fix is superseded.** `-3644`'s segment resolution is still what gives a
  clean `t = 0`, and I-2 is only meaningful starting from one.

---

## Amendment, 2026-08-07 — the term is on the device. I-1 holds; **I-2 failed and its criterion was wrong**

CPU-Warp, ρ = 8 and 24, level 2, segment-resolved shells, filamin, strain 1e-3, CG at 1e-3 pN.

### The gates

| # | outcome |
|---|---|
| **I-1** device == NumPy | **PASS** — energy `5.4e-13` / `1.5e-13`, force `3.7e-12` / `5.5e-12` relative |
| **I-2** the clean shell stays clean | **FAILED** — 1 and 13 segment pairs, not 0 |
| **I-3** residency | **PASS** — `bulk_readbacks = 0`, and **2** pair-list builds, reported not hidden |
| **I-4** the list is not stale | **PASS** — 0.0171 / 0.0368 µm moved against a skin/2 of 0.05 |
| **I-5** the solver still converges | **PASS** — 13,171 / 22,329 iterations against 9,692 / 30,844 without |
| **I-6** what excluded volume costs | **measured** — `−0.093%` (ρ = 8), `+0.099%` (ρ = 24) |
| **I-7** the retained control | **PASS** — the term is off by default and the old refusal still fires |

### I-2 failed, and the criterion is the reason — the eighth time on this project

**A WCA equilibrium has contacts inside `r_c` by construction.** `r_c = 2^(1/6)σ = 1.1225 σ`, the
potential is purely repulsive with its zero at `r_c`, so two filaments pressed together by crosslink
tension **must** settle at some `r < r_c`. Asking for *zero pairs inside `r_c`* is asking that no
filament touch any other, which a dense shell cannot satisfy and should not.

The physics is in the **separation**, not the count, and `ALEPH-PORT-3644`'s amendment had already
named that distinction *before* this measurement — *"one is a shell that is not a shell; the other is
a shell in mild contact"*:

| ρ | steric on device | pairs inside `r_c` | **closest approach** | |
|---:|---|---:|---:|---|
| 8 | off | 2 | **0.9726 σ** | inside σ — interpenetrating |
| 8 | **on** | 1 | **1.1000 σ** | outside σ — touching |
| 24 | off | 14 | **0.9444 σ** | inside σ — interpenetrating |
| 24 | **on** | 13 | **1.0109 σ** | outside σ — touching |

**The count barely moves and the separation moves across σ.** Without the term the solver compresses
filaments 3–6% into each other; with it they rest at or beyond contact. That is exactly what the port
was for, and I-2 as written cannot see it.

**I-2 is recorded FAILED.** A gate re-run against a criterion chosen after seeing the result is not a
gate — `ALEPH-PORT-3645`'s amendment records this lane declining the same move four hours earlier and
the reasoning has not changed. The corrected criterion — *closest approach ≥ σ* — belongs in its own
entry, pre-registered, and the separations above are what a reader needs either way.

### I-6, and it is small

| ρ | `K_A` relaxed, term off | term on | change |
|---:|---:|---:|---:|
| 8 | 175.2694 | 175.1056 | **−0.093%** |
| 24 | 1436.107 | 1437.524 | **+0.099%** |

**One part in a thousand, and the sign is not consistent.** Read together with `-3644`'s corrected
G-6 — construction-level interpenetration worth `−0.84%` against a run-to-run spread of `±0.5%` —
the picture is consistent and worth stating plainly:

> **Excluded volume barely moves this cortex's areal modulus.** At ρ ≤ 24 the shell's stiffness is
> set by its crosslinks and its bending, and contact is a correction of order 0.1%. What contact
> changes is the *structure* — whether filaments are 6% inside each other or resting against each
> other — and a modulus is simply not a sensitive probe of that.

That is a result about the observable as much as about the physics, and it is the same lesson as
`-3642` §0: the areal modulus is blind to things a cortex reader cares about.

### Cost

Wall clock roughly doubles on CPU (16.1 s against 7.9 s at ρ = 8; 133.6 s against 65.5 s at ρ = 24)
for 8,161 and 73,069 skin-list pairs. **Two** pair-list rebuilds per relaxation, not per iteration,
so `ALEPH-PORT-3637`'s residency invariant is intact as written and `pair_list_builds` reports the
host reads that are outside it.

### And the driver had no test at all

`grep -rl ResidentCortexRelaxation tests/` returned **nothing** before today. The module that
produces every relaxed cortical modulus in this repository was uncontrolled, which is why a force
law could be added to it without anything going red first.
`tests/runtime/test_resident_cortex_steric.py` is seven controls, including one that fails if the
launch is missing from `assemble_forces` even though the term would still compute correctly beside
it — the shape of defect this lane has already met once, in an exporter that silently ignored a
placement and whose only tell was two runs returning bit-identical numbers.

---

## Amendment, 2026-08-07 — **I-1 FAILS above ρ ≈ 24, and the diagnosis is exact**

Job 51, `4090-1`, ρ = 100:

| ρ | device pairs | E device | E NumPy | rel ΔE | **rel ΔF** |
|---:|---:|---:|---:|---:|---:|
| 8 | 8,161 | 2.4802537530e-06 | 2.4802537530e-06 | 5.4e-13 | 3.7e-12 |
| 24 | 73,069 | 3.8751932242e-05 | 3.8751932242e-05 | 1.5e-13 | 5.5e-12 |
| **48** | 293,949 | 1.5583361381e-04 | 1.5566864697e-04 | **1.1e-03** | **4.4e-01** |
| **100** | 1,284,469 | 7.02078009e-04 | 7.02284197e-04 | **2.9e-04** | **4.3e-01** |

And at ρ = 100 the term **did nothing**: 196 segment pairs after relaxation with it on, 196 with it
off.

### It is not the card and it is not the skin

Reproduced on **CPU-Warp** at ρ = 48 with the identical numbers, so it is not a CUDA `atomic_add`
ordering effect. And running the device with `steric_skin_um = 0.0`, which makes its candidate set
*identical* to the NumPy law's (91,428 both), gives **exactly the same** disagreement — so it is not
the skin either.

### It is one contact, and here it is

At ρ = 48 there are **782** contacts inside `r_c`. The NumPy law's position-quantised deduplication
merges **one** of them; this entry's topological rule drops **zero**. That single contact is the
whole 44% force disagreement:

```
L = [8244 8245]  R = [43878 43879]   s = 0.998466   t = 0.537753   d = 1.1109 sigma
L = [8245 8246]  R = [43878 43879]   s = 0.000167   t = 0.536905   d = 1.1113 sigma
                                     contact points 0.165 nm apart, quantum 0.35 nm
```

**One crossing at node 8245, reported by both of that node's segments.** The filament is very
slightly bent, so `s` is `0.998466` rather than `1.0` and `0.000167` rather than `0.0`.

**The kernel's rule tests `s >= 1.0` exactly.** On a bent filament an exactly-endpoint contact is a
measure-zero condition, so the rule fires almost never and the shared node is pushed twice. It
worked at ρ = 8 and 24 by luck — no such near-endpoint pair happened to be in contact there.

**The NumPy rule is right and this entry's is wrong.** §2a claimed *"the two rules agree wherever
both apply"*, and that claim was checked only on exactly-coincident test configurations, which is
the one case where a bent filament does not occur.

### I-1 is recorded FAILED and the term is not fit above ρ ≈ 24

`steric_on_device` stays **off by default**, so nothing committed is affected and I-7 still holds.
What is on the record is that the device law is at parity to `1e-12` at ρ = 8 and 24 and **not at
parity at 48 and 100**, and that a relaxation at native density with it on is not currently a
measurement of anything.

### Three fixes, and the choice is a real one

1. **Widen the endpoint test** to `s >= 1 - tol`. Cheapest, and it introduces a fitted constant into
   a construction that has been param-free since `ALEPH-PORT-3636`. It would also wrongly drop
   genuine near-endpoint contacts that have no duplicate.
2. **Local ownership.** A thread recomputes its left segment's *predecessor*'s closest approach to
   the same right segment; if that contact is within the quantum, the predecessor owns it and this
   thread returns. Fully local, deterministic, ~3× the kernel work — but it keeps the *first in
   index order* where NumPy keeps the *closest*, and the two duplicates above differ in distance by
   2.4e-6 µm, so exact parity to `1e-10` is still not guaranteed.
3. **Pair the filaments, not the segments.** Build the list over *filament* pairs and let one thread
   scan that pair's segment products for the single closest approach. One crossing is then one
   contact **by construction** and no deduplication rule is needed anywhere — but two filaments that
   genuinely cross twice would get one contact instead of two.

**3 is the one that removes the class of defect rather than the instance**, and its cost is
bounded: at three nodes a filament pair is a 2 × 2 scan. Its risk is stated above and is measurable.
That is `ALEPH-PORT-3647`, with its gates written before its code, and this entry does not pre-empt
the choice further than naming which one it would defend.
