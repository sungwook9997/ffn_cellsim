# ALEPH-PORT-3644 — the builder cannot see the overlaps it is supposed to remove

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3644` |
| Lane | `46143f30` Lane W2 (cortex) |
| Status | `PROPOSED` |
| Written | `2026-08-06` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **INTERNAL FIX.** The segment geometry `ALEPH-PORT-3643` already landed, applied to the construction instead of the force. Nothing is read from a provider. |
| Aleph target | `aleph/vertical/cortex_surface_coupling.py` |
| Depends on | `ALEPH-PORT-3643` (`segment_closest_approach`, `refresh_steric_segment_pairs`) |
| Blocks | re-measuring every cortex modulus; the node count becoming a parameter |
| Exists because | `ALEPH-PORT-3643` gate F-7: a shell the resolution graded **clean** carries **286** interpenetrations at ρ = 24. |

---

## 0. The measurement that forces this

`ALEPH-PORT-3643` amendment, F-7, on shells `resolve_cross_filament_overlaps` ran to completion:

| ρ | filaments | node pairs inside `r_c` | **segment pairs inside `r_c`** | \|F\| node | \|F\| segment |
|---:|---:|---:|---:|---:|---:|
| 8 | 2,493 | **0** | **17** | 0 | 3.56e+17 |
| 24 | 7,480 | **0** | **286** | 0 | 7.71e+32 |

`r⁻¹²` reaching 1e32 pN means the closest approach is a small fraction of σ. These filaments are not
brushing past each other; they are **through** each other, and the builder returned the shell as
settled.

## 1. Why: the resolution and the law shared one blind spot

`resolve_cross_filament_overlaps` sweeps `cKDTree(moved).query_pairs(cutoff)` over **nodes** and
pushes node *centres* apart. `ALEPH-PORT-3643` §1 already established what that misses — at a 150 nm
node spacing and σ = 7 nm the filament is a dotted line, and a crossing that lands between nodes is
invisible.

The consequence is worse than a missed force, because the resolution's **stopping criterion is the
same query**:

```
loop:  pairs = node pairs inside r_c
       if none: return "settled"        <- and "none" was never "no overlap"
```

So the shell it produces passes its own gate **by construction**. Every extra sweep drives the node
pairs to zero and does nothing to a crossing between nodes, and the loop then reports success. A
test that shares the defect of the thing it tests cannot fail, and this one did not.

**This is the same defect at a second site, not a new one.** It is a separate entry because the fix
is a separate function and because F-7 measured it after 3643's gates were written.

## 2. The fix

`segment_closest_approach` already returns, for a segment pair, the parameters `s, t ∈ [0,1]` of the
closest approach and the distance. The resolution becomes: for every cross-filament **segment** pair
closer than `r_c`, push the two closest-approach *points* apart to `1.01·r_c` and distribute each
point's half of the push to its two endpoints by the interpolation weights `(1−s, s)` and `(1−t, t)`.

Three properties follow, and the third is why this is the right shape rather than a bigger hammer:

* **It is a strict generalisation.** When the closest approach falls on an endpoint (`s = 0`), the
  full push lands on that node and the rule reduces exactly to the current one. Nothing that the
  node-pair resolution separates correctly is separated differently.
* **The target stays param-free.** `1.01·r_c`, `r_c = 2^(1/6)·σ`, σ is the sourced F-actin diameter.
  No rate, no threshold, no fitted constant — the same property that made `ALEPH-PORT-3636` a
  construction rather than a tuning, preserved.
* **The stopping criterion stops sharing the defect.** Convergence is tested on the **segment** count,
  which is the quantity F-7 measured. A shell that settles under this rule has zero *segment* pairs,
  which is strictly stronger than what any shell in this repository has been certified for.

### A cost this does not hide

Pushing four endpoints apart **stretches the segments**. The node-pair rule has exactly the same
property and it has never been measured; the extension does not make it worse, but it does make it
worth measuring, so G-5 measures it. If the construction injects large stretch energy, the relaxation
inherits a starting state further from its minimum, and the honest place for that is a number.

## 3. The re-measurement, and why it belongs to this entry

**Every cortical modulus this lane has published was measured on an interpenetrating shell.** How
much that is worth is not established by F-7 — a 1e32 pN force says how *close* two rods are, not how
much stress the structure carries once they are apart.

The only way to know is to build the same shell both ways and compare, which is G-6, and it is the
number the whole segment-steric arc has been heading for. It carries **no threshold**: the delta is
the result.

## 4. The default, stated in advance so the switch is not silent

`StericLaw.SEGMENT_PAIR` shipped **opt-in**, so `ALEPH-PORT-3643` moved no committed number. This
entry proposes the opposite for the resolution, and says why before the measurement rather than
after:

> **If G-1..G-5 pass, `OverlapResolution.SEGMENT_PAIR` becomes the default**, with G-6's delta in the
> commit message. The node-pair resolution is not a slower approximation of the right answer; it is
> *known to certify shells that interpenetrate*. Leaving a known-wrong construction as the default
> because changing it moves numbers is how a defect becomes a baseline.

`OverlapResolution.NODE_PAIR` stays reachable, because every published cortex number was produced
under it and reproducing them has to remain possible.

## 5. The gates — written before the code

| # | Gate | Threshold |
|---|---|---|
| **G-1** | **It removes what F-7 found.** The ρ = 8 and ρ = 24 shells, resolved under the new rule, re-counted with `cross_filament_segment_contacts`. | **0** segment pairs, where the node rule left 17 and 286 |
| **G-2** | **The node count stays the same and so does the ownership.** Resolution moves nodes; it may not add, drop or re-assign one. | shapes and `filament_of_node` identical |
| **G-3** | **It reduces to the current rule when the closest approach is at an endpoint.** Two filaments crossing at a shared node position. | displacements agree to `1e-12` µm |
| **G-4** | **It terminates at native density.** ρ = 100, 31,165 filaments, 93,495 nodes. | settles within `max_sweeps`, sweep count reported |
| **G-5** | **How much does the construction stretch the filaments?** Segment length before and after, max and RMS. | **measured, no threshold** |
| **G-6** | **What did the invisible interpenetration cost?** `K_A` affine and relaxed on the same seed and density, node-resolved vs segment-resolved, at ρ = 8 and ρ = 24. | **measured, no threshold** |
| **G-7** | **The retained control survives.** With `OverlapResolution.NODE_PAIR` requested, `K_A` = 1081.819088 at ρ = 0.5, level 2. | bitwise |

**G-1 is the port's reason and can fail honestly.** Pushing two segments apart at their crossing can
pivot them into a third, and a rule that oscillates is a real possible outcome at native density —
two nearly-parallel filaments have a shallow closest-approach minimum and the push direction is
poorly conditioned there. If it does not settle, **that is a result and gets reported**, and the
answer would be a different formulation (a length-preserving projection, or resolution during
placement rather than after it) rather than a raised sweep cap.

**G-6 carries no threshold deliberately.** This lane has written a threshold against a guess **six**
times on 2026-08-06 — `ALEPH-PORT-3641`'s amendment records the fifth — and every one was wrong. The
size of this effect is exactly the kind of thing that has been mispredicted, in both directions:
interpenetrating filaments could be *stiffening* the shell (a 1e32 pN contact is a very strong
constraint if the relaxation ever saw it) or *softening* it (crossings that pass through carry no
load at all). Both stories are plausible and neither is evidence.

## 6. What is NOT claimed

- **Not that this makes the shell physical.** It makes it non-interpenetrating. Filament orientation
  statistics, length distribution and the absence of a real polymerisation history are untouched.
- **Not that the published moduli are wrong by a known amount** — until G-6 runs, they are wrong by
  an unmeasured amount, which is a different and more honest statement.
- **Not that the segment steric law is on by default.** `ALEPH-PORT-3643` left the *force* opt-in and
  this entry does not change that. A shell that is genuinely non-interpenetrating is the
  precondition for switching it, not a consequence.
- **Not that the node count is addressed.** Still unblocked, still unchanged. Refining the mesh under
  a segment-aware resolution is now safe in a way it was not this morning.
- **Not a performance claim.** The segment pair search is `O(candidates)` per sweep with a search
  radius of `r_c + longest segment`, which at native density is dominated by the segment length, not
  by σ. G-4 reports the sweep count and the wall clock; neither is a gate.

---

## Amendment, 2026-08-07 — every gate, and the default did **not** move

### The gates

| # | outcome |
|---|---|
| **G-1** removes what F-7 found | **PASS** — 17 → **0** (ρ = 8) and 286 → **0** (ρ = 24), reproducing F-7's counts exactly on the node rule first |
| **G-2** node count and ownership untouched | **PASS** — `(7479, 3)` and `(22440, 3)`, `filament_of_node` identical |
| **G-3** reduces to the node rule at an endpoint | **PASS** — max difference **exactly `0`**, both rules settling in 2 sweeps |
| **G-4** terminates at native density | **PASS** — ρ = 100, 31,165 filaments, 93,495 nodes: **14 sweeps, 0 segment pairs, 1.7 s** |
| **G-5** how much it stretches the filaments | **measured** — rms `3.3e-4`, max `2.5e-2` of rest length at ρ = 100 |
| **G-6** what the interpenetration cost | **measured, and the answer is "nothing the affine estimator can see"** |
| **G-7** the retained control | **PASS** — `K_A` = 1081.819088 |

### G-1 and G-4, the two that mattered

| ρ | filaments | rule | sweeps | node pairs | **segment pairs** | max \|F\| under the segment law |
|---:|---:|---|---:|---:|---:|---:|
| 8 | 2,493 | node | 2 | 0 | **17** | 4.87e+17 pN |
| 8 | 2,493 | **segment** | 2 | 0 | **0** | **0** |
| 24 | 7,480 | node | 2 | 0 | **286** | 7.82e+32 pN |
| 24 | 7,480 | **segment** | 6 | 0 | **0** | **0** |
| 48 | 14,959 | node | 2 | 0 | **1,200** | 3.32e+34 pN |
| 48 | 14,959 | **segment** | 14 | 0 | **0** | **0** |
| **100** | **31,165** | node | 2 | 0 | **5,159** | **1.91e+39 pN** |
| **100** | **31,165** | **segment** | **14** | 0 | **0** | **0** |

The node rule terminates in **2 sweeps at every density** — it runs out of things it can see, not out of
overlaps. The segment rule takes 14 and leaves nothing.

### G-6: the affine modulus cannot see interpenetration either

| ρ | `K_A` affine, node-resolved | `K_A` affine, segment-resolved | relative change |
|---:|---:|---:|---:|
| 8 | 14274.711744 | 14274.6903542 | **−1.5e-6** |
| 24 | 42827.7413963 | 42827.7270104 | **−3.4e-7** |
| 100 | 178410.5 | 178409.8 | **−4e-6** |

Removing 5,159 interpenetrations moves the affine areal modulus by **four parts in a million**. The
numbers are not bit-identical — the crosslink count moves 2,492 → 2,491 at ρ = 24 and the nodes move
up to 12.1 nm — so the resolution is reaching the network; the estimator is simply blind to this.

**That is a third blind spot in the same estimator and it is consistent with the other two.**
`ALEPH-PORT-3642` §0 records that the affine bound reports 14,270 pN/µm for a shell that is 2,436
disconnected fragments. It is blind to connectivity, and it is blind to overlap, for the same
reason: it evaluates a *virial at fixed affine displacement* and never asks whether the
configuration is one a structure could hold.

**So G-6 is not yet answered.** What the interpenetration cost is a question about the **relaxed**
modulus, which needs the GPU driver and a PI card citation. What is now on the record is the
narrower and still useful fact that *the affine number does not move*, and that anyone reading an
affine `K_A` as evidence about a shell's steric state is reading something that is not in it.

### G-5: the stretch, and it is where the conflict comes from

| ρ | rule | rms `ΔL/L₀` | max `ΔL/L₀` | max node move |
|---:|---|---:|---:|---:|
| 24 | node | 2.27e-4 | 1.48e-2 | 2.51e-3 µm |
| 24 | segment | 3.28e-4 | 1.46e-2 | 1.08e-2 µm |
| 100 | node | — | 1.48e-2 | 3.53e-3 µm |
| 100 | segment | — | 2.47e-2 | **1.21e-2 µm** |

## The default did not move, and §4 said it would

**§4 of this entry pre-committed:** *"If G-1..G-5 pass, `OverlapResolution.SEGMENT_PAIR` becomes the
default."* **They passed. The switch was not made.**

It fails **`ALEPH-PORT-3636` gate G-C** — *"max node displacement < σ = 0.007 µm (the relaxation does
not move filaments off-shell)"* — at **1.21e-2 µm**, 1.7× the threshold.

That gate belongs to another entry, and **a lane that moves a threshold so its own new code can pass
is the failure this project has recorded more than any other.** `CLAUDE.md` §2.4 and §2.5 are the
rules; this is reported, not resolved.

### The evidence a decision would need, measured rather than argued

G-C states a threshold *and* a reason, and the two have come apart. "Off-shell" is a radial claim,
so it was measured directly at ρ = 100, thickness parameter 0.2 µm:

| rule | max \|move\| | max \|Δ radial\| | radial fraction | shell radius span, before → after |
|---|---:|---:|---:|---|
| node | 3.53e-3 µm | 3.02e-3 µm | 0.496 | 0.298490 → **0.298490** |
| segment | 1.21e-2 µm | 1.19e-2 µm | **0.903** | 0.298490 → **0.301171** |

The segment rule's move is **90% radial** — so it is the kind of motion G-C was written to catch —
and it is **6.0% of the shell thickness**, changing the shell's total radial span by **0.9%**.

**The threshold fails; the reason does not.** Which of the two is the gate is a question for whoever
owns `ALEPH-PORT-3636`, and offering the measurement is the whole of this lane's part in it.

### Status

`OverlapResolution.SEGMENT_PAIR` is implemented, gated and **opt-in**. `NODE_PAIR` remains the
default, so no committed number moves and
`docs/results/2026-08-06-overlap-free-native-cortex/verify.py` still reproduces every claim.

**What is on the record is that the default is now known-wrong rather than suspected-wrong**, with
5,159 measured interpenetrations at native density behind that, and that the thing standing between
the fix and the tree is one threshold in a neighbouring entry.

---

## Amendment, 2026-08-07 — G-6 relaxed, on `4090-1`, and the answer is not the one this entry expected

PI citation: *"gpu 관련하여서 12시간 연장 허가 지금 함"* — 2026-08-06T15:59:18.249Z,
`~/.claude/projects/-Users-sw1-Project-Aleph/46143f30-d1ab-4f8a-9c66-b63204a3b7a4.jsonl#8a6b369e-f305-4a05-98c4-f0cca542db7d`,
narrowed to `4090-1` by a later instruction in the same session. Slurm job 43.
Filamin (60 nm, 8.2e5 pN/µm), strain 1e-3, CG at 1e-3 pN, budget 200,000.

### G-6, the number the whole arc was heading for

| ρ | filaments | `K_A` relaxed, node-resolved | `K_A` relaxed, segment-resolved | change |
|---:|---:|---:|---:|---:|
| 8 | 2,493 | 175.5439 | 176.2796 | **+0.419%** |
| 24 | 7,480 | 1456.309 | 1437.768 | **−1.273%** |
| 100 | 31,165 | 34776.42 | 35744.79 | **+2.785%** *(provisional — see below)* |

**Removing 5,159 interpenetrations changes the relaxed areal modulus by under 3%, and the sign is
not even consistent across density.**

### Why, and it is the finding rather than a disappointment

`ALEPH-PORT-3643` measured the steric force at those contacts as **1.9e39 pN**. A force that large
changing a modulus by 2.8% is not a small effect — it is *no effect*, and the reason is that

> **no code path in this engine ever evaluated it.**

* The **affine** estimator is a virial over axial, bending and crosslink terms at fixed affine
  displacement. It moved by four parts in a million (this entry's earlier amendment).
* The **relaxed** modulus comes from `ResidentCortexRelaxation`, and `ALEPH-PORT-3637` **does not put
  the steric term on the device at all**. `assemble_forces` has no contact term.

So the 1.9e39 pN was computed here, by a diagnostic written for this port, and by nothing else. The
interpenetration was invisible to the *construction* (that is `-3644`), invisible to the *affine
bound*, and invisible to the *solver*. **Three independent blind spots, and the shell passed all
three.**

That reframes what `-3643` and `-3644` bought. They did not correct a modulus. They **made a defect
measurable that the engine had no instrument for**, which is the precondition for the fix and is not
the fix.

### The second question, which was not a gate and should have been

Each shell was counted at build and again on the coordinates the solver returned:

| ρ | rule | segment pairs @ build | @ affine strain | **@ relaxed** | closest, relaxed | max \|F\| |
|---:|---|---:|---:|---:|---:|---:|
| 8 | node | 17 | 17 | 18 | 2.76e-4 µm | 4.7e+17 |
| 8 | **segment** | **0** | 0 | **2** | 6.81e-3 µm | 2.4e-1 |
| 24 | node | 286 | 286 | 286 | 1.60e-5 µm | 5.3e+33 |
| 24 | **segment** | **0** | 0 | **13** | 6.61e-3 µm | 4.0e-1 |
| 100 | node | 5,159 | 5,150 | 5,189 | 6.47e-6 µm | 7.4e+38 |
| 100 | **segment** | **0** | 0 | **123** | 6.53e-3 µm | 4.9e-1 |

**The relaxation creates interpenetration.** A shell built with zero segment overlaps ends a single
relaxation with 123 of them at native density. The solver cannot see contact, so nothing stops it
walking two filaments through each other while it lowers the energy of the terms it *can* see.

**The two defects are not the same size and this entry will not blur them.** The construction defect
left filaments essentially *coincident* — closest approach 6.5e-6 µm, which is `σ/1080`. The solver
defect leaves them **7% inside σ** — closest 6.5e-3 µm against σ = 7e-3 — with a contact force under
0.5 pN. One is a shell that is not a shell; the other is a shell in mild contact. Both are wrong.
Only the first was catastrophic.

**This is `ALEPH-PORT-3646`**: the steric law has to go on the device. It is the larger remaining
defect and it is why `-3645` §5 says the construction fix is a precondition, not a result.

### The ρ = 100 number is provisional, and the reason is in the solver

| ρ | rule | `K_A` relaxed | iterations | L∞ | **rms** |
|---:|---|---:|---:|---:|---:|
| 24 | segment | 1437.77 | 22,274 | 0.0349 | 0.00047 |
| 100 | node | 34776.4 | 36,204 | 0.823 | 0.00391 |
| 100 | **segment** | 35744.8 | **1,167** | 7.69 | **0.0378** |

CG stops when its line search finds no further descent, so **1,167 is not a budget limit — it is the
solver giving up**, thirty times earlier and at ten times the residual of the run it is being
compared against. **A `K_A` measured there is not comparable to one measured at 36,204**, and the
+2.785% rests on exactly that comparison.

`scripts/_ws_cortex_convergence.py` is running on `4090-1` (job 44) over three initial steps and
three tolerances to establish whether that number is a property of the shell or of where the line
search happened to start. **Until it lands, +2.785% is not a result and is marked as such here
rather than quoted onward.**

This is the sixth time on this project that a number turned out to be about the estimator rather
than the thing; it is the first time it was caught in the same session by a column that was added
because the previous five had been.

---

## Amendment, 2026-08-07 — **RETRACTION.** "A third blind spot in the same estimator" is wrong

The amendment above concluded, from G-6's four-parts-per-million affine change, that

> *"the affine estimator is structurally blind to interpenetration, for the same reason it is blind
> to connectivity — three independent blind spots, and the shell passed all three."*

**That is retracted.** Two independent adversarial verifiers refuted it, and this lane reproduced
their measurement before accepting it.

### The measurement had the overlap-sensitive term switched off

`CortexFilamentNetwork.internal_forces()` is `axial + bending + crosslink + **steric**`, and
`steric_term()` branches on `self.steric_law`. Both gate scripts do the same three lines — set
`SEGMENT_PAIR`, read the diagnostic, **set it back to `NODE_PAIR`** — and *then* call
`areal_modulus`. So the affine virial was evaluated under the one law whose own docstring says it
reads exactly 0 pN at a gap of exactly 0 through a between-node crossing.

Re-running the **same** `areal_modulus` on the **same** networks at the **same** coordinates,
changing only `steric_law`:

| ρ | `K_A`, steric = `NODE_PAIR` | `K_A`, steric = `SEGMENT_PAIR` |
|---:|---:|---:|
| 8 | 14,274.71 | **1.73e+12** |
| 24 | 42,827.74 | **2.08e+26** |
| 100 | 178,410.55 | **1.48e+32** |

**The affine estimator sees interpenetration perfectly well.** It was handed a blindfold.

### And a control with the confounders at exactly zero

One verifier rigidly translated a single crosslink-free filament on the clean ρ = 8 shell so that one
crossing lands at `s = t = 0.5` — between nodes — at a closest approach of `0.001 σ`. Rigid
translation with no crosslink makes `ΔE_axial = ΔE_bending = ΔE_crosslink = 0.00e+00` *exactly*, so
the steric term is the only channel:

| | |
|---|---|
| `K_A` with steric = `NODE_PAIR` | 14274.6903542 → 14275.0250462 — **+2.3e-5**, and that is the centroid shift |
| `K_A` with steric = `SEGMENT_PAIR` | 14274.6903542 → **2.40e+31** — 27 orders of magnitude |
| the same crossing placed **at a node** | `NODE_PAIR` also goes to 2.40e+31 |

The last row is decisive. **The blindness is not a property of the affine estimator and not even of
the node-pair law in general — it is a property of where the crossing lands relative to the mesh**,
which is `ALEPH-PORT-3643`'s already-recorded defect, appearing at a third *site* rather than being
a third *defect*.

### What survives

The arithmetic. `14274.711744 → 14274.6903542`, `42827.7413963 → 42827.7270104`,
`178410.547 → 178409.827` all reproduce exactly, as do G-1..G-5 and G-7. The supported statement is
narrower and has a different subject:

> The **node-pair steric law** cannot see a between-node interpenetration — which `-3643` already
> established — and because `internal_forces()` sums the steric term, an affine `K_A` evaluated under
> that law inherits the blindness. G-6 therefore differenced two moduli from which the only
> interpenetration-sensitive term had been removed, and the residual few-ppm motion is the response
> of the axial, bending and crosslink terms to a sub-nanometre node displacement.

**G-6 as run had no discriminating power**, and that was designed into the measurement rather than
discovered by it.

### The relaxed half is unaffected and stands

`ALEPH-PORT-3637` puts **no steric term of any kind** on the device, so the relaxed comparison is
genuinely between two shells the solver could not tell apart on contact. The 0 → 2 / 13 / **123**
segment pairs created *during* relaxation stand, and remain the reason `ALEPH-PORT-3646` is the next
entry.

### Why this is recorded rather than edited away

Six adversarial verifiers were run against three claims from this session; **all six returned
refuted**, and two of the three claims were materially wrong. The amendment above was written with
the numbers in front of it and drew a conclusion the numbers did not support. Deleting it would hide
that the reasoning failed while the arithmetic held, which is the failure mode worth being able to
see.

---

## Amendment, 2026-08-07 — **CORRECTION.** G-6's ρ = 100 delta was `+2.785%`. It is `−0.84%`.

The amendment above marked `+2.785%` provisional because the segment-resolved run stopped at 1,167
CG iterations with an RMS residual ten times the run it was compared against. Two nine-run
convergence sweeps on `4090-1` (jobs 48 and 49, three initial steps × three tolerances each) now
say what that number actually was.

| resolution | spread over 9 runs | best-converged `K_A` | its RMS residual |
|---|---:|---:|---:|
| node-resolved | **0.572%** | 34755.68 | 0.001269 |
| segment-resolved | **0.435%** | 34464.42 | 0.001264 |

The two best runs are converged to the **same** RMS to three digits, so they are comparable in a way
the original pair was not.

> **Corrected G-6 at ρ = 100: `34464.42 / 34755.68 − 1` = **−0.838%**.**
> The published `+2.785%` was the wrong sign and 3.3× the magnitude, and it was a property of where
> one line search happened to stop.

### And a second result, which is the more useful one

**The relaxed areal modulus is reproducible to about ±0.5%, not to six significant figures.**
Across nine runs at ρ = 100 the answer moves by 0.435% (segment) and 0.572% (node) purely from the
solver's starting step and tolerance. Every relaxed `K_A` this lane has written to six or seven
digits — including in this ledger — carries **three**.

That makes the corrected delta of −0.84% **comparable to its own reproducibility**. The honest
statement is:

> Removing 5,159 construction-level interpenetrations changes the relaxed areal modulus at ρ = 100
> by **−0.8%, against a run-to-run spread of ±0.5%**. It is marginally resolvable at best, and it
> is not the several-fold effect a 1.9e39 pN contact force might suggest — because
> `ALEPH-PORT-3637` puts no steric term on the device, so the solver never evaluated that force.

### What did NOT change

The node-resolved shells keep **5,188–5,191** segment pairs after relaxation in all nine runs, at a
closest approach of **6.47e-6 µm** — `σ/1080`, essentially coincident — regardless of tolerance or
step. The interpenetration is a property of the construction, and the solver neither creates nor
removes it there.
