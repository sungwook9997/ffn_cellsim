# ALEPH-PORT-3647 — pair the filaments, not the segments, and the duplication stops existing

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3647` |
| Lane | `46143f30` Lane W2 (cortex) |
| Status | `PROPOSED` |
| Written | `2026-08-07` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **RE-DERIVATION.** A change of what the steric neighbour list enumerates. Nothing is read from a provider. |
| Aleph target | `aleph/vertical/cortex_filaments.py`, `aleph/runtime/cortex_relax_kernels.py`, `aleph/runtime/resident_cortex.py` |
| Depends on | `ALEPH-PORT-3643` (the law), `-3646` (whose I-1 this fixes) |
| Exists because | Three deduplication rules have now been written for one defect, and two of them were wrong. |

---

## 0. The defect has been fixed twice and is still there

A rod split into segments reports **one crossing at a shared node as two contacts**, once per
adjacent segment. Three rules have been written against it:

| where | rule | outcome |
|---|---|---|
| `segment_steric_energy_and_forces` (`-3643`) | quantise the contact position at `0.05 σ`, keep the closest | works |
| `_segment_separation_step` (`-3644`) | *none* | **4× over-push**, ratio exactly `4.000000000`; fixed in `-3645` by copying the first rule |
| `k_segment_steric` (`-3646`) | topological: drop a contact at `s = 1` when the segment has a successor | **fails on a bent filament**, where `s = 0.998466` |

The third failure is the informative one. It is not that the rule was badly implemented; it is that
**an exact-endpoint test is a measure-zero condition** and a real filament is never exactly straight.
`-3646` §2a asserted *"the two rules agree wherever both apply"* and that was checked only on
exactly-coincident test geometry — the one case where the failure cannot occur.

Three rules, two wrong, one defect. **The rules are the problem.**

## 1. The measurement that decides the design, taken before the code

If two filaments routinely touched in more than one place, "one contact per filament pair" would be
the wrong law and this entry would not be written. Measured on shells at 0.98 compression, which is
what puts contacts in at all:

| ρ | contacts inside `r_c` | filament pairs | pairs with **> 1** contact | max per pair | 2nd/1st distance |
|---:|---:|---:|---:|---:|---:|
| 8 | 13 | 13 | **0** | 1 | — |
| 24 | 202 | 202 | **0** | 1 | — |
| 48 | 782 | 781 | **1** | 2 | **1.00031** |
| 100 | 3,210 | 3,194 | **16** | 2 | **1.00196** |

Two things, and the second is the one that settles it.

* **Almost every filament pair carries exactly one contact.** 3,194 pairs hold 3,210 contacts at
  native density.
* **Every multi-contact pair is a duplicate, not a second crossing.** The second contact sits at
  `1.002×` the distance of the first. A genuinely separate crossing between two 300 nm rods would be
  at a completely different separation; `1.002×` is the same crossing seen from the segment on the
  other side of a shared node.

So collapsing a filament pair to its single closest approach **removes exactly the duplicates and
nothing else**, at every density measured.

## 2. What changes

The neighbour list enumerates **filament pairs** instead of segment pairs, and one thread scans that
pair's segment products for the single closest approach and applies the WCA there.

* **No deduplication rule anywhere.** One crossing is one contact by construction, which is why this
  removes the class of defect rather than the instance.
* **Bounded work.** At three nodes a filament pair is a `2 × 2` scan. It grows as `(n−1)²` with
  refinement, which is the same growth the segment-pair list already had and is bounded by the
  filament, not by the shell.
* **Cheaper list.** 31,165 filaments instead of 62,330 segments, and the pair count falls with it.

`_segment_separation_step`'s host-side resolution keeps its position-quantised rule for now: it is
measured correct, `-3645` gated it, and changing two things at once would make the next failure
ambiguous.

### The cost, stated plainly

**Two filaments that genuinely cross twice will get one contact.** That is a real approximation and
the table in §1 is the whole argument that it is currently harmless — 0 such pairs at ρ ≤ 24, and at
ρ = 100 the 16 candidates are all duplicates. It is a **density- and length-dependent** claim, not a
theorem: longer filaments, a finer mesh, or a denser shell could each break it, and J-2 is the gate
that keeps measuring it rather than assuming it stays true.

## 3. The gates — written before the code

| # | Gate | Threshold |
|---|---|---|
| **J-1** | **Device equals NumPy**, both filament-pair, energy and per-node force, at ρ = 8, 24, 48, 100. | `≤ 1e-10` relative — the threshold `-3646` I-1 failed at 48 and 100 |
| **J-2** | **How often does one filament pair carry two genuine contacts?** Second-closest over closest, at every density, and after a relaxation as well as at build. | **measured, no threshold.** A second contact at `> 1.1×` the first is a real crossing this law would lose, and the count of those is the result |
| **J-3** | **Refinement invariance survives.** `-3643` F-2's fixed crossing at `n` = 3, 5, 9, 17, 33. | varies `< 5%`, as F-2 required |
| **J-4** | **Two crossing rods still repel**, `-3643` F-1's configuration with the crossing between nodes. | `\|F\| > 0` and rising as the gap closes |
| **J-5** | **The shell stops interpenetrating.** Closest approach after one relaxation at ρ = 8, 24, 100, with the term on and off. | **on: `≥ σ`.** Off it is 0.94–0.97 σ. This is the criterion `-3646` I-2 should have carried and did not |
| **J-6** | **The solver still converges.** Iterations and RMS residual, term on against off. | within the `±0.5%` spread `-3645`'s W-3 measured |
| **J-7** | **Every committed number is untouched.** `steric_on_device` off by default; `StericLaw.NODE_PAIR` still the default; `K_A` = 1081.819088. | bitwise |

**J-5 is written as a separation and not as a count, deliberately.** `-3646` I-2 asked for *zero
pairs inside `r_c`* and failed on a shell that was behaving correctly, because a WCA equilibrium has
contacts inside `r_c` by construction — `r_c = 1.1225 σ` and the potential's zero is there. That is
the eighth threshold this lane has written against the wrong estimator, and this one is written
against the quantity `-3644` already identified as the physical one.

**J-2 carries no threshold and can make this entry wrong.** If a relaxation drives two filaments into
a genuine second crossing — which the build-time table cannot rule out, because it was measured at
`t = 0` — then one contact per pair loses a force, and the answer is a per-pair *list* of contacts
rather than a single one. That is a bigger kernel and it would be a different entry.

## 4. What is NOT claimed

- **Not that this makes the modulus move.** `-3646` I-6 measured excluded volume at ~0.1% of the
  areal modulus at ρ ≤ 24 and `-3644`'s corrected G-6 at −0.84% for construction interpenetration.
  This makes the force *correct*, which is a different property from making it *large*.
- **Not that the host-side resolution is fixed by this.** `_segment_separation_step` keeps the
  quantised rule.
- **Not that the segment-pair law is withdrawn.** `StericLaw.SEGMENT_PAIR` stays, because `-3643`'s
  F-1..F-5 were measured on it and reproducing them has to remain possible.
- **Not that ε = 50 pN/µm is sourced.** Row 23 of the axis registry, class `axis`, marked UNSOURCED,
  and still unsourced after this.

---

## Amendment, 2026-08-07 — **J-1 passes everywhere `-3646` failed.** J-5 fails, and the reason is a parameter

`4090-1`, job 54. PI citation `8a6b369e-f305-4a05-98c4-f0cca542db7d`, narrowed to one card.

### J-1, which is the whole reason this entry exists

| ρ | filament pairs | E device | E NumPy | rel ΔE | **rel ΔF** |
|---:|---:|---:|---:|---:|---:|
| 8 | 5,263 | 2.480253753e-06 | 2.480253753e-06 | 1.7e-13 | 8.2e-13 |
| 24 | 47,803 | 3.875193224e-05 | 3.875193224e-05 | 5.9e-14 | 4.2e-12 |
| **48** | 191,939 | 1.556686470e-04 | 1.556686470e-04 | **2.6e-13** | **5.2e-12** |
| **100** | 834,855 | 6.999367663e-04 | 6.999367663e-04 | **1.5e-14** | **6.1e-12** |

**PASS at every density**, including the two where `-3646`'s segment-pair kernel was **44% wrong**.
Removing the deduplication rule removed the defect, which is what §0 predicted and is the only part
of this entry that was a prediction.

The pair list is also **35% smaller** at native density — 834,855 filament pairs against 1,284,469
segment pairs — because a filament is a coarser object than its segments.

### J-5, and the failure is informative

| ρ | steric on device | `K_A` relaxed | iterations | rms | **closest / σ** |
|---:|---|---:|---:|---:|---:|
| 8 | off | 175.4536 | 8,281 | 0.00063 | 0.9726 |
| 8 | **on** | 175.383 | 10,312 | 0.00083 | **1.0993** ✓ |
| 24 | off | 1436.5 | 30,560 | 0.00082 | 0.9443 |
| 24 | **on** | 1450.879 | 4,760 | 0.00463 | **0.9696** ✗ |
| 48 | off | 6440.468 | 61,848 | 0.00181 | **0.1870** |
| 48 | **on** | 6442.826 | 58,217 | 0.00159 | **0.8954** ✗ |
| 100 | off | 34447.07 | 50,961 | 0.00123 | **0.0737** |
| 100 | **on** | 34465.5 | 43,012 | 0.00233 | **1.0097** ✓ |

**J-5 is recorded FAILED**: it required `≥ σ` at every density and got it at two of four.

**But look at what it does to the shells that were worst.** Without a contact force the relaxation
drives filaments to **0.187 σ** at ρ = 48 and **0.0737 σ** at ρ = 100 — that is 0.5 nm apart on a
7 nm diameter, which is through each other, and it is *the solver doing it*, not the construction.
With the term: **4.79×** and **13.7×** further apart.

> The two failures are at **0.9696 σ** and **0.8954 σ** — a contact compressed 3% and 10%, not an
> interpenetration. The distinction is `-3644`'s and it is the one that matters here.

### Why it does not reach σ, and it is not the law

`ALEPH-PORT-3643`'s WCA carries `ε = 50 pN/µm`, which is **row 23 of the axis registry, class
`axis`, marked UNSOURCED**. The crosslinks in this shell carry `8.2e5 pN/µm` — **16,400× stiffer**.
A crosslink under tension can therefore compress a steric contact well inside σ and the energy still
goes down, because the barrier it is pushing against is four orders too soft to stop it.

**So J-5 is measuring the parameter, not the law**, and the evidence that it is the parameter is
that J-1 passes to `1e-14` — the force is exactly the one `-3643` specified, evaluated exactly
right. What is unsourced is how hard a filament is.

That is a finding this arc has been walking toward since `-3643` and it is now unavoidable: **the
cortex's excluded volume has no sourced strength**, and every structural conclusion above depends on
it. Sourcing `ε` is the next thing worth doing and it is a literature question, not a code question.

### J-6, and it agrees with everything else

| ρ | `K_A` off → on | change |
|---:|---|---:|
| 8 | 175.4536 → 175.383 | **−0.040%** |
| 24 | 1436.5 → 1450.879 | +1.001% *(the on-run stopped at 4,760 iterations against 30,560; confounded)* |
| 48 | 6440.468 → 6442.826 | **+0.037%** |
| 100 | 34447.07 → 34465.5 | **+0.054%** |

**Under a tenth of a percent wherever the two runs are comparably converged**, against a `±0.5%`
run-to-run spread. Consistent with `-3646` I-6 and with `-3644`'s corrected G-6.

> **The consolidated result of the whole `-3643`…`-3647` arc:** excluded volume changes this
> cortex's areal modulus by **less than 0.1%** and changes how far its filaments are from each other
> by up to **13.7×**. The modulus is not a probe of contact. Anyone reading `K_A` as evidence that a
> shell is physically sound is reading something that is not in it.

---

## Amendment 2, 2026-08-07 — J-5's failure was the solver, and the distribution says so

`ALEPH-PORT-3649` found that the CG line search compared an energy measured *without* a retraction
against energies measured *with* one, and stopped as soon as per-step progress fell below that fixed
`5e-13` penalty. Every J-5 number was taken under that.

Re-measured with the fix, and reporting the **whole contact distribution** rather than the minimum
J-5 gates on:

| ρ | contacts inside `r_c` | min | p5 | median | p95 | max | **< σ** | < 0.95 σ | < 0.9 σ |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 1 | 1.1001 | 1.1001 | 1.1001 | 1.1001 | 1.1001 | **0** | 0 | 0 |
| 24 | 11 | **1.0108** | 1.0426 | 1.1153 | 1.1208 | 1.1214 | **0** | 0 | 0 |

**Not one contact is inside σ at either density.** The whole distribution sits between `1.01 σ` and
`1.1214 σ`, and `r_c/σ = 1.1225` — so every contact is between the filament surface and the cutoff,
which is exactly and only what a correct WCA equilibrium looks like.

### What this changes

`ALEPH-PORT-3647` recorded J-5 **FAILED** on the strength of `0.9696 σ` at ρ = 24 and `0.8954 σ` at
ρ = 48. The ρ = 24 figure was already shown to be convergence (study C reached 1.0108 at 18,649
iterations); with the line-search fix it is **1.0108 with the entire distribution above σ**.

**ρ = 48 and ρ = 100 have not been re-measured with the fix** — `slurmctld` went down on the
workstation at about 13:5x KST and no GPU allocation has been possible since, and those two densities
are not CPU-feasible in the remaining window. **J-5 therefore stays recorded FAILED**, because the
density it failed at has not been retested. It is *suspected* to have been the same defect, and
suspicion is not a gate outcome.

### And J-5's criterion was the wrong estimator, which is separate from whether it passes

`closest approach` is a **minimum over every contact in the shell**. At ρ = 24 that one number was
1.0108 while the median was 1.1153 — the min sits 4.5% below the median, so a single marginal pair
decides the gate. This lane has now been wrong nine times gating on an extremum, and the table above
is the shape a criterion should have.

**Not re-gated.** The corrected criterion belongs in its own entry, pre-registered, and J-5 stands
as written until ρ = 48 is remeasured.
