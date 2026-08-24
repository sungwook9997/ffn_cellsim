# ALEPH-PORT-3645 — a threshold and the reason it was written for have come apart

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3645` |
| Lane | `46143f30` Lane W2 (cortex) |
| Status | `PROPOSED` |
| Written | `2026-08-07` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **DEFAULT CHANGE + GATE SUPERSESSION.** No new physics. One enum default moves and one gate's criterion is replaced by the quantity it was standing in for. |
| Aleph target | `aleph/vertical/cortex_surface_coupling.py`, `docs/results/2026-08-06-overlap-free-native-cortex/` |
| Depends on | `ALEPH-PORT-3644` (all seven gates), `ALEPH-PORT-3636` (whose gate this supersedes) |
| Exists because | `3644` §4 pre-committed to a default change if its gates passed. They passed, and it was held back for a reason that turned out to be wrong. |

---

## 0. What is actually being decided

`ALEPH-PORT-3644` measured that the committed cortical construction certifies shells carrying
**5,159 interpenetrations at ρ = 100**, at a WCA force of 1.9e39 pN, and that
`OverlapResolution.SEGMENT_PAIR` removes all of them in 14 sweeps. It shipped **opt-in** anyway.

The stated reason was that switching fails `ALEPH-PORT-3636` gate G-C — *"the relaxation does not
move filaments off the shell. Max node displacement `< σ = 0.007 µm`"* — at `1.21e-2` µm.

**That reason rested on a premise this entry checked and found false.** `3644`'s amendment says the
gate "belongs to another entry" and invokes `CLAUDE.md` §2.4. `ALEPH-PORT-3636`'s header reads
**`Lane | 46143f30 Lane W2 (lead)`** — the *same lane*, one day earlier. §2.4 forbids a lane moving a
guard *to unblock itself*; it does not forbid a lane replacing its own proxy with the quantity the
proxy was written to stand for, in the open, with the old threshold preserved.

So the question is not "may this be changed" but "is the change right", and that is a measurement.

## 1. The threshold and the reason are two different claims

G-C states both, in one row:

```
| G-C | The relaxation does not move filaments OFF THE SHELL.  <- the reason
       Max node displacement.                                 <- the proxy
       < sigma = 0.007 um                                     <- the threshold
```

**"Off-shell" is a radial statement about the shell's extent.** "Max node displacement" is a scalar
over all nodes in all directions, dominated by whichever single node moved most. The two coincide
when displacements are tiny and isotropic, which is what the node-pair rule produces, and that is
why the proxy was serviceable when it was written.

>  **CORRECTION, 2026-08-07 — the span column below is WRONG.** It was measured about the node
>  centroid, which under `FilamentPlacement.POISSON` sits **4.91e-2 µm** from the construction
>  centre — 4× the displacement being resolved. The corrected numbers are in the amendment, and
>  they are what H-2 was decided on. The *radial-fraction* column survives and is conservative.

Measured at ρ = 100 with the declared shell thickness `T = 0.2 µm`:

| rule | max \|move\| | max \|Δradius\| | radial fraction | shell radial span, before → after |
|---|---:|---:|---:|---|
| node | 3.53e-3 µm | 3.02e-3 µm | 0.496 | 0.298490 → 0.298490 |
| segment | 1.21e-2 µm | 1.19e-2 µm | 0.903 | 0.298490 → **0.301171** |

The segment rule's motion is **90% radial**, which is exactly the kind of motion G-C was written to
catch — so this is not a case of the gate firing on something irrelevant. It is a case of the gate's
**magnitude** being calibrated on a rule that could not see most of the overlaps.

**The shell's radial extent grew by 2.681 nm on a 298.490 nm band — 0.9%, and 0.38 σ.**

## 2. The replacement, and why it needs no new constant

A gate that swaps one invented number for another has bought nothing. The criterion below is built
from `σ`, which is the F-actin diameter the card already sources and the same quantity that defines
`r_c` and therefore the whole construction:

> **The shell's radial span may not grow by more than σ.**

The argument is short enough to check. The resolution's entire job is to put two filament *surfaces*
`r_c = 2^(1/6)σ` apart. Making that room can push the outermost filament outward — by at most about
one filament diameter, because beyond that it is no longer making room for a contact, it is
translating the shell. So `Δspan ≤ σ` is the boundary between "settled" and "displaced", and it
carries no fitted constant.

It also measures the right object. **A node moving 12 nm inside an existing 298 nm band is not
off-shell; the band growing is.** G-C's proxy could not tell those apart and this does.

### What is NOT claimed about the replacement

- **Not that σ is the only defensible boundary.** It is *a* boundary derived from a sourced quantity
  rather than chosen to let this change through, and the measured value (0.38 σ) is not near it, so
  the conclusion does not depend on where in `[0.38σ, ∞)` the line is drawn. Had the measurement
  been 0.95 σ this entry would say so and the gate would be worth less.
- **Not that G-C was wrong.** It was right about the rule it was written for, and it caught a real
  difference. It is superseded, not withdrawn, and its number stays in `3636`.

## 3. The default moves

`OverlapResolution.SEGMENT_PAIR` becomes the default for `resolve_cross_filament_overlaps` and for
`build_radial_cortex_network`. `NODE_PAIR` stays reachable, because every published cortex number was
produced under it and reproducing them has to remain possible.

**Numbers will move.** `3644` G-6 measured the affine `K_A` change as four parts in a million, so the
movement is small there; the relaxed modulus is being measured on `4090-1` as this is written and its
delta is not yet known. Whatever it is, it goes in the commit.

**Against leaving the default alone:** the node-pair construction is not a cheaper approximation of
the right answer. It is *known* to return shells with thousands of interpenetrations while reporting
zero. Leaving a known-wrong default because changing it moves numbers is how a defect becomes a
baseline, and this project has a proposal on file saying so.

## 4. The gates — written before the code

| # | Gate | Threshold |
|---|---|---|
| **H-1** | **The default is the segment rule and it works through the builder.** A shell built with `overlap_relaxation_sweeps > 0` and no explicit resolution, at ρ = 8, 24, 48, 100. | **0** segment pairs at every density |
| **H-2** | **The shell stays a shell.** Radial span before and after resolution, at every density. | `Δspan ≤ σ = 0.007 µm` |
| **H-3** | **Every published number is still reachable.** `OverlapResolution.NODE_PAIR` requested explicitly reproduces `3636`'s ρ = 8 / 24 / 100 shells and `K_A`. | bitwise |
| **H-4** | **The retained control is untouched.** `K_A` = 1081.819088 at ρ = 0.5, level 2, through `build_vertical`. | bitwise |
| **H-5** | **G-C's own quantity, still reported.** Max node displacement under the new default. | **measured, no threshold** — it is no longer a criterion and pretending otherwise would hide the change |
| **H-6** | **The committed reproduction still runs.** `docs/results/2026-08-06-overlap-free-native-cortex/verify.py`. | exit 0, with any moved number updated **in the README as an amendment**, never silently |
| **H-7** | **What the default change costs the observable.** Affine and relaxed `K_A` before and after, at ρ = 8, 24, 100. | **measured, no threshold** |

**H-2 can fail honestly**, and at densities this entry has not yet measured it is a real risk: the
span growth is a max over the whole shell, and at ρ = 48 or a thickness of 0 it may be larger than at
ρ = 100. If it exceeds σ anywhere the default does **not** move and this entry reports that instead.

**H-6 is the one that catches dishonesty rather than error.** A default change that quietly alters a
committed result's numbers and leaves its README saying the old ones is worse than not making the
change.

## 5. What is NOT claimed

- **Not that the cortex is now overlap-free in any run.** This fixes the *construction*. The GPU
  relaxation has **no device-side steric term** (`ALEPH-PORT-3637`), so a shell built clean can gain
  contacts during a solve — measured at ρ = 8: **0 → 2** across one relaxation. That is
  `ALEPH-PORT-3646` and it is the larger defect.
- **Not that the moduli are now right.** They are now measured on a structure that is not
  interpenetrating at `t = 0`, which is a precondition and not a result.
- **Not that G-C's replacement is itself validated.** `Δspan ≤ σ` is derived, not tested against a
  case where it *should* fail. A negative control — a deliberately runaway resolution that the gate
  must catch — is worth having and is not in this entry.

---

## Amendment, 2026-08-07 — **H-2 FAILED. The default did not move.** And two defects found on the way

Six adversarial verifiers were run against this session's three headline claims, each prompted to
*refute* rather than confirm and each re-measuring rather than reading. **All six returned
`refuted`.** Every refutation was then reproduced by this lane directly before anything was acted on.

### The gates

| # | outcome |
|---|---|
| **H-1** default leaves zero segment pairs | **PASS** — ρ = 8, 24, 48, 100 |
| **H-2** the shell stays a shell, `Δspan ≤ σ` | **FAILED at ρ = 100: 1.1170 σ** |
| **H-3** `NODE_PAIR` reproduces the pre-`3645` shells | **PASS** — 17 / 286 / 1200 / 5159, exactly |
| **H-4** the retained control | **PASS** — `K_A` = 1081.819088 |
| **H-5** G-C's own quantity | **measured** — 1.17 σ / 1.55 σ / 1.55 σ / 1.73 σ |
| **H-6** the committed reproduction | **PASS** — exit 0 |
| **H-8** *(new)* one crossing, one contact | **PASS** — ratio 1.000000000, was **4.000000000** |

**§4 said: "if H-2 fails anywhere the default does not move and this entry reports that instead."**
It failed. `NODE_PAIR` remains the default. That is what a pre-registered gate is for, and changing
the criterion now — after seeing it fail — to one that passes is the single thing pre-registration
exists to prevent.

### Defect 1 — the span was measured about the wrong centre, and it was this lane's own error

`build_radial_cortex_network` places every filament radially about **`surface.mean(axis=0)`**
(`cortex_surface_coupling.py:1222`). The measurement used `positions.mean(axis=0)` — the **node**
centroid. Under `FilamentPlacement.POISSON` those are not the same point, and the difference is not
small relative to what is being measured:

| | |
|---|---|
| `icosphere_arrays(2, 5.0)` centroid | `5.4e-17` µm — exactly centred, every `\|v\| = 5.000000000000` |
| node-cloud centroid, ρ = 100 | **4.91e-2 µm** from it |
| node-cloud centroid, ρ = 24 | **7.47e-2 µm** from it |
| the displacement being resolved | 1.2e-2 µm — **4× to 6× smaller than the error** |

| ρ | centre | span before | span after | Δspan | Δ/σ |
|---:|---|---:|---:|---:|---:|
| 100 | node centroid *(wrong)* | 0.298490 | 0.301171 | 2.68e-3 | 0.383 |
| 100 | **surface centroid** | **0.202246** | **0.210066** | **7.82e-3** | **1.117** |
| 24 | node centroid *(wrong)* | 0.342279 | 0.342279 | 0.00e+00 | 0.000 |
| 24 | **surface centroid** | **0.202210** | 0.202526 | 3.17e-4 | 0.045 |

**A 0.2 µm shell was reported as having a 0.298490 µm radial span, to six decimal places, and
nothing caught it.** Measured correctly it is **0.202246**, which matches the band the builder
analytically produces — `[R − inset − T, R − inset]` plus a tangent segment's outward bulge
`√(r² + half²) − r` = 0.002259 — to **1.3e-5 µm**. Forty-eight per cent of the published span was
noise in the centre estimate.

This is the seventh time on this project a number turned out to be about the *estimator* rather than
the thing. It is the first time the estimator's own error exceeded the effect it was resolving.

### Defect 2 — a crossing at a node was pushed **four times**

`_segment_contacts` registers one contact per **segment pair**, and `np.add.at` sums all of them. A
crossing whose closest approach lands on a node is inside `r_c` of **both** of that node's segments
on **both** filaments, so it registers **four times**:

```
pair [0 1] x [3 4]   s=1.0 t=1.0
pair [0 1] x [4 5]   s=1.0 t=0.0
pair [1 2] x [3 4]   s=0.0 t=1.0
pair [1 2] x [4 5]   s=0.0 t=0.0     <- one crossing, four pushes
```

Measured on two 3-node filaments crossing perpendicular at their middle nodes: the shared node moved
`8.8716e-3` µm where the node rule moves it `2.2179e-3` — **ratio exactly `4.000000000`** — landing
the gap at 0.02124 µm against a target of 0.00794, an overshoot of **2.68×**.

**This broke the one property `ALEPH-PORT-3644` gate G-3 claimed and could not test.** G-3 asserted
the segment rule *"reduces exactly to the node rule at an endpoint"* and measured a difference of
exactly 0 — on a configuration (two filaments meeting end-to-end) that happens to register **once**.
The gate passed on a case that could not fail. The per-pair identity is real; it does not lift to a
rule-level identity, and G-3's evidence was structurally incapable of noticing.

`segment_steric_energy_and_forces` had the identical duplication and already fixed it by
deduplicating on the contact's own position. **This is that fix at the second site**, which is where
it should have gone at the same time. `H-8` now tests it, and `H-8b` tests that a genuine
between-node contact survives deduplication (2 registered → 1 kept, not 0).

**It is not why H-2 failed.** Measured after the fix, `max |move|` is unchanged at every density —
8.17e-3 / 1.08e-2 / 1.08e-2 / 1.21e-2 — because a crossing landing exactly on a node is measure-zero
in a Poisson placement. It is a correctness fix, not the cause.

### What H-2's failure does and does not license

`Δspan` is a **max-minus-min over two nodes**, and the two nodes are not even the same before and
after: at ρ = 100 under the segment rule the argmax and argmin both change identity. Every
distributional measure of the same change is 45× to 180× smaller — radius std +0.021%, p1–p99
+0.074%, p5–p95 +0.062%, against the span's +3.87%.

**So H-2 may itself be the wrong estimator, and that is recorded here and deliberately not acted
on.** `ALEPH-PORT-3641`'s amendment records this lane failing exactly this way with an L∞ residual —
a max statistic dominated by one node, mistaken for a convergence criterion. Writing H-2 against
another max statistic was avoidable and was not avoided.

But a gate re-run with a friendlier statistic *after* it has failed is not a gate. If `Δspan` is the
wrong criterion, the argument for a better one has to be made **before** the measurement, in its own
entry, and the default stays where it is until then.

### Status

- `OverlapResolution.NODE_PAIR` remains the default. `SEGMENT_PAIR` remains opt-in.
- `_deduplicate_contacts` is landed and gated.
- **No committed number moves.** `docs/results/2026-08-06-overlap-free-native-cortex/verify.py`
  exits 0 and reproduces every claim.
- The `Δspan` growth is not a defect in the resolution. Two crowded filaments *should* separate
  radially in a shell; whether the construction should hold them inside the declared thickness, or
  the declared thickness should be an outcome rather than an input, is a real question and is not
  this entry's.
