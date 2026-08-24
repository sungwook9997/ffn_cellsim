# The balance gate accepts a visibly diverging run — measured, not argued

**Status:** MEASUREMENT. Native, `cuda:0` (RTX 4090), 4,361,496 nodes. Two runs, one control.
Artifacts: `aleph/outputs/ac/world_phase4/phase4_first.json` (relax = 0) and `phase4_relax.json`
(relax = 1e-4 µm/pN). Driver `aleph/scripts/world_phase4_native.py`.

**What is claimed here:** that `ledger.py`'s force-balance predicate cannot distinguish a static
configuration from one whose residual is growing four orders of magnitude, at native population.
**What is NOT claimed:** any magnitude. The coefficients are declared test points and the residual is
therefore not a physical residual. That is exactly why this run can say something about the GATE
without saying anything about the cell.

---

## 1. The claim, before this run

`STATE.md`, the PHASE 4 module and several 2026-08-20 commits assert that the acceptance predicate
**cannot fail for physics**. The argument: `ledger.py:239-259` accepts when
`|reaction + traction|² ≤ tol²`, and `tol²` comes from `assemble_balance_tolerance_ratio(n_terms)` —
the **Higham summation bound**, `n_terms · eps64`. `sf_motor_slice.py:65` says the predicate tests
**adjoint closure**, not convergence.

That was an argument from the source. It had never been measured.

## 2. What was measured

Ten steps, membrane surface tension + discrete Helfrich bending bound over arena-addressed ranges,
`|Σ F|` over every live node against the summation bound.

| step | residual, relax = 0 [pN] | residual, relax = 1e-4 [pN] | tolerance [pN] |
|---:|---:|---:|---:|
| 0 | 2.055e-13 | 1.804e-13 | 2.683e-06 |
| 2 | 2.057e-13 | 2.479e-11 | 2.683e-06 |
| 4 | 2.187e-13 | 1.368e-10 | 2.683e-06 |
| 6 | 1.453e-13 | 1.809e-09 | 2.683e-06 |
| 7 | 2.150e-13 | **7.053e-09** | 2.683e-06 |
| 9 | 3.458e-13 | 1.724e-09 | 2.683e-06 |

**Every one of the twenty steps is inside the tolerance.** The static arm by seven orders of
magnitude; the moving arm, at its worst step, by **380×**.

## 3. What that means, and what it does not

**The static arm is the mechanism, stated plainly.** With `relax = 0` nothing moves and the residual
sits at ~2e-13 pN, which is not a small physical force — it is **float64 summation error on an
internally closed force system**. Every force here appears with both signs, so `Σ F` is zero by
construction and what remains is rounding. **The gate is measuring the accumulator, not the physics.**

**The moving arm is the finding.** With `relax = 1e-4` the geometry moves, the force field changes,
and the residual climbs **1.8e-13 → 7.1e-9 over seven steps — four orders of magnitude**. That is a
run visibly leaving its starting configuration. The gate would have returned `ACCEPTED` at every
step, and been correct to: adjoint closure never broke. Nothing about the run's *behaviour* is
visible to it.

⚠ **This is not a bug in `ledger.py`.** The predicate does what its own docstring says. The defect is
that a convergence gate was resting on it — recorded on (e) 1, and now measured.

⚠ **The residual growth is NOT a physics result.** `relax` is a declared numerical relaxation, not a
drag law; the arena's per-node `mobility` is still zero because a mobility belongs to a drag law
nobody has bound. The growth shows the wiring is live and the gate is blind to it. It shows nothing
about whether a real cell would diverge.

### 3.1 A second, smaller observation, recorded because it will matter later

The static arm's residual is **not bit-identical across steps** — 2.055e-13, 2.185e-13, 1.453e-13 —
while nothing moves and the same kernels run on the same data. That is atomic-add ordering in the
force accumulation, which is nondeterministic on GPU by construction. It is far below anything that
matters here, but **any future acceptance criterion built on a residual DIFFERENCE has to clear this
floor first**, and the floor is ~1e-13 pN at this population. A criterion comparing residuals between
steps would be reading accumulator noise unless its threshold sits well above it.

## 3.2 ⚠ THE SHARPEST FORM: adding physics makes the gate WEAKER

A third arm, run after the first two. The cortex was bound — `laws.network_warp.link_spring_kernel`
over **4,128,840** axial segments and `laws.forces_warp.cytosim_bending_kernel` over **4,060,026**
bending triples, at the built chord rather than the requested step. **8.2 million force terms added
to a run that previously had only the membrane surface.**

| | membrane only | + cortex |
|---|---:|---:|
| force terms | 409,600 faces + hinges | **+ 8,188,866** |
| residual, step 0 | 2.055e-13 pN | 2.228e-13 pN |
| residual, step 9 | 3.458e-13 pN | 2.117e-13 pN |
| **tolerance** | **2.683e-06 pN** | **2.298e-04 pN** |
| s/step | 0.1475 | 0.1618 |

**The residual did not move. The tolerance grew 86×.**

That is the mechanism in its clearest form. `tol` is the Higham bound and scales with the summation's
own magnitude — with `n_terms` and with `Σ|F|` — while the residual stays at machine epsilon because
the added forces are INTERNAL and cancel in pairs. So:

> **Every force added to this model makes its acceptance gate more permissive.**

A gate that loosens as the physics gets richer is not a weak gate; it is a gate pointed at the wrong
quantity. And it is the opposite of the intuition a reader brings — that a model with more physics in
it is being held to a higher standard.

⚠ **Control, and it is a strong one:** the sign gate returns `-1.1504499137865293e-03` pN with the
cortex bound against `-1.1504499137865295e-03` without it — agreeing to **15 significant figures**.
The membrane forces are bit-for-bit what they were, so the added cortex is genuinely additive and the
residual comparison is between the same run plus terms, not two different runs.

## 4. What this changes

Nothing about (e) 1's status: it stays undecidable until an observable exists to write a criterion
against. What changes is that the reason is now a measurement rather than a reading of the source,
and the measurement is at **full native population**, not a slice.

**And it puts a number on the ceiling.** Any predicate reading `|Σ F|` on an internally closed system
is reading rounding. To see behaviour, a criterion has to read something the system does not conserve
by construction — a state observable, which is what (e) 1's stationarity was reaching for and what
cortical tension γ is the candidate for. `world/observe_gamma.py` landed 2026-08-21 to emit it.

## 5. Provenance

* build: `world/step.py`, `world/laws_bind.py`, `world/geometry.py` at the commits of 2026-08-21
* device: `cuda:0`, RTX 4090, inside Slurm job 72
* sign gate: PASS both arms, mean radial force `-1.1504e-03` pN under a positive tension — a closed
  surface pulled inward, asserted before any step so a sign flip aborts rather than producing a
  plausible run. The two arms agree to 15 digits, which is the control on the gate itself.
* cost: **0.1475 s/step** static, **0.2694 s/step** moving, at 4,361,496 nodes.


---

## 7. The general form: a gate is only sound when its quantity and its tolerance scale TOGETHER

Assembled from three sessions' findings, and it is sharper than any of them alone.

§3.2 measured that the balance gate's tolerance grows 86× when 8.2 M force terms are added while its
residual does not move. Session B then noted that its own contact contract is IMMUNE to that, and
gave the reason — not luck, and not the unrepresentable count:

> Every quantity in `check_complementarity()` is a **max-norm and not a sum**:
> `penetration = max(0, −min Φ)`, `adhesion = max(0, −min γ)`, `product = max|Φ·γ|`, and the
> tolerance is `gap_tol·max|γ| + force_tol·max|Φ|`. All four are **intensive in the constraint
> count**. A step with a million active contacts is judged by the same numbers as one with ten.

Ratcheted by replicating the active set 2048× and asserting **not one float moves**, with the
converse also pinned: one bad pair in 100,001 is still a violation, not a dilution.

**Putting the three together gives the rule.**

| gate | the QUANTITY | its TOLERANCE | scale together? |
|---|---|---|---|
| balance (`ledger.py`) | `max` over nodes → stays at machine epsilon | **sum** over terms → grows with the model | ⚠ **no** |
| contact (`check_complementarity`) | `max` over constraints → intensive | `max` over constraints → intensive | **yes** |
| γ stationarity | sum over plane crossings, averaged over 64 orientations → quiet, no floor | (still the PI's) | — |

> **A gate is unsound when its quantity and its tolerance have different extensivity.** The balance
> gate's defect is not that it uses a Higham bound — a Higham bound is the right tolerance for a
> summation. It is that the bound is EXTENSIVE in the term count while the residual it bounds is
> INTENSIVE, so every force added to the model buys slack the physics did not pay for.

⚠ **Two consequences worth carrying into (e) 1's wording.**

1. **Session A's kinetic feedback is the extensive case at its worst.** A binding event ADDS a force
   term, so a residual-based tolerance grows with the very events being judged — and unlike the static
   arm, the term count is not constant between steps. **If (e) 1 stays residual-based, the tolerance
   must be pinned to declared CAPACITY, never to the live term count.**
2. **The property is a norm choice, and a norm choice can be reverted silently.** Session B's ratchet
   exists because someone later saying *"let's make the residual a sum"* would remove the immunity
   without touching a threshold. **Extensivity is not visible in a diff** — it is visible in a test
   that replicates the population and asserts the number does not move.

⚠ **Not claimed:** that a max-norm is always right, or that γ's construction makes it a good gate. The
claim is about the PAIRING. A sum with a sum-scaled tolerance is sound; so is a max with a max-scaled
one. This engine shipped a max against a sum.


---

## 8. ⚠ SUPERSEDED — see `THE_GATE_IS_EXACT_ABOUT_THE_WRONG_THING_2026-08-21.md`

**This section is wrong in two ways and is kept only so the correction has something to point at.**
Its baseline rows are a DIFFERENT CELL (4,361,496 nodes / 0 populations) from its seven-population
row (4,558,554 / 11), so the table is three builds rather than a series; and its "+9.8% step cost"
is a difference between two single timings, quoted as a measurement, from distributions later
measured to be ±7% wide. The replacement measures five runs per configuration on ONE cell, finds the
tolerance ratio to be **626x with a spread of exactly 0.000%**, and finds the residual it bounds
swinging by up to **1037%** across identical runs.

## 8 (as written, superseded). The strongest form, measured last: seven populations bound

Run 2026-08-21 with every bindable strand population bound at once — cortex, microtubule,
intermediate filament, filopodium, lamellipodium, stress fibre and transverse arc.

| | membrane only | + cortex | **all seven** |
|---|---:|---:|---:|
| force terms | 409,600 faces | + 8.19 M | **8,489,866** |
| residual [pN] | 2.055e-13 | 2.228e-13 | **4.631e-12** |
| **tolerance [pN]** | **2.683e-06** | 2.298e-04 | **1.755e-03** |
| margin | 1.3e+07 | 1.0e+09 | **3.8e+08** |
| s/step | 0.1475 | 0.1618 | **0.1620** |

**The tolerance is 654× what it was with the membrane alone.** The gate is now **380 million times
looser** than the residual it bounds.

⚠ **And the cost of the physics is 10%.** Twenty-one times more force terms — 409,600 to 8.49 M — for
a 9.8% longer step. So the trade the extensivity defect makes is the worst possible one: **the model
gets an order of magnitude richer at almost no compute cost, and its acceptance gate gets three
orders of magnitude weaker for free.** Nothing about the run gets harder to pass as it gets more
complete.

**This is the same measurement as §3.2, extended, and it is reported because the scale changes how it
reads.** At 86× a reader can wonder whether the tolerance is merely generous. At 654× on a run where
adding the physics cost 10% of the wall clock, the pairing is the whole story.

⚠ Still no magnitude: seven declared test points, none sourced for this cell.
