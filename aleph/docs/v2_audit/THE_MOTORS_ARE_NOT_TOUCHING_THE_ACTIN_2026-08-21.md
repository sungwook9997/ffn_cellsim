# The motors are 0.250 µm away from the actin they pull, the check that would see it asks a different question, and no PI ruling can fix it

**2026-08-21, ~06:00 KST.** Read out of the TIER 0 curve by the lane that is blocked on the eight
crossbridge constants, then verified independently here from source. **This is upstream of PI queue
item 11 and changes what has to be decided first.**

## 1. The arithmetic

| | source | shell [µm] |
|---|---|---|
| NMII minifilaments | `world/geometry.py:135` — `r_cell − tip_clearance(0.25) − 0.30 = 6.95`, thickness 0.20 | **[6.850, 7.050]** |
| cortical actin | `world/build/__init__.py:197` — `r_shell = r_cell − 0.5·h_cortex = 7.40`, thickness `h_cortex = 0.20` | **[7.300, 7.500]** |

**The two shells do not overlap.** Outermost NMII node → innermost cortex node = **exactly 0.250 µm.**

That single number predicts every feature of the TIER 0 curve, which was measured before anyone
looked for it:

| observed | what a 0.250 µm gap predicts |
|---|---|
| 0 stationed at reach ≤ 0.20 | below the gap — nothing reaches anything ✓ |
| first non-zero at 0.30 | just clears 0.25 ✓ |
| all 442 by 0.60 | gap + both shell thicknesses ✓ |
| `no_cortex_in_reach` beats `no_antiparallel` by 1–2 orders | a rigid radial offset hits **every** minifilament identically and is indifferent to polarity ✓ |

## 2. ⚠ My "the threshold sits on the head offset" was a coincidence

I recorded that the count is still 0 at reach = 0.200 µm, *exactly* the built `head_offset_um`, and
called it striking. **The observation is correct and the causal reading was wrong.** The head arm is
**tangential** (`build/nmii.py`: `ey = cross(c, ex)`), so it cannot reduce a *radial* separation at
all: a 0.200 µm tangential offset costs **2.7 nm** radially (`offset²/2r`). The causal number is
**0.250**, and it is not a property of the motor — it is the difference between two shells.

The lane that found it also rejected its own first explanation (discretisation: `cortex_seg_um = 0.05`
caps midpoint effects at 25 nm, which cannot make 200 nm). **Two candidate causes were killed by
arithmetic before the surviving one was reported.**

## 3. The check that would have caught it asks a different question — a seventh sibling

```python
# world/geometry.py:278
assert fp.nmii_radius_um + fp.nmii_thickness_um / 2 + 0.20 < fp.r_cell_um
```

The comment says *"7.40 was the incumbent's shell radius and it puts the head bulge through the
membrane."* **That is true**, and the correction moved NMII out of the population it exists to pull.

Then `phase1_inside.json` reported `n_outside: 0` and **passed.**

> ⚠ **"Inside the membrane" and "inside the cortex" are different questions, and the motors were only
> ever asked the first.**

The guard compares against `r_cell`. `build/nmii.py`'s own shell-excursion guard compares against
**its own** radius and thickness. **Each is correct within its own scope, and no check owns the
relationship between the two shells.** A constant chosen to fix one containment problem created a
coupling problem, and every gate in the path was looking elsewhere.

⚠ And the guard's `+ 0.20` term appears to treat the **tangential** head arm as a **radial** one —
which, by §2's arithmetic, over-reserves the radial budget by roughly seventy-fold (0.200 µm reserved,
2.7 nm actually used). ⚠ That is an inference from the arithmetic, not a claim about intent, and
`world/geometry.py` is not the reporting lane's file.

## 4. ⚠ The consequence for PI queue item 11 — the ruling cannot land here

The `hand_kmc` NMIIA capture-radius proxy is **0.210 µm**. The separation is **0.250 µm**.

> **The capture radius is smaller than the distance. In this build, no value of any of the eight
> constants produces a single bond.**

**TIER 1 cannot run on this cell whatever the PI decides.** Placement is upstream of the ruling.

⚠ **And note what TIER 0's pre-registered conclusion did and did not buy.** I wrote before the run:
*"if `n_stationed` is 0, no value of the eight constants makes this family contractile."* It was not 0
— but **why** it was not 0 is the answer: because the sweep raised `reach_um` above 0.25, and reach is
a **build-time geometric argument, not a binding radius**. **"A station exists" and "a head can bind"
are different statements in this build**, and only the second is what item 11 is about. The
pre-registered conclusion was sound and the number that satisfied it came from the axis, not the cell.

## 5. Polarity is not the wall, and it is not idle either

`n_minus_row_first_choice` saturates at **192 of 442**, not 442. At reaches where the anti-parallel
requirement causes no failures it is still moving the placement in **250 of 442 (57%)** of cases. PI
decision 4 was necessary; it is simply not what is blocking today. The column exists to separate
*"causes no failures"* from *"costs nothing"*, and it did.

## 6. What landed on the instrument, and why it matters more than this run

`station_census` now returns `nearest_cortex_node_um` (min / p05 / median / p95 / max) for **every**
minifilament, **independent of reach**.

`n_no_cortex_in_reach` only ever said *"further than this reach"*. It never said how far. **That is why
this run had to sweep nine reaches and read the distance back out of where the count moved** — the
instrument was holding the distance the whole time. It is one kd-query.

No recorded number changes and the question is not redefined: a test asserts the distances are
**identical** at reach 0.05 and at reach 1000.0, so they cannot depend on what was asked. Percentiles
rather than a mean, because the failure being diagnosed is **a sub-population sitting at a non-zero
floor**, and the mean of a bimodal distribution hides exactly that. **The next TIER 0 run shows
`min ≈ 0.25` in one call, with no ladder.**

## 7. What has to happen, and whose it is

Putting NMII inside the cortical shell is the fix, and it belongs to `world/geometry.py` and
`world/build/nmii.py` — **not** to the lane that found this, which owns neither. Whether both
constraints are simultaneously satisfiable is arithmetic: the minifilaments must sit in [7.30, 7.50]
while the head bulge stays inside the membrane. Per §2 a tangential arm costs 2.7 nm radially, and
`build/nmii.py`'s own excursion calculation already says so — so **7.40 appears to leave the heads
3 nm inside the membrane**, and the 0.20 subtraction that moved them to 6.95 looks like a tangential
arm charged as a radial one. ⚠ **That must be confirmed by the owner of those files before anything
moves**, and it is PI queue item 14.

---

## 8. Addendum — tonight the builder would still have built 26,520 of them (`da75726e`)

Item 14 fixes the placement. **This is the other half**, and it is the half that matters after the
placement is fixed, because it is what stops the same defect being built again.

⚠ **TIER 0 *counted* stations; it did not stop anything.** Called with `reach_um = 0.6`,
`build_nmii_cortex_crossbridge` would have happily built **26,520 crossbridges** — 442 minifilaments
× 30 heads × 2 sides — straight across the 0.250 µm gap. It would have succeeded. It would have
reported `n_outside: 0`. **26,520 rows that read as motors**, in exactly the way tonight's other six
defects read as verdicts.

It now refuses. Three properties, and each is why it is a **guard** and not a threshold:

1. **It invents no number.** It uses the caller's own declared `kinetics.capture_radius_um`, a value
   the caller has to declare anyway.
2. **It is the kernel's own condition.** `propose_crossbridge_transitions_kernel` tests
   `L <= capture_radius`; this tests the same comparison against the **built** configuration, not a
   paraphrase of it.
3. **It is decisive at t₀ for a structural reason.** `CrossbridgeState` is born all-free by design;
   nothing moves until something binds, and nothing binds beyond its capture radius. **So if every
   bond is out of range as built, the population is inert for EVERY value of the eight constants** —
   not for today's guesses, for all of them.

A *partial* over-range set is **reported and not refused** — `n_beyond_capture_at_t0` rides in the
census — because some heads starting out of range is a placement fact and not necessarily an error.

**And the refusal names the right owner.** A message that reads as a kinetics problem sends the
reader hunting for a rate constant; the cause is 0.250 µm between two builders' shells. So the text
says *"a PLACEMENT result, not a kinetics one"*, says *"inert for EVERY value of the eight PI-GAPs"*,
and points at `station_census().nearest_cortex_node_um`, which answers **without building anything**.
Three tests hold those three properties, because the wrong search is the expensive part.

### Where the guard lives, and why that is the argument

The same move `bond.py` made for the ERM count, aimed at something different: **the relationship
between two populations.** That is precisely what neither builder's guard could see —

    geometry.py:278      head bulge   vs  r_cell          (the MEMBRANE)
    build/nmii.py        minifilament vs  its OWN shell

**Each is correct inside its own scope. Nothing owned the relationship.** It is now owned by the
**consumer** — the code that has to couple the two is where coupling-possibility is checked. Putting
it in either builder would require that builder to know about the other, which is the coupling
`families/__init__.py` exists to prevent by keeping them in separate files.
