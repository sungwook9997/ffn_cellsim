# The thirteenth was mine, it agreed with what I expected, and it had a kill command behind it

**2026-08-21, ~06:55 KST.** The other twelve defects this night were caught before they cost anything.
This one nearly destroyed a 1.3-hour run in flight.

## 1. What I was trying to do

The 9,000-step τ run emits no progress of any kind. It was sized from **4.06 s/step**; a later
isolated timing of the same two calls said **6.14 s/step**, under which it would be killed at ~6,672
steps and — because snapshots accumulate in RAM and the record is written after the loop — **write
nothing at all.** I needed the actual rate and the process would not tell me.

So I inferred it from outside. Each step reads back **576.0 MB**; at 4 KiB pages that is 140,625
minor faults, so `faults / 140,625` should give steps.

## 2. What it said, and why I believed it

```
1,080,681 minor faults / 47 s  ÷  140,625 faults/step  =  7.68 steps  →  6.12 s/step
```

**6.12, sitting beside the 6.14 I already had.** I read the agreement as confirmation and moved
straight to consequences: the run is doomed, it will save nothing, kill it and relaunch.

## 3. ⚠ The check I "verified" it with cannot fail

```
steps  = faults / 140,625
check:   faults / steps == 140,625   ✓
```

**That is the same division twice.** Dividing by a constant and confirming that multiplying back
returns the original holds for every input, including every wrong one. I wrote it into a message as a
self-consistency check — *"140,625 measured vs 140,625 predicted"* — four hours into a night spent
cataloguing checks whose verdict does not depend on their subject.

## 4. What it actually was

The author of `observe_gamma.py` asked one question: **have you checked transparent huge pages?**

```
/proc/<pid>/smaps_rollup    AnonHugePages:  1,978,368 kB
/sys/kernel/mm/transparent_hugepage/enabled   [always] madvise never
```

Calibrated on the same kernel — twelve cycles of allocate, touch and free two 288 MB float64 arrays:

| | faults per 576 MB cycle |
|---|---:|
| **measured** | **329** |
| my assumption (4 KiB) | 140,625 |
| **error** | **427×** |

With the measured conversion the rate comes out at **14 ms/step**, which is faster than the device
path's own readback and therefore impossible. **So the faults are not the readback at all** — they are
NumPy temporaries inside the estimator, and the method does not count steps. **I know nothing new
about the rate.** Both candidates stand exactly where they were.

## 5. Why this one is worse than the other twelve

Every other defect tonight was found by *using* an instrument and getting a refusal, a crash, or a
number that would not reconcile. **This one produced a number that reconciled beautifully with a prior
I already held, and the reconciliation was the reason I stopped checking.**

> ⚠ **Agreement with an expectation is not verification. It is the condition under which verification
> stops happening.**

And the consequence was not a document. It was a `kill` on a run that had cost 1.3 hours of a
twelve-hour grant, taken on a conversion factor I had assumed rather than measured, defended by a
tautology.

## 6. What stopped it

**A colleague asking for the one measurement I had skipped**, in a message that also said *"512배 차이라
잘못 고르면 진행도가 완전히 틀립니다"* — a 512× spread, choose wrong and the progress figure is
entirely wrong. It was 427×.

⚠ **Nothing in my own process would have caught it.** The number was plausible, it agreed with an
independent-looking measurement, and my check was vacuous. **The thing that worked was sending the
raw numbers to someone else rather than the conclusion** — which is the same reason the device gate's
weak fixture was caught this morning, and it is the only defence in this file that generalises.

## 7. What was not done

⚠ **The run continues.** Nothing about it has changed, and there is no evidence it is off the rate it
was sized from. The candidate that says it fits — 4.06 s/step — was measured **solo**, in a 300-step
run of this exact configuration; the candidate that says it does not was measured **beside this very
run**, and the same conditions expanded a GPU kernel by 1.26× on a card our job owns outright, which
is host-side contention. **That hypothesis has a test that costs no card: re-measure solo when τ ends.**
It is recorded as a hypothesis with its test, not as the answer.

---

## 8. The standing conclusion, and it is not the one these documents are written as if it were

Every audit file this night produced is implicitly addressed to someone who does not yet know the
rule. That framing is wrong, and the sharpest evidence came from the author of `observe_gamma.py`
correcting my praise of their own catch:

> ⚠ *"I said 'pass an explicit centre to both sides' — and then, an hour later, wrote a check that did
> not hold the geometry fixed. I broke the rule immediately after stating it. It is not a knowledge
> problem, it is an application problem. **Knowing does not stop it.**"*

Check that against the night's list. **Every one of these was committed by someone who could have
stated the rule that forbids it, and several within an hour of stating it:**

* the balance-gate table comparing two different cells — committed by the author of the document
  about comparing two different cells;
* τ run I sized from a rate measured on a different cell — **the same error, five hours after
  correcting it**, because `s/step` had stopped feeling like a quantity with a cell attached;
* a page-fault conversion "verified" by dividing by a constant and multiplying back — written four
  hours into cataloguing checks whose verdict does not depend on their subject;
* a device gate whose fixture put ≤1 element per thread — written by the author of the two-stage
  reduction the fixture was meant to exercise;
* a centre-equivalence check that compared two centres — written **one hour** after its author told
  me to stop comparing two centres;
* `assert ... or True` inside the Sanity Gate written to prevent exactly that.

⚠ **So the defence cannot be "state the rule better".** The rules were stated, in writing, by the
people who then broke them. What actually caught these, in every single case, was one of two things:

1. **Executing the thing** — starting the app, running the gate, re-running the measurement. Not
   reading it.
2. **Handing the raw numbers to someone else.** Not the conclusion: the numbers. The stale fixture was
   caught because a PASS line was forwarded verbatim with its element counts; the page-fault error was
   caught because the faults-per-second went out before the rate did.

**Both are properties of the process, not of the person.** That is what these documents should be
recommending, and mostly they are not — they are recommending vigilance, which is what everyone
already had.
