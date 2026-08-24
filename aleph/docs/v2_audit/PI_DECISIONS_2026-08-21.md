# PI decisions, 2026-08-21 00:28 KST — twenty items, sixteen closed

**Status:** RATIFIED. Transcribed from the PI in session, transcript
`d380041d-06a6-49ea-9daf-1e4e57e76bf0`, 2026-08-21 00:28 KST. An agent may write this citation and may
never write a grant; this file is the citation, and the decisions below are the PI's, not this session's.

**What it unblocks.** `PHASE3_DECISION_QUEUE_2026-08-20.md` closed 2026-08-20 with **14 modules, 0
families, 82 questions**. Sixteen of the twenty items above it are now answered. This file is what a
session reads before opening a family module.

⚠ **Two of these overturn this session's own recommendation.** They are marked. A recommendation that
survives is worth less than one that was overruled with a reason, and the reasons are recorded so the
next session does not re-argue them.

---

## TIER 1 — the ones that unblock several lanes

### 1. The three cardless CONTACT edges are NOT bonds — they are constraints ✅ (a)

`membrane_cortex_contact`, `nucleus_cortex_contact`, `lamellipodium_membrane_contact`.

**Ruling:** a non-penetration field, outside `BondFamily`. `world/` gets a destination for it.

**Why it is the right shape:** the missing chemistry card is not a gap in these edges, it is their
IDENTITY — contact has no binding chemistry. And `BondCount` would demand *how many*, while a contact
count is a function of CONFIGURATION, not a population: it changes every step. aLENS places the same
physics at unilateral complementarity `0 <= Phi_u ⊥ gamma_u >= 0`, not at a bond.

**Closes:** D2, D6, and the third module that shares the shape.

#### ⚠ 1a. CORRECTION — the enumeration in this section was wrong. The RULING is not.

Reported by session B on 2026-08-21 while landing `world/contact.py`, and verified here in the frozen
source rather than relayed. **This session mis-transcribed which edges are cardless.**

`lamellipodium_membrane_contact` is **not** cardless. `contracts.py:791-796`:

```python
"lamellipodium_membrane_contact", ConnectorFamily.CONTACT,
"lamellipodium", "membrane", True, True,          # kinetics, commit_on_accept
chemistry_card="actin_membrane_brownian_ratchet_contact",
```

The cardless CONTACT edges are `membrane_cortex_contact` (`:981`), `nucleus_cortex_contact` (`:960`)
and `membrane_ecm_contact` (`:945`) — all three `False, False`.

⚠ **And `membrane_ecm_contact` is one of the eight EXTERNAL connectors this port plan excludes.** So
inside the cell there are **two** cardless CONTACT edges, not three. The count in the heading of §4 item
1 — *"one ruling closes three modules"* — was arithmetic on a miscount.

**The ruling stands anyway, and for a better reason than the enumeration gave it.** The unilateral gap
law is shared by FIVE declared edges (`connector_joints.py:7-13` calls them *"the same unilateral gap
law"* in the source's own words), and `DUPLICATE_CRITERION` already says a shared force law is not a
duplicate. So `lamellipodium_membrane_contact` has **two shapes at once**:

* **a constraint — unconditionally.** A barbed end may push the membrane and may never pull it inward,
  under any reading of the ratchet. Declared 2026-08-21.
* **a bond — only under a TETHERED ratchet reading.** `MogilnerOster2003_BJ` is the literature that
  would decide it and its ingested chunk count is **zero**, so the existing refusal stands verbatim.

The PI's 3-way grouping and the source's 3/2 split are therefore both correct: they are counting
different things. **This is a report that the edge is two, not a request to amend the decision.**

#### 1b. What session B settled that the decision did not ask

**The duplicate criterion for cardless things is the GAP FUNCTION.**

> A constraint's identity is its gap function — which primitives, measured under which metric. The
> surface pair is SCOPE, not identity, and the force law is downstream and decides nothing.

Same rule as `bond.py`'s, aimed at a different object. A bond's physics lives in its chemistry, so its
identity is the card; a constraint has no chemistry and **its physics IS its geometry**, so its identity
is Φ. Both reject the same thing — the declared component pair, which is bookkeeping on one side and
scope on the other, and identity on neither.

**Consequence: the five merge into ONE `ContactFamily` with several scopes.**

**And the count is not unanswered — it is UNREPRESENTABLE.** `ContactScope`/`ContactFamily` have no
field a count could occupy, and a `dataclasses.fields()` check fails the moment one appears. The
per-vertex trap each module recorded is now *unreachable* rather than merely not-yet-stepped-in.

**There is no stiffness in the base form.** Following aLENS, `0 ≤ Φ ⊥ γ ≥ 0` with γ a multiplier — a
stiffness does not exist. `regularisation_stiffness=None` is that shape; a number makes it a penalty
APPROXIMATION and `provenance_row()` records which was in the loop. **Zero is refused**: it would read
as configured while switching the constraint off.

**Still open, and not closed by this decision:** `d0`, the contact distance, is registered in no scope.
The port source's test `0.030 um` is a FIXTURE, and it is the first value someone needing a number would
promote, so the tests block it. And `membrane_ecm_contact` is a fourth scope of the same law with no
module — recorded so covering three of three does not read as complete.

### 2. ⚠ NMII crossbridge opens the kinetics machinery — NOT ERM ✅ (b), OVERRULING this session

This session recommended ERM on the arithmetic: three of its four `BondCount` answers already exist, so
one parameter would have opened it, whereas a processive motor needs a rate law, a duty ratio and
catch-slip detachment all written for the first time.

**The PI chose the motor anyway.** Recorded without re-argument. It is now buildable because decision 6
landed in the same sitting and fixed the head count the builder refuses without.

**Consequence to hold on to:** the first kinetic family is therefore the HARDEST one, so `bond.py`'s
deferred rate law is being designed against the most demanding consumer rather than the easiest. That is
a defensible reason to choose it and it is how the interface should be judged.

### 3. ERM density — a declared axis, band and test point, marked TEMPORARY ✅ (a)

**Not** the running 231.8/µm². That number is `subdiv 7 = 163,842` vertices over `4π(7.5)² = 706.9 µm²` —
**a mesh artefact wearing a physiological label**, the original ERM defect and the reason the whole class
has a name.

⚠ **The three L0 runs of 2026-08-20 ran `subdiv 3` = 0.9/µm²** — two and a half orders below even the
diagnostic-only proxy.

`temporary` is the PI's word and it belongs in the artifact: this axis expires when a sourced areal
density arrives, and any result standing on it carries that.

### 4. `StrandPopulation` carries PER-STRAND polarity ✅ (a)

Today it is one `int` for the whole population, documented as *"+1 if the barbed end is the last node of
each filament"* — so every filament runs the same way and **an anti-parallel pair is not merely
unverifiable, it is unrepresentable.** A bipolar minifilament straddles two ANTI-parallel filaments; in a
uniform-polarity population every station is a spring wearing a motor's name.

**Not free:** polarity selects `tip_node` (`strand.py:202`), which every barbed-end law addresses, so a
per-strand array propagates into the growth kernels. Owners: A2 (the population record) and G (which
stations are legitimate).

**Unblocks:** the cortex contraction channel — wall 1 of 2. See decision 5 for wall 2.

---

## TIER 2

### 5. A head-resolved population is NOT fed to the lumped channel ✅ (a)

`laws/myosin_linear.minifilament_kernel` is one link between two actin anchors; no head appears. G's 442
are 34 nodes each and exist FOR head resolution. Feeding them in demotes them to what
`engine/nmii_actuator.py:492` refuses as *"aggregate or two-anchor NMII mechanics is diagnostic-only"*.

If a cortex contraction slice is wanted, the honest route is a SEPARATE lumped diagnostic built from the
cortex builder's own filament pairs and labelled as one. The head-resolved population survives undemoted.

### 6. 30 heads/side and 2.0 pN ratified as a declared axis ✅ (a) — ⚠ WITH A CONDITION

**PI note: "band length must be rethought."** The band this session proposed, [10, 30], is NOT ratified.
Only the test point 30 and the axis are.

That condition is well aimed. The band's two ends do not have the same standing:

* **10** comes from `KB-DRAFT-3.B-14`, whose citations field says of itself *"Hill 1938; Kovács 2003;
  Stam-Hocky 2015; Billington 2013; Veigel 2002; Stachowiak 2009 — ALL absent from SE/KU/corpus."*
* **30** comes from `KB-3.18`, `verified/High` — but verified as *"Cortex molecular composition"*, and
  its nine citations are composition REVIEWS. Neither Billington nor Stam is among them.
* `source_audit`: `Billington2013` = **CHECK**, `StamAlbertsGardelMunro2015` = **CHECK**.

**So the band's endpoints are a draft row's self-declared-absent citations and a review-level number read
out of scope.** A band is a claim about what values are admissible, and neither end supports one.

**Open, carried forward:** what the band actually is. Until then the axis has a test point and no band,
and that is stated rather than papered over.

**Measured tonight, so the fork is priced rather than argued:** at H=10 the arena is 4,762,742 nodes and
at H=30 it is 4,780,422 — the head population forks 2.8x and the arena total moves **0.37%**. This is not
a memory decision. What it sizes is PHASE 3's bond capacity array and the KMC event count.

### 7. `KB-DRAFT-7-05` is RE-SCOPED, not retired ✅ (b) — PI: *"Hill 1938 too old"*

The row calls a Hill force-velocity the sanctioned oracle; the engine runs LINEAR by PI ratification
2026-07-07 (`laws/myosin_linear.py:7-10`). `myosin_linear.py` noted *"CLAUDE.md to be annotated"* and the
KB row was never updated, so engine and knowledge base say opposite things and the KB row carries
`confidence: High`.

**Ruling:** scope it. Hill 1938 is correct FOR MUSCLE. Retiring the row deletes something true; leaving
it is a contradiction. Writing the scope is the fix — and this is KB failure type 1 (a row verified for
one purpose, quoted for another) being closed the way type 1 should be closed.

### 8. `sf_arc` is a DIFFERENT STRUCTURE, not a rename ✅ (a)

`sf_arc` is a transverse/dorsal ARC; what PHASE 1 built is a VENTRAL FA-to-FA bundle. Renaming the
contract to match the build would attach an existing name to a structure nobody built.

**Consequence:** an `sf_arc` population needs its own builder. `test_families_census`'s
`xfail(strict=True)` comes off only when it exists.

**Separately: `nucleus` -> `nuclear_envelope` IS a rename** and is unified on the latter.

### 9. LINC isoform split + SUN-site shared budget -> PHASE 4 ✅

Deferred by the PI to PHASE 4 rather than answered now. D3/D4/D5 stay blocked and that is deliberate.

### 10. Perinuclear cap — declared as a band ✅

### 11. Leading-edge width — 5 µm test point ✅

The sourced band is 5–20 µm (`KB-3.7`, ~100 filaments per µm of leading edge); the lower end is the test
point. The current 200 filaments implies a 2 µm edge — **2.5–10x under band.**

⚠ Not cosmetic: tip load is shared among the barbed ends at the edge, so an under-count raises
per-filament load and **the ratchet reports that as stall**. A run at a placeholder count can produce a
number that reads as "the motor is weak".

---

## TIER 3 — approved as recommended, EXCEPT 12

### 12. ⛔ NOT decided — stays open

**Whether a conclusion written before `world/` existed may be read as a conclusion about `world/`.**
This session recommended *no by default, with symbol-level re-confirmation as the promotion path*. The PI
did not take it. The item stays open and no memo conclusion may be promoted to the canonical tree on this
session's say-so.

Mechanical path-rekeying (~330 citations) is NOT blocked by this — chronology settles which file each
memo meant, and `git log --follow` returns it.

### 13–20. Approved ✅

| | item | ruling |
|---|---|---|
| 13 | `docs/session_prompts/**` owned by no lane | add to the `contracts` lane |
| 14 | `.gitignore` forbids the figures `CLAUDE.md` requires committed | narrow exception for the PHASE 1 render |
| 15 | `nvidia-smi --accounting-mode=1` | PI enables it; `exact_peak_gpu_bytes` is otherwise permanently null |
| 16 | `STRUCTURE.md` names directories that do not exist | repair, delegated to this lane |
| 17 | master-plan IDENTITY claim | remains deferred |
| 18 | (e) 1 stationarity amendment | remains deferred — see below |
| 19 | blebbistatin `f_stall` sweep | remains deferred |
| 20 | render-commit policy | with 14 |

⚠ **18 deferred has a consequence the PI asked to work around, not around it.** PHASE 4 is *"the step,
AND WHAT ACCEPTS IT"*, and the acceptance predicate is undefined: `balance_ok` is a Higham rounding bound
and cannot fail for physics; `descent_ratio` reads a post-rollback value and cannot succeed. So tonight's
PHASE 4 separates **"the step ran"** from **"the step was accepted"**, and emits
`ACCEPTANCE_UNDEFINED` rather than a verdict.

**And the PI directed the blocker itself be attacked:** emit cortical tension γ on the resting path. (e) 1
is undecidable today because a stationarity gate cannot be written against an observable the run does not
produce, and γ is that observable. Emitting it turns (e) 1 from undecidable into answerable.

---

## GPU

**4090-1, 12-hour extension approved.** ⚠ Not yet submitted: `gpu-submit` refuses while job 72 holds the
card's token (`selected GPU token is not free`). Job 72 has ~4.5 h left at the time of writing; the
extension goes in when it frees. No other card is taken — the grant names 4090-1.

## Scope directed for tonight

> **"페이즈 3 다 끝내고 기계+텐션 방출하고 페이즈 4까지 완성"** — finish PHASE 3, build the step
> machinery and emit tension, complete PHASE 4.
