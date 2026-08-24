# Four new lanes — session prompts, 2026-07-29

Written because the four areas below are more than one session can hold, and because a session that
starts without the day's findings will repeat them. Each prompt is self-contained: copy one into a
fresh session verbatim.

**Every lane starts the same way** — the boot protocol in `CLAUDE.md` §Session boot, then `STATE.md`,
then `STATE_NONQUOTABLE.md`. The anchors written on 2026-07-29 that every lane needs:
`#native2026-07-29` (a census scalar that reads 1.0 for two runs 1.95x apart in cost),
`#popkinetics2026-07-29` (which components can have an emergent population and which cannot),
`#gpuspeed2026-07-29` and `#gpu2026-07-29` (the second GPU host and what it is worth).

**The one reframe that reorganises all four.** Population counts — MT, IF, filopodia, FA — are
**OUTPUTS the platform infers**, not inputs a PI signs. `CLAUDE.md` §What this is says so
("Parameters are OUTPUTS; literature supplies priors"), the cell says so (MT count moves with dynamic
instability; filopodia form and retract), and the engine already refuses the wrong framing in code:
`aleph/engine/population.py` — "count changes only via the free-list, **never by resampling to a
target**". What needs a literature prior is the **rate**, never the count. A lane that asks "what
number should this be?" has taken a wrong turn; the question is "what events move it, and what is the
prior on their rates?"

**The defect class every lane will meet.** Four guards were found on 2026-07-29 that each printed a
pass while covering less than their name promised: `fraction_of_native` (cortex only),
`t0` (the first SAMPLED step, not step 0), `balance_ok_d` (adjoint wiring, not convergence), and a
kernel-cache check that asked about directory names when the arch tag is a file component. Zero of
them were arithmetic bugs. Before trusting any scalar named like a guarantee, read what it computes.

---

## Lane 1 — `unknowns-register`

```
export FFN_SESSION=unknowns-register
```

Owns: `aleph/docs/v2_audit/**` (already `free` in `ownership.yaml`). Touches no runtime code.

**The problem.** `PI_GAP_EVIDENCE_CARDS_2026-07-25.md` holds 33 cards and they are not the same kind
of thing. IF subunit-exchange turnover is unknown because the literature is paywalled and the corpus
genuinely lacks it (searched 2026-07-29: vimentin/keratin x exchange/turnover/half-life/FRAP over the
whole content layer returns 9 candidate chunks and ZERO extractable rates). MT count is "unknown"
only in the sense that it is what the platform exists to infer. Mixing them cost this project a wrong
turn on 2026-07-29, when the Lead told the PI that MT counts needed a signature.

**Deliverable.** One register, every entry classified:

| class | meaning | what closes it |
|---|---|---|
| `LITERATURE-ABSENT` | a rate/constant nobody has published, or paywalled | acquiring the paper |
| `INFERENCE-TARGET` | an output the sweep recovers with uncertainty | a prior + identifiable observable |
| `DESIGN-DECISION` | a modelling choice with no external truth | PI ratification |
| `MECHANISM-MISSING` | the slot exists but no event can move it | code |

Cover the 33 cards **and** every UNSET slot in code (`grep -rn "PI-GAP\|PI_GAP\|UNSET\|NO default"`
under `aleph/ac/`). Two worked examples to calibrate against: `DynamicInstabilityRateCard` is
`INFERENCE-TARGET` with its prior already encoded (`rusan_2001_llcpk_proxy`, `ratified=False`) —
mechanism bound, awaiting data, not a signature. `IntermediateFilamentTurnoverCard` is
`LITERATURE-ABSENT` with no candidate at all.

**Then the column that makes it useful:** for each `INFERENCE-TARGET`, which observable constrains it,
and by which measurement method. `aleph/docs/v2_audit/CELLTYPE_MODALITY_MATRIX_2026-07-29.md` is the
input. Targets with no constraining observable are the ones stage 4 cannot recover, and saying so is
the point — `ROADMAP.md` makes the progress denominator "identifiable parameters, MEASURED".

**Do not** create `ValidationGate` or `ModelContract` rows: `CLAUDE.md` says those are PI-authored.

---

## Lane 2 — `measurement-modality`

```
export FFN_SESSION=measurement-modality
```

Owns: `aleph/virtual_cell/observation_operator.py`, a new `aleph/engine/probe_*.py`,
`aleph/scripts/probe_*.py`. **Split the `virtual-cell` lane's glob before starting** — it currently
claims `aleph/virtual_cell/**` and the two would collide.

**The idea, in the PI's words.** Cortical tension measured by different protocols differs enormously.
That spread is not noise to be averaged away — it is compartment resolution. Implement each real
experimental protocol as a virtual measurement, apply them all to ONE cell state, and require the
*pattern of differences* to match what the literature reports. Where it does not match localises which
compartment is wrong, because each protocol touches a different one.

**The target is already ratified and quantitative.** `KB-6.1.6` (PI-ratified) classifies three
MDA-MB-231 "cortical tension" numbers as **distinct quantities by protocol**: tether/membrane
~0.05 mN/m (transient), subsecond cortical prestress ~10 mN/m @0.01 s, AFM cortical-shell surface
tension ~0.30 mN/m (seconds). A 200x spread, on one cell line, recorded as physics rather than
disagreement. Supporting rows: `KB-1.14` (AFM ~1 um tip is 10-100x stiffer than bulk rheology),
`KB-6.1.1` ("Anchor RATIOS not absolute Pa"), `KB-6.3.1` (AFM 3.2 kPa vs micropipette 0.1-0.5 kPa on
the same cells), `KB-3.B1.1` (the tether formula `T = F^2/(8*B*pi^2)`), `KB-6.1.7` (MCF-7 parallel-plate
AFM gamma_eff = 0.27 mN/m, n=27, IQR reported).

**What exists and what does not.** `sandbox_inverse.measure_modality_contrast` already asks the right
question — does an added modality lift the Fisher information along the direction the baseline was
worst at, with the direction fixed BEFORE seeing the augmented result. It is wired to imaging
modalities only. Mechanical protocols appear nowhere in the codebase (searched: aspiration,
microplate, optical stretcher, indentation — zero hits), and `observation_operator.py` declares AFM
force-indentation and traction in `UNIMPLEMENTED_OPERATORS`. Note the discipline there before adding
anything: `validate_operator_registration` refuses an operator that claims a BLOCKED modality, added
after an external review smuggled a hardcoded `5.0` through as a traction reading.

**Suggested order, from the matrix rather than from one pair.** The non-adherent lineage first —
HL60 / Jurkat / K562 (`KB-6.2.4`) is "the cleanest cortex-only ref (no FA/clutch)", carries TWO methods
(AFM + microwell), reports HL60 855+/-670 vs Jurkat 48+/-35 Pa (~18x), and cortex is the only
compartment this engine currently runs at native. Then AFM sharp-tip vs colloidal on the breast triplet
(the only same-family pair differing solely in probe size). Then MRS, whose value is that ONE run
yields eta (separates, d'=3.70) and G (does not, d'=0.97) — a target that is structurally immune to
outcome-fitting, because reproducing "one axis separates and the other does not" cannot be tuned.

**Verify each protocol against a closed-form oracle first** (`aleph/validation/oracles/`, read-only).
A protocol that cannot reproduce its own analytic limit will not be believed on a cell.

---

## Lane 3 — `nn-interpret`

```
export FFN_SESSION=nn-interpret
```

Owns: a new `aleph/inference/**` and `aleph/tests/ac/inference/**`. Keep out of
`aleph/scripts/**` (held whole by `scripts-legacy`) — agree a `nn_*.py` prefix carve-out first.

**Read this before assuming you are blocked.** `virtual_cell/surrogate.py` opens with "Training a
neural surrogate on simulation data is LOCKED: no accepted native trace exists, so any surrogate would
be trained on force-accepted trajectories." That lock is correct and it is **narrow**: it covers
learning from THIS project's simulation output. A network trained on EXPERIMENTAL data does not touch
an accepted trajectory and is not what that sentence forbids. The distinction is not written down
anywhere yet — writing it down is part of this lane's first deliverable, so the next reader does not
stop at the word "LOCKED".

**The data that exists.** `aleph/references/260313_{Bare,Pre,Lam4}.csv` — 4,554 rows total across
three conditions, time-resolved from segmentation masks at 60-minute intervals, columns
`Frame, Time_min, Area_um2, Perimeter_px, Solidity, Circularity, A0_um2, A_over_A0,
EffectiveRadius_um, Series`. `A_over_A0` is exactly `ROADMAP.md`'s stage-0 adherent observable
("spreading area A, top-down silhouette, the PI-ratified metric").

**Establish provenance before training on it.** `EffectiveRadius_um` runs 227-301 um and `Area_um2`
16-28 x10^4 — far too large for a single cell, so these are very likely spheroids or colonies. If so
they are a Layer-2 collective observable, not stage-0 single-cell spreading, and the lane's target
changes accordingly. `Bare`/`Pre`/`Lam4` read as substrate coatings (laminin?) but that is a guess.
Resolve it from `aleph/references/*.md` or by asking the PI. **Do not train until the object being
measured is named** — this is the same failure as computing a modulus from a state that did not
contain what the operator needed.

---

## Lane 4 — `nn-structure`

```
export FFN_SESSION=nn-structure
```

Owns: `aleph/inference/structure/**`. Coordinate with Lane 3 on the shared package root.

**Goal.** A network that infers compartment/connector structure from experimental data alone — which
compartments must be present, and how coupled, to produce the observations.

**Start by proving the data is insufficient, because it currently is.** One spreading time series
across three conditions cannot distinguish 14 compartments. The honest first deliverable is not a model
but a **requirement**: given the composed graph (14 components / 36 connectors across 13 families,
which `ac_composed_world_dump.py` builds and validates today), what observations would identify it, and
which of those exist? `CELLTYPE_MODALITY_MATRIX_2026-07-29.md` shows 22 cell types, of which 19 carry
only ONE measurement method — so cross-modality identification is possible today for the breast triplet
and, weakly, the non-adherent trio, and nowhere else.

**The lever this lane has that others do not.** Structure is constrained by what a model CANNOT
express. `virtual_cell/archetypes.py` already records the pattern: an `active rod/cable graph` has zero
in-plane shear modulus **by construction**, so an RBC observable is not merely unmeasured but
uncomputable under that architecture. Enumerating those impossibilities across the 22-type panel is a
structural constraint that needs no training data at all, and it is the right first move.

---

## What none of the four are blocked on

Stages 1-5 cannot move until something converges: `SOLVE_COUPLED` has zero production tasks and
`composed_native.py::solve_candidate` is one force-assembly pass with no iteration. **All four lanes
above are outside that block** — Lane 1 is documents and KB, Lanes 3/4 train on experimental data, and
Lane 2 can build and oracle-verify protocols before any converged trajectory exists. That is the
argument for opening them now rather than waiting.

## Still PI decisions

1. **`ROADMAP.md`'s K=2.** The matrix supports keeping MCF7 + MDA-MB-231 as the stage-5 library (the
   only lineage with multiple measurement methods) while adopting the other 20 types as a
   VALIDATION panel. Not the same as widening K, and the two should not be conflated.
2. **Which mechanical modality to build first** — the suggestion above is from the matrix, not a ruling.
3. **The gate amendment** (convergence vs stationarity), still held, with its stage-1 sentence still
   defective (`balance_ok_d` cannot fail for non-convergence).
4. **Boot budget 46,000 -> 48,000 B** — now 30 bytes from the ceiling.
