# ALEPH-PORT-1705 — focal_adhesion: a geometry-less clutch graph, in series, with no state owner upstream

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1705` |
| Lane | `L17 registry — actomyosin load-path compartments` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` — **verdict: RE-DERIVE.** The strongest RE-DERIVE of the five: the upstream state owner does not exist there at all. |

---

## 1. Aleph API

```python
from aleph.state.census_actomyosin import FOCAL_ADHESION, NO_FILAMENT_OWNERSHIP_RULE
```

Aleph target file: `aleph/state/census_actomyosin.py`. This entry authorises the `focal_adhesion`
contract entry. No clutch kinetics, no bond law, no maturation rate is authorised.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/solid/adhesion_clutch.py`, `ffn_sim/ac/engine/load_path.py`, `ffn_sim/ff/fa_clutch_warp.py`, `ffn_sim/ff/fa_maturation.py`, `ffn_sim/ff/fa_ecm.py` |
| Source symbol(s) | `AdhesionClutchCompartment`, `attwood_slip_off_rate`, `LoadPathJointRuntime`, `ActorSegmentRuntime`, `clutch_spring_kernel`, `clutch_catchslip_kmc_kernel`, `resolve_clutch`, `talin_unfold_rate`, `fa_growth_rate`, `k_int_effective`, `attach_clutches_to_ecm` |
| Read from | **`git show be0e5876:<path>`** — the immutable commit. |
| Working tree == commit? | `yes` — `git diff be0e5876 -- ffn_sim/ac/solid/adhesion_clutch.py`, `-- ffn_sim/ac/engine/load_path.py`, `-- ffn_sim/ff/fa_clutch_warp.py`, `-- ffn_sim/ff/fa_maturation.py` are all empty. |

| Path | `sha256:` | Lines |
|---|---|---|
| `ffn_sim/ac/solid/adhesion_clutch.py` | `sha256:80eb4ea0fd90c9dc` | 118 |
| `ffn_sim/ac/engine/load_path.py` | `sha256:b117c9ac503fb97c` | 785 |
| `ffn_sim/ff/fa_clutch_warp.py` | `sha256:2db8028d978d2a1e` | 108 |
| `ffn_sim/ff/fa_maturation.py` | `sha256:04aeb73dd88c2cef` | 134 |

## 3. Why source-derived porting beats clean-room

**It does not, and this is the clearest RE-DERIVE of the five.** The reason is not code quality. It is
that the assets are in the wrong places relative to each other, in three specific ways.

**There is no focal-adhesion state owner in that engine at all.** Its own audit says so and my search
concurs: the component is declared in a contract table and has no owner class anywhere. So there is no
implementation of this compartment to port. What exists is three disconnected fragments:

1. A **compartment class of 118 lines** whose accumulate path launches a two-sided Hookean spring — and
   which owns a catch–slip off-rate helper that its accumulate path **never calls**. So the "stochastic
   clutch" in that class is deterministic: the rate law is present as a function and absent as
   behaviour. That is the exact declared-but-unevaluated shape this project was restarted over,
   appearing at function granularity rather than at connector granularity.
2. A **series-joint runtime of 785 lines** with four real kernels — which is a *linear spring* whose own
   docstring concedes that candidate transitions arrive from a caller, so it has no clutch kinetics
   either. Test-only.
3. The **only real clutch kinetics and the only maturation code in the repository** — a catch–slip
   kinetic Monte Carlo kernel, and a talin-unfolding-plus-vinculin maturation module with a Hill growth
   term — both of which live in a different subtree that the engine imports **zero times**. Real code,
   with real closed-form tests, on no path that produced any result.

So the honest description is: the good kinetics is unreachable, the reachable code has no kinetics, and
the compartment has no owner. Porting any one fragment would import the disconnection along with it,
and porting the set would be porting an architecture — which master plan §5.3 forbids outright.

**Two further blockers, both fatal on their own.** The α2β1–collagen rate law has **no kinetics anywhere
in that tree**, by its own gap card; and the catch–slip parameterisation used as a default has **no
resolvable citation**, by its own audit. A clutch with no bond law is not a clutch, and a bond law
whose provenance does not resolve cannot be inherited under PLAN §0.2.

Nothing crossed. One *argument* is worth recording, and it is in §4.

## 4. Physical or mathematical law represented

No kinetics was ported. The law derived here is the one that makes "not two independent springs" a
quantitative statement rather than a style note, and it is the reason the registry declares the two
halves a single mechanical group.

**Series compliances add; parallel stiffnesses add.** For the path
`actin/SF → talin/integrin state → collagen ligand`, let the actin-side element have stiffness `k₁` and
the ligand-side element `k₂`. The two are **in series**: the same load `F` passes through both, and the
total extension is the sum of the two extensions. So

```
δ_total = F/k₁ + F/k₂        ⇒        k_series = (1/k₁ + 1/k₂)⁻¹ = k₁k₂/(k₁+k₂)
```

Evaluating them as two independent springs each attached to ground computes `k₁ + k₂` instead. The ratio
of the wrong answer to the right one is

```
(k₁ + k₂) / (k₁k₂/(k₁+k₂)) = (k₁ + k₂)² / (k₁k₂) ≥ 4
```

with equality only at `k₁ = k₂`. So the error is **never smaller than a factor of four**, and it grows
without bound as the stiffnesses separate — at `k₂ = 100k₁` it is 102×. Two consequences that matter more
than the number:

- The wrong model is **always stiffer**, never softer. So the error has a fixed sign, and a fitted
  prefactor absorbing it would have to be less than one quarter — which means the misfit is not
  absorbable without a prefactor so far from unity that it would itself be a finding.
- The **softer element dominates** the true series stiffness (`k_series → k₁` as `k₂ → ∞`). So the
  compliant half of the path sets the compliance, and that is precisely what a clutch model is *for*:
  the whole point of a molecular clutch is that its compliance and its rupture, not the actin's, set the
  traction. Evaluating in parallel throws away the mechanism the model was built to represent.

**A geometry-less representation forbids a geometric claim, and forbids it structurally.** There is no
adhesion mesh, no adhesion volume, no plaque shape. So adhesion area is not a small-error quantity
here — it is not a quantity at all, and any number reported as an adhesion area would be manufactured
from something else. `FOCAL_ADHESION.unsupported_claims` names area, aspect ratio and spatial growth for
that reason.

**An index into another owner's array is not a reference.** The actin-side state must reference an actin
owner's material point *through a connector*. In the reference implementation the actin side is a global
integer index into a whole-cell position array, and the accumulate path indexes straight into it. That is
ownership spelled differently: the adhesion writes into arrays it does not own, with no adjoint pair and
no connector to remove in a control. `FOCAL_ADHESION.mechanical_role` states that the reference is not
ownership, and `validate_group` requires the no-filament-ownership rule on the entry.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| actin-side stiffness | pN/µm | N/m | > 0, **unsourced** |
| ligand-side stiffness | pN/µm | N/m | > 0, **unsourced** |
| series stiffness | pN/µm | N/m | derived, `≤ min(k₁, k₂)` |
| bond binding rate | 1/s | 1/s | ≥ 0, **unsourced, no law exists** |
| load-dependent off-rate | 1/s | 1/s | ≥ 0, **unsourced** |
| occupancy | 1 | 1 | [0, 1], an **output** |
| clutch count | count | 1 | ≥ 0, an **output** |

Singular cases:

- **An unbound clutch.** Legal, and it carries **exactly zero**, not a small residual. This is the same
  unilateral-exactness argument as `ALEPH-PORT-1103` I1: the inactive branch has identically zero
  derivative, so the control is an equality against zero and a tolerance would let a small wrong
  traction through in exactly the place that matters.
- **One side of the series unbound.** The series stiffness is zero, not `k` of the bound half. A
  parallel evaluation would report the bound half's stiffness — a load path that is broken reporting a
  load. That is the failure mode §4 quantifies.
- **A zero stiffness on either side.** `k_series = 0` exactly, by the formula, with no division by zero
  if the reciprocal form is avoided. A naive `1/k₁ + 1/k₂` implementation would raise; the product form
  `k₁k₂/(k₁+k₂)` does not, and that is the numerically conditioned form to use.
- **An adhesion area requested.** Refused. It does not exist in this representation.
- **A target adhesion count requested.** Refused; the count is an output.

Invariants:

- **I1.** The representation is declared geometry-less and names the absence of a mesh.
  `tests/state/test_census_actomyosin.py::test_focal_adhesion_is_geometry_less`.
- **I2.** The series path is declared verbatim, and both wrong readings — passive solid block, and
  independent spring — are named as excluded.
  `test_focal_adhesion_declares_the_series_path_and_refuses_the_two_wrong_readings`.
- **I3.** Both sides of the clutch plus occupancy, maturation and accepted binding topology are owned;
  no actin filament population is.
  `test_focal_adhesion_owns_both_sides_of_the_clutch_but_no_actin`.
- **I4.** The entry states that a reference to an actin material point is not ownership. Same test.
- **I5.** Its count is declared an output with four rate priors named as owed.
  `test_each_dynamic_population_records_that_its_count_is_an_output`.
- **I6.** No state key it owns is claimed by another owner.
  `test_no_state_key_is_claimed_by_two_owners`.

## 6. Source evidence class and known retractions

Looked in that repository's audit directory, its gap-evidence cards, and its test tree.

- Its component audit rates `focal_adhesion` at its **lowest wiring rung — declared-only** — and states
  outright that no state-owner class exists anywhere in that engine, and that the series-joint runtime is
  a linear spring, test-only, with no clutch kinetic Monte Carlo.
- Its blocking gap is recorded as the α2β1–collagen rate law: the declared chemistry card **has no
  kinetics anywhere in the tree**.
- Its gap cards record the compartment and its clutch connector as `SEAMED` — series joint kernel-bound,
  catch-bond rate law for the collagen ligand **absent**.
- Its own audit records **no resolvable citation** for the catch–slip parameterisation used as a default,
  under both the motor and the adhesion use.

Reachability: the compartment class is **dead** outside its own tests. The series-joint runtime is
**test-only**. The clutch kinetic Monte Carlo and the maturation module are on **no engine path at all**
— the engine imports the maturation module zero times.

Tests: small and honest. Four tests on the compartment class, two of them gated on a device being
present; three on the clutch kernels including a kinetic-Monte-Carlo-against-analytic comparison; real
closed-form checks on the maturation module. These are the *right* tests. They are all attached to code
that no result came through, which is the whole finding.

**Retraction in substance, found:** the traction magnitude is recorded `INVALID` and the compartment
`SEAMED` in that project's own gap cards, and the bond law's citation does not resolve. Nothing numeric
was inherited, so there is nothing here that could be retracted numerically.

## 7. Independent oracle or derivation

For this entry's declarative content: the series identity in §4 is exact algebra, and the `≥ 4` bound is
a two-line consequence of the arithmetic–geometric mean inequality. Re-derived above rather than quoted,
and it needs no reference.

For the kinetics a later lane would write, Aleph's available oracles are all closed-form and none is a
comparison against that repository:

- **The series identity itself**, measured: apply a known load across an assembled two-element joint and
  read the total extension. It must equal `F/k₁ + F/k₂` exactly, and the parallel-evaluation bug fails
  this by at least 4×, which is enormous compared with any round-off floor. That is the negative control
  a future lane owes.
- **The exact zero on the unbound branch**, asserted with `==`.
- **Bell-law endpoints**: the off-rate reduces to the bare rate at zero load, exactly.
- **Detailed balance at zero load** fixes the ratio of binding to unbinding rates exactly, independent
  of both magnitudes — so the *ratio* is testable even while both rates are unsourced, which is the
  most useful thing available given the gaps.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_actomyosin.py::test_focal_adhesion_is_geometry_less` | The representation declares itself geometry-less and names the absent mesh. |
| Positive | `tests/state/test_census_actomyosin.py::test_focal_adhesion_declares_the_series_path_and_refuses_the_two_wrong_readings` | The series path text is present, and both wrong readings are named as excluded. |
| Positive | `tests/state/test_census_actomyosin.py::test_focal_adhesion_owns_both_sides_of_the_clutch_but_no_actin` | Five clutch state blocks are owned, no filament population is, and the reference-is-not-ownership statement is present. |
| Positive | `tests/state/test_census_actomyosin.py::test_non_filament_participants_say_so` | Both non-filament participants carry the no-ownership rule. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_a_participant_that_drops_the_no_filament_rule_is_refused` | Replacing `focal_adhesion.mechanical_role` with the plausible-looking "A series clutch." is refused. That mutation is realistic: it is *true*, it is shorter, and it silently drops the statement that makes an actin index not be ownership. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_dropping_the_count_is_an_output_rule_is_refused` | Replacing `focal_adhesion.approximation` with "Geometry-less clutches." is refused; the count rule is mandatory on the entry. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_an_entry_with_no_stated_limits_is_refused` | An entry with no `unsupported_claims` is refused — the five disclaimers here include the parallel-spring one, which is the whole of §4. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_two_owners_sharing_a_state_key_is_refused` | The shared-state-key refusal, which is what would fire if a future adhesion lane declared an actin block as its own. |

## 10. Numerical and precision envelope

No arithmetic crosses this boundary. Controls are exact string and set comparisons, and no tolerance is
stated because the properties are discrete.

Two precision points are recorded for the lane that writes the mechanics, both derived above rather than
inherited. First, **use the product form** `k₁k₂/(k₁+k₂)` and not the reciprocal sum: the product form is
well-conditioned as either stiffness goes to zero and needs no special case, while the reciprocal form
raises. Second, **the unbound branch must be exactly zero and asserted with `==`**, not with a
tolerance: the derivative is identically zero on the whole branch rather than merely small near the
gate, so `==` is available, and the only reason to prefer a tolerance would be to let something through.
Both points are the same lesson `ALEPH-PORT-1103` recorded for the unilateral tether, arrived at
independently for a series joint.

## 11. Production-backend residency and transfer

Host-side declarations; no device residency, no transfer, not in the step loop.

The residency note for a future adhesion lane, and it is the sharp one: because the joint is
**geometry-less**, its state is a small graph rather than a mesh, so the natural residency is a compact
device-side clutch table — and the actin side is a *connector endpoint reference*, which must be
resolved through the connector rather than by indexing the actin owner's array directly. Indexing
directly is cheaper and is exactly the shortcut that turns a reference into ownership. Recording it here
means the cost of doing it properly is budgeted rather than discovered under pressure.

## 12. Comments and docstrings to discard

Read, and **not** carried: every compartment, runtime, kernel and rate-function name; every module path;
the evidence-rung vocabulary; the gap-card identifiers; the two sourced clutch constants and the
attribution on them; the catch–slip default parameterisation and its unresolvable citation; the
force-cap value in the parked subtree.

The docstring that concedes candidate transitions arrive from a caller is the one worth naming: it is
*honest*, and it is also unreadable without that repository open, because "the caller" is a structure
Aleph does not have. What replaces it: §4 above, which derives why the two halves are one group from the
series algebra and puts a lower bound of 4× on the error of getting it wrong — a statement anyone can
check against Aleph's own code with nothing else on disk.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** `PROPOSED`. The contract entry and its controls pass (`tests/state/test_census_actomyosin.py`, 88 passed). The substantive content is the series derivation and its `≥ 4×` bound in §4, and the three-way disconnection finding in §3; neither reviewed by a human. |
| Reviewer | agent-proposed, **unratified**. |
| Rollback | Delete the `FOCAL_ADHESION` entry, `NO_FILAMENT_OWNERSHIP_RULE`'s second consumer, and the tests. `validate_group` then fails on the missing participant. Nothing outside the module depends on them. |

## 14. Honest limits

- **Unratified.**
- **No adhesion mechanics exists in Aleph.** This is a declaration and nothing more.
- **Every rate is `UNSOURCED`, and one of them has no law at all.** The α2β1–collagen bond law is not
  merely unsourced in Aleph — no formulation of it was found in the reference either. So a future lane
  faces a literature question, not a porting question, and this entry cannot reduce that.
- **The `≥ 4×` bound is exact algebra about a modelling error, not a measurement of anything.** It says
  what the parallel-evaluation bug costs; it does not say that anyone made it. I did not measure the
  reference's series joint against its own oracle.
- **The three-way disconnection finding is `AUDIT_READ` plus a search.** I read the files and searched
  for imports; I did not execute that tree. That the engine imports the maturation module zero times is
  a grep result, which is strong for absence and says nothing about whether some dynamic path exists.
- **The geometry-less choice is the registry's and I did not evaluate it.** A geometry-less adhesion
  cannot produce an adhesion area, and adhesion area is a very common experimental observable. Whether
  that is an acceptable scope limit is a decision, not a fact, and it belongs to the PI. Recorded here as
  a question rather than settled by this lane.
