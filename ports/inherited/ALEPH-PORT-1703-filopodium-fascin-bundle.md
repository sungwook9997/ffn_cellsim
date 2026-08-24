# ALEPH-PORT-1703 — filopodium: a bundled strut whose root is a connector, not a weld

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1703` |
| Lane | `L17 registry — actomyosin load-path compartments` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` — **verdict: RE-DERIVE.** One *discipline* crosses; no code, no identifier, no constant. |

---

## 1. Aleph API

```python
from aleph.state.census_actomyosin import FILOPODIUM
```

Aleph target file: `aleph/state/census_actomyosin.py`. This entry authorises the `filopodium` contract
entry. No bundle mechanics and no crosslinker stiffness is authorised, because none exists.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/protrusion.py`, `ffn_sim/ac/weave/regions.py` |
| Source symbol(s) | `FilopodiumFascinBundleMechanics`, `FilopodiumStateOwner`, `_fascin_bundle_kernel`, `fascin_bundle_pair_reference`, the three fascin geometry/stiffness module constants, `_filopodium_arch` |
| Read from | **`git show be0e5876:<path>`** — the immutable commit. |
| Working tree == commit? | `yes` — `git diff be0e5876 -- ffn_sim/ac/engine/protrusion.py` and `-- ffn_sim/ac/weave/regions.py` are both empty. |

| Path | `sha256:` | Lines |
|---|---|---|
| `ffn_sim/ac/engine/protrusion.py` | `sha256:98300c76d61b3a6d` | 1107 |
| `ffn_sim/ac/weave/regions.py` | `sha256:091ac7ea8d1a0645` | 598 |

A second implementation exists in a parked subtree of that repository, with its own kernels and a hard
force cap. It is out of that project's own import closure and was not audited beyond confirming it is
parked; nothing from it crosses and nothing here depends on it.

## 3. Why source-derived porting beats clean-room

**RE-DERIVE.** The fascin crosslink there is a Hookean pair force, which is the least portable thing
imaginable — any competent author writes it in four lines, and `ports/TEMPLATE.md` §3 is explicit that a
convention anyone would reach independently is a reason to write it clean-room, not to port it.

What *is* worth carrying is a **discipline**, and it is the best single thing I found in that
repository:

> Two of the three fascin constants — the interfilament spacing and the axial repeat — carry real
> citations and are used. The third, the crosslink stiffness, is a slot with **no default at all**, and
> the mechanics delegate **raises on construction** unless a positive finite value is supplied. Beside
> it, a comment forbids substituting a different crosslinker's stiffness, which is the specific wrong
> answer available.

That is the correct handling of a missing parameter and it is rarer than it should be. The instinct
under schedule pressure is to reach for a nearby crosslinker's value with a comment saying
"provisional", and then the provisional value is in every result. The alternative — refuse at
construction — makes the gap load-bearing: nothing can run until someone decides, so the decision
cannot be deferred by accident.

Aleph's counterpart is that `FILOPODIUM.approximation` names the crosslinker stiffness as owed and
states that a consumer needing one must refuse. That is a weaker guarantee than a constructor that
raises, because a contract cannot refuse on behalf of code that does not exist yet, and this entry says
so rather than claiming parity. **What crosses is the discipline. No constant, no identifier, no
comment.**

## 4. Physical or mathematical law represented

No mechanics was ported. The declarative law is about the root, and it is the one that distinguishes a
filopodium from a relabelled cortical filament:

**A root is a connector, and that is what makes detachment possible.** Suppose instead the filopodial
bundle's basal nodes *were* cortical nodes — shared, not connected. Then:

1. There is no state whose change constitutes detachment. Detachment would have to be a topology edit
   on the cortex's own arrays, performed by something that does not own them.
2. The load at the root is not transferred, it is *already* internal, so no adjoint pair exists and
   force closure at the root is a tautology.
3. Removing the connector to isolate the root load path — the standard negative control — is
   impossible, because there is no connector to remove.

So "the root is a connector" is not stylistic. It is the condition under which the root has state,
carries a falsifiable transfer, and can be switched off in a control. The same argument applies to the
tip contact and to the nascent adhesions along the shaft.

**A bundle is a rod with an effective rigidity, and the effective rigidity is where the crosslinker
hides.** A bundle of `n` filaments crosslinked densely enough to shear-lock behaves as one rod with
rigidity scaling like `n²` times the single-filament value; crosslinked loosely, like `n` times. The
crossover is set by the crosslinker stiffness, which is the unsourced parameter. So the missing
constant is not a detail of one force term — it selects between two scalings that differ by a factor of
`n`, and for a bundle of ten filaments that is a factor of ten in the buckling load. That is why
`FILOPODIUM.unsupported_claims` disclaims buckling load and critical length specifically.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| bundle node position | µm | m | finite |
| interfilament spacing | µm | m | > 0 |
| crosslinker stiffness | pN/µm | N/m | > 0, **unsourced, no value carried** |
| tip elongation velocity | µm/s | m/s | ≥ 0, load-dependent |
| root attachment rate | 1/s | 1/s | ≥ 0, **unsourced** |

Singular cases:

- **A filopodium with no root attachment.** Legal: a retracting or newly nucleated bundle. It carries
  no root load, and the transmitted load must be exactly zero rather than small.
- **A bundle of one filament.** Legal and degenerate: the crosslinker topology is empty and the
  effective rigidity is the single-filament value. A consumer that assumes `n ≥ 2` must say so.
- **A crosslinker stiffness requested with no source.** Refused, never defaulted, and never borrowed
  from another crosslinker.
- **A target filopodium count requested.** Refused; the count is an output.

Invariants:

- **I1.** `filopodium` records the disjointness rule verbatim and owns a name-prefixed filament
  population block. `tests/state/test_census_actomyosin.py::test_each_filament_owner_records_the_disjointness_rule`,
  `test_each_filament_owner_names_its_own_filament_population_block`.
- **I2.** It owns a root attachment state, so detachment is a state change rather than a topology edit
  on somebody else's arrays. Asserted through the owned-state block list in
  `test_no_state_key_is_claimed_by_two_owners`.
- **I3.** Its count is declared an output and at least three rate priors are named as owed.
  `test_each_dynamic_population_records_that_its_count_is_an_output`,
  `test_each_dynamic_population_names_the_rate_priors_evidence_is_owed_on`.
- **I4.** It disclaims buckling load and lifetime.
  `test_no_entry_claims_a_magnitude_it_has_no_source_for`.
- **I5.** No state key it owns is claimed by another owner.
  `test_no_state_key_is_claimed_by_two_owners`.

## 6. Source evidence class and known retractions

Looked in that repository's audit directory, roadmap, and gap-evidence cards.

- `AUDIT_AC_ENGINE_2026-07-25.md:38` rates `filopodium` **`KERNEL_BOUND`, zero CUDA execution,
  "declared-only in practice"**, and `:89` lists it among components never executed. `:50` records the
  filopodial motor connector as `CONTRACTED` with **no bind target, because neither owner is ever
  constructed**.
- `AC_ENGINE_COMPLETION_ROADMAP_2026-07-23.md:53` marks filopodium **`SEAMED`** with activation blocked
  on the crosslinker stiffness.
- Its own compartment-validation track lists "the crosslinker stiffness of `None` must raise" as a gate
  still to be written — so the raise exists in the code and the *test* of the raise did not, at the time
  that document was written.
- Its composed-world artifact hardcodes a filopodium count with a comment saying the population size is
  a declared gap and is deliberately not asserted. That is a placeholder, honestly labelled. It is not
  resampling, and I want to be precise about that: I searched for resampling, rejection sampling and
  count-rescaling machinery for filopodium and adhesion counts and **found none**.

Reachability: the owner, the mechanics delegate and the kernel are **dead code** outside tests. The
kernel is never launched outside its own test module.

Tests: the fascin block is inside a 726-line plumbing module. Two of its filopodium tests are *guard*
tests asserting the raise fires, which is the right kind of test but is not physics. Exactly one is
numeric — an equal-and-opposite check on the host reference pair.

Retractions: none found for the two sourced fascin geometry constants; I looked in the audit directory
and the gap cards. Neither was inherited. The stiffness has nothing to retract because it has no value.

## 7. Independent oracle or derivation

For this entry's declarative content: the same exact ownership identity as `ALEPH-PORT-1701` §7, plus
the exact statement that a detached root transmits exactly zero. Both are discrete, so both are
asserted with `==`.

For the mechanics a later lane would write, Aleph has an oracle the reference did not use: **Euler
buckling with a known exponent.** A slender clamped-free rod of rigidity `B` and length `L` buckles at
`F_c = π²B/(4L²)`, so a bundle's effective rigidity can be *measured* from a buckling sweep and
compared against the `n` versus `n²` scalings in §4 — a convergence-rate check with a known exponent,
which `ports/TEMPLATE.md` §7 names as a legitimate independent oracle. That would decide the
crosslinker regime without needing the crosslinker constant, which is the useful part: it turns the
missing parameter into a measurable rather than a blocker.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_census_actomyosin.py::test_each_filament_owner_records_the_disjointness_rule` | The rule is present verbatim on `filopodium`. |
| Positive | `tests/state/test_census_actomyosin.py::test_each_filament_owner_names_its_own_filament_population_block` | It owns `filopodium_filament_population`, name-prefixed so a collision is detectable. |
| Positive | `tests/state/test_census_actomyosin.py::test_each_dynamic_population_names_the_rate_priors_evidence_is_owed_on` | Four priors are named as owed, including the crosslinker stiffness. |
| Positive | `tests/state/test_census_actomyosin.py::test_no_entry_claims_a_magnitude_it_has_no_source_for` | It disclaims buckling load, critical length and lifetime. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_two_owners_sharing_a_state_key_is_refused` | Giving `filopodium` the lamellipodial filament-population key is refused by name — the declarative form of "the filopodium is a slice of somebody else's array". |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_an_explicit_owner_that_owns_nothing_is_refused` | `filopodium` with an empty `owned_state` is refused. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_a_filament_owner_that_drops_the_disjointness_rule_is_refused` | A filament owner whose `mechanical_role` omits the rule is refused. |
| Negative (must fail) | `tests/state/test_census_actomyosin.py::test_resampling_to_a_target_count_is_named_as_excluded` | The entry must name resampling, rejection and rescaling as excluded; stripping any of the three fails. |

## 10. Numerical and precision envelope

No arithmetic crosses. Controls are exact string and set comparisons, asserted with `==` and `in`; no
tolerance is stated because the properties are discrete and a tolerance would weaken a check that
should not be weakened.

For the record, since it bears on a later lane: a Hookean bundle crosslink at physiological
interfilament spacing is a **stiff** short-range spring, and the explicit stability limit
`dt < 2/(M·k)` scales inversely with it. So the missing crosslinker stiffness is not only a physics gap
— it also sets the step size, which means it sets the cost. A "provisional" value would have silently
set both.

## 11. Production-backend residency and transfer

Host-side declarations, no device residency, no transfer, not in the step loop. When filopodial
mechanics is written, the residency note this entry leaves is that the crosslink is a short-range pair
interaction within one owner's array, so it needs no cross-ownership transfer at all — unlike the root,
which needs a two-array adjoint scatter for exactly the reasons recorded in `ALEPH-PORT-1702` §3.

## 12. Comments and docstrings to discard

Read, and **not** carried: the delegate, owner, kernel and reference-function names; the constant
identifiers; the two citations attached to the geometry constants; the module paths; the gap-card
labels; the status vocabulary; the parked subtree's force-cap value.

The forbidding comment beside the empty stiffness slot is the one piece of prose I was tempted by, and
it is exactly the kind of thing PLAN §0.2 forbids carrying — it references that project's own
parameter names and its own gap identifiers, so it would be unreadable without its repository open.
What replaces it: §3 and §4 above, which state the discipline and the `n` versus `n²` consequence in
Aleph's terms, plus `FILOPODIUM.approximation`, which says the slot stays empty and a consumer must
refuse.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** `PROPOSED`. The contract entry and its controls pass (`tests/state/test_census_actomyosin.py`, 88 passed). The substantive content is the root-is-a-connector derivation in §4 and the discipline in §3, neither reviewed by a human. |
| Reviewer | agent-proposed, **unratified**. |
| Rollback | Delete the `FILOPODIUM` entry and its tests. `validate_group` then fails on a missing filament owner. Nothing outside the module depends on it. |

## 14. Honest limits

- **Unratified.**
- **No filopodial mechanics exists in Aleph.** This is a declaration.
- **Aleph's refusal is weaker than the reference's.** That project refuses at construction; Aleph
  currently refuses in prose, because there is no constructor to refuse in. When a mechanics lane
  lands, it should raise on an unsupplied crosslinker stiffness, and this entry is where that
  obligation is written down.
- **The `n` versus `n²` scaling argument is textbook shear-lag reasoning, restated, and it is not
  measured.** The factor-of-`n` claim is an order-of-magnitude statement about which regime matters,
  not a derived prefactor.
- **The Euler-buckling oracle in §7 is proposed, not written.** No test implements it.
- **"Dead code outside tests" is `AUDIT_READ`.** I searched that tree and read its own audit; I did not
  execute anything there.
- **The parked duplicate implementation was not audited.** It has real kernels and a hard force cap. If
  a filopodium vertical is authorised it deserves its own entry rather than inheriting this verdict by
  silence.
