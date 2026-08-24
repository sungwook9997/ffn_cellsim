# ALEPH-PORT-2501 — sf_arc as an explicit active rod/cable graph

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-2501` |
| Lane | `L25 physics — sf_arc mechanics` |
| Status | `PROPOSED` |
| Written | `2026-07-30` — **after the code, which is a deviation from PLAN §0.2.5, recorded in §9** |
| Port class | `RE-DERIVED` |
| Verdict | **RE-DERIVE.** Zero lines, zero identifiers, zero constants taken. Verified by symbol search, not asserted. |
| Authorises | `aleph/vertical/sf_arc.py` |
| Controls | `tests/vertical/test_sf_arc_controls.py` — 34 tests, mutation-checked |

---

## 1. Aleph API

```python
from aleph.vertical.sf_arc import (
    SFArcGraph, Bundle, Subpopulation, MaterialCard, ElementType,
    InternalJoint, AttachmentSite, AttachmentRole, MaterialPoint, MaterialPointSink,
    axial_energy_and_forces, bending_energy_and_forces, joint_energy_and_forces,
    assert_no_cable_carries_compression, assert_state_keys_disjoint_from,
    owned_state_keys, default_material_cards, straight_bundle_nodes,
)
```

This authorises the **mechanics** of the `sf_arc` owner. The *contract* — one `E` owner and its four
`I` subpopulations — was authorised separately by `ALEPH-PORT-1702` and is not re-litigated here. The
material-card numbers are placeholders and no claim is made about their values; see §8.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY, never modified) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Nearest source paths | `ffn_sim/ac/weave/stress_fiber.py` (462 lines), `ffn_sim/ac/weave/regions.py`, `ffn_sim/ff/architecture_spec.py` |
| Lines taken | **0** |

**The verdict is checked rather than claimed.** Every public identifier this module defines was
searched for across the whole reference tree:

| Identifier | Reference files containing it |
|---|---|
| `SFArcGraph`, `MaterialCard`, `InternalJoint`, `AttachmentSite`, `MaterialPointSink` | 0 |
| `axial_energy_and_forces`, `bending_energy_and_forces`, `straight_bundle_nodes` | 0 |

## 3. Why RE-DERIVE, and what the reference actually is

The reference has no `sf_arc` **owner**. `ac/weave/stress_fiber.py` is a *load-path extraction and
wiring* module: `extract_sf_bundles(cell)` reads bundles back out of an existing woven cell,
`stress_fiber_load_path` walks a path that already exists, and `wire_stress_fibers` attaches motor
reactions to it. It is a reader over someone else's state, not a state owner with its own arrays,
its own energy and its own transaction.

That is a different object from what the docs specify, so there was nothing to port even in
principle. The design here — one owner, four subpopulations as *labels* rather than as separate
owners, private node arrays, attachment by material coordinate — comes from the registered contract
in `aleph/state/census_actomyosin.py`, not from the reference.

**One argument crosses, and no code.** The audit of the reference cortex owner
(`ac/engine/cortex_population.py`) found a **zero-energy axial mode per filament**: bending was
bound, an axial law was not, and a second-difference bending energy vanishes for *any* uniform node
spacing including collapse to a point. Inextensibility was deferred to a constraint projector the
owner never referenced, so the reference's own buckling test passed only because the test called
that projector directly.

This module does not have that defect, and the claim is now measured rather than reasoned:
`TestTheAxialLawIsActuallyBound` uniformly contracts a straight rod by 20% and requires the energy to
notice. It also asserts that the *bending* term is blind to that contraction — confirming that the
axial law is the only thing standing between this owner and the reference's failure, and that it is
actually bound.

## 4. What the controls establish

`tests/vertical/test_sf_arc_controls.py`, 34 tests, all passing.

| Claim | How it is checked |
|---|---|
| `F = -grad E` for the axial, bending and joint terms **separately** | central finite difference, observed convergence order > 1.6, relative residual < 1e-6 |
| The terms sum to the reported total, and the total gradient matches | independent FD on `potential_energy_pn_um` |
| A cable never carries compression | the slack branch returns **exactly** `0.0`, testable with `==` |
| The isolated graph injects no net force **or moment** | residual over *constituent* scale < 1e-12 |
| The active channel is internally self-equilibrated | net active force over constituent scale < 1e-12 |
| Contraction does positive active work; extension negative | signed work, both directions |
| Axial stiffness is resolution-independent | identical end-to-end strain at 3/5/9/17 nodes, spread < 1e-12 |
| A material point resolves to the arclength it names | interpolated position, and **moment** conservation on scatter |
| Rollback restores positions bit-identically and discards proposed topology | `np.array_equal`, not `allclose` |
| A role may not attach to a subpopulation the registry withholds | LINC on a ventral fibre is refused |

**The FD step is chosen above the round-off floor.** Steps are 4e-4, 2e-4, 1e-4 µm, all well above
`eps^(1/3)`. `PLAN.md` §6.1 records the afternoon this was got wrong: below the floor the measured
error is cancellation rather than truncation, the apparent order goes negative, and a correct
gradient looks broken.

### 4a. The controls were checked for the ability to fail

A passing test proves nothing unless it can fail. Four mutants were introduced into the shipped
module and each was caught:

| Mutant | Tests failed |
|---|---|
| bending energy coefficient `0.5 → 0.4` | 2 |
| bending gradient divides by the wrong segment length (`lv → lu`) | 2 |
| bending middle-node force sign flipped | 5 |
| rollback restores `positions * (1 + 2^-40)` instead of exactly | 1 |

The source was restored and verified identical to `HEAD` after each. The module additionally ships
two deliberate breaks as flags — `flip_axial_sign` and `allow_cable_compression` — so the negative
controls drive the *shipped* code path rather than a copy of it. Both defaults are `False` and no
production caller sets either.

## 5. Units, domains, singular cases, invariants

**Units.** Length µm, force pN, energy pN·µm, time s, throughout and without conversion. Axial modulus
`axial_modulus_pn` is a stiffness times a length (pN), so the per-segment stiffness is
`axial_modulus_pn / L0` and the same number means the same material at any discretisation — this is
the resolution-independence the control measures, not an assertion. Bending rigidity is pN·µm².

**Domains.** Every rest length must be `> 0` and every Voronoi length must be `> 0`; both are checked
and raise rather than clamp. A material coordinate must lie within its bundle's rest length; a site
beyond the end is refused rather than clamped, because a clamped attachment silently slides to the tip
and reports a plausible force from the wrong place.

**Singular cases.** A collapsed segment (zero length) raises, since the unit tangent is undefined
there. A straight bending triple gives **exactly** zero energy and exactly zero force, not a small
number — the second term is `|t₂−t₁|²` and both tangents are identical. A slack cable gives **exactly**
zero, because `dE/dL` is identically zero on the whole slack branch rather than merely small near the
gate, which is what makes it testable with `==`.

**Invariants, each with a control.**

| Invariant | Control |
|---|---|
| The internal force has zero resultant and zero moment | residual / constituent scale < 1e-12 |
| The active force has zero resultant | residual / constituent scale < 1e-12 |
| Bundle node ranges tile `[0, N)` with no gap or overlap | `assert_filament_ownership_is_unique` |
| Per-subpopulation node counts partition the nodes | census sums to `positions.shape[0]` |
| Every state key is `sf_arc_`-prefixed and unique | `owned_state_keys` |
| Rollback restores state bit-identically | `np.array_equal`, not `allclose` |
| Bending energy is invariant under uniform contraction | asserted at `abs=1e-9` |
| End-to-end axial stiffness is independent of node count | spread < 1e-12 across 3/5/9/17 nodes |

## 6. Positive control

The **positive control** is the finite-difference gradient agreement, run per energy term separately
rather than on the sum: a summed check is exactly the one that a missing or doubled term survives,
because three wrong gradients can still add up correctly. Each term must first show `energy > 0.0` in
the test configuration — otherwise the comparison is between two zeros and passes for a term that does
nothing at all. Measured: relative residual < 1e-6 with an observed convergence order above 1.6 against
a theoretical 2.

A second positive control is analytic rather than numerical: for a straight rod the bending energy must
be exactly zero, and in the small-angle regime doubling the transverse offset must quadruple the
bending energy. Measured **3.9988** and **3.9952** at offsets of 0.01→0.02→0.04 µm — approaching 4 from below, which is the tangent law's saturation still faintly present at these amplitudes and is the direction it must approach from.

## 7. Deliberately failing negative control

Three, and each one must fail:

1. **`flip_axial_sign=True`** negates the axial force so a stretched cable pushes. The gradient test
   must then disagree by more than 50% of the force scale. If it does not, the positive control is not
   sensitive to the axial force at all and proves nothing about it.
2. **`allow_cable_compression=True`** lets a one-sided element carry compression. The graph is squeezed
   to half its rest length and `compressive_load_in_cables_pn()` must report a non-zero load, after
   which `assert_no_cable_carries_compression` must raise. The intermediate assertion is deliberate:
   `PLAN.md` §6.1 records a reporter that returned `0.0` from an unpopulated cache and made a
   deliberately broken strut look innocent, so the detector is required to *see* the defect before it
   is asked to refuse it.
3. The **same configuration on the default path must be clean**, or the detector fires on everything
   and its refusal carries no information.

Both breaks are shipped flags rather than hand-edited copies, so the control exercises the real code
path. A negative control run against a duplicate of the code only proves the duplicate is broken.

## 8. What is NOT established

- **No material-card value is evidence of anything.** `axial_modulus_pn`, `bending_rigidity_pn_um2`
  and `active_tension_pn` are placeholders chosen to make the algebra sharp. Quantitative status is
  `BLOCKED`; citation status is `UNSOURCED`.
- **No connector is wired.** This module provides the endpoint machinery — `AttachmentRole`,
  `MaterialPointSink` — that nine declared connectors need. It wires none of them. The wired count
  contributed by this entry is **0**.
- **Not integrated into `aleph/vertical/assembly.py`.** Deliberately. It has controls now, but it has
  never been stepped inside a world, and the first vertical's Laplace/FD results were obtained on the
  elastic shell cortex and do not transfer.
- **No dissipation, no thermal channel, no crosslink kinetics, no turnover.** The active tension is a
  stand-in for an NMII connector this owner does not have yet.
- **The bending law is a unit-tangent law**, `E = (κ/2h)|t₂−t₁|²`, not a second-difference law. It is
  quadratic in the turning angle only at small angles — the measured energy ratio at a 0.4 µm offset
  on a unit baseline is 3.59, not 4. This is the law saturating, not an error, but the two laws
  disagree about a sharply bent filament and anything downstream must know which one it has.

## 9. Process deviation, recorded rather than smoothed over

**This entry was written after the code, which PLAN §0.2.5 forbids.** The L25 lane wrote the module
and then died on an API error before writing either the ledger entry or the tests. The module cited
this file by name in its docstring, and the file did not exist — the same shape as the L11 failure in
`PLAN.md` §6.1, where docstrings cited tests that had never been written.

For four hours the module sat in the repository at 1,328 lines, parsing and importing, with no tests
and no ledger entry. Nothing consumed it, which is the only reason this cost nothing. The rule it
broke exists precisely because a module that parses reads as finished.

The tests, the mutation check and this entry were written afterwards by the coordinator session. The
code was **read line by line and not trusted**; what §4 records is what was measured, and §5 records
what was not.
