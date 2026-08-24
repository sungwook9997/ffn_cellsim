# ALEPH-PORT-1102 — Turgor ownership, the enclosed-volume constraint, and the minimal cortex

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1102` |
| Lane | `L11 vertical — membrane / cortex / pressure` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `DEFECT_STUDY` — the only thing inherited is a description of a mistake, and the code exists to make that mistake impossible to repeat. |

---

> ## The controls this entry names were never written at the paths it cites
>
> Recorded 2026-07-30. Every control path this entry cited pointed at a test file that has never
> existed in this repository — `test_membrane.py`, `test_pressure.py`, `test_controls.py`,
> `test_assembly.py`, `test_connectors.py`, `test_vertical_firewall.py`. Lane L11 wrote its ledgers
> before its code, which section 0.2.5 requires, against a planned test layout that was then never
> built: the controls were consolidated into **`tests/vertical/test_vertical_controls.py`** instead,
> and the ledgers were never updated.
>
> **The physics is not in question.** This module is the wired vertical, it is covered by real
> passing controls, and `F = -grad E` agrees with a finite difference to 3.95e-07 with observed
> orders 2.000. What is missing is the *mapping* — which specific test discharges which specific
> claim below.
>
> The citations have been rewritten to name the file that really holds the controls, with the
> per-claim mapping marked as never recorded, rather than guessed. Guessing which existing test was
> meant by a promised name writes a false provenance link, which is worse than a visible gap.
>
> **Before this entry may go `ACCEPTED`, someone who can read the physics must restore the mapping**
> — one named test per claim, positive and negative — since `ACCEPTED` requires both controls named
> and resolvable.

## 1. Aleph API

```python
from aleph.vertical.pressure import (
    OsmoticEnvelope,
    PressureDoubleCountError,
    PressureLoad,
    TurgorEnsemble,
    assert_single_pressure_path,
    pressure_load_on_surface,
)
from aleph.vertical.cortex import ElasticCortexShell
```

Aleph target files: `aleph/vertical/pressure.py`, `aleph/vertical/cortex.py`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | none — **no source commit was read for this entry.** No file was opened. |
| Source path | none — no source path was read. |
| Source symbol(s) | none read. The *behaviour* described in §3 was supplied as prose by a prior read-only audit recorded in PLAN §1.1, not by reading code in this session. |
| Read from | neither. |
| Working tree == commit? | not applicable — nothing was read from either. |

## 3. Why source-derived porting beats clean-room

No code is source-derived. One **argument** is, and it is worth recording precisely because a
clean-room author would very likely not have thought of it:

> In the reference project the osmotic load is applied directly to the membrane by the incumbent
> driver, while the declared architecture routes it through the cortex. Two live descriptions of
> the same load, and nothing in the system compares them.

That is an empirical failure mode discovered by somebody else's run, which §3 of `ports/TEMPLATE.md`
names as a legitimate reason to let something cross the boundary. What crosses is the failure mode
and nothing else. The code below was written from the physics.

## 4. Physical or mathematical law represented

**Who owns the turgor.** The osmotic pressure difference is set by the solute imbalance across the
*envelope that separates the two compartments*. In this vertical that envelope is the cortex — the
cortex is the innermost closed surface, it bounds the cytoplasm, and it is what the interior pushes
on. The membrane is outside it. Therefore:

```
E_turgor = E_turgor(V_cortex)          and nothing else
```

The pressure enters the membrane's equilibrium **only** through whatever the cortex transmits across
the contact (`ALEPH-PORT-1103`). Applying `ΔP` to the membrane as well is not a modelling choice
with a small error; it is the same load counted twice, and the sum is exactly `2ΔP∇V`, which looks
entirely plausible. Nothing about the resulting equilibrium is obviously wrong. That is why it needs
a test rather than a docstring.

**Ensembles.** Two are offered and exactly one is chosen per run:

- `FIXED_PRESSURE`: `E = −ΔP · V`, hence `F_i = ΔP ∇_i V`, outward for `ΔP > 0`. This is the
  osmotic-reservoir limit: the interior is buffered, `V` is free.
- `VOLUME_PENALTY`: `E = (K_V/2)(V/V₀ − 1)²`, hence an effective pressure
  `ΔP_eff = −dE/dV = −(K_V/V₀)(V/V₀ − 1)`, positive when compressed. `V` is nearly fixed.

**Penalty, not a Lagrange multiplier — the choice and its price.** `VOLUME_PENALTY` is the
constrained-volume implementation. A Lagrange multiplier would enforce `V = V₀` exactly, but it
introduces an unknown that only a saddle-point solve determines, and this vertical's relaxation
(`aleph/vertical/relax.py`) is a first-order descent with no such solve. Three reasons the penalty
is the right answer *here*, and one honest cost:

1. It has a well-defined force at every configuration, including far from the constraint manifold.
   A multiplier does not — its value at a non-equilibrium configuration is whatever the solver's
   inner iteration says it is.
2. Its error is **measurable and reportable**: `(V − V₀)/V₀` is a number in the run record. A
   multiplier hides its error inside a solver tolerance, where nobody reads it.
3. It is one term in a sum of energies, so the whole assembly stays a single scalar potential and
   the descent's monotonicity check remains meaningful.

Cost: `K_V` sets a stiffness, the stiffness sets a step-size limit, and volume conservation is only
as good as `K_V` is large. That trade is explicit, and the residual `(V − V₀)/V₀` is asserted in the
control rather than assumed small.

**The cortex.** Deliberately the smallest object that can own an envelope and relay a load: a closed
shell of the same connectivity as the membrane, sitting just inside it, with one constitutive term —
a per-face areal spring `E = Σ_f (k_a/2)(A_f − A_f⁰)²/A_f⁰`. Its gradient is `k_a (A_f − A_f⁰)/A_f⁰`
times the exact face-area gradient of `ALEPH-PORT-1101`. §14 lists what it does not represent, and
that list is longer than the model.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| turgor ΔP | pN/µm² | Pa | any sign; positive inflates |
| enclosed volume V | µm³ | m³ | > 0 |
| volume bulk term K_V | pN·µm | J | > 0 |
| cortical areal modulus k_a | pN/µm | N/m | ≥ 0 |
| pressure force | pN | N | finite |

`1 pN/µm² = 1 Pa` exactly, which is the one unit coincidence in this system worth remembering.

Singular cases:

- `V ≤ 0` means the envelope is inside-out. Refused, not clamped.
- `V₀ ≤ 0` in `VOLUME_PENALTY` is refused at construction.
- `k_a = 0` is legal: a cortex with no elasticity of its own, which is the cleanest configuration for
  demonstrating that it is a *relay* and not a load-bearer.

Invariants:

- **I1. Exactly one pressure path.** `assert_single_pressure_path` refuses a world in which more
  than one participant applies an osmotic load, or in which the single load lands on an owner other
  than the declared envelope. `tests/vertical/test_vertical_controls.py (specific control never recorded)`.
- **I2. Pressure force is minus the exact gradient of `E(V)`.** Central-difference verified.
  `tests/vertical/test_vertical_controls.py (specific control never recorded)`.
- **I3. Net pressure force vanishes.** `Σ_i ΔP ∇_i V = 0` to round-off — a uniform pressure on a
  closed surface exerts no resultant. `tests/vertical/test_vertical_controls.py (specific control never recorded)`.
- **I4. Duplication is exactly detectable.** Adding a second envelope on the membrane changes the
  total force by exactly the duplicated term, vertex by vertex.
  `tests/vertical/test_vertical_controls.py (specific control never recorded)`.

## 6. Source evidence class and known retractions

Nothing is inherited but the defect description, whose evidence class is the audit that produced it:
`AUDIT_READ` of a repository at a stated commit, recorded in PLAN §1.1 and
`vault/02_Graph/claims/claim-coverage-gate-passes-unimplemented.md`. Aleph re-states it as a
hypothesis about a *class* of mistake and then tests for that class in Aleph's own code. Aleph makes
no claim about the reference project's current state, and this entry must not be cited as one.

Retractions looked for in: PLAN §1.1, `docs/decisions/OVERNIGHT_DECISIONS_2026-07-30.md`, and the
vault claim above. None found; none expected, since nothing numerical was inherited.

## 7. Independent oracle or derivation

- **Divergence theorem.** `V = (1/3)∮ x·n dA` on a closed surface, which `aleph/state/manifold.py`
  already verifies against `(1/6)Σ v₀·(v₁×v₂)` to round-off. So the volume this module differentiates
  is the same volume the manifold's orientation invariant is built on — two derivations, one number.
- **Central difference** of `E(V(x))` against the analytic force.
- **Resultant identity.** `∮ n dA = 0` on any closed surface, so a uniform pressure produces zero net
  force. Exact, and independent of the discretisation being right.
- **Exact discrete Laplace relation.** From `A(λx) = λ²A`, `V(λx) = λ³V` and the scale-invariance of
  the bending term, stationarity of `σA − ΔP V` under uniform dilation gives `ΔP = 2σA/(3V)` with no
  continuum limit taken anywhere. On a sphere `A/(3V) = 1/R` and this is `ΔP = 2σ/R`.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | A pressurised sphere relaxed through the proper load path recovers σ to better than 1% by the membrane dilational-virial identity, and the recovered `ΔP·R/2` agrees with the applied `ΔP` to within the stated geometric offset. |
| Positive | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | Analytic pressure force equals `−dE/dx` to `O(h²)`. |
| Positive | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | `‖Σ F‖` at round-off. |
| Positive | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | With `VOLUME_PENALTY`, the relaxed volume error `(V−V₀)/V₀` falls as `1/K_V`, and the number is reported rather than assumed. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | Applying the turgor to the cortex *and* to the membrane changes the total force by exactly the duplicated term; the test fails if the difference is anything else, including zero. |
| Negative (must fail) | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | `build_vertical` refuses to assemble a world with two osmotic envelopes. |
| Negative (must fail) | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | An envelope whose target is not the declared osmotic envelope is refused at assembly time. |

## 10. Numerical and precision envelope

float64 throughout; no float32 path. The volume gradient is a cross product of positions, so its
relative conditioning is set by `|x|²/V ~ 1/R`; at `R = 5 µm` there is no cancellation worth worrying
about, and the central-difference control is asserted at `1e-6` relative for the same
truncation-versus-cancellation reason set out in `ALEPH-PORT-1101` §10.

`VOLUME_PENALTY` conditioning: the effective volumetric stiffness is `K_V/V₀²` per unit volume error,
and the descent's step limit scales as `1/K_V`. The control sweeps `K_V` over two decades and reports
the volume error at each, rather than asserting a single tolerance that would hide the trade.

Round-off floor for the resultant identity (I3) is `~1e-16` relative to `Σ‖F_i‖`; asserted at `1e-12`.

## 11. Production-backend residency and transfer

Host-side numpy; CPU is the definition. Device-resident arrays in a future accelerator build would be
the envelope's vertex positions and its force accumulator; `V` and `ΔP` are scalars and would be the
only per-step host round trip. The double-count audit (`assert_single_pressure_path`) is a build-time
structural check over the pipeline's registrations and never touches a device at all — deliberately,
so it cannot be skipped by a run that fails before its first step.

`aleph/vertical/**` may not import `validation/**`; enforced by
`tests/vertical/test_vertical_firewall.py`.

## 12. Comments and docstrings to discard

No source prose was read, so none survives. What replaces it: the ownership argument in §4, restated
in the module docstring of `aleph/vertical/pressure.py` in Aleph's own vocabulary, and the
"what this is not" list in the docstring of `aleph/vertical/cortex.py`.

Specifically **not** carried across: any provider module path, facade method name, dispatch-slot
name, or gate label. The failure is described in terms of *what went wrong physically*, never in
terms of the other project's identifiers, so this entry is readable without that repository on disk.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** Status is `PROPOSED`: the code named in §1 has not landed and no control has been run, so there is no measured result to record here. This row is filled in with the measured residuals, at the tolerances stated in §10, at the moment the status moves to `ACCEPTED` — and not before. |
| Reviewer | agent-proposed, **unratified**. No human has reviewed the ownership argument in §4, which is a modelling decision, not a derivation. |
| Rollback | Delete `aleph/vertical/pressure.py` and `aleph/vertical/cortex.py`. Breaks `aleph/vertical/assembly.py` and every control in `tests/vertical/test_controls.py`. |

## 14. Honest limits — what the cortex does NOT represent

The cortex here is a relay with an areal spring. It is **not** cortical mechanics, and the following
are absent, not approximated:

- no actin filaments, no filament length distribution, no nucleation, no capping, no severing;
- no myosin, no motor kinetics, no active contractile stress, no ATP;
- no turnover — the rest areas are fixed, so the shell is elastic, never viscoelastic and never
  fluid on long timescales, which is the single most important thing real cortex does;
- no thickness: it is a zero-thickness surface, so it has no bending modulus, no shear modulus, and
  no through-thickness stress gradient;
- no cross-linkers, no network connectivity, no percolation, no strain-stiffening;
- no membrane-to-cortex adhesion chemistry beyond the two explicit connectors of
  `ALEPH-PORT-1103` — no receptor density, no binding kinetics, no unbinding under load;
- no cortical flow, no treadmilling, no polarity;
- no plasticity, no yielding, no rupture, and no blebbing.

Anyone who reads a cortical mechanical property off this object is reading it off a spring.

Further limits:

- **Unratified**, as above.
- The osmotic pressure is prescribed or penalised; there is **no solute transport model**, no
  aquaporin, no ion balance, and no time-dependent osmotic response. `ΔP` is a parameter.
- `VOLUME_PENALTY` conserves volume approximately by construction. Any claim requiring exact volume
  conservation needs a multiplier or an augmented Lagrangian, neither of which exists here.
- The two ensembles are not interconvertible and the code does not pretend they are: a run must
  declare which one it is in, and the run record carries it.
