# ALEPH-PORT-1103 — Tensile-only tether, compressive-only contact, and the assembled load path

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1103` |
| Lane | `L11 vertical — membrane / cortex / pressure` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `DEFECT_STUDY` — one inherited defect description; no inherited code. |

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
from aleph.vertical.connectors import (
    ConnectorGeometryError,
    CortexMembraneContact,
    ErmTether,
    UnilateralSpring,
    site_pair_forces,
)
from aleph.vertical.assembly import (
    LoadPathError,
    VerticalWorld,
    build_vertical,
    transmitted_load_on_membrane,
)
from aleph.vertical.relax import (
    MonotonePotentialDescent,
    RelaxationReport,
    relax_to_equilibrium,
)
```

Aleph target files: `aleph/vertical/connectors.py`, `aleph/vertical/assembly.py`,
`aleph/vertical/relax.py`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | none — **no source commit was read for this entry.** No file was opened. |
| Source path | none read. |
| Source symbol(s) | none read. Two *behaviours* were supplied as prose by a prior read-only audit recorded in PLAN §1.1. |
| Read from | neither. |
| Working tree == commit? | not applicable — nothing was read from either. |

## 3. Why source-derived porting beats clean-room

It does not, for the code. It does for one **argument**, which is the reason this lane exists at all:

> A tensile-only molecular tether was implemented and a compressive contact was declared but never
> evaluated, while a structural coverage gate passed for it — because the gate proved that the
> connector had been *assigned* to something callable, not that anything ever called it. Meanwhile
> the incumbent driver sidestepped both by pre-stretching the tether's rest lengths so it held the
> turgor at `t = 0`.

That is a case analysis which is easy to get subtly wrong and which was discovered by somebody
else's run. `ports/TEMPLATE.md` §3 names exactly that as a legitimate reason to let something cross.
What crosses is the case analysis. No code, no identifiers, no constants.

The Aleph consequence, stated as physics rather than as a bug report:

- A molecular linker is a **tether**, not a strut. Slack, compressed, or broken ⇒ **exactly zero**
  force, not a small one and not a negative tension.
- Therefore a **separate, non-adhesive compressive** element is required to carry an outward load
  across the cortex/membrane gap. Without it the turgor has no path to the membrane.
- A pre-stretched rest length that makes the tether hold the load at `t = 0` is a third thing again:
  it makes the *initial condition* carry the physics, so the model appears to work and its declared
  load path is fiction. This vertical's default configuration deliberately leaves the tether slack,
  so that if the contact is removed the transmitted load is exactly zero and the failure is loud.

## 4. Physical or mathematical law represented

**Unilateral spring.** Both connectors are one object with one gate. For a scalar separation `s`
along the line of centres and a rest separation `s₀`:

```
tether  (tensile-only)      U(s) = (k/2) max(0, s − s₀)²
contact (compressive-only)  U(s) = (k/2) max(0, s₀ − s)²
```

`U` is `C¹` at `s = s₀` (both value and derivative vanish there) and `C⁰` but not `C²` — the
stiffness jumps. That is the correct regularity for a unilateral constraint and it is why the sign
test is a test of an *exact zero*, not of a small number: `dU/ds` is identically zero on the whole
inactive branch, not merely small near the gate.

**Central force, hence closing moments.** Both act along `u = (x_a − x_b)/‖x_a − x_b‖`. A pair of
equal and opposite forces along the line of centres has zero net moment about any point. Choosing a
non-central direction — the cortex surface normal, say, which is tempting because it is the physical
direction of a contact traction — would leave a net couple on the two-body system with nothing to
absorb it. That is why the line of centres is used and the resulting `central_force=True` claim is
verified numerically rather than declared.

**Step-average force by the discrete gradient.** The force magnitude reported over a step from `s₀`
to `s₁` is

```
F̄ = −(U(s₁) − U(s₀)) / (s₁ − s₀)        (and −U′(s) when s₁ = s₀)
```

This is the discrete-gradient (Gonzalez) rule. Three properties, all of which matter here:

1. It closes the energy books **exactly**, by construction: the endpoint work equals `−ΔU` to
   round-off, for any step, including one that crosses the engagement gate mid-step. A midpoint rule
   does not, and would leave an energy residual at exactly the moment the physics is most delicate.
2. On a quadratic branch it reduces to the ordinary midpoint force `−k(s̄ − s₀)`, so nothing exotic
   happens in the common case.
3. It preserves the exact zero: if both ends of the step are on the inactive branch, `U` is zero at
   both, so `F̄ = 0` exactly. The unilateral property survives the time integration, which a
   naive average of the two endpoint forces would also do — but a *midpoint-geometry* rule would
   not, since the midpoint of a slack and a taut configuration can be taut.

**Distributed coupling, one adjoint pair.** The connectors act at every one of the `V` site pairs,
but `Connector.accumulate` takes exactly two endpoint handles and returns exactly one `AdjointPair`.
The resolution is not to weaken the protocol: each handle stands for a whole body's site set, and
the pair reports the **resultant** of the whole distributed coupling —

```
force_a = Σ_i f_i^(a)      force_b = Σ_i f_i^(b)      point_a = point_b = c
moment_a = Σ_i (x_i^(a) − c) × f_i^(a)      moment_b = Σ_i (x_i^(b) − c) × f_i^(b)
```

with `c` the common centroid of all sites. The moment residual then reduces to `moment_a + moment_b`,
which is the true net moment of the distributed coupling, assembled from two independent per-site
sums. It is a real check, not an identity: it vanishes only because every per-site force is central.

The two sides are accumulated **independently** — each from its own directed unit vector — so force
closure remains a check. `force_b := −force_a` is never written anywhere in this module, and the
control asserts that the two sums cancel rather than assuming it.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| separation `s` | µm | m | > 0 |
| rest separation `s₀` | µm | m | > 0 |
| stiffness `k` | pN/µm | N/m | ≥ 0 |
| force | pN | N | finite |
| stored energy | pN·µm | J | ≥ 0 |
| mobility | µm/(pN·s) | m/(N·s) | > 0 |

Singular cases:

- Coincident sites (`s = 0`) have no line of centres. Refused with `ConnectorGeometryError`, never
  regularised — a contact pair that has collapsed to a point is a geometric failure, and a softened
  denominator would turn it into a large finite force pointing in an arbitrary direction.
- A broken tether (`bound = False`) returns exactly zero regardless of `s`, and stays broken until
  something rebinds it. Rupture is state, not a threshold re-evaluated each step.
- `k = 0` is legal and returns exactly zero; it is how a connector is switched off without being
  removed from the schedule, which is the honest way to isolate a load path in a control.

Invariants:

- **I1. Compressed tether carries exactly zero.** Not `< tol`. Exactly `0.0`.
  `tests/vertical/test_connectors.py::test_compressed_tether_carries_exactly_zero`.
- **I2. Stretched contact carries exactly zero.** The non-adhesive half of the same statement.
  `tests/vertical/test_connectors_solid.py::test_both_contacts_push_when_close_and_carry_exactly_zero_when_separated`.
- **I3. Force closure.** `‖F_a + F_b‖` at round-off relative to `‖F_a‖ + ‖F_b‖`, from two
  independent sums. `tests/vertical/test_vertical_controls.py (specific control never recorded)`.
- **I4. Moment closure.** `‖M_a + M_b‖` at round-off relative to the moment scale, same test.
- **I5. Energy closure.** `W + ΔU + D − A = 0` at round-off for a step with non-zero endpoint
  displacement, including a step that crosses the engagement gate.
  `tests/vertical/test_connectors.py::test_energy_books_close_across_the_engagement_gate`.
- **I6. Turgor reaches the membrane only through the contact.** With the contact excluded, the load
  transmitted to the membrane is exactly zero.
  `tests/vertical/test_vertical_controls.py (specific control never recorded)`.
- **I7. A registered-but-unbound connector cannot pass.** The pipeline's schedule check and the
  coverage gate both refuse it.
  `tests/vertical/test_vertical_controls.py (specific control never recorded)`.

## 6. Source evidence class and known retractions

The inherited defect description carries evidence class `AUDIT_READ` at the commit recorded in
PLAN §1.1 and in `vault/02_Graph/claims/claim-coverage-gate-passes-unimplemented.md`. Aleph restates
it as a hypothesis about a class of mistake, tests for that class in Aleph's own code, and makes no
claim about the reference project's present state. Looked for retractions in PLAN §1.1,
`docs/decisions/OVERNIGHT_DECISIONS_2026-07-30.md`, and the vault claim; none found. No numerical
value was inherited, so there is nothing that could be retracted numerically.

## 7. Independent oracle or derivation

- **The unilateral condition itself.** `dU/ds ≡ 0` on the inactive branch is an identity, so the
  control is an equality against exact zero and not a tolerance.
- **Newton's third law and the two-force-member theorem**, evaluated by
  `aleph/runtime/ledger.py` from two independently accumulated sums.
- **The discrete-gradient energy identity** `W = −ΔU`, exact by construction of `F̄`, and asserted
  including across the gate crossing — which is precisely where a midpoint rule fails.
- **Global dilational virial.** At equilibrium, for the membrane alone,
  `Σ_i F_i^(internal) · x_i + Σ_i F_i^(transmitted) · x_i = 0`. With `Σ_i F_i^(tension) · x_i =
  −2σA` and `Σ_i F_i^(bending) · x_i = 0` (exact scale invariance, `ALEPH-PORT-1101` I5), this gives
  `σ = W_transmitted / (2A)` with no continuum limit taken. This is the identity the Laplace control
  reads σ off, and it is derived here rather than fitted.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | Force and moment residuals of the distributed pair are at round-off relative to their own scales, from independent accumulations. |
| Positive | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | A taut tether pulls, with magnitude `k(s − s₀)`, along the line of centres. |
| Positive | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | A compressed contact pushes, with magnitude `k(s₀ − s)`, along the line of centres. |
| Positive | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | With the proper load path assembled, the load reaching the membrane equals the turgor the cortex carries, to the cortex's own elastic residual, which is reported. |
| Positive | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | The whole assembly's net force is at round-off with no external load. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_connectors.py::test_compressed_tether_carries_exactly_zero` | A compressed tether returns exactly `0.0`, not a negative tension. Fails the instant the tether becomes a strut. |
| Negative (must fail) | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | The deliberately-broken bidirectional variant *does* carry compression, and the audit detects it. If this stops detecting, I1 is unguarded. |
| Negative (must fail) | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | With the contact excluded — coverage still green, because the exclusion is declared — the transmitted load is exactly zero and the membrane loses volume. |
| Negative (must fail) | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | The exact structural failure this lane exists for: a connector that is registered, that has an owner, and that no slot ever evaluates. |
| Negative (must fail) | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | A connector bound to a slot that launches nothing fails the witness check. |

## 10. Numerical and precision envelope

float64 throughout. Force, moment and energy residuals are checked by `aleph/runtime/ledger.py`
against a mixed absolute/relative band; this lane uses `Tolerance(relative=1e-9, absolute=1e-12)`,
the ledger default, because every residual here is a cancellation of `O(V)` float64 terms and the
accumulated round-off is `O(√V·ε) ≈ 1e-15` at `V = 162`. Six decades of margin.

The discrete-gradient force switches to `−U′(s)` when `|s₁ − s₀| < 1e-13 µm`, which is where the
difference quotient's cancellation error would exceed the value it computes. The switch is exact at
the crossover in the sense that both branches agree to `O(Δs)` there, and the control exercises both.

The unilateral zeros are exact zeros, so they are asserted with `==`, not with a tolerance. That is
deliberate: a tolerance would let a small wrong force through, and a small wrong force in this exact
place is what the whole lane is about.

Relaxation: first-order descent with adaptive step. Stability limit is `dt < 2/(M·k_max)`; the
controller finds it by rejection rather than by assuming it, so a stiff contact slows the run instead
of exploding it. Convergence is reported as `max‖F_i‖` and is never assumed.

## 11. Production-backend residency and transfer

Host-side numpy; CPU is the definition. In a future accelerator build the device-resident arrays are
the two bodies' positions and force accumulators plus the per-site rest lengths and bound flags; the
per-step host transfer is the six scalars of the `AdjointPair` resultant and the energy scalars —
which is exactly the point of reporting a resultant rather than a per-site pair list.

`aleph/vertical/**` may not import `validation/**`. Enforced by
`tests/firewall/test_scope_firewall.py::test_runtime_never_imports_the_analytic_oracles`.

## 12. Comments and docstrings to discard

No source prose was read, so none survives. The replacement is §3 and §4 above, restated in the
module docstrings of `aleph/vertical/connectors.py` and `aleph/vertical/assembly.py` in terms of
physics only.

Deliberately **not** carried into `aleph/**`: the other project's connector names, facade method
names, dispatch-slot names, gate labels, module paths, or repository name. The failure is described
as "a declared load path that is never evaluated", which is a statement anyone can check against
Aleph's own code without that repository on disk. `tests/ports/test_port_discipline.py` scans for
leaks of all of them.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** Status is `PROPOSED`: the code named in §1 has not landed and no control has been run, so there is no measured result to record here. This row is filled in with the measured residuals, at the tolerances stated in §10, at the moment the status moves to `ACCEPTED` — and not before. |
| Reviewer | agent-proposed, **unratified**. No human has reviewed the load-path argument in §3–§4, which is the substantive claim of this lane. |
| Rollback | Delete `aleph/vertical/connectors.py`, `assembly.py`, `relax.py` and their tests. Breaks every control in `tests/vertical/test_controls.py`. Nothing outside `aleph/vertical/` depends on them. |

## 14. Honest limits

- **Unratified.** The load-path argument is the whole point of the lane and no human has signed it.
- **The separation is the line-of-centres distance, not a signed normal offset.** This is what makes
  the pair central and its moments cancel, and it is also a real limitation: if a membrane site
  passes *through* the cortex the distance becomes positive again on the far side and the connectors
  respond as though it were outside. This vertical operates in the regime where the membrane stays
  outside the cortex; the code refuses only the coincident case, not the inverted one. A signed
  formulation would need a torque-carrying third body and is out of scope here.
- **The tether has no unbinding kinetics.** `bound` is a flag that nothing in this vertical ever
  flips. There is no force-dependent rupture, no rebinding rate, no receptor density, and therefore
  nothing here is evidence about tether lifetime or about blebbing.
- **The contact has no friction and no adhesion.** It is normal, unilateral, and elastic. There is no
  tangential traction at all, so the membrane can slide over the cortex at zero cost. Real
  membrane-cortex contact does not.
- **One site pair per vertex index.** The correspondence between membrane and cortex sites is fixed
  by construction (the two shells share a connectivity), not discovered by a proximity search. That
  is a strong simplification and it is the reason no contact-candidate search exists in this lane.
  A deforming cell would need one.
- **`relax.py` is scaffolding.** It is a first-order descent whose only job is to bring a surface
  close enough to equilibrium that the controls can be read. It is not the production solver, it has
  no inertia, no thermal forcing, no implicit treatment of the stiff contact, and its acceptance
  predicate (`MonotonePotentialDescent`) is **not** a proposal for `ALEPH-DQ-104` — that decision is
  reserved for the PI and this predicate would be a bad answer to it, because monotone energy descent
  is a property of a minimiser, not of a physical trajectory.
- **No dynamics, no timescales, no thermal statistics.** Everything here is quasi-static. Nothing in
  this entry is evidence about relaxation rates, viscosity, or the fluctuation spectrum.
