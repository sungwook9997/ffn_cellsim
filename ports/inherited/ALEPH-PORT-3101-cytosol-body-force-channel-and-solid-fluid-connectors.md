# ALEPH-PORT-3101 — a body-force channel in the Darcy flux law, and the five solid-to-fluid connectors

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3101` |
| Lane | lead session `f70de564` — connector layer |
| Status | `PROPOSED` |
| Written | `2026-07-30` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Decision this implements | `docs/decisions/PROPOSAL-solid-to-fluid-is-two-books-not-one.md` |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
# added to aleph.vertical.cytosol
CytosolField.body_force_pn_per_um3        # read-only (dof_count, dim) view
CytosolField.clear_body_force             # zero it, as clear_reaction_source does
CytosolField.body_force_resultant_pn      # the integral, for measuring force_b
CytosolField.internal_face_body_force     # the face projection the flux law consumes
TransferStencil.scatter_force_pn          # (dim,) [pN] through the SAME weight vector

# new module aleph.vertical.connectors_fluid
FluidCouplingError
ImmersedDragConnector
build_sf_cytosol_transfer
build_nmii_cytosol_transfer
```

**Two connectors, not five.** An earlier draft of this section named a `MovingBoundaryConnector` and
three further builders. They do not exist, and the entry is corrected rather than the names left
standing: `tests/ports/test_port_discipline.py`'s ORPHAN-C signal exists precisely to catch a ledger
naming an Aleph target that was never written, and it is the one signal in that file that cannot be
gamed. `ALEPH-PORT-2801` §12 records the same class of defect in `ecm.py`'s `__all__` — five exported
names that did not exist — and the rule it drew is the one applied here: **an absent name is a
visible gap and a stub is an invisible one.** What blocks the other three is in §14.

`aleph.vertical.wiring.CONNECTOR_PROTOCOL` changes from `("evaluate_sites",)` to
`("evaluate_sites", "accumulate")`. **This is a strengthening**; §14 records the check that it costs
nothing today.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | **none named and none read** |
| Source symbol(s) | **none** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened for this work. |
| Working tree == commit? | not applicable — nothing was read |

`PLAN.md` §7 records that the reference project's `SOLVE_COUPLED` has a declaration and zero
production tasks, so it has never converged a solid and a fluid together. There is nothing there to
port even if the policy allowed it.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** Darcy's law with a body force,

```
q = -(k/mu) (grad p - b)
```

is the textbook form — it is how gravity enters Darcy flow — and the discretisation is forced by the
operators `aleph/vertical/cytosol.py` already has. Per `ports/TEMPLATE.md` §3 that is a reason to
write clean-room.

## 4. Physical or mathematical law represented

**The claim in one line: a solid-to-fluid connector keeps two conservation books, and only the
momentum book is a force pair.**

The mass book is the scalar volumetric source `TransferStencil.scatter` already writes into
`_reaction_source`, which enters the content balance at `cytosol.py:1408`. It exists even for a body
moving with zero relative velocity and therefore zero drag, so it is **not** the reaction to the
drag.

The momentum book is what is missing. Today `internal_face_flux_um_per_s` is

```
q_f = -mobility * (p_hi - p_lo) / dx                       # cytosol.py:1187
```

with no body-force term, so an immersed body has nowhere to deliver its drag reaction and `force_b`
could not be reported honestly. Adding `b` [pN/µm³]:

```
q_f = -mobility * ( (p_hi - p_lo)/dx  -  b_f )
b_f = 0.5 * (b_lo + b_hi) . n_f          # n_f the face normal, lo -> hi
```

Its consequence for the monolithic solve. The content balance is
`V dc/dt + V div q = V source`, and with `div q = mobility (Lap p) - mobility (div b)` the extra term
is known, so it moves to the right-hand side:

```
b_rhs[n_u:] += V * mobility * div_b
(div_b)_lo += b_f / dx ,   (div_b)_hi -= b_f / dx        # internal faces only
```

Internal faces only, matching the Laplacian, which is assembled from exactly those faces
(`cytosol.py:1039-1050`). Using a different face set for the two terms would break the telescoping
that makes the mass balance an identity rather than a tolerance.

**Why the scatter is exact.** `TransferStencil` holds one weight vector, renormalised to sum to one
over the fluid cells in the support (`_multilinear_weights`, `cytosol.py:1517`). Scattering a force
through those same weights therefore deposits a total equal to what was handed in — the same
transpose property that already makes the power the solid loses equal the power the fluid gains. The
new `scatter_force_pn` uses `self._weights`, not a second stencil; a body that takes drag through one
stencil and returns it through another is a momentum source whose error scales with its speed, and
that is stated on the class already.

**Why `force_b` is measured and not assigned.** `AdjointPair`'s docstring
(`aleph/runtime/participant.py:310`) forbids `force_b := -force_a`, because that makes the
force-closure check a tautology. So the connector reads the fluid's body-force array **before and
after** its scatter and reports the measured delta, integrated over cell volume. It equals `-force_a`
to round-off when everything is right, and it does **not** when the weights fail to normalise, when
the stencil is stale, or when the scatter is partial. That is the difference between a check and an
identity.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI | Domain |
|---|---|---|---|
| `b`, pore-fluid body force | pN/µm³ | 1e6 N/m³ | finite |
| `mobility_um4_per_pn_s` | µm⁴/(pN·s) | — | `> 0` |
| scattered force | pN | 1e-12 N | finite |
| drag coefficient | pN·s/µm | — | `>= 0` |

Singular and boundary cases:

- **Non-finite scattered force**: refused. A non-finite force is a physics finding, not a value to
  pass along.
- **Wrong-shape scattered force**: refused; `(dim,)` exactly.
- **A transfer point with no fluid cell in its support**: already refused by `_multilinear_weights`,
  and that refusal is what stops a body returning its drag to nobody.
- **An owner declaring both drag channels**: still refused by `declare_drag_channel`. A connector
  using this channel declares `RESOLVED_DARCY`.
- **A stencil issued before a re-classification**: `_install_classification` reallocates the sinks,
  so a stencil held across one points at a dead array. Pre-existing and **not fixed here**; recorded
  in §14.

Invariants, each with the test that will assert it:

- **I1.** A scattered force is recovered exactly: `V * sum_cells b == force handed in`.
- **I2.** The scatter is the transpose of `interpolate_vector` — `<v, scatter(f)> == <interpolate_vector(v), f>`.
- **I3.** A uniform body force with the pressure free carries a Darcy flux equal to `mobility * |b|`
  on every internal face, and the discrete solve reproduces it.
- **I4.** Zero body force reproduces the previous flux **bit-for-bit** — the term is additive and
  inert when unused.
- **I5.** `force_a + force_b` closes against the constituent scale, with `force_b` read from the
  fluid array.
- **I6.** The body force is a **per-step accumulator, not owned state**, and is cleared rather than
  rolled back — exactly as `_reaction_source` is. `CytosolField.snapshot` (`cytosol.py:1674`) holds
  nine members and `_reaction_source` is not among them, deliberately: a quantity connectors rebuild
  every step has no pre-candidate value to restore. Asserted as *absence* — the control checks that
  the snapshot tuple's length is unchanged, so a later author who adds the body force to it fails a
  test rather than silently making a rejected step replay last step's forces.

## 6. Source evidence class and known retractions

No claim is made about the reference implementation. `PLAN.md` §7's finding that `SOLVE_COUPLED` has
zero production tasks is quoted from Aleph's own audit record and not re-verified here.

The `cytosol` contract's `citation_status` is `UNSOURCED` and nothing here changes it.
`ALEPH-PORT-1803`'s rule that no porosity, permeability or viscosity value may appear in code is
untouched — every magnitude in the controls is chosen to make the algebra sharp.

## 7. Independent oracle or derivation

1. **The uniform-body-force closed form.** With `b` uniform and the pressure unconstrained, the
   steady Darcy discharge is exactly `q = mobility * b` on every internal face, independent of the
   grid. A rescaled body-force term cannot reproduce it.
2. **The discrete transpose identity**, checked as an inner-product equality rather than by
   inspecting the weights.
3. **Exact recovery of the scattered total**, which is an oracle for the renormalisation and is
   independent of every constant in the flux law.

## 8. Positive control

**Measured this session on this tree** (`/Users/sw1/miniconda3/envs/aleph/bin/python`, macOS/darwin,
float64). Every number below came from a run in this session; none is quoted from memory.

| Control | Location | Asserts | Measured |
|---|---|---|---|
| Positive | `tests/vertical/test_connectors_fluid.py::TestTheBodyForceChannelIsExact::test_the_scattered_force_is_recovered_exactly` | I1 | scattered `(3.0, −1.25, 0.75)` pN, recovered with max error **4.44e-16** |
| Positive | `tests/vertical/test_connectors_fluid.py::TestTheBodyForceChannelIsExact::test_the_scatter_is_the_transpose_of_the_read` | I2 | `<v, scatter(f)>` = −1.0110964458948895, `<interpolate_vector(v), f>` = −1.0110964458948892, difference **2.22e-16** |
| Positive | `tests/vertical/test_connectors_fluid.py::TestTheFluxLawAgainstItsClosedForm::test_a_uniform_body_force_carries_the_predicted_darcy_flux` | I3 | `b = 7.0` pN/µm³, mobility 2.5 → predicted **17.5** µm/s; measured min = max = **17.5** on every x-face, and **exactly 0.0** on every face normal to y or z |
| Positive | `tests/vertical/test_connectors_fluid.py::TestTheFluxLawAgainstItsClosedForm::test_zero_body_force_reproduces_the_previous_flux_bit_for_bit` | I4 | `np.array_equal` **True** on a random pressure field |
| Positive | `tests/vertical/test_connectors_fluid.py::TestTheFluxLawAgainstItsClosedForm::test_the_body_force_reaches_the_solve_and_changes_the_pressure` | the channel is wired into the monolithic solve, not only the flux accessor | pressure difference **> 1e-9** pN/µm² after one backward-Euler step |
| Positive | `tests/vertical/test_connectors_fluid.py::TestTheConnectorDeliversAndTheLedgerCanSeeIt::test_the_pair_closes_against_the_constituent_scale` | I5, judged by the real `ForceWorkLedger` | verdict `ok` |

**Suite for this entry: 23 passed** (`tests/vertical/test_connectors_fluid.py`, 0.17 s).
`tests/vertical/test_cytosol_controls.py` re-run after the additions: **47 passed** — the pre-existing
owner is undisturbed, which is what I4's bit-identity is there to guarantee.

The I3 result is the one worth reading twice. The cross-axis flux is **exactly `0.0`**, not small: a
body force along `x` must drive nothing through a face normal to `y`, and a face-projection bug that
picked up the wrong component would show up here and nowhere else — a total-flux check would average
it away.

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_connectors_fluid.py::TestTheNegativeControlsFail::test_an_unnormalised_scatter_loses_momentum_and_the_closure_sees_it` | With the weight renormalisation defeated, the measured `force_b` no longer cancels `force_a` and the ledger's closure check fails |
| Negative (must fail) | `tests/vertical/test_connectors_fluid.py::TestTheNegativeControlsFail::test_a_connector_that_scatters_nothing_is_caught_by_the_closure` | A connector that computes a drag and delivers none reports `force_b == 0` and fails closure, rather than reporting a comfortable pair |

The second control is the one that matters. `PLAN.md` §6.1 records a negative control that could not
fail because a cold cache returned `0.0` for both the correct element and the broken one. Here the
failing quantity is read from the **fluid's** array, so "nothing was delivered" and "delivery
balanced" are different numbers.

## 9a. A mutant no measurement in this entry can kill — added 2026-07-31

**Reported by the `connectors_transfer` lane (`ALEPH-PORT-3304`), which found it in its own module
and stated that it applies unchanged here.** Verified against this module before recording.

The forbidden shortcut is `force_b := -force_a`, which `AdjointPair`'s docstring names as "the single
most effective way to make a broken connector look correct". §4 of this entry argues at length that
`force_b` is *measured* rather than assigned, and that argument is sound for a **lost or partial
delivery** — the negative controls in §9 both fail as intended.

It does **not** hold against the source-level shortcut, and the reason is a proof rather than a gap
in the controls:

```
force_a = -zeta * v        force_b = +zeta * v
```

For an exactly linear drag law `-(zeta*v)` and `zeta*(-v)` are **bitwise identical** in IEEE-754.
Measured here over 2000 random draws of `zeta` and `v`: **zero differences**. Reading `force_b` back
from the fluid's own array does not rescue it either, because the values scattered are the same
values. So no number this module can produce distinguishes the correct implementation from the
forbidden one.

`tests/vertical/test_connectors_fluid.py::TestTheNegativeControlsFail::test_the_two_sides_are_computed_independently_in_the_source`
therefore reads the source, and **its own docstring labels it as textual rather than physical** so a
later reader does not bank it as a measurement. It is weaker than a numerical control and is the only
thing standing between this connector and a shortcut nothing else here can see.

**The general form, worth carrying beyond this entry:** a Newton-pair closure check is only a real
check when the two sides come from computations that can *disagree*. For a linear law evaluated from
one state vector they cannot, and `site_pair_forces` — the house pattern every distributed connector
here reuses — has the same property. The `ALEPH-PORT-3303` lane recorded the same thing independently
from the other direction, measuring a closure residual of exactly `0.0` and filing it as a limitation
rather than a result. Two lanes, two modules, one conclusion.

## 10. Numerical and precision envelope

float64 throughout; no reduced-precision path. The body-force term is additive in the flux law, so
with `b = 0` the computed flux is bit-identical to the previous expression — I4 asserts `==`, not
`allclose`, because an additive zero is exact in IEEE-754 and anything looser would hide a term that
is quietly always active.

Closure (I5) is asserted against the **constituent** scale, never the resultant. `PLAN.md` §2.5
records a resultant-based tolerance rejecting 40,000 consecutive steps of a correct relaxation, and
the immersed case has the same shape: a body in the middle of the fluid delivers a large force at
every stencil cell whose sum is the thing that has to cancel.

## 11. Production-backend residency and transfer

Host, numpy, float64. **No GPU work of any kind was run by this lane, and no authorization was
sought or held.** The scatter is an `np.add.at` over at most `2^dim` cells and is a natural device
kernel; the monolithic solve is dense and is the part that will not scale, which is a pre-existing
property of this module and not something this entry changes.

## 12. Comments and docstrings to discard

Nothing to discard: no reference prose entered, because no reference file was read.

One piece of **Aleph's own** prose is corrected rather than discarded:
`docs/design/CONNECTOR_WIRING_MAP.md` §3a states that a solid-to-fluid coupling's adjoint pair "is a
vector force on the solid against a scalar source in the fluid, which is a correct Biot coupling and
is not two force blocks." The second half is right and the first half conflates the two books; the
proposal in `docs/decisions/` states the correction and §3a is updated to point at it.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** 23 controls pass (`tests/vertical/test_connectors_fluid.py`, 2026-07-30, 0.17 s) and the pre-existing owner is undisturbed (`test_cytosol_controls.py`, 47 passed). Status stays `PROPOSED`: the mutation testing this project requires of an `ACCEPTED` entry has **not** been run against this module, and two of the five connectors this entry names are not wired — see §14. |
| Reviewer | Agent-proposed. Unratified. No PI review. |
| Rollback | Revert the `cytosol.py` additions (three named members, all additive), delete `aleph/vertical/connectors_fluid.py` and its tests, and remove the five `Binding` rows. Restore `CONNECTOR_PROTOCOL` to `("evaluate_sites",)`. Nothing else imports any of it. |

## 14. Honest limits — what this entry does NOT establish

**Nothing here is evidence that a coupled solid–fluid system converges.** It gives momentum a place
to go. `SOLVE_COUPLED` — a whole cell converging as one — is what the reference project declared and
never ran, and this project has not run it either.

**The momentum reaches the skeleton through the pressure field, not directly.** In the quasi-static
Darcy limit the pore fluid has no momentum reservoir: a body force enters the flux law, changes the
pressure, and the pressure acts on the skeleton through the coupling blocks the module already has.
So the pair `force_a + force_b = 0` is a statement about what the *connector* exchanged, and the
subsequent transmission to the skeleton is the solve's business. **A control that the skeleton
ultimately receives the momentum is not written**, and would need a steady-state solve with a known
analytic answer.

**No magnitude here is a cytoplasm number.** Evidence rung `ANALYTIC_ORACLE`, quantitative status
`BLOCKED`. Drag coefficients in the controls are chosen to make the algebra sharp.

**Kinetics are absent, not implied.** None of the five connectors attaches, detaches, or turns over.

**The stale-stencil hole is pre-existing and stays open.** `_install_classification` reallocates
`_reaction_source` (and now the body-force array), so a `TransferStencil` held across a
re-classification writes into a dead array and loses its scatter silently. `TransferStencil` stamps
no topology epoch, so `endpoints.py`'s `EPOCH_NOT_STAMPED` applies and the staleness guard cannot
fire for it. Reported rather than fixed: it is a change to an owner's transaction discipline and
belongs with that owner's lane.

**Three of the five connectors this entry set out to unblock are NOT wired, and the count it
contributes is two.**

| Connector | State | Why |
|---|---|---|
| `sf_cytosol_transfer` | **wired** | `ImmersedDragConnector`, 23 controls |
| `nmii_cytosol_transfer` | **wired** | same class, different parameters |
| `membrane_cytosol_boundary` | not wired | needs a moving-boundary connector over `BoundaryFaceStencil`, which is a **different geometry** — a prescribed normal flux and a per-face traction, not an immersed point. The framing error in §3a blocked it, but removing that error does not write it. |
| `nucleus_cytosol_boundary` | not wired | same, on the interior boundary of the same domain |
| `surface_porous_transfer` | not wired | §4 below |

The two boundary connectors are the smaller remaining job of the two kinds, and `BoundaryFaceStencil`
already carries both halves (`prescribe_normal_flux`, `pressure_traction_pn`, and a per-face
`outward_normals`). Their sign convention is the one thing the registry flags as easy to get
backwards: an inward-facing normal would make the nucleus inflate under compression.

**`surface_porous_transfer` is wired against the cortex's porous-quadrature endpoint, which
`cortex_filaments.py:296` deliberately withholds.** Whether a filament graph presents a surface to
the fluid at all is the open cortex question in
`PROPOSAL-cortex-is-a-filament-graph-not-the-envelope.md`. If that connector cannot be given a real
endpoint it is **not wired**, and it stays out of `WIRED` rather than being counted — see §14 of that
proposal and `CONNECTOR_WIRING_MAP.md` §4.
