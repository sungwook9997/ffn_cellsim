# ALEPH-PORT-2002 — cytosol immersed-transfer connector contracts (group D, 6 continuous edges)

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-2002` |
| Lane | `L20 cross-system and cytosol connector contracts` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |

---

## 1. Aleph API

```python
from aleph.state.connectors_crosssystem_cytosol import (
    CYTOSOL_IMMERSED_TRANSFER,      # the six group-D contracts, in manifest order
    LANE_CONNECTORS,                # group C then group D, fourteen in all
    PRIOR_ART_NOTES,
    PRIOR_ART_REALIZATION,
    PriorArtRealization,
    SCHEMA_MISMATCH,
    SCHEMA_SOURCE,
    ConnectorContract,
    ConnectorFamily,
    build_contract,
    validate_contract,
)
```

`CYTOSKELETAL_CROSS_SYSTEM`, `LINC_CONNECTORS`, `connector`, `connectors_touching` and
`ConnectorContractError` are authorised by `ALEPH-PORT-2001`. Nothing outside these two lists is
covered.

The six contracts: `sf_cytosol_transfer`, `mt_cytosol_transfer`, `if_cytosol_transfer`,
`lamellipodium_cytosol_transfer`, `filopodium_cytosol_transfer`, `nmii_cytosol_transfer`. All six are
`family=C`, `commit_on_accept=False`, `bidirectional=True`, `adjoint_required=True`,
`internal_to=None`, `implemented=False`, and all six name the adjoint reaction in `mechanism`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/runtime.py`, `ffn_sim/ac/engine/contracts.py`, `ffn_sim/ac/engine/dispatch.py`, `ffn_sim/ac/engine/stress_fiber.py`, `ffn_sim/ac/engine/microtubule_rig.py`, `ffn_sim/ac/engine/intermediate_filament_rig.py`, `ffn_sim/ac/engine/protrusion.py`, `ffn_sim/ac/engine/nmii_actuator.py`, `ffn_sim/ac/engine/interior_column_slice.py`, `ffn_sim/ac/engine/composed_native.py`, `ffn_sim/ac/engine/monomer_flux.py` |
| Source symbol(s) | `ImmersedTransferConnector`, `CytosolFieldEndpoint`, `IntermediateFilamentFluidTransfer`, `ProtrusionGraphConnector`, `ProtrusionStepBindings`, `ProtrusionActors.accumulate_mechanics`, `_REQUIRED_CONNECTORS` (protrusion and NMII), `NMIIActuator.accumulate_candidate`, `canonical_facade_claims`, `accumulate_transfer` (all definitions), `_excluded_facade_noop` |
| Read from | **working tree** |
| Working tree == commit? | **yes for every file read except `microtubule_rig.py`**, which is one of the 31 uncommitted changes. `git diff be0e5876 -- ffn_sim/ac/engine/microtubule_rig.py` was inspected before the `mt_cytosol_transfer` verdict was written. Every other file cited above is clean at this commit. |

> The commit-versus-working-tree row exists because at `be0e5876` this tree carries 31 uncommitted
> changes. `HEAD` *is* `be0e5876`, so those are deltas on top of the cited commit rather than a
> different revision — but a reader reproducing this on a clean checkout will see a slightly different
> `microtubule_rig.py`, and only that file.

**File digests of every source file read for this entry** (`shasum -a 256` of the working-tree copy,
2026-07-30). The commit names a revision; the digest names the bytes actually read, which is the only
provenance that survives a dirty tree:

```
sha256:c0905caa0a2b9277c7fd902811fb4591c2737a8c8d760d3febb714afd8894fe1  ffn_sim/ac/engine/runtime.py
sha256:acb92dbd918b1f70a824b192e5f28940cb97ae74edbfa6b8a51ee17755257e72  ffn_sim/ac/engine/contracts.py
sha256:325307d38f7df0989c1bf71f421f5fb8ea9cf5341f3dd739e2ad890f70d2311a  ffn_sim/ac/engine/dispatch.py
sha256:209982c9def60d7602d815a402dc74dc2ec60879a4bd6826eda5603fbfd3a990  ffn_sim/ac/engine/stress_fiber.py
sha256:8728a9060ea7a5eadc7803c176cd300c81fde41dc55bcd1e6ea04a13186c0d0e  ffn_sim/ac/engine/microtubule_rig.py
sha256:8560fe7dd8b499821732df15a5a83d9327cc93b3e7813ebc7945e5d4a338362a  ffn_sim/ac/engine/intermediate_filament_rig.py
sha256:98300c76d61b3a6d572d91fa838b27fe127a50bb43eb19040fd8e81baba8fc63  ffn_sim/ac/engine/protrusion.py
sha256:1523ed040a259971a55cc34281b6b6c763f7b9a299fa0901e03ebf5d79a2ceae  ffn_sim/ac/engine/nmii_actuator.py
sha256:f9cacfc80f3bab393c6fd165d069b3ce438b266684f031b8cc2dd7e7cdc604d8  ffn_sim/ac/engine/interior_column_slice.py
sha256:d3115b37bdb2d3110841661f557435c5ee8b5a14ae71c5874c4d8b39aeb9fc63  ffn_sim/ac/engine/composed_native.py
sha256:3adebe166574634b73c5037ce9328d69e9ffd3c50422765188ae23be736723b7  ffn_sim/ac/engine/monomer_flux.py
```

`microtubule_rig.py` is the one file above that differs from `be0e5876`; its digest is the working-tree
digest, not the committed one. Every other digest is both.

## 3. Why source-derived porting beats clean-room

**It does not, and no code was ported.** The six declarations were re-derived from Aleph's own
manuscript extract `docs/manuscripts/extracted/APPENDIX_A_mechanics_registry.txt` (connector legend
plus group D). Endpoint owners, endpoint roles and mechanisms come from there; the three-part
decomposition of an immersed transfer, the momentum-source argument, the work-conjugacy argument, the
unbound-motor argument and the monomer-reservoir argument were written here.

The archive was read to answer the audit question in §6 — which of the six has an implementation — and
the answer is: none of them.

One thing was considered and rejected: the `ImmersedTransferConnector` protocol shape itself
(`accumulate_transfer(solid, cytosol)` plus transaction and ledger hooks). It is a reasonable
interface. It is also a two-argument method signature, which is not a liftable asset, and the audit
found a concrete reason Aleph should *not* copy the surrounding arrangement — see the protrusion defect
in §6, where the analogous protrusion interface omits the fluid endpoint and thereby makes two of these
six connectors structurally unimplementable through their own declared seam. Inheriting an interface
family that already failed that way once would be inheriting the mistake.

## 4. Physical or mathematical law represented

There *is* a law here, unlike group C, and stating it is what makes the `adjoint_required` field
mean something.

Let a solid quadrature point sit at `x_p` inside a fluid whose velocity field is `u`. Interpolation to
the point is a linear operator `S` acting on the fluid degrees of freedom:

    u_p = S u

Drag on the solid is a function of the slip:

    f_p = -zeta (v_p - u_p)

with `zeta` the local drag tensor (anisotropic for a filament: different along and across the axis).
The reaction delivered to the fluid must be the **transpose** of the same operator applied to the same
force:

    g = -S^T f_p

Two things follow, and they are the reason the transpose is not one scatter choice among several.

**Momentum.** With `S` a partition of unity in the sense `1^T S = 1^T` (every unit of interpolation
weight is accounted for), the total force on the pair is

    sum(f_p) + 1^T g = sum(f_p) - 1^T S^T f_p = sum(f_p) - sum(f_p) = 0

identically, at any configuration. Use any other scatter and the residual is nonzero, grows with `|v_p
- u_p|`, and therefore *looks like a velocity-dependent constitutive law* rather than an error. That is
the specific reason a momentum leak here is dangerous: it is plausible.

**Power.** The power the solid loses is `v_p . f_p`; the power the fluid gains is `u^T g = u^T (-S^T
f_p) = -(S u)^T f_p = -u_p . f_p`. Their sum is `(v_p - u_p) . f_p = -zeta |v_p - u_p|^2 <= 0` for
positive-semidefinite `zeta`. The pair dissipates and never generates. Any scatter other than `S^T`
breaks this and can make the pair a power source at some configuration, which no energy audit averaged
over a step will reliably catch.

So `adjoint_required=True` on all six is not a style preference: it is the condition under which the
force ledger and the work ledger are both meaningful. This is why the field is refused as `False`
rather than warned about.

**The second law on the two protrusion edges.** `lamellipodium_cytosol_transfer` and
`filopodium_cytosol_transfer` additionally carry monomer flux. If a barbed end elongates at rate `r`
consuming monomer at concentration `c(x_b)`, the sink on the field is `-r` delivered at `x_b` through
the same interpolation transpose, and the total monomer in the closed cell is conserved up to
polymerised mass:

    d/dt ( integral c dV + N_monomers_in_filaments ) = 0

Declaring drag alone leaves `c` unread and unmodified, which is formally the limit `c -> infinity`: an
infinite reservoir. Growth then has no monomer cost and no transport limit, and protrusion velocity
becomes a free parameter instead of a predicted quantity.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| (none) | — | — | This entry declares no numeric quantity. No drag coefficient, viscosity, permeability, polymerisation rate or monomer concentration is asserted. |

The law in §4 is stated in symbols precisely because no constant crosses this boundary. When these
connectors are implemented, `zeta` will need a sourced value with a literature citation and its own
ledger entry; declaring one here would be an unsourced constant with nothing to check it.

Singular and boundary cases, each with the behaviour Aleph requires:

- **A quadrature point outside the fluid domain.** Not a case these declarations can resolve; it is a
  requirement on any implementation, recorded as such — the interpolation stencil must refuse rather
  than clamp, because a clamped stencil silently violates `1^T S = 1^T` and breaks the momentum
  identity of §4.
- **Zero slip.** `f_p = 0` and `g = 0` exactly, not approximately: the expression is linear in the
  slip, so the zero is exact on the whole branch. That is testable with `==` when implemented.
- **`commit_on_accept=True` on a `C` connector.** Refused. Continuous transfer owns no discrete state,
  so the claim would give rollback something to restore that never existed.
- **`adjoint_required=False`.** Refused. §4 is why.
- **`bidirectional=False`.** Refused at construction and again by the lane gate.
- **A group-D connector whose `mechanism` omits the adjoint reaction.** Caught by test, not by the type
  — see I3 below. The declaration would still construct, which is exactly why the test exists.

Invariants that must hold, each with the test that asserts it:

- **I1.** Exactly six contracts, in manifest order —
  `test_group_d_is_the_six_named_connectors_in_order`.
- **I2.** Every group-D connector is `family=C` and sets `commit_on_accept=False` —
  `test_every_group_d_connector_is_continuous`, `test_no_group_d_connector_claims_a_commit`.
- **I3.** Every group-D `mechanism` names the adjoint reaction, states it is mandatory and not
  optional, and names the momentum-source consequence —
  `test_every_group_d_mechanism_names_the_adjoint_reaction`,
  `test_every_group_d_mechanism_states_the_reaction_is_not_optional`.
- **I4.** Every group-D `mechanism` names drag and velocity interpolation, i.e. all three parts of an
  immersed transfer — `test_every_group_d_mechanism_names_drag_and_velocity_interpolation`.
- **I5.** Every fluid-side endpoint is `cytosol` and every fluid-side role names a stencil —
  `test_every_group_d_fluid_endpoint_is_the_cytosol`,
  `test_every_group_d_fluid_role_is_a_transfer_stencil`.
- **I6.** The two protrusion edges carry the barbed-end sink and the G-actin field in both the
  mechanism and both endpoint roles —
  `test_the_protrusion_pair_carries_the_barbed_end_sink_as_well`.
- **I7.** The other four declare **no** monomer channel, so the dual role is not spread by copy-paste
  to connectors that do not have it —
  `test_the_four_non_protrusion_transfers_claim_no_monomer_channel`.
- **I8.** `nmii_cytosol_transfer` states that unbound-head minifilaments are included, and states the
  free-rigid-body and duty-ratio consequences —
  `test_nmii_transfer_covers_minifilaments_with_unbound_heads`.
- **I9.** Nothing claims to be implemented —
  `test_nothing_in_this_lane_claims_to_be_implemented`.

## 6. Source evidence class and known retractions

**None of the six has a concrete implementation anywhere in the archived tree.** Method:
`grep -rn "def accumulate_transfer"` over the whole repository yields seven definitions — one Protocol
declaration in `runtime.py`, one narrowed Protocol re-declaration in `intermediate_filament_rig.py`,
**one concrete implementation** in `interior_column_slice.py`, and four test doubles. The single
concrete implementation belongs to `surface_porous_transfer` (`cortex <-> cytosol`), a group-A edge,
and its constructor refuses any other name outright. So the exactly-one concrete immersed transfer in
the tree is not one of these six.

| connector | verdict | what it rests on |
|---|---|---|
| `sf_cytosol_transfer` | **SEAMED** | Typed slot `StressFiberStepBindings.cytosol_transfer: ImmersedTransferConnector`, alongside a declared `CytosolFieldEndpoint`. Satisfied only by a test double. |
| `mt_cytosol_transfer` | **SEAMED** | Typed slot `MicrotubuleStepBindings.cytosol_transfer`, validated at construction for name/endpoints, plus a declared `cytosol` field endpoint. Satisfied only by a test double. |
| `if_cytosol_transfer` | **SEAMED** | The best-shaped of the six: the IF rig narrows the protocol to a filament-specific `IntermediateFilamentFluidTransfer` sub-interface with its own field endpoint. Still satisfied only by a test double. |
| `lamellipodium_cytosol_transfer` | **SEAMED, AND THE SEAM CANNOT CARRY THE PHYSICS** | See the defect below. |
| `filopodium_cytosol_transfer` | **SEAMED, AND THE SEAM CANNOT CARRY THE PHYSICS** | Same defect. |
| `nmii_cytosol_transfer` | **NAME ONLY, AND IT PASSES A COVERAGE GATE** | See the second defect below. |

**Summary: 0 of 6 implemented. 4 typed seams, 2 of which are structurally inadequate, and 1 name.**
This is the largest unimplemented block among the fourteen this lane declares, and it is worse than
"unimplemented" in two specific ways.

**Defect 1 — the protrusion immersed-transfer seam has no fluid endpoint.**
`ProtrusionStepBindings` *requires* a runtime for `lamellipodium_cytosol_transfer` and
`filopodium_cytosol_transfer`: `_REQUIRED_CONNECTORS` lists both, and construction raises on a missing
or extra name. But the slot is typed as the facade's generic `ProtrusionGraphConnector`, whose
accumulation method is `accumulate_actor(actor: ProtrusionActorView) -> None` — it receives the
protrusion's own geometry view and **nothing else**. There is no `CytosolFieldEndpoint` argument, no
fluid view, and no field parameter of any kind, in that method or in `ProtrusionActors.accumulate_mechanics`
which calls it. So no object satisfying the required interface can read the porous-fluid velocity or
scatter a reaction into it. The declared contract's endpoint roles say "active filament quadrature and
barbed-end sink" against "porous drag and G-actin transport stencil", and neither half is reachable
through the seam that is supposed to carry it. Nor does a G-actin channel exist elsewhere: the only
monomer-flux machinery in the tree is a separate mechanism keyed on a connector name
(`"cortex_gactin_flux"`) that is not among the 36 registered connectors at all.

This is a stronger failure than an unfilled seam. An unfilled seam is work not done; this is work that
cannot be done without changing the interface, while the bindings assert it is required and the
architecture validator confirms the contract is bidirectional with adjoint transfer. Every structural
check passes and the physics has nowhere to go.

**Defect 2 — `nmii_cytosol_transfer` satisfies a dispatch coverage gate while being a name.**
`canonical_facade_claims()` assigns the edge to `NMIIActuator.accumulate_candidate`, with an in-code
comment explaining that it was added because an all-heads-unbound minifilament would otherwise be a
free rigid body — the physics is understood and written down. But `accumulate_candidate` iterates a
fixed dictionary of four motor edges (`nmii_sf_motor`, `nmii_cortex_motor`,
`nmii_lamellipodium_motor`, `nmii_filopodium_motor`) and dispatches those; `nmii_cytosol_transfer` is
not in it. The module contains **no occurrence of the string `cytosol`** anywhere. So the coverage
claim is satisfied structurally — the named method exists and is called exactly once — and the edge is
never evaluated by anything.

This is a second instance of the defect an earlier audit recorded for `membrane_cortex_contact`, which
was assigned to a surface facade method that had no slot for it. The gate is blind to both for the
same reason: it checks that a claimed edge is assigned to a callable that runs exactly once, not that
the callable does anything with that edge. **A coverage gate over names cannot distinguish a
connector from a string.** Aleph's runtime lane already answers this — coverage there requires
scheduled, reached, *and* evaluated-something with a witness count — and these two findings are the
independent confirmation that the third condition is the load-bearing one.

**Third finding, weaker but worth recording — the facades themselves are never built.** Five of the
six group-D edges are claimed by `StressFiberActor`, `MicrotubuleRig`, `IntermediateFilamentRig` and
`ProtrusionActors`. None of those classes is constructed anywhere outside `ffn_sim/tests/`. In the one
composed native path, `compose_native_cell_world` substitutes `_excluded_facade_noop` for the
microtubule, IF, protrusion, surface and fluid facade types, and the docstring is honest about it. So
even if the seams were filled, nothing in production would call them.

**Retractions and self-labelled gaps (looked in `STATE.md`, `STATE_NONQUOTABLE.md`,
`outputs/tag_kb/run_audit_report.md`, `outputs/ac/cell_assembled/composed_world_v2.json`, and the
module docstrings):** `composed_world_v2.json` self-labels `evidence: CENSUS-WIRED`,
`quantitative_claim_status: BLOCKED`. The interior-column connector's docstring is candid that it
implements only porous drag and pressure load and no separate no-flux kernel, and that its
snapshot/rollback hooks are structural no-ops. No retraction was found for any of the six here —
because there is nothing to retract: no claim was ever made about them.

**Reachability and test quality:** all six are unreachable. The tests that touch them are
architecture-validation and binding-shape tests exercising the plumbing, plus four test doubles whose
purpose is to be counted. Not one exercises a drag law or an adjoint scatter.

## 7. Independent oracle or derivation

* **Aleph's manuscript extract** fixes the six names, their endpoint owners, their endpoint roles and
  their mechanisms, and it is Aleph material. Every one is pinned by test, so a declaration drifting
  from the manuscript fails.
* **The derivation in §4** is the independent law, written from the interpolation operator rather than
  from any source. It is the oracle any future implementation must be checked against, and it supplies
  two exact checks — `sum(f_solid) + sum(f_fluid) = 0` identically, and non-positive pair power — which
  hold at *every* configuration and therefore do not require a converged state to test. That property
  matters here: an identity that holds only at equilibrium cannot separate a wiring error from
  incomplete convergence, which is exactly the trap this project has already met once, where a max-norm
  residual gate could not see a signed cancellation.
* **Internal consistency** is machine-checked: family determines commit semantics, the fluid endpoint
  is `cytosol` for all six, and the dual-role assertions are two-sided — the protrusion pair must
  declare the monomer channel and the other four must not, so a copy-paste in either direction fails.

Agreement with the archive is **not** the check and could not be: there is nothing there to agree with.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_connectors_crosssystem_cytosol.py::TestAllFourteenAreConstructible::test_group_d_is_the_six_named_connectors_in_order` | The six names appear exactly once each, in manifest order. |
| Positive | `…::TestGroupDNamesTheAdjointReaction::test_every_group_d_mechanism_names_the_adjoint_reaction` | All six name the adjoint reaction. |
| Positive | `…::TestGroupDNamesTheAdjointReaction::test_every_group_d_mechanism_states_the_reaction_is_not_optional` | All six state it is mandatory and name the momentum-source consequence. |
| Positive | `…::TestGroupDNamesTheAdjointReaction::test_the_protrusion_pair_carries_the_barbed_end_sink_as_well` | Both protrusion edges declare the barbed-end sink and the G-actin field in the mechanism *and* in both endpoint roles. |
| Positive | `…::TestGroupDNamesTheAdjointReaction::test_nmii_transfer_covers_minifilaments_with_unbound_heads` | The NMII edge includes unbound-head minifilaments and states both consequences. |
| Positive | `…::TestBidirectionalCannotBeFalse::test_all_fourteen_pass_the_lane_gate` | `validate_contract` accepts all fourteen — so a gate that rejected everything cannot pass as correct. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/state/test_connectors_crosssystem_cytosol.py::TestBidirectionalCannotBeFalse::test_the_lane_gate_rejects_adjoint_required_false` | A **continuous** contract with `adjoint_required=False` is refused. This is the group-D-specific negative control: drag is where dropping the reaction looks harmless. |
| Negative (must fail) | `…::TestValidationCanFail::test_a_continuous_connector_claiming_a_commit_is_rejected` | A `C` contract with `commit_on_accept=True` is refused. |
| Negative (must fail) | `…::TestBidirectionalCannotBeFalse::test_the_contract_type_itself_rejects_bidirectional_false` | `bidirectional=False` refused at construction. |
| Negative (must fail) | `…::TestBidirectionalCannotBeFalse::test_the_lane_gate_also_rejects_bidirectional_false` | Refused again by the lane gate, probed on a duck-typed stub so construction cannot mask the gate. |
| Negative (must fail) | `…::TestGroupDNamesTheAdjointReaction::test_the_four_non_protrusion_transfers_claim_no_monomer_channel` | The inverse control on I6: a monomer channel spread to a connector that has none fails. Without it, satisfying I6 by adding the words everywhere would pass. |

`test_a_valid_probe_declaration_constructs` establishes that the probe used for the rejections is well
formed, so each rejection differs from an accepted declaration by exactly one field.

## 10. Numerical and precision envelope

No arithmetic in this entry. Every assertion is an exact comparison of strings, booleans, enum members
and tuple order, so there is no working precision, no accumulation precision, no conditioning and no
tolerance. Tuple order is compared with `==`, so a connector inserted in the wrong position fails —
order is part of the manifest, not a presentation detail.

The derivation in §4 does carry a precision consequence, recorded now so it is not discovered later:
`sum(f_solid) + sum(f_fluid) = 0` is exact in real arithmetic and, in floating point, is a sum of
signed terms. It is therefore vulnerable to the same trap this project already documented for the
Laplace-law check — a max-norm residual cannot see a signed cancellation, and a small absolute value
can mean cancellation rather than correctness. Any future implementation of these six must gate on the
signed accounting identity with a stated tolerance relative to the *magnitude of the terms being
summed*, not on a max norm of the residual force.

Outside the declared vocabulary the behaviour is refusal, not degradation: unknown names raise
`KeyError`, out-of-vocabulary families raise at construction, `commit_on_accept=True` on a `C`
connector raises `ConnectorContractError`.

## 11. Production-backend residency and transfer

**This code never runs on the production backend.** Frozen dataclasses of strings, booleans and enum
members, host-resident, read at build time. No device allocation, no transfer, no host round-trip per
step, no per-step cost.

That is what keeps this module importable with no GPU and no `warp` present — the condition every test
in this lane runs under, and the condition gate R1 requires.

The residency question the *implementations* will face is recorded here because it is the reason these
six are hard and the group-C edges are not: an immersed transfer touches two owners' device arrays plus
a field, so the interpolation stencil, its transpose, the solid quadrature and the fluid degrees of
freedom must be co-resident or explicitly staged. A per-step host round-trip would dominate the cost.
This entry authorises none of that and makes no claim about it.

## 12. Comments and docstrings to discard

Nothing from the provider survives. Specifically discarded:

* **Provider module paths and symbol references** — every `ffn_sim.ac.engine.*` path and every
  cross-reference into that tree. They live in §2 and §6 of this entry and nowhere in `aleph/**`.
* **Provider status vocabulary as code** — `SEAMED`, `KERNEL_BOUND`, `PI-GAP`, `CENSUS-WIRED`,
  `Card-5 seam`, `T10`, gate labels. Replaced by Aleph's `PriorArtRealization` enum, whose members are
  defined by what was found rather than by an inherited label.
* **Dated authority claims** — "PI-ratified 2026-07-25", "PI D5-A, 2026-07-28". The provider's decision
  record, not Aleph's.
* **The `per_node_pressure_traction` / `aggregate_or_lumped` fidelity-marker attribute pattern.**
  Read and not carried: a boolean asserting its own fidelity is a claim with no check behind it, which
  is the shape this project exists to eliminate. Where Aleph needs that distinction it belongs in a
  test, not in an attribute.
* **`chemistry_card` strings** and the `"cortex_gactin_flux"` connector name, which is not in the
  registered 36 and would import an ontology Aleph has not declared.
* **Sourced constants.** No viscosity, drag coefficient, permeability or polymerisation rate crosses
  this boundary. In particular the cytoplasm-viscosity datum that appears in the archive is *not*
  inherited; any future use needs its own ledger entry and its own literature citation.

Replaced by: the derivation in §4 written from the interpolation operator; per-contract prose stating
why each connector's omission would be a specific physical error; and `PRIOR_ART_NOTES`, which states
the audit basis in Aleph's vocabulary and points here for the paths.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** `PROPOSED`. As of 2026-07-30, `tests/state/test_connectors_crosssystem_cytosol.py` passes 61/61 across both entries. That is evidence the declarations are self-consistent and that the validation can fail; it is not evidence that any transfer works, because none is implemented. |
| Reviewer | **Agent-proposed, unratified.** Written by lane L20 during autonomous operation. No PI review. |
| Rollback | Delete `aleph/state/connectors_crosssystem_cytosol.py` and `tests/state/test_connectors_crosssystem_cytosol.py`. Nothing breaks: nothing imports either as of this writing. The manifest loses six of thirty-six connector declarations, which blocks a state-schema completeness check and therefore the representation lanes' axis-order work — the dependency `ALEPH-PD-002` identified. |

## 14. Honest limits

What this entry does **not** establish:

* **No physics.** Not one of the six transfers a force or a monomer in Aleph. `implemented=False` is
  the status. A green suite means the declarations are consistent; consistency is not correctness.
* **The derivation in §4 is unexercised.** It is derived here and no test runs it, because there is
  nothing to run it against. Carried as `UNVERIFIED`. In particular the claim that `1^T S = 1^T` holds
  for whatever stencil Aleph eventually chooses is an assumption about a stencil that does not exist
  yet, and the momentum identity depends on it.
* **No drag law is chosen.** `zeta` is written as a tensor because filament drag is anisotropic, and
  which anisotropic form (slender-body, Stokes cylinder, an effective porous-medium law) is undecided.
  Nothing here constrains that choice, and the choice will matter more than the connector declaration
  does.
* **The audit in §6 is a reading of a call graph, not an execution.** Done by grep and by following
  constructor call sites on a CPU-only Mac. A runtime-registered implementation would be invisible to
  it; none was found, but absence by grep is weaker than absence by running. No GPU job was run — this
  session runs zero.
* **The protrusion-seam defect is a claim about a method signature**, which is the strongest form of
  this kind of claim available without running the code, and it is still a static claim. It says no
  object satisfying the declared interface can reach the fluid; it does not prove that no
  implementation exists somewhere that reaches the fluid by another route. None was found.
* **`microtubule_rig.py` is modified in the working tree.** The diff against `be0e5876` was inspected
  before the `mt_cytosol_transfer` verdict was written, but a reader reproducing this on a clean
  checkout will see a different version of that one file.
* **The G-actin field itself is not declared by this lane.** These two contracts name a G-actin
  transport stencil as an endpoint role; whether Aleph's cytosol owner actually declares a monomer
  field is another lane's question. If it does not, these two contracts name a role that has no owner,
  and no test here can catch that — it needs a cross-manifest check that does not exist yet. **This is
  the most likely way these six declarations are currently wrong.**
* **`composite_group` is `None` for all six** and `internal_to` is `None` for all six. Both follow from
  the manuscript extract. If either changes, the tests assert the old answer and will not catch it.
* **Nothing here is evidence about `ALEPH-DQ-104`.** Group D sets `commit_on_accept=False`, which
  presumes an acceptance predicate exists to be irrelevant to; which predicate is the PI's decision.
