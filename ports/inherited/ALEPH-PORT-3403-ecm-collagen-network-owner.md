# ALEPH-PORT-3403 — the `ecm` state owner: `CollagenNetwork` over the 2801 element library

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3403` |
| Lane | `L34 ecm state owner` |
| Status | `PROPOSED` |
| Written | `2026-07-31` (**before** the code, per PLAN §0.2.5 — see §14 for what the code then changed) |
| Port class | `RE-DERIVED` |

---

## 1. Aleph API

The exact public surface this entry authorises, all of it in `aleph/vertical/ecm_network.py`.
Nothing outside this list, and no other file, is covered.

```python
from aleph.vertical.ecm_network import (
    COMPONENT_NAME,
    CollagenNetwork,
    CollagenSiteEndpoint,
    ECMEndpointRole,
    ECMRoleError,
    ECMSiteType,
    ECMStaleHandleError,
    EpochStampedMaterialPointSink,
    FarFieldClamp,
    OWNED_STATE_KEYS,
    PUBLISHABLE_ENDPOINT_ROLES,
    ECM_ENDPOINT_ROLES,
    SheetGeometry,
    WorldBoundaryFrame,
    assert_far_field_frame_owns_no_state,
    assert_state_keys_disjoint_from,
    build_planar_collagen_sheet,
    owned_state_keys,
)
```

**This entry is the state owner and nothing else.** The three energy laws, the material cards, the
fibre/crosslink/ligand data types and the barycentric sink live in `aleph/vertical/ecm.py` under
`ALEPH-PORT-2801`, are gradient-checked to order 2.000 there, and are **imported, not re-derived**.
This module allocates the state, runs the accepted-step transaction, owns the topology epoch,
resolves ligands against a live network, remodels, evolves damage, clamps the far field and
publishes the endpoints. `aleph/vertical/ecm.py` is not modified by this entry.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — cited only because `ALEPH-PORT-2801` cites it; not resolved or read by this lane |
| Source path | `ffn_sim/ac/engine/ecm_mechanics.py`, `ffn_sim/ac/engine/ecm_world.py` — **named, not read** |
| Source symbol(s) | **none** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened by this lane. |
| Working tree == commit? | not applicable — nothing was read, so there is no revision to have mis-cited |

The Aleph inputs that *were* read, and which are the real provenance of this module: the registered
`ecm` contract in `aleph/state/census_environment_surface.py`; the four registered connector
contracts naming `ecm` in `aleph/state/connectors_surface_traction.py`; the design sections of
`aleph/vertical/ecm.py`'s docstring, which describe a `CollagenNetwork` that was never written and
which this module implements; `docs/design/ENDPOINT_CONTRACT.md` I1–I3; and
`aleph/vertical/cortex_filaments.py` as the canonical owner shape.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** An accepted-step transaction (`snapshot`/`commit`/
`rollback`), a monotone counter bumped on committed topology change, a barycentric material-point
resolver and a Dirichlet clamp on a boundary node set are all constructions any competent author
re-derives directly; per `ports/TEMPLATE.md` §3 that is a reason to write clean-room, not to port.
The one non-obvious design decision here — that ligand *material coordinates* are frozen at
construction so that rewriting rest lengths cannot move an adhesion — is stated in Aleph's own
`ecm.py` docstring and in the registered contract's `mechanical_role`, and is re-derived below.

## 4. Physical or mathematical law represented

**Mostly not a law. This is an ownership and transaction argument**, and saying so is the honest
answer §4 asks for. The three energies are `ALEPH-PORT-2801`'s and are unchanged. What is new:

**The accepted-step transaction.** A candidate step may move nodes, propose crosslinks, propose an
unbind, propose remodelling and propose damage. None of it is real until `commit()`. `rollback()`
restores the state bit-identically and discards every proposal. The reason is mechanical rather
than stylistic: a bond formed inside a candidate step that is then rejected leaves the network
stiffer than its own accepted history says it is, and no force check can see it, because the force
check is computed *from* the topology.

**The topology epoch.** `ecm_topology_epoch` advances on a committed crosslink bind, a committed
unbind, and committed remodelling — exactly the four events the registered contract names
("crosslink binding, unbinding, remapping and remodeling"). It does **not** advance on damage:
damage moves no node and rewires nothing, and bumping the epoch for it would invalidate every
outstanding handle on every step in which any segment yields.

**Ligands are material coordinates, and the invariance is now derivable rather than asserted.** Let
`s` be a ligand's rest arclength in the reference material frame fixed at construction, and let the
per-node table `S[i]` be that same frame. Resolution is `k = searchsorted(S, s) - 1`,
`w = (s - S[k]) / (S[k+1] - S[k])`. Remodelling rewrites `segment_rest_um` (the stress-free
lengths that the axial energy divides by) and damage rewrites the effective modulus; **neither
appears in the two lines above**. So the resolved `(node_a, node_b, w)` is bit-identical across a
remodelling commit. That is a *proof about which arrays the resolver reads*, and the way to test it
is to ship a flag that makes the resolver read the other array — `remodel_material_coordinates`,
which rewrites `S` too — and require the healthy branch to be bit-identical while the broken branch
moves the ligand.

**The far field is a clamp, not a spring.** `world_boundary` is registered `[B]`: rigid, no
compliance, no degrees of freedom, no snapshot. So the far field is a Dirichlet constraint on a
declared boundary node set: those nodes do not move in `SOLVE`, and the reaction is *reported*:

`R = − Σ_{i ∈ clamped} f_i`

with `f` the total assembled force (internal + connector). The identity that makes this checkable:
the internal forces sum to zero over all nodes (every term is a translation-invariant pair), so for
an external load `L` applied at free nodes, at equilibrium of the free nodes `R = −Σ L` exactly.
Giving the anchor a stiffness instead would hand the frame a compliance the contract forbids and
would make the far-field reaction a function of a spring constant nobody can measure.

**A 2D sheet is a declared geometry, not an emergent one.** `build_planar_collagen_sheet` lays warp
fibres at one offset along the plane normal and weft fibres at another, separated by exactly the
crosslink rest length, and crosslinks them at every intersection. The offset is not cosmetic: two
fibres crossing at the same point give a crosslink of zero length, whose direction is undefined and
which `crosslink_energy_and_forces` refuses. `SheetGeometry` records the plane normal axis, the
in-plane extent and the sheet thickness, so the dimensionality is a field that can be asserted
rather than a property of the numbers somebody typed.

## 4a. Mutation testing of the controls

Three mutants applied to `aleph/vertical/ecm_network.py`, the controls run against each, the source
restored and verified byte-identically. Filled in from the run; see §14 for what it found.

| # | Mutation | Tests killed | Survived? |
|---|---|---|---|
| M1 | `commit()`: bump `topology_epoch` on *every* commit rather than only on a committed topology change | **3** | no |
| M2 | `accumulate(SOLVE)`: drop `step[self.clamped_node_mask] = 0.0`, so the far field is declared and not applied | **2** | no |
| M3 | `propose_damage_from_strain`: `np.minimum` for `np.maximum`, so damage heals | **4** | no |
| M4 | `crosslink_arrays`: write the **a**-side interpolation weight into the b side | **0 — SURVIVED** | see below |

Source restored and verified **byte-identical** after each: the SHA-256 of the restored file equals
`308b367a6c2ac20401483043f14079493d873aade7ac48b0518d16d29b016484`, the hash of the pre-mutation
copy, and the file is otherwise untouched.

**The exercise found a real gap, which is the reason to run it.** M4 was added as a fourth probe on
the suspicion that it would survive, and it did: **all 90 controls passed with the b-side endpoint of
every crosslink written into the flat arrays with the wrong interpolation weight.** The reason is
specific and worth stating, because it generalises to every owner in this repository that flattens
resolved material points into element arrays:

* the error is *internally consistent*. The same wrong weight builds the endpoint position and
  scatters the gradient, so `F = −grad E` still holds **exactly** and every finite-difference control
  passes at order 2;
* the four scattered forces still sum to zero and still carry no net moment about the origin, so
  force and moment closure pass **exactly**;
* the pin against the element library passes, because the pin hands
  `crosslink_energy_and_forces` the same wrong arrays the owner uses.

What is actually wrong is that the load arrives at the wrong place along the fibre — a purely
geometric fact, invisible to every conservation identity. The control written to kill it is
`::TestTheSheetIsTwoDimensionalAndSaysSo::test_the_freshly_woven_sheet_is_stress_free`: a freshly
built weave joins two points that are *at the intersection*, exactly one sheet thickness apart, so it
must store **zero** crosslink energy at construction. With M4 applied that control fails and the
other 91 still pass, which is the shape of the finding. Also added:
`::test_a_resolved_material_point_sits_at_the_arclength_it_names`, which measures the resolved
position against a ruler rather than against the resolver's own weights.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| node position | µm | 1e-6 m | finite |
| segment rest length | µm | 1e-6 m | `> 0`, strictly |
| material coordinate `s` | µm | 1e-6 m | `[0, L_rest]` of its fibre, at construction |
| damage `d` | dimensionless | — | `[0, max_damage]`, `max_damage < 1` |
| remodelling fraction | dimensionless | — | `(0, 1]` |
| force | pN | 1e-12 N | finite |
| energy | pN·µm | 1e-18 J | `>= 0` |
| `topology_epoch` | count | — | monotone non-decreasing within an accepted history |
| mobility | µm/(pN·s) | — | `> 0` |

Singular and boundary cases, each with the behaviour Aleph requires:

- **`rollback()` with no snapshot**: refused with `ECMTopologyError`. A rollback that silently does
  nothing is indistinguishable from one that worked.
- **A crosslink joining a fibre to itself**: refused at proposal with `ECMTopologyError` (the
  element library refuses it at construction; the owner refuses it again at the point of proposal so
  the message names the network).
- **A duplicate crosslink id**, live or pending: refused. Two links with one id cannot be unbound.
- **`propose_unbind` of an unknown or already-pending id**: refused with `ECMTopologyError`.
- **A material coordinate outside its fibre**: refused, not clamped. Clamping an attachment to the
  end of a fibre moves a load path without telling anybody.
- **An endpoint role this owner does not publish**: refused with `ECMRoleError`.
- **A handle used after a committed topology change**: refused with `ECMStaleHandleError`.
- **A handle block whose sinks were issued at different epochs**: refused, not reconciled.
- **A scattered force of the wrong shape or a non-finite one**: refused.
- **Fibre node ranges that do not tile `[0, N)`**: refused with `ECMOwnershipError`.
- **A far-field frame object carrying state**: refused with `WorldBoundaryStateError`.
- **A remodelling fraction outside `(0, 1]`**: refused. Zero is a no-op dressed as an event and a
  fraction above one drives the rest length past the current length, which is not relaxation.
- **A sheet whose nodes do not lie within the declared thickness of the declared plane**: refused.

Invariants, each with the test that asserts it (all in
`tests/vertical/test_ecm_network_controls.py` unless stated):

- **I1.** `F = −grad E` per term on the owner's own assembled state, order > 1.6 —
  `::TestEachTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order`
- **I2.** The owner's terms are bit-identically the element library's functions —
  `::TestTheOwnerDrivesTheElementLibraryAndDoesNotRestateIt::test_each_owner_term_is_bit_identical_to_the_element_library_call`
- **I3.** Force and moment closure against the **constituent** scale, plus translation invariance,
  determinism and purity — `::TestClosureAndInvarianceOnTheIsolatedSheet`
- **I4.** `rollback()` restores all seven state blocks bit-identically —
  `::TestTheAcceptedStepTransaction::test_rollback_restores_every_state_block_bit_identically`
- **I5.** The epoch advances only on a committed topology change —
  `::TestTheAcceptedStepTransaction::test_the_epoch_advances_only_on_a_committed_topology_change`
- **I6.** A ligand resolves bit-identically across remodelling —
  `::TestLigandsSurviveRemodelling::test_the_resolved_ligand_is_bit_identical_across_a_remodelling_commit`
- **I7.** Declared state keys equal the registered contract, and all seven are actually allocated —
  `::TestTheDeclarationsAgreeWithTheCode::test_the_owned_state_keys_are_exactly_the_registered_contract`,
  `::TestTheDeclarationsAgreeWithTheCode::test_every_declared_state_block_is_actually_allocated`
- **I8.** The far-field reaction equals minus the applied load at equilibrium —
  `::TestTheFarFieldIsAClampAndNotASpring::test_the_reaction_is_minus_the_load_the_constraint_removed`
- **I9.** No state key of this owner is declared by any other registered owner —
  `tests/state/test_census_key_ownership.py::test_no_state_key_has_two_registered_owners`

## 6. Source evidence class and known retractions

No claim is made about the reference implementation's ECM: none of it was read by this lane, and
`ALEPH-PORT-2801` §6 records the same for the element library. Where I looked: nowhere inside
`/Users/sw1/ffn_cellsim`.

Aleph-side retractions that **do** bear on this entry, and which is why they are listed:

- `ALEPH-PORT-2801` §14 is a self-retraction: it records that `ecm.py` declares seven state blocks
  and allocates none, that `__all__` exported five symbols that never existed (including
  `CollagenNetwork` and `FarFieldAnchor`), and that its docstring described them in the present
  tense. This entry exists to close that gap; §14 below states exactly how much of it is closed.
- `ALEPH-PORT-2901` (`focal_adhesion`) is `REJECTED`, and its archived code is marked do-not-copy.
  Nothing here reads it.
- The registered `ecm` contract carries `citation_status: UNSOURCED` and this entry does not
  improve it. `default_collagen_card()` returns placeholder magnitudes and says so.

## 7. Independent oracle or derivation

Four, none of which is the reference implementation:

1. **Central finite differences** of the owner's own assembled energy, per term. Independent of
   every constant in the law.
2. **The element library itself, called directly.** The owner's `axial_term()` must be
   bit-identical to `ecm.axial_energy_and_forces` on the owner's published arrays. This is the pin
   that stops the owner and the library drifting into two models of one matrix — the failure
   `ALEPH-PORT-2801` §14 names when it says an element library with no owner "obliges every
   consumer to reassemble the network by hand".
3. **Global equilibrium at the clamp.** `R = −Σ L` is an identity of statics, not of this code, and
   it is checked at a relaxed equilibrium rather than at construction.
4. **Bit-identity of the ligand resolution across remodelling**, which is a statement about which
   arrays the resolver reads and is checked with `==` on the node pair and the weight rather than
   with a tolerance.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_ecm_network_controls.py::TestEachTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order` | Per term, `energy > 0` first, then max&#124;F+∇E&#124;/max&#124;F&#124; at 4e-4/2e-4/1e-4 µm with observed order > 1.6 |
| Positive | `tests/vertical/test_ecm_network_controls.py::TestTheOwnerDrivesTheElementLibraryAndDoesNotRestateIt::test_each_owner_term_is_bit_identical_to_the_element_library_call` | `np.array_equal` between each owner term and the `aleph.vertical.ecm` function on the owner's own arrays |
| Positive | `tests/vertical/test_ecm_network_controls.py::TestTheAcceptedStepTransaction::test_rollback_restores_every_state_block_bit_identically` | All seven registered blocks restored with `==` after a step that moved nodes, bound a crosslink, remodelled and damaged |
| Positive | `tests/vertical/test_ecm_network_controls.py::TestTheFarFieldIsAClampAndNotASpring::test_the_reaction_is_minus_the_load_the_constraint_removed` | At a relaxed equilibrium under an applied load, `R + ΣL` is below 1e-9 of the constituent force scale |
| Positive | `tests/vertical/test_ecm_network_controls.py::TestLigandsSurviveRemodelling::test_the_resolved_ligand_is_bit_identical_across_a_remodelling_commit` | `(node_a, node_b, weight)` identical with `==` across a commit that rewrote every rest length |
| Positive | `tests/vertical/test_ecm_network_controls.py::TestTheSheetIsTwoDimensionalAndSaysSo::test_the_freshly_woven_sheet_is_stress_free` | A freshly built weave stores < 1e-20 pN·µm of crosslink energy and every crosslink spans exactly one sheet thickness to 1e-12 µm — the control the mutation exercise produced |

**Measured, this session, on this tree** (`/Users/sw1/miniconda3/envs/aleph/bin/python`, numpy
2.5.1). Every number below was produced by a run in this session; none is quoted from memory. The
configuration is the woven 3x3 sheet of `_loaded_sheet`: 30 nodes, 24 segments, 18 bending stencils,
9 crosslinks, perturbed off every symmetry from a pinned seed, with damage committed on 10 of 24
segments (max 0.1649) and held fixed.

Finite-difference agreement, steps 4e-4 / 2e-4 / 1e-4 µm:

| Term | Energy (pN·µm) | Force scale (pN) | rel. max&#124;F+∇E&#124; at 1e-4 | observed order |
|---|---|---|---|---|
| axial | 33.214900 | 125.749 | 2.293e-09 | **2.000** |
| bending | 2.771134 | 7.046 | 6.786e-09 | **2.000** |
| crosslink | 1.861768 | 3.166 | 4.046e-08 | **2.000** |
| total | 37.847803 | — | 1.959e-09 | — |

Closure on the isolated sheet (resultant ÷ constituent scale):

| Term | force | moment |
|---|---|---|
| axial | 4.069e-19 | 1.309e-17 |
| bending | 8.242e-18 | 2.631e-18 |
| crosslink | 1.725e-17 | 1.175e-17 |
| total | **1.748e-17** | — |

The far field, after 8,000 overdamped steps at `dt = 1e-3` under a point load `(2, −1, 0)` pN
delivered through a published endpoint:

| Quantity | Measured |
|---|---|
| free-node residual | 2.706e-10 pN |
| free-node residual ÷ constituent scale (6.290 pN) | 4.302e-11 |
| reaction `R` | `(−2.000000000, +1.000000000, −6.668e-11)` pN |
| &#124;R + ΣL&#124; ÷ constituent scale | **2.500e-11** |

## 9. Deliberately failing negative control

Every break is a flag on the **shipped** owner, default `False`, so the control drives the real code
path rather than a hand-edited copy. Each is asserted in **both** directions: the broken branch
misbehaves and the healthy branch is exactly inert.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_ecm_network_controls.py::TestTheNegativeControlsFail::test_a_flipped_axial_sign_is_caught_by_the_gradient_test` | `flip_axial_sign=True` drives the FD disagreement above 0.5 relative; healthy branch below 1e-7 |
| Negative (must fail) | `tests/vertical/test_ecm_network_controls.py::TestTheNegativeControlsFail::test_a_crosslink_committed_without_an_accepted_step_survives_a_rollback` | `commit_topology_without_accept=True` leaves a bond the rejected step never accepted; healthy branch rolls back to exactly the committed count |
| Negative (must fail) | `tests/vertical/test_ecm_network_controls.py::TestTheNegativeControlsFail::test_a_stale_handle_scatters_silently_when_the_epoch_guard_is_off` | `allow_stale_handles=True` scatters into a network that has been rewired; healthy branch raises `ECMStaleHandleError` |
| Negative (must fail) | `tests/vertical/test_ecm_network_controls.py::TestLigandsSurviveRemodelling::test_remodelling_the_material_frame_moves_the_ligand_which_is_the_defect` | `remodel_material_coordinates=True` moves the resolved ligand; healthy branch is bit-identical |
| Negative (must fail) | `tests/vertical/test_ecm_network_controls.py::TestTheMissingTermControl::test_omitting_the_crosslink_law_makes_a_stretched_network_free` | `omit_crosslink_law=True` gives **exactly** zero energy difference on a configuration that only crosslinks resist; healthy branch is strictly positive |
| Negative (must fail) | `tests/vertical/test_ecm_network_controls.py::TestTheNegativeControlsFail::test_a_cable_carrying_compression_is_visible_at_the_network_level` | `allow_cable_compression=True` makes a slack sheet store energy; healthy branch stores exactly `0.0` |

## 10. Numerical and precision envelope

Float64 throughout; no reduced-precision path exists.

Finite-difference steps 4e-4 / 2e-4 / 1e-4 µm on coordinates of order 1–8 µm. The central-difference
floor is `eps^(1/3) ≈ 6e-6` **relative**; the smallest step used is well above it. Below the floor
cancellation dominates, the apparent order goes negative, and a correct gradient is
indistinguishable from a wrong one.

Tolerances, and why each is not looser:

- Per-term FD agreement asserted at `1e-6` relative. Measured worst case **4.0e-08** (crosslink),
  best **2.3e-09** (axial) — a factor of 25 of headroom. The flipped-sign break sits at **2.0000**,
  seven orders above the gate, so the negative control cannot pass by accident.
- Observed FD orders asserted at `> 1.6`; measured **2.000** for all three terms.
- Force and moment closure asserted at `1e-13` of the **constituent** scale (`Σ|f_i|`,
  `Σ|r_i||f_i|`), never the resultant. Measured worst residual **1.7e-17** (force) and **1.3e-17**
  (moment). The gap is deliberate: closure is an identity of the algebra and its residual is set by
  summation order rather than by the physics. A resultant-based tolerance would get *stricter the
  more correct the physics is* and once rejected 40,000 consecutive steps of a good relaxation.
- The far-field identity asserted at `1e-9` of the constituent force scale, and the free-node
  residual at `1e-8`; measured **2.5e-11** and **4.3e-11**. The gap here is *not* deliberate slack —
  the identity is exact in the algebra, but the equilibrium it is evaluated at is reached by
  relaxation and is not exact, so the tolerance is set by how long that relaxation is run.
- Explicit overdamped stability: the update is stable while `dt·μ·k_node < 2`, and the stiffest node
  carries two segments of `600 / 1.5 = 400` pN/µm plus a crosslink, giving a bound near `2.5e-3` s.
  Measured: at `dt = 1e-3` the residual falls geometrically (9.2e-04 at 1,000 steps, 7.6e-06 at
  3,000, 2.7e-08 at 6,000, 6.3e-12 at 12,000); at `dt = 2e-3` it does not fall at all. The control
  runs 8,000 steps, and the margin is measured rather than assumed.
- Bit-identity claims (rollback, ligand resolution, purity, determinism) are asserted with
  `np.array_equal` and `==`, never `allclose`. There is no tolerance to loosen.

Outside the envelope: every singular case in §5 refuses with a typed exception rather than
degrading. The exception, and it is inherited: `default_collagen_card()`'s magnitudes are
placeholders, so this module will happily integrate a matrix whose stiffness is not collagen's.

## 11. Production-backend residency and transfer

Host, numpy, float64. The owner takes its `Backend` from `StepContext` and uses it for the
accumulate/scale/max-abs operations in `accumulate`, exactly as `cortex_filaments.py` does, so the
witness counter advances and a participant that did nothing cannot pretend otherwise. The element
laws themselves are numpy and are `O(S)`/`O(T)`/`O(C)` gather–scatter kernels with no host
round-trip per step. **No GPU work of any kind was run by this lane**, and no GPU authorization was
sought or held.

## 12. Comments and docstrings to discard

Nothing to discard from the reference: no reference file was read.

What was deliberately **not** carried over from Aleph's own prose, which is the interesting part:

- `ecm.py`'s design sections name a class `FarFieldAnchor`. This module ships `FarFieldClamp`
  instead, because `ecm_far_field_anchor` is census item ② and it is **undecided whether it is a
  connector at all**; a class named after the registry row would read as that row being
  implemented. The rename is a deviation from `ALEPH-PORT-2801`'s design prose and is recorded here
  rather than made silently.
- `ecm.py`'s `IMPLEMENTATION STATUS` block says the owner does not exist. It is **left exactly as
  it is**: `ecm.py` belongs to another lane and this entry does not modify it. So the tree now
  carries a stale sentence in `ecm.py`, and that is named in §14 rather than fixed across an
  ownership boundary.
- No `ALEPH-PD-*` file, no `decided_by` field, and no claim of PI ratification appears anywhere in
  this work.

`tests/vertical/test_ecm_network_controls.py::TestTheDeclarationsAgreeWithTheCode::test_no_reference_project_vocabulary_survives`
asserts the module text carries none of `ffn_cellsim`, `ffn_sim`, `dcm_contact`, `surface_body`, and
`::test_the_module_imports_nothing_from_validation` asserts the `aleph/** → validation/**` firewall.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** 92 controls pass (`tests/vertical/test_ecm_network_controls.py`, 2026-07-31, 1.95 s), together with the 60 element-library controls of `ALEPH-PORT-2801` unchanged and `tests/state/` green. Four mutants applied, killing 3 / 2 / 4 / **0**; the fourth survived and is the finding in §4a; the control written to kill it now does. Source restored byte-identically after every one. Status stays `PROPOSED`: no PI review, no ratification, the magnitudes are still `UNSOURCED`, and this owner has never been stepped inside `aleph/vertical/assembly.py`. |
| Reviewer | Agent-proposed (lane L34). Unratified. No PI review. |
| Rollback | Delete `aleph/vertical/ecm_network.py` and `tests/vertical/test_ecm_network_controls.py`. Nothing breaks: no module in the tree imports either, no `Binding` was added to `aleph/vertical/wiring.py`, and `aleph/vertical/ecm.py` is byte-unchanged. |

## 14. Honest limits — what this entry does NOT establish

**Wired connector count contributed: zero.** `aleph/vertical/wiring.py` is untouched, deliberately —
a `Binding` row is added in the same commit that wires a connector, and wiring is the lead session's
lane. What this module changes is that the `ecm` side of `membrane_ecm_contact`,
`integrin_collagen_clutch` and `ecm_crosslink` is now *reachable*: it publishes a
`CoupledSiteEndpoint`-shaped handle and an `EndpointHandle` factory. Reachable is not wired, and the
count in `wiring.py` is still the honest one.

**`ecm_far_field_anchor` is implemented as a boundary condition and NOT as a connector.** Census
item ② is undecided — whether that registry row is a connector at all, and `HANDOFF.md` records that
finalizing it would require moving two `tests/firewall/` assertions, which `CLAUDE.md` §2 rule 4
forbids doing to unblock oneself. So `PUBLISHABLE_ENDPOINT_ROLES` **withholds** the far-field role,
`FarFieldClamp` is a Dirichlet constraint on this owner's own nodes, and nothing here presents it as
a load path between two owners. If the PI decides it *is* a connector, this clamp is what its `ecm`
side would be built from, and the withholding becomes wrong in exactly one line.

**`world_boundary` is a frame with no state and no behaviour.** `WorldBoundaryFrame` holds a pose,
owns nothing, and `assert_far_field_frame_owns_no_state` refuses any object presented as the frame
that carries state. That is the whole of it. It is **not** a `world_boundary` owner module, it is
not registered as a participant, and it receives nothing: the reaction is reported by the clamp and
is deliberately not stored, because a frame that accumulates a reaction has state and a reference
frame with state is a body.

**Damage evolution now exists but its envelope is not physics.** `FibreCard`'s
`damage_onset_strain` / `damage_failure_strain` linear softening envelope is *implemented*,
irreversible and monotone, and committed only on an accepted step. Nothing about the envelope's
*shape* is sourced: linear softening between two strains is a modelling choice made here, no
collagen measurement stands behind either strain, and no rate dependence, no fatigue and no healing
exist. The mechanics is tested; the material is `UNSOURCED`.

**Remodelling is one law and it is the simplest one.** `segment_rest_um += fraction × (L − L0)`, a
uniform relaxation of the stress-free length toward the current length, applied on accept. It is
irreversible and it is not derived from a protease model, a fibre-recruitment model or a
viscoelastic one — `ALEPH-PD-004` excludes matrix chemistry, so remodelling here is mechanical state
and any timescale it produces belongs to `fraction` and the step controller, not to biology.

**No fluid, and therefore no rate dependence that means anything.** Per `ALEPH-PD-004` there is no
interstitial phase; any rate dependence this network shows comes from the step controller.

**No steric exclusion between fibres.** Two fibres of this sheet can pass through each other at zero
cost. The cortex owner has a WCA term; this one does not, and the sheet's crosslinks are what hold
its topology. Absent, not lumped.

**No collision with the cell, no adhesion kinetics, no traction measurement.** This owner publishes
where an adhesion may attach. Whether one attaches, with what lifetime, and what traction it
transmits belongs to `focal_adhesion` and the clutch connector.

**The 2D sheet has a thickness and that is deliberate.** Warp and weft are offset along the plane
normal by exactly the crosslink rest length, because a zero-length crosslink has no direction and is
refused. So "planar" means "planar to within the declared `sheet_thickness_um`", it is a declared
field rather than an implied property, and a control asserts it. A genuinely zero-thickness sheet
would need a crosslink law with a rest length of zero, which is a different element.

**This owner has never been stepped inside a world.** It is not in `aleph/vertical/assembly.py`, it
has never run beside another owner, and `accumulate` has only ever been called by its own controls.
Everything above is a statement about one participant in isolation.

**A stale sentence now exists in `aleph/vertical/ecm.py`.** Its `IMPLEMENTATION STATUS` block says
"the `ecm` state owner does not exist yet", which stopped being true when this module landed.
`ecm.py` is another lane's file and this entry does not edit it. Reported, not fixed.

**Physical magnitudes are not established.** Evidence rung `ANALYTIC_ORACLE`, quantitative status
`BLOCKED`. No modulus, rigidity, crosslink stiffness, damage threshold, remodelling rate or sheet
dimension here is a collagen number.

**A control suite is not a simulation, and 92 green controls are 92 statements about one participant
evaluated by its own test file.** Nothing here says the matrix is right; it says the owner computes
what it declares it computes, keeps what the registry says it keeps, and refuses what it says it
refuses.

**Not run:** any GPU job (zero, and no authorization was sought or held), and the whole-repo suite.
Only `tests/vertical/test_ecm_network_controls.py`, `tests/vertical/test_ecm_controls.py`,
`tests/state/` and `tests/ports/` were exercised, because other sessions are writing this tree
concurrently and `CLAUDE.md` §1 says to report a foreign breakage rather than conflate it with mine.

**Two guards are red because this work exists in the working tree and not in the index, and both
are guards doing their job.** Neither was edited to unblock this lane (`CLAUDE.md` §2 rule 4), and
both are repaired by the commit that lands the module rather than by anything here:

- `tests/vertical/test_package_exports.py::test_every_module_present_is_listed` — `ecm_network` is
  a module in `aleph/vertical/` and is not in that package's `__all__`. **The landing commit must
  add it there**, and must do so in the same commit, because the sibling control
  `::test_every_listed_name_is_a_module_git_actually_tracks` fails if the name is listed while the
  file is untracked. `aleph/vertical/__init__.py` is a shared file and outside this lane's
  allocation, so it was not edited here.
- `tests/ports/test_index_names_only_tracked_files.py` — see below.

**One test is red because of this entry, and it is red for a reason that resolves on commit.**
`tests/ports/test_index_names_only_tracked_files.py` — a guard added by a concurrent lane during
this session — fails because `ports/ledger/INDEX.md`, which that lane regenerated, now names this
entry while the entry is still untracked. It goes green when the entry is committed. `INDEX.md` is
a shared generated file that this lane does not own, so it was not regenerated here; the three
`tests/ports/` failures that were present when this work started (`test_no_orphan_ports`,
`test_named_controls_resolve_to_real_tests`, `test_index_is_not_stale`, all from the
`ALEPH-PORT-3404` / `aleph/scenarios/` lane) had cleared by the time it finished.
