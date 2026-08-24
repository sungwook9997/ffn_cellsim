# ALEPH-PORT-3202 — the microtubule owner as an explicit dynamic rod graph

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3202` |
| Lane | `L32 microtubule / internal frame` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (**before** the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
from aleph.vertical.microtubule import (
    MICROTUBULE_ENDPOINT_ROLES,
    MTOC_STATE_KEYS,
    OWNED_STATE_KEYS,
    PUBLISHABLE_ENDPOINT_ROLES,
    ROLE_SITE_KINDS,
    CaptureState,
    MaterialPoint,
    MaterialPointSink,
    MicrotubuleAttachmentRole,
    MicrotubuleAttachmentSite,
    MicrotubuleCard,
    MicrotubuleNetwork,
    MicrotubuleOwnershipError,
    MicrotubuleRateUnavailableError,
    MicrotubuleRoleError,
    MicrotubuleStaleHandleError,
    MicrotubuleTopologyError,
    PlusEndPhase,
    Rod,
    SiteKind,
    assert_state_keys_disjoint_from,
    assert_the_mtoc_is_never_an_endpoint,
    assert_the_rod_can_carry_compression,
    axial_energy_and_forces,
    bending_energy_and_forces,
    default_microtubule_card,
    discrete_euler_critical_load_pn,
    minus_end_anchor_energy_and_forces,
    owned_state_keys,
    straight_rod_nodes,
)
```

This entry authorises **a state owner**, not an element library. `MicrotubuleNetwork` allocates all
seven registered state blocks, participates in the accepted-step transaction, and publishes endpoint
handles. It wires **zero** connectors: `aleph/vertical/wiring.py` is not touched by this lane and no
`Binding` row is added. §14 records exactly what is therefore unestablished.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — the revision `ALEPH-PORT-1801` audited |
| Source path | `ffn_sim/ac/engine/microtubule_rig.py`, `ffn_sim/ac/solid/microtubule.py`, `ffn_sim/ff/microtubule.py` — **named, not read** |
| Source symbol(s) | **none** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened by this lane, at either revision. |
| Working tree == commit? | not applicable to this lane — nothing was read, so there is no revision to have mis-cited. `ALEPH-PORT-1801` §2 records that the two revisions differ by one comment block in `microtubule_rig.py` and that the other two paths are clean. |

The paths are recorded so a later auditor can establish independently whether Aleph's re-derivation
agrees or disagrees with the reference. Listing them is **not** a claim that they were consulted.
This lane's only facts about the reference are **inherited from `ALEPH-PORT-1801`**, which did read
both revisions, and they are used here as an Aleph-local document rather than re-derived:

- the reference's plus-end kinetics are a real stochastic behaviour with commit/rollback;
- its **bending** path is self-labelled `STATIC_BENDING_ADAPTER_DYNAMIC_TOPOLOGY_PENDING` — a static
  bending adapter over a frozen bead-chain aster, with the dynamic-topology case declared pending;
- its transition rates are unset by design, as priors to be inferred.

The middle one is the fact that shapes this module, and it is stated as inherited rather than as
measured here: bending and dynamic topology are **not consistent with each other** in that tree.
This module's whole difficulty is that they must be, because a rod whose node count changes and
whose bending stencil is rebuilt on the change is exactly the case that adapter defers.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** `ALEPH-PORT-1801` §3 already reached that verdict for this
component and this entry inherits it rather than re-arguing it. Restated for the mechanics, which
1801 did not authorise:

- A two-sided segment spring with a rest-length normalisation, a tangent-difference bending stencil,
  and a Hookean central-force anchor between a rod base and a hub point are laws any competent author
  re-derives in an afternoon. Per `ports/TEMPLATE.md` §3 that is a reason to write clean-room.
- The candidate/commit/rollback split for a topology change is genuinely good and Aleph already owns
  it — `ALEPH-PORT-305`. Taking it now would import an architecture Aleph has re-derived better.
- The one thing worth taking is a *defect study*, not code: a bending adapter validated only at a
  frozen topology, shipped in a rig whose topology is dynamic. That is recorded in §4 as the
  configuration this module's controls are built to detect, and it cost zero lines to inherit.

**Recorded as prior art, unported** (both already named by 1801, repeated because they are now
load-bearing in code rather than in a contract): a `topology_generation` counter so a length change
is an observable committed event rather than a quiet array resize; and keeping transition rates out
of the component.

## 4. Physical or mathematical law represented

Three energy terms with their exact analytic gradients, plus one topology discipline that is not a
law and is named as such. Units µm / pN / pN·µm throughout.

### 4.1 Axial — two-sided, because this is the compression-capable member

`E = Σ_s (k_s / 2) (L_s − L0_s)² / L0_s`, giving segment tension `T = k (L − L0) / L0` which pulls
the two ends together. Dividing by the rest length makes the modulus resolution-independent: a rod
cut into twice as many segments has twice as many springs each twice as stiff, and the same
end-to-end stiffness.

**Two-sided is the physics claim of this component and not a default.** The registered contract calls
the microtubule "the compression-capable member of the internal frame", and the `mt_nucleus_linc`
contract says it is "the only LINC path that can deliver a compressive nuclear load". A one-sided
(tension-only) axial law would leave that sentence with no implementation: a compressed rod would
store exactly `0.0` and deliver exactly `0.0`, and the nucleus would receive no compressive load
through the one edge declared to carry it. The break flag `tension_only_axial` reproduces that
configuration on the shipped code path so a control can drive it in both directions.

### 4.2 Bending — the term the buckling oracle interrogates

On unit tangents `t1`, `t2` of the two segments meeting at interior node `i`:
`E_i = (κ_i / 2h_i) |t2 − t1|²`, with `h_i` the **rest** Voronoi length `(L0_prev + L0_next)/2`.

Because `h` is a rest length, a uniformly stretched or contracted straight rod has exactly zero
bending energy and exactly zero bending force, so bending and stretching are decoupled and a buckling
measurement is a statement about `κ` alone. Small-angle limit: `|t2 − t1| = 2 sin(θ/2) → θ → C h`, so
`E_i → (κ/2) C² h` and the sum tends to `∫ (κ/2) C² ds` — Euler–Bernoulli with coefficient `κ`.

**This is the same discrete convention `aleph/vertical/cortex_filaments.py` and `aleph/vertical/ecm.py`
carry, and the agreement is deliberate and asserted rather than assumed.** Aleph should have one
bending convention, not three that disagree about a sharply bent rod. The stencil here is written
independently (this module imports nothing from either) and is pinned to `sf_arc`'s numerically by a
control.

### 4.3 The buckling oracle, re-derived for this stencil

A rod with bending rigidity that can carry compression is a buckling problem, and a critical load is
an **eigenvalue of the whole operator**, so an accidental rescaling of one term cannot reproduce it
the way it can reproduce a single energy value. `ALEPH-PORT-2801` §4 derives this for the ECM
stencil; the derivation below is carried out again here for this module's stencil, and the
conclusion — **the conventions match** — is a result of the derivation rather than an assumption
imported with the formula.

Take a straight pinned–pinned chain along `x` with rest spacing `h`, `N` nodes, `L = (N−1)h`, and a
small transverse displacement `y_i`.

*Bending.* To first order in `y`, `t2 − t1 = [(y_{i+1} − y_i)/h − (y_i − y_{i−1})/h] ŷ`, so

    E_bend = Σ_i (κ / 2h) · (y_{i+1} − 2y_i + y_{i−1})² / h²  =  (κ / 2h³) Σ_i (δ²y)_i².

For the fundamental mode `y_j = A sin(π x_j / L) = A sin(π j h / L)` the second difference is exact:
`(δ²y)_i = −4 sin²(πh/2L) · y_i`. Hence `E_bend = (8κ/h³) sin⁴(πh/2L) Σ y².`

*Work of the compressive end load `F`.* Each segment shortens in projection by
`√(h² + Δy²) − h ≈ Δy²/2h`, so the load does work `(F/2h) Σ_s (y_{s+1} − y_s)²`. For the same mode
`Σ_s (y_{s+1} − y_s)² = 4 sin²(πh/2L) Σ y²`, giving `W = (2F/h) sin²(πh/2L) Σ y².`

*Equating* the two:

    F_c = (4κ / h²) · sin²(πh / 2L).

Expanding `sin²(x) = x²(1 − x²/3 + …)` at `x = πh/2L` gives `F_c → π²κ/L²` **from below**, with
relative deficit `−π²h²/(12 L²)`, i.e. second order in `h`. Both numbers are asserted: checking only
the continuum limit lets a wrong `h` convention hide inside the discretisation error, and checking
only the discrete formula lets a wrong continuum coefficient survive.

`discrete_euler_critical_load_pn(kappa, length_um, node_count)` is the closed form, shipped in the
module so a consumer can compare against it without re-deriving, and exercised by the controls
against a measured eigenvalue.

**Stated as required: this module's discrete convention is the same as `ALEPH-PORT-2801` §4's**, and
the identical critical-load formula is the evidence, not the premise. If a future change to either
stencil breaks that agreement, the two entries will disagree about `F_c` at the same `κ, L, h`, which
is a detectable failure rather than a silent divergence.

### 4.4 The minus-end anchor, and the zero it can only report

`E = Σ_r (k_r / 2) (|p_base,r − X| − a_r)²` between rod `r`'s minus-end node and the MTOC hub point
`X`. Central force: the force on the base is along `X − p_base` and the reaction on the hub is its
exact negative, so the pair injects neither net force nor net moment into the world.

**`mtoc_torque` can only contain zeros, and this module says so in the value rather than in a
comment.** The reaction on the hub acts *at the hub point*, so its moment about the hub is
`(X − X) × f = 0` identically — for any central `U`, not only a Hookean one, at any extension.
`ALEPH-PORT-1801` §4(b) derives this from the contract side and the registered `centrosome_mtoc`
entry lists "MTOC torque balance" among its unsupported claims. Here the consequence is
implemented honestly: `mtoc_torque()` returns an exact `(3,)` zero vector and its docstring names the
moment-free construction as the reason. A control asserts it is exactly zero **and** that the
accompanying `mtoc_force()` is not, so "the hub is inert" cannot pass as "the hub is moment-free".

### 4.5 The topology discipline, which is not a law

The graph is dynamic in the strong sense: **the number of active nodes on a rod is state**. A rod
carries a build-time `capacity` and a live `active_node_count`; segments, bending triples and rest
lengths are derived from the active prefix and rebuilt when it changes.

A length change is proposed during a candidate step and applied only by `commit()`, which bumps
`topology_generation`. That counter is the whole reason the registered contract declares it: it makes
a length change a *committed topology event* rather than a quiet array resize. An endpoint handle
records the generation it was issued at and refuses to scatter afterwards, so a connector cannot
write into a rod whose node set has moved underneath it. The break flag `allow_quiet_topology_resize`
applies the resize without bumping the counter — the reference's permanently-inert-epoch
configuration — so the control that detects it drives the shipped path.

`plus_end_phase` is **read**, not stored and ignored: a rod may add a node only while `GROWING` and
remove one only while `SHRINKING`. A phase flip does not bump `topology_generation`, because the node
set does not move; it increments `catastrophe_rescue_state`.

### 4.6 What is deliberately absent: every rate

The registered contract is explicit — "Any transition rate — catastrophe, rescue, capture — is an
inference target of the run, not a constant this contract supplies." So this module ships the state
and the committed-transition machinery and **no rate value at all**. There is no default, no
placeholder and no zero. `MicrotubuleNetwork.catastrophe_rate_per_s()` and `.rescue_rate_per_s()`
raise `MicrotubuleRateUnavailableError`, naming the contract clause, so a consumer that needs a rate
**refuses** rather than silently running at a number nobody sourced. A zero default would be worse
than an exception: a microtubule that never catastrophes is a specific and flattering model, and it
would arrive without anybody choosing it.

## 4a. Mutation testing of the controls

Three mutants applied to `aleph/vertical/microtubule.py`, the control suite run against each, then
the source restored and verified. Recorded here because a control suite no mutant can break is
decoration. Every number below was produced by a run in this session.

| # | Mutation | Tests killed | Notes |
|---|---|---|---|
| M1 | `bending`: rest Voronoi `h` → **current** segment lengths (`kappa / h` → `kappa / ((|u|+|v|)/2)`) | **4** | per-term FD (bending), total-force FD, the bent-rod scale-invariance control, and the cross-owner agreement with `sf_arc`'s stencil |
| M2 | `axial`: drop the `/ rest` normalisation (`k Δ²` for `k Δ²/L0`) | **5** | all three buckling-load controls, the buckling convergence control, and the resolution-independence closed form. `F = −grad E` stays exactly true, so **no finite difference sees it at all** |
| M3 | `commit`: bump `topology_generation` on **every** commit, not only on a committed length change | **2** | the generation-advances/stale-handle control and the quiet-resize negative control |

All three killed, none survived. Source restored and verified **byte-identical** — SHA-256
`3a5cf8d0f56aea2f8a10f5ccb3d9dee7a5737e3a292e68290f7ebaf1abb5a677` before and after all three, on
the final revision of the module.

**The exercise found a real gap, which is the reason to run it.** M2 was originally killed by only
**four** tests, and the surviving one was
`::TestTheBendingLawAgainstEulerBuckling::test_the_critical_load_matches_the_discrete_prediction[5]`.
The reason is exact rather than statistical: the buckling rod was `L = 4 µm`, so the coarsest chain
had `h = 4/4 = 1.0` **exactly**, and dividing by a rest length of 1.0 is the identity. That
parametrization was therefore blind to a modulus that had stopped being resolution-independent,
while reading in the report as a passing buckling control.

Fixed by moving the buckling rod to `L = 5 µm`, whose three spacings are 1.25, 0.625 and 0.3125 —
none of them 1 — after which M2 kills all three parametrizations. The reason is written into the
class attribute's own comment rather than left in this ledger, because the next person to round `5.0`
down to `4.0` for tidiness will be reading the test and not this file.

M2 is the counterpart lesson to M1 in the other direction: dropping the rest-length normalisation
leaves `F = −grad E` **exactly** true, so no finite-difference check can see it, and every one of its
five kills came from the two independent closed forms. A gradient control and an oracle are not
substitutes for each other.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| node position, hub position | µm | 1e-6 m | finite |
| `axial_modulus_pn` | pN | 1e-12 N | `> 0` (stiffness × length) |
| `bending_rigidity_pn_um2` | pN·µm² | 1e-24 N·m² | `> 0` — a microtubule with zero rigidity is a cable and cannot buckle |
| `minus_end_anchor_stiffness_pn_per_um` | pN/µm | 1e-6 N/m | `>= 0` |
| `minus_end_anchor_rest_um` | µm | 1e-6 m | `> 0` |
| `rest_spacing_um` | µm | 1e-6 m | `> 0` (the length one growth event adds) |
| energy | pN·µm | 1e-18 J | `>= 0` |
| `topology_generation` | — (count) | — | monotone non-decreasing |
| catastrophe / rescue rate | — | — | **not defined here**; refused, §4.6 |

Singular and boundary cases, each with the behaviour Aleph requires. **Refuse, never clamp:**

- **Collapsed segment** (`L = 0`): refused with `ValueError`. The direction is genuinely undefined.
- **A rod base coincident with the hub**: refused. The anchor direction is undefined and a silent
  epsilon returns a finite force whose direction came from the last subtraction's round-off.
- **`rest_length <= 0`** on a segment: refused; the normalisation would divide by zero.
- **A rod with fewer than two active nodes**: refused. A rod with one node has no segment, no
  tangent, and no plus end distinguishable from its minus end.
- **Growth past `capacity`**: refused with `MicrotubuleTopologyError`. The array is not silently
  reallocated, because a reallocation is exactly the quiet resize `topology_generation` exists to
  make impossible.
- **Shrinkage below two active nodes**: refused with `MicrotubuleTopologyError`.
- **Growth while not `GROWING`, or shrinkage while not `SHRINKING`**: refused with
  `MicrotubuleTopologyError`. This is what makes `plus_end_phase` read rather than inert.
- **A material coordinate outside the rod's *active* arclength**: refused, not clamped. Clamping an
  attachment to the tip moves a load path without telling anybody, and on a shrinking rod it would do
  so every step.
- **A `PLUS_END` site requested for `mt_sf_spectraplakin`**: refused with `MicrotubuleRoleError`. The
  registry's role for that edge is a *lattice site*.
- **A length change of more than one node**: refused. A rod that jumps several lattice sites in one
  accepted step has a topology nothing observed.
- **A plus-end site carrying a material coordinate**, or a lattice/quadrature site carrying none:
  refused at site construction. The plus end is not addressable by coordinate — that is what makes
  it the plus end rather than a lattice point that happens to be last.
- **A quadrature site with no quadrature length**, or a non-quadrature site carrying one: refused.
  Without the length a drag integral over the rod is resolution-dependent, which makes a
  reorientation time scale a statement about the node count.
- **A capture recorded against a site whose role is not `mt_cortex_capture`**, or against a site on
  another rod, or against an unpublished site: refused with `MicrotubuleRoleError`.
- **A rate requested**: refused with `MicrotubuleRateUnavailableError`, §4.6.
- **`PROPOSE_KINETICS`**: refused with `ValueError`. That is the phase a dynamic-instability driver
  would occupy and this component holds no rate to drive one with; accepting the phase and returning
  an empty receipt would be indistinguishable, in a coverage table, from a participant that ran.
- **A handle used after a committed length change**: refused with `MicrotubuleStaleHandleError`. A
  *phase* change deliberately does **not** stale a handle: it moves no node, and bumping the counter
  for it would invalidate every outstanding handle for no reason.
- **`rollback()` with no snapshot**: refused. Restoring nothing is not the same as restoring.
- **Empty element set**: returns exactly `(0.0, zeros_like(positions))`, not a crash.

Invariants, each with the test that asserts it (all in `tests/vertical/test_microtubule_controls.py`
unless stated):

- **I1.** `F = −grad E` per term, order 2 —
  `::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order`
- **I2.** Force and moment closure per term against the **constituent** scale —
  `::TestForceAndMomentClosureOnTheIsolatedNetwork::test_every_term_closes_on_its_own`,
  `::test_every_term_carries_no_net_moment`
- **I3.** Translation invariance of the total energy, hub included —
  `::TestForceAndMomentClosureOnTheIsolatedNetwork::test_translating_the_whole_network_and_its_hub_changes_no_energy`
- **I4.** Resolution-independent axial modulus —
  `::TestTheRodIsCompressionCapable::test_end_to_end_stiffness_does_not_depend_on_node_count`
- **I5.** No term mutates its inputs (`np.array_equal`, not `allclose`) —
  `::TestTheTermsArePureFunctions::test_no_term_mutates_its_arguments`
- **I6.** Declared state keys equal the registered contract —
  `::TestOwnershipAndTheStateContract::test_the_declared_keys_are_exactly_the_registered_ones`
- **I7.** No state key of this owner is declared by any other registered owner —
  `tests/state/test_census_key_ownership.py::test_no_state_key_has_two_registered_owners`, and
  `::TestOwnershipAndTheStateContract::test_the_declared_keys_collide_with_no_other_registered_owner`
- **I8.** A committed length change bumps `topology_generation` and staled every outstanding handle —
  `::TestTheDynamicTopologyIsCommittedAndNotQuiet::test_a_committed_length_change_advances_the_generation_and_stales_a_handle`
- **I9.** Rollback restores bit-identically —
  `::TestTheDynamicTopologyIsCommittedAndNotQuiet::test_a_rejected_step_advances_nothing`

## 6. Source evidence class and known retractions

No claim is made here about the reference implementation's microtubule code, because this lane opened
none of it. Everything this entry says about it is quoted from `ALEPH-PORT-1801`, which did read both
revisions and states its own evidence class:

- **Self-declared status, and it is low for exactly the path this module needs.**
  `STATIC_ADAPTER_STATUS = "STATIC_BENDING_ADAPTER_DYNAMIC_TOPOLOGY_PENDING"` is exported from the
  reference's public surface at both revisions. A module that exports `"..._PENDING"` as part of its
  API is telling the reader that path is incomplete.
- **A retraction found by 1801, and it is uncommitted.** A comment block in the working tree retracts
  the committed revision's claim that the dynamic-instability transitions were deferred. A reader of
  `be0e5876` alone gets the retracted claim. It does not affect this module, which inherits neither
  the kinetics nor the rates.
- **Where this lane looked: nowhere inside the reference** — deliberately, so this derivation could
  not be contaminated by it. Whether Aleph's re-derivation agrees or disagrees with
  `ffn_sim/ac/engine/microtubule_rig.py` is **unknown here** and this entry does not guess.

The registered `microtubule` contract carries `citation_status: UNSOURCED`, and nothing here changes
that. `default_microtubule_card()` returns placeholder magnitudes and says so in its docstring.

## 7. Independent oracle or derivation

Four, none of which is the reference implementation:

1. **Central finite differences** of the module's own energies, per term — an oracle for the gradient
   claim that is independent of every constant in the law.
2. **Euler buckling**, `F_c = (4κ/h²) sin²(πh/2L)` discrete and `π²κ/L²` continuum, derived in §4.3
   for this stencil. This is the oracle the brief for this lane names, and it is the one that catches
   what a finite difference cannot: a critical load is an eigenvalue of the whole operator, so a
   rescaling of a single term cannot reproduce it. Mutant M2 is the demonstration — dropping the
   rest-length normalisation leaves `F = −grad E` exactly true and is invisible to every finite
   difference, and every one of its kills came from this oracle.
3. **The end-load closed form** `ΔL = F L0 / k` for a compressed chain, giving resolution
   independence of the axial modulus. Independent of the bending law entirely.
4. **The moment-free identity** of a central-force anchor to a single hub point, §4.4 — an exact
   algebraic statement, checked against the implementation's own reported `mtoc_torque`.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_microtubule_controls.py::TestEachEnergyTermIsTheGradientItClaimsToBe::test_the_force_matches_the_finite_difference_at_second_order` | Per term, `energy > 0` first, then max &#124;F + ∇E&#124; / max&#124;F&#124; at steps 4e-4/2e-4/1e-4 µm with observed order > 1.6 |
| Positive | `tests/vertical/test_microtubule_controls.py::TestTheBendingLawAgainstEulerBuckling::test_the_critical_load_matches_the_discrete_prediction` | Measured buckling load equals `(4κ/h²)sin²(πh/2L)` to 2e-4 relative at N = 5, 9, 17 |
| Positive | `tests/vertical/test_microtubule_controls.py::TestTheRodIsCompressionCapable::test_a_compressed_rod_stores_energy_and_pushes_back` | A rod shortened below its rest length stores strictly positive axial energy and delivers a strictly outward force at its ends |

**Measured this session, on this tree**, with `/Users/sw1/miniconda3/envs/aleph/bin/python`. Every
number in this section was produced by a run in this session; none is quoted from memory.

Finite-difference agreement, steps 4e-4 / 2e-4 / 1e-4 µm:

| Term | Energy (pN·µm) | Force scale (pN) | max&#124;F+∇E&#124; at 1e-4 (pN) | relative | observed order |
|---|---|---|---|---|---|
| axial | 60.530801 | 305.20 | 6.283e-06 | 2.059e-08 | **2.000** |
| bending | 28.651534 | 91.03 | 4.791e-06 | 5.263e-08 | **2.000** |
| minus_end_anchor | 22.152942 | 43.44 | 1.220e-06 | 2.807e-08 | **2.000** |
| total | 111.335277 | 290.79 | 5.276e-06 | 1.814e-08 | **2.000** |

Euler buckling, `L = 5 µm`, `κ = 20 pN·µm²`, axial modulus `5e5 pN`, continuum
`π²κ/L² = 7.895684 pN`:

| N | h (µm) | measured F_c (pN) | discrete prediction (pN) | vs discrete | vs continuum |
|---|---|---|---|---|---|
| 5 | 1.2500 | 7.498277 | 7.498066 | 2.803e-05 | 5.0332 % |
| 9 | 0.6250 | 7.794964 | 7.794736 | 2.921e-05 | 1.2756 % |
| 17 | 0.3125 | 7.870581 | 7.870349 | 2.951e-05 | 0.3179 % |

Observed order of convergence to the continuum value: **1.992**, and every measured load approaches
it **from below**, as §4.3 requires. The residual ~2.9e-05 disagreement with the discrete formula is
stable across all three resolutions — it does not shrink with `h`, so it is not a discretisation
error. It is the compression of the base state that the linearised formula neglects, of order
`F_c/k = 1.6e-05`; the factor of about two between the two is not accounted for here and is named
rather than fitted.

Closure on the isolated network (resultant ÷ constituent scale). The anchor rows include the hub
reaction, which is the other half of that pair and lives outside the node array:

| Term | force | moment |
|---|---|---|
| axial | 0.00e+00 | 9.46e-18 |
| bending | 7.09e-18 | 3.49e-17 |
| minus_end_anchor | 0.00e+00 | 3.62e-17 |
| total | **1.21e-17** | **2.08e-17** |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_microtubule_controls.py::TestTheRodIsCompressionCapable::test_the_same_rod_is_exactly_inert_when_the_axial_law_is_made_tension_only` | With `tension_only_axial=True` a rod compressed to 0.85 of rest stores **exactly** `0.0` and every force is **exactly** `0.0`; with the shipped law the same rod stores > 0 and pushes |
| Negative (must fail) | `tests/vertical/test_microtubule_controls.py::TestTheNegativeControlsFail::test_a_flipped_axial_sign_is_caught_by_the_gradient_test` | With `flip_axial_sign=True` the FD disagreement rises above **0.5 relative**, while the energy is bit-identical to the healthy branch |
| Negative (must fail) | `tests/vertical/test_microtubule_controls.py::TestTheNegativeControlsFail::test_a_quiet_topology_resize_lets_a_stale_handle_scatter` | With `allow_quiet_topology_resize=True` a handle held across a committed length change scatters silently into the resized rod; with the shipped code the same handle raises `MicrotubuleStaleHandleError` |
| Negative (must fail) | `tests/vertical/test_microtubule_controls.py::TestTheBendingLawAgainstEulerBuckling::test_a_rod_with_the_bending_law_omitted_has_no_positive_critical_load` | With `omit_bending_law=True` the transverse Hessian of the compressed chain has no positive generalised eigenvalue, i.e. `F_c = 0`: the rod buckles under any load |

All four breaks are flags on the **shipped** owner and shipped functions, defaulting to `False`, so
each control drives the real code path rather than a hand-edited copy. Each is asserted in **both**
directions on purpose: asserting only that the healthy branch is zero would pass for a function that
returns zero unconditionally, which is the shape of the `PLAN.md` §6.1 defect where a cold cache read
made a deliberately broken strut look innocent.

## 10. Numerical and precision envelope

Working and accumulation precision: float64 throughout; no reduced-precision path exists.

Finite-difference steps 4e-4 / 2e-4 / 1e-4 µm on coordinates of order 1–6 µm. The floor for a central
difference is `eps^(1/3) ≈ 6e-6` **relative**; the smallest step used is 17× above it. Below that
floor cancellation in the energy evaluation dominates, the apparent order goes negative, and a
correct gradient is indistinguishable from a wrong one.

Tolerances, and why they are not looser:

- Per-term relative FD agreement asserted at `1e-7` against a measured worst case of `5.3e-08`
  (bending) — only a factor of **2** of headroom, which is the tightest tolerance here and is stated
  as such. It is not flaky: the geometry is fixed, there is no RNG, and the FD errors are
  bit-reproducible run to run. It is tight enough that mutant M1 cannot pass, and loosening it to
  `1e-6` would still catch M1 — the tightness buys margin against a *smaller* future gradient error,
  not against this one.
- Closure asserted at `1e-13` of the **constituent** scale (`Σ|f_i|`, `Σ|r_i||f_i|`) against a
  measured worst case of `3.6e-17`. Never against the resultant: a resultant-based tolerance gets
  *stricter the more correct the physics is*, because a correct Newton pair cancels, and `PLAN.md`
  records that defect rejecting forty thousand consecutive steps while reporting a plausible tension.
- The buckling load asserted at `2e-4` relative against the discrete prediction, which is ~7× the
  systematic `2.9e-05` base-state offset documented in §8. Asserted tighter, it would fail for a
  reason that is understood and is not a defect; the offset is systematic, so a tolerance below it
  would fail at every resolution at once and look like a physics error.
- The convergence-to-continuum control asserted at `4e-3` relative on the finest chain against a
  measured `3.18e-03`. That is 1.26× of headroom and it is the second-tightest number here — it is
  bounded below by the true `O(h²)` discretisation error at `N = 17`, so the only way to widen it is
  to refine the chain, which costs quadratically in the Hessian build.
- Exact zeros (`== 0.0`, the number) for the tension-only branch, the `mtoc_torque` vector, and the
  active-power channel. These are algebraic identities of a bound gate, not small numbers, and
  asserting them with a tolerance would let a nearly-zero leak pass.

Outside the envelope: every singular case in §5 refuses with an exception rather than degrading
silently. There is no input range over which this module returns a plausible number it cannot stand
behind — **with one exception that is stated rather than hidden**: the material card is unsourced
(§6), so every magnitude it produces is dimensionally right and physically unbacked.

## 11. Production-backend residency and transfer

Host, numpy, float64. `accumulate` is the only method that touches `aleph/runtime/`'s `Backend`
protocol, through `ctx.backend`'s `scale` / `add` / `max_abs`. No device residency, no GPU work of
any kind was run by this lane, and no GPU authorization was sought or held.

When this does run on a production backend, the axial and bending terms are `O(S)` and `O(T)`
gather/scatter kernels with no host round-trip per step. The **topology rebuild is different and is
the one thing worth flagging now**: `commit()` rebuilds the segment and triple tables when a length
change lands, which is a host-side ragged operation whose cost is `O(active nodes)` and which cannot
be fused into a step kernel. That is a consequence of the topology being genuinely dynamic and is not
a defect, but a backend that assumes a static incidence table will need `topology_generation` as its
invalidation signal — which is a second reason the counter exists.

## 12. Comments and docstrings to discard

Nothing to discard: no reference prose entered this module, because no reference file was read by
this lane. The absence is asserted rather than asserted-in-prose —
`::TestTheDeclarationsAgreeWithTheCode::test_no_reference_project_vocabulary_survives` checks the
module text carries none of `ffn_cellsim`, `ffn_sim`, `microtubule_rig`, `static_bending_adapter`,
`dcm_contact`, `surface_body`; and `::test_the_module_imports_nothing_from_validation` asserts the
`aleph/** → validation/**` firewall.

Three pieces of vocabulary from the reference are **named in this ledger and forbidden in the
package**, which is the split `tests/ports/test_port_discipline.py` §1 requires:
`STATIC_ADAPTER_STATUS`, `MicrotubuleRodBuildSpec`, `DynamicInstabilityRateCard`. The third is the
sharpest: Aleph has no rate card, deliberately (§4.6), and a symbol by that name appearing here later
would be a rate arriving without anybody choosing it.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** `111 passed in 0.42s` (`tests/vertical/test_microtubule_controls.py`, 2026-07-31, `/Users/sw1/miniconda3/envs/aleph/bin/python`). Three mutants applied, killing 4 / 5 / 2 tests respectively; none survived; source restored and verified byte-identical by SHA-256 (`ec3e8fef…`). Status stays `PROPOSED`: no PI review, the material card is unsourced, and the module is wired into no assembled vertical. |
| Reviewer | Agent-proposed (lane L32). Unratified. No PI review. |
| Rollback | Delete `aleph/vertical/microtubule.py` and `tests/vertical/test_microtubule_controls.py`. Nothing breaks: no module in the tree imports either, `aleph/vertical/__init__.py` does not export it, and no connector reaches it. |

## 14. Honest limits — what this entry does NOT establish

**Wired connector count contributed by this module: zero.** Four connectors name `microtubule` as an
endpoint — `mt_cortex_capture`, `mt_nucleus_linc`, `mt_sf_spectraplakin`, `mt_cytosol_transfer` — and
this module publishes the microtubule-side endpoint of all four and wires none of them. No `Binding`
row is added to `aleph/vertical/wiring.py`; that is the lead session's lane and a row is added in the
same commit that wires the edge, never before. What is established is that a connector *could* reach
this owner; that no force has ever crossed one of these edges is unestablished and is the next lane's
work.

**No role was withheld, and that is a claim worth stating rather than a default.**
`cortex_filaments.py` withholds `surface_porous_transfer` because a filament graph has no shell to
put a quadrature point on, and withholding is the correct move when the representation genuinely
cannot honour a role. Here all four registered roles are honourable by a rod graph: a plus end is the
last active node, a lattice site is a material coordinate, and a quadrature point is a lattice site
carrying a quadrature length. `PUBLISHABLE_ENDPOINT_ROLES` therefore equals the full role set. **This
is not the same as saying all four can be wired**: `mt_cytosol_transfer` is blocked by an undecided
protocol question — a solid-to-fluid pair is a vector force against a scalar Darcy source and
`wiring.py`'s `CONNECTOR_PROTOCOL` cannot express it — recorded in `docs/design/ENDPOINT_CONTRACT.md`
§2c and in `docs/decisions/PROPOSAL-solid-to-fluid-is-two-books-not-one.md`. That blockage is in the
protocol and not in this representation, which is why the honest move here is to publish and say so
rather than to withhold and imply the rod graph is the problem.

**The quadrature endpoint carries two attributes the declared protocol has no room for.**
`quadrature_length_um` and `tangent_unit` are needed for a resolution-independent, anisotropic drag —
the `mt_cytosol_transfer` contract says a microtubule's drag "is strongly anisotropic and sets the
rate at which the filament can reorient at all" — and `CoupledSiteEndpoint` declares neither. This
module publishes them as extra attributes on the sink, which a structural protocol permits and a
connector written strictly against the protocol will not see. **That is a gap in the endpoint
contract, reported rather than worked around**, and it is the same gap the PROPOSAL above is about.

**Dynamic instability is state and machinery, not a model.** No catastrophe rate, no rescue rate, no
capture rate, no load dependence of any of them, and no stochastic driver. Nothing in this module
decides *when* a transition happens; a caller proposes one and the module commits it. That is
deliberate (§4.6) and it means the following are unestablished here: that any transition sequence
this module can represent is a physical one; that the phase machine's state space is the right one;
that `growth_shrinkage_state`'s counters are the quantities an inference procedure will want. The
first inference lane that consumes them may well find they are not.

**The MTOC is a hub and its torque is an honest zero, not a computed one.** `mtoc_torque()` returns
exact zeros because the anchor is a central force to a single point, and §4.4 proves it must. What is
**not** established is anything about a real centrosome: no centriole pair, no pericentriolar
material, no nucleation geometry, no orientation dynamics. The rotational half of `mtoc_pose` is
carried and never driven; a control asserts it does not move rather than pretending it might.
Anything requiring MTOC torque needs an anchor law with a finite moment arm, which is a different
model and not a refinement of this one.

**Physical magnitudes are not established.** Evidence rung `ANALYTIC_ORACLE`, quantitative status
`BLOCKED`. No modulus, rigidity, anchor stiffness or spacing here is a microtubule number; they are
chosen to make the algebra sharp. `citation_status` for the `microtubule` contract remains
`UNSOURCED` and this entry does not improve it. In particular the buckling controls establish that
the *discrete operator* has the critical load its own energy implies — they establish nothing about
whether a real microtubule at that length buckles at that load.

**Never stepped inside a world.** The module is not in `aleph/vertical/assembly.py` and has never
been advanced by the runtime. `accumulate` is exercised phase by phase by a control; a sequence of
accepted and rejected steps against other participants is not. Specifically unestablished: that the
overdamped update in `SOLVE` is stable at any particular `dt`, and that the topology rebuild
interleaves correctly with another owner's commit.

**The reference implementation is unaudited by this lane.** Everything in §2 and §6 is inherited from
`ALEPH-PORT-1801`. Whether this re-derivation agrees or disagrees with the reference's rod rig is
unknown, and this entry deliberately does not guess.

**Not run:** any GPU job (zero, and no authorization sought or held), and the full suite. Only
`tests/vertical/test_microtubule_controls.py`, `tests/state/` and `tests/ports/` were exercised, with
other sessions writing concurrently — `CLAUDE.md` §1 says to report a foreign breakage rather than fix
it, which means not conflating one with mine. Two failures in `tests/ports/` predate this lane and
are reported in the lane summary, not repaired here.

**`ports/ledger/INDEX.md` is not regenerated by this lane.** It is a shared file and adding this entry
makes it stale. That is reported to the lead session rather than fixed here, because regenerating it
would sweep another session's concurrently-added entry into this lane's change.
