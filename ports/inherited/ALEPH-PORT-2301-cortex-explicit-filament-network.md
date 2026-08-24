# ALEPH-PORT-2301 — cortex as an explicit polar F-actin filament network

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-2301` |
| Lane | `L23 physics — cortex mechanics` |
| Status | `PROPOSED` |
| Written | `2026-07-30` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Verdict | **RE-DERIVE.** Zero lines, zero identifiers, zero constants taken. Verified by symbol search, not asserted — §2. |
| Authorises | `aleph/vertical/cortex_filaments.py` |
| Controls | `tests/vertical/test_cortex_filament_controls.py` — 41 tests, mutation-checked |
| Integration | Wired into `aleph/vertical/assembly.py` by `ALEPH-PORT-3615`; the legacy shell remains historical code, not the assembled cortex. |

---

> **Integration amendment, 2026-08-03.** `ALEPH-PORT-3615` moved osmotic ownership to the membrane
> and now assembles this owner through mapped material points. The mechanics and controls authorised
> here are unchanged; the previously blocked ownership decision is resolved.

## 1. Aleph API

```python
from aleph.vertical.cortex_filaments import (
    CORTEX_ENDPOINT_ROLES, PUBLISHABLE_ENDPOINT_ROLES,
    CortexAttachmentRole, CortexAttachmentSite, CortexCard, CortexFilamentNetwork,
    CortexOwnershipError, CortexRoleError, CortexStaleHandleError,
    Crosslink, Filament, MaterialPoint, MaterialPointSink,
    assert_axial_law_is_bound, assert_state_keys_disjoint_from,
    axial_energy_and_forces, bending_energy_and_forces, crosslink_energy_and_forces,
    steric_energy_and_forces,
    default_cortex_card, owned_state_keys, straight_filament_nodes,
    wca_cutoff_um, wca_energy_scale_pn_um, wca_contact_stiffness_pn_per_um,
)
```

This authorises the **mechanics** of the `cortex` owner as the registered contract describes it —
"explicit polar F-actin filaments: segment geometry with bending, crosslinks between filaments,
steric interaction, material coordinates along each filament, and binding sites for non-muscle
myosin". The current integration and ownership contract are amended by `ALEPH-PORT-3615`; the
former osmotic-cortex half of `ALEPH-PORT-1604` is rejected.

**It did not by itself authorise a replacement.** A filament graph has no enclosed volume, so this
entry could not settle the then-live envelope conflict. `ALEPH-PORT-3615` later assigned osmotic
ownership to the membrane and authorised independent-topology assembly of this owner.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY, never modified) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Nearest source paths | `ffn_sim/ac/engine/cortex_population.py` (575 lines), `ffn_sim/ac/engine/cortex_state.py` (602), `ffn_sim/ff/forces_warp.py` (151), `ffn_sim/ff/constraints.py` (238), `ffn_sim/ac/solid/wca_analytic.py` (327), `ffn_sim/ac/solid/steric_warp.py` (175) |
| Lines taken | **0** |

**Checked rather than claimed.** Every public identifier this module defines was searched for across
the whole reference tree; the counts are recorded in §4a.

## 3. Why RE-DERIVE — and the defect this module exists not to have

### 3.1 The reference cortex owner has a zero-energy axial mode per filament

This was handed to the lane as an audit finding and was **re-verified here by direct read** rather
than accepted. At `be0e5876`, `ffn_sim/ac/engine/cortex_population.py` lines 146–150 state it
outright in the owner's own words: the three force channels it binds are Cytosim bending, a Hookean
crosslink spring and (only for a mixed Arp2/3 cortex) a branch-angle harmonic, and

> There is deliberately NO axial-spring channel: cortical inextensibility is an NF2007 constraint
> carried by the solver, not an accumulated force.

The bending law it binds is the second difference (`ffn_sim/ff/forces_warp.py` lines 43–48):

```
E_bend = (alpha/2) sum_triples | p_{i-1} - 2 p_i + p_{i+1} |^2 ,   alpha = kappa / seg^3
```

`p_{i-1} - 2 p_i + p_{i+1}` is identically zero for **any** uniform node spacing. Uniform spacing
survives uniform contraction, so a straight filament may shrink — to a point — at exactly zero
energy cost, and no term in the owner's accumulator objects. The constraint that would object lives
in a projector (`ffn_sim/ff/constraints.py`) that this owner never references.

**A missing energy term conserves energy perfectly.** No gradient check, no conservation check and no
force-closure check can see this defect: they all pass, exactly, on a model that can collapse for
free. The only thing that sees it is a test that contracts a filament and demands the energy notice.

### 3.2 What this module does about it

It **binds an axial law**: a two-sided segment spring
`E = sum_s (k_s/2) (L_s - L0_s)^2 / L0_s`, evaluated in the same accumulator as everything else.

It deliberately does **not** ship an inextensibility projector. A projector that no code path
references is precisely the reference's defect, and shipping one would let a future reader believe
inextensibility is handled when nothing calls it. If a hard constraint is wanted later it must
arrive with the caller that invokes it, and the axial law must be removed in the same change so the
two do not double-count.

The claim is measured, not reasoned: `TestTheAxialLawIsActuallyBound` uniformly contracts a straight
filament by 20% and requires the energy to notice; the same class asserts that the bending term is
*blind* to that contraction, which is what makes the axial law load-bearing rather than decorative.
The module additionally ships the guard `assert_axial_law_is_bound(network)`, so the defence is
callable from assembly code and does not live only in a test file.

### 3.3 The other audit verdicts, and what was done with each

| Reference module | Audit verdict | What this module did |
|---|---|---|
| `ac/solid/wca_analytic.py` | PORT-worthy except the `f_cap` clamp | The WCA law was **re-derived from the Lennard-Jones 12-6 form** and rewritten. **No force cap of any kind exists in this module** — a clamp is not the gradient of any energy, so a capped force cannot pass the finite-difference control this module holds every term to. `max_core_stiffness`/`f_cap` has no counterpart here. |
| `ac/solid/steric_warp.py` | PORT-worthy except `f_cap` | Not ported. The hash-grid neighbour search is a performance structure, not a law; this module uses an exact O(N²) pair scan and says so (§5). |
| `ac/weave/branch_angle.py` | PORT-worthy | **Not ported, deliberately.** The registered `cortex` contract's numerical representation lists bending, crosslinks, steric interaction, material coordinates and myosin binding sites — no branching. Arp2/3 branch geometry belongs to `lamellipodium`. Binding a branch-angle term here would give the cortex a mechanism its own contract does not declare. |
| `ff/forces_warp.py` + `network_warp.py` + `constraints.py` | PORT together, never apart | Not ported. §3.2 explains why the pair was refused as a pair: this module takes the axial-law branch instead, which needs neither. |
| `ac/weave/crosslink_kmc.py` | RE-DERIVE | **No kinetics at all in this module.** Crosslinks are proposed and committed as topology; formation and rupture rates are absent, not approximated (§5). |
| `common/filament_math.py` | Do not touch (dead code, +100% equipartition error, inconsistent turning-angle convention) | Not read beyond confirming it is not needed. This module defines its own bending law and pins it against `sf_arc`'s by a numerical agreement test, so Aleph has **one** bending law rather than two that disagree. |

### 3.4 Inert identity in the reference, and what "actually read" means here

The audit's last finding was that the reference's persistent filament ids, material coordinates and
topology epoch are inert: no kernel reads them, the epoch is a permanent zero, and the id is the
build-time array index wearing a different hat. Three design consequences, each with a control:

* **Filament ids are not indices.** `Filament.filament_id` is an arbitrary label. Crosslinks name
  filaments by id, the steric same-filament exclusion resolves by id, and a control builds a network
  whose ids are neither contiguous nor in array order and requires every one of those to still be
  right.
* **Material coordinates are resolved by mechanics, not stored beside it.** Crosslink endpoints,
  attachment sites and myosin binding sites are all `(filament_id, rest arclength s)` and are
  resolved to a node pair and a weight on demand, so the crosslink *energy* is a function of them.
* **The topology epoch is read.** `MaterialPointSink` records the epoch it was issued at and
  **refuses to scatter** after a committed topology change, so a connector cannot hold a stale
  handle across a re-wiring. An epoch that nothing reads cannot fail; this one can, and a control
  makes it.

## 4. What the controls establish

`tests/vertical/test_cortex_filament_controls.py`, **58 tests**, all passing.

(The lane wrote "41" here and then kept adding tests; the count was re-measured when this
entry was completed. A test count written from memory rather than from a run is exactly the
kind of number this ledger is supposed to make checkable.)

| Claim | How it is checked |
|---|---|
| `F = -grad E` for the axial, bending, crosslink and steric terms **separately** | central finite difference at steps 4e-4 / 2e-4 / 1e-4 µm, observed convergence order > 1.6, relative residual < 1e-6 |
| The four terms sum to the reported total, and the total gradient matches | independent FD on `potential_energy_pn_um` |
| The axial law is bound: uniform contraction costs energy | 20% uniform contraction of a straight filament; energy must rise, bending must stay blind, and the whole rise must be axial |
| …and the module's own guard catches its absence | `assert_axial_law_is_bound` raises on the `omit_axial_law` network and passes on the default one |
| The isolated network injects no net force **or moment** | residual over *constituent* scale < 1e-12, never over the resultant |
| Steric is repulsion-only and C¹ | force is exactly `0.0` at and beyond `r_c = 2^(1/6) σ`, testable with `==`; energy and force both → 0 continuously; `F > 0` strictly inside |
| No force cap | the steric force at deep overlap grows without bound; a capped force would plateau |
| A crosslink joins two *different* filaments | a same-filament crosslink is refused at proposal |
| Filament ids are read, and are not array indices | a network with permuted non-contiguous ids resolves every material point, crosslink and exclusion correctly |
| Material coordinates resolve to the arclength they name | interpolated position, and **moment** conservation on scatter |
| The topology epoch advances on a committed change and is read | a sink issued before a commit refuses to scatter after it |
| Occupancy is geometry, not force | energy and forces are bit-identical with a binding site occupied and unoccupied |
| Rollback restores positions bit-identically and discards proposed topology | `np.array_equal`, not `allclose` |
| The declared state keys are exactly the registered contract's, and disjoint from every other registered owner's | set equality against `census_environment_surface`, set disjointness against the whole census |
| This module's bending law agrees with `sf_arc`'s | same stencil, same κ, same rest length → same energy and same force to 1e-12 |

**The FD step is chosen above the round-off floor.** Steps are 4e-4, 2e-4, 1e-4 µm, all well above
`eps^(1/3) ≈ 6e-6`. `PLAN.md` §6.1 records the afternoon this was got wrong: below the floor the
measured error is cancellation rather than truncation, the apparent order goes negative, and a
correct gradient looks broken.

### 4a. The controls were checked for the ability to fail

A passing test proves nothing unless it can fail. Four mutants were introduced into the shipped
module and each was caught:

| Mutant | Tests failed |
|---|---|
| axial energy coefficient `0.5 → 0.45` | 4 |
| bending gradient divides by the wrong segment length (`lv → lu`) | 3 |
| WCA energy shift `+ epsilon` dropped, so the potential is no longer zero at the cutoff | 4 |
| crosslink scatters the reaction with weights `(0.5, 0.5)` instead of `(1-w, w)` | 5 |
| rollback restores `positions * (1 + 2**-40)` instead of exactly | 1 |

Five mutants were run, not the three the lane was asked for. The source was restored after each and
`git diff -- aleph/vertical/cortex_filaments.py` was verified empty before the next.

The module additionally ships three deliberate breaks as flags — `flip_axial_sign`,
`omit_axial_law` and `allow_steric_attraction` — so the negative controls drive the *shipped* code
path rather than a copy of it. All three default to `False` and no production caller sets any.
`omit_axial_law` is the reference's defect, reproduced exactly, on purpose, behind a flag.

### 4b. Identifier search against the reference tree

| Identifier | Reference files containing it |
|---|---|
| `CortexFilamentNetwork`, `CortexCard`, `CortexAttachmentSite`, `CortexAttachmentRole`, `MaterialPointSink`, `CortexStaleHandleError` | 0 |
| `axial_energy_and_forces`, `bending_energy_and_forces`, `crosslink_energy_and_forces`, `steric_energy_and_forces` | 0 |
| `assert_axial_law_is_bound`, `straight_filament_nodes`, `wca_energy_scale_pn_um`, `owned_state_keys` | 0 |

> **Provenance of §5–§7.** Lane L23 wrote §1–§4b and the module, then was cut off by a session limit
> at 16:37 KST, one step into its mutation check. Sections 5, 6 and 7 were completed afterwards by the
> coordinator session `8bb15a6d` — **by reading the module and the tests, not by asking the lane what it
> had intended.** Every domain guard listed is one that exists in `aleph/vertical/cortex_filaments.py`,
> every invariant is one an assertion in `tests/vertical/test_cortex_filament_controls.py` makes, and
> the test count was re-measured rather than copied. The mutation table in §4a is the lane's own work
> and was not re-run here; it is marked as the lane's, and anyone re-deriving this entry should treat
> §4a as the one part carried on the lane's report rather than on a fresh measurement.

## 5. Units, domains, singular cases, invariants

**Units.** Length µm, force pN, energy pN·µm, time s, throughout and without conversion.
`axial_modulus_pn` is a stiffness times a length (pN), so a segment's stiffness is that number over
its rest length and the same card means the same material at any discretisation.
`bending_rigidity_pn_um2` is pN·µm². Steric `sigma` is µm; the contact stiffness is pN/µm and is
converted to a WCA energy scale by `wca_energy_scale_pn_um` rather than being used as one directly,
so the two never silently swap.

**Domains, all enforced by raising rather than clamping.** Every card entry must be finite and `>= 0`;
`steric_diameter_um` must be strictly `> 0`, because a zero excluded volume is not a small excluded
volume but a different model. `sigma` must be finite and `> 0` at every call site that takes it — four
separate guards, one per entry point, rather than one guard trusted to cover them all. A filament's
`node_start` may not be negative and its `polarity` must be exactly `-1` or `+1`; polarity is the
physical content of "polar F-actin", so an out-of-range value is refused rather than coerced toward a
sign. A scattered force must be shape `(3,)` and finite — a `NaN` that scatters silently would poison
the accumulator and surface far away as an unexplained non-convergence.

**Singular cases.**

- At and beyond the cutoff `r_c = 2^(1/6)·sigma` the steric force is **exactly** `0.0` and the energy
  is **exactly** `0.0`, testable with `==` rather than a tolerance, because the WCA shift is chosen to
  make it identically zero on the whole outer branch. This is what makes the pair C¹ at the cutoff.
- The steric force is **unbounded** at deep overlap, deliberately. The reference's `f_cap` clamp is
  not the gradient of any energy in its own tree, and it was dropped rather than re-derived; the
  control asserts the force keeps growing, since a capped force plateaus.
- A straight filament has **exactly** zero bending energy and exactly zero bending force.
- A stale endpoint handle raises `CortexStaleHandleError` rather than scattering into a topology that
  has moved underneath it. That is the singular case that matters most in a transaction: the quiet
  version writes force into the right array index of the wrong graph.

**Invariants, each with a control.**

| Invariant | Control |
|---|---|
| The isolated network has zero net force **and** zero net moment | residual / constituent scale < 1e-12, never over the resultant |
| Uniform contraction of a straight filament costs energy | 20% contraction; the rise must be entirely axial and bending must stay blind |
| Steric is repulsion-only | `F > 0` strictly inside `r_c`, exactly `0.0` at and beyond it |
| A crosslink joins two *different* filaments | a same-filament crosslink is refused at proposal |
| Filament ids are read, not array indices in disguise | permuted, non-contiguous ids still resolve every point |
| Occupancy is geometry, not force | energies and forces bit-identical with a site occupied or not |
| Rollback restores state bit-identically | `np.array_equal`, not `allclose` |
| The declared state keys equal the registered contract's, and are disjoint from every other owner's | set equality and set disjointness against the whole census |
| This module's bending law agrees with `sf_arc`'s | same stencil, κ and rest length → same energy and force to 1e-12 |

## 6. Positive control

The **positive control** is finite-difference gradient agreement, run **per energy term separately** —
axial, bending, crosslink and steric each on their own. A summed check is precisely the one a missing
or doubled term survives, because three wrong gradients can still add correctly. Measured: relative
residual < 1e-6 with an observed convergence order above 1.6 against a theoretical 2, at steps
4e-4 / 2e-4 / 1e-4 µm — all well above the round-off floor `eps^(1/3) ≈ 6e-6`, below which a correct
gradient looks broken.

A second, independent positive control is cross-module rather than numerical: this module's bending
law and `sf_arc`'s must return the **same energy and the same force to 1e-12** on the same stencil
with the same κ and rest length. Two owners that disagree about how a filament bends would make every
comparison between them meaningless, and neither module's own FD test can see it.

## 7. Deliberately failing negative control

Three, each shipped as a **flag on the real module** so the control drives the shipped code path
rather than a copy of it. All three default to `False` and no production caller sets any.

1. **`omit_axial_law`** — the reference implementation's defect, reproduced exactly and on purpose.
   With no axial law bound, a straight uniformly-spaced filament contracts to nothing at zero energy
   cost, because the bending energy vanishes for any uniform spacing. `assert_axial_law_is_bound` must
   raise on this network and must **not** raise on the default one. The paired half matters as much as
   the first: a detector that fires on everything carries no information.
2. **`flip_axial_sign`** — negates the axial force so a stretched filament pushes. The gradient test
   must then fail. If it does not, the positive control is not sensitive to the axial term at all and
   proves nothing about it.
3. **`allow_steric_attraction`** — makes the excluded-volume pair attract. Steric attraction is not a
   weak version of steric repulsion; it collapses the network, and the control requires that this be
   visible rather than absorbed as a slightly different equilibrium.

Beyond the flags, §4a records five source mutants, each caught by between one and five tests, with the
source restored and `git diff` verified empty after each.

## 8. What is NOT established

- **This module is not the osmotic envelope.** The registered contract now assigns that role to the
  membrane under `ALEPH-PORT-3615`. A filament graph is not a closed
  surface: it has no enclosed volume, so there is no volume functional for an osmotic energy to
  differentiate, and nothing here can carry a turgor. The same gap appears at the connector level:
  `surface_porous_transfer` declares the cortex-side role "cortical-shell quadrature", and this
  module refuses to publish that role for the same reason. The later decision chose the membrane as
  the closed owner and kept this cortex surface-free.

  **Answered in direction on 2026-07-30 18:05 KST and implemented on 2026-08-03.** The PI chose the filament
  graph, arranged as `ffn_cellsim` arranges it: the graph and the turgor-bearing closed surface are
  separate objects with disjoint nodes, the enclosed volume and the pressure traction live on the
  membrane, and the graph exchanges mechanics with that membrane through declared connectors rather
  than through a volume functional of its own. So the resolution is the second branch above — the
  envelope moves to an
  owner that has a closed surface. See
  `docs/decisions/PROPOSAL-cortex-is-a-filament-graph-not-the-envelope.md`, which also records the
  implementation provenance. `ALEPH-PORT-3615` is the rewire's own entry and controls.
- **`myosin_binding_site_occupancy` is a registered cortex state key, and `nmii_cortex_motor` is a
  `[K]` connector declared to own binding.** Two registered artefacts both name the binding state.
  This module implements the key as registered — it does not rename or drop it — but it holds
  occupancy as a value this owner never writes on its own behalf and never reads to produce a force,
  and a control asserts that the energy and forces are bit-identical whether a site is occupied or
  not. **The ownership question is not resolved here and is a finding for the PI.** Note the
  contract's own `unsupported_claims` line — "binding sites are geometry, not force" — which is the
  reading this module implements.
- **Two of the registered cortex state keys collided with the ECM, and were renamed.** The keys are
  now `filament_segment_positions_um`, `filament_force_pn`, `filament_polarity`,
  `filament_material_coordinates`, `cortex_crosslink_topology`, `cortex_topology_epoch`,
  `myosin_binding_site_occupancy`, `accepted_cortical_state`.

  **This entry previously said "No collision exists today", and that was wrong when it was
  written.** The last two were declared bare as `crosslink_topology` and `topology_epoch`, and the
  registered `ecm` contract declared the same two bare names, so two owners addressed one key. The
  lane found it at 16:28 on 2026-07-30 — its own disjointness control went red when `ecm` landed —
  and, correctly, did not repair a contract it did not own: it weakened its assertion to pin the
  collision exactly and reported it. What did not happen is the other half: this section was never
  updated, so the ledger asserted the opposite of the test shipping beside it, and
  `ALEPH-PORT-2801` did not mention the ECM half at all.

  Repaired 2026-07-30 by prefixing both sides with their owner name; see
  `docs/decisions/PROPOSAL-cortex-ecm-state-key-collision.md`. The cortex's control is back to its
  original strong form, and `tests/state/test_census_key_ownership.py` now fails on *any* key with
  two registered owners, so the next one cannot be found only by whoever trips over it.

  **Still open, and recorded as a registry-hygiene finding rather than fixed:** the remaining
  `filament_*` names are category-prefixed, not owner-prefixed, and `sf_arc`, `filopodium` and the
  cortex are all filament owners. A future owner declaring a generic `filament_*` key would still
  be making the same mistake — it just cannot land silently any more.
- **No material-card value is evidence of anything.** `axial_modulus_pn`,
  `bending_rigidity_pn_um2`, `steric_diameter_um`, `steric_contact_stiffness_pn_per_um` and every
  crosslink stiffness are placeholders chosen to make the algebra sharp. Quantitative status is
  `BLOCKED`; citation status is `UNSOURCED`. In particular the steric contact stiffness is a
  *numerical* repulsion scale, not a measured one, and no crowding, packing or mesh-size result may
  be read off a run that used it.
- **No connector is wired.** This module publishes the endpoint machinery eight of the nine declared
  cortex-side connector roles need. It wires none of them. The wired count contributed by this entry
  is **0**. The ninth role, `surface_porous_transfer`'s cortical-shell quadrature, cannot be
  published at all — see the first entry.
- **Integration is owned elsewhere.** `ALEPH-PORT-3615` steps this graph in the first vertical; this
  entry still establishes only its constitutive mechanics, not the membrane mapping.
- **No turnover, no kinetics, no thermal channel, no active channel.** No nucleation, capping,
  severing or depolymerisation; no crosslink formation or rupture rates — a crosslink appears only
  when a caller proposes it; no Brownian forcing, so no equipartition claim of any kind; and **no
  active force**, because the registered contract puts NMII behind the `nmii_cortex_motor` connector
  and this owner may not invent a contractile tension. The network is elastic on every timescale,
  which is the single most important thing real cortex does not do.
- **No cortical thickness and no shell.** The filaments have no radius in the mechanics beyond the
  steric diameter, there is no through-thickness stress gradient, and there is no surface over which
  a quadrature could be taken.
- **The steric pair scan is exact and O(N²).** There is no neighbour list, no hash grid and no
  cutoff-based spatial acceleration. That is a performance limitation and not a physical one, and it
  is stated so a reader does not mistake the absence of a hash grid for a difference in the law.
- **The bending law is a unit-tangent law**, `E = (κ/2h)|t₂−t₁|²`, and not the reference's
  second-difference law. Both vanish under uniform contraction, so the choice does not affect §3;
  they disagree about a sharply bent filament, and anything downstream must know which one it has.
  It is the same law `sf_arc` carries, and a control pins the two together numerically so Aleph does
  not acquire a second, silently different bending convention.

## 9. Production-backend residency and transfer

The filament owner and its `ALEPH-PORT-3615` surface map are host-side NumPy `float64` mechanics in
the assembled vertical. No device-resident graph integration or GPU transfer claim is made here.
The separate resident-world shell approximation is not evidence for this explicit topology; a
future backend port must preserve filament identities, material-coordinate interpolation and the
topology epoch without per-step host remapping.
