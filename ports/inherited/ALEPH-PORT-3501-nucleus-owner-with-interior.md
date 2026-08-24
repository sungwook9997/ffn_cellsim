# ALEPH-PORT-3501 — the nucleus gets an interior: a volume constraint and a chromatin network

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3501` |
| Lane | `46143f30` nucleus lane |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `SOURCE-DERIVED` |
| Exists because | The registered `nucleus` contract declares `chromatin_network_state` and `enclosed_volume`, and no implementation carries either. The nucleus as built resists areal stretch and not compression at constant area. |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list, and no other file, is
covered by it.

```python
from aleph.vertical.nucleus_interior import (
    ChromatinCard,
    NucleoplasmCard,
    chromatin_energy_and_forces,
    enclosed_volume_um3,          # RE-EXPORTED from aleph.vertical.membrane, not defined here
    reference_vertex_areas_um2,
    volume_energy_and_forces,
)
from aleph.vertical.nucleus import (
    NucleusOwner,
    NuclearEnvelopeBoundary,
    build_spherical_nucleus,
    owned_state_keys,
)
```

`aleph/vertical/nucleus_lamina.py` is modified in exactly one place — an appended
`owned_state_keys()` declaring the registered contract's nine keys — and its two energy laws
(`lamina_areal_energy_and_forces`, `curvature_energy_and_forces`, `ALEPH-PORT-3001`) are **imported
and consumed unchanged**. This entry re-derives none of them.

`aleph/vertical/membrane.py` and `aleph/vertical/cortex.py` each gain an `owned_state_keys()` too.
That is a declaration read off the registry, not physics, and it is covered here only because it
lands in the same commit series.

**The volume law is NOT authorised by this entry, because this entry does not implement it.**
`aleph/vertical/membrane.py` has shipped `enclosed_volume_um3` and its exact gradient
`volume_gradient` since `ALEPH-PORT-1101`, and `nucleus_interior` imports both.
`volume_energy_and_forces` contributes only the constitutive part — the strain measure, the energy,
and the chain rule `−dE/dV · ∂V/∂x`. The first version of this module carried its own copy; §3.1
records why that was wrong and what it measured.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/cell/compartments.py` (the nucleoplasm volume kernels), `ffn_sim/ac/nucleus/chromatin_analytic.py`, `ffn_sim/ac/nucleus/lamina_analytic.py` |
| Source symbol(s) | `_nucleoplasm_pressure_kernel`, `nucleoplasm_volume_force_arr_kernel`, `wlc_tension`, `wlc_energy`, `wlc_small_strain_stiffness` |
| Read from | **working tree**, 2026-07-31 |
| Working tree == commit? | **yes**, for all three paths — `git status --porcelain` on those paths returned empty at `be0e5876`. |

The `Read from`/`Working tree` rows are the ones `ports/TEMPLATE.md` says are not bureaucracy, and
they are answered by measurement rather than by assumption: at this same commit the source tree has
carried uncommitted changes before, so "clean" is a thing to check per path and not per repository.

Aleph inputs that were also read, and which are **the larger part of this module's provenance** —
larger than §3's first draft admitted:

- the registered `nucleus`, `chromatin`, `nucleoplasm` and `nuclear_envelope` contracts in
  `aleph/state/full_census.py`;
- `aleph/vertical/nucleus_lamina.py` (`ALEPH-PORT-3001`) for the surface laws, the `LaminaCard`
  signature and the reference-areas convention;
- **`aleph/vertical/membrane.py` (`ALEPH-PORT-1101`) — `enclosed_volume_um3` and `volume_gradient`,
  which this module now imports rather than re-implements**, and `HelfrichMembrane` as the canonical
  owner shape;
- **`aleph/vertical/pressure.py` (`ALEPH-PORT-1102`) — `TurgorEnsemble.VOLUME_PENALTY`, which is the
  same quadratic volume functional on a different compartment, and `assert_single_pressure_path` /
  `osmotic_sources_among` / `PressureDoubleCountError`, which are the audit this owner must be
  visible to.** See §3.1 and §5's I9.
- `aleph/state/manifold.py` for `icosphere_arrays`, `TriangulatedSurface` and the outward-winding
  guarantee; `docs/design/ENDPOINT_CONTRACT.md` I1–I3.

## 3. Why source-derived porting beats clean-room

**Partly it does and partly it does not, and the split is per-term.**

**The volume constraint: NOT source-derived after all, and the first draft of this section was
wrong.** It argued that the vertex gradient `∂V/∂x_a = (b × c)/6`, with its three cyclic
assignments, is the kind of easily-mis-rotated case analysis §3 calls a legitimate reason to port —
and that the port consisted of re-deriving it and finding that `ffn_cellsim`'s kernel agrees.

**The identical construction was already in this repository, controlled, and one import away.**
`aleph/vertical/membrane.py:356` `volume_gradient` is that exact expression, and
`aleph/vertical/pressure.py:97` already implements `E = (K_V/2)(V/V₀−1)²` over it as
`TurgorEnsemble.VOLUME_PENALTY`. Resting a `SOURCE-DERIVED` justification on agreement with an
external kernel, when the same construction was in-repo with its own controls and its own ledger
entry, is a provenance error: it credits the wrong source and it hides a duplication. §3.1 records
what that duplication measured before it was removed.

What survives as genuinely source-derived is narrow and is stated as such: the source's kernel was
read, and its `p_vol = −K(V−V₀)/V₀` with `f = p_vol ∂V/∂x` **corroborates** the sign and factor
convention arrived at independently in §4.2. Corroboration is worth recording and it is not a port.
Were this entry re-classified today the honest class would be `RE-DERIVED` for the volume half; it
is left `SOURCE-DERIVED` because the chromatin half took one structural idea from
`ac/nucleus/chromatin_analytic.py` and because a class is not an agent's to revise after review —
that is flagged in §14 rather than changed here.

### 3.1 What the duplicate measured, before it was removed

There were **three** `enclosed_volume_um3` in the package. The membrane's uses
`v₀ · (v₁ × v₂)` summed per face; this module's first version used `(a × b) · c` with a fused
`einsum`. Both are correct — the scalar triple product is cyclic — and they disagreed:

| Comparison | Measured |
|---|---|
| membrane vs nucleus form, level-2 sphere at rest | `9.9e-14` µm³ |
| worst over 200 randomly deformed meshes | `2.0e-13` µm³ |
| this module's scatter vs `−dE/dV · membrane.volume_gradient` | `7.1e-15` pN |

Meanwhile `membrane.enclosed_volume_um3`'s own docstring promised that the volume it differentiates
and the volume the manifold validates are *"one number, not two."* The duplicate made that false,
and it would have stayed false silently, because nothing compares two correct implementations.

Removed by importing both functions. `enclosed_volume_um3` is re-exported so the API in §1 is
unchanged, and `test_the_volume_law_has_exactly_one_implementation_in_the_package` asserts the
**identity of the function object** rather than the equality of its output — two implementations
that agree today still drift; one implementation cannot.

### 3.2 The two `K_V` are not interchangeable, and a reader who spots the duplication will get it wrong

```
this entry     E = (K_V / 2) · V_ref · (V/V_ref − 1)²      K_V in pN/µm²
pressure.py    E = (K_V / 2) ·         (V/V₀   − 1)²      K_V in pN·µm
```

so `K_pressure = K_nucleus · V_ref`, a factor of ~1.1e2 µm³ on a 3 µm nucleus. Both spellings are
defensible — one is a modulus, the other is already the compartment's integrated stiffness — and
now that the two live one import apart, the trap is copying a number across without the factor.
Stated here, in `volume_energy_and_forces`'s docstring, and exercised by
`TestTheVolumeTermIsCountableByTheSinglePathAudit`, whose double-count control performs exactly this
conversion in order to make the two loads bit-identical.

**The chromatin term: NOT ported, and this entry says so plainly.** The source's chromatin is a
Marko–Siggia worm-like-chain population (`wlc_tension`/`wlc_energy`) whose small-strain limit
`k = 3 k_BT / (2 ℓ_p L)` feeds a three-regime areal law on the envelope. This module implements a
**quadratic radial network** instead. It takes from the source exactly one structural idea — that
the interior's stiffness is carried by strands anchored between the nuclear interior and the
envelope's own vertices — and nothing else: no WLC form, no persistence length, no `k_BT`
convention, no magnitude. Under §3's own test that is a reason to write clean-room, and it was.

**The one design divergence, stated because a future reader would otherwise assume it was an
oversight.** In the source, chromatin lives *inside* the lamina's areal law as the small-strain
branch of a three-regime split (`ac/nucleus/lamina_analytic.py`). Here it is a **separate radial
term**. Both are parallel wirings on the envelope's own vertices, so neither achieves a true series
compliance — theirs hides that inside one function and this one states it. Separating it is what
lets the control ask *"did the lamina stay blind to this?"*, a question that cannot be asked of a
single fused function. That is a preference for testability, not a claim that their encoding is
wrong.

**`PLAN.md` §0.2 applied in full.** Every original comment and docstring is discarded (§12), the law
is re-derived below, and `PLAN.md` §7's prior was taken seriously: of ten candidates an earlier
audit rated most liftable, one survived being run. The source was read to check a derivation, not to
supply one.

## 4. Physical or mathematical law represented

### 4.1 The enclosed volume, and why it references no origin

For a closed, oriented, outward-wound triangulated surface, the divergence theorem with the field
`F = x/3` (`div F = 1`) gives

```
V = ∫_Ω dV = (1/3) ∮_∂Ω x · n dA = (1/6) Σ_f (a_f × b_f) · c_f
```

the last step because a triangle with vertices `(a, b, c)` has `n dA = (b − a) × (c − a) / 2` and
`x · n dA` integrates over the face to `(a × b) · c / 2`, times the `1/3`.

**Translation invariance is a theorem here, not a hope.** Replacing `x → x + t` changes the surface
integrand by `t · n dA`, and `∮ n dA = 0` for any closed surface, so `V` is unchanged. Written in
the triple-product form this is the statement that the sum of the face area-vectors vanishes. It is
asserted by a control rather than assumed, because a volume built by orienting normals *against the
coordinate origin* — which is what `PLAN.md` §7 records the reference project's `surface_manifold`
doing, flipping 106 of 320 normals on a shell translated by `3R` — would fail exactly this.

### 4.2 The nucleoplasm volume energy

The `nucleoplasm` contract says the volume constraint *"is what makes the nucleus resist compression
at constant area"*. The minimal potential with that property, and with a linear pressure response
about the rest state, is a quadratic penalty on the volumetric strain `ε_V = V/V_ref − 1`:

```
E_vol = (K_V / 2) · V_ref · ε_V²                                     [pN·µm]
```

with `K_V` in pN/µm². The induced pressure is `p = −∂E/∂V = −K_V ε_V` [pN/µm²], positive when
`V < V_ref` — the nucleus pushes back on compression. The vertex force is

```
f_a = −∂E/∂x_a = −K_V ε_V · ∂V/∂x_a ,   ∂V/∂x_a = Σ_{f ∋ a} (x_b × x_c)/6
```

with `(a, b, c)` in the face's own winding order, so a face contributes `(b×c)/6` to its first
corner, `(c×a)/6` to its second and `(a×b)/6` to its third. Setting `K_V = 0` gives exactly zero
energy and exactly zero force — the inert model is a **point in the same family**, not a separate
code path, which is what makes the negative control in §9 honest.

This is the same law the source's device kernels carry (`p_vol = −K(V − V₀)/V₀`, `f = p_vol ∂V/∂x`),
re-derived above and then compared.

### 4.3 The chromatin network

The `chromatin` contract says it is *"the nucleus's interior stiffness… why a nucleus resists
compression rather than only areal stretch."* Encode it as one radial strand per envelope vertex,
anchored at the nuclear centroid, rest length `r₀`:

```
E_chr = (G / 2) Σ_v  w_v · ((r_v − r₀)/r₀)² ,     r_v = |x_v − c|
```

`G` in pN/µm², `w_v` in µm², so `E_chr` carries pN·µm. `w_v` is the vertex's share of the
**reference** surface area (one third of each incident face's rest area), which is what makes the
response a property of the surface rather than of the vertex count: refining the mesh doubles the
number of strands and halves each one's weight.

```
f_v = −∂E/∂x_v = −(G w_v ε_v / r₀) · (x_v − c)/r_v ,   ε_v = (r_v − r₀)/r₀
```

**The weights are material, and that is a derivation and not a convenience.** If `w_v` were
recomputed from the *current* positions, the true gradient would carry a second term,

```
dE/dx = (G/2) Σ_v (dw_v/dx) ε_v²  +  G Σ_v w_v ε_v (dε_v/dx)
```

and differentiating only the strain drops the first. The nucleus contract already anticipates the
fix by declaring `face_reference_areas` as owned state: reference areas are part of the design, and
`lamina_areal_energy_and_forces` already takes them as a parameter. The convention this entry
converged on for chromatin is the one that module's author had already reached.

The centroid is a **parameter** of the function rather than something computed inside it, for the
same reason: recomputing it makes the energy depend on a perturbation twice, and the gradient check
would fail for a reason that has nothing to do with the law.

### 4.4 What this is not

Two potentials evaluated on the **same** vertex positions add, so their stiffnesses add: this
chromatin term is wired in **parallel** with the lamina, whatever the biology calls the arrangement.
A true series compliance requires the chromatin to carry its own degrees of freedom — an interior
node set — which this entry does not introduce. §14 states it and
`docs/decisions/PROPOSAL-chromatin-series-needs-interior-degrees-of-freedom.md` raises it.

### 4.5 The centroid the owner uses, and the correction it therefore owes

`chromatin_energy_and_forces` differentiates at a **fixed** anchor. That is the right shape for the
function — a finite-difference probe must be able to hold it still — and it is the wrong physics for
a body floating in space, because a nucleus whose strand anchor did not follow it would resist being
carried. `NucleusOwner` therefore uses the live centroid `c = (1/N) Σ_k x_k` and corrects for it.
With `E = Σ_v φ(|x_v − c|)`:

```
dE/dx_k = φ'(r_k) u_k + (dc/dx_k)ᵀ (dE/dc)
        = φ'(r_k) u_k − (1/N) Σ_v φ'(r_v) u_v
```

since `dc/dx_k = (1/N) I` and `dE/dc = −Σ_v φ'(r_v) u_v`. The first term is exactly what the
fixed-anchor function returns, so the entire correction is

```
f = f_fixed − mean(f_fixed)
```

— subtract the mean force, which is manifestly the projection that leaves the body with no net
internal force. **This is not a formality**, and the numbers below are labelled by which quantity
they describe, because two different peaks are in play and conflating them overstates the case:

| Quantity | Uncorrected | Corrected |
|---|---|---|
| net internal force of the **chromatin term alone** | `2.9777` pN | `3.33e-15` pN |
| peak single-vertex force of the **chromatin term alone** | `3.4230` pN | — |
| …so the spurious net is this fraction of the **chromatin** peak | **87.0%** | — |
| peak single-vertex force of the **assembled owner** (all four terms) | `376.31` pN | — |
| …so the spurious net is this fraction of the **assembled** peak | **0.79%** | — |
| net internal force of the **assembled owner** | — | `4.09e-13` pN, and see below |

Both framings are true and they answer different questions. 87% is the right number for "is this
term's own force nearly all translation?" — it is. 0.79% is the right number for "how large is the
defect against what the nucleus actually feels?" — it is small but not negligible, and it is a
*systematic* drift of the whole body rather than noise, so it does not average away over a
trajectory.

**The corrected net is a cancellation residual and must not be quoted to three significant
figures.** This row previously read `4.26e-14`, which **reproduces under no code path this module
has ever shipped** and whose provenance is unknown — it was not the pre-fix-round-1 scatter, and it
is not any combination of the two volume functions with the two scatters. Re-measured
deterministically:

```
net component-max, by (volume function, scatter)
  own (a x b).c      + own np.add.at          3.55e-13
  own (a x b).c      + membrane.volume_gradient 3.59e-13
  membrane form      + own np.add.at          4.19e-13
  membrane form      + membrane.volume_gradient 4.01e-13     <- shipped
  o.internal_forces() as the control calls it  4.0856e-13    <- the reproducible figure
```

The peak single-vertex force is `376.31` pN and the net is `~4e-13` pN: **fifteen orders of
cancellation.** Perturbing `dE_dV` by a *single ULP* moves the net between `2.9e-13` and `4.0e-13`.
So the low digits are a property of the summation path, not of the physics, and the only stable
claim is the ratio to the peak — which is what the control asserts (`< 1e-9 × peak`, i.e.
`3.76e-07` pN) and what this entry should have quoted. Recorded in
`test_the_internal_forces_sum_to_zero`'s own docstring so the next reader meets it before quoting
it again.

Every other number in this table was re-measured and reproduces exactly.

**Two controls pin the correction**, not three: `test_the_internal_forces_sum_to_zero` and
`test_the_internal_force_is_minus_the_gradient_of_the_owners_own_energy` (order 2 against the
owner's *total* energy). `test_a_rigid_translation_costs_the_owner_nothing` does **not** — measured,
it passes with the correction removed, because the correction is force-only while the energy uses
the live centroid either way. It pins the other necessary half: that the anchor follows the body at
all.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `positions`, `centroid` | µm | m | finite; `(N, 3)` float64 |
| `enclosed_volume_um3` | µm³ | m³ | `> 0` for an outward-wound closed surface |
| `bulk_modulus_pn_per_um2` (`K_V`) | pN/µm² | Pa | `≥ 0`; `0` is the inert model |
| `reference_volume_um3` (`V_ref`) | µm³ | m³ | `> 0`, refused otherwise |
| `shear_modulus_pn_per_um2` (`G`) | pN/µm² | Pa | `≥ 0`; `0` is the inert model |
| `reference_radius_um` (`r₀`) | µm | m | `> 0`, refused otherwise |
| `reference_vertex_area_um2` (`w_v`) | µm² | m² | `≥ 0`, one entry per vertex |
| `drag_pn_s_per_um` | pN·s/µm | N·s/m | `≥ 0`; consumed by the rate response, never by the energy |
| energies | pN·µm | J | finite |
| forces | pN | N | finite |

Singular and boundary cases, each with the behaviour Aleph requires:

- **A vertex coincident with the centroid.** `r_v = 0` makes the radial direction undefined.
  `chromatin_energy_and_forces` **raises** `ValueError`. It does not regularise: a strand with no
  direction has no force, and returning zero would be a silent answer to an ill-posed question.
- **`K_V = 0` / `G = 0`.** Exactly zero energy, exactly zero force, by early return. Not "small".
- **A non-positive `V_ref` or `r₀`.** Refused at card construction, so the refusal happens where the
  material is declared rather than deep inside an energy evaluation.
- **A negative modulus.** Refused at card construction; a negative stiffness is a different physical
  claim, not a parameter choice.
- **An inward-wound surface** would give a negative volume. Not silently `abs()`-ed: the sphere
  control asserts a *two-sided* band that pins the sign, and a failure there is a finding about
  `aleph.state.manifold`'s winding.

Invariants that must hold, each with the test that asserts it:

- **I1. Translation invariance of the volume.** `enclosed_volume_um3` is unchanged by a rigid
  translation to `rel=1e-12` — `tests/vertical/test_nucleus_interior.py::test_it_is_translation_invariant`.
  Exact in exact arithmetic and *conditioned* in float64: the invariance holds to `rel=1e-9` out to
  a `1e3 µm` offset and degrades past it (§10.4), asserted by
  `tests/vertical/test_nucleus_interior.py::test_the_translation_invariance_has_a_measured_range`.
- **I2. `F = −∇E`, at second order.** Central differences converge at order 2.00 over a measured
  band for both terms —
  `tests/vertical/test_nucleus_interior.py::test_the_finite_difference_agrees_and_the_order_is_two`
  and `::test_the_chromatin_finite_difference_agrees_and_the_order_is_two`.
- **I3. The discrete volume converges to the exact sphere from below.**
  `tests/vertical/test_nucleus_interior.py::test_a_sphere_recovers_four_thirds_pi_r_cubed`.
- **I4. Mesh-resolution independence of the chromatin response.** Refining the icosphere by one
  level changes the energy of a fixed uniform radial strain by less than the surface-area
  discretisation error, because the weights are areas and not counts —
  `tests/vertical/test_nucleus_interior.py::test_the_chromatin_response_is_mesh_resolution_independent`.
- **I5. The two terms are functions of one coordinate set.** `tests/vertical/test_nucleus_interior.py::test_the_two_terms_still_share_one_coordinate_set`.
  **This is a forward guard, not evidence for §4.4**, and the first draft of this entry had it
  wrong. Additivity of the second difference is a property of a *linear operator*: measured on
  `sum(x⁴)` and `sum(|x|³)` — two functions with nothing to do with each other or with a
  nucleus, and contributing 83% and 17% of the total stiffness so neither is a passenger —
  the same probe is additive to `rel = 7.9e-12`. It therefore cannot fail for the reason §4.4 cares
  about, and it survived all of M1–M10. **The parallel-not-series finding rests on the algebra in
  §4.4**, which needs no probe. What this invariant is *for* is the future: it goes red the day a
  chromatin term acquires its own interior degrees of freedom, which is exactly the change
  `PROPOSAL-chromatin-series-needs-interior-degrees-of-freedom.md` asks the PI to authorise.
- **I6. The owner's total is the exact sum of its four terms**, compared with `==` and not
  `pytest.approx` — `tests/vertical/test_nucleus_owner.py::test_the_owner_reports_the_sum_of_its_four_terms_exactly`.
- **I7. Rollback is bit-identical** —
  `tests/vertical/test_nucleus_owner.py::test_a_rejected_candidate_restores_the_owner_bit_identically`.
- **I8. The assembled nucleus is an isolated body.** Its internal forces sum to zero to `1e-9` of
  the largest single force and a rigid translation costs it exactly nothing —
  `tests/vertical/test_nucleus_owner.py::test_the_internal_forces_sum_to_zero` and
  `::test_a_rigid_translation_costs_the_owner_nothing`. This is **not** automatic: see §4.5.
- **I9. The nucleus's volume load is countable by the single-path audit.**
  `NucleusOwner.applies_osmotic_load` is `True` and `pressure_target_owner` is its own name, so
  `aleph.vertical.pressure.assert_single_pressure_path` sees it and refuses a second osmotic load on
  the same envelope — `tests/vertical/test_nucleus_owner.py::TestTheVolumeTermIsCountableByTheSinglePathAudit`.
  **Three of that class's four controls call the guard**; the fourth,
  `test_the_doubled_load_is_exactly_the_duplicated_term`, does not — it establishes that the second
  load is bit-for-bit the term the owner already carries, which is what makes the double count
  silent, and it is not killed by M8. An earlier draft called all four "controls that drive the
  guard". That is the same over-claim this entry corrected once already by measurement (§4.5's
  "three controls pin the correction", which M9 showed to be two), and it is corrected the same
  way: M8 kills three, and three is the number.

## 6. Source evidence class and known retractions

**Where I looked:** `/Users/sw1/ffn_cellsim/STATE.md` §(c) (the pushed blocklist, 17 active
entries), `STATE.md` §(b), and the call graph of each symbol named in §2, by grep over the source
tree.

- **No magnitude from the source enters this entry or this code.** Every number in §10 and in the
  controls was measured in Aleph, on Aleph's own meshes, with cards this lane chose. The blocklist
  check is therefore vacuous by construction — and it was still performed, because "we quoted
  nothing" is a claim that has to be true and not merely intended. §(c) entry 4 ("any `CONNECTED`
  rung; first composed native cell — nothing SOLVES, so no magnitude is physical") is the entry that
  would bite if a whole-cell number were lifted; none is.
- **Reachability.** `_nucleoplasm_pressure_kernel` / `nucleoplasm_volume_force_arr_kernel` are
  reachable from the source's production cell step (`ac/cell/driver.py` calls
  `cell.nucleus.accumulate(pos, f)`), which is the one place the source is unambiguously ahead of
  Aleph on the exact thing this entry builds. `ac/nucleus/chromatin_analytic.py` is referenced by
  **tests and documents only** — no production caller anywhere in that tree — so its evidence class
  is "host-side reference code with unit tests", not "runs in a cell".
- **Whether the source's tests exercise the law or the plumbing:** its WLC tests check `dE/dx = f`
  by finite difference, which is a real check of the law. Not inherited, because the law is not
  ported.
- **No retraction found** attached to the nucleoplasm volume kernels. §(c) carries no nucleus-,
  chromatin- or lamina-specific entry.
- **Evidence class of everything this entry produces: `ANALYTIC_ORACLE`.** Nothing here is a
  measurement of a nucleus, and §14 says what that forbids.

## 7. Independent oracle or derivation

Two, and the second is the point of the entry.

1. **`F = −∇E` by central finite difference, with an observed convergence order.** The step sizes
   are **measured, not inherited**: the error curve is plotted before the band is chosen. For the
   volume energy the clean `h²` band is `3e-3 … 2e-1` µm — three orders of magnitude *above* the
   `eps^(1/3) ≈ 6e-6` µm that `PLAN.md` §6.1's rule of thumb would suggest — because the energy is a
   degree-6 polynomial along a straight line whose third derivative is small against an energy of
   order `1e2`, so truncation falls under the round-off floor early. Below the band the apparent
   order reads `3.24`, then `7.39`, then `−6.12`; an error that **grows** as `h` shrinks is
   round-off, and it is distinguishable from a wrong gradient, which **plateaus**. Both failure
   signatures are recorded in §10 because telling them apart is the whole skill.
2. **The volume-constraint oracle, corrected by measurement.** The first form of this oracle said:
   squeeze the sphere at nearly constant area, and the lamina will be nearly blind while the volume
   term carries all of it. **That is false, and the reason is structural.**
   `lamina_areal_energy_and_forces` charges **per face**, so a deformation that holds the *total*
   area fixed while redistributing it among faces is one the lamina notices loudly. Measured at
   exactly constant total area (solved, `z·0.80, x·1.102189, y·1.102189`, `ΔA/A = +5.6e-07`,
   `ΔV/V = −2.81%`): `E_lamina = 268.8`, `E_volume = 10.8`, `E_chromatin = 51.7` pN·µm. The lamina
   carries **81%**, not a negligible share.

   The oracle that survives is exact rather than approximate and is stronger. **The lamina's areal
   energy is a function of the face areas and of nothing else**, so it scores any two configurations
   with identical face areas identically. A reflection through a coordinate plane is an isometry of
   every triangle, so it preserves every face area *bit-identically* — and it negates the enclosed
   volume. So there exist two configurations that a nucleus-without-an-interior scores exactly the
   same and that differ by the whole volume of the nucleus. Bound, the volume term charges
   `2 K_V V = 54635.10` pN·µm for the pair, checked against the closed form and not against a
   recorded number. No tolerance is asked to carry the argument.

Neither oracle is agreement with `ffn_cellsim`, and no control in this entry compares a number
against one from there.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_nucleus_interior.py::test_a_sphere_recovers_four_thirds_pi_r_cubed` | A level-4 icosphere of `R = 3 µm` recovers `4πR³/3` from **below**, with relative error under 1%. The two-sided band pins the sign, so an inward winding fails rather than passing on magnitude. |
| Positive | `tests/vertical/test_nucleus_interior.py::test_it_is_translation_invariant` | The enclosed volume is invariant under a `(9, −4, 2.5) µm` translation to `rel=1e-12` — the origin-referencing failure mode, tested for. |
| Positive | `tests/vertical/test_nucleus_interior.py::test_the_finite_difference_agrees_and_the_order_is_two` | The volume force is `−∇E` at observed order `2.0000` across `h = 2e-2 … 5e-3 µm`, absolute error `< 1e-4` at the coarse end. |
| Positive | `tests/vertical/test_nucleus_interior.py::test_the_chromatin_finite_difference_agrees_and_the_order_is_two` | The chromatin force is `−∇E` at observed order `2.0000` across `h = 1e-1 … 3.125e-3 µm`, with the centroid and the weights held fixed across the step. |
| Positive | `tests/vertical/test_nucleus_interior.py::test_compressing_a_nucleus_at_constant_area_costs_energy` | A prolate squeeze that moves the volume and not the area costs more than 1 pN·µm, against a rest energy that is exactly zero. |
| Forward guard, **not** evidence | `tests/vertical/test_nucleus_interior.py::test_the_two_terms_still_share_one_coordinate_set` | The radial stiffnesses add to `rel=1e-9`. Listed here for completeness and explicitly **demoted**: additivity follows from the linearity of the second difference and holds for any two energies on shared coordinates (`7.9e-12` on two unrelated ones), so it is not evidence for §4.4. It fires when a future chromatin term acquires its own degrees of freedom. See I5. |
| Positive | `tests/vertical/test_nucleus_interior.py::test_the_volume_law_has_exactly_one_implementation_in_the_package` | `nucleus_interior.enclosed_volume_um3 is membrane.enclosed_volume_um3` — the identity of the function object, not the equality of its output. §3.1. |
| Positive | `tests/vertical/test_nucleus_interior.py::test_the_volume_force_is_the_membranes_gradient_scaled_by_the_constitutive_factor` | The scatter is `np.array_equal` to `−dE/dV · membrane.volume_gradient(...)`, so this module has grown no second geometry path. |
| Positive | `tests/vertical/test_nucleus_interior.py::test_the_volume_falls_as_the_cube_of_a_uniform_radial_scaling` | §7 oracle 2's **exponent**: `V(s) = s³ V(1)` at fixed face count, to `rel=1e-12` over `s ∈ {0.5 … 2.0}` (measured worst `1.6e-15`). `test_a_sphere_recovers_four_thirds_pi_r_cubed` pins one radius; this pins the power. |
| Positive | `tests/vertical/test_nucleus_owner.py::test_the_reported_volume_falls_as_the_cube_of_a_uniform_scaling` | The same law through the owner's `enclosed_volume` state key. Replaces a control that compared the method against the function it delegates to, which could only fail if `return` stopped working. |
| Positive | `tests/vertical/test_nucleus_owner.py::test_the_doubled_load_is_exactly_the_duplicated_term` | With units converted per §3.2, an `OsmoticEnvelope(VOLUME_PENALTY)` on the nuclear envelope produces a load `np.array_equal` to the owner's own `volume_term()` force — the double count is bit-for-bit the same term, which is why it is silent. |
| Positive | `tests/vertical/test_nucleus_interior.py::test_the_chromatin_response_is_mesh_resolution_independent` | Refining the icosphere by one level (4× the faces) changes the chromatin energy of a fixed uniform radial strain by under 1%, because the weights are reference areas and not vertex counts. |
| Positive | `tests/vertical/test_nucleus_owner.py::test_the_owner_reports_the_sum_of_its_four_terms_exactly` | `potential_energy_pn_um()` equals the sum of the four term energies with `==`, so a term that is computed and never added is caught. |
| Positive | `tests/vertical/test_nucleus_owner.py::test_the_owner_declares_exactly_the_registered_contract_keys` | `owned_state_keys()` equals the registered `nucleus` contract's `owned_state`, in order. |
| Positive | `tests/vertical/test_nucleus_interior.py::test_the_translation_invariance_has_a_measured_range` | The volume is invariant to `rel=1e-9` out to a `1e3 µm` offset. Written **because a mutant survived** — see §9.1 M3. |
| Positive | `tests/vertical/test_nucleus_owner.py::test_every_term_is_live_in_the_configuration_the_controls_use` | The vacuity guard for the sum above: all four term energies exceed `1e-9` in the configuration the owner controls use, so a sum of four zeros cannot pass for a sum of four terms. |
| Positive | `tests/vertical/test_nucleus_owner.py::test_the_internal_forces_sum_to_zero` | The assembled nucleus exerts no net force on itself: the vector sum of the internal forces is below `1e-9` of the largest single vertex force. **Measured on this control's own configuration, which is the ASSEMBLED owner: peak `376.31` pN, net `4.0856e-13` pN, tolerance `3.76e-07` pN.** Two corrections live in this row. An early draft quoted `3.3e-15` against `3.423`, which are the CHROMATIN TERM IN ISOLATION (§4.5) and made the control look ~100× tighter than it is. Its replacement, `4.26e-14`, reproduced under no code path at all; the figure above is the reproducible one. **The net is a 15-order cancellation residual whose digits depend on the summation path — §4.5 — so the stable claim is the ratio, not the value.** |
| Positive | `tests/vertical/test_nucleus_owner.py::test_a_rigid_translation_costs_the_owner_nothing` | The total energy is invariant under a rigid translation to `rel=1e-12`. |
| Positive | `tests/vertical/test_nucleus_owner.py::test_the_internal_force_is_minus_the_gradient_of_the_owners_own_energy` | Order 2 over `h = 1e-2 … 2.5e-3 µm` against the owner's **total** energy, so a term whose gradient disagrees with its own energy is caught at assembly and not only in isolation. |
| Positive | `tests/vertical/test_nucleus_owner.py::test_the_declared_drag_sets_the_mobility` | `drag_pn_s_per_um` reaches the overdamped update as `1/drag`, so the card's rate number is wired rather than decorative. Asserts no timescale — see §14. |
| Positive | `tests/vertical/test_nucleus_owner.py::test_a_handle_outlives_its_topology_and_is_refused` | `ENDPOINT_CONTRACT.md` I3: a boundary handle issued against one topology generation refuses to be read after a bump. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_nucleus_interior.py::test_a_nucleus_with_no_interior_is_free_to_be_compressed` | **The control that would have caught this gap**, in the exact form §7 arrived at. A reflection through the xy-plane leaves every face area bit-identical (asserted with `np.array_equal`, not a tolerance), so the lamina scores the nucleus and its own inversion *identically* — it is blind to its own volume. The enclosed volume meanwhile negates exactly. Bind the interior and the pair costs `2 K_V V = 54635.10` pN·µm, checked against the closed form. If the interior ever stops being bound, the second half fails; if the lamina ever silently acquires a volume response, the first half fails. |
| Negative (must fail) | `tests/vertical/test_nucleus_interior.py::test_the_interior_charges_for_a_compression_the_lamina_charges_nothing_for` | The ordinary-deformation companion, and it records the correction of §7: at *exactly* constant total area the nucleus loses 2.81% of its volume, a nucleus with no interior charges exactly `0.0` for that loss, and the bound interior charges `62.55` pN·µm. |
| Negative (must fail) | `tests/vertical/test_nucleus_interior.py::test_the_error_curve_falls_rather_than_plateauing` | The chromatin gradient's error must fall by more than 100× over a 32× reduction in `h`. A plateau is a wrong analytic gradient; a rise is round-off; the two are different and this control asserts the *shape* of the curve rather than one point on it. |
| Negative (must fail) | `tests/vertical/test_nucleus_owner.py::test_the_reduced_coordinates_key_is_declared_and_refused` | `reduced_coordinates()` raises `NotImplementedError`. A silently absent modal reduction is indistinguishable from one that reduces to the identity. |
| Negative (must fail) | `tests/vertical/test_nucleus_owner.py::test_a_rollback_without_a_snapshot_is_refused` | An un-advance with nothing to un-advance to raises rather than leaving the state where the rejected candidate put it. |
| Negative (must fail) | `tests/vertical/test_nucleus_interior.py::test_a_zero_bulk_modulus_makes_the_resistance_vanish_exactly` | `K_V = 0` gives energy exactly `0.0` and a force array with exactly zero non-zeros — driving the shipped path, not a copy of it. |
| Negative (must fail) | `tests/vertical/test_nucleus_interior.py::test_a_zero_shear_modulus_makes_the_chromatin_vanish_exactly` | The same statement for `G = 0`. |
| Negative (must fail) | `tests/vertical/test_nucleus_interior.py::test_a_vertex_on_the_centroid_is_refused` | A vertex coincident with the centroid raises `ValueError` rather than returning a regularised force for an undefined direction. |
| Negative (must fail) | `tests/vertical/test_nucleus_interior.py::test_a_card_with_an_impossible_material_is_refused` | A negative modulus, a non-positive reference volume and a non-positive reference radius are each refused at card construction. |
| Negative (must fail) | `tests/vertical/test_nucleus_owner.py::test_an_owner_scheduled_in_a_phase_it_declares_no_work_for_refuses` | The owner raises rather than silently doing nothing in a phase it has no work in. |

### 9.1 Mutation table — measured against the shipped source

**Ten** mutants, re-run after fix round 2 under `PYTHONDONTWRITEBYTECODE=1`, against the whole of
`tests/vertical/test_nucleus_interior.py` and `tests/vertical/test_nucleus_owner.py` (**47
controls**, baseline green). Each mutant is applied to the **shipped** source, reverted, and the
file confirmed byte-identical to `HEAD` with `git diff --exit-code`; every target file is also
checked clean *before* each mutant, so a concurrent lane's edit can neither be mistaken for a result
nor destroyed by the study. **No test file was modified at any point** — mutating a control measures
the mutation and not the suite.

> **Methodological requirement, learned the expensive way and recorded so nobody repeats it:
> `PYTHONDONTWRITEBYTECODE=1` is mandatory for a mutation study.** CPython's bytecode cache
> validates a `.pyc` against the source's **mtime in whole seconds plus its size in bytes**. A
> mutant that preserves byte length — which most single-token mutations do — and is applied and
> reverted inside the same wall-clock second therefore leaves a `.pyc` that the interpreter
> considers *valid*, and the next run silently executes the wrong bytecode. The symptom is
> phantom kills on a git-clean tree, which looks exactly like a working study. An independent
> re-review of this entry hit precisely that on its first pass. The numbers below were re-measured
> under `PYTHONDONTWRITEBYTECODE=1` and are unchanged, but "unchanged" is a measurement and not an
> assumption, which is why the re-run happened.

**M2 and M3 moved to `aleph/vertical/membrane.py`**, because §3.1 moved the volume law there. A
mutation of the volume law is now a mutation of that file, and the question these two answer is
whether *this entry's* controls still catch a corrupted volume once the code is somebody else's.
They do.

**M8 and M9 are new in this round**, covering the code fix round 1 added: the osmotic marker and the
self-force correction. New code gets mutants or it is not covered.

**Result: 9 killed, 1 equivalent.** The one survivor is a finding and is written up below rather
than dressed as a pass.

| # | File | Mutation | Verdict | Controls that died |
|---|---|---|---|---|
| M1 | `nucleus_interior` | drop the `0.5` from the volume energy | **killed** (4) | `test_the_finite_difference_agrees_and_the_order_is_two`, `test_the_internal_force_is_minus_the_gradient_of_the_owners_own_energy`, `test_a_nucleus_with_no_interior_is_free_to_be_compressed`, `test_the_interior_charges_for_a_compression_the_lamina_charges_nothing_for` |
| M2 | `membrane` | `enclosed_volume_um3`: swap the two crossed vertices | **killed** (33) | `test_a_sphere_recovers_four_thirds_pi_r_cubed` first — the volume changes sign and the two-sided band is what sees it. The owner then fails to construct at all, hence 33. |
| M3 | `membrane` | `enclosed_volume_um3`: subtract the centroid first | **SURVIVED — equivalent** | none, of 47. §9.2. |
| M4 | `nucleus_interior` | chromatin weights `1.0` instead of the reference vertex area | **killed** (3) | `test_the_chromatin_response_is_mesh_resolution_independent` (**written for this mutant**), `test_a_uniform_radial_contraction_costs_energy`, `test_the_interior_charges_for_a_compression_the_lamina_charges_nothing_for` |
| M5 | `nucleus_interior` | return `+dE_dV · ∂V/∂x` instead of `−` | **killed** (4) | `test_the_finite_difference_agrees_and_the_order_is_two`, `test_the_internal_force_is_minus_the_gradient_of_the_owners_own_energy`, `test_the_volume_force_is_the_membranes_gradient_scaled_by_the_constitutive_factor`, `test_the_doubled_load_is_exactly_the_duplicated_term` |
| M6 | `nucleus_interior` | recompute the chromatin centroid internally | **killed** (3) | `test_the_chromatin_finite_difference_agrees_and_the_order_is_two`, `test_the_error_curve_falls_rather_than_plateauing`, `test_a_vertex_on_the_centroid_is_refused` |
| M7 | `nucleus_interior` | recompute the vertex weights from the **current** positions | **killed** (5), **by the right failure mode** | `test_the_chromatin_finite_difference_agrees_and_the_order_is_two`, `test_the_error_curve_falls_rather_than_plateauing`, `test_a_uniform_radial_contraction_costs_energy`, `test_the_interior_charges_for_a_compression_the_lamina_charges_nothing_for`, `test_the_internal_force_is_minus_the_gradient_of_the_owners_own_energy` |
| M8 | `nucleus` | drop `applies_osmotic_load`, hiding the volume term from the pressure audit | **killed** (3) | `test_the_audit_refuses_a_second_osmotic_load_on_the_nuclear_envelope`, `test_the_audit_accepts_the_nucleus_carrying_its_volume_term_alone`, `test_the_audit_is_single_compartment_and_two_compartments_need_two_calls` |
| M9 | `nucleus` | drop the mean-subtraction self-force correction | **killed** (2) | `test_the_internal_forces_sum_to_zero`, `test_the_internal_force_is_minus_the_gradient_of_the_owners_own_energy` |
| M10 | `membrane` | `enclosed_volume_um3` := area × mean radius — an impostor that is exact for a sphere, translation invariant, and degree-three homogeneous | **killed** (5) | **`test_a_linear_map_scales_the_volume_by_its_determinant`** (added for it), `test_a_nucleus_with_no_interior_is_free_to_be_compressed`, `test_the_interior_charges_for_a_compression_the_lamina_charges_nothing_for`, `test_the_finite_difference_agrees_and_the_order_is_two`, `test_the_internal_force_is_minus_the_gradient_of_the_owners_own_energy` |

**M10 is the one that changed a verdict.** Before fix round 2 this impostor **survived the entire
`R³` oracle** — it passes `test_the_volume_falls_as_the_cube_of_a_uniform_radial_scaling` at every
scale factor to twelve digits, and it passes `test_a_sphere_recovers_four_thirds_pi_r_cubed`,
because a uniform scaling about the centroid cannot distinguish the enclosed volume from *any*
translation-invariant degree-three-homogeneous formula. The cube law pins homogeneity; it does not
pin the volume. The determinant control does, and M10 is the evidence that it does rather than the
claim that it does.

**M9 confirms the count in §4.5 rather than asserting it.** Exactly two controls died, and
`test_a_rigid_translation_costs_the_owner_nothing` was **not** among them — which is the measurement
behind saying two controls pin the correction and not three.

**M8 is the one that matters for Track C.** It is the state this entry shipped in before fix round
1: the owner applies a volume load and the single-path audit cannot see it. Three controls now die
on it, and they die because they *drive the guard* — they build the doubled world and require
`PressureDoubleCountError` — rather than asserting that a boolean attribute has a particular value.

**M7's failure mode was confirmed, not assumed.** The error curve under the mutant was re-measured
while it was applied, and it is a plateau and not a rise:

```
h          error        order
1.0e-01    2.343113e-03    --
5.0e-02    2.531455e-03  -0.1115
2.5e-02    2.578542e-03  -0.0266
1.25e-02   2.590314e-03  -0.0066
6.25e-03   2.593257e-03  -0.0016
3.125e-03  2.593993e-03  -0.0004
1.0e-03    2.594213e-03  -0.0001
1.0e-04    2.594238e-03  -0.0000
```

The error converges to `2.5942e-03` and stops. Round-off *grows* as `h` shrinks (see §10.1, which
reaches order `−6.12`); a wrong analytic gradient does not move at all. A control that killed this
for the wrong reason would not catch its next variant, which is why
`test_the_error_curve_falls_rather_than_plateauing` asserts the *shape* of the curve.

### 9.2 M3 survived, and it is an equivalent mutant — the proof, and what was written instead

**The build plan expected `test_it_is_translation_invariant` to kill M3. It cannot, and the reason
is the identity that makes the volume translation-invariant in the first place.** For a closed,
consistently wound surface `Σ_f n_f dA = 0`, so

```
(1/6) Σ_f ((a−c₀) × (b−c₀)) · (c−c₀)  =  (1/6) Σ_f (a × b) · c    for every c₀
```

Subtracting the centroid is therefore **another way of being translation invariant**, which makes
the invariance test the one control structurally incapable of distinguishing the two forms. M3
changes no value on this function's documented domain, and all 47 controls passing under it is the
correct outcome rather than a hole.

Where the two forms genuinely differ is **conditioning**, and there the shipped form is the worse of
the two: the triple products grow as the cube of the offset while their sum stays of order `1e2`.
Measured relative error against the origin value (§10.4). So rather than manufacture a control that
could only assert an implementation detail, this lane wrote
`test_the_translation_invariance_has_a_measured_range`, which pins the range over which the shipped
form's invariance actually holds. That is a real strengthening of the suite and it records a limit
nobody had written down; it is **not** a kill, and this entry does not claim one.

## 10. Numerical and precision envelope

**Working and accumulation precision:** float64 throughout, host-side NumPy. No reduction is done in
a lower precision than the array it reduces; `enclosed_volume_um3` accumulates the face triple
products with `np.einsum`, whose pairwise summation is what keeps the level-4 sphere's relative
error at the discretisation level rather than at the accumulation level.

### 10.1 The volume gradient's error curve, measured

Level-2 icosphere, `R = 3 µm`, jittered by `N(0, 0.05 µm)`, `K_V = 250 pN/µm²`,
`V_ref = 0.9 · V(x)`:

```
h          error       order
2.0e-01    1.346e-03     --
1.0e-01    3.366e-04   2.0000
5.0e-02    8.414e-05   2.0000
2.5e-02    2.103e-05   2.0000
1.25e-02   5.259e-06   2.0000
6.25e-03   1.315e-06   2.0000
3.125e-03  3.289e-07   1.9990
1.0e-03    3.475e-08   3.2425   <- round-off floor begins
4.0e-04    2.077e-10   7.3865   <- meaningless
1.0e-04    1.449e-08  -6.1245   <- the error GROWS as h shrinks
```

The asserted band is `2e-2 … 5e-3`, in the middle of the clean region, with the order asserted in
`(1.95, 2.05)` and the coarse-end absolute error under `1e-4`. A looser tolerance would admit M1 and
M5; a tighter one would be asserting the round-off floor.

### 10.2 The chromatin gradient's error curve, measured

Same mesh, `G = 120 pN/µm²`, `r₀ = 3 µm`, weights and centroid held fixed: order `2.0000` across
`1e-1 … 3.125e-3`.

### 10.3 What a wrong gradient looks like, so it is not mistaken for round-off

The first version of `chromatin_energy_and_forces` computed the weights from live positions and
differentiated only the strain. Its error curve:

```
h          error        order
1.0e-01    5.248e-03      --
5.0e-02    5.171e-03    0.0213
2.5e-02    5.152e-03    0.0054
1.25e-02   5.147e-03    0.0013
3.125e-03  5.146e-03    0.0001
```

**The error converges to a non-zero constant.** Not a floor that rises, not a slope of 2 — a
plateau. Round-off *grows* as `h` shrinks; a wrong analytic gradient *does not move*. `PLAN.md` §6.1
records the project losing a night to the other half of this confusion, and the rule that transfers
is: never choose a finite-difference step from a formula without plotting the curve first. The
optimum is a property of the function, not of `eps`.

### 10.4 The volume's translation invariance is exact in theory and conditioned in float64

Relative error of `enclosed_volume_um3` against its value at the origin, for a level-2 sphere
translated along `(1, −0.4, 0.25)`, comparing the shipped form with the centroid-shifted form of
§9.2's M3:

```
offset [um]   shipped      centroid-shifted
9.0e+00       8.9e-16      4.4e-16
1.0e+02       5.0e-13      4.4e-16
1.0e+03       5.2e-11      2.7e-15
1.0e+04       3.2e-07      6.1e-14
1.0e+05       3.3e-04      4.2e-13
1.0e+06       4.5e-01      3.0e-12
```

`test_the_translation_invariance_has_a_measured_range` asserts `rel < 1e-9` at `1e2` and `1e3` µm,
which brackets any configuration a cell simulation occupies. Past `1e4` µm the divergence-theorem
form needs its origin moved to the centroid; that is recorded here rather than fixed, because
changing it now would be changing shipped numerics for a regime nothing in this project uses.

**Range of validity and behaviour outside it.** The volume and chromatin laws are polynomial and
smooth for any finite configuration, so there is no strain beyond which the *expression* degrades;
what degrades is the *modelling claim*, and §14 carries that. The refusals in §5 are hard errors
rather than silent degradation.

## 11. Production-backend residency and transfer

Host, CPU, NumPy float64, in the same residency as every other `aleph/vertical/` law. The energy
functions are pure — arrays in, `(float, ndarray)` out, nothing cached, nothing mutated in place
except the accumulator each one allocates itself — so they are safe to call from a finite-difference
probe, which is precisely how the controls use them.

`NucleusOwner` holds the only mutable state: `positions`, `forces`, one snapshot buffer, and the
reference quantities fixed at construction. No device transfer, no host round-trip, **no GPU** —
nothing in this entry was run on one and nothing in it needs one. If a device backend is added, the
volume reduction is the term that needs attention: it is a global scalar reduction over faces inside
a force evaluation, which is the shape that costs a synchronisation.

## 12. Comments and docstrings to discard

Discarded from the source, verified absent by grep over both new modules:

- The device/gate vocabulary attached to the volume kernels — the zero-round-trip gate label, the
  `wp.kernel` / `wp.tid()` / `wp.atomic_add` framing, the "one thread" note, and the sibling-kernel
  cross-references. Aleph's version is a host reduction and a `np.add.at` scatter, and it says so.
- The source's provenance prose on the WLC module: its `k_BT` value and 310 K convention, its
  persistence-length parameterisation, its internal parameter-set labels and the identifier for the
  parameter file the small-strain stiffness feeds. None of it is used, because the WLC is not
  ported.
- Every source module path, class name and internal document reference. `tests/ports/test_port_discipline.py`
  enforces this mechanically over `aleph/**`; naming the source is required *here* and forbidden
  *there*.

**What replaces them:** a module docstring that states why the interior is a separate module from
the lamina in terms of what it lets a control ask; the derivation in §4 restated in the functions'
own docstrings; and — the part with no counterpart in the source — the measured plateau of §10.3,
written into `chromatin_energy_and_forces` where the next author will meet it.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | `2026-07-31` — `tests/vertical/test_nucleus_interior.py` (21) and `tests/vertical/test_nucleus_owner.py` (26) pass, 47 controls; gradient order `2.0000` for both interior terms over the bands in §10 and order 2 for the owner's total; **9 of 10 mutants killed and M3 proved equivalent** (§9.2), M4 by a control written for it and M7 by the plateau signature and not merely by a failure. |
| Reviewer | **Agent-proposed and unratified.** Written and measured by session `46143f30`; no PI review has happened and this entry carries no `decided_by` field. |
| Status, and why not the one the build plan asked for | The plan's Task A5 step 4 says to set `LANDED`. `ports/regenerate_index.py:64` maps `LANDED → ACCEPTED`, and `ACCEPTED` is a ratification — `ports/TEMPLATE.md` §13 asks for a reviewer and says agent-proposed entries "remain unratified". An agent stamping its own port `ACCEPTED` is the same move as writing `decided_by: PI`. So this entry stays at `PROPOSED`, matching all twelve recent lane entries in `ports/ledger/`. The divergence from the plan is deliberate and is reported rather than taken. |
| Rollback | Delete `aleph/vertical/nucleus_interior.py`, `aleph/vertical/nucleus.py` and their two test files, and remove `owned_state_keys()` from `nucleus_lamina.py` / `membrane.py` / `cortex.py`. What breaks: `aleph.scenarios.audit` returns those three owners to "declares no state keys at all", and the nucleus returns to resisting areal stretch only. Nothing else imports these modules — no connector is wired to the owner and `COMPARTMENT_MODULES` is unchanged (see §14). |

## 14. Honest limits

- **`K_V`, `G`, `r₀` and `V_ref` are a material card, not a measurement.** Every value used in the
  controls was chosen to make a control sharp. **No stiffness computed from them may be reported as
  a property of a real nucleus**, in any document, at any rung. The cards are `UNSOURCED`; nothing
  in this entry attempts a literature reading and none has happened.
- **This is a parallel wiring, not a series compliance.** Both terms are evaluated on the envelope's
  own vertices, so their energies add and their stiffnesses add. **The support for that is the
  algebra of §4.4, not a measurement** — the stiffness probe that appeared to measure it is additive
  for any two energies on shared coordinates and so cannot fail for this reason (I5). The census's
  "in series" describes the
  *physical* arrangement, interior behind surface; encoding it as two additive potentials on shared
  coordinates is a **modelling choice**. A true series response needs interior degrees of freedom,
  which is a different and larger port. Raised in
  `docs/decisions/PROPOSAL-chromatin-series-needs-interior-degrees-of-freedom.md`; it is a question
  for the PI and this entry does not answer it.
- **`reduced_coordinates` is declared and NOT implemented.** The contract's modal reduction is not
  in this entry at all. `owned_state_keys()` declares the contract's nine keys because that is what
  the audit compares against the registry; `NucleusOwner` carries eight of them and refuses to
  pretend about the ninth.
- **The rate dependence is declared and only half-carried.** `drag_pn_s_per_um` sits on the card and
  reaches the owner's overdamped update as a mobility; the nucleoplasm contract's *"relaxation takes
  time rather than happening instantly"* is therefore representable but has **no control asserting a
  timescale**, because a relaxation parameter in this vertical is not evidence about real
  timescales. `UNVERIFIED`.
- **Task A6 step 3 is deferred, deliberately.** `aleph/scenarios/audit.py`'s
  `COMPARTMENT_MODULES["nucleus"]` still points at `aleph.vertical.nucleus_lamina`, because the
  scenario lane holds that file with uncommitted work (`docs/ACTIVE_SESSIONS.md`, `CLAUDE.md` §1).
  So `owned_state_keys()` stays in `nucleus_lamina.py` where the audit reads it, and `nucleus.py`
  re-exports the same tuple from that one definition rather than owning a second copy. Until the
  repoint happens, **the audit reports the lamina module and not the owner**, which is a smaller
  claim than this entry could make.
- **No connector is wired to this owner.** `nucleus_cytosol_boundary` remains unwired.
  `NucleusOwner.boundary_endpoint()` publishes the solid half — face set, outward normals, epoch
  stamp — and there is no `Binding` row in `wiring.py`, per `ENDPOINT_CONTRACT.md` §3 rule 5. The
  connector is Track C's.
- **Nothing here has been exercised in a whole-cell run.** The controls are unit-level and
  `ANALYTIC_ORACLE` rung. `aleph.scenarios.audit`'s own note applies: schedulable is not correct and
  exercised is not correct either.
- **The volume's translation invariance is conditioned, and past `1e4 µm` from the origin it stops
  holding to any useful tolerance.** Measured in §10.4, asserted only over the range a cell
  occupies, and **not fixed** — the better-conditioned form exists (§9.2) and adopting it would
  change shipped numerics for a regime nothing here uses.
- **One mutant is equivalent rather than killed.** §9.2 M3. The suite's kill rate is
  9/10 and M3 is unkillable by any behavioural control on the documented domain. Reporting
  it as 10/10 would have been the easy sentence and it would have been false.
- **`assert_single_pressure_path` is a single-compartment audit and a cell with a nucleus has two.**
  The guard asks "is there exactly one osmotic load, and does it land on the declared envelope?" A
  world holding both the cortex turgor and this owner has two legitimate loads on two compartments,
  and an unpartitioned call over it refuses. The usage that works today, guard unmodified, is one
  call per compartment over candidates partitioned by `pressure_target_owner`; both behaviours are
  pinned by `test_the_audit_is_single_compartment_and_two_compartments_need_two_calls`. **Whether
  the guard should learn about compartments is a PI question this entry does not answer.** It was
  not modified: `CLAUDE.md` §2 rule 4, and it belongs to `ALEPH-PORT-1102`.
- **The port class is arguably wrong for the volume half and was not changed.** §3 now records that
  the volume construction was in-repo all along, which makes that half `RE-DERIVED` rather than
  `SOURCE-DERIVED`. The entry keeps its declared class because the chromatin half genuinely took a
  structural idea from the source and because re-classifying an entry after review is not an agent's
  call. Flagged for the reviewer.
- **`NucleusOwner` now carries `triangles` and `volume_um3()`, which makes the double-count hazard
  reachable.** Before them an `OsmoticEnvelope` could not be attached to this owner at all — but
  only by `AttributeError`, which is accidental protection one `@property` away from vanishing and
  which no reader would have recognised as load-bearing. The real protection is
  `applies_osmotic_load` plus the audit, and it is now exercised rather than theoretical (I9).
- **The comparison with the source was structural, not numerical.** The volume gradient's *form* was
  compared against the source kernel; no number was compared against any number of theirs, and
  agreement with `ffn_cellsim` is not evidence for anything in §8 or §9.
