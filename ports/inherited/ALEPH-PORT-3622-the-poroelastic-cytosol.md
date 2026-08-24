# ALEPH-PORT-3622 — the poroelastic cytosol: the diffusivity the coupling was never asked for

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3622` |
| Lane | `46143f30` — the poroelastic-cytosol diffusivity lane |
| Status | `PROPOSED` |
| Written | `2026-08-04` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` — nothing was ported; the derivation in §4 is carried out here against Aleph's own discrete operators |
| Aleph target | `aleph/vertical/cytosol.py` (append-only), `tests/vertical/test_cytosol_poroelastic_diffusivity.py` (new) |
| Supersedes | Nothing. **Extends `ALEPH-PORT-2401`**, which is where the Biot coupling already lives. |

---

## 0. The brief's premise was false, and that is the first finding

This lane was briefed as follows: *"the porous SOLID skeleton has no elasticity inside the cytosol
… our `mobility` has no E in it … nothing represents the meshwork INSIDE the cytosolic volume."*
It asked for the solid phase of a two-phase model to be built.

**The solid phase already exists and has since 2026-07-30.** It was checked before a line was
written, by enumeration rather than by reading the brief:

| Claim in the brief | Measured state of the tree |
|---|---|
| "the porous solid skeleton has no elasticity inside the cytosol" | `PoroelasticCard.oedometric_modulus_pn_per_um2` is a **required** field with no default. `CytosolField.assemble_coupled_system` puts `-(V·H/dt)·DᵀD` in the solid block. |
| "nothing represents the meshwork inside the cytosolic volume" | A staggered skeleton displacement lives on the coupling-axis faces; `solid_dof_count` counts it; `dilatation_operator` maps it to cell dilatation. |
| "an elastic stress that couples to the pressure field … compressing the cytosol raises pressure through the SKELETON" | `fluid_content()` is `ζ = S·p + α·div u`. The two off-diagonal blocks are exact transposes (`coupling_blocks`), the system is solved monolithically in one `np.linalg.solve`, and `effective_biot_coefficient` exists specifically because an earlier version's break left `α·θ` on the right-hand side. |
| "a negative control that would catch a decoupled skeleton … plant it and kill it" | `CytosolField.decouple_solid` is exactly that, shipped as a flag on the owner so the control drives the production code path. |

`tests/vertical/test_cytosol_controls.py` is **47 controls, all green**, measured this session.

Per `CLAUDE.md` §2.5 a contradiction is a finding and not a licence, so the lane was **rescoped**
rather than either (a) rebuilding a working coupling on top of itself or (b) declaring the work
done and landing nothing. What follows is the one part of the brief that names something the tree
genuinely does not have.

## 1. Aleph API

The gap is narrow and exact. `PoroelasticCard` reports **one** diffusivity:

```python
@property
def consolidation_coefficient_um2_per_s(self) -> float:
    return self.mobility_um4_per_pn_s / self.storage_um2_per_pn        # c_v = mobility / S
```

and its own docstring names the trap it leaves open — *"it is **not** the coupled problem's time
constant: the coupled system's response also carries the skeleton's stiffness through `alpha` and
the oedometric modulus. Reading a consolidation time off this alone is the decoupled reading of a
coupled model."* The module states the trap and then **provides no way out of it.** Enumerated,
not pattern-matched: `PoroelasticCard` has exactly six public members and `CytosolField` sixty-eight,
and no member of either returns a diffusivity that contains the modulus.

So the brief's one true sentence is *"our `mobility` has no E in it"* — true of the **reported
number**, false of the **solved system**. This entry closes that, append-only:

```python
# EDITED — appended to aleph/vertical/cytosol.py; no existing line is changed
darcy_mobility_um4_per_pn_s(...)                     # module function: k/mu from pore geometry
PoroelasticCard.skeleton_storage_um2_per_pn          # alpha^2 / H_oed
PoroelasticCard.poroelastic_diffusivity_um2_per_s    # D_p, the coupled diffusivity
CytosolField.poroelastic_diffusivity_um2_per_s       # D_p as the ASSEMBLED SYSTEM sees it
```

```python
# CONSUMED, READ-ONLY — every one of these is another lane's landed work
CytosolField.assemble_coupled_system, .step, .dilatation, .effective_biot_coefficient,
             .decouple_solid, .cell_centres_um, .set_pressure, .pressure
assert_uniform_pressure_carries_no_darcy_flux, PoroelasticCard.consolidation_coefficient_um2_per_s
```

**`tests/vertical/test_cytosol_controls.py` is read and driven, never edited.** It is
`ALEPH-PORT-2401`'s acceptance; adding this lane's assertions to it would launder new work through
another entry's controls. The new controls live in their own module.

**No default card, and no coefficient literal, enters `aleph/**`.** `ALEPH-PORT-1803` forbids it
while the organelle homogenization map is undocumented, and the existing control asserting that
`default_poroelastic_card` does not exist stays green. `darcy_mobility_um4_per_pn_s` is
keyword-only with **every argument required**: it is a conversion, not a source of numbers. Every
biological value in §10 lives in the test module and in this file, where it is labelled.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` — cited for the record only |
| Source path and symbol | None. No `ffn_sim/` path was opened for this entry. |
| What was taken | **Nothing.** No line, no identifier, no constant, no coefficient. |
| Literature consulted | Moeendarbary et al. 2013 *Nat. Mater.* 12:253–261; Mogilner & Manhart 2018 *Annu. Rev. Fluid Mech.* 50:347–370; Strychalski 2021 3D bleb model. Cited for their claims in §6; none of their numbers is hard-coded in `aleph/**`. |

## 3. Why source-derived porting beats clean-room

**It does not, here — nothing was ported and the entry is `RE-DERIVED` from scratch.** The result
in §4 is four lines of algebra on operators Aleph already owns (`D`, `Dᵀ`, `Lap`), and a port would
have supplied a formula without the discrete statement that makes it checkable on *our* grid. The
literature's `D_p ~ E·ξ²/μ` is a scaling written for a continuum; §4.4 shows it is the
incompressible corner of what the discrete system actually runs at, which is a stronger claim than
copying the scaling would have produced and a weaker number than the scaling implies.

## 4. Physical or mathematical law represented

### 4.1 What the assembled system already is

From `assemble_coupled_system`, multiplying the solid row by `-dt/V` and using
`θ = D·u + clamp`, the two rows of the landed system are exactly:

```
solid    Dᵀ (H θ  -  α p)  =  f_ext / V                      (1)
fluid    [ S (p¹ - p⁰) + α (θ¹ - θ⁰) ] / dt  =  - mob · Lap p¹   (2)
```

with `Lap = GᵀG` positive semi-definite (so `-Lap` is the diffusion operator), `H` the drained
oedometric modulus, `α` the Biot coefficient, `S` the storativity, `mob = k/μ` the Darcy mobility.
Equation (1) is the quasi-static effective-stress balance; (2) is the fluid-content balance. This
is the coupling the brief asked to be built, and it is already there.

### 4.2 The elimination the module never performs

`Dᵀ` maps cell fields to coupling-axis face fields, `(Dᵀs)_f = (s_lo - s_hi)/dx`, so its null space
is **the fields constant along the coupling axis**. With no external solid load, (1) says

```
H θ - α p  =  c ,      c constant along the coupling axis                    (3)
```

Every column of `D` sums to exactly zero, so `Σθ = 0` on each axis line when the confining clamp is
zero; summing (3) over a line therefore fixes `c = -α·mean(p)`, and for a zero-mean pressure mode
`c = 0`. Hence

```
θ  =  (α / H) · p                                                            (4)
```

— *the skeleton's dilatation is slaved to the pressure by its own stiffness.* Substituting (4) into
(2):

```
( S + α²/H ) · ∂p/∂t  =  mob · ∇²p                                           (5)
```

which is a diffusion equation whose diffusivity is **not** `mob/S`:

```
        mobility                    k / μ
D_p = ------------------  =  -------------------                             (6)
       S  +  α² / H          S  +  α² / H_oed
```

`α²/H_oed` is the storage the **skeleton** contributes — the volume the meshwork gives up per unit
pressure because it is elastic rather than rigid. Setting `α = 0` (the `decouple_solid` break)
returns `D_p = mob/S = c_v` exactly, which is the sense in which the existing reported number is
the decoupled reading.

### 4.3 Units, in ours

`mob` [µm⁴/(pN·s)] ÷ ( `S` [µm²/pN] + `α²`[1] / `H` [pN/µm²] ) = µm⁴/(pN·s) · pN/µm² = **µm²/s**.
Both terms of the denominator are µm²/pN, so they are addable, which is the dimensional check that
the skeleton term is a *storage* and not a modulus in disguise.

### 4.4 Recovering the literature's scaling, and what it costs

Moeendarbary et al. write the poroelastic diffusion constant as `D_p ~ E·ξ²/μ`. Take (6) to the
incompressible corner — solid grains and pore fluid both far stiffer than the meshwork, so `S → 0`
and `α → 1` — and use the pore-scale permeability scaling `k = a·ξ²`:

```
D_p  →  mob · H  =  (a ξ² / μ) · H_oed   ≈   E ξ² / μ                        (7)
```

So the literature's formula is the `S → 0, α → 1` corner of (6), with `E` playing the oedometric
modulus and `a` a pore-shape factor of order one. This matters in two directions:

* it makes `E` *derivable in our parameters* rather than asserted — which is what the brief asked
  for — and
* it says the scaling **over**-estimates `D_p` whenever the storativity is not negligible, by the
  factor `(1 + S·H/α²)`. At the §10 configuration that factor is `1.10`, so the literature formula
  is 10 % high there. A number that agrees with a scaling law to within 10 % is agreeing with the
  scaling law's corner, not with a measurement, and §14.2 says so.

`darcy_mobility_um4_per_pn_s(pore_size_um=ξ, fluid_viscosity_pn_s_per_um2=μ, pore_shape_factor=a)`
returns `a·ξ²/μ`. It has no defaults; `a` is required precisely because it is the unsourced number
in (7) and hiding it in a default would make `k = ξ²` read as a measurement.

## 5. Units, domains, singular cases, invariants

| Symbol | Meaning | Unit | Domain |
|---|---|---|---|
| `mob` | Darcy mobility `k/μ` | µm⁴/(pN·s) | `> 0` |
| `S` | storativity `1/M` | µm²/pN | `> 0` |
| `H_oed` | drained oedometric modulus | pN/µm² (`= 1 Pa`) | `> 0` |
| `α` | Biot coefficient | — | `(0, 1]` |
| `α²/H_oed` | **skeleton storage** | µm²/pN | `> 0`, `= 0` only under the break |
| `D_p` | coupled poroelastic diffusivity | µm²/s | `0 < D_p < c_v` |

**Singular cases, all reachable and all already refused by `PoroelasticCard.__post_init__`:**
`H_oed → 0` (a skeleton with no stiffness) would send `α²/H → ∞` and `D_p → 0`; the constructor
refuses `H_oed ≤ 0`, so the limit is approached and never attained. `S → 0` is refused for the same
reason, which means `D_p = mob·H/α²` is a limit this module cannot construct and (7) is therefore an
asymptote rather than a configuration. No new refusal is added by this entry — every argument that
could produce a division by zero is already validated by the landed constructor, and a second
validation would be a second place to keep in step.

**Invariants this entry adds, each machine-checked:**

1. `D_p < c_v` **strictly**, for every constructible card. The skeleton can only ever slow the
   pressure down, because `α²/H > 0` on the whole constructible domain.
2. Under `decouple_solid=True`, `CytosolField.poroelastic_diffusivity_um2_per_s == c_v` **exactly**
   — asserted with `==`, not a tolerance, because `effective_biot_coefficient` returns a hard `0.0`
   and `0.0*0.0/H` is exactly `0.0`. An exact zero is a claim (`PLAN.md`), so it is claimed exactly.
3. The **measured** decay rate of a pressure mode under `CytosolField.step` equals the **reported**
   `D_p`. This is the one that has teeth; §9 is about it.

The four invariants the brief names as load-bearing are `ALEPH-PORT-2401`'s and are **re-asserted
here on this entry's own configurations** rather than assumed to carry: a uniform pressure field
carries exactly zero constituent Darcy flux, the prescribed boundary flux defaults to exactly `0.0`,
the divergence-theorem residual is exactly `0.0`, and every one of those reads the constituent scale
and never the resultant.

## 6. Source evidence class and known retractions

Evidence rung **`ANALYTIC_ORACLE`**; quantitative status **`BLOCKED`**, unchanged from
`ALEPH-PORT-2401` and for its reason: `ALEPH-PORT-1803` holds that this field's porosity,
permeability and viscosity cannot be interpreted biologically while the organelle homogenization map
is undocumented. What is established here is a **relation between coefficients and solver
behaviour**, not a measurement of cytoplasm.

Literature claims, cited and not re-derived:

* **Moeendarbary et al. 2013** — cytoplasm behaves as a poroelastic material; deformation rate is
  limited by intracellular water redistribution through the meshwork; reported `D_p` of order
  1–100 µm²/s by microindentation. **Used for**: the plausibility band in §10 and the scaling in
  §4.4. **Not used for**: any number in `aleph/**`.
* **Mogilner & Manhart 2018** — the model hierarchy from cytoplasmic flow through two-phase
  poroviscous to poroelastic. **Used for**: placing this module in the hierarchy (§14.1).
* **Strychalski 2021** — a viscous-fluid cytoplasm gives broader, faster blebs; a poroelastic one
  gives consistent morphology, because the elastic response relieves pressure and limits bleb size.
  **Used for**: naming, in §14.3, the observable this entry does *not* reach.

No retraction is known against any of the three. **None of them is evidence that this module's
discretisation is correct**, and none is cited as such.

## 7. Independent oracle or derivation

The oracle is (5)–(6), derived in §4 from the assembled operators and **not** from the code that
computes `D_p`. The two are independent in the way that matters: the property is four arithmetic
operations on card fields, and the check is a time-stepped solve of a 64-cell coupled saddle-point
system through `CytosolField.step`. A mistake in the algebra of (6) cannot also be present in
`np.linalg.solve`.

The measurement, for a zero-mean cosine mode `p ∝ cos(π x / L)` on the coupling axis, uses the
**discrete** eigenvalue of the assembled no-flux Laplacian rather than the continuum `k²`:

```
λ = ( 2 - 2 cos(π dx / L) ) / dx²
```

and backward Euler gives the exact per-step ratio `p^{n+1} = p^n / (1 + D_p λ dt)`, so `D_p` is
recovered as `((ratio)^(-1/steps) - 1) / (λ dt)` with **no discretisation error at all** — the mode
is an exact eigenvector of the discrete operator, so the recovered number is the solver's diffusivity
to round-off and not to truncation. Using `k²` instead would have introduced an `O(dx²)` error and
made the test a convergence study wearing an identity's clothes.

**One transient is real and is not error.** The skeleton is quasi-static, so a run started from
`u = 0` with a pressure already imposed is off the manifold (4), and the *first* step alone carries
an instantaneous drop of `S/(S + α²/H)` as the skeleton takes up its share. At the §8 card that is a
factor of two, and measuring across it recovers `D_p` values of order 100 rather than 2. It is
skipped by warming up five steps before the measurement window opens, and it is documented here
because a later reader who removes the warm-up will see a wrong answer with no error message.

## 8. Positive control

`tests/vertical/test_cytosol_poroelastic_diffusivity.py` — **15 test functions, 27 cases**, of which
12 functions are positive controls and 3 are the negative controls of §9.

* `test_the_measured_decay_rate_is_the_reported_poroelastic_diffusivity` — the load-bearing one.
  Five cards spanning `α ∈ [0.5, 1.0]` and `H ∈ [0.8, 300]` pN/µm², measured through `step` against
  §7's oracle. Gate `1e-11`; **observed worst case `2.479e-13`**, per-card below.
* `test_the_measurement_would_reject_the_decoupled_formula` — the same five, against `mob/S`.
  Asserts the measurement *separates* them, so a card too stiff to discriminate fails rather than
  passes quietly.
* `test_the_reported_diffusivity_matches_the_hand_derivation` — (6) written out independently, on
  all six cards.
* `test_the_skeleton_strictly_slows_the_pressure_down` — `D_p < c_v` strictly on every card.
* `test_the_sharp_card_halves_the_diffusivity_exactly` — `α²/H == S` exactly, so `D_p == c_v/2`
  exactly. Asserted with `==`; the arithmetic is exact in binary and a tolerance would hide a
  changed card.
* `test_the_literature_configuration_lands_in_the_reported_band` /
  `test_the_decoupled_reading_would_have_been_an_order_of_magnitude_fast` /
  `test_the_mobility_is_derived_from_pore_geometry_and_not_written_down` — §10's configuration,
  its 11× decoupled contrast, and the `E ξ²/μ` reduction with its 1.10 overestimate factor.
* `test_the_four_load_bearing_invariants_survive_this_entry` — the brief's four, re-asserted on this
  entry's two cards with `==` where exactness is claimed, and each guarded against being vacuous
  (a boundary-face count, and a non-uniform field for the divergence identity).
* `test_there_is_still_no_default_poroelastic_card`,
  `test_the_mobility_helper_refuses_a_missing_argument`,
  `test_the_mobility_helper_refuses_a_degenerate_value` — `ALEPH-PORT-1803`, still holding.

Measured, this session:

| `α` | `H` [pN/µm²] | `c_v` | `D_p` reported | `D_p` measured | rel. err |
|---|---|---|---|---|---|
| 1.0 | 2.0 | 4.000000 | 2.0000000000 | 2.0000000000 | `1.13e-14` |
| 0.8 | 300.0 | 4.000000 | 3.9830058417 | 3.9830058417 | `2.21e-13` |
| 0.9 | 3.0 | 4.000000 | 2.5974025974 | 2.5974025974 | `2.48e-13` |
| 1.0 | 1.0 | 4.000000 | 1.3333333333 | 1.3333333333 | `1.15e-14` |
| 0.5 | 0.8 | 4.000000 | 2.4615384615 | 2.4615384615 | `1.47e-13` |

## 9. Deliberately failing negative control

The brief names the failure precisely: *"an elastic term that is computed and never reaches the
pressure — every existing test passes, the model is unchanged, and `D_p` is a number in a
docstring."* That mutant is **constructible in this tree without editing a line**, because
`decouple_solid=True` is exactly a skeleton that is computed and never reaches the pressure.

`test_a_decoupled_skeleton_makes_the_cards_diffusivity_a_lie` drives it:

| Quantity, on a field with `decouple_solid=True` | Value | Why |
|---|---|---|
| `card.poroelastic_diffusivity_um2_per_s` | `2.0` | the card does not know about the break — **this is the number in the docstring** |
| `field.poroelastic_diffusivity_um2_per_s` | `4.0` | routes through `effective_biot_coefficient`, which is `0.0` |
| **measured** decay rate | `4.0` | what the solver actually does |

The control asserts the measurement agrees with the **field** and disagrees with the **card** by a
factor of two. It fails if `D_p` is put only on the card; it fails if the field property reads
`card.biot_coefficient` instead of `effective_biot_coefficient`; and it fails if the elastic term is
ever detached from the assembly while the reported number goes on including it. That is the
"computed but never reaches the pressure" mutant, planted and killed.

`test_the_measurement_would_reject_the_decoupled_formula` is the second half, and it exists because
the first is not sufficient. It asserts that the measured rate on a **coupled** field is *not*
`mob/S` — separated by a factor of two on the sharp card — so a `D_p` that silently dropped the
`α²/H` term would be caught even though such a `D_p` is a perfectly well-typed float that every
`ALEPH-PORT-2401` control still passes.

**Both controls can fail, and were seen to fail.** Before the field property was routed through
`effective_biot_coefficient` the first control reported `2.0` against a measured `4.0`; the run is
in §13.

## 10. Numerical and precision envelope

All float64. No GPU, no float32, no device residency (§11).

**The configuration `D_p` is reported at**, and every number in it lives in the test module and in
this file, never in `aleph/**`:

| Parameter | Value | Unit | Where it comes from |
|---|---|---|---|
| pore size `ξ` | `0.02` | µm (20 nm) | Moeendarbary's stated cytoskeletal mesh range, 10–100 nm |
| cytosol viscosity `μ` | `5.0e-3` | pN·s/µm² (5 mPa·s) | a few times water, as that paper's discussion assumes |
| pore shape factor `a` | `1.0` | — | `k = ξ²`, the scaling's own convention, **required not defaulted** |
| ⇒ mobility `k/μ` | `0.08` | µm⁴/(pN·s) | `darcy_mobility_um4_per_pn_s`, computed |
| oedometric modulus `H` | `1.0e3` | pN/µm² (1 kPa) | cytoplasmic elastic modulus, order 1 kPa |
| Biot `α` | `1.0` | — | incompressible constituents |
| storativity `S` | `1.0e-4` | µm²/pN | `M = 10 kPa` |
| **⇒ `D_p`** | **`72.7`** | **µm²/s** | equation (6) |
| for contrast, `c_v = mob/S` | `800.0` | µm²/s | the decoupled reading — **11× too fast** |

`D_p = 72.7 µm²/s` sits inside Moeendarbary's reported 1–100 µm²/s. **This is a plausibility check
and not a validation**, for three reasons stated here so no later reader has to reconstruct them:
the parameters were chosen to land in the band rather than measured; `D_p ∝ ξ²`, so the 10–100 nm
mesh range alone spans a factor of 100 and covers the entire reported band by itself; and the
skeleton term dominates the denominator here by 10:1, which is what makes `E` the controlling
parameter — exactly as the literature says — but is a property of the chosen `S`.

The decoupled number is the point of the entry: at this configuration the module would have
reported `800 µm²/s`, an order of magnitude fast and outside the literature band, and nothing in the
tree would have contradicted it.

**Precision.** The decay measurement recovers `D_p` to a relative error of `≤ 2.5e-13` across five
cards — round-off of a dense `126 × 126` solve repeated 45 times, not truncation, because the mode
is an exact eigenvector (§7). The gate is set at `1e-11`, roughly forty times the observed worst
case: tight enough that the factor-of-two mutant in §9 is rejected by eleven orders of magnitude,
loose enough not to be a BLAS-version tripwire. The exact claims — the `== c_v` collapse under the
break, the `== 0.0` boundary flux, the `== 0.0` divergence residual — are asserted with `==` and
carry no tolerance at all.

## 11. Production-backend residency and transfer

**None, and deliberately.** This entry adds no array, no kernel and no state: `D_p` is a scalar
property of five floats already resident wherever the card is, and it is computed on the host at the
moment it is asked for. There is nothing to transfer to a device and no `warp` kernel is written.
`CytosolField`'s dense assembly is CPU-only in the landed module and this entry does not change
that. No GPU run was made and none was requested — `CLAUDE.md` §3 makes asking a standing
requirement, and there was nothing here to ask about.

## 12. Comments and docstrings to discard

None to discard: no provider file was opened, so there was no provider prose to strip. All new
prose is written for this entry. Two conventions of the surrounding module are followed because they
are load-bearing rather than stylistic — the reported quantity names the thing it is *not*
(`consolidation_coefficient_um2_per_s` already warns it is the decoupled reading; the new
`poroelastic_diffusivity_um2_per_s` says which one it is and what breaks it), and every deliberate
break is a flag on the shipped object rather than a sabotaged copy.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending PI review.** The code and controls are written and green; the status stays `PROPOSED` because the physical claim in §10 is a plausibility check against literature and the PI's decision `#69f82385` asked for the poroelastic model to be *completed*, which §14 says this does not do. |
| Controls | `tests/vertical/test_cytosol_poroelastic_diffusivity.py` — 15 test functions / 27 cases, all green. `tests/vertical/test_cytosol_controls.py` — `ALEPH-PORT-2401`'s 47 controls, re-run **unedited**, all green. Directly affected neighbours (`test_package_exports`, `test_connectors_fluid`, `test_connectors_transfer`, `test_owner_state_keys_match_contracts`): **164 passed** together with the two above. `tests/firewall/` 26 passed, `tests/ports/` 21 passed after `INDEX.md` was regenerated. |
| Negative controls seen to fail | Yes — **two mutants, both planted in the shipped file and both killed, then reverted and the file diffed back to byte-identical.** *Mutant 1*, `CytosolField.poroelastic_diffusivity_um2_per_s` reading `card.biot_coefficient` instead of `effective_biot_coefficient`: 2 failed / 25 passed, reporting `D_p = 2 um^2/s but measurably diffuses at 4 um^2/s`. *Mutant 2*, the `α²/H` term zeroed so the modulus never reaches the reported number — the brief's "computed and never reaches the pressure": 6 failed / 21 passed, **including the `H = 300` card where the effect is only 0.4 %**, which is the evidence that the gate discriminates in the least favourable regime and not only on the sharp card. |
| Reviewer | Unassigned — no agent may sign this off. |
| Rollback | Delete the four appended members, the three `__all__` entries and the new test module. No existing line of `aleph/vertical/cytosol.py` is modified, so rollback cannot disturb `ALEPH-PORT-2401`. |

## 14. Honest limits

**14.1 This does not complete the poroelastic model, and the brief's own framing is why.** The PI's
decision `#69f82385` selected *"poroelastic 완성 — 세포질 골격에 탄성 추가"*. The elasticity was
already there; what this entry adds is the **number that says how fast the coupled system relaxes**,
which is the observable the literature reports and the one the tree could not previously state. That
is a real gap closed and it is smaller than "completing the model".

**14.2 The `72.7 µm²/s` is not a measurement of cytoplasm.** §10 lists the three reasons. The
strongest of them: `D_p ∝ ξ²` and the pore size is known only to a factor of ten, so the model's own
uncertainty on this number is a factor of a hundred — wider than the literature band it agrees with.
Agreement of that kind is a constraint that the model is not absurd, and nothing more.

**14.3 The uniaxial-strain restriction still stands and it bounds `D_p`'s meaning.** The skeleton
deforms along one declared axis with no shear and no deviatoric stress (`ALEPH-PORT-2401` §14).
Equation (3) uses the null space of `Dᵀ` on that axis, so `D_p` as derived here is the **confined**
poroelastic diffusivity. Moeendarbary's microindentation is not a confined geometry; the two numbers
are the same order and are not the same quantity. **Strychalski's bleb result is therefore out of
reach**: relieving pressure to limit a bleb needs the unconfined response and a moving boundary, and
this owner has neither.

**14.4 Zero external solid load is an assumption of the derivation, not of the code.** (3) holds
when `f_ext = 0`. `assemble_coupled_system` accepts `external_solid_load_pn`, and under a non-zero
one `H θ - α p` is no longer constant along the axis and the single-mode decay rate is no longer
`D_p λ`. The property returns the same number regardless. **A caller applying a solid load and
reading `D_p` as a time constant is reading a number that does not apply**, and nothing in the code
refuses it. This is the sharpest limitation in this entry.

**14.5 `UNVERIFIED` on a device.** Nothing here has been executed anywhere but this host's CPU.

**14.6 The brief asked for a solid phase and got a diffusivity, which is a smaller deliverable.**
Recorded plainly rather than dressed up: the premise check in §0 consumed the part of the lane that
would have built new physics, and what remained was genuinely absent but narrow. A reviewer who
wants the larger thing should read §0 first and decide whether the premise or the scope was wrong.
