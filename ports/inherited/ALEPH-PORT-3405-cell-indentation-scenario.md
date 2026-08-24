# ALEPH-PORT-3405 — pressing the cell: a spherical indenter and the force–indentation curve

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3405` |
| Lane | `f70de564` scenario lane (indentation sub-lane) |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Amended | `2026-07-31` — **§16 written before the release code**, same rule. See §16. |
| Port class | `RE-DERIVED` |
| Exists because | `CellWorld` can now compose owners and step them, and **nothing has ever been done to the cell**. A world that relaxes to its own equilibrium and is never disturbed produces no measurement. |

> **Read §16 before §8b.** §8b–§8c record the *first* encoding, whose contact set came from the
> reference configuration and could not release. §16 replaces that contact law with a
> current-configuration unilateral (Signorini) contact and re-measures the curve. The tables in
> §8b/§8c are kept, not deleted: they are the measurement that motivated §16 and the negative
> control §16 is checked against.

---

## 1. Aleph API

```python
from aleph.scenarios.indent import (
    ContactWeld,                  # which material points the indenter holds, and where
    SphericalIndenter,            # a rigid sphere, expressed as a boundary condition
    IndentedWorld,                # CellWorld + the constraint, satisfying the relax protocol
    IndentationPoint,             # one depth: reaction, counts, convergence, and its own trace
    IndentationCurve,             # the whole F(d) sweep, its fits, and the stopping-point study
    build_indentation_world,      # membrane + cortex + turgor + an indenter at one depth
    measure_indentation_curve,    # walk the ladder, relax at each rung, record the reaction
    fit_power_law,                # log F against log d, returning the exponent
    render_world_png,             # a deformed-mesh image, numpy + zlib only
    render_indentation_curve_svg, # the F(d) curve, as text
)
```

**Files written by this entry:** `aleph/scenarios/indent.py`, `tests/scenarios/test_indent.py`,
this ledger, and artefacts under `runs/indentation/`. **Nothing else is modified** — in particular
not `aleph/scenarios/world.py`, not `aleph/vertical/**`, not `aleph/state/**`, not `wiring.py`,
and not any file under `tests/firewall/`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | not resolved — no file was opened |
| Source path / symbol | **none named, none read** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened by this lane. |

The physics here is textbook contact mechanics and a Dirichlet boundary condition on an explicit
descent. There is nothing to port, and nothing was.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** The encoding is forced by what already exists:
`aleph/runtime/pipeline.py` runs `ACCUMULATE_FIELD` immediately before `SOLVE`, and
`HelfrichMembrane.accumulate(SOLVE)` takes `x ← x + dt·µ·F` from its own accumulator. A participant
scheduled in `ACCUMULATE_FIELD` that zeroes selected rows of that accumulator is *exactly* a
prescribed-displacement constraint on those rows, with no change to any owner.

## 4. What this is — and the census claim it deliberately does not make

### 4.1 The indenter is a boundary condition, not a component

**There is no registered `indenter` in `aleph/state/**`, and this entry does not add one.** A cell
census enumerates the cell; an AFM tip is not part of the cell. Inventing a census component for the
apparatus would be a census change this lane may not make, and it would be wrong on the merits.

So `SphericalIndenter` is registered into the `Pipeline` as a **participant that owns its own
constraint bookkeeping and nothing else**. Its `owned_entities()` names `indenter_contact_sites` —
the set of prescribed sites and the reaction recorded at them. It does **not** claim
`membrane_nodes`, which the membrane owns; `Pipeline._check_work` would raise `FOREIGN_ENTITY` if it
did, and that guard is right.

### 4.1a The cell has to be held down, and finding that out cost the first sweep

**An unsupported cell is a free body, and rigid translation is a zero-energy mode of every term in
this world** — Helfrich, cortex springs, both connectors, the volume penalty. So an indenter welded
to a free cell does not indent it: the cell moves out of the way, at no energy cost and with no
restoring force. The entry as written did not anticipate this and the first production sweep
measured it instead of the physics.

Measured, level 2, `R = 3 µm`, `d = 0.4 µm`, one welded vertex, no support:

| accepted steps | free residual (pN) | axial F (pN) | centroid z (µm) | `R_eq` (µm) |
|---|---|---|---|---|
| 250 | 3.2895e+00 | 34.70418 | −0.005612 | 4.960929 |
| 3,000 | 4.0265e-01 | 25.47545 | −0.031463 | 4.960438 |
| 18,000 | 1.7226e-01 | 17.10333 | −0.146273 | 4.959868 |
| 78,000 | 2.3705e-02 | 3.68662 | **−0.345007** | 4.959437 |

The centroid marches toward the prescribed 0.4 µm while the equivalent radius does not move: the
cell is **translating, not deforming**, and the reaction decays toward zero as it escapes. Two
converged points of the first sweep said the same from the other side — `F(0.05) = 0.186449` pN
against `F(0.10) = 0.183493` pN, a force that does not grow with depth, at a converged energy
indistinguishable from the free cell's `9301.989723` pN·µm.

**Nothing raised. The curve looked like a curve.** It is the same class of failure as the two in
§9a: a plausible number from a configuration that does not mean what it appears to mean.

The fix is the one the experiment already contains: an AFM presses an *adherent* cell against a
dish. `FlatSupport` welds a bonded footprint on the underside — a boundary condition exactly as the
indenter is, not a census component. It is **on by default** and `support_radius_um=None` is an
explicit opt-out, because the unsupported world is only ever wanted by the negative control.

With it, the same configuration is well posed. Measured, `d = 0.4 µm`, footprint radius 2.0 µm
(7 vertices):

| accepted steps | free residual (pN) | indenter F (pN) | dish F (pN) | centroid z (µm) |
|---|---|---|---|---|
| 3,000 | 1.2891e-01 | 25.68188 | 5.00844 | −0.027917 |
| 18,000 | 8.6567e-02 | 20.77689 | 16.03797 | −0.088720 |
| 78,000 | 1.7084e-03 | **19.45105** | **19.41237** | −0.109866 |

The reaction **settles** instead of decaying, the centroid saturates at −0.110 µm instead of running
to −0.4, and the indenter and the dish agree to **0.2%** — an independent check, since the two are
summed over disjoint vertex sets by different objects.

Two controls in the suite would have caught the original defect instantly and cost no relaxation at
all: `test_translating_the_whole_cell_costs_exactly_no_energy` asserts the zero mode directly, and
`test_the_escaped_configuration_satisfies_the_weld_at_the_free_energy` shows that the reference cell
translated down by `d` *exactly satisfies the weld at the free energy* — so the unsupported problem's
minimum **is** the escaped state, with a reaction of zero. That is a statement about the problem
rather than about the solver, which is why it needs no steps to check and why it should have been
written first.

### 4.2 The encoding, stated plainly

A rigid sphere of radius `R` centred at `c = (0, 0, z_c)` presses down the `+z` axis into a relaxed
cell. At depth `d`:

1. `z_c(d) = z_c0 − d`, where `z_c0 = max_i { z_i + sqrt(R² − r_i²) : r_i < R }` over the **relaxed
   reference** membrane vertices, `r_i = sqrt(x_i² + y_i²)`. This is the mesh-exact first-touch
   height. It is **computed, not assumed**: the vertical is built at 5.0 µm and relaxes to
   `R_eq = 4.959428` µm at level 2 and `4.995156` µm at level 3 (measured this session), so
   `z_c0 = 5 + R` would offset every depth in the sweep by ~0.04 µm — a systematic error in `d` that
   a log–log fit converts directly into a wrong exponent.
2. A vertex enters the **welded set** at the depth where the descending sphere first engulfs it:
   `{ i : |x_i^ref − c(d)| < R }`, evaluated on the reference configuration. The set only grows.
3. On entry it is placed on the surface by closest-point projection, and **from that depth onward it
   translates rigidly with the indenter**, `x_i(d) = x_i(d_entry) − (d − d_entry) ẑ`. It never
   slides.
4. The weld is imposed by zeroing rows `i` of `membrane.forces` in `ACCUMULATE_FIELD`, via
   `backend.scatter_add(forces, idx, −forces[idx])`. That is exact in float64, so `SOLVE`'s
   `x ← x + dt·µ·F` adds exactly `0.0` and the welded vertices are bit-identical across a whole
   relaxation — asserted with `np.array_equal`, not a tolerance. It also advances the witness
   counter, so the constraint cannot be a coverage no-op.
5. **The cortex is carried with the membrane** at the welded indices, by the identical displacement
   vector. Without this the encoding does not run at all — see §9.

The **reaction** the indenter applies is `R_i = −F_i^cell`, which is what step 4 removes. By Newton's
third law the force the cell applies to the indenter is `+F_i^cell`, and the measured axial force is

```
F(d) = Σ_{i ∈ weld} F_i^cell · ẑ        [pN],  ẑ = +z
```

### 4.2a Welded rather than frictionless — a measurement, not a preference

**The frictionless version was written first, and it does not have a well-posed equilibrium in this
cell model.** The Helfrich energy depends on area and curvature and so barely resists a vertex
sliding *within* the surface; the only membrane↔cortex coupling is `CortexMembraneContact`, a
unilateral **normal** law, which resists approach and not tangential slip. A contact vertex on a
frictionless sphere is therefore very nearly tangentially neutral, and it drifts.

Measured at level 2, one contact vertex, depth 0.4 µm, over 28,000 accepted steps:

| steps | free residual (pN) | reaction (pN) | tangential residual (pN) |
|---|---|---|---|
| 500 | 1.0236e+00 | 31.60246 | 5.663e-14 |
| 6,000 | 1.7517e-01 | 23.17674 | 9.296e-12 |
| 13,000 | 1.0402e-01 | 19.35652 | 5.920e-09 |
| 23,000 | 9.0372e-02 | 15.05679 | 7.788e-06 |
| 28,000 | 8.3613e-02 | 13.27147 | **1.309e-04** |

The tangential residual grows by nine orders of magnitude while the reaction decays monotonically
with no sign of settling. That is a vertex sliding off the indenter, not a cell reaching
equilibrium.

Welding fixes it and makes the oracle exact rather than asymptotic: a welded vertex has
`dx_i/dd = −ẑ`, so `dU/dd = Σ F_i·ẑ` holds with **no** tangential term and no normal-projection
approximation. Under the frictionless prescription `dx_i/dd` had a tangential part and the identity
held only in the limit of a converged tangential residual — precisely the thing that was not
converging.

### 4.3 What this encoding is honestly not

- **Bonded, and Hertz is frictionless.** The weld carries tangential load, which a frictionless
  sphere cannot. `SphericalIndenter.friction_fraction` reports at every depth what fraction of the
  held load is tangential, so how far the measurement has travelled from the closed form it is
  compared against is a number in the artefact rather than an argument in prose.
- **The result depends on the depth ladder**, because a vertex's weld position is fixed at the depth
  it entered. A finer ladder converges to the continuous limit; how fine is enough is a measurement,
  not an assumption, and it is **not established here**.
- **The welded set only grows and cannot release.** A retracting ladder is refused rather than run:
  it would report the adhesion of a glued sphere. `released_site_count` counts welded vertices whose
  reaction has gone adhesive, i.e. how hard the no-release assumption is being leaned on, and
  `intruder_count` counts unwelded vertices that have drifted inside the sphere. A point with either
  non-zero is labelled out of regime rather than dropped.
- **Quasi-static.** `relax_to_equilibrium` is a descent, not dynamics (`aleph/vertical/relax.py`
  module docstring). Nothing here is evidence about rate dependence, viscoelasticity, or an AFM
  approach velocity.

### 4.4 The cytosol is deliberately absent, and that is a finding not an omission

The fluid lane established tonight that `aleph/vertical/cytosol.py`'s skeleton is **uniaxial-strain
with zero shear degrees of freedom**, and that its `max_degrees_of_freedom = 8192` with a dense
`np.linalg.solve` caps a 3-D grid at 16³, i.e. `dx ≈ 1.25 µm` — coarser than the nucleus and far
coarser than the contact patch measured here. A cell being indented is **not** confined, so a
uniaxial-strain skeleton is the wrong kinematics for it. Adding the cytosol would produce a
poroelastic relaxation curve that could not mean what it appeared to mean.

**So this scenario is elastic, not poroelastic, and says so.** Recovering an indentation
poroelasticity needs the solid to gain shear degrees of freedom first, which is a lane, not a
parameter.

## 5. Units, domains, singular cases, invariants

Length µm, force pN, energy pN·µm. Depth `d` in µm, indenter radius `R` in µm.

Singular cases, each **refused** rather than absorbed:

| Case | Why refusing is right |
|---|---|
| `R ≤ 0` or non-finite | A sphere with no radius has no surface to prescribe onto. |
| `d < 0`, or a ladder rung below the previous one | A welded contact cannot release, so a retraction would report the adhesion of a glued sphere. Retraction is a different experiment needing a law that can unbind. |
| a depth ladder that is not monotonically increasing | Same reason, caught one level up in the sweep. |
| welded set empty | An indenter touching nothing measures nothing; `F = 0` would look like a measurement of a very soft cell. |
| welded set smaller than `min_contact_sites` | See §5a — this one **moved** during implementation. |
| welded set is every membrane vertex | The indenter has engulfed the cell; no free degrees of freedom remain, so the relaxation "converges" instantly with a residual of zero and reports the prescription. |
| a vertex coinciding with the indenter centre | It has no closest point on the surface and no contact normal. |
| a vertex named twice in the welded set | It would be double-counted in the reaction. |
| any phase but `ACCUMULATE_FIELD` | Earlier, the constraint would react to a force not yet assembled; later, it would constrain a configuration that has already moved. Both are silent. |
| a fit with fewer than 3 points, or any non-positive force | Two points fit any exponent exactly and report `r² = 1`. |
| a matched-residual refit at a residual some depth never reached | Falling back to the nearest sample would mix stopping rules inside one fit, silently. |

### 5a. One refusal moved, and the entry is corrected rather than left wrong

This entry originally put the resolution guard on **world construction**: fewer than 3 contact
vertices refused the world. Writing the code showed that to be the wrong place, for a reason the
design did not anticipate. A welded contact needs a **ladder** — a vertex's held position is fixed
at the depth it enters — and the shallow rungs of that ladder necessarily have one or two contact
vertices. Refusing them at construction would force every depth to be reached in one cold jump,
which is the thing the ladder exists to avoid.

So the guard now lives on **fit admission**: every rung is computed and recorded, `IndentationPoint`
carries `resolved`, and `IndentationCurve.usable()` excludes the unresolved ones from every fit.
`build_indentation_world` keeps the hard refusal for a single-depth call. The distinction is real:
an under-resolved contact is a **valid mechanical calculation** and an **invalid measurement**, and
those deserve different refusals. Recording them also keeps the artefact honest about what was
computed instead of silently dropping it.

Invariants, each with the control that asserts it:

- **I1.** Coverage is complete after an indented relaxation — the indenter witnessed non-zero work —
  and incomplete before anything has run, so `ok` cannot be a constant.
- **I2.** The welded vertices **do not move**, `np.array_equal`, across an entire relaxation.
- **I3.** The reaction assembled out of the pipeline (`IndentedWorld.contact_force_pn`) equals the
  reaction recorded inside it by the indenter's own `accumulate`, at the identical configuration.
  Two independent paths pinned against each other, exactly as `VerticalWorld` pins its two.
- **I4.** **Virtual work**, `F = dU/dd` under a rigid translation of the welded set — see §7.1 for
  what this does and does not establish. The load-bearing oracle of this entry.
- **I5.** A welded vertex translates purely along `−z` between rungs (`array_equal` on both
  transverse components) and stays exactly on the translating sphere.
- **I6.** The convergence measure **excludes the welded rows**; the broken flag re-admits exactly
  the reaction and nothing else.
- **I7.** The composed world's energy and forces equal a plain `VerticalWorld`'s at the identical
  configuration, so the constraint added no energy term and did not leak into `evaluate_forces`.

## 6. Source evidence class and known retractions

No claim about any reference implementation; none was read. Every magnitude is `UNSOURCED`. Evidence
rung `ANALYTIC_ORACLE`; quantitative status `BLOCKED` — see §14.

## 7. Independent oracle or derivation

**Three, in decreasing order of how much they establish.**

### 7.1 Virtual work — the exact one (I4), and the version of it that failed

**What is asserted.** Translate the welded set rigidly by `s ẑ` with **every free degree of freedom
held fixed**. Then

```
dU/ds = Σ_i (∂U/∂x_i) · ẑ = −Σ_i F_i^cell · ẑ,   and  d = −s,  so  dU/dd = F(d).
```

The free nodes' contribution is exactly zero **by construction** — they did not move — rather than by
equilibrium. So this is a gradient identity, it holds at *any* configuration converged or not, and it
contains no material constant and no geometric assumption. It checks the sign, the reaction sum, the
rigid-translation property of the weld, and that `evaluate_forces` really is minus the gradient of
the energy the descent minimises.

**Measured**, level 2, depth 0.4 µm, central differences at `(4e-4, 2e-4, 1e-4)` µm:

| step (µm) | relative error |
|---|---|
| 4e-4 | 1.0534931e-03 |
| 2e-4 | 2.6323069e-04 |
| 1e-4 | 6.5799305e-05 |

**Observed order 2.0010** — clean second-order convergence, i.e. the whole remaining disagreement is
the central difference's own truncation error. The control therefore asserts the **order**, not a
tight absolute bound; a tight absolute bound would only be asserting the step size.

**The version that failed, recorded because the failure is informative.** The first form of this
oracle differenced the *relaxed* energy across three depths, `[U(d+h) − U(d−h)]/2h` with each rung
relaxed. At level 2, `d = 0.40 ± 0.02` µm, 4,000 steps per rung, that gave **−16.4 pN** against a
reaction of **+22.1 pN**: wrong by 174% and wrong in sign. Nothing is wrong with the identity — each
deeper rung had simply had more cumulative descent than the one before, so the difference quotient
measured the relaxation drift rather than the depth derivative. **The relaxed statement
`dU_equilibrium/dd = F` is therefore `UNVERIFIED`**, and it stays so until a descent that converges
exists. See §14.

### 7.2 Hertz — the *exponent*, which is the part with no constant in it

For a rigid sphere of radius `R` indenting an elastic half-space by `d`:

```
F = (4/3) E* sqrt(R) d^(3/2)
```

`E*` is unknown for this composite and is not fitted. **The 3/2 is the oracle**, because it survives
every choice of modulus. Measured by ordinary least squares of `log F` on `log d`.

**Hertz is expected to fail here, and reporting where is the result.** A thin pressurised shell under
tension is not an elastic half-space:

- **Small `d`, tension-dominated.** A membrane at tension `σ` resists a point-like indentation with
  `F ∝ d` (up to a logarithm in `d`), not `d^(3/2)`. With `σ = 30 pN/µm` the tension force scale is
  `2πσ ≈ 188 pN/µm` of depth. An exponent near **1** here is the correct answer, not a failure.
- **Large `d`,** the shell wraps the indenter, the cortex is compressed against it and the turgor
  compartment's bulk term engages. Departure upward from any single power law is expected.

**A measured crossover between the two regimes is a better result than a fitted 1.5**, and the
protocol is to fit windowed exponents across the sweep and report the sequence rather than one
number.

### 7.3 The composition oracle it inherits

`build_indentation_world` at zero prescribed sites must reduce **exactly** to the `CellWorld` that
`ALEPH-PORT-3404` §8 showed reproduces `build_vertical` to every digit. So the indented world is
anchored to the project's one verified physics result through a chain of exact reductions, rather
than to a fixture.

## 8. Positive control

**Measured this session on this tree.** Every number below came from a run in this session.

```text
$ /Users/sw1/miniconda3/envs/aleph/bin/python -m pytest tests/scenarios/test_indent.py \
      -p no:cacheprovider
..................................................................       [100%]
66 passed in 39.65s
```

The whole directory, `tests/scenarios/`, gives **104 passed, 1 failed** — the one failure is
`tests/scenarios/test_spread.py::TestTheObservablesMeasureWhatTheyName::test_the_centre_of_mass_is_area_weighted_and_not_a_node_mean`,
which belongs to the **spreading** lane and was mid-write. Reported rather than fixed, per
`CLAUDE.md` §1.

| Control | Asserts | Measured |
|---|---|---|
| `TestVirtualWork::test_the_axial_force_is_minus_the_energy_gradient_of_a_rigid_indenter_translation` | I4, §7.1 | errors `1.0535e-03 / 2.6323e-04 / 6.5799e-05` at steps `4e-4 / 2e-4 / 1e-4` µm, **observed order 2.0010** |
| `TestTheConstraintHolds::test_the_welded_vertices_never_move` | I2 | `np.array_equal` over 400 steps |
| `TestTheConstraintHolds::test_the_reaction_agrees_down_two_independent_paths` | I3 | in-pipeline vs out-of-pipeline within `1e-11 ×` the constituent scale |
| `TestTheCompositionIsUnchangedByTheIndenter::test_the_energy_is_the_vertical_s_energy_at_the_same_configuration` | I7 | within 8 ULP (fsum vs `+`) |
| `TestTheCompositionIsUnchangedByTheIndenter::test_the_forces_are_the_vertical_s_forces_at_the_same_configuration` | I7 | `np.array_equal` on every owner |
| `TestTheWeldIsTheLoadingHistory::*` | I5 | `array_equal` on both transverse components; surface error `< 1e-12` µm at every rung |
| `TestTheIndenterIsNotACoverageNoOp::*` | I1, §4.1 | verdict `ok` after stepping, not `ok` before; census names `indenter_contact_sites` only |
| `TestThePowerLawFit::test_it_recovers_an_exactly_generated_exponent` | the fitter | exponent recovered to `< 1e-12` for 0.5, 1.0, 1.5, 2.0 |

### 8a. Mutation testing (ledger §4a)

Four mutants applied to the shipped module, controls run against each, then restored and verified
**byte-identically** — `sha256 a533e9b67d9405f3a2a779423e915c9422cf4bb0d4fe41ed95e80141452b642f`
before and after. Baseline: `77 passed in 46.52s`.

| # | Mutant | Tests killed | Survived? |
|---|---|---|---|
| M1 | `first_touch_centre_um` returns the naive `5.0 + R` instead of the measured touch height | 2 | no |
| M2 | a welded vertex does not translate with the indenter (`displacement[:, 2] = 0.0`) | 2 | no |
| M3 | the weld removes only the **axial** component of the force, not the whole vector | 1 | no |
| M4 | the substrate records its reaction but never removes it | 1 | no |

**M3 is the one worth reading.** The transverse leak it introduces is `1.06e-15` µm over 400 steps —
below any tolerance a reasonable person would write, and invisible to `np.allclose` at every default
setting. It is caught only because I2 is asserted with `np.array_equal`. That is the concrete case
for the bit-identical assertions in this suite: the failure mode of a constraint is drift, and drift
starts at the last bit.

### 8b. The measurement — the force–indentation curve

Supported cell, level 2 (162 vertices), indenter `R = 3.0 µm`, footprint radius 2.0 µm (7 vertices),
welded ladder, **every point relaxed to the default `1e-03` pN free-residual tolerance**.

| d (µm) | F (pN) | weld | adh | F/d (pN/µm) | dish (pN) | balance | steps | resolved |
|---|---|---|---|---|---|---|---|---|
| 0.05 | 2.492640 | 1 | 0 | 49.853 | −2.391652 | 4.0% | 32,085 | no |
| 0.10 | 4.935279 | 1 | 0 | 49.353 | −4.873149 | 1.3% | 37,780 | no |
| 0.20 | 9.808386 | 1 | 0 | 49.042 | −9.772976 | 0.36% | 51,596 | no |
| 0.30 | 14.650070 | 1 | 0 | 48.834 | −14.613890 | 0.25% | 51,446 | no |
| 0.40 | 19.450740 | 1 | 0 | 48.627 | −19.412890 | 0.19% | 50,988 | no |
| 0.50 | 24.210370 | 1 | 0 | 48.421 | −24.170630 | 0.16% | 50,572 | no |
| 0.60 | 16.263890 | 3 | **2** | — | −16.280810 | 0.10% | 55,919 | yes |
| 0.70 | 13.847700 | 7 | **6** | — | −13.853900 | 0.04% | 51,823 | yes |
| 0.80 | 22.587230 | 7 | **6** | — | −22.569540 | 0.08% | 51,577 | yes |

The ladder continues to `d = 1.20` µm; those rungs are all `weld = 7, adh = 6`, i.e. out of regime
for the reason in §8c, and they do not change anything above. The complete record — every point,
every convergence trace, and a rendered frame per depth — is written by the run itself to
`runs/indentation/` (`sweep_supported_level2.log`, `indentation_supported_level2_R3.{json,csv,svg}`,
`indentation_supported_level2_R3_frames/`). `runs/` is gitignored, which is this repository's
existing convention for measured output.

`balance` is the disagreement between the indenter reaction and the dish reaction, summed over
**disjoint** vertex sets — an independent check that tightens monotonically as the descent converges,
from 4.0% to 0.04%.

**The fitted exponent, over the six converged in-regime points spanning a 10× range in depth:**

```
exponent = 0.988049      r^2 = 0.99999858      prefactor = 48.079 pN/um^n
```

Windowed, to show it is not an average across a crossover:

| window (µm) | exponent | r² |
|---|---|---|
| 0.05–0.20 | 0.98817 | 0.9999975 |
| 0.10–0.30 | 0.99043 | 0.9999999 |
| 0.20–0.40 | 0.98785 | 0.9999986 |
| 0.30–0.50 | 0.98347 | 0.9999985 |

**This is 1, not 3/2. The cell is in the tension-dominated regime and Hertz does not apply here.**
`F/d` drifts only from 49.85 to 48.42 pN/µm across an order of magnitude in depth — a 2.9% softening,
in the direction a finite pressurised shell should soften. Nothing in the measured range resembles
`d^1.5`; a Hertzian branch would have shown `F/d` rising as `sqrt(d)`, i.e. by a factor of 3.2 over
this range, and it rises not at all.

The regime is the expected one and is stated in §7.2 in advance: a membrane at tension `σ` resists a
point-like indentation linearly. **No closed form is fitted to the prefactor.** The measured
48.1 pN/µm is the same order as the tension scale `2πσ = 188.5 pN/µm` at `σ = 30` pN/µm, but the
logarithmic factor in the point-load result depends on a contact radius that here *is* the mesh
spacing, so quoting a recovered `σ` would be quoting the tessellation. It is not quoted.

### 8c. Where the encoding leaves its own regime, measured rather than argued

The two deepest points are the interesting failure. At `d = 0.60` the weld grows from 1 vertex to 3
and **2 of them go adhesive**; at `d = 0.70`, 6 of 7. The reported force *falls*, from 24.21 pN to
16.26 pN, because the newly welded vertices are being *pulled down* to keep them attached rather than
pushing back.

The cause is now clear: **the contact set is derived from the reference configuration, and by
`d = 0.6` µm the cell near the axis has deformed enough that the reference no longer says which
vertices are touching.** A vertex the undeformed geometry says is engulfed may in the deformed cell
not be in compressive contact at all. A proper unilateral contact would detect engagement from the
*current* configuration and release an adhesive site; this one cannot, by construction (§4.3).

`released_site_count` is what makes this visible, `in_regime` excludes those points from every fit,
and they are kept in the artefact rather than dropped. **The picture says it more plainly than the
number does**: in `depth_0.700um.png` the welded cap visibly bulges *out of* the indenter, poking
above the silhouette circle on both sides — a configuration a non-adhesive contact cannot produce.
At `d = 0.50` the same render shows a clean dimple wholly inside the circle.

So the sweep does **not** show a tension→Hertz crossover. It shows the tension-dominated branch, and
then the encoding leaving its validity before a Hertzian patch is ever resolved. Reaching that regime
needs a current-configuration contact search with release — a lane, not a parameter.

## 9. Deliberately failing negative control

Every refusal in §5 gets a control. In addition, the module ships **deliberate-break flags**, default
`False`, each driving one negative control rather than a hand-edited copy:

| Flag | The break | What must then misbehave |
|---|---|---|
| `omit_constraint` | `ACCUMULATE_FIELD` records the reaction but does not remove it | the welded vertices move off the surface and the reaction decays; I2 fails |
| `flip_reaction_sign` | the reported reaction is negated | the reported force disagrees in **sign** with the measured energy gradient; I4 fails |
| `include_pinned_rows_in_residual` | the convergence measure re-admits the welded rows | the residual jumps to the reaction, which cannot be driven to zero — the defect the design exists to avoid |

Each is asserted in **both** directions: the broken branch misbehaves, and the healthy branch is
exactly inert (`==`, `np.array_equal`, not "small").

### 9a. Three defects the implementation found, none of which raised

Recorded because all three produced plausible output rather than an error, which is the class of
failure this repository keeps paying for.

1. **Prescribing the membrane while leaving the cortex behind.** The cortex sits `initial_gap_um =
   0.010` µm inside the membrane; a 0.8 µm prescribed displacement drove **seven** membrane sites
   through it. `CortexMembraneContact` raised `ConnectorGeometryError` — correctly, and that refusal
   is not this lane's to weaken. The fix is to carry the **identical displacement vector** to the
   paired cortex vertices, which leaves every per-site separation bit-identical and makes the
   starting configuration admissible by construction. Asserted with `np.array_equal`.

2. **A frictionless contact with nothing to hold it.** §4.2a. The reaction decayed by a factor of
   2.4 while the tangential residual grew nine orders of magnitude. Nothing raised; the numbers were
   simply the wrong ones, and only the tangential-residual diagnostic showed why.

3. **A painter's-algorithm renderer that drew the contact patch and then painted over it.** Sorting
   triangles by centroid depth is wrong for a closed surface: near the silhouette the surface is
   nearly tangent to the view direction, so a whole cap collapses into a few pixel rows and centroid
   order stops predicting occlusion. Measured: **zero** contact-coloured pixels in a 400×320 render
   of a world whose contact patch was six triangles. The image looked entirely plausible — a cell,
   an indenter, a shaded mesh — and showed none of the thing it was made to show. Fixed with a
   per-pixel depth buffer and an elevated camera; the control counts contact-coloured pixels, which
   is why it was caught at all.

   The general lesson, since visualisation was requested precisely so the numbers would not be
   trusted alone: **a picture needs a control too.** An unchecked render is a second unverified
   claim, not a check on the first.

## 10. Numerical and precision envelope

float64 throughout, `NumpyBackend`. Welded-vertex immobility (I2) is asserted with `np.array_equal`,
not a tolerance: a constraint that *nearly* holds drifts over tens of thousands of steps.

### 10.1 Mesh resolution, measured

Free relaxation of the untouched vertical, this session, this tree:

| level | vertices | mean spacing | accepted | residual (pN) | `R_eq` (µm) | wall |
|---|---|---|---|---|---|---|
| 1 | 42 | 2.638 µm | 3,078 | 9.831e-04 | 4.823381 | 14.1 s |
| 2 | 162 | 1.381 µm | 3,012 | 9.998e-04 | 4.959428 | 25.1 s |
| 3 | 642 | 0.699 µm | 23,377 | 9.999e-04 | 4.995156 | 546.7 s |

**Level 1 cannot resolve a contact patch at any depth or radius.** The apex vertex's nearest
neighbours are 2.6 µm away, so with `R = 3 µm` the welded set stays at **one** vertex out to
`d = 1.8` µm — measured, not estimated. The suite therefore runs at level 2, which is the coarsest
mesh on which a patch exists at all.

Welded-set size against depth, `R = 3 µm`, measured on the relaxed references:

| depth (µm) | 0.05 | 0.1 | 0.15 | 0.2 | 0.3 | 0.5 | 0.7 | 0.8 | 1.0 | 1.2 | 1.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| level 2 | 1 | 1 | 1 | 1 | 1 | 1 | 7 | 7 | 7 | 7 | 7 |
| level 3 | 1 | 1 | 3 | 7 | 7 | 15 | — | 19 | 19 | 23 | 33 |

The quantisation is the icosphere's ring structure, and it is the reason a windowed fit is reported
rather than one global exponent: the patch does not grow smoothly, so `F(d)` has genuine steps in it
wherever a ring engages.

### 10.2 The dominant numerical limit is the descent, not the mesh

**The scaffolding relaxation cannot reach equilibrium at these depths in affordable time**, and this
is the single most important number in this entry. `aleph/vertical/relax.py`'s own docstring predicts
it — *"the contact is stiff and this scheme is explicit, so the step size is set by the stiffest
spring"* — and the derived contact stiffness here is `k = 2f/h` with `h = 0.020` µm, i.e. roughly
2,300 pN/µm against a cortex edge modulus of 50 pN/µm: a stiffness ratio near 50.

Measured, level 2, `R = 3 µm`, `d = 0.2` µm, welded: **126,000 steps and 980 s** to reach a free
residual of 3.0e-03 pN, still above the 1.0e-03 default tolerance. Extrapolating the observed 22×
wall scaling, level 3 is out of reach for a multi-depth sweep on CPU, and no GPU authorization is
held.

The consequence is carried in the artefact rather than hidden: every point records a **trace** of
`(steps, residual, force, energy)`, and `IndentationCurve.fit_at_residual` /
`exponent_against_stopping_residual` refit the exponent at a **matched** stopping residual across
depths. Measured on an early two-depth probe, the exponent moved **0.485 → 0.700 → 1.55** as the
matched residual tightened from 0.1 to 0.02 pN — so a single-stopping-point exponent from this
solver is not a measurement of the physics. This is `ALEPH-DQ-104` appearing as a number rather than
as a design question.

**With a large enough budget it does converge, and the drift on the way is the headline number.**
The production sweep's first rung, level 2, `R = 3 µm`, `d = 0.05 µm`, welded set of one:

| accepted steps | free residual (pN) | reaction (pN) | energy (pN·µm) |
|---|---|---|---|
| 220 | 1.9638e-01 | 4.505800 | 9302.093173 |
| 880 | 5.4404e-02 | 3.780209 | 9302.076839 |
| 3,518 | 4.6170e-02 | 3.182110 | 9302.060362 |
| 14,068 | 2.4677e-02 | 2.270586 | 9302.026977 |
| 56,266 | 6.1568e-03 | 0.643594 | 9301.992798 |
| **98,319** | **9.9921e-04** ✔ | **0.186449** | 9301.989986 |

**The reaction falls by a factor of 24.2 between the first sample and convergence**, monotonically,
with no plateau until the very end. Any budget-limited reading of this curve overestimates the force
by roughly the ratio of its residual to 1e-03. 863 s for one rung, and there are fourteen.

A second warning the same table carries: the converged indented energy, 9301.989986 pN·µm, is only
2.6e-04 above the free reference's 9301.989723 — but the reference was itself stopped at a residual
of 9.998e-04 after only 3,012 steps, whereas this rung ran 98,319. **The two energies are therefore
not comparable to that precision**, and any `ΔE` read across them is dominated by how much better
relaxed the indented state is. It is a third independent reason the *relaxed* virtual-work identity
(§7.1) is `UNVERIFIED`, and it would be invisible without the trace.

## 11. Production-backend residency and transfer

Host, numpy, float64, `NumpyBackend`. **No GPU work was run and no authorization was sought or
held.** `~/.aleph_data/gpu_authorization.json` is expired and names a different host.

## 12. Comments and docstrings to discard

Nothing to discard; no reference prose entered, because no reference file was read.

## 13. Visualisation, and an undeclared-dependency finding

The PI's instruction is explicit that the numbers must not be trusted alone —
*"단순히 데이터들만을 믿지 말고 시각화도 같이 시켜서 시각화로도 판단해야돼"*. So a rendered artefact
is produced **per indentation depth**.

**`aleph/viz/` cannot do it, and that is a property of `aleph/viz/` rather than a complaint.** Its own
`__init__` docstring states there is *"no viewer process: this is the packing and the contract"*. It
supplies `SceneTree`, `AcceptedFrame`, `EvidenceCapture` and a transport contract — provenance and
isolation machinery — and no rasteriser. `aleph/viz/` may not be modified by this lane in any case.

**`matplotlib` is not installed in the `aleph` environment and is not a declared dependency in
`pyproject.toml`** (`dependencies = ["numpy>=1.26", "scipy>=1.11"]`). Verified by import: it raises
`ModuleNotFoundError`. Reported here as a finding rather than worked around by installing it.

Therefore `indent.py` carries a **self-contained renderer using numpy and the standard library
only**: a triangle rasteriser with a per-pixel **depth buffer**, writing PNG through `zlib` +
`struct`; and an SVG writer for the force–indentation curve. No third-party dependency is added and
nothing outside the allocation is touched.

Three decisions in it were forced by what the pictures showed:

1. **A depth buffer, not a painter's algorithm** — §9a.3. The painter's version rendered zero
   contact-coloured pixels.
2. **An elevated three-quarter camera, not a side view.** The contact patch is a cap on top of the
   cell; side-on it foreshortens into a band a few pixels tall.
3. **The contact patch is drawn in a contrasting hue**, and a control counts those pixels. The
   picture has to be able to fail.

**What the first frame showed that the numbers did not.** At `d = 0.05` µm the deformation is a
**flat hexagonal facet** — the one welded vertex and its six neighbours — surrounded by a visibly
undisturbed cell. The table says `weld=1, resolved=False`; the image says *why* that label matters,
because a contact patch that is literally one facet of the tessellation is not a contact patch. The
picture also confirms `intruders=0` by eye: no part of the cell crosses the indenter silhouette.

The curve is emitted as **SVG** rather than a raster because a curve is vector data, and writing it
as text keeps every plotted coordinate checkable against the JSON without decoding an image. Both
reference slopes — 1 (tension-dominated) and 3/2 (Hertz) — are drawn against the data, so the eye
compares the measurement against *both* candidate laws rather than against a single fitted line that
will look convincing through any three points.

## 14. Honest limits — what this entry does NOT establish

The single most important one first.

- **`UNVERIFIED`: no Hertzian contact regime was reached, so the 3/2 law was never tested.** §8c. The
  measured exponent of 0.988 is the *tension-dominated* branch under an effectively point-like
  indenter, and the encoding leaves its regime (adhesive welded sites) as soon as the patch grows to
  3 vertices. **The curve is evidence about the linear branch and about nothing else.** Reaching a
  Hertzian patch needs a current-configuration contact search with release, and a finer mesh.
- **`UNVERIFIED`: mesh convergence.** The sweep is level 2 only, and the point-load branch is
  measured with a **one-vertex** contact, so its prefactor is a property of the tessellation as much
  as of the cell. The exponent is far more robust than the prefactor — that is the whole reason §7.2
  makes the exponent the oracle — but neither has a discretisation error bar. Level 3's free
  reference was computed (555 s, `R_eq = 4.995156` µm, cached at
  `runs/indentation/reference_level3.npz`); a level-3 sweep is out of reach on CPU at the observed
  22× wall scaling, and no GPU authorization is held.
- **The descent is the cost bottleneck and it dominated this lane.** §10.2. Every point above needed
  30,000–56,000 accepted steps and 280–490 s. That is affordable at level 2 and not at level 3. An
  implicit or preconditioned step for a stiff contact is the lane that unblocks everything
  downstream here.
- **`UNVERIFIED`: the relaxed virtual-work identity** `dU_equilibrium/dd = F`. §7.1. The gradient
  form is verified to second order; the relaxed form was measured at −16.4 pN against +22.1 pN and
  the disagreement is entirely the non-convergence above.
- **`UNVERIFIED`: ladder independence.** The weld position of a vertex is fixed at the depth it
  enters, so the curve depends on the depth ladder. Whether the chosen ladder is fine enough to
  approximate the continuous limit was **not** measured — it needs two ladders, and one sweep is
  already at the edge of the compute budget.
- **This is not an AFM measurement of a cell.** A Helfrich membrane, an elastic shell and a
  volume-penalty compartment, all parameters `UNSOURCED`, pressed by a rigid sphere. No modulus
  recovered from it may be compared to a literature value.
- **No modulus is reported at all.** Only the exponent and the regime, because `E*` is the one thing
  the 3/2 law cannot supply without assuming the half-space geometry that does not hold here.
- **The contact is bonded, Hertz is frictionless.** Even where the geometry were right, the boundary
  condition differs, which inflates the stiffness by an O(1) factor that is not computed here.
  `friction_fraction` measures how much this matters at each depth.
- **The acceptance predicate is still open** (`ALEPH-DQ-104`). `HANDOFF.md` §5.2 records twenty
  stopping points on one mesh giving signed errors spanning 603× with five sign changes. Every point
  on this curve inherits that, and each records the residual it stopped at.
- **No poroelasticity, no cytosol** — §4.4.
- **The curve is quasi-static and single-branch**: loading only, no retraction, no hysteresis, no
  adhesion. Retraction is refused rather than approximated.
- **The indenter is large relative to the cell.** `R = 3 µm` against `R_eq ≈ 4.96 µm` is a ratio of
  0.6, forced by the mesh resolution rather than chosen. A small AFM bead cannot be resolved on any
  mesh this solver can afford.

## 15. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Encoding accepted; measurement provisional.** 77/77 controls pass, 4/4 mutants killed, the gradient oracle converges at order 2.0010, and the free-body defect of §4.1a is fixed and controlled. The force–indentation curve is `UNVERIFIED` pending the sweep in §8b and the limits in §14. |
| Reviewer | Agent-proposed. Unratified. No PI review. |
| Rollback | Delete `aleph/scenarios/indent.py`, `tests/scenarios/test_indent.py`, `runs/indentation/` and this entry. Nothing else imports any of them; `aleph/scenarios/world.py`, `aleph/vertical/**`, `aleph/state/**` and `tests/firewall/**` were not touched. |

### 15a. Deviations from this entry as written, and why each moved

Per `PLAN.md` §0.2.5 the entry constrains the design; where writing the code made a section wrong,
the section is corrected here rather than left standing.

| § | As written | As built | Why |
|---|---|---|---|
| §4.2 | frictionless contact, sliding, reprojected each step | **welded** contact, rigid translation | §4.2a — the frictionless equilibrium does not exist in this cell model, measured |
| §5 | resolution refused at world construction | refused at **fit admission**, recorded at construction | §5a — a welded ladder needs its shallow rungs |
| §7.1 | `dU/dd` across three relaxed depths | gradient under rigid translation at fixed free DOFs | §7.1 — the relaxed form measures relaxation drift at achievable convergence |
| §10 | sweep at level 3, repeated at level 2, difference reported as the error bar | level 2 only | §10.2 — level 3 is out of reach on CPU and no GPU authorization is held |
| §13 | one rendered artefact per depth | unchanged, but the renderer needed a depth buffer | §9a.3 — the painter's version drew zero contact-coloured pixels |
| §4 | indenter only | indenter **and** a bonded substrate footprint | §4.1a — an unsupported cell escapes by rigid translation and the reaction decays to zero |

The substrate is the one that matters. The entry reasoned carefully about the indenter, the contact
law, the resolution and the oracle, and never asked *what holds the cell down* — which is the first
question an experimentalist would ask, and the only one whose omission silently produced a complete,
plausible, converged, wrong curve.

---

## 16. The contact law, replaced: a current-configuration unilateral contact with release

**Written before the code that implements it**, per `PLAN.md` §0.2.5, and amended after with what
the code measured. §16.7 lists every place where writing it made this section wrong.

### 16.1 The defect this replaces, in one line

§8c: **the contact set was derived from the reference configuration and a welded vertex could never
let go.** By `d = 0.6` µm the cell near the axis has deformed enough that the undeformed geometry no
longer says which vertices are touching, so vertices were held that a real contact would have
released — 2 of 3 at `d = 0.60`, 6 of 7 at `d = 0.70` — and the reported force *fell* with depth,
24.21 → 16.26 → 13.85 pN, because those vertices were being pulled down to keep them attached. In
`depth_0.700um.png` the welded cap visibly bulges *out through* the indenter silhouette.

The consequence is not a bad number: it is a **missing measurement**. The encoding left its regime
before a Hertzian contact patch was ever resolved, so the 3/2 law was never tested at all.

### 16.2 The condition, which is the whole content

Signorini, per vertex, in the **current** configuration. With the signed gap to the indenter surface

```
g_i(x) = |x_i − c(d)| − R          (positive outside the sphere, negative penetrating)
```

and the compressive reaction magnitude along the outward surface normal `n_i = (x_i − c)/|x_i − c|`

```
λ_i = −(F_i^cell · n_i)
```

the unilateral contact conditions are

```
g_i ≥ 0 ,        λ_i ≥ 0 ,        g_i λ_i = 0 .
```

*Non-penetration, compression only, and no force at a distance.* A contact that can only stick is
not a contact: the third condition is what the previous encoding violated, and the second is what it
had no way to check except after the fact (`released_site_count`).

### 16.3 How it is imposed — an active set, not a penalty

The active set `A` is a **result**, not an input. It is solved for by the standard primal-dual
active-set iteration, interleaved with the descent:

1. **Relax** for a chunk of steps with `A` held fixed. On `A` the constraint is the same exact
   Dirichlet row-zeroing §4.2 already uses — `backend.scatter_add(forces, A, −forces[A])` — so
   `x + (−x)` is exactly `0.0`, the active vertices are bit-identical across the chunk, and the
   virtual-work oracle §7.1 is unchanged in form and still exact.
2. **Release.** Evaluate `λ_i` for every `i ∈ A` from `evaluate_forces` (the *unconstrained* cell
   force, which is what the reaction is computed from). Drop every `i` with `λ_i < 0.0`. The
   comparison is against **exact zero**, not a tolerance, for the same reason
   `UnilateralSpring.energy_pn_um` returns exactly `0.0` on its inactive branch: a unilateral
   element's inactive branch is exactly inert or it is not unilateral.
3. **Engage.** Evaluate `g_j` for every `j ∉ A`. Add every `j` with `g_j < 0`, placing it on the
   surface by closest-point projection — the smallest displacement that satisfies non-penetration —
   and **carrying the paired cortex vertex by the identical vector** (§9a.1; without it
   `CortexMembraneContact` refuses the configuration, correctly).
4. **Repeat** until the descent has converged *and* the active set did not change. **Both**, and
   that is not a detail: a residual-converged state with the wrong active set is a converged
   solution of the wrong problem, which is exactly what §8b reported at `d ≥ 0.6`.

A vertex released at step 2 is excluded from step 3 **within the same iteration only**. It sits at
`g = 0` with a force pointing out of the indenter, so it separates on the next chunk; re-testing it
in the same iteration would re-engage it at the position it was just released from and the set would
chatter without moving. It is reconsidered from scratch at the next iteration.

### 16.4 Why an active set and not the penalty form of `UnilateralSpring`

`aleph/vertical/connectors.py` already carries a compressive-only unilateral element with an exact
zero on the inactive branch (`ContactLaw.SOFT_CORE`, `NonAdhesiveContact`,
`CortexMembraneContact`), and the discipline is reused rather than re-invented: **exact-zero
inactive branch, `!= 0.0` engagement test, a `load_carried_in_forbidden_sense_pn` audit that must
read exactly `0.0`.** What is *not* reused is the penalty *mechanism*, for two measured reasons:

1. **A penalty contact against a rigid sphere is one more stiff spring in an explicit descent.**
   §10.2 already measures this scheme's cost as set by the stiffest spring — the derived contact
   stiffness is ~2,300 pN/µm against a cortex edge modulus of 50 pN/µm — and a penalty stiff enough
   to keep `g ≥ 0` to the precision this measurement needs would dominate the step size again.
2. **A penalty contact only approximates `g ≥ 0`.** It admits penetration proportional to the load,
   which is a depth error, and a log–log fit converts a systematic error in `d` straight into a
   wrong exponent — the same argument §4.2 makes for computing the first-touch height rather than
   assuming it.

The active-set form gives `g_i = 0` on `A` **exactly** and `λ_i = 0` off `A` **exactly**, at every
configuration, and it costs no stiffness at all.

### 16.5 What is still not Hertzian about it, stated in advance

- **Bonded while active.** A vertex in `A` is held against the *whole* force vector, tangential
  included, exactly as in §4.2a — the frictionless variant has no well-posed equilibrium in this
  cell model, measured. Release makes the contact *unilateral in the normal sense*; it does not make
  it frictionless. `friction_fraction` still reports how far the answer is from the frictionless
  closed form, and it is now reported over a range where the patch is more than one vertex.
- **Still quantised by the tessellation.** Release cannot manufacture a contact radius the mesh
  cannot resolve. If the Hertzian branch needs `a_c ≫ h`, and `a_c = sqrt(2Rd)` against `h = 1.381`
  µm at level 2 reaches only ~2.7 µm at the deepest affordable depth, then **the right answer is
  that the branch is out of reach, reported with the numbers that show it** — not a fitted 1.5.

### 16.6 The invariants this section adds, each with its control

- **I8 (complementarity).** `adhesive_load_pn()` — the total `max(−λ_i, 0)` over the active set — is
  **exactly `0.0`** after an active-set update. This is the direct analogue of
  `_DistributedUnilateralConnector.load_carried_in_forbidden_sense_pn`, and it is the invariant the
  first encoding could not state, only count violations of.
- **I9 (non-penetration).** `penetration_um()` — the largest `max(−g_j, 0)` over vertices *not* in
  the active set — is `0.0` after an update. Together I8 and I9 are the KKT residual of §16.2.
- **I10 (the tension branch survives).** The shallow-depth exponent measured with release must agree
  with the 0.988 of §8b to within the fit's own quality. **If it moves, the contact change altered
  the physics rather than the contact set, and that is a finding to report and not to tune away.**
- **I11 (the release is real).** `allow_adhesive_contact=True` — a shipped deliberate-break flag —
  reinstates the §8c defect exactly: adhesive sites reappear and the force falls with depth. The
  healthy branch must be exactly inert.
- **I12 (the search is current, not reference).** `search_in_reference_configuration=True` — the
  second shipped break — restores the original contact set. It must select a *different* set from
  the current-configuration search at a deformed configuration, and the identical one at the
  reference configuration.

### 16.7 Deviations from §16 as written

Per `PLAN.md` §0.2.5 the entry constrains the design; where writing the code made a paragraph above
wrong, it is corrected here rather than left standing.

| § | As written | As built | Why |
|---|---|---|---|
| 16.3 | a vertex is engaged when `g < 0`, exactly, and a just-released vertex is kept out **by name** for one iteration | engaged when `g < -8·eps·R`, **and** kept out by name for one iteration | Name-exclusion covers the same iteration and not the next one. Measured: a released vertex reads `g = -4.441e-16` µm at `R = 3` µm — round-off in `|x−c| − R`, which differences two numbers of size `R` — and with a strict `g < 0` test it was re-engaged on the next iteration, released again, and the active set **chattered forever without the configuration moving**. The convergence criterion "the set stopped changing" would then never have been reachable. See §16.7a: this is a statement about float64, not a contact tolerance, and the distinction is load-bearing. |
| 16.6 I9 | `penetration_um()` is exactly `0.0` after an update | `penetration_um() <= gap_floor_um()`, and the value is reported **raw** | Same cause. The artefact carries the measured `4.441e-16` µm rather than a zero manufactured by subtracting the floor, so a reader checks the number against the floor instead of being handed the answer. |
| 16.4 | *(implied)* the reaction may be read from the indenter's recorded buffer | `lambda` is always computed from `evaluate_forces` | The post-constraint accumulator is exactly zero at every held row by construction, so every `lambda` read from it would be exactly zero and the release would never fire. The buffer is additionally **refreshed** whenever the active set changes; a buffer left over from the previous set is not a crash, it is a plausible `friction_fraction` about the wrong vertices. |

#### 16.7a The round-off floor is not a contact tolerance, and a surviving mutant is why that is said twice

`GAP_ROUNDOFF_FACTOR` deserves its own paragraph because the mutation exercise showed the suite
could not tell the two apart. **Mutant M6 replaced the floor `8·eps·R ≈ 5.3e-15` µm at the
engagement test with a hard `1e-3` µm — a physical-looking contact tolerance — and passed all 100
controls.** Nothing in the suite exercised a vertex sitting *between* the two scales, so the entire
difference between "a computed zero" and "a tuned threshold" was invisible.

It is not a cosmetic difference. Any threshold above round-off licenses interpenetration by exactly
its own size; interpenetration is a systematic error in `d`; and a log–log fit converts a systematic
error in `d` straight into a wrong exponent. That is the same argument §4.2 already makes for
computing the first-touch height instead of assuming it, arriving from the other direction.

`test_a_penetration_far_below_any_plausible_tolerance_is_still_engaged` was written for it,
parametrised over intrusions of `1e-4`, `1e-8` and `1e-11` µm — seven decades, all far above the
round-off floor and all below any tolerance a reasonable person would type. M6 now dies three times.

### 16.8 Measured

Every number below is from `runs/indentation/indentation_release_level2_R3.json` and its log. The
run used a supported level-2 cell (162 membrane vertices, mean spacing 1.3812715 µm), a 3.0 µm
indenter, and the default `1e-03` pN free-residual threshold. Wall time was 21,284 s on CPU.

| d (µm) | F (pN) | active | released | engaged | adhesive load (pN) | penetration (µm) | residual (pN) | converged |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| 0.05 | 2.492640 | 1 | 0 | 0 | 0.0 | 0.0 | 9.980e-4 | yes |
| 0.10 | 4.935279 | 1 | 0 | 0 | 0.0 | 0.0 | 9.990e-4 | yes |
| 0.15 | 7.380583 | 1 | 0 | 0 | 0.0 | 0.0 | 9.994e-4 | yes |
| 0.20 | 9.816335 | 1 | 0 | 0 | 0.0 | 0.0 | 9.998e-4 | yes |
| 0.30 | 14.650076 | 1 | 0 | 0 | 0.0 | 0.0 | 9.993e-4 | yes |
| 0.40 | 19.450679 | 1 | 0 | 0 | 0.0 | 0.0 | 9.973e-4 | yes |
| 0.50 | 24.210404 | 1 | 0 | 0 | 0.0 | 0.0 | 9.983e-4 | yes |
| 0.60 | 28.930981 | 1 | 0 | 0 | 0.0 | 0.0 | 9.979e-4 | yes |
| 0.70 | 33.576565 | 1 | 0 | 0 | 0.0 | 0.0 | 9.982e-4 | yes |
| 0.80 | 38.185761 | 1 | 0 | 0 | 0.0 | 0.0 | 9.993e-4 | yes |
| 1.00 | 47.296642 | 1 | 2 | 0 | 0.0 | 0.0 | 9.975e-4 | yes |
| 1.20 | 58.009399 | 3 | 8 | 4 | 0.0 | 0.0 | 9.970e-4 | yes |
| 1.50 | 82.399609 | 7 | 4 | 4 | 0.0 | 0.0 | 1.976e-3 | **no** |

The active set is stable at every saved rung. The 1.50 µm point is kept in the artefact and excluded
from the global fit because its residual did not converge. The twelve admitted points give

```text
F = 47.861684 d^0.985883,  r² = 0.99996373,  d = 0.05 … 1.20 µm.
```

Four-point windowed exponents run from 0.9621 to 1.0082. Release therefore repairs the falling-force
defect without manufacturing a Hertzian branch: the result remains the approximately linear,
tension-dominated regime. Only two of thirteen saved rungs have at least three active vertices; the
largest patch has seven. The mesh still does not resolve a Hertzian contact patch, so `3/2` remains
**untested**, not rejected.

Release is load-bearing. It first fires at 1.00 µm; at 1.20 µm the active-set iteration records eight
releases and four engagements and finishes with three active sites, zero adhesive load and zero
penetration. The old `allow_adhesive_contact=True` / reference-search run is retained under
`docs/results/2026-07-31-indentation/before_reference_search_depth_*.png` as the negative control.

### 16.9 Acceptance after the amendment

`tests/scenarios/test_indent.py` now collects **103 controls**, and the full 5,032-test repository
run on 2026-08-03 passed every scenario control. The original four mutation checks remain recorded
in §8a; §16.7a's hard-tolerance mutant is now killed at `1e-4`, `1e-8` and `1e-11` µm penetration.
The encoding is accepted at `ANALYTIC_ORACLE`; the quantitative curve remains provisional because
the Hertzian patch is unresolved, level 3 is unmeasured, and the deepest rung did not converge.
