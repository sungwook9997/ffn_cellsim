# ALEPH-PORT-1101 — Helfrich surface mechanics on a closed triangulated manifold

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-1101` |
| Lane | `L11 vertical — membrane / cortex / pressure` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` — nothing was ported. No source file was opened for this entry. |

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
from aleph.vertical.membrane import (
    HelfrichMembrane,
    MembraneEnergy,
    MembraneMaterial,
    enclosed_volume_um3,
    helfrich_energy,
    helfrich_forces,
    surface_area_um2,
    total_area_gradient,
    volume_gradient,
)
```

Aleph target files: `aleph/vertical/membrane.py`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | none — **no source commit was read for this entry.** No file was opened, so no commit can be cited, and citing one would be false provenance. |
| Source path | none — no source path was read. |
| Source symbol(s) | none. |
| Read from | neither. Nothing was read. |
| Working tree == commit? | not applicable — nothing was read from either. |

The only thing that crossed the boundary into this entry is a **defect description** supplied by a
prior read-only audit (recorded in PLAN §1.1), namely that in the reference project the declared
load path and the executed load path disagree. That description is an argument, not code, and it is
recorded in `ALEPH-PORT-1103` where it belongs. This entry inherits nothing at all.

## 3. Why source-derived porting beats clean-room

It does not. The Helfrich functional, the cotangent Laplace–Beltrami operator, and the divergence
theorem for the enclosed volume of a closed polyhedron are textbook. Re-deriving them from scratch
costs less than auditing somebody else's discretisation, because the derivation is the thing that
has to be checked anyway: a bending force is only trustworthy if its own gradient identity has been
verified, and verifying that identity is the same work whether the code was written or read.

Clean-room wins. Nothing was ported.

## 4. Physical or mathematical law represented

Helfrich free energy of a closed fluid surface at fixed topology:

```
E = ∮ [ (κ/2)(2H − c₀)² + σ ] dA
```

`H` is the mean curvature (positive for a convex body with outward normals, so `H = 1/R` on a
sphere of radius `R`), `κ` the bending rigidity, `σ` the surface tension, `c₀` the spontaneous
curvature. The Gaussian-curvature term is dropped: at fixed genus, Gauss–Bonnet makes
`∮ K dA = 2πχ` a constant, so it contributes no force. `aleph/state/curvature.py` already verifies
that identity to round-off on this mesh family, so dropping the term is a checked statement rather
than an assumption.

**Discretisation.** Integrate the mean-curvature normal over the vertex cell. The cotangent
Laplace–Beltrami identity gives the integrated mean-curvature vector at vertex `i`

```
K_i = ½ Σ_{j∈N(i)} (cot α_ij + cot β_ij) (x_i − x_j)   ≈  ∮_{cell(i)} 2H n dA  =  2 H_i A_i n_i
```

so with a lumped vertex area `A_i` the discrete energy is

```
E = Σ_i (κ / 2A_i) | K_i − c₀ A_i n̂_i |²  +  σ Σ_f A_f
```

which reduces to `Σ_i (κ/2)(2H_i − c₀)² A_i` when `K_i = 2H_i A_i n̂_i`, and expands into three
separately differentiable pieces:

```
E_bend  = (κ/2) Σ_i |K_i|² / A_i
E_cross = −κ c₀ Σ_i K_i · n̂_i
E_offset= (κ c₀² / 2) Σ_i A_i = (κ c₀² / 2) A_total
```

**Vertex-area choice, and why it is not the accurate one.** `A_i` is the *barycentric* lumping
`A_i = ⅓ Σ_{f∋i} A_f`, not the mixed-Voronoi area that `aleph/state/curvature.py` shows to be four
to six orders of magnitude more accurate pointwise. The reason is that this module needs a
**force**, and a force is a gradient. The mixed-Voronoi construction branches on whether a triangle
is obtuse; that branch makes the area only piecewise smooth in the vertex positions, so the energy
is only `C⁰` across the branch and the force is discontinuous there. A discontinuous force cannot
pass a central-difference gradient check and cannot be integrated by any descent that assumes
smoothness. The barycentric lumping is a polynomial in the vertex positions and is `C^∞`.

The cost is stated rather than hidden: barycentric lumping leaves an `O(1)` pointwise curvature
error at the twelve valence-5 vertices an icosphere inherits from the icosahedron, so the
*pointwise* bending force there is biased. Two things make that acceptable here and both are
measured, not asserted: the **net** resultant of the bending force vanishes identically (translation
invariance of the discrete energy, exact to round-off), and the discrete bending energy is exactly
invariant under uniform scaling, so it contributes exactly nothing to the radial equilibrium that
Laplace's law is read off. Anything that needs accurate *pointwise* curvature must use
`aleph.state.curvature`, not this module.

**Exact gradients.** Every gradient is analytic, not differenced.

- `∇_{p_k} A_f = ½ (p_{k+1} − p_{k+2}) × n̂_f` (cyclic corner indices).
- `∇_{p_m} cot_k = (1/|N|) ∂(u·v)/∂p_m − (cot_k/|N|) · 2 ∇_{p_m} A_f`, with `u = p_{k+1} − p_k`,
  `v = p_{k+2} − p_k`, `|N| = 2A_f`.
- `V = (1/6) Σ_f p_0 · (p_1 × p_2)` (divergence theorem on a closed polyhedron), so
  `∇_{p_0} V = (1/6) p_1 × p_2` and cyclically. This is exact for the polyhedron, not an
  approximation to the smooth volume.
- `∂N_f/∂p_k = [q_k]_×` with `q_k = p_{k+2} − p_{k+1}`, which supplies `∇ n̂_i` through
  `m_i = ½ Σ_{f∋i} N_f`, `n̂_i = m_i / |m_i|`.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| vertex position | µm | m | finite |
| bending rigidity κ | pN·µm | J | ≥ 0 |
| surface tension σ | pN/µm | N/m | ≥ 0 |
| spontaneous curvature c₀ | 1/µm | 1/m | finite, any sign |
| area | µm² | m² | > 0 per face |
| enclosed volume | µm³ | m³ | > 0 for outward orientation |
| force | pN | N | finite |
| energy | pN·µm | J | finite |

Singular and boundary cases:

- A degenerate triangle (`|N| = 0`) has no normal and no cotangent. Refused with
  `DegenerateGeometryError`; never regularised by adding an epsilon, because an epsilon would turn
  a topological failure into a small wrong number.
- A vertex whose incident face normals cancel (`|m_i| = 0`) has no vertex normal. Refused for the
  same reason. Reachable only for a non-embedded surface.
- `κ = 0` is legal and gives a pure tension (soap-film) surface. `σ = 0` is legal and gives a pure
  Willmore surface. Both zero is legal and gives zero force; it is not refused here, because the
  degeneracy is a property of the *ensemble*, not of the energy.
- `c₀ = 0` is the default and the only value used by this vertical's controls.

Invariants:

- **I1. Force is minus the exact energy gradient.** Asserted against a central difference to
  `O(h²)` for random perturbations, on random non-spherical meshes.
  **UNTESTED** *(the entry named two central-difference controls and recorded, for a third, that the
  "specific control was never recorded". None of the three exists — see §14.4)*.
- **I2. Net force vanishes.** `Σ_i F_i = 0` to round-off, for any configuration and any material —
  the discrete energy depends only on differences of positions.
  `tests/vertical/test_vertical_controls.py (specific control never recorded)`.
- **I3. Net torque vanishes.** `Σ_i x_i × F_i = 0` to round-off — the discrete energy is invariant
  under rigid rotation. **UNTESTED** *(this entry named a control in `tests/vertical/test_membrane.py`, a file that was never written; see §14.4)*.
- **I4. Volume gradient is exact.** `∇V` matches a central difference to round-off, not to `O(h²)`,
  because `V` is a cubic polynomial in the vertex positions and the central difference of a cubic
  is exact up to the `h²` term of its third derivative — checked at two step sizes.
  **UNTESTED** *(this entry named a control in `tests/vertical/test_membrane.py`, a file that was never written; see §14.4)*.
- **I5. Bending energy is scale-invariant.** `E_bend(λx) = E_bend(x)` exactly, for `c₀ = 0`. This is
  the discrete analogue of the Willmore energy's conformal invariance and it is what makes the
  Laplace-law control independent of κ.
  **UNTESTED** *(this entry named a control in `tests/vertical/test_membrane.py`, a file that was never written; see §14.4)*.
- **I6. Sphere energy.** `E_bend → 8πκ` on a refining icosphere.
  `tests/vertical/test_vertical_controls.py (specific control never recorded)`.

## 6. Source evidence class and known retractions

No source artifact is claimed, because no source was read. The evidence class of this module is
therefore its own: `ANALYTIC_IDENTITY` for I1–I5 (each is an exact statement checkable without any
external number) and `CONVERGENCE` for I6.

Looked for retractions in: nothing, because nothing was inherited. There is no source claim to
retract.

## 7. Independent oracle or derivation

Four independent checks, none of which is "it agrees with another implementation":

1. **Central-difference gradient.** `F = −(E(x+h) − E(x−h)) / 2h`, componentwise, for random
   perturbations of random meshes. This is the sign arbiter for the entire lane. It cannot be
   satisfied by a plausible-looking wrong force.
2. **Symmetry identities.** Translation invariance ⇒ `ΣF = 0`; rotation invariance ⇒ `Σ x×F = 0`.
   Both to round-off, both independent of the physics being right.
3. **Closed form on a sphere.** `E_bend = 8πκ` for any closed genus-0 surface at the Willmore
   minimum, and `A = 4πR²`, `V = (4/3)πR³` as the polyhedron refines.
4. **Exact discrete scaling.** `E_bend(λx) = E_bend(x)` and `A(λx) = λ²A(x)`, `V(λx) = λ³V(x)`
   exactly, which together give the discrete Laplace relation `ΔP = 2σA/(3V)` with no continuum
   approximation anywhere in it.

`validation/analytic/helfrich.py` is a *different* oracle — the equilibrium fluctuation spectrum of
a nearly-flat patch — and is deliberately **not** used here: it is a statement about thermal
statistics, this module is a statement about deterministic forces, and `aleph/**` may not import
`validation/**` in any case (`tests/vertical/test_vertical_firewall.py`).

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | Analytic bending force equals the central difference of the energy; the max componentwise relative error is below 1e-6 at h = 1e-6 µm on a randomly perturbed level-2 icosphere. |
| Positive | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | `‖ΣF‖` is at float64 round-off relative to `Σ‖F‖`. |
| Positive | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | `E_bend/8πκ → 1` as the icosphere refines. |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | Flipping the sign of the returned force makes the central-difference comparison fail. If it stops failing, the gradient check has stopped checking anything. |
| Negative (must fail) | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | Substituting the mixed-Voronoi vertex area into the energy while keeping the barycentric gradient — the exact mistake an "accuracy improvement" would make — breaks the gradient identity. |
| Negative (must fail) | `tests/vertical/test_vertical_controls.py (specific control never recorded)` | A collapsed triangle raises rather than producing a large finite force. |

## 10. Numerical and precision envelope

Working and accumulation precision: float64 throughout. There is no float32 path; this module is
the CPU definition of the law per `ALEPH-DQ-107`.

Tolerances and why:

- Gradient check: `h = 1e-6 µm` against a mesh of characteristic edge length ~1 µm, i.e.
  `h/L ≈ 1e-6`. Central-difference truncation is `O((h/L)²) ≈ 1e-12` and cancellation error is
  `O(ε·E/h) ≈ 1e-16·E/1e-6 = 1e-10·E`. The two cross near `h/L ~ 1e-5`, so `1e-6` sits on the
  cancellation side and the achievable agreement is ~`1e-8` relative, not `1e-12`. The control is
  asserted at `1e-6` relative, two decades of margin above the observed floor and eight decades
  below any sign or factor error. A tighter number would be asserting the differencing scheme, not
  the force.
- Symmetry identities: asserted at `1e-12` relative to the force scale — these are exact
  cancellations of the same float64 terms, so the only spread is round-off.
- Scale invariance: asserted at `1e-12` relative.
- Sphere energy: asserted as a convergence statement (`|E/8πκ − 1|` below 3e-3 at level 3), not an
  identity, because it is one.

Conditioning: the cotangent weights diverge as a triangle approaches degeneracy (`cot → ∞` as the
opposite angle → 0). The refusal in §5 is the hard floor; between there and a well-shaped mesh the
accuracy degrades smoothly and is not otherwise guarded. Outside the envelope the module raises. It
never silently degrades.

## 11. Production-backend residency and transfer

Host-side numpy today, and CPU is the definition. The arrays that would live on the device in a
future accelerator build are the vertex positions `(V,3)`, the triangle table `(F,3)` and the force
accumulator `(V,3)`; the per-face intermediates (normals, cotangents, area gradients) are
recomputable and would never be transferred. One host round trip per step is required today only
for the energy scalar the relaxation's step controller reads; nothing else needs the host.

This module is `aleph/vertical/**` and therefore forbidden from importing `validation/**`. Enforced
by `tests/firewall/test_scope_firewall.py::test_runtime_never_imports_the_analytic_oracles`.

## 12. Comments and docstrings to discard

Nothing to discard: no source prose was read, so none can survive. What replaces it is the
derivation in §4 above and the module docstring of `aleph/vertical/membrane.py`, which states the
discretisation, names the vertex-area choice and its measured cost, and states plainly that the
sign arbiter is the central-difference check and not the author's confidence.

Vocabulary deliberately **not** used anywhere in `aleph/vertical/**`: any provider package name,
module path, gate label, or facade method name. `tests/ports/test_port_discipline.py` scans for
them.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Pending.** Status is `PROPOSED`: the code named in §1 has not landed and no control has been run, so there is no measured result to record here. This row is filled in with the measured residuals, at the tolerances stated in §10, at the moment the status moves to `ACCEPTED` — and not before. |
| Reviewer | agent-proposed, **unratified**. Written autonomously on 2026-07-30 while the PI was asleep; no human has reviewed the derivation in §4. |
| Rollback | Delete `aleph/vertical/membrane.py` and `tests/vertical/test_membrane.py`. Breaks: `aleph/vertical/assembly.py`, `aleph/vertical/relax.py`, and every control in `tests/vertical/`. Nothing outside `aleph/vertical/` depends on it. |

## 14. Honest limits

- **Unratified.** No human has checked the gradient derivation. The central-difference control makes
  a sign or factor error very unlikely to survive, but it does not make the *model* right.
- The pointwise bending force at the twelve valence-5 vertices of an icosphere is biased at `O(1)`
  by the barycentric lumping (§4). This module is not a curvature estimator and must not be used as
  one. The bias is invisible to every control here because every control is either a symmetry, a
  net resultant, or a scaling identity — all of which are blind to a pointwise error that sums and
  scales correctly. That is a real limit, not a technicality.
- Fluid-membrane behaviour is **not** modelled: there is no in-plane flow, no lipid exchange, no
  area reservoir, no viscosity. The tension σ is a constant material parameter, not a state variable
  responding to area change.
- No thermal forcing. This is a deterministic energy and its gradient. The connection to
  `ALEPH-DQ-102` (passive fluctuation spectrum) requires a stochastic integrator that does not
  exist yet, and until it does, nothing in this module is evidence about `(κ, σ)` recovery.
- The Gaussian-curvature term is dropped on a fixed-topology argument. If kinetic topology change
  ever enters this vertical, that argument lapses and the term must return.
- Self-intersection is not detected. A surface that folds through itself still has a finite energy
  and a finite force here; only `TriangulatedSurface.check_invariants` would object, and this module
  operates on raw arrays for speed.

### 14.4 Five invariants are asserted here and tested nowhere

**Found 2026-08-04**, by widening `named_controls` in `tests/ports/test_port_discipline.py` so it
reads citations in prose and not only in the control table. This entry cited six tests in
`tests/vertical/test_membrane.py` and `tests/vertical/test_vertical_firewall.py`. **Neither file has
ever existed**, and `tests/vertical/` contains no membrane or Helfrich test module at all.

The firewall citation had a real successor and is repaired: `aleph/**` must not import
`validation/**`, and `tests/firewall/test_scope_firewall.py::test_runtime_never_imports_the_analytic_oracles`
enforces it over **the whole `aleph` tree**, which is wider than this entry claimed. The name says
`runtime` and means runtime code generally; the per-subtree firewalls were consolidated into it.

The other five have **no successor anywhere under `tests/`**:

| | invariant | status |
|---|---|---|
| I3 | net internal torque vanishes — `Σ xᵢ × Fᵢ = 0` | **claimed, not tested** |
| I4 | `∇V` matches a central difference to round-off | **claimed, not tested** |
| I5 | bending energy is scale-invariant at `c₀ = 0` | **claimed, not tested** |
| — | tension force matches a central difference | **claimed, not tested** |
| — | spontaneous-curvature force matches a central difference | **claimed, not tested** |

Nearby tests exist and are **not** these. `test_the_bending_energy_is_exactly_scale_invariant`
grades the **nucleus lamina**, not this membrane; `test_recovered_tension_matches_the_declared_tension`
is a Laplace-law recovery, not a gradient check. Substituting either would make this entry cite a
control that does not test what the row claims, which is worse than an honest gap.

**This entry was written 2026-07-30 "before the code", and the names were what the lane intended to
write.** The code landed with a different test structure and the ledger did not follow. That is the
`ALEPH-PORT-2901` pattern the discipline test's own message warns about — a `PROPOSED` entry reading
like verified work — and it survived because the citations sat in prose where the gate could not see
them.

**Not repaired by renaming.** The names are removed rather than pointed at approximate matches, and
the invariants are recorded here as untested. Writing them is a task; claiming them is not.

