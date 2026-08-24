# ALEPH-PORT-3620 — the moving-boundary coupling, between two endpoints that do not share a geometry

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3620` |
| Lane | `46143f30` — CUDA track |
| Status | `PROPOSED` |
| Written | `2026-08-04` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Exists because | `ALEPH-PORT-3617` §1 says it in as many words: it ports *"that operator and its two gradings — **not the connectors**"*. The operator exists and is graded on an RTX 4090 (`-3619`). `membrane_cytosol_boundary` and `nucleus_cytosol_boundary` are still `physics_absent`, and this entry is the piece between the operator and them. |
| Blocks | `wiring.py` 30/36 → 32/36, **if and only if §14 is satisfied** |

---

## 1. Aleph API

New surface, all of it under `aleph/vertical/`:

```
aleph/vertical/connectors_boundary.py     MovingBoundaryCoupling  (new module)
aleph/vertical/wiring.py                  two Binding rows, in the same commit as the class
```

**No change to `aleph/runtime/`.** The law is `-3617`'s and is called, not re-implemented.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository identity | `/Users/sw1/ffn_cellsim` — **READ ONLY**, and read only through the audit already committed at `.superpowers/sdd/BUILD_PLAN-2026-07-31-engine-gaps/source-audit-bc.md` |
| Source path | `ac/cell/fsi_coupling.py` (the two-directional coupling driver), `ac/fluid/ibm_reference.py` (its host oracle) |

## 3. Why source-derived porting beats clean-room

The reference project has run this coupling for years and its **failures** are the valuable part —
a clean-room implementation would rediscover them. Two are already inherited into `-3617` and
graded there: the boundary divergence that is a factor of two low behind interior-only gates, and
the host oracle that carries no mask while the production path it certifies runs the masked form.

## 4. The law, and the thing that is actually hard

### 4.1 The two directions

    solid -> fluid   the surface's own normal motion becomes a prescribed normal flux
                     (the KINEMATIC condition: the fluid may not pass through a moving wall)
    fluid -> solid   the fluid's pressure traction becomes a force on the surface's vertices

Both endpoints already exist and neither is the gap:

| half | endpoint | published by | supplies |
|---|---|---|---|
| fluid | `BoundaryFaceStencil` | `CytosolField.boundary_endpoint` | `outward_normals()`, `prescribe_normal_flux(q)`, `pressure_traction_pn()` |
| solid | `NuclearEnvelopeBoundary` | `NucleusOwner.boundary_endpoint` | `face_centroids_um()`, `face_areas_um2()`, `outward_normals()`, `normal_velocity_um_per_s()` |

### 4.2 **The two face sets are not the same faces, and this is the whole entry**

Measured in the source, not assumed:

- `BoundaryFaceStencil.outward_normals` is built from `faces.boundary_axis` and
  `faces.boundary_sign`. Every fluid normal is exactly `±e_x`, `±e_y` or `±e_z`. These are
  **axis-aligned Cartesian grid faces**.
- `NuclearEnvelopeBoundary.outward_normals` comes from the triangulated envelope. These are
  **Lagrangian triangles with arbitrary normals**, and the surface moves.

So the connector is not a relabelling between two views of one surface. **There is no face
correspondence to find**, and any implementation that pairs face `i` with face `i` is wrong in a way
that would look right on a sphere-in-a-cube fixture.

`-3617`'s masked immersed-boundary pair is the operator for exactly this. **But it maps Lagrangian
nodes to cell CENTRES over a 4³ block, and the fluid endpoint here speaks in FACES.** That
mismatch is real, it is not resolved by this section, and §14.1 records it as the open question
rather than papering over it.

### 4.3 The identity the coupling must satisfy

Whatever the geometric bridge turns out to be, the pair must remain **adjoint** in the same sense
`-3617` §5 establishes:

    the work the solid does on the fluid  ==  the work the fluid does on the solid,
    to the sign, over one exchange

If it does not, the coupled system has a numerical energy source that no owner's own energy books
can see — each side balances and the pair does not. That is the failure this entry exists to make
impossible, and it is why §9's negative control is an energy control rather than a force comparison.

## 5. Units, domains, singular cases, invariants

| quantity | unit |
|---|---|
| prescribed normal flux | µm/s, **positive out of the fluid** |
| pressure traction | pN per face |
| face area | µm² |
| normal velocity | µm/s, **positive along the surface's outward normal** |

**The two sign conventions above point in opposite senses and that is the first place this can go
wrong.** `prescribe_normal_flux` is *positive out of the fluid*; `normal_velocity_um_per_s` is
*positive along the solid's outward normal*, which the nucleus docstring says is "a volume being
taken away from" the fluid. They agree, but only after reading both docstrings, and a control
asserts the sign rather than a comment claiming it.

**Singular cases:**

- **A surface at rest.** Every prescribed flux must be **exactly `0.0`**, not small. `-3617`'s
  cytosol docstring makes that exactness a census statement: the cross-boundary volumetric flux is
  zero, so any nuclear volume change comes from the internal volume constraint and never from fluid
  crossing the surface. A coupling that leaked `1e-18` would make that claim false.
- **A uniform pressure field.** The resultant traction on a closed surface is zero by construction,
  so a control that reads the resultant passes on a coupling that transfers nothing. Per-face, per-
  vertex, always.
- **A surface partly outside the fluid domain.** Undefined here and refused rather than clamped.

## 6. Source evidence class and known retractions

**Evidence class: SOURCE-DERIVED** for the coupling's shape, **RE-DERIVED** for its arithmetic.

**Retraction carried forward and repeated because it is easy to re-acquire:** an earlier revision of
`BUILD_PLAN-2026-07-31-engine-gaps.md` said the reference stubs the `fluid ← solid` direction to
zero and carries one importer. **All three parts were false.** There are seven importers and they
run every outer step. Do not build against the stub story.

## 7. Independent oracle or derivation

The **work identity** of §4.3, computed on the host from the two endpoints' own published
quantities, with no appeal to the connector's internals. It is independent in the sense that
matters: a connector that computed both sides of its own identity would be certifying itself.

## 8. Positive control

`tests/vertical/test_connectors_boundary.py` — **to be written with the class, in the same commit.**

- a surface at rest prescribes **exactly `0.0`** on every fluid face
- a surface advancing uniformly prescribes a flux whose sign is *out of the fluid*
- the work identity of §4.3 holds to its stated tolerance
- the traction reaches **vertices**, and the per-face constituent scale is non-zero when the
  resultant is zero

## 9. Deliberately failing negative control

At least four mutants, each aimed at a failure the obvious control would miss:

| mutant | why the obvious control misses it |
|---|---|
| `WRONG_FLUX_SIGN` | the resultant work still balances in magnitude; only the sign of the exchange is wrong, and an energy *magnitude* check passes |
| `WRONG_FACE_PAIRING` (index `i` ↔ index `i`) | **passes on any fixture where the two face sets happen to be ordered alike** — which is why §8's fixture must not be one |
| `WRONG_RESULTANT_TRACTION` (net instead of per-face) | exactly zero net on a closed surface, so a resultant control cannot see it |
| `WRONG_UNSCALED_FLUX` (velocity used as flux without the area weight) | dimensionally plausible and correct at unit area, which a unit-sphere fixture supplies |

## 10. Numerical and precision envelope

To be measured. The prescribed-flux-at-rest control is an **exact zero** and takes no tolerance.
The work identity takes a relative tolerance to be stated when measured, not chosen in advance.

## 11. Production-backend residency and transfer

Deferred and stated as deferred. `-3617`'s kernels are device-resident; this connector's first
version runs on the host, because the geometric bridge in §4.2 has to be settled before it is worth
putting on a device. **No claim of device residency is made for this entry.**

## 12. Comments and docstrings to discard

Every comment and docstring in `ac/cell/fsi_coupling.py` and `ac/fluid/ibm_reference.py`, per
`PLAN.md` §0.2. The mathematics is re-derived and the prose is new.

## 13. Acceptance

Accepted when §8 and §9 pass **and** `wiring.py` reports 32/36 with both `Binding` rows resolving.
Until then this entry is `PROPOSED` and the connectors stay `physics_absent`.

## 14. Honest limits

> ### ⚠ AMENDED 2026-08-04 15:1x — §4.2 was one level too low
>
> §4.2 asks how to transfer between two face sets. **That question presumes the fluid's face set is
> already determined, and it is not.** The fluid's boundary faces are derived:
>
> ```python
> faces = np.flatnonzero(self.faces.boundary_neighbour == CellClass.NUCLEUS_INCLUSION)
> ```
>
> — from `cell_class`, a **per-cell voxel classification** of which control volumes are nucleus.
> **Nothing under `aleph/` computes it from the envelope.** It arrives as a constructor argument and
> is never recomputed, while `CellClass`'s own docstring says the opposite is required: *"the
> membrane and the nuclear envelope move, so which volumes are interior fluid changes during the
> run. A reclassification is therefore a topology change that commits with the step."* The
> declaration is ahead of the code, which is the same shape `wiring.py` reports for the connectors.
>
> So the work divides in two, and the second cannot start before the first:
>
> | | question | state |
> |---|---|---|
> | **1** | which cells does the moving envelope occupy — i.e. **what is `cell_class`** | **no code at all** |
> | **2** | how do flux and traction cross between those faces and the triangles | `-3617`'s operator, aimed at cell centres |
>
> **This retracts my own recommendation.** I proposed option A — extend `-3617`'s operator to face
> centres — having looked only at question 2. A does not touch question 1, and neither does B.

### 14.1 The classification must change. `PI decision, 2026-08-04T06:16:34.433Z`

> *"시나리오가 고정 마스크를 주는 방식은 말도 안된다고 생각함. 매 스텝까지는 아니더라도 변화하게
> 해야되는 거 아님?"*
> `~/.claude/projects/-Users-sw1-Project-Aleph/46143f30-….jsonl#4cf0cbef-f036-4391-b57d-a5f7835054a6`

A scenario-supplied fixed mask is refused. The open question is no longer *whether* the
classification changes but **when**, and that is a design question with a sound answer rather than a
taste.

> ### ⚠ REFUTED BY MEASUREMENT, 2026-08-04 15:2x — the bound below does not work
>
> §14.2.5 said the margin had never been computed and that it was a measurement to take **before**
> building. It was taken, and it says **do not build this.**

**The proposal was: reclassify on a bound, not on a schedule.**

The surface is affine between its vertices, so no point of it moves further in a step than its
furthest vertex. Accumulate `D = Σ max_v |Δx_v|` since the last classification, record the margin
`m = min over near-surface cell centres of |distance to the surface|`, and reclassify when `D ≥ m`.
Sound by construction: a class change requires crossing a cell centre, and the surface cannot have
reached the nearest one while `D < m`.

**Measured on a real envelope — icosphere, radius 3 µm — against realistic grids:**

| mesh | `dx` | near-surface cells | `m` |
|---|---|---|---|
| level 2 (162 v, 320 f) | 1.00 | 448 | 0.0018 µm = **0.0018 dx** |
| level 2 | 0.50 | 1 846 | **0.000000** |
| level 2 | 0.25 | 7 290 | **0.000000** |
| level 3 (642 v, 1280 f) | 1.00 | 448 | 0.034 µm = 0.034 dx |
| level 3 | 0.50 | 1 846 | **0.000000** |

**`m` is zero at every realistic spacing, so the trigger fires on the first step and the bound
degenerates to "every step".** Distance computed point-to-*triangle*, not point-to-vertex — the
vertex form overestimates and would have flattered the design.

**Why, and it is obvious afterwards:** `m` is a **minimum over thousands of samples** of a quantity
spread roughly uniformly over `[0, dx/2]`. With 7 290 near-surface cells the smallest is
indistinguishable from zero. The flaw is the `min`: it is set by whichever single cell happens to lie
nearest the surface, and I never asked whether one cell flipping is a problem worth an early
reclassification.

### 14.1a So the question changes: is reclassifying every step actually expensive?

Measured, host NumPy, exact solid-angle winding over a closed surface — no ray casting, no
tie-breaking:

| grid | cells | × faces | work | wall |
|---|---|---|---|---|
| 15³ | 3 375 | 320 | 1.1 M | **205 ms** |
| 29³ | 24 389 | 320 | 7.8 M | **938 ms** |

≈ **5–8 M triangle-point evaluations per second on the host.** Too slow to run every step: a step is
milliseconds and this is hundreds.

Restricting to cells the surface could reach buys **2–3×**, not the two orders needed —
1 846/3 375 (55 %) and 7 290/24 389 (30 %).

**But the shape of the work is `cells × faces` with a per-cell reduction, which is the shape of
every kernel already in `law_kernels.py`.** 7.8 M evaluations is small for an RTX 4090. So the
answer this measurement points at is **not a cleverer trigger — it is putting the classification on
a kernel and reclassifying every step because it is then cheap.** That is also the track this
project is already on, which is a reason to believe it rather than a reason to be suspicious of it.

**UNVERIFIED:** no kernel has been written and no device figure exists. "Small for a 4090" is an
expectation, and §14.2 keeps it as one.

**One thing the measurement corrected on the way.** I compared the winding classification against
`|centre| < radius` and called the analytic sphere "exact". It is not: the icosphere is a
**polyhedron**, and the 96 of 24 389 cells (0.39 %) where they disagree lie between the polyhedron
and the sphere. The winding test was right about the surface it was given and my reference was the
wrong surface.

### 14.1b The volume is defined twice and nothing reconciles the two

Raised by the PI, 2026-08-04 15:2x, and it sits **above** everything else in this entry.

| side | definition | behaviour |
|---|---|---|
| solid | `enclosed_volume_um3(vertices, triangles)` — the divergence theorem over the closed surface. Registered as the nucleus's **owned state** and used by its volume constraint. | **continuous** in the vertex positions |
| fluid | `cell_volume_um3 = dx**dim` — the volume of **one grid control volume**. Nothing counts `INTERIOR_FLUID` cells to form a total; searched and absent. | **piecewise constant**, changing only when the classification changes |

So the envelope's volume varies smoothly while the fluid's voxel count does not move at all — and
then, when one control volume flips, it jumps by `dx³` at once. At `dx = 0.5 µm` that is
**0.125 µm³ against a 3 µm nucleus's 113 µm³ — roughly 0.1 % in a single step**, from no fluid
motion whatever.

**Why this is not a rounding detail.** `CellClass.NUCLEUS_INCLUSION`'s own docstring makes an
explicit modelling claim: the envelope is impermeable, so *"the relative flux across a fluid/nucleus
face is exactly zero rather than small — the census registers that as a modelling statement with
teeth"*. A reclassification moves `dx³` of volume between the fluid and the nucleus **with no flux
at all**, which contradicts the claim it is supposed to protect. Volume would be teleporting.

**This reframes §14.1.** I asked *when to reclassify* and treated it as a cost question. The
question underneath is **how much the two volume definitions are allowed to disagree, and where the
difference is booked** — and no reclassification schedule, however clever, answers it.

**Two questions for the PI, neither of them mine to settle:**

1. **Which volume is authoritative?** Is the closed-surface volume the truth and the voxel count an
   approximation of it, or is the grid the fluid's ledger and the surface volume a separate quantity
   that need not agree? — **ANSWERED, see §14.1b-R.**
2. **Where is the `dx³` jump booked at a reclassification?** Absorbed into the nucleus's volume
   constraint, carried as an explicit correction flux, or declared as an accepted discretisation
   error with a stated bound. **Leaving it unbooked makes the "exactly zero" statement false**, and
   that statement is load-bearing in the census. — **ANSWERED by the same ruling: the third of these
   three, and §14.1b-M states the bound.**

~~**UNVERIFIED:** the 0.1 % figure is arithmetic from `dx` and a nominal radius, not a measurement of
a running configuration.~~ **Measured 2026-08-04, §14.1b-M.** The arithmetic survived — `0.110584 %`
against the icosphere's own volume where the estimate said "roughly 0.1 %" — but it was the wrong
quantity to have been reassured by, and the measurement is what shows that. §14.1b-M.

### 14.1b-R RULED: the solid surface volume is authoritative

**This records a PI ruling and carries no PI authority itself. No `decided_by` field, deliberately;
the coordinates are what make the quote checkable by a session that did not receive it.**

| Field | Value |
|---|---|
| Transcript | `~/.claude/projects/-Users-sw1-Project-Aleph/46143f30-d1ab-4f8a-9c66-b63204a3b7a4.jsonl` |
| Record uuid | `#f29d6731-83a3-40cd-9b90-36dfa26c5427` |
| Timestamp | `2026-08-04T08:21:25.021Z` |
| Verbatim | *"고체 표면 부피가 정본"* |

Read out: **the solid surface's volume is the authoritative text.** So `enclosed_volume_um3(vertices,
triangles)` — already the nucleus's registered owned state and already what its volume constraint
acts on — **is** the nucleus's volume, at every instant and continuously in the vertex positions.
The count of `NUCLEUS_INCLUSION` control volumes is a **rendering of that volume onto the grid**, not
a second ledger of it.

**The consequence, which is the whole reason the ruling matters.** A reclassification transfers
nothing. The nucleus's volume did not change by `dx³`; the *grid's approximation* of the nucleus's
volume changed by `dx³`. That is a **discretisation error with a stated bound** — §14.1b-M states
it — and not a flux, not a source, not a quantity that has to be conserved against anything.

**And so the census statement is true again, exactly as written.** `CellClass.NUCLEUS_INCLUSION`'s
docstring says *"the relative flux across a fluid/nucleus face is exactly zero rather than small"*.
That claim is about a **face**. A reclassification is not a face event: no fluid crosses any face
when a cell changes class, and every fluid/nucleus face still carries the prescribed normal flux
that `FaceTable` defaults to exactly `0.0`. The statement was only ever endangered by reading the
voxel count as a volume ledger — by treating `dx³` as something that *went* somewhere. The ruling
says the voxel count is not a ledger, so nothing went anywhere, and **"exactly zero" is exactly zero
and not a tolerance**.

**What the ruling does not settle, named rather than assumed.** It is about *volume*. A cell that
leaves `INTERIOR_FLUID` also takes a **pressure degree of freedom** with it — `CytosolField`
carries one `_pressure` entry per fluid cell and a storage term that scales with `cell_volume_um3` —
and whether *that* is likewise a pure approximation artefact or a second unbooked quantity is a
different question. It is not measured here and is added to §14.2.

### 14.1b-M The bound, measured

`scripts/measure_control_volume_fraction.py`, host NumPy, CPU. The authoritative volume is read from
`aleph.vertical.membrane.enclosed_volume_um3` on the icosphere the nucleus is actually built from —
not from `4/3 π R³` — so this bounds the object the ruling is about. Classification is exact
solid-angle winding (Van Oosterom–Strackee), the same rule §14.1a used. A voxel count depends on
where the grid sits relative to the body, so every count below is a range over **12 grid offsets**:
five shifted equally on all three axes and seven shifted independently per axis, because a shift
common to all three is a symmetric special case and would have understated the spread.

The authoritative volume is itself level-dependent, and reporting it as "113 µm³" hid that:

| icosphere | faces | `V_surface` [µm³] | vs. the smooth ball `4/3 π 3³` |
|---|---|---|---|
| level 2 | 320 | 109.270206359 | −3.384 % |
| level 3 | 1 280 | 112.124002062 | −0.861 % |
| level 4 | 5 120 | 112.852951596 | −0.216 % |
| level 5 | 20 480 | **113.036173631** | −0.054 % |

The bound is taken on level 5. **A 3 µm nucleus at level 2 — the mesh §14.1a timed — is already
3.4 % smaller than the sphere it is drawn as, which is thirty times the single-cell jump this
section was worried about.**

**The answer to the question §14.1b asked** — what fraction of the nucleus's volume is one control
volume:

| `dx` [µm] | one control volume [µm³] | as a fraction of `V_surface` |
|---|---|---|
| 0.50 | 0.125 | `1.105840688e-03` = **0.110584 %** |
| 0.25 | 0.015625 | `1.382300860e-04` = **0.013823 %** |

Halving `dx` divides one reclassification's share by exactly 8, as `dx³` requires. **The estimate at
the head of this section was right**: "roughly 0.1 %" against a measured 0.110584 %.

**And that is the number this section should not have stopped at.** Two more, from the same run:

| `dx` [µm] | standing gap `(N·dx³ − V)/V`, measured over 12 offsets | rigorous cut-cell bound |
|---|---|---|
| 0.50 | −1.027 % … +0.963 % (`N` = 895…913) | **74.534 % … 77.851 %** |
| 0.25 | −1.539 % … +0.341 % (`N` = 7 123…7 259) | **36.935 % … 37.654 %** |

* The **standing gap** is how far the two definitions already disagree *at rest*, before anything
  moves. It is **an order of magnitude larger than one reclassification** at both spacings, and
  halving `dx` did not shrink it — the spread is 1.99 % wide at `dx = 0.5` and 1.88 % wide at
  `dx = 0.25`. A single-jump figure of 0.11 % therefore bounds nothing about the agreement of the
  two volumes; it bounds only the size of one step of the error.
* The **cut-cell bound** counts the cells the surface actually passes through — for a sphere,
  exactly those whose cube has its nearest point inside and its farthest point outside. Every other
  cell is unambiguously in or out, so the entire disagreement lives in these cells and their total
  volume is the rigorous bound. At `dx = 0.5` **three quarters of the nucleus's volume sits in cells
  the surface cuts.** Twelve cells across a diameter is not a discretisation of a sphere; it is
  mostly boundary layer. **The ≈1 % agreement measured above is cancellation, not accuracy**, and a
  bound derived from it would be a bound derived from luck.

**So the honest statement of the discretisation error the ruling now licenses is: `0.11 %` (`dx =
0.5`) or `0.014 %` (`dx = 0.25`) per reclassification event, against a standing disagreement of
about `1 %` and a rigorous worst case of `75 %` and `37 %` respectively.** Only the first of those
three is small. Recording the first alone would have been the same mistake as recording "roughly
0.1 %" and moving on.

**What this does not measure.** No running configuration was stepped: this is a static nucleus at
its reference radius, and the accumulation of many reclassification events over a trajectory is not
bounded here — only the size of one, and the worst case they could reach. A drift measurement needs
the coupling to exist, which is what the rest of this entry is for.

### 14.2 Everything else still open

1. **Question 2 remains as §4.2 states it**, and the choice between extending `-3617`'s operator,
   publishing a cell-centred fluid endpoint, or an explicit triangle↔face intersection is unmade.
   It should be decided **after** question 1, because the classification determines what the face
   set even is. **§14.1b's question 1 is now answered (§14.1b-R), so this precondition is met and
   this item is unblocked.**
2. **`membrane_cytosol_boundary` may not be servable by the same class.** Only the nucleus endpoint
   was read. Serving both from one class is an expectation, not a measurement.
3. **No device claim.** §11.
4. **The work identity's tolerance is unstated** because it is unmeasured. A number chosen before
   the measurement is a number the measurement cannot contradict.
5. ~~The margin `m` has never been computed.~~ **Computed 2026-08-04, and it refuted the design.**
   See the block at the head of §14.1. The predicted failure mode — *"`m` is routinely a fraction of
   a step's motion, so the bound degenerates to every step and buys nothing"* — is exactly what
   happened, except that `m` is not a fraction of the motion but **zero**.
6. **The kernel classification is UNVERIFIED.** §14.1a argues it from the shape of the work and one
   host throughput figure. No kernel exists, no device measurement has been taken, and "small for a
   4090" is an expectation. It is the next thing to measure, not the next thing to assume.
7. **The classification rule itself is unchosen.** Solid-angle winding was used to take the timing
   because it is exact for a closed surface and vectorises; that is not a decision that it is the
   rule. A volume-fraction rule would change both the cost and what "a cell changed class" means.
   **§14.1b-M sharpens the stake:** at `dx = 0.5` three quarters of the nucleus's volume is in cells
   the surface cuts, and a volume-fraction rule is exactly a rule about those cells.
8. **The reclassified cell's pressure degree of freedom is unbooked, and §14.1b-R does not book it.**
   The ruling is about volume. `CytosolField` carries one `_pressure` entry per fluid cell and a
   storage term scaling with `cell_volume_um3`; a cell leaving `INTERIOR_FLUID` takes that entry
   with it. Whether that is a second approximation artefact or an actual unbooked quantity is
   **unmeasured** — it is named here rather than assumed benign, because assuming it benign is the
   same move §14.1b was written to catch.
9. **The `dx = 0.5` grid may simply be too coarse for a 3 µm nucleus**, and nothing in this entry has
   asked that question. §14.1b-M measures a rigorous cut-cell bound of 74.5–77.9 % there. That is not
   a reason to reject the discretisation on its own — the error may still cancel in the quantities
   the coupling actually transfers — but "may still cancel" is a hypothesis and no measurement here
   supports it.
