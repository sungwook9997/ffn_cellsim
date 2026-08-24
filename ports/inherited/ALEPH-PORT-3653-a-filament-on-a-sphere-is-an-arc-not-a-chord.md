# ALEPH-PORT-3653 — a filament on a sphere is an arc, not a chord

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3653` |
| Lane | `a40e55a2 S-OBSERVE handover, scope extended by PI reassignment` |
| Status | `PROPOSED` |
| Written | `2026-08-08` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Exists because | `ALEPH-PORT-3652` moved the filament from 10 nm to a sourced 120 nm, and the zero-thickness shell control went red. The PI, 2026-08-08 19:1x KST: *"노드를 꼭 그렇게 지정해놔야한다면 처음부터 구면 곡률에 맞게 필라멘트 곡률 넣어서 배열해놔도 되니깐"* — and that is the repair. |

---

> **NOT LANDED. Written, implemented, measured, and then reverted on the measurement.** The code
> below was built and run on 2026-08-08; it makes the shell exact to `1e-15 µm` and **moves six
> controls**, one of which is a statement about the physics rather than about a tolerance. Six
> controls need six adjudications and this session did not have them, so the change is recorded
> here at full detail — one step from landing — rather than merged with its reds open. See §7a.

## 1. Aleph API

`aleph/vertical/cortex_surface_coupling.py` — **proposed, currently reverted**:

- `_arc_nodes(midpoint, tangent, half, arc) -> np.ndarray` — filament node positions on the great
  circle through `midpoint`, replacing the tangent chord at both construction sites (the retained
  one-per-vertex control at `:1319` and the areal-density population at `:1365`).

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none — no provider file was read.** |
| Source commit | **not applicable**; nothing read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source path | **not applicable** |
| Source symbol(s) | **not applicable**; the Aleph symbol is named in §1 |
| Read from | **plane geometry**, plus this repository's own Cytosim harness — see §7 |
| Working tree == commit? | **not applicable** — no provider revision is cited |

## 3. Why source-derived porting beats clean-room

It does not; nothing is source-derived. Stated explicitly rather than left blank.

## 4. Physical or mathematical law represented

A filament of rest length `L` lying **on** a sphere of radius `R` subtends an angle `L/R`, and every
one of its points is at radius `R`. A filament laid along the **tangent** is a chord: its midpoint
is at `R` and its ends are at `√(R² + (L/2)²)`, i.e.

```
    bulge  ≈  L² / (8R)
```

which is a construction artefact, not a physical displacement. It is also **quadratic in L**, which
is why it was invisible for as long as the filament was wrong:

| filament | bulge at R = 4.98 µm |
|---:|---:|
| 0.010 µm (shipped until 2026-08-08) | 2.5e-06 µm = **0.0025 nm** |
| 0.060 µm | 9.0e-05 µm = 0.090 nm |
| **0.120 µm (sourced, `-3652`)** | **3.6e-04 µm = 0.36 nm** |

The zero-thickness control allows `1e-5 µm`. The sourced filament exceeds it by **36×**; the old
one cleared it by 4×. **The chord was always wrong and only the wrong filament length hid it.**

**Rest arclength is preserved exactly.** The subtended angle is `s/R`, so the arclength between two
nodes is `R·Δθ = Δs`. Nothing that reads a material coordinate changes meaning.

## 5. Units, domains, singular cases, invariants

- µm throughout. `R = ‖midpoint‖ > 0` — a midpoint at the origin has no radial direction, and the
  callers place midpoints on a shell, so this cannot arise from them.
- **I1.** Every node of every filament sits at `‖midpoint‖` to float64 round-off.
- **I2.** Node-to-node rest arclength is unchanged from the chord construction by construction.
- **I3.** As `L/R → 0` the arc converges to the chord, so this is a refinement of the old behaviour
  and not a different model.

## 6. Source evidence class and known retractions

Geometry. Nothing to retract. The bulge formula is stated exactly rather than to leading order in
the code; `L²/(8R)` appears in this entry only as the reader's sanity check.

## 7. Independent oracle or derivation

**Two, and the second is this repository's own external oracle.**

**(a) The shell closes.** `‖node − centre‖` must be one number for every node. Measured after:

| | radial spread | tolerance |
|---|---:|---:|
| `subdivision_level=1`, 126 nodes | **8.882e-16 µm** | 1e-5 |
| `subdivision_level=2`, 486 nodes | **3.553e-15 µm** | 1e-5 |

That is float64 round-off — **ten orders of magnitude inside the tolerance**, where the chord at the
old 10 nm filament sat at 2.5e-6, only 4× inside it.

**(b) It is the configuration Cytosim already grades Aleph against.**
`validation/independent/cytosim.py::arc_configuration` builds *"`node_count` points on a circular
arc"*, and its docstring gives the reason in its own words: *"one of this lane's surveys measured
that a naive 'L-bend' gives an end-node force ratio of exactly 1.0 when the kink is at an interior
node and 1.5 only when it is at node 1 or `n-2`, so a comparison that does not state where the bend
is measures nothing. **A uniform arc has curvature everywhere and has no such ambiguity.**"

So `ALEPH-PORT-3639`/`-3642` graded Aleph's bending law on an arc and agreed with Cytosim to
**0.016%** — while the cortex those laws run on was built from chords. This entry makes the
construction the shape the oracle certified.

## 7a. What it moves, measured — and why it was reverted rather than landed

Implemented, `tests/vertical` run twice to separate this change from `ALEPH-PORT-3652`'s default:

| change | reds | which |
|---|---:|---|
| filament default `0.005 → 0.060` alone | **3** | `test_descent_tolerance` ×2, `test_numerical_invariance` ×1 |
| **arc placement alone** | **6** | `test_cortex_crosslink_topology` ×4, `test_cortex_population_density` ×1, `test_vertical_controls::test_a_bidirectional_tether_is_caught_carrying_compression` ×1 |
| both together | 9 | the union; they do not interact |

**One of the six is not a tolerance, and it is the reason this is a decision rather than a patch.**
`test_the_population_is_built_at_rest_under_both_placements` fails, and it is **right to fail**.
`segment_rest_um` is computed from the actual node separations, so an arc is at **axial** rest —
the axial term is not what moved. What moved is **bending**: a chord through three nodes is
straight and carries zero curvature, while an arc carries the surface's curvature and therefore a
non-zero bending energy at construction.

**That energy is physical.** A cortical filament lying on a curved cell surface *is* bent by that
curvature, and the corresponding prestress is real, not an artefact. But it means the cortex is no
longer *"built at rest"*, which is a property several controls and at least one published number
depend on. Whether that property should be kept is a modelling judgement about what the
construction is supposed to represent, and it is not a judgement this session was in a position to
make six times in a row.

**So the change is reverted in code and preserved here in full.** Everything needed to land it is in
§1 and §4; what is missing is six adjudications, each of which asks the same question — *is this
control's premise still true once a filament follows the surface it lies on?*

## 8. Positive control

**One exists and is what caught this. The rest are described, deliberately not named**, because the
change is reverted (see the note under the title) and **naming a test nobody has written is the
failure `PLAN` §0.2.4 exists over** — `tests/ports::test_named_controls_resolve_to_real_tests`
caught this entry doing exactly that, on 2026-08-09, and the names are removed rather than the guard
argued with.

- **Exists:**
  `tests/vertical/test_cortex_population_density.py::test_without_a_thickness_every_cortical_node_sits_at_one_radius`
  — the control that found the defect. With the arc it passes at `1e-15 µm`; with the chord at the
  sourced filament it fails at `3.6e-4`, and with the chord at the old 10 nm filament it passed at
  `2.5e-6`, four times inside a tolerance it now clears by ten orders.
- **To be written when this lands:** an invariant that node-to-node rest arclength is unchanged
  between the chord and arc constructions, which is what makes this a refinement rather than a new
  model.

## 9. Deliberately failing negative control

Also described rather than named, for the same reason. When this lands it needs two:

- the chord bulge computed directly at the sourced filament, asserted to exceed the shell control's
  tolerance by more than 30× — so the control that caught this is shown able to catch it;
- the same computation at 0.010 µm clearing that tolerance, which is **why the defect survived**
  until `ALEPH-PORT-3652` lengthened the filament.

Both are arithmetic on `L²/8R` and neither needs a build, so neither is expensive; they are absent
because the change is absent, not because they are hard.

## 10. Numerical and precision envelope

`cos`/`sin` on an angle of order `L/R ≈ 0.024` rad at the sourced filament — far from any branch or
cancellation. The measured residual radial spread is `~1e-15 µm` on a 4.98 µm radius, i.e. `~2e-16`
relative, which is one ULP.

## 11. Production-backend residency and transfer

**None.** Build-time host geometry. It changes node *positions*, not array shapes, dtypes, or the
transfer. The existing host/device parity controls remain the check, and they compare positions, so
they will see this change and must be re-run — which is what §13 requires.

## 12. Comments and docstrings to discard

Nothing copied. The Cytosim docstring in §7(b) is quoted for support and cited, not inherited.

## 13. Acceptance

`PROPOSED`. `AUDITED` when §8 and §9 pass **and** the host/device parity controls have been re-run,
because this moves every cortical node position and those controls compare positions bitwise.

## 14. Honest limits

- **It moves every cortical node.** Every cortex number in this repository was taken on chords. The
  displacement is sub-nanometre, but "small" is not "zero" and the parity controls compare bits.
- **It does not make the filament flexible.** The nodes are placed on an arc; the filament is still
  a polyline between them and the segments are still straight. This fixes where the nodes *start*,
  not what shape the filament can take.
- **The shell is a sphere here.** `_arc_nodes` uses `‖midpoint‖` as the local radius, which is exact
  for a sphere and an approximation for any other surface. Nothing in this repository builds a
  cortex on a non-spherical surface today; when something does, this is where it will be wrong.
- **The repair route was scoped on 2026-08-09 and is larger than it looks.** The one red that is a
  statement about physics — the cortex no longer being *"built at rest"* — dissolves if the filament
  bending law gains an intrinsic curvature, so an arc is at rest **in its curved state**. But
  `bending_energy_and_forces`'s own docstring records that it is *"pinned numerically by a control"*
  to `aleph.vertical.sf_arc`'s law, *"so Aleph has one bending convention rather than two that
  disagree"* — `tests/vertical/test_cortex_filament_controls.py::test_this_bending_law_agrees_with_the_sf_arc_owners_bending_law`.
  **`sf_arc.py` is not this session's path**, and the term also has to reach
  `aleph/runtime/cortex_relax_kernels.py` with its parity re-run. An *optional* term defaulting to
  zero keeps the pin, but then the arc still cannot land without opting in, and opting in needs the
  device kernel — so the optional form does not actually unblock this entry. **Reported and not
  attempted.**
- **It does not address the placement of the other six compartments**
  (`PROPOSAL-only-the-cortex-has-a-resolution.md` §4). Those are blobs, and a blob's filaments being
  arcs does not un-blob it.
