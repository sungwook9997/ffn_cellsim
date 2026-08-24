# ALEPH-PORT-3630 — the fluid grid has no origin, and the fixtures are not in one frame

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3630` |
| Lane | `46143f30` Lane W2 (lead) |
| Status | `PROPOSED` |
| Written | `2026-08-04` — **before the code**, per `CLAUDE.md` §3 |
| Port class | `NOT A PORT` — nothing is taken from `/Users/sw1/ffn_cellsim`; no law changes and no physics is derived |
| Aleph target | `aleph/vertical/cytosol.py` (one new defaulted constructor argument), `aleph/scenarios/whole_cell.py` (the fixture placement) |
| Exists because | The nine `cytosol` coupling rows cannot be wired honestly. The fluid grid covers one eighth of the cell, in a corner, and `CytosolField` has no way to be placed anywhere else. |
| What this entry does **not** claim | That the resulting geometry is anatomical. It is not; `FIXTURE_SOURCE` still reads *"no primary source was read"*. This makes the fixture **self-consistent**, which is a different and lower claim than making it right. |

---

## 1. The measurement

`build_whole_cell()` places its compartments in two unrelated frames:

| owner | bounding box [µm] |
|---|---|
| `membrane` | `[-5, -5, -5] .. [5, 5, 5]` |
| `cortex` | `[-4.98, …] .. [4.98, …]` |
| `nucleus` | `[-3, -3, -3] .. [3, 3, 3]` |
| `sf_arc` | `[0, 0, 0] .. [6, 6.5, 1.1]` |
| `ecm` | `[0, 0, 0] .. [6, 6, 0.02]` |
| `microtubule` | `[0, 0, -0.2] .. [3.26, 4.07, 1.07]` |
| `lamellipodium`, `filopodium`, `intermediate_filament` | first octant, ≲ 3 µm |
| **`cytosol` grid** | **`(8,8,8) × dx 0.5` = `[0, 4]³`** |

The membrane is a 10 µm sphere about the origin. **The fluid occupies a 4 µm cube in one octant of
it** — roughly one eighth of the cell's bounding box, and none of the other seven.

Measured consequence, independently reproduced while mapping the transfer rows:

- **6 of `sf_arc`'s 11 nodes are refused** by `CytosolField.transfer_endpoint`. The refusals are
  exact — `(6.0, 1.0, 0.4)`, `(5.5, 3.0, 1.1)`, `(1.0, 6.0, 0.3)`, `(3.0, 6.25, 0.6)`,
  `(5.0, 6.5, 0.9)` — and the owner is right to refuse them: *"an immersed body there would read
  nothing and return its drag to nobody."*
- **`nmii` sits at `x ∈ [-0.15, 0.15]`, which is outside the domain**, and survives only because
  `_multilinear_weights` drops the out-of-range index and renormalises what is left. It is accepted
  by the clipping rule rather than by being in the fluid — a distinction that would never appear in
  any output.

## 2. Why this is a fixture defect and not a physics decision

`FIXTURE_SOURCE` reads *"UNVERIFIED whole-cell assembly fixture; no primary source was read"*, and
the module docstring says the geometry *"is where the constructors put it rather than anywhere
anatomical"*. Nothing about the current placement is a claim. Two constructors defaulted to the
origin and two others defaulted to the first octant, and no test ever compared them.

**No PI decision is required to make the fluid contain the cell.** What would require one is
choosing *anatomical* positions, and this entry does not do that.

## 3. What is actually missing in the code

`CytosolField` has no origin. Its axes are built at `cytosol.py:408` as

    axes = [(np.arange(n) + 0.5) * dx for n in shape]

so the domain is always `[0, n·dx]` on every axis. A grid that can only start at the origin cannot
contain a body centred there, at any resolution: growing `shape` extends it in `+x, +y, +z` only.

The point-to-cell map has the same assumption in two more places — `cytosol.py:2102`
(`base = floor(point/dx − 0.5)`) and `:2113` (`centre = (index + 0.5)·dx`).

## 4. The change

**One defaulted constructor argument**, `origin_um`, default `(0, 0, 0)`, threaded through those
three expressions:

    axes   = (arange(n) + 0.5) * dx + origin[axis]
    base   = floor((point - origin) / dx - 0.5)
    centre = (index + 0.5) * dx + origin

**Every existing call is bit-identical**, because adding `0.0` to a float is exact. That is the
whole argument for why this may touch a module that is the parity oracle for 109 kernels: the
default path's arithmetic does not change, it gains a term that is exactly zero.

The scenario then places the fluid on the cell: `shape=(24,24,24)`, `dx=0.5`,
`origin_um=(-6, -6, -6)` → `[-6, 6]³`, which contains the 10 µm membrane with a 1 µm margin, and the
nucleus inclusion moves with it.

## 5. Cost

13,824 cells against 512. That is affordable **only because `ALEPH-PORT-3625` landed today**: the
dense `np.linalg.solve` it replaced was `O(n³)` time and `O(n²)` memory, and 13,824 cells of coupled
Biot unknowns would have been a matrix of tens of gigabytes. The matrix-free Krylov solve makes the
domain a resolution choice rather than a wall.

Recorded because the ordering matters: **this entry was not possible yesterday.**

## 6. Units, domains, singular cases, invariants

Units: length µm. `origin_um` is a position, not an offset in cells, so a caller cannot express a
half-cell shift by accident.

**Invariants:**

| invariant | how it holds |
|---|---|
| default path bit-identical | `+ 0.0` is exact; asserted by comparing every existing law case before and after |
| a point maps to the same cell under a shifted grid and an equally shifted point | translation covariance; asserted directly |
| the nucleus inclusion stays where the nucleus is | the inclusion mask is expressed in cell indices, so it moves with the origin by construction, and the control measures the resulting box against `nucleus.positions` |

**Singular case:** a point exactly on the domain boundary. The existing clipping rule accepts it by
dropping out-of-range indices and renormalising. That behaviour is *not* changed here — but it is
the reason `nmii` looked fine while sitting outside the fluid, so this entry adds a control that
**names which sites are accepted by clipping rather than by containment**, so the next one cannot
hide.

## 7. Positive control, negative control

**None of these is written, and this section no longer names them as though they were.** Amended
2026-08-06 by lane `d95461e6`: the three controls below were specified here as test identifiers
while the port itself was never executed, and `test_named_controls_resolve_to_real_tests` flagged
two of them as citing nothing. The guard's own docstring states the obligation this entry was
failing: *"do not write down the name of a test you have not written"* — and it is explicit that
nothing forces an entry to name controls before it has them. So the specifications are kept in full
and the names are withdrawn until the code exists to carry them. **Nothing else in this entry was
changed, and no judgement is offered here on whether the port should proceed** — that is the
reviewer's, and §11 still records the status as `PROPOSED` with the PI as reviewer.

**Positive (specified, unwritten).** Every registered cytosol law case is built both with the
default origin and with no origin at all, and the state arrays must be `np.array_equal`. Not close:
the added term is exactly zero, so any difference at all means the threading is wrong.

**Negative (specified, unwritten).** A point beyond the domain must raise, not be silently
renormalised onto the nearest interior cell. And a whole-cell control that every transfer site is
inside the fluid **by containment rather than by clipping** — which is the control the measurement
in §1 says was missing, and which §1 predicts fails today.

## 8. Numerical and precision envelope

Host `float64`. The one exactness claim is the default path (`+ 0.0`), asserted with
`np.array_equal` rather than a tolerance. Nothing else here is a precision question: an origin is a
translation, and every law in the module is translation-invariant in the sense that it depends on
differences of positions and on cell indices, never on absolute coordinates.

## 9. Production-backend residency and transfer

Unchanged. The origin is host state used to build host index arrays; no device path is touched, and
no kernel takes a position in absolute coordinates.

## 10. Source identity, why not clean-room, discarded prose

No source repository, no path, no symbol, no commit — nothing was ported and
`/Users/sw1/ffn_cellsim` was not opened for this entry. Clean-room does not apply where there is no
room to clean: the change is three arithmetic expressions in Aleph's own module. No foreign prose
exists to discard.

## 11. Independent oracle

None, and that is stated rather than glossed. The invariants in §6 are self-consistency: they say
the grid is where it claims to be, not that it is where a cell's cytoplasm is. **The only external
check in this world remains the exterior drag against `6πµa`**, and it does not touch the interior
fluid.

## 12. Acceptance, reviewer, rollback

Status `PROPOSED`. Reviewer: the PI. Rollback is the default value: pass no `origin_um` and every
call returns to the current arithmetic exactly. The scenario's grid size is a separate line and
reverts independently.

## 13. Open for the PI

**The fixture is now self-consistent and still not anatomical.** A 24³ grid at `dx = 0.5` resolves
the cytosol at 500 nm, which is coarser than the cortical mesh it is supposed to percolate through.
Whether that matters depends on what the first experiments measure, and that is a scope question
rather than this entry's.
