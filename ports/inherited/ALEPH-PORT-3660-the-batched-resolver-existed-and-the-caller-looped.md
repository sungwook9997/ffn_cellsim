# ALEPH-PORT-3660 — the batched resolver existed and the caller looped anyway

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3660` |
| Lane | `b4fad06b` |
| Status | `AUDITED` |
| Written | `2026-08-09` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Follows | `ALEPH-PORT-3628`, which built `resolve_material_points`, and `ALEPH-PORT-3659`, which is the same defect one file over |
| Found by | re-profiling after `-3659`: `_cytosol_resolver.resolve` became **52 %** of the step |

---

## 1. Aleph API

- **New:** `CortexFilamentNetwork.endpoint_sinks(site_ids) -> tuple[MaterialPointSink, ...]` — the
  plural of `endpoint_sink`, resolving every site in one `resolve_material_points` call.
- **Changed:** `aleph/scenarios/whole_cell_couplings.py::_cytosol_resolver.resolve` uses it **when
  the owner publishes it**, and otherwise keeps the generator. Every other owner is untouched.

**Nothing else moves.** The returned objects are the same class, in the same order, carrying the
same `MaterialPoint` field values, so every downstream consumer — `_ImmersedSolidBlock`, the
velocity loop, the quadrature loop — sees what it saw.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none.** Nothing read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | **not applicable** — no provider commit is behind this |
| Source path and symbol | **not applicable** — no `ffn_sim` path and no provider symbol was read |
| Read from | this tree's own profile, and `ALEPH-PORT-3628`'s `resolve_material_points` |

## 3. Why source-derived porting beats clean-room

It does not; nothing is source-derived. **RE-DERIVED** from this repository's own measurement.

## 4. Physical or mathematical law represented

**None.** This is the same call, batched. `ALEPH-PORT-3628` §4 established that
`resolve_material_points` *"performs the same IEEE-754 double operations on the same operands as the
singular form, so the results are equal as bit patterns rather than to a tolerance"*, and that claim
is what this entry consumes rather than re-derives.

## 5. Units, domains, singular cases, invariants

- **I1 — bit-identity.** `endpoint_sinks(ids)` and `tuple(endpoint_sink(i) for i in ids)` produce
  material points whose `node_a`, `node_b`, `weight` and `position` are equal as bit patterns.
- **I2 — refusal order is preserved.** A loop stops at the first bad site; so must the batch, with
  the same exception type and the same site named. `-3628` already specifies and tests this.
- **I3 — the epoch is the same.** Each sink records `topology_epoch` at issue; batching must not
  make one sink older than another, and it does not, because the batch is issued at one instant.
- **I4 — owners that do not publish `endpoint_sinks` are untouched**, including `microtubule`, whose
  sinks carry `tangent_unit`. **A generic batched path that dropped that attribute would silently
  give the anisotropic drag an arbitrary axis** — so the caller checks for the method rather than
  assuming every owner has one.

## 6. Source evidence class and known retractions

Not a literature claim. Evidence class: **profile and equality test on this tree.**

## 7. Independent oracle or derivation

**The singular loop is the oracle**, and it stays in the tree as the fallback path, so the control
is a comparison against live code rather than against a frozen copy.

**Measured before the code, `cProfile`, 10 steps of the sourced whole cell after `-3659`:**

| | cumulative / 10 steps | share |
|---|---:|---|
| whole relaxation | 11.76 s | |
| **`_cytosol_resolver.resolve`** | **6.15 s** | **52 %** |
| `endpoint_sink` | 3.01 s | 311,716 calls |
| … `resolve_material_point` inside it | 2.68 s (1.17 s `tottime`) | 311,716 calls |

**31,172 Python-level site resolutions per step**, each building a `MaterialPoint` and a
`MaterialPointSink` and running its own `np.searchsorted` on a two-element table.

**Prediction, and it is deliberately more modest than `-3659`'s was.** `-3659` predicted 1.7–1.9×
from a 51 % profile share and measured **1.41×**, because `cProfile` inflates Python-call-heavy code
relative to array code. Applying that correction to an identical situation: **predicted 1.25–1.45×
on the step**, and the objects still have to be constructed, so the whole 52 % is not recoverable.
**A gain below 1.10× refutes the diagnosis.**

## 8. Positive control

`tests/vertical/test_endpoint_sinks_batch.py`:

- `test_the_batch_matches_the_loop_bit_for_bit` — I1, `np.array_equal` on positions and weights and
  `==` on the integer node indices, over every site of a real cortex.
- `test_the_order_is_the_order_of_the_ids` — a batch that returns the right set in the wrong order
  would pass I1's aggregate comparison if it were written carelessly; this pins the correspondence.
- `test_the_sinks_scatter_the_same_force` — the claim that reaches the physics: the same forces
  scattered through both sets of sinks give a bit-identical force array.

## 9. Deliberately failing negative control

- `test_an_unknown_site_is_refused_the_same_way` — I2, the batch raises what the loop raises.
- `test_the_batch_and_the_loop_can_disagree` — a **vacuity control**: move the nodes between the two
  resolutions and assert the positions then differ, so the equality above is not comparing a cached
  object with itself.

## 10. Numerical and precision envelope

**No tolerance.** `np.array_equal` throughout, per `-3628` §4.

## 11. Production-backend residency and transfer

**None.** Host-side resolution of indices and weights; the scatter path, its dtype and its ordering
are unchanged.

## 12. Comments and docstrings to discard

Nothing copied. The generator in `_cytosol_resolver` keeps its comment about `transfer_endpoints`
having already been batched for the same reason — that comment is the record that this call site was
optimised once and this half of it was missed.

## 13. Acceptance

`AUDITED`. §8 and §9 pass — `tests/vertical/test_endpoint_sinks_batch.py`, **6 passed**.

**Measured, and it landed under §7's range. Recorded as measured.**

| | s/step at 10 | s/step at 25 | step speedup | cumulative from before `-3659` |
|---|---:|---:|---:|---:|
| before `-3659` | 1.3291 | 1.3205 | — | — |
| after `-3659` | 0.9280 | 0.9351 | 1.41× | 1.41× |
| **after `-3660`** | **0.8037** | **0.8190** | **1.16×** | **1.65×** |

**Predicted 1.25–1.45×, measured 1.16×, refutation line 1.10×.** It clears the line, so the
diagnosis stands: the isolated resolution really is **2.27× faster** — 31,168 sites in 0.0911 s
against 0.2071 s, measured directly and bit-identical — but that call is a smaller share of the step
than `cProfile` showed. **Twice now the profiler's share over-predicted the gain, in the same
direction and by a similar factor.** The correction applied here (`-3659`'s 1.7–1.9 → 1.41) was
itself not enough; a third estimate on this engine should discount a Python-heavy profile share by
about **2×**, not 1.3×.

**The physics did not move.** Engine clock `1.045848e-07`, `1.950294e-06`, `9.823787e-06` s at 1, 10
and 25 steps — bit-identical to the two runs before it.

## 14. Honest limits

- **Still lever L1, still bounded.** `ENGINE_PERFORMANCE_TARGET.md` §3 measures that ≈20× is needed
  and that the method, not the wall time per step, is where most of it must come from.
- **Only one owner gets the batched path.** `intermediate_filament`, `microtubule` and the
  protrusions keep the loop. They have 12–64 sites against the cortex's 31,165, so the loop costs
  nothing there — but that means this entry's gain **shrinks as the other compartments get their
  sourced element counts** (defects D2/D3), and the loop will have to be batched generically then.
- **The object construction remains.** One `MaterialPoint` and one `MaterialPointSink` per site are
  still built in Python. Removing those needs the connector interface to accept a batch handle,
  which is a contract change and not this entry.
