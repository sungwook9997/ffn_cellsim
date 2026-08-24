# ALEPH-PORT-3629 — construction memory and sweep readiness

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3629` |
| Lane | `codex-019fccbd` — construction memory and sweep readiness |
| Status | `PROPOSED` |
| Written | `2026-08-04` — **before production code**, per `PLAN.md` §0.2.5 and `CLAUDE.md` §3 |
| Port class | `NOT A PORT` — no source code, identifier, constant, or prose is taken from `/Users/sw1/ffn_cellsim` |
| Initial target | `aleph/vertical/medium.py`; the PI explicitly reassigned this older `PI (direct)` path in the current task after the ownership conflict was reported |
| Controls | New focused construction-memory and bit-identity controls under `tests/vertical/` |

---

## 0. Question and acceptance contract

This entry asks where the measured construction-time memory amplification comes from and permits
only substitutions that preserve the complete relaxation result bit-for-bit. The initial supplied
observation is 2,710 owner nodes with about 0.13 MB of position/force arrays but about 0.48 GB peak
process memory. Those figures are **UNVERIFIED in this entry until re-measured**.

Acceptance is numerical rather than rhetorical:

1. attribute retained and peak construction memory by allocation traceback and report MB and share;
2. report before/after `(nodes, process memory, ms/step)` for subdivision levels 2, 3, and 4;
3. compare positions after equal relaxation recipes with `numpy.array_equal`, not a tolerance;
4. attempt subdivision level 5 and report success or the exact failure;
5. run focused controls and then the repository suite, attributing foreign failures without editing
   their owners' files.

## 1. Source identity

`/Users/sw1/ffn_cellsim` is read-only. Its array-based crosslink storage may be inspected only as an
architectural comparison. This change will be derived from Aleph's own allocation traces and public
behaviour. Lines taken: **0**.

## 2. Authority and concurrent ownership

Lane W2 `46143f30` owns all of `aleph/scenarios/**`; those files are driven read-only and the
scenario-side call change is reported to W2 rather than written here. S5
`1afc4d25` owns the material-point batch work in `cortex_filaments.py`, `ecm_network.py`, and
`sf_arc.py`; those paths are also read-only here. S2/S3 runtime and cytosol paths, every guard and
firewall test, `aleph/represent/**`, and `aleph/learn/**` are excluded.

## 3. Physical invariant

No physical law, parameter, evaluation order, or floating-point operation in a law may change.
Construction representation may change only when its public results and the complete relaxed owner
positions remain bit-identical. A memory improvement that changes one position bit is rejected.

## 4. Measurement plan written before implementation

- run the supplied fresh-process level 2/3/4 benchmark;
- use `tracemalloc` snapshots around `build_whole_cell`, grouped by traceback and source line;
- use process RSS/peak RSS alongside Python-traced bytes, since native NumPy allocations and
  allocator high-water marks are not fully represented by `tracemalloc`;
- census large arrays and object populations, including crosslinks, mesh neighbourhoods, and any
  candidate all-pairs allocation;
- profile construction time separately from relaxation time;
- preserve a baseline process/code checkout or captured arrays sufficient for post-change exact
  comparison.

## 5. Evidence, controls, results, and rollback

**Pending measurement.** This section will be completed with allocation totals, before/after tables,
bit-identity results, tests, and an explicit rollback path before the lane closes.

---

## Sections completed by Lane W2 (`46143f30`), 2026-08-04 22:xx

**Written by a different lane from the one that did the work**, and that is stated rather than
hidden. `codex-019fccbd` closed with these mandatory fields unfilled. Every claim below is read off
the code and the tests in the same commit.

### Source repository identity, source path and symbol

**None.** Port class is `NOT A PORT`; the header already records that no source code, identifier,
constant or prose was taken from `/Users/sw1/ffn_cellsim`. No source path, symbol or commit exists to
cite, and the field is answered rather than skipped so an absent citation cannot be read as an
unrecorded one.

### Why source-derived porting beats clean-room

**It does not.** The restriction/prolongation pair is standard multigrid structure and was written
against this repository's own `SurfaceQuadrature`. Nothing was derived from a reference
implementation.

### Independent oracle or derivation

`test_a_10242_node_sphere_uses_a_constant_operator_and_improves_the_drag_error` grades against
**Stokes' law**, `6*pi*mu*a`, which is analytic and outside this repository. That makes it the second
externally-checkable number in the whole cell, and the *only* oracle in this entry that is not an
internal identity.

The other controls grade against identities rather than references, which is weaker and is worth
naming as such: `test_the_restriction_and_its_adjoint_preserve_virtual_work_force_and_moment` is
self-consistency, not correctness. A restriction that halved every force and a prolongation that
doubled it would satisfy virtual work exactly.

### Positive control

`test_no_budget_is_bit_identical_to_a_non_binding_budget`. With a budget that does not bind, the
restricted path must reproduce the unrestricted one **bit for bit** — so the machinery is proved to
be transparent when it is not needed, before any claim is made about what it does when it is.

`test_the_restriction_and_its_adjoint_preserve_virtual_work_force_and_moment` adds the property that
makes the coarsening admissible at all: restriction and prolongation are adjoint, so no work, no net
force and **no net moment** is created or destroyed by the change of representation. Moment is the
one that would otherwise fail silently — a non-adjoint pair conserves force while torquing the cell.

### Deliberately failing negative control

`test_a_budget_caps_only_the_dense_operator_not_the_membrane_interface` is the control against the
mistake this entry is most likely to make: capping the *interface* rather than the *operator* would
silently coarsen what the membrane sees, and every force would remain finite while the surface the
fluid resolved stopped being the surface the cell has.

`test_a_proposed_budget_is_transactional` drives the rejected-candidate path: a budget change inside
a step that is then rolled back must leave no trace, or the medium resolves a surface the cell never
had.

**Named as thin:** neither is a mutation of the shipped arithmetic. This entry ships no
deliberate-break flag, so nothing here measures what happens when the restriction is *wrong* rather
than merely mis-scoped. That is a gap in this entry and it is recorded rather than papered over.

### Production-backend residency and transfer

Host `float64` throughout. The restriction is dense linear algebra on host NumPy arrays and no part
of it crosses to a device. What the entry changes is the *size* of that host work: the dense mobility
is `O(M^2)` in the operator node count rather than `O(V^2)` in the membrane's, so the cost stops
tracking the mesh. No device residency is claimed and none is needed for the result.

### Comments and docstrings to discard

**None.** Nothing was copied; all prose is new and written in Aleph's vocabulary.
