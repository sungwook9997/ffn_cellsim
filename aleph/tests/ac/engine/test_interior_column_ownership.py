"""Card-5: membrane / cortex / nucleus must ADDRESS disjoint blocks, not the whole global array.

``STATE.md`` (c) 4 — *"Arrays are still shared"* — is the reason no ``CONNECTED`` rung and no "first
composed native cell" may be quoted.  In ``scripts/ac_gate_b_interior_column_native.py`` the membrane,
cortex and nucleus owners are handed the SAME ``cell.pos_d``/``cell.f_d``, so each of the three claims
every node in the cell, and the charter's line — *co-location in a shared array is NEVER a connection* —
is violated by construction.

This module is the host-side half of the definition of done.  It fixes two things the CUDA lane cannot
be trusted to expose on its own:

* **the check itself.**  :func:`assert_component_state_disjoint` compared device POINTERS for equality.
  A Warp slice view (``arr[a:b]``) is a real array with an OFFSET pointer, so two views that overlap —
  the exact defect a careless split introduces — passed the guard.  ``test_overlapping_views_*`` is the
  positive control: it fails against the pointer-identity version of the check.
* **the block arithmetic.**  An offset into the wrong array still indexes something, silently.
  :func:`assert_disjoint_blocks` is where a wrong ``node_off`` becomes an exception instead of a number.

What this file does NOT establish: that the native driver's owners are disjoint at full population.
That is a CUDA measurement with a run record beside its figures, and it is the only thing that can move
(c) 4 — by the PI, not by this test.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pytest

from aleph.engine.cortex_state import assert_component_state_disjoint
from aleph.engine.interior_column_slice import assert_disjoint_blocks

# ── CUDA metadata doubles.  A Warp vec3d array is rank-1 with a 24-byte stride; these carry exactly the
# metadata `_byte_span` reads (ptr / shape / strides / device) and evaluate no physics on the host. ────
_VEC3D_STRIDE = 24


@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:0"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _View:
    """A rank-1 device array or a slice view of one: the metadata, and nothing else."""

    ptr: int
    n: int
    device: _FakeDevice = field(default_factory=_FakeDevice)

    @property
    def shape(self) -> tuple[int, ...]:
        return (self.n,)

    @property
    def strides(self) -> tuple[int, ...]:
        return (_VEC3D_STRIDE,)


@dataclass(frozen=True, slots=True)
class _Owner:
    name: str
    position_d: object
    force_d: object


def _slice(base_ptr: int, node_off: int, n_nodes: int) -> _View:
    """The view ``base[node_off : node_off + n_nodes]`` would have, with Warp's pointer arithmetic."""
    return _View(base_ptr + node_off * _VEC3D_STRIDE, n_nodes)


# Two buffers, laid out as the built cell lays them out: actin [0, n_actin), then nucleus, then membrane.
_POS, _FORCE = 0x1000_0000, 0x2000_0000
_N_ACTIN, _NUC_OFF, _NUC_N, _MEM_OFF, _MEM_N = 494_802, 494_802, 642, 495_444, 40_962
_N_TOTAL = _MEM_OFF + _MEM_N


def _blocks() -> dict[str, tuple[int, int]]:
    return {
        "cortex": (0, _N_ACTIN),
        "nucleus": (_NUC_OFF, _NUC_N),
        "membrane": (_MEM_OFF, _MEM_N),
    }


def _owners_from(blocks: dict[str, tuple[int, int]]) -> list[_Owner]:
    return [
        _Owner(name, _slice(_POS, off, n), _slice(_FORCE, off, n))
        for name, (off, n) in blocks.items()
    ]


# ── the arrangement on disk today ────────────────────────────────────────────────────────────────────
def test_the_shared_global_array_arrangement_is_refused() -> None:
    """Three owners each handed the whole ``cell.pos_d``/``cell.f_d`` — what the native driver builds.

    Caught by pointer identity alone, so this passed before Card-5 too.  It is here because it is the
    arrangement the change is against, and a guard nobody exercises against the real defect is how the
    membrane/cortex pairwise check in ``SurfaceBody`` came to have never executed outside tests.
    """
    whole_pos, whole_force = _View(_POS, _N_TOTAL), _View(_FORCE, _N_TOTAL)
    shared = [_Owner(name, whole_pos, whole_force) for name in ("membrane", "cortex", "nucleus")]
    with pytest.raises(ValueError, match="NEVER a connection"):
        assert_component_state_disjoint(shared)


# ── the positive control: this is what the pointer-identity check could not see ───────────────────────
def test_overlapping_views_of_one_buffer_are_refused() -> None:
    """Two views that OVERLAP are one shared allocation wearing two pointers.

    The failure mode a split introduces: a block length taken from the wrong compartment, so the
    membrane's view runs into the nucleus's.  Both views have distinct ``ptr`` values, so a check that
    compares pointers for equality accepts this — which is why the guard now spans bytes.
    """
    overlapping = dict(_blocks())
    overlapping["nucleus"] = (_NUC_OFF, _NUC_N + 1)  # one node too long: it reaches into the membrane
    overlapping["membrane"] = (_NUC_OFF + _NUC_N, _MEM_N)
    with pytest.raises(ValueError, match="overlaps"):
        assert_component_state_disjoint(_owners_from(overlapping))


def test_disjoint_views_of_one_buffer_are_accepted() -> None:
    """Adjacent, non-overlapping blocks of one buffer are exclusive ADDRESSING and must pass.

    They are not private allocations, and nothing here claims they are — that distinction is the whole
    content of (c) 4 and it is carried in the docstrings of the builders, not asserted away here.
    """
    assert_component_state_disjoint(_owners_from(_blocks()))


def test_a_components_force_may_not_be_its_own_position() -> None:
    """Position and force must come from DIFFERENT buffers, or the force is a displacement.

    Caught by the same-allocation branch rather than the dedicated position/force one, because a component
    handed one array twice hits the composition-wide check first.  Asserted on the behaviour, not on which
    branch produces it.
    """
    same = _slice(_POS, 0, _N_ACTIN)
    with pytest.raises(ValueError, match="NEVER a connection"):
        assert_component_state_disjoint([_Owner("cortex", same, same)])


# ── the block arithmetic, where a wrong offset silently indexes something ─────────────────────────────
def test_blocks_that_leave_the_global_array_are_refused() -> None:
    """A block running past ``n_total`` reads another allocation entirely."""
    with pytest.raises(ValueError, match="outside the global array"):
        assert_disjoint_blocks({"membrane": (_MEM_OFF, _MEM_N + 1)}, n_total=_N_TOTAL)
    with pytest.raises(ValueError, match="outside the global array"):
        assert_disjoint_blocks({"membrane": (-1, _MEM_N)}, n_total=_N_TOTAL)


def test_overlapping_blocks_are_refused_and_name_both_components() -> None:
    """The error must name both sides: which two compartments were given the same nodes."""
    with pytest.raises(ValueError, match="nucleus.*membrane|membrane.*nucleus"):
        assert_disjoint_blocks(
            {"nucleus": (_NUC_OFF, _NUC_N + 8), "membrane": (_NUC_OFF + _NUC_N, _MEM_N)},
            n_total=_N_TOTAL,
        )


def test_an_empty_block_is_refused() -> None:
    """A zero-length compartment is a build defect, not an empty component."""
    with pytest.raises(ValueError, match="at least one node"):
        assert_disjoint_blocks({"membrane": (_MEM_OFF, 0)}, n_total=_N_TOTAL)


def test_the_built_cell_layout_is_disjoint() -> None:
    """The layout the driver derives from the built cell: actin, then nucleus, then membrane."""
    assert_disjoint_blocks(_blocks(), n_total=_N_TOTAL)


def test_native_driver_cannot_force_accept_a_nonconverged_inner_solve() -> None:
    """C-2 guard: the transaction predicate must come from the inner solve, never ``ones``.

    This is intentionally a source-level gate.  The dev host has no CUDA, while the failure it protects
    against occurred only in the native driver: three identical non-converged candidates advanced the clock
    because the script supplied ``wp.ones`` to :class:`CellTransaction`.  Importing the script cannot exercise
    that path without a full native allocation, but its predicate wiring is completely visible in source.
    """
    driver = Path(__file__).resolve().parents[3] / "scripts" / "ac_gate_b_interior_column_native.py"
    source = driver.read_text(encoding="utf-8")
    assert "accepted_ones" not in source
    assert "wp.copy(accepted_d, last_inner_report.converged_d)" in source
    assert "C-2 BLOCKED" in source


def test_native_driver_compares_existing_solvers_under_the_same_device_gate() -> None:
    """Numerical-operator selection must not create a second acceptance definition."""
    driver = Path(__file__).resolve().parents[3] / "scripts" / "ac_gate_b_interior_column_native.py"
    source = driver.read_text(encoding="utf-8")
    assert '"--inner-solver"' in source
    assert "inner_solver=args.inner_solver" in source
    assert "implicit_line_search_steps=args.line_search_steps" in source
    assert "line_search_objective=args.line_search_objective" in source
    assert "C-2 acceptance always retains max projected force" in source
    assert "max_inner_retries=args.max_inner_retries" in source
    assert "implicit_coarse_modes=args.coarse_modes" in source
    assert "implicit_coarse_iterations=args.coarse_iterations" in source
    assert "wp.copy(accepted_d, last_inner_report.converged_d)" in source


def test_c2_native_driver_defaults_to_the_grid_converged_membrane() -> None:
    """C-2 must not silently regress to the measured subdivision-6 membrane artifact."""
    driver = Path(__file__).resolve().parents[3] / "scripts" / "ac_gate_b_interior_column_native.py"
    source = driver.read_text(encoding="utf-8")
    assert '"--membrane-subdiv", type=int, default=8' in source
    assert "grid-truncation control" in source
    assert "cannot support acceptance" in source


def test_c2_diagnostic_is_capture_only_and_written_after_transaction_resolution() -> None:
    """Instrumentation may expose a rejected candidate but must never turn it into accepted state."""
    driver = Path(__file__).resolve().parents[3] / "scripts" / "ac_gate_b_interior_column_native.py"
    source = driver.read_text(encoding="utf-8")
    assert '"--diagnostic-out"' in source
    assert "capture_candidate=args.diagnostic_out is not None" in source
    assert source.index("slc.step(solve") < source.index("_write_c2_diagnostic(", source.index("slc.step(solve"))
    assert "C2_REJECTED_ROLLED_BACK" in source
    assert "rejected_candidate_never_advanced_biological_time" in source
    assert '"accepted_trial_counts": trial_counts' in source
    assert "accelerated_by_family_and_scale" in source
    assert "_c2_diagnostic_step_path(" in source
    assert '"physical_step_index": int(physical_step_index)' in source
