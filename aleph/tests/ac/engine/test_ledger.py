"""Structural + CUDA gates for the typed global cell-engine ledger.

The CPU-green gates inject metadata-only array doubles (authoritative physics never executes) and prove the
ledger owns disjoint CUDA accumulators, exposes the far-field / force / work / mass / topology slots a
``LedgerContributor`` pushes into, and validates every contribution's metadata before any launch.  The
CUDA_UNIT gate (GPU lane only) proves the reduction kernels sum correctly and that the global force-balance
gate — reactions + cell-traction resultant → 0 — is genuinely CHECKED, not assumed.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pytest
import warp as wp

import aleph.engine.ledger as ledger_module
from aleph.engine.ledger import (
    GlobalCellLedger,
    assemble_balance_tolerance_ratio,
    make_global_cell_ledger,
)


def _cuda_available() -> bool:
    """Return whether a real Warp CUDA device is present (GPU lane only)."""
    try:
        wp.init()
        return bool(wp.get_device().is_cuda)
    except Exception:  # pragma: no cover - depends on hardware
        return False


@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:0"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    """Metadata-only CUDA array double; authoritative physics never executes here."""

    ptr: int
    shape: tuple[int, ...]
    dtype: object
    device: _FakeDevice = _FakeDevice()


def _ledger(*, ptr_base: int = 100, device: _FakeDevice | None = None, n_topology: int = 3) -> GlobalCellLedger:
    if device is None:
        device = _FakeDevice()
    return GlobalCellLedger(
        reaction_resultant_d=_FakeArray(ptr_base, (1,), wp.vec3d, device),
        traction_resultant_d=_FakeArray(ptr_base + 1, (1,), wp.vec3d, device),
        force_resultant_d=_FakeArray(ptr_base + 2, (1,), wp.vec3d, device),
        work_d=_FakeArray(ptr_base + 3, (1,), wp.float64, device),
        mass_d=_FakeArray(ptr_base + 4, (1,), wp.float64, device),
        topology_counts_d=_FakeArray(ptr_base + 5, (n_topology,), wp.int32, device),
        balance_residual_sq_d=_FakeArray(ptr_base + 6, (1,), wp.float64, device),
        balance_ok_d=_FakeArray(ptr_base + 7, (1,), wp.int32, device),
    )


def test_ledger_owns_disjoint_cuda_accumulators_and_slots() -> None:
    ledger = _ledger()

    assert ledger.device == "cuda:0"
    assert ledger.n_topology_slots == 3
    # every slot a LedgerContributor pushes into is present and callable.
    for slot in (
        "add_far_field_reaction",
        "add_cell_traction",
        "add_force_resultant",
        "add_work",
        "add_mass",
        "add_topology_delta",
        "assemble_balance",
        "reset",
    ):
        assert callable(getattr(ledger, slot))


def test_ledger_rejects_host_free_wrong_dtype_and_aliased_accumulators() -> None:
    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    with pytest.raises(ValueError, match="CUDA device array"):
        _ledger(device=cpu)

    base = _ledger()
    with pytest.raises(TypeError, match="dtype"):
        replace(base, work_d=replace(base.work_d, dtype=wp.float32))
    with pytest.raises(ValueError, match="single-entry resultant"):
        replace(base, reaction_resultant_d=replace(base.reaction_resultant_d, shape=(2,)))
    with pytest.raises(ValueError, match="one accepted scalar"):
        replace(base, work_d=replace(base.work_d, shape=(3,)))
    with pytest.raises(ValueError, match="one device gate flag"):
        replace(base, balance_ok_d=replace(base.balance_ok_d, shape=(2,)))
    with pytest.raises(ValueError, match="distinct storage"):
        replace(base, mass_d=replace(base.work_d))
    other = _FakeDevice(alias="cuda:1")
    with pytest.raises(ValueError, match="share one CUDA device"):
        replace(base, mass_d=replace(base.mass_d, device=other))


def test_contribution_metadata_is_validated_before_any_launch() -> None:
    ledger = _ledger()
    good_vec = _FakeArray(900, (5,), wp.vec3d)
    good_scalar = _FakeArray(901, (1,), wp.float64)

    # A wrong-dtype / wrong-rank contribution raises in Python, before wp.launch — CPU-green off-CUDA.
    with pytest.raises(TypeError, match="dtype"):
        ledger.add_far_field_reaction(_FakeArray(902, (5,), wp.float64), good_scalar)
    with pytest.raises(TypeError, match="dtype"):
        ledger.add_cell_traction(good_scalar)
    with pytest.raises(ValueError, match="topology-slot count"):
        ledger.add_topology_delta(_FakeArray(903, (2,), wp.int32))  # ledger has 3 slots
    with pytest.raises(ValueError, match="single-entry squared-tolerance"):
        ledger.assemble_balance(_FakeArray(904, (2,), wp.float64))
    # cross-device contribution is rejected before launch too.
    with pytest.raises(ValueError, match="share one CUDA device"):
        ledger.add_force_resultant(_FakeArray(905, (5,), wp.vec3d, _FakeDevice(alias="cuda:7")))
    del good_vec  # only the metadata is exercised; no launch happens on the CUDA-free host.


def test_make_global_cell_ledger_validates_topology_slot_count() -> None:
    with pytest.raises(ValueError, match="positive integer"):
        make_global_cell_ledger(n_topology_slots=0, device="cuda:0")
    with pytest.raises(ValueError, match="positive integer"):
        make_global_cell_ledger(n_topology_slots=True, device="cuda:0")


@pytest.mark.skipif(not _cuda_available(), reason="CUDA_UNIT ledger reduction/balance gate runs on the GPU lane")
def test_reaction_plus_traction_balance_is_checked_not_assumed() -> None:  # pragma: no cover - GPU lane
    import numpy as np

    device = "cuda:0"
    ledger = make_global_cell_ledger(n_topology_slots=2, device=device)

    # Far-field reaction resultant R = Σ reaction; cell→ECM traction resultant T = Σ traction.
    reaction = wp.array(
        np.array([[1.0, -2.0, 3.0], [0.5, 0.5, 0.5]], dtype=np.float64), dtype=wp.vec3d, device=device
    )  # Σ = [1.5, -1.5, 3.5]
    work = wp.array(np.array([0.25], dtype=np.float64), dtype=wp.float64, device=device)
    traction = wp.array(
        np.array([[-1.5, 1.5, -3.5]], dtype=np.float64), dtype=wp.vec3d, device=device
    )  # Σ = -R exactly → balanced
    tol_sq = wp.array(np.array([1.0e-18], dtype=np.float64), dtype=wp.float64, device=device)

    ledger.add_far_field_reaction(reaction, work)
    ledger.add_cell_traction(traction)
    ledger.assemble_balance(tol_sq)

    assert np.allclose(ledger.reaction_resultant_d.numpy()[0], [1.5, -1.5, 3.5])
    assert np.allclose(ledger.traction_resultant_d.numpy()[0], [-1.5, 1.5, -3.5])
    assert float(ledger.work_d.numpy()[0]) == pytest.approx(0.25)
    assert float(ledger.balance_residual_sq_d.numpy()[0]) == pytest.approx(0.0, abs=1e-12)
    assert int(ledger.balance_ok_d.numpy()[0]) == 1  # reactions + traction sum to zero → gate passes

    # An imbalance (traction no longer equal-and-opposite) trips the gate: the check is real, not assumed.
    ledger.assemble_balance  # noqa: B018 - keep reference for readers
    extra = wp.array(np.array([[10.0, 0.0, 0.0]], dtype=np.float64), dtype=wp.vec3d, device=device)
    ledger.add_cell_traction(extra)
    ledger.assemble_balance(tol_sq)
    assert float(ledger.balance_residual_sq_d.numpy()[0]) > 1.0
    assert int(ledger.balance_ok_d.numpy()[0]) == 0

    # reset() zeroes every accumulator on the device.
    ledger.reset()
    assert np.allclose(ledger.reaction_resultant_d.numpy()[0], [0.0, 0.0, 0.0])
    assert int(ledger.balance_ok_d.numpy()[0]) == 0


@pytest.mark.skipif(not _cuda_available(), reason="CUDA_UNIT boundary→ledger forward runs on the GPU lane")
def test_far_field_boundary_forwards_committed_reaction_into_the_concrete_ledger() -> None:  # pragma: no cover
    import numpy as np

    from aleph.engine.ecm_world import FarFieldDirichletBoundaryRuntime

    device = "cuda:0"
    committed_reaction = wp.array(
        np.array([[2.0, 0.0, -1.0], [-0.5, 1.0, 0.0]], dtype=np.float64), dtype=wp.vec3d, device=device
    )
    runtime = FarFieldDirichletBoundaryRuntime(
        pinned_index_d=wp.array(np.array([0, 1], dtype=np.int32), dtype=wp.int32, device=device),
        target_position_d=wp.zeros(2, dtype=wp.vec3d, device=device),
        reaction_d=wp.zeros(2, dtype=wp.vec3d, device=device),
        committed_reaction_d=committed_reaction,
        committed_target_d=wp.zeros(2, dtype=wp.vec3d, device=device),
        boundary_work_d=wp.array(np.array([0.75], dtype=np.float64), dtype=wp.float64, device=device),
    )
    ledger = make_global_cell_ledger(device=device)

    # isinstance(ledger, GlobalCellLedger) → the concrete typed slot, not the duck-typed fallback.
    runtime.accumulate_ledger(ledger)

    assert isinstance(ledger, GlobalCellLedger)
    assert np.allclose(ledger.reaction_resultant_d.numpy()[0], [1.5, 1.0, -1.0])
    assert float(ledger.work_d.numpy()[0]) == pytest.approx(0.75)


# ── D8: the assemble_balance tolerance is DERIVED, and takes no supplied magnitude ────────────────


def test_balance_tolerance_is_a_float64_summation_bound_not_a_chosen_number() -> None:
    """γ_n is the Higham forward bound for summing n terms — a property of float64, not of the cell."""
    u = float(np.finfo(np.float64).eps) / 2.0

    # a single term accumulates no round-off, so it demands exact cancellation
    assert assemble_balance_tolerance_ratio(1) == 0.0
    # two terms: exactly one rounding
    assert assemble_balance_tolerance_ratio(2) == pytest.approx(u / (1.0 - u), rel=1e-12)
    # closed form holds at the native population
    for n in (100, 494_802, 1_413_720):
        slack = (n - 1) * u
        assert assemble_balance_tolerance_ratio(n) == pytest.approx(slack / (1.0 - slack), rel=1e-12)

    # monotone in n: more accumulation can only loosen the round-off floor
    ratios = [assemble_balance_tolerance_ratio(n) for n in (1, 2, 100, 494_802, 1_413_720)]
    assert ratios == sorted(ratios)

    # it stays a ROUND-OFF floor, not a physics threshold: at the native population with a 2e4 pN
    # force scale the tolerance is ~1e-6 pN, five orders under the 0.21 pN gate. A run cannot pass
    # this by being sloppy — only by actually balancing.
    assert assemble_balance_tolerance_ratio(494_802) * 2.0e4 < 1.0e-5


def test_balance_tolerance_rejects_counts_that_make_the_bound_meaningless() -> None:
    with pytest.raises(ValueError, match="at least one accumulated force contribution"):
        assemble_balance_tolerance_ratio(0)
    with pytest.raises(ValueError, match="at least one accumulated force contribution"):
        assemble_balance_tolerance_ratio(-5)
    with pytest.raises(ValueError, match="vacuous"):
        assemble_balance_tolerance_ratio(10**17)


def test_derive_balance_tolerance_validates_before_any_launch() -> None:
    """The derivation slot exists and rejects bad metadata on the host, never on the device."""
    ledger = _ledger()
    assert callable(ledger.derive_balance_tolerance)

    good = _FakeArray(900, (1,), wp.float64)
    with pytest.raises(ValueError, match="single-entry"):              # wrong shape
        ledger.derive_balance_tolerance(_FakeArray(901, (2,), wp.float64), good)
    with pytest.raises(TypeError, match="dtype"):                      # wrong dtype
        ledger.derive_balance_tolerance(_FakeArray(902, (1,), wp.int32), good)
    with pytest.raises(ValueError):                                    # foreign device
        ledger.derive_balance_tolerance(
            _FakeArray(903, (1,), wp.float64, _FakeDevice(alias="cuda:1")), good
        )


def test_the_gate_takes_no_supplied_force_magnitude() -> None:
    """The whole point of D8: the only free input is a COUNT.

    A tolerance that accepted a force scale would have relocated the magic number instead of removing
    it, so the scale is read on the device from the two resultants the ledger already holds.
    """
    import inspect

    params = inspect.signature(assemble_balance_tolerance_ratio).parameters
    assert list(params) == ["n_terms"]
    src = inspect.getsource(GlobalCellLedger.derive_balance_tolerance)
    assert "self.reaction_resultant_d" in src and "self.traction_resultant_d" in src


def test_add_body_force_validates_the_channel_and_the_array_before_any_launch() -> None:
    """The two-body reading of the gate: one body per channel, and neither computed from the other."""
    ledger = _ledger()
    assert callable(ledger.add_body_force)

    with pytest.raises(ValueError, match="reaction' or 'traction'"):
        ledger.add_body_force(_FakeArray(910, (4,), wp.vec3d), side="both")
    with pytest.raises(TypeError, match="dtype"):
        ledger.add_body_force(_FakeArray(911, (4,), wp.float64), side="reaction")
    with pytest.raises(ValueError):                                    # foreign device
        ledger.add_body_force(
            _FakeArray(912, (4,), wp.vec3d, _FakeDevice(alias="cuda:1")), side="traction"
        )


def test_balance_scale_is_the_sum_of_magnitudes_not_the_magnitude_of_the_sum() -> None:
    """The Higham bound's scale is Σ|x_i|, and for a DIPOLE that is nothing like |Σx_i|.

    This is the defect that would have made the gate reject every physically correct step of a motor
    lane: a bipolar minifilament's forces very nearly cancel, so a tolerance scaled by the RESULTANT
    would sit orders of magnitude below the round-off it exists to admit.
    """
    import inspect

    source = inspect.getsource(GlobalCellLedger.add_body_force)
    assert "_reduce_vec3_magnitude_kernel" in source

    # the fallback is preserved for channels that accumulate no scale, so an existing ledger built
    # without the accumulator behaves exactly as before
    ledger = _ledger()
    assert ledger.balance_scale_d is None

    # …and the derivation kernel prefers the accumulated scale when there is one
    kernel_source = inspect.getsource(ledger_module._derive_balance_tolerance_kernel)
    assert "accumulated_magnitude" in kernel_source
    assert "wp.length(reaction_resultant[0]) + wp.length(traction_resultant[0])" in kernel_source
