"""Runtime protocols shared by independently implemented cell-engine components."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import warp as wp


def _device_storage_key(array: object) -> tuple[object, ...]:
    """Return device storage identity without copying device data."""
    ptr = getattr(array, "ptr", None)
    if ptr is None:
        return ("object", id(array))
    return ("device-ptr", str(getattr(array, "device", None)), int(ptr))


@dataclass(frozen=True, slots=True)
class CytosolFieldEndpoint:
    """Public pressure/residual view shared by immersed solid--cytosol transfers."""

    component: str
    pressure_d: wp.array
    residual_d: wp.array

    def __post_init__(self) -> None:
        if self.component != "cytosol":
            raise ValueError("fluid endpoint must be the registered cytosol component")
        for label, array in (("pressure_d", self.pressure_d), ("residual_d", self.residual_d)):
            device = getattr(array, "device", None)
            if not bool(getattr(device, "is_cuda", False)):
                raise ValueError(f"cytosol.{label} must be a Warp CUDA device array")
            if getattr(array, "dtype", None) != wp.float64:
                raise TypeError(f"cytosol.{label} must have dtype wp.float64")
            shape = getattr(array, "shape", None)
            if not isinstance(shape, tuple) or len(shape) != 3 or any(int(size) <= 0 for size in shape):
                raise ValueError(f"cytosol.{label} must be a non-empty three-dimensional array")
        if self.pressure_d.shape != self.residual_d.shape:
            raise ValueError("cytosol pressure and residual arrays must have identical shape")
        if str(self.pressure_d.device) != str(self.residual_d.device):
            raise ValueError("cytosol pressure and residual arrays must share one CUDA device")
        if _device_storage_key(self.pressure_d) == _device_storage_key(self.residual_d):
            raise ValueError("cytosol pressure and residual arrays must not alias")


@runtime_checkable
class MechanicsContributor(Protocol):
    """A component or connector that adds force to state it is authorized to address."""

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Add the current candidate's mechanical force without committing kinetics."""


@runtime_checkable
class TransactionParticipant(Protocol):
    """Device-resident state that participates in the global physical-step transaction.

    Reversible mechanics/field state must snapshot and roll back even when a contract has
    ``commit_on_accept=False``.  That flag means there is no separate irreversible kinetic
    commit; it never exempts a runtime from candidate-step restoration.
    """

    def snapshot_candidate(self) -> None:
        """Capture the accepted state before candidate mutation."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore the snapshot under the scheduler-owned rejected predicate."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit kinetics/remodelling only under the final accepted predicate."""


@runtime_checkable
class LedgerContributor(Protocol):
    """A component that exposes device ledger terms to the global acceptance gate."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Add force/work/mass/topology terms without deciding acceptance."""


@runtime_checkable
class EventProposer(Protocol):
    """A kinetic component or connector that proposes candidate events (no force, no host state).

    This is the fifth scheduler method (alongside ``snapshot_candidate`` / ``rollback`` /
    ``commit_irreversible`` / ``accumulate_ledger``).  It generalises the pre-existing
    ``propose_candidate`` (LINC / MT) into one signature.  An event raises a local *rate* only; the
    force it implies materialises later inside ``commit_irreversible`` from whatever actually binds
    under the accepted predicate.
    """

    def propose_events(
        self,
        rates: object,
        dt_phys: float,
        rng_seed: int,
        neighbors: object,
    ) -> None:
        """Write this participant's device candidate-event buffer only.

        Reads this participant's LOCAL rate fields (set by :class:`~aleph.engine.cell_state.CellState`)
        plus neighbour candidate lists; writes its own candidate-event buffer.  MUST NOT write any
        force/pos array, mutate another participant, or read authoritative host state.
        """


@runtime_checkable
class EventClock(Protocol):
    """Device-resident biological clock + reproducible RNG epoch, advanced only on an accepted step.

    ``base_seed`` is the run's fixed host seed; the per-epoch draw derives on-device from
    ``base_seed`` and the device-resident accepted-step index, so a rejected candidate re-draws
    identically.  ``advance`` is device-gated on the acceptance predicate — inner mechanical
    iterations neither advance the clock nor the epoch, and no acceptance value is read on the host.
    """

    base_seed: int

    def advance(self, accepted: wp.array, dt_phys: float) -> None:
        """Advance biological time and the RNG epoch on the device iff ``accepted`` selects it."""


@runtime_checkable
class ImmersedTransferConnector(TransactionParticipant, LedgerContributor, Protocol):
    """Graph-owned adjoint solid--cytosol transfer with reversible device stencil state."""

    name: str
    component_a: str
    component_b: str

    def accumulate_transfer(self, solid: object, cytosol: CytosolFieldEndpoint) -> None:
        """Gather field load and scatter the work-conjugate reaction to a public solid view."""
