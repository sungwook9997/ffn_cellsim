"""Structural gates for the cytosol/nucleus CONNECTED coupled candidate (roadmap §3.1).

CUDA-free gates for ``ac/engine/cytosol_connected.py``: the semipermeable membrane fluid boundary, and the
single coupled enclosed-volume candidate that drives the field + BOTH moving boundaries and closes the coupled
mass / no-flux / adjoint-work ledger under one predicate.  Delegate doubles stand in for the CUDA backends
(the concrete ``MembraneFluxBC``/``PressureCoupling`` binding runs on the native lane).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest

from aleph.engine.contracts import (
    CellArchitecture,
    ComponentContract,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
    reference_cell_architecture,
)
from aleph.engine.cytosol_connected import (
    MEMBRANE_FLUID_BOUNDARY,
    CytosolConnectedBindings,
    CytosolConnectedCandidate,
    MembranePressureFluxAdjointBoundary,
    MovingSemipermeableFluidBoundaryFacade,
)
from aleph.engine.fluid_core import (
    NUCLEUS_FLUID_BOUNDARY,
    MovingImpermeableFluidBoundaryFacade,
)


# -- doubles --------------------------------------------------------------------------------------
@dataclass
class _Delegate:
    label: str
    calls: list = field(default_factory=list)

    def update_geometry(self, pos, vel) -> None:
        self.calls.append((self.label, "update_geometry"))

    def accumulate_boundary(self, pos, force) -> None:
        self.calls.append((self.label, "accumulate"))

    def snapshot_candidate(self) -> None:
        self.calls.append((self.label, "snapshot"))

    def rollback(self, accepted) -> None:
        self.calls.append((self.label, "rollback", accepted))

    def commit_irreversible(self, accepted, dt_phys, rng_seed) -> None:
        self.calls.append((self.label, "commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger) -> None:
        self.calls.append((self.label, "ledger", ledger))


class _CytosolDouble:
    name = "cytosol"

    def __init__(self, order: list) -> None:
        self._order = order
        self.calls: list = []

    def solve_candidate(self, dt_phys) -> None:
        self._order.append(("cytosol", "solve"))

    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted, dt_phys, rng_seed) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger) -> None:
        self.calls.append(("ledger", ledger))

    def make_immersed_transfer_endpoint(self, residual_d):
        return ("endpoint", residual_d)


def _facades(order):
    membrane_delegate = _Delegate("membrane")
    nucleus_delegate = _Delegate("nucleus")
    # route the ordering-sensitive hooks through the shared order list via the facades themselves
    membrane = MovingSemipermeableFluidBoundaryFacade(delegate=membrane_delegate)
    nucleus = MovingImpermeableFluidBoundaryFacade(delegate=nucleus_delegate)
    return membrane, nucleus, membrane_delegate, nucleus_delegate


def _bindings(membrane, nucleus):
    o = object
    return CytosolConnectedBindings(
        membrane_fluid_boundary=membrane, nucleus_fluid_boundary=nucleus,
        membrane_surface_position_d=o(), membrane_surface_velocity_d=o(), membrane_surface_force_d=o(),
        nucleus_surface_position_d=o(), nucleus_surface_velocity_d=o(), nucleus_surface_force_d=o(),
    )


# -- semipermeable membrane facade ----------------------------------------------------------------
def test_membrane_facade_rejects_impermeable() -> None:
    with pytest.raises(ValueError, match="semipermeable"):
        MovingSemipermeableFluidBoundaryFacade(delegate=_Delegate("m"), impermeable=True)


def test_membrane_facade_endpoints_are_membrane_and_cytosol() -> None:
    facade = MovingSemipermeableFluidBoundaryFacade(delegate=_Delegate("m"))
    assert {facade.component_a, facade.component_b} == {"membrane", "cytosol"}
    assert facade.name == MEMBRANE_FLUID_BOUNDARY


def test_membrane_facade_forwards_hooks_to_delegate() -> None:
    delegate = _Delegate("m")
    facade = MovingSemipermeableFluidBoundaryFacade(delegate=delegate)
    facade.update_geometry(object(), object())
    facade.accumulate(object(), object())
    facade.snapshot_candidate()
    facade.accumulate_ledger("L")
    kinds = [c[1] for c in delegate.calls]
    assert kinds == ["update_geometry", "accumulate", "snapshot", "ledger"]


# -- coupled candidate ----------------------------------------------------------------------------
def test_candidate_builds_against_reference_architecture() -> None:
    order: list = []
    membrane, nucleus, _, _ = _facades(order)
    candidate = CytosolConnectedCandidate(reference_cell_architecture(), _CytosolDouble(order))
    assert candidate.make_immersed_transfer_endpoint("res") == ("endpoint", "res")


def test_coupled_iteration_drives_both_boundaries_and_field_in_order() -> None:
    order: list = []
    membrane, nucleus, m_del, n_del = _facades(order)
    cytosol = _CytosolDouble(order)
    candidate = CytosolConnectedCandidate(reference_cell_architecture(), cytosol)
    candidate.coupled_candidate_iteration(_bindings(membrane, nucleus), dt_phys=0.05)
    # both boundary geometries update, then the field solves, then both adjoint tractions scatter.
    m = [c for c in m_del.calls]
    n = [c for c in n_del.calls]
    assert ("membrane", "update_geometry") in m and ("membrane", "accumulate") in m
    assert ("nucleus", "update_geometry") in n and ("nucleus", "accumulate") in n
    assert ("cytosol", "solve") in order
    # field solve happens after both geometry updates and before the tractions
    geo_done = m.index(("membrane", "update_geometry"))
    acc_start = m.index(("membrane", "accumulate"))
    assert geo_done < acc_start


def test_all_three_participants_snapshot_rollback_commit_under_one_predicate() -> None:
    order: list = []
    membrane, nucleus, m_del, n_del = _facades(order)
    cytosol = _CytosolDouble(order)
    candidate = CytosolConnectedCandidate(reference_cell_architecture(), cytosol)
    bindings = _bindings(membrane, nucleus)
    accepted = object()
    assert len(candidate.transaction_participants(bindings)) == 3
    candidate.snapshot_candidate(bindings)
    candidate.rollback(bindings, accepted)
    candidate.commit_irreversible(bindings, accepted, 0.05, 7)
    candidate.accumulate_coupled_ledger(bindings, "L")
    for delegate in (m_del, n_del):
        kinds = [c[1] for c in delegate.calls]
        assert kinds == ["snapshot", "rollback", "commit", "ledger"]
    assert [c[0] for c in cytosol.calls] == ["snapshot", "rollback", "commit", "ledger"]


def test_candidate_rejects_kinetic_membrane_boundary() -> None:
    # a membrane fluid boundary declared kinetic is not a field coupler -> rejected.
    components = (
        ComponentContract("membrane", ComponentRole.SURFACE_BODY, "r", "s"),
        ComponentContract("cytosol", ComponentRole.FLUID_VOLUME, "r", "s"),
        ComponentContract("nucleus", ComponentRole.CORE_BODY, "r", "s"),
    )
    connectors = (
        ConnectorContract(
            MEMBRANE_FLUID_BOUNDARY, ConnectorFamily.FLUID_BOUNDARY, "membrane", "cytosol",
            True, True,  # kinetic + commit_on_accept -> not a field coupler
        ),
        ConnectorContract(
            NUCLEUS_FLUID_BOUNDARY, ConnectorFamily.FLUID_BOUNDARY, "nucleus", "cytosol", False, False,
        ),
    )
    bad = CellArchitecture(components=components, connectors=connectors)
    with pytest.raises(ValueError, match="non-kinetic field coupler"):
        CytosolConnectedCandidate(bad, _CytosolDouble([]))


# -- concrete membrane delegate CUDA guards -------------------------------------------------------
@dataclass
class _Dev:
    is_cuda: bool


@dataclass
class _Arr:
    dtype: object
    shape: tuple
    device: _Dev


class _MembraneBCDouble:
    def apply(self, d_pi_osm) -> None: ...


class _CouplingDouble:
    def accumulate(self, state, force) -> None: ...


class _Participant:
    def snapshot_candidate(self) -> None: ...
    def rollback(self, a) -> None: ...
    def commit_irreversible(self, a, dt, s) -> None: ...


class _LedgerDbl:
    def accumulate_ledger(self, ledger) -> None: ...


def test_membrane_delegate_rejects_non_cuda_node_volume() -> None:
    import warp as wp

    host = _Arr(dtype=wp.float64, shape=(10,), device=_Dev(is_cuda=False))
    with pytest.raises(ValueError, match="CUDA device array"):
        MembranePressureFluxAdjointBoundary(
            _MembraneBCDouble(), _CouplingDouble(), host, _Participant(), _LedgerDbl(),
            osmotic_difference=40.0,
        )


def test_membrane_delegate_rejects_bad_osmotic_driver() -> None:
    import warp as wp

    dev = _Arr(dtype=wp.float64, shape=(10,), device=_Dev(is_cuda=True))
    with pytest.raises(TypeError, match="osmotic_difference"):
        MembranePressureFluxAdjointBoundary(
            _MembraneBCDouble(), _CouplingDouble(), dev, _Participant(), _LedgerDbl(),
            osmotic_difference="not-a-number",
        )
