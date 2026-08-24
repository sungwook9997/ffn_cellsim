"""Structural gates for the ``immersed_transfer`` connector family — one runtime, seven declared edges.

CUDA-free: the pressure coupling is a recording double, so these gates check the BINDING CONTRACT (which
edges the family admits, where its endpoints come from, what it refuses) and never physics.  The
``-alpha*V*grad(p)`` scatter itself is the injected coupling's, and it runs first on the GPU.

The gate that matters most here is :func:`test_endpoints_come_from_the_contract_not_the_caller`.  The
runtime this generalises used to validate a hard-coded ``{cortex, cytosol}`` pair, which could only ever
protect the one edge it named and would have had to be copied six times.  Reading the endpoints from the
:class:`~aleph.engine.contracts.ConnectorContract` makes a caller unable to mis-state them, because it is
never asked.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest
import warp as wp

from aleph.engine.contracts import ConnectorFamily, reference_cell_architecture
from aleph.engine.immersed_transfer import (
    IMMERSED_TRANSFER_EDGES,
    ImmersedPorousTransfer,
    build_immersed_porous_transfer,
    immersed_transfer_edges,
)

CYTOSOL = "cytosol"


@dataclass(frozen=True, slots=True)
class _FakeArray:
    """Metadata-only array double.  ``dtype`` is real because ``SurfaceQuadratureState`` validates it."""

    tag: str
    dtype: object = wp.vec3d


@dataclass(slots=True)
class _RecordingPressureCoupling:
    """Stands in for the real ``PressureCoupling``; records the (state, force) it was asked to scatter."""

    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, state: object, force: object) -> None:
        self.calls.append((state, force))


@dataclass(frozen=True, slots=True)
class _SolidView:
    position_d: _FakeArray
    force_d: _FakeArray


def _transfer(name: str) -> tuple[ImmersedPorousTransfer, _RecordingPressureCoupling]:
    pc = _RecordingPressureCoupling()
    return (
        build_immersed_porous_transfer(
            name=name, pressure_coupling=pc, node_volume_d=_FakeArray("node_vol", wp.float64),
        ),
        pc,
    )


def test_the_family_has_seven_declared_edges_and_the_list_is_read_from_the_contract() -> None:
    """A new solid that declares a cytosol transfer is covered without editing the family module."""
    arch = reference_cell_architecture()
    declared = tuple(
        c.name for c in arch.connectors if c.family is ConnectorFamily.IMMERSED_TRANSFER
    )
    assert immersed_transfer_edges() == declared
    assert IMMERSED_TRANSFER_EDGES == declared
    assert len(declared) == 7


def test_every_declared_edge_in_the_family_constructs_from_the_one_runtime() -> None:
    """Seven edges, one implementation — the point of the family."""
    for name in IMMERSED_TRANSFER_EDGES:
        transfer, _ = _transfer(name)
        assert transfer.name == name


def test_endpoints_come_from_the_contract_not_the_caller() -> None:
    """The constructor is never told which components it joins, so it cannot be told wrong."""
    arch = reference_cell_architecture()
    by_name = {c.name: c for c in arch.connectors}
    for name in IMMERSED_TRANSFER_EDGES:
        transfer, _ = _transfer(name)
        contract = by_name[name]
        assert (transfer.component_a, transfer.component_b) == (contract.component_a, contract.component_b)
        assert transfer.component_b == CYTOSOL, "every immersed transfer ends at the cytosol field"
        assert transfer.bidirectional is contract.bidirectional
        assert transfer.adjoint_transfer_required is contract.adjoint_transfer_required


def test_an_edge_from_another_family_is_refused() -> None:
    """An immersed transfer may not stand in for a motor or a LINC edge."""
    with pytest.raises(ValueError, match="not 'immersed_transfer'"):
        _transfer("nmii_cortex_motor")


def test_an_undeclared_edge_name_is_refused() -> None:
    with pytest.raises(ValueError, match="not declared in the architecture"):
        _transfer("cytosol_transfer_that_does_not_exist")


def test_a_coupling_without_accumulate_is_refused_at_construction() -> None:
    with pytest.raises(TypeError, match="PressureCoupling"):
        ImmersedPorousTransfer(
            name="sf_cytosol_transfer", pressure_coupling=object(), node_volume_d=_FakeArray("nv"),
        )


def test_transfer_scatters_onto_the_solid_owned_force_array_with_its_own_node_volumes() -> None:
    """The traction lands on the immersed component's OWN force array, never on a shared one."""
    transfer, pc = _transfer("sf_cytosol_transfer")
    solid = _SolidView(position_d=_FakeArray("sf_pos"), force_d=_FakeArray("sf_force"))
    transfer.accumulate_transfer(solid, cytosol=_FakeArray("cytosol_endpoint"))

    assert len(pc.calls) == 1
    state, force = pc.calls[0]
    assert state.node_pos.tag == "sf_pos"
    assert state.node_volume.tag == "node_vol"
    assert force.tag == "sf_force"


def test_mechanics_alias_and_transfer_hook_scatter_the_same_thing() -> None:
    transfer, pc = _transfer("mt_cytosol_transfer")
    pos, force = _FakeArray("mt_pos"), _FakeArray("mt_force")
    transfer.accumulate(pos, force)
    transfer.accumulate_transfer(_SolidView(position_d=pos, force_d=force), cytosol=None)

    assert [(s.node_pos.tag, f.tag) for s, f in pc.calls] == [("mt_pos", "mt_force")] * 2


def test_transaction_api_is_complete_and_owns_no_reversible_state() -> None:
    """All three hooks (a partial API is a ``TypeError`` at world construction) and all three are no-ops."""
    transfer, pc = _transfer("if_cytosol_transfer")
    transfer.snapshot_candidate()
    transfer.rollback(_FakeArray("accepted"))
    transfer.commit_irreversible(_FakeArray("accepted"), dt_phys=0.05, rng_seed=1)
    assert pc.calls == [], "a non-kinetic field coupler advances and restores nothing of its own"


def test_fidelity_markers_say_per_node_not_lumped() -> None:
    transfer, _ = _transfer("filopodium_cytosol_transfer")
    assert transfer.per_node_pressure_traction is True
    assert transfer.aggregate_or_lumped is False
