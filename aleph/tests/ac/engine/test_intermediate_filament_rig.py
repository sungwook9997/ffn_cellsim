"""Structural gates for the nonlinear intermediate-filament rig seam."""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import pytest
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ConnectorContract,
    ConnectorFamily,
    reference_cell_architecture,
)
from aleph.engine.intermediate_filament_rig import (
    IF_CYTOSOL_TRANSFER,
    IF_KERATIN_WALL_EA_PN,
    IF_LINEAR_TANGENT_REGIME,
    IF_LP_KERATIN_UM,
    IF_LP_VIMENTIN_UM,
    IF_NUCLEUS_LINC,
    IF_SF_PLECTIN,
    IF_VIMENTIN_WALL_EA_PN,
    IF_WLC_STRAIN_STIFFENING_REGIME,
    KBT_37C_PN_UM,
    LEGACY_REFERENCE_STATUS,
    BackboneTurnoverKinetics,
    CytosolFieldEndpoint,
    IntermediateFilamentTurnoverCard,
    IntermediateFilamentRig,
    IntermediateFilamentRigStateOwner,
    IntermediateFilamentRigView,
    IntermediateFilamentStepBindings,
    LegacyIntermediateFilamentCageAdapter,
    LinearBackboneCableAdapter,
    NonlinearCableCard,
    NonlinearWlcCableAdapter,
    StructuralEndpoint,
    backbone_turnover_reference,
    keratin_nonlinear_cable_card,
    keratin_vimentin_backbone_stiffness,
    keratin_vimentin_persistence_length,
    vimentin_nonlinear_cable_card,
)
from aleph.laws.network_warp import link_spring_kernel, wlc_spring_kernel, xl_turnover_kernel


@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:0"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    """CUDA metadata double; no test evaluates authoritative physics on the host."""

    ptr: int
    shape: tuple[int, ...] = (8,)
    dtype: object = wp.vec3d
    device: _FakeDevice = _FakeDevice()


@dataclass(slots=True)
class _MechanicsSpy:
    n_nodes: int = 8
    component_local: bool = True
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        self.calls.append((pos, force))


@dataclass(slots=True)
class _TransactionSpy:
    covered: tuple[object, ...]
    component_name: str = "intermediate_filament"
    calls: list[tuple[object, ...]] = field(default_factory=list)

    def owned_arrays(self) -> tuple[wp.array, ...]:
        return self.covered  # type: ignore[return-value]

    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))


@dataclass(slots=True)
class _LedgerSpy:
    calls: list[object] = field(default_factory=list)

    def accumulate_ledger(self, ledger: object) -> None:
        self.calls.append(ledger)


@dataclass(slots=True)
class _GraphConnectorSpy:
    name: str
    component_a: str
    component_b: str
    mechanics_calls: list[tuple[IntermediateFilamentRigView, StructuralEndpoint]] = field(
        default_factory=list
    )
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate_rig(
        self,
        rig: IntermediateFilamentRigView,
        endpoint: StructuralEndpoint,
    ) -> None:
        self.mechanics_calls.append((rig, endpoint))

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


@dataclass(slots=True)
class _FluidTransferSpy:
    name: str = IF_CYTOSOL_TRANSFER
    component_a: str = "intermediate_filament"
    component_b: str = "cytosol"
    mechanics_calls: list[tuple[IntermediateFilamentRigView, CytosolFieldEndpoint]] = field(
        default_factory=list
    )
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate_transfer(
        self,
        rig: IntermediateFilamentRigView,
        cytosol: CytosolFieldEndpoint,
    ) -> None:
        self.mechanics_calls.append((rig, cytosol))

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


def _state_arrays() -> dict[str, _FakeArray]:
    return {
        "position_d": _FakeArray(100),
        "force_d": _FakeArray(101),
        "fiber_offset_d": _FakeArray(102, shape=(3,), dtype=wp.int32),
        "retained_count_d": _FakeArray(103, shape=(2,), dtype=wp.int32),
        "constitutive_state_d": _FakeArray(104, shape=(7,), dtype=wp.float64),
        "turnover_state_d": _FakeArray(105, shape=(2,), dtype=wp.int32),
        "material_epoch_d": _FakeArray(106, shape=(2,), dtype=wp.int32),
        "internal_crosslink_pair_d": _FakeArray(107, shape=(4,), dtype=wp.vec2i),
        "internal_crosslink_state_d": _FakeArray(108, shape=(4,), dtype=wp.int32),
        "mapping_epoch_d": _FakeArray(109, shape=(1,), dtype=wp.int32),
        "topology_epoch_d": _FakeArray(110, shape=(1,), dtype=wp.int32),
    }


def _state() -> tuple[
    IntermediateFilamentRigStateOwner,
    _MechanicsSpy,
    _TransactionSpy,
    _LedgerSpy,
]:
    arrays = _state_arrays()
    mechanics = _MechanicsSpy()
    transaction = _TransactionSpy(
        covered=(
            arrays["position_d"],
            arrays["retained_count_d"],
            arrays["constitutive_state_d"],
            arrays["turnover_state_d"],
            arrays["material_epoch_d"],
            arrays["internal_crosslink_pair_d"],
            arrays["internal_crosslink_state_d"],
            arrays["mapping_epoch_d"],
            arrays["topology_epoch_d"],
        )
    )
    ledger = _LedgerSpy()
    state = IntermediateFilamentRigStateOwner(
        n_filaments=2,
        mechanics=mechanics,
        transaction=transaction,
        ledger=ledger,
        **arrays,  # type: ignore[arg-type]
    )
    return state, mechanics, transaction, ledger


def _bindings() -> tuple[
    IntermediateFilamentStepBindings,
    _GraphConnectorSpy,
    _GraphConnectorSpy,
    _FluidTransferSpy,
]:
    nucleus = StructuralEndpoint("nucleus", _FakeArray(200), _FakeArray(201))
    sf_arc = StructuralEndpoint("sf_arc", _FakeArray(300), _FakeArray(301))
    cytosol = CytosolFieldEndpoint(
        "cytosol",
        _FakeArray(400, shape=(4, 5, 6), dtype=wp.float64),
        _FakeArray(401, shape=(4, 5, 6), dtype=wp.float64),
    )
    linc = _GraphConnectorSpy(IF_NUCLEUS_LINC, "intermediate_filament", "nucleus")
    plectin = _GraphConnectorSpy(IF_SF_PLECTIN, "intermediate_filament", "sf_arc")
    transfer = _FluidTransferSpy()
    return (
        IntermediateFilamentStepBindings(nucleus, sf_arc, cytosol, linc, plectin, transfer),
        linc,
        plectin,
        transfer,
    )


def _architecture_without(connector_name: str) -> CellArchitecture:
    base = reference_cell_architecture()
    return CellArchitecture(
        base.components,
        tuple(connector for connector in base.connectors if connector.name != connector_name),
    )


def _architecture_replacing(replacement: ConnectorContract) -> CellArchitecture:
    base = reference_cell_architecture()
    return CellArchitecture(
        base.components,
        tuple(
            replacement if connector.name == replacement.name else connector
            for connector in base.connectors
        ),
    )


def test_reference_architecture_registers_if_linc_plectin_and_fluid_transfer() -> None:
    architecture = reference_cell_architecture()
    state, *_ = _state()
    rig = IntermediateFilamentRig(architecture, state)
    contracts = {connector.name: connector for connector in architecture.connectors}

    assert rig.architecture is architecture
    assert contracts[IF_NUCLEUS_LINC].family is ConnectorFamily.LINC
    assert contracts[IF_SF_PLECTIN].family is ConnectorFamily.PLECTIN
    assert contracts[IF_NUCLEUS_LINC].kinetics and contracts[IF_NUCLEUS_LINC].commit_on_accept
    assert contracts[IF_SF_PLECTIN].kinetics and contracts[IF_SF_PLECTIN].commit_on_accept
    assert contracts[IF_CYTOSOL_TRANSFER].family is ConnectorFamily.IMMERSED_TRANSFER
    assert not contracts[IF_CYTOSOL_TRANSFER].kinetics
    assert not contracts[IF_CYTOSOL_TRANSFER].commit_on_accept


@pytest.mark.parametrize("missing", [IF_NUCLEUS_LINC, IF_SF_PLECTIN, IF_CYTOSOL_TRANSFER])
def test_rig_refuses_to_hide_a_required_graph_edge(missing: str) -> None:
    state, *_ = _state()
    with pytest.raises(ValueError, match=missing):
        IntermediateFilamentRig(_architecture_without(missing), state)


def test_architecture_enforces_connector_family_and_kinetic_semantics() -> None:
    state, *_ = _state()
    wrong_linc = ConnectorContract(
        IF_NUCLEUS_LINC,
        ConnectorFamily.PLECTIN,
        "intermediate_filament",
        "nucleus",
        True,
        True,
    )
    with pytest.raises(ValueError, match="ConnectorFamily.LINC"):
        IntermediateFilamentRig(_architecture_replacing(wrong_linc), state)

    kinetic_fluid = ConnectorContract(
        IF_CYTOSOL_TRANSFER,
        ConnectorFamily.IMMERSED_TRANSFER,
        "intermediate_filament",
        "cytosol",
        True,
        True,
    )
    with pytest.raises(ValueError, match="must be non-kinetic"):
        IntermediateFilamentRig(_architecture_replacing(kinetic_fluid), state)


def test_state_exposes_only_material_coordinate_geometry_to_connectors() -> None:
    state, *_ = _state()
    view = state.geometry()

    assert view.position_d is state.position_d
    assert view.retained_count_d is state.retained_count_d
    assert view.material_epoch_d is state.material_epoch_d
    assert view.mapping_epoch_d is state.mapping_epoch_d
    assert view.topology_epoch_d is state.topology_epoch_d
    assert view.n_filaments == 2
    assert not hasattr(view, "constitutive_state_d")
    assert not hasattr(view, "internal_crosslink_state_d")


def test_state_rejects_cpu_cross_device_and_aliased_arrays() -> None:
    state, *_ = _state()
    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    with pytest.raises(ValueError, match="CUDA device array"):
        replace(state, position_d=_FakeArray(999, device=cpu))

    other_gpu = _FakeDevice(alias="cuda:1")
    with pytest.raises(ValueError, match="share one CUDA device"):
        replace(state, force_d=_FakeArray(998, device=other_gpu))

    with pytest.raises(ValueError, match="must not alias"):
        replace(
            state,
            internal_crosslink_state_d=_FakeArray(
                state.internal_crosslink_pair_d.ptr,
                shape=(4,),
                dtype=wp.int32,
            ),
        )


def test_transaction_covers_nonlinear_reduction_turnover_crosslinks_and_epochs() -> None:
    state, _, transaction, _ = _state()
    missing_topology_epoch = _TransactionSpy(covered=transaction.covered[:-1])
    with pytest.raises(ValueError, match="transaction must cover"):
        replace(state, transaction=missing_topology_epoch)

    wrong_owner = _TransactionSpy(covered=transaction.covered, component_name="global_monolith")
    with pytest.raises(ValueError, match="owned by intermediate_filament"):
        replace(state, transaction=wrong_owner)


def test_facade_has_no_private_linc_plectin_or_fluid_state() -> None:
    state, *_ = _state()
    rig = IntermediateFilamentRig(reference_cell_architecture(), state)
    bindings, linc, plectin, transfer = _bindings()

    assert set(type(rig).__slots__) == {"architecture", "state"}
    assert not hasattr(rig, "linc_state")
    assert bindings.nucleus_linc is linc
    assert bindings.sf_plectin is plectin
    assert bindings.cytosol_transfer is transfer


def test_mechanics_wires_public_view_to_external_nucleus_sf_and_cytosol() -> None:
    state, mechanics, *_ = _state()
    rig = IntermediateFilamentRig(reference_cell_architecture(), state)
    bindings, linc, plectin, transfer = _bindings()
    rig.accumulate_mechanics(bindings)

    assert mechanics.calls == [(state.position_d, state.force_d)]
    linc_view, nucleus_endpoint = linc.mechanics_calls[0]
    plectin_view, sf_endpoint = plectin.mechanics_calls[0]
    fluid_view, fluid_endpoint = transfer.mechanics_calls[0]
    assert linc_view.position_d is state.position_d
    assert plectin_view.mapping_epoch_d is state.mapping_epoch_d
    assert fluid_view.topology_epoch_d is state.topology_epoch_d
    assert nucleus_endpoint is bindings.nucleus
    assert sf_endpoint is bindings.sf_arc
    assert fluid_endpoint is bindings.cytosol


def test_transaction_and_ledger_include_reversible_fluid_stencil_cache() -> None:
    state, _, transaction, state_ledger = _state()
    rig = IntermediateFilamentRig(reference_cell_architecture(), state)
    bindings, linc, plectin, transfer = _bindings()
    accepted_d = object()
    ledger = object()

    rig.snapshot_candidate(bindings)
    rig.rollback(bindings, accepted_d)
    rig.commit_irreversible(bindings, accepted_d, dt_phys=0.04, rng_seed=91)
    rig.accumulate_ledger(bindings, ledger)

    expected = [
        ("snapshot",),
        ("rollback", accepted_d),
        ("commit", accepted_d, 0.04, 91),
    ]
    assert transaction.calls == expected
    assert linc.transaction_calls == expected
    assert plectin.transaction_calls == expected
    assert transfer.transaction_calls == expected
    assert state_ledger.calls == [ledger]
    assert linc.ledger_calls == plectin.ledger_calls == transfer.ledger_calls == [ledger]


def test_bindings_reject_misdirected_or_private_connectors() -> None:
    bindings, linc, plectin, transfer = _bindings()
    with pytest.raises(ValueError, match=IF_NUCLEUS_LINC):
        replace(bindings, nucleus_linc=replace(linc, name="embedded_linc_spring"))
    with pytest.raises(ValueError, match="incorrect component endpoints"):
        replace(bindings, sf_plectin=replace(plectin, component_b="nucleus"))
    with pytest.raises(ValueError, match=IF_CYTOSOL_TRANSFER):
        replace(bindings, cytosol_transfer=replace(transfer, name="private_drag"))


def test_external_endpoints_validate_cuda_shape_dtype_and_storage() -> None:
    with pytest.raises(ValueError, match="must not alias"):
        StructuralEndpoint("nucleus", _FakeArray(1), _FakeArray(1))
    with pytest.raises(TypeError, match="float64"):
        CytosolFieldEndpoint(
            "cytosol",
            _FakeArray(2, shape=(2, 2, 2), dtype=wp.int32),
            _FakeArray(3, shape=(2, 2, 2), dtype=wp.float64),
        )
    with pytest.raises(ValueError, match="shape"):
        CytosolFieldEndpoint(
            "cytosol",
            _FakeArray(4, shape=(2, 2, 2), dtype=wp.float64),
            _FakeArray(5, shape=(2, 2, 3), dtype=wp.float64),
        )


@dataclass(slots=True)
class _LegacyCageSpy:
    n_if: int = 8
    n_bonds: int = 13
    n_backbone: int = 7
    n_linc: int = 2
    n_anchor: int = 2
    n_xl: int = 2
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        self.calls.append((pos, force))


def test_legacy_cage_adapter_is_reference_only_and_truthfully_labelled() -> None:
    legacy = _LegacyCageSpy()
    adapter = LegacyIntermediateFilamentCageAdapter(legacy)
    combined_pos = _FakeArray(1)
    combined_force = _FakeArray(2)
    adapter.accumulate(combined_pos, combined_force)

    assert not adapter.component_local
    assert adapter.status == LEGACY_REFERENCE_STATUS
    assert legacy.calls == [(combined_pos, combined_force)]

    with pytest.raises(ValueError, match="non-empty cage"):
        LegacyIntermediateFilamentCageAdapter(_LegacyCageSpy(n_if=0))
    with pytest.raises(ValueError, match="exactly partition"):
        LegacyIntermediateFilamentCageAdapter(_LegacyCageSpy(n_backbone=6))


def test_component_state_rejects_monolithic_legacy_adapter_and_count_mismatch() -> None:
    state, *_ = _state()
    adapter = LegacyIntermediateFilamentCageAdapter(_LegacyCageSpy())
    with pytest.raises(ValueError, match="monolithic legacy IF mechanics"):
        replace(state, mechanics=adapter)

    with pytest.raises(ValueError, match="expects 7 nodes"):
        replace(state, mechanics=_MechanicsSpy(n_nodes=7))


# ── Component-local backbone kernel binding (SEAMED → KERNEL_BOUND) ────────────────────────────────


def _backbone_segment(ptr: int, n_seg: int = 5) -> _FakeArray:
    return _FakeArray(ptr, shape=(n_seg, 2), dtype=wp.int32)


def _backbone_scalar(ptr: int, n_seg: int = 5) -> _FakeArray:
    return _FakeArray(ptr, shape=(n_seg,), dtype=wp.float64)


def _linear_backbone_adapter(n_nodes: int = 8) -> LinearBackboneCableAdapter:
    return LinearBackboneCableAdapter(
        n_nodes=n_nodes,
        segment_d=_backbone_segment(700),
        stiffness_d=_backbone_scalar(701),
        rest_d=_backbone_scalar(702),
        device="cuda:0",
        cell_type="keratin",
    )


def _sourced_card() -> NonlinearCableCard:
    return NonlinearCableCard(
        persistence_length_um=IF_LP_VIMENTIN_UM,
        axial_stiffness_pn=470.0,          # caller-supplied EA (PI GAP) — a plausible positive placeholder
        thermal_energy_pn_um=KBT_37C_PN_UM,
        crossover_ratio=0.9,               # caller-supplied x_max (PI GAP) in (0, 1)
    )


def _wlc_backbone_adapter(n_nodes: int = 8) -> NonlinearWlcCableAdapter:
    return NonlinearWlcCableAdapter(
        n_nodes=n_nodes,
        segment_d=_backbone_segment(710),
        contour_length_d=_backbone_scalar(711),
        card=_sourced_card(),
        device="cuda:0",
        cell_type="vimentin",
    )


def test_linear_backbone_adapter_binds_reused_link_spring_kernel() -> None:
    adapter = _linear_backbone_adapter()
    assert adapter.bound_kernel is link_spring_kernel
    assert adapter.component_local is True
    assert adapter.regime == IF_LINEAR_TANGENT_REGIME
    assert adapter.n_nodes == 8
    # a component-local backbone carries no graph/crosslink connector-bond family counts
    for attr in ("n_linc", "n_anchor", "n_xl", "n_bonds"):
        assert not hasattr(adapter, attr)


def test_wlc_backbone_adapter_binds_reused_wlc_kernel() -> None:
    adapter = _wlc_backbone_adapter()
    assert adapter.bound_kernel is wlc_spring_kernel
    assert adapter.component_local is True
    assert adapter.regime == IF_WLC_STRAIN_STIFFENING_REGIME


def test_state_owner_installs_component_local_backbone_kernel_backend() -> None:
    # the monolithic cage is rejected; a real-kernel-bound component-local backbone backend is accepted.
    state, *_ = _state()
    installed = replace(state, mechanics=_linear_backbone_adapter(n_nodes=8))
    assert installed.mechanics.bound_kernel is link_spring_kernel
    assert installed.mechanics.regime == IF_LINEAR_TANGENT_REGIME

    installed_wlc = replace(state, mechanics=_wlc_backbone_adapter(n_nodes=8))
    assert installed_wlc.mechanics.bound_kernel is wlc_spring_kernel

    with pytest.raises(ValueError, match="expects 5 nodes|expects 8 nodes"):
        replace(state, mechanics=_linear_backbone_adapter(n_nodes=5))


def test_backbone_adapter_rejects_global_offset_cpu_and_aliased_arrays() -> None:
    with pytest.raises(ValueError, match="node_offset == 0"):
        replace(_linear_backbone_adapter(), node_offset=3)

    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    with pytest.raises(ValueError, match="CUDA device array"):
        replace(_linear_backbone_adapter(), segment_d=_FakeArray(720, shape=(5, 2), dtype=wp.int32, device=cpu))

    with pytest.raises(TypeError, match="int32"):
        replace(_linear_backbone_adapter(), segment_d=_FakeArray(721, shape=(5, 2), dtype=wp.float64))

    with pytest.raises(ValueError, match="shape \\(n_segments, 2\\)"):
        replace(_linear_backbone_adapter(), segment_d=_FakeArray(722, shape=(5, 3), dtype=wp.int32))

    with pytest.raises(ValueError, match="must not alias"):
        replace(_linear_backbone_adapter(), rest_d=_backbone_scalar(701))  # aliases stiffness ptr

    with pytest.raises(ValueError, match="cell_type"):
        replace(_linear_backbone_adapter(), cell_type="actin")


def test_nonlinear_card_forbids_fabricated_or_out_of_band_constants() -> None:
    # Lp band, positive EA/kBT, and 0 < x_max < 1 are all enforced; no defaults are provided.
    with pytest.raises(ValueError, match="band"):
        replace(_sourced_card(), persistence_length_um=5.0)
    with pytest.raises(ValueError, match="EA"):
        replace(_sourced_card(), axial_stiffness_pn=0.0)
    with pytest.raises(ValueError, match="kBT"):
        replace(_sourced_card(), thermal_energy_pn_um=0.0)
    with pytest.raises(ValueError, match="x_max"):
        replace(_sourced_card(), crossover_ratio=1.0)


def test_sourced_helpers_reuse_ff_resolver_and_band_checked_lp() -> None:
    k_keratin = keratin_vimentin_backbone_stiffness("keratin", n_fil=60, l_seg_um=0.5)
    k_vimentin = keratin_vimentin_backbone_stiffness("vimentin", n_fil=60, l_seg_um=0.5)
    # vimentin is the stiffer network (9 vs 6 MPa) → larger derived k_bb; both strictly positive.
    assert 0.0 < k_keratin < k_vimentin
    assert keratin_vimentin_persistence_length("keratin") == IF_LP_KERATIN_UM
    assert keratin_vimentin_persistence_length("vimentin") == IF_LP_VIMENTIN_UM
    # physiological kBT at 37 °C ≈ 4.28e-3 pN·µm (derived physical constant, not a tuned magic number).
    assert 4.2e-3 < KBT_37C_PN_UM < 4.35e-3
    with pytest.raises(ValueError, match="cell_type"):
        keratin_vimentin_backbone_stiffness("actin", n_fil=60, l_seg_um=0.5)


# ── IF turnover kinetics + sourced/NOT-FOUND WLC cards (PI-GAP slots) ─────────────────────────────


def test_turnover_card_is_an_unset_slot_with_no_default() -> None:
    card = IntermediateFilamentTurnoverCard(
        k_off0_per_s=0.1, x_beta_over_kt_per_pn=0.05, provenance="PI card"
    )
    assert card.ratified is False
    with pytest.raises(ValueError, match="finite and nonnegative"):
        IntermediateFilamentTurnoverCard(
            k_off0_per_s=-1.0, x_beta_over_kt_per_pn=0.05, provenance="x"
        )
    with pytest.raises(ValueError, match="provenance"):
        IntermediateFilamentTurnoverCard(
            k_off0_per_s=0.1, x_beta_over_kt_per_pn=0.05, provenance=""
        )


def _turnover_kinetics() -> BackboneTurnoverKinetics:
    return BackboneTurnoverKinetics(
        n_nodes=8,
        segment_d=_backbone_segment(800),
        stiffness_d=_backbone_scalar(801),
        rest_d=_backbone_scalar(802),
        card=IntermediateFilamentTurnoverCard(
            k_off0_per_s=0.1, x_beta_over_kt_per_pn=0.05, provenance="PI card"
        ),
        device="cuda:0",
    )


def test_backbone_turnover_binds_reused_xl_turnover_kernel() -> None:
    kinetics = _turnover_kinetics()
    assert kinetics.bound_kernel is xl_turnover_kernel
    # a segment table that is not (S, 2) node pairs is rejected (component-local backbone-only contract).
    with pytest.raises(ValueError, match=r"\(n_segments, 2\)"):
        BackboneTurnoverKinetics(
            n_nodes=8,
            segment_d=_FakeArray(810, shape=(5, 3), dtype=wp.int32),
            stiffness_d=_backbone_scalar(811),
            rest_d=_backbone_scalar(812),
            card=IntermediateFilamentTurnoverCard(
                k_off0_per_s=0.1, x_beta_over_kt_per_pn=0.05, provenance="PI card"
            ),
            device="cuda:0",
        )
    with pytest.raises(ValueError, match="must not alias"):
        BackboneTurnoverKinetics(
            n_nodes=8,
            segment_d=_backbone_segment(820),
            stiffness_d=_backbone_scalar(821),
            rest_d=_backbone_scalar(821),  # aliases stiffness ptr
            card=IntermediateFilamentTurnoverCard(
                k_off0_per_s=0.1, x_beta_over_kt_per_pn=0.05, provenance="PI card"
            ),
            device="cuda:0",
        )


def test_backbone_turnover_reference_creeps_rest_toward_current_length() -> None:
    # a stretched segment (L > r0) relaxes its rest length TOWARD the current length (stress relaxation).
    crept = backbone_turnover_reference(
        length_um=1.2, rest_um=1.0, stiffness_pn_per_um=50.0,
        k_off0_per_s=1.0, x_beta_over_kt_per_pn=0.1, dt_real=0.5,
    )
    assert 1.0 < crept < 1.2
    # zero off-rate → no turnover, rest length unchanged.
    frozen = backbone_turnover_reference(
        length_um=1.2, rest_um=1.0, stiffness_pn_per_um=50.0,
        k_off0_per_s=0.0, x_beta_over_kt_per_pn=0.1, dt_real=0.5,
    )
    assert frozen == pytest.approx(1.0)
    # higher force → higher Bell off-rate → faster creep toward L.
    faster = backbone_turnover_reference(
        length_um=1.4, rest_um=1.0, stiffness_pn_per_um=50.0,
        k_off0_per_s=1.0, x_beta_over_kt_per_pn=0.1, dt_real=0.5,
    )
    assert (faster - 1.0) / (1.4 - 1.0) > (crept - 1.0) / (1.2 - 1.0)
    with pytest.raises(ValueError, match="dt_real must be positive"):
        backbone_turnover_reference(
            length_um=1.2, rest_um=1.0, stiffness_pn_per_um=50.0,
            k_off0_per_s=1.0, x_beta_over_kt_per_pn=0.1, dt_real=0.0,
        )


def test_vimentin_card_fills_sourced_anchors_and_leaves_crossover_a_slot() -> None:
    card = vimentin_nonlinear_cable_card(crossover_ratio=0.9)
    assert card.persistence_length_um == IF_LP_VIMENTIN_UM
    assert card.axial_stiffness_pn == IF_VIMENTIN_WALL_EA_PN   # FOUND vimentin wall EA (Qin 2009)
    assert card.thermal_energy_pn_um == KBT_37C_PN_UM
    assert card.crossover_ratio == 0.9                          # caller-supplied modeling slot
    # crossover is validated by NonlinearCableCard: an out-of-range x/Lc is rejected, not silently accepted.
    with pytest.raises(ValueError, match="crossover_ratio"):
        vimentin_nonlinear_cable_card(crossover_ratio=1.5)


def test_keratin_card_is_blocked_as_a_pi_gap() -> None:
    # keratin (primary MCF7 IF) EA/x_max are NOT FOUND: the factory refuses to invent the wall stiffness.
    assert IF_KERATIN_WALL_EA_PN is None
    with pytest.raises(NotImplementedError, match="keratin K8/K18 axial EA/x_max NOT FOUND"):
        keratin_nonlinear_cable_card(crossover_ratio=0.9)
