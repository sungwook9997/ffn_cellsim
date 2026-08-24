"""Ownership and connector seam for the nonlinear intermediate-filament rig.

The landed IF cage is a useful Hookean reference, but one global bond array mixes backbone, LINC, cortical
anchor, and crosslink state.  This module keeps the production seam component-local:

* IF owns CUDA cable geometry, nonlinear history, reduction mapping, turnover, and internal crosslinks;
* nucleus LINC, SF plectin, and cytosol immersed transfer remain graph-owned bindings;
* every mutable authoritative array is explicitly covered by the IF transaction delegate; and
* the legacy monolithic cage can be called as a reference adapter but is rejected as a component backend.

No method here reads authoritative state back to the host or advances a private physical clock.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
)
from aleph.engine.runtime import (
    CytosolFieldEndpoint,
    ImmersedTransferConnector,
    LedgerContributor,
    MechanicsContributor,
    TransactionParticipant,
)
from aleph.components.solid.intermediate_filament import IF_E_KERATIN_PA, IF_E_VIMENTIN_PA
from aleph.laws.intermediate_filaments import LP_BAND_UM, resolve_intermediate_filaments
from aleph.laws.network_warp import link_spring_kernel, wlc_spring_kernel, xl_turnover_kernel

INTERMEDIATE_FILAMENT_COMPONENT = "intermediate_filament"
NUCLEUS_COMPONENT = "nucleus"
SF_COMPONENT = "sf_arc"
CYTOSOL_COMPONENT = "cytosol"
IF_NUCLEUS_LINC = "if_nucleus_linc"
IF_SF_PLECTIN = "if_sf_plectin"
IF_CYTOSOL_TRANSFER = "if_cytosol_transfer"
LEGACY_REFERENCE_STATUS = "MONOLITHIC_HOOKEAN_REFERENCE_NONLINEAR_TURNOVER_PENDING"

# ── Component-local backbone constitutive regimes (bound to reused ff Warp kernels) ──────────────────
# The linear tangent limit is the sourced, immediately-installable backbone law (k_bb = E_if·A_if/l_seg is
# DERIVED in ff.resolve_intermediate_filaments — no magic number).  The WLC strain-stiffening regime is the
# production direction, but its α→β unfolding crossover (x_max) and enthalpic wall (EA) for IF are a PI GAP
# (Kreplak 2005 10.1016/j.jmb.2005.09.092; Block 2017 10.1103/PhysRevLett.118.048101 give the qualitative
# three regimes, not a ratified card) — this module NEVER defaults them; a caller must supply a sourced card.
IF_LINEAR_TANGENT_REGIME = "LINEAR_TANGENT_REFERENCE"
IF_WLC_STRAIN_STIFFENING_REGIME = "WLC_STRAIN_STIFFENING_UNFOLDING_CARD_REQUIRED"

# Physiological thermal energy k_B·T at 37 °C (310.15 K).  Physical constant, grid-invariant, derivable —
# passes the Magic-Number Block.  1 pN·µm = 1e-18 J, so kBT = 1.380649e-23 J/K · 310.15 K / 1e-18 ≈ 4.28e-3.
KBT_37C_PN_UM = 1.380649e-23 * 310.15 / 1.0e-18

# Single-filament persistence length presets (µm), literature-sourced and band-checked against LP_BAND_UM.
# Keratin K8/K18 (epithelial / MCF7) softer; vimentin (mesenchymal / MDA-MB-231) stiffer — the EMT contrast.
IF_LP_KERATIN_UM = 0.30    # Lichtenstern 2012 10.1016/j.jsb.2011.11.003
IF_LP_VIMENTIN_UM = 0.50   # Lin 2010 10.1016/j.jmb.2010.04.054

# ── Sourced TENSILE (NOT bending) axial-WLC constitutive anchors — vimentin FOUND, keratin NOT FOUND ──
# PI_GAP_LITERATURE_SOURCING_2026-07-23 §3.  Use the AXIAL tensile numbers: bending-derived E (300-900 MPa)
# is ~100× the axial small-strain E (anisotropy artifact) and MUST NOT be used for a tensile cable.  These are
# provenance constants for building a NonlinearCableCard candidate — they are NEVER auto-wired as adapter
# defaults; a caller supplies the card explicitly (so a fabricated crossover cannot slip in).
IF_VIMENTIN_WALL_EA_PN = 13000.0   # ~13 nN post-α→β-unfolding enthalpic EA (hollow A=71.5 nm²),
#                                    Qin/Kreplak/Buehler 2009 PLoS ONE 4:e7294 (toe EA ~0.5-0.9 nN)
IF_VIMENTIN_XMAX_STRAIN = 2.5      # ~200-300% rupture/plateau strain, Block 2018 sciadv.aat1161 /
#                                    Forsting 2019 nanolett.9b02972 — a MODELING NOTE, not the x/Lc crossover
IF_KERATIN_WALL_EA_PN: float | None = None   # NOT FOUND — Lorenz 2019 PRL 123:188102 / 2023 Matter 6:2019
#                                              (K8/K18 EA/x_max) paywalled → PI-GAP slot, no default

_CELL_TYPE_E_PA = {"keratin": IF_E_KERATIN_PA, "vimentin": IF_E_VIMENTIN_PA}
_CELL_TYPE_LP_UM = {"keratin": IF_LP_KERATIN_UM, "vimentin": IF_LP_VIMENTIN_UM}


def _storage_key(array: object) -> tuple[object, ...]:
    """Return storage identity without reading device data."""
    ptr = getattr(array, "ptr", None)
    if ptr is None:
        return ("object", id(array))
    return ("device-ptr", str(getattr(array, "device", None)), int(ptr))


def _validate_device_array(
    array: object,
    *,
    label: str,
    dtype: object,
    ndim: int = 1,
    shape: tuple[int, ...] | None = None,
) -> None:
    """Validate CUDA residence, dtype, dimensionality, and optional exact shape."""
    device = getattr(array, "device", None)
    if not bool(getattr(device, "is_cuda", False)):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != dtype:
        raise TypeError(f"{label} must have dtype {dtype}")
    actual_shape = getattr(array, "shape", None)
    if not isinstance(actual_shape, tuple) or len(actual_shape) != ndim:
        raise ValueError(f"{label} must be {ndim}-dimensional")
    if any(int(size) <= 0 for size in actual_shape):
        raise ValueError(f"{label} must be non-empty")
    if shape is not None and actual_shape != shape:
        raise ValueError(f"{label} must have shape {shape}, got {actual_shape}")


def _require_same_device(arrays: tuple[object, ...], *, label: str) -> None:
    if len({str(getattr(array, "device", None)) for array in arrays}) != 1:
        raise ValueError(f"{label} arrays must share one CUDA device")


@runtime_checkable
class IntermediateFilamentTransaction(TransactionParticipant, Protocol):
    """Transaction delegate proving coverage of authoritative IF state."""

    component_name: str

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Return arrays covered by snapshot/rollback without a device read."""


@dataclass(frozen=True, slots=True)
class IntermediateFilamentRigView:
    """Non-owning material-coordinate view exposed to registered graph connectors."""

    position_d: wp.array
    force_d: wp.array
    fiber_offset_d: wp.array
    retained_count_d: wp.array
    material_epoch_d: wp.array
    mapping_epoch_d: wp.array
    topology_epoch_d: wp.array
    n_filaments: int


@dataclass(frozen=True, slots=True)
class StructuralEndpoint:
    """Live nucleus or SF mechanics endpoint owned by another component."""

    component: str
    position_d: wp.array
    force_d: wp.array

    def __post_init__(self) -> None:
        if self.component not in {NUCLEUS_COMPONENT, SF_COMPONENT}:
            raise ValueError("IF structural endpoint must be nucleus or sf_arc")
        _validate_device_array(
            self.position_d,
            label=f"{self.component}.position_d",
            dtype=wp.vec3d,
        )
        _validate_device_array(
            self.force_d,
            label=f"{self.component}.force_d",
            dtype=wp.vec3d,
            shape=self.position_d.shape,
        )
        _require_same_device((self.position_d, self.force_d), label=self.component)
        if _storage_key(self.position_d) == _storage_key(self.force_d):
            raise ValueError(f"{self.component} position and force arrays must not alias")


@runtime_checkable
class IntermediateFilamentGraphConnector(Protocol):
    """Graph-owned kinetic LINC or plectin joint."""

    name: str
    component_a: str
    component_b: str

    def accumulate_rig(
        self,
        rig: IntermediateFilamentRigView,
        endpoint: StructuralEndpoint,
    ) -> None:
        """Scatter equal-and-opposite joint loads through public device views."""

    def snapshot_candidate(self) -> None:
        """Snapshot graph-owned endpoint and kinetic state."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore connector state under rejection."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit bind/unbind kinetics only after global acceptance."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute joint force/work/topology terms."""


@runtime_checkable
class IntermediateFilamentFluidTransfer(ImmersedTransferConnector, Protocol):
    """Graph-owned adjoint interpolation/scatter and reversible stencil cache."""

    name: str
    component_a: str
    component_b: str

    def accumulate_transfer(
        self,
        rig: IntermediateFilamentRigView,
        cytosol: CytosolFieldEndpoint,
    ) -> None:
        """Gather porous-fluid load and scatter the equal work-conjugate reaction."""

    def snapshot_candidate(self) -> None:
        """Snapshot mapping/stencil caches before candidate geometry changes."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore cached mappings under rejection."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Finalize accepted non-kinetic caches; no private kinetics are advanced."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute fluid/solid transfer work and residual terms."""


@dataclass(frozen=True, slots=True)
class IntermediateFilamentStepBindings:
    """External graph endpoints and connector state supplied for one candidate step."""

    nucleus: StructuralEndpoint
    sf_arc: StructuralEndpoint
    cytosol: CytosolFieldEndpoint
    nucleus_linc: IntermediateFilamentGraphConnector
    sf_plectin: IntermediateFilamentGraphConnector
    cytosol_transfer: IntermediateFilamentFluidTransfer

    def __post_init__(self) -> None:
        if self.nucleus.component != NUCLEUS_COMPONENT:
            raise ValueError("nucleus endpoint must be the registered nucleus component")
        if self.sf_arc.component != SF_COMPONENT:
            raise ValueError("SF endpoint must be the registered sf_arc component")
        _validate_runtime_connector(
            self.nucleus_linc,
            expected_name=IF_NUCLEUS_LINC,
            endpoints={INTERMEDIATE_FILAMENT_COMPONENT, NUCLEUS_COMPONENT},
        )
        _validate_runtime_connector(
            self.sf_plectin,
            expected_name=IF_SF_PLECTIN,
            endpoints={INTERMEDIATE_FILAMENT_COMPONENT, SF_COMPONENT},
        )
        _validate_runtime_connector(
            self.cytosol_transfer,
            expected_name=IF_CYTOSOL_TRANSFER,
            endpoints={INTERMEDIATE_FILAMENT_COMPONENT, CYTOSOL_COMPONENT},
        )


def _validate_runtime_connector(
    connector: object,
    *,
    expected_name: str,
    endpoints: set[str],
) -> None:
    if getattr(connector, "name", None) != expected_name:
        raise ValueError(f"expected graph connector {expected_name!r}")
    if {
        getattr(connector, "component_a", None),
        getattr(connector, "component_b", None),
    } != endpoints:
        raise ValueError(f"{expected_name} has incorrect component endpoints")


@dataclass(frozen=True, slots=True)
class IntermediateFilamentRigStateOwner:
    """Authoritative component-local CUDA state for a nonlinear IF cable graph."""

    n_filaments: int
    position_d: wp.array
    force_d: wp.array
    fiber_offset_d: wp.array
    retained_count_d: wp.array
    constitutive_state_d: wp.array
    turnover_state_d: wp.array
    material_epoch_d: wp.array
    internal_crosslink_pair_d: wp.array
    internal_crosslink_state_d: wp.array
    mapping_epoch_d: wp.array
    topology_epoch_d: wp.array
    mechanics: MechanicsContributor
    transaction: IntermediateFilamentTransaction
    ledger: LedgerContributor

    def __post_init__(self) -> None:
        if self.n_filaments <= 0:
            raise ValueError("n_filaments must be positive")
        _validate_device_array(
            self.position_d,
            label="intermediate_filament.position_d",
            dtype=wp.vec3d,
        )
        n_nodes = int(self.position_d.shape[0])
        _validate_device_array(
            self.force_d,
            label="intermediate_filament.force_d",
            dtype=wp.vec3d,
            shape=(n_nodes,),
        )
        _validate_device_array(
            self.fiber_offset_d,
            label="intermediate_filament.fiber_offset_d",
            dtype=wp.int32,
            shape=(self.n_filaments + 1,),
        )
        for label, array in (
            ("retained_count_d", self.retained_count_d),
            ("turnover_state_d", self.turnover_state_d),
            ("material_epoch_d", self.material_epoch_d),
        ):
            _validate_device_array(
                array,
                label=f"intermediate_filament.{label}",
                dtype=wp.int32,
                shape=(self.n_filaments,),
            )
        _validate_device_array(
            self.constitutive_state_d,
            label="intermediate_filament.constitutive_state_d",
            dtype=wp.float64,
        )
        _validate_device_array(
            self.internal_crosslink_pair_d,
            label="intermediate_filament.internal_crosslink_pair_d",
            dtype=wp.vec2i,
        )
        n_crosslink_slots = int(self.internal_crosslink_pair_d.shape[0])
        _validate_device_array(
            self.internal_crosslink_state_d,
            label="intermediate_filament.internal_crosslink_state_d",
            dtype=wp.int32,
            shape=(n_crosslink_slots,),
        )
        for label, array in (
            ("mapping_epoch_d", self.mapping_epoch_d),
            ("topology_epoch_d", self.topology_epoch_d),
        ):
            _validate_device_array(
                array,
                label=f"intermediate_filament.{label}",
                dtype=wp.int32,
                shape=(1,),
            )
        all_arrays = (
            self.position_d,
            self.force_d,
            self.fiber_offset_d,
            self.retained_count_d,
            self.constitutive_state_d,
            self.turnover_state_d,
            self.material_epoch_d,
            self.internal_crosslink_pair_d,
            self.internal_crosslink_state_d,
            self.mapping_epoch_d,
            self.topology_epoch_d,
        )
        _require_same_device(all_arrays, label="intermediate-filament")
        if len({_storage_key(array) for array in all_arrays}) != len(all_arrays):
            raise ValueError("IF state, force, topology, and history arrays must not alias")
        if getattr(self.mechanics, "component_local", True) is False:
            raise ValueError("monolithic legacy IF mechanics cannot serve a component-local state owner")
        n_expected = getattr(self.mechanics, "n_nodes", None)
        if n_expected is not None and int(n_expected) != n_nodes:
            raise ValueError(f"IF mechanics expects {int(n_expected)} nodes but state owns {n_nodes}")
        if self.transaction.component_name != INTERMEDIATE_FILAMENT_COMPONENT:
            raise ValueError("IF transaction delegate must be owned by intermediate_filament")
        mutable = (
            self.position_d,
            self.retained_count_d,
            self.constitutive_state_d,
            self.turnover_state_d,
            self.material_epoch_d,
            self.internal_crosslink_pair_d,
            self.internal_crosslink_state_d,
            self.mapping_epoch_d,
            self.topology_epoch_d,
        )
        covered = {_storage_key(array) for array in self.transaction.owned_arrays()}
        if {_storage_key(array) for array in mutable} - covered:
            raise ValueError(
                "IF transaction must cover geometry, nonlinear history, reduction, turnover, "
                "crosslinks, and epochs"
            )

    def geometry(self) -> IntermediateFilamentRigView:
        """Return the non-owning material-coordinate connector view."""
        return IntermediateFilamentRigView(
            position_d=self.position_d,
            force_d=self.force_d,
            fiber_offset_d=self.fiber_offset_d,
            retained_count_d=self.retained_count_d,
            material_epoch_d=self.material_epoch_d,
            mapping_epoch_d=self.mapping_epoch_d,
            topology_epoch_d=self.topology_epoch_d,
            n_filaments=self.n_filaments,
        )

    def accumulate(self) -> None:
        """Launch component-local cable mechanics."""
        self.mechanics.accumulate(self.position_d, self.force_d)

    def snapshot_candidate(self) -> None:
        """Snapshot all mutable authoritative IF state."""
        self.transaction.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        """Restore IF state under the scheduler-owned rejected predicate."""
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit turnover/crosslinks/history only under final acceptance."""
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute IF mechanics/topology/history terms."""
        self.ledger.accumulate_ledger(ledger)


def _find_connector(architecture: CellArchitecture, name: str) -> ConnectorContract:
    for connector in architecture.connectors:
        if connector.name == name:
            return connector
    raise ValueError(f"architecture must register {name!r}")


def _validate_build_connector(
    architecture: CellArchitecture,
    *,
    name: str,
    family: ConnectorFamily,
    endpoints: set[str],
    kinetic: bool,
) -> None:
    connector = _find_connector(architecture, name)
    if connector.family is not family:
        raise ValueError(f"{name} must use ConnectorFamily.{family.name}")
    if {connector.component_a, connector.component_b} != endpoints:
        raise ValueError(f"{name} has incorrect component endpoints")
    if connector.kinetics is not kinetic or connector.commit_on_accept is not kinetic:
        qualifier = "kinetic and accepted-step committed" if kinetic else "non-kinetic"
        raise ValueError(f"{name} must be {qualifier}")
    if not connector.bidirectional or not connector.adjoint_transfer_required:
        raise ValueError(f"{name} must be bidirectional with adjoint transfer")


def _validate_architecture(architecture: CellArchitecture) -> None:
    component = architecture.component(INTERMEDIATE_FILAMENT_COMPONENT)
    if component.role is not ComponentRole.STRUCTURAL_RIG:
        raise ValueError("intermediate_filament must have ComponentRole.STRUCTURAL_RIG")
    if not component.owns_geometry or not component.dynamically_evolving:
        raise ValueError("intermediate_filament must own geometry and be dynamically evolving")
    _validate_build_connector(
        architecture,
        name=IF_NUCLEUS_LINC,
        family=ConnectorFamily.LINC,
        endpoints={INTERMEDIATE_FILAMENT_COMPONENT, NUCLEUS_COMPONENT},
        kinetic=True,
    )
    _validate_build_connector(
        architecture,
        name=IF_SF_PLECTIN,
        family=ConnectorFamily.PLECTIN,
        endpoints={INTERMEDIATE_FILAMENT_COMPONENT, SF_COMPONENT},
        kinetic=True,
    )
    _validate_build_connector(
        architecture,
        name=IF_CYTOSOL_TRANSFER,
        family=ConnectorFamily.IMMERSED_TRANSFER,
        endpoints={INTERMEDIATE_FILAMENT_COMPONENT, CYTOSOL_COMPONENT},
        kinetic=False,
    )


@dataclass(frozen=True, slots=True)
class IntermediateFilamentRig:
    """Component facade; all cross-component state is supplied by step bindings."""

    architecture: CellArchitecture
    state: IntermediateFilamentRigStateOwner

    def __post_init__(self) -> None:
        _validate_architecture(self.architecture)

    def accumulate_mechanics(self, bindings: IntermediateFilamentStepBindings) -> None:
        """Launch IF mechanics followed by each registered connector contribution."""
        self.state.accumulate()
        view = self.state.geometry()
        bindings.nucleus_linc.accumulate_rig(view, bindings.nucleus)
        bindings.sf_plectin.accumulate_rig(view, bindings.sf_arc)
        bindings.cytosol_transfer.accumulate_transfer(view, bindings.cytosol)

    def transaction_participants(
        self,
        bindings: IntermediateFilamentStepBindings,
    ) -> tuple[TransactionParticipant, ...]:
        """Return state, kinetic joints, and reversible fluid-transfer cache in fixed order."""
        return (
            self.state,
            bindings.nucleus_linc,
            bindings.sf_plectin,
            bindings.cytosol_transfer,
        )

    def snapshot_candidate(self, bindings: IntermediateFilamentStepBindings) -> None:
        """Snapshot every IF-connected transaction participant."""
        for participant in self.transaction_participants(bindings):
            participant.snapshot_candidate()

    def rollback(self, bindings: IntermediateFilamentStepBindings, accepted: wp.array) -> None:
        """Forward one scheduler-owned predicate to every participant."""
        for participant in self.transaction_participants(bindings):
            participant.rollback(accepted)

    def commit_irreversible(
        self,
        bindings: IntermediateFilamentStepBindings,
        accepted: wp.array,
        dt_phys: float,
        rng_seed: int,
    ) -> None:
        """Commit topology/history/joint state once after global acceptance."""
        for participant in self.transaction_participants(bindings):
            participant.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, bindings: IntermediateFilamentStepBindings, ledger: object) -> None:
        """Collect component and connector ledger terms without deciding acceptance."""
        self.state.accumulate_ledger(ledger)
        bindings.nucleus_linc.accumulate_ledger(ledger)
        bindings.sf_plectin.accumulate_ledger(ledger)
        bindings.cytosol_transfer.accumulate_ledger(ledger)


class _LegacyIntermediateFilamentCage(Protocol):
    n_if: int
    n_bonds: int
    n_backbone: int
    n_linc: int
    n_anchor: int
    n_xl: int

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch the legacy global-index Hookean cage kernel."""


@dataclass(frozen=True, slots=True)
class LegacyIntermediateFilamentCageAdapter:
    """Reference-only adapter for the monolithic Hookean IF cage.

    The delegate is intentionally not component-local: its global bond array embeds external LINC and cortex
    anchors plus internal crosslinks. :class:`IntermediateFilamentRigStateOwner` rejects it as a backend.
    """

    compartment: _LegacyIntermediateFilamentCage

    def __post_init__(self) -> None:
        counts = (
            int(self.compartment.n_backbone),
            int(self.compartment.n_linc),
            int(self.compartment.n_anchor),
            int(self.compartment.n_xl),
        )
        if int(self.compartment.n_if) <= 0 or int(self.compartment.n_bonds) <= 0:
            raise ValueError("legacy IF adapter requires a non-empty cage")
        if any(count < 0 for count in counts) or sum(counts) != int(self.compartment.n_bonds):
            raise ValueError("legacy IF bond-family counts must exactly partition n_bonds")

    @property
    def component_local(self) -> bool:
        """Declare that global-index connector bonds prevent component-local use."""
        return False

    @property
    def status(self) -> str:
        """Prevent the Hookean reference from being reported as the completed IF model."""
        return LEGACY_REFERENCE_STATUS

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Delegate one legacy reference force evaluation on its combined global arrays."""
        self.compartment.accumulate(pos, force)


def keratin_vimentin_backbone_stiffness(cell_type: str, *, n_fil: int, l_seg_um: float) -> float:
    """Return the DERIVED backbone axial stiffness k_bb [pN/µm] for a keratin/vimentin cage.

    Pure reuse of :func:`ff.resolve_intermediate_filaments` (k_bb = E_if·A_if/l_seg); introduces no new
    constant and range-checks E against the sourced band.  Used by construction to fill the linear adapter's
    per-segment stiffness; it is not a tuning knob.
    """
    if cell_type not in _CELL_TYPE_E_PA:
        raise ValueError("cell_type must be 'keratin' or 'vimentin'")
    resolved = resolve_intermediate_filaments(
        n_fil=n_fil, E_if_Pa=_CELL_TYPE_E_PA[cell_type], l_seg_um=l_seg_um
    )
    return float(resolved.k_bb)


def keratin_vimentin_persistence_length(cell_type: str) -> float:
    """Return the sourced single-filament persistence length Lp [µm], band-checked against LP_BAND_UM."""
    if cell_type not in _CELL_TYPE_LP_UM:
        raise ValueError("cell_type must be 'keratin' or 'vimentin'")
    lp = _CELL_TYPE_LP_UM[cell_type]
    if not (LP_BAND_UM[0] <= lp <= LP_BAND_UM[1]):
        raise ValueError(f"Lp {lp} µm outside sourced band {LP_BAND_UM}")
    return float(lp)


def _validate_local_backbone_segments(
    segment_d: object,
    *,
    n_nodes: int,
    node_offset: int,
) -> int:
    """Validate a component-local backbone-only segment table; return the segment count.

    A component-local backbone carries only IF↔IF backbone bonds with node indices in ``[0, n_nodes)``.
    Global-offset tables (which embed LINC/anchor/crosslink connector bonds) are rejected: index residence
    cannot be device-read here, so ``node_offset == 0`` is required as the local-indexing provenance flag,
    mirroring the microtubule bending adapter.
    """
    if node_offset != 0:
        raise ValueError(
            "component-local IF backbone adapter requires node_offset == 0 (global-index cages carry "
            "graph-owned LINC/anchor/crosslink bonds and are rejected)"
        )
    _validate_device_array(
        segment_d,
        label="if_backbone.segment_d",
        dtype=wp.int32,
        ndim=2,
    )
    if int(segment_d.shape[1]) != 2:
        raise ValueError("if_backbone.segment_d must have shape (n_segments, 2)")
    return int(segment_d.shape[0])


@dataclass(frozen=True, slots=True)
class LinearBackboneCableAdapter:
    """Component-local backbone-only Hookean cable bound to the reused ``ff.link_spring_kernel``.

    This is the sourced LINEAR TANGENT LIMIT of the IF constitutive law (derived k_bb, no magic number).  It
    carries only backbone segments with LOCAL node indices; LINC, cortex anchor, and internal crosslinks are
    NOT here (graph- and component-crosslink-owned).  Unlike :class:`LegacyIntermediateFilamentCageAdapter`
    it advertises ``component_local == True`` and a matching ``n_nodes``, so
    :class:`IntermediateFilamentRigStateOwner` accepts it as a production backbone mechanics backend.
    """

    n_nodes: int
    segment_d: wp.array
    stiffness_d: wp.array
    rest_d: wp.array
    device: str
    cell_type: str
    node_offset: int = 0

    def __post_init__(self) -> None:
        if self.n_nodes <= 0:
            raise ValueError("n_nodes must be positive")
        if self.cell_type not in _CELL_TYPE_E_PA:
            raise ValueError("cell_type must be 'keratin' or 'vimentin'")
        n_seg = _validate_local_backbone_segments(
            self.segment_d, n_nodes=self.n_nodes, node_offset=self.node_offset
        )
        _validate_device_array(
            self.stiffness_d,
            label="if_backbone.stiffness_d",
            dtype=wp.float64,
            shape=(n_seg,),
        )
        _validate_device_array(
            self.rest_d,
            label="if_backbone.rest_d",
            dtype=wp.float64,
            shape=(n_seg,),
        )
        _require_same_device(
            (self.segment_d, self.stiffness_d, self.rest_d), label="if_backbone"
        )
        if len(
            {_storage_key(a) for a in (self.segment_d, self.stiffness_d, self.rest_d)}
        ) != 3:
            raise ValueError("IF backbone segment/stiffness/rest arrays must not alias")

    @property
    def component_local(self) -> bool:
        """Declare backbone-local ownership so the state owner installs it as a production backend."""
        return True

    @property
    def regime(self) -> str:
        """Name the sourced constitutive regime bound by this adapter."""
        return IF_LINEAR_TANGENT_REGIME

    @property
    def bound_kernel(self) -> object:
        """Expose the reused Warp kernel this adapter launches (no reimplemented physics)."""
        return link_spring_kernel

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch the reused Hookean backbone kernel over component-local segments."""
        n_seg = int(self.segment_d.shape[0])
        if n_seg:
            wp.launch(
                link_spring_kernel,
                dim=n_seg,
                inputs=[pos, self.segment_d, self.stiffness_d, self.rest_d],
                outputs=[force],
                device=self.device,
            )


@dataclass(frozen=True, slots=True)
class NonlinearCableCard:
    """Sourced three-regime WLC constitutive card for one IF material.

    ``persistence_length_um`` and ``thermal_energy_pn_um`` are sourced/derivable, but ``axial_stiffness_pn``
    (EA enthalpic wall) and ``crossover_ratio`` (x_max entropic→enthalpic transition for α→β unfolding) are a
    PI GAP for intermediate filaments and MUST be supplied explicitly — this card provides no defaults, so a
    fabricated crossover cannot slip in.
    """

    persistence_length_um: float
    axial_stiffness_pn: float
    thermal_energy_pn_um: float
    crossover_ratio: float

    def __post_init__(self) -> None:
        if not (LP_BAND_UM[0] <= self.persistence_length_um <= LP_BAND_UM[1]):
            raise ValueError(f"persistence_length_um outside sourced band {LP_BAND_UM}")
        if self.axial_stiffness_pn <= 0.0:
            raise ValueError("axial_stiffness_pn (EA) must be positive [pN]")
        if self.thermal_energy_pn_um <= 0.0:
            raise ValueError("thermal_energy_pn_um (kBT) must be positive [pN·µm]")
        if not (0.0 < self.crossover_ratio < 1.0):
            raise ValueError("crossover_ratio (x_max) must lie in (0, 1)")


@dataclass(frozen=True, slots=True)
class NonlinearWlcCableAdapter:
    """Component-local backbone-only extensible-WLC cable bound to the reused ``ff.wlc_spring_kernel``.

    This is the production strain-stiffening direction (soft entropic coil → finite-extensibility enthalpic
    wall).  It is component-local exactly like the linear adapter, but requires an explicit
    :class:`NonlinearCableCard`; the IF α→β unfolding plateau/crossover is unresolved, so the caller — not
    this module — supplies the sourced card.  Per-segment contour length replaces rest length.
    """

    n_nodes: int
    segment_d: wp.array
    contour_length_d: wp.array
    card: NonlinearCableCard
    device: str
    cell_type: str
    node_offset: int = 0

    def __post_init__(self) -> None:
        if self.n_nodes <= 0:
            raise ValueError("n_nodes must be positive")
        if self.cell_type not in _CELL_TYPE_E_PA:
            raise ValueError("cell_type must be 'keratin' or 'vimentin'")
        n_seg = _validate_local_backbone_segments(
            self.segment_d, n_nodes=self.n_nodes, node_offset=self.node_offset
        )
        _validate_device_array(
            self.contour_length_d,
            label="if_backbone.contour_length_d",
            dtype=wp.float64,
            shape=(n_seg,),
        )
        _require_same_device((self.segment_d, self.contour_length_d), label="if_backbone")
        if _storage_key(self.segment_d) == _storage_key(self.contour_length_d):
            raise ValueError("IF backbone segment/contour arrays must not alias")

    @property
    def component_local(self) -> bool:
        """Declare backbone-local ownership so the state owner installs it as a production backend."""
        return True

    @property
    def regime(self) -> str:
        """Name the production constitutive regime; the label records the unresolved-card requirement."""
        return IF_WLC_STRAIN_STIFFENING_REGIME

    @property
    def bound_kernel(self) -> object:
        """Expose the reused Warp WLC kernel this adapter launches (no reimplemented physics)."""
        return wlc_spring_kernel

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch the reused extensible-WLC backbone kernel over component-local segments."""
        n_seg = int(self.segment_d.shape[0])
        if n_seg:
            wp.launch(
                wlc_spring_kernel,
                dim=n_seg,
                inputs=[
                    pos,
                    self.segment_d,
                    self.contour_length_d,
                    wp.float64(self.card.persistence_length_um),
                    wp.float64(self.card.axial_stiffness_pn),
                    wp.float64(self.card.thermal_energy_pn_um),
                    wp.float64(self.card.crossover_ratio),
                ],
                outputs=[force],
                device=self.device,
            )


def vimentin_nonlinear_cable_card(*, crossover_ratio: float) -> NonlinearCableCard:
    """Build a vimentin WLC card from the SOURCED tensile anchors, leaving the x/Lc crossover as a caller slot.

    ``persistence_length_um`` (Lin 2010) and ``axial_stiffness_pn`` (Qin 2009 post-unfolding EA) are FOUND and
    filled here; ``thermal_energy_pn_um`` is the derived 37 °C kBT.  ``crossover_ratio`` (the entropic→enthalpic
    x_max as a fraction of contour length) is NOT determined by the sourced ~200-300% *strain* number — that is
    a rest-relative rupture strain, not the WLC x/Lc crossover — so the mapping stays a caller-supplied modeling
    choice rather than being fabricated from an incompatible quantity.
    """
    return NonlinearCableCard(
        persistence_length_um=IF_LP_VIMENTIN_UM,
        axial_stiffness_pn=IF_VIMENTIN_WALL_EA_PN,
        thermal_energy_pn_um=KBT_37C_PN_UM,
        crossover_ratio=crossover_ratio,
    )


def keratin_nonlinear_cable_card(*, crossover_ratio: float) -> NonlinearCableCard:
    """Keratin (primary MCF7-epithelial IF) WLC card — BLOCKED: EA/x_max are a PI-GAP.

    Raises unconditionally: the definitive K8/K18 optical-tweezer EA/x_max (Lorenz 2019/2023) are paywalled and
    NOT FOUND, so there is no sourced ``axial_stiffness_pn`` to fill.  A caller that has a PI-ratified keratin
    card must construct :class:`NonlinearCableCard` directly; this function refuses to invent the wall stiffness.
    """
    raise NotImplementedError(
        "keratin K8/K18 axial EA/x_max NOT FOUND (Lorenz 2019/2023 paywalled) — PI-GAP; "
        "supply a ratified NonlinearCableCard explicitly"
    )


@dataclass(frozen=True, slots=True)
class IntermediateFilamentTurnoverCard:
    """Subunit-exchange turnover kinetics for the IF backbone — an UNSET PI-GAP slot.

    IF turnover lets the network dissipate stored strain by subunit exchange (rest length creeps toward the
    current length; each segment becomes a Maxwell element).  Both rates have NO default: no IF-specific
    ``k_off0`` / Bell force scale is sourced in ``PI_GAP_LITERATURE_SOURCING_2026-07-23``, so a caller supplies
    a ratified card.  Semantics reuse the established viscoelastic-turnover primitive (NF2007 §10.1 Bell slip,
    ``ff.network_warp.xl_turnover_kernel``) exactly — this card is only the two unset numbers.
    """

    k_off0_per_s: float
    x_beta_over_kt_per_pn: float
    provenance: str
    ratified: bool = False

    def __post_init__(self) -> None:
        for label, value in (
            ("k_off0_per_s", self.k_off0_per_s),
            ("x_beta_over_kt_per_pn", self.x_beta_over_kt_per_pn),
        ):
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{label} must be finite and nonnegative")
        if not self.provenance.strip():
            raise ValueError("provenance must record the rate source (or an explicit PI-GAP marker)")


@dataclass(frozen=True, slots=True)
class BackboneTurnoverKinetics:
    """Bind the reused ``xl_turnover_kernel`` to a component-local backbone's segment/stiffness/rest arrays.

    Turnover is an accepted-step irreversible kinetic: it creeps ``rest_d`` IN PLACE toward the current segment
    length under a Bell-slip off-rate, so the caller must launch it only under an accepted predicate and cover
    ``rest_d`` in the IF transaction snapshot/rollback.  The physics kernel is REUSED, not re-derived; this
    adapter only supplies the actor-owned arrays and the UNSET-slot :class:`IntermediateFilamentTurnoverCard`.
    """

    n_nodes: int
    segment_d: wp.array
    stiffness_d: wp.array
    rest_d: wp.array
    card: IntermediateFilamentTurnoverCard
    device: str

    def __post_init__(self) -> None:
        if self.n_nodes <= 0:
            raise ValueError("n_nodes must be positive")
        n_seg = _validate_local_backbone_segments(self.segment_d, n_nodes=self.n_nodes, node_offset=0)
        for label, array in (("stiffness_d", self.stiffness_d), ("rest_d", self.rest_d)):
            _validate_device_array(
                array,
                label=f"if_turnover.{label}",
                dtype=wp.float64,
                shape=(n_seg,),
            )
        _require_same_device((self.segment_d, self.stiffness_d, self.rest_d), label="if_turnover")
        if len({_storage_key(a) for a in (self.segment_d, self.stiffness_d, self.rest_d)}) != 3:
            raise ValueError("IF turnover segment/stiffness/rest arrays must not alias")

    @property
    def bound_kernel(self) -> object:
        """Expose the reused Warp kernel this adapter launches (no reimplemented physics)."""
        return xl_turnover_kernel

    def apply_turnover(self, pos: wp.array, dt_real: float) -> None:
        """Launch one reused viscoelastic-turnover pass over component-local backbone segments.

        ``dt_real`` is the outer physical-time increment; ``k_off0·dt_real`` is the per-event turnover fraction.
        Mutates ``rest_d`` in place — invoke only on an accepted step.
        """
        if not (dt_real > 0.0):
            raise ValueError("dt_real must be positive")
        n_seg = int(self.segment_d.shape[0])
        if n_seg:
            wp.launch(
                xl_turnover_kernel,
                dim=n_seg,
                inputs=[
                    pos,
                    self.segment_d,
                    self.stiffness_d,
                    self.rest_d,
                    wp.float64(self.card.k_off0_per_s * dt_real),
                    wp.float64(self.card.x_beta_over_kt_per_pn),
                ],
                device=self.device,
            )


def backbone_turnover_reference(
    *,
    length_um: float,
    rest_um: float,
    stiffness_pn_per_um: float,
    k_off0_per_s: float,
    x_beta_over_kt_per_pn: float,
    dt_real: float,
) -> float:
    """Pure-Python mirror of one ``xl_turnover_kernel`` event: return the crept rest length [µm].

    CPU-green ground truth for :meth:`BackboneTurnoverKinetics.apply_turnover`.  Force ``F=k(L-r0)``; Bell-slip
    per-event off-rate ``k_off·dt = k_off0·dt·exp(|F|·xβ)``; ensemble fraction ``1-exp(-k_off·dt)`` rebinds
    force-free at the current geometry, so ``r0`` creeps toward ``L``.
    """
    if dt_real <= 0.0:
        raise ValueError("dt_real must be positive")
    if stiffness_pn_per_um < 0.0 or k_off0_per_s < 0.0 or x_beta_over_kt_per_pn < 0.0:
        raise ValueError("stiffness, k_off0, and Bell force scale must be nonnegative")
    force_pn = stiffness_pn_per_um * (length_um - rest_um)
    koff_dt = k_off0_per_s * dt_real * np.exp(abs(force_pn) * x_beta_over_kt_per_pn)
    frac = 1.0 - np.exp(-koff_dt)
    return float(rest_um + frac * (length_um - rest_um))


__all__ = [
    "BackboneTurnoverKinetics",
    "CYTOSOL_COMPONENT",
    "CytosolFieldEndpoint",
    "IF_CYTOSOL_TRANSFER",
    "IF_KERATIN_WALL_EA_PN",
    "IF_LINEAR_TANGENT_REGIME",
    "IF_LP_KERATIN_UM",
    "IF_LP_VIMENTIN_UM",
    "IF_NUCLEUS_LINC",
    "IF_SF_PLECTIN",
    "IF_VIMENTIN_WALL_EA_PN",
    "IF_VIMENTIN_XMAX_STRAIN",
    "IF_WLC_STRAIN_STIFFENING_REGIME",
    "INTERMEDIATE_FILAMENT_COMPONENT",
    "IntermediateFilamentFluidTransfer",
    "IntermediateFilamentGraphConnector",
    "IntermediateFilamentRig",
    "IntermediateFilamentRigStateOwner",
    "IntermediateFilamentRigView",
    "IntermediateFilamentStepBindings",
    "IntermediateFilamentTransaction",
    "IntermediateFilamentTurnoverCard",
    "KBT_37C_PN_UM",
    "LEGACY_REFERENCE_STATUS",
    "LegacyIntermediateFilamentCageAdapter",
    "LinearBackboneCableAdapter",
    "NUCLEUS_COMPONENT",
    "NonlinearCableCard",
    "NonlinearWlcCableAdapter",
    "SF_COMPONENT",
    "StructuralEndpoint",
    "backbone_turnover_reference",
    "keratin_nonlinear_cable_card",
    "keratin_vimentin_backbone_stiffness",
    "keratin_vimentin_persistence_length",
    "vimentin_nonlinear_cable_card",
]
