"""Component/connector architecture for the dynamically connected Active Cell engine."""

from aleph.engine.actor import CellActor
from aleph.engine.composition import (
    ComposedCellRuntimes,
    ComposedCellWorld,
    FacadeDispatch,
    build_composed_cell_world,
    canonical_facade_phase,
)
from aleph.engine.composed_native import (
    ComposedNativeCell,
    SFArcStateOwner,
    build_native_composed_cell_world,
    build_sf_arc_state_owner,
    compose_native_cell_world,
)
from aleph.engine.contracts import (
    EVIDENCE_ORDER,
    LADDER_FLOOR,
    OSMOTIC_ENVELOPE,
    OSMOTIC_LOADED_SURFACE,
    CellArchitecture,
    ComponentContract,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
    ConnectorScope,
    EvidenceLabel,
    EvidenceRung,
    GateVerdict,
    QuantitativeClaim,
    VoidCeiling,
    classify_gate,
    reference_cell_architecture,
    rung_rank,
)
from aleph.engine.dump_state import (
    ComposedWorldDump,
    ComposedWorldEntryDump,
    dump_composed_world_state,
)
from aleph.engine.cell_state import (
    CellState,
    Distribution,
    Rate,
    assert_cellstate_writes_no_force,
    mcf7_reference_state,
)
from aleph.engine.events import WarpEventClock, make_event_clock
from aleph.engine.ledger import GlobalCellLedger, make_global_cell_ledger
from aleph.engine.population import PopulationLedger, assert_disjoint_populations
from aleph.engine.runtime import (
    EventClock,
    EventProposer,
    LedgerContributor,
    MechanicsContributor,
    TransactionParticipant,
)
from aleph.engine.transaction import CellTransaction
from aleph.engine.world import CellWorldTransaction

__all__ = [
    "assert_cellstate_writes_no_force",
    "assert_disjoint_populations",
    "build_composed_cell_world",
    "build_native_composed_cell_world",
    "build_sf_arc_state_owner",
    "canonical_facade_phase",
    "CellActor",
    "CellArchitecture",
    "CellState",
    "CellTransaction",
    "CellWorldTransaction",
    "classify_gate",
    "ComponentContract",
    "ComponentRole",
    "compose_native_cell_world",
    "ComposedCellRuntimes",
    "ComposedCellWorld",
    "ComposedNativeCell",
    "ComposedWorldDump",
    "ComposedWorldEntryDump",
    "ConnectorContract",
    "ConnectorFamily",
    "ConnectorScope",
    "Distribution",
    "dump_composed_world_state",
    "EnsembleRunner",
    "EnsembleSummary",
    "EventClock",
    "EventLog",
    "EventProposer",
    "EventRecord",
    "EVIDENCE_ORDER",
    "EvidenceLabel",
    "EvidenceRung",
    "FacadeDispatch",
    "GateVerdict",
    "GlobalCellLedger",
    "LADDER_FLOOR",
    "LedgerContributor",
    "make_event_clock",
    "make_global_cell_ledger",
    "mcf7_reference_state",
    "MechanicsContributor",
    "mix_ensembles",
    "ObservableDistribution",
    "OSMOTIC_ENVELOPE",
    "OSMOTIC_LOADED_SURFACE",
    "PopulationLedger",
    "QuantitativeClaim",
    "Rate",
    "RealizationResult",
    "reference_cell_architecture",
    "rung_rank",
    "SFArcStateOwner",
    "TransactionParticipant",
    "VoidCeiling",
    "WarpEventClock",
]


#: Re-exported from `ac.engine.ensemble` on demand.  Eagerly importing it cost every native run 566
#: lines for an API no simulation path touches — `ac/engine/__init__` is imported by anything that
#: reaches the engine at all, so a re-export here is a tax on the whole tree (2026-07-28).
_LAZY_ENSEMBLE = ['EnsembleRunner', 'EnsembleSummary', 'EventLog', 'EventRecord', 'ObservableDistribution', 'RealizationResult', 'mix_ensembles']


def __getattr__(name: str):
    """Resolve the deferred `ensemble` re-exports (PEP 562), so `from ac.engine import X` still works."""
    if name in _LAZY_ENSEMBLE:
        import importlib

        return getattr(importlib.import_module("aleph.engine.ensemble"), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
