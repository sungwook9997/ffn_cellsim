"""CUDA-resident extracellular-matrix topology primitives.

This package owns collagen geometry/topology only.  Constitutive mechanics, crosslink
chemistry, FA clutch state, and coupled-core scheduling remain separate owners.
"""

from aleph.components.ecm.device_schema import ECMTopologyState, PopulationSlot
from aleph.components.ecm.mikado_topology import (
    MikadoInitConfig,
    MikadoTopologyBuilder,
    SegmentContactCandidates,
)
from aleph.components.ecm.remodel_runtime import ECMRemodelRuntime
from aleph.components.ecm.remodel_schema import (
    ECMEndpointRemapView,
    ECMRemodelPlan,
    ECMRemodelProposal,
    FiberRemodelEvent,
    ProposalStatus,
    RemapKind,
    RemodelEvent,
)

__all__ = [
    "ECMTopologyState",
    "ECMEndpointRemapView",
    "ECMRemodelPlan",
    "ECMRemodelProposal",
    "ECMRemodelRuntime",
    "FiberRemodelEvent",
    "MikadoInitConfig",
    "MikadoTopologyBuilder",
    "PopulationSlot",
    "ProposalStatus",
    "RemapKind",
    "RemodelEvent",
    "SegmentContactCandidates",
]
