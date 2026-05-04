"""V2 image-constrained cell-resolved mechanobiology package.

The v2 package is intentionally independent of Taichi/GPU at import time. It
starts with data contracts and biological state schemas, then grows toward
single-cell and cell-resolved spheroid simulation.
"""

from acs.v2.active_contour import (
    ActiveContourParameters,
    ActiveContourParametersError,
    ActiveContourState,
    compute_rate_max,
    ellipse_polygon_vertices,
    regular_polygon_vertices,
)
from acs.v2.active_contour_harness import (
    HarnessRun,
    StepDiagnostics,
    compute_finite_n_residual,
    equilibrium_radius_from_cubic,
    equilibrium_radius_from_quartic,
    run_active_contour_test,
)
from acs.v2.dynamics.active_contour import (
    ActiveContourStepError,
    compute_area_forces,
    compute_cortex_forces,
    energy_area,
    energy_cortex,
    step,
)
from acs.v2.dynamics.ecm_open_loop import (
    ECMOpenLoopError,
    accumulate_prescribed_traction,
)
from acs.v2.cell_cluster import CellClusterState
from acs.v2.data_contract import (
    ArtifactKind,
    ImagingDatasetSpec,
    MetricSpec,
    V2DataContract,
)
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState
from acs.v2.junction import JunctionState
from acs.v2.measurement_boundary import (
    MeasurementBoundary,
    MeasurementBoundaryError,
)
from acs.v2.metrics import (
    MetricRegistry,
    MetricRegistryError,
    RegisteredMetric,
    boundary_perimeter_um,
    boundary_projected_area_um2,
    default_registry,
)
from acs.v2.output.frame_dump import (
    FrameDumpReadResult,
    read_frame,
    write_frame,
)
from acs.v2.protrusion import ProtrusionEvent
from acs.v2.single_cell import SingleCellState
from acs.v2.viz.stub3d import render_frame_html, render_frame_png

__all__ = [
    "ActiveContourParameters",
    "ActiveContourParametersError",
    "ActiveContourState",
    "ActiveContourStepError",
    "ArtifactKind",
    "CellClusterState",
    "ECMSubstrateState",
    "ECMOpenLoopError",
    "FocalAdhesionState",
    "FrameDumpReadResult",
    "HarnessRun",
    "ImagingDatasetSpec",
    "JunctionState",
    "MeasurementBoundary",
    "MeasurementBoundaryError",
    "MetricRegistry",
    "MetricRegistryError",
    "MetricSpec",
    "ProtrusionEvent",
    "RegisteredMetric",
    "SingleCellState",
    "StepDiagnostics",
    "V2DataContract",
    "accumulate_prescribed_traction",
    "boundary_perimeter_um",
    "boundary_projected_area_um2",
    "compute_area_forces",
    "compute_cortex_forces",
    "compute_finite_n_residual",
    "compute_rate_max",
    "default_registry",
    "ellipse_polygon_vertices",
    "energy_area",
    "energy_cortex",
    "equilibrium_radius_from_cubic",
    "equilibrium_radius_from_quartic",
    "read_frame",
    "regular_polygon_vertices",
    "run_active_contour_test",
    "step",
    "render_frame_html",
    "render_frame_png",
    "write_frame",
]
