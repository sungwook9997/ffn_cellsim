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
from acs.v2.dynamics.closed_loop_phase_d import (
    FAToECMResponseResult,
    PhaseDNoOpStepResult,
    step_ecm_to_fa_bias,
    step_fa_to_ecm_response,
    step_phase_d_no_op,
)
from acs.v2.dynamics.closed_loop_phase_e import (
    PhaseEStepResult,
    step_closed_loop_phase_e_v1,
)
from acs.v2.dynamics.closed_loop_phase_e_sweep import (
    PhaseEV1Item5SweepConfig,
    PhaseEV1Item5SweepMetadata,
    PhaseEV1Item5SweepResult,
    run_phase_e_v1_sensitivity_sweep,
)
from acs.v2.dynamics.ecm_constitutive_response import (
    ECMConstitutiveResponseError,
    ECMOrientationResponseDiagnostics,
    ECMOrientationResponseResult,
    K_ORIENT_PER_S,
    TAU_ALIGN_RANGE_S,
    TRACTION_REF_NN_PER_UM2,
    step_ecm_orientation_response,
)
from acs.v2.dynamics.ecm_lyapunov_metric import (
    ECMOrientationLyapunovMetricDiagnostics,
    ECMOrientationLyapunovMetricResult,
    MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL,
    MAX_SQ_FROBENIUS_DIFF_PER_CELL,
    compute_ecm_orientation_lyapunov_metric,
)
from acs.v2.dynamics.ecm_open_loop import (
    ECMOpenLoopError,
    accumulate_prescribed_traction,
    apply_prescribed_density_rate,
    apply_prescribed_orientation_rate,
    apply_prescribed_stiffness_rate,
)
from acs.v2.dynamics.ecm_to_fa_bias import (
    ECMSampledAtFAs,
    ECMToFABiasResult,
    FAToECMBiasError,
    RATE_NAMES,
    compute_ecm_to_fa_bias_active,
    compute_ecm_to_fa_bias_neutral,
    sample_ecm_at_fa_positions,
)
from acs.v2.dynamics.fa_to_ecm_scattering import (
    FAToECMScatteringError,
    scatter_fa_traction_to_ecm_bilinear,
)
from acs.v2.dynamics.focal_adhesion import (
    FocalAdhesionDynamicsDiagnostics,
    FocalAdhesionDynamicsError,
    FocalAdhesionDynamicsParameters,
    FocalAdhesionDynamicsResult,
    compute_radial_tangential_decomposition,
    step_focal_adhesions_static,
)
from acs.v2.dynamics.protrusion_coupled_focal_adhesion import (
    ProtrusionCoupledDynamicsDiagnostics,
    ProtrusionCoupledDynamicsResult,
    ProtrusionStateMultipliers,
    step_protrusion_coupled_focal_adhesions,
)
from acs.v2.ecm_open_loop_harness import (
    EcmOlRun,
    EcmOlScenario,
    EcmOlStepDiagnostics,
    make_default_ecm,
    run_ecm_ol_scenario,
)
from acs.v2.cell_cluster import CellClusterState
from acs.v2.data_contract import (
    ArtifactLayoutEntry,
    ArtifactKind,
    ImagingDatasetSpec,
    MetricSpec,
    SegmentationProvenance,
    V2DataContract,
    canonical_artifact_layout,
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
    csv_scalar_metric,
    csv_spheroid_a_over_a0,
    csv_spheroid_area_um2,
    csv_spheroid_effective_radius_um,
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
    "ArtifactLayoutEntry",
    "ArtifactKind",
    "CellClusterState",
    "ECMConstitutiveResponseError",
    "ECMOpenLoopError",
    "ECMOrientationLyapunovMetricDiagnostics",
    "ECMOrientationLyapunovMetricResult",
    "ECMOrientationResponseDiagnostics",
    "ECMOrientationResponseResult",
    "ECMSampledAtFAs",
    "ECMSubstrateState",
    "ECMToFABiasResult",
    "EcmOlRun",
    "EcmOlScenario",
    "EcmOlStepDiagnostics",
    "FocalAdhesionDynamicsError",
    "FocalAdhesionDynamicsParameters",
    "FAToECMBiasError",
    "FAToECMResponseResult",
    "FAToECMScatteringError",
    "FocalAdhesionDynamicsDiagnostics",
    "FocalAdhesionDynamicsResult",
    "FocalAdhesionState",
    "FrameDumpReadResult",
    "HarnessRun",
    "ImagingDatasetSpec",
    "JunctionState",
    "K_ORIENT_PER_S",
    "MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL",
    "MAX_SQ_FROBENIUS_DIFF_PER_CELL",
    "MeasurementBoundary",
    "MeasurementBoundaryError",
    "MetricRegistry",
    "MetricRegistryError",
    "MetricSpec",
    "PhaseDNoOpStepResult",
    "PhaseEStepResult",
    "PhaseEV1Item5SweepConfig",
    "PhaseEV1Item5SweepMetadata",
    "PhaseEV1Item5SweepResult",
    "ProtrusionCoupledDynamicsDiagnostics",
    "ProtrusionCoupledDynamicsResult",
    "ProtrusionEvent",
    "ProtrusionStateMultipliers",
    "RATE_NAMES",
    "RegisteredMetric",
    "SegmentationProvenance",
    "SingleCellState",
    "StepDiagnostics",
    "TAU_ALIGN_RANGE_S",
    "TRACTION_REF_NN_PER_UM2",
    "V2DataContract",
    "accumulate_prescribed_traction",
    "apply_prescribed_density_rate",
    "apply_prescribed_orientation_rate",
    "apply_prescribed_stiffness_rate",
    "boundary_perimeter_um",
    "boundary_projected_area_um2",
    "canonical_artifact_layout",
    "compute_area_forces",
    "compute_cortex_forces",
    "compute_ecm_orientation_lyapunov_metric",
    "compute_ecm_to_fa_bias_active",
    "compute_ecm_to_fa_bias_neutral",
    "compute_finite_n_residual",
    "compute_radial_tangential_decomposition",
    "compute_rate_max",
    "csv_scalar_metric",
    "csv_spheroid_a_over_a0",
    "csv_spheroid_area_um2",
    "csv_spheroid_effective_radius_um",
    "default_registry",
    "ellipse_polygon_vertices",
    "energy_area",
    "energy_cortex",
    "equilibrium_radius_from_cubic",
    "make_default_ecm",
    "run_ecm_ol_scenario",
    "equilibrium_radius_from_quartic",
    "read_frame",
    "regular_polygon_vertices",
    "run_active_contour_test",
    "run_phase_e_v1_sensitivity_sweep",
    "step",
    "step_closed_loop_phase_e_v1",
    "step_ecm_orientation_response",
    "step_ecm_to_fa_bias",
    "step_fa_to_ecm_response",
    "step_focal_adhesions_static",
    "step_phase_d_no_op",
    "step_protrusion_coupled_focal_adhesions",
    "render_frame_html",
    "render_frame_png",
    "sample_ecm_at_fa_positions",
    "scatter_fa_traction_to_ecm_bilinear",
    "write_frame",
]
