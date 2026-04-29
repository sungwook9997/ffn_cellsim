"""acs.visualization — plotting + figure-generation utilities for the
ActiveCellSim project's layer-by-layer narrative.

Functions:
- `aggregate_narrative_table()`: scan results/ directory for all stage
  pilot artifacts; return tidy DataFrame for plotting.
- `plot_layer_by_layer_R_drift()`: 8+ layer-by-layer R drift bar chart.
- `plot_phenotype_comparison()`: 3-phenotype Bare/Pre/Lam4 R drift +
  A/A₀_topdown comparison (Stage 1b vs Phase 4-v2).
- `plot_AA0_trajectory()`: A/A₀_topdown(t) trajectories overlaid (multi-
  stage + production when available).
- `plot_zeta_response_curve()`: 13-point ζ instability curve.
- `plot_3d_snapshots()`: per-stage final-frame 3D particle scatter
  (matplotlib for headless SSH compatibility).
- `generate_summary()`: writes results/visualization/SUMMARY.md with
  inline figure references + commentary.

PI Visualization directive 2026-04-29 (parallel with Production Lam4).
"""

from acs.visualization.plots import (
    aggregate_narrative_table,
    generate_summary,
    plot_3d_snapshots,
    plot_AA0_trajectory,
    plot_layer_by_layer_R_drift,
    plot_layer_x_phenotype_summary,
    plot_phenotype_comparison,
    plot_pi_comparison,
    plot_production_bucket_classification,
    plot_production_trajectory,
    plot_zeta_response_curve,
)

__all__ = [
    "aggregate_narrative_table",
    "generate_summary",
    "plot_3d_snapshots",
    "plot_AA0_trajectory",
    "plot_layer_by_layer_R_drift",
    "plot_layer_x_phenotype_summary",
    "plot_phenotype_comparison",
    "plot_pi_comparison",
    "plot_production_bucket_classification",
    "plot_production_trajectory",
    "plot_zeta_response_curve",
]
