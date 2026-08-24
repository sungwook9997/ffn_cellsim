"""Pure summary gates for the post-loop WCA contact-graph diagnostic."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.incumbent.contact_graph_diagnostics import (
    ContactGraphDiagnostics,
    summarize_wca_contact_graph,
)


def test_contact_graph_summary_counts_symmetric_pairs_once_and_reports_hotspot() -> None:
    diagnostics = ContactGraphDiagnostics(
        degree=np.array([2, 1, 1, 0], dtype=np.int32),
        capped_degree=np.array([1, 1, 0, 0], dtype=np.int32),
        nearest_distance_um=np.array([0.005, 0.006, 0.007, np.inf]),
        max_pair_force_pn=np.array([1000.0, 1000.0, 20.0, 0.0]),
        radial_curvature_sum_pn_per_um=np.array([50.0, 0.0, 25.0, 0.0]),
        nearest_neighbor=np.array([1, 0, 0, -1], dtype=np.int32),
        strongest_neighbor=np.array([1, 0, 0, -1], dtype=np.int32),
    )
    summary = summarize_wca_contact_graph(diagnostics, sigma_um=0.01, top_node_index=0)

    assert summary["undirected_contact_count"] == 2
    assert summary["undirected_capped_contact_count"] == 1
    assert summary["contact_node_count"] == 3
    assert summary["nearest_distance_over_sigma_quantiles"][-1] == pytest.approx(0.7)
    assert summary["projected_residual_top_node"] == {
        "index": 0,
        "contact_degree": 2,
        "capped_contact_degree": 1,
        "nearest_distance_over_sigma": 0.5,
        "max_pair_force_pN": 1000.0,
        "radial_curvature_sum_pN_per_um": 50.0,
        "nearest_neighbor": 1,
        "strongest_neighbor": 1,
    }


def test_contact_graph_summary_rejects_asymmetric_directed_counts() -> None:
    diagnostics = ContactGraphDiagnostics(
        degree=np.array([1, 0], dtype=np.int32),
        capped_degree=np.zeros(2, dtype=np.int32),
        nearest_distance_um=np.array([0.005, np.inf]),
        max_pair_force_pn=np.array([1.0, 0.0]),
        radial_curvature_sum_pn_per_um=np.array([1.0, 0.0]),
        nearest_neighbor=np.array([1, -1], dtype=np.int32),
        strongest_neighbor=np.array([1, -1], dtype=np.int32),
    )
    with pytest.raises(RuntimeError, match="odd directed contact count"):
        summarize_wca_contact_graph(diagnostics, sigma_um=0.01, top_node_index=0)
