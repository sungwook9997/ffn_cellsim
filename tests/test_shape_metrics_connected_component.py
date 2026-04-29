from __future__ import annotations

import numpy as np

from acs.analysis.shape_metrics import (
    top_down_connected_component_area,
    top_down_projection_area,
)


def test_connected_component_area_detects_single_outlier_hull_inflation():
    theta = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
    core = np.column_stack([
        np.cos(theta),
        np.sin(theta),
        np.zeros_like(theta),
    ])
    outlier = np.array([[8.0, 0.0, 0.0]])
    x = np.vstack([core, outlier])

    all_hull = top_down_projection_area(x)
    core_metric = top_down_connected_component_area(x, link_radius=0.6)

    assert core_metric["n_component"] == 32
    assert core_metric["fraction"] == 32 / 33
    assert all_hull / core_metric["area"] > 2.0


def test_connected_component_area_tracks_single_connected_sheet():
    grid_x, grid_y = np.meshgrid(np.linspace(-1.0, 1.0, 8), np.linspace(-1.0, 1.0, 8))
    x = np.column_stack([
        grid_x.ravel(),
        grid_y.ravel(),
        np.zeros(grid_x.size),
    ])

    metric = top_down_connected_component_area(x, link_radius=0.45)

    assert metric["n_component"] == x.shape[0]
    assert metric["fraction"] == 1.0
    assert np.isfinite(metric["area"])
