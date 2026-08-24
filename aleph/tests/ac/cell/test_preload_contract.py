"""Pure-host sanity gates for the membrane--ERM--cortex preload capacity contract."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.compartments import density_resolved_erm_pairs
from aleph.components.incumbent.preload_contract import evaluate_preload_capacity, triangle_surface_area


def test_current_continuum_force_substitution_is_diagnostic_arithmetic_only() -> None:
    capacity = evaluate_preload_capacity(
        pressure_pa=40.0,
        radius_um=7.5,
        membrane_area_um2=703.4902162858084,
        membrane_tension_pn_per_um=10.0,
        erm_tether_count=642,
        erm_rupture_force_pn=11.43470677829702,
    )
    assert capacity.required_erm_count_min == 2297
    assert capacity.required_erm_density_min_per_um2 == pytest.approx(3.2649139201533086)
    assert capacity.total_pressure_capacity_upper_bound_pa == pytest.approx(13.101896015590846)
    assert capacity.pressure_capacity_ratio == pytest.approx(0.32754740038977115)
    assert not capacity.mechanically_feasible_upper_bound


def test_membrane_tension_alone_closes_zero_residual_boundary() -> None:
    capacity = evaluate_preload_capacity(
        pressure_pa=2.0,
        radius_um=10.0,
        membrane_area_um2=100.0,
        membrane_tension_pn_per_um=10.0,
        erm_tether_count=0,
        erm_rupture_force_pn=5.0,
    )
    assert capacity.membrane_laplace_pressure_pa == 2.0
    assert capacity.residual_pressure_for_erm_pa == 0.0
    assert capacity.required_erm_count_min == 0
    assert capacity.mechanically_feasible_upper_bound


def test_required_density_is_grid_invariant_and_count_scales_with_area() -> None:
    coarse = evaluate_preload_capacity(
        pressure_pa=30.0, radius_um=10.0, membrane_area_um2=100.0,
        membrane_tension_pn_per_um=0.0, erm_tether_count=400, erm_rupture_force_pn=10.0,
    )
    fine = evaluate_preload_capacity(
        pressure_pa=30.0, radius_um=10.0, membrane_area_um2=400.0,
        membrane_tension_pn_per_um=0.0, erm_tether_count=1600, erm_rupture_force_pn=10.0,
    )
    assert coarse.required_erm_density_min_per_um2 == fine.required_erm_density_min_per_um2
    assert fine.required_erm_count_min == 4 * coarse.required_erm_count_min
    assert fine.pressure_capacity_ratio == coarse.pressure_capacity_ratio


def test_triangle_area_and_invalid_inputs() -> None:
    verts = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0]])
    assert triangle_surface_area(verts, np.array([[0, 1, 2]])) == 3.0
    with pytest.raises(ValueError, match="positive"):
        evaluate_preload_capacity(
            pressure_pa=1.0, radius_um=0.0, membrane_area_um2=1.0,
            membrane_tension_pn_per_um=0.0, erm_tether_count=0, erm_rupture_force_pn=1.0,
        )


def test_density_resolved_pairs_have_exact_count_and_unique_cortex_endpoints() -> None:
    membrane = np.array([
        [0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
    ])
    cortex = np.repeat(0.9 * membrane, 2, axis=0)
    mem_idx, cortex_idx = density_resolved_erm_pairs(
        membrane, cortex, surface_area_um2=4.0, density_per_um2=1.0, reach_um=0.2)
    assert mem_idx.shape == cortex_idx.shape == (4,)
    assert np.unique(cortex_idx).size == 4
    assert np.all(np.linalg.norm(membrane[mem_idx] - cortex[cortex_idx], axis=1) <= 0.2)


def test_density_builder_refuses_more_linkers_than_mechanistic_cortex_nodes() -> None:
    membrane = np.array([[0.0, 0.0, 1.0]])
    cortex = np.array([[0.0, 0.0, 0.9]])
    with pytest.raises(ValueError, match="increase mechanistic cortex resolution"):
        density_resolved_erm_pairs(
            membrane, cortex, surface_area_um2=1.0, density_per_um2=2.0, reach_um=0.2)


def test_mcf7_production_provenance_cannot_exist_without_a_density() -> None:
    with pytest.raises(ValueError, match="MCF7 production ERM provenance"):
        build_cell(CellConfig(erm_density_mcf7_production=True))
