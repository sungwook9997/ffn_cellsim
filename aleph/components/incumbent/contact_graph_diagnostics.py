"""Post-solve CUDA diagnostics for the live WCA contact graph.

The mechanical solver deliberately keeps the physical-time loop device resident.  This module is therefore
diagnostic-only: it is called after a candidate solve has ended and reports the topology that the next
contact-cluster preconditioner must resolve.  It does not alter positions, forces, kinetics, or acceptance.

Sanity Gate:
    * Dimensions: distance is um, pair force is pN, and radial curvature is pN/um.
    * Boundary cases: inactive and contact-free nodes report zero degree, zero force/curvature, and an infinite
      nearest-contact distance; same-filament pairs are excluded exactly as in the production WCA kernel.
    * Conservation/counting: the own-row traversal reports directed degree.  Every admissible pair is visited
      twice, so the host summary divides the summed degree by two and requires an even directed count.
    * Numerical: the force-cap predicate and WCA formula match ``wca_steric_kernel``; capped contacts report
      zero radial derivative, matching the analytic implicit operator.
    * Residency: neighbor discovery and per-node statistics run on Warp CUDA.  Host reads occur only after the
      quasi-static candidate and physical-time hot loop have ended.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

_TWO_POW_1_6 = wp.constant(wp.float64(2.0 ** (1.0 / 6.0)))
_INFINITY = wp.constant(wp.float64(float("inf")))


@wp.kernel
def wca_contact_graph_stats_kernel(
    grid: wp.uint64,
    query_points: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    fiber_id: wp.array(dtype=wp.int32),
    active_nodes: wp.array(dtype=wp.int32),
    radius: wp.float32,
    sigma: wp.float64,
    epsilon: wp.float64,
    force_cap: wp.float64,
    degree: wp.array(dtype=wp.int32),
    capped_degree: wp.array(dtype=wp.int32),
    nearest_distance: wp.array(dtype=wp.float64),
    max_pair_force: wp.array(dtype=wp.float64),
    radial_curvature_sum: wp.array(dtype=wp.float64),
    nearest_neighbor: wp.array(dtype=wp.int32),
    strongest_neighbor: wp.array(dtype=wp.int32),
) -> None:
    """Measure each live node's WCA neighbors with the production pair predicates."""
    i = wp.tid()
    zero = wp.float64(0.0)
    degree[i] = wp.int32(0)
    capped_degree[i] = wp.int32(0)
    nearest_distance[i] = _INFINITY
    max_pair_force[i] = zero
    radial_curvature_sum[i] = zero
    nearest_neighbor[i] = wp.int32(-1)
    strongest_neighbor[i] = wp.int32(-1)
    if active_nodes[i] < wp.int32(0):
        return

    pi = pos[i]
    fi = fiber_id[i]
    cutoff = _TWO_POW_1_6 * sigma
    query = wp.hash_grid_query(grid, query_points[i], radius)
    j = wp.int32(0)
    while wp.hash_grid_query_next(query, j):
        if j != i and active_nodes[j] >= wp.int32(0) and fiber_id[j] != fi:
            delta = pi - pos[j]
            length = wp.length(delta)
            if length > wp.float64(1.0e-12) and length < cutoff:
                ratio = sigma / length
                sr6 = ratio * ratio * ratio
                sr6 = sr6 * sr6
                force_magnitude = (wp.float64(24.0) * epsilon / length) * (
                    wp.float64(2.0) * sr6 * sr6 - sr6)
                degree[i] = degree[i] + wp.int32(1)
                if length < nearest_distance[i]:
                    nearest_distance[i] = length
                    nearest_neighbor[i] = j
                reported_force = force_magnitude
                if force_cap > zero and force_magnitude > force_cap:
                    reported_force = force_cap
                if reported_force > max_pair_force[i]:
                    max_pair_force[i] = reported_force
                    strongest_neighbor[i] = j
                if force_cap > zero and force_magnitude > force_cap:
                    capped_degree[i] = capped_degree[i] + wp.int32(1)
                else:
                    radial_curvature_sum[i] = radial_curvature_sum[i] + (
                        wp.float64(24.0) * epsilon / (length * length)
                    ) * (wp.float64(26.0) * sr6 * sr6 - wp.float64(7.0) * sr6)


@dataclass(frozen=True)
class ContactGraphDiagnostics:
    """Host-side post-loop copy of per-node contact-graph statistics."""

    degree: np.ndarray
    capped_degree: np.ndarray
    nearest_distance_um: np.ndarray
    max_pair_force_pn: np.ndarray
    radial_curvature_sum_pn_per_um: np.ndarray
    nearest_neighbor: np.ndarray
    strongest_neighbor: np.ndarray


def measure_wca_contact_graph(steric: object, pos: wp.array) -> ContactGraphDiagnostics:
    """Measure the current production WCA graph on CUDA and copy it after synchronization."""
    n = int(steric.n)
    device = steric.device
    degree_d = wp.zeros(n, dtype=wp.int32, device=device)
    capped_degree_d = wp.zeros(n, dtype=wp.int32, device=device)
    nearest_d = wp.empty(n, dtype=wp.float64, device=device)
    max_force_d = wp.zeros(n, dtype=wp.float64, device=device)
    radial_d = wp.zeros(n, dtype=wp.float64, device=device)
    nearest_neighbor_d = wp.empty(n, dtype=wp.int32, device=device)
    strongest_neighbor_d = wp.empty(n, dtype=wp.int32, device=device)

    # Rebuild the same live grid at the final diagnostic position.  This is outside the solve loop.
    steric._accumulate_pos(pos, wp.zeros(n, dtype=wp.vec3d, device=device))
    wp.launch(
        wca_contact_graph_stats_kernel,
        dim=n,
        inputs=[
            steric.grid.id,
            steric._qpts,
            pos,
            steric.fiber_id,
            steric.active,
            wp.float32(steric.r_c),
            wp.float64(steric.sigma),
            wp.float64(steric.epsilon),
            wp.float64(steric.f_cap),
            degree_d,
            capped_degree_d,
            nearest_d,
            max_force_d,
            radial_d,
            nearest_neighbor_d,
            strongest_neighbor_d,
        ],
        device=device,
    )
    wp.synchronize_device(device)
    return ContactGraphDiagnostics(
        degree=degree_d.numpy(),
        capped_degree=capped_degree_d.numpy(),
        nearest_distance_um=nearest_d.numpy(),
        max_pair_force_pn=max_force_d.numpy(),
        radial_curvature_sum_pn_per_um=radial_d.numpy(),
        nearest_neighbor=nearest_neighbor_d.numpy(),
        strongest_neighbor=strongest_neighbor_d.numpy(),
    )


def summarize_wca_contact_graph(
    diagnostics: ContactGraphDiagnostics,
    *,
    sigma_um: float,
    top_node_index: int,
) -> dict[str, object]:
    """Return exact pair counts, distribution summaries, and the residual hotspot neighborhood."""
    degree = np.asarray(diagnostics.degree, dtype=np.int64)
    capped = np.asarray(diagnostics.capped_degree, dtype=np.int64)
    directed = int(degree.sum())
    directed_capped = int(capped.sum())
    if directed % 2 or directed_capped % 2:
        raise RuntimeError("symmetric WCA graph produced an odd directed contact count")
    contact_nodes = degree > 0
    finite_distance = np.isfinite(diagnostics.nearest_distance_um)

    def _quantiles(values: np.ndarray) -> list[float]:
        if values.size == 0:
            return []
        return np.quantile(values, [0.0, 0.5, 0.9, 0.99, 1.0]).tolist()

    top = int(top_node_index)
    top_valid = 0 <= top < degree.size
    return {
        "undirected_contact_count": directed // 2,
        "undirected_capped_contact_count": directed_capped // 2,
        "contact_node_count": int(np.count_nonzero(contact_nodes)),
        "degree_quantiles": _quantiles(degree[contact_nodes]),
        "capped_degree_quantiles": _quantiles(capped[contact_nodes]),
        "nearest_distance_over_sigma_quantiles": _quantiles(
            diagnostics.nearest_distance_um[finite_distance] / float(sigma_um)),
        "max_pair_force_pN_quantiles": _quantiles(
            diagnostics.max_pair_force_pn[contact_nodes]),
        "radial_curvature_sum_pN_per_um_quantiles": _quantiles(
            diagnostics.radial_curvature_sum_pn_per_um[contact_nodes]),
        "projected_residual_top_node": {
            "index": top,
            "contact_degree": int(degree[top]) if top_valid else None,
            "capped_contact_degree": int(capped[top]) if top_valid else None,
            "nearest_distance_over_sigma": (
                float(diagnostics.nearest_distance_um[top] / sigma_um)
                if top_valid and np.isfinite(diagnostics.nearest_distance_um[top]) else None
            ),
            "max_pair_force_pN": (
                float(diagnostics.max_pair_force_pn[top]) if top_valid else None
            ),
            "radial_curvature_sum_pN_per_um": (
                float(diagnostics.radial_curvature_sum_pn_per_um[top]) if top_valid else None
            ),
            "nearest_neighbor": int(diagnostics.nearest_neighbor[top]) if top_valid else None,
            "strongest_neighbor": int(diagnostics.strongest_neighbor[top]) if top_valid else None,
        },
    }
