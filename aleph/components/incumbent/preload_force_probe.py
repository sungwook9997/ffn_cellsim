#!/usr/bin/env python3
"""Diagnose a rejected physiological preload candidate by CUDA force-family decomposition.

This is a gate-side diagnostic, not a simulation runtime.  It runs one native quasi-static candidate on CUDA,
captures that candidate with a device-to-device copy immediately before the transactional rollback, and performs
all host reads only after the inner loop has ended.  No force family is removed or rescaled.

Sanity Gate:
    * Dimensions: every reported force is pN; ERM extension is um and ``k_erm * extension`` is pN.
    * Boundary cases: an absent family reports zero; unbound ERM states contribute neither extension nor force.
    * Conservation: each family reports its vector sum as well as its maximum nodal magnitude.
    * Numerical: the candidate is the exact rejected device state, not a host reconstruction.
    * Sign sense: positive ERM extension denotes membrane-cortex separation beyond the formation length.
    * Measurement protocol: synchronization and NumPy reads happen after the physical/inner hot loop only.
"""

from __future__ import annotations

import argparse
import json
import time
from collections.abc import Callable
from pathlib import Path

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.contact_graph_diagnostics import (
    measure_wca_contact_graph,
    summarize_wca_contact_graph,
)
from aleph.components.incumbent.driver import make_inner_solve
from aleph.components.incumbent.erm_tether import erm_tether_force_kernel
from aleph.world.gpu_memory import query_process_peak
from aleph.components.incumbent.inner_mechanics import project_constraint_forces_kernel
from aleph.laws.forces_warp import cytosim_bending_kernel
from aleph.laws.membrane_surface import helfrich_bending_kernel, membrane_area_kernel
from aleph.laws.network_warp import _zero, link_spring_kernel


def _node_family(cell: object, index: int) -> str:
    """Name the contiguous population block containing a global node index."""
    if index < cell.n_actin:
        return "actin"
    if cell.nucleus is not None and index < cell.nucleus.node_off:
        return "myosin"
    if cell.nucleus is not None and index < cell.nucleus.node_off + cell.nucleus.n_verts:
        return "nucleus"
    if cell.membrane is not None and index < cell.membrane.node_off + cell.membrane.n_verts:
        return "membrane"
    return "other"


def _force_stats(cell: object, launch: Callable[[], None]) -> dict[str, object]:
    """Launch one additive family into the shared force buffer and summarize after synchronization."""
    wp.launch(_zero, dim=cell.n_total, inputs=[cell.f_d], device=cell.device)
    launch()
    wp.synchronize_device(cell.device)
    force = cell.f_d.numpy()
    magnitude = np.linalg.norm(force, axis=1)
    top = int(np.argmax(magnitude)) if magnitude.size else -1
    return {
        "max_nodal_force_pN": float(magnitude[top]) if top >= 0 else 0.0,
        "rms_nodal_force_pN": float(np.sqrt(np.mean(magnitude * magnitude))) if magnitude.size else 0.0,
        "net_force_pN": np.sum(force, axis=0).tolist(),
        "top_node_index": top,
        "top_node_family": _node_family(cell, top) if top >= 0 else "none",
    }


def _array_stats(cell: object, force: np.ndarray) -> dict[str, object]:
    """Summarize a host diagnostic copy of a device force array."""
    magnitude = np.linalg.norm(force, axis=1)
    top = int(np.argmax(magnitude)) if magnitude.size else -1
    return {
        "max_nodal_force_pN": float(magnitude[top]) if top >= 0 else 0.0,
        "rms_nodal_force_pN": float(np.sqrt(np.mean(magnitude * magnitude))) if magnitude.size else 0.0,
        "net_force_pN": np.sum(force, axis=0).tolist(),
        "top_node_index": top,
        "top_node_family": _node_family(cell, top) if top >= 0 else "none",
    }


def run_probe(args: argparse.Namespace) -> dict[str, object]:
    """Run and summarize one full-native preload candidate."""
    wp.init()
    build_t0 = time.perf_counter()
    cell = build_cell(CellConfig(
        n_filaments=args.n_filaments,
        with_myosin=True,
        with_steric=not args.no_steric,
        with_nucleus=True,
        with_membrane=True,
        with_pressure=True,
        erm_density_per_um2=args.erm_density,
        erm_density_source=args.erm_density_source,
        device=args.device,
    ))
    wp.synchronize_device(cell.device)
    build_wall_s = time.perf_counter() - build_t0
    setup_t0 = time.perf_counter()
    inner = make_inner_solve(
        cell,
        n_inner=args.inner_iters,
        reshape_every=args.reshape_every,
        max_inner_retries=args.max_inner_retries,
        capture_candidate=True,
        inner_solver=args.inner_solver,
        implicit_cg_max_iterations=args.implicit_cg_max_iterations,
        implicit_line_search_steps=args.implicit_line_search_steps,
        implicit_coarse_iterations=args.implicit_coarse_iterations,
        implicit_coarse_modes=args.implicit_coarse_modes,
        rkc_stages=args.rkc_stages,
    )
    wp.synchronize_device(cell.device)
    accelerator_setup_wall_s = time.perf_counter() - setup_t0
    solve_t0 = time.perf_counter()
    report = inner(args.dt_phys)
    wp.synchronize_device(cell.device)
    solve_wall_s = time.perf_counter() - solve_t0

    candidate = inner.candidate_pos_d
    if candidate is None:
        raise RuntimeError("diagnostic candidate capture was not allocated")
    wp.copy(cell.pos_d, candidate)
    wp.synchronize_device(cell.device)

    families: dict[str, dict[str, object]] = {}
    families["actin_bending"] = _force_stats(
        cell,
        lambda: wp.launch(
            cytosim_bending_kernel,
            dim=cell.n_tri,
            inputs=[cell.pos_d, cell.tri_d, cell.alpha_d, cell.f_d],
            device=cell.device,
        ),
    )
    families["actin_crosslinks"] = _force_stats(
        cell,
        lambda: wp.launch(
            link_spring_kernel,
            dim=cell.n_xl,
            inputs=[cell.pos_d, cell.xl_d, cell.kxl_d, cell.r0xl_d, cell.f_d],
            device=cell.device,
        ),
    )
    families["myosin"] = _force_stats(
        cell,
        lambda: cell.myosin.accumulate(cell.pos_d, cell.f_d),
    )
    families["nucleus"] = _force_stats(
        cell,
        lambda: cell.nucleus.accumulate(cell.pos_d, cell.f_d),
    )
    membrane = cell.membrane
    families["membrane_bending"] = _force_stats(
        cell,
        lambda: wp.launch(
            helfrich_bending_kernel,
            dim=membrane.n_hinges,
            inputs=[cell.pos_d, membrane.hinges_d, wp.float64(membrane.kappa_tilde)],
            outputs=[cell.f_d],
            device=cell.device,
        ),
    )
    families["membrane_area_tension"] = _force_stats(
        cell,
        lambda: wp.launch(
            membrane_area_kernel,
            dim=membrane.n_faces,
            inputs=[cell.pos_d, membrane.faces_d, wp.float64(membrane.gamma_mem)],
            outputs=[cell.f_d],
            device=cell.device,
        ),
    )
    families["erm_tethers"] = _force_stats(
        cell,
        lambda: wp.launch(
            erm_tether_force_kernel,
            dim=membrane.n_erm,
            inputs=[
                cell.pos_d,
                membrane.erm_m_d,
                membrane.erm_c_d,
                membrane.erm_bound_d,
                wp.float64(membrane.k_erm),
                membrane.erm_rest_d,
                wp.float64(membrane.f_rupt),
            ],
            outputs=[cell.f_d],
            device=cell.device,
        ),
    )
    families["membrane_pressure"] = _force_stats(
        cell,
        lambda: cell.membrane_pressure.accumulate(cell.pos_d, cell.f_d),
    )
    families["steric"] = _force_stats(
        cell,
        lambda: cell.steric.accumulate(cell.state, cell.f_d) if cell.steric is not None else None,
    )
    families["bulk_pressure"] = _force_stats(
        cell,
        lambda: cell.pressure.accumulate(cell.state, cell.f_d),
    )

    # Reassemble the exact total after the family probes altered the shared force buffer.
    from aleph.components.incumbent.driver import _accumulate_all

    families["total"] = _force_stats(cell, lambda: _accumulate_all(cell, cell.pos_d, cell.f_d))

    # Apply the same exact NF2007 tangent-space projector as the runtime to the final candidate force. This
    # separates large raw constraint reactions from the residual that can still move an inextensible fiber.
    projected_force_d = wp.empty_like(cell.f_d)
    wp.copy(projected_force_d, cell.f_d)
    projector_finite_d = wp.ones(1, dtype=wp.int32, device=cell.device)
    wp.launch(
        project_constraint_forces_kernel,
        dim=cell.n_fibers,
        inputs=[
            cell.pos_d,
            cell.f_d,
            cell.foff_d,
            cell.soff_d,
            projected_force_d,
            wp.zeros(int(cell.srest_d.shape[0]), dtype=wp.float64, device=cell.device),
            wp.zeros(int(cell.srest_d.shape[0]), dtype=wp.float64, device=cell.device),
            projector_finite_d,
        ],
        device=cell.device,
    )
    wp.synchronize_device(cell.device)
    projected_force = projected_force_d.numpy()
    projected_stats = _array_stats(cell, projected_force)
    projected_stats["finite"] = bool(projector_finite_d.numpy()[0])
    fiber_offset = cell.foff_d.numpy().astype(np.int64)
    fiber_force = np.add.reduceat(projected_force[:cell.n_actin], fiber_offset[:-1], axis=0)
    fiber_force_magnitude = np.linalg.norm(fiber_force, axis=1)
    top_fiber = int(np.argmax(fiber_force_magnitude))
    projected_top_node = int(projected_stats["top_node_index"])
    top_node_fiber = (
        int(np.searchsorted(fiber_offset, projected_top_node, side="right") - 1)
        if 0 <= projected_top_node < cell.n_actin else None
    )
    projected_stats["fiber_translation_residual"] = {
        "max_fiber_sum_force_pN": float(fiber_force_magnitude[top_fiber]),
        "rms_fiber_sum_force_pN": float(np.sqrt(np.mean(fiber_force_magnitude**2))),
        "top_fiber_index": top_fiber,
        "top_node_fiber_index": top_node_fiber,
        "top_node_fiber_sum_force_pN": (
            fiber_force[top_node_fiber].tolist() if top_node_fiber is not None else None
        ),
        "top_node_fiber_sum_force_magnitude_pN": (
            float(fiber_force_magnitude[top_node_fiber]) if top_node_fiber is not None else None
        ),
    }

    contact_graph = None
    if cell.steric is not None:
        contact_diagnostics = measure_wca_contact_graph(cell.steric, cell.pos_d)
        contact_graph = summarize_wca_contact_graph(
            contact_diagnostics,
            sigma_um=float(cell.steric.sigma),
            top_node_index=int(projected_stats["top_node_index"]),
        )
        xlink_degree = np.zeros(cell.n_total, dtype=np.int64)
        crosslinks = cell.xl_d.numpy().astype(np.int64)
        if crosslinks.size:
            np.add.at(xlink_degree, crosslinks[:, 0], 1)
            np.add.at(xlink_degree, crosslinks[:, 1], 1)
        top_node = int(projected_stats["top_node_index"])
        contact_graph["projected_residual_top_node"]["crosslink_degree"] = (
            int(xlink_degree[top_node]) if 0 <= top_node < xlink_degree.size else None
        )
        if 0 <= top_node < xlink_degree.size:
            incident = crosslinks[(crosslinks[:, 0] == top_node) | (crosslinks[:, 1] == top_node)]
            xlink_neighbors = np.where(incident[:, 0] == top_node, incident[:, 1], incident[:, 0])
            contact_graph["projected_residual_top_node"]["crosslink_neighbors"] = (
                xlink_neighbors.astype(int).tolist()
            )
            contact_graph["projected_residual_top_node"]["strongest_contact_is_crosslinked"] = bool(
                int(contact_graph["projected_residual_top_node"]["strongest_neighbor"])
                in set(xlink_neighbors.tolist())
            )

    pos = cell.pos_d.numpy()
    mem_idx = membrane.erm_m_d.numpy().astype(np.int64)
    cortex_idx = membrane.erm_c_d.numpy().astype(np.int64)
    rest = membrane.erm_rest_d.numpy()
    bound = membrane.erm_bound_d.numpy().astype(bool)
    extension = np.linalg.norm(pos[mem_idx] - pos[cortex_idx], axis=1) - rest
    tension = membrane.k_erm * np.maximum(extension, 0.0)
    active_tension = tension[bound]
    mem_pos = pos[membrane.node_off:membrane.node_off + membrane.n_verts]
    mem_radius = np.linalg.norm(mem_pos, axis=1)
    cortex_radius = np.linalg.norm(pos[:cell.n_actin], axis=1)

    preprojection = inner.candidate_preprojection_pos_d
    if preprojection is None:
        raise RuntimeError("diagnostic preprojection candidate capture was not allocated")
    pre_pos = preprojection.numpy()
    projection_delta = np.linalg.norm(pos - pre_pos, axis=1)
    wp.copy(cell.pos_d, preprojection)
    preprojection_families = {
        "actin_crosslinks": _force_stats(
            cell,
            lambda: wp.launch(
                link_spring_kernel,
                dim=cell.n_xl,
                inputs=[cell.pos_d, cell.xl_d, cell.kxl_d, cell.r0xl_d, cell.f_d],
                device=cell.device,
            ),
        ),
        "erm_tethers": _force_stats(
            cell,
            lambda: wp.launch(
                erm_tether_force_kernel,
                dim=membrane.n_erm,
                inputs=[
                    cell.pos_d,
                    membrane.erm_m_d,
                    membrane.erm_c_d,
                    membrane.erm_bound_d,
                    wp.float64(membrane.k_erm),
                    membrane.erm_rest_d,
                    wp.float64(membrane.f_rupt),
                ],
                outputs=[cell.f_d],
                device=cell.device,
            ),
        ),
        "total": _force_stats(cell, lambda: _accumulate_all(cell, cell.pos_d, cell.f_d)),
    }
    history_iteration = inner.history_iteration_d.numpy()
    history_displacement = inner.history_displacement_d.numpy()
    history_constraint = inner.history_constraint_d.numpy()
    history_projected_force = inner.history_projected_force_d.numpy()
    reported_iterations = int(report.iters_d.numpy()[0])
    history_valid = (history_iteration > 0) & (history_iteration <= reported_iterations)
    tolerance_um = float(np.sqrt(np.finfo(np.float64).eps) * cell.convergence_length_um)
    dt_mu = float(report.dt_mu_d.numpy()[0])
    convergence_force = float(inner.convergence_force_d.numpy()[0])
    wp.synchronize_device(cell.device)
    try:
        from warp._src.context import runtime as warp_runtime

        warp_device = wp.get_device(cell.device)
        cell.ledger["warp_mempool_used_current_after_solver_bytes"] = int(
            warp_runtime.core.wp_cuda_device_get_mempool_used_mem_current(warp_device.ordinal))
        cell.ledger["warp_mempool_used_high_after_solver_bytes"] = int(
            warp_runtime.core.wp_cuda_device_get_mempool_used_mem_high(warp_device.ordinal))
    except Exception as exc:  # noqa: BLE001 - the native ledger must retain an unavailable reason
        cell.ledger["warp_mempool_after_solver_error"] = f"{type(exc).__name__}: {exc}"
    cell.ledger.update(query_process_peak(cell.device).ledger_fields())

    return {
        "status": "DIAGNOSTIC_NOT_PRODUCTION_RATIFICATION",
        "device": cell.device,
        "config": {
            "n_filaments": cell.cfg.n_filaments,
            "with_steric": cell.steric is not None,
            "erm_density_per_um2": args.erm_density,
            "erm_density_source": args.erm_density_source,
            "inner_iters": args.inner_iters,
            "max_inner_retries": args.max_inner_retries,
            "reshape_every": args.reshape_every,
            "dt_phys_s": args.dt_phys,
            "inner_solver": args.inner_solver,
            "implicit_cg_max_iterations": args.implicit_cg_max_iterations,
            "implicit_line_search_steps": args.implicit_line_search_steps,
            "implicit_coarse_iterations": args.implicit_coarse_iterations,
            "implicit_coarse_modes": args.implicit_coarse_modes,
            "rkc_stages": args.rkc_stages,
        },
        "solver": {
            "converged": bool(report.converged_d.numpy()[0]),
            "iterations": reported_iterations,
            "attempts": int(report.attempts_d.numpy()[0]),
            "residual_candidate_pN": float(report.residual_d.numpy()[0]),
            "max_displacement_um": float(report.max_displacement_d.numpy()[0]),
            "constraint_residual_um": float(report.constraint_residual_d.numpy()[0]),
            "dt_mu": dt_mu,
            "tolerance_um": tolerance_um,
            "displacement_to_tolerance_ratio": float(report.max_displacement_d.numpy()[0]) / tolerance_um,
            "projected_force_convergence_pN": convergence_force,
            "projected_force_tolerance_pN": tolerance_um / dt_mu,
            "projected_force_to_tolerance_ratio": convergence_force * dt_mu / tolerance_um,
            "projected_force_candidate": projected_stats,
            "convergence_history": {
                "iteration": history_iteration[history_valid].tolist(),
                "max_displacement_um": history_displacement[history_valid].tolist(),
                "constraint_residual_um": history_constraint[history_valid].tolist(),
                "projected_force_pN": history_projected_force[history_valid].tolist(),
            },
            "implicit_regularization_pN_per_um": (
                float(inner.implicit_regularization_d.numpy()[0])
                if inner.implicit_regularization_d is not None else None
            ),
            "implicit_last_cg_converged": (
                bool(inner.implicit_cg.converged.numpy()[0]) if inner.implicit_cg is not None else None
            ),
            "implicit_last_cg_iterations": (
                int(inner.implicit_cg.iterations.numpy()[0]) if inner.implicit_cg is not None else None
            ),
            "implicit_last_cg_relative_residual_norm": (
                float(np.sqrt(inner.implicit_cg.rr.numpy()[0] / inner.implicit_cg.rr_initial.numpy()[0]))
                if inner.implicit_cg is not None and inner.implicit_cg.rr_initial.numpy()[0] > 0.0 else None
            ),
            "implicit_guard_accept_count": (
                int(inner.implicit_accept_count_d.numpy()[0])
                if inner.implicit_accept_count_d is not None else None
            ),
            "implicit_last_candidate_residual_pN": (
                float(inner.implicit_residual_d.numpy()[0])
                if inner.implicit_residual_d is not None else None
            ),
            "explicit_last_candidate_residual_pN": (
                float(inner.explicit_residual_d.numpy()[0])
                if inner.explicit_residual_d is not None else None
            ),
            "implicit_best_candidate_residual_pN": (
                float(inner.best_residual_d.numpy()[0])
                if inner.best_residual_d is not None else None
            ),
            "accelerator_guard_accept_count": (
                int(inner.implicit_accept_count_d.numpy()[0])
                if inner.implicit_accept_count_d is not None else None
            ),
            "explicit_candidate_accept_count": (
                int(inner.explicit_accept_count_d.numpy()[0])
                if inner.explicit_accept_count_d is not None else None
            ),
            "stationary_candidate_accept_count": (
                int(inner.stationary_accept_count_d.numpy()[0])
                if inner.stationary_accept_count_d is not None else None
            ),
            "accelerator_last_trial_residual_pN": (
                float(inner.implicit_residual_d.numpy()[0])
                if inner.implicit_residual_d is not None else None
            ),
            "accelerator_best_residual_pN": (
                float(inner.best_residual_d.numpy()[0])
                if inner.best_residual_d is not None else None
            ),
            "anderson_last_coefficient": (
                float(inner.anderson.coefficient.numpy()[0]) if inner.anderson is not None else None
            ),
            "implicit_line_search_scales": list(inner.implicit_line_search_scales),
            "implicit_scale_accept_counts": (
                inner.implicit_scale_accept_counts_d.numpy().tolist()
                if inner.implicit_scale_accept_counts_d is not None else None
            ),
            "accelerator_candidate_families": list(inner.accelerator_candidate_families),
        },
        "timing": {
            "build_wall_s": build_wall_s,
            "accelerator_setup_wall_s": accelerator_setup_wall_s,
            "inner_solve_wall_s": solve_wall_s,
            "inner_solve_wall_s_per_reported_iteration": (
                solve_wall_s / reported_iterations if reported_iterations > 0 else None
            ),
        },
        "ledger": cell.ledger,
        "force_families": families,
        "contact_graph_candidate": contact_graph,
        "last_projection": {
            "max_node_displacement_um": float(projection_delta.max()),
            "rms_node_displacement_um": float(np.sqrt(np.mean(projection_delta * projection_delta))),
            "top_node_index": int(np.argmax(projection_delta)),
            "top_node_family": _node_family(cell, int(np.argmax(projection_delta))),
            "force_families_before_projection": preprojection_families,
        },
        "erm_candidate": {
            "explicit_count": membrane.n_erm,
            "bound_count": int(np.count_nonzero(bound)),
            "extension_um_quantiles": np.quantile(extension[bound], [0.0, 0.5, 0.9, 0.99, 1.0]).tolist(),
            "compressed_count": int(np.count_nonzero(extension[bound] < 0.0)),
            "tensile_count": int(np.count_nonzero(extension[bound] > 0.0)),
            "tension_pN_quantiles": np.quantile(active_tension, [0.0, 0.5, 0.9, 0.99, 1.0]).tolist(),
            "at_or_above_rupture_count": int(np.count_nonzero(active_tension >= membrane.f_rupt)),
        },
        "geometry": {
            "membrane_radius_um": {
                "min": float(mem_radius.min()),
                "mean": float(mem_radius.mean()),
                "max": float(mem_radius.max()),
            },
            "cortex_node_radius_um": {
                "min": float(cortex_radius.min()),
                "mean": float(cortex_radius.mean()),
                "max": float(cortex_radius.max()),
            },
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-filaments", type=int, default=70686)
    parser.add_argument("--erm-density", type=float, required=True)
    parser.add_argument("--erm-density-source", required=True)
    parser.add_argument("--inner-iters", type=int, default=6000)
    parser.add_argument("--max-inner-retries", type=int, default=0)
    parser.add_argument("--reshape-every", type=int, default=20)
    parser.add_argument("--dt-phys", type=float, default=1.0e-4)
    parser.add_argument("--device", default=None)
    parser.add_argument("--no-steric", action="store_true")
    parser.add_argument(
        "--inner-solver",
        choices=("explicit", "analytic_implicit", "block_descent", "anderson", "rkc1",
                 "contact_schwarz", "fiber_contact_schwarz", "tournament", "contact_tournament",
                 "cluster_tournament", "rigid_contact_schwarz", "rigid_cluster_tournament"),
        default="explicit",
    )
    parser.add_argument("--implicit-cg-max-iterations", type=int, default=32)
    parser.add_argument("--implicit-line-search-steps", type=int, default=4)
    parser.add_argument("--implicit-coarse-iterations", type=int, default=8)
    parser.add_argument("--implicit-coarse-modes", type=int, default=0)
    parser.add_argument("--rkc-stages", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run_probe(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
