#!/usr/bin/env python
r"""Full-native C-2 force-family audit at the suspended-cell resting baseline.

This is instrumentation, not a second model and not a parameter sweep.  C-2 currently reports only one
whole-cell scalar, so an operator fix would otherwise be guesswork: the scalar does not say which component
owns the maximum, which force family supplies it, or whether the large term is cancelled before the
inextensibility projection.  This driver builds the SAME full-native configuration and runs the SAME resting
preloads as ``ac_gate_b_interior_column_native.py``, then evaluates each existing force channel separately.

It writes raw and NF2007-projected total-force statistics per disjoint component block, plus every channel's
vector at the node that owns the projected maximum.  It changes no rest length beyond the production-facing
preload functions, launches no physical step, commits no biology and declares no gate threshold.

Sanity Gate:
    * dimensions: every recorded vector and norm is force [pN]; resultants are sums of nodal forces [pN].
    * conservation: the vector sum of separately assembled channels is compared against the ordinary all-
      channel assembly; the relative difference is recorded rather than assumed.
    * boundary/sign: component blocks are derived from the built offsets and asserted disjoint; projection
      changes only actin-fiber constraint directions and leaves nucleus/membrane nodes copied through.
    * numerical: all reductions are descriptive.  No tolerance is introduced and no observed value can pass
      or fail a gate.
    * runtime: CUDA is mandatory.  Host readback occurs only after each force evaluation; there is no physical
      or inner-mechanical loop in this diagnostic.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import CellConfig, build_cell
from aleph.components.incumbent.driver import (
    OMITTABLE_CHANNELS,
    _accumulate_all,
    _preload_cortex_pretension,
    _preload_erm_resting_balance,
    make_inner_solve,
)
from aleph.components.incumbent.inner_mechanics import project_constraint_forces_kernel
from aleph.engine.interior_column_slice import assert_disjoint_blocks


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--filaments", type=int, default=70686)
    parser.add_argument("--membrane-subdiv", type=int, default=6, dest="membrane_subdiv")
    parser.add_argument("--nucleus-subdiv", type=int, default=3, dest="nucleus_subdiv")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--cortex-prestrain", type=float, default=0.0, dest="cortex_prestrain")
    parser.add_argument("--declared-commit", type=str, required=True, dest="declared_commit")
    parser.add_argument("--candidate-solver", type=str, default=None, dest="candidate_solver")
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--outer", type=int, default=40)
    parser.add_argument("--max-inner-retries", type=int, default=0, dest="max_inner_retries")
    parser.add_argument("--line-search-steps", type=int, default=4, dest="line_search_steps")
    parser.add_argument("--line-search-objective", choices=("max_force", "l2_squared"),
                        default="max_force", dest="line_search_objective")
    parser.add_argument("--coarse-modes", type=int, default=0, dest="coarse_modes")
    parser.add_argument("--coarse-iterations", type=int, default=8, dest="coarse_iterations")
    return parser


def _block_stats(force: np.ndarray, blocks: dict[str, tuple[int, int]]) -> dict[str, dict]:
    stats: dict[str, dict] = {}
    for name, (offset, count) in blocks.items():
        values = force[offset:offset + count]
        magnitude = np.linalg.norm(values, axis=1)
        local = int(np.argmax(magnitude)) if magnitude.size else -1
        stats[name] = {
            "node_count": int(count),
            "max_norm_pn": float(magnitude[local]) if local >= 0 else 0.0,
            "max_global_node": int(offset + local) if local >= 0 else None,
            "resultant_pn": [float(v) for v in values.sum(axis=0)],
            "sum_norm_pn": float(magnitude.sum()),
        }
    return stats


def _project(cell, position_d: wp.array, force_d: wp.array) -> np.ndarray:
    projected_d = wp.empty_like(force_d)
    wp.copy(projected_d, force_d)
    diag_d = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=cell.device)
    rhs_d = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=cell.device)
    finite_d = wp.ones(1, dtype=wp.int32, device=cell.device)
    wp.launch(
        project_constraint_forces_kernel,
        dim=cell.n_fibers,
        inputs=[position_d, force_d, cell.foff_d, cell.soff_d, projected_d, diag_d, rhs_d, finite_d],
        device=cell.device,
    )
    wp.synchronize_device(cell.device)
    if int(finite_d.numpy()[0]) != 1:
        raise RuntimeError("NF2007 force projection reported a non-positive pivot")
    return projected_d.numpy()


def main() -> None:
    args = _parser().parse_args()
    wp.init()
    device = wp.get_device()
    if not device.is_cuda:
        raise RuntimeError(f"I0-A: C-2 force audit requires CUDA; resolved {device}")

    cfg = CellConfig(
        n_filaments=args.filaments,
        with_myosin=False,
        overlap_free_cortex=True,
        with_membrane=True,
        membrane_subdivisions=args.membrane_subdiv,
        with_nucleus=True,
        nucleus_subdivisions=args.nucleus_subdiv,
        with_pressure=True,
        erm_radial_pairing=True,
        seed=args.seed,
    )
    cell = build_cell(cfg)
    blocks = {
        "cortex": (0, int(cell.n_actin)),
        "nucleus": (int(cell.nucleus.node_off), int(cell.nucleus.n_verts)),
        "membrane": (int(cell.membrane.node_off), int(cell.membrane.n_verts)),
    }
    assert_disjoint_blocks(blocks, n_total=int(cell.n_total))

    preload = {
        "cortex": _preload_cortex_pretension(cell, args.cortex_prestrain),
        "erm": _preload_erm_resting_balance(cell),
    }
    audited_position_d = cell.pos_d
    candidate = None
    if args.candidate_solver is not None:
        inner_solve = make_inner_solve(
            cell,
            n_inner=args.outer,
            max_inner_retries=args.max_inner_retries,
            capture_candidate=True,
            inner_solver=args.candidate_solver,
            implicit_line_search_steps=args.line_search_steps,
            line_search_objective=args.line_search_objective,
            implicit_coarse_modes=args.coarse_modes,
            implicit_coarse_iterations=args.coarse_iterations,
        )
        report = inner_solve(float(args.dt))
        inner_converged = bool(int(report.converged_d.numpy()[0]))
        if inner_converged:
            raise RuntimeError(
                "candidate force audit is for a rejected C-2 geometry; this solve converged, so use the "
                "accepted-state validation path instead"
            )
        audited_position_d = inner_solve.candidate_pos_d
        candidate = {
            "solver": str(args.candidate_solver),
            "dt_phys_s": float(args.dt),
            "iteration_budget": int(inner_solve.iteration_budget),
            "coarse_modes": int(args.coarse_modes),
            "coarse_iterations": int(args.coarse_iterations),
            "line_search_objective": str(args.line_search_objective),
            "reported_iterations": int(report.iters_d.numpy()[0]),
            "inner_converged": inner_converged,
            "candidate_residual_pn": float(report.residual_d.numpy()[0]),
            "authoritative_geometry_rolled_back": True,
            "audited_geometry": "device_capture_immediately_before_rejected_candidate_rollback",
        }
    # Most incumbent force channels consume the explicit ``pos`` argument, while steric and pressure consume
    # ``cell.state``.  Audit through a shallow cell view that points BOTH routes at the same captured geometry;
    # never swap the authoritative cell's position references, even transiently.
    audit_cell = copy.copy(cell)
    audit_cell.pos_d = audited_position_d
    audit_cell.state = SimpleNamespace(
        pos=audited_position_d,
        node_pos=audited_position_d,
        node_volume=cell.state.node_volume,
    )
    work_d = wp.zeros(audit_cell.n_total, dtype=wp.vec3d, device=audit_cell.device)
    channel_arrays: dict[str, np.ndarray] = {}
    channel_stats: dict[str, dict] = {}
    channels = sorted(OMITTABLE_CHANNELS)
    for channel in channels:
        _accumulate_all(audit_cell, audited_position_d, work_d, omit=OMITTABLE_CHANNELS - {channel})
        wp.synchronize_device(cell.device)
        values = work_d.numpy()
        channel_arrays[channel] = values
        channel_stats[channel] = _block_stats(values, blocks)

    _accumulate_all(audit_cell, audited_position_d, work_d)
    wp.synchronize_device(cell.device)
    total = work_d.numpy()
    projected = _project(audit_cell, audited_position_d, work_d)
    summed = np.sum(np.stack([channel_arrays[name] for name in channels], axis=0), axis=0)
    scale = max(float(np.linalg.norm(total)), np.finfo(float).tiny)
    assembly_relative_error = float(np.linalg.norm(summed - total) / scale)

    projected_norm = np.linalg.norm(projected, axis=1)
    max_node = int(np.argmax(projected_norm))
    owner = next(name for name, (off, count) in blocks.items() if off <= max_node < off + count)
    pos = audited_position_d.numpy()
    cortex_pos = pos[:cell.n_actin]
    cortex_center = cortex_pos.mean(axis=0)
    radial = cortex_pos - cortex_center
    radius = np.linalg.norm(radial, axis=1)
    radial_unit = np.divide(
        radial,
        radius[:, None],
        out=np.zeros_like(radial),
        where=radius[:, None] > np.finfo(float).tiny,
    )
    outward_radial_load = float(np.einsum("ij,ij->", total[:cell.n_actin], radial_unit))
    links = cell.xl_d.numpy().astype(np.int64)
    link_vec = pos[links[:, 1]] - pos[links[:, 0]]
    link_length_sq = np.einsum("ij,ij->i", link_vec, link_vec)
    virial_denominator = float(np.dot(cell.kxl_d.numpy(), link_length_sq))
    mean_radius = float(radius.mean())
    derived_prestrain = (
        outward_radial_load * mean_radius / virial_denominator
        if virial_denominator > np.finfo(float).tiny else None
    )
    at_max = {
        name: {
            "vector_pn": [float(v) for v in channel_arrays[name][max_node]],
            "norm_pn": float(np.linalg.norm(channel_arrays[name][max_node])),
        }
        for name in channels
    }

    record = {
        "schema": "aleph-c2-force-audit@1",
        "kind": "diagnostic",
        "quantitative_claim": "BLOCKED",
        "declared_commit": args.declared_commit,
        "device": str(device),
        "configuration": {
            "filaments": int(args.filaments),
            "membrane_subdiv": int(args.membrane_subdiv),
            "nucleus_subdiv": int(args.nucleus_subdiv),
            "seed": int(args.seed),
            "cortex_prestrain": float(args.cortex_prestrain),
            "with_myosin": False,
            "erm_radial_pairing": True,
        },
        "census": {
            "actin_filaments": int(cell.n_fibers),
            "actin_nodes": int(cell.n_actin),
            "total_nodes": int(cell.n_total),
            "crosslinks": int(cell.n_xl),
            "blocks": {name: {"offset": off, "count": count} for name, (off, count) in blocks.items()},
        },
        "preload": preload,
        "candidate": candidate,
        "assembly_relative_error": assembly_relative_error,
        "raw_total": _block_stats(total, blocks),
        "projected_total": _block_stats(projected, blocks),
        "projected_max": {
            "global_node": max_node,
            "owner": owner,
            "norm_pn": float(projected_norm[max_node]),
            "raw_total_vector_pn": [float(v) for v in total[max_node]],
            "projected_total_vector_pn": [float(v) for v in projected[max_node]],
            "channel_contributions": at_max,
        },
        "laplace_virial_prestrain_derivation": {
            "identity": "sum(k*l^2)*epsilon = mean_radius*sum(F_cortex dot radial_unit)",
            "outward_radial_load_pn": outward_radial_load,
            "mean_radius_um": mean_radius,
            "crosslink_virial_denominator_pn_um": virial_denominator,
            "derived_dimensionless_prestrain": derived_prestrain,
            "status": "DIAGNOSTIC_DERIVATION_NOT_APPLIED",
        },
        "channels": channel_stats,
        "may_not_be_quoted_for": [
            "a force magnitude", "a convergence rate", "an operator property", "cortical tension",
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "out": str(args.out),
        "projected_max": record["projected_max"],
        "projected_blocks": record["projected_total"],
        "assembly_relative_error": assembly_relative_error,
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
