#!/usr/bin/env python3
"""Native foundation gates for live geometry, Biot conservation, pressure preload, and residency.

This runner deliberately distinguishes a component law from the whole physiological contract.  It can prove
the live-mesh classifier, conservative closed-domain updates, true-area membrane flux, device Terzaghi/Green,
pressure surface work, and zero-D2H hot path.  It still cannot call NG-3 complete
unless the full 70,686-filament baseline reaches the registered mechanics convergence criterion.

Every Warp launch is CUDA-only.  NumPy is used only after synchronization as an independent acceptance oracle
or artifact writer, never as simulation state inside a physical-time step.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import (
    BIOT_ALPHA,
    BIOT_MOBILITY,
    BIOT_STORAGE_S,
    PI_0_PA,
    CellConfig,
    build_cell,
)
from aleph.components.incumbent.driver import make_inner_solve
from aleph.components.incumbent.erm_tether import erm_tether_force_kernel
from aleph.components.incumbent.fsi_coupling import SolidDilatationCoupling
from aleph.components.incumbent.membrane_pressure import signed_volume, uniform_pressure_force_reference
from aleph.components.fluid.biot_substrate import BiotSubstrate
from aleph.components.fluid.consolidation_analytic import (
    terzaghi_degree_of_consolidation,
    terzaghi_excess_pressure,
    time_factor,
)
from aleph.components.fluid.field_grid import FLUID, NUCLEUS, OUTSIDE, FieldGrid
from aleph.components.fluid.scheduler import DeviceInnerSolveReport, PhysicalScheduler
from aleph.components.fluid.surface_trace import membrane_surface_pressure_trace

__all__ = [
    "gate_live_domain_parity",
    "gate_ng2_conservation",
    "gate_ng3_pressure_traction",
    "gate_ng6_device_residency",
    "gate_device_retry_multistep",
    "gate_outer_rejection_transaction",
    "run_gates",
]

REL_TOL = 1.0e-10
CONTENT_REL_TOL = 1.0e-11
_FLUID = wp.constant(1)


@wp.kernel
def _affine_surface_kernel(
    pos: wp.array(dtype=wp.vec3d),
    node_off: wp.int32,
    centre: wp.vec3d,
    scale: wp.vec3d,
) -> None:
    """Apply a prescribed volume-preserving diagonal affine map to one live surface on CUDA."""
    i = wp.tid()
    p = pos[node_off + i] - centre
    pos[node_off + i] = centre + wp.vec3d(p[0] * scale[0], p[1] * scale[1], p[2] * scale[2])


@wp.kernel
def _force_first_erm_above_rupture_kernel(
    pos: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    rest: wp.array(dtype=wp.float64),
    k_erm: wp.float64,
    f_rupt: wp.float64,
) -> None:
    """Prescribe one super-threshold ERM extension for the commit/rollback gate."""
    mi = membrane_idx[0]
    ci = cortex_idx[0]
    delta = pos[mi] - pos[ci]
    length = wp.length(delta)
    direction = wp.vec3d(1.0, 0.0, 0.0)
    if length > wp.float64(1.0e-12):
        direction = delta / length
    pos[mi] = pos[ci] + direction * (rest[0] + wp.float64(2.0) * f_rupt / k_erm)


@wp.kernel
def _content_change_kernel(
    p: wp.array3d(dtype=wp.float64),
    p_before: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    storage: wp.float64,
    cell_volume: wp.float64,
    delta: wp.array(dtype=wp.float64),
) -> None:
    """Reduce content change cellwise, avoiding cancellation from subtracting two large total reductions."""
    i, j, k = wp.tid()
    if mask[i, j, k] == _FLUID:
        wp.atomic_add(delta, 0, storage * (p[i, j, k] - p_before[i, j, k]) * cell_volume)


@wp.kernel
def _surface_trace_probe_kernel(
    pos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    origin: wp.vec3d,
    dx: wp.float64,
    pressure_trace: wp.array(dtype=wp.float64),
    valid: wp.array(dtype=wp.int32),
) -> None:
    """Expose the production one-sided surface trace for an independent affine-field oracle."""
    t = wp.tid()
    x0 = pos[faces[t, 0]]
    x1 = pos[faces[t, 1]]
    x2 = pos[faces[t, 2]]
    area2 = wp.cross(x1 - x0, x2 - x0)
    length = wp.length(area2)
    if length <= wp.float64(0.0):
        pressure_trace[t] = wp.float64(0.0)
        valid[t] = 0
        return
    centroid = (x0 + x1 + x2) / wp.float64(3.0)
    trace = membrane_surface_pressure_trace(p, mask, origin, dx, centroid, area2 / length)
    pressure_trace[t] = trace[0]
    valid[t] = wp.int32(trace[1])


@wp.kernel
def _node_work_kernel(
    force: wp.array(dtype=wp.vec3d),
    velocity: wp.array(dtype=wp.vec3d),
    work: wp.array(dtype=wp.float64),
) -> None:
    """Reduce the Lagrangian pressure-force power on CUDA."""
    i = wp.tid()
    wp.atomic_add(work, 0, wp.dot(force[i], velocity[i]))


@wp.kernel
def _grid_pressure_work_kernel(
    gradp: wp.array3d(dtype=wp.vec3d),
    spread_momentum: wp.array3d(dtype=wp.vec3d),
    alpha: wp.float64,
    work: wp.array(dtype=wp.float64),
) -> None:
    """Reduce the Eulerian power paired with ``-alpha grad(p)`` on CUDA."""
    i, j, k = wp.tid()
    wp.atomic_add(work, 0, -alpha * wp.dot(gradp[i, j, k], spread_momentum[i, j, k]))


def _fixture(device: str | None, *, n_filaments: int = 32):
    """Build the production membrane/nucleus/grid topology with a cheap diagnostic cortex."""
    return build_cell(CellConfig(
        n_filaments=n_filaments,
        with_myosin=False,
        with_steric=False,
        with_nucleus=True,
        with_membrane=True,
        with_pressure=True,
        device=device,
    ))


def _grid_points(grid: object) -> np.ndarray:
    idx = np.indices(grid.shape, dtype=np.float64).reshape(3, -1).T
    return np.asarray(grid.origin, dtype=np.float64)[None, :] + float(grid.dx) * idx


def _mesh_inside(points: np.ndarray, verts: np.ndarray, faces: np.ndarray, *, chunk: int = 256) -> np.ndarray:
    """Independent float64 signed-solid-angle occupancy, vectorized in bounded chunks."""
    x = np.asarray(verts, dtype=np.float64)
    f = np.asarray(faces, dtype=np.int64)
    v0, v1, v2 = x[f[:, 0]], x[f[:, 1]], x[f[:, 2]]
    lo, hi = x.min(axis=0), x.max(axis=0)
    candidates = np.flatnonzero(np.all((points >= lo) & (points <= hi), axis=1))
    inside = np.zeros(points.shape[0], dtype=bool)
    for start in range(0, candidates.size, chunk):
        ids = candidates[start:start + chunk]
        q = points[ids, None, :]
        a, b, c = v0[None, :, :] - q, v1[None, :, :] - q, v2[None, :, :] - q
        la = np.linalg.norm(a, axis=2)
        lb = np.linalg.norm(b, axis=2)
        lc = np.linalg.norm(c, axis=2)
        num = np.einsum("qfi,qfi->qf", a, np.cross(b, c))
        den = (
            la * lb * lc
            + np.einsum("qfi,qfi->qf", a, b) * lc
            + np.einsum("qfi,qfi->qf", a, c) * lb
            + np.einsum("qfi,qfi->qf", b, c) * la
        )
        with np.errstate(invalid="ignore"):
            winding = 2.0 * np.sum(np.arctan2(num, den), axis=1) / (4.0 * np.pi)
        inside[ids] = winding >= 0.5
    return inside


def _provider_state(cell: object) -> dict[str, int]:
    membrane = cell.membrane_mask_provider
    nucleus_live = cell.nucleus_mask_provider._live  # gate-only topology identity; never runtime mutation
    return {
        "membrane_mesh_id": int(membrane.mesh.id),
        "membrane_points_ptr": int(membrane.points_d.ptr),
        "nucleus_mesh_id": int(nucleus_live.mesh.id),
        "nucleus_points_ptr": int(nucleus_live.points_d.ptr),
    }


def _classify_against_oracle(cell: object, label: str) -> tuple[dict[str, Any], np.ndarray]:
    membrane_provider = cell.membrane_mask_provider
    nucleus_provider = cell.nucleus_mask_provider
    membrane_provider.reset_diagnostics()
    nucleus_provider._live.reset_diagnostics()  # explicit gate window
    cell.domain.classify()
    wp.synchronize_device(cell.device)

    mask = cell.grid.mask.numpy()
    pos = cell.pos_d.numpy()
    points = _grid_points(cell.grid)
    mem = cell.membrane
    nuc = cell.nucleus
    mem_verts = pos[mem.node_off:mem.node_off + mem.n_verts]
    nuc_verts = pos[nuc.node_off:nuc.node_off + nuc.n_verts]
    mem_inside = _mesh_inside(points, mem_verts, mem.mesh.faces)
    nuc_inside = _mesh_inside(points, nuc_verts, nuc.mesh.faces)
    expected = np.full(points.shape[0], OUTSIDE, dtype=np.int32)
    expected[mem_inside] = FLUID
    expected[mem_inside & nuc_inside] = NUCLEUS
    expected = expected.reshape(cell.grid.shape)

    mem_ties = membrane_provider.surface_tie_mask_d.numpy().astype(bool)
    nuc_ties = nucleus_provider.surface_tie_mask_d.numpy().astype(bool)
    mem_fail = membrane_provider.query_failure_mask_d.numpy().astype(bool)
    nuc_fail = nucleus_provider.query_failure_mask_d.numpy().astype(bool)
    excluded = mem_ties | nuc_ties | mem_fail | nuc_fail
    mismatch = mask != expected
    valid_mismatches = int(np.count_nonzero(mismatch & ~excluded))
    query_failures = int(mem_fail.sum() + nuc_fail.sum())
    report = {
        "label": label,
        "query_failures": query_failures,
        "membrane_surface_ties": int(mem_ties.sum()),
        "nucleus_surface_ties": int(nuc_ties.sum()),
        "valid_cell_mismatches": valid_mismatches,
        "tie_cell_mismatches_reported": int(np.count_nonzero(mismatch & excluded)),
        "fluid_cells": int(np.count_nonzero(mask == FLUID)),
        "nucleus_cells": int(np.count_nonzero(mask == NUCLEUS)),
        "outside_cells": int(np.count_nonzero(mask == OUTSIDE)),
        "pass": query_failures == 0 and valid_mismatches == 0,
    }
    return report, mismatch.astype(np.int8)


def _deform_boundaries(cell: object) -> None:
    # Fixed rational structural perturbation, independent of any gate outcome; determinant is exactly one.
    sx, sy = 51.0 / 50.0, 99.0 / 100.0
    sz = 1.0 / (sx * sy)
    scale = wp.vec3d(sx, sy, sz)
    for comp in (cell.membrane, cell.nucleus):
        verts = cell.pos_d.numpy()[comp.node_off:comp.node_off + comp.n_verts]
        centre = verts.mean(axis=0)
        wp.launch(
            _affine_surface_kernel,
            dim=comp.n_verts,
            inputs=[cell.pos_d, wp.int32(comp.node_off), wp.vec3d(*centre.tolist()), scale],
            device=cell.device,
        )


def gate_live_domain_parity(device: str | None) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Compare live CUDA winding classification with an independent float64 oracle at two geometries."""
    cell = _fixture(device)
    identity_before = _provider_state(cell)
    rest, rest_mismatch = _classify_against_oracle(cell, "rest")
    _deform_boundaries(cell)
    deformed, deformed_mismatch = _classify_against_oracle(cell, "volume_preserving_affine")
    identity_after = _provider_state(cell)
    identity_stable = identity_before == identity_after
    report = {
        "rest": rest,
        "deformed": deformed,
        "allocation_identity_before": identity_before,
        "allocation_identity_after": identity_after,
        "allocation_identity_stable": identity_stable,
        "pass": bool(rest["pass"] and deformed["pass"] and identity_stable),
    }
    mid = cell.grid.shape[2] // 2
    return report, {
        "live_mismatch_rest_midplane": rest_mismatch[:, :, mid],
        "live_mismatch_deformed_midplane": deformed_mismatch[:, :, mid],
    }


def _device_content(grid: object, storage: float) -> float:
    wp.synchronize_device(grid.device)
    return float(grid.total_content_device(storage).numpy()[0])


def _device_content_change(grid: object, p_before: wp.array, storage: float) -> float:
    """Read a post-step device reduction of ``S*sum(p-p_before)*dV`` for a fixed mask."""
    with wp.ScopedDevice(grid.device):
        delta = wp.zeros(1, dtype=wp.float64, device=grid.device)
        wp.launch(
            _content_change_kernel,
            dim=grid.shape,
            inputs=[grid.p, p_before, grid.mask, wp.float64(storage), wp.float64(grid.dx**3)],
            outputs=[delta],
            device=grid.device,
        )
    return float(delta.numpy()[0])


def _relative(error: float, *scales: float) -> float:
    return abs(error) / max(*(abs(v) for v in scales), np.finfo(np.float64).tiny)


def _gate_device_terzaghi(device: str | None) -> dict[str, Any]:
    """Run a doubly-drained 1-D slab through the CUDA Dirichlet operator and compare the closed form."""
    n = 80
    half_h = 5.0
    length = 2.0 * half_h
    dx = length / n
    grid = FieldGrid((n, 1, 1), dx, device=device)
    substrate = BiotSubstrate(grid, mobility=BIOT_MOBILITY, storage_S=BIOT_STORAGE_S, alpha=BIOT_ALPHA)
    grid.set_pressure(np.ones(grid.shape, dtype=np.float64))
    target_s = 0.15
    # The slab has one active spatial direction, so use its exact explicit diffusion limit rather than the
    # conservative 3-D FieldGrid helper.  A factor one-half is the same pre-registered CPU oracle safety.
    dt_limit = BIOT_STORAGE_S * dx * dx / (2.0 * BIOT_MOBILITY)
    n_steps = int(np.ceil(target_s / (0.5 * dt_limit)))
    dt = target_s / n_steps
    for _ in range(n_steps):
        substrate.step_dirichlet_x(dt, 0.0, 0.0)
    wp.synchronize_device(grid.device)
    pressure = grid.p.numpy()[:, 0, 0]
    tv = float(time_factor(substrate.c_v, target_s, half_h))
    z_over_h = (np.arange(n, dtype=np.float64) + 0.5) * dx / half_h
    oracle = terzaghi_excess_pressure(z_over_h, tv, u0=1.0)
    u_numeric = 1.0 - float(pressure.mean())
    u_oracle = float(terzaghi_degree_of_consolidation(np.array([tv]))[0])
    u_abs = abs(u_numeric - u_oracle)
    profile_rmse = float(np.sqrt(np.mean((pressure - oracle) ** 2)))
    passed = bool(u_abs <= 3.0e-3 and profile_rmse <= 5.0e-3)
    return {
        "n_cells": n,
        "n_steps": n_steps,
        "time_factor": tv,
        "degree_numeric": u_numeric,
        "degree_oracle": u_oracle,
        "degree_absolute_error": u_abs,
        "profile_rmse": profile_rmse,
        "pass": passed,
    }


def _gate_device_green(device: str | None) -> dict[str, Any]:
    """Run the registered compact Gaussian packet and read the CUDA spread slope ``2*d*c_v``."""
    n = 61
    dx = 0.6
    grid = FieldGrid((n, n, n), dx, device=device)
    substrate = BiotSubstrate(grid, mobility=BIOT_MOBILITY, storage_S=BIOT_STORAGE_S, alpha=BIOT_ALPHA)
    xs = (np.arange(n, dtype=np.float64) - n // 2) * dx
    xx, yy, zz = np.meshgrid(xs, xs, xs, indexing="ij")
    r2 = xx**2 + yy**2 + zz**2
    sigma0 = 2.0 * dx
    pressure = np.exp(-r2 / (2.0 * sigma0**2))
    grid.set_pressure(np.ascontiguousarray(pressure))

    def mean_sq_radius(field: np.ndarray) -> float:
        weights = np.clip(field, 0.0, None)
        return float(np.sum(r2 * weights) / np.sum(weights))

    r2_start = mean_sq_radius(pressure)
    dt = substrate.cfl_dt(safety=0.5)
    n_steps = 60
    for _ in range(n_steps):
        substrate.step(dt)
    wp.synchronize_device(grid.device)
    pressure_end = grid.p.numpy()
    r2_end = mean_sq_radius(pressure_end)
    measured = (r2_end - r2_start) / (n_steps * dt)
    expected = 2.0 * 3.0 * substrate.c_v
    relative_error = abs(measured - expected) / expected
    return {
        "n_cells_per_axis": n,
        "n_steps": n_steps,
        "measured_d_r2_dt_um2_per_s": measured,
        "expected_2dcv_um2_per_s": expected,
        "relative_error": relative_error,
        "pass": bool(relative_error <= 0.05),
    }


def gate_ng2_conservation(device: str | None) -> dict[str, Any]:
    """Run the NG-2 component suite and the canonical accepted-transaction device ledger."""
    cell = _fixture(device)
    g, sub = cell.grid, cell.substrate
    mask = g.mask.numpy()
    fluid = mask == FLUID
    shape = g.shape
    dt = 0.5 * sub.cfl_dt()
    zeros = np.zeros(shape, dtype=np.float64)

    # Constant state is bit-preserving on a frozen impermeable live domain.
    p_uniform = np.full(shape, PI_0_PA, dtype=np.float64)
    g.set_pressure(p_uniform)
    g.s_water.assign(zeros)
    g.s_membrane.zero_()
    g.div_vs.assign(zeros)
    for _ in range(4):
        sub.step(dt)
    p_after_constant = g.p.numpy()
    constant_max_abs = float(np.max(np.abs(p_after_constant[fluid] - PI_0_PA)))
    constant_pass = bool(np.array_equal(p_after_constant[fluid], p_uniform[fluid]))

    # Smooth non-uniform pressure, no source/flux: device total content must close.
    idx = np.indices(shape, dtype=np.float64)
    p0 = PI_0_PA + np.sin((idx[0] + 1.0) / shape[0]) + np.cos((idx[1] + 1.0) / shape[1])
    g.set_pressure(np.ascontiguousarray(p0))
    before = _device_content(g, sub.storage_S)
    for _ in range(4):
        sub.step(dt)
    after = _device_content(g, sub.storage_S)
    impermeable_rel = _relative(after - before, before, after)
    impermeable_pass = impermeable_rel <= CONTENT_REL_TOL

    # One-step source + solid-dilatation balance.
    g.set_pressure(p_uniform)
    source = np.zeros(shape, dtype=np.float64)
    dilatation = np.zeros(shape, dtype=np.float64)
    source[fluid] = np.sin((idx[0][fluid] + 1.0) / shape[0]) * 1.0e-2
    dilatation[fluid] = np.cos((idx[2][fluid] + 1.0) / shape[2]) * 1.0e-2
    g.s_water.assign(np.ascontiguousarray(source))
    g.s_membrane.zero_()
    g.div_vs.assign(np.ascontiguousarray(dilatation))
    before = _device_content(g, sub.storage_S)
    sub.step(dt)
    after = _device_content(g, sub.storage_S)
    source_delta = after - before
    expected_delta = dt * float(np.sum((source - BIOT_ALPHA * dilatation)[fluid])) * g.dx**3
    balance_error = source_delta - expected_delta
    balance_rel = _relative(balance_error, source_delta, expected_delta)
    balance_pass = balance_rel <= REL_TOL

    # Nonzero Kedem-Katchalsky flux on the true live-triangle area.  Doubling the resting osmotic difference is
    # the natural factor-two boundary case (zero at one Pi, positive at two Pi), not a fitted physiological value.
    g.set_pressure(p_uniform)
    g.s_water.assign(zeros)
    g.div_vs.assign(zeros)
    d_pi_flux = 2.0 * PI_0_PA
    cell.membrane_bc.reset_diagnostics()
    cell.membrane_bc.apply(d_pi_flux)
    membrane_source = g.s_membrane.numpy()
    source_integral = float(np.sum(membrane_source[fluid])) * g.dx**3
    integrated_surface_flux = float(cell.membrane_bc.integrated_flux_d.numpy()[0])
    flux_unresolved = int(cell.membrane_bc.unresolved_faces_d.numpy()[0])
    mem = cell.membrane
    pos = cell.pos_d.numpy()[mem.node_off:mem.node_off + mem.n_verts]
    area2 = np.cross(pos[mem.mesh.faces[:, 1]] - pos[mem.mesh.faces[:, 0]],
                     pos[mem.mesh.faces[:, 2]] - pos[mem.mesh.faces[:, 0]])
    live_area = 0.5 * float(np.linalg.norm(area2, axis=1).sum())
    analytic_surface_flux = cell.membrane_bc.L_p * (d_pi_flux - PI_0_PA) * live_area
    p_before_flux_d = wp.clone(g.p)
    sub.step(dt)
    content_flux_delta = _device_content_change(g, p_before_flux_d, sub.storage_S)
    flux_surface_rel = _relative(
        integrated_surface_flux - analytic_surface_flux, integrated_surface_flux, analytic_surface_flux)
    source_surface_rel = _relative(
        source_integral - integrated_surface_flux, source_integral, integrated_surface_flux)
    content_surface_rel = _relative(
        content_flux_delta - dt * integrated_surface_flux, content_flux_delta, dt * integrated_surface_flux)
    membrane_flux_pass = bool(
        cell.membrane_bc.discretization == "LIVE_TRIANGLE_CONSERVATIVE"
        and flux_unresolved == 0
        and flux_surface_rel <= REL_TOL
        and source_surface_rel <= REL_TOL
        and content_surface_rel <= REL_TOL
    )

    # Moving-domain closure from uniform state; remap accounting must equal the discrete content change.
    g.set_pressure(p_uniform)
    g.p_bar.assign(p_uniform)
    g.s_water.assign(zeros)
    g.div_vs.assign(zeros)
    before = _device_content(g, sub.storage_S)
    cell.domain.bind_storage(sub.storage_S)
    _deform_boundaries(cell)
    moving_d = cell.domain.remap()
    moving = float(moving_d.numpy()[0])
    after = _device_content(g, sub.storage_S)
    moving_error = (after - before) - moving
    moving_rel = _relative(moving_error, after - before, moving)
    p_moved = g.p.numpy()
    moved_constant_max_abs = float(np.max(np.abs(p_moved[g.mask.numpy() == FLUID] - PI_0_PA)))
    moving_pass = moving_rel <= REL_TOL and moved_constant_max_abs <= 1.0e-12

    terzaghi = _gate_device_terzaghi(device)
    green = _gate_device_green(device)
    coupled = _gate_canonical_coupled_ledger(device)
    implemented_pass = bool(
        constant_pass and impermeable_pass and balance_pass and membrane_flux_pass and moving_pass)
    component_pass = bool(implemented_pass and terzaghi["pass"] and green["pass"])
    contract_pass = bool(component_pass and coupled["pass"])
    return {
        "constant_state": {"max_abs_Pa": constant_max_abs, "pass": constant_pass},
        "impermeable": {"relative_content_drift": impermeable_rel, "pass": impermeable_pass},
        "source_and_dilatation": {
            "delta_content": source_delta,
            "expected_delta": expected_delta,
            "balance_error": balance_error,
            "relative_error": balance_rel,
            "pass": balance_pass,
        },
        "moving_domain": {
            "moving_face_delta": moving,
            "balance_error": moving_error,
            "relative_error": moving_rel,
            "constant_state_max_abs_Pa": moved_constant_max_abs,
            "pass": moving_pass,
        },
        "membrane_content_equals_flux": {
            "discretization": cell.membrane_bc.discretization,
            "live_triangle_area_um2": live_area,
            "unresolved_faces": flux_unresolved,
            "integrated_surface_flux_um3_per_s": integrated_surface_flux,
            "analytic_surface_flux_um3_per_s": analytic_surface_flux,
            "grid_source_integral_um3_per_s": source_integral,
            "content_delta": content_flux_delta,
            "expected_content_delta": dt * integrated_surface_flux,
            "surface_quadrature_relative_error": flux_surface_rel,
            "grid_deposition_relative_error": source_surface_rel,
            "content_balance_relative_error": content_surface_rel,
            "pass": membrane_flux_pass,
        },
        "device_terzaghi_dirichlet": terzaghi,
        "device_green_impulse": green,
        "implemented_subgates_pass": implemented_pass,
        "component_suite_status": "PASS" if component_pass else "FAIL",
        "canonical_coupled_outer_step": coupled,
        "status": "PASS" if contract_pass else ("FAIL_COMPONENT" if not component_pass else "FAIL_COUPLED"),
    }


def _gate_canonical_coupled_ledger(device: str | None) -> dict[str, Any]:
    """Close all nonzero NG-2 terms in one accepted canonical scheduler transaction."""
    cell = _fixture(device)
    g = cell.grid
    mask = g.mask.numpy()
    fluid = mask == FLUID
    pressure = np.full(g.shape, PI_0_PA, dtype=np.float64)
    interior = np.zeros(g.shape, dtype=np.float64)
    # Fixed boundary-case loads, not fitted production parameters: both independent channels are nonzero and
    # have opposite signs in the registered Biot balance, which makes a sign swap immediately visible.
    interior[fluid] = 1.0e-2
    g.set_pressure(pressure)
    g.p_bar.assign(pressure)
    g.s_water.assign(interior)
    g.s_membrane.zero_()
    # Seed div(v_s)^n through the production FSI map, not by assigning its Eulerian output. Grid-point nodes
    # give the four-point Peskin kernel exact affine first-order reproduction over the enclosed cell domain.
    fsi_pos = _grid_points(g)
    fsi_centre = np.asarray(g.origin) + 0.5 * (np.asarray(g.shape) - 1.0) * g.dx
    dilation_rate = 1.0 / 8.0
    fsi_velocity = dilation_rate * (fsi_pos - fsi_centre[None, :])
    fsi_volume = np.full(fsi_pos.shape[0], g.dx**3, dtype=np.float64)
    fsi_active = np.ones(fsi_pos.shape[0], dtype=np.int32)
    fsi_coupling = SolidDilatationCoupling(g)
    fsi_coupling.update(
        wp.array(np.ascontiguousarray(fsi_pos), dtype=wp.vec3d, device=cell.device),
        wp.array(np.ascontiguousarray(fsi_velocity), dtype=wp.vec3d, device=cell.device),
        wp.array(fsi_volume, dtype=wp.float64, device=cell.device),
        wp.array(fsi_active, dtype=wp.int32, device=cell.device),
    )
    cell.membrane_bc.reset_diagnostics()
    cell.membrane_mask_provider.reset_diagnostics()
    cell.nucleus_mask_provider._live.reset_diagnostics()

    pos = cell.pos_d.numpy()
    surfaces: list[tuple[object, wp.vec3d]] = []
    for comp in (cell.membrane, cell.nucleus):
        centre = pos[comp.node_off:comp.node_off + comp.n_verts].mean(axis=0)
        surfaces.append((comp, wp.vec3d(*centre.tolist())))
    with wp.ScopedDevice(cell.device):
        zero_f64 = wp.zeros(1, dtype=wp.float64, device=cell.device)
        zero_i32 = wp.zeros(1, dtype=wp.int32, device=cell.device)
        one_i32 = wp.ones(1, dtype=wp.int32, device=cell.device)

    def accepted_affine_inner(_dt_phys: float, _fluid_valid_d: wp.array) -> DeviceInnerSolveReport:
        # A rational 2% isotropic expansion makes the moving-domain term nonzero without introducing a fitted
        # material parameter. It is an acceptance oracle executed entirely on CUDA inside the scheduler.
        scale = wp.vec3d(51.0 / 50.0, 51.0 / 50.0, 51.0 / 50.0)
        for comp, centre in surfaces:
            wp.launch(
                _affine_surface_kernel,
                dim=comp.n_verts,
                inputs=[cell.pos_d, wp.int32(comp.node_off), centre, scale],
                device=cell.device,
            )
        return DeviceInnerSolveReport(
            residual_d=zero_f64,
            iters_d=zero_i32,
            attempts_d=one_i32,
            converged_d=one_i32,
            max_displacement_d=zero_f64,
            constraint_residual_d=zero_f64,
            finite_d=one_i32,
            dt_mu_d=zero_f64,
        )

    accepted_affine_inner.failure_counters = (  # type: ignore[attr-defined]
        cell.membrane_bc.unresolved_faces_d,
        cell.membrane_mask_provider.query_failures_d,
        cell.nucleus_mask_provider.query_failures_d,
    )
    scheduler = PhysicalScheduler(
        substrate=cell.substrate,
        domain=cell.domain,
        membrane_bc=cell.membrane_bc,
        inner_solve=accepted_affine_inner,
        osmotic_difference=lambda _t: 2.0 * PI_0_PA,
    )
    report = scheduler.outer_step(1.0e-3).readback()
    scale = max(
        abs(report.content_delta),
        abs(report.membrane_surface_flux_integral),
        abs(report.membrane_grid_source_integral),
        abs(report.interior_source_integral),
        abs(report.solid_dilatation_integral),
        abs(report.moving_face_delta),
        np.finfo(np.float64).tiny,
    )
    solver_rel = abs(report.solver_conservation_error) / scale
    physical_rel = abs(report.physical_conservation_error) / scale
    deposition_rel = abs(report.membrane_deposition_error) / scale
    nonzero_terms = {
        "membrane_surface": bool(
            abs(report.membrane_surface_flux_integral) > np.finfo(np.float64).tiny),
        "membrane_grid": bool(abs(report.membrane_grid_source_integral) > np.finfo(np.float64).tiny),
        "interior": bool(abs(report.interior_source_integral) > np.finfo(np.float64).tiny),
        "solid_dilatation": bool(abs(report.solid_dilatation_integral) > np.finfo(np.float64).tiny),
        "moving_face": bool(abs(report.moving_face_delta) > np.finfo(np.float64).tiny),
    }
    passed = bool(
        report.outer_accepted
        and not report.outer_rolled_back
        and report.t == report.dt_phys
        and all(nonzero_terms.values())
        and solver_rel <= REL_TOL
        and physical_rel <= REL_TOL
        and deposition_rel <= REL_TOL
    )
    return {
        "time_index_contract": (
            "ledger consumes div(v_s)^n during Biot subcycles; post-remap div(v_s)^(n+1) is next-step state"
        ),
        "outer_accepted": report.outer_accepted,
        "outer_rolled_back": report.outer_rolled_back,
        "accepted_time_s": report.t,
        "n_biot_subcycles": report.n_biot_subcycles,
        "content_start": report.content_start,
        "content_end": report.total_content,
        "content_delta": report.content_delta,
        "membrane_surface_flux_integral": report.membrane_surface_flux_integral,
        "membrane_grid_source_integral": report.membrane_grid_source_integral,
        "interior_source_integral": report.interior_source_integral,
        "solid_dilatation_integral": report.solid_dilatation_integral,
        "moving_face_delta": report.moving_face_delta,
        "nonzero_terms": nonzero_terms,
        "membrane_deposition_error": report.membrane_deposition_error,
        "solver_conservation_error": report.solver_conservation_error,
        "physical_conservation_error": report.physical_conservation_error,
        "membrane_deposition_relative_error": deposition_rel,
        "solver_conservation_relative_error": solver_rel,
        "physical_conservation_relative_error": physical_rel,
        "status": "PASS" if passed else "FAIL",
        "pass": passed,
    }


def _fsi_conservation_contract(device: str | None) -> dict[str, Any]:
    """Run CUDA rigid-motion, affine-dilation, and partition-of-unity FSI transfer oracles."""
    n = 10
    dx = 0.5
    origin = (-0.5 * (n - 1) * dx,) * 3
    grid = FieldGrid((n, n, n), dx, origin=origin, device=device)
    coupling = SolidDilatationCoupling(grid)
    idx = np.indices(grid.shape, dtype=np.float64).reshape(3, -1).T
    pos = np.asarray(origin, dtype=np.float64)[None, :] + dx * idx
    volume = np.full(pos.shape[0], dx**3, dtype=np.float64)
    active = np.ones(pos.shape[0], dtype=np.int32)

    def upload_velocity(velocity: np.ndarray) -> tuple[wp.array, wp.array, wp.array, wp.array]:
        return (
            wp.array(np.ascontiguousarray(pos), dtype=wp.vec3d, device=grid.device),
            wp.array(np.ascontiguousarray(velocity), dtype=wp.vec3d, device=grid.device),
            wp.array(volume, dtype=wp.float64, device=grid.device),
            wp.array(active, dtype=wp.int32, device=grid.device),
        )

    rigid_v = np.array([0.7, -0.3, 0.2], dtype=np.float64)
    rigid = np.broadcast_to(rigid_v, pos.shape).copy()
    coupling.update(*upload_velocity(rigid))
    wp.synchronize_device(grid.device)
    rigid_grid = coupling._vs_grid.numpy()  # gate-only observation of production CUDA scratch
    rigid_div = grid.div_vs.numpy()
    rigid_velocity_error = float(np.max(np.abs(rigid_grid - rigid_v)))
    rigid_div_error = float(np.max(np.abs(rigid_div)))
    rigid_scale = max(float(np.max(np.abs(rigid_v))) / dx, np.finfo(np.float64).tiny)
    rigid_pass = bool(rigid_velocity_error / max(np.max(np.abs(rigid_v)), 1.0) <= REL_TOL
                      and rigid_div_error / rigid_scale <= REL_TOL)

    dilation_rate = 0.125
    affine = dilation_rate * pos
    coupling.update(*upload_velocity(affine))
    wp.synchronize_device(grid.device)
    affine_grid = coupling._vs_grid.numpy()
    affine_div = grid.div_vs.numpy()
    velocity_core = (slice(2, -2),) * 3
    divergence_core = (slice(3, -3),) * 3
    affine_expected = affine.reshape(grid.shape + (3,))
    affine_velocity_error = float(np.max(np.abs(affine_grid[velocity_core] - affine_expected[velocity_core])))
    expected_divergence = 3.0 * dilation_rate
    affine_div_error = float(np.max(np.abs(affine_div[divergence_core] - expected_divergence)))
    measured_integral = float(np.sum(affine_div[divergence_core])) * dx**3
    expected_integral = expected_divergence * float(np.prod(affine_div[divergence_core].shape)) * dx**3
    affine_integral_rel = _relative(measured_integral - expected_integral, measured_integral, expected_integral)
    affine_pass = bool(
        affine_velocity_error / max(float(np.max(np.abs(affine))), 1.0) <= REL_TOL
        and affine_div_error / expected_divergence <= REL_TOL
        and affine_integral_rel <= REL_TOL
    )

    interior_selector = np.all((idx >= 2) & (idx <= n - 3), axis=1)
    pos_i = pos[interior_selector]
    vel_i = np.column_stack((
        0.2 + 0.03 * pos_i[:, 0],
        -0.1 + 0.02 * pos_i[:, 1],
        0.4 - 0.01 * pos_i[:, 2],
    ))
    vol_i = np.full(pos_i.shape[0], dx**3, dtype=np.float64)
    active_i = np.ones(pos_i.shape[0], dtype=np.int32)
    coupling.update(
        wp.array(np.ascontiguousarray(pos_i), dtype=wp.vec3d, device=grid.device),
        wp.array(np.ascontiguousarray(vel_i), dtype=wp.vec3d, device=grid.device),
        wp.array(vol_i, dtype=wp.float64, device=grid.device),
        wp.array(active_i, dtype=wp.int32, device=grid.device),
    )
    wp.synchronize_device(grid.device)
    spread_weight = float(np.sum(coupling._weight.numpy()))
    spread_momentum = np.sum(coupling._momentum.numpy(), axis=(0, 1, 2))
    node_weight = float(np.sum(vol_i))
    node_momentum = np.sum(vel_i * vol_i[:, None], axis=0)
    weight_rel = _relative(spread_weight - node_weight, spread_weight, node_weight)
    momentum_rel = float(np.linalg.norm(spread_momentum - node_momentum) /
                         max(np.linalg.norm(spread_momentum), np.linalg.norm(node_momentum),
                             np.finfo(np.float64).tiny))
    spread_pass = bool(weight_rel <= REL_TOL and momentum_rel <= REL_TOL)
    passed = bool(rigid_pass and affine_pass and spread_pass)
    return {
        "rigid_translation_zero_divergence": {
            "velocity_max_abs_error_um_per_s": rigid_velocity_error,
            "divergence_max_abs_error_per_s": rigid_div_error,
            "pass": rigid_pass,
        },
        "affine_dilation_integral_identity": {
            "dilation_rate_per_s": dilation_rate,
            "velocity_core_max_abs_error_um_per_s": affine_velocity_error,
            "divergence_core_max_abs_error_per_s": affine_div_error,
            "measured_integral_um3_per_s": measured_integral,
            "expected_integral_um3_per_s": expected_integral,
            "relative_error": affine_integral_rel,
            "pass": affine_pass,
        },
        "spread_partition_conservation": {
            "node_weight_um3": node_weight,
            "grid_weight_um3": spread_weight,
            "node_momentum_um4_per_s": node_momentum.tolist(),
            "grid_momentum_um4_per_s": spread_momentum.tolist(),
            "weight_relative_error": weight_rel,
            "momentum_relative_error": momentum_rel,
            "pass": spread_pass,
        },
        "status": "PASS" if passed else "FAIL",
        "pass": passed,
    }


def _ng5_closed_system_contract(device: str | None) -> dict[str, Any]:
    """Gate pressure-force closure and the actual CUDA Peskin work transpose at native cortex density.

    Sanity Gate:
        * Dimensions: bulk ``-alpha grad(p) V`` and surface ``p n dA`` are both pN; their sum is an
          internal-force residual.  Nodal and Eulerian powers are both pN um/s because the raw spread stores
          ``sum(V_n v_n phi)`` rather than a density.
        * Boundary cases: the uniform closed-surface case is already certified by NG-3.  A nonzero affine
          pressure gradient makes both the bulk and surface totals nonzero and opposite in a closed system.
        * Conservation/work: ``sum_n f_n dot v_n`` must equal
          ``sum_grid (-alpha grad(p)) dot sum_n(V_n v_n phi)`` to round-off for the production kernels.
        * Numerical: the fixed structural tolerance is ``REL_TOL``; no physiological magnitude is fitted.
        * Residency: force interpolation, velocity spread, and both work reductions execute on CUDA.  Host
          reads occur only after synchronization for this acceptance report.
    """
    cell = build_cell(CellConfig(
        n_filaments=70686,
        with_myosin=True,
        with_steric=True,
        with_nucleus=True,
        with_membrane=True,
        with_pressure=True,
        device=device,
    ))
    if cell.pressure is None or cell.membrane_pressure is None or cell.grid is None or cell.membrane is None:
        raise RuntimeError("NG-5 requires the complete pressure/membrane production topology")

    gradient = np.array([1.0, -0.5, 0.25], dtype=np.float64)
    points = _grid_points(cell.grid).reshape((*cell.grid.shape, 3))
    p_linear = PI_0_PA + np.einsum("...i,i->...", points, gradient)
    cell.grid.set_pressure(np.ascontiguousarray(p_linear))

    cell.f_d.zero_()
    cell.pressure.accumulate(cell.state, cell.f_d)
    wp.synchronize_device(cell.device)
    bulk_force = cell.f_d.numpy().copy()

    cell.membrane_pressure.reset_diagnostics()
    cell.membrane_pressure.accumulate(cell.pos_d, cell.f_d)
    wp.synchronize_device(cell.device)
    combined_force = cell.f_d.numpy().copy()
    unresolved = int(cell.membrane_pressure.unresolved_faces_d.numpy()[0])
    surface_force = combined_force - bulk_force

    bulk_total = np.sum(bulk_force, axis=0)
    surface_total = np.sum(surface_force, axis=0)
    combined_total = np.sum(combined_force, axis=0)
    node_volume = cell.node_volume_d.numpy()
    node_volume_sum = float(np.sum(node_volume))
    pos = cell.pos_d.numpy()
    mem = cell.membrane
    membrane_volume = signed_volume(
        pos[mem.node_off:mem.node_off + mem.n_verts],
        mem.mesh.faces,
    )
    expected_bulk = -BIOT_ALPHA * node_volume_sum * gradient
    expected_surface = membrane_volume * gradient
    force_scale = max(
        float(np.linalg.norm(bulk_total)),
        float(np.linalg.norm(surface_total)),
        np.finfo(np.float64).tiny,
    )
    net_force_rel = float(np.linalg.norm(combined_total) / force_scale)
    bulk_oracle_rel = float(np.linalg.norm(bulk_total - expected_bulk) /
                            max(np.linalg.norm(expected_bulk), np.finfo(np.float64).tiny))
    surface_oracle_rel = float(np.linalg.norm(surface_total - expected_surface) /
                               max(np.linalg.norm(expected_surface), np.finfo(np.float64).tiny))
    force_pass = bool(
        unresolved == 0
        and np.linalg.norm(bulk_total) > 0.0
        and np.linalg.norm(surface_total) > 0.0
        and bulk_oracle_rel <= REL_TOL
        and surface_oracle_rel <= REL_TOL
        and net_force_rel <= REL_TOL
    )

    # Choose the dissipative response direction v_n = f_n / V_n on the actin skeleton.  This makes the work
    # sign non-vacuous (strictly positive) while leaving the transpose identity sensitive to every Peskin
    # weight used independently by PressureCoupling and SolidDilatationCoupling.
    velocity = np.zeros_like(bulk_force)
    active_volume = node_volume > 0.0
    velocity[active_volume] = bulk_force[active_volume] / node_volume[active_volume, None]
    velocity_d = wp.array(np.ascontiguousarray(velocity), dtype=wp.vec3d, device=cell.device)
    fsi = SolidDilatationCoupling(cell.grid)
    fsi.update(cell.pos_d, velocity_d, cell.node_volume_d, cell.solid_active_d)
    node_work_d = wp.zeros(1, dtype=wp.float64, device=cell.device)
    grid_work_d = wp.zeros(1, dtype=wp.float64, device=cell.device)
    wp.launch(
        _node_work_kernel,
        dim=cell.n_total,
        inputs=[cell.f_d, velocity_d],
        outputs=[node_work_d],
        device=cell.device,
    )
    wp.launch(
        _grid_pressure_work_kernel,
        dim=cell.grid.shape,
        inputs=[cell.pressure._gradp, fsi._momentum, wp.float64(BIOT_ALPHA)],
        outputs=[grid_work_d],
        device=cell.device,
    )
    wp.synchronize_device(cell.device)
    node_work = float(node_work_d.numpy()[0])
    grid_work = float(grid_work_d.numpy()[0])
    work_rel = _relative(node_work - grid_work, node_work, grid_work)
    work_pass = bool(node_work > 0.0 and grid_work > 0.0 and work_rel <= REL_TOL)

    passed = bool(force_pass and work_pass)
    return {
        "native_cortical_filaments": cell.n_fibers,
        "native_actin_nodes": cell.n_actin,
        "uniform_pressure_closed_surface": "COVERED_BY_NG3_COMPONENT",
        "spatial_affine_pressure_bulk_plus_surface_net_force": {
            "gradient_Pa_per_um": gradient.tolist(),
            "bulk_total_pN": bulk_total.tolist(),
            "surface_total_pN": surface_total.tolist(),
            "combined_total_pN": combined_total.tolist(),
            "node_control_volume_um3": node_volume_sum,
            "membrane_discrete_volume_um3": membrane_volume,
            "bulk_oracle_relative_error": bulk_oracle_rel,
            "surface_oracle_relative_error": surface_oracle_rel,
            "net_force_relative_error": net_force_rel,
            "unresolved_membrane_faces": unresolved,
            "pass": force_pass,
        },
        "fsi_transfer_adjoint_work_identity": {
            "node_work_pN_um_per_s": node_work,
            "grid_work_pN_um_per_s": grid_work,
            "relative_error": work_rel,
            "nonnegative_dissipation_pass": bool(node_work > 0.0 and grid_work > 0.0),
            "pass": work_pass,
        },
        "status": "PASS" if passed else "FAIL_OPEN",
        "pass": passed,
    }


def _pressure_case(cell: object, pressure: float) -> tuple[np.ndarray, int]:
    cell.grid.set_pressure(np.full(cell.grid.shape, pressure, dtype=np.float64))
    cell.f_d.zero_()
    cell.membrane_pressure.reset_diagnostics()
    cell.membrane_pressure.accumulate(cell.pos_d, cell.f_d)
    wp.synchronize_device(cell.device)
    force = cell.f_d.numpy()
    unresolved = int(cell.membrane_pressure.unresolved_faces_d.numpy()[0])
    mem = cell.membrane
    return force[mem.node_off:mem.node_off + mem.n_verts], unresolved


def gate_ng3_pressure_traction(
    device: str | None, native_resting: dict[str, Any] | None,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Certify the uniform pressure surface law, then defer whole-t0 equilibrium to the native report."""
    cell = _fixture(device)
    pos = cell.pos_d.numpy()
    mem = cell.membrane
    verts = pos[mem.node_off:mem.node_off + mem.n_verts]
    ref = uniform_pressure_force_reference(verts, mem.mesh.faces, PI_0_PA)
    f_zero, unresolved_zero = _pressure_case(cell, 0.0)
    f_pos, unresolved_pos = _pressure_case(cell, PI_0_PA)
    f_neg, unresolved_neg = _pressure_case(cell, -PI_0_PA)
    scale = max(float(np.linalg.norm(ref, axis=1).max()), np.finfo(np.float64).tiny)
    parity_rel = float(np.linalg.norm(f_pos - ref, axis=1).max() / scale)
    reversal_rel = float(np.linalg.norm(f_neg + f_pos, axis=1).max() / scale)
    closure_rel = float(np.linalg.norm(f_pos.sum(axis=0)) / np.linalg.norm(f_pos, axis=1).sum())
    virial = float(np.sum(f_pos * verts))
    volume6 = float(np.einsum(
        "ij,ij->i", verts[mem.mesh.faces[:, 0]],
        np.cross(verts[mem.mesh.faces[:, 1]], verts[mem.mesh.faces[:, 2]])).sum())
    virial_ref = 3.0 * PI_0_PA * volume6 / 6.0
    virial_rel = _relative(virial - virial_ref, virial, virial_ref)
    uniform_pass = bool(
        np.array_equal(f_zero, np.zeros_like(f_zero))
        and unresolved_zero == unresolved_pos == unresolved_neg == 0
        and parity_rel < REL_TOL and reversal_rel < REL_TOL
        and closure_rel < REL_TOL and virial_rel < REL_TOL
    )

    # Affine trace gate: a continuum linear pressure field must be reproduced at the actual surface centroid,
    # not merely for the uniform special case.
    gradient = np.array([1.0, -0.5, 0.25], dtype=np.float64)
    points = _grid_points(cell.grid).reshape((*cell.grid.shape, 3))
    p_linear = PI_0_PA + np.einsum("...i,i->...", points, gradient)
    cell.grid.set_pressure(np.ascontiguousarray(p_linear))
    cell.f_d.zero_()
    cell.membrane_pressure.reset_diagnostics()
    cell.membrane_pressure.accumulate(cell.pos_d, cell.f_d)
    wp.synchronize_device(cell.device)
    f_linear = cell.f_d.numpy()[mem.node_off:mem.node_off + mem.n_verts]
    unresolved_linear = int(cell.membrane_pressure.unresolved_faces_d.numpy()[0])
    faces = mem.mesh.faces
    area2_linear = np.cross(verts[faces[:, 1]] - verts[faces[:, 0]],
                            verts[faces[:, 2]] - verts[faces[:, 0]])
    centroids = (verts[faces[:, 0]] + verts[faces[:, 1]] + verts[faces[:, 2]]) / 3.0
    p_trace = PI_0_PA + centroids @ gradient
    trace_device_d = wp.empty(mem.n_faces, dtype=wp.float64, device=cell.device)
    trace_valid_d = wp.empty(mem.n_faces, dtype=wp.int32, device=cell.device)
    wp.launch(
        _surface_trace_probe_kernel,
        dim=mem.n_faces,
        inputs=[cell.pos_d, mem.faces_d, cell.grid.p, cell.grid.mask,
                wp.vec3d(*cell.grid.origin.tolist()), wp.float64(cell.grid.dx)],
        outputs=[trace_device_d, trace_valid_d],
        device=cell.device,
    )
    wp.synchronize_device(cell.device)
    trace_device = trace_device_d.numpy()
    trace_valid = trace_valid_d.numpy()
    trace_scale = max(float(np.abs(p_trace).max()), np.finfo(np.float64).tiny)
    pressure_trace_rel = float(np.max(np.abs(trace_device - p_trace)) / trace_scale)
    per_vertex = p_trace[:, None] * area2_linear / 6.0
    f_linear_ref = np.zeros_like(verts)
    np.add.at(f_linear_ref, faces[:, 0], per_vertex)
    np.add.at(f_linear_ref, faces[:, 1], per_vertex)
    np.add.at(f_linear_ref, faces[:, 2], per_vertex)
    linear_scale = max(float(np.linalg.norm(f_linear_ref, axis=1).max()), np.finfo(np.float64).tiny)
    linear_rel = float(np.linalg.norm(f_linear - f_linear_ref, axis=1).max() / linear_scale)
    linear_pass = bool(
        unresolved_linear == 0 and np.all(trace_valid == 1)
        and pressure_trace_rel < REL_TOL and linear_rel < REL_TOL)
    component_pass = bool(uniform_pass and linear_pass)
    native_root = native_resting or {}
    native_report = native_root.get("report", native_root)
    native_ledger = native_root.get("ledger", native_report.get("ledger", {}))
    native_config = native_root.get("config", {})
    native_converged = bool(native_report.get("inner_converged", False))
    native_population = int(native_ledger.get("N_unique_active_fibers", 0))
    config_pass = bool(
        int(native_config.get("n_filaments", 0)) == 70686
        and all(bool(native_config.get(key, False)) for key in (
            "with_myosin", "with_steric", "with_nucleus", "with_membrane", "with_pressure",
        ))
    )
    population_pass = bool(
        native_population == 70686
        and int(native_ledger.get("n_actin_nodes", 0)) == 494802
        and int(native_ledger.get("n_total_nodes", 0)) >= 494802
    )
    diagnostics_pass = bool(
        native_converged
        and bool(native_report.get("stable", False))
        and bool(native_report.get("outer_accepted", False))
        and not bool(native_report.get("outer_rolled_back", True))
        and bool(native_report.get("inner_finite", False))
        and bool(native_report.get("fluid_finite", False))
        and bool(native_report.get("pos_finite", False))
        and bool(native_report.get("force_finite", False))
        and int(native_report.get("membrane_flux_unresolved_faces", -1)) == 0
        and int(native_report.get("membrane_pressure_unresolved_faces", -1)) == 0
        and int(native_report.get("membrane_mask_query_failures", -1)) == 0
        and int(native_report.get("nucleus_mask_query_failures", -1)) == 0
    )
    diagnostic_capacity_feasible = bool(
        native_ledger.get("preload_mechanically_feasible_upper_bound", False))
    preload_capacity_feasible = bool(
        diagnostic_capacity_feasible
        and native_ledger.get("preload_capacity_force_basis") == "SOURCE_GROUNDED_SINGLE_ERM_FORCE"
    )
    preload_density_sourced = bool(
        native_ledger.get("erm_pairing_mode") == "EXPLICIT_DENSITY_RESOLVED"
        and native_ledger.get("erm_density_source") not in (None, "", "ABSENT_FROM_CONTRACT_GRAPH")
        and native_ledger.get("erm_density_mcf7_production") is True
    )
    preload_kinetics_sourced = bool(
        native_ledger.get("erm_kinetics_mode") == "BELL_SLIP_ON_OFF"
        and native_ledger.get("erm_kinetics_source") not in (None, "", "ABSENT_FROM_CONTRACT_GRAPH")
        and all(
            isinstance(native_ledger.get(key), (int, float)) and float(native_ledger[key]) > 0.0
            for key in (
                "erm_bell_k_on_s",
                "erm_bell_k_off0_s",
                "erm_bell_force_pn",
                "erm_capture_radius_um",
            )
        )
        and native_ledger.get("erm_rebind_rest_policy") == "formation_length"
    )
    preload_contract_pass = bool(
        preload_capacity_feasible and preload_density_sourced and preload_kinetics_sourced)
    physiological_pass = bool(
        config_pass and population_pass and diagnostics_pass and preload_contract_pass)
    return {
        "surface_law": {
            "zero_is_bit_zero": bool(np.array_equal(f_zero, np.zeros_like(f_zero))),
            "unresolved_faces": max(unresolved_zero, unresolved_pos, unresolved_neg),
            "device_oracle_relative_error": parity_rel,
            "sign_reversal_relative_error": reversal_rel,
            "closed_surface_force_relative_error": closure_rel,
            "virial_relative_error": virial_rel,
            "pass": uniform_pass,
        },
        "young_laplace_required_tension_pN_per_um": PI_0_PA * 7.5 / 2.0,
        "spatial_linear_field_trace": {
            "gradient_Pa_per_um": gradient.tolist(),
            "unresolved_faces": unresolved_linear,
            "invalid_trace_faces": int(np.count_nonzero(trace_valid != 1)),
            "pressure_trace_relative_error": pressure_trace_rel,
            "device_oracle_relative_error": linear_rel,
            "pass": linear_pass,
            "status": "PASS" if linear_pass else "FAIL_OPEN",
        },
        "native_inner_converged": native_converged,
        "native_cortical_filaments": native_population,
        "native_config_pass": config_pass,
        "native_population_pass": population_pass,
        "native_diagnostics_pass": diagnostics_pass,
        "preload_capacity": {
            "membrane_area_um2": native_ledger.get("preload_membrane_area_um2"),
            "membrane_laplace_pressure_Pa": native_ledger.get("preload_membrane_laplace_pressure_pa"),
            "erm_tether_count": native_ledger.get("preload_erm_tether_count"),
            "erm_density_per_um2": native_ledger.get("preload_erm_density_per_um2"),
            "diagnostic_membrane_tether_force_pN": native_ledger.get("preload_erm_rupture_force_pn"),
            "pressure_capacity_upper_bound_Pa": native_ledger.get(
                "preload_total_pressure_capacity_upper_bound_pa"),
            "pressure_capacity_ratio": native_ledger.get("preload_pressure_capacity_ratio"),
            "required_erm_density_min_per_um2": native_ledger.get(
                "preload_required_erm_density_min_per_um2"),
            "required_erm_count_min": native_ledger.get("preload_required_erm_count_min"),
            "force_basis": native_ledger.get("preload_capacity_force_basis", "MISSING"),
            "diagnostic_membrane_tether_force_capacity_feasible": diagnostic_capacity_feasible,
            "mechanically_feasible_upper_bound": preload_capacity_feasible,
            "density_source": native_ledger.get("erm_density_source", "MISSING"),
            "mcf7_production_provenance": bool(
                native_ledger.get("erm_density_mcf7_production", False)),
            "density_sourced": preload_density_sourced,
            "kinetics_mode": native_ledger.get("erm_kinetics_mode", "MISSING"),
            "kinetics_source": native_ledger.get("erm_kinetics_source", "MISSING"),
            "bell_k_on_s": native_ledger.get("erm_bell_k_on_s"),
            "bell_k_off0_s": native_ledger.get("erm_bell_k_off0_s"),
            "bell_force_pn": native_ledger.get("erm_bell_force_pn"),
            "capture_radius_um": native_ledger.get("erm_capture_radius_um"),
            "rebind_rest_policy": native_ledger.get("erm_rebind_rest_policy", "MISSING"),
            "kinetics_sourced": preload_kinetics_sourced,
            "pass": preload_contract_pass,
        },
        "physiological_t0_pass": physiological_pass,
        "status": "PASS" if component_pass and physiological_pass else "FAIL_OPEN",
    }, {
        "pressure_force_device_norm": np.linalg.norm(f_pos, axis=1),
        "pressure_force_oracle_norm": np.linalg.norm(ref, axis=1),
        "pressure_linear_device_norm": np.linalg.norm(f_linear, axis=1),
        "pressure_linear_oracle_norm": np.linalg.norm(f_linear_ref, axis=1),
    }


def gate_ng6_device_residency(
    device: str | None, native_resting: dict[str, Any] | None,
) -> dict[str, Any]:
    """Use Warp's directional memcpy profiler with a deliberate DtoH positive control."""
    cell = _fixture(device)
    inner = make_inner_solve(cell, n_inner=1, reshape_every=1)
    warmup_scheduler = PhysicalScheduler(
        substrate=cell.substrate,
        domain=cell.domain,
        membrane_bc=cell.membrane_bc,
        inner_solve=inner,
        osmotic_difference=lambda _t: PI_0_PA,
    )
    # Compile/allocate outside the measured window.
    warmup_scheduler.outer_step(1.0e-3)
    # Use a fresh instance so the measured attempt starts at the same accepted t=0 state after the warm-up
    # candidate was rolled back; this isolates compilation/allocation from the residency window.
    scheduler = PhysicalScheduler(
        substrate=cell.substrate,
        domain=cell.domain,
        membrane_bc=cell.membrane_bc,
        inner_solve=inner,
        osmotic_difference=lambda _t: PI_0_PA,
    )
    wp.synchronize_device(cell.device)
    wp.timing_begin(wp.TIMING_MEMCPY)
    report_d = scheduler.outer_step(1.0e-3)
    records = wp.timing_end(synchronize=True)
    hot_dtoh = [record for record in records if record.name == "memcpy DtoH"]
    wp.timing_begin(wp.TIMING_MEMCPY)
    report_d.inner_residual_d.numpy()
    positive_records = wp.timing_end(synchronize=True)
    positive_dtoh = [record for record in positive_records if record.name == "memcpy DtoH"]
    small_pass = len(hot_dtoh) == 0 and len(positive_dtoh) > 0

    native_report = (native_resting or {}).get("report", native_resting or {})
    native_available = "hot_loop_dtoh_copies" in native_report
    native_pass = bool(
        native_available
        and int(native_report["hot_loop_dtoh_copies"]) == 0
        and int(native_report.get("profiler_positive_control_dtoh_copies", 0)) > 0
    )
    return {
        "small_fixture": {
            "hot_loop_dtoh_copies": len(hot_dtoh),
            "hot_loop_dtoh_elapsed_ms": float(sum(record.elapsed for record in hot_dtoh)),
            "positive_control_dtoh_copies": len(positive_dtoh),
            "pass": small_pass,
        },
        "full_native": {
            "available": native_available,
            "hot_loop_dtoh_copies": native_report.get("hot_loop_dtoh_copies"),
            "positive_control_dtoh_copies": native_report.get("profiler_positive_control_dtoh_copies"),
            "pass": native_pass,
        },
        "status": "PASS" if small_pass and native_pass else "INCOMPLETE",
    }


def gate_device_retry_multistep(device: str | None) -> dict[str, Any]:
    """Prove same-time retry, convergence no-op, and outer fail-stop control on CUDA.

    The retry budget is numerical only.  A force-free candidate must converge in its first chunk and leave all
    later fixed launches as no-ops.  An unconverged coupled candidate must consume every configured chunk,
    roll back, and device-disable the remaining pre-scheduled outer steps.  A control without a mechanical
    failure must commit every physical step and reach exactly the prescribed ramp time.
    """
    force_free_cfg = CellConfig(
        n_filaments=32,
        with_myosin=False,
        with_steric=False,
        with_nucleus=False,
        with_membrane=False,
        with_pressure=False,
        device=device,
    )
    force_free_control = build_cell(force_free_cfg)
    control_inner = make_inner_solve(
        force_free_control, n_inner=1, reshape_every=1, max_inner_retries=0)
    control_d = control_inner(1.0e-3)
    wp.synchronize_device(force_free_control.device)
    control_converged = bool(control_d.converged_d.numpy()[0])
    control_pos = force_free_control.pos_d.numpy().copy()

    force_free = build_cell(force_free_cfg)
    short_circuit_inner = make_inner_solve(
        force_free, n_inner=1, reshape_every=1, max_inner_retries=2)
    short_d = short_circuit_inner(1.0e-3)
    wp.synchronize_device(force_free.device)
    short_converged = bool(short_d.converged_d.numpy()[0])
    short_attempts = int(short_d.attempts_d.numpy()[0])
    short_iters = int(short_d.iters_d.numpy()[0])
    short_exact = bool(np.array_equal(force_free.pos_d.numpy(), control_pos))
    short_pass = bool(
        control_converged and short_converged and short_attempts == 1 and short_iters == 1 and short_exact)

    # A retry is a compute-budget partition, not another numerical operator. Three iterations in one chunk and
    # one iteration in each of three chunks must therefore agree within the registered float64 position
    # resolution; independent CUDA atomic accumulation need not be bit-identical across two fresh builds.
    continuous = _fixture(device)
    continuous_inner = make_inner_solve(
        continuous, n_inner=3, reshape_every=1, max_inner_retries=0, capture_candidate=True)
    continuous_inner(1.0e-3)
    partitioned = _fixture(device)
    partitioned_inner = make_inner_solve(
        partitioned, n_inner=1, reshape_every=1, max_inner_retries=2, capture_candidate=True)
    partitioned_inner(1.0e-3)
    wp.synchronize_device(partitioned.device)
    continuous_candidate = continuous_inner.candidate_pos_d.numpy()
    partitioned_candidate = partitioned_inner.candidate_pos_d.numpy()
    trajectory_delta_um = float(np.linalg.norm(
        continuous_candidate - partitioned_candidate, axis=1).max())
    trajectory_tolerance_um = float(
        np.sqrt(np.finfo(np.float64).eps) * continuous.convergence_length_um)
    trajectory_pass = bool(trajectory_delta_um <= trajectory_tolerance_um)

    rejected = _fixture(device)
    rejected_grid = rejected.grid
    rejected_before = {
        "pos": rejected.pos_d.numpy().copy(),
        "p": rejected_grid.p.numpy().copy(),
        "p_new": rejected_grid.p_new.numpy().copy(),
        "mask": rejected_grid.mask.numpy().copy(),
        "div_vs": rejected_grid.div_vs.numpy().copy(),
        "s_membrane": rejected_grid.s_membrane.numpy().copy(),
        "s_total": rejected_grid.s_total.numpy().copy(),
    }
    rejected_scheduler = PhysicalScheduler(
        substrate=rejected.substrate,
        domain=rejected.domain,
        membrane_bc=rejected.membrane_bc,
        inner_solve=make_inner_solve(
            rejected, n_inner=1, reshape_every=1, max_inner_retries=2),
        osmotic_difference=lambda _t: PI_0_PA,
    )
    rejected_reports = [report.readback() for report in rejected_scheduler.run(3.0e-3, 3)]
    rejected_after = {
        "pos": rejected.pos_d.numpy(),
        "p": rejected_grid.p.numpy(),
        "p_new": rejected_grid.p_new.numpy(),
        "mask": rejected_grid.mask.numpy(),
        "div_vs": rejected_grid.div_vs.numpy(),
        "s_membrane": rejected_grid.s_membrane.numpy(),
        "s_total": rejected_grid.s_total.numpy(),
    }
    rejected_exact = {
        name: bool(np.array_equal(rejected_before[name], rejected_after[name])) for name in rejected_before
    }
    first, *disabled = rejected_reports
    retry_exhaustion_pass = bool(
        first.outer_step_enabled
        and not first.outer_accepted
        and first.outer_rolled_back
        and first.inner_attempts == 3
        and first.inner_iters == 3
        and first.t == 0.0
        and all(not report.outer_step_enabled for report in disabled)
        and all(not report.outer_accepted and report.outer_rolled_back for report in disabled)
        and all(report.inner_attempts == 0 and report.inner_iters == 0 for report in disabled)
        and all(report.t == 0.0 for report in rejected_reports)
        and all(rejected_exact.values())
    )

    accepted = _fixture(device)
    accepted_scheduler = PhysicalScheduler(
        substrate=accepted.substrate,
        domain=accepted.domain,
        membrane_bc=accepted.membrane_bc,
        osmotic_difference=lambda _t: PI_0_PA,
    )
    accepted_reports = [report.readback() for report in accepted_scheduler.run(3.0e-3, 3)]
    expected_times = np.arange(1, 4, dtype=np.float64) * 1.0e-3
    accepted_times = np.asarray([report.t for report in accepted_reports])
    accepted_pass = bool(
        all(report.outer_step_enabled and report.outer_accepted and not report.outer_rolled_back
            for report in accepted_reports)
        and np.allclose(accepted_times, expected_times, rtol=0.0, atol=np.finfo(np.float64).eps)
    )

    passed = bool(short_pass and trajectory_pass and retry_exhaustion_pass and accepted_pass)
    return {
        "converged_first_chunk": {
            "configured_attempts": 3,
            "executed_attempts": short_attempts,
            "first_converged_iteration": short_iters,
            "matches_one_chunk_control_bit_exact": short_exact,
            "pass": short_pass,
        },
        "chunk_partition_trajectory": {
            "continuous_iterations": 3,
            "partitioned_iterations": [1, 1, 1],
            "candidate_max_difference_um": trajectory_delta_um,
            "fixed_numerical_tolerance_um": trajectory_tolerance_um,
            "pass": trajectory_pass,
        },
        "retry_exhaustion_fail_stop": {
            "configured_attempts": 3,
            "first_step_executed_attempts": first.inner_attempts,
            "first_step_iterations": first.inner_iters,
            "step_enabled": [report.outer_step_enabled for report in rejected_reports],
            "outer_accepted": [report.outer_accepted for report in rejected_reports],
            "accepted_time_s": [report.t for report in rejected_reports],
            "bit_exact_restoration": rejected_exact,
            "pass": retry_exhaustion_pass,
        },
        "accepted_multistep_control": {
            "step_enabled": [report.outer_step_enabled for report in accepted_reports],
            "outer_accepted": [report.outer_accepted for report in accepted_reports],
            "accepted_time_s": accepted_times.tolist(),
            "prescribed_ramp_time_s": 3.0e-3,
            "pass": accepted_pass,
        },
        "status": "PASS" if passed else "FAIL",
    }


def gate_outer_rejection_transaction(device: str | None) -> dict[str, Any]:
    """Prove an unconverged coupled attempt commits no geometry, field, source, mask, or physical time."""
    cell = _fixture(device)
    g = cell.grid
    mask_h = g.mask.numpy()
    idx = np.indices(g.shape, dtype=np.float64)
    interior = np.zeros(g.shape, dtype=np.float64)
    interior[mask_h == FLUID] = np.sin((idx[0][mask_h == FLUID] + 1.0) / g.shape[0]) * 1.0e-2
    membrane_sentinel = np.full(g.shape, 1.0e-3, dtype=np.float64)
    total_sentinel = np.full(g.shape, 2.0e-3, dtype=np.float64)
    g.s_water.assign(np.ascontiguousarray(interior))
    g.s_membrane.assign(membrane_sentinel)
    g.s_total.assign(total_sentinel)
    before = {
        "pos": cell.pos_d.numpy().copy(),
        "p": g.p.numpy().copy(),
        "p_new": g.p_new.numpy().copy(),
        "mask": mask_h.copy(),
        "div_vs": g.div_vs.numpy().copy(),
        "s_membrane": g.s_membrane.numpy().copy(),
        "s_total": g.s_total.numpy().copy(),
    }
    scheduler = PhysicalScheduler(
        substrate=cell.substrate,
        domain=cell.domain,
        membrane_bc=cell.membrane_bc,
        inner_solve=make_inner_solve(cell, n_inner=1, reshape_every=1),
        osmotic_difference=lambda _t: PI_0_PA,
    )
    report = scheduler.outer_step(1.0e-3).readback()
    after = {
        "pos": cell.pos_d.numpy(),
        "p": g.p.numpy(),
        "p_new": g.p_new.numpy(),
        "mask": g.mask.numpy(),
        "div_vs": g.div_vs.numpy(),
        "s_membrane": g.s_membrane.numpy(),
        "s_total": g.s_total.numpy(),
    }
    exact = {name: bool(np.array_equal(before[name], after[name])) for name in before}
    unconverged_pass = bool(
        not report.outer_accepted
        and report.outer_rolled_back
        and report.t == 0.0
        and report.moving_face_delta == 0.0
        and report.content_delta == 0.0
        and report.membrane_surface_flux_integral == 0.0
        and report.membrane_grid_source_integral == 0.0
        and report.interior_source_integral == 0.0
        and report.solid_dilatation_integral == 0.0
        and report.membrane_deposition_error == 0.0
        and report.solver_conservation_error == 0.0
        and report.physical_conservation_error == 0.0
        and all(exact.values())
    )

    # A non-finite interior source must reject even without a mechanical solver that could incidentally expose
    # it through force interpolation.  Pressure is restored from a finite device snapshot.
    nan_cell = _fixture(device)
    nan_grid = nan_cell.grid
    nan_source = np.zeros(nan_grid.shape, dtype=np.float64)
    first_fluid = tuple(np.argwhere(nan_grid.mask.numpy() == FLUID)[0])
    nan_source[first_fluid] = np.nan
    nan_grid.s_water.assign(nan_source)
    nan_pos_before = nan_cell.pos_d.numpy().copy()
    nan_p_before = nan_grid.p.numpy().copy()
    nan_p_new_before = nan_grid.p_new.numpy().copy()
    nan_scheduler = PhysicalScheduler(
        substrate=nan_cell.substrate,
        domain=nan_cell.domain,
        membrane_bc=nan_cell.membrane_bc,
        inner_solve=make_inner_solve(nan_cell, n_inner=1, reshape_every=1),
        osmotic_difference=lambda _t: PI_0_PA,
    )
    nan_report = nan_scheduler.outer_step(1.0e-3).readback()
    nonfinite_pass = bool(
        not nan_report.fluid_finite
        and not nan_report.outer_accepted
        and nan_report.outer_rolled_back
        and nan_report.t == 0.0
        and nan_report.inner_iters == 0
        and np.array_equal(nan_cell.pos_d.numpy(), nan_pos_before)
        and np.array_equal(nan_grid.p.numpy(), nan_p_before)
        and np.array_equal(nan_grid.p_new.numpy(), nan_p_new_before)
    )

    # Invalid prescribed time is rejected before the transaction-attempt latch or any D2D snapshot/mutation.
    dt_cell = _fixture(device)
    dt_p_before = dt_cell.grid.p.numpy().copy()
    dt_scheduler = PhysicalScheduler(substrate=dt_cell.substrate, domain=dt_cell.domain)
    invalid_dt_raised = False
    try:
        dt_scheduler.outer_step(0.0)
    except ValueError:
        invalid_dt_raised = True
    invalid_dt_pass = bool(
        invalid_dt_raised and dt_scheduler._t == 0.0
        and int(dt_scheduler._run_active_d.numpy()[0]) == 1  # gate-only post-call inspection
        and np.array_equal(dt_cell.grid.p.numpy(), dt_p_before)
    )

    # ERM rupture is irreversible biology and must not be driven by dimensionless inner iterations.  Force one
    # tether to exactly twice its derived rupture force, reject the mechanical candidate, and require bound
    # state plus geometry to remain bit-exact.  A separate explicit outer-boundary commit must then rupture it.
    erm_cell = _fixture(device)
    erm = erm_cell.membrane
    if erm is None or erm.n_erm == 0:
        raise RuntimeError("ERM transaction gate requires at least one explicit membrane-cortex tether")
    wp.launch(
        _force_first_erm_above_rupture_kernel,
        dim=1,
        inputs=[erm_cell.pos_d, erm.erm_m_d, erm.erm_c_d, erm.erm_rest_d,
                wp.float64(erm.k_erm), wp.float64(erm.f_rupt)],
        device=erm_cell.device,
    )
    wp.synchronize_device(erm_cell.device)
    erm_pos_before = erm_cell.pos_d.numpy().copy()
    erm_bound_before = erm.erm_bound_d.numpy().copy()
    mi = int(erm.erm_m_d.numpy()[0])
    ci = int(erm.erm_c_d.numpy()[0])
    forced_tension = erm.k_erm * (
        float(np.linalg.norm(erm_pos_before[mi] - erm_pos_before[ci])) - float(erm.erm_rest_d.numpy()[0]))
    force_probe_d = wp.zeros(erm_cell.n_total, dtype=wp.vec3d, device=erm_cell.device)
    wp.launch(
        erm_tether_force_kernel,
        dim=1,
        inputs=[erm_cell.pos_d, erm.erm_m_d, erm.erm_c_d, erm.erm_bound_d,
                wp.float64(erm.k_erm), erm.erm_rest_d, wp.float64(erm.f_rupt)],
        outputs=[force_probe_d],
        device=erm_cell.device,
    )
    wp.synchronize_device(erm_cell.device)
    force_probe = force_probe_d.numpy()
    candidate_force_zero = bool(np.array_equal(force_probe, np.zeros_like(force_probe)))
    force_evaluation_state_exact = bool(np.array_equal(erm.erm_bound_d.numpy(), erm_bound_before))
    erm_scheduler = PhysicalScheduler(
        substrate=erm_cell.substrate,
        domain=erm_cell.domain,
        membrane_bc=erm_cell.membrane_bc,
        inner_solve=make_inner_solve(erm_cell, n_inner=1, reshape_every=1),
        osmotic_difference=lambda _t: PI_0_PA,
    )
    erm_report = erm_scheduler.outer_step(1.0e-3).readback()
    erm_bound_after_rejection = erm.erm_bound_d.numpy().copy()
    erm_pos_after_rejection = erm_cell.pos_d.numpy().copy()
    rejected_bound_exact = bool(np.array_equal(erm_bound_after_rejection, erm_bound_before))
    rejected_geometry_exact = bool(np.array_equal(erm_pos_after_rejection, erm_pos_before))
    rejected_erm_pass = bool(
        forced_tension > erm.f_rupt
        and candidate_force_zero
        and force_evaluation_state_exact
        and not erm_report.outer_accepted
        and erm_report.outer_rolled_back
        and rejected_bound_exact
        and rejected_geometry_exact
    )
    accepted_d = wp.ones(1, dtype=wp.int32, device=erm_cell.device)
    erm.rupture_step(erm_cell.pos_d, accepted_d)
    wp.synchronize_device(erm_cell.device)
    explicit_commit_pass = bool(
        int(np.count_nonzero(erm_bound_before)) - int(np.count_nonzero(erm.erm_bound_d.numpy())) >= 1)
    erm_transaction_pass = bool(rejected_erm_pass and explicit_commit_pass)

    passed = bool(unconverged_pass and nonfinite_pass and invalid_dt_pass and erm_transaction_pass)
    return {
        "unconverged_candidate": {
            "outer_accepted": report.outer_accepted,
            "outer_rolled_back": report.outer_rolled_back,
            "accepted_time_s": report.t,
            "attempted_time_s": report.attempted_t,
            "moving_face_delta": report.moving_face_delta,
            "committed_content_delta": report.content_delta,
            "committed_membrane_surface_flux_integral": report.membrane_surface_flux_integral,
            "committed_membrane_grid_source_integral": report.membrane_grid_source_integral,
            "committed_interior_source_integral": report.interior_source_integral,
            "committed_solid_dilatation_integral": report.solid_dilatation_integral,
            "committed_solver_conservation_error": report.solver_conservation_error,
            "committed_physical_conservation_error": report.physical_conservation_error,
            "bit_exact_restoration": exact,
            "pass": unconverged_pass,
        },
        "nonfinite_fluid": {
            "fluid_finite": nan_report.fluid_finite,
            "outer_accepted": nan_report.outer_accepted,
            "outer_rolled_back": nan_report.outer_rolled_back,
            "accepted_time_s": nan_report.t,
            "inner_updates_executed": nan_report.inner_iters,
            "geometry_unchanged_bit_exact": bool(np.array_equal(nan_cell.pos_d.numpy(), nan_pos_before)),
            "pressure_restored_bit_exact": bool(np.array_equal(nan_grid.p.numpy(), nan_p_before)),
            "pressure_scratch_restored_bit_exact": bool(
                np.array_equal(nan_grid.p_new.numpy(), nan_p_new_before)),
            "pass": nonfinite_pass,
        },
        "invalid_dt_no_mutation": {
            "raised_before_attempt": invalid_dt_raised,
            "pressure_unchanged": bool(np.array_equal(dt_cell.grid.p.numpy(), dt_p_before)),
            "pass": invalid_dt_pass,
        },
        "erm_commit_only_rupture": {
            "forced_tension_to_threshold_ratio": forced_tension / erm.f_rupt,
            "super_threshold_candidate_force_zero": candidate_force_zero,
            "force_evaluation_bound_state_bit_exact": force_evaluation_state_exact,
            "candidate_outer_accepted": erm_report.outer_accepted,
            "candidate_outer_rolled_back": erm_report.outer_rolled_back,
            "rejected_bound_state_bit_exact": rejected_bound_exact,
            "rejected_geometry_bit_exact": rejected_geometry_exact,
            "explicit_commit_ruptured_tethers": int(np.count_nonzero(erm_bound_before))
            - int(np.count_nonzero(erm.erm_bound_d.numpy())),
            "pass": erm_transaction_pass,
        },
        "pass": passed,
    }


def _native_measurement_contract(native_resting: dict[str, Any] | None) -> dict[str, Any]:
    """Audit NG-9 timing/population/peak-byte evidence without relabelling partial counters as exact."""
    root = native_resting or {}
    report = root.get("report", root)
    ledger = root.get("ledger", report.get("ledger", {}))
    available = bool(root)
    population_pass = bool(
        int(ledger.get("N_unique_active_fibers", 0)) == 70686
        and int(ledger.get("n_actin_nodes", 0)) == 494802
        and int(ledger.get("n_total_nodes", 0)) >= 494802
    )
    outer_timing_pass = bool(float(report.get("outer_step_wall_s", 0.0)) > 0.0)
    exact_peak_bytes = ledger.get("gpu_bytes_peak_exact")
    exact_peak_pass = bool(
        ledger.get("gpu_memory_accounting_status") == "EXACT_WHOLE_DEVICE_PEAK"
        and isinstance(exact_peak_bytes, int)
        and exact_peak_bytes > 0
    )
    status = "PASS" if available and population_pass and outer_timing_pass and exact_peak_pass else "INCOMPLETE"
    return {
        "available": available,
        "population_pass": population_pass,
        "outer_step_wall_time_pass": outer_timing_pass,
        "outer_step_wall_s": report.get("outer_step_wall_s"),
        "exact_peak_gpu_bytes_pass": exact_peak_pass,
        "memory_accounting_status": ledger.get("gpu_memory_accounting_status", "MISSING"),
        "memory_accounting_probe_status": ledger.get("gpu_memory_accounting_probe_status", "MISSING"),
        "exact_peak_gpu_bytes": exact_peak_bytes,
        "mempool_high_water_bytes": ledger.get("warp_mempool_used_high_bytes"),
        "reason": (
            "Warp mempool high-water and endpoint device usage do not include/prove every CUDA allocation's "
            "whole-step peak; inspect memory_accounting_probe_status for the NVML exact-counter blocker"
        ) if not exact_peak_pass else "exact NVML process-lifetime peak was recorded in bytes",
        "status": status,
    }


def run_gates(device: str | None, native_resting: dict[str, Any] | None) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    live, live_arrays = gate_live_domain_parity(device)
    ng2 = gate_ng2_conservation(device)
    ng3, pressure_arrays = gate_ng3_pressure_traction(device, native_resting)
    ng6 = gate_ng6_device_residency(device, native_resting)
    retry_multistep = gate_device_retry_multistep(device)
    rejection = gate_outer_rejection_transaction(device)
    native_measurement = _native_measurement_contract(native_resting)
    fsi_conservation = _fsi_conservation_contract(device)
    ng5_closed_system = _ng5_closed_system_contract(device)
    results = {
        "schema": "ffn-ac-foundation-gates-v2",
        "device": str(wp.get_device(device)),
        "live_domain": live,
        "ng2": ng2,
        "ng3": ng3,
        "ng6": ng6,
        "device_retry_multistep": retry_multistep,
        "outer_rejection_transaction": rejection,
        "fsi_conservation": fsi_conservation,
        "ng5_closed_system_force": ng5_closed_system,
        "ng9_native_measurement": native_measurement,
    }
    complete = (
        live["pass"] and rejection["pass"] and retry_multistep["status"] == "PASS"
        and ng2["status"] == "PASS" and ng3["status"] == "PASS"
        and ng6["status"] == "PASS" and fsi_conservation["status"] == "PASS"
        and ng5_closed_system["status"] == "PASS" and native_measurement["status"] == "PASS"
    )
    results["overall_status"] = "PASS" if complete else "INCOMPLETE_OR_FAIL"
    return results, {**live_arrays, **pressure_arrays}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default=None)
    parser.add_argument("--native-resting-json", default="")
    parser.add_argument("--json", required=True)
    parser.add_argument("--npz", default="")
    args = parser.parse_args()
    wp.init()
    dev = wp.get_device(args.device)
    if not dev.is_cuda:
        raise RuntimeError("foundation native gates require CUDA; CPU simulation is forbidden")
    native = None
    if args.native_resting_json:
        native = json.loads(Path(args.native_resting_json).read_text(encoding="utf-8"))
    results, arrays = run_gates(args.device, native)
    out = Path(args.json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")
    if args.npz:
        npz = Path(args.npz)
        npz.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(npz, **arrays)
    print(json.dumps(results, indent=2, sort_keys=True))
    return 0 if results["overall_status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
