"""GPU-resident outer physical clock + inner mechanical solve orchestration (I1a/P5).

One physical-time step advances the conservative Biot field at its own CFL, relaxes the solid mechanics at the
new fluid state, and then remaps the fluid domain to the newly moved membrane/nuclear meshes.  The membrane
Kedem-Katchalsky source is rebuilt before *every* Biot subcycle because it depends on the evolving local
pressure.  Remapping after mechanics removes the historical one-outer-step geometry lag.

The hot path returns device scalars.  Explicit :meth:`DeviceOuterStepReport.readback` is the only conversion to
Python numbers and is called by a runner after (never from inside) the physical-time loop.  This separation is
the enforceable NG-6 contract: launch control and prescribed time live on the host, while every authoritative
field, boundary, residual and conservation diagnostic remains CUDA-resident.

A mechanical candidate is transactional.  Both pressure buffers, classification, solid dilatation and derived
membrane source channels are snapshotted device-to-device before the fluid subcycles.  If mechanics does not
converge, or if any registered device failure counter is non-zero, the candidate geometry and all snapshotted
fields are restored on CUDA and accepted physical time does not advance.  A persistent device run-active latch
implements fail-stop multi-step control: after the first rejected transaction, every remaining pre-scheduled
outer step is disabled and rolled back without a host read.  Numerical mechanics retry chunks occur inside the
same outer transaction and therefore never masquerade as physical time.

Sanity Gate:
    * Time: ``dt_phys`` is physical seconds; inner iterations are never reinterpreted as time.
    * CFL: ``n_sub = ceil(dt_phys / dt_cfl)`` and ``dt_sub = dt_phys/n_sub``.
    * Ordering: current domain -> pressure-dependent fluid subcycles -> mechanics -> live-domain remap.
    * Conservation: one device-resident ledger integrates the true membrane surface flux, its deposited grid
      source, the independent interior source, the ``-alpha*div(v_s)`` source actually consumed by every Biot
      subcycle, and the post-mechanics moving-face transport.  The newly spread ``div(v_s)`` after remap is
      explicitly the next staggered step's input; it is never back-dated into the current content balance.
    * Rejection: a failed candidate restores every authoritative field and leaves accepted device time unchanged.
      It also deactivates all later pre-scheduled steps on device; the host never branches on acceptance.
    * Residency: ``outer_step`` and ``run`` contain no ``.numpy()``, host scalar extraction, or state-dependent
      early exit.  Post-loop report materialization is explicit.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field

import warp as wp

from aleph.components.fluid.biot_substrate import BiotSubstrate
from aleph.components.fluid.boundary import MembraneFluxBC
from aleph.components.fluid.domain import Domain

__all__ = [
    "DeviceInnerSolveReport",
    "DeviceOuterStepReport",
    "OuterStepReport",
    "PhysicalScheduler",
]


@wp.kernel
def _conditional_restore_f64_3d_kernel(
    state: wp.array3d(dtype=wp.float64),
    start: wp.array3d(dtype=wp.float64),
    accepted: wp.array(dtype=wp.int32),
) -> None:
    """Restore one float64 field when the outer candidate was rejected."""
    i, j, k = wp.tid()
    if accepted[0] == 0:
        state[i, j, k] = start[i, j, k]


@wp.kernel
def _conditional_restore_i32_3d_kernel(
    state: wp.array3d(dtype=wp.int32),
    start: wp.array3d(dtype=wp.int32),
    accepted: wp.array(dtype=wp.int32),
) -> None:
    """Restore one int32 field when the outer candidate was rejected."""
    i, j, k = wp.tid()
    if accepted[0] == 0:
        state[i, j, k] = start[i, j, k]


@wp.kernel
def _invalidate_if_counter_nonzero_kernel(
    accepted: wp.array(dtype=wp.int32),
    failure_counter: wp.array(dtype=wp.int32),
) -> None:
    """Latch a device diagnostic failure into the outer acceptance flag."""
    if failure_counter[0] != 0:
        accepted[0] = 0


@wp.kernel
def _finite_f64_3d_kernel(state: wp.array3d(dtype=wp.float64), finite: wp.array(dtype=wp.int32)) -> None:
    """Latch non-finite field state without a host readback."""
    i, j, k = wp.tid()
    if not wp.isfinite(state[i, j, k]):
        finite[0] = 0


@wp.kernel
def _invalidate_if_zero_kernel(accepted: wp.array(dtype=wp.int32), valid: wp.array(dtype=wp.int32)) -> None:
    """Fold a device validity latch into outer-step acceptance."""
    if valid[0] == 0:
        accepted[0] = 0


@wp.kernel
def _initialize_outer_acceptance_kernel(
    accepted: wp.array(dtype=wp.int32),
    step_enabled: wp.array(dtype=wp.int32),
    inner_converged: wp.array(dtype=wp.int32),
) -> None:
    """Start one outer acceptance latch from run activity and mechanics convergence."""
    accepted[0] = step_enabled[0] * inner_converged[0]


@wp.kernel
def _deactivate_run_if_rejected_kernel(
    run_active: wp.array(dtype=wp.int32), accepted: wp.array(dtype=wp.int32),
) -> None:
    """Fail-stop the remaining fixed-schedule outer loop after the first rejection."""
    if accepted[0] == 0:
        run_active[0] = 0


@wp.kernel
def _zero_if_rejected_kernel(value: wp.array(dtype=wp.float64), accepted: wp.array(dtype=wp.int32)) -> None:
    """Zero a candidate-only scalar when its outer step is rejected."""
    if accepted[0] == 0:
        value[0] = wp.float64(0.0)


@wp.kernel
def _write_rejection_kernel(accepted: wp.array(dtype=wp.int32), rejected: wp.array(dtype=wp.int32)) -> None:
    """Materialize ``1-accepted`` without a host branch."""
    rejected[0] = 1 - accepted[0]


@wp.kernel
def _advance_time_if_accepted_kernel(
    accepted_time: wp.array(dtype=wp.float64), dt_phys: wp.float64, accepted: wp.array(dtype=wp.int32),
) -> None:
    """Advance the authoritative physical clock only for an accepted coupled state."""
    if accepted[0] != 0:
        accepted_time[0] += dt_phys


@wp.kernel
def _accumulate_grid_conservation_kernel(
    mask: wp.array3d(dtype=wp.int32),
    interior_source: wp.array3d(dtype=wp.float64),
    membrane_source: wp.array3d(dtype=wp.float64),
    div_vs: wp.array3d(dtype=wp.float64),
    dt_sub: wp.float64,
    alpha: wp.float64,
    cell_volume: wp.float64,
    interior_integral: wp.array(dtype=wp.float64),
    membrane_grid_integral: wp.array(dtype=wp.float64),
    dilatation_integral: wp.array(dtype=wp.float64),
) -> None:
    """Integrate every grid source consumed by one Biot subcycle into the outer NG-2 ledger."""
    i, j, k = wp.tid()
    if mask[i, j, k] == 1:
        scale = dt_sub * cell_volume
        wp.atomic_add(interior_integral, 0, interior_source[i, j, k] * scale)
        wp.atomic_add(membrane_grid_integral, 0, membrane_source[i, j, k] * scale)
        wp.atomic_add(dilatation_integral, 0, -alpha * div_vs[i, j, k] * scale)


@wp.kernel
def _accumulate_scalar_rate_kernel(
    integral: wp.array(dtype=wp.float64), rate: wp.array(dtype=wp.float64), dt_sub: wp.float64,
) -> None:
    """Accumulate a device scalar rate over one prescribed subcycle without a host readback."""
    integral[0] += rate[0] * dt_sub


@wp.kernel
def _finalize_conservation_kernel(
    content_start: wp.array(dtype=wp.float64),
    content_end: wp.array(dtype=wp.float64),
    membrane_surface: wp.array(dtype=wp.float64),
    membrane_grid: wp.array(dtype=wp.float64),
    interior: wp.array(dtype=wp.float64),
    dilatation: wp.array(dtype=wp.float64),
    moving_face: wp.array(dtype=wp.float64),
    content_delta: wp.array(dtype=wp.float64),
    surface_deposition_error: wp.array(dtype=wp.float64),
    solver_balance_error: wp.array(dtype=wp.float64),
    physical_balance_error: wp.array(dtype=wp.float64),
) -> None:
    """Close the committed outer-step ledger against both grid and true-surface membrane fluxes."""
    delta = content_end[0] - content_start[0]
    shared = interior[0] + dilatation[0] + moving_face[0]
    content_delta[0] = delta
    surface_deposition_error[0] = membrane_grid[0] - membrane_surface[0]
    solver_balance_error[0] = delta - (membrane_grid[0] + shared)
    physical_balance_error[0] = delta - (membrane_surface[0] + shared)


@dataclass
class DeviceInnerSolveReport:
    """Device-resident mechanics diagnostics returned by the injected inner solve."""

    residual_d: wp.array
    iters_d: wp.array
    attempts_d: wp.array
    converged_d: wp.array
    max_displacement_d: wp.array
    constraint_residual_d: wp.array
    finite_d: wp.array
    dt_mu_d: wp.array


@dataclass
class OuterStepReport:
    """Materialized post-loop diagnostics from one physical step."""

    t: float
    attempted_t: float
    dt_phys: float
    n_biot_subcycles: int
    total_content: float
    content_start: float
    content_delta: float
    membrane_surface_flux_integral: float
    membrane_grid_source_integral: float
    interior_source_integral: float
    solid_dilatation_integral: float
    moving_face_delta: float
    membrane_deposition_error: float
    solver_conservation_error: float
    physical_conservation_error: float
    inner_residual: float
    inner_iters: int
    inner_attempts: int
    inner_converged: bool
    inner_max_displacement: float
    inner_constraint_residual: float
    inner_finite: bool
    inner_dt_mu: float
    outer_step_enabled: bool
    outer_accepted: bool
    outer_rolled_back: bool
    fluid_finite: bool


@dataclass
class DeviceOuterStepReport:
    """One outer-step report whose authoritative diagnostics remain on CUDA."""

    device: str
    accepted_time_d: wp.array
    attempted_t: float
    dt_phys: float
    n_biot_subcycles: int
    total_content_d: wp.array
    content_start_d: wp.array
    content_delta_d: wp.array
    membrane_surface_flux_integral_d: wp.array
    membrane_grid_source_integral_d: wp.array
    interior_source_integral_d: wp.array
    solid_dilatation_integral_d: wp.array
    moving_face_delta_d: wp.array
    membrane_deposition_error_d: wp.array
    solver_conservation_error_d: wp.array
    physical_conservation_error_d: wp.array
    inner_residual_d: wp.array
    inner_iters_d: wp.array
    inner_attempts_d: wp.array
    inner_converged_d: wp.array
    inner_max_displacement_d: wp.array
    inner_constraint_residual_d: wp.array
    inner_finite_d: wp.array
    inner_dt_mu_d: wp.array
    outer_step_enabled_d: wp.array
    outer_accepted_d: wp.array
    outer_rolled_back_d: wp.array
    fluid_finite_d: wp.array

    def readback(self) -> OuterStepReport:
        """Materialize diagnostics after the physical-time loop (the explicit NG-6 readback boundary)."""
        wp.synchronize_device(self.device)
        return OuterStepReport(
            t=float(self.accepted_time_d.numpy()[0]),
            attempted_t=self.attempted_t,
            dt_phys=self.dt_phys,
            n_biot_subcycles=self.n_biot_subcycles,
            total_content=float(self.total_content_d.numpy()[0]),
            content_start=float(self.content_start_d.numpy()[0]),
            content_delta=float(self.content_delta_d.numpy()[0]),
            membrane_surface_flux_integral=float(self.membrane_surface_flux_integral_d.numpy()[0]),
            membrane_grid_source_integral=float(self.membrane_grid_source_integral_d.numpy()[0]),
            interior_source_integral=float(self.interior_source_integral_d.numpy()[0]),
            solid_dilatation_integral=float(self.solid_dilatation_integral_d.numpy()[0]),
            moving_face_delta=float(self.moving_face_delta_d.numpy()[0]),
            membrane_deposition_error=float(self.membrane_deposition_error_d.numpy()[0]),
            solver_conservation_error=float(self.solver_conservation_error_d.numpy()[0]),
            physical_conservation_error=float(self.physical_conservation_error_d.numpy()[0]),
            inner_residual=float(self.inner_residual_d.numpy()[0]),
            inner_iters=int(self.inner_iters_d.numpy()[0]),
            inner_attempts=int(self.inner_attempts_d.numpy()[0]),
            inner_converged=bool(self.inner_converged_d.numpy()[0]),
            inner_max_displacement=float(self.inner_max_displacement_d.numpy()[0]),
            inner_constraint_residual=float(self.inner_constraint_residual_d.numpy()[0]),
            inner_finite=bool(self.inner_finite_d.numpy()[0]),
            inner_dt_mu=float(self.inner_dt_mu_d.numpy()[0]),
            outer_step_enabled=bool(self.outer_step_enabled_d.numpy()[0]),
            outer_accepted=bool(self.outer_accepted_d.numpy()[0]),
            outer_rolled_back=bool(self.outer_rolled_back_d.numpy()[0]),
            fluid_finite=bool(self.fluid_finite_d.numpy()[0]),
        )


@dataclass
class PhysicalScheduler:
    """Drive the GPU-resident outer physical clock over a Biot substrate and live domain."""

    substrate: BiotSubstrate
    domain: Domain
    membrane_bc: MembraneFluxBC | None = None
    inner_solve: Callable[[float, wp.array], DeviceInnerSolveReport] | None = None
    osmotic_difference: Callable[[float], float] = field(default=lambda t: 0.0)
    _t: float = 0.0

    def __post_init__(self) -> None:
        self.domain.bind_storage(self.substrate.storage_S)
        g = self.substrate.grid
        with wp.ScopedDevice(self.substrate.grid.device):
            self._zero_f64_d = wp.zeros(1, dtype=wp.float64)
            self._zero_i32_d = wp.zeros(1, dtype=wp.int32)
            self._accepted_time_d = wp.full(1, float(self._t), dtype=wp.float64)
            self._run_active_d = wp.ones(1, dtype=wp.int32)
            self._attempt_enabled_d = wp.ones(1, dtype=wp.int32)
            self._outer_accepted_d = wp.ones(1, dtype=wp.int32)
            self._rejected_d = wp.zeros(1, dtype=wp.int32)
            self._fluid_finite_d = wp.ones(1, dtype=wp.int32)
            self._p_start_d = wp.empty(g.shape, dtype=wp.float64)
            self._p_new_start_d = wp.empty(g.shape, dtype=wp.float64)
            self._mask_start_d = wp.empty(g.shape, dtype=wp.int32)
            self._div_start_d = wp.empty(g.shape, dtype=wp.float64)
            self._s_membrane_start_d = wp.empty(g.shape, dtype=wp.float64)
            self._s_total_start_d = wp.empty(g.shape, dtype=wp.float64)
            self._content_start_d = wp.zeros(1, dtype=wp.float64)
            self._content_end_d = wp.zeros(1, dtype=wp.float64)
            self._content_delta_d = wp.zeros(1, dtype=wp.float64)
            self._membrane_surface_integral_d = wp.zeros(1, dtype=wp.float64)
            self._membrane_grid_integral_d = wp.zeros(1, dtype=wp.float64)
            self._interior_integral_d = wp.zeros(1, dtype=wp.float64)
            self._dilatation_integral_d = wp.zeros(1, dtype=wp.float64)
            self._membrane_deposition_error_d = wp.zeros(1, dtype=wp.float64)
            self._solver_conservation_error_d = wp.zeros(1, dtype=wp.float64)
            self._physical_conservation_error_d = wp.zeros(1, dtype=wp.float64)

    @staticmethod
    def _snapshot_scalar(a: wp.array) -> wp.array:
        """Retain one scalar for a report via D2D clone (never a D2H read)."""
        return wp.clone(a)

    def outer_step(self, dt_phys: float, *, biot_cfl_safety: float = 0.9) -> DeviceOuterStepReport:
        """Advance one physical step without transferring authoritative state to the host."""
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        if not math.isfinite(biot_cfl_safety) or not 0.0 < biot_cfl_safety <= 1.0:
            raise ValueError("biot_cfl_safety must be finite and in (0, 1]")
        sub = self.substrate
        g = sub.grid

        # Fixed host launch schedule, device decision: a prior rejection leaves ``run_active=0``.  The current
        # enable bit is snapshotted D2D for reporting, and also initializes the fluid-valid latch so a disabled
        # step cannot move geometry or become accepted even though its prescribed kernels are still launched.
        wp.copy(self._attempt_enabled_d, self._run_active_d)
        wp.copy(self._fluid_finite_d, self._attempt_enabled_d)

        # Transaction start: every authoritative field that the coupled attempt may mutate is retained D2D.
        wp.copy(self._p_start_d, g.p)
        wp.copy(self._p_new_start_d, g.p_new)
        wp.copy(self._mask_start_d, g.mask)
        wp.copy(self._div_start_d, g.div_vs)
        wp.copy(self._s_membrane_start_d, g.s_membrane)
        wp.copy(self._s_total_start_d, g.s_total)
        wp.copy(self._content_start_d, g.total_content_device(sub.storage_S))
        for scalar_d in (
            self._content_delta_d,
            self._membrane_surface_integral_d,
            self._membrane_grid_integral_d,
            self._interior_integral_d,
            self._dilatation_integral_d,
            self._membrane_deposition_error_d,
            self._solver_conservation_error_d,
            self._physical_conservation_error_d,
        ):
            scalar_d.zero_()

        dt_cfl = sub.cfl_dt(biot_cfl_safety)
        n_sub = max(1, int(dt_phys / dt_cfl) + (1 if dt_phys % dt_cfl > 0 else 0))
        dt_sub = dt_phys / n_sub
        for substep in range(n_sub):
            if self.membrane_bc is not None:
                # Rebuild only the pressure-dependent membrane channel. Interior sources have an independent
                # owner and survive every subcycle.
                self.membrane_bc.apply(self.osmotic_difference(self._t + substep * dt_sub))
                wp.launch(
                    _accumulate_scalar_rate_kernel,
                    dim=1,
                    inputs=[self._membrane_surface_integral_d, self.membrane_bc.integrated_flux_d,
                            wp.float64(dt_sub)],
                    device=g.device,
                )
            wp.launch(
                _accumulate_grid_conservation_kernel,
                dim=g.shape,
                inputs=[g.mask, g.s_water, g.s_membrane, g.div_vs, wp.float64(dt_sub), wp.float64(sub.alpha),
                        wp.float64(g.dx**g.dim)],
                outputs=[self._interior_integral_d, self._membrane_grid_integral_d,
                         self._dilatation_integral_d],
                device=g.device,
            )
            sub.step(dt_sub)

        # Audit fluid inputs before mechanics.  The device latch is passed into the inner solver so every
        # position/projection update becomes a no-op when a Biot/source field is non-finite; no D2H branch is
        # needed and invalid pressure never reaches authoritative geometry.
        for field_d in (g.p, g.p_new, g.p_bar, g.div_vs, g.s_water, g.s_membrane, g.s_total):
            wp.launch(_finite_f64_3d_kernel, dim=g.shape, inputs=[field_d, self._fluid_finite_d], device=g.device)

        inner = DeviceInnerSolveReport(
            residual_d=self._zero_f64_d,
            iters_d=self._zero_i32_d,
            attempts_d=self._zero_i32_d,
            converged_d=self._attempt_enabled_d,
            max_displacement_d=self._zero_f64_d,
            constraint_residual_d=self._zero_f64_d,
            finite_d=self._attempt_enabled_d,
            dt_mu_d=self._zero_f64_d,
        )
        if self.inner_solve is not None:
            inner = self.inner_solve(dt_phys, self._fluid_finite_d)

        # Re-audit after mechanics before acceptance.  This catches any future inner hook that writes a fluid
        # channel; the pre-mechanics audit above is what prevents invalid input from moving geometry.
        for field_d in (g.p, g.p_new, g.p_bar, g.div_vs, g.s_water, g.s_membrane, g.s_total):
            wp.launch(_finite_f64_3d_kernel, dim=g.shape, inputs=[field_d, self._fluid_finite_d], device=g.device)
        wp.launch(
            _initialize_outer_acceptance_kernel,
            dim=1,
            inputs=[self._outer_accepted_d, self._attempt_enabled_d, inner.converged_d],
            device=g.device,
        )
        wp.launch(_invalidate_if_zero_kernel, dim=1,
                  inputs=[self._outer_accepted_d, self._fluid_finite_d], device=g.device)

        failure_counters = getattr(self.inner_solve, "failure_counters", ()) if self.inner_solve is not None else ()
        for counter in failure_counters:
            wp.launch(_invalidate_if_counter_nonzero_kernel, dim=1, inputs=[self._outer_accepted_d, counter],
                      device=g.device)
        rollback = getattr(self.inner_solve, "rollback", None) if self.inner_solve is not None else None
        if rollback is not None:
            rollback(self._outer_accepted_d)

        # Mechanics has just moved both closed surfaces.  Classify/remap NOW so the next fluid step sees the
        # current geometry and this step records its moving-face content (no one-step lag).
        moving_d = self.domain.remap()
        # Live-mesh queries happen during remap.  Latch any new device failure, restore the mechanical candidate,
        # then refresh the mesh acceleration structures against the now-authoritative (possibly restored) nodes.
        for counter in failure_counters:
            wp.launch(_invalidate_if_counter_nonzero_kernel, dim=1, inputs=[self._outer_accepted_d, counter],
                      device=g.device)
        if rollback is not None:
            rollback(self._outer_accepted_d)
        self.domain.classify()
        if self.inner_solve is not None:
            post_remap = getattr(self.inner_solve, "post_remap", None)
            if post_remap is not None:
                # The new membrane-minus-nucleus mask is now authoritative. Spread accepted solid motion on
                # this geometry so the next explicit staggered Biot step receives div(v_s) without a mask lag.
                post_remap(dt_phys, self._outer_accepted_d)

        for field_d in (g.p, g.p_new, g.p_bar, g.div_vs, g.s_water, g.s_membrane, g.s_total):
            wp.launch(_finite_f64_3d_kernel, dim=g.shape, inputs=[field_d, self._fluid_finite_d], device=g.device)
        wp.launch(_invalidate_if_zero_kernel, dim=1,
                  inputs=[self._outer_accepted_d, self._fluid_finite_d], device=g.device)
        if rollback is not None:
            rollback(self._outer_accepted_d)
        # If the post-remap FSI spread invalidated the transaction, return mesh accelerators to the restored
        # geometry before restoring the Eulerian mask snapshot below.
        self.domain.classify()

        # A rejected outer candidate leaves no pressure, domain, source, dilatation, geometry, time, or moving-
        # face trace behind.  These fixed launches contain only device-side predicates, hence no D2H branch.
        wp.launch(_conditional_restore_f64_3d_kernel, dim=g.shape,
                  inputs=[g.p, self._p_start_d, self._outer_accepted_d], device=g.device)
        wp.launch(_conditional_restore_f64_3d_kernel, dim=g.shape,
                  inputs=[g.p_new, self._p_new_start_d, self._outer_accepted_d], device=g.device)
        wp.launch(_conditional_restore_i32_3d_kernel, dim=g.shape,
                  inputs=[g.mask, self._mask_start_d, self._outer_accepted_d], device=g.device)
        wp.launch(_conditional_restore_f64_3d_kernel, dim=g.shape,
                  inputs=[g.div_vs, self._div_start_d, self._outer_accepted_d], device=g.device)
        wp.launch(_conditional_restore_f64_3d_kernel, dim=g.shape,
                  inputs=[g.s_membrane, self._s_membrane_start_d, self._outer_accepted_d], device=g.device)
        wp.launch(_conditional_restore_f64_3d_kernel, dim=g.shape,
                  inputs=[g.s_total, self._s_total_start_d, self._outer_accepted_d], device=g.device)
        wp.launch(_zero_if_rejected_kernel, dim=1, inputs=[moving_d, self._outer_accepted_d], device=g.device)
        for scalar_d in (
            self._membrane_surface_integral_d,
            self._membrane_grid_integral_d,
            self._interior_integral_d,
            self._dilatation_integral_d,
        ):
            wp.launch(_zero_if_rejected_kernel, dim=1,
                      inputs=[scalar_d, self._outer_accepted_d], device=g.device)
        wp.launch(_write_rejection_kernel, dim=1,
                  inputs=[self._outer_accepted_d, self._rejected_d], device=g.device)
        wp.launch(_advance_time_if_accepted_kernel, dim=1,
                  inputs=[self._accepted_time_d, wp.float64(dt_phys), self._outer_accepted_d], device=g.device)
        wp.launch(_deactivate_run_if_rejected_kernel, dim=1,
                  inputs=[self._run_active_d, self._outer_accepted_d], device=g.device)
        commit_irreversible = (
            getattr(self.inner_solve, "commit_irreversible", None) if self.inner_solve is not None else None)
        if commit_irreversible is not None:
            # Fixed host launch schedule; the irreversible update itself is predicated by the authoritative
            # device acceptance scalar.  No acceptance value crosses to Python inside the physical-time step.
            commit_irreversible(self._outer_accepted_d)
        content_d = g.total_content_device(sub.storage_S)
        wp.copy(self._content_end_d, content_d)
        wp.launch(
            _finalize_conservation_kernel,
            dim=1,
            inputs=[self._content_start_d, self._content_end_d, self._membrane_surface_integral_d,
                    self._membrane_grid_integral_d, self._interior_integral_d, self._dilatation_integral_d,
                    moving_d],
            outputs=[self._content_delta_d, self._membrane_deposition_error_d,
                     self._solver_conservation_error_d, self._physical_conservation_error_d],
            device=g.device,
        )

        self._t += dt_phys
        return DeviceOuterStepReport(
            device=g.device,
            accepted_time_d=self._snapshot_scalar(self._accepted_time_d),
            attempted_t=self._t,
            dt_phys=dt_phys,
            n_biot_subcycles=n_sub,
            total_content_d=self._snapshot_scalar(content_d),
            content_start_d=self._snapshot_scalar(self._content_start_d),
            content_delta_d=self._snapshot_scalar(self._content_delta_d),
            membrane_surface_flux_integral_d=self._snapshot_scalar(self._membrane_surface_integral_d),
            membrane_grid_source_integral_d=self._snapshot_scalar(self._membrane_grid_integral_d),
            interior_source_integral_d=self._snapshot_scalar(self._interior_integral_d),
            solid_dilatation_integral_d=self._snapshot_scalar(self._dilatation_integral_d),
            moving_face_delta_d=self._snapshot_scalar(moving_d),
            membrane_deposition_error_d=self._snapshot_scalar(self._membrane_deposition_error_d),
            solver_conservation_error_d=self._snapshot_scalar(self._solver_conservation_error_d),
            physical_conservation_error_d=self._snapshot_scalar(self._physical_conservation_error_d),
            inner_residual_d=self._snapshot_scalar(inner.residual_d),
            inner_iters_d=self._snapshot_scalar(inner.iters_d),
            inner_attempts_d=self._snapshot_scalar(inner.attempts_d),
            inner_converged_d=self._snapshot_scalar(inner.converged_d),
            inner_max_displacement_d=self._snapshot_scalar(inner.max_displacement_d),
            inner_constraint_residual_d=self._snapshot_scalar(inner.constraint_residual_d),
            inner_finite_d=self._snapshot_scalar(inner.finite_d),
            inner_dt_mu_d=self._snapshot_scalar(inner.dt_mu_d),
            outer_step_enabled_d=self._snapshot_scalar(self._attempt_enabled_d),
            outer_accepted_d=self._snapshot_scalar(self._outer_accepted_d),
            outer_rolled_back_d=self._snapshot_scalar(self._rejected_d),
            fluid_finite_d=self._snapshot_scalar(self._fluid_finite_d),
        )

    def run(self, t_ramp: float, n_phys: int, **kwargs: object) -> list[DeviceOuterStepReport]:
        """Launch a fixed multi-step schedule whose fail-stop decision remains device-resident."""
        if not math.isfinite(t_ramp) or t_ramp <= 0.0:
            raise ValueError("t_ramp must be finite and positive")
        if isinstance(n_phys, bool) or not isinstance(n_phys, int) or n_phys <= 0:
            raise ValueError("n_phys must be a positive integer")
        dt_phys = t_ramp / n_phys
        return [self.outer_step(dt_phys, **kwargs) for _ in range(n_phys)]  # type: ignore[arg-type]
