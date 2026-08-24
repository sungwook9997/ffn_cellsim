r"""The ``extracellular_medium`` component runtime — exterior Stokes drag on the cell's outer surface.

KU-3.5 (whole-cell load path).  PI decision D5-A, 2026-07-28: this is T10's implementation of the 14th
component and its ``membrane_medium_traction`` connector, which until now were CONTRACTED and nothing more.

WHAT WAS ACTUALLY BROKEN.  The cell's six whole-body rigid modes carry zero internal stiffness by
construction -- translating or rotating the entire cell costs no elastic energy -- so every implicit path
in the tree holds them with a numerical term: ``aI`` with ``a = 1/(dt/gamma)`` in
:mod:`aleph.engine.sf_implicit` and ``compute_regularization_kernel`` in ``ac/cell/``.  That term is a
mobility the CFL bound already computed; it is not a force law, it has no relation to the cell's size or
shape, and it can therefore neither be checked against a measurement nor LOAD a rigid mode.  An exterior
medium is the physically correct object in that slot, and it is the reason the plan rates the gap
BLOCKS_CRAWL rather than cosmetic: a cell cannot translate through nothing.

WHAT THIS MODULE ADDS, PRECISELY.  A device-resident regularised-Stokeslet single-layer operator over the
cell's outer surface quadrature, and its inverse applied matrix-free by conjugate gradients:

    ``u = M F``            (mobility: point forces on the surface -> surface velocities)
    ``f_medium = -M^-1 v`` (resistance: prescribed surface velocity -> force the medium exerts back)

with ``v = (x_candidate - x_committed) / dt_phys`` the candidate surface velocity.  The tangent of that
force is ``M^-1 / dt_phys``, symmetric positive definite on the WHOLE space including the six rigid modes
-- so the operator is non-singular there for a physical reason, with the sphere limits ``6 pi mu a`` and
``8 pi mu a^3`` as exact checks, rather than because a regulariser was added.  Closed forms and the
grid-invariance instrument live in :mod:`aleph.engine.medium_stokes_analytic`.

WHY NOT DARCY, AND WHY NOT A GRID.  Outside the cell ``k -> inf``, so the Brinkman screening length
``l_B = sqrt(k) -> inf`` and Darcy is wrong there rather than approximate.  What is left is Stokes, and the
single-layer integral representation IS the exterior Stokes solution -- the only approximations are the
surface quadrature and the blob width, both refinable and both measured by
``medium_stokes_analytic.sphere_drag_convergence``.  A grid solve would add a truncated far field and an
immersed-boundary blob on top, and at any affordable ``dx`` could not resolve the cell--substrate
lubrication film that crawl traction lives in.

TRAP #4 IS NOT OPTIONAL, AND THE AUDIT CHANGED WHAT IT MEANS HERE.  The plan requires that when a medium
starts dragging, the per-node drag ``gamma_node = 6 pi eta R / Nc`` be removed in the SAME change, because
cytoplasm viscosity applied in both particle drag and a background fluid is a 2x error.  The measured
situation (``docs/v2_audit/T10_GAMMA_NODE_AUDIT_2026-07-28.md``) is that no such term has ever existed in
``ac/``: both live sites are in ``ff/`` (``network_warp.py``'s whole-cell compression driver and
``implicit_ff.py``'s modal COM drag) and their only callers are four legacy ``scripts/ff_*.py`` drivers.
So the obligation for this lane is not a deletion but an INVARIANT, and
:func:`assert_single_dissipation_owner` is the machine-checkable form of it: the exterior medium must be
the only velocity-proportional drag any component addressing this surface applies.  A guard is what keeps
the 2x error from arriving later by import, which a one-time deletion would not.

WHAT THIS MODULE DELIBERATELY DOES NOT DO.  (1) No wall.  ``MediumWallMode.HALF_SPACE_BLAKE`` is declared
and raises: the basal no-slip image system is exact and closed-form, but it is the asymmetric world
boundary (T10 checklist item 3) and is a separate landing, so an unimplemented mode fails loudly instead of
silently returning free-space drag under a wall's name.  (2) No DISPATCH claim -- the
``membrane_medium_traction`` canonical facade claim belongs to
:class:`~aleph.engine.surface_body.SurfaceBody`, whose lane owns the membrane surface quadrature.  The
connector RUNTIME that claim binds is :class:`MembraneMediumTraction`, at the foot of this module
(added 2026-08-09): the component here computes the traction, and that facade is the edge identity it
is bound under.  (3) No magnitude.  ``mu_medium`` has no KnowledgeClaim and no SourceEvidence in the knowledge base (queried
2026-07-28) and is therefore a PI-GAP: the settings object REQUIRES it with no default, and every gate is
a shape/exactness gate that holds at any positive viscosity.

engine units: length um, force pN, velocity um/s, viscosity Pa*s (numerically pN*s/um^2, ``ff/units.py``).
Runtime: NVIDIA Warp on CUDA only (I0-A).  The module is CPU-importable -- kernels JIT lazily -- while
:func:`build_exterior_stokes_medium` allocates device memory and is a gbook lane.

Sanity Gate (self-tested in ``tests/ac/engine/test_medium_exterior.py``):
    * dimensional: ``[M] = um/(pN*s)``, so ``[M^-1 v] = pN`` and the translation resistance carries
      pN*s/um, matching ``6 pi mu a``.
    * boundary cases: zero surface velocity gives exactly zero medium force (no spurious source); a
      single quadrature point reduces to the isolated self-mobility ``1/(4 pi mu eps)``.
    * conservation/invariant: the assembled mobility is EXACTLY symmetric (the kernel is even in the
      separation vector, so the two blocks are the same IEEE products) and positive definite, hence the
      six-mode grand resistance is SPD -- the rigid modes are loaded, not masked.
    * sign-sense: the medium force opposes surface motion, so ``sum(v . f_medium) < 0`` strictly for any
      non-zero velocity field; a positive hydrodynamic power is a sign error, never a physiological one.
    * numerical: float64 throughout; CG residual reduction is reported, never asserted from a host read
      inside the loop, and the operator's condition number is measured rather than assumed.
    * measurement-protocol: the resistance is reduced about a stated centre, and the reported drag is a
      property of the surface quadrature at a stated ``epsilon/h`` -- refining the surface is the lever,
      not lowering ``epsilon``.
    * grid-invariance (Magic-Number Block): ``epsilon`` is the mesh's own mean nearest-neighbour spacing,
      so it falls under refinement and is never selected against an outcome; the measured convergence
      order in ``h`` is the evidence.
"""

from __future__ import annotations

import ast
import math
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import numpy as np
import warp as wp

# Generic device-scalar CG arithmetic, reused rather than rebuilt.  These kernels carry no cortex physics:
# they are the preconditioned-CG scalar recurrences over `wp.vec3d` arrays, and the operator is supplied by
# the caller.  Reusing them keeps one CG convergence/finiteness latch semantics across the engine.
from aleph.components.incumbent.implicit_mechanics import (
    _cg_alpha_kernel,
    _cg_begin_kernel,
    _cg_beta_kernel,
    _cg_direction_kernel,
    _cg_finalize_kernel,
    _cg_update_kernel,
    _dot_kernel,
)
from aleph.engine.contracts import ConnectorFamily, reference_cell_architecture
from aleph.engine.ledger import GlobalCellLedger
from aleph.engine.medium_stokes_analytic import mean_nearest_neighbour_spacing

__all__ = [
    "EXTERIOR_MEDIUM_COMPONENT",
    "MEMBRANE_MEDIUM_CONNECTOR",
    "ExteriorStokesMedium",
    "ExteriorStokesMediumSettings",
    "MediumWallMode",
    "MembraneMediumTraction",
    "assert_single_dissipation_owner",
    "build_exterior_stokes_medium",
    "derive_blob_epsilon",
]

EXTERIOR_MEDIUM_COMPONENT = "extracellular_medium"
MEMBRANE_MEDIUM_CONNECTOR = "membrane_medium_traction"

#: Sites carrying a per-node velocity-proportional drag that must never coexist with the exterior medium
#: on one surface (PI framework trap #4).  Measured live sites as of 2026-07-28; see the audit document.
FORBIDDEN_DRAG_SYMBOLS: tuple[str, ...] = (
    "physical_node_gammas",
    "simulate_whole_cell_compression_on_device",
    "implicit_step_current",
    "ff_implicit_step_gpu",
)


class MediumWallMode(StrEnum):
    """Which exterior Green's function the medium uses.

    ``FREE_SPACE`` is the unbounded Oseen/Cortez kernel.  ``HALF_SPACE_BLAKE`` is the no-slip planar-wall
    image system that the asymmetric world boundary (basal 2D collagen, free face media) needs; it is
    declared so the seam exists and raises so it cannot silently return free-space drag.
    """

    FREE_SPACE = "free_space"
    HALF_SPACE_BLAKE = "half_space_blake"


@wp.kernel
def _regularized_stokeslet_matvec_kernel(
    position: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),
    epsilon: wp.float64,
    inv_eight_pi_mu: wp.float64,
    velocity: wp.array(dtype=wp.vec3d),
) -> None:
    r"""Apply the free-space regularised-Stokeslet mobility, ``velocity = M force``.

    One thread per target quadrature point, summing over every source point including itself (the
    regularised kernel is finite at zero separation, which is why no self-term special case appears).
    The evaluation is a direct O(N^2) sum with no neighbour list: the exterior kernel decays as ``1/r``
    and truncating it would silently break both the Stokes far field and the operator's symmetry.
    """
    m = wp.tid()
    eps2 = epsilon * epsilon
    target = position[m]
    accumulated = wp.vec3d(0.0, 0.0, 0.0)
    for n in range(force.shape[0]):
        dx = target - position[n]
        f = force[n]
        r2 = wp.dot(dx, dx)
        soft = r2 + eps2
        denom = soft * wp.sqrt(soft)
        isotropic = (r2 + wp.float64(2.0) * eps2) / denom
        projected = wp.dot(dx, f) / denom
        accumulated = accumulated + isotropic * f + projected * dx
    velocity[m] = inv_eight_pi_mu * accumulated


@wp.kernel
def _gather_positions_kernel(
    position: wp.array(dtype=wp.vec3d),
    surface_index: wp.array(dtype=wp.int32),
    gathered: wp.array(dtype=wp.vec3d),
) -> None:
    """Copy the addressed surface positions into the medium's own contiguous quadrature array."""
    t = wp.tid()
    gathered[t] = position[surface_index[t]]


@wp.kernel
def _surface_velocity_kernel(
    position: wp.array(dtype=wp.vec3d),
    committed_position: wp.array(dtype=wp.vec3d),
    surface_index: wp.array(dtype=wp.int32),
    inv_dt_phys: wp.array(dtype=wp.float64),
    velocity: wp.array(dtype=wp.vec3d),
) -> None:
    """Gather the candidate surface velocity ``(x_candidate - x_committed)/dt_phys``.

    ``inv_dt_phys`` is a device scalar written by the scheduler-facing setter, not a host read of
    authoritative state: ``dt_phys`` is the outer clock the scheduler already owns on the host.
    """
    t = wp.tid()
    k = surface_index[t]
    velocity[t] = inv_dt_phys[0] * (position[k] - committed_position[t])


@wp.kernel
def _scatter_medium_force_kernel(
    solved_force: wp.array(dtype=wp.vec3d),
    surface_index: wp.array(dtype=wp.int32),
    force: wp.array(dtype=wp.vec3d),
    traction: wp.array(dtype=wp.vec3d),
    reaction: wp.array(dtype=wp.vec3d),
) -> None:
    r"""Scatter ``f_medium = -M^-1 v`` onto the surface and record its adjoint reaction.

    ``solved_force`` holds ``M^-1 v`` -- the force the FLUID would have to receive to move at the surface
    velocity -- so the force on the SOLID is its negative.  The equal-and-opposite half is recorded
    explicitly in ``reaction`` rather than dropped: the medium is an unbounded reservoir with no node to
    receive it, so without this the connector would be one-sided and Newton's third law would hold only by
    assertion.  No work is accumulated here, deliberately -- inner mechanical iterations are not physical
    time, so a per-candidate work sum would integrate the solver's path instead of the step's.
    """
    t = wp.tid()
    k = surface_index[t]
    f_medium = -solved_force[t]
    traction[t] = f_medium
    reaction[t] = solved_force[t]
    wp.atomic_add(force, k, f_medium)


@wp.kernel
def _medium_rollback_kernel(
    accepted: wp.array(dtype=wp.int32),
    committed_traction: wp.array(dtype=wp.vec3d),
    traction: wp.array(dtype=wp.vec3d),
) -> None:
    """Restore the committed traction under a rejected outer step; accepted state is never mutated."""
    t = wp.tid()
    if accepted[0] == 0:
        traction[t] = committed_traction[t]


@wp.kernel
def _medium_commit_kernel(
    accepted: wp.array(dtype=wp.int32),
    quadrature_position: wp.array(dtype=wp.vec3d),
    traction: wp.array(dtype=wp.vec3d),
    committed_position: wp.array(dtype=wp.vec3d),
    committed_traction: wp.array(dtype=wp.vec3d),
    dissipated_work: wp.array(dtype=wp.float64),
) -> None:
    r"""Advance the committed reference and traction, and bank the step's dissipated work, on acceptance.

    The committed position is what the next step's velocity is measured against, so committing it is
    exactly the act of advancing the medium's history; a rejected candidate leaves it untouched and the
    re-attempted step therefore sees an identical velocity.  Work is banked HERE and only here, as
    ``f_medium . (x_accepted - x_committed)`` in pN*um -- reaction-conjugate work over the accepted
    interval, the same protocol the ECM far-field anchor uses, and dimensionally a work rather than the
    power a per-iteration sum would have produced.  It is strictly negative for any real motion, which is
    the sign-sense gate: a medium that adds energy is a sign error.  No private clock is advanced.
    """
    t = wp.tid()
    if accepted[0] != 0:
        current = quadrature_position[t]
        wp.atomic_add(dissipated_work, 0, wp.dot(traction[t], current - committed_position[t]))
        committed_position[t] = current
        committed_traction[t] = traction[t]


@dataclass(frozen=True, slots=True)
class ExteriorStokesMediumSettings:
    """Physical and numerical settings of the exterior medium.

    Attributes:
        viscosity_pa_s: Medium dynamic viscosity in Pa*s.  REQUIRED with no default: this is a PI-GAP
            (no KnowledgeClaim, no SourceEvidence as of 2026-07-28), and a default would be a convenience
            value entering production silently, which the physiological-baseline rule forbids.
        blob_epsilon_um: Regularisation length in um.  Derive it from the surface with
            :func:`derive_blob_epsilon` rather than choosing it; see the module Magic-Number note.
        wall_mode: Which exterior Green's function to use.
        cg_max_iterations: Iteration cap for the matrix-free resistance solve.
        cg_relative_tolerance: Relative residual reduction the CG latch reports convergence at.
    """

    viscosity_pa_s: float
    blob_epsilon_um: float
    wall_mode: MediumWallMode = MediumWallMode.FREE_SPACE
    cg_max_iterations: int = 200
    cg_relative_tolerance: float = 1.0e-10

    def __post_init__(self) -> None:
        for label, value in (
            ("viscosity_pa_s", self.viscosity_pa_s),
            ("blob_epsilon_um", self.blob_epsilon_um),
            ("cg_relative_tolerance", self.cg_relative_tolerance),
        ):
            if not (isinstance(value, float) and math.isfinite(value) and value > 0.0):
                raise ValueError(f"{label} must be a positive finite float, got {value!r}")
        if isinstance(self.cg_max_iterations, bool) or not isinstance(self.cg_max_iterations, int):
            raise TypeError("cg_max_iterations must be an int")
        if self.cg_max_iterations < 1:
            raise ValueError("cg_max_iterations must be at least one")
        if self.wall_mode is MediumWallMode.HALF_SPACE_BLAKE:
            raise NotImplementedError(
                "HALF_SPACE_BLAKE is declared but not implemented: the basal no-slip image system is "
                "T10 checklist item 3 (the asymmetric world boundary) and is a separate landing. "
                "Returning free-space drag under a wall's name would be a silent 2x-class error."
            )
        if self.wall_mode is not MediumWallMode.FREE_SPACE:
            raise ValueError(f"unsupported wall mode {self.wall_mode!r}")

    @property
    def inv_eight_pi_viscosity(self) -> float:
        """The kernel prefactor ``1/(8 pi mu)`` in ``um/(pN*s)``."""
        return 1.0 / (8.0 * math.pi * self.viscosity_pa_s)


def derive_blob_epsilon(positions: np.ndarray, *, epsilon_ratio: float = 1.0) -> float:
    """Return the DERIVED blob length: ``epsilon_ratio`` times the quadrature's own spacing.

    The default ratio of one means the blob width equals the point spacing, so the discrete single layer
    exactly resolves its own blob.  Holding the ratio fixed while refining the surface is what makes the
    scheme grid-invariant, and is why this is a derivation rather than a tuning knob: the ratio is never
    adjusted to move a result, and ``medium_stokes_analytic.sphere_drag_convergence`` reports the order in
    ``h`` that justifies it.

    Args:
        positions: Host copy of the surface quadrature, shape ``(N, 3)``, in um.
        epsilon_ratio: Fixed ``epsilon / h``.

    Returns:
        Blob regularisation length in um.
    """
    if not (isinstance(epsilon_ratio, float) and math.isfinite(epsilon_ratio) and epsilon_ratio > 0.0):
        raise ValueError(f"epsilon_ratio must be a positive finite float, got {epsilon_ratio!r}")
    return epsilon_ratio * mean_nearest_neighbour_spacing(positions)


def assert_single_dissipation_owner(
    module_paths: tuple[Path, ...] | list[Path],
    *,
    forbidden: tuple[str, ...] = FORBIDDEN_DRAG_SYMBOLS,
) -> None:
    """Enforce PI framework trap #4 by AST: no medium-bearing module may also import a per-node drag.

    The obligation the plan states as "remove ``gamma_node`` in the same change" is, for this lane, an
    invariant rather than a deletion: the audit found no such term in ``ac/`` at all, so what has to be
    prevented is its ARRIVAL alongside the medium.  A guard does that permanently; a one-time deletion
    would not.

    The check is an AST walk over import statements, not a source regex, so a mention inside a docstring
    or a comment -- of which this module deliberately has several -- cannot trip it, and an aliased import
    cannot evade it.

    Args:
        module_paths: Python source files to inspect.
        forbidden: Symbol names that carry a per-node or whole-cell velocity-proportional drag.

    Raises:
        ValueError: If any inspected module imports one of the forbidden drag symbols.
    """
    offences: list[str] = []
    for path in module_paths:
        tree = ast.parse(Path(path).read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                for alias in node.names:
                    if alias.name in forbidden:
                        offences.append(f"{path}:{node.lineno} imports {alias.name}")
            elif isinstance(node, ast.Attribute) and node.attr in forbidden:
                offences.append(f"{path}:{node.lineno} references attribute {node.attr}")
    if offences:
        raise ValueError(
            "PI framework trap #4: a module bound to the exterior medium also reaches a per-node drag, "
            "which double-counts viscosity (a 2x error, not a modelling choice): " + "; ".join(offences)
        )


def _device_is_cuda(array: object) -> bool:
    """Return whether an array-like object declares a CUDA device."""
    return bool(getattr(getattr(array, "device", None), "is_cuda", False))


def _validate_device_array(
    array: object, *, label: str, dtype: object, shape: tuple[int, ...]
) -> None:
    """Validate dtype, shape and CUDA residency of one authoritative medium array."""
    if getattr(array, "dtype", None) != dtype:
        raise TypeError(f"{label} must have dtype {dtype}")
    if tuple(getattr(array, "shape", ())) != shape:
        raise ValueError(f"{label} must have shape {shape}")
    if not _device_is_cuda(array):
        raise ValueError(f"{label} must be a Warp CUDA device array (I0-A: GPU is the only runtime)")


@dataclass(frozen=True, slots=True)
class ExteriorStokesMedium:
    """The exterior-medium component and the CUDA arrays it exclusively owns.

    The component owns no geometry (``owns_geometry=False`` in its contract): it binds the surface
    quadrature the membrane owns, addressing it through ``surface_index_d``.  What it does own is the
    traction it applies, the committed surface reference its velocity is measured against, the adjoint
    reaction the unbounded reservoir has no node to receive, the dissipated power, and the CG scratch.

    Every method is device-side.  Nothing here reads an authoritative device value on the host, advances a
    private physical clock, or mutates committed state outside the scheduler-owned acceptance scalar; the
    two host-readable properties are labelled out-of-loop reporting only.
    """

    settings: ExteriorStokesMediumSettings
    surface_index_d: wp.array
    committed_position_d: wp.array
    traction_d: wp.array
    committed_traction_d: wp.array
    reaction_d: wp.array
    dissipated_work_d: wp.array
    inv_dt_phys_d: wp.array
    velocity_d: wp.array
    quadrature_position_d: wp.array
    solved_force_d: wp.array
    _residual_d: wp.array
    _direction_d: wp.array
    _action_d: wp.array
    _rr_d: wp.array
    _rr_initial_d: wp.array
    _rz_d: wp.array
    _rz_old_d: wp.array
    _p_ap_d: wp.array
    _alpha_d: wp.array
    _beta_d: wp.array
    _active_d: wp.array
    _converged_d: wp.array
    _finite_d: wp.array
    _iterations_d: wp.array

    def __post_init__(self) -> None:
        n = int(self.surface_index_d.shape[0])
        if n < 1:
            raise ValueError("an exterior medium needs at least one surface quadrature point")
        _validate_device_array(
            self.surface_index_d, label="medium.surface_index_d", dtype=wp.int32, shape=(n,)
        )
        vector_arrays = (
            "committed_position_d",
            "traction_d",
            "committed_traction_d",
            "reaction_d",
            "velocity_d",
            "quadrature_position_d",
            "solved_force_d",
            "_residual_d",
            "_direction_d",
            "_action_d",
        )
        for label in vector_arrays:
            _validate_device_array(
                getattr(self, label), label=f"medium.{label}", dtype=wp.vec3d, shape=(n,)
            )
        scalar_arrays = (
            "dissipated_work_d",
            "inv_dt_phys_d",
            "_rr_d",
            "_rr_initial_d",
            "_rz_d",
            "_rz_old_d",
            "_p_ap_d",
            "_alpha_d",
            "_beta_d",
        )
        for label in scalar_arrays:
            _validate_device_array(
                getattr(self, label), label=f"medium.{label}", dtype=wp.float64, shape=(1,)
            )
        for label in ("_active_d", "_converged_d", "_finite_d", "_iterations_d"):
            _validate_device_array(
                getattr(self, label), label=f"medium.{label}", dtype=wp.int32, shape=(1,)
            )
        owned = tuple(
            getattr(self, label)
            for label in vector_arrays + scalar_arrays + ("surface_index_d", "_active_d",
                                                          "_converged_d", "_finite_d", "_iterations_d")
        )
        devices = {str(array.device) for array in owned}
        if len(devices) != 1:
            raise ValueError(f"every exterior-medium array must share one CUDA device, got {devices}")
        pointers = [int(array.ptr) for array in owned if getattr(array, "ptr", None) is not None]
        if len(set(pointers)) != len(pointers):
            raise ValueError("exterior-medium arrays must not alias one another")

    @property
    def device(self) -> str:
        """CUDA device shared by every medium array."""
        return str(self.surface_index_d.device)

    @property
    def n_quadrature(self) -> int:
        """Number of surface quadrature points, from host-visible metadata only."""
        return int(self.surface_index_d.shape[0])

    def set_step(self, dt_phys: float) -> None:
        """Publish the outer clock's step to the device as ``1/dt_phys``.

        ``dt_phys`` is the scheduler's host-owned physical step, so writing it to the device is a
        configuration push and not a readback of authoritative state.

        Args:
            dt_phys: Outer physical step in seconds; must be positive and finite.

        Raises:
            ValueError: If ``dt_phys`` is not positive and finite.
        """
        if not (math.isfinite(dt_phys) and dt_phys > 0.0):
            raise ValueError(f"dt_phys must be positive-finite, got {dt_phys!r}")
        self.inv_dt_phys_d.fill_(wp.float64(1.0 / float(dt_phys)))

    def mobility_apply(self, force_d: wp.array, velocity_d: wp.array) -> None:
        """Apply the mobility ``velocity = M force`` matrix-free on the device.

        Args:
            force_d: Point forces on the quadrature, ``wp.vec3d`` of length ``n_quadrature``, in pN.
            velocity_d: Destination surface velocities in um/s; may not alias ``force_d``.
        """
        if getattr(force_d, "ptr", None) is not None and force_d.ptr == velocity_d.ptr:
            raise ValueError("mobility_apply needs distinct input and output arrays")
        wp.launch(
            _regularized_stokeslet_matvec_kernel,
            dim=self.n_quadrature,
            inputs=[
                self.quadrature_position_d,
                force_d,
                wp.float64(self.settings.blob_epsilon_um),
                wp.float64(self.settings.inv_eight_pi_viscosity),
                velocity_d,
            ],
            device=self.device,
        )

    def _solve_resistance(self) -> None:
        r"""Solve ``M F = v`` for the point forces by matrix-free CG, leaving ``F`` in ``solved_force_d``.

        ``M`` is symmetric positive definite, so unpreconditioned CG is the right method and its
        convergence latch is the same device-scalar recurrence the rest of the engine uses.  The mobility
        has a constant diagonal ``1/(4 pi mu eps)`` -- every quadrature point has the same self-mobility --
        so a Jacobi preconditioner is exactly a uniform rescaling and would not change the iteration count;
        that is why none is applied, rather than as an omission.
        """
        d = self.device
        n = self.n_quadrature
        self.solved_force_d.zero_()
        wp.copy(self._residual_d, self.velocity_d)
        wp.copy(self._direction_d, self.velocity_d)
        self._finite_d.fill_(wp.int32(1))
        self._rr_d.zero_()
        wp.launch(_dot_kernel, dim=n, inputs=[self._residual_d, self._residual_d, self._rr_d], device=d)
        wp.launch(
            _cg_begin_kernel,
            dim=1,
            inputs=[
                self._rr_d,
                self._rr_d,
                self._rr_initial_d,
                self._rz_old_d,
                self._active_d,
                self._converged_d,
                self._finite_d,
                self._iterations_d,
            ],
            device=d,
        )
        eps64 = wp.float64(self.settings.cg_relative_tolerance ** 2)
        for iteration in range(1, self.settings.cg_max_iterations + 1):
            self.mobility_apply(self._direction_d, self._action_d)
            self._p_ap_d.zero_()
            wp.launch(_dot_kernel, dim=n, inputs=[self._direction_d, self._action_d, self._p_ap_d], device=d)
            wp.launch(
                _cg_alpha_kernel,
                dim=1,
                inputs=[self._rz_old_d, self._p_ap_d, self._active_d, self._finite_d, self._alpha_d],
                device=d,
            )
            wp.launch(
                _cg_update_kernel,
                dim=n,
                inputs=[
                    self.solved_force_d,
                    self._residual_d,
                    self._direction_d,
                    self._action_d,
                    self._alpha_d,
                    self._active_d,
                ],
                device=d,
            )
            self._rz_d.zero_()
            wp.launch(_dot_kernel, dim=n, inputs=[self._residual_d, self._residual_d, self._rz_d], device=d)
            wp.launch(
                _cg_finalize_kernel,
                dim=1,
                inputs=[
                    self._rz_d,
                    self._rr_initial_d,
                    eps64,
                    wp.int32(iteration),
                    self._active_d,
                    self._converged_d,
                    self._finite_d,
                    self._iterations_d,
                ],
                device=d,
            )
            wp.launch(
                _cg_beta_kernel,
                dim=1,
                inputs=[self._rz_d, self._rz_old_d, self._active_d, self._finite_d, self._beta_d],
                device=d,
            )
            wp.launch(
                _cg_direction_kernel,
                dim=n,
                inputs=[self._residual_d, self._beta_d, self._active_d, self._direction_d],
                device=d,
            )

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        r"""Add the exterior medium's traction to the candidate surface force.

        Gathers the candidate surface velocity against the committed reference, solves the resistance
        problem, and scatters ``f_medium = -M^-1 v`` onto the caller's force accumulator while recording
        the adjoint reaction and the dissipated power.  Also refreshes the medium's own copy of the
        quadrature positions from the caller's live array, so the kernel evaluates the Green's function on
        the CURRENT surface rather than a stale one -- the surface deforms, and an exterior operator built
        once would silently freeze the cell's shape.

        Args:
            pos: The membrane/surface owner's live position array, ``wp.vec3d``.
            force: The owner's force accumulator, ``wp.vec3d``, same length as ``pos``.
        """
        d = self.device
        n = self.n_quadrature
        wp.launch(
            _gather_positions_kernel,
            dim=n,
            inputs=[pos, self.surface_index_d, self.quadrature_position_d],
            device=d,
        )
        wp.launch(
            _surface_velocity_kernel,
            dim=n,
            inputs=[pos, self.committed_position_d, self.surface_index_d, self.inv_dt_phys_d,
                    self.velocity_d],
            device=d,
        )
        self._solve_resistance()
        wp.launch(
            _scatter_medium_force_kernel,
            dim=n,
            inputs=[
                self.solved_force_d,
                self.surface_index_d,
                force,
                self.traction_d,
                self.reaction_d,
            ],
            device=d,
        )

    def snapshot_candidate(self) -> None:
        """Reset the candidate traction to its last accepted value by a device copy.

        ``dissipated_work_d`` is deliberately NOT reset: it is a cumulative accepted-step ledger, banked
        only inside :meth:`commit_irreversible`, so a rejected candidate never touched it.
        """
        wp.copy(self.traction_d, self.committed_traction_d)

    def rollback(self, accepted: wp.array) -> None:
        """Restore the committed traction under a rejected outer-step predicate."""
        wp.launch(
            _medium_rollback_kernel,
            dim=self.n_quadrature,
            inputs=[accepted, self.committed_traction_d, self.traction_d],
            device=self.device,
        )

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Advance the committed surface reference and traction only under the accepted predicate.

        Args:
            accepted: Device acceptance scalar owned by the scheduler.
            dt_phys: Outer physical step in seconds, for validation only; no private clock is advanced.
            rng_seed: Run seed; the medium is deterministic and draws no random numbers.
        """
        if not (math.isfinite(dt_phys) and dt_phys > 0.0):
            raise ValueError(f"dt_phys must be positive-finite, got {dt_phys!r}")
        if isinstance(rng_seed, bool) or not isinstance(rng_seed, int) or rng_seed < 0:
            raise ValueError("rng_seed must be a nonnegative integer")
        wp.launch(
            _medium_commit_kernel,
            dim=self.n_quadrature,
            inputs=[
                accepted,
                self.quadrature_position_d,
                self.traction_d,
                self.committed_position_d,
                self.committed_traction_d,
                self.dissipated_work_d,
            ],
            device=self.device,
        )

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward the medium reaction resultant and its dissipated work to the global cell ledger.

        The reaction is the reservoir's half of the ``membrane_medium_traction`` adjoint pair, so it enters
        the same boundary-reaction channel the ECM far-field anchor uses: an exterior boundary receives a
        resultant that the global force balance must see rather than discard.  A ledger sink lacking the
        hook is tolerated as a schema gap, and no authoritative device value is read either way.
        """
        if isinstance(ledger, GlobalCellLedger):
            ledger.add_far_field_reaction(self.reaction_d, self.dissipated_work_d)
            return
        recorder = getattr(ledger, "add_far_field_reaction", None)
        if callable(recorder):
            recorder(self.reaction_d, self.dissipated_work_d)

    def converged_flag(self) -> bool:
        """Host read of the CG convergence latch, for OUT-OF-LOOP reporting only."""
        return bool(int(self._converged_d.numpy()[0]))

    def iteration_count(self) -> int:
        """Host read of the CG iteration count reached, for OUT-OF-LOOP reporting only."""
        return int(self._iterations_d.numpy()[0])


def build_exterior_stokes_medium(
    *,
    settings: ExteriorStokesMediumSettings,
    surface_index: np.ndarray,
    initial_position: np.ndarray,
    device: str,
) -> ExteriorStokesMedium:
    """Allocate the exterior medium's device state against an existing surface quadrature.

    Args:
        settings: Validated physical/numerical settings.
        surface_index: Indices into the surface owner's node array, shape ``(N,)``.
        initial_position: Committed surface reference positions, shape ``(N, 3)``, in um.  This is the
            state the first step's velocity is measured against, so it must be the accepted resting
            configuration and not an arbitrary zero.
        device: Warp CUDA device string; a CPU device is rejected by I0-A.

    Returns:
        The allocated :class:`ExteriorStokesMedium`.

    Raises:
        ValueError: If the index/position shapes disagree or the device is not CUDA.
    """
    index = np.ascontiguousarray(np.asarray(surface_index, dtype=np.int32).reshape(-1))
    reference = np.ascontiguousarray(np.asarray(initial_position, dtype=np.float64))
    if reference.ndim != 2 or reference.shape != (index.shape[0], 3):
        raise ValueError("initial_position must have shape (N, 3) matching surface_index")
    if index.shape[0] < 1:
        raise ValueError("an exterior medium needs at least one surface quadrature point")
    if np.any(index < 0):
        raise ValueError("surface_index entries must be non-negative node indices")
    if np.unique(index).shape[0] != index.shape[0]:
        raise ValueError("surface_index must address each surface node at most once")
    n = int(index.shape[0])

    def zeros_vec() -> wp.array:
        return wp.zeros(n, dtype=wp.vec3d, device=device)

    def zeros_scalar(dtype: object) -> wp.array:
        return wp.zeros(1, dtype=dtype, device=device)

    return ExteriorStokesMedium(
        settings=settings,
        surface_index_d=wp.array(index, dtype=wp.int32, device=device),
        committed_position_d=wp.array(reference, dtype=wp.vec3d, device=device),
        traction_d=zeros_vec(),
        committed_traction_d=zeros_vec(),
        reaction_d=zeros_vec(),
        dissipated_work_d=zeros_scalar(wp.float64),
        inv_dt_phys_d=zeros_scalar(wp.float64),
        velocity_d=zeros_vec(),
        quadrature_position_d=wp.array(reference, dtype=wp.vec3d, device=device),
        solved_force_d=zeros_vec(),
        _residual_d=zeros_vec(),
        _direction_d=zeros_vec(),
        _action_d=zeros_vec(),
        _rr_d=zeros_scalar(wp.float64),
        _rr_initial_d=zeros_scalar(wp.float64),
        _rz_d=zeros_scalar(wp.float64),
        _rz_old_d=zeros_scalar(wp.float64),
        _p_ap_d=zeros_scalar(wp.float64),
        _alpha_d=zeros_scalar(wp.float64),
        _beta_d=zeros_scalar(wp.float64),
        _active_d=zeros_scalar(wp.int32),
        _converged_d=zeros_scalar(wp.int32),
        _finite_d=zeros_scalar(wp.int32),
        _iterations_d=zeros_scalar(wp.int32),
    )


@dataclass(frozen=True, slots=True)
class MembraneMediumTraction:
    """The ``membrane_medium_traction`` CONNECTOR — the edge, not the medium that computes it.

    WHAT THIS CLOSES, and what it does not.  ``ExteriorStokesMedium`` above owns the exterior solve, the
    committed surface reference and the reaction reservoir, and it already exposes ``accumulate`` plus the
    three transaction hooks and ``accumulate_ledger``.  What it does NOT carry is a connector identity —
    no ``name``, no ``component_a``/``component_b`` — because it is the 14th COMPONENT.  So the 36th
    connector stayed unbindable for want of an identity, not for want of physics, and
    ``dispatch.py`` said of the claim that it "does not assert that any medium traction is evaluated".
    This facade is what makes that sentence false: bound here, the traction is evaluated.

    WHY A NON-MOVING FRAME STILL NEEDS AN ADJOINT (PI, 2026-08-09, ratified for ``ecm_far_field_anchor``
    and applied identically here).  ``extracellular_medium`` declares ``dynamically_evolving=False``,
    which governs POSITION: the far field does not move.  ``adjoint_transfer_required=True`` governs
    FORCE: the reaction is ACCUMULATED into the reservoir rather than discarded.  A clamp would discard
    it, and the balance ledger's two channels would then agree by OMISSION — which is precisely the
    failure a one-sided scatter is supposed to be caught by.  ``accumulate_ledger`` forwards the medium
    reaction, which the medium itself already calls "the reservoir's half of the adjoint pair".

    Structurally the twin of :class:`~aleph.engine.ecm_world.ECMBoundaryAnchorFacade`, and deliberately
    so: two environment boundaries with two different bookkeeping conventions is how one of them ends up
    silently discarding its reaction.
    """

    delegate: object
    name: str = MEMBRANE_MEDIUM_CONNECTOR
    component_a: str = "membrane"
    component_b: str = EXTERIOR_MEDIUM_COMPONENT
    bidirectional: bool = True
    adjoint_transfer_required: bool = True

    def __post_init__(self) -> None:
        contract = next(
            (
                c
                for c in reference_cell_architecture().connectors
                if c.name == MEMBRANE_MEDIUM_CONNECTOR
            ),
            None,
        )
        if contract is None:
            raise ValueError(f"architecture must declare {MEMBRANE_MEDIUM_CONNECTOR!r}")
        if contract.family is not ConnectorFamily.ENVIRONMENT_BOUNDARY:
            raise ValueError("membrane_medium_traction must use ENVIRONMENT_BOUNDARY")
        if {self.component_a, self.component_b} != {contract.component_a, contract.component_b}:
            raise ValueError(
                "membrane medium traction endpoints must be membrane and extracellular_medium"
            )
        if contract.kinetics or contract.commit_on_accept:
            raise ValueError("membrane_medium_traction is a non-kinetic boundary constraint")
        if not (self.bidirectional and self.adjoint_transfer_required):
            raise ValueError("membrane medium traction must expose reaction and adjoint work")
        # Duck-typed rather than `isinstance(..., ExteriorStokesMedium)`, matching
        # `ECMBoundaryAnchorFacade`'s Protocol delegate: the requirement is the exterior solve's API,
        # and a nominal check would additionally forbid driving the structural gate with a recording
        # double — which is how every other connector in this tree is gated CUDA-free.
        missing = tuple(
            hook
            for hook in (
                "accumulate", "snapshot_candidate", "rollback",
                "commit_irreversible", "accumulate_ledger",
            )
            if not callable(getattr(self.delegate, hook, None))
        )
        if missing:
            raise TypeError(
                f"membrane medium traction delegate has an incomplete exterior-medium API; "
                f"missing {missing}"
            )

    def accumulate_boundary(self, pos: wp.array, force: wp.array) -> None:
        """Scatter the exterior traction onto the membrane's OWN force array (caller owns zeroing)."""
        self.delegate.accumulate(pos, force)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Mechanics-hook alias for :meth:`accumulate_boundary`."""
        self.delegate.accumulate(pos, force)

    def snapshot_candidate(self) -> None:
        """Forward the transaction snapshot hook."""
        self.delegate.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        """Forward the scheduler-owned predicate without host inspection."""
        self.delegate.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Advance the committed surface reference and traction under the accepted predicate."""
        self.delegate.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward the reservoir reaction and dissipated boundary work."""
        self.delegate.accumulate_ledger(ledger)
