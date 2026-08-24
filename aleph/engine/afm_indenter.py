r"""CUDA runtime for the experiment-only spherical AFM indenter.

The apparatus reuses :func:`aleph.engine.connector_joints._unilateral_contact_force_kernel` **without a
second contact law**.  A sphere is represented exactly for node contact by one indenter-centre point: each
membrane vertex is a degenerate segment ``(i, i)``, the centre is the degenerate segment ``(0, 0)``, and the
connector's zero-force separation is ``sphere_radius + surface_contact_gap``.  The existing gap-gated kernel
therefore pushes an overlapping membrane vertex radially out of the sphere and scatters the equal-and-opposite
force onto the centre.  No triangle broad phase or heuristic pair search is needed: every membrane vertex has
one fixed-capacity pair and separated vertices are structurally force-free.

The centre is prescribed apparatus state, not a massive simulation particle.  A candidate depth is installed
by a CUDA kernel, snapshotted before the physical step and restored on rejection.  ``reaction[0]`` is the
upward force on the centre (the cantilever readout), matching the established ECM-indenter convention;
``reaction[1]`` is the active-contact count.  Reaction is reduced on-device after every force pass.

Sanity Gate:
    * dimensions: centre/radius/gap/depth [um], stiffness [pN/um], force/reaction[0] [pN].
    * sign: for a centre above the cell, an overlapping vertex receives a downward force and the centre an
      equal upward reaction.  At/beyond radius+gap both are exactly zero.
    * conservation: the reused connector kernel scatters ``+f`` and ``-f``; the centre force is the negative
      sum of membrane contact forces to round-off.
    * boundary: radius, gap and stiffness are required positive finite inputs; a non-CUDA device raises;
      rejection restores the centre exactly and accepted state keeps the prescribed depth.
    * runtime: pair/contact/reaction arrays remain device-resident.  Host reads reaction only between accepted
      physical steps; no per-inner-iteration readback exists.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import warp as wp

from aleph.engine.connector_joints import _unilateral_contact_force_kernel
from aleph.engine.load_path import JointState

__all__ = [
    "AFMIndenterForceHook",
    "AFMIndenterRuntime",
    "AFMIndenterRuntimeSpec",
    "sphere_contact_reference",
]


@wp.kernel
def _set_indenter_depth_kernel(
    centre: wp.array(dtype=wp.vec3d),
    initial_centre: wp.array(dtype=wp.vec3d),
    depth_um: wp.float64,
) -> None:
    """Prescribe downward indentation along global ``-z`` from the immutable initial centre."""
    c = initial_centre[0]
    centre[0] = wp.vec3d(c[0], c[1], c[2] - depth_um)


@wp.kernel
def _conditional_restore_centre_kernel(
    centre: wp.array(dtype=wp.vec3d),
    snapshot: wp.array(dtype=wp.vec3d),
    accepted: wp.array(dtype=wp.int32),
) -> None:
    """Restore apparatus candidate state iff the coupled physical step rejected."""
    if accepted[0] == 0:
        centre[0] = snapshot[0]


@wp.kernel
def _indenter_reaction_kernel(
    centre_force: wp.array(dtype=wp.vec3d),
    load: wp.array(dtype=wp.float64),
    reaction: wp.array(dtype=wp.float64),
) -> None:
    """Reduce the predeclared cantilever channel and active-contact count."""
    t = wp.tid()
    if t == 0:
        reaction[0] = centre_force[0][2]
    if load[t] > wp.float64(0.0):
        wp.atomic_add(reaction, 1, wp.float64(1.0))


def _positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive; got {value!r}")
    return value


@dataclass(frozen=True, slots=True)
class AFMIndenterRuntimeSpec:
    """Sourced apparatus/contact values needed to allocate the runtime."""

    radius_um: float
    surface_contact_gap_um: float
    stiffness_pn_per_um: float
    provenance: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "radius_um", _positive(self.radius_um, "indenter radius [um]"))
        object.__setattr__(
            self, "surface_contact_gap_um", _positive(self.surface_contact_gap_um, "surface contact gap [um]")
        )
        object.__setattr__(
            self, "stiffness_pn_per_um", _positive(self.stiffness_pn_per_um, "contact stiffness [pN/um]")
        )
        if not str(self.provenance).strip():
            raise ValueError("indenter/contact provenance is required")

    @property
    def centre_contact_distance_um(self) -> float:
        """The exact ``rest`` value consumed by the reused point-pair kernel."""
        return self.radius_um + self.surface_contact_gap_um


def sphere_contact_reference(
    centre_to_vertex_um: float,
    spec: AFMIndenterRuntimeSpec,
) -> float:
    """Signed membrane-vertex force magnitude for the exact centre-point representation."""
    gap = float(centre_to_vertex_um) - spec.centre_contact_distance_um
    return spec.stiffness_pn_per_um * gap if gap < 0.0 else 0.0


class AFMIndenterRuntime:
    """Fixed-capacity membrane↔sphere connector plus accepted-gated apparatus state."""

    name = "indenter_membrane_contact"
    component_a = "afm_indenter"
    component_b = "membrane"

    def __init__(
        self,
        *,
        membrane_position_d: wp.array,
        membrane_force_d: wp.array,
        initial_centre_um: tuple[float, float, float],
        spec: AFMIndenterRuntimeSpec,
        device: str,
    ) -> None:
        resolved = wp.get_device(device)
        if not resolved.is_cuda:
            raise RuntimeError(f"I0-A: AFM indenter runtime requires CUDA; resolved {resolved}")
        if membrane_position_d.shape != membrane_force_d.shape or len(membrane_position_d.shape) != 1:
            raise ValueError("membrane position/force arrays must be matching one-dimensional vec3 arrays")
        self.device = str(resolved)
        self.membrane_position_d = membrane_position_d
        self.membrane_force_d = membrane_force_d
        self.spec = spec
        self.n_pairs = int(membrane_position_d.shape[0])
        if self.n_pairs <= 0:
            raise ValueError("AFM contact requires at least one membrane vertex")
        centre_np = np.asarray([initial_centre_um], dtype=np.float64)
        if centre_np.shape != (1, 3) or not np.isfinite(centre_np).all():
            raise ValueError("initial indenter centre must be three finite coordinates [um]")
        indices = np.arange(self.n_pairs, dtype=np.int32)
        membrane_segments = np.column_stack([indices, indices]).astype(np.int32, copy=False)
        centre_segments = np.zeros((1, 2), dtype=np.int32)
        engaged = int(JointState.ACTIN_ENGAGED)
        with wp.ScopedDevice(self.device):
            self.centre_d = wp.array(centre_np, dtype=wp.vec3d, device=self.device)
            self.initial_centre_d = wp.array(centre_np, dtype=wp.vec3d, device=self.device)
            self.snapshot_centre_d = wp.empty_like(self.centre_d)
            self.centre_force_d = wp.zeros(1, dtype=wp.vec3d, device=self.device)
            self.membrane_segments_d = wp.array(
                membrane_segments, dtype=wp.int32, ndim=2, device=self.device)
            self.centre_segments_d = wp.array(centre_segments, dtype=wp.int32, ndim=2, device=self.device)
            self.active_d = wp.ones(self.n_pairs, dtype=wp.int32, device=self.device)
            self.state_d = wp.full(self.n_pairs, engaged, dtype=wp.int32, device=self.device)
            self.element_a_d = wp.array(indices, dtype=wp.int32, device=self.device)
            self.element_b_d = wp.zeros(self.n_pairs, dtype=wp.int32, device=self.device)
            self.u_a_d = wp.zeros(self.n_pairs, dtype=wp.float64, device=self.device)
            self.u_b_d = wp.zeros(self.n_pairs, dtype=wp.float64, device=self.device)
            self.stiffness_d = wp.full(
                self.n_pairs, spec.stiffness_pn_per_um, dtype=wp.float64, device=self.device)
            self.rest_d = wp.full(
                self.n_pairs, spec.centre_contact_distance_um, dtype=wp.float64, device=self.device)
            self.load_d = wp.zeros(self.n_pairs, dtype=wp.float64, device=self.device)
            self.energy_d = wp.zeros(self.n_pairs, dtype=wp.float64, device=self.device)
            self.force_on_membrane_d = wp.zeros(self.n_pairs, dtype=wp.vec3d, device=self.device)
            self.reaction_d = wp.zeros(2, dtype=wp.float64, device=self.device)

    def snapshot_candidate(self) -> None:
        """Save the accepted apparatus centre before installing a candidate depth."""
        wp.copy(self.snapshot_centre_d, self.centre_d)

    def set_candidate_depth(self, depth_um: float) -> None:
        """Install a prescribed nonnegative depth without reading simulation state."""
        depth = float(depth_um)
        if not math.isfinite(depth) or depth < 0.0:
            raise ValueError("candidate indentation depth must be finite and nonnegative")
        wp.launch(
            _set_indenter_depth_kernel,
            dim=1,
            inputs=[self.centre_d, self.initial_centre_d, wp.float64(depth)],
            device=self.device,
        )

    def accumulate_on(self, membrane_position_d: wp.array, membrane_force_d: wp.array) -> None:
        """Add contact to the exact candidate membrane arrays supplied by an inner solve."""
        if (
            membrane_position_d.shape != membrane_force_d.shape
            or len(membrane_position_d.shape) != 1
            or int(membrane_position_d.shape[0]) != self.n_pairs
        ):
            raise ValueError(
                "candidate membrane position/force arrays must be matching one-dimensional vec3 arrays "
                f"of the fixed contact capacity {self.n_pairs}"
            )
        self.centre_force_d.zero_()
        wp.launch(
            _unilateral_contact_force_kernel,
            dim=self.n_pairs,
            inputs=[
                membrane_position_d,
                membrane_force_d,
                self.membrane_segments_d,
                self.centre_d,
                self.centre_force_d,
                self.centre_segments_d,
                self.active_d,
                self.state_d,
                self.element_a_d,
                self.element_b_d,
                self.u_a_d,
                self.u_b_d,
                self.stiffness_d,
                self.rest_d,
                self.load_d,
                self.energy_d,
                self.force_on_membrane_d,
            ],
            device=self.device,
        )
        self.reaction_d.zero_()
        wp.launch(
            _indenter_reaction_kernel,
            dim=self.n_pairs,
            inputs=[self.centre_force_d, self.load_d, self.reaction_d],
            device=self.device,
        )

    def accumulate(self) -> None:
        """Add contact to the construction-time membrane views (standalone connector API)."""
        self.accumulate_on(self.membrane_position_d, self.membrane_force_d)

    def rollback(self, accepted_d: wp.array) -> None:
        """Restore rejected prescribed motion under the same device predicate as the cell."""
        wp.launch(
            _conditional_restore_centre_kernel,
            dim=1,
            inputs=[self.centre_d, self.snapshot_centre_d, accepted_d],
            device=self.device,
        )

    def commit_irreversible(self, accepted_d: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Complete the transaction API; prescribed apparatus motion has no irreversible event."""


@dataclass(frozen=True, slots=True)
class AFMIndenterForceHook:
    """Map global inner-solve candidate arrays onto the membrane block without a host readback."""

    runtime: AFMIndenterRuntime
    membrane_node_offset: int
    total_node_count: int

    def __post_init__(self) -> None:
        offset = int(self.membrane_node_offset)
        total = int(self.total_node_count)
        if offset < 0 or total <= 0 or offset + self.runtime.n_pairs > total:
            raise ValueError(
                "AFM membrane block must lie inside the global candidate arrays: "
                f"offset={offset}, n_membrane={self.runtime.n_pairs}, n_total={total}"
            )
        object.__setattr__(self, "membrane_node_offset", offset)
        object.__setattr__(self, "total_node_count", total)

    @property
    def stiffness_bound_pn_per_um(self) -> float:
        """Per-node contact tangent added to the inner solver's explicit stability scale."""
        return self.runtime.spec.stiffness_pn_per_um

    def __call__(self, global_position_d: wp.array, global_force_d: wp.array) -> None:
        """Launch the sixth contact edge on the current trial geometry and force accumulator."""
        if (
            int(global_position_d.shape[0]) != self.total_node_count
            or global_position_d.shape != global_force_d.shape
        ):
            raise ValueError("AFM force hook received arrays outside its declared global node layout")
        start = self.membrane_node_offset
        stop = start + self.runtime.n_pairs
        self.runtime.accumulate_on(
            global_position_d[start:stop],
            global_force_d[start:stop],
        )
