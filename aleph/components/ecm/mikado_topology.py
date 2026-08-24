"""Warp-CUDA construction and material-point contacts for a 3-D Mikado network.

Each collagen fiber is an explicit ordered chain of model points, consecutive axial
segments, and consecutive bending triples.  Random centers/orientations and all topology
arrays are generated directly on CUDA.  A device HashGrid accelerates segment-midpoint
queries, followed by exact segment--segment closest-point geometry.  The result is a
candidate buffer addressed by ``(segment, u, generation)`` on both fibers.

The candidates are not crosslink bonds.  The graph-owned collagen crosslink connector
later selects/associates them under accepted-step kinetics.  This separation prevents a
topology builder from silently creating static chemistry or double-owning bond state.

Sanity Gate:
    * Units: box, contour, segment, capture, material coordinates, rest lengths, and
      separation are all micrometres; ``u`` is dimensionless.
    * Geometry: node order is monotone along one clipped straight rod; segment and bend
      indices never cross fibers.  The midpoint query radius is
      ``max_segment_length + capture`` and therefore cannot prune an admissible pair.
    * Boundary: at least one far-field face is required and only tags topology.  Reactions
      and boundary work remain outside this module.
    * Contact: only distinct-fiber pairs with ``segment_b > segment_a`` are emitted, so
      every geometric segment pair appears once.  Exact closest points, not node proxies,
      determine capture.
    * Residency: generation, query, prefix scan, and candidate fill execute on CUDA.
      One initialization-only scalar read sizes an exact candidate buffer; no method is
      called from the physical-time loop and no geometry is downloaded.
    * Constitutive hold: this module does not import the conflicting collagen modulus
      bands (30--100 Pa versus 5--100 Pa) and cannot close a material-production gate.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import IntFlag

import warp as wp

from aleph.components.ecm.device_schema import ECMTopologyState, require_cuda_device

__all__ = ["MikadoInitConfig", "MikadoTopologyBuilder", "SegmentContactCandidates"]


class _Face(IntFlag):
    X_LO = 1
    X_HI = 2
    Y_LO = 4
    Y_HI = 8
    Z_LO = 16
    Z_HI = 32


_FACE_BY_NAME = {
    "x_lo": _Face.X_LO,
    "x_hi": _Face.X_HI,
    "y_lo": _Face.Y_LO,
    "y_hi": _Face.Y_HI,
    "z_lo": _Face.Z_LO,
    "z_hi": _Face.Z_HI,
}


@dataclass(frozen=True, slots=True)
class MikadoInitConfig:
    """Physically specified initialization inputs for one 3-D collagen network.

    ``n_fibers`` must come from the concentration/geometry population calculation;
    this class never lowers it to fit memory.  ``pin_margin_um`` and ``pin_faces`` are
    required physiological boundary inputs rather than a relaxed/free default.
    """

    box_lo_um: tuple[float, float, float]
    box_hi_um: tuple[float, float, float]
    n_fibers: int
    fiber_length_um: float
    target_segment_um: float
    crosslink_capture_um: float
    pin_faces: tuple[str, ...]
    pin_margin_um: float
    rng_seed: int
    max_refinement_level: int
    persistent_id_base: int = 0

    def __post_init__(self) -> None:
        if len(self.box_lo_um) != 3 or len(self.box_hi_um) != 3:
            raise ValueError("Mikado topology requires exactly three box coordinates")
        values = (*self.box_lo_um, *self.box_hi_um)
        if any(isinstance(value, bool) or not math.isfinite(value) for value in values):
            raise ValueError("box coordinates must be finite numbers")
        if any(hi <= lo for lo, hi in zip(self.box_lo_um, self.box_hi_um, strict=True)):
            raise ValueError("every box_hi coordinate must exceed box_lo")
        if (
            isinstance(self.n_fibers, bool)
            or not isinstance(self.n_fibers, int)
            or self.n_fibers <= 0
        ):
            raise ValueError("n_fibers must be a positive integer from physiological density")
        for label, value in (
            ("fiber_length_um", self.fiber_length_um),
            ("target_segment_um", self.target_segment_um),
            ("crosslink_capture_um", self.crosslink_capture_um),
            ("pin_margin_um", self.pin_margin_um),
        ):
            if isinstance(value, bool) or not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{label} must be finite and positive")
        unknown = set(self.pin_faces) - set(_FACE_BY_NAME)
        if unknown:
            raise ValueError(f"unknown far-field faces: {sorted(unknown)}")
        if not self.pin_faces:
            raise ValueError("production ECM cannot be initialized as a free network")
        if len(set(self.pin_faces)) != len(self.pin_faces):
            raise ValueError("pin_faces must not contain duplicates")
        extents = tuple(hi - lo for lo, hi in zip(self.box_lo_um, self.box_hi_um, strict=True))
        if self.pin_margin_um > min(extents):
            raise ValueError("pin_margin_um cannot exceed the smallest box extent")
        if (
            isinstance(self.rng_seed, bool)
            or not isinstance(self.rng_seed, int)
            or not -(2**31) <= self.rng_seed < 2**31
        ):
            raise ValueError("rng_seed must fit the Warp int32 RNG contract")
        if (
            isinstance(self.max_refinement_level, bool)
            or not isinstance(self.max_refinement_level, int)
            or self.max_refinement_level < 0
        ):
            raise ValueError("max_refinement_level must be a nonnegative integer")
        if isinstance(self.persistent_id_base, bool) or not isinstance(
            self.persistent_id_base, int
        ):
            raise ValueError("persistent_id_base must be an integer")
        last_id = self.persistent_id_base + self.initial_identity_count - 1
        if not -(2**63) <= self.persistent_id_base <= last_id < 2**63:
            raise ValueError("initial persistent topology IDs must fit int64 without wrapping")
        capacities = (
            self.n_fibers * self.node_slots_per_fiber,
            self.n_fibers * self.segment_slots_per_fiber,
            self.n_fibers * self.bend_slots_per_fiber,
        )
        if max(capacities) >= 2**31:
            raise ValueError("refinement topology exceeds int32 device indexing")

    @property
    def nodes_per_fiber(self) -> int:
        """Smallest chain whose uncut segment length does not exceed the target."""
        return max(2, math.ceil(self.fiber_length_um / self.target_segment_um) + 1)

    @property
    def max_segment_length_um(self) -> float:
        """Conservative segment bound before line--box clipping."""
        return self.fiber_length_um / (self.nodes_per_fiber - 1)

    @property
    def segment_slots_per_fiber(self) -> int:
        """Worst-case segment slots under binary refinement to the declared level."""
        return (self.nodes_per_fiber - 1) * (1 << self.max_refinement_level)

    @property
    def node_slots_per_fiber(self) -> int:
        """Worst-case ordered node slots under binary refinement."""
        return self.segment_slots_per_fiber + 1

    @property
    def bend_slots_per_fiber(self) -> int:
        """Worst-case consecutive bending-triple slots under binary refinement."""
        return max(0, self.node_slots_per_fiber - 2)

    @property
    def initial_identity_count(self) -> int:
        """Persistent IDs consumed by initial fibers, nodes, and segments."""
        return self.n_fibers * (2 * self.nodes_per_fiber)

    @property
    def pin_face_mask(self) -> int:
        """Bit mask consumed by the CUDA boundary-tag kernel."""
        mask = _Face(0)
        for face in self.pin_faces:
            mask |= _FACE_BY_NAME[face]
        return int(mask)

    @property
    def grid_dimensions(self) -> tuple[int, int, int]:
        """HashGrid dimensions derived from the conservative contact query radius."""
        radius = self.max_segment_length_um + self.crosslink_capture_um
        return tuple(
            max(1, math.ceil((hi - lo) / radius))
            for lo, hi in zip(self.box_lo_um, self.box_hi_um, strict=True)
        )


@wp.func
def _clipped_rod_interval(
    center: wp.vec3d,
    direction: wp.vec3d,
    box_lo: wp.vec3d,
    box_hi: wp.vec3d,
    half_length: wp.float64,
):
    """Intersect ``center + t direction`` with the box and desired contour interval."""
    t_lo = -half_length
    t_hi = half_length
    tiny = wp.float64(1.0e-15)
    for axis in range(3):
        velocity = direction[axis]
        if wp.abs(velocity) > tiny:
            ta = (box_lo[axis] - center[axis]) / velocity
            tb = (box_hi[axis] - center[axis]) / velocity
            axis_lo = wp.min(ta, tb)
            axis_hi = wp.max(ta, tb)
            t_lo = wp.max(t_lo, axis_lo)
            t_hi = wp.min(t_hi, axis_hi)
    return t_lo, t_hi


@wp.func
def _node_boundary_mask(
    point: wp.vec3d,
    box_lo: wp.vec3d,
    box_hi: wp.vec3d,
    selected_faces: wp.int32,
    margin: wp.float64,
) -> wp.int32:
    """Return selected box faces whose boundary band contains ``point``."""
    mask = wp.int32(0)
    if selected_faces & wp.int32(1) and point[0] - box_lo[0] <= margin:
        mask = mask | wp.int32(1)
    if selected_faces & wp.int32(2) and box_hi[0] - point[0] <= margin:
        mask = mask | wp.int32(2)
    if selected_faces & wp.int32(4) and point[1] - box_lo[1] <= margin:
        mask = mask | wp.int32(4)
    if selected_faces & wp.int32(8) and box_hi[1] - point[1] <= margin:
        mask = mask | wp.int32(8)
    if selected_faces & wp.int32(16) and point[2] - box_lo[2] <= margin:
        mask = mask | wp.int32(16)
    if selected_faces & wp.int32(32) and box_hi[2] - point[2] <= margin:
        mask = mask | wp.int32(32)
    return mask


@wp.kernel
def _initialize_mikado_kernel(
    box_lo: wp.vec3d,
    box_hi: wp.vec3d,
    fiber_length: wp.float64,
    n_fibers: wp.int32,
    nodes_per_fiber: wp.int32,
    node_slots_per_fiber: wp.int32,
    segment_slots_per_fiber: wp.int32,
    bend_slots_per_fiber: wp.int32,
    seed: wp.int32,
    persistent_id_base: wp.int64,
    selected_faces: wp.int32,
    pin_margin: wp.float64,
    position: wp.array(dtype=wp.vec3d),
    reference_position: wp.array(dtype=wp.vec3d),
    node_active: wp.array(dtype=wp.int32),
    node_fiber: wp.array(dtype=wp.int32),
    node_local_index: wp.array(dtype=wp.int32),
    node_generation: wp.array(dtype=wp.int32),
    persistent_node_id: wp.array(dtype=wp.int64),
    boundary_face_mask: wp.array(dtype=wp.int32),
    fiber_offsets: wp.array(dtype=wp.int32),
    fiber_active: wp.array(dtype=wp.int32),
    persistent_fiber_id: wp.array(dtype=wp.int64),
    fiber_generation: wp.array(dtype=wp.int32),
    segments: wp.array2d(dtype=wp.int32),
    segment_fiber: wp.array(dtype=wp.int32),
    segment_active: wp.array(dtype=wp.int32),
    segment_generation: wp.array(dtype=wp.int32),
    persistent_segment_id: wp.array(dtype=wp.int64),
    refinement_root_segment_id: wp.array(dtype=wp.int64),
    refinement_path: wp.array(dtype=wp.int64),
    segment_prev: wp.array(dtype=wp.int32),
    segment_next: wp.array(dtype=wp.int32),
    segment_bend: wp.array(dtype=wp.int32),
    segment_material_s0: wp.array(dtype=wp.float64),
    segment_material_s1: wp.array(dtype=wp.float64),
    segment_rest_length: wp.array(dtype=wp.float64),
    bend_triples: wp.array2d(dtype=wp.int32),
    bend_fiber: wp.array(dtype=wp.int32),
    bend_active: wp.array(dtype=wp.int32),
    bend_left_segment: wp.array(dtype=wp.int32),
    bend_right_segment: wp.array(dtype=wp.int32),
) -> None:
    """Generate one disjoint straight-chain fiber per CUDA thread."""
    fiber = wp.tid()
    state = wp.rand_init(seed, fiber)
    center = wp.vec3d(
        box_lo[0] + wp.float64(wp.randf(state)) * (box_hi[0] - box_lo[0]),
        box_lo[1] + wp.float64(wp.randf(state)) * (box_hi[1] - box_lo[1]),
        box_lo[2] + wp.float64(wp.randf(state)) * (box_hi[2] - box_lo[2]),
    )
    sampled = wp.sample_unit_sphere(state)
    direction = wp.vec3d(wp.float64(sampled[0]), wp.float64(sampled[1]), wp.float64(sampled[2]))
    direction = wp.normalize(direction)
    t_lo, t_hi = _clipped_rod_interval(
        center, direction, box_lo, box_hi, wp.float64(0.5) * fiber_length
    )
    first = center + t_lo * direction
    last = center + t_hi * direction
    clipped_length = wp.length(last - first)
    denominator = wp.float64(nodes_per_fiber - wp.int32(1))
    node_base = fiber * node_slots_per_fiber
    segment_base = fiber * segment_slots_per_fiber
    bend_base = fiber * bend_slots_per_fiber

    fiber_offsets[fiber] = node_base
    fiber_active[fiber] = wp.int32(1)
    persistent_fiber_id[fiber] = persistent_id_base + wp.int64(fiber)
    fiber_generation[fiber] = wp.int32(0)

    for local in range(nodes_per_fiber):
        node = node_base + local
        fraction = wp.float64(local) / denominator
        point = first + fraction * (last - first)
        position[node] = point
        reference_position[node] = point
        node_active[node] = wp.int32(1)
        node_fiber[node] = fiber
        node_local_index[node] = local
        node_generation[node] = wp.int32(0)
        packed_node = fiber * nodes_per_fiber + local
        persistent_node_id[node] = persistent_id_base + wp.int64(n_fibers + packed_node)
        boundary_face_mask[node] = _node_boundary_mask(
            point, box_lo, box_hi, selected_faces, pin_margin
        )

    if fiber == n_fibers - wp.int32(1):
        fiber_offsets[n_fibers] = node_base + node_slots_per_fiber
    for local in range(nodes_per_fiber - wp.int32(1)):
        segment = segment_base + local
        node_a = node_base + local
        node_b = node_a + wp.int32(1)
        s0 = clipped_length * wp.float64(local) / denominator
        s1 = clipped_length * wp.float64(local + wp.int32(1)) / denominator
        segments[segment, 0] = node_a
        segments[segment, 1] = node_b
        segment_fiber[segment] = fiber
        segment_active[segment] = wp.int32(1)
        segment_generation[segment] = wp.int32(0)
        packed_segment = fiber * (nodes_per_fiber - wp.int32(1)) + local
        segment_id = persistent_id_base + wp.int64(
            n_fibers + n_fibers * nodes_per_fiber + packed_segment
        )
        persistent_segment_id[segment] = segment_id
        refinement_root_segment_id[segment] = segment_id
        refinement_path[segment] = wp.int64(0)
        if local > wp.int32(0):
            segment_prev[segment] = segment - wp.int32(1)
        if local + wp.int32(1) < nodes_per_fiber - wp.int32(1):
            segment_next[segment] = segment + wp.int32(1)
            segment_bend[segment] = bend_base + local
        segment_material_s0[segment] = s0
        segment_material_s1[segment] = s1
        segment_rest_length[segment] = s1 - s0

    for local in range(wp.max(nodes_per_fiber - wp.int32(2), wp.int32(0))):
        bend = bend_base + local
        node = node_base + local
        bend_triples[bend, 0] = node
        bend_triples[bend, 1] = node + wp.int32(1)
        bend_triples[bend, 2] = node + wp.int32(2)
        bend_fiber[bend] = fiber
        bend_active[bend] = wp.int32(1)
        bend_left_segment[bend] = segment_base + local
        bend_right_segment[bend] = segment_base + local + wp.int32(1)


@wp.kernel
def _initialize_population_kernel(
    n_fibers: wp.int32,
    n_nodes: wp.int32,
    n_segments: wp.int32,
    n_bends: wp.int32,
    population: wp.array(dtype=wp.int32),
) -> None:
    """Initialize device-resident active counts in one CUDA thread."""
    if wp.tid() == 0:
        population[wp.int32(0)] = n_fibers
        population[wp.int32(1)] = n_nodes
        population[wp.int32(2)] = n_segments
        population[wp.int32(3)] = n_bends


@wp.kernel
def _initialize_node_free_list_kernel(
    active_slots_per_fiber: wp.int32,
    capacity_slots_per_fiber: wp.int32,
    free_list: wp.array(dtype=wp.int32),
) -> None:
    """Pack every dormant node slot into the device free stack."""
    slot = wp.tid()
    local = slot % capacity_slots_per_fiber
    if local >= active_slots_per_fiber:
        fiber = slot // capacity_slots_per_fiber
        dormant_per_fiber = capacity_slots_per_fiber - active_slots_per_fiber
        free_index = fiber * dormant_per_fiber + local - active_slots_per_fiber
        free_list[free_index] = slot


@wp.kernel
def _initialize_element_free_list_kernel(
    active_slots_per_fiber: wp.int32,
    capacity_slots_per_fiber: wp.int32,
    free_list: wp.array(dtype=wp.int32),
) -> None:
    """Pack dormant segment or bend slots into its device free stack."""
    slot = wp.tid()
    local = slot % capacity_slots_per_fiber
    if local >= active_slots_per_fiber:
        fiber = slot // capacity_slots_per_fiber
        dormant_per_fiber = capacity_slots_per_fiber - active_slots_per_fiber
        free_index = fiber * dormant_per_fiber + local - active_slots_per_fiber
        free_list[free_index] = slot


@wp.kernel
def _initialize_free_population_kernel(
    free_nodes: wp.int32,
    free_segments: wp.int32,
    free_bends: wp.int32,
    free_population: wp.array(dtype=wp.int32),
) -> None:
    """Publish valid free-stack extents on the device."""
    if wp.tid() == 0:
        free_population[0] = free_nodes
        free_population[1] = free_segments
        free_population[2] = free_bends


@wp.kernel
def _initialize_identity_cursor_kernel(
    next_persistent_id: wp.int64,
    cursor: wp.array(dtype=wp.int64),
) -> None:
    """Set the device-only identity allocator after all initial identities."""
    if wp.tid() == 0:
        cursor[0] = next_persistent_id


@wp.func
def _closest_segment_points(
    p0: wp.vec3d,
    p1: wp.vec3d,
    q0: wp.vec3d,
    q1: wp.vec3d,
):
    """Return clamped material coordinates and closest points on two 3-D segments."""
    d1 = p1 - p0
    d2 = q1 - q0
    r = p0 - q0
    a = wp.dot(d1, d1)
    e = wp.dot(d2, d2)
    f = wp.dot(d2, r)
    eps = wp.float64(1.0e-24)
    s = wp.float64(0.0)
    t = wp.float64(0.0)
    if a <= eps and e <= eps:
        return s, t, p0, q0, wp.length(p0 - q0)
    if a <= eps:
        t = wp.clamp(f / e, wp.float64(0.0), wp.float64(1.0))
    else:
        c = wp.dot(d1, r)
        if e <= eps:
            s = wp.clamp(-c / a, wp.float64(0.0), wp.float64(1.0))
        else:
            b = wp.dot(d1, d2)
            denominator = a * e - b * b
            if denominator > eps:
                s = wp.clamp(
                    (b * f - c * e) / denominator,
                    wp.float64(0.0),
                    wp.float64(1.0),
                )
            t = (b * s + f) / e
            if t < wp.float64(0.0):
                t = wp.float64(0.0)
                s = wp.clamp(-c / a, wp.float64(0.0), wp.float64(1.0))
            elif t > wp.float64(1.0):
                t = wp.float64(1.0)
                s = wp.clamp((b - c) / a, wp.float64(0.0), wp.float64(1.0))
    point_p = p0 + s * d1
    point_q = q0 + t * d2
    return s, t, point_p, point_q, wp.length(point_p - point_q)


@wp.kernel
def _segment_midpoints_kernel(
    position: wp.array(dtype=wp.vec3d),
    segments: wp.array2d(dtype=wp.int32),
    active_segments_per_fiber: wp.int32,
    segment_slots_per_fiber: wp.int32,
    midpoint: wp.array(dtype=wp.vec3d),
) -> None:
    """Build exact-f64 midpoint query points for the device HashGrid."""
    packed_segment = wp.tid()
    fiber = packed_segment // active_segments_per_fiber
    local = packed_segment % active_segments_per_fiber
    segment = fiber * segment_slots_per_fiber + local
    midpoint[packed_segment] = wp.float64(0.5) * (
        position[segments[segment, 0]] + position[segments[segment, 1]]
    )


@wp.kernel
def _count_segment_contacts_kernel(
    grid: wp.uint64,
    midpoint: wp.array(dtype=wp.vec3d),
    position: wp.array(dtype=wp.vec3d),
    segments: wp.array2d(dtype=wp.int32),
    segment_fiber: wp.array(dtype=wp.int32),
    segment_active: wp.array(dtype=wp.int32),
    active_segments_per_fiber: wp.int32,
    segment_slots_per_fiber: wp.int32,
    query_radius: wp.float64,
    capture_radius: wp.float64,
    counts: wp.array(dtype=wp.int32),
) -> None:
    """Count each canonical, exact segment-contact pair on CUDA."""
    packed_a = wp.tid()
    fiber_a = packed_a // active_segments_per_fiber
    local_a = packed_a % active_segments_per_fiber
    segment_a = fiber_a * segment_slots_per_fiber + local_a
    count = wp.int32(0)
    if segment_active[segment_a] > wp.int32(0):
        query = wp.hash_grid_query(grid, midpoint[packed_a], query_radius)
        packed_b = wp.int32(0)
        while wp.hash_grid_query_next(query, packed_b):
            fiber_b = packed_b // active_segments_per_fiber
            local_b = packed_b % active_segments_per_fiber
            segment_b = fiber_b * segment_slots_per_fiber + local_b
            if (
                packed_b > packed_a
                and segment_active[segment_b] > wp.int32(0)
                and segment_fiber[segment_b] != segment_fiber[segment_a]
            ):
                _u_a, _u_b, _point_a, _point_b, distance = _closest_segment_points(
                    position[segments[segment_a, 0]],
                    position[segments[segment_a, 1]],
                    position[segments[segment_b, 0]],
                    position[segments[segment_b, 1]],
                )
                if distance <= capture_radius:
                    count = count + wp.int32(1)
    counts[packed_a] = count


@wp.kernel
def _contact_total_kernel(
    counts: wp.array(dtype=wp.int32),
    offsets: wp.array(dtype=wp.int32),
    total: wp.array(dtype=wp.int32),
) -> None:
    """Extract the final exclusive-scan extent into one device scalar."""
    if wp.tid() == 0:
        last = wp.int32(counts.shape[0] - 1)
        total[0] = offsets[last] + counts[last]


@wp.kernel
def _fill_segment_contacts_kernel(
    grid: wp.uint64,
    midpoint: wp.array(dtype=wp.vec3d),
    position: wp.array(dtype=wp.vec3d),
    segments: wp.array2d(dtype=wp.int32),
    segment_fiber: wp.array(dtype=wp.int32),
    segment_active: wp.array(dtype=wp.int32),
    segment_generation: wp.array(dtype=wp.int32),
    persistent_segment_id: wp.array(dtype=wp.int64),
    offsets: wp.array(dtype=wp.int32),
    active_segments_per_fiber: wp.int32,
    segment_slots_per_fiber: wp.int32,
    query_radius: wp.float64,
    capture_radius: wp.float64,
    out_segment_a: wp.array(dtype=wp.int32),
    out_segment_b: wp.array(dtype=wp.int32),
    out_fiber_a: wp.array(dtype=wp.int32),
    out_fiber_b: wp.array(dtype=wp.int32),
    out_u_a: wp.array(dtype=wp.float64),
    out_u_b: wp.array(dtype=wp.float64),
    out_generation_a: wp.array(dtype=wp.int32),
    out_generation_b: wp.array(dtype=wp.int32),
    out_persistent_segment_a: wp.array(dtype=wp.int64),
    out_persistent_segment_b: wp.array(dtype=wp.int64),
    out_separation: wp.array(dtype=wp.float64),
) -> None:
    """Fill the exact candidate buffer using the previously scanned row offsets."""
    packed_a = wp.tid()
    fiber_a = packed_a // active_segments_per_fiber
    local_a = packed_a % active_segments_per_fiber
    segment_a = fiber_a * segment_slots_per_fiber + local_a
    if segment_active[segment_a] <= wp.int32(0):
        return
    cursor = offsets[packed_a]
    query = wp.hash_grid_query(grid, midpoint[packed_a], query_radius)
    packed_b = wp.int32(0)
    while wp.hash_grid_query_next(query, packed_b):
        fiber_b = packed_b // active_segments_per_fiber
        local_b = packed_b % active_segments_per_fiber
        segment_b = fiber_b * segment_slots_per_fiber + local_b
        if (
            packed_b > packed_a
            and segment_active[segment_b] > wp.int32(0)
            and segment_fiber[segment_b] != segment_fiber[segment_a]
        ):
            u_a, u_b, _point_a, _point_b, distance = _closest_segment_points(
                position[segments[segment_a, 0]],
                position[segments[segment_a, 1]],
                position[segments[segment_b, 0]],
                position[segments[segment_b, 1]],
            )
            if distance <= capture_radius:
                out_segment_a[cursor] = segment_a
                out_segment_b[cursor] = segment_b
                out_fiber_a[cursor] = segment_fiber[segment_a]
                out_fiber_b[cursor] = segment_fiber[segment_b]
                out_u_a[cursor] = u_a
                out_u_b[cursor] = u_b
                out_generation_a[cursor] = segment_generation[segment_a]
                out_generation_b[cursor] = segment_generation[segment_b]
                out_persistent_segment_a[cursor] = persistent_segment_id[segment_a]
                out_persistent_segment_b[cursor] = persistent_segment_id[segment_b]
                out_separation[cursor] = distance
                cursor = cursor + wp.int32(1)


@dataclass(frozen=True, slots=True)
class SegmentContactCandidates:
    """Graph-input candidates; no association state or crosslink chemistry is stored."""

    device: str
    count_d: wp.array
    segment_a_d: wp.array
    segment_b_d: wp.array
    fiber_a_d: wp.array
    fiber_b_d: wp.array
    u_a_d: wp.array
    u_b_d: wp.array
    generation_a_d: wp.array
    generation_b_d: wp.array
    persistent_segment_a_d: wp.array
    persistent_segment_b_d: wp.array
    separation_d: wp.array
    initialization_count: int

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Return candidate/scratch ownership for device-byte accounting."""
        return (
            self.count_d,
            self.segment_a_d,
            self.segment_b_d,
            self.fiber_a_d,
            self.fiber_b_d,
            self.u_a_d,
            self.u_b_d,
            self.generation_a_d,
            self.generation_b_d,
            self.persistent_segment_a_d,
            self.persistent_segment_b_d,
            self.separation_d,
        )

    def allocated_array_bytes(self) -> int:
        """Return exact Warp-array capacity bytes (HashGrid storage is separate)."""
        return sum(int(array.capacity) for array in self.owned_arrays())


class MikadoTopologyBuilder:
    """One-shot CUDA builder for initial collagen geometry and capture candidates."""

    def __init__(self, config: MikadoInitConfig, *, device: str | None = None) -> None:
        self.config = config
        self.device = require_cuda_device(device)

    def initialize(self) -> ECMTopologyState:
        """Allocate and populate the explicit fiber/segment/bend topology on CUDA."""
        cfg = self.config
        topology = ECMTopologyState.allocate(
            n_fibers=cfg.n_fibers,
            initial_nodes_per_fiber=cfg.nodes_per_fiber,
            node_slots_per_fiber=cfg.node_slots_per_fiber,
            device=self.device,
        )
        with wp.ScopedDevice(self.device):
            wp.launch(
                _initialize_mikado_kernel,
                dim=cfg.n_fibers,
                inputs=[
                    wp.vec3d(*cfg.box_lo_um),
                    wp.vec3d(*cfg.box_hi_um),
                    wp.float64(cfg.fiber_length_um),
                    wp.int32(cfg.n_fibers),
                    wp.int32(cfg.nodes_per_fiber),
                    wp.int32(cfg.node_slots_per_fiber),
                    wp.int32(cfg.segment_slots_per_fiber),
                    wp.int32(cfg.bend_slots_per_fiber),
                    wp.int32(cfg.rng_seed),
                    wp.int64(cfg.persistent_id_base),
                    wp.int32(cfg.pin_face_mask),
                    wp.float64(cfg.pin_margin_um),
                ],
                outputs=[
                    topology.position_d,
                    topology.reference_position_d,
                    topology.node_active_d,
                    topology.node_fiber_d,
                    topology.node_local_index_d,
                    topology.node_generation_d,
                    topology.persistent_node_id_d,
                    topology.boundary_face_mask_d,
                    topology.fiber_offsets_d,
                    topology.fiber_active_d,
                    topology.persistent_fiber_id_d,
                    topology.fiber_generation_d,
                    topology.segments_d,
                    topology.segment_fiber_d,
                    topology.segment_active_d,
                    topology.segment_generation_d,
                    topology.persistent_segment_id_d,
                    topology.refinement_root_segment_id_d,
                    topology.refinement_path_d,
                    topology.segment_prev_d,
                    topology.segment_next_d,
                    topology.segment_bend_d,
                    topology.segment_material_s0_d,
                    topology.segment_material_s1_d,
                    topology.segment_rest_length_d,
                    topology.bend_triples_d,
                    topology.bend_fiber_d,
                    topology.bend_active_d,
                    topology.bend_left_segment_d,
                    topology.bend_right_segment_d,
                ],
                device=self.device,
            )
            wp.launch(
                _initialize_population_kernel,
                dim=1,
                inputs=[
                    wp.int32(topology.n_fiber_capacity),
                    wp.int32(topology.n_initial_nodes),
                    wp.int32(topology.n_initial_segments),
                    wp.int32(topology.n_initial_bends),
                ],
                outputs=[topology.active_population_d],
                device=self.device,
            )
            wp.launch(
                _initialize_node_free_list_kernel,
                dim=topology.n_node_capacity,
                inputs=[
                    wp.int32(cfg.nodes_per_fiber),
                    wp.int32(cfg.node_slots_per_fiber),
                ],
                outputs=[topology.node_free_list_d],
                device=self.device,
            )
            wp.launch(
                _initialize_element_free_list_kernel,
                dim=topology.n_segment_capacity,
                inputs=[
                    wp.int32(cfg.nodes_per_fiber - 1),
                    wp.int32(cfg.segment_slots_per_fiber),
                ],
                outputs=[topology.segment_free_list_d],
                device=self.device,
            )
            if topology.n_bend_capacity:
                wp.launch(
                    _initialize_element_free_list_kernel,
                    dim=topology.n_bend_capacity,
                    inputs=[
                        wp.int32(max(0, cfg.nodes_per_fiber - 2)),
                        wp.int32(cfg.bend_slots_per_fiber),
                    ],
                    outputs=[topology.bend_free_list_d],
                    device=self.device,
                )
            wp.launch(
                _initialize_free_population_kernel,
                dim=1,
                inputs=[
                    wp.int32(topology.n_node_capacity - topology.n_initial_nodes),
                    wp.int32(topology.n_segment_capacity - topology.n_initial_segments),
                    wp.int32(topology.n_bend_capacity - topology.n_initial_bends),
                ],
                outputs=[topology.free_population_d],
                device=self.device,
            )
            wp.launch(
                _initialize_identity_cursor_kernel,
                dim=1,
                inputs=[wp.int64(cfg.persistent_id_base + cfg.initial_identity_count)],
                outputs=[topology.next_persistent_id_d],
                device=self.device,
            )
        return topology

    def initial_contact_candidates(
        self,
        topology: ECMTopologyState,
    ) -> SegmentContactCandidates:
        """Build exact initial material-point capture candidates on CUDA.

        The only device-to-host transfer is the final integer extent of the initialization
        prefix scan.  It sizes the immutable candidate allocation before physical time
        starts; neither topology geometry nor a per-step decision is downloaded.
        """
        if topology.device != self.device:
            raise ValueError("topology and Mikado builder must use the same CUDA device")
        if topology.n_fiber_capacity != self.config.n_fibers:
            raise ValueError("topology fiber capacity does not match the initialization config")
        cfg = self.config
        n_segments = topology.n_initial_segments
        active_segments_per_fiber = cfg.nodes_per_fiber - 1
        dims = cfg.grid_dimensions
        query_radius = cfg.max_segment_length_um + cfg.crosslink_capture_um
        with wp.ScopedDevice(self.device):
            midpoint_d = wp.zeros(n_segments, dtype=wp.vec3d)
            counts_d = wp.zeros(n_segments, dtype=wp.int32)
            offsets_d = wp.zeros(n_segments, dtype=wp.int32)
            total_d = wp.zeros(1, dtype=wp.int32)
            grid = wp.HashGrid(*dims, device=self.device, dtype=wp.float64)
            wp.launch(
                _segment_midpoints_kernel,
                dim=n_segments,
                inputs=[
                    topology.position_d,
                    topology.segments_d,
                    wp.int32(active_segments_per_fiber),
                    wp.int32(cfg.segment_slots_per_fiber),
                ],
                outputs=[midpoint_d],
                device=self.device,
            )
            grid.build(midpoint_d, query_radius)
            wp.launch(
                _count_segment_contacts_kernel,
                dim=n_segments,
                inputs=[
                    grid.id,
                    midpoint_d,
                    topology.position_d,
                    topology.segments_d,
                    topology.segment_fiber_d,
                    topology.segment_active_d,
                    wp.int32(active_segments_per_fiber),
                    wp.int32(cfg.segment_slots_per_fiber),
                    wp.float64(query_radius),
                    wp.float64(cfg.crosslink_capture_um),
                ],
                outputs=[counts_d],
                device=self.device,
            )
            wp.utils.array_scan(counts_d, offsets_d, inclusive=False)
            wp.launch(
                _contact_total_kernel,
                dim=1,
                inputs=[counts_d, offsets_d],
                outputs=[total_d],
                device=self.device,
            )
            candidate_count = int(total_d.numpy()[0])
            if candidate_count < 0:
                raise RuntimeError("segment contact count overflowed int32")
            count_d = total_d
            segment_a_d = wp.zeros(candidate_count, dtype=wp.int32)
            segment_b_d = wp.zeros(candidate_count, dtype=wp.int32)
            fiber_a_d = wp.zeros(candidate_count, dtype=wp.int32)
            fiber_b_d = wp.zeros(candidate_count, dtype=wp.int32)
            u_a_d = wp.zeros(candidate_count, dtype=wp.float64)
            u_b_d = wp.zeros(candidate_count, dtype=wp.float64)
            generation_a_d = wp.zeros(candidate_count, dtype=wp.int32)
            generation_b_d = wp.zeros(candidate_count, dtype=wp.int32)
            persistent_segment_a_d = wp.zeros(candidate_count, dtype=wp.int64)
            persistent_segment_b_d = wp.zeros(candidate_count, dtype=wp.int64)
            separation_d = wp.zeros(candidate_count, dtype=wp.float64)
            if candidate_count:
                wp.launch(
                    _fill_segment_contacts_kernel,
                    dim=n_segments,
                    inputs=[
                        grid.id,
                        midpoint_d,
                        topology.position_d,
                        topology.segments_d,
                        topology.segment_fiber_d,
                        topology.segment_active_d,
                        topology.segment_generation_d,
                        topology.persistent_segment_id_d,
                        offsets_d,
                        wp.int32(active_segments_per_fiber),
                        wp.int32(cfg.segment_slots_per_fiber),
                        wp.float64(query_radius),
                        wp.float64(cfg.crosslink_capture_um),
                    ],
                    outputs=[
                        segment_a_d,
                        segment_b_d,
                        fiber_a_d,
                        fiber_b_d,
                        u_a_d,
                        u_b_d,
                        generation_a_d,
                        generation_b_d,
                        persistent_segment_a_d,
                        persistent_segment_b_d,
                        separation_d,
                    ],
                    device=self.device,
                )
        return SegmentContactCandidates(
            device=self.device,
            count_d=count_d,
            segment_a_d=segment_a_d,
            segment_b_d=segment_b_d,
            fiber_a_d=fiber_a_d,
            fiber_b_d=fiber_b_d,
            u_a_d=u_a_d,
            u_b_d=u_b_d,
            generation_a_d=generation_a_d,
            generation_b_d=generation_b_d,
            persistent_segment_a_d=persistent_segment_a_d,
            persistent_segment_b_d=persistent_segment_b_d,
            separation_d=separation_d,
            initialization_count=candidate_count,
        )
