"""The actin cortex: one strand POPULATION, claimed once and built by kernels on the device.

THE CLAIM MODEL IS THE DESIGN DECISION THIS FILE CARRIES.  ``strand.build_strand`` claims three ranges
per filament — nodes, segments, angles — which is correct for one filament and wrong for seventy
thousand: the native cortex would arrive as **212,058 claims**, and every arena operation that walks
the claim list (``population_of``, ``assert_partitioned``, ``census``, ``assert_disjoint_spans``) is
linear or quadratic in that number.  ``population_of`` alone becomes a 212,058-entry scan per query.

The arena already says what the unit is: *"a population is a NAMED, CONTIGUOUS, HALF-OPEN ID RANGE"*.
So this module claims **three ranges for the whole cortex**, and the per-filament addressing that
``build_strand`` got from separate claims comes from arithmetic instead::

    filament f, node k     ->  nodes.lo + f * nodes_per_strand + k
    filament f, segment k  ->  segments.lo + f * (nodes_per_strand - 1) + k

which is exact because every filament in a population has the SAME node count.  That is not a
restriction smuggled in for convenience — it is what a population IS here.  The cortex is not one
label-blind network: formin-nucleated long filaments and Arp2/3-nucleated short ones are two different
lengths and therefore two different populations, each uniform, coupled by connectors and never by a
shared range.  ``CLAUDE.md``'s "one filament, one component" already required exactly that split.

CONSTRUCTION IS A KERNEL, AND THE FILAMENT SEEDS NEVER EXIST ON THE HOST.  PI 2026-08-20 settled that
the arena is built by Warp CUDA kernels rather than assembled in NumPy and uploaded.  Taken literally
that would still allow a host array of 70,686 filament orientations to be uploaded; this module does
not do that either.  Each filament's centre direction, tangent and shell radius are drawn INSIDE the
kernel from ``wp.rand_init(seed, filament)``, so every node of a filament re-derives the same frame and
no per-filament array is ever allocated on either side.  The redundancy is ~61 extra draws per filament
in an arithmetic-bound kernel, against one fewer host allocation and one fewer transfer.  The build is
reproducible from ``(seed, n_filaments)`` alone.

THE REST LENGTH IS MEASURED FROM WHAT WAS BUILT, NOT FROM WHAT WAS ASKED FOR.  A filament is laid as a
great-circle arc on its own shell radius, so consecutive nodes are separated by a CHORD slightly
shorter than the arc step.  Rather than record the requested ``seg_um`` and let a law inherit the
difference, ``seg_rest_um`` is the built chord, computed on the device from the positions that exist.
The population therefore starts axially unstrained by construction, and ``strand.py``'s rule — *"the
realised length is what a downstream law must read"* — is satisfied without anyone having to remember
it.  ``seg_arc_um`` stays on the REQUESTED step, because the arc coordinate is a material coordinate
along the rest contour and is what a motor binds at; the two are different quantities.

THE POPULATION COMES FROM DENSITY x GEOMETRY.  ``n_filaments`` is never passed in.  It is
``round(areal_density_um2 x 4 pi radius_um^2)`` at the radius the shell is actually placed at, so the
count and the placement cannot disagree.  ⚠ The incumbent's 70,686 was taken at the CELL radius 7.5 µm
while the shell is placed ~h_cortex/2 inside it; deriving both from one radius is what this signature
enforces, and it changes the count.  That difference is a fact about the incumbent, reported rather
than reproduced.

WHAT IS DELIBERATELY ABSENT.  No bending modulus, no axial law, no crosslink, no motor, no steric
term, no overlap resolution.  PHASE 1 is geometry.  In particular the incumbent's build-time
interpenetration fix (``cortex_overlap_mode``) is NOT ported: it is a relaxation under a WCA potential,
which is physics, and it belongs to the phase that binds ``laws/``.  The radial dispersion here is the
h_cortex shell thickness and nothing more — a filament sits at its own radius inside the layer, which
is the physiological arrangement, not a repair for one.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``radius_um``, ``thickness_um``, ``contour_um``, ``seg_um`` and the arc coordinate
    are lengths [µm]; ``areal_density_um2`` is [µm⁻²]; ``seed`` and every index are dimensionless.
  * boundary — every length has NO default and the build refuses without it; non-finite or
    non-positive lengths and densities are refused; a discretisation giving fewer than 3 nodes is
    refused, because a filament with no bending triple is a different physical object rather than a
    coarser one; a shell thickness that would reach a non-positive radius is refused; a tangent drawn
    parallel to its own centre direction falls back to a fixed axis instead of normalising a zero.
  * conservation/invariant — three claims, contiguous, one population; every segment and angle index
    is inside this population's node range by arithmetic, and ``_demo`` checks the arithmetic against
    ``strand.build_strand`` on a small case; the arena's ``assert_partitioned`` covers the rest.
  * CFL/precision — no integration and no force.  float64 for every position and length.  The RNG
    draws are float32 because they choose a random DIRECTION, which carries no accuracy requirement;
    they are cast and re-normalised in float64, so a node sits on its shell radius to float64.
  * sign sense — ``polarity`` is ``+1`` when the barbed end is the last node of a filament, matching
    ``strand.Strand``; the tip is derived from it so the two cannot disagree.  Nothing is accumulated,
    so there is no force sign to check yet.
  * measurement protocol — construction launches kernels and reads nothing back per element.  The only
    host readbacks are the summary statistics in :meth:`StrandPopulation.record`, which are build-time
    and optional.

engine units: length µm, areal density µm⁻².  Runtime: NVIDIA Warp on CUDA.  CPU-importable; the
builder requires a CUDA arena and raises otherwise.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from dataclasses import fields as _fields

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.world.arena import Claim, Kind, WorldArena

__all__ = ["StrandPopulation", "build_strand_population", "build_cortex", "filaments_for_density"]


def filaments_for_density(areal_density_um2: float, radius_um: float) -> int:
    """Filament count from a sourced areal density over the shell it is placed on.

    ``round(rho x 4 pi R^2)``.  Both arguments are the SAME radius the filaments are laid at, which is
    the point: a count derived at one radius and a placement at another is two geometries wearing one
    name.

    Raises:
        ValueError: on a non-positive density or radius, or a combination giving fewer than one filament.
    """
    if not (math.isfinite(areal_density_um2) and areal_density_um2 > 0.0):
        raise ValueError(f"areal_density_um2 must be finite and positive; got {areal_density_um2!r}")
    if not (math.isfinite(radius_um) and radius_um > 0.0):
        raise ValueError(f"radius_um must be finite and positive; got {radius_um!r}")
    n = int(round(areal_density_um2 * 4.0 * math.pi * radius_um * radius_um))
    if n < 1:
        raise ValueError(
            f"areal density {areal_density_um2} um^-2 over 4*pi*{radius_um}^2 um^2 gives {n} filaments"
        )
    return n


@dataclass(frozen=True, slots=True)
class StrandPopulation:
    """One filament population: THREE claims for all of it, and device-resident topology.

    Attributes:
        population: the name the ranges are attributed to, e.g. ``"cortex"``.
        nodes / segments / angles: the three claims — one each, for the whole population.
        n_strands: filaments in the population, DERIVED from density x geometry.
        nodes_per_strand: uniform across the population; this uniformity is what lets a filament's
            addresses be arithmetic rather than a stored offset table.
        seg_node: ``(S, 2)`` device int32 segment endpoints, GLOBAL arena node ids.
        seg_rest_um: ``(S,)`` device float64 rest lengths [µm] — the BUILT chord, not the requested step.
        seg_arc_um: ``(S, 2)`` device float64 rest-configuration arc interval ``[s0, s1)`` [µm], on the
            requested step: this is the material coordinate a motor binds at.
        angle_idx: ``(A, 3)`` device int32 bending triples, GLOBAL ids.
        radius_um / thickness_um: the shell the filaments sit in [µm]; each filament takes its own
            radius inside ``radius +/- thickness/2``.
        centre_um: the shell centre [µm] — recorded, because a population whose position is not in its
            own record cannot be checked against another population's.
        areal_density_um2: the sourced density the count came from [µm⁻²].
        density_provenance: where that density came from.  Free text, carried into the record so a
            PI-GAP label travels inside the artifact rather than beside it.
        contour_um_requested / contour_um_realised: asked-for and built contour [µm].  These agree:
            a contour is divided into a whole number of segments, so it is the STEP that rounds.
        seg_um_requested / seg_um_realised: the discretisation step asked for and the one the rounding
            produced [µm].  Distinct again from the built CHORD, which is shorter than the arc step by
            the curvature of the shell and is what ``seg_rest_um`` holds.
        seed: the only other thing the build depends on.  ``(seed, n_strands)`` reproduces it exactly.
        polarity: ``+1`` if the barbed end is the last node of each filament.
        topology_bytes: device bytes this population's own arrays hold, EXCLUDING the arena's shared
            node arrays.
    """

    population: str
    nodes: Claim
    segments: Claim
    angles: Claim
    n_strands: int
    nodes_per_strand: int
    seg_node: wp.array = field(repr=False)
    seg_rest_um: wp.array = field(repr=False)
    seg_arc_um: wp.array = field(repr=False)
    angle_idx: wp.array = field(repr=False)
    radius_um: float
    thickness_um: float
    centre_um: tuple[float, float, float]
    areal_density_um2: float
    density_provenance: str
    contour_um_requested: float
    contour_um_realised: float
    seg_um_requested: float
    seg_um_realised: float
    seed: int
    polarity: int
    topology_bytes: int
    #: Per-filament polarity, or ``None`` when the population is genuinely uniform.  PI decision 4,
    #: 2026-08-21.  Host ``int8``; ``n_strands`` entries, so 70,686 filaments cost 70 kB and no kernel
    #: has to read it to answer a host-side question.
    polarity_per_strand: npt.NDArray[np.int8] | None = None
    #: The same array on the device, for the laws that address a barbed end inside a kernel.
    polarity_d: wp.array | None = field(default=None, repr=False)
    #: Where the per-strand polarity came from.  Free text, carried into the record so the basis
    #: travels INSIDE the artifact rather than beside it — the same rule ``density_provenance`` follows.
    polarity_basis: str = "scalar; population is uniform by construction"

    @property
    def segments_per_strand(self) -> int:
        """Segments per filament, ``N - 1``."""
        return self.nodes_per_strand - 1

    def node_id(self, strand: int, k: int) -> int:
        """GLOBAL arena node id of node ``k`` of filament ``strand``.

        The replacement for a per-filament claim: the address is arithmetic over one contiguous range.
        """
        if not 0 <= strand < self.n_strands:
            raise IndexError(f"strand {strand} outside [0, {self.n_strands})")
        if not 0 <= k < self.nodes_per_strand:
            raise IndexError(f"node {k} outside [0, {self.nodes_per_strand})")
        return self.nodes.lo + strand * self.nodes_per_strand + k

    def tip_node(self, strand: int) -> int:
        """GLOBAL id of filament ``strand``'s barbed end.

        Reads :attr:`polarity_per_strand` when the population has one and falls back to the scalar
        :attr:`polarity` only when it does not.

        ⚠ **Why the array exists at all.** Until 2026-08-21 this was the scalar for every filament, and
        the docstring read *"+1 if the barbed end is the last node of EACH filament"* — so a whole
        cortex ran one way and **an anti-parallel pair was not merely unverifiable, it was
        unrepresentable.** A bipolar minifilament straddles two ANTI-parallel filaments; in a
        uniform-polarity population every station it could take is a spring wearing a motor's name, and
        `plan_contraction` returned ``bindable: True`` without ever reading this field. That is the ERM
        defect class with a hidden DEGREE OF FREEDOM in place of a hidden count — and it is worse,
        because a filled field reads as an answer while a missing one at least looks blank.
        """
        if self.polarity_per_strand is not None:
            if not 0 <= strand < self.n_strands:
                raise IndexError(f"strand {strand} outside [0, {self.n_strands})")
            sign = int(self.polarity_per_strand[strand])
        else:
            sign = self.polarity
        return self.node_id(strand, self.nodes_per_strand - 1 if sign > 0 else 0)

    @property
    def polarity_is_mixed(self) -> bool:
        """Whether this population can hold an anti-parallel pair at all.

        A station is only a motor if its two filaments oppose. This is the question a connector has to
        ask BEFORE it places one, and before 2026-08-21 nothing could ask it.
        """
        return self.polarity_per_strand is not None and bool(
            (self.polarity_per_strand > 0).any() and (self.polarity_per_strand < 0).any())

    def antiparallel(self, a: int, b: int) -> bool:
        """Whether filaments ``a`` and ``b`` run in opposite directions.

        Raises:
            ValueError: if the population carries no per-strand polarity. **Refusing is the point** —
                a uniform population cannot answer this, and answering ``False`` would read as "these
                two happen to be parallel" rather than "this population cannot say".
        """
        if self.polarity_per_strand is None:
            raise ValueError(
                f"population {self.population!r} carries a single scalar polarity, so it cannot say "
                "whether two filaments oppose. Build it with per-strand polarity (PI decision 4, "
                "2026-08-21) before asking; a bipolar motor station is not placeable without it."
            )
        return int(self.polarity_per_strand[a]) != int(self.polarity_per_strand[b])

    def record(self, *, with_stats: bool = True) -> dict[str, object]:
        """Host-side summary for a run record.

        ``with_stats`` reads the segment rest lengths back to report their spread — build-time only,
        and skippable, because it is a full-population D2H transfer and this module refuses to pretend
        those are free.
        """
        out: dict[str, object] = {
            "population": self.population,
            "kind": "strand_population",
            "n_strands": self.n_strands,
            "nodes_per_strand": self.nodes_per_strand,
            "n_nodes": self.nodes.count,
            "n_segments": self.segments.count,
            "n_angle3": self.angles.count,
            "radius_um": self.radius_um,
            "thickness_um": self.thickness_um,
            "centre_um": list(self.centre_um),
            "areal_density_um2": self.areal_density_um2,
            "density_provenance": self.density_provenance,
            "contour_um_requested": self.contour_um_requested,
            "contour_um_realised": self.contour_um_realised,
            "seg_um_requested": self.seg_um_requested,
            "seg_um_realised_arc": self.seg_um_realised,
            "seed": self.seed,
            "polarity": self.polarity,
            "polarity_per_strand": (
                None if self.polarity_per_strand is None else {
                    "n_plus": int((self.polarity_per_strand > 0).sum()),
                    "n_minus": int((self.polarity_per_strand < 0).sum()),
                    "mixed": self.polarity_is_mixed,
                    "basis": self.polarity_basis,
                }),
            "topology_bytes": self.topology_bytes,
            "claims": {
                "node": [self.nodes.lo, self.nodes.hi],
                "segment": [self.segments.lo, self.segments.hi],
                "angle3": [self.angles.lo, self.angles.hi],
            },
        }
        if with_stats:
            rest = self.seg_rest_um.numpy()
            out["seg_chord_um_mean"] = float(rest.mean())
            out["seg_chord_um_min"] = float(rest.min())
            out["seg_chord_um_max"] = float(rest.max())
        return out


# ── kernels ─────────────────────────────────────────────────────────────────────────────────────
@wp.kernel
def _strand_positions(
    seed: wp.int32,
    n_per: wp.int32,
    seg: wp.float64,
    radius: wp.float64,
    thickness: wp.float64,
    centre: wp.vec3d,
    node_lo: wp.int32,
    range_id: wp.int32,
    pos: wp.array(dtype=wp.vec3d),
    strand_id: wp.array(dtype=wp.int32),
    rng_id: wp.array(dtype=wp.int32),
):
    """Lay every node of every filament as a great-circle arc on its own shell radius.

    One thread per NODE, not per filament: the frame is re-derived from ``wp.rand_init(seed, f)`` so
    all ``n_per`` threads of a filament agree without any of them writing a shared value.  That is what
    removes the per-filament host array the PI's construction ruling was aimed at.

    The arc is centred on the filament's own centre direction, so ``s`` runs symmetrically about it and
    the population is not biased toward one end.  ``theta = s / r`` uses the filament's OWN radius, so
    every filament realises the same arc length regardless of where in the shell thickness it sits.
    """
    i = wp.tid()
    f = i // n_per
    k = i - f * n_per

    state = wp.rand_init(seed, f)
    c32 = wp.sample_unit_sphere_surface(state)
    c = wp.normalize(wp.vec3d(wp.float64(c32[0]), wp.float64(c32[1]), wp.float64(c32[2])))
    t32 = wp.sample_unit_sphere_surface(state)
    t = wp.vec3d(wp.float64(t32[0]), wp.float64(t32[1]), wp.float64(t32[2]))
    t = t - wp.dot(t, c) * c
    # A tangent drawn parallel to its own centre direction has no tangent-plane part. It cannot
    # happen at any rate worth naming, and normalising a zero would put a NaN into the arena, so the
    # degenerate draw takes a fixed axis instead of a random one.
    if wp.length(t) < wp.float64(1.0e-9):
        t = wp.vec3d(wp.float64(1.0), wp.float64(0.0), wp.float64(0.0))
        t = t - wp.dot(t, c) * c
        if wp.length(t) < wp.float64(1.0e-9):
            t = wp.vec3d(wp.float64(0.0), wp.float64(1.0), wp.float64(0.0))
            t = t - wp.dot(t, c) * c
    t = wp.normalize(t)

    r = radius + thickness * (wp.float64(wp.randf(state)) - wp.float64(0.5))
    s = (wp.float64(k) - wp.float64(n_per - 1) * wp.float64(0.5)) * seg
    th = s / r
    pos[node_lo + i] = (c * wp.cos(th) + t * wp.sin(th)) * r + centre
    strand_id[node_lo + i] = f
    rng_id[node_lo + i] = range_id


@wp.kernel
def _strand_segments(
    n_per: wp.int32,
    n_seg_per: wp.int32,
    seg: wp.float64,
    node_lo: wp.int32,
    pos: wp.array(dtype=wp.vec3d),
    seg_node: wp.array(dtype=wp.int32, ndim=2),
    seg_rest: wp.array(dtype=wp.float64),
    seg_arc: wp.array(dtype=wp.float64, ndim=2),
):
    """Segment endpoints, rest length and material arc interval, one thread per segment.

    ``seg_rest`` is the distance between the two nodes that were actually built, so the population
    starts axially unstrained and no law inherits the arc-versus-chord difference.  ``seg_arc`` stays
    on the requested step because it is a coordinate along the material, not a distance in space.
    """
    j = wp.tid()
    f = j // n_seg_per
    k = j - f * n_seg_per
    i0 = node_lo + f * n_per + k
    seg_node[j, 0] = i0
    seg_node[j, 1] = i0 + 1
    seg_rest[j] = wp.length(pos[i0 + 1] - pos[i0])
    seg_arc[j, 0] = wp.float64(k) * seg
    seg_arc[j, 1] = wp.float64(k + 1) * seg


@wp.kernel
def _strand_angles(
    n_per: wp.int32,
    n_ang_per: wp.int32,
    node_lo: wp.int32,
    angle_idx: wp.array(dtype=wp.int32, ndim=2),
):
    """Bending triples, one thread per triple.  ``N - 2`` per filament, and never across filaments."""
    j = wp.tid()
    f = j // n_ang_per
    k = j - f * n_ang_per
    i0 = node_lo + f * n_per + k
    angle_idx[j, 0] = i0
    angle_idx[j, 1] = i0 + 1
    angle_idx[j, 2] = i0 + 2


# ── the builder ─────────────────────────────────────────────────────────────────────────────────
def build_strand_population(
    arena: WorldArena,
    population: str,
    *,
    n_strands: int,
    radius_um: float,
    thickness_um: float,
    contour_um: float,
    seg_um: float,
    areal_density_um2: float,
    density_provenance: str,
    centre_um: tuple[float, float, float] = (0.0, 0.0, 0.0),
    seed: int = 0,
    polarity: int = +1,
    polarity_mix: float | None = None,
    polarity_basis: str | None = None,
) -> StrandPopulation:
    """Claim ONE range per kind and build ``n_strands`` shell filaments into ``arena`` on the device.

    Callers normally reach this through :func:`build_cortex`, which derives ``n_strands`` rather than
    accepting it.  It is separate so a second uniform shell population — an Arp2/3 short-filament
    cortex, for instance — can be built without either one being able to type its own count.

    Args:
        arena: the world to claim from.  Must hold a CUDA allocation.
        population: the name the ranges are attributed to.
        n_strands: filaments to build.  Derived by the caller from density x geometry.
        radius_um: mid-shell radius the filaments are laid on [µm].
        thickness_um: shell thickness [µm]; each filament takes its own radius in
            ``radius +/- thickness/2``, which is the layer being a layer rather than a sphere.
        contour_um: rest contour length of each filament [µm].
        seg_um: requested discretisation step [µm].
        areal_density_um2: the density ``n_strands`` came from — recorded, not re-used.
        density_provenance: where that density came from.
        centre_um: shell centre [µm].
        seed: the RNG seed.  ``(seed, n_strands)`` reproduces the build exactly.
        polarity: ``+1`` if the barbed end is the last node of each filament.

    Returns:
        The built :class:`StrandPopulation`.

    Raises:
        ValueError: on a non-positive count or length, a thickness reaching a non-positive radius, a
            polarity that is not ±1, or a discretisation giving fewer than 3 nodes per filament.
        RuntimeError: if the arena holds no device allocation.
    """
    if n_strands < 1:
        raise ValueError(f"n_strands must be at least 1; got {n_strands}")
    if polarity not in (+1, -1):
        raise ValueError(f"polarity must be +1 or -1; got {polarity!r}")
    for name, value in (("radius_um", radius_um), ("contour_um", contour_um), ("seg_um", seg_um)):
        if not (math.isfinite(value) and value > 0.0):
            raise ValueError(f"{name} must be finite and positive; got {value!r}")
    if not (math.isfinite(thickness_um) and thickness_um >= 0.0):
        raise ValueError(f"thickness_um must be finite and nonnegative; got {thickness_um!r}")
    if radius_um - 0.5 * thickness_um <= 0.0:
        raise ValueError(
            f"a shell of thickness {thickness_um} um about radius {radius_um} um reaches the centre; "
            "the inner face of the layer would be at or through the origin."
        )

    n_seg_per = int(round(contour_um / seg_um))
    n_per = n_seg_per + 1
    if n_per < 3:
        raise ValueError(
            f"contour_um={contour_um} at seg_um={seg_um} gives {n_per} nodes per filament; a strand "
            "needs at least 3 so it carries one bending triple. A filament with no bending stiffness "
            "is a different physical object, not a coarser one."
        )
    n_ang_per = n_per - 2
    if not arena.node_arrays:
        raise RuntimeError(
            f"build_strand_population({population}) needs a CUDA arena — this one was built with "
            "device=None, the bookkeeping-only mode. Construction is a Warp kernel (PI 2026-08-20); "
            "there is no host build path to fall back to."
        )
    contour_realised = n_seg_per * (contour_um / n_seg_per)  # exact; the step is what rounds, not the contour

    n_nodes = n_strands * n_per
    n_segments = n_strands * n_seg_per
    n_angles = n_strands * n_ang_per

    nodes = arena.claim(population, Kind.NODE, n_nodes)
    segments = arena.claim(population, Kind.SEGMENT, n_segments)
    angles = arena.claim(population, Kind.ANGLE3, n_angles)
    range_id = len(arena.claims(kind=Kind.NODE)) - 1

    device = arena.node_arrays["position"].device
    seg_step = contour_um / n_seg_per
    with wp.ScopedDevice(device):
        wp.launch(
            _strand_positions,
            dim=n_nodes,
            inputs=[
                int(seed), n_per, wp.float64(seg_step), wp.float64(radius_um),
                wp.float64(thickness_um), wp.vec3d(*(float(x) for x in centre_um)),
                nodes.lo, range_id,
                arena.node_arrays["position"], arena.node_arrays["strand_id"],
                arena.node_arrays["range_id"],
            ],
        )
        seg_node = wp.zeros((n_segments, 2), dtype=wp.int32)
        seg_rest = wp.zeros(n_segments, dtype=wp.float64)
        seg_arc = wp.zeros((n_segments, 2), dtype=wp.float64)
        wp.launch(
            _strand_segments,
            dim=n_segments,
            inputs=[
                n_per, n_seg_per, wp.float64(seg_step), nodes.lo,
                arena.node_arrays["position"], seg_node, seg_rest, seg_arc,
            ],
        )
        angle_idx = wp.zeros((n_angles, 3), dtype=wp.int32)
        wp.launch(_strand_angles, dim=n_angles, inputs=[n_per, n_ang_per, nodes.lo, angle_idx])

    # ponytail: seg_node, seg_arc and angle_idx are MATERIALISED, and all three are derivable from
    # (j, n_per, n_seg_per, nodes.lo) by the same arithmetic node_id() uses — at native that is ~148 MB
    # of the cortex's ~181 MB of topology. They are stored because PHASE 2's strategy is to BIND
    # laws/ unchanged, and those kernels take index ARRAYS; computing the indices in-kernel would mean
    # forking every law. Upgrade path if the layout ever binds on memory: pass (n_per, n_seg_per, lo)
    # and index arithmetically, one law at a time, measuring each.
    topology_bytes = int(seg_node.size * 4 + seg_rest.size * 8 + seg_arc.size * 8 + angle_idx.size * 4)
    # ── per-filament polarity, PI decision 4 (2026-08-21) ────────────────────────────────────────
    # ⚠ Until this landed, `polarity` was ONE int for the whole population, documented as "+1 if the
    # barbed end is the last node of EACH filament". So every filament ran the same way and an
    # ANTI-PARALLEL PAIR WAS UNREPRESENTABLE — not unverifiable, unrepresentable. A bipolar NMII
    # minifilament straddles two anti-parallel filaments; in a uniform population every station it
    # could take is a spring wearing a motor's name, and `plan_contraction` returned bindable: True
    # without reading the field at all.
    #
    # `polarity_mix` is the fraction pointing +1 and it has NO DEFAULT: a caller either states the
    # architecture or gets the old uniform behaviour explicitly. 0.5 is not a convenience — it is the
    # ABSENCE OF RECTIFICATION, which is what a cortex is. The contrast is written into this repo
    # already: `build/stress_fiber.py` quotes architecture_spec's *"this rectification is REQUIRED,
    # random mixed polarity does NOT contract"* for a stress fibre. A cortex is the other case, and
    # that is why cortical NMII generates tension at all.
    pol_host: npt.NDArray[np.int8] | None = None
    pol_d = None
    basis = polarity_basis or "scalar; population is uniform by construction"
    if polarity_mix is not None:
        if not 0.0 <= polarity_mix <= 1.0:
            raise ValueError(f"polarity_mix is a fraction in [0, 1]; got {polarity_mix!r}")
        rng = np.random.default_rng(int(seed))
        pol_host = np.where(rng.random(n_strands) < polarity_mix, 1, -1).astype(np.int8)
        pol_d = wp.array(pol_host.astype(np.int32), dtype=wp.int32)
        if polarity_basis is None:
            raise ValueError(
                "polarity_mix was given without polarity_basis. A mixed-polarity population is an "
                "architectural claim about the structure, so the artifact has to say where the claim "
                "came from — the same rule density_provenance follows, and for the same reason: the "
                "ERM count was a mesh artefact because nobody had to say."
            )

    return StrandPopulation(
        polarity_per_strand=pol_host, polarity_d=pol_d, polarity_basis=basis,
        population=population, nodes=nodes, segments=segments, angles=angles,
        n_strands=n_strands, nodes_per_strand=n_per,
        seg_node=seg_node, seg_rest_um=seg_rest, seg_arc_um=seg_arc, angle_idx=angle_idx,
        radius_um=float(radius_um), thickness_um=float(thickness_um),
        centre_um=tuple(float(x) for x in centre_um),
        areal_density_um2=float(areal_density_um2), density_provenance=density_provenance,
        contour_um_requested=float(contour_um), contour_um_realised=float(contour_realised),
        seg_um_requested=float(seg_um), seg_um_realised=float(seg_step),
        seed=int(seed), polarity=int(polarity),
        topology_bytes=topology_bytes,
    )


def build_cortex(
    arena: WorldArena,
    *,
    radius_um: float | None = None,
    thickness_um: float | None = None,
    areal_density_um2: float | None = None,
    contour_um: float | None = None,
    seg_um: float | None = None,
    density_provenance: str = "",
    centre_um: tuple[float, float, float] = (0.0, 0.0, 0.0),
    seed: int = 0,
    population: str = "cortex",
    polarity_mix: float | None = 0.5,
) -> StrandPopulation:
    """The cortical actin shell, with its filament count DERIVED from density x geometry.

    Every argument that carries physiology has NO default and the build refuses without it.  That is
    the ratified pattern (``strand.build_strand``, ``MembraneAreaCard.capacity``) and it matters most
    here: ``seg_um = 0.5`` is the constant the provenance audit ranks second most load-bearing in the
    engine with a source of *"none (0 KB rows)"*, and a default in this signature would put it back
    under a new name.

    Args:
        arena: the world to claim from.
        radius_um: MID-SHELL radius [µm].  **REQUIRED.**  The cortex is a layer just inside the
            membrane, so this is the cell radius minus about half the cortex thickness — not the cell
            radius.  It is also the radius the filament count is derived over, so the count and the
            placement are the same geometry.
        thickness_um: cortex thickness [µm].  **REQUIRED.**
        areal_density_um2: filament areal density [µm⁻²].  **REQUIRED.**
        contour_um: filament rest contour length [µm].  **REQUIRED.**
        seg_um: discretisation step [µm].  **REQUIRED.**
        density_provenance: where ``areal_density_um2`` came from.
        centre_um: cell centre [µm].
        seed: RNG seed for the filament frames.
        population: the claim name.

    Returns:
        The built :class:`StrandPopulation`.

    Raises:
        ValueError: if any required value is missing, or on any refusal from
            :func:`build_strand_population` / :func:`filaments_for_density`.
    """
    missing = [
        name for name, value in (
            ("radius_um", radius_um), ("thickness_um", thickness_um),
            ("areal_density_um2", areal_density_um2), ("contour_um", contour_um), ("seg_um", seg_um),
        ) if value is None
    ]
    if missing:
        raise ValueError(
            f"{missing} have no defaults and must be declared. Every value labelled physiological is "
            "physiological for some particular cell and is not certain, so it arrives as a swept axis "
            "with a scope — never as a literal chosen here. seg_um in particular: the audit ranks the "
            "incumbent's 0.5 um second most load-bearing with a source of 'none (0 KB rows)'."
        )
    return build_strand_population(
        arena, population,
        n_strands=filaments_for_density(areal_density_um2, radius_um),
        radius_um=radius_um, thickness_um=thickness_um, contour_um=contour_um, seg_um=seg_um,
        areal_density_um2=areal_density_um2, density_provenance=density_provenance,
        centre_um=centre_um, seed=seed,
        polarity_mix=polarity_mix,
        polarity_basis=(
            None if polarity_mix is None else
            "UNRECTIFIED — a cortex has no mechanism that aligns filament polarity, so the mix is "
            "isotropic. This is not a convenience default: the contrast is written into this repo "
            "already, at build/stress_fiber.py, which quotes architecture_spec's 'this rectification "
            "is REQUIRED, random mixed polarity does NOT contract' for a VENTRAL FIBRE. The cortex is "
            "the other case, and an unrectified mix is the premise under which cortical NMII "
            "generates tension at all. PI decision 4, 2026-08-21."
        ),
    )


# ── self-check ──────────────────────────────────────────────────────────────────────────────────
def _demo() -> None:
    """Self-check: the derivation, the address arithmetic, and every refusal.

    Everything a device would compute is checked at the level a device is not needed for — the counts,
    the addressing that replaced the per-strand claim, and the refusals.  The kernels themselves are
    checked against ``strand.build_strand`` on the GPU host, where they can run.
    """
    from aleph.world.strand import build_strand

    # The count is derived, and it reproduces the incumbent's 70,686 at the CELL radius.
    assert filaments_for_density(100.0, 7.5) == 70_686, "100 um^-2 x 4 pi 7.5^2"
    # ...and does NOT at the shell radius, which is where the filaments actually are. Reported, not hidden.
    assert filaments_for_density(100.0, 7.4) == 68_813
    assert filaments_for_density(100.0, 7.4) < filaments_for_density(100.0, 7.5)

    # ── PI decision 4: a population must be able to say whether two filaments oppose ──────────────
    import numpy as _np
    uniform = StrandPopulation(
        population="u", nodes=Claim("u", Kind.NODE, 0, 40), segments=Claim("u", Kind.SEGMENT, 0, 36),
        angles=Claim("u", Kind.ANGLE3, 0, 32), n_strands=4, nodes_per_strand=10,
        seg_node=None, seg_rest_um=None, seg_arc_um=None, angle_idx=None,
        radius_um=7.4, thickness_um=0.2, centre_um=(0.0, 0.0, 0.0), areal_density_um2=100.0,
        density_provenance="self-check", contour_um_requested=1.0, contour_um_realised=1.0,
        seg_um_requested=0.25, seg_um_realised=0.25, seed=0, polarity=+1, topology_bytes=0)
    mixed = StrandPopulation(
        **{**{f.name: getattr(uniform, f.name) for f in _fields(uniform)},
           "polarity_per_strand": _np.array([1, -1, 1, -1], dtype=_np.int8),
           "polarity_basis": "self-check"})

    # A uniform population REFUSES the question rather than answering False. Answering False would
    # read as "these two happen to be parallel"; the truth is "this population cannot say".
    assert not uniform.polarity_is_mixed
    try:
        uniform.antiparallel(0, 1)
    except ValueError as exc:
        assert "cannot say" in str(exc)
    else:
        raise AssertionError("a uniform population must refuse antiparallel(), not answer it")

    assert mixed.polarity_is_mixed
    assert mixed.antiparallel(0, 1) and not mixed.antiparallel(0, 2)
    # and the tip follows the per-strand sign, which is the whole point
    assert mixed.tip_node(0) == mixed.node_id(0, 9), "+1 puts the barbed end last"
    assert mixed.tip_node(1) == mixed.node_id(1, 0), "-1 puts it first"
    # a population that is all one way is NOT mixed even with an array present
    all_plus = StrandPopulation(
        **{**{f.name: getattr(uniform, f.name) for f in _fields(uniform)},
           "polarity_per_strand": _np.ones(4, dtype=_np.int8), "polarity_basis": "self-check"})
    assert not all_plus.polarity_is_mixed, "an array of one sign holds no anti-parallel pair either"

    # Nodes per filament follow N = L/seg + 1, which is not proportional to 1/seg.
    for seg, expect in ((0.5, 7), (0.05, 61)):
        assert int(round(3.0 / seg)) + 1 == expect
    assert 70_686 * 7 == 494_802, "the incumbent's node count, reproduced from the same arithmetic"

    # The address arithmetic that replaced the per-strand claim, against the arena's own bookkeeping.
    arena = WorldArena(capacity={Kind.NODE: 1_000, Kind.SEGMENT: 1_000, Kind.ANGLE3: 1_000})
    n_per, n_str = 5, 7
    nodes = arena.claim("cortex", Kind.NODE, n_str * n_per)
    segments = arena.claim("cortex", Kind.SEGMENT, n_str * (n_per - 1))
    angles = arena.claim("cortex", Kind.ANGLE3, n_str * (n_per - 2))
    pop = StrandPopulation(
        population="cortex", nodes=nodes, segments=segments, angles=angles,
        n_strands=n_str, nodes_per_strand=n_per,
        seg_node=None, seg_rest_um=None, seg_arc_um=None, angle_idx=None,
        radius_um=7.4, thickness_um=0.2, centre_um=(0.0, 0.0, 0.0),
        areal_density_um2=100.0, density_provenance="demo",
        contour_um_requested=1.0, contour_um_realised=1.0, seg_um_requested=0.25,
        seg_um_realised=0.25, seed=0, polarity=+1, topology_bytes=0,
    )
    assert pop.node_id(0, 0) == nodes.lo
    assert pop.node_id(n_str - 1, n_per - 1) == nodes.hi - 1, "the last node is the end of the range"
    assert pop.tip_node(3) == pop.node_id(3, n_per - 1), "polarity +1 puts the barbed end last"
    assert pop.segments_per_strand == n_per - 1
    # Consecutive filaments do not overlap and leave no gap — the claim is one contiguous range.
    assert pop.node_id(1, 0) == pop.node_id(0, n_per - 1) + 1
    for bad in ((n_str, 0), (0, n_per), (-1, 0)):
        try:
            pop.node_id(*bad)
        except IndexError:
            pass
        else:  # pragma: no cover
            raise AssertionError(f"node_id{bad} must refuse")
    arena.assert_partitioned()

    # ONE claim per kind for the whole population, against build_strand's three PER STRAND.
    assert len(arena.claims(population="cortex")) == 3
    per_strand = WorldArena(capacity={Kind.NODE: 1_000, Kind.SEGMENT: 1_000, Kind.ANGLE3: 1_000})
    for i in range(n_str):
        build_strand(per_strand, "cortex", start=(0.0, float(i), 0.0), direction=(1.0, 0.0, 0.0),
                     contour_um=1.0, seg_um=0.25)
    assert len(per_strand.claims(population="cortex")) == 3 * n_str
    # Same nodes, same segments, same triples — only the claim COUNT differs.
    assert per_strand.n_live(Kind.NODE) == arena.n_live(Kind.NODE)
    assert per_strand.n_live(Kind.SEGMENT) == arena.n_live(Kind.SEGMENT)
    assert per_strand.n_live(Kind.ANGLE3) == arena.n_live(Kind.ANGLE3)
    # At native this is the difference the phase exists to remove.
    assert 3 * 70_686 == 212_058, "the per-strand claim count the plan names"

    # Every physiological length refuses to default.
    live = WorldArena(capacity={Kind.NODE: 100, Kind.SEGMENT: 100, Kind.ANGLE3: 100})
    try:
        build_cortex(live, radius_um=7.4, thickness_um=0.2, areal_density_um2=100.0, contour_um=3.0)
    except ValueError as exc:
        assert "seg_um" in str(exc) and "no defaults" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a missing seg_um must refuse")
    try:
        build_cortex(live)
    except ValueError as exc:
        assert all(k in str(exc) for k in ("radius_um", "thickness_um", "seg_um"))
    else:  # pragma: no cover
        raise AssertionError("build_cortex() must refuse")

    # A bookkeeping-only arena refuses to build: there is no host construction path.
    try:
        build_cortex(live, radius_um=7.4, thickness_um=0.2, areal_density_um2=1.0e-3,
                     contour_um=3.0, seg_um=1.0)
    except RuntimeError as exc:
        assert "needs a CUDA arena" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device=None arena must refuse to build")

    # A discretisation with no bending triple is a different object, and a shell through the origin
    # is not a shell. Both refuse before anything is claimed.
    device_free = WorldArena(capacity={Kind.NODE: 100})
    for kwargs, needle in (
        ({"contour_um": 1.0, "seg_um": 1.0}, "at least 3"),
        ({"contour_um": 3.0, "seg_um": 0.05, "thickness_um": 20.0}, "reaches the centre"),
    ):
        call = {"radius_um": 7.4, "thickness_um": 0.2, "areal_density_um2": 100.0} | kwargs
        try:
            build_cortex(device_free, **call)
        except (ValueError, RuntimeError) as exc:
            assert needle in str(exc) or "needs a CUDA arena" in str(exc)

    n_fil = filaments_for_density(100.0, 7.4)
    print(
        f"cortex self-check OK — density 100 um^-2 at R=7.4 um -> {n_fil:,} filaments; "
        f"seg 0.05 um over 3.0 um -> 61 nodes each -> {n_fil * 61:,} nodes, "
        f"{3} claims (was {3 * n_fil:,})"
    )


if __name__ == "__main__":
    _demo()
