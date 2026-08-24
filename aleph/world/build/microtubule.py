r"""The microtubule aster, built into the arena by a CUDA kernel that computes its own geometry.

WHAT IS DIFFERENT FROM THE BUILDER THIS REPLACES.  ``components/solid/microtubule.py`` builds the aster
on the HOST — ``_fibonacci_directions`` in NumPy, a Python list comprehension per arm, then
``build_fiber_network`` — and uploads the result.  At the native cortex segment length that pattern
costs a full host mirror of every population before the device sees one node, and the mirror is exactly
what ``CLAUDE.md`` forbids keeping.  Here the arm direction, the bead spacing and the node position are
COMPUTED IN THE KERNEL from the node's own index.  Nothing is uploaded: the only host→device traffic is
a handful of scalars in the launch record.

ONE CLAIM PER POPULATION, NOT PER STRAND.  ``build_strand`` claims three ranges per filament, so a
native cortex would stand at 212,058 claims and every ``population_of`` query would walk them.  A
population whose strands are all the same length does not need that: it is ONE contiguous node range,
and the strand a node belongs to is ``(index - lo) // nodes_per_strand``.  This module claims four
ranges — NODE, STRAND, SEGMENT, ANGLE3 — for the whole aster, and that arithmetic is the topology.

**THE TOPOLOGY IS ARITHMETIC, SO THERE IS NO TOPOLOGY ARRAY.**  Segment ``s`` of strand ``k`` joins
nodes ``lo + k*n_per + i`` and ``+1``; angle triple ``a`` is the same with three.  A law kernel
reconstructs that from two integers, which is why this module claims SEGMENT and ANGLE3 ranges — the
COUNTS are real and the arena must reserve them — while allocating no index arrays for them.  A
population whose strands differ in length would need explicit arrays; none of the five compartments
this session owns is such a population, and the moment one is, that is when the arrays get claimed.

THE MTOC IS A SEPARATE NODE.  Every arm's inner end is one segment out from it and the hub is a BOND
(PHASE 3), never a shared node.  Arms sharing a node would be a weld, and the incumbent already learned
this the other way round: ``laws/microtubule.merge_aster_into_cortex`` records that the MTOC is held by
finite-rest crosslinks "NOT a zero-length weld (which ``implicit_ff`` drops at L<1e-9)".

⚠ **THE COUNT IS A PI-GAP AND THIS MODULE REFUSES TO INVENT ONE.**  ``n_mt`` has no default here and
none is reachable: the incumbent's own docstring records "the MT COUNT per MCF7 cell is MCF7-absent —
proxy ~250–600 MTs/epithelial cell (order-of-magnitude); ``n_mt`` here is a provisional geometry, not a
sourced production count", and the value it ships is 40.  So the caller supplies a
:class:`~aleph.world.bond.BondCount`, which cannot be constructed without a scope and a provenance, and
the artifact carries ``PI_GAP`` where that is the truth.  ``seg_um`` likewise has no default: the
sourced 0.05 µm is a CORTICAL mesh figure and says nothing about the discretisation of a 25 nm tube.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``centre``, ``reach_R_um``, ``seg_um`` [µm]; the count is dimensionless; nothing here
    carries a force, a stiffness or a time.
  * boundary — an arm shorter than 3 beads carries no bending triple and is refused, the same rule and
    the same message as ``build_strand``; a non-positive or non-finite length is refused; a count that
    resolves to zero strands is refused, since an empty range reads downstream as "present but silent".
  * conservation/invariant — the four claims are contiguous and the kernel writes only inside the NODE
    claim: every global index it forms is ``lo + t`` for ``t`` in ``[0, n_nodes)`` and the launch
    dimension IS ``n_nodes``, so out-of-range addressing is unreachable rather than checked.
    ``arena.assert_partitioned()`` is asserted after the build.
  * CFL/precision — no integration; float64 throughout, matching the arena's ``vec3d`` arrays.  The
    Fibonacci direction is evaluated in float64 so two arms 500 apart are still distinguishable.
  * sign sense — the aster radiates OUTWARD: bead ``i`` sits at ``(i+1)·seg_um`` from the MTOC, so the
    tip is the last node of each strand and the polarity convention matches ``Strand.polarity = +1``.
  * measurement protocol — the build launches one kernel of ``dim=n_nodes`` plus one of ``dim=1`` for
    the MTOC, and reads nothing back.  Node counts are host bookkeeping; any byte figure must come from
    the driver, not from this module's arithmetic.

engine units: length µm.  Runtime: NVIDIA Warp on CUDA; the build refuses a device-less arena.
"""

from __future__ import annotations

import warp as wp

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount
# The refusal message for a missing length lives in ONE place on purpose: `test_layer_directions`'
# own argument, applied to prose — a second spelling of "has no default" is a second thing to keep
# true, and the tests that assert on it would pass while the copy drifted.
from aleph.world.strand import _require_length

__all__ = ["plan_microtubules", "build_microtubules", "lay_aster_kernel", "lay_point_kernel"]

#: ``π(1+√5)`` — the golden-angle azimuth step of the Fibonacci sphere, matching
#: ``laws/microtubule._fibonacci_directions`` exactly so the two builders place arm ``k`` in the same
#: direction. A different spiral would be a different aster, and the A/B against the incumbent in
#: PHASE 2 compares forces node by node.
GOLDEN_AZIMUTH = wp.constant(wp.float64(10.166407384630519))


@wp.kernel
def lay_aster_kernel(
    lo: wp.int32, n_per: wp.int32, n_strand: wp.int32, strand_base: wp.int32, range_ordinal: wp.int32,
    cx: wp.float64, cy: wp.float64, cz0: wp.float64, seg_um: wp.float64,
    position: wp.array(dtype=wp.vec3d), strand_id: wp.array(dtype=wp.int32),
    range_id: wp.array(dtype=wp.int32),
):
    """Write one aster arm bead per thread, deriving its arm direction from its own index.

    The thread id ``t`` runs over the arm beads only (the MTOC is written separately). ``k = t / n_per``
    is the arm and ``i`` the bead along it, using C integer division rather than ``//`` or ``%`` so the
    arithmetic is the same in every Warp version this has to build under.
    """
    t = wp.tid()
    k = t / n_per
    i = t - k * n_per

    kd = wp.float64(k) + wp.float64(0.5)
    cosp = wp.float64(1.0) - wp.float64(2.0) * kd / wp.float64(n_strand)
    sinp = wp.sqrt(wp.max(wp.float64(0.0), wp.float64(1.0) - cosp * cosp))
    az = GOLDEN_AZIMUTH * kd
    dx = sinp * wp.cos(az)
    dy = sinp * wp.sin(az)

    r = wp.float64(i + 1) * seg_um
    g = lo + t
    position[g] = wp.vec3d(cx + r * dx, cy + r * dy, cz0 + r * cosp)
    strand_id[g] = strand_base + k
    range_id[g] = range_ordinal


@wp.kernel
def lay_point_kernel(
    g0: wp.int32, strand_tag: wp.int32, range_ordinal: wp.int32,
    px: wp.float64, py: wp.float64, pz: wp.float64,
    position: wp.array(dtype=wp.vec3d), strand_id: wp.array(dtype=wp.int32),
    range_id: wp.array(dtype=wp.int32),
):
    """Write one named point — here the MTOC, which belongs to no arm and so carries ``strand_tag``.

    A one-thread launch rather than a branch inside :func:`lay_aster_kernel`: the hub is one node out of
    hundreds of thousands, and a per-thread predicate that is false everywhere but once costs the whole
    launch its uniformity for no gain.
    """
    if wp.tid() == 0:
        position[g0] = wp.vec3d(px, py, pz)
        strand_id[g0] = strand_tag
        range_id[g0] = range_ordinal


def plan_microtubules(*, count: BondCount, reach_R_um: float | None = None,
                      seg_um: float | None = None) -> dict[str, object]:
    """Resolve the aster's inventory from its declared count and geometry, touching no device.

    Separated from :func:`build_microtubules` so the arithmetic and every refusal can be exercised
    without a card, which is the only part of this module a CPU machine is allowed to run.

    Args:
        count: how many microtubules, for which cell, on whose authority.  Basis ``explicit`` (a count
            per cell) or ``volumetric`` (a density per µm³ of cytoplasm) — with ``volumetric`` the
            support is supplied at build time.  **REQUIRED**, and it cannot be built without a scope
            and a provenance, which is what keeps a proxy from being recorded as a measurement.
        reach_R_um: the radius the arm tips reach [µm].  **REQUIRED** — no default.  Tips short of the
            cortex leave the aster mechanically disengaged, which is the defect
            ``INTERNAL_DISPLACEMENT_DIAGNOSIS_2026-07-16`` recorded.
        seg_um: bead spacing along the tube [µm].  **REQUIRED** — no default.  The sourced 0.05 µm is a
            CORTICAL mesh size and is not a statement about discretising a microtubule.

    Returns:
        The inventory: ``n_strands``, ``nodes_per_strand``, ``n_nodes`` (arms + the MTOC), ``n_segments``,
        ``n_angles``, the realised segment length, and the count's provenance row.

    Raises:
        ValueError: on a missing or non-positive length, or a discretisation giving fewer than 3 beads
            per arm — an arm with no bending triple is a different physical object, not a coarser one.
    """
    reach = _require_length("reach_R_um", reach_R_um)
    seg = _require_length("seg_um", seg_um)
    n_per = int(round(reach / seg))
    if n_per < 3:
        raise ValueError(
            f"reach_R_um={reach} at seg_um={seg} gives {n_per} beads per arm; an arm needs at least 3 so "
            "it carries one bending triple. A tube with no bending stiffness is a different physical "
            "object, not a coarser one — and bending is the ONLY force a microtubule contributes."
        )
    return {
        "kind": "microtubule",
        "count_basis": count.basis,
        "nodes_per_strand": n_per,
        "seg_um_requested": seg,
        "seg_um_realised": reach / n_per,
        "reach_R_um": reach,
        "provenance": {"n_microtubules": _provenance_row(count)},
    }


def _provenance_row(count: BondCount) -> dict[str, object]:
    """The count with its scope and authority attached, in the shape ``BondFamily`` uses.

    Together on purpose, and for the reason ``bond.py`` gives: a count separated from its scope is how
    a provisional geometry acquires a physiological label.
    """
    return {
        "basis": count.basis, "value": count.value, "support_unit": count.support_unit,
        "scope": count.scope, "source_class": str(count.source_class), "provenance": count.provenance,
    }


def build_microtubules(
    arena: WorldArena,
    *,
    count: BondCount,
    support: float = 1.0,
    centre: tuple[float, float, float],
    reach_R_um: float | None = None,
    seg_um: float | None = None,
    population: str = "microtubule",
) -> dict[str, object]:
    """Claim the aster's ranges and lay its geometry down with a CUDA kernel.

    Args:
        arena: the world to claim from.  Must hold a CUDA allocation.
        count: see :func:`plan_microtubules`.
        support: the MEASURED support the count resolves against — 1.0 for an ``explicit`` count, the
            cytoplasmic volume [µm³] for a ``volumetric`` one.  Measured, never analytic.
        centre: the MTOC position [µm].
        reach_R_um: see :func:`plan_microtubules`.  **REQUIRED**.
        seg_um: see :func:`plan_microtubules`.  **REQUIRED**.
        population: the name the four claims are attributed to.

    Returns:
        The census of what was built: the inventory from :func:`plan_microtubules` plus the claims.

    Raises:
        RuntimeError: if ``arena`` holds no device allocation.  There is no CPU construction path here
            and there must never be one — a geometry laid down on the host is a host mirror, and the
            residency rule exists because the mirror is what gets stepped by accident.
        ValueError: on a refusal from :func:`plan_microtubules`, or a count resolving to zero strands.
    """
    inventory = plan_microtubules(count=count, reach_R_um=reach_R_um, seg_um=seg_um)
    n_strand = count.resolve(support)
    if n_strand < 1:
        raise ValueError(
            f"{population}: the count resolves to {n_strand} microtubules against support={support}. "
            "A population that is 'present but has nothing' is indistinguishable downstream from one "
            "that failed to build, and the two are different facts."
        )
    n_per = int(inventory["nodes_per_strand"])
    if not arena.node_arrays:
        raise RuntimeError(
            f"{population}: this arena holds no device allocation, and construction here is a CUDA "
            "kernel rather than a host build that is uploaded. Warp CUDA is the only runtime; build the "
            f"arena with a CUDA device. The inventory it would have built: {n_strand} tubes x {n_per} "
            f"beads + 1 MTOC = {n_strand * n_per + 1} nodes."
        )

    n_arm_nodes = n_strand * n_per
    nodes = arena.claim(population, Kind.NODE, n_arm_nodes + 1)   # +1: the MTOC is its own node
    strands = arena.claim(population, Kind.STRAND, n_strand)
    segments = arena.claim(population, Kind.SEGMENT, n_strand * (n_per - 1))
    angles = arena.claim(population, Kind.ANGLE3, n_strand * (n_per - 2))
    ordinal = len(arena.claims()) - 4

    seg_realised = float(inventory["seg_um_realised"])
    wp.launch(
        lay_aster_kernel, dim=n_arm_nodes,
        inputs=[nodes.lo, n_per, n_strand, strands.lo, ordinal,
                float(centre[0]), float(centre[1]), float(centre[2]), seg_realised],
        outputs=[arena.node_arrays["position"], arena.node_arrays["strand_id"],
                 arena.node_arrays["range_id"]],
        device=arena.device,
    )
    wp.launch(
        lay_point_kernel, dim=1,
        inputs=[nodes.lo + n_arm_nodes, -1, ordinal,
                float(centre[0]), float(centre[1]), float(centre[2])],
        outputs=[arena.node_arrays["position"], arena.node_arrays["strand_id"],
                 arena.node_arrays["range_id"]],
        device=arena.device,
    )
    arena.assert_partitioned()

    return inventory | {
        "population": population, "n_strands": n_strand, "n_nodes": nodes.count,
        "n_segments": segments.count, "n_angles": angles.count, "mtoc_node": nodes.lo + n_arm_nodes,
        "claims": {"node": (nodes.lo, nodes.hi), "strand": (strands.lo, strands.hi),
                   "segment": (segments.lo, segments.hi), "angle3": (angles.lo, angles.hi)},
        "support": float(support),
    }


def _codegen_check() -> int:
    """Type-check this module's kernels by emitting their CUDA source, with no device and no execution.

    The dev machine has no card, so the alternative is shipping kernel syntax that
    nobody has put in front of a compiler; this generates the CUDA C++ Warp would compile, which type-
    checks every expression, and produces NO number — so it is not a CPU result under any reading.
    """
    from warp._src.context import ModuleBuilder  # private, deliberately: there is no public codegen API

    module = wp.get_module(__name__)
    return len(ModuleBuilder(module, module.options).codegen("cuda"))


def _demo() -> None:
    """Self-check: the inventory arithmetic and every refusal, none of which needs a device."""
    # The count the incumbent ships is 40, against a proxy band of 250-600. Neither is a measurement of
    # this cell, so the class is what says so — and the class is what the artifact will carry.
    gap = BondCount(basis="explicit", value=500.0, scope="MCF7 epithelial, interphase — PROXY, the "
                    "MT count for MCF7 is absent from the KB",
                    source_class="PI_GAP",
                    provenance="components/solid/microtubule.py docstring: '~250-600 MTs/epithelial "
                               "cell (order-of-magnitude)'; the shipped n_mt=40 is 6-15x under it")

    inv = plan_microtubules(count=gap, reach_R_um=7.4, seg_um=0.05)
    assert inv["nodes_per_strand"] == 148, inv["nodes_per_strand"]
    assert abs(inv["seg_um_realised"] - 0.05) < 1e-15
    assert inv["provenance"]["n_microtubules"]["source_class"] == "PI_GAP"

    # 7.4 at 0.5 (the incumbent's spacing) does not divide evenly: 15 beads of 0.49333 um.
    coarse = plan_microtubules(count=gap, reach_R_um=7.4, seg_um=0.5)
    assert coarse["nodes_per_strand"] == 15
    assert coarse["seg_um_realised"] != coarse["seg_um_requested"]

    # Both lengths refuse to default, in the one spelling the codebase has for it.
    for kwargs in ({"reach_R_um": 7.4}, {"seg_um": 0.05}, {}):
        try:
            plan_microtubules(count=gap, **kwargs)
        except ValueError as exc:
            assert "no default" in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"plan_microtubules{kwargs} must refuse")

    # Too coarse to carry a bending triple is a different object, not a coarser one.
    try:
        plan_microtubules(count=gap, reach_R_um=1.0, seg_um=0.5)
    except ValueError as exc:
        assert "at least 3" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a 2-bead arm must refuse")

    # A count with no scope cannot be constructed at all — the refusal is upstream of this module.
    try:
        BondCount(basis="explicit", value=40.0, scope="", source_class="PI_GAP", provenance="x")
    except ValueError as exc:
        assert "scope is required" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a scopeless count must refuse")

    # And a device-less arena refuses the BUILD while still reporting what it would have built.
    arena = WorldArena(capacity={Kind.NODE: 10, Kind.STRAND: 10, Kind.SEGMENT: 10, Kind.ANGLE3: 10})
    try:
        build_microtubules(arena, count=gap, centre=(0.0, 0.0, 0.0), reach_R_um=7.4, seg_um=0.05)
    except RuntimeError as exc:
        assert "no device allocation" in str(exc) and "74001 nodes" in str(exc), str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device-less build must refuse")
    assert arena.n_live(Kind.NODE) == 0, "a refused build claims nothing"

    chars = _codegen_check()
    assert chars > 0

    print(f"microtubule self-check OK — CUDA codegen {chars} chars; {500} tubes x {inv['nodes_per_strand']} beads + 1 MTOC = "
          f"{500 * int(inv['nodes_per_strand']) + 1} nodes at seg {inv['seg_um_realised']} um")


if __name__ == "__main__":
    _demo()
