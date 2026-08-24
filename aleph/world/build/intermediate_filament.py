r"""The perinuclear intermediate-filament cage: radial spokes, laid down by a CUDA kernel.

WHAT THIS COMPARTMENT IS FOR.  The IF cage is the nucleus↔cortex load path.  Without it the nucleus is
mechanically DECOUPLED from the cell surface — the defect ``INTERNAL_DISPLACEMENT_DIAGNOSIS_2026-07-16``
recorded and the reason ``components/solid/intermediate_filament.py`` exists.  Geometrically it is the
simplest population this session owns: ``n_spokes`` chains running radially across the cytoplasmic gap,
on Fibonacci directions, from just outside the nuclear envelope to just inside the cortex.

WHAT IS DIFFERENT FROM THE BUILDER THIS REPLACES.  The incumbent builds the cage on the host and then
runs a ``cKDTree`` to find, for each spoke end, the nearest nucleus bead and the nearest cortex node.
Those nearest-neighbour searches are the LINC and cortex-anchor BONDS — PHASE 3 — and they need the two
other populations to exist first.  This module builds the spokes and nothing else: no LINC, no anchor,
no plectin crosslink, no stiffness.  Splitting it that way is not tidiness.  ``k_linc`` carries the note
*"8 pN is a TENSION not a stiffness"* in ``compartments.py``, so the bond that search feeds is blocked
on a PI decision, and building the geometry behind that decision would have coupled the two.

⚠ **THREE PI-GAPS, AND ALL THREE ARE REFUSED HERE RATHER THAN DEFAULTED.**

  * ``n_spokes`` — the incumbent ships ``n_fil = 60`` and its own docstring calls the native cage count
    a "provisional PI GAP".  It arrives as a :class:`~aleph.world.bond.BondCount`, which cannot be
    built without a scope and a provenance.
  * ``gap_um`` — the incumbent hardcodes ``r0, r1 = R_nuc + 0.15, R_cortex - 0.15``.  0.15 µm is a
    CONVENIENCE with no source, and it sets how much of the cytoplasm the spoke actually spans, so it
    is declared here rather than carried.
  * ``seg_um`` — no default, for the reason ``build_strand`` gives.

WHAT THE KERNEL DOES AND DOES NOT KNOW.  It computes a spoke direction from the thread's own index by
the same Fibonacci spiral ``laws/microtubule._fibonacci_directions`` uses, so spoke ``k`` points where
the incumbent's spoke ``k`` points and PHASE 2's node-by-node A/B is meaningful.  It does not know the
nucleus mesh, the cortex mesh, or which bead is nearest to what.

TOPOLOGY IS ARITHMETIC.  One claim per population per kind; segment ``s`` of spoke ``k`` joins
``lo + k*n_per + i`` and ``+1``.  See :mod:`aleph.world.build.microtubule` for the full argument.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``centre``, ``R_nuc_um``, ``R_cortex_um``, ``gap_um``, ``seg_um`` [µm]; the spoke
    count is dimensionless.  No stiffness and no force appear in this module at all.
  * boundary — a span that is non-positive (the gap eats the cytoplasm, or the nucleus is larger than
    the cortex) is refused with both radii named; a spoke of fewer than 3 beads is refused, the same
    rule as ``build_strand``; a count resolving to zero spokes is refused.
  * conservation/invariant — the kernel writes only ``lo + t`` for ``t`` in ``[0, n_nodes)`` and the
    launch dimension IS ``n_nodes``, so leaving the claim is unreachable.  ``assert_partitioned`` after.
  * CFL/precision — no integration; float64 throughout.  The realised segment length is the span
    divided by the whole number of segments and is the value reported, never the requested one.
  * sign sense — the spoke runs OUTWARD: bead 0 is the inner (nuclear) end and the last bead is the
    outer (cortical) end, so ``inner = lo + k*n_per`` and ``outer = inner + n_per - 1``, which is what
    PHASE 3's LINC and anchor families will address.
  * measurement protocol — one kernel launch of ``dim=n_nodes``; nothing is read back.

engine units: length µm.  Runtime: NVIDIA Warp on CUDA; the build refuses a device-less arena.
"""

from __future__ import annotations

import warp as wp

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount
from aleph.world.strand import _require_length

__all__ = ["plan_intermediate_filaments", "build_intermediate_filaments", "lay_spokes_kernel"]

#: ``π(1+√5)``, the golden-angle azimuth — identical to the aster's, so spoke ``k`` and arm ``k`` share
#: a direction table and the two populations interleave rather than stacking on one axis.
GOLDEN_AZIMUTH = wp.constant(wp.float64(10.166407384630519))


@wp.kernel
def lay_spokes_kernel(
    lo: wp.int32, n_per: wp.int32, n_strand: wp.int32, strand_base: wp.int32, range_ordinal: wp.int32,
    cx: wp.float64, cy: wp.float64, cz0: wp.float64, r_inner: wp.float64, seg_um: wp.float64,
    position: wp.array(dtype=wp.vec3d), strand_id: wp.array(dtype=wp.int32),
    range_id: wp.array(dtype=wp.int32),
):
    """Write one spoke bead per thread at radius ``r_inner + i*seg_um`` along Fibonacci direction ``k``."""
    t = wp.tid()
    k = t / n_per
    i = t - k * n_per

    kd = wp.float64(k) + wp.float64(0.5)
    cosp = wp.float64(1.0) - wp.float64(2.0) * kd / wp.float64(n_strand)
    sinp = wp.sqrt(wp.max(wp.float64(0.0), wp.float64(1.0) - cosp * cosp))
    az = GOLDEN_AZIMUTH * kd

    r = r_inner + wp.float64(i) * seg_um
    g = lo + t
    position[g] = wp.vec3d(cx + r * sinp * wp.cos(az), cy + r * sinp * wp.sin(az), cz0 + r * cosp)
    strand_id[g] = strand_base + k
    range_id[g] = range_ordinal


def plan_intermediate_filaments(
    *, count: BondCount, R_nuc_um: float | None = None, R_cortex_um: float | None = None,
    gap_um: float | None = None, seg_um: float | None = None,
) -> dict[str, object]:
    """Resolve the cage's inventory from its declared count and geometry, touching no device.

    Args:
        count: how many spokes, for which cell, on whose authority.  **REQUIRED** — the incumbent's
            ``n_fil = 60`` is labelled a provisional PI GAP by its own docstring, so this arrives with
            a source class that says which it is.
        R_nuc_um: nuclear radius [µm].  **REQUIRED** — no default.
        R_cortex_um: cortical radius [µm].  **REQUIRED** — no default.
        gap_um: clearance held at BOTH ends, so the spoke spans ``R_cortex - R_nuc - 2*gap`` [µm].
            **REQUIRED** — the incumbent's 0.15 is an unsourced convenience and it decides how much of
            the cytoplasm the load path actually crosses.
        seg_um: bead spacing along the spoke [µm].  **REQUIRED** — no default.

    Returns:
        ``n_strands`` is absent — it needs the support — but every per-spoke figure is here:
        ``nodes_per_strand``, ``span_um``, the realised segment length, and the provenance row.

    Raises:
        ValueError: on a missing or non-positive length, a non-positive span, or a spoke of fewer than
            3 beads.
    """
    r_nuc = _require_length("R_nuc_um", R_nuc_um)
    r_cor = _require_length("R_cortex_um", R_cortex_um)
    gap = _require_length("gap_um", gap_um)
    seg = _require_length("seg_um", seg_um)

    span = r_cor - r_nuc - 2.0 * gap
    if span <= 0.0:
        raise ValueError(
            f"the cytoplasmic span is {span} µm: R_cortex_um={r_cor} minus R_nuc_um={r_nuc} minus twice "
            f"gap_um={gap} leaves nothing for the spoke to cross. The cage is the nucleus-to-cortex load "
            "path; with no span there is no path, which is the decoupled nucleus this compartment exists "
            "to fix."
        )
    n_seg = int(round(span / seg))
    n_per = n_seg + 1
    if n_per < 3:
        raise ValueError(
            f"a span of {span} µm at seg_um={seg} gives {n_per} beads per spoke; a spoke needs at least "
            "3 so it carries one bending triple. A filament with no bending stiffness is a different "
            "physical object, not a coarser one."
        )
    return {
        "kind": "intermediate_filament",
        "count_basis": count.basis,
        "nodes_per_strand": n_per,
        "span_um": span,
        "r_inner_um": r_nuc + gap,
        "seg_um_requested": seg,
        "seg_um_realised": span / n_seg,
        "provenance": {"n_spokes": {
            "basis": count.basis, "value": count.value, "support_unit": count.support_unit,
            "scope": count.scope, "source_class": str(count.source_class),
            "provenance": count.provenance,
        }},
    }


def build_intermediate_filaments(
    arena: WorldArena,
    *,
    count: BondCount,
    support: float = 1.0,
    centre: tuple[float, float, float],
    R_nuc_um: float | None = None,
    R_cortex_um: float | None = None,
    gap_um: float | None = None,
    seg_um: float | None = None,
    population: str = "intermediate_filament",
) -> dict[str, object]:
    """Claim the cage's ranges and lay its spokes down with a CUDA kernel.

    Args:
        arena: the world to claim from.  Must hold a CUDA allocation.
        count: see :func:`plan_intermediate_filaments`.
        support: the MEASURED support the count resolves against — 1.0 for an ``explicit`` count, the
            nuclear surface area [µm²] for an ``areal`` one.  Measured, never analytic: the area of the
            mesh that exists, not of the sphere it approximates.
        centre: the nuclear centre [µm]; the spokes are radial about it.
        R_nuc_um / R_cortex_um / gap_um / seg_um: see :func:`plan_intermediate_filaments`.
        population: the name the four claims are attributed to.

    Returns:
        The census of what was built.

    Raises:
        RuntimeError: if ``arena`` holds no device allocation.
        ValueError: on any refusal from :func:`plan_intermediate_filaments`, or a count resolving to
            zero spokes.
    """
    inventory = plan_intermediate_filaments(
        count=count, R_nuc_um=R_nuc_um, R_cortex_um=R_cortex_um, gap_um=gap_um, seg_um=seg_um)
    n_strand = count.resolve(support)
    if n_strand < 1:
        raise ValueError(
            f"{population}: the count resolves to {n_strand} spokes against support={support}. A cage "
            "with no spokes is not a weak load path, it is no load path — and it would be recorded as "
            "'present' either way, which is the fact an empty claim destroys."
        )
    n_per = int(inventory["nodes_per_strand"])
    if not arena.node_arrays:
        raise RuntimeError(
            f"{population}: this arena holds no device allocation, and construction here is a CUDA "
            "kernel rather than a host build that is uploaded. Warp CUDA is the only runtime; build the "
            f"arena with a CUDA device. The inventory it would have built: {n_strand} spokes x {n_per} "
            f"beads = {n_strand * n_per} nodes."
        )

    n_nodes = n_strand * n_per
    nodes = arena.claim(population, Kind.NODE, n_nodes)
    strands = arena.claim(population, Kind.STRAND, n_strand)
    segments = arena.claim(population, Kind.SEGMENT, n_strand * (n_per - 1))
    angles = arena.claim(population, Kind.ANGLE3, n_strand * (n_per - 2))
    ordinal = len(arena.claims()) - 4

    wp.launch(
        lay_spokes_kernel, dim=n_nodes,
        inputs=[nodes.lo, n_per, n_strand, strands.lo, ordinal,
                float(centre[0]), float(centre[1]), float(centre[2]),
                float(inventory["r_inner_um"]), float(inventory["seg_um_realised"])],
        outputs=[arena.node_arrays["position"], arena.node_arrays["strand_id"],
                 arena.node_arrays["range_id"]],
        device=arena.device,
    )
    arena.assert_partitioned()

    return inventory | {
        "population": population, "n_strands": n_strand, "n_nodes": nodes.count,
        "n_segments": segments.count, "n_angles": angles.count,
        # The two ends PHASE 3 will address, as arithmetic rather than as a stored index list.
        "inner_end_stride": n_per, "inner_end_lo": nodes.lo, "outer_end_offset": n_per - 1,
        "claims": {"node": (nodes.lo, nodes.hi), "strand": (strands.lo, strands.hi),
                   "segment": (segments.lo, segments.hi), "angle3": (angles.lo, angles.hi)},
        "support": float(support),
    }


def _codegen_check() -> int:
    """Type-check this module's kernels by emitting their CUDA source, with no device and no execution.

    The dev machine has no card, so the alternative is shipping kernel syntax that nobody has put in
    front of a compiler; this generates the CUDA C++ Warp would compile, which type-checks every
    expression, and produces NO number — so it is not a CPU result under any reading.
    """
    from warp._src.context import ModuleBuilder  # private, deliberately: there is no public codegen API

    module = wp.get_module(__name__)
    return len(ModuleBuilder(module, module.options).codegen("cuda"))


def _demo() -> None:
    """Self-check: the inventory arithmetic and every refusal, none of which needs a device."""
    gap_count = BondCount(
        basis="explicit", value=60.0,
        scope="MCF7 epithelial, interphase — keratin K8/K18 cage; the native cage count is not sourced",
        source_class="PI_GAP",
        provenance="components/solid/intermediate_filament.py docstring: 'GAP magnitudes (native cage "
                   "count n_fil, ...) are provisional PI GAPs'")

    inv = plan_intermediate_filaments(count=gap_count, R_nuc_um=5.0, R_cortex_um=7.5, gap_um=0.15,
                                      seg_um=0.05)
    assert abs(float(inv["span_um"]) - 2.2) < 1e-12, inv["span_um"]
    assert inv["nodes_per_strand"] == 45, inv["nodes_per_strand"]
    assert abs(float(inv["r_inner_um"]) - 5.15) < 1e-12

    # A gap that eats the cytoplasm names both radii rather than building a spoke of length zero.
    try:
        plan_intermediate_filaments(count=gap_count, R_nuc_um=5.0, R_cortex_um=7.5, gap_um=2.0,
                                    seg_um=0.05)
    except ValueError as exc:
        assert "no path" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a non-positive span must refuse")

    # Every length refuses to default, in the one spelling the codebase has for it.
    base = {"R_nuc_um": 5.0, "R_cortex_um": 7.5, "gap_um": 0.15, "seg_um": 0.05}
    for missing in base:
        kwargs = {k: v for k, v in base.items() if k != missing}
        try:
            plan_intermediate_filaments(count=gap_count, **kwargs)
        except ValueError as exc:
            assert "no default" in str(exc) and missing in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"a missing {missing} must refuse")

    # Too coarse to carry a bending triple is a different object, not a coarser one.
    try:
        plan_intermediate_filaments(count=gap_count, R_nuc_um=5.0, R_cortex_um=7.5, gap_um=0.15,
                                    seg_um=2.0)
    except ValueError as exc:
        assert "at least 3" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a 2-bead spoke must refuse")

    arena = WorldArena(capacity={Kind.NODE: 10, Kind.STRAND: 10, Kind.SEGMENT: 10, Kind.ANGLE3: 10})
    try:
        build_intermediate_filaments(arena, count=gap_count, centre=(0.0, 0.0, 0.0), **base)
    except RuntimeError as exc:
        assert "no device allocation" in str(exc) and "2700 nodes" in str(exc), str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device-less build must refuse")
    assert arena.n_live(Kind.NODE) == 0, "a refused build claims nothing"

    chars = _codegen_check()
    assert chars > 0

    print(f"intermediate_filament self-check OK — CUDA codegen {chars} chars; 60 spokes x "
          f"{inv['nodes_per_strand']} beads = {60 * int(inv['nodes_per_strand'])} nodes across a "
          f"{inv['span_um']:.2f} µm span")


if __name__ == "__main__":
    _demo()
