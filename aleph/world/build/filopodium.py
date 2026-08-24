r"""Filopodia: parallel formin bundles rooted on the cell surface, laid down by a CUDA kernel.

A FILOPODIUM IS A BUNDLE, AND THAT IS THE WHOLE STRUCTURAL CLAIM.  ``laws/architecture_spec.FILOPODIUM``
states it as one row of the unified weaving table — a tight PARALLEL formin bundle, 10–30 filaments,
uniform polarity (all barbed ends at the tip), no motor in the core, crosslinked at ~7–8 nm spacing.
Geometry here is exactly that and nothing more: ``n_filaments`` straight chains per finger, hex-packed
in the plane perpendicular to the finger, rooted at the cell surface and running outward.

THE HEX PACKING IS COMPUTED FROM THE THREAD INDEX, NOT FROM A TABLE.  Row ``j/m``, column ``j - row*m``,
odd rows offset by half a spacing, row pitch ``spacing·√3/2``.  So the bundle's cross-section is a real
hexagonal lattice and it costs no host array — which matters because the alternative pattern (build the
lattice in NumPy, upload it) is the host mirror this port exists to remove.

⚠ **TWO COUNTS, AND THE PER-CELL ONE HAS NO SOURCE AT ALL.**  ``architecture_spec`` gives filaments per
bundle (20, band 10–30, and the bundler is a stand-in: *"fascin kinetics are PI-gated"*).  It gives
nothing whatever for **how many filopodia a cell has** — that number is absent from the KB, not merely
unratified.  Both arrive as :class:`~aleph.world.bond.BondCount`, which cannot be constructed without a
scope and a provenance, so the artifact records which of the two is which.

⚠ **AND THE PLACEMENT IS A PHYSIOLOGICAL STATEMENT, SO IT IS DECLARED.**  Filopodia on a polarised
migrating cell sit at the leading edge; on a rounded cell they are spread.  ``cap_deg`` is the half-angle
of the spherical cap the fingers are distributed over, about ``axis`` — 180° is the whole sphere and 30°
is a leading edge.  It has no default, because "spread evenly over the sphere" is not a neutral choice,
it is the unpolarised cell, and picking it silently would be picking a cell state.

WHAT IS NOT HERE.  The fascin crosslinks between neighbouring filaments in a bundle, the nascent
adhesion clutch at the tip, and the membrane the finger pushes into are all BONDS — PHASE 3 — and
``filopodium_nascent_clutch`` is one of the eight EXTERNAL connectors this plan does not touch at all.

TOPOLOGY IS ARITHMETIC.  One claim per population per kind; filament ``(f, j)`` is strand
``f*n_filaments + j`` and its bead ``i`` is node ``lo + (f*n_filaments + j)*n_per + i``.  See
:mod:`aleph.world.build.microtubule` for the argument.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``centre``, ``root_R_um``, ``contour_um``, ``seg_um``, ``spacing_um`` [µm];
    ``cap_deg`` is an angle [deg] and is converted once, on the host; the counts are dimensionless.
  * boundary — a filament of fewer than 3 beads is refused, the same rule as ``build_strand``; a
    ``cap_deg`` outside (0, 180] is refused; either count resolving to zero is refused; a degenerate
    ``axis`` is refused before it reaches the kernel, where a normalise would return NaN silently.
  * conservation/invariant — the kernel writes only ``lo + t`` for ``t`` in ``[0, n_nodes)`` and the
    launch dimension IS ``n_nodes``.  The lattice offset is bounded by ``spacing·m``, so two fingers
    cannot interpenetrate unless the caller's own spacing says they do.  ``assert_partitioned`` after.
  * CFL/precision — no integration; float64 throughout.  The realised segment length is reported.
  * sign sense — polarity is UNIFORM and outward: bead 0 is the root, the last bead is the tip, so
    ``Strand.polarity = +1`` and the barbed end is the last node.  This is not a convention choice —
    a filopodium with mixed polarity does not elongate.
  * measurement protocol — one kernel launch of ``dim=n_nodes``; nothing is read back.

engine units: length µm, angle deg on the interface and rad inside.  Runtime: NVIDIA Warp on CUDA; the
build refuses a device-less arena.
"""

from __future__ import annotations

import math

import warp as wp

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount
from aleph.world.strand import _require_length

__all__ = ["plan_filopodia", "build_filopodia", "lay_bundles_kernel", "first_tangent"]

GOLDEN_AZIMUTH = wp.constant(wp.float64(10.166407384630519))
#: ``√3/2`` — the row pitch of a hexagonal lattice at unit spacing.  Written as a constant rather than
#: computed per thread because it is the one number in the packing that is geometry rather than index.
HEX_ROW_PITCH = wp.constant(wp.float64(0.8660254037844386))


@wp.func
def first_tangent(d: wp.vec3d) -> wp.vec3d:
    """A unit vector perpendicular to ``d``, chosen so the cross product never degenerates.

    The seed axis is swapped when ``d`` is nearly parallel to x, which is the whole trick: a fixed seed
    gives a zero-length cross product for one direction in the population, and that filament's entire
    bundle would collapse onto its own axis — a defect that shows up as a rendering artefact long before
    anyone reads it as NaN.
    """
    a = wp.vec3d(wp.float64(1.0), wp.float64(0.0), wp.float64(0.0))
    if wp.abs(d[0]) > wp.float64(0.9):
        a = wp.vec3d(wp.float64(0.0), wp.float64(1.0), wp.float64(0.0))
    return wp.normalize(wp.cross(d, a))


@wp.kernel
def lay_bundles_kernel(
    lo: wp.int32, n_per: wp.int32, n_fil: wp.int32, n_finger: wp.int32, hex_cols: wp.int32,
    hex_rows: wp.int32, strand_base: wp.int32, range_ordinal: wp.int32,
    cx: wp.float64, cy: wp.float64, cz0: wp.float64,
    axx: wp.float64, axy: wp.float64, axz: wp.float64, cos_cap: wp.float64,
    root_R_um: wp.float64, spacing_um: wp.float64, seg_um: wp.float64,
    position: wp.array(dtype=wp.vec3d), strand_id: wp.array(dtype=wp.int32),
    range_id: wp.array(dtype=wp.int32),
):
    """Write one bundle bead per thread: finger direction, hex offset and bead position, all from ``tid``."""
    t = wp.tid()
    per_finger = n_fil * n_per
    f = t / per_finger
    rem = t - f * per_finger
    j = rem / n_per
    i = rem - j * n_per

    # The finger direction: a Fibonacci spiral restricted to a spherical cap about `axis`.
    axis = wp.vec3d(axx, axy, axz)
    b1 = first_tangent(axis)
    b2 = wp.cross(axis, b1)
    fd = wp.float64(f) + wp.float64(0.5)
    cosp = wp.float64(1.0) - (wp.float64(1.0) - cos_cap) * fd / wp.float64(n_finger)
    sinp = wp.sqrt(wp.max(wp.float64(0.0), wp.float64(1.0) - cosp * cosp))
    az = GOLDEN_AZIMUTH * fd
    d = wp.normalize(sinp * (wp.cos(az) * b1 + wp.sin(az) * b2) + cosp * axis)

    # The hex lattice in the plane perpendicular to the finger.
    u1 = first_tangent(d)
    u2 = wp.cross(d, u1)
    row = j / hex_cols
    col = j - row * hex_cols
    stagger = wp.float64(row - wp.int32(2) * (row / wp.int32(2))) * wp.float64(0.5)
    ox = (wp.float64(col) - wp.float64(0.5) * wp.float64(hex_cols - 1) + stagger) * spacing_um
    oy = (wp.float64(row) - wp.float64(0.5) * wp.float64(hex_rows - 1)) * spacing_um * HEX_ROW_PITCH

    root = wp.vec3d(cx, cy, cz0) + root_R_um * d
    p = root + ox * u1 + oy * u2 + wp.float64(i) * seg_um * d

    g = lo + t
    position[g] = p
    strand_id[g] = strand_base + f * n_fil + j
    range_id[g] = range_ordinal


def plan_filopodia(
    *, count: BondCount, bundle_count: BondCount, support: float = 1.0,
    contour_um: float | None = None, seg_um: float | None = None, spacing_um: float | None = None,
    cap_deg: float | None = None, root_R_um: float | None = None,
) -> dict[str, object]:
    """Resolve the filopodial inventory from its two declared counts and its geometry, no device touched.

    Args:
        count: how many filopodia this cell has, for which cell, on whose authority.  **REQUIRED**, and
            ⚠ no value for it exists in the KB — it is absent, not merely unratified.
        bundle_count: filaments per bundle.  **REQUIRED**.  ``architecture_spec`` gives 20 in a band of
            10–30, with the bundler itself a PI-gated stand-in.
        support: the MEASURED support ``count`` resolves against — 1.0 for an ``explicit`` count, the
            cell's surface area [µm²] for an ``areal`` one.  Measured, never analytic.
        contour_um: filament contour length [µm].  **REQUIRED** — no default.
        seg_um: bead spacing [µm].  **REQUIRED** — no default.
        spacing_um: centre-to-centre filament spacing in the bundle [µm].  **REQUIRED** — the ~7–8 nm
            fascin figure is geometric but it is not this engine's to assume.
        cap_deg: half-angle [deg] of the cap the fingers are spread over, about ``axis``.  **REQUIRED**
            — it encodes whether the cell is polarised, which is a cell state, not a build detail.
        root_R_um: radius at which the fingers root on the cell surface [µm].  **REQUIRED**.

    Returns:
        The inventory, including both provenance rows.

    Raises:
        ValueError: on a missing or non-positive length, a ``cap_deg`` outside (0, 180], a filament of
            fewer than 3 beads, or either count resolving to zero.
    """
    contour = _require_length("contour_um", contour_um)
    seg = _require_length("seg_um", seg_um)
    spacing = _require_length("spacing_um", spacing_um)
    root_r = _require_length("root_R_um", root_R_um)
    cap = _require_length("cap_deg", cap_deg)
    if cap > 180.0:
        raise ValueError(
            f"cap_deg={cap} exceeds 180°, which is the whole sphere. A cap wider than the sphere is not "
            "a more spread-out cell, it is an arithmetic error in the caller's polarisation."
        )

    n_seg = int(round(contour / seg))
    n_per = n_seg + 1
    if n_per < 3:
        raise ValueError(
            f"contour_um={contour} at seg_um={seg} gives {n_per} beads per filament; a filament needs at "
            "least 3 so it carries one bending triple. A filament with no bending stiffness is a "
            "different physical object, not a coarser one."
        )
    n_finger = count.resolve(support)
    n_fil = bundle_count.resolve(1.0)
    for label, value in (("filopodia", n_finger), ("filaments per bundle", n_fil)):
        if value < 1:
            raise ValueError(
                f"the count resolves to {value} {label}. A population that is 'present but has nothing' "
                "is indistinguishable downstream from one that failed to build, and the two are "
                "different facts."
            )
    hex_cols = int(math.ceil(math.sqrt(n_fil)))
    return {
        "kind": "filopodium",
        "n_strands": n_finger * n_fil,
        "n_fingers": n_finger,
        "filaments_per_bundle": n_fil,
        "nodes_per_strand": n_per,
        "hex_cols": hex_cols,
        "hex_rows": int(math.ceil(n_fil / hex_cols)),
        "spacing_um": spacing,
        "root_R_um": root_r,
        "cap_deg": cap,
        "cos_cap": math.cos(math.radians(cap)),
        "seg_um_requested": seg,
        "seg_um_realised": contour / n_seg,
        "contour_um_realised": contour,
        "provenance": {
            "n_filopodia": _row(count), "filaments_per_bundle": _row(bundle_count),
        },
    }


def _row(count: BondCount) -> dict[str, object]:
    """The count with its scope and authority attached, in the shape ``BondFamily.provenance_row`` uses."""
    return {
        "basis": count.basis, "value": count.value, "support_unit": count.support_unit,
        "scope": count.scope, "source_class": str(count.source_class), "provenance": count.provenance,
    }


def build_filopodia(
    arena: WorldArena,
    *,
    count: BondCount,
    bundle_count: BondCount,
    support: float = 1.0,
    centre: tuple[float, float, float],
    axis: tuple[float, float, float],
    contour_um: float | None = None,
    seg_um: float | None = None,
    spacing_um: float | None = None,
    cap_deg: float | None = None,
    root_R_um: float | None = None,
    population: str = "filopodium",
) -> dict[str, object]:
    """Claim the filopodial ranges and lay every bundle down with one CUDA kernel.

    Args:
        arena: the world to claim from.  Must hold a CUDA allocation.
        centre: the cell centre the fingers root about [µm].
        axis: the polarisation axis the cap is centred on; normalised here.
        count / bundle_count / support / contour_um / seg_um / spacing_um / cap_deg / root_R_um: see
            :func:`plan_filopodia`.
        population: the name the four claims are attributed to.

    Returns:
        The census of what was built.

    Raises:
        RuntimeError: if ``arena`` holds no device allocation.
        ValueError: on any refusal from :func:`plan_filopodia`, or a degenerate ``axis``.
    """
    inventory = plan_filopodia(
        count=count, bundle_count=bundle_count, support=support, contour_um=contour_um, seg_um=seg_um,
        spacing_um=spacing_um, cap_deg=cap_deg, root_R_um=root_R_um)

    ax = [float(v) for v in axis]
    norm = math.sqrt(sum(v * v for v in ax))
    if not math.isfinite(norm) or norm <= 0.0:
        raise ValueError(
            f"axis is degenerate; got {axis!r}. Normalising it in the kernel would return NaN for every "
            "finger and the population would stand as a cloud of NaN positions that claims fine."
        )
    ax = [v / norm for v in ax]

    n_strand = int(inventory["n_strands"])
    n_per = int(inventory["nodes_per_strand"])
    n_nodes = n_strand * n_per
    if not arena.node_arrays:
        raise RuntimeError(
            f"{population}: this arena holds no device allocation, and construction here is a CUDA "
            "kernel rather than a host build that is uploaded. Warp CUDA is the only runtime; build the "
            f"arena with a CUDA device. The inventory it would have built: {inventory['n_fingers']} "
            f"fingers x {inventory['filaments_per_bundle']} filaments x {n_per} beads = {n_nodes} nodes."
        )

    nodes = arena.claim(population, Kind.NODE, n_nodes)
    strands = arena.claim(population, Kind.STRAND, n_strand)
    segments = arena.claim(population, Kind.SEGMENT, n_strand * (n_per - 1))
    angles = arena.claim(population, Kind.ANGLE3, n_strand * (n_per - 2))
    ordinal = len(arena.claims()) - 4

    wp.launch(
        lay_bundles_kernel, dim=n_nodes,
        inputs=[nodes.lo, n_per, int(inventory["filaments_per_bundle"]), int(inventory["n_fingers"]),
                int(inventory["hex_cols"]), int(inventory["hex_rows"]), strands.lo, ordinal,
                float(centre[0]), float(centre[1]), float(centre[2]), ax[0], ax[1], ax[2],
                float(inventory["cos_cap"]), float(inventory["root_R_um"]),
                float(inventory["spacing_um"]), float(inventory["seg_um_realised"])],
        outputs=[arena.node_arrays["position"], arena.node_arrays["strand_id"],
                 arena.node_arrays["range_id"]],
        device=arena.device,
    )
    arena.assert_partitioned()

    return inventory | {
        "population": population, "n_nodes": nodes.count, "n_segments": segments.count,
        "n_angles": angles.count, "axis": tuple(ax),
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
    """Self-check: the inventory arithmetic, the hex lattice shape, and every refusal."""
    per_cell = BondCount(
        basis="explicit", value=50.0,
        scope="MCF7 epithelial, spread on glass — NO VALUE EXISTS for filopodia per cell",
        source_class="PI_GAP",
        provenance="absent from the contract graph; architecture_spec.FILOPODIUM gives filaments per "
                   "BUNDLE only and says nothing about how many bundles a cell carries")
    per_bundle = BondCount(
        basis="explicit", value=20.0, scope="filopodial core bundle, generic — band 10-30",
        source_class="UNRATIFIED_PROXY",
        provenance="laws/architecture_spec.FILOPODIUM: n_filaments=20, 'Parallel formin bundle (10-30 "
                   "filaments, ~7-8 nm crosslink spacing). fascin = PI-gated.'")

    inv = plan_filopodia(count=per_cell, bundle_count=per_bundle, contour_um=3.0, seg_um=0.05,
                         spacing_um=0.008, cap_deg=60.0, root_R_um=7.5)
    assert inv["nodes_per_strand"] == 61, inv["nodes_per_strand"]
    assert inv["n_strands"] == 1_000 and inv["n_fingers"] == 50
    # 20 filaments tile a 5-wide lattice exactly 4 rows deep — the packing covers the bundle with no
    # empty row, which is what makes the cross-section a bundle rather than a ragged patch.
    assert (inv["hex_cols"], inv["hex_rows"]) == (5, 4), (inv["hex_cols"], inv["hex_rows"])
    assert abs(float(inv["cos_cap"]) - 0.5) < 1e-12, "a 60° cap has cos = 1/2"

    # Every length refuses to default.
    base = {"contour_um": 3.0, "seg_um": 0.05, "spacing_um": 0.008, "cap_deg": 60.0, "root_R_um": 7.5}
    for missing in base:
        kwargs = {k: v for k, v in base.items() if k != missing}
        try:
            plan_filopodia(count=per_cell, bundle_count=per_bundle, **kwargs)
        except ValueError as exc:
            assert "no default" in str(exc) and missing in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"a missing {missing} must refuse")

    # A cap wider than the sphere is an arithmetic error, not a more spread-out cell.
    try:
        plan_filopodia(count=per_cell, bundle_count=per_bundle, **(base | {"cap_deg": 200.0}))
    except ValueError as exc:
        assert "whole sphere" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("cap_deg > 180 must refuse")

    # Too coarse to carry a bending triple is a different object, not a coarser one.
    try:
        plan_filopodia(count=per_cell, bundle_count=per_bundle, **(base | {"seg_um": 3.0}))
    except ValueError as exc:
        assert "at least 3" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a 2-bead filament must refuse")

    arena = WorldArena(capacity={Kind.NODE: 10, Kind.STRAND: 10, Kind.SEGMENT: 10, Kind.ANGLE3: 10})
    try:
        build_filopodia(arena, count=per_cell, bundle_count=per_bundle, centre=(0.0, 0.0, 0.0),
                        axis=(0.0, 0.0, 1.0), **base)
    except RuntimeError as exc:
        assert "no device allocation" in str(exc) and "61000 nodes" in str(exc), str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device-less build must refuse")

    # A degenerate axis is refused on the host, where it is still a number rather than a NaN cloud.
    try:
        build_filopodia(arena, count=per_cell, bundle_count=per_bundle, centre=(0.0, 0.0, 0.0),
                        axis=(0.0, 0.0, 0.0), **base)
    except ValueError as exc:
        assert "degenerate" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a zero axis must refuse")
    assert arena.n_live(Kind.NODE) == 0, "a refused build claims nothing"

    chars = _codegen_check()
    assert chars > 0

    print(f"filopodium self-check OK — CUDA codegen {chars} chars; {inv['n_fingers']} fingers x "
          f"{inv['filaments_per_bundle']} filaments x {inv['nodes_per_strand']} beads = "
          f"{inv['n_strands'] * int(inv['nodes_per_strand'])} nodes")


if __name__ == "__main__":
    _demo()
