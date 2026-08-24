r"""Ventral stress fibres: FA-to-FA actomyosin bundles on the basal plane, laid down by a CUDA kernel.

WHAT STANDS HERE AND WHAT DOES NOT.  A ventral stress fibre is a bundle spanning two focal adhesions,
with SARCOMERIC organisation along it: α-actinin Z-bodies at a periodic spacing, NMII bands
anti-registered between them, and formin-nucleated filaments whose barbed ends point at the two FA ends
so the bundle CONTRACTS (``architecture_spec.STRESS_FIBER``: *"this rectification is REQUIRED, random
mixed polarity does NOT contract"*).  Every one of those is a BOND or a law — PHASE 3 and PHASE 2.
What PHASE 1 builds is the bundle geometry: ``n_fibres`` parallel bundles across the basal footprint,
each of ``n_filaments`` hex-packed chains along the fibre axis.

⚠ **AND THE ENGINE HAS ALREADY BEEN BITTEN BY BUILDING THE STATIONS BEFORE THE POPULATION.**
``STATE.md`` (e) item 5 records that SF geometry REFUSES a motor station — ``N_stations: 0`` — so
``nmii_sf_motor`` is unbindable *for a model reason*, with the standing instruction "do not force it".
This module therefore does not place a station at all, and the sarcomere period is carried as a
DECLARED axis in the census rather than materialised as geometry that would decide the question early.

⚠ **THREE PI-GAPS, ALL REFUSED.**  ``architecture_spec``'s own note is explicit: *"sarcomere_um +
N_filaments PI-GATED (no KB point value)"*, with filaments per fibre a band of 7–30 and the shipped
value 20.  And **how many stress fibres a cell has is absent entirely** — the spec describes ONE fibre.
Both counts arrive as :class:`~aleph.world.bond.BondCount`, which cannot be built without a scope and a
provenance; ``sarcomere_um`` is required and recorded, never used to place anything.

THE ONE INVARIANT THIS POPULATION HAS THAT THE OTHERS DO NOT.  Two neighbouring fibres must not
interpenetrate: the bundle's own cross-section is ``hex_cols x spacing_um`` wide and the fibres are
``pitch_um`` apart, so ``pitch_um`` must exceed the bundle width.  It is asserted rather than assumed,
because the two numbers come from different sources — the packing is ~12 nm and the pitch is a cell-
scale figure, and nothing else in the build would notice them crossing.

TOPOLOGY IS ARITHMETIC.  One claim per population per kind; filament ``(f, j)`` is strand
``f*n_filaments + j``.  See :mod:`aleph.world.build.microtubule` for the argument.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``origin``, ``length_um``, ``pitch_um``, ``spacing_um``, ``seg_um``,
    ``sarcomere_um`` [µm]; the counts are dimensionless.  No stiffness, no motor, no force.
  * boundary — a filament of fewer than 3 beads is refused; either count resolving to zero is refused;
    a ``pitch_um`` narrower than the bundle cross-section is refused as interpenetration; an ``axis``
    parallel to the substrate ``normal`` is refused, since a VENTRAL fibre lies in the basal plane.
  * conservation/invariant — the kernel writes only ``lo + t`` for ``t`` in ``[0, n_nodes)`` and the
    launch dimension IS ``n_nodes``.  ``assert_partitioned`` is asserted after the build.
  * CFL/precision — no integration; float64 throughout.  The realised segment length is reported, and
    it is the one a bending coefficient ``α = κ/seg³`` must read.
  * sign sense — the fibres are laid symmetrically about ``origin`` in the lateral direction: fibre
    ``f`` sits at ``(f - (N-1)/2)·pitch``, so the footprint is centred rather than growing off to one
    side.  Polarity is NOT set here: the graded antiparallel rectification is a PHASE 2 property of the
    law that reads the strand, and writing a polarity now would fix it before its law exists.
  * measurement protocol — one kernel launch of ``dim=n_nodes``; nothing is read back.

engine units: length µm.  Runtime: NVIDIA Warp on CUDA; the build refuses a device-less arena.
"""

from __future__ import annotations

import math

import warp as wp

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount
from aleph.world.strand import _require_length

__all__ = ["plan_stress_fibers", "build_stress_fibers", "lay_ventral_bundles_kernel"]

#: ``√3/2`` — the row pitch of a hexagonal lattice at unit spacing.
HEX_ROW_PITCH = wp.constant(wp.float64(0.8660254037844386))


@wp.kernel
def lay_ventral_bundles_kernel(
    lo: wp.int32, n_per: wp.int32, n_fil: wp.int32, n_fiber: wp.int32, hex_cols: wp.int32,
    hex_rows: wp.int32, strand_base: wp.int32, range_ordinal: wp.int32,
    ox: wp.float64, oy: wp.float64, oz: wp.float64,
    ax: wp.float64, ay: wp.float64, az: wp.float64,
    lx: wp.float64, ly: wp.float64, lz: wp.float64,
    nx: wp.float64, ny: wp.float64, nz: wp.float64,
    pitch_um: wp.float64, spacing_um: wp.float64, seg_um: wp.float64,
    position: wp.array(dtype=wp.vec3d), strand_id: wp.array(dtype=wp.int32),
    range_id: wp.array(dtype=wp.int32),
):
    """Write one bundle bead per thread: fibre offset, hex offset and axial position, all from ``tid``.

    The frame ``(axis, lateral, normal)`` is orthonormal and the SAME for every fibre — ventral fibres
    are parallel by construction — so unlike the filopodial kernel this one builds no per-strand frame.
    """
    t = wp.tid()
    per_fiber = n_fil * n_per
    f = t / per_fiber
    rem = t - f * per_fiber
    j = rem / n_per
    i = rem - j * n_per

    axis = wp.vec3d(ax, ay, az)
    lat = wp.vec3d(lx, ly, lz)
    nrm = wp.vec3d(nx, ny, nz)

    row = j / hex_cols
    col = j - row * hex_cols
    stagger = wp.float64(row - wp.int32(2) * (row / wp.int32(2))) * wp.float64(0.5)
    hex_l = (wp.float64(col) - wp.float64(0.5) * wp.float64(hex_cols - 1) + stagger) * spacing_um
    hex_n = (wp.float64(row) - wp.float64(0.5) * wp.float64(hex_rows - 1)) * spacing_um * HEX_ROW_PITCH
    fiber_l = (wp.float64(f) - wp.float64(0.5) * wp.float64(n_fiber - 1)) * pitch_um

    g = lo + t
    position[g] = (wp.vec3d(ox, oy, oz) + (fiber_l + hex_l) * lat + hex_n * nrm
                   + wp.float64(i) * seg_um * axis)
    strand_id[g] = strand_base + f * n_fil + j
    range_id[g] = range_ordinal


def plan_stress_fibers(
    *, count: BondCount, bundle_count: BondCount, support: float = 1.0,
    length_um: float | None = None, seg_um: float | None = None, pitch_um: float | None = None,
    spacing_um: float | None = None, sarcomere_um: float | None = None,
) -> dict[str, object]:
    """Resolve the ventral-fibre inventory from its two declared counts and geometry, no device touched.

    Args:
        count: how many ventral stress fibres this cell has.  **REQUIRED**, and ⚠ no value for it
            exists — ``architecture_spec.STRESS_FIBER`` describes a SINGLE fibre.
        bundle_count: filaments per fibre.  **REQUIRED**.  Band 7–30, shipped 20, and the spec's own
            note says *"N_filaments PI-GATED (no KB point value)"*.
        support: the MEASURED support ``count`` resolves against — 1.0 for an ``explicit`` count, the
            basal contact area [µm²] for an ``areal`` one.  Measured, never analytic.
        length_um: FA-to-FA fibre length [µm].  **REQUIRED** — no default.
        seg_um: bead spacing [µm].  **REQUIRED** — no default.
        pitch_um: centre-to-centre spacing between neighbouring fibres [µm].  **REQUIRED**.
        spacing_um: centre-to-centre spacing between filaments inside one bundle [µm].  **REQUIRED** —
            the ~12 nm actin packing is itself flagged PI-gated in the spec.
        sarcomere_um: the Z-body / NMII band period [µm].  **REQUIRED**, and deliberately RECORDED
            rather than used: it is PI-gated, and placing stations from it now would answer the
            ``N_stations: 0`` question with a build detail instead of a decision.

    Returns:
        The inventory, with both provenance rows and the declared sarcomere period.

    Raises:
        ValueError: on a missing or non-positive length, a filament of fewer than 3 beads, either count
            resolving to zero, or a pitch narrower than the bundle's own cross-section.
    """
    length = _require_length("length_um", length_um)
    seg = _require_length("seg_um", seg_um)
    pitch = _require_length("pitch_um", pitch_um)
    spacing = _require_length("spacing_um", spacing_um)
    sarcomere = _require_length("sarcomere_um", sarcomere_um)

    n_seg = int(round(length / seg))
    n_per = n_seg + 1
    if n_per < 3:
        raise ValueError(
            f"length_um={length} at seg_um={seg} gives {n_per} beads per filament; a filament needs at "
            "least 3 so it carries one bending triple. A filament with no bending stiffness is a "
            "different physical object, not a coarser one."
        )
    n_fiber = count.resolve(support)
    n_fil = bundle_count.resolve(1.0)
    for label, value in (("stress fibres", n_fiber), ("filaments per fibre", n_fil)):
        if value < 1:
            raise ValueError(
                f"the count resolves to {value} {label}. A population that is 'present but has nothing' "
                "is indistinguishable downstream from one that failed to build."
            )
    hex_cols = int(math.ceil(math.sqrt(n_fil)))
    hex_rows = int(math.ceil(n_fil / hex_cols))
    bundle_width = hex_cols * spacing
    if n_fiber > 1 and pitch <= bundle_width:
        raise ValueError(
            f"pitch_um={pitch} is not wider than the bundle's own cross-section {bundle_width} µm "
            f"({hex_cols} columns at spacing_um={spacing}), so neighbouring fibres interpenetrate. The "
            "packing and the pitch come from different sources — nothing else in this build would "
            "notice them crossing, which is why it is asserted here."
        )
    return {
        "kind": "stress_fiber",
        "n_strands": n_fiber * n_fil,
        "n_fibers": n_fiber,
        "filaments_per_fiber": n_fil,
        "nodes_per_strand": n_per,
        "hex_cols": hex_cols, "hex_rows": hex_rows,
        "bundle_width_um": bundle_width,
        "pitch_um": pitch, "spacing_um": spacing,
        "seg_um_requested": seg,
        "seg_um_realised": length / n_seg,
        "length_um_realised": length,
        # Declared, recorded, and NOT used to place anything — see the module docstring and STATE.md (e) 5.
        "sarcomere_um_declared": sarcomere,
        "n_stations_placed": 0,
        "provenance": {"n_stress_fibers": _row(count), "filaments_per_fiber": _row(bundle_count)},
    }


def _row(count: BondCount) -> dict[str, object]:
    """The count with its scope and authority attached, in the shape ``BondFamily.provenance_row`` uses."""
    return {
        "basis": count.basis, "value": count.value, "support_unit": count.support_unit,
        "scope": count.scope, "source_class": str(count.source_class), "provenance": count.provenance,
    }


def build_stress_fibers(
    arena: WorldArena,
    *,
    count: BondCount,
    bundle_count: BondCount,
    support: float = 1.0,
    origin: tuple[float, float, float],
    axis: tuple[float, float, float],
    normal: tuple[float, float, float],
    length_um: float | None = None,
    seg_um: float | None = None,
    pitch_um: float | None = None,
    spacing_um: float | None = None,
    sarcomere_um: float | None = None,
    population: str = "stress_fiber",
) -> dict[str, object]:
    """Claim the ventral-fibre ranges and lay every bundle down with one CUDA kernel.

    Args:
        arena: the world to claim from.  Must hold a CUDA allocation.
        origin: the basal-plane point the centre fibre starts at [µm] — its first FA end.
        axis: the fibre direction, in the basal plane; normalised here.
        normal: the substrate normal.  The lateral direction is ``normal x axis``.
        count / bundle_count / support / length_um / seg_um / pitch_um / spacing_um / sarcomere_um: see
            :func:`plan_stress_fibers`.
        population: the name the four claims are attributed to.

    Returns:
        The census of what was built.

    Raises:
        RuntimeError: if ``arena`` holds no device allocation.
        ValueError: on any refusal from :func:`plan_stress_fibers`, or an ``axis`` parallel to
            ``normal`` — a ventral fibre lies IN the basal plane.
    """
    inventory = plan_stress_fibers(
        count=count, bundle_count=bundle_count, support=support, length_um=length_um, seg_um=seg_um,
        pitch_um=pitch_um, spacing_um=spacing_um, sarcomere_um=sarcomere_um)

    a = _unit("axis", axis)
    n = _unit("normal", normal)
    lat = (n[1] * a[2] - n[2] * a[1], n[2] * a[0] - n[0] * a[2], n[0] * a[1] - n[1] * a[0])
    lat_norm = math.sqrt(sum(v * v for v in lat))
    if lat_norm <= 1e-9:
        raise ValueError(
            f"axis {axis!r} is parallel to normal {normal!r}. A VENTRAL stress fibre lies in the basal "
            "plane between two focal adhesions; one along the substrate normal is a different structure "
            "with a different name, not a rotated ventral fibre."
        )
    lat = tuple(v / lat_norm for v in lat)
    nrm = (a[1] * lat[2] - a[2] * lat[1], a[2] * lat[0] - a[0] * lat[2], a[0] * lat[1] - a[1] * lat[0])

    n_strand = int(inventory["n_strands"])
    n_per = int(inventory["nodes_per_strand"])
    n_nodes = n_strand * n_per
    if not arena.node_arrays:
        raise RuntimeError(
            f"{population}: this arena holds no device allocation, and construction here is a CUDA "
            "kernel rather than a host build that is uploaded. Warp CUDA is the only runtime; build the "
            f"arena with a CUDA device. The inventory it would have built: {inventory['n_fibers']} "
            f"fibres x {inventory['filaments_per_fiber']} filaments x {n_per} beads = {n_nodes} nodes."
        )

    nodes = arena.claim(population, Kind.NODE, n_nodes)
    strands = arena.claim(population, Kind.STRAND, n_strand)
    segments = arena.claim(population, Kind.SEGMENT, n_strand * (n_per - 1))
    angles = arena.claim(population, Kind.ANGLE3, n_strand * (n_per - 2))
    ordinal = len(arena.claims()) - 4

    wp.launch(
        lay_ventral_bundles_kernel, dim=n_nodes,
        inputs=[nodes.lo, n_per, int(inventory["filaments_per_fiber"]), int(inventory["n_fibers"]),
                int(inventory["hex_cols"]), int(inventory["hex_rows"]), strands.lo, ordinal,
                float(origin[0]), float(origin[1]), float(origin[2]), a[0], a[1], a[2],
                lat[0], lat[1], lat[2], nrm[0], nrm[1], nrm[2],
                float(inventory["pitch_um"]), float(inventory["spacing_um"]),
                float(inventory["seg_um_realised"])],
        outputs=[arena.node_arrays["position"], arena.node_arrays["strand_id"],
                 arena.node_arrays["range_id"]],
        device=arena.device,
    )
    arena.assert_partitioned()

    return inventory | {
        "population": population, "n_nodes": nodes.count, "n_segments": segments.count,
        "n_angles": angles.count, "axis": a, "lateral": lat, "normal": nrm,
        "claims": {"node": (nodes.lo, nodes.hi), "strand": (strands.lo, strands.hi),
                   "segment": (segments.lo, segments.hi), "angle3": (angles.lo, angles.hi)},
        "support": float(support),
    }


def _unit(name: str, v: tuple[float, float, float]) -> tuple[float, float, float]:
    """Normalise a direction on the host, or refuse — a NaN direction claims cleanly and renders as dust."""
    x = [float(c) for c in v]
    norm = math.sqrt(sum(c * c for c in x))
    if not math.isfinite(norm) or norm <= 0.0:
        raise ValueError(f"{name} is degenerate; got {v!r}")
    return tuple(c / norm for c in x)


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
    """Self-check: the inventory arithmetic, the interpenetration guard, and every refusal."""
    per_cell = BondCount(
        basis="explicit", value=50.0,
        scope="MCF7 epithelial on glass — NO VALUE EXISTS for ventral stress fibres per cell",
        source_class="PI_GAP",
        provenance="laws/architecture_spec.STRESS_FIBER describes ONE fibre; the per-cell number is "
                   "absent from the contract graph")
    per_fiber = BondCount(
        basis="explicit", value=20.0, scope="ventral SF bundle, generic — band 7-30",
        source_class="PI_GAP",
        provenance="laws/architecture_spec.STRESS_FIBER notes: 'sarcomere_um + N_filaments PI-GATED "
                   "(no KB point value)'")

    base = {"length_um": 10.0, "seg_um": 0.05, "pitch_um": 1.0, "spacing_um": 0.012,
            "sarcomere_um": 1.0}
    inv = plan_stress_fibers(count=per_cell, bundle_count=per_fiber, **base)
    assert inv["nodes_per_strand"] == 201, inv["nodes_per_strand"]
    assert inv["n_strands"] == 1_000 and inv["n_fibers"] == 50
    assert inv["n_stations_placed"] == 0, "no motor station is placed; STATE.md (e) 5 says do not force it"
    assert abs(float(inv["bundle_width_um"]) - 0.06) < 1e-12, inv["bundle_width_um"]

    # The interpenetration guard: at a 12 nm packing the bundle is 60 nm wide, so a pitch of 50 nm is
    # two fibres occupying the same filaments. Neither number alone looks wrong.
    try:
        plan_stress_fibers(count=per_cell, bundle_count=per_fiber, **(base | {"pitch_um": 0.05}))
    except ValueError as exc:
        assert "interpenetrate" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a pitch inside the bundle width must refuse")

    for missing in base:
        kwargs = {k: v for k, v in base.items() if k != missing}
        try:
            plan_stress_fibers(count=per_cell, bundle_count=per_fiber, **kwargs)
        except ValueError as exc:
            assert "no default" in str(exc) and missing in str(exc), str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"a missing {missing} must refuse")

    try:
        plan_stress_fibers(count=per_cell, bundle_count=per_fiber, **(base | {"seg_um": 10.0}))
    except ValueError as exc:
        assert "at least 3" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a 2-bead filament must refuse")

    arena = WorldArena(capacity={Kind.NODE: 10, Kind.STRAND: 10, Kind.SEGMENT: 10, Kind.ANGLE3: 10})
    try:
        build_stress_fibers(arena, count=per_cell, bundle_count=per_fiber, origin=(0.0, 0.0, 0.0),
                            axis=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0), **base)
    except RuntimeError as exc:
        assert "no device allocation" in str(exc) and "201000 nodes" in str(exc), str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device-less build must refuse")

    try:
        build_stress_fibers(arena, count=per_cell, bundle_count=per_fiber, origin=(0.0, 0.0, 0.0),
                            axis=(0.0, 0.0, 1.0), normal=(0.0, 0.0, 1.0), **base)
    except ValueError as exc:
        assert "different structure" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a fibre along the substrate normal must refuse")
    assert arena.n_live(Kind.NODE) == 0, "a refused build claims nothing"

    chars = _codegen_check()
    assert chars > 0

    print(f"stress_fiber self-check OK — CUDA codegen {chars} chars; {inv['n_fibers']} fibres x "
          f"{inv['filaments_per_fiber']} filaments x {inv['nodes_per_strand']} beads = "
          f"{inv['n_strands'] * int(inv['nodes_per_strand'])} nodes; bundle "
          f"{inv['bundle_width_um'] * 1000:.0f} nm wide at pitch {inv['pitch_um']} µm")


if __name__ == "__main__":
    _demo()
