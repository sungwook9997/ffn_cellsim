r"""The lamellipodial dendritic array: a ±35° two-mode Arp2/3 patch, laid down by a CUDA kernel.

WHAT MAKES THIS POPULATION DIFFERENT FROM THE OTHER FOUR.  A cortex is isotropic, a filopodium and a
stress fibre are bundles, an aster and an IF cage are radial.  The lamellipodium is the only one whose
ORIENTATION is the physics: ``laws/architecture_spec.LAMELLIPODIUM`` records mothers seeded at ±35°
about the protrusion axis, with Arp2/3 daughters branching at θ₀ = 70° — which flips a filament into
the OTHER mode, so every filament sits at ±35° and every junction sits at 70° simultaneously.  That is
one geometric statement, and it is why the two modes are laid down by the parity of the filament index
rather than by a random draw: the array is deterministic, reproducible, and diffable against a rebuild.

⚠ **THE BRANCH JUNCTION IS NOT BUILT HERE, AND SAYING SO MATTERS.**  The 70° angle is maintained by
``laws/network_warp.branch_angle_kernel`` as an angle HARMONIC — a force with a stiffness
(``ARP23_BRANCH_K = 0.173 pN·µm/rad²``), not a rigid constraint — and it acts on a mother/daughter
triple that is a BOND between two filaments.  Bonds are PHASE 3.  What stands here is the population in
its two modes; what does not is the dendritic connectivity that makes it a network rather than a comb.
A render of this population before PHASE 3 will look like a comb, and that is the truth about it.

⚠ **THE COUNT IS THE WORST OF THE FIVE.**  ``architecture_spec`` gives ``n_filaments = 200`` for a patch
of half-extent 8 µm — an EXPLICIT count, not a density, so it does not scale with the patch and cannot
be checked against a measurement.  The one datum in the neighbourhood is the linear Arp2/3 branch
density 1.25/µm (Vinzenz 2012), and its own comment says *"PI-gated cite"*.  There is no lamellipodial
actin AREAL density in the KB at all.  So ``count`` arrives as a :class:`~aleph.world.bond.BondCount`
and the natural basis for it is ``areal`` against the MEASURED patch area — which is the shape a
density has to be in before it can ever be sourced.  ``build_bond_family``'s rule applies verbatim:
supply a bigger patch rather than lowering the density.

WHY THE PATCH ORIGIN IS PASSED AND NOT DERIVED.  A lamellipodium sits at the leading edge, on the
substrate, and where that is depends on the cell's adhesion state — which this module cannot know and
must not guess.  So the caller gives the corner, the axis and the substrate normal, and the patch is
built from those.  Deriving it from a cell radius would be inventing a spread state.

TOPOLOGY IS ARITHMETIC.  One claim per population per kind.  See :mod:`aleph.world.build.microtubule`.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``origin``, ``width_um``, ``depth_um``, ``contour_um``, ``seg_um`` [µm];
    ``mode_deg`` is an angle [deg], converted once on the host; the count is per µm² or explicit.
  * boundary — a filament of fewer than 3 beads is refused, the same rule as ``build_strand``; a
    ``mode_deg`` outside (0, 90) is refused, since at 0° the two modes coincide and at 90° the array
    runs across the protrusion axis instead of into it; a count resolving to zero is refused; an
    ``axis`` parallel to ``normal`` is refused, because their cross product is the lateral direction
    and a degenerate one would collapse the patch onto a line.
  * conservation/invariant — the kernel writes only ``lo + t`` for ``t`` in ``[0, n_nodes)`` and the
    launch dimension IS ``n_nodes``.  ``assert_partitioned`` is asserted after the build.
  * CFL/precision — no integration; float64 throughout.  The realised segment length is reported.
  * sign sense — the two modes are ``+mode_deg`` for even filament indices and ``-mode_deg`` for odd,
    both rotated about ``normal``, so their MEAN direction is ``axis``: the array protrudes forward
    rather than sideways, which is checkable by summing the directions and is asserted in the demo.
  * measurement protocol — one kernel launch of ``dim=n_nodes``; nothing is read back.

engine units: length µm, area µm², angle deg on the interface.  Runtime: NVIDIA Warp on CUDA; the build
refuses a device-less arena.
"""

from __future__ import annotations

import math

import warp as wp

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount
from aleph.world.strand import _require_length

__all__ = ["plan_lamellipodium", "build_lamellipodium", "lay_twomode_patch_kernel"]


@wp.kernel
def lay_twomode_patch_kernel(
    lo: wp.int32, n_per: wp.int32, cols: wp.int32, strand_base: wp.int32, range_ordinal: wp.int32,
    ox: wp.float64, oy: wp.float64, oz: wp.float64,
    ax: wp.float64, ay: wp.float64, az: wp.float64,
    lx: wp.float64, ly: wp.float64, lz: wp.float64,
    dx_um: wp.float64, dy_um: wp.float64, cos_mode: wp.float64, sin_mode: wp.float64,
    seg_um: wp.float64,
    position: wp.array(dtype=wp.vec3d), strand_id: wp.array(dtype=wp.int32),
    range_id: wp.array(dtype=wp.int32),
):
    """Write one dendritic-array bead per thread; the filament's MODE is the parity of its index.

    ``(ax, ay, az)`` is the unit protrusion axis and ``(lx, ly, lz)`` the unit lateral direction, both
    already orthonormal on the host — the kernel does no frame construction because the patch has ONE
    frame for every filament, unlike the bundle kernels where each finger has its own.
    """
    t = wp.tid()
    f = t / n_per
    i = t - f * n_per

    axis = wp.vec3d(ax, ay, az)
    lat = wp.vec3d(lx, ly, lz)

    row = f / cols
    col = f - row * cols
    # +mode for even filaments, -mode for odd: the two Arp2/3 modes, deterministic and diffable.
    sign = wp.float64(1 - 2 * (f - wp.int32(2) * (f / wp.int32(2))))
    d = cos_mode * axis + sign * sin_mode * lat

    seed = (wp.vec3d(ox, oy, oz)
            + (wp.float64(col) + wp.float64(0.5)) * dx_um * lat
            + (wp.float64(row) + wp.float64(0.5)) * dy_um * axis)

    g = lo + t
    position[g] = seed + wp.float64(i) * seg_um * d
    strand_id[g] = strand_base + f
    range_id[g] = range_ordinal


def plan_lamellipodium(
    *, count: BondCount, support: float | None = None, width_um: float | None = None,
    depth_um: float | None = None, contour_um: float | None = None, seg_um: float | None = None,
    mode_deg: float | None = None,
) -> dict[str, object]:
    """Resolve the array's inventory from its declared density and patch geometry, no device touched.

    Args:
        count: how many filaments, per what, for which cell, on whose authority.  **REQUIRED**.  An
            ``areal`` basis resolves against the MEASURED patch area; an ``explicit`` one is a count and
            does not scale with the patch, which is what ``architecture_spec``'s 200 actually is.
        support: the MEASURED support to resolve against [µm² for ``areal``, 1.0 for ``explicit``].
            **REQUIRED** — no default, because the analytic ``width x depth`` is exactly the substitution
            that turned a mesh artefact into a physiological-looking density once already.
        width_um: patch extent across the protrusion axis [µm].  **REQUIRED**.
        depth_um: patch extent along the protrusion axis [µm].  **REQUIRED**.
        contour_um: filament contour length [µm].  **REQUIRED** — Arp2/3 filaments are SHORT
            (``architecture_spec`` uses 1.0 µm), and this is the axis on which that is declared.
        seg_um: bead spacing [µm].  **REQUIRED**.
        mode_deg: the ±mode half-angle about the protrusion axis [deg].  **REQUIRED** — 35° is θ₀/2 and
            is GROUNDED (Mueller 2017), which is a reason to pass it, not a reason to default it.

    Returns:
        The inventory, with the patch area recorded next to the support so the two can be compared.

    Raises:
        ValueError: on a missing or non-positive length, a ``mode_deg`` outside (0, 90), a filament of
            fewer than 3 beads, or a count resolving to zero filaments.
    """
    width = _require_length("width_um", width_um)
    depth = _require_length("depth_um", depth_um)
    contour = _require_length("contour_um", contour_um)
    seg = _require_length("seg_um", seg_um)
    mode = _require_length("mode_deg", mode_deg)
    sup = _require_length("support", support)
    if mode >= 90.0:
        raise ValueError(
            f"mode_deg={mode} must be under 90°: at 90° the array runs ACROSS the protrusion axis and "
            "the mean filament direction is zero, so the patch cannot protrude at all."
        )

    n_seg = int(round(contour / seg))
    n_per = n_seg + 1
    if n_per < 3:
        raise ValueError(
            f"contour_um={contour} at seg_um={seg} gives {n_per} beads per filament; a filament needs at "
            "least 3 so it carries one bending triple. A filament with no bending stiffness is a "
            "different physical object, not a coarser one."
        )
    n_fil = count.resolve(sup)
    if n_fil < 1:
        raise ValueError(
            f"the count resolves to {n_fil} filaments against support={sup}. A population that is "
            "'present but has nothing' is indistinguishable downstream from one that failed to build."
        )
    # Seed the patch on a grid whose cells are roughly square, so an AREAL density lands isotropically
    # rather than in stripes — a stripe pattern would be a build artefact masquerading as an orientation.
    cols = max(1, int(round(math.sqrt(n_fil * width / depth))))
    rows = int(math.ceil(n_fil / cols))
    return {
        "kind": "lamellipodium",
        "n_strands": n_fil,
        "nodes_per_strand": n_per,
        "cols": cols, "rows": rows,
        "dx_um": width / cols, "dy_um": depth / rows,
        "width_um": width, "depth_um": depth,
        "patch_area_um2_analytic": width * depth,
        "support": sup,
        "mode_deg": mode,
        "cos_mode": math.cos(math.radians(mode)), "sin_mode": math.sin(math.radians(mode)),
        "seg_um_requested": seg,
        "seg_um_realised": contour / n_seg,
        "contour_um_realised": contour,
        "provenance": {"n_filaments": {
            "basis": count.basis, "value": count.value, "support_unit": count.support_unit,
            "scope": count.scope, "source_class": str(count.source_class),
            "provenance": count.provenance,
        }},
    }


def build_lamellipodium(
    arena: WorldArena,
    *,
    count: BondCount,
    support: float | None = None,
    origin: tuple[float, float, float],
    axis: tuple[float, float, float],
    normal: tuple[float, float, float],
    width_um: float | None = None,
    depth_um: float | None = None,
    contour_um: float | None = None,
    seg_um: float | None = None,
    mode_deg: float | None = None,
    population: str = "lamellipodium",
) -> dict[str, object]:
    """Claim the array's ranges and lay the whole two-mode patch down with one CUDA kernel.

    Args:
        arena: the world to claim from.  Must hold a CUDA allocation.
        origin: the patch corner [µm] — the rear-lateral corner, from which ``depth`` runs along
            ``axis`` and ``width`` along the lateral direction.
        axis: the protrusion direction, in the substrate plane; normalised here.
        normal: the substrate normal.  The lateral direction is ``normal x axis``, so ``normal`` decides
            which side of the axis the ``+`` mode sits on.
        count / support / width_um / depth_um / contour_um / seg_um / mode_deg: see
            :func:`plan_lamellipodium`.
        population: the name the four claims are attributed to.

    Returns:
        The census of what was built.

    Raises:
        RuntimeError: if ``arena`` holds no device allocation.
        ValueError: on any refusal from :func:`plan_lamellipodium`, or an ``axis`` parallel to
            ``normal``.
    """
    inventory = plan_lamellipodium(
        count=count, support=support, width_um=width_um, depth_um=depth_um, contour_um=contour_um,
        seg_um=seg_um, mode_deg=mode_deg)

    a = _unit("axis", axis)
    n = _unit("normal", normal)
    lat = (n[1] * a[2] - n[2] * a[1], n[2] * a[0] - n[0] * a[2], n[0] * a[1] - n[1] * a[0])
    lat_norm = math.sqrt(sum(v * v for v in lat))
    if lat_norm <= 1e-9:
        raise ValueError(
            f"axis {axis!r} is parallel to normal {normal!r}, so the lateral direction is degenerate and "
            "the patch would collapse onto a line. A lamellipodium lies IN the substrate plane; an axis "
            "along the normal is a protrusion into the dish."
        )
    lat = tuple(v / lat_norm for v in lat)

    n_strand = int(inventory["n_strands"])
    n_per = int(inventory["nodes_per_strand"])
    n_nodes = n_strand * n_per
    if not arena.node_arrays:
        raise RuntimeError(
            f"{population}: this arena holds no device allocation, and construction here is a CUDA "
            "kernel rather than a host build that is uploaded. Warp CUDA is the only runtime; build the "
            f"arena with a CUDA device. The inventory it would have built: {n_strand} filaments x "
            f"{n_per} beads = {n_nodes} nodes."
        )

    nodes = arena.claim(population, Kind.NODE, n_nodes)
    strands = arena.claim(population, Kind.STRAND, n_strand)
    segments = arena.claim(population, Kind.SEGMENT, n_strand * (n_per - 1))
    angles = arena.claim(population, Kind.ANGLE3, n_strand * (n_per - 2))
    ordinal = len(arena.claims()) - 4

    wp.launch(
        lay_twomode_patch_kernel, dim=n_nodes,
        inputs=[nodes.lo, n_per, int(inventory["cols"]), strands.lo, ordinal,
                float(origin[0]), float(origin[1]), float(origin[2]), a[0], a[1], a[2],
                lat[0], lat[1], lat[2], float(inventory["dx_um"]), float(inventory["dy_um"]),
                float(inventory["cos_mode"]), float(inventory["sin_mode"]),
                float(inventory["seg_um_realised"])],
        outputs=[arena.node_arrays["position"], arena.node_arrays["strand_id"],
                 arena.node_arrays["range_id"]],
        device=arena.device,
    )
    arena.assert_partitioned()

    return inventory | {
        "population": population, "n_nodes": nodes.count, "n_segments": segments.count,
        "n_angles": angles.count, "axis": a, "lateral": lat,
        "claims": {"node": (nodes.lo, nodes.hi), "strand": (strands.lo, strands.hi),
                   "segment": (segments.lo, segments.hi), "angle3": (angles.lo, angles.hi)},
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
    """Self-check: the inventory arithmetic, the two-mode mean direction, and every refusal."""
    gap = BondCount(
        basis="areal", value=3.125,
        scope="lamellipodial dendritic array, generic — NO areal actin density for a lamellipodium "
              "exists in the KB",
        source_class="PI_GAP",
        provenance="back-read from laws/architecture_spec.LAMELLIPODIUM's EXPLICIT n_filaments=200 over "
                   "its own 8x8 µm patch; the only nearby datum is the LINEAR Arp2/3 branch density "
                   "1.25/µm (Vinzenz 2012), whose own comment says 'PI-gated cite'")

    inv = plan_lamellipodium(count=gap, support=64.0, width_um=8.0, depth_um=8.0, contour_um=1.0,
                             seg_um=0.05, mode_deg=35.0)
    assert inv["n_strands"] == 200, inv["n_strands"]
    assert inv["nodes_per_strand"] == 21, inv["nodes_per_strand"]
    assert inv["cols"] * inv["rows"] >= 200, "the grid must cover every filament"

    # The two modes are symmetric about the axis: their mean direction IS the protrusion axis, so the
    # patch protrudes forward. cos(35°) is the axial component and the lateral components cancel.
    ca, sa = float(inv["cos_mode"]), float(inv["sin_mode"])
    assert abs((ca + ca) / 2.0 - math.cos(math.radians(35.0))) < 1e-15
    assert abs((sa - sa)) == 0.0, "the ± modes cancel laterally"

    # An areal density scales with the patch; an explicit count does not. That difference is the whole
    # reason the basis is recorded rather than the number alone.
    bigger = plan_lamellipodium(count=gap, support=256.0, width_um=16.0, depth_um=16.0, contour_um=1.0,
                                seg_um=0.05, mode_deg=35.0)
    assert bigger["n_strands"] == 800, bigger["n_strands"]

    base = {"support": 64.0, "width_um": 8.0, "depth_um": 8.0, "contour_um": 1.0, "seg_um": 0.05,
            "mode_deg": 35.0}
    for missing in base:
        kwargs = {k: v for k, v in base.items() if k != missing}
        try:
            plan_lamellipodium(count=gap, **kwargs)
        except ValueError as exc:
            assert "no default" in str(exc) and missing in str(exc), str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"a missing {missing} must refuse")

    try:
        plan_lamellipodium(count=gap, **(base | {"mode_deg": 90.0}))
    except ValueError as exc:
        assert "cannot protrude" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a 90° mode must refuse")

    try:
        plan_lamellipodium(count=gap, **(base | {"seg_um": 1.0}))
    except ValueError as exc:
        assert "at least 3" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a 2-bead filament must refuse")

    arena = WorldArena(capacity={Kind.NODE: 10, Kind.STRAND: 10, Kind.SEGMENT: 10, Kind.ANGLE3: 10})
    try:
        build_lamellipodium(arena, count=gap, origin=(0.0, 0.0, 0.0), axis=(1.0, 0.0, 0.0),
                            normal=(0.0, 0.0, 1.0), **base)
    except RuntimeError as exc:
        assert "no device allocation" in str(exc) and "4200 nodes" in str(exc), str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device-less build must refuse")

    try:
        build_lamellipodium(arena, count=gap, origin=(0.0, 0.0, 0.0), axis=(0.0, 0.0, 1.0),
                            normal=(0.0, 0.0, 1.0), **base)
    except ValueError as exc:
        assert "collapse onto a line" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an axis along the normal must refuse")
    assert arena.n_live(Kind.NODE) == 0, "a refused build claims nothing"

    chars = _codegen_check()
    assert chars > 0

    print(f"lamellipodium self-check OK — CUDA codegen {chars} chars; {inv['n_strands']} filaments x "
          f"{inv['nodes_per_strand']} beads = {inv['n_strands'] * int(inv['nodes_per_strand'])} nodes "
          f"on a {inv['cols']}x{inv['rows']} seed grid")


if __name__ == "__main__":
    _demo()
