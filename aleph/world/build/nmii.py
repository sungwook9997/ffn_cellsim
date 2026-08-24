r"""NMII motor heads: the population that is missing, and the bond that is not built here.

WHY THE ACTIVE LOAD PATH IS INERT.  ``lamellipodium``, ``filopodium`` and ``sf_arc`` are declared
``ComponentRole.ACTIVE_LOAD_PATH`` (``engine/contracts.py``), but an actin bundle is not itself
active — the activity is NMII, and in the arena the NMII population is **zero**.  The frozen
incumbent carries 442 bipolar minifilaments (``components/motor/force_budget_ledger.py``:
``CURRENT_N_MF``), so every head that would load those bundles exists only in the port source.  This
module builds the population; it does not make it pull.

⚠ **THIS IS A POPULATION AND A COUPLING, AND ONLY THE POPULATION LANDS HERE.**
  1. minifilament backbone + heads as a NODE population — this module;
  2. head → actin is a **BOND**, and ``world/bond.py`` says in as many words that kinetics, attach /
     detach, free lists and snapshot twins "arrive with the first family that has kinetics".  The NMII
     crossbridge IS that first family.  Choosing its shape is therefore a gate-contract decision and
     is **surfaced to the PI, not invented here**.  ``build_nmii`` claims no ``Kind.BOND`` and reports
     ``n_bonds: 0`` with the reason attached, so the absence is a recorded refusal rather than an
     oversight.  The head↔backbone arm is NOT that bond: an arm never breaks, so it is a SEGMENT —
     structure — exactly as ``bond.py`` draws the line.

TOPOLOGY IS PORTED BY READING, NOT BY IMPORT.  ``components/motor/minifilament_topology.py`` already
holds the head-resolved Stam-Hocky bipolar layout and its invariants, but ``world/`` may not import
``components/`` (``tests/architecture/test_layer_directions.py``), so the index arithmetic is restated
here and its self-check asserts the same invariants that oracle asserts:

    one minifilament = ``[bb_0 .. bb_{n_bb-1}, head+_0 .. head+_{H-1}, head-_0 .. head-_{H-1}]``
    n_particles     = n_bb + 2H
    segments        = (n_bb - 1) backbone chain + 2H head arms
    angle3          = (n_bb - 2) backbone bending + 2H head-arm orientation
    anchor(head+_i) = (i * half) // H,  anchor(head-_i) = n_bb - 1 - (i * half) // H,
                      half = max(n_bb // 2, 1)          — the axial bipolar split

⚠ **FOUR PI-GAPS.  EVERY ONE IS A REQUIRED ARGUMENT WITH NO DEFAULT, AND THE REFUSAL NAMES BOTH SIDES.**
  * ``n_heads_per_side`` — **10** (AFINES) vs **28–30** (Billington 2013).  The 2026-08-20 source audit
    found "30 motors per filament end" in Melli 2018 *eLife* and Tripathi 2021 *JoVE*, which SUPPORTS
    the ``claim_b`` side; it does not close the gap, and this module does not close it either.  It is
    the single largest lever on head count: at 442 minifilaments the choice is 8,840 vs 24,752 heads.
  * ``backbone_length_um`` — 0.301 µm (Billington 2013) is the archived figure and is PI-GAP.
  * ``head_offset_um`` — 0.200 µm is ARCHIVED, sourced to nothing.
  * ``n_bb`` — 14 is the archived HOOMD figure and is PI-GAP.

  ``F_stall_head`` (0.5 vs 2.0 pN) is **deliberately not an argument here.**  The 2026-08-20 source
  audit found that *neither value is a measurement*.  It is a property of the crossbridge BOND's force
  law, and this module builds no force law, so accepting it would be storing an unratified magnitude in
  a geometry builder.  It arrives with item 2 above, at the PI.

WHERE THE MINIFILAMENTS GO, AND WHY ``e_y`` IS TANGENT.  Each minifilament is placed on the cortical
shell with a random radial direction, a random backbone axis in the tangent plane, and the head-offset
direction ``e_y = ĉ × ê_x`` — also tangent.  Both head rows therefore stay in the shell instead of one
row being pushed out through the membrane and the other in through the cytosol.  This is a PLACEMENT,
not a straddle: which anti-parallel actin pair each side binds is decided by the bond (item 2), and
writing a straddle now would answer that question with a build detail.  ``minifilament_topology.
straddle_frame`` is the port target for when it is answered.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``backbone_length_um``, ``head_offset_um``, ``radius_um``, ``thickness_um``,
    ``centre_um`` [µm]; the count is areal [µm⁻²] resolved against a MEASURED shell area [µm²]; the
    bead and head counts are dimensionless.  No stiffness, no rate, no force anywhere in this module.
  * boundary — ``n_bb < 2`` is refused (a backbone is a chain); ``n_heads_per_side < 1`` is refused (a
    bipolar filament with no heads is a rod); a count resolving to zero is refused, because a
    population that is present but empty is indistinguishable downstream from one that failed to
    build; a node claim exceeding arena capacity raises through ``arena.claim``.
  * conservation/invariant — the kernel writes only ``lo + t`` for ``t`` in ``[0, n_nodes)`` and the
    launch dimension IS ``n_nodes``; ``assert_partitioned`` is asserted after the build.  Bipolar
    symmetry is structural: ``H`` heads at ``+offset`` and ``H`` at ``−offset`` about the SAME axis, so
    the net transverse offset is exactly zero by construction rather than to a tolerance.
  * CFL/precision — no integration; float64 throughout.  The backbone bead spacing
    ``L_bb/(n_bb−1)`` is reported, because that is the rest length the backbone bond law must read.
  * sign sense — ``+`` heads sit on ``+e_y`` and ``−`` heads on ``−e_y``, anchored to the two opposite
    halves of the backbone.  That opposition IS the bipolarity; a same-side placement would build a
    filament that translates rather than contracts.  No polarity or walk direction is written: those
    belong to the crossbridge, not to the geometry.
  * measurement protocol — one kernel launch of ``dim = n_nodes``; nothing is read back.

engine units: length µm.  Runtime: NVIDIA Warp on CUDA; the build refuses a device-less arena.
"""

from __future__ import annotations

import warp as wp

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount
from aleph.world.strand import _require_length

__all__ = ["NMII_PI_GAPS", "plan_nmii", "build_nmii", "lay_minifilaments_kernel"]

#: The unresolved magnitude behind each required argument, quoted in the refusal so a caller cannot
#: pick a side without reading which sides there are.  Keyed by argument name.
NMII_PI_GAPS: dict[str, str] = {
    "n_bb": "backbone beads: 14 is the archived HOOMD figure (cortex/myosin.py, deleted at 1a9ded66) "
            "and is an I0-B3 PI-GAP — no Contract-Graph point value.",
    "n_heads_per_side": "heads per side H: claim_a=10 (AFINES) vs claim_b=28 (Billington 2013), "
                        "params_i0b3.yaml:68,72. The 2026-08-20 source audit found '30 motors per "
                        "filament end' in Melli 2018 eLife and Tripathi 2021 JoVE, which SUPPORTS "
                        "claim_b — it does not close the gap. Declare a side with its scope.",
    "backbone_length_um": "backbone contour L_bb: 0.301 µm (Billington 2013) is the archived figure "
                          "and is an I0-B3 PI-GAP.",
    "head_offset_um": "perpendicular head↔backbone offset r0_head: 0.200 µm is ARCHIVED and sourced "
                      "to nothing. PI-GAP.",
}

#: What ``build_nmii`` refuses to build, and why the refusal is the correct outcome rather than a gap
#: in the implementation.  Emitted in the census so the absence travels with the artifact.
_BOND_REFUSAL = (
    "head→actin crossbridge. This is the FIRST bond family with kinetics, and world/bond.py defers "
    "attach/detach, the free list and the snapshot twins to exactly that family. Its shape is a "
    "gate-contract decision (CLAUDE.md §Working model) and F_stall_head — 0.5 vs 2.0 pN — was found by "
    "the 2026-08-20 source audit to be a measurement on NEITHER side. Surfaced to the PI, not built."
)


@wp.kernel
def lay_minifilaments_kernel(
    seed: wp.int32, lo: wp.int32, n_per: wp.int32, n_bb: wp.int32, n_side: wp.int32,
    strand_base: wp.int32, range_ordinal: wp.int32,
    radius: wp.float64, thickness: wp.float64, centre: wp.vec3d,
    bb_seg_um: wp.float64, offset_um: wp.float64,
    position: wp.array(dtype=wp.vec3d), strand_id: wp.array(dtype=wp.int32),
    range_id: wp.array(dtype=wp.int32),
):
    """Write one bead — backbone or head — per thread, with the placement frame re-derived from ``tid``.

    One thread per NODE, not per minifilament: the frame comes from ``wp.rand_init(seed, m)`` so all
    ``n_per`` threads of one minifilament agree on it without any of them writing a shared value, and
    no per-minifilament host array is ever built.  Same idiom as
    :func:`aleph.world.build.cortex._strand_positions`, for the same construction ruling.

    A head is placed by ANCHOR BEAD rather than by its own axial coordinate, which is what keeps the
    arm perpendicular and of length ``offset_um`` exactly: the head-arm angle spring then rests
    force-free at ``π/2`` and the backbone bending spring at ``π``, so the built minifilament starts
    unstrained in every internal law rather than in the axial one alone.
    """
    t = wp.tid()
    m = t / n_per
    k = t - m * n_per

    state = wp.rand_init(seed, m)
    c32 = wp.sample_unit_sphere_surface(state)
    c = wp.normalize(wp.vec3d(wp.float64(c32[0]), wp.float64(c32[1]), wp.float64(c32[2])))
    a32 = wp.sample_unit_sphere_surface(state)
    ex = wp.vec3d(wp.float64(a32[0]), wp.float64(a32[1]), wp.float64(a32[2]))
    ex = ex - wp.dot(ex, c) * c
    # A backbone axis drawn parallel to its own radial direction has no tangent-plane part. It cannot
    # happen at any rate worth naming, and normalising a zero would put a NaN into the arena, so the
    # degenerate draw takes a fixed axis instead of a random one.
    if wp.length(ex) < wp.float64(1.0e-9):
        ex = wp.vec3d(wp.float64(1.0), wp.float64(0.0), wp.float64(0.0))
        ex = ex - wp.dot(ex, c) * c
        if wp.length(ex) < wp.float64(1.0e-9):
            ex = wp.vec3d(wp.float64(0.0), wp.float64(1.0), wp.float64(0.0))
            ex = ex - wp.dot(ex, c) * c
    ex = wp.normalize(ex)
    ey = wp.cross(c, ex)  # unit already: c ⟂ ex and both are unit. Tangent, so both head rows stay in the shell.

    # `thickness` here is the INSET draw span, not the shell thickness — see build_nmii's launch.
    r = radius + thickness * (wp.float64(wp.randf(state)) - wp.float64(0.5))
    origin = c * r + centre

    half = wp.max(n_bb / wp.int32(2), wp.int32(1))
    bead = k
    side = wp.float64(0.0)
    if k >= n_bb:
        h = k - n_bb
        if h < n_side:
            bead = (h * half) / n_side
            side = wp.float64(1.0)
        else:
            hh = h - n_side
            bead = n_bb - wp.int32(1) - (hh * half) / n_side
            side = wp.float64(-1.0)

    axial = (wp.float64(bead) - wp.float64(0.5) * wp.float64(n_bb - 1)) * bb_seg_um
    g = lo + t
    position[g] = origin + axial * ex + side * offset_um * ey
    strand_id[g] = strand_base + m
    range_id[g] = range_ordinal


def _require_gap(name: str, value: object) -> float:
    """Return a declared PI-GAP value, or refuse with both sides of the gap named.

    ``_require_length``'s refusal says a value must be declared.  This one also says WHAT is unresolved
    about it, because for these four the caller's next question is "declared as what?" and the answer
    is not a number this repository is allowed to choose.
    """
    if value is None:
        raise ValueError(f"{name} has no default and must be declared. PI-GAP — {NMII_PI_GAPS[name]}")
    return float(value)


def _require_int(name: str, value: object, minimum: int, why: str) -> int:
    """Return a declared count at or above ``minimum``, or refuse with the structural reason."""
    n = int(_require_gap(name, value))
    if n < minimum:
        raise ValueError(f"{name} must be >= {minimum}: {why} Got {n}.")
    return n


def plan_nmii(
    *, count: BondCount, support: float, n_bb: int | None = None,
    n_heads_per_side: int | None = None, backbone_length_um: float | None = None,
    head_offset_um: float | None = None, radius_um: float | None = None,
    thickness_um: float | None = None,
) -> dict[str, object]:
    """Resolve the minifilament inventory from its declared count and topology, no device touched.

    Args:
        count: how many bipolar minifilaments this cell has.  **REQUIRED**, and ⚠ the only DIRECT
            cortical NMII density in the corpus is Nie 2015, 0.625/µm²; the 16–21/µm² "physiological"
            band is a back-calculated GAP target, not a measurement (``force_budget_ledger``).  Which
            of those a run uses is the caller's declaration, carried with its scope.
        support: the MEASURED support ``count`` resolves against — the cortical shell area [µm²] for
            an ``areal`` count, 1.0 for an ``explicit`` one.  Measured, never analytic.
        n_bb: backbone beads per minifilament.  **REQUIRED**, PI-GAP — see :data:`NMII_PI_GAPS`.
        n_heads_per_side: heads ``H`` on each anti-parallel side.  **REQUIRED**, PI-GAP, and the
            largest single lever on the head count.
        backbone_length_um: backbone contour length [µm].  **REQUIRED**, PI-GAP.
        head_offset_um: perpendicular head↔backbone offset [µm].  **REQUIRED**, PI-GAP.
        radius_um: mid-shell radius the minifilaments are placed on [µm].  **REQUIRED**.
        thickness_um: shell thickness [µm]; each minifilament takes its own radius in
            ``radius ± thickness/2``, which is the layer being a layer rather than a sphere.
            **REQUIRED**.

    Returns:
        The inventory, with the provenance row, the derived backbone rest length, and the bond
        refusal.

    Raises:
        ValueError: on any undeclared PI-GAP, a backbone of fewer than 2 beads, a side with no heads,
            a count resolving to zero, or a geometry whose beads leave the shell.
    """
    n_bb_v = _require_int("n_bb", n_bb, 2, "a backbone is a chain, and a chain needs two ends.")
    h = _require_int("n_heads_per_side", n_heads_per_side, 1,
                     "a bipolar minifilament with no heads on a side is a rod, not a motor.")
    l_bb = _require_gap("backbone_length_um", backbone_length_um)
    offset = _require_gap("head_offset_um", head_offset_um)
    for name, v in (("backbone_length_um", l_bb), ("head_offset_um", offset)):
        if not (v > 0.0):
            raise ValueError(f"{name} must be finite and positive; got {v!r}")
    radius = _require_length("radius_um", radius_um)
    thickness = _require_length("thickness_um", thickness_um)

    n_mf = count.resolve(support)
    if n_mf < 1:
        raise ValueError(
            f"the count resolves to {n_mf} minifilaments. A population that is 'present but has "
            "nothing' is indistinguishable downstream from one that failed to build — and an NMII "
            "population of zero is exactly the state the arena is already in."
        )

    # The backbone is a straight CHORD on a curved shell (beads bow inward by the sagitta) and the
    # heads sit a tangent offset out (bulging outward by offset²/2r). Both are second order and tiny
    # at native, but they come from different sources and nothing else in the build would notice them
    # crossing the shell — the same argument the stress-fibre pitch guard is made of.
    excursion = l_bb * l_bb / (8.0 * radius) + offset * offset / (2.0 * radius)
    if excursion > 0.5 * thickness:
        raise ValueError(
            f"a minifilament of L_bb={l_bb} µm with a {offset} µm head offset cannot be INSET into its "
            f"shell: the chord sagitta plus the head bulge is {excursion:.4g} µm against a "
            f"half-thickness of {0.5 * thickness:.4g} µm at radius {radius} µm, so there is no span "
            "left to draw a filament radius from. Beads outside the cortical layer are in the "
            "membrane or in the cytosol, which are other components' claims."
        )
    # ⚠ THE DRAW IS INSET BY THE EXCURSION, AND UNTIL 2026-08-22 IT WAS NOT.
    #
    # The guard above tested FIT — "does a minifilament fit inside the shell" — and the kernel then
    # drew each filament's radius from the FULL thickness, so a filament drawn at the outer edge put
    # its beads that excursion PAST the shell whatever the guard returned. The guard could not fail
    # for the thing it appeared to protect, which is the family `geometry.py`'s own self-check names
    # by pointing at `balance_ok` and `descent_ratio`.
    #
    # Measured before the fix, over a ten-point radius ladder at full native: the outermost node sat a
    # CONSTANT ~3.3 nm past `radius + thickness/2` at every radius, against this formula's 4.2 nm — so
    # the magnitude was never wrong and insetting by it is conservative. At radius 7.40 that 3.3 nm was
    # 414 of 32,708 nodes outside the membrane and `assert_inside_membrane` refused the build.
    #
    # With the inset, the guard above becomes exactly "is there room to inset", which is a condition
    # that can actually bind — the same line, now load-bearing.
    draw_span = max(thickness - 2.0 * excursion, 0.0)

    n_per = n_bb_v + 2 * h
    return {
        "kind": "nmii",
        "n_minifilaments": n_mf,
        "n_bb": n_bb_v,
        "n_heads_per_side": h,
        "n_heads": n_mf * 2 * h,
        "nodes_per_minifilament": n_per,
        "n_strands": n_mf,
        "segments_per_minifilament": (n_bb_v - 1) + 2 * h,
        "angles_per_minifilament": max(n_bb_v - 2, 0) + 2 * h,
        "backbone_length_um": l_bb,
        # The rest length the backbone bond law must read — reported, never chosen from a target.
        "backbone_seg_um": l_bb / (n_bb_v - 1),
        "head_offset_um": offset,
        "radius_um": radius, "thickness_um": thickness,
        "shell_excursion_um": excursion,
        #: The span a filament radius is actually drawn from — the shell thickness less the
        #: excursion at each edge. Recorded so a reader can see the inset rather than infer it.
        "draw_span_um": draw_span,
        "n_bonds": 0,
        "bond_refusal": _BOND_REFUSAL,
        "pi_gaps": dict(NMII_PI_GAPS),
        "provenance": {"n_minifilaments": {
            "basis": count.basis, "value": count.value, "support_unit": count.support_unit,
            "scope": count.scope, "source_class": str(count.source_class),
            "provenance": count.provenance}},
    }


def build_nmii(
    arena: WorldArena,
    *,
    count: BondCount,
    support: float,
    centre_um: tuple[float, float, float] = (0.0, 0.0, 0.0),
    n_bb: int | None = None,
    n_heads_per_side: int | None = None,
    backbone_length_um: float | None = None,
    head_offset_um: float | None = None,
    radius_um: float | None = None,
    thickness_um: float | None = None,
    seed: int = 0,
    population: str = "nmii",
) -> dict[str, object]:
    """Claim the NMII ranges and lay every bipolar minifilament down with one CUDA kernel.

    Args:
        arena: the world to claim from.  Must hold a CUDA allocation.
        centre_um: the shell centre [µm] — recorded, because a population whose position is not in its
            record cannot be placed against another population afterwards.
        seed: the device RNG seed the placement frames are drawn from.  Recorded, so the build is
            reproducible from the artifact alone.
        population: the name the four claims are attributed to.
        count / support / n_bb / n_heads_per_side / backbone_length_um / head_offset_um / radius_um /
            thickness_um: see :func:`plan_nmii`.

    Returns:
        The census of what was built — including ``n_bonds: 0`` and the reason, which is the point.

    Raises:
        RuntimeError: if ``arena`` holds no device allocation.
        ValueError: on any refusal from :func:`plan_nmii`.
    """
    inventory = plan_nmii(
        count=count, support=support, n_bb=n_bb, n_heads_per_side=n_heads_per_side,
        backbone_length_um=backbone_length_um, head_offset_um=head_offset_um,
        radius_um=radius_um, thickness_um=thickness_um)

    n_mf = int(inventory["n_minifilaments"])
    n_per = int(inventory["nodes_per_minifilament"])
    n_nodes = n_mf * n_per
    if not arena.node_arrays:
        raise RuntimeError(
            f"{population}: this arena holds no device allocation, and construction here is a CUDA "
            "kernel rather than a host build that is uploaded. Warp CUDA is the only runtime; build "
            f"the arena with a CUDA device. The inventory it would have built: {n_mf} minifilaments x "
            f"{n_per} beads = {n_nodes} nodes carrying {inventory['n_heads']} heads."
        )

    nodes = arena.claim(population, Kind.NODE, n_nodes)
    strands = arena.claim(population, Kind.STRAND, n_mf)
    segments = arena.claim(population, Kind.SEGMENT, n_mf * int(inventory["segments_per_minifilament"]))
    angles = arena.claim(population, Kind.ANGLE3, n_mf * int(inventory["angles_per_minifilament"]))
    ordinal = len(arena.claims()) - 4

    wp.launch(
        lay_minifilaments_kernel, dim=n_nodes,
        inputs=[int(seed), nodes.lo, n_per, int(inventory["n_bb"]),
                int(inventory["n_heads_per_side"]), strands.lo, ordinal,
                # ⚠ draw_span, NOT thickness. The kernel draws a filament radius uniformly over
                # what it is given; giving it the full thickness is what put beads outside the shell.
                float(inventory["radius_um"]), float(inventory["draw_span_um"]),
                wp.vec3d(float(centre_um[0]), float(centre_um[1]), float(centre_um[2])),
                float(inventory["backbone_seg_um"]), float(inventory["head_offset_um"])],
        outputs=[arena.node_arrays["position"], arena.node_arrays["strand_id"],
                 arena.node_arrays["range_id"]],
        device=arena.device,
    )
    arena.assert_partitioned()

    return inventory | {
        "population": population, "n_nodes": nodes.count, "n_segments": segments.count,
        "n_angles": angles.count, "centre_um": tuple(float(c) for c in centre_um), "seed": int(seed),
        "claims": {"node": (nodes.lo, nodes.hi), "strand": (strands.lo, strands.hi),
                   "segment": (segments.lo, segments.hi), "angle3": (angles.lo, angles.hi)},
        "support": float(support),
    }


def _codegen_check() -> int:
    """Type-check this module's kernel by emitting its CUDA source, with no device and no execution.

    The dev machine has no card, so the alternative is shipping kernel syntax that nobody has put in
    front of a compiler; this generates the CUDA C++ Warp would compile, which type-checks every
    expression, and produces NO number — so it is not a CPU result under any reading.
    """
    from warp._src.context import ModuleBuilder  # private, deliberately: there is no public codegen API

    module = wp.get_module(__name__)
    return len(ModuleBuilder(module, module.options).codegen("cuda"))


def _anchor_beads(n_bb: int, h: int) -> tuple[list[int], list[int]]:
    """The host mirror of the kernel's anchor arithmetic — the one thing here that can be wrong.

    Restated rather than imported: ``world/`` may not import ``components/``
    (``tests/architecture/test_layer_directions.py``), so the port is a reading of
    ``minifilament_topology.head_backbone_bonds`` and this function is what makes the reading checkable.
    """
    half = max(n_bb // 2, 1)
    return ([(i * half) // h for i in range(h)],
            [n_bb - 1 - (i * half) // h for i in range(h)])


def _demo() -> None:
    """Self-check: the count invariants, the bipolar split, and every refusal — no device touched."""
    gap = dict(scope="MCF7 cortical shell — placeholder scope for the self-check",
               source_class="PI_GAP", provenance="self-check only; not a run configuration")
    native = dict(n_bb=14, n_heads_per_side=10, backbone_length_um=0.301, head_offset_um=0.200,
                  radius_um=7.40, thickness_um=0.20)

    # Every PI-GAP is required, and the refusal names both sides of the gap.
    for missing, expect in (("n_bb", "archived HOOMD"), ("n_heads_per_side", "claim_b"),
                            ("backbone_length_um", "Billington"), ("head_offset_um", "ARCHIVED")):
        try:
            plan_nmii(count=BondCount(basis="explicit", value=442.0, **gap), support=1.0,
                      **{**native, missing: None})
        except ValueError as exc:
            assert "PI-GAP" in str(exc) and expect in str(exc), (missing, exc)
        else:  # pragma: no cover
            raise AssertionError(f"{missing} must be refused, not defaulted")

    # A backbone of one bead is not a chain; a side with no heads is not a motor.
    for kwargs, expect in ((dict(n_bb=1), "a chain needs two ends"),
                           (dict(n_heads_per_side=0), "is a rod, not a motor")):
        try:
            plan_nmii(count=BondCount(basis="explicit", value=442.0, **gap), support=1.0,
                      **{**native, **kwargs})
        except ValueError as exc:
            assert expect in str(exc), (kwargs, exc)
        else:  # pragma: no cover
            raise AssertionError(f"plan_nmii({kwargs}) must refuse")

    # The incumbent's population, reproduced: 442 minifilaments, 34 beads, 8,840 heads.
    inv = plan_nmii(count=BondCount(basis="explicit", value=442.0, **gap), support=1.0, **native)
    assert inv["n_minifilaments"] == 442 and inv["nodes_per_minifilament"] == 34
    assert inv["n_heads"] == 8840, "442 x 2 x 10 — components/motor/force_budget_ledger.CURRENT_N_HEADS"
    assert inv["segments_per_minifilament"] == 13 + 20 and inv["angles_per_minifilament"] == 12 + 20
    assert abs(inv["backbone_seg_um"] - 0.301 / 13) < 1e-15
    assert inv["n_bonds"] == 0 and "PI" in str(inv["bond_refusal"])

    # claim_b is the lever, and the audit's Melli/Tripathi finding is on that side.
    inv_b = plan_nmii(count=BondCount(basis="explicit", value=442.0, **gap), support=1.0,
                      **{**native, "n_heads_per_side": 28})
    assert inv_b["n_heads"] == 24752, "the H gap is a 2.8x lever on the head population"

    # The only DIRECT cortical density resolves against a MEASURED shell area, not an analytic one.
    nie = BondCount(basis="areal", value=0.625, scope="cortical NMII minifilaments, Nie 2015",
                    source_class="SOURCED",
                    provenance="Nie 2015 — the only DIRECT cortical NMII areal density in the corpus; "
                               "the 16-21/um^2 band is a back-calculated GAP target, not a measurement")
    # 0.625/um^2 x 706.86 um^2 = 442 — the incumbent's CURRENT_N_MF is that density over that area,
    # so the population this module builds IS the one the frozen engine carries, not a new number.
    # 706.86 um^2 is the R=7.5 um MEMBRANE area the incumbent's density was taken against, while the
    # minifilaments sit on the 7.40 um cortex shell; force_budget_ledger flags that 1.4% mismatch
    # itself. A run supplies its own MEASURED shell area and gets a slightly smaller count.
    assert plan_nmii(count=nie, support=706.86, **native)["n_minifilaments"] == 442

    # A density low enough to build nothing is refused here even though BondCount allows it: an empty
    # NMII population is the state the arena is already in, and it must not look like a build.
    try:
        plan_nmii(count=BondCount(basis="areal", value=1e-9, **gap), support=1.0, **native)
    except ValueError as exc:
        assert "resolves to 0 minifilaments" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a count resolving to zero must refuse")

    # Geometry that leaves the shell is refused. The guard is deliberately loose: the offset is
    # TANGENT, so it costs only offset^2/2r radially — at native it spends 0.4% of the half-thickness,
    # and it takes ~1.2 um before an arm reaches out of a 0.2 um layer. It is a nonsense catcher.
    assert inv["shell_excursion_um"] < 0.05 * 0.5 * native["thickness_um"]

    # ⚠ THE INSET, PINNED AS AN INVARIANT RATHER THAN A COMMENT. The draw span plus an excursion at
    # each edge must not exceed the shell — which is the promise the guard above only APPEARED to make
    # until 2026-08-22, when the kernel was drawing from the full thickness and putting beads a
    # measured ~3.3 nm past the shell at every radius on a ten-point ladder.
    half_reach = 0.5 * inv["draw_span_um"] + inv["shell_excursion_um"]
    assert half_reach <= 0.5 * inv["thickness_um"] + 1e-12, (half_reach, inv["thickness_um"])
    # And the span is the whole shell less the two insets, not something smaller chosen for comfort.
    assert abs(inv["draw_span_um"] - (inv["thickness_um"] - 2.0 * inv["shell_excursion_um"])) < 1e-12
    # A shell with no room to inset must RAISE, not silently draw from a negative span.
    try:
        plan_nmii(count=BondCount(basis="explicit", value=442.0, **gap), support=1.0,
                  **{**native, "thickness_um": 2.0 * inv["shell_excursion_um"] * 0.9})
    except ValueError as exc:
        assert "cannot be INSET" in str(exc), exc
    else:
        raise AssertionError("a shell too thin to inset must be refused")
    try:
        plan_nmii(count=BondCount(basis="explicit", value=442.0, **gap), support=1.0,
                  **{**native, "head_offset_um": 2.0})
    except ValueError as exc:
        assert "cannot be INSET" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a minifilament that leaves the cortical layer must refuse")

    # The bipolar split: H anchors on each half, opposite ends, mirror images of one another.
    plus, minus = _anchor_beads(14, 10)
    assert len(plus) == len(minus) == 10
    assert max(plus) < 7 <= min(minus), "the two sides must anchor to opposite halves of the backbone"
    assert [13 - b for b in plus] == minus, "the split is mirror-symmetric, so the two clamps match"
    assert _anchor_beads(2, 1) == ([0], [1]), "the minimum backbone still splits"

    print(f"nmii self-check OK — {inv['n_minifilaments']} minifilaments, {inv['n_heads']} heads, "
          f"{inv['n_bonds']} bonds (refused, see census); codegen {_codegen_check()} chars")


if __name__ == "__main__":
    _demo()
