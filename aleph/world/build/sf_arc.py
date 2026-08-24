r"""Transverse arcs: curved, FA-FREE actomyosin bundles in the lamella, laid down by a CUDA kernel.

**WHY THIS FILE EXISTS AT ALL — PI decision 8, 2026-08-21** (`docs/v2_audit/PI_DECISIONS_2026-08-21.md`
§8): ``sf_arc`` is a DIFFERENT STRUCTURE, not another name for what PHASE 1 built.  Six contract
endpoints address ``sf_arc``; ``build/stress_fiber.py`` claims population ``stress_fiber``; and the
temptation all week was to rename one to the other.  The PI refused: *"Renaming the contract to match
the build would attach an existing name to a structure nobody built."*  So this is the builder for the
structure that was missing, and ``test_families_census``'s ``xfail(strict=True)`` comes off when it
exists — **by Lead, not here**.

WHAT MAKES IT A DIFFERENT STRUCTURE.  Four differences, none cosmetic, each visible in the geometry
this module writes:

1. **No focal adhesion at either end.**  A ventral fibre spans two FAs and that is its definition.  A
   transverse arc is FA-FREE: it is born by Tm4/myosin-II condensation of the lamellipodial network
   plus Arp2/3-filament annealing (`KB-DRAFT-7-07`, Tojkander 2011 *Curr Biol*, Burnette 2011 *Nat Cell
   Biol*) and its ends terminate on DORSAL fibres, not on the substrate.
2. **Curved, not straight.**  It runs PARALLEL to the leading edge, concave toward it — a circular arc
   concentric with the contact rim, which is what this kernel writes.
3. **Lifted off the basal plane.**  Ventral fibres lie ON the substrate; arcs sit in the lamella, above
   it, and flow rearward.  ``lift_um`` is where that goes and it has no default.
4. **The count is set by a lamellar SPACING, not by an FA pattern.**  Different question, different
   gap — see below.

⚠ **AND A THIRD STRUCTURE IS NAMED HERE PRECISELY SO IT IS NOT SILENTLY FOLDED IN.**  The DORSAL stress
fibre — formin-nucleated at an FA, rising out of the basal plane, radial rather than circumferential,
and **non-contractile** — is neither a ventral fibre nor a transverse arc.  PI decision 8 writes the
population as *"transverse/dorsal ARC"*; this module builds the TRANSVERSE arc and reports
``dorsal_sf_built: False`` in its own census, because building a radial FA-anchored bundle under the
name ``sf_arc`` would repeat, one level down, exactly the conflation decision 8 was about.  **Open for
the PI: does population ``sf_arc`` owe a dorsal builder too, or is that a fourth population?**

THE COUNT, AND WHY IT IS STILL A GAP.  :func:`arc_count_from_footprint` derives *how many arcs the
lamellar annulus admits* — ``floor((r_outer - r_min) / pitch) + 1``, the same shape as
``world.geometry.CellFootprint.n_stress_fibres`` — so the count is never typed at a call site.
⚠ That is a DERIVATION, not a source: it is only as good as the arc PITCH, and **there is no arc
spacing in the Contract-Graph.**  Searched 2026-08-21: the vault's four arc-bearing rows
(`KB-DRAFT-7-07`, Tojkander2011, Tojkander2012, Vallenius2013, Burnette2011) give the assembly PATHWAY
and no spacing, no count and no arc-per-cell figure.  So the derived count carries the pitch's OWN
source class — ``PI_GAP`` today — and becomes ``DERIVED`` by itself on the day a pitch is registered.
Filaments per arc is the same PI-GAP band 7-30 that ``stress_fiber.py`` records.  ⚠ 2026-08-20's
lesson is in force: **50 fibres at a typed 0.5 µm pitch put 95.5% of that population outside the
cell.**  Nothing here is typed at a call site and every placement comes from
:func:`aleph.world.geometry.footprint`.

⚠ **THE ARC LENGTH INHERITS A NUMBER THE PI HAS ALREADY RULED IS UNDER BAND.**  An arc's chord is the
LEADING-EDGE width, taken from the footprint's ``lamellipodium_width_um``.  PI decision 11 rules the
sourced leading edge is 5-20 µm (`KB-3.7`) and what is built is **2.5-10x under** it.  The arcs
this module builds are therefore short by the same factor, and they get longer the moment the
footprint does.  Recorded rather than corrected: the leading edge is not this file's to set.

NO STATION AND NO POLARITY, BOTH DELIBERATE.  ``n_stations_placed`` is 0 for the reason
``stress_fiber.py`` gives — ``STATE.md`` (e) 5, *"do not force it"*.  And no ``strand_polarity`` is
written: PI decision 4 makes per-strand polarity a ``StrandPopulation`` property owned by A2/G, and
``active/contraction.py`` already asserts that ``sf_arc`` has no polarity *"by decision"*, refusing an
anti-parallel witness on that ground.  Writing one here would answer a ratified open question with a
build detail.

TOPOLOGY IS ARITHMETIC.  One claim per population per kind; filament ``(a, j)`` of arc ``a`` is strand
``a*n_filaments + j``.  See :mod:`aleph.world.build.microtubule` for the argument.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``pitch_um``, ``lift_um``, ``seg_um``, ``spacing_um`` and every radius [µm]; the
    half-angle and ``dphi`` [rad]; the counts dimensionless.  No stiffness, no motor, no force.
  * boundary — a filament of fewer than 3 beads is refused; either count resolving to zero is refused;
    a ``pitch_um`` narrower than the bundle cross-section is refused as interpenetration; an innermost
    radius that reaches the cell axis is refused; a leading edge wider than the arc's own diameter is
    refused rather than clamped into ``asin``; a ``lift_um`` that carries the arc past the cell's
    mid-plane is refused, since a lamellar arc is not an equatorial one; and the worst-case NODE — the
    outermost filament of the outermost arc, at the top of its hex cross-section — is checked against
    the built membrane radius BEFORE anything is claimed.
  * conservation/invariant — the kernel writes only ``lo + t`` for ``t`` in ``[0, n_nodes)`` and the
    launch dimension IS ``n_nodes``.  ``assert_partitioned`` is asserted after the build.
  * CFL/precision — no integration; float64 throughout.  Segment length varies BY ARC — the angular
    step is shared, so ``seg = r * dphi`` and the outer arc is the longest.  The plan sizes the bead
    count from the OUTER arc so every realised segment is at or below the requested ``seg_um``, and
    both extremes are reported: a bending coefficient ``α = κ/seg³`` must read the per-arc value.
  * sign sense — arcs are concentric about the footprint's disc centre and centred on the +x
    protrusion axis, so the population is symmetric about the protrusion axis rather than growing off
    to one side; ``a = 0`` is the OUTERMOST (leading) arc and increasing ``a`` moves rearward, which is
    the direction of retrograde flow.
  * measurement protocol — one kernel launch of ``dim=n_nodes``; nothing is read back.
    :func:`host_arc_positions` is a numpy twin of the kernel, used ONLY so the PI membrane rule can be
    asserted on a machine with no card — the same role ``active/contraction._twin`` plays.

engine units: length µm, angle rad.  Runtime: NVIDIA Warp on CUDA; the build refuses a device-less arena.
"""

from __future__ import annotations

import math

import numpy as np
import warp as wp

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount, SourceClass
from aleph.world.geometry import CellFootprint
from aleph.world.strand import _require_length

__all__ = [
    "ARC_PITCH_GAP", "arc_count_from_footprint", "plan_sf_arcs", "build_sf_arcs",
    "host_arc_positions", "lay_transverse_arcs_kernel",
]

#: ``√3/2`` — the row pitch of a hexagonal lattice at unit spacing.  Same lattice as ``stress_fiber``.
HEX_ROW_PITCH = wp.constant(wp.float64(0.8660254037844386))

#: Why the arc spacing is a gap, in the shape a ``BondCount`` provenance takes.  Quoted by callers so
#: the PI queue reads the same sentence the code refuses on.
ARC_PITCH_GAP = (
    "transverse-arc lamellar SPACING is ABSENT_FROM_CONTRACT_GRAPH. Searched 2026-08-21: KB-DRAFT-7-07 "
    "and its four arc sources (Tojkander2011_CurrBiol, Tojkander2012_JCellSci, Vallenius2013_OpenBiol, "
    "Burnette2011_NatCellBiol) carry the assembly PATHWAY (Tm4/myosin-II condensation + Arp2/3 "
    "annealing) and no spacing, no per-cell count and no arc width. Registering one is a ModelContract "
    "change and is PI-authored."
)


@wp.kernel
def lay_transverse_arcs_kernel(
    lo: wp.int32, n_per: wp.int32, n_fil: wp.int32, hex_cols: wp.int32, hex_rows: wp.int32,
    strand_base: wp.int32, range_ordinal: wp.int32,
    cx: wp.float64, cy: wp.float64, z_arc: wp.float64,
    r_outer: wp.float64, pitch_um: wp.float64, spacing_um: wp.float64,
    phi0: wp.float64, dphi: wp.float64,
    position: wp.array(dtype=wp.vec3d), strand_id: wp.array(dtype=wp.int32),
    range_id: wp.array(dtype=wp.int32),
):
    """Write one arc bead per thread: arc index, hex offset and angular position, all from ``tid``.

    The bundle cross-section is carried in the LOCAL ``(radial, z)`` frame of each bead, which is what
    makes this kernel different from the ventral one: there the frame is shared by the whole
    population, here it turns with ``phi``.  ``a = 0`` is the outermost arc.
    """
    t = wp.tid()
    per_arc = n_fil * n_per
    a = t / per_arc
    rem = t - a * per_arc
    j = rem / n_per
    i = rem - j * n_per

    row = j / hex_cols
    col = j - row * hex_cols
    stagger = wp.float64(row - wp.int32(2) * (row / wp.int32(2))) * wp.float64(0.5)
    hex_r = (wp.float64(col) - wp.float64(0.5) * wp.float64(hex_cols - 1) + stagger) * spacing_um
    hex_z = (wp.float64(row) - wp.float64(0.5) * wp.float64(hex_rows - 1)) * spacing_um * HEX_ROW_PITCH

    r = r_outer - wp.float64(a) * pitch_um + hex_r
    phi = phi0 + wp.float64(i) * dphi

    g = lo + t
    position[g] = wp.vec3d(cx + r * wp.cos(phi), cy + r * wp.sin(phi), z_arc + hex_z)
    strand_id[g] = strand_base + a * n_fil + j
    range_id[g] = range_ordinal


def _envelope(fp: CellFootprint, *, lift_um: float, spacing_um: float, n_fil: int) -> dict[str, float]:
    """The arc placement envelope, DERIVED from the built cell — never from a call-site constant.

    Args:
        fp: the footprint of the cell that was actually built.
        lift_um: how far above the basal plane the arc plane sits [µm].
        spacing_um: filament centre-to-centre spacing inside one bundle [µm].
        n_fil: filaments per arc, already resolved.

    Returns:
        ``z_arc``, the outermost admissible arc radius, the bundle's own half-extents, and the
        half-angle its chord subtends.

    Raises:
        ValueError: if the lift carries the arc plane to or past the cell's mid-plane, or if the
            leading edge is wider than the arc that has to span it.
    """
    if not (0.0 < lift_um < abs(fp.z_basal_um)):
        raise ValueError(
            f"lift_um={lift_um} does not place an arc in the LAMELLA: the basal plane is at "
            f"z={fp.z_basal_um} µm, so a positive lift below {abs(fp.z_basal_um)} µm is the band. A "
            "lift at or past the cell's mid-plane is an equatorial ring, which is a different "
            "structure with a different name — the error PI decision 8 exists to prevent."
        )
    hex_cols = int(math.ceil(math.sqrt(n_fil)))
    hex_rows = int(math.ceil(n_fil / hex_cols))
    half_w = 0.5 * hex_cols * spacing_um
    half_h = 0.5 * hex_rows * spacing_um * 0.8660254037844386
    z_arc = fp.z_basal_um + lift_um
    # The cell is WIDER at the arc plane than at the substrate — that is the whole point of lifting it.
    # `clearance` is the footprint's own margin and is reused rather than re-invented; the bundle's own
    # half-width is then subtracted, because the envelope places a CENTRELINE and the kernel writes a
    # cross-section around it.
    r_span = math.sqrt(fp.r_cell_um ** 2 - z_arc ** 2)
    r_outer = fp.clearance * r_span - half_w
    half_chord = 0.5 * fp.lamellipodium_width_um
    if r_outer <= 0.0 or half_chord >= r_outer:
        raise ValueError(
            f"the leading edge is {fp.lamellipodium_width_um} µm wide and the outermost arc admitted at "
            f"z={z_arc:.4g} µm has radius {r_outer:.4g} µm, so no arc concentric with the rim can span "
            "it. An arc is a chord of the cell, not a straight bundle laid across it."
        )
    return {
        "z_arc_um": z_arc, "r_outer_um": r_outer, "r_span_um": r_span,
        "bundle_half_width_um": half_w, "bundle_half_height_um": half_h,
        "hex_cols": float(hex_cols), "hex_rows": float(hex_rows),
        "half_angle_rad": math.asin(half_chord / r_outer),
        "chord_um": fp.lamellipodium_width_um,
    }


def arc_count_from_footprint(
    fp: CellFootprint,
    *,
    pitch_um: float,
    spacing_um: float,
    lift_um: float,
    n_filaments_per_arc: int,
    pitch_scope: str,
    pitch_provenance: str,
    pitch_source_class: SourceClass | str,
) -> BondCount:
    """How many arcs the lamellar annulus ADMITS at ``pitch_um``, as a :class:`BondCount`.

    ``floor((r_outer - r_min) / pitch) + 1``, the same shape as ``CellFootprint.n_stress_fibres``:
    the pitch is held and the count yields, which is the direction that keeps a sourced quantity
    sourced.  ``r_min`` is the innermost radius a bundle can occupy without reaching the cell axis.

    ⚠ **The result is only as good as the pitch**, so the pitch's own source class is carried through:
    a ``PI_GAP`` pitch gives a ``PI_GAP`` count, and the same call returns ``DERIVED`` unchanged on the
    day a spacing is registered.  A derivation does not launder its inputs.

    Args:
        fp: the footprint of the built cell.
        pitch_um: centre-to-centre radial spacing between neighbouring arcs [µm].
        spacing_um: filament spacing inside one bundle [µm].
        lift_um: arc-plane height above the basal plane [µm].
        n_filaments_per_arc: filaments per arc, already resolved.
        pitch_scope: the cell, state and assay the PITCH is true for.
        pitch_provenance: the citation or the derivation for the pitch — :data:`ARC_PITCH_GAP` today.
        pitch_source_class: where the pitch came from.

    Returns:
        The arc count, with the derivation and the pitch's provenance travelling with it.
    """
    pitch = _require_length("pitch_um", pitch_um)
    env = _envelope(fp, lift_um=lift_um, spacing_um=spacing_um, n_fil=n_filaments_per_arc)
    span = env["r_outer_um"] - env["bundle_half_width_um"]
    n = int(span / pitch) + 1
    cls = SourceClass(pitch_source_class)
    # A derivation from a sourced input is DERIVED; a derivation from a blank is still that blank.
    out = SourceClass.DERIVED if cls in (SourceClass.SOURCED, SourceClass.DERIVED) else cls
    return BondCount(
        basis="explicit", value=float(n), scope=pitch_scope,
        source_class=out,
        provenance=(
            f"DERIVED: floor((r_outer - r_min)/pitch) + 1 = floor(({env['r_outer_um']:.4f} - "
            f"{env['bundle_half_width_um']:.4f})/{pitch}) + 1 = {n}, over the lamellar annulus of the "
            f"BUILT cell (R={fp.r_cell_um} µm, arc plane z={env['z_arc_um']:.4f} µm). The derivation "
            f"inherits the PITCH's source class ({cls}), which is: {pitch_provenance}"
        ),
    )


def plan_sf_arcs(
    fp: CellFootprint,
    *,
    count: BondCount,
    bundle_count: BondCount,
    support: float = 1.0,
    pitch_um: float | None = None,
    lift_um: float | None = None,
    seg_um: float | None = None,
    spacing_um: float | None = None,
) -> dict[str, object]:
    """Resolve the transverse-arc inventory against the built cell, no device touched.

    Args:
        fp: the footprint of the cell that was built.  Every placement comes from here.
        count: how many transverse arcs this cell has.  **REQUIRED** — see
            :func:`arc_count_from_footprint`, and never a literal at the call site.
        bundle_count: filaments per arc.  **REQUIRED**.  Band 7-30 and PI-GATED, the same gap
            ``stress_fiber.py`` records for the ventral bundle.
        support: the MEASURED support ``count`` resolves against; 1.0 for an ``explicit`` count.
        pitch_um: centre-to-centre radial spacing between neighbouring arcs [µm].  **REQUIRED**.
        lift_um: how far the arc plane sits above the basal plane [µm].  **REQUIRED**.
        seg_um: bead spacing along the arc [µm].  **REQUIRED**; it is an upper bound here, since the
            realised segment is ``r·dphi`` and shrinks on the inner arcs.
        spacing_um: filament centre-to-centre spacing inside one bundle [µm].  **REQUIRED**.

    Returns:
        The inventory, with both provenance rows, the realised segment extremes and the envelope it
        was derived from.

    Raises:
        ValueError: on a missing or non-positive length, either count resolving to zero, a filament of
            fewer than 3 beads, a pitch narrower than the bundle cross-section, an innermost arc that
            reaches the cell axis, or a worst-case node outside the built membrane.
    """
    pitch = _require_length("pitch_um", pitch_um)
    lift = _require_length("lift_um", lift_um)
    seg = _require_length("seg_um", seg_um)
    spacing = _require_length("spacing_um", spacing_um)

    n_arcs = count.resolve(support)
    n_fil = bundle_count.resolve(1.0)
    for label, value in (("transverse arcs", n_arcs), ("filaments per arc", n_fil)):
        if value < 1:
            raise ValueError(
                f"the count resolves to {value} {label}. A population that is 'present but has nothing' "
                "is indistinguishable downstream from one that failed to build."
            )

    env = _envelope(fp, lift_um=lift, spacing_um=spacing, n_fil=n_fil)
    r_outer, half_w, half_h = env["r_outer_um"], env["bundle_half_width_um"], env["bundle_half_height_um"]

    if n_arcs > 1 and pitch <= 2.0 * half_w:
        raise ValueError(
            f"pitch_um={pitch} is not wider than the bundle's own cross-section {2.0 * half_w} µm "
            f"({int(env['hex_cols'])} columns at spacing_um={spacing}), so neighbouring arcs "
            "interpenetrate. The packing and the spacing come from different sources — nothing else in "
            "this build would notice them crossing, which is why it is asserted here."
        )
    r_inner = r_outer - (n_arcs - 1) * pitch
    if r_inner <= half_w:
        raise ValueError(
            f"{n_arcs} arcs at pitch {pitch} µm put the innermost centreline at r={r_inner:.4g} µm, "
            f"inside its own bundle half-width {half_w:.4g} µm — the arc has collapsed onto the cell "
            "axis and is no longer concentric with anything. Fewer arcs or a tighter pitch, but not a "
            "clamp: an arc folded through the axis renders as a rosette and reads as a build."
        )

    arc_len_outer = 2.0 * env["half_angle_rad"] * r_outer
    n_seg = int(math.ceil(arc_len_outer / seg))
    n_per = n_seg + 1
    if n_per < 3:
        raise ValueError(
            f"the outer arc is {arc_len_outer:.4g} µm long at seg_um={seg}, giving {n_per} beads per "
            "filament; a filament needs at least 3 so it carries one bending triple. A filament with "
            "no bending stiffness is a different physical object, not a coarser one."
        )
    dphi = 2.0 * env["half_angle_rad"] / n_seg

    # THE PI RULE, 2026-08-20 — checked on the WORST NODE before anything is claimed. The outermost
    # filament of the outermost arc, at the top of its hex cross-section, is the only candidate.
    worst_r = math.hypot(r_outer + half_w, abs(env["z_arc_um"]) + half_h)
    if worst_r >= fp.r_cell_um:
        raise ValueError(
            f"the outermost arc's outermost filament reaches r={worst_r:.4g} µm against a built "
            f"membrane of {fp.r_cell_um} µm. Every node lies inside the built membrane (PI 2026-08-20); "
            "`assert_partitioned` passing does NOT cover this, because partition is a property of ID "
            "ranges and says nothing about where a node is."
        )

    return {
        "kind": "sf_arc",
        "structure": "transverse arc — FA-FREE, curved, lamellar. NOT a ventral FA-to-FA bundle "
                     "(build/stress_fiber.py) and NOT a dorsal FA-anchored fibre. PI decision 8, "
                     "2026-08-21.",
        "dorsal_sf_built": False,
        "n_strands": n_arcs * n_fil,
        "n_arcs": n_arcs,
        "filaments_per_arc": n_fil,
        "nodes_per_strand": n_per,
        "hex_cols": int(env["hex_cols"]), "hex_rows": int(env["hex_rows"]),
        "bundle_width_um": 2.0 * half_w, "bundle_height_um": 2.0 * half_h,
        "pitch_um": pitch, "spacing_um": spacing, "lift_um": lift,
        "z_arc_um": env["z_arc_um"],
        "r_outer_um": r_outer, "r_inner_um": r_inner,
        "half_angle_rad": env["half_angle_rad"], "dphi_rad": dphi, "phi0_rad": -env["half_angle_rad"],
        "chord_um": env["chord_um"],
        "chord_note": "the arc chord IS the leading-edge width from the footprint. PI decision 11 "
                      "rules the sourced leading edge is 5-20 um (KB-3.7) and the built one is "
                      "2.5-10x under it, so these arcs are short by the same factor.",
        "seg_um_requested": seg,
        "seg_um_realised_outer": r_outer * dphi,
        "seg_um_realised_inner": r_inner * dphi,
        "arc_length_um_outer": arc_len_outer,
        "worst_node_radius_um": worst_r,
        "n_stations_placed": 0,
        "polarity_note": "ABSENT, and it is a GAP rather than a default. A transverse arc is "
                         "CONTRACTILE (Tm4/myosin-II condensation, KB-DRAFT-7-07), so an arc "
                         "population has to be able to hold an anti-parallel pair before NMII binds "
                         "one. The per-strand machinery landed 2026-08-21 (832e4313, PI decision 4) on "
                         "build/cortex.py's StrandPopulation, which this builder does not use; wiring "
                         "it is an architectural claim about arc rectification and owes a "
                         "polarity_basis. Surfaced, not guessed.",
        "provenance": {"n_transverse_arcs": _row(count), "filaments_per_arc": _row(bundle_count)},
    }


def _row(count: BondCount) -> dict[str, object]:
    """The count with its scope and authority attached, in the shape ``BondFamily.provenance_row`` uses."""
    return {
        "basis": count.basis, "value": count.value, "support_unit": count.support_unit,
        "scope": count.scope, "source_class": str(count.source_class), "provenance": count.provenance,
    }


def host_arc_positions(inventory: dict[str, object], centre_xy: tuple[float, float]) -> np.ndarray:
    """Numpy twin of :func:`lay_transverse_arcs_kernel`, for the membrane rule on a card-less machine.

    It is not a second implementation of a law — there is no law here, only placement — and it exists
    for the reason ``active/contraction._twin`` exists: the PI rule *every node lies inside the built
    membrane* is worth asserting on the dev machine, and the alternative is asserting it for the first
    time on the GPU host.

    Args:
        inventory: what :func:`plan_sf_arcs` returned.
        centre_xy: the footprint disc centre the arcs are concentric with [µm].

    Returns:
        ``(n_nodes, 3)`` positions in the same order the kernel writes them.
    """
    n_arcs, n_fil = int(inventory["n_arcs"]), int(inventory["filaments_per_arc"])
    n_per, cols, rows = (int(inventory[k]) for k in ("nodes_per_strand", "hex_cols", "hex_rows"))
    spacing = float(inventory["spacing_um"])

    t = np.arange(n_arcs * n_fil * n_per)
    a, rem = np.divmod(t, n_fil * n_per)
    j, i = np.divmod(rem, n_per)
    row, col = np.divmod(j, cols)
    hex_r = (col - 0.5 * (cols - 1) + 0.5 * (row % 2)) * spacing
    hex_z = (row - 0.5 * (rows - 1)) * spacing * 0.8660254037844386

    r = float(inventory["r_outer_um"]) - a * float(inventory["pitch_um"]) + hex_r
    phi = float(inventory["phi0_rad"]) + i * float(inventory["dphi_rad"])
    return np.stack([centre_xy[0] + r * np.cos(phi), centre_xy[1] + r * np.sin(phi),
                     float(inventory["z_arc_um"]) + hex_z], axis=1)


def build_sf_arcs(
    arena: WorldArena,
    fp: CellFootprint,
    *,
    count: BondCount,
    bundle_count: BondCount,
    support: float = 1.0,
    centre_xy: tuple[float, float] = (0.0, 0.0),
    pitch_um: float | None = None,
    lift_um: float | None = None,
    seg_um: float | None = None,
    spacing_um: float | None = None,
    population: str = "sf_arc",
) -> dict[str, object]:
    """Claim the transverse-arc ranges and lay every arc down with one CUDA kernel.

    Args:
        arena: the world to claim from.  Must hold a CUDA allocation.
        fp: the footprint of the built cell; every placement is derived from it.
        centre_xy: the basal disc centre the arcs are concentric with [µm].  The footprint's own
            origins are built about ``(0, 0)``, which is the default.
        count / bundle_count / support / pitch_um / lift_um / seg_um / spacing_um: see
            :func:`plan_sf_arcs`.
        population: the name the four claims are attributed to.

    Returns:
        The census of what was built.

    Raises:
        RuntimeError: if ``arena`` holds no device allocation.
        ValueError: on any refusal from :func:`plan_sf_arcs`.
    """
    inventory = plan_sf_arcs(fp, count=count, bundle_count=bundle_count, support=support,
                             pitch_um=pitch_um, lift_um=lift_um, seg_um=seg_um, spacing_um=spacing_um)

    n_strand = int(inventory["n_strands"])
    n_per = int(inventory["nodes_per_strand"])
    n_nodes = n_strand * n_per
    if not arena.node_arrays:
        raise RuntimeError(
            f"{population}: this arena holds no device allocation, and construction here is a CUDA "
            "kernel rather than a host build that is uploaded. Warp CUDA is the only runtime; build the "
            f"arena with a CUDA device. The inventory it would have built: {inventory['n_arcs']} arcs "
            f"x {inventory['filaments_per_arc']} filaments x {n_per} beads = {n_nodes} nodes."
        )

    nodes = arena.claim(population, Kind.NODE, n_nodes)
    strands = arena.claim(population, Kind.STRAND, n_strand)
    segments = arena.claim(population, Kind.SEGMENT, n_strand * (n_per - 1))
    angles = arena.claim(population, Kind.ANGLE3, n_strand * (n_per - 2))
    ordinal = len(arena.claims()) - 4

    wp.launch(
        lay_transverse_arcs_kernel, dim=n_nodes,
        inputs=[nodes.lo, n_per, int(inventory["filaments_per_arc"]), int(inventory["hex_cols"]),
                int(inventory["hex_rows"]), strands.lo, ordinal,
                float(centre_xy[0]), float(centre_xy[1]), float(inventory["z_arc_um"]),
                float(inventory["r_outer_um"]), float(inventory["pitch_um"]),
                float(inventory["spacing_um"]), float(inventory["phi0_rad"]),
                float(inventory["dphi_rad"])],
        outputs=[arena.node_arrays["position"], arena.node_arrays["strand_id"],
                 arena.node_arrays["range_id"]],
        device=arena.device,
    )
    arena.assert_partitioned()

    return inventory | {
        "population": population, "n_nodes": nodes.count, "n_segments": segments.count,
        "n_angles": angles.count, "centre_xy_um": (float(centre_xy[0]), float(centre_xy[1])),
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


#: The two axes a caller must declare, quoted by ``_demo`` and by the tests so the gap has one wording.
_PITCH_SCOPE = "MCF7 epithelial on glass — NO VALUE EXISTS for transverse-arc lamellar spacing"


def _demo() -> None:
    """Self-check: the derivation, the PI membrane rule on the built nodes, and every refusal."""
    from aleph.world.geometry import footprint

    fp = footprint(7.5, -7.0)
    per_arc = BondCount(
        basis="explicit", value=20.0, scope="transverse arc bundle, generic — band 7-30",
        source_class=SourceClass.PI_GAP,
        provenance="laws/architecture_spec.STRESS_FIBER notes 'N_filaments PI-GATED (no KB point "
                   "value)'; no arc-specific bundle count exists either")
    axes = {"pitch_um": 0.6, "lift_um": 0.4, "seg_um": 0.05, "spacing_um": 0.012}

    # 1. THE COUNT IS DERIVED FROM THE BUILT CELL, and it inherits the pitch's class rather than
    #    laundering it. Nothing is typed at a call site — that is the 2026-08-20 defect.
    count = arc_count_from_footprint(
        fp, pitch_um=axes["pitch_um"], spacing_um=axes["spacing_um"], lift_um=axes["lift_um"],
        n_filaments_per_arc=20, pitch_scope=_PITCH_SCOPE, pitch_provenance=ARC_PITCH_GAP,
        pitch_source_class=SourceClass.PI_GAP)
    assert count.source_class is SourceClass.PI_GAP, "a PI_GAP pitch may not yield a DERIVED count"
    assert "ABSENT_FROM_CONTRACT_GRAPH" in count.provenance and "floor(" in count.provenance
    assert count.value >= 1.0

    inv = plan_sf_arcs(fp, count=count, bundle_count=per_arc, **axes)
    assert inv["kind"] == "sf_arc" and inv["dorsal_sf_built"] is False
    assert inv["n_stations_placed"] == 0, "no station is placed; STATE.md (e) 5 says do not force it"
    assert "strand_polarity" not in inv, "polarity is PI decision 4's, not this builder's"
    assert inv["n_arcs"] == int(count.value) and inv["n_strands"] == inv["n_arcs"] * 20
    # The angular step is shared, so the realised segment SHRINKS inward and never exceeds the request.
    assert inv["seg_um_realised_inner"] < inv["seg_um_realised_outer"] <= axes["seg_um"] + 1e-12

    # 2. THE PI RULE, ON THE ACTUAL NODES — every one inside the built membrane, checked by the
    #    checker that owns the rule rather than by a repeat of the arithmetic that placed them.
    from aleph.world.geometry import assert_inside_membrane

    pos = host_arc_positions(inv, (0.0, 0.0))
    assert pos.shape == (inv["n_strands"] * int(inv["nodes_per_strand"]), 3)

    class _P:
        def __init__(self, lo: int, hi: int) -> None:
            self.nodes = type("R", (), {"lo": lo, "hi": hi})()

    rep = assert_inside_membrane(pos, {"sf_arc": _P(0, pos.shape[0])}, fp.r_cell_um)
    assert rep["sf_arc"]["n_outside"] == 0
    assert rep["sf_arc"]["r_max_um"] <= inv["worst_node_radius_um"] + 1e-9

    # The arcs are lifted OFF the substrate and curved — the two things that make this not a ventral
    # fibre are properties of the built nodes, not of the docstring.
    assert pos[:, 2].min() > fp.z_basal_um, "an arc is lifted OFF the substrate; a ventral fibre is on it"
    assert abs(float(np.mean(pos[:, 2])) - inv["z_arc_um"]) < 1e-9
    outer = pos[:int(inv["nodes_per_strand"])]
    assert np.ptp(np.hypot(outer[:, 0], outer[:, 1])) < 1e-9, "an arc is concentric, so r is constant"
    assert np.ptp(outer[:, 1]) > 1e-3, "and it is CURVED — a straight bundle would be a ventral fibre"

    # 3. Every refusal.
    for missing in axes:
        try:
            plan_sf_arcs(fp, count=count, bundle_count=per_arc,
                         **{k: v for k, v in axes.items() if k != missing})
        except ValueError as exc:
            assert "no default" in str(exc) and missing in str(exc), str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"a missing {missing} must refuse")

    for bad, frag in (
        ({"pitch_um": 0.05}, "interpenetrate"),          # inside the bundle's own cross-section
        ({"lift_um": 7.0}, "equatorial ring"),           # at the mid-plane, a different structure
        ({"seg_um": 100.0}, "at least 3"),               # a filament with no bending triple
    ):
        try:
            plan_sf_arcs(fp, count=count, bundle_count=per_arc, **(axes | bad))
        except ValueError as exc:
            assert frag in str(exc), (bad, exc)
        else:  # pragma: no cover
            raise AssertionError(f"{bad} must refuse")

    # Too many arcs at this pitch fold the innermost one through the cell axis.
    many = BondCount(basis="explicit", value=count.value * 4, scope=_PITCH_SCOPE,
                     source_class=SourceClass.PI_GAP, provenance="deliberately over-subscribed")
    try:
        plan_sf_arcs(fp, count=many, bundle_count=per_arc, **axes)
    except ValueError as exc:
        assert "collapsed onto the cell axis" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an over-subscribed annulus must refuse")

    # A device-less arena refuses, and reports the inventory it would have built.
    arena = WorldArena(capacity={Kind.NODE: 10, Kind.STRAND: 10, Kind.SEGMENT: 10, Kind.ANGLE3: 10})
    try:
        build_sf_arcs(arena, fp, count=count, bundle_count=per_arc, **axes)
    except RuntimeError as exc:
        assert "no device allocation" in str(exc) and f"{inv['n_arcs']} arcs" in str(exc), str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device-less build must refuse")
    assert arena.n_live(Kind.NODE) == 0, "a refused build claims nothing"

    chars = _codegen_check()
    assert chars > 0

    print(f"sf_arc self-check OK — CUDA codegen {chars} chars; {inv['n_arcs']} transverse arcs x "
          f"{inv['filaments_per_arc']} filaments x {inv['nodes_per_strand']} beads = "
          f"{inv['n_strands'] * int(inv['nodes_per_strand'])} nodes; arc plane z={inv['z_arc_um']:.3f} "
          f"µm, r {inv['r_inner_um']:.3f}-{inv['r_outer_um']:.3f} µm, chord {inv['chord_um']:.3f} µm "
          f"(PI decision 11: 2.5-10x under band), worst node r={inv['worst_node_radius_um']:.3f} µm")


if __name__ == "__main__":
    _demo()
