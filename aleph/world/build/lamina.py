r"""The lamin meshwork: the nucleus's load-bearing shell, at ITS OWN spacing — which is sourced.

**THE SPACING THE ENVELOPE COULD NOT FIND IS IN THE CORPUS.**  :mod:`aleph.world.build.envelope`
carries :data:`~aleph.world.build.envelope.ENVELOPE_MESH_PI_GAP`, whose sentence is that "no lamin
meshwork spacing is registered in this repository", so it borrowed the cortical 50-100 nm band and
labelled the borrowing an analogy.  That sentence is now false.  ``Stephens2017_MBoC`` — the paper the
KB's own chromatin and lamina claims already cite — carries a parameter table on p7 giving three
lamina numbers together, attributed to Shimi et al. 2015:

======================================  ==========  ================
quantity                                simulation  experimental
======================================  ==========  ================
nodes in lamina, ``N``                        1 000  1 000 - 3 000
node connectivity, ``z``                      ~4.5   4
lamin filament length ``(4πR²/N)^½`` [µm]      0.6    0.4
======================================  ==========  ================

**The lamin meshwork is 0.4 µm — four to eight times COARSER than the cortex, not equal to it.**  So
the analogy does not become a derivation; it inverts.  Applying this repository's own rule (a shell's
resolution follows the physical structure that holds it) to the nuclear envelope gives icosphere level
**4**, 2 562 vertices, mean edge 0.384 µm — where the envelope currently stands at level 6 and 40 962.
Two independent columns of that one table agree at the nuclear radius: the count derived from the
spacing (2 043) lands inside the independently reported node band (1 000-3 000), and the spacing that
count realises (0.400 µm) lands on the reported filament length.  That agreement is what makes it a
derivation.  Reconciling it with the envelope's level is a PI decision and is NOT taken here.

⚠ **THE CITATION IS SECOND-HAND AND IS LABELLED AS SUCH.**  Shimi et al. 2015 is in neither the PDF
corpus nor ``source_evidence``; what is in the corpus is Stephens' READING of it.  That is a real
grade below a sourced value and :data:`LAMIN_MESHWORK_SOURCE` says so in the sentence it carries into
the artifact.  It is nonetheless the only one of the three declared envelope resolutions with a
physical source behind it — the other two are the KB gap card's ``subdivisions = 3`` and the
envelope's own analogy.

**THIS IS A DISJOINT POPULATION, NOT A RELABELLED ENVELOPE.**  The charter's rule is that components
own state and couple only through connectors, "never shared nodes and never a permanent weld", and the
lamina and the envelope are two different networks — a lamin intermediate-filament meshwork and a
lipid bilayer surface — at two different spacings.  Sharing the envelope's vertices would be
co-location, which "is NEVER a connection".  So this module claims its own NODE range and the LINC and
lamina-envelope couplings are BONDS, later.

**NODES ONLY, AND THAT IS THE POINT.**  A closed triangulation would hand the meshwork connectivity
``z`` = 6, against the sourced 4 — 1.5x the load-bearing filaments at the same node count, which is a
stiffness error dressed as a mesh choice.  The topology is a BOND question and bonds are PHASE 3, so
laying one down now would bake the wrong answer into the arena layout before anyone chose it.  What
this module builds is the node population at the sourced areal spacing; :func:`plan_lamina` reports
``n_filaments_implied`` at the sourced ``z`` so the PHASE 3 claim can be sized, and reports nothing as
built that is not built.

WHAT IS DELIBERATELY ABSENT.  No areal tension, no bending, no rupture criterion, no LINC tether, no
chromatin tether, and no lamin-A/C versus lamin-B split — those are laws carrying parameters that
``components/incumbent/compartments.py`` itself marks "I0-B2 PI GAP, TEST values".  Also absent: any
standoff between the lamina and the inner nuclear membrane.  The lamina is built AT the nuclear
radius, and that is not a physiological choice — the reported lamina thickness is tens of nanometres,
an order of magnitude below this population's own 0.4 µm spacing, so the two surfaces are the same
surface at this resolution and a nonzero offset would be a number invented under the resolution that
could resolve it.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``radius_um``, ``filament_um`` and the realised spacing [µm]; ``N`` and ``z`` are
    dimensionless; the areal spacing is ``(area / count)^½`` and so carries µm, which is the source's
    own definition rather than a re-derivation.
  * boundary — radius and filament length are refused without a value and refused non-positive; a
    count derived outside the independently sourced band :data:`LAMIN_NODE_BAND` is REFUSED rather
    than built, which is what rejects a caller who reaches for the cortical spacing by habit (0.1 µm
    at the nuclear radius overshoots the band by ~11x); an empty claim is refused by the arena.
  * conservation/invariant — one contiguous NODE claim; every global index the kernel forms is
    ``lo + t`` for ``t`` in ``[0, N)`` and the launch dimension IS ``N``, so out-of-range addressing is
    unreachable rather than checked.  ``assert_partitioned`` is asserted after the build.
  * CFL/precision — no integration; float64 throughout, matching the arena's ``vec3d`` arrays.  The
    golden-angle azimuth is evaluated in float64 so two nodes 2 000 apart remain distinguishable.
  * sign sense — the only signed quantity is the radius, asserted positive; nodes are placed on the
    sphere of that radius about ``centre_um`` and nothing is accumulated.
  * measurement protocol — one kernel of ``dim=N``; nothing is read back.  The realised spacing is the
    source's areal formula, NOT a measured nearest-neighbour distance, and is reported under a name
    that says so.

engine units: length µm.  Runtime: NVIDIA Warp on CUDA; the build refuses a device-less arena.
"""

from __future__ import annotations

import math

import warp as wp

from aleph.world.arena import Kind, WorldArena
# The same golden-angle azimuth the aster uses. ONE definition of the spiral in this package: two
# populations placed by two spellings of the same construction is a difference nobody would find.
from aleph.world.build.microtubule import GOLDEN_AZIMUTH
from aleph.world.strand import _require_length

__all__ = [
    "LAMIN_CONNECTIVITY",
    "LAMIN_FILAMENT_UM",
    "LAMIN_MESHWORK_SOURCE",
    "LAMIN_NODE_BAND",
    "build_lamina",
    "lay_shell_kernel",
    "plan_lamina",
]

#: Lamin filament length [µm] — the meshwork edge, i.e. the areal spacing ``(4πR²/N)^½``.  The
#: EXPERIMENTAL column; the simulation column of the same table is 0.6 µm, which is that paper's
#: rounding of 0.56 at its own ``N`` = 1 000 and ``R`` = 5 µm.
LAMIN_FILAMENT_UM = 0.4

#: Independently reported node count of the lamin meshwork, from the same table.  Used as a REFUSAL
#: band on the count derived from :data:`LAMIN_FILAMENT_UM`, not as a source for the count itself —
#: the point of two numbers from one measurement is that either can check the other.
LAMIN_NODE_BAND = (1000, 3000)

#: Reported node connectivity of the meshwork.  NOT used to build anything here: it sizes the PHASE 3
#: bond claim and it is the number a closed triangulation would get wrong by 1.5x.
LAMIN_CONNECTIVITY = 4.0

#: The provenance sentence a caller carries, verbatim, into the run record.  It is a sentence rather
#: than a citation key because the GRADE is the part a later reader has to see: this is Stephens'
#: reading of Shimi, and Shimi is not in this repository.
LAMIN_MESHWORK_SOURCE = (
    "SOURCED_SECONDARY: lamin meshwork spacing 0.4 um (experimental column), node count 1000-3000, "
    "connectivity z=4 — Stephens et al. 2017 Mol Biol Cell 28:1984 (10.1091/mbc.E16-09-0653) Table 1, "
    "p7, attributing all three to Shimi et al. 2015. ⚠ Shimi 2015 is in NEITHER the PDF corpus NOR "
    "source_evidence, so this is a SECOND-HAND citation and not a primary one. It nonetheless "
    "supersedes the cortical-mesh analogy in envelope.ENVELOPE_MESH_PI_GAP, which is 4-8x too FINE: "
    "the lamin meshwork is coarser than the cortex, not equal to it."
)


@wp.kernel
def lay_shell_kernel(
    lo: wp.int32, n: wp.int32, range_ordinal: wp.int32,
    cx: wp.float64, cy: wp.float64, cz: wp.float64, radius: wp.float64,
    position: wp.array(dtype=wp.vec3d), strand_id: wp.array(dtype=wp.int32),
    range_id: wp.array(dtype=wp.int32),
):
    """Place one meshwork node per thread on the sphere, from its own index and nothing else.

    The golden-angle spiral spreads ``n`` points at near-uniform areal density, which is the
    distribution the source's own ``(4πR²/N)^½`` presumes when it converts a count into a length.  A
    triangulated sphere would place them just as uniformly and would additionally assert a
    connectivity nobody chose — see the module docstring.
    """
    t = wp.tid()
    td = wp.float64(t) + wp.float64(0.5)
    cosp = wp.float64(1.0) - wp.float64(2.0) * td / wp.float64(n)
    sinp = wp.sqrt(wp.max(wp.float64(0.0), wp.float64(1.0) - cosp * cosp))
    az = GOLDEN_AZIMUTH * td

    g = lo + t
    position[g] = wp.vec3d(
        cx + radius * sinp * wp.cos(az),
        cy + radius * sinp * wp.sin(az),
        cz + radius * cosp,
    )
    strand_id[g] = -1          # a meshwork node belongs to no strand, as the MTOC does not
    range_id[g] = range_ordinal


def plan_lamina(*, radius_um: float | None = None, filament_um: float | None = None) -> dict[str, object]:
    """Derive the meshwork's inventory from its radius and filament length, touching no device.

    The count is ``N = 4πR² / L²`` — the source's own ``L = (4πR²/N)^½`` inverted, so this is that
    statement read the other way and not a second model of the same thing.  The derived ``N`` is then
    checked against :data:`LAMIN_NODE_BAND`, which the SAME table reports independently.

    Args:
        radius_um: the radius the meshwork sits on [µm].  **REQUIRED, no default.**  Derive it with
            :func:`~aleph.world.build.envelope.nuclear_radius_um` so the lamina and the envelope move
            together instead of drifting apart under a swept cell radius.
        filament_um: lamin filament length = meshwork areal spacing [µm].  **REQUIRED, no default.**
            :data:`LAMIN_FILAMENT_UM` is the sourced value; passing a number here rather than reading
            a default is what keeps the spacing a declared axis.

    Returns:
        The inventory: ``n_nodes``, the realised areal spacing, the implied PHASE 3 filament count at
        the sourced connectivity, and the band the count was checked against.

    Raises:
        ValueError: on a missing or non-positive length, or on a count falling outside
            :data:`LAMIN_NODE_BAND` — the two sourced numbers disagreeing is a fact about the inputs,
            and widening the band to accept it would be choosing a constant to make a check pass.
    """
    radius = _require_length("radius_um", radius_um)
    filament = _require_length("filament_um", filament_um)

    area_um2 = 4.0 * math.pi * radius * radius
    n_nodes = int(round(area_um2 / (filament * filament)))
    lo, hi = LAMIN_NODE_BAND
    if not (lo <= n_nodes <= hi):
        raise ValueError(
            f"filament_um={filament} at radius_um={radius} gives N = 4piR^2/L^2 = {n_nodes} lamin "
            f"meshwork nodes, outside the independently reported band {lo}-{hi}. Both numbers come "
            "from ONE table (Stephens 2017 MBoC Table 1 p7, after Shimi 2015) and are supposed to "
            "agree; that they do not means the spacing passed here is not this structure's. The "
            f"cortical mesh band lands here — 0.100 um gives {int(round(area_um2 / 0.01))} nodes — "
            "which is the substitution envelope.ENVELOPE_MESH_PI_GAP makes and this refuses. Do NOT "
            "widen the band: it is a contract written before the run."
        )
    return {
        "kind": "lamin_meshwork",
        "n_nodes": n_nodes,
        "radius_um": radius,
        "shell_area_um2": area_um2,
        "filament_um_requested": filament,
        # The source's areal definition, not a measured nearest-neighbour distance — named to say so.
        "areal_spacing_um_realised": math.sqrt(area_um2 / n_nodes),
        "connectivity_sourced": LAMIN_CONNECTIVITY,
        "n_filaments_implied": int(round(0.5 * n_nodes * LAMIN_CONNECTIVITY)),
        "node_band_sourced": list(LAMIN_NODE_BAND),
        "standoff_um": 0.0,
        # ⚠ MACHINE-READABLE, because the decision above is already argued in this module's docstring
        # and CONSUMERS DO NOT READ DOCSTRINGS. On 2026-08-21 `world_export_cell.py` died with
        # `KeyError: 'n_strands'` the first time the full cell was exported, and the Lead filed the
        # node-only state as an unrecorded defect — it is neither unrecorded nor a defect. What was
        # missing is the LAYER: the census returned `n_filaments_implied` (4,086) and nothing else,
        # and a count named "implied" still reads as a count of things that exist. The day it lands
        # in a stiffness calculation nothing stops it.
        #
        # So the absence is now a field, in the output of the module that did not build it, where a
        # consumer must step over it rather than miss it.
        "n_segments": 0,
        "topology_is_built": False,
        "topology_note": (
            "NODES ONLY, BY DECISION — see this module's docstring. A closed triangulation gives "
            "z = 6 against the SOURCED z = 4, which is 1.5x the load-bearing filaments at the same "
            "node count: a stiffness error dressed as a mesh choice. Topology is a BOND question and "
            "bonds are PHASE 3, so laying one down now bakes the wrong answer into the arena layout "
            "before anyone chooses it. n_filaments_implied is WHAT PHASE 3 WOULD CLAIM, never what "
            "exists. A consumer that needs edges must REFUSE, not interpolate."),
    }


def build_lamina(
    arena: WorldArena,
    *,
    radius_um: float | None = None,
    filament_um: float | None = None,
    spacing_provenance: str | None = None,
    centre_um: tuple[float, float, float] = (0.0, 0.0, 0.0),
    population: str = "lamina",
) -> dict[str, object]:
    """Claim the meshwork's NODE range and lay its geometry down with a CUDA kernel.

    Args:
        arena: the world to claim from.  Must hold a CUDA allocation.
        radius_um: see :func:`plan_lamina`.  **REQUIRED.**
        filament_um: see :func:`plan_lamina`.  **REQUIRED.**
        spacing_provenance: where ``filament_um`` came from.  **REQUIRED, and it may not be empty** —
            the same rule :func:`~aleph.world.build.envelope.build_envelope` applies to its own
            spacing, for the same reason: the GRADE of this number is second-hand, and a grade that is
            not written into the artifact is a grade that gets quoted as a primary source later. Pass
            :data:`LAMIN_MESHWORK_SOURCE` verbatim unless a primary citation has arrived.
        centre_um: nucleus centre [µm].  A coordinate origin, so it may default; concentric with the
            envelope unless a caller says otherwise.
        population: the claim name.  A label, never a type.

    Returns:
        The census: the inventory from :func:`plan_lamina`, the claim, and the provenance.

    Raises:
        RuntimeError: if ``arena`` holds no device allocation.  There is no host construction path and
            there must never be one.
        ValueError: on any refusal from :func:`plan_lamina`, or an empty provenance.
    """
    inventory = plan_lamina(radius_um=radius_um, filament_um=filament_um)
    if not spacing_provenance or not spacing_provenance.strip():
        raise ValueError(
            "spacing_provenance must be a non-empty statement of where filament_um came from. The "
            "lamin meshwork spacing IS sourced, but SECOND-HAND (Stephens 2017 Table 1 reading Shimi "
            "2015, which this repository does not hold), and second-hand is a grade a later reader "
            "has to be able to see. Pass LAMIN_MESHWORK_SOURCE while that stands."
        )
    n_nodes = int(inventory["n_nodes"])
    if not arena.node_arrays:
        raise RuntimeError(
            f"{population}: this arena holds no device allocation, and construction here is a CUDA "
            "kernel rather than a host build that is uploaded. Warp CUDA is the only runtime; build "
            f"the arena with a CUDA device. The inventory it would have built: {n_nodes} meshwork "
            f"nodes at {inventory['areal_spacing_um_realised']:.4f} um."
        )

    nodes = arena.claim(population, Kind.NODE, n_nodes)
    # The NODE-claim ordinal, as membrane.py and cortex.py both tag it. Nothing reads range_id
    # semantically today; it is written so that a later reader can, and so two builders in one arena
    # do not tag their nodes with the same number.
    range_ordinal = len(arena.claims(kind=Kind.NODE)) - 1

    wp.launch(
        lay_shell_kernel, dim=n_nodes,
        inputs=[nodes.lo, n_nodes, range_ordinal,
                float(centre_um[0]), float(centre_um[1]), float(centre_um[2]),
                float(inventory["radius_um"])],
        outputs=[arena.node_arrays["position"], arena.node_arrays["strand_id"],
                 arena.node_arrays["range_id"]],
        device=arena.device,
    )
    arena.assert_partitioned()

    return inventory | {
        "population": population,
        "centre_um": tuple(float(x) for x in centre_um),
        "claims": {"node": (nodes.lo, nodes.hi)},
        "spacing_provenance": spacing_provenance,
    }


def _codegen_check() -> int:
    """Type-check this module's kernel by emitting its CUDA source, with no device and no execution.

    The dev machine has no card, so the alternative is shipping kernel syntax nobody has put in front
    of a compiler.  This produces NO number, so it is not a CPU result under any reading.
    """
    from warp._src.context import ModuleBuilder  # private, deliberately: there is no public codegen API

    module = wp.get_module(__name__)
    return len(ModuleBuilder(module, module.options).codegen("cuda"))


# ── self-check ──────────────────────────────────────────────────────────────────────────────────
def _demo() -> None:
    """Self-check: the two sourced numbers agreeing, and every refusal.  No device needed."""
    from aleph.world.build.envelope import nuclear_radius_um

    r_nuc = nuclear_radius_um()                      # 0.68 x 7.5 = 5.1 um, the value A1 built on
    assert math.isclose(r_nuc, 5.1)

    inv = plan_lamina(radius_um=r_nuc, filament_um=LAMIN_FILAMENT_UM)

    # THE CROSS-CHECK. Two columns of one table, taken independently: the count DERIVED from the
    # spacing must land inside the count REPORTED, and the spacing that count realises must land back
    # on the spacing reported. Neither is used to produce the other.
    assert inv["n_nodes"] == 2043, inv["n_nodes"]
    assert LAMIN_NODE_BAND[0] <= inv["n_nodes"] <= LAMIN_NODE_BAND[1]
    assert abs(inv["areal_spacing_um_realised"] - LAMIN_FILAMENT_UM) < 1e-4
    assert inv["n_filaments_implied"] == 4086, inv["n_filaments_implied"]

    # It moves with the cell, because the radius does — a swept R_cell cannot leave the meshwork at a
    # spacing that was right for a different nucleus.
    assert plan_lamina(radius_um=nuclear_radius_um("ifc") * 0.795,
                       filament_um=LAMIN_FILAMENT_UM)["n_nodes"] > inv["n_nodes"]

    # The refusal that matters: the CORTICAL band, which is what the envelope currently substitutes.
    try:
        plan_lamina(radius_um=r_nuc, filament_um=0.100)
    except ValueError as exc:
        assert "outside the independently reported band" in str(exc)
        assert "32685" in str(exc), str(exc)          # 10.9x the top of the band, quantified in situ
    else:  # pragma: no cover
        raise AssertionError("the cortical mesh spacing must refuse at the nuclear radius")

    # Both lengths refuse to default, in this package's one spelling of that refusal.
    for kwargs in ({"radius_um": r_nuc}, {"filament_um": LAMIN_FILAMENT_UM}, {}):
        try:
            plan_lamina(**kwargs)
        except ValueError as exc:
            assert "no default" in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"plan_lamina{kwargs} must refuse")

    # Provenance refuses to default AND refuses to be blank; the grade is second-hand and must travel.
    arena = WorldArena(capacity={Kind.NODE: 4000})
    for prov in (None, "", "   "):
        try:
            build_lamina(arena, radius_um=r_nuc, filament_um=LAMIN_FILAMENT_UM,
                         spacing_provenance=prov)
        except ValueError as exc:
            assert "non-empty statement" in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"spacing_provenance={prov!r} must refuse")

    # With one, it gets past every refusal and stops at the missing device, not before.
    try:
        build_lamina(arena, radius_um=r_nuc, filament_um=LAMIN_FILAMENT_UM,
                     spacing_provenance=LAMIN_MESHWORK_SOURCE)
    except RuntimeError as exc:
        assert "no device allocation" in str(exc) and "2043 meshwork nodes" in str(exc), str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device-less build must refuse")
    assert arena.n_live(Kind.NODE) == 0, "a refused build claims nothing"

    assert "SOURCED_SECONDARY" in LAMIN_MESHWORK_SOURCE and "Shimi" in LAMIN_MESHWORK_SOURCE

    chars = _codegen_check()
    assert chars > 0

    print(
        f"lamina self-check OK — CUDA codegen {chars} chars; R_nuc {r_nuc:.2f} um, filament "
        f"{LAMIN_FILAMENT_UM} um -> N = {inv['n_nodes']:,} nodes (sourced band "
        f"{LAMIN_NODE_BAND[0]:,}-{LAMIN_NODE_BAND[1]:,}), realised areal spacing "
        f"{inv['areal_spacing_um_realised']:.4f} um, {inv['n_filaments_implied']:,} filaments implied "
        f"at the sourced z={LAMIN_CONNECTIVITY} (PHASE 3, not claimed here)"
    )


if __name__ == "__main__":
    _demo()
