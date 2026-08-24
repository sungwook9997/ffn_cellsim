r"""Chromatin as a crosslinked polymer of Mbp domains — a NODE population, not an areal modulus.

**WHY THIS IS A POPULATION AND NOT A NUMBER ON THE ENVELOPE.**  ``components/incumbent`` has no
chromatin geometry at all: chromatin enters as ``k_chrom = 15.0 pN/µm`` — one of the constants that
file itself marks *"I0-B2 PI GAP, TEST values"* — folded into the per-FACE
``lamina_areal_tension_kernel``, and ``components/nucleus/chromatin_analytic.py`` is a host-side
worm-like-chain ORACLE with no runtime object behind it.  That is a lumped representation of a
mechanism, which is the pattern ``CLAUDE.md`` inverts: *"Mechanistic over lumped, at every design
decision."*  The choice was put to the Lead with both candidates and the node population was ratified
2026-08-20; the deciding argument was that **the source the incumbent borrows its chromatin
parameters FROM models chromatin as a polymer of discrete subunits**, so lumping it while quoting the
same table uses half of one source.

WHAT THE SOURCE SPECIFIES.  ``Stephens2017_MBoC`` Table 1, p7 — the same table
:mod:`aleph.world.build.lamina` reads for the meshwork spacing — carries the chromatin interior
complete, with an interpretation column:

=========================================  =======  ==============================================
quantity                                     value  interpretation
=========================================  =======  ==============================================
subunits in polymer, ``Np``                    552  1 subunit ≈ 1-10 Mbp
subunit diameter, ``σp`` [µm]                  0.6  Mbp domain size
intersubunit spring ``kp`` [nN/µm]             1.6  domain Young's modulus ≈ 1 kPa
crosslinks, ``Nc``                              55  20 % of domains linked to a distant domain
links to shell, ``Ns``                          40  ~50 % of subunits in the outer nucleus
=========================================  =======  ==============================================

Two of those check out and one does not, and this module reports all three rather than the two that
agree.  ``Nc`` reproduces EXACTLY from its own interpretation — 20 % of 552 domains, two per
crosslink, is 55.2 — so it is DERIVED here rather than typed.  ``Ns`` does not: a subunit-thick outer
shell at the nuclear radius holds ~173 of the 552 subunits, half of which is ~86, against the 40
stated.  **That reading of ``Ns`` does not hold**, so ``Ns`` is carried as the source's declared
number with the discrepancy attached, and it is a PHASE 3 bond count in any case.

⚠ **``Np`` IS A PI-GAP AND THIS MODULE REFUSES TO INVENT ONE.**  552 is a model value for a 10 µm
HeLa nucleus, not a measurement of MCF7, and the derivation that would replace it —
``Np = genome size / Mbp domain size`` — is blocked: no genome size for MCF7, which is hypertriploid
with a rearranged karyotype, is registered in this repository.  So the count arrives as a
:class:`~aleph.world.bond.BondCount`, which cannot be constructed without a scope and an authority,
and the artifact carries ``PI_GAP`` where that is the truth.  ``σp`` is likewise required, for the
smaller reason that a domain size is a physiological axis.

⚠ **THE OCCUPANCY THE SOURCE'S OWN NUMBERS IMPLY IS 11 %, AND NOTHING HERE CAN CHECK IT.**  552
domains of 0.6 µm inside a 5.1 µm nucleus fill ``Np (σp/2R)³`` = 11.2 % of it.  Interphase chromatin
occupancy is not registered in this repository — a corpus search for chromatin volume fraction and
for chromosome territories returns nothing — so this module REPORTS the fraction and does not judge
it.  It is a coarse-graining of a genome into Mbp blobs, not a packing measurement, and reading it as
one would be exactly the promotion ``CLAUDE.md`` lists.

**THE CONFIGURATION IS ONE SEEDED DRAW, AND IT SAYS SO.**  No measured chromatin spatial organisation
exists in this repository, so the rest configuration is a confined random walk of step ``σp`` from the
nuclear centre — the cortex builder's own answer to the same problem, where the seed is a declared
NUMERICAL axis and a conclusion that depends on it is a conclusion about one draw.  The walk starts at
the physiological operating point in the sense that matters here: **every backbone bond is exactly one
subunit diameter long at rest**, so no pre-strain is baked into the layout by the placement.  The
boundary is reflecting rather than clamping, which leaves the radial distribution unbiased instead of
piling subunits onto the shell.  Self-avoidance is deliberately absent: at 11 % occupancy the source's
own model uses soft potentials and its domains interpenetrate, so a hard-sphere walk would be a
different physical object.

WHAT IS DELIBERATELY ABSENT.  No intersubunit spring (``kp`` is a law), no crosslink bonds, no
chromatin-lamina tether, no heterochromatin/euchromatin split, no bending — the source connects its
subunits with extensible springs of zero bending modulus, so no ANGLE3 range is claimed and a later
bending law would be a change of model, not a missing array.  ``Nc`` and ``Ns`` are BONDS and bonds
are PHASE 3; their counts are reported so the PHASE 3 claim can be sized, and nothing is recorded as
built that is not built.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``subunit_um``, ``radius_um``, the walk step and the reflecting radius [µm]; counts
    and the volume fraction are dimensionless; the fraction is ``Np (σp/2R)³``, which is a ratio of
    volumes and carries no unit.
  * boundary — count, diameter and radius are refused without a value and refused non-positive; a
    subunit that does not fit inside the nucleus is refused; a volume fraction at or above 1 is
    refused, since more chromatin than nucleus is not a coarser model but a wrong one; a count
    resolving below 2 is refused, since a polymer with no bond is not a polymer.
  * conservation/invariant — three contiguous claims (NODE, STRAND, SEGMENT); the walk writes indices
    ``lo + i`` for ``i`` in ``[0, Np)`` only, and the loop bound IS ``Np``; the backbone topology is
    the arithmetic ``(i, i+1)`` and so no index array exists to disagree with it.
    ``assert_partitioned`` is asserted after the build.
  * CFL/precision — no integration; float64 positions.  ``wp.randf`` is float32 and is widened at the
    point of use, as ``cortex.py`` and ``nmii.py`` do, so the seed reproduces bit-for-bit.
  * sign sense — nothing is accumulated; the one signed quantity is the reflecting radius
    ``R - σp/2``, asserted positive, which is what keeps every domain wholly inside the envelope.
  * measurement protocol — one kernel of ``dim=1`` walking the chain on the device.  A walk is
    sequential and 552 steps is microseconds; the alternative is a host mirror, which is the thing the
    residency rule exists to prevent.  Nothing is read back.

engine units: length µm.  Runtime: NVIDIA Warp on CUDA; the build refuses a device-less arena.
"""

from __future__ import annotations

import math

import warp as wp

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount
from aleph.world.strand import _require_length

__all__ = [
    "CHROMATIN_SOURCE",
    "CROSSLINKED_DOMAIN_FRACTION",
    "SHELL_LINKS_DECLARED",
    "SUBUNIT_UM",
    "build_chromatin",
    "lay_chain_kernel",
    "plan_chromatin",
]

#: Chromatin subunit diameter [µm] — one Mbp-scale domain, the source's own interpretation of its
#: ``σp``.  A domain size, not a discretisation: halving it is a different coarse-graining of the
#: genome and changes what a subunit IS.
SUBUNIT_UM = 0.6

#: Fraction of domains physically linked to a distant domain.  The interpretation column of ``Nc``,
#: from which ``Nc`` reproduces exactly — which is why the fraction is the constant here and the
#: crosslink count is derived from it.
CROSSLINKED_DOMAIN_FRACTION = 0.20

#: Links to the shell as the source DECLARES them.  Carried, not derived: its own interpretation
#: ("~50 % of subunits in the outer nucleus") gives ~86 at the nuclear radius, not 40.
SHELL_LINKS_DECLARED = 40

#: The provenance sentence a caller carries verbatim into the run record.  Stephens 2017 is in the PDF
#: corpus, so unlike the lamina spacing this is a primary reading — of a MODEL, for a different cell.
CHROMATIN_SOURCE = (
    "Stephens et al. 2017 Mol Biol Cell 28:1984 (10.1091/mbc.E16-09-0653) Table 1, p7 — crosslinked "
    "polymer interior: Np=552 subunits, subunit diameter 0.6 um interpreted as an Mbp domain, 20% of "
    "domains crosslinked to a distant domain (Nc=55), Ns=40 links to the shell. ⚠ These are that "
    "paper's MODEL values for a 2R=10 um HeLa nucleus, NOT measurements of MCF7. The count in "
    "particular is a PI_GAP: deriving it needs MCF7 genome size / Mbp domain size, and no genome size "
    "for this hypertriploid line is registered in this repository."
)


@wp.kernel
def lay_chain_kernel(
    lo: wp.int32, n: wp.int32, strand_tag: wp.int32, range_ordinal: wp.int32, seed: wp.int32,
    step_um: wp.float64, r_reflect: wp.float64, cx: wp.float64, cy: wp.float64, cz: wp.float64,
    position: wp.array(dtype=wp.vec3d), strand_id: wp.array(dtype=wp.int32),
    range_id: wp.array(dtype=wp.int32),
):
    """Walk the whole chain in ONE thread, so every backbone bond is exactly one subunit long.

    Sequential by nature — subunit ``i+1`` is one step from subunit ``i`` and from nothing else — so
    this is a ``dim=1`` launch rather than a parallel placement that would have to be stitched into a
    chain afterwards.  552 steps on the device costs microseconds; the parallel alternative costs a
    host mirror, which is what may not exist.

    The boundary REFLECTS: a step leaving the ball is folded back through the sphere at the same
    direction, ``r -> 2 r_reflect - r``.  Clamping would pile subunits onto the shell and quietly
    change the radial distribution, which is the one distributional property this placement has.
    """
    if wp.tid() != 0:
        return

    state = wp.rand_init(seed, 0)
    x = wp.float64(0.0)
    y = wp.float64(0.0)
    z = wp.float64(0.0)

    for i in range(n):
        g = lo + i
        position[g] = wp.vec3d(cx + x, cy + y, cz + z)
        strand_id[g] = strand_tag
        range_id[g] = range_ordinal

        # Isotropic direction: cos(theta) uniform on [-1, 1], azimuth uniform on [0, 2pi).
        u = wp.float64(2.0) * wp.float64(wp.randf(state)) - wp.float64(1.0)
        az = wp.float64(6.283185307179586) * wp.float64(wp.randf(state))
        s = wp.sqrt(wp.max(wp.float64(0.0), wp.float64(1.0) - u * u))
        x = x + step_um * s * wp.cos(az)
        y = y + step_um * s * wp.sin(az)
        z = z + step_um * u

        r = wp.sqrt(x * x + y * y + z * z)
        if r > r_reflect:
            scale = (wp.float64(2.0) * r_reflect - r) / r
            x = x * scale
            y = y * scale
            z = z * scale


def plan_chromatin(
    *, count: BondCount, subunit_um: float | None = None, radius_um: float | None = None
) -> dict[str, object]:
    """Resolve the polymer's inventory from its declared count and the nucleus it lives in.

    Separated from :func:`build_chromatin` so the arithmetic and every refusal can be exercised
    without a card — the only part of this module a CPU machine may run.

    Args:
        count: how many subunits, for which cell, on whose authority.  **REQUIRED**, and it cannot be
            constructed without a scope and a provenance.  Basis ``explicit`` (a count per nucleus) or
            ``volumetric`` (a density per µm³ of nucleoplasm), with the support supplied at build time.
        subunit_um: chromatin domain diameter [µm].  **REQUIRED, no default.**  :data:`SUBUNIT_UM` is
            the sourced value; it is also the backbone rest length and the walk step, so it is one
            number governing three things and may not be a builder's private choice.
        radius_um: the nuclear radius the polymer is confined to [µm].  **REQUIRED, no default.**
            Derive it with :func:`~aleph.world.build.envelope.nuclear_radius_um` so the interior
            cannot outgrow the envelope that will contain it.

    Returns:
        The inventory: ``n_subunits``, the reflecting radius, the volume fraction, the derived
        crosslink count and the declared shell-link count with its discrepancy, and the count's
        provenance row.

    Raises:
        ValueError: on a missing or non-positive length; on a subunit that does not fit inside the
            nucleus; on a volume fraction at or above 1; or on a count resolving below 2, since a
            polymer with no bond is not a coarser polymer but a different object.
    """
    sigma = _require_length("subunit_um", subunit_um)
    radius = _require_length("radius_um", radius_um)

    r_reflect = radius - 0.5 * sigma
    if r_reflect <= 0.0:
        raise ValueError(
            f"subunit_um={sigma} does not fit inside radius_um={radius}: a domain of that diameter "
            "centred anywhere in the nucleus already crosses the envelope. The subunit is a physical "
            "domain size, so this is a statement about the two values, not a discretisation to relax."
        )
    return {
        "kind": "chromatin_polymer",
        "count_basis": count.basis,
        "subunit_um": sigma,
        "radius_um": radius,
        "r_reflect_um": r_reflect,
        "step_um": sigma,                       # backbone rest length == walk step == domain diameter
        "crosslinked_domain_fraction": CROSSLINKED_DOMAIN_FRACTION,
        "n_shell_links_declared": SHELL_LINKS_DECLARED,
        "provenance": {"n_subunits": _provenance_row(count), "source": CHROMATIN_SOURCE},
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


def _resolved(inventory: dict[str, object], n_subunits: int) -> dict[str, object]:
    """The count-dependent half of the inventory, once a support has resolved the count.

    Split out because everything here needs ``n_subunits`` and nothing here needs a device, so the
    volume-fraction refusal fires in :func:`plan_chromatin`'s regime rather than after a claim.
    """
    sigma = float(inventory["subunit_um"])
    radius = float(inventory["radius_um"])
    if n_subunits < 2:
        raise ValueError(
            f"the count resolves to {n_subunits} chromatin subunits. A polymer with no bond is not a "
            "coarser polymer, it is a different object — and a population that is 'present but has "
            "nothing' is indistinguishable downstream from one that failed to build."
        )
    phi = n_subunits * (sigma / (2.0 * radius)) ** 3
    if phi >= 1.0:
        raise ValueError(
            f"{n_subunits} subunits of {sigma} um fill {phi:.3f} of a {radius} um nucleus. More "
            "chromatin than nucleus is not a coarser model; check the count against the domain size "
            "before raising either."
        )
    # Nc reproduces EXACTLY from its own interpretation column: 20% of the domains, two per crosslink.
    n_crosslinks = int(round(0.5 * CROSSLINKED_DOMAIN_FRACTION * n_subunits))
    # The Ns interpretation, evaluated at THIS geometry, so the discrepancy travels as a number.
    outer_fraction = 1.0 - ((radius - sigma) / radius) ** 3 if radius > sigma else 1.0
    return {
        "n_subunits": n_subunits,
        "n_backbone_segments": n_subunits - 1,
        "volume_fraction": phi,
        "n_crosslinks_derived": n_crosslinks,
        "n_shell_links_implied_by_interpretation": int(round(0.5 * outer_fraction * n_subunits)),
    }


def build_chromatin(
    arena: WorldArena,
    *,
    count: BondCount,
    support: float = 1.0,
    subunit_um: float | None = None,
    radius_um: float | None = None,
    centre_um: tuple[float, float, float] = (0.0, 0.0, 0.0),
    seed: int = 0,
    population: str = "chromatin",
) -> dict[str, object]:
    """Claim the polymer's ranges and walk its rest configuration with a CUDA kernel.

    Args:
        arena: the world to claim from.  Must hold a CUDA allocation.
        count: see :func:`plan_chromatin`.
        support: the MEASURED support the count resolves against — 1.0 for an ``explicit`` count, the
            nucleoplasm volume [µm³] for a ``volumetric`` one.  Measured, never analytic: the
            envelope's discrete enclosed volume, not ``4/3 πR³``.
        subunit_um: see :func:`plan_chromatin`.  **REQUIRED.**
        radius_um: see :func:`plan_chromatin`.  **REQUIRED.**
        centre_um: nucleus centre [µm], and the first subunit's position.  A coordinate origin.
        seed: the walk's RNG seed.  NUMERICAL, not physiological — it selects one realisation of a
            configuration nothing has measured, and a conclusion that depends on it is a conclusion
            about one draw.  It may default because 0 is as arbitrary as any other and the number is
            recorded either way.
        population: the claim name.

    Returns:
        The census: the inventory, the resolved counts, the three claims and the seed.

    Raises:
        RuntimeError: if ``arena`` holds no device allocation.  There is no host construction path.
        ValueError: on any refusal from :func:`plan_chromatin`.
    """
    inventory = plan_chromatin(count=count, subunit_um=subunit_um, radius_um=radius_um)
    resolved = _resolved(inventory, count.resolve(support))
    n_sub = int(resolved["n_subunits"])
    if not arena.node_arrays:
        raise RuntimeError(
            f"{population}: this arena holds no device allocation, and construction here is a CUDA "
            "kernel rather than a host build that is uploaded. Warp CUDA is the only runtime; build "
            f"the arena with a CUDA device. The inventory it would have built: {n_sub} subunits + "
            f"{n_sub - 1} backbone segments at {inventory['step_um']} um."
        )

    nodes = arena.claim(population, Kind.NODE, n_sub)
    strands = arena.claim(population, Kind.STRAND, 1)
    segments = arena.claim(population, Kind.SEGMENT, n_sub - 1)
    range_ordinal = len(arena.claims(kind=Kind.NODE)) - 1

    wp.launch(
        lay_chain_kernel, dim=1,
        inputs=[nodes.lo, n_sub, strands.lo, range_ordinal, int(seed),
                float(inventory["step_um"]), float(inventory["r_reflect_um"]),
                float(centre_um[0]), float(centre_um[1]), float(centre_um[2])],
        outputs=[arena.node_arrays["position"], arena.node_arrays["strand_id"],
                 arena.node_arrays["range_id"]],
        device=arena.device,
    )
    arena.assert_partitioned()

    return inventory | resolved | {
        "population": population,
        # ⚠ A RESTATEMENT OF THIS BUILDER'S OWN CLAIM, not new information: it claims one STRAND and
        # `n_sub - 1` SEGMENTS over `n_sub` nodes, so the chain is n_strands=1 x nodes_per_strand=n_sub
        # by construction three lines above. Recorded because a consumer that cannot see the topology
        # behaves as if there is none: `world_export_cell.py` died with `KeyError: 'n_strands'` on
        # 2026-08-21 the first time the full twelve-population cell was exported, and before that had
        # silently never included this population at all. The segments were always claimed; nothing
        # said so in a form a reader could use.
        "n_strands": 1,
        "nodes_per_strand": int(n_sub),
        "centre_um": tuple(float(x) for x in centre_um),
        "seed": int(seed),
        "support": float(support),
        "claims": {"node": (nodes.lo, nodes.hi), "strand": (strands.lo, strands.hi),
                   "segment": (segments.lo, segments.hi)},
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
    """Self-check: the one derived number that reproduces, the one that does not, and the refusals."""
    from aleph.world.build.envelope import nuclear_radius_um

    r_nuc = nuclear_radius_um()                      # 0.68 x 7.5 = 5.1 um, the value A1 built on
    gap = BondCount(
        basis="explicit", value=552.0,
        scope="HeLa MODEL value at 2R=10 um, interphase — NOT a measurement of MCF7",
        source_class="PI_GAP",
        provenance=CHROMATIN_SOURCE,
    )

    inv = plan_chromatin(count=gap, subunit_um=SUBUNIT_UM, radius_um=r_nuc)
    res = _resolved(inv, gap.resolve(1.0))

    assert res["n_subunits"] == 552 and res["n_backbone_segments"] == 551
    assert math.isclose(inv["r_reflect_um"], 4.8), inv["r_reflect_um"]

    # THE ONE THAT REPRODUCES: Nc from its own interpretation column, to the unit.
    assert res["n_crosslinks_derived"] == 55, res["n_crosslinks_derived"]

    # THE ONE THAT DOES NOT: the Ns interpretation gives ~86 at this geometry, against the stated 40.
    assert res["n_shell_links_implied_by_interpretation"] == 86, res["n_shell_links_implied_by_interpretation"]
    assert SHELL_LINKS_DECLARED == 40
    assert res["n_shell_links_implied_by_interpretation"] != SHELL_LINKS_DECLARED, (
        "if these ever agree, the discrepancy this module reports has been closed and the docstring "
        "that describes it is stale"
    )

    # The occupancy the source's own numbers imply. Reported, never judged — nothing here can check it.
    assert abs(res["volume_fraction"] - 0.1124) < 1e-4, res["volume_fraction"]

    # A domain that does not fit inside the nucleus is a statement about two values, not a mesh.
    try:
        plan_chromatin(count=gap, subunit_um=2.0 * r_nuc, radius_um=r_nuc)
    except ValueError as exc:
        assert "does not fit inside" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a subunit larger than the nucleus must refuse")

    # More chromatin than nucleus is not a coarser model.
    try:
        _resolved(inv, 10_000)
    except ValueError as exc:
        assert "More chromatin than nucleus" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a volume fraction above 1 must refuse")

    # A polymer with no bond is a different object.
    try:
        _resolved(inv, 1)
    except ValueError as exc:
        assert "not a coarser polymer" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a 1-subunit polymer must refuse")

    # Both lengths refuse to default, in this package's one spelling of that refusal.
    for kwargs in ({"subunit_um": SUBUNIT_UM}, {"radius_um": r_nuc}, {}):
        try:
            plan_chromatin(count=gap, **kwargs)
        except ValueError as exc:
            assert "no default" in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"plan_chromatin{kwargs} must refuse")

    # A count with no scope cannot be constructed at all — the refusal is upstream of this module.
    try:
        BondCount(basis="explicit", value=552.0, scope="", source_class="PI_GAP", provenance="x")
    except ValueError as exc:
        assert "scope is required" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a scopeless count must refuse")

    # And a device-less arena refuses the BUILD while still reporting what it would have built.
    arena = WorldArena(capacity={Kind.NODE: 1000, Kind.STRAND: 4, Kind.SEGMENT: 1000})
    try:
        build_chromatin(arena, count=gap, subunit_um=SUBUNIT_UM, radius_um=r_nuc)
    except RuntimeError as exc:
        assert "no device allocation" in str(exc) and "552 subunits" in str(exc), str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device-less build must refuse")
    assert arena.n_live(Kind.NODE) == 0, "a refused build claims nothing"

    chars = _codegen_check()
    assert chars > 0

    print(
        f"chromatin self-check OK — CUDA codegen {chars} chars; R_nuc {r_nuc:.2f} um, subunit "
        f"{SUBUNIT_UM} um -> {res['n_subunits']} subunits / {res['n_backbone_segments']} backbone "
        f"segments, occupancy {res['volume_fraction'] * 100:.1f}% (unchecked — no measured value in "
        f"the KB); Nc derived {res['n_crosslinks_derived']} == 55 sourced; Ns interpretation "
        f"{res['n_shell_links_implied_by_interpretation']} vs {SHELL_LINKS_DECLARED} declared — "
        "REPORTED, not reconciled. Both are PHASE 3 bonds, not claimed here."
    )


if __name__ == "__main__":
    _demo()
