r"""``filopodium_membrane_tip`` — the bundled barbed-end tip pushing into the plasma membrane.

**THIS MODULE BUILDS NOTHING, AND THAT IS ITS PRODUCT.** Three separate reasons, each sufficient on its
own, and they are independent so closing one does not close the others:

  1. it is a **DUPLICATE** of ``lamellipodium_membrane_contact`` — same physics, argued below;
  2. it is a **KINETIC** family, and :mod:`aleph.world.bond` deliberately has no kinetics yet;
  3. its **count is unanswered** — nine questions, listed in :data:`SPEC` ``.blocked_by``.

---------------------------------------------------------------------------------------------------
1. WHY IT IS A DUPLICATE — the argument is from physics, and the port source states it as a test
---------------------------------------------------------------------------------------------------

The criterion is :data:`~aleph.world.families.DUPLICATE_CRITERION`, converged by three sessions on
2026-08-20: **a family's identity is its CHEMISTRY CARD; a shared runtime class is not a duplicate and
a shared force law is not a duplicate.** This edge meets it on the card, which is the decisive test,
and the shared law and runtime are corroboration that is not on its own sufficient:

  * **Same chemistry card — the criterion.** ``contracts.py:795`` and ``:841`` both write
    ``chemistry_card="actin_membrane_brownian_ratchet_contact"`` — the same string, twice, on the two
    edges. :class:`~aleph.world.bond.BondFamily` names the card *"the identity that actually
    distinguishes families"* and says two edges sharing one **are** one family.

**Two controls, and they are inside this connector's own neighbourhood, so folding by card is
demonstrably not folding by name:**

  * ``filopodium_cortex_root`` is the OTHER edge on this very population and it is NOT folded: same
    ``filopodium`` component, and D13 records that it carries ``formin_fascin_root_coupling`` while its
    own ``TRANSIENT_ACTIN`` sibling carries ``transient_actin_crosslink``. Sharing a component folds
    nothing.
  * The three steric contacts — ``membrane_cortex_contact``, ``membrane_ecm_contact``,
    ``nucleus_cortex_contact`` — share this edge's FORCE LAW exactly and are NOT folded, because their
    cards differ and ``kinetics`` is False. Sharing a law folds nothing either. That is the criterion's
    second clause, met by a case rather than by assertion.

The corroboration, from the frozen port source arguing it against itself:

  * **Same force law.** ``engine/connector_joints.py`` opens by saying all five ``CONTACT`` edges are
    *"the same unilateral gap law"* and that the module therefore adds *"exactly ONE kernel"*. Both
    edges get it; ``test_a_contact_edge_gets_the_unilateral_law_and_nothing_else_does`` asserts so.
  * **Same kinetic law.** ``test_the_two_contact_physics_are_told_apart_by_kinetics_not_by_the
    _chemistry_card`` resolves every contact edge to a ``JointKind`` and asserts, literally,
    ``ratchet == {"lamellipodium_membrane_contact", "filopodium_membrane_tip"}`` — the two of them and
    nothing else share ``JointKind.BROWNIAN_RATCHET_CONTACT``. The other three are ``STERIC_CONTACT``.
  * **Same runtime.** Both are ``ContractJointConnector`` in the committed device-run census, both
    ``NOT_BOUND``.
  * **The only difference is which population the actin node sits in** — ``"growing barbed-end material
    point"`` vs ``"bundled barbed-end tip"``. That is precisely what ``bond.py`` refuses to store: a
    bond holds two GLOBAL node indices and the component pair is read off the arena's ID ranges at
    query time. A merged family reports ``{("filopodium", "membrane"): n, ("lamellipodium",
    "membrane"): m}`` from :meth:`~aleph.world.bond.BondFamily.component_pairs` with no second
    declaration anywhere.

**What could have refuted it, and did not.** A different card would have — that is the criterion, and
D9's ``if_sf_plectin`` is the case where it came out the other way on exactly this shape. A parameter
difference would not be enough — ``rest_um``
and ``stiffness_pn_per_um`` are per-bond arrays, so a filopodial tip and a lamellipodial edge may carry
different contact distances inside ONE family. Membrane tension resisting a finger is likewise not a
counter-example: that force belongs to the membrane component's own law, not to this connector. What
WOULD refute it is a different force law or a different rate law, and the port source's own tests say
there is neither.

**Consequence for whoever owns the survivor** (``lamellipodium_membrane_contact``, session D14): a
:class:`~aleph.world.bond.BondFamily` carries ONE ``BondCount`` and ONE support. So the lamellipodial
and filopodial contacts must be counted in a COMMON basis. That is a real constraint the merge creates
and it is listed in ``blocked_by``, not resolved here.

---------------------------------------------------------------------------------------------------
2. WHY IT IS ALSO THE KINETIC FIRST — escalate, do not build
---------------------------------------------------------------------------------------------------

``bond.py`` §"WHAT IS DELIBERATELY ABSENT" defers attach/detach, the free list and snapshot twins to
*"the first family that has kinetics"*. This edge is exactly that: ``kinetics=True`` and
``commit_on_accept=True`` in the contract, a ``BROWNIAN_RATCHET_CONTACT`` rate law, and
``remap_on_accept=True``. ``connector_joints.py`` is emphatic that its ``propose_events`` **raises**
without an injected rate law, because *"a kinetic connector that silently proposes nothing is
indistinguishable from a permanent weld"* — the one thing the charter says a connector may never be.

So a family built here would either be a weld wearing a ratchet's name, or it would require inventing
the very shape ``bond.py`` refused to guess. Both are forbidden. **PI decision.**

---------------------------------------------------------------------------------------------------
3. THE ERM CHECK — was the count a mesh number wearing a physiological label?
---------------------------------------------------------------------------------------------------

Asked, and the answer is **no, but the mesh has its hand on a different lever.**

Not the ERM defect: the ERM count was resolved against the membrane VERTEX SET, so the discretisation
set the population. Here the population is anchored on the actin side — barbed ends of filopodial
filaments — which is a filament count, not a mesh count. Refining the membrane does not create tips.

But the mesh decides whether that count is REPRESENTABLE, which is the same defect with the roles
swapped. Arithmetic over the committed PHASE-1 record (``world_phase1/phase1a_native.json``,
``phase1_native.json``) — geometry only, no physics claimed:

  * membrane is subdivision 7: 163,842 vertices, 327,680 faces, mean edge **0.0706 µm** over
    706.8 µm² — i.e. 231.8 vertices/µm², which is the very number ``bond.py`` names as the mesh
    number that wore the label "physiological 235/µm²";
  * a filopodial bundle is 20 filaments hex-packed at 0.008 µm: 5 columns × 4 rows, so its
    cross-section is **0.032 × 0.021 µm** — about a THIRD of one membrane triangle.

So all 20 barbed ends of one finger fall inside one or two triangles and would bond to the SAME
membrane vertex. "One bond per barbed end" and "one bond per bundle tip" are therefore not two
discretisations of one model: at this mesh they differ by a factor of 20 in the stiffness delivered to
a single vertex, and choosing the first without saying so is a load-sharing decision made silently by
the mesh. It is a decision, so it is declared (``blocked_by``) rather than taken.

:func:`bundle_fits_inside_one_membrane_triangle` is that arithmetic, kept runnable so the finding
fails loudly when the subdivision or the fascin spacing changes rather than going stale in prose. It is
robust to the unsourced input: with 20 filaments the bundle stays inside one triangle for any spacing
below ~0.0176 µm, and the fascin figure is 0.007–0.008 µm.

⚠ Both filopodial counts are themselves ``PI_GAP`` placeholders — PHASE 1 passed 50 fingers and 20
filaments per bundle through ``gap()``, whose scope field reads *"NO VALUE EXISTS"*. So even the one
basis that looks derivable (one contact per barbed end → 1,000) is a product of two numbers with no
source. ``world/build/filopodium.py`` says the per-cell filopodium count is *absent from the KB, not
merely unratified*.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — nothing is constructed, so no force, length or stiffness is produced here. The one
    computation is a length comparison, µm on both sides, and it is asserted to be so.
  * boundary — the refusal is unconditional: :func:`build_filopodium_membrane_tip` raises for every
    input, including a fully-specified one, because the reasons are the duplication and the kinetics
    and neither is an argument the caller can supply.
  * conservation/invariant — no bond is claimed, so the arena is untouched and ``assert_partitioned``
    is unaffected. This module allocates nothing.
  * CFL/precision — no integration; the arithmetic is float64 host arithmetic on committed record
    values.
  * sign sense — recorded for the family that will eventually exist, NOT implemented here: the
    unilateral law is non-zero only when the gap is below the contact distance and then pushes the tip
    and the membrane APART. A contact that pulls is an adhesion nobody declared, and it is a separate
    edge (``filopodium_nascent_fa``) when it is wanted.
  * measurement protocol — host only; no device is touched, nothing is uploaded, nothing is read back.

engine units: length µm, stiffness pN/µm.  Runtime: pure host bookkeeping, CPU-importable.
"""

from __future__ import annotations

import math

from aleph.world.families import ConnectorGapError, FamilySpec

__all__ = [
    "CONNECTOR",
    "DUPLICATE_OF_ENTRY",
    "SPEC",
    "build_filopodium_membrane_tip",
    "bundle_fits_inside_one_membrane_triangle",
]

CONNECTOR = "filopodium_membrane_tip"

#: The survivor of the duplication, for Lead to fold into ``families.DUPLICATE_OF``. It is exposed
#: here rather than written into ``families/__init__.py`` because fourteen sessions share that file
#: and one module per connector is the only isolation this package has.
DUPLICATE_OF_ENTRY: dict[str, str] = {CONNECTOR: "lamellipodium_membrane_contact"}

#: The unanswered questions, one per name, in the order Lead should queue them: the two rulings first,
#: then the count, then the parameters that have no home anywhere in the tree.
_BLOCKED_BY: tuple[str, ...] = (
    "DUPLICATE_OF:lamellipodium_membrane_contact — same JointKind.BROWNIAN_RATCHET_CONTACT, same "
    "chemistry card, same unilateral kernel, same runtime. This module builds nothing; the merged "
    "family belongs to that connector's owner.",
    "KINETIC_FAMILY_FIRST — bond.py defers attach/detach, the free list and snapshot twins to 'the "
    "first family that has kinetics'. This edge is kinetics=True + commit_on_accept=True with a "
    "ratchet rate law, so it is that family. PI decides the shape; it may not be invented here.",
    "n_filopodia_per_cell — ABSENT from the Contract-Graph, not merely unratified "
    "(world/build/filopodium.py). PHASE 1 passed 50 as an explicit PI_GAP placeholder.",
    "filaments_per_bundle — architecture_spec gives 20 in a 10-30 band and calls the bundler itself a "
    "PI-gated stand-in. PHASE 1 passed 20 as an explicit PI_GAP placeholder.",
    "contact_count_model — one bond per barbed end, or one per bundle tip? At membrane subdivision 7 "
    "a whole 20-filament bundle sits inside one triangle, so the two differ by 20x in the stiffness "
    "delivered to a single membrane vertex. A load-sharing decision, not a discretisation detail.",
    "capacity_or_instantaneous — the contract sets remap_on_accept=True, i.e. the pair set is meant to "
    "be re-discovered by a broad phase on each accepted step, while BondCount resolves ONE fixed "
    "number. Is the count the capacity of contactable tips, or the bound population? No broad phase "
    "exists in either tree.",
    "contact_distance_um — no value anywhere; connector_joints.py states on purpose that it defines "
    "no contact distance, so the gap at which the force vanishes has never been sourced.",
    "stiffness_pn_per_um — same: every magnitude was to arrive on a spec, and no spec was ever built.",
    "common_basis_with_lamellipodium — a BondFamily carries ONE BondCount and ONE support, so merging "
    "the two contacts forces one basis across a bundled tip population and a branched leading edge. "
    "Which basis is a consequence of the merge and is not settled by either connector alone.",
)

#: What this module declares about itself. ``basis`` records the CANDIDATE basis — one contact per
#: barbed end — and the fifth entry of ``blocked_by`` is the open question of whether that is right.
SPEC = FamilySpec(
    connector=CONNECTOR,
    populations=("filopodium", "membrane"),
    basis="per_filament",
    blocked_by=_BLOCKED_BY,
    source="none — see blocked_by; both filopodial counts are PI_GAP placeholders and neither the "
           "contact distance nor the stiffness exists in any tree.",
)


def bundle_fits_inside_one_membrane_triangle(
    *, n_filaments: int, spacing_um: float, membrane_edge_um: float
) -> bool:
    """Whether a hex-packed bundle's cross-section is narrower than one membrane triangle edge.

    The hex packing is ``build_filopodium``'s own: ``ceil(sqrt(n))`` columns, the rest in rows, odd
    rows staggered by half a spacing and rows pitched by ``sqrt(3)/2``. The width is the larger extent,
    so this is the generous comparison — the bundle is smaller still along the other axis.

    Args:
        n_filaments: filaments in the bundle.
        spacing_um: centre-to-centre filament spacing [µm].
        membrane_edge_um: mean membrane mesh edge length [µm].

    Returns:
        True when every barbed end of one finger lands within a single triangle, i.e. when the mesh
        cannot tell one bond per barbed end from one bond per bundle tip.

    Raises:
        ValueError: on a non-positive count, spacing or edge length.
    """
    if n_filaments < 1 or spacing_um <= 0.0 or membrane_edge_um <= 0.0:
        raise ValueError(
            f"n_filaments={n_filaments}, spacing_um={spacing_um}, membrane_edge_um={membrane_edge_um}: "
            "a bundle needs at least one filament and both lengths must be positive"
        )
    cols = math.ceil(math.sqrt(n_filaments))
    rows = math.ceil(n_filaments / cols)
    width_um = (cols - 1) * spacing_um
    height_um = (rows - 1) * spacing_um * (3.0 ** 0.5) / 2.0
    return max(width_um, height_um) < membrane_edge_um


def build_filopodium_membrane_tip(*_args: object, **_kwargs: object) -> None:
    """Refuse, with the nine unanswered questions by name.

    Present so a caller that reaches for this family gets the refusal and its reasons rather than an
    ``AttributeError``, which would read as "not written yet" instead of "must not be written".

    Raises:
        ConnectorGapError: always.
    """
    raise ConnectorGapError(
        CONNECTOR,
        list(_BLOCKED_BY),
        note=(
            "This connector is the SAME PHYSICS as lamellipodium_membrane_contact and must not build a "
            "second family; it is also the first kinetic family, which bond.py deferred on purpose. "
            "Both are PI decisions, and neither is an argument a caller can supply."
        ),
    )


def _demo() -> None:
    """Self-check: the module refuses, names why, and the mesh arithmetic is re-run rather than quoted."""
    assert not SPEC.buildable, "a module that builds nothing must not report itself buildable"
    assert SPEC.connector == CONNECTOR
    assert SPEC.populations == ("filopodium", "membrane")
    assert len(SPEC.blocked_by) == 9, "nine questions is the product; losing one loses a decision"
    assert DUPLICATE_OF_ENTRY == {CONNECTOR: "lamellipodium_membrane_contact"}

    # The refusal is unconditional — a fully-specified call is refused for the same reasons.
    for args, kwargs in (((), {}), ((None,), {"count": 1000, "stiffness_pn_per_um": 100.0})):
        try:
            build_filopodium_membrane_tip(*args, **kwargs)
        except ConnectorGapError as exc:
            assert exc.connector == CONNECTOR
            assert len(exc.missing) == 9
            assert "DUPLICATE_OF" in exc.missing[0] and "KINETIC_FAMILY_FIRST" in exc.missing[1]
        else:  # pragma: no cover
            raise AssertionError("this family must never be built")

    # PHASE 1 as committed: membrane subdiv 7 mean edge 0.0706 um (phase1a_native.json), bundle of 20
    # at the 8 nm fascin spacing (world_phase1_native.py). The whole finger sits in one triangle.
    assert bundle_fits_inside_one_membrane_triangle(
        n_filaments=20, spacing_um=0.008, membrane_edge_um=0.0706)
    # ...and it stays true across the 10-30 band and the 7-8 nm spacing, so the finding is not an
    # artefact of the two PI_GAP placeholders.
    for n in (10, 20, 30):
        for s in (0.007, 0.008):
            assert bundle_fits_inside_one_membrane_triangle(
                n_filaments=n, spacing_um=s, membrane_edge_um=0.0706), (n, s)
    # It is a real comparison and not a tautology: coarsen the bundle and it stops holding.
    assert not bundle_fits_inside_one_membrane_triangle(
        n_filaments=20, spacing_um=0.030, membrane_edge_um=0.0706)
    try:
        bundle_fits_inside_one_membrane_triangle(n_filaments=0, spacing_um=0.008, membrane_edge_um=0.07)
    except ValueError as exc:
        assert "at least one filament" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an empty bundle must refuse")

    print(f"{CONNECTOR}: BLOCKED — {len(SPEC.blocked_by)} unanswered; duplicate of "
          f"{DUPLICATE_OF_ENTRY[CONNECTOR]}; kinetic-first escalated to PI")


if __name__ == "__main__":
    _demo()
