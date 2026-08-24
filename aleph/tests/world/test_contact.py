"""The contact type, and the one property no single connector module can check for itself.

WHAT ``_demo()`` ALREADY COVERS, so this file does not repeat it. :mod:`aleph.world.contact` and each
of the three connector modules run their own self-checks, and ``test_families_census`` invokes every
``_demo`` in isolation. Those cover the refusals, the field-level validation and the per-module
declarations.

WHAT ONLY THIS FILE CAN CHECK. :data:`~aleph.world.contact.CONTACT_IDENTITY` is a CROSS-MODULE claim —
*the three edges are ONE constraint family with three scopes* — and a module that sees only its own
scope cannot test it. That is the same reason ``test_families_census`` exists for ``DUPLICATE_OF``:
three duplicate verdicts landed within one hour, each written by a session that could not see the
others' files. This is that check for the constraint side.

⚠ AND IT IS A RATCHET, NOT A DESCRIPTION. Every assertion here fails on the way the type is most
likely to be misused later — a count field appearing, a scope acquiring a private gap law, a stiffness
appearing without being declared as a regularisation, a fixture ``d0`` being pinned as a magnitude. The
`0.030 um` in the port source's tests is exactly such a fixture and is the value most likely to be
promoted by someone who needs a number to make a run start.

Runtime: host-only. Nothing here touches a device, and that is deliberate — a constraint declaration
is host bookkeeping, and a test that needed CUDA to check a dataclass would be testing the wrong thing.
"""

from __future__ import annotations

import dataclasses

import pytest

from aleph.world.contact import (
    CONTACT_IDENTITY,
    GAP_LAW_POINT_SURFACE,
    ContactDistance,
    ContactFamily,
    ContactScope,
    check_complementarity,
)
from aleph.world.families import lamellipodium_membrane_contact as d12
from aleph.world.families import membrane_cortex_contact as d2
from aleph.world.families import nucleus_cortex_contact as d6
from aleph.world.bond import SourceClass

#: The three modules PI decision 1 names. ``membrane_ecm_contact`` is the fourth scope of the same law
#: and has no module yet — recorded here rather than silently omitted, because a census that quietly
#: covers three of four reads as complete.
RULED = (d2, d6, d12)


def test_the_three_edges_are_one_family_and_not_three() -> None:
    """CONTACT_IDENTITY's cross-module claim, which no single module can make for itself."""
    scopes = tuple(m.CONTACT for m in RULED)
    family = ContactFamily(gap_law=GAP_LAW_POINT_SURFACE, scopes=scopes)

    assert len(family.scopes) == 3
    assert {s.connector for s in family.scopes} == {
        "membrane_cortex_contact", "nucleus_cortex_contact", "lamellipodium_membrane_contact"
    }
    # ONE gap law across all three. A module that needed its own would be declaring a second family,
    # and the merge this decision rests on would be a name rather than a physics statement.
    assert family.gap_law == GAP_LAW_POINT_SURFACE
    assert "gap function" in CONTACT_IDENTITY and "SCOPE, not identity" in CONTACT_IDENTITY

    # All three are blocked, all three on a contact distance, and none of them on a count.
    assert family.blocked().keys() == {s.connector for s in scopes}
    for connector, reasons in family.blocked().items():
        joined = " ".join(reasons).lower()
        assert "d0" in joined, f"{connector}: blocked on something other than the contact distance"
        for count_word in ("how many", "against what", "areal density", "per vertex"):
            assert count_word not in joined, (
                f"{connector}: a COUNT question is in the live queue for a constraint. It was "
                "dissolved by PI decision 1, not deferred — the active set is an output of the broad "
                "phase and changes every step."
            )


def test_a_count_has_nowhere_to_live() -> None:
    """The ERM defect is kept out by the SHAPE of the type, not by a convention anyone must remember."""
    fields = {f.name for f in dataclasses.fields(ContactScope)}
    fields |= {f.name for f in dataclasses.fields(ContactFamily)}
    forbidden = {"count", "n_contacts", "n_pairs", "pairs", "capacity", "density", "support"}
    assert not forbidden & fields, (
        f"a count-shaped field has appeared on the contact type: {sorted(forbidden & fields)}. "
        "A contact count is a function of the CONFIGURATION; storing one is how a broad-phase number "
        "acquires a physiological label, which is the defect bond.py is shaped around."
    )


def test_declaring_a_constraint_claims_no_arena_capacity() -> None:
    """A constraint is an address resolved from live geometry (arena.py:76-78), not an allocation."""
    from aleph.world.arena import Kind, WorldArena

    arena = WorldArena(capacity={Kind.NODE: 16, Kind.SEGMENT: 16, Kind.ANGLE3: 16, Kind.BOND: 16})
    ContactFamily(gap_law=GAP_LAW_POINT_SURFACE, scopes=tuple(m.CONTACT for m in RULED))
    assert arena.n_live(Kind.BOND) == 0
    arena.assert_partitioned()


def test_no_scope_is_enforceable_and_no_d0_has_been_pinned() -> None:
    """The ruling gave these edges a destination. It did not give them a magnitude."""
    for m in RULED:
        assert not m.CONTACT.enforceable, (
            f"{m.CONTACT.connector}: a contact distance has been supplied. If it is sourced this test "
            "is what should change; if it is the 0.030 um fixture from the port source's tests, it is "
            "a fixture wearing a physiological label."
        )


def test_every_scope_names_the_wrong_value_it_is_most_likely_to_be_given() -> None:
    """A KB sweep (2026-08-21, 12,810 chunks) found no d0 for any scope — and a NEAR-MISS for each.

    ⚠ That is the dangerous shape of this gap. Each scope has, sitting next to it in the corpus, a
    quantity of the right dimension and the wrong meaning:

    * membrane/cortex — the cortex **thickness**, 100-200 nm, sourced in five papers. How thick the
      layer is, not how far it sits from the bilayer.
    * envelope/cortex — ``Fang2016_PhysRevE``'s repulsion, whose gap is a RATIO of change against the
      INITIAL configuration, so the build sets the reference. No absolute d0 exists in it at all, and
      its ``source_audit`` verdict is CHECK.
    * lamellipodium/membrane — the ratchet's ``delta = 2.7 nm``, a monomer half-length: the step a
      polymerisation event advances by, not a separation at which a force vanishes.

    Each is dimensionally correct, well sourced, and wrong. Naming it in the blocker is what stops it
    being rediscovered as an answer — an empty blocker would leave the near-miss as the only number
    in reach.
    """
    traps = {
        "membrane_cortex_contact": "thickness",
        "nucleus_cortex_contact": "initial configuration",
        "lamellipodium_membrane_contact": "2.7 nm",
    }
    for m in RULED:
        joined = " ".join(m.CONTACT.blocked_by)
        trap = traps[m.CONTACT.connector]
        assert trap.lower() in joined.lower(), (
            f"{m.CONTACT.connector}: the blocker no longer names the near-miss value ({trap!r}). "
            "The gap is not merely unanswered — a wrong number of the right dimension is in reach, "
            "and the blocker is what keeps it from being read as the answer."
        )


def test_the_bond_census_still_reads_these_as_no_family() -> None:
    """A destination for the constraint must not read as a bond becoming buildable."""
    for m in RULED:
        assert not m.SPEC.buildable, f"{m.SPEC.connector}: ruled not a bond, yet the SPEC reads buildable"
        assert m.SPEC.blocked_by


def test_the_ratchet_edge_carries_two_shapes_not_one() -> None:
    """D12's constraint is unconditional; its bond waits on the elastic/tethered reading.

    ⚠ This is where PI decision 1's wording and ``contracts.py:791`` disagree — the decision calls
    this edge cardless and it is not. The disagreement is resolved by the edge being two things, and
    that resolution is only true while the two sets of questions stay separate.
    """
    assert not set(d12.CONTACT.blocked_by) & set(d12.BLOCKED_BY)
    joined = " ".join(d12.BLOCKED_BY)
    assert "is-this-a-bond" in joined, "the bond question must survive the constraint ruling"
    assert "elastic ratchet" in joined and "tethered ratchet" in joined


def test_complementarity_is_scored_and_the_three_failures_are_distinguishable() -> None:
    """0 <= Phi ⊥ gamma >= 0 — and which of the three clauses broke, because they mean different things."""
    ok = check_complementarity([0.5, 0.0], [0.0, 12.0], gap_tol_um=1e-9, force_tol_pn=1e-9)
    assert ok["penetration_um"] == 0.0 and ok["adhesion_pn"] == 0.0 and ok["product_pn_um"] == 0.0

    # An empty active set is a physical answer: nothing is touching.
    assert check_complementarity([], [], gap_tol_um=0.0, force_tol_pn=0.0)["n_constraints"] == 0.0

    for phi, gamma, expect in (
        ([-0.02], [5.0], "non-penetration violated"),
        ([0.10], [-5.0], "unilaterality violated"),
        ([0.10], [5.0], "complementarity violated"),
    ):
        with pytest.raises(AssertionError, match=expect):
            check_complementarity(phi, gamma, gap_tol_um=1e-9, force_tol_pn=1e-9)


def test_the_product_tolerance_is_derived_and_not_a_free_threshold() -> None:
    """A chosen product tolerance would be a threshold edited after seeing a residual.

    The admissible product is ``gap_tol*max|gamma| + force_tol*max|Phi|`` and nothing else — the
    first-order error the two supplied tolerances already admit. So a gap sitting exactly at its own
    tolerance while loaded PASSES, and one at twice it FAILS, with no third number involved.
    """
    at_tolerance = check_complementarity([-1e-3], [10.0], gap_tol_um=1e-3, force_tol_pn=0.0)
    assert at_tolerance["product_pn_um"] == pytest.approx(at_tolerance["product_tol_pn_um"])

    with pytest.raises(AssertionError, match="non-penetration violated"):
        check_complementarity([-2e-3], [10.0], gap_tol_um=1e-3, force_tol_pn=0.0)

    # And the signature has no place to put one, which is the durable half of this check.
    import inspect

    params = set(inspect.signature(check_complementarity).parameters)
    assert params == {"phi_um", "gamma_pn", "gap_tol_um", "force_tol_pn"}


def test_no_verdict_rests_on_a_threshold_the_row_does_not_carry() -> None:
    """A PASS must be falsifiable from the artifact alone, and two of the three once were not.

    ⚠ Caught by Lead reading the function rather than its description, 2026-08-21. The row recorded
    ``penetration_um`` and ``adhesion_pn`` but not the tolerances they were compared against, so a
    PASS at ``gap_tol_um=1e-9`` and a PASS at ``1e-3`` produced identical artifacts. A provenance
    stamp does not close that: it fixes WHICH CODE RAN and the tolerance is an ARGUMENT.

    The check is derived from the signature rather than from a written list, so a fourth tolerance
    added later fails here until it is put in the row too.
    """
    import inspect

    row = check_complementarity([0.5], [0.0], gap_tol_um=1e-9, force_tol_pn=1e-9)
    tol_params = {p for p in inspect.signature(check_complementarity).parameters if "_tol" in p}
    assert tol_params, "the signature no longer names its tolerances; this check has gone vacuous"
    missing = tol_params - row.keys()
    assert not missing, (
        f"these thresholds judge a verdict but are absent from the artifact row: {sorted(missing)}. "
        "A row that carries a residual without the tolerance it was compared against turns a PASS "
        "into an unfalsifiable claim — worse than a wrong number, which can at least be checked."
    )

    # And the recorded thresholds are the ones actually used, not defaults echoed back.
    other = check_complementarity([0.5], [0.0], gap_tol_um=2.5e-4, force_tol_pn=7.0)
    assert other["gap_tol_um"] == 2.5e-4 and other["force_tol_pn"] == 7.0


def test_the_gate_does_not_loosen_with_the_number_of_constraints() -> None:
    """⚠ Session A, 2026-08-21: a Higham-style tolerance grows with the TERM COUNT, so a gate is
    loosened by the very events it is meant to judge — measured at 8.2M added force terms, residual
    unchanged, tolerance 86x larger.

    A contact set is the worst possible case for that feedback, because its size is a function of the
    configuration and therefore of the step being judged. The defence here is structural and is worth
    stating precisely rather than as a lucky side effect: **every quantity in this contract is a
    max-norm, never a sum.** ``penetration_um``, ``adhesion_pn`` and ``product_pn_um`` are maxima, and
    the derived ``product_tol_pn_um`` is built from ``max|gamma|`` and ``max|Phi|`` — all four are
    INTENSIVE in the constraint count.

    So a step that activates a million contacts is judged by the same numbers as one that activates
    ten. Duplicating the entire active set changes nothing at all, which this asserts on the nose. A
    residual that summed, or a tolerance that scaled with ``C``, would fail here.
    """
    phi = [0.4, -1e-4, 0.0, 0.9]
    gamma = [0.0, 3.0, 7.5, 0.0]
    small = check_complementarity(phi, gamma, gap_tol_um=1e-3, force_tol_pn=1e-9)
    large = check_complementarity(phi * 2048, gamma * 2048, gap_tol_um=1e-3, force_tol_pn=1e-9)

    assert large["n_constraints"] == small["n_constraints"] * 2048
    # Every key except the count itself, derived rather than listed — a quantity added later is
    # covered without anyone remembering to add it here.
    for key in small.keys() - {"n_constraints"}:
        assert large[key] == small[key], (
            f"{key} moved when the active set grew 2048x with identical physics. Every quantity in "
            "this contract must be intensive in the constraint count, or activating contacts would "
            "loosen the gate that judges the step that activated them."
        )

    # And the violation threshold does not move either: a single bad pair among a million is still a
    # violation, not a dilution.
    with pytest.raises(AssertionError, match="non-penetration violated"):
        check_complementarity([0.4] * 100_000 + [-2e-3], [0.0] * 100_000 + [3.0],
                              gap_tol_um=1e-3, force_tol_pn=1e-9)


def test_a_stiffness_cannot_appear_without_declaring_itself_a_regularisation() -> None:
    """None is the physics (gamma is a multiplier); a number is an approximation of it, and says so."""
    scopes = tuple(m.CONTACT for m in RULED)
    assert ContactFamily(gap_law=GAP_LAW_POINT_SURFACE, scopes=scopes).is_complementarity

    penalised = ContactFamily(
        gap_law=GAP_LAW_POINT_SURFACE, scopes=scopes, regularisation_stiffness_pn_per_um=1.0e3)
    assert not penalised.is_complementarity
    assert penalised.provenance_row()["formulation"] == "regularised_penalty"

    # Zero is neither form and would disable the constraint while looking configured.
    with pytest.raises(ValueError, match="silently disable"):
        ContactFamily(gap_law=GAP_LAW_POINT_SURFACE, scopes=scopes,
                      regularisation_stiffness_pn_per_um=0.0)


def test_the_same_surface_pair_twice_is_the_duplication_the_shape_exposes() -> None:
    """Under CONTACT_IDENTITY the pair is scope, so a repeated pair is one constraint declared twice."""
    scopes = tuple(m.CONTACT for m in RULED)
    with pytest.raises(ValueError, match="declared twice"):
        ContactFamily(gap_law=GAP_LAW_POINT_SURFACE, scopes=scopes + (scopes[0],))


def test_a_contact_distance_refuses_to_be_a_bare_number() -> None:
    """The same four questions BondCount asks, minus the ones a constraint cannot be asked."""
    good = ContactDistance(0.03, "MCF7 interphase adherent 37C", SourceClass.PI_GAP, "no datum yet")
    assert good.value_um == 0.03

    for kwargs, expect in (
        ({"value_um": 0.0}, "absence of a contact"),
        ({"value_um": float("inf")}, "finite"),
        ({"scope": "   "}, "scope is required"),
        ({"provenance": ""}, "provenance is required"),
    ):
        with pytest.raises(ValueError, match=expect):
            dataclasses.replace(good, **kwargs)


def test_a_scope_must_be_answered_or_blocked_and_never_both() -> None:
    """A partial answer that reads as a whole one is how a fixture becomes a magnitude."""
    with pytest.raises(ValueError, match="refusal with no content"):
        ContactScope(connector="x", point_population="a", surface_population="b")

    d0 = ContactDistance(0.03, "fixture", SourceClass.CONVENIENCE, "test_connector_joints.py")
    with pytest.raises(ValueError, match="Say which"):
        ContactScope(connector="x", point_population="a", surface_population="b", d0=d0,
                     blocked_by=("still open",))

    with pytest.raises(ValueError, match="cannot be kept out of itself"):
        ContactScope(connector="x", point_population="a", surface_population="a",
                     blocked_by=("open",))


def test_the_surface_side_is_named_and_every_declared_one_owns_faces() -> None:
    """⚠ The lamina question (PI queue item 17), answered against the type rather than in prose.

    A non-penetration constraint is about a REGION a point may not enter, and ``surface.py`` records
    that a FACE is *"the only primitive whose measure is an AREA"*. So the surface side must own
    FACEs. ``build/lamina.py`` claims ``Kind.NODE`` and nothing else — 2,043 points whose
    ``n_filaments_implied`` is derived from a sourced connectivity rather than built — so the lamina
    could not be one; ``build/envelope.py`` claims ``NODE`` + ``FACE`` + ``ANGLE4`` and can.

    ⚠ Before roles existed this was not merely unchecked, it was UNSAYABLE — and the three declared
    scopes had already drifted apart on order without anything noticing: two were written
    surface-first, one point-first. That drift is the evidence the role naming was not decoration.

    This test pins the roles that ARE declared. It cannot pin the FACE requirement itself, because
    that needs a population's primitive census, which belongs to the builders and not to a
    declaration type — recorded in the class docstring rather than faked here.
    """
    surfaces = {m.CONTACT.surface_population for m in RULED}
    assert surfaces == {"membrane", "nuclear_envelope"}, (
        f"a scope names {sorted(surfaces)} as a surface. Every surface side must own FACEs; a "
        "NODE-only population (the nuclear lamina is the live example) bounds no region and cannot "
        "be one. See PI queue item 17."
    )
    for m in RULED:
        s = m.CONTACT
        assert s.point_population != s.surface_population
        assert s.populations == (s.point_population, s.surface_population), "role order, not declaration order"
