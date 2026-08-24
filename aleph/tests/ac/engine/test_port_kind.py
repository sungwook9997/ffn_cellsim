"""`PortKind` — the invariant that makes "how many connectors" a question with an answer.

The declaration count mixes four kinds of object, so it answers nothing cleanly. These tests pin the
derivation, and they pin the counts, because the counts are what a reader will quote.
"""

from __future__ import annotations

import pytest

from aleph.engine.contracts import (
    CellArchitecture,
    ComponentContract,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
    ConnectorScope,
    PortKind,
    reference_cell_architecture,
)

ARCH = reference_cell_architecture()


def test_every_declared_connector_classifies() -> None:
    """No connector may fall outside the four kinds — an unclassified edge would be uncountable."""
    kinds = {c.name: ARCH.port_kind(c.name) for c in ARCH.connectors}
    assert len(kinds) == len(ARCH.connectors)
    assert all(isinstance(k, PortKind) for k in kinds.values())


def test_the_counts_are_what_the_repository_quotes() -> None:
    """38 declarations = 28 power ports + 6 series halves + 2 reservoir ports + 2 internal elements.

    The six series halves are THREE springs (a mechanical_group binds one composite runtime), so the
    number of places two bodies exchange force is 31 — not 38.
    """
    kinds = [ARCH.port_kind(c.name) for c in ARCH.connectors]
    assert len(kinds) == 38
    assert kinds.count(PortKind.POWER_PORT) == 28
    assert kinds.count(PortKind.SERIES_HALF) == 6
    assert kinds.count(PortKind.RESERVOIR_PORT) == 2
    assert kinds.count(PortKind.INTERNAL_ELEMENT) == 2
    assert len(ARCH.power_ports()) == 28

    groups = {c.mechanical_group for c in ARCH.connectors if c.mechanical_group}
    assert len(groups) == 3, "six series halves must be exactly three composite runtimes"


def test_internal_elements_couple_nothing() -> None:
    """An edge whose endpoints are one component is a constitutive element, not a port."""
    internal = [c for c in ARCH.connectors if ARCH.port_kind(c.name) is PortKind.INTERNAL_ELEMENT]
    assert {c.name for c in internal} == {"dorsal_arc_crosslink", "ecm_crosslink"}
    assert all(c.component_a == c.component_b for c in internal)


def test_reservoir_ports_are_the_ones_whose_far_side_does_not_move() -> None:
    """The discriminator is `dynamically_evolving`, NOT the ENVIRONMENT role.

    `ecm` is ENVIRONMENT and takes part in real couplings, so classifying by role would wrongly
    demote every ECM edge to a boundary condition.
    """
    reservoir = {c.name for c in ARCH.connectors if ARCH.port_kind(c.name) is PortKind.RESERVOIR_PORT}
    assert reservoir == {"ecm_far_field_anchor", "membrane_medium_traction"}

    ecm = ARCH.component("ecm")
    assert ecm.role is ComponentRole.ENVIRONMENT and ecm.dynamically_evolving
    assert "membrane_ecm_contact" not in reservoir


def test_power_ports_are_between_two_evolving_bodies() -> None:
    """The balance invariant is only meaningful where both sides own force and both sides move."""
    for name in ARCH.power_ports():
        c = next(x for x in ARCH.connectors if x.name == name)
        assert c.component_a != c.component_b
        assert c.mechanical_group is None
        assert ARCH.component(c.component_a).dynamically_evolving
        assert ARCH.component(c.component_b).dynamically_evolving


def test_an_ambiguous_declaration_raises_instead_of_being_classified_silently() -> None:
    """Two rules firing at once must fail loudly — silent precedence is how a miscount survives."""
    arch = CellArchitecture(
        components=(
            ComponentContract("a", ComponentRole.ACTIVE_LOAD_PATH, "rod", "graph"),
            ComponentContract(
                "frame", ComponentRole.ENVIRONMENT, "frame", "none",
                owns_geometry=False, dynamically_evolving=False,
            ),
        ),
        connectors=(
            ConnectorContract(
                "bad", ConnectorFamily.CONTACT, "a", "a", False, False,
                scope=ConnectorScope.INTERNAL, mechanical_group="some_series",
            ),
        ),
    )
    with pytest.raises(ValueError, match="inconsistent"):
        arch.port_kind("bad")


def test_an_unknown_name_raises() -> None:
    with pytest.raises(KeyError):
        ARCH.port_kind("no_such_connector")


def test_a_cut_isolates_its_connector_for_all_but_one_pair() -> None:
    """The balance gate's unit is a CUT. It names a connector only when the cut isolates one.

    27 of 28 power ports are their endpoints' sole coupling, so the cut IS the connector there. The
    one exception is deliberate and must stay visible rather than be attributed to either edge.
    """
    cuts = ARCH.balance_cuts()
    assert sum(len(v) for v in cuts.values()) == 28
    assert len(cuts) == 27

    shared = {k: v for k, v in cuts.items() if len(v) > 1}
    assert list(shared) == [("cortex", "membrane")]
    assert set(shared[("cortex", "membrane")]) == {"membrane_erm_cortex", "membrane_cortex_contact"}

    attributable = ARCH.attributable_cuts()
    assert len(attributable) == 26
    assert "membrane_erm_cortex" not in attributable.values()


def test_cuts_carry_only_power_ports() -> None:
    """A reservoir port has no peer body to cancel against; an internal element has no cut at all."""
    listed = {name for names in ARCH.balance_cuts().values() for name in names}
    for excluded in ("ecm_far_field_anchor", "membrane_medium_traction",
                     "dorsal_arc_crosslink", "ecm_crosslink", "fa_actin_anchor"):
        assert excluded not in listed
