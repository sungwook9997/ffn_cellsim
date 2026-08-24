"""Runtime-binding gates for the component-first Cell actor."""

from __future__ import annotations

import pytest

from aleph.engine import CellActor, reference_cell_architecture


class _Runtime:
    def accumulate(self) -> None:
        pass


def test_actor_rejects_undeclared_or_duplicate_component() -> None:
    actor = CellActor(reference_cell_architecture())
    actor.bind_component("membrane", object())
    with pytest.raises(ValueError, match="already bound"):
        actor.bind_component("membrane", object())
    with pytest.raises(KeyError):
        actor.bind_component("unknown", object())


def test_geometryless_fa_semantic_edges_share_one_composite_runtime() -> None:
    actor = CellActor(reference_cell_architecture())
    composite_clutch = _Runtime()
    actor.bind_connector("fa_actin_anchor", composite_clutch)
    actor.bind_connector("integrin_collagen_clutch", composite_clutch)
    assert actor.connector_runtime("fa_actin_anchor") is actor.connector_runtime(
        "integrin_collagen_clutch"
    )
    assert actor.runtimes_with("accumulate") == (composite_clutch,)


def test_double_stiffness_from_two_fa_runtime_objects_is_rejected() -> None:
    actor = CellActor(reference_cell_architecture())
    actor.bind_connector("fa_actin_anchor", _Runtime())
    with pytest.raises(ValueError, match="one composite runtime"):
        actor.bind_connector("integrin_collagen_clutch", _Runtime())


def test_missing_bindings_are_explicit() -> None:
    actor = CellActor(reference_cell_architecture())
    actor.bind_component("membrane", object())
    missing = actor.missing_bindings()
    assert "membrane" not in missing["components"]
    assert "cortex" in missing["components"]
    assert "membrane_erm_cortex" in missing["connectors"]
    with pytest.raises(ValueError, match="not fully bound"):
        actor.assert_fully_bound()


def test_complete_names_with_inert_objects_are_not_a_production_cell() -> None:
    architecture = reference_cell_architecture()
    actor = CellActor(architecture)
    inert = object()
    for component in architecture.components:
        actor.bind_component(component.name, inert)
    for connector in architecture.connectors:
        actor.bind_connector(connector.name, inert)
    assert actor.missing_bindings() == {"components": (), "connectors": ()}
    with pytest.raises(TypeError, match="incomplete runtime API"):
        actor.assert_fully_bound()
