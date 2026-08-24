"""Coverage and exact-once gates for the whole-cell candidate pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import pytest

from aleph.engine.contracts import CellArchitecture, reference_cell_architecture
from aleph.engine.dispatch import (
    CandidatePhase,
    CandidateTask,
    CellCandidatePipeline,
    canonical_facade_claims,
    validate_canonical_facade_claims,
)


@dataclass(slots=True)
class _RuntimeSpy:
    name: str
    calls: list[str] = field(default_factory=list)

    def run(self) -> None:
        self.calls.append(self.name)

    def run_with(self, marker: str, *, suffix: str) -> None:
        self.calls.append(f"{marker}{suffix}")


def _complete_tasks() -> tuple[CandidateTask, ...]:
    architecture = reference_cell_architecture()
    tasks: list[CandidateTask] = []
    for component in architecture.components:
        if component.owns_geometry:
            runtime = _RuntimeSpy(f"component:{component.name}")
            tasks.append(
                CandidateTask(
                    f"component:{component.name}",
                    CandidatePhase.ASSEMBLE_INTERNAL,
                    runtime,
                    "run",
                    component_claims=(component.name,),
                )
            )
    for connector in architecture.connectors:
        runtime = _RuntimeSpy(f"connector:{connector.name}")
        tasks.append(
            CandidateTask(
                f"connector:{connector.name}",
                CandidatePhase.ASSEMBLE_CONNECTORS,
                runtime,
                "run",
                connector_claims=(connector.name,),
            )
        )
    return tuple(tasks)


def test_reference_graph_requires_all_36_connector_dispatch_claims() -> None:
    architecture = reference_cell_architecture()
    pipeline = CellCandidatePipeline(architecture, _complete_tasks())

    # Deliberate tripwire: 35 -> 36 on PI D5-A (membrane_medium_traction, declared only).
    assert len(architecture.connectors) == 38
    assert {pipeline.connector_owner(edge.name) for edge in architecture.connectors} == {
        f"connector:{edge.name}" for edge in architecture.connectors
    }
    assert set(pipeline._component_owner) == {
        component.name for component in architecture.components if component.owns_geometry
    }


def test_real_canonical_facades_claim_all_common_edges_exactly_once() -> None:
    architecture = reference_cell_architecture()
    validate_canonical_facade_claims(architecture)
    claims = canonical_facade_claims()

    assert len(claims) == 8
    assert sum(len(claim.connector_claims) for claim in claims) == 38
    assert sum(len(claim.component_claims) for claim in claims) == 11


def test_missing_component_or_connector_claim_rejects_build() -> None:
    architecture = reference_cell_architecture()
    tasks = _complete_tasks()
    with pytest.raises(ValueError, match="missing_connectors"):
        CellCandidatePipeline(architecture, tasks[:-1])

    first_component = next(
        index for index, task in enumerate(tasks) if task.component_claims
    )
    without_component = (*tasks[:first_component], *tasks[first_component + 1 :])
    with pytest.raises(ValueError, match="missing_components"):
        CellCandidatePipeline(architecture, without_component)


def test_duplicate_and_unknown_claims_reject_build() -> None:
    architecture = reference_cell_architecture()
    tasks = _complete_tasks()
    connector_task = next(task for task in tasks if task.connector_claims)
    duplicate_runtime = _RuntimeSpy("duplicate")
    duplicate = CandidateTask(
        "duplicate-edge",
        CandidatePhase.ASSEMBLE_FIELDS,
        duplicate_runtime,
        "run",
        connector_claims=connector_task.connector_claims,
    )
    with pytest.raises(ValueError, match="claimed by both"):
        CellCandidatePipeline(architecture, (*tasks, duplicate))

    unknown = replace(
        connector_task,
        name="unknown-edge",
        runtime=_RuntimeSpy("unknown"),
        connector_claims=("not_in_the_graph",),
    )
    replaced = tuple(unknown if task is connector_task else task for task in tasks)
    with pytest.raises(ValueError, match="unknown connectors"):
        CellCandidatePipeline(architecture, replaced)


def test_same_runtime_cannot_dispatch_twice_in_one_phase() -> None:
    architecture = reference_cell_architecture()
    tasks = list(_complete_tasks())
    first = next(index for index, task in enumerate(tasks) if task.connector_claims)
    second = next(index for index in range(first + 1, len(tasks)) if tasks[index].connector_claims)
    tasks[second] = replace(tasks[second], runtime=tasks[first].runtime)
    with pytest.raises(ValueError, match="scheduled more than once"):
        CellCandidatePipeline(architecture, tuple(tasks))


def test_composite_semantic_edges_share_one_mechanics_call() -> None:
    architecture = reference_cell_architecture()
    tasks = list(_complete_tasks())
    group = architecture.mechanical_group("alpha2beta1_collagen_series")
    group_names = {connector.name for connector in group}
    tasks = [task for task in tasks if not (set(task.connector_claims) & group_names)]
    runtime = _RuntimeSpy("composite")
    tasks.append(
        CandidateTask(
            "fa-series-composite",
            CandidatePhase.ASSEMBLE_CONNECTORS,
            runtime,
            "run",
            connector_claims=tuple(sorted(group_names)),
        )
    )
    pipeline = CellCandidatePipeline(architecture, tuple(tasks))
    pipeline.run_phase(CandidatePhase.ASSEMBLE_CONNECTORS)

    assert runtime.calls == ["composite"]
    assert {pipeline.connector_owner(name) for name in group_names} == {"fa-series-composite"}


def test_phase_order_and_argument_forwarding_are_deterministic() -> None:
    architecture = reference_cell_architecture()
    tasks = list(_complete_tasks())
    connector_index = next(index for index, task in enumerate(tasks) if task.connector_claims)
    component_index = next(index for index, task in enumerate(tasks) if task.component_claims)
    order: list[str] = []

    @dataclass(slots=True)
    class _OrderRuntime:
        label: str

        def run(self) -> None:
            order.append(self.label)

        def run_with(self, marker: str, *, suffix: str) -> None:
            order.append(f"{marker}{suffix}")

    tasks[connector_index] = replace(
        tasks[connector_index],
        phase=CandidatePhase.DISCOVER_GEOMETRY,
        runtime=_OrderRuntime("discover"),
    )
    tasks[component_index] = replace(
        tasks[component_index],
        phase=CandidatePhase.SOLVE_COUPLED,
        runtime=_OrderRuntime("unused"),
        method_name="run_with",
        args=("solve",),
        kwargs=(("suffix", "-done"),),
    )
    pipeline = CellCandidatePipeline(architecture, tuple(tasks))
    pipeline.run_candidate()

    assert order == ["discover", "solve-done"]


def test_task_rejects_missing_method_and_repeated_keyword() -> None:
    runtime = _RuntimeSpy("x")
    with pytest.raises(TypeError, match="no callable"):
        CandidateTask("bad", CandidatePhase.ASSEMBLE_INTERNAL, runtime, "missing")
    with pytest.raises(ValueError, match="repeats a keyword"):
        CandidateTask(
            "bad-kwargs",
            CandidatePhase.ASSEMBLE_INTERNAL,
            runtime,
            "run_with",
            kwargs=(("suffix", "a"), ("suffix", "b")),
        )


# --- the connector runtime census: what makes "38/38" a measurement -------------------------------


def test_every_declared_connector_names_a_runtime_exactly_once() -> None:
    """The census must cover the architecture with no hole and no stale entry, in either direction."""
    from aleph.engine.dispatch import connector_runtime_census, validate_connector_runtime_census

    architecture = reference_cell_architecture()
    validate_connector_runtime_census(architecture)
    census = connector_runtime_census()
    assert len(census) == len(architecture.connectors)
    assert len({binding.connector for binding in census}) == len(census)


def test_every_named_runtime_symbol_imports_and_is_not_a_bare_protocol() -> None:
    """A Protocol is a SEAM, not a runtime, and a census that cannot tell them apart reports a lie.

    This gate is the reason the census exists rather than a hand count.  Writing it immediately caught
    ``ecm_world.CompositeECMClutch`` — the mature FA series' declared interface — standing in for its
    actual implementation, ``CompositeLoadPathECMClutchAdapter``.  A hand count would have recorded the
    edge as served either way.
    """
    import importlib

    from aleph.engine.dispatch import connector_runtime_census

    for binding in connector_runtime_census():
        module = importlib.import_module(binding.module)
        runtime = getattr(module, binding.symbol, None)
        assert runtime is not None, f"{binding.connector}: {binding.module}.{binding.symbol} is gone"
        assert isinstance(runtime, type), f"{binding.connector}: {binding.symbol} is not a class"
        assert not getattr(runtime, "_is_protocol", False), (
            f"{binding.connector}: {binding.symbol} is a Protocol, i.e. a declared seam with no "
            "implementation behind it — the census may not count it as a runtime"
        )


def test_every_named_runtime_carries_the_transaction_and_mechanics_api() -> None:
    """A connector runtime the world transaction cannot drive is not a runtime.

    ``CellActor.assert_fully_bound`` enforces exactly this at build; asserting it statically means the
    census cannot claim an edge is served by a class that would be rejected the moment it was bound.
    """
    import importlib

    from aleph.engine.dispatch import connector_runtime_census

    mechanics_hooks = (
        "accumulate", "accumulate_actor", "accumulate_boundary", "accumulate_candidate",
        "accumulate_contact", "accumulate_internal", "accumulate_motor", "accumulate_pair",
        "accumulate_rig", "accumulate_transfer", "solve_candidate",
    )
    for binding in connector_runtime_census():
        runtime = getattr(importlib.import_module(binding.module), binding.symbol)
        for hook in ("snapshot_candidate", "rollback", "commit_irreversible", "accumulate_ledger"):
            assert callable(getattr(runtime, hook, None)), f"{binding.symbol} is missing {hook}"
        assert any(callable(getattr(runtime, hook, None)) for hook in mechanics_hooks), (
            f"{binding.symbol} has no mechanics/coupling hook"
        )


def test_a_census_hole_is_reported_as_a_hole() -> None:
    """A negative control: the coverage check must be able to fail."""
    from aleph.engine.contracts import ConnectorContract, ConnectorFamily
    from aleph.engine.dispatch import validate_connector_runtime_census

    architecture = reference_cell_architecture()
    extended = CellArchitecture(
        components=architecture.components,
        connectors=(
            *architecture.connectors,
            ConnectorContract(
                "cortex_ecm_invented", ConnectorFamily.CONTACT, "cortex", "ecm", False, False,
            ),
        ),
    )
    with pytest.raises(ValueError, match="without_a_runtime=\\['cortex_ecm_invented'\\]"):
        validate_connector_runtime_census(extended)
