"""Runtime registry for a component-first Active Cell actor.

The registry is intentionally agnostic to component physics.  It proves that runtime objects are
bound to declared ownership/connectivity contracts and deduplicates composite mechanical joints.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from aleph.engine.contracts import CellArchitecture, ConnectorContract

_TRANSACTION_HOOKS = ("snapshot_candidate", "rollback", "commit_irreversible")
_EVENT_HOOK = "propose_events"
_MECHANICS_HOOKS = (
    "accumulate",
    "accumulate_actor",
    "accumulate_boundary",
    "accumulate_candidate",
    "accumulate_internal",
    "accumulate_motor",
    "accumulate_pair",
    "accumulate_rig",
    "accumulate_transfer",
    "solve_candidate",
)


@dataclass(slots=True)
class CellActor:
    """Bind runtime components and connectors to a validated :class:`CellArchitecture`.

    Two semantic connector edges in one ``mechanical_group`` must bind to the same runtime
    joint object.  This is the build-time guard that prevents a geometry-less focal adhesion
    from becoming two springs in series by accident.
    """

    architecture: CellArchitecture
    _components: dict[str, object] = field(default_factory=dict, init=False)
    _connectors: dict[str, object] = field(default_factory=dict, init=False)

    def bind_component(self, name: str, runtime: object) -> None:
        """Bind one state-owning runtime object to a declared component name."""
        self.architecture.component(name)
        if name in self._components:
            raise ValueError(f"component {name!r} is already bound")
        self._components[name] = runtime

    def bind_connector(self, name: str, runtime: object) -> None:
        """Bind one runtime joint, enforcing composite mechanical-group identity."""
        contract = self.connector_contract(name)
        if name in self._connectors:
            raise ValueError(f"connector {name!r} is already bound")
        if contract.mechanical_group is not None:
            for peer in self.architecture.mechanical_group(contract.mechanical_group):
                existing = self._connectors.get(peer.name)
                if existing is not None and existing is not runtime:
                    raise ValueError(
                        f"mechanical group {contract.mechanical_group!r} must bind one composite runtime"
                    )
        self._connectors[name] = runtime

    def connector_contract(self, name: str) -> ConnectorContract:
        """Return a named connector declaration or raise ``KeyError``."""
        for connector in self.architecture.connectors:
            if connector.name == name:
                return connector
        raise KeyError(name)

    def component_runtime(self, name: str) -> object:
        """Return one bound component runtime."""
        return self._components[name]

    def connector_runtime(self, name: str) -> object:
        """Return one bound connector runtime."""
        return self._connectors[name]

    def missing_bindings(self) -> dict[str, tuple[str, ...]]:
        """Report declarations that have no runtime object yet."""
        components = tuple(
            component.name for component in self.architecture.components
            if component.name not in self._components
        )
        connectors = tuple(
            connector.name for connector in self.architecture.connectors
            if connector.name not in self._connectors
        )
        return {"components": components, "connectors": connectors}

    def assert_fully_bound(self) -> None:
        """Reject a production build with any name or runtime-capability hole."""
        missing = self.missing_bindings()
        if missing["components"] or missing["connectors"]:
            raise ValueError(
                "cell actor is not fully bound; "
                f"components={missing['components']}, connectors={missing['connectors']}"
            )
        for contract in self.architecture.components:
            runtime = self._components[contract.name]
            if contract.dynamically_evolving:
                self._require_hooks(
                    runtime,
                    (*_TRANSACTION_HOOKS, "accumulate_ledger"),
                    label=f"component {contract.name!r}",
                )
            if contract.has_events:
                self._require_hooks(
                    runtime, (_EVENT_HOOK,), label=f"component {contract.name!r}"
                )
            if contract.owns_geometry and not any(
                callable(getattr(runtime, hook, None)) for hook in _MECHANICS_HOOKS
            ):
                raise TypeError(
                    f"component {contract.name!r} has no registered mechanics/field hook"
                )
        for contract in self.architecture.connectors:
            runtime = self._connectors[contract.name]
            self._require_hooks(
                runtime,
                (*_TRANSACTION_HOOKS, "accumulate_ledger"),
                label=f"connector {contract.name!r}",
            )
            if contract.kinetics:
                self._require_hooks(
                    runtime, (_EVENT_HOOK,), label=f"connector {contract.name!r}"
                )
            if not any(callable(getattr(runtime, hook, None)) for hook in _MECHANICS_HOOKS):
                raise TypeError(f"connector {contract.name!r} has no mechanics/coupling hook")
            if contract.mechanical_group is not None:
                group = getattr(runtime, "mechanical_group", contract.mechanical_group)
                if group != contract.mechanical_group:
                    raise ValueError(
                        f"connector {contract.name!r} runtime has the wrong mechanical group"
                    )
                continue
            if getattr(runtime, "name", None) != contract.name:
                raise ValueError(f"connector {contract.name!r} runtime has the wrong name")
            endpoints = {
                getattr(runtime, "component_a", None),
                getattr(runtime, "component_b", None),
            }
            if endpoints != {contract.component_a, contract.component_b}:
                raise ValueError(f"connector {contract.name!r} runtime has the wrong endpoints")

    @staticmethod
    def _require_hooks(runtime: object, hooks: tuple[str, ...], *, label: str) -> None:
        missing = tuple(hook for hook in hooks if not callable(getattr(runtime, hook, None)))
        if missing:
            raise TypeError(f"{label} has an incomplete runtime API; missing {missing}")

    def unique_runtime_objects(self) -> tuple[object, ...]:
        """Return each component/joint object once, preserving deterministic registration order."""
        unique: list[object] = []
        identities: set[int] = set()
        for runtime in (*self._components.values(), *self._connectors.values()):
            identity = id(runtime)
            if identity not in identities:
                identities.add(identity)
                unique.append(runtime)
        return tuple(unique)

    def runtimes_with(self, method_name: str) -> tuple[Any, ...]:
        """Return unique runtime objects implementing a named engine hook."""
        return tuple(
            runtime for runtime in self.unique_runtime_objects()
            if callable(getattr(runtime, method_name, None))
        )
