"""Endpoint/owner validator for manifests and populations against the CANONICAL engine grammar.

The virtual-cell layer can now name components and connectors freely
(:class:`~aleph.virtual_cell.contracts.CellStateManifest`), and it can lay out a heterotypic
population (:mod:`aleph.virtual_cell.population`).  Neither of those checks the names against the
architecture the engine actually declares, so a manifest may currently reference a component that does
not exist, wire a connector to the wrong endpoint, or -- the project's oldest defect -- assume that two
things sitting in the same array are mechanically coupled.  This module is that check.

The census is READ FROM CODE, never transcribed.  ``aleph/engine/contracts.py`` declares the
components and connectors; ``aleph/engine/dispatch.py`` declares, in
``canonical_facade_claims()``, which single dispatch facade OWNS each component and connector edge.
Both are loaded here without importing the ``aleph.engine`` package, because that package's
``__init__`` imports the fluid/NMII/surface rigs and therefore initialises Warp: ``contracts.py`` is a
leaf module and is executed from its file path, and the ownership manifest is extracted from
``dispatch.py`` with :mod:`ast` (its function body imports the eight facade classes, which would pull
Warp in again).  When either source cannot be located or parsed, :class:`CensusUnavailableError` is
raised and nothing is validated -- there is deliberately NO fallback list, because a list copied out of
a status document is exactly the drift this module exists to catch.  A caller with its own census may
inject one through :meth:`EngineGrammar.from_declarations`.

**The engine's ownership rule, as measured rather than assumed.**  ``ConnectorContract`` carries no
owner field at all; ownership is the exact-once dispatch claim enforced by
``dispatch.validate_canonical_facade_claims``.  A facade owns a set of components and a set of
connectors, and the expected relation -- the owning facade also owns one of the connector's endpoint
components -- holds for 35 of the 36 declared connectors.  The exception is ``fa_actin_anchor``
(``sf_arc`` <-> ``focal_adhesion``), owned by ``ECMWorld``, which owns neither endpoint.  That is not a
defect: ``focal_adhesion`` is declared ``owns_geometry=False`` so it has no facade of its own, and
``fa_actin_anchor`` shares the ``alpha2beta1_collagen_series`` mechanical group with
``integrin_collagen_clutch``, so the composite series joint is one runtime object owned once by the
facade that owns the ECM end.  The rule encoded here is therefore the measured one -- *the owner must
own an endpoint UNLESS the connector belongs to a mechanical group* -- and the exemption is reported as
an ``EXPANSION`` finding rather than silently dropped, so a future census change cannot hide behind it.

Nothing here runs physics, allocates device memory, or grants an evidence rung.  A clean report means
the wiring is expressible in the declared grammar; it says nothing about whether the physics is right.

Sanity Gate:
  * dimensional: this module carries no physical units.  The only quantities it compares are
    co-location extents in [um], and they are compared against each other, never against a tolerance
    chosen to make a check pass -- an overlap is an overlap at any scale.
  * boundary: an empty component set, a self-referential connector instance, a duplicate cell ID, and
    a connector whose two endpoints resolve to the same reference are each rejected explicitly rather
    than passing through as a degenerate no-op.
  * conservation-invariant: the invariant enforced is topological, not energetic -- every declared
    connector edge has exactly one owner and exactly two endpoints, and a composite mechanical group
    resolves to ONE owner, so a series joint can never be counted as two springs.
  * measurement-protocol: co-location is measured as shared array identity OR overlapping axis-aligned
    extent, and it is never evidence of coupling.  Two components that overlap and have no declared
    connector are reported ``CO_LOCATED_WITHOUT_CONNECTOR``; the check has no threshold to tune.
  * numerical: extents are compared with closed-interval arithmetic on finite floats only; a
    non-finite bound is rejected at construction, so no comparison can silently succeed on a NaN.
"""

from __future__ import annotations

import ast
import importlib.util
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from types import MappingProxyType

__all__ = [
    "CellJunction",
    "CellWiring",
    "CensusUnavailableError",
    "ComponentDeclaration",
    "ComponentExtent",
    "ConnectorDeclaration",
    "ConnectorInstance",
    "EngineGrammar",
    "GrammarReport",
    "GrammarViolation",
    "NAMESPACE_SEPARATOR",
    "PopulationWiring",
    "SHARED_NAMESPACE",
    "Severity",
    "SharedComponent",
    "ViolationCode",
    "check_colocation",
    "load_engine_grammar",
    "namespaced_ref",
    "parse_ref",
    "shared_ref",
    "validate_manifest",
    "validate_population",
]

#: Separator between a cell ID and the canonical component name it owns.
NAMESPACE_SEPARATOR = "::"

#: Reserved namespace for an object that is referenced by many cells rather than owned by one.
SHARED_NAMESPACE = "shared"

_ENGINE_CONTRACTS_MODULE = "_ffn_virtual_cell_engine_census_contracts"


class CensusUnavailableError(RuntimeError):
    """The canonical component/connector census could not be read from the engine source."""


def _nonempty(value: object, *, what: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{what} must be a non-empty string")
    return value.strip()


def _unique(values: Sequence[str], *, what: str) -> tuple[str, ...]:
    if isinstance(values, str) or not isinstance(values, Sequence):
        raise TypeError(f"{what} must be a sequence of strings")
    normalized = tuple(_nonempty(value, what=f"{what} entry") for value in values)
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{what} must not contain duplicates")
    return normalized


class Severity(StrEnum):
    """Whether a finding refuses the wiring or opens a programme expansion.

    ``REFUSED`` means the wiring is not expressible in the declared grammar and must not be built.
    ``EXPANSION`` means the wiring is structurally consistent but names something the canonical
    architecture does not declare -- a missing component or missing law -- which is a project to open,
    not a bug to patch.  Keeping them apart is what stops "the grammar has no cell-cell junction"
    from being reported with the same weight as "this connector has no such endpoint".
    """

    REFUSED = "REFUSED"
    EXPANSION = "EXPANSION"


class ViolationCode(StrEnum):
    """Typed grammar findings; the string value is what an artifact records."""

    UNKNOWN_COMPONENT = "unknown-component"
    UNKNOWN_CONNECTOR = "unknown-connector"
    ENDPOINT_NOT_PRESENT = "endpoint-not-present"
    ENDPOINT_MISMATCH = "endpoint-mismatch"
    OWNER_MISSING = "owner-missing"
    OWNER_NOT_ENDPOINT = "owner-not-endpoint"
    OWNER_NOT_ENDPOINT_COMPOSITE = "owner-not-endpoint-composite"
    ONE_SIDED_CONNECTOR = "one-sided-connector"
    ADJOINT_NOT_DECLARED = "adjoint-not-declared"
    UNCONNECTED_COMPONENT = "unconnected-component"
    CO_LOCATED_WITHOUT_CONNECTOR = "co-located-without-connector"
    COMPONENT_NOT_NAMESPACED = "component-not-namespaced"
    DUPLICATE_CELL_ID = "duplicate-cell-id"
    COMPONENT_NAME_COLLISION = "component-name-collision"
    CROSS_CELL_INTRACELLULAR_CONNECTOR = "cross-cell-intracellular-connector"
    PER_CELL_COMPONENT_SHARED = "per-cell-component-shared"
    SHARED_COMPONENT_DUPLICATED = "shared-component-duplicated"
    SHARED_COMPONENT_UNDECLARED = "shared-component-undeclared"
    JUNCTION_LAW_NOT_IN_GRAMMAR = "junction-law-not-in-grammar"
    JUNCTION_ENDPOINTS_SAME_CELL = "junction-endpoints-same-cell"


#: Which programme expansion each finding routes to (shared lane brief, "failure routing").
_EXPANSION: Mapping[ViolationCode, str] = MappingProxyType(
    {
        ViolationCode.UNKNOWN_COMPONENT: "archetype failure -> missing component / missing law project",
        ViolationCode.UNKNOWN_CONNECTOR: "archetype failure -> missing component / missing law project",
        ViolationCode.ENDPOINT_NOT_PRESENT: "wiring failure -> component set is incomplete for this archetype",
        ViolationCode.ENDPOINT_MISMATCH: "wiring failure -> connector instance contradicts the declared edge",
        ViolationCode.OWNER_MISSING: "wiring failure -> exact-once dispatch ownership must be assigned",
        ViolationCode.OWNER_NOT_ENDPOINT: "wiring failure -> exact-once dispatch ownership must be assigned",
        ViolationCode.OWNER_NOT_ENDPOINT_COMPOSITE: "census note -> composite mechanical group owns the edge",
        ViolationCode.ONE_SIDED_CONNECTOR: "wiring failure -> Newton's 3rd law requires a bidirectional edge",
        ViolationCode.ADJOINT_NOT_DECLARED: "wiring failure -> adjoint closure must be declared per connector",
        ViolationCode.UNCONNECTED_COMPONENT: "wiring failure -> declare a connector or drop the component",
        ViolationCode.CO_LOCATED_WITHOUT_CONNECTOR: "wiring failure -> declare a connector; co-location is never connection",
        ViolationCode.COMPONENT_NOT_NAMESPACED: "wiring failure -> namespace per-cell state by cell ID",
        ViolationCode.DUPLICATE_CELL_ID: "wiring failure -> cell identity must be unique in a population",
        ViolationCode.COMPONENT_NAME_COLLISION: "wiring failure -> namespace per-cell state by cell ID",
        ViolationCode.CROSS_CELL_INTRACELLULAR_CONNECTOR: "archetype failure -> missing component / missing law project",
        ViolationCode.PER_CELL_COMPONENT_SHARED: "wiring failure -> per-cell state cannot be one shared object",
        ViolationCode.SHARED_COMPONENT_DUPLICATED: "wiring failure -> one environment object referenced by many cells",
        ViolationCode.SHARED_COMPONENT_UNDECLARED: "wiring failure -> declare the shared environment object",
        ViolationCode.JUNCTION_LAW_NOT_IN_GRAMMAR: "archetype failure -> missing component / missing law project",
        ViolationCode.JUNCTION_ENDPOINTS_SAME_CELL: "wiring failure -> a junction joins two distinct cells",
    }
)


@dataclass(frozen=True, slots=True)
class GrammarViolation:
    """One typed finding, its subject, and the expansion it opens."""

    code: ViolationCode
    subject: str
    detail: str
    severity: Severity = Severity.REFUSED

    def __post_init__(self) -> None:
        if not isinstance(self.code, ViolationCode):
            raise TypeError("code must be a ViolationCode")
        if not isinstance(self.severity, Severity):
            raise TypeError("severity must be a Severity")
        object.__setattr__(self, "subject", _nonempty(self.subject, what="violation subject"))
        object.__setattr__(self, "detail", _nonempty(self.detail, what="violation detail"))

    @property
    def expansion(self) -> str:
        """Return the programme expansion this finding routes to."""
        return _EXPANSION[self.code]

    def to_dict(self) -> dict[str, str]:
        """Return a JSON-able record of this finding."""
        return {
            "code": self.code.value,
            "subject": self.subject,
            "detail": self.detail,
            "severity": self.severity.value,
            "expansion": self.expansion,
        }


@dataclass(frozen=True, slots=True)
class GrammarReport:
    """The findings of one validation pass plus the checks that produced them."""

    subject_id: str
    findings: tuple[GrammarViolation, ...]
    checks_run: tuple[str, ...]
    census_source: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "subject_id", _nonempty(self.subject_id, what="subject_id"))
        findings = tuple(self.findings)
        if any(not isinstance(finding, GrammarViolation) for finding in findings):
            raise TypeError("findings must contain GrammarViolation values")
        object.__setattr__(self, "findings", findings)
        object.__setattr__(self, "checks_run", _unique(self.checks_run, what="checks_run"))
        object.__setattr__(
            self, "census_source", _nonempty(self.census_source, what="census_source")
        )

    @property
    def refusals(self) -> tuple[GrammarViolation, ...]:
        """Return only the findings that refuse the wiring."""
        return tuple(f for f in self.findings if f.severity is Severity.REFUSED)

    @property
    def expansions(self) -> tuple[GrammarViolation, ...]:
        """Return only the findings that open a programme expansion."""
        return tuple(f for f in self.findings if f.severity is Severity.EXPANSION)

    @property
    def ok(self) -> bool:
        """True when nothing refuses the wiring; expansion findings do not make it false."""
        return not self.refusals

    @property
    def codes(self) -> frozenset[ViolationCode]:
        """Return the distinct finding codes present."""
        return frozenset(finding.code for finding in self.findings)

    def raise_for_refusals(self) -> None:
        """Raise ``ValueError`` naming every refusal, or return when there are none."""
        if self.ok:
            return
        detail = "; ".join(f"{f.code.value}[{f.subject}]: {f.detail}" for f in self.refusals)
        raise ValueError(f"{self.subject_id} is not expressible in the canonical grammar: {detail}")

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-able record of this report."""
        return {
            "subject_id": self.subject_id,
            "census_source": self.census_source,
            "checks_run": list(self.checks_run),
            "ok": self.ok,
            "findings": [finding.to_dict() for finding in self.findings],
        }


@dataclass(frozen=True, slots=True)
class ComponentDeclaration:
    """One canonical state-owning component, as declared by the engine."""

    name: str
    role: str
    representation: str
    solver: str
    owns_geometry: bool
    dynamically_evolving: bool

    def __post_init__(self) -> None:
        for name in ("name", "role", "representation", "solver"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        for name in ("owns_geometry", "dynamically_evolving"):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be bool")


@dataclass(frozen=True, slots=True)
class ConnectorDeclaration:
    """One canonical connector edge, its endpoints, and the facade that owns its dispatch."""

    name: str
    family: str
    endpoint_a: str
    endpoint_b: str
    scope: str
    kinetics: bool
    commit_on_accept: bool
    bidirectional: bool
    adjoint_transfer_required: bool
    mechanical_group: str | None
    owner: str

    def __post_init__(self) -> None:
        for name in ("name", "family", "endpoint_a", "endpoint_b", "scope", "owner"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        for name in (
            "kinetics",
            "commit_on_accept",
            "bidirectional",
            "adjoint_transfer_required",
        ):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be bool")
        if self.mechanical_group is not None:
            object.__setattr__(
                self, "mechanical_group", _nonempty(self.mechanical_group, what="mechanical_group")
            )

    @property
    def endpoints(self) -> frozenset[str]:
        """Return the connector's two declared endpoint component names."""
        return frozenset({self.endpoint_a, self.endpoint_b})

    @property
    def internal(self) -> bool:
        """True when both endpoints are the same component (an intra-component edge)."""
        return self.endpoint_a == self.endpoint_b


@dataclass(frozen=True, slots=True)
class EngineGrammar:
    """The canonical component/connector census plus per-edge dispatch ownership."""

    components: tuple[ComponentDeclaration, ...]
    connectors: tuple[ConnectorDeclaration, ...]
    facade_components: Mapping[str, tuple[str, ...]]
    source: str

    def __post_init__(self) -> None:
        components = tuple(self.components)
        if not components or any(
            not isinstance(component, ComponentDeclaration) for component in components
        ):
            raise ValueError("grammar needs at least one ComponentDeclaration")
        names = tuple(component.name for component in components)
        if len(set(names)) != len(names):
            raise ValueError("component names must be unique")
        object.__setattr__(self, "components", components)

        connectors = tuple(self.connectors)
        if any(not isinstance(connector, ConnectorDeclaration) for connector in connectors):
            raise TypeError("connectors must contain ConnectorDeclaration values")
        connector_names = tuple(connector.name for connector in connectors)
        if len(set(connector_names)) != len(connector_names):
            raise ValueError("connector names must be unique")
        known = set(names)
        for connector in connectors:
            missing = connector.endpoints - known
            if missing:
                raise ValueError(
                    f"connector {connector.name!r} references undeclared components "
                    f"{sorted(missing)}"
                )
        object.__setattr__(self, "connectors", connectors)

        if not isinstance(self.facade_components, Mapping):
            raise TypeError("facade_components must be a mapping")
        frozen = {
            _nonempty(facade, what="facade name"): _unique(
                tuple(claims), what=f"facade {facade!r} component claims"
            )
            for facade, claims in self.facade_components.items()
        }
        object.__setattr__(self, "facade_components", MappingProxyType(frozen))
        object.__setattr__(self, "source", _nonempty(self.source, what="source"))

    @classmethod
    def from_declarations(
        cls,
        components: Sequence[ComponentDeclaration],
        connectors: Sequence[ConnectorDeclaration],
        facade_components: Mapping[str, Sequence[str]],
        *,
        source: str,
    ) -> EngineGrammar:
        """Build a grammar from an explicitly supplied census.

        This is the escape hatch for a caller that has its own authoritative census (a different
        engine, or a test fixture).  It is deliberately explicit: there is no default census, so no
        code path can quietly substitute a transcribed list for the one the engine declares.

        Args:
            components: The declared state-owning components.
            connectors: The declared connector edges, each already carrying its owner.
            facade_components: Owner name -> the component names that owner claims.
            source: A human-readable provenance string recorded on every report.

        Returns:
            The validated :class:`EngineGrammar`.
        """
        return cls(
            components=tuple(components),
            connectors=tuple(connectors),
            facade_components={key: tuple(value) for key, value in facade_components.items()},
            source=source,
        )

    @property
    def component_names(self) -> frozenset[str]:
        """Return every declared component name."""
        return frozenset(component.name for component in self.components)

    @property
    def connector_names(self) -> frozenset[str]:
        """Return every declared connector name."""
        return frozenset(connector.name for connector in self.connectors)

    @property
    def environment_components(self) -> frozenset[str]:
        """Return components whose role is ``environment`` -- the shared, never-per-cell ones."""
        return frozenset(
            component.name for component in self.components if component.role == "environment"
        )

    def component(self, name: str) -> ComponentDeclaration:
        """Return one declared component or raise ``KeyError``."""
        for component in self.components:
            if component.name == name:
                return component
        raise KeyError(name)

    def connector(self, name: str) -> ConnectorDeclaration:
        """Return one declared connector or raise ``KeyError``."""
        for connector in self.connectors:
            if connector.name == name:
                return connector
        raise KeyError(name)

    def connectors_touching(self, component_name: str) -> tuple[ConnectorDeclaration, ...]:
        """Return every declared connector with ``component_name`` at either endpoint."""
        return tuple(
            connector for connector in self.connectors if component_name in connector.endpoints
        )

    def induced_connectors(
        self, component_names: Iterable[str]
    ) -> tuple[ConnectorDeclaration, ...]:
        """Return the connectors whose BOTH endpoints lie inside ``component_names``.

        Args:
            component_names: The component subset an archetype declares.

        Returns:
            The induced subgraph's connectors, in declaration order.

        Raises:
            KeyError: If any name is not declared.
        """
        selected = {_nonempty(name, what="component name") for name in component_names}
        for name in selected:
            self.component(name)
        return tuple(connector for connector in self.connectors if connector.endpoints <= selected)

    def isolated_components(self, component_names: Iterable[str]) -> tuple[str, ...]:
        """Return components in the subset that no induced connector touches.

        Co-location in a component list is never a connection, so a component that survives the
        induced subgraph with no edge is UNCONNECTED and is reported as such.
        """
        selected = {_nonempty(name, what="component name") for name in component_names}
        induced = self.induced_connectors(selected)
        touched: set[str] = set()
        for connector in induced:
            touched |= connector.endpoints
        return tuple(sorted(selected - touched))

    def owner_of(self, connector_name: str) -> str:
        """Return the single facade that owns one connector's dispatch."""
        return self.connector(connector_name).owner


def _load_engine_contracts_module(path: Path) -> object:
    """Execute the engine's leaf ``contracts.py`` from its file path, without importing its package.

    ``aleph.engine.__init__`` imports the fluid, NMII and surface rigs, which import Warp.  The
    virtual-cell layer must stay CPU-only and Warp-free at import time, and ``contracts.py`` itself
    imports nothing but ``dataclasses`` and ``enum``, so executing it directly is both safe and
    faithful -- the census comes from running the engine's own declaration, not from a copy of it.
    """
    cached = sys.modules.get(_ENGINE_CONTRACTS_MODULE)
    if cached is not None and getattr(cached, "__file__", None) == str(path):
        return cached
    spec = importlib.util.spec_from_file_location(_ENGINE_CONTRACTS_MODULE, path)
    if spec is None or spec.loader is None:
        raise CensusUnavailableError(f"cannot load the engine component census from {path}")
    module = importlib.util.module_from_spec(spec)
    # Registered before execution: ``dataclasses(slots=True)`` resolves the defining module by name.
    sys.modules[_ENGINE_CONTRACTS_MODULE] = module
    try:
        spec.loader.exec_module(module)
    except Exception as error:  # pragma: no cover - only on a broken engine tree
        sys.modules.pop(_ENGINE_CONTRACTS_MODULE, None)
        raise CensusUnavailableError(
            f"engine census at {path} failed to execute: {error}"
        ) from error
    return module


def _extract_facade_claims(path: Path) -> tuple[dict[str, tuple[str, ...]], dict[str, str]]:
    """Extract ``(facade -> components, connector -> facade)`` from ``dispatch.py`` with :mod:`ast`.

    ``canonical_facade_claims()`` imports the eight facade classes inside its body, so calling it
    would initialise Warp.  The manifest is a literal tuple of ``CanonicalFacadeClaim(...)`` calls, so
    it is read statically instead.  Anything that does not parse into that shape raises
    :class:`CensusUnavailableError` -- an unreadable ownership manifest must stop validation, never
    fall back to an assumed one.
    """
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except (OSError, SyntaxError) as error:
        raise CensusUnavailableError(
            f"cannot parse the engine dispatch manifest at {path}"
        ) from error

    function = next(
        (
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef) and node.name == "canonical_facade_claims"
        ),
        None,
    )
    if function is None:
        raise CensusUnavailableError(
            f"{path} declares no canonical_facade_claims(); connector ownership cannot be read"
        )
    returned = next(
        (node.value for node in ast.walk(function) if isinstance(node, ast.Return)), None
    )
    if not isinstance(returned, ast.Tuple):
        raise CensusUnavailableError(
            "canonical_facade_claims() no longer returns a literal tuple of claims"
        )

    facade_components: dict[str, tuple[str, ...]] = {}
    connector_owner: dict[str, str] = {}
    for element in returned.elts:
        if not isinstance(element, ast.Call):
            raise CensusUnavailableError("canonical_facade_claims() holds a non-call element")
        arguments = {index: value for index, value in enumerate(element.args)} | {
            keyword.arg: keyword.value for keyword in element.keywords if keyword.arg
        }
        owner_node = arguments.get(0, arguments.get("owner_type"))
        components_node = arguments.get(2, arguments.get("component_claims"))
        connectors_node = arguments.get(3, arguments.get("connector_claims"))
        if not isinstance(owner_node, ast.Name):
            raise CensusUnavailableError("a facade claim has no statically readable owner type")
        try:
            components = tuple(ast.literal_eval(components_node))
            connectors = tuple(ast.literal_eval(connectors_node))
        except (ValueError, TypeError) as error:
            raise CensusUnavailableError(
                f"facade {owner_node.id} has non-literal component/connector claims"
            ) from error
        if owner_node.id in facade_components:
            raise CensusUnavailableError(f"facade {owner_node.id} is claimed twice")
        facade_components[owner_node.id] = components
        for connector in connectors:
            if connector in connector_owner:
                raise CensusUnavailableError(
                    f"connector {connector!r} has duplicate canonical facade owners"
                )
            connector_owner[connector] = owner_node.id
    if not facade_components or not connector_owner:
        raise CensusUnavailableError("the engine dispatch manifest is empty")
    return facade_components, connector_owner


def load_engine_grammar(
    *,
    contracts_path: Path | str | None = None,
    dispatch_path: Path | str | None = None,
) -> EngineGrammar:
    """Read the canonical component/connector grammar out of the engine source.

    Args:
        contracts_path: Override for ``aleph/engine/contracts.py``.
        dispatch_path: Override for ``aleph/engine/dispatch.py``.

    Returns:
        The :class:`EngineGrammar` the engine actually declares, with per-connector ownership.

    Raises:
        CensusUnavailableError: If either source is missing, unparsable, or no longer exposes the
            declaration this reader depends on.  There is no fallback census by design.
    """
    engine_root = Path(__file__).resolve().parents[1] / "engine"
    contracts = Path(contracts_path) if contracts_path is not None else engine_root / "contracts.py"
    dispatch = Path(dispatch_path) if dispatch_path is not None else engine_root / "dispatch.py"
    for path in (contracts, dispatch):
        if not path.is_file():
            raise CensusUnavailableError(f"canonical engine census source not found: {path}")

    module = _load_engine_contracts_module(contracts)
    builder = getattr(module, "reference_cell_architecture", None)
    if not callable(builder):
        raise CensusUnavailableError(
            f"{contracts} exposes no reference_cell_architecture(); the census cannot be read"
        )
    architecture = builder()
    facade_components, connector_owner = _extract_facade_claims(dispatch)

    components = tuple(
        ComponentDeclaration(
            name=contract.name,
            role=str(contract.role),
            representation=contract.representation,
            solver=contract.solver,
            owns_geometry=bool(contract.owns_geometry),
            dynamically_evolving=bool(contract.dynamically_evolving),
        )
        for contract in architecture.components
    )
    unowned = sorted({contract.name for contract in architecture.connectors} - set(connector_owner))
    if unowned:
        raise CensusUnavailableError(
            "the dispatch manifest does not own every declared connector: " + ", ".join(unowned)
        )
    connectors = tuple(
        ConnectorDeclaration(
            name=contract.name,
            family=str(contract.family),
            endpoint_a=contract.component_a,
            endpoint_b=contract.component_b,
            scope=str(contract.scope),
            kinetics=bool(contract.kinetics),
            commit_on_accept=bool(contract.commit_on_accept),
            bidirectional=bool(contract.bidirectional),
            adjoint_transfer_required=bool(contract.adjoint_transfer_required),
            mechanical_group=contract.mechanical_group,
            owner=connector_owner[contract.name],
        )
        for contract in architecture.connectors
    )
    source = f"{contracts.as_posix()} + {dispatch.as_posix()}"
    return EngineGrammar(
        components=components,
        connectors=connectors,
        facade_components=facade_components,
        source=source,
    )


def namespaced_ref(cell_id: str, component_name: str) -> str:
    """Return the population-unique reference for one cell's own component state."""
    cell = _nonempty(cell_id, what="cell_id")
    component = _nonempty(component_name, what="component_name")
    if NAMESPACE_SEPARATOR in cell or NAMESPACE_SEPARATOR in component:
        raise ValueError(f"cell_id and component_name must not contain {NAMESPACE_SEPARATOR!r}")
    if cell == SHARED_NAMESPACE:
        raise ValueError(f"{SHARED_NAMESPACE!r} is reserved for shared objects")
    return f"{cell}{NAMESPACE_SEPARATOR}{component}"


def shared_ref(object_id: str) -> str:
    """Return the reference for ONE object that many cells reference rather than own."""
    identity = _nonempty(object_id, what="object_id")
    if NAMESPACE_SEPARATOR in identity:
        raise ValueError(f"object_id must not contain {NAMESPACE_SEPARATOR!r}")
    return f"{SHARED_NAMESPACE}{NAMESPACE_SEPARATOR}{identity}"


def parse_ref(reference: str) -> tuple[str | None, str]:
    """Split a reference into ``(cell_id_or_None, name)``.

    Args:
        reference: A reference produced by :func:`namespaced_ref` or :func:`shared_ref`, or a bare
            component name (which is exactly the un-namespaced case the validator refuses).

    Returns:
        ``(cell_id, component_name)`` for a per-cell reference, ``(None, object_id)`` for a shared
        reference, and ``(None, name)`` for a bare name -- callers distinguish the last two by
        testing for the separator themselves.
    """
    value = _nonempty(reference, what="reference")
    if NAMESPACE_SEPARATOR not in value:
        return None, value
    owner, _, name = value.partition(NAMESPACE_SEPARATOR)
    if owner == SHARED_NAMESPACE:
        return None, name
    return owner, name


@dataclass(frozen=True, slots=True)
class ComponentExtent:
    """One component's array identity and axis-aligned extent, for the co-location check.

    Neither field is evidence of coupling.  They exist so the validator can find the pairs that LOOK
    coupled -- same array, or overlapping coordinates -- and prove that a declared connector exists.
    """

    ref: str
    lower_um: tuple[float, float, float]
    upper_um: tuple[float, float, float]
    array_id: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "ref", _nonempty(self.ref, what="extent ref"))
        for name in ("lower_um", "upper_um"):
            values = tuple(getattr(self, name))
            if len(values) != 3 or any(
                isinstance(value, bool) or not isinstance(value, (int, float)) for value in values
            ):
                raise ValueError(f"{name} must be a (x, y, z) tuple of finite numbers")
            floats = tuple(float(value) for value in values)
            if any(value != value or value in (float("inf"), float("-inf")) for value in floats):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, floats)
        if any(hi < lo for lo, hi in zip(self.lower_um, self.upper_um, strict=True)):
            raise ValueError("upper_um must be componentwise >= lower_um")
        if self.array_id is not None:
            object.__setattr__(self, "array_id", _nonempty(self.array_id, what="array_id"))

    def overlaps(self, other: ComponentExtent) -> bool:
        """True when the two extents share an array or their closed extents intersect."""
        if self.array_id is not None and self.array_id == other.array_id:
            return True
        return all(
            self.lower_um[axis] <= other.upper_um[axis]
            and other.lower_um[axis] <= self.upper_um[axis]
            for axis in range(3)
        )


@dataclass(frozen=True, slots=True)
class ConnectorInstance:
    """One instantiated connector edge between two resolved component references."""

    instance_id: str
    connector_name: str
    endpoint_a_ref: str
    endpoint_b_ref: str
    owner: str
    bidirectional: bool = True
    adjoint_declared: bool = True

    def __post_init__(self) -> None:
        for name in (
            "instance_id",
            "connector_name",
            "endpoint_a_ref",
            "endpoint_b_ref",
            "owner",
        ):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        for name in ("bidirectional", "adjoint_declared"):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be bool")

    @property
    def endpoint_refs(self) -> frozenset[str]:
        """Return the instance's two resolved endpoint references."""
        return frozenset({self.endpoint_a_ref, self.endpoint_b_ref})


@dataclass(frozen=True, slots=True)
class CellJunction:
    """One explicit cell-cell junction; the ONLY sanctioned route between two cells' state."""

    junction_id: str
    law_id: str
    cell_a: str
    cell_b: str
    endpoint_a_ref: str
    endpoint_b_ref: str

    def __post_init__(self) -> None:
        for name in (
            "junction_id",
            "law_id",
            "cell_a",
            "cell_b",
            "endpoint_a_ref",
            "endpoint_b_ref",
        ):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))

    @property
    def endpoint_refs(self) -> frozenset[str]:
        """Return the junction's two resolved endpoint references."""
        return frozenset({self.endpoint_a_ref, self.endpoint_b_ref})


@dataclass(frozen=True, slots=True)
class CellWiring:
    """One cell instance's owned components, referenced shared components, and connector names."""

    cell_id: str
    archetype_id: str
    owned_components: tuple[str, ...]
    shared_components: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for name in ("cell_id", "archetype_id"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        if NAMESPACE_SEPARATOR in self.cell_id:
            raise ValueError(f"cell_id must not contain {NAMESPACE_SEPARATOR!r}")
        object.__setattr__(
            self, "owned_components", _unique(self.owned_components, what="owned_components")
        )
        object.__setattr__(
            self, "shared_components", _unique(self.shared_components, what="shared_components")
        )
        overlap = set(self.owned_components) & set(self.shared_components)
        if overlap:
            raise ValueError(
                "a component cannot be both owned and shared: " + ", ".join(sorted(overlap))
            )

    def ref(self, component_name: str) -> str:
        """Return this cell's namespaced reference for one component it owns."""
        if component_name not in self.owned_components:
            raise KeyError(f"cell {self.cell_id!r} does not own component {component_name!r}")
        return namespaced_ref(self.cell_id, component_name)


@dataclass(frozen=True, slots=True)
class SharedComponent:
    """ONE environment object referenced by many cells -- never duplicated per cell."""

    name: str
    object_id: str
    referencing_cell_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        for name in ("name", "object_id"):
            object.__setattr__(self, name, _nonempty(getattr(self, name), what=name))
        object.__setattr__(
            self,
            "referencing_cell_ids",
            _unique(self.referencing_cell_ids, what="referencing_cell_ids"),
        )
        if not self.referencing_cell_ids:
            raise ValueError(f"shared component {self.name!r} is referenced by no cell")

    @property
    def ref(self) -> str:
        """Return the single reference every referencing cell must point at."""
        return shared_ref(self.object_id)


@dataclass(frozen=True, slots=True)
class PopulationWiring:
    """A heterotypic population's cells, shared objects, connector instances, and junctions."""

    population_id: str
    cells: tuple[CellWiring, ...]
    shared: tuple[SharedComponent, ...] = ()
    connectors: tuple[ConnectorInstance, ...] = ()
    junctions: tuple[CellJunction, ...] = ()
    extents: tuple[ComponentExtent, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "population_id", _nonempty(self.population_id, what="population_id")
        )
        for name, kind in (
            ("cells", CellWiring),
            ("shared", SharedComponent),
            ("connectors", ConnectorInstance),
            ("junctions", CellJunction),
            ("extents", ComponentExtent),
        ):
            values = tuple(getattr(self, name))
            if any(not isinstance(value, kind) for value in values):
                raise TypeError(f"{name} must contain {kind.__name__} values")
            object.__setattr__(self, name, values)
        if not self.cells:
            raise ValueError("a population needs at least one cell")


def _check_connector_declaration(
    connector: ConnectorDeclaration, subject: str
) -> list[GrammarViolation]:
    findings: list[GrammarViolation] = []
    if not connector.bidirectional:
        findings.append(
            GrammarViolation(
                ViolationCode.ONE_SIDED_CONNECTOR,
                subject,
                f"connector {connector.name!r} is declared one-way; a mechanical connection must be "
                "bidirectional (Newton's 3rd law)",
            )
        )
    if not connector.adjoint_transfer_required:
        findings.append(
            GrammarViolation(
                ViolationCode.ADJOINT_NOT_DECLARED,
                subject,
                f"connector {connector.name!r} declares no adjoint closure; a one-sided declaration "
                "is refused",
            )
        )
    return findings


def _check_owner(
    connector: ConnectorDeclaration, grammar: EngineGrammar, subject: str
) -> list[GrammarViolation]:
    """Check the measured ownership rule for one connector.

    The rule read out of ``dispatch.validate_canonical_facade_claims`` is exact-once facade
    ownership.  The additional expectation -- the owner also owns an endpoint -- holds for every
    declared connector except the composite ``alpha2beta1_collagen_series`` joint, so a mechanical
    group exempts it and the exemption is reported as an ``EXPANSION`` note rather than dropped.
    """
    findings: list[GrammarViolation] = []
    claimed = grammar.facade_components.get(connector.owner)
    if claimed is None:
        findings.append(
            GrammarViolation(
                ViolationCode.OWNER_MISSING,
                subject,
                f"connector {connector.name!r} names owner {connector.owner!r}, which claims no "
                "components in the canonical dispatch manifest",
            )
        )
        return findings
    if connector.endpoints & set(claimed):
        return findings
    if connector.mechanical_group is not None:
        findings.append(
            GrammarViolation(
                ViolationCode.OWNER_NOT_ENDPOINT_COMPOSITE,
                subject,
                f"connector {connector.name!r} is owned by {connector.owner!r}, which owns neither "
                f"endpoint {sorted(connector.endpoints)}; permitted because the edge belongs to "
                f"mechanical group {connector.mechanical_group!r}, so the composite joint is one "
                "runtime object owned once",
                severity=Severity.EXPANSION,
            )
        )
        return findings
    findings.append(
        GrammarViolation(
            ViolationCode.OWNER_NOT_ENDPOINT,
            subject,
            f"connector {connector.name!r} is owned by {connector.owner!r}, which owns none of its "
            f"endpoints {sorted(connector.endpoints)} and belongs to no mechanical group",
        )
    )
    return findings


def validate_manifest(manifest: object, grammar: EngineGrammar) -> GrammarReport:
    """Validate one cell manifest's component/connector names against the canonical grammar.

    Args:
        manifest: Any object exposing ``manifest_id``, ``components`` and ``connectors`` -- a
            :class:`~aleph.virtual_cell.contracts.CellStateManifest` satisfies this without this
            module importing it.
        grammar: The canonical grammar, normally from :func:`load_engine_grammar`.

    Returns:
        A :class:`GrammarReport`.  ``ok`` is true when nothing refuses the manifest; composite-owner
        notes are ``EXPANSION`` findings and do not make it false.

    Raises:
        TypeError: If ``manifest`` does not expose the three required attributes.
    """
    if not isinstance(grammar, EngineGrammar):
        raise TypeError("grammar must be an EngineGrammar")
    for attribute in ("manifest_id", "components", "connectors"):
        if not hasattr(manifest, attribute):
            raise TypeError(f"manifest must expose {attribute!r}")
    subject_id = _nonempty(manifest.manifest_id, what="manifest_id")
    components = tuple(manifest.components)
    connectors = tuple(manifest.connectors)

    findings: list[GrammarViolation] = []
    known_components = grammar.component_names
    for name in components:
        if name not in known_components:
            findings.append(
                GrammarViolation(
                    ViolationCode.UNKNOWN_COMPONENT,
                    name,
                    f"{name!r} is not declared in the canonical census ({grammar.source})",
                )
            )
    present = set(components) & known_components

    touched: set[str] = set()
    for name in connectors:
        try:
            connector = grammar.connector(name)
        except KeyError:
            findings.append(
                GrammarViolation(
                    ViolationCode.UNKNOWN_CONNECTOR,
                    name,
                    f"{name!r} is not declared in the canonical census ({grammar.source})",
                )
            )
            continue
        missing_endpoints = sorted(connector.endpoints - present)
        if missing_endpoints:
            findings.append(
                GrammarViolation(
                    ViolationCode.ENDPOINT_NOT_PRESENT,
                    name,
                    f"connector {name!r} needs endpoint components {missing_endpoints}, which this "
                    "manifest does not declare",
                )
            )
            continue
        touched |= connector.endpoints
        findings.extend(_check_connector_declaration(connector, name))
        findings.extend(_check_owner(connector, grammar, name))

    for name in sorted(present - touched):
        findings.append(
            GrammarViolation(
                ViolationCode.UNCONNECTED_COMPONENT,
                name,
                f"component {name!r} is declared but no declared connector touches it; sharing a "
                "manifest is not a mechanical connection",
            )
        )

    return GrammarReport(
        subject_id=subject_id,
        findings=tuple(findings),
        checks_run=(
            "component-declared",
            "connector-declared",
            "connector-endpoints-present",
            "connector-bidirectional",
            "connector-adjoint-declared",
            "connector-single-owner",
            "component-connected",
        ),
        census_source=grammar.source,
    )


def check_colocation(
    extents: Sequence[ComponentExtent],
    connectors: Sequence[ConnectorInstance],
    junctions: Sequence[CellJunction] = (),
) -> tuple[GrammarViolation, ...]:
    """Report every co-located component pair that no declared connector or junction joins.

    Co-location -- a shared array, or an overlapping coordinate range -- is NEVER a connection.  This
    is the project's oldest invariant and the reason components own their own state.

    Args:
        extents: One extent per component reference; duplicate references are rejected.
        connectors: The declared connector instances.
        junctions: The declared cell-cell junctions.

    Returns:
        One ``CO_LOCATED_WITHOUT_CONNECTOR`` finding per unjoined co-located pair, sorted.

    Raises:
        ValueError: If two extents share a reference.
    """
    ordered = tuple(extents)
    refs = [extent.ref for extent in ordered]
    if len(set(refs)) != len(refs):
        raise ValueError("each component reference may declare only one extent")
    joined = {frozenset(instance.endpoint_refs) for instance in connectors}
    joined |= {frozenset(junction.endpoint_refs) for junction in junctions}

    findings: list[GrammarViolation] = []
    for index, first in enumerate(ordered):
        for second in ordered[index + 1 :]:
            if not first.overlaps(second):
                continue
            if frozenset({first.ref, second.ref}) in joined:
                continue
            pair = tuple(sorted((first.ref, second.ref)))
            shared_array = first.array_id is not None and first.array_id == second.array_id
            reason = (
                f"share array {first.array_id!r}" if shared_array else "have overlapping extents"
            )
            findings.append(
                GrammarViolation(
                    ViolationCode.CO_LOCATED_WITHOUT_CONNECTOR,
                    f"{pair[0]} | {pair[1]}",
                    f"{pair[0]} and {pair[1]} {reason} but no connector or junction joins them: "
                    "UNCONNECTED -- co-location is never connection",
                )
            )
    return tuple(sorted(findings, key=lambda finding: finding.subject))


def validate_population(wiring: PopulationWiring, grammar: EngineGrammar) -> GrammarReport:
    """Validate a heterotypic population's namespacing, sharing, and cross-cell wiring.

    The checks, in the order they run:

    1. cell IDs are unique, so two instances can never resolve to the same component reference;
    2. every per-cell component is namespaced by cell ID, and no per-cell component is declared as a
       shared object;
    3. every ``environment``-role component is declared SHARED -- one object referenced by many cells
       -- and never owned per cell in a multi-cell population;
    4. every connector instance names a declared connector, resolves to references that exist, and
       matches the declared edge's endpoint components;
    5. a connector instance whose two endpoints resolve to DIFFERENT cells is refused: cell-cell
       coupling goes through an explicit junction, never through two cells' cortices;
    6. every junction joins two distinct cells, and its law is checked against the canonical census.

    Args:
        wiring: The population wiring to validate.
        grammar: The canonical grammar, normally from :func:`load_engine_grammar`.

    Returns:
        A :class:`GrammarReport` for the population.
    """
    if not isinstance(wiring, PopulationWiring):
        raise TypeError("wiring must be a PopulationWiring")
    if not isinstance(grammar, EngineGrammar):
        raise TypeError("grammar must be an EngineGrammar")

    findings: list[GrammarViolation] = []
    cell_ids = [cell.cell_id for cell in wiring.cells]
    duplicates = sorted({cell_id for cell_id in cell_ids if cell_ids.count(cell_id) > 1})
    for cell_id in duplicates:
        findings.append(
            GrammarViolation(
                ViolationCode.DUPLICATE_CELL_ID,
                cell_id,
                f"cell ID {cell_id!r} appears more than once; namespacing cannot separate them",
            )
        )
    known_cells = set(cell_ids)
    environment = grammar.environment_components
    multicell = len(known_cells) > 1

    owned_ref_to_cell: dict[str, str] = {}
    ref_component: dict[str, str] = {}
    for cell in wiring.cells:
        for name in cell.owned_components:
            if name not in grammar.component_names:
                findings.append(
                    GrammarViolation(
                        ViolationCode.UNKNOWN_COMPONENT,
                        f"{cell.cell_id}:{name}",
                        f"{name!r} is not declared in the canonical census ({grammar.source})",
                    )
                )
                continue
            if multicell and name in environment:
                findings.append(
                    GrammarViolation(
                        ViolationCode.SHARED_COMPONENT_DUPLICATED,
                        f"{cell.cell_id}:{name}",
                        f"environment component {name!r} is owned per cell; in a population it must "
                        "be ONE shared object referenced by every cell, not duplicated",
                    )
                )
            reference = namespaced_ref(cell.cell_id, name)
            if reference in owned_ref_to_cell:
                findings.append(
                    GrammarViolation(
                        ViolationCode.COMPONENT_NAME_COLLISION,
                        reference,
                        f"reference {reference!r} is claimed by more than one cell instance",
                    )
                )
            owned_ref_to_cell[reference] = cell.cell_id
            ref_component[reference] = name
        for name in cell.shared_components:
            if name not in grammar.component_names:
                findings.append(
                    GrammarViolation(
                        ViolationCode.UNKNOWN_COMPONENT,
                        f"{cell.cell_id}:{name}",
                        f"{name!r} is not declared in the canonical census ({grammar.source})",
                    )
                )
            elif name not in environment:
                findings.append(
                    GrammarViolation(
                        ViolationCode.PER_CELL_COMPONENT_SHARED,
                        f"{cell.cell_id}:{name}",
                        f"component {name!r} owns per-cell state (role "
                        f"{grammar.component(name).role!r}) and cannot be one object shared by "
                        "several cells",
                    )
                )

    shared_refs: dict[str, str] = {}
    for shared in wiring.shared:
        if shared.name not in grammar.component_names:
            findings.append(
                GrammarViolation(
                    ViolationCode.UNKNOWN_COMPONENT,
                    shared.object_id,
                    f"{shared.name!r} is not declared in the canonical census ({grammar.source})",
                )
            )
            continue
        unknown = sorted(set(shared.referencing_cell_ids) - known_cells)
        if unknown:
            findings.append(
                GrammarViolation(
                    ViolationCode.SHARED_COMPONENT_UNDECLARED,
                    shared.object_id,
                    f"shared {shared.name!r} is referenced by unknown cells {unknown}",
                )
            )
        shared_refs[shared.ref] = shared.name
        ref_component[shared.ref] = shared.name

    declared_shared_names = {shared.name for shared in wiring.shared}
    for cell in wiring.cells:
        for name in cell.shared_components:
            if name in grammar.component_names and name not in declared_shared_names:
                findings.append(
                    GrammarViolation(
                        ViolationCode.SHARED_COMPONENT_UNDECLARED,
                        f"{cell.cell_id}:{name}",
                        f"cell {cell.cell_id!r} references shared {name!r} but the population "
                        "declares no such shared object",
                    )
                )

    known_refs = set(owned_ref_to_cell) | set(shared_refs)
    for instance in wiring.connectors:
        subject = instance.instance_id
        try:
            connector = grammar.connector(instance.connector_name)
        except KeyError:
            findings.append(
                GrammarViolation(
                    ViolationCode.UNKNOWN_CONNECTOR,
                    subject,
                    f"{instance.connector_name!r} is not declared in the canonical census "
                    f"({grammar.source})",
                )
            )
            continue
        findings.extend(_check_connector_declaration(connector, subject))
        findings.extend(_check_owner(connector, grammar, subject))
        if not instance.bidirectional:
            findings.append(
                GrammarViolation(
                    ViolationCode.ONE_SIDED_CONNECTOR,
                    subject,
                    "connector instance is declared one-way; a mechanical connection must be "
                    "bidirectional",
                )
            )
        if not instance.adjoint_declared:
            findings.append(
                GrammarViolation(
                    ViolationCode.ADJOINT_NOT_DECLARED,
                    subject,
                    "connector instance declares no adjoint closure; a one-sided declaration is "
                    "refused",
                )
            )
        if instance.owner != connector.owner:
            findings.append(
                GrammarViolation(
                    ViolationCode.OWNER_MISSING,
                    subject,
                    f"instance names owner {instance.owner!r} but the canonical manifest gives "
                    f"{connector.owner!r}; ownership is exact-once and is not per instance",
                )
            )

        endpoint_cells: list[str | None] = []
        endpoint_components: list[str | None] = []
        for reference in (instance.endpoint_a_ref, instance.endpoint_b_ref):
            if NAMESPACE_SEPARATOR not in reference:
                findings.append(
                    GrammarViolation(
                        ViolationCode.COMPONENT_NOT_NAMESPACED,
                        subject,
                        f"endpoint {reference!r} is a bare component name; per-cell state must be "
                        f"namespaced {'<cell>' + NAMESPACE_SEPARATOR + '<component>'} and a shared "
                        f"object must be {SHARED_NAMESPACE + NAMESPACE_SEPARATOR + '<object>'}",
                    )
                )
                endpoint_cells.append(None)
                endpoint_components.append(None)
                continue
            if reference not in known_refs:
                findings.append(
                    GrammarViolation(
                        ViolationCode.ENDPOINT_NOT_PRESENT,
                        subject,
                        f"endpoint {reference!r} resolves to nothing this population declares",
                    )
                )
                endpoint_cells.append(None)
                endpoint_components.append(None)
                continue
            endpoint_cells.append(owned_ref_to_cell.get(reference))
            endpoint_components.append(ref_component[reference])

        if all(component is not None for component in endpoint_components):
            declared = connector.endpoints
            instantiated = set(endpoint_components)
            if instantiated != declared and not (connector.internal and instantiated == declared):
                findings.append(
                    GrammarViolation(
                        ViolationCode.ENDPOINT_MISMATCH,
                        subject,
                        f"instance joins components {sorted(instantiated)} but connector "
                        f"{connector.name!r} declares {sorted(declared)}",
                    )
                )
        first, second = endpoint_cells
        if first is not None and second is not None and first != second:
            findings.append(
                GrammarViolation(
                    ViolationCode.CROSS_CELL_INTRACELLULAR_CONNECTOR,
                    subject,
                    f"intracellular connector {connector.name!r} is wired between cells {first!r} "
                    f"and {second!r}; cell-cell coupling must go through an explicit junction, "
                    "never directly between two cells' components",
                )
            )

    for junction in wiring.junctions:
        if junction.cell_a == junction.cell_b:
            findings.append(
                GrammarViolation(
                    ViolationCode.JUNCTION_ENDPOINTS_SAME_CELL,
                    junction.junction_id,
                    f"junction {junction.junction_id!r} names cell {junction.cell_a!r} twice",
                )
            )
        for reference in junction.endpoint_refs:
            if reference not in known_refs:
                findings.append(
                    GrammarViolation(
                        ViolationCode.ENDPOINT_NOT_PRESENT,
                        junction.junction_id,
                        f"junction endpoint {reference!r} resolves to nothing this population "
                        "declares",
                    )
                )
        if junction.law_id not in grammar.connector_names:
            findings.append(
                GrammarViolation(
                    ViolationCode.JUNCTION_LAW_NOT_IN_GRAMMAR,
                    junction.junction_id,
                    f"junction law {junction.law_id!r} is not a declared connector; the canonical "
                    f"census ({grammar.source}) declares NO cell-cell junction connector at all, so "
                    "every multicell coupling in this layer is a missing-law project, not a wiring "
                    "detail",
                    severity=Severity.EXPANSION,
                )
            )

    findings.extend(check_colocation(wiring.extents, wiring.connectors, wiring.junctions))

    return GrammarReport(
        subject_id=wiring.population_id,
        findings=tuple(findings),
        checks_run=(
            "cell-id-unique",
            "per-cell-component-namespaced",
            "shared-environment-not-duplicated",
            "connector-declared",
            "connector-endpoints-resolve",
            "connector-endpoint-components-match",
            "connector-single-owner",
            "connector-bidirectional",
            "connector-adjoint-declared",
            "cross-cell-connector-refused",
            "junction-law-declared",
            "co-location-is-not-connection",
        ),
        census_source=grammar.source,
    )
