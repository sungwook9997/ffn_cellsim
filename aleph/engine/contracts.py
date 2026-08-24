"""Declarative component and connector contracts for the FF-AC cell engine.

This module describes *ownership and connectivity*.  It deliberately contains no physics
implementation and never treats co-location in an array as a mechanical connection.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class ComponentRole(StrEnum):
    """Numerical role played by one independently owned cell-engine component."""

    SURFACE_BODY = "surface_body"
    FLUID_VOLUME = "fluid_volume"
    CORE_BODY = "core_body"
    ACTIVE_LOAD_PATH = "active_load_path"
    STRUCTURAL_RIG = "structural_rig"
    ENVIRONMENT = "environment"


class ConnectorFamily(StrEnum):
    """First-class dynamic joint families in the initial cell composition."""

    ERM = "erm"
    FA_CLUTCH = "fa_clutch"
    ACTIN_ANCHOR = "actin_anchor"
    TRANSIENT_ACTIN = "transient_actin"
    LINC = "linc"
    PLECTIN = "plectin"
    SPECTRAPLAKIN = "spectraplakin"
    MOTOR = "motor"
    IMMERSED_TRANSFER = "immersed_transfer"
    CHEMICAL_FLUX = "chemical_flux"
    FLUID_BOUNDARY = "fluid_boundary"
    FIBER_CROSSLINK = "fiber_crosslink"
    ENVIRONMENT_BOUNDARY = "environment_boundary"
    CONTACT = "contact"


class ConnectorScope(StrEnum):
    """Whether a connector joins actors or entities owned by one actor."""

    INTER_COMPONENT = "inter_component"
    INTERNAL = "internal"


class PortKind(StrEnum):
    """WHAT a declared connector is, physically — so that counting them means something.

    The count "38 connectors" mixes four different kinds of object, which is why it answers no
    question cleanly: it is not the number of couplings, nor of runtime objects, nor of force laws.
    It is the number of DECLARATIONS.  This enum names the kinds so a reader (and a gate) can ask for
    the one they meant.

    The vocabulary is the standard one for interconnected physical systems — bond graphs (Paynter),
    port-Hamiltonian systems (van der Schaft), Modelica ``connector``.  There a port carries a pair of
    conjugate variables, a *potential* (equal across the connection) and a *flow* (summing to zero),
    whose product is power; a connection is a power-conserving interconnection.  This repository
    already implements exactly that check for ONE cut — :meth:`GlobalCellLedger.assemble_balance`
    reduces two INDEPENDENTLY sourced force channels and requires them to cancel within a derived
    tolerance — and ``ledger.py`` already states the generalisation: *"any cut of the composed cell
    into two sub-bodies instantiates it."*  What was missing is the declaration of WHICH edges are
    supposed to satisfy it.

    ``POWER_PORT``
        A genuine coupling between two state-owning components.  Both sides own force arrays, and the
        flows must cancel.  This is the only kind for which the balance invariant is meaningful, and
        the only kind :meth:`CellArchitecture.power_ports` returns.
    ``RESERVOIR_PORT``
        The far side is not a body that moves — a far-field reference frame or an exterior medium
        (``dynamically_evolving=False``).  The reaction does not cancel against a peer; it accumulates
        in the environment reservoir.  ``ledger.py`` already keeps that as its own channel.
    ``SERIES_HALF``
        One of two semantic edges that a ``mechanical_group`` requires to bind ONE composite runtime.
        Two such declarations are ONE spring — counting both as couplings double-counts.
    ``INTERNAL_ELEMENT``
        Both endpoints are inside a single component (``scope=INTERNAL``).  It couples nothing, so it
        is a constitutive element of that component and not a connector in the port sense at all.

    DERIVED, never declared.  Every rule below reads fields that are already on the contracts, so this
    cannot drift away from the declaration it describes — which a 38-row hand-typed field would.
    """

    POWER_PORT = "power_port"
    RESERVOIR_PORT = "reservoir_port"
    SERIES_HALF = "series_half"
    INTERNAL_ELEMENT = "internal_element"


class EvidenceRung(StrEnum):
    """How far a component/connector has been wired and executed — a STRUCTURAL axis only.

    A rung answers "does this run, and against what": is the seam declared, are real Warp kernels
    bound behind it, has it executed on CUDA, is its connector actually dispatched inside an
    accepted-step transaction.  It answers **nothing** about whether a magnitude the run produced
    may be quoted; that is :class:`QuantitativeClaim`, a separate and orthogonal field.

    Splitting the two is PI decision **D1, option B (2026-07-28)**.  The 8-state ladder had no
    honest slot for "mechanically connected, quantitatively BLOCKED", which under the 2026-07-25
    reframe (mechanical connectedness is the gate, magnitudes explicitly are not) is the SUCCESS
    case, not a shortfall.  Inserting a ninth rung would have forced a permanent argument about
    whether it sorts above or below ``CONNECTED``.  It does not sort at all — it is a second axis,
    and the ``nmii_sf_motor`` artifact had already been written that way by hand
    (``MECHANISM PASS`` + ``quantitative_claim_status: BLOCKED``) before any vocabulary existed.

    The first two members sit BELOW the ladder because they prove wiring rather than physics:
    ``SCHEMA-PROOF`` (topology only) and ``CENSUS-WIRED`` (population census wired, no physics
    gate).  :data:`LADDER_FLOOR` is the boundary.

    This enum is the single definition.  It previously existed only as a literal tuple inside
    ``scripts/ac_architecture_dashboard.py`` — a visualization script — which is why every driver
    stamped its rung as a hand-typed free string and one of them stamped a rung it had not earned.
    """

    SCHEMA_PROOF = "SCHEMA-PROOF"
    CENSUS_WIRED = "CENSUS-WIRED"
    CONTRACTED = "CONTRACTED"
    SEAMED = "SEAMED"
    KERNEL_BOUND = "KERNEL_BOUND"
    CUDA_UNIT = "CUDA_UNIT"
    CONNECTED = "CONNECTED"
    NATIVE = "NATIVE"
    OPTIMISED = "OPTIMISED"
    PRODUCTION = "PRODUCTION"


class QuantitativeClaim(StrEnum):
    """Whether the MAGNITUDES a run produced may be quoted — orthogonal to :class:`EvidenceRung`.

    * ``BLOCKED`` — the mechanism is demonstrated but no number may leave the artifact: an
      unsourced parameter (PI-GAP), an unconverged inner solve, or a residual/signal ratio too
      large to separate physics from solver transient.  Under the current gate this is the
      expected and acceptable state, and it is the default.
    * ``OPEN`` — nothing structurally forbids quoting a magnitude, but none has been validated.
    * ``CONFIRMED`` — a magnitude is quotable, with the artifact naming what validated it.
    """

    BLOCKED = "BLOCKED"
    OPEN = "OPEN"
    CONFIRMED = "CONFIRMED"


class GateVerdict(StrEnum):
    """The outcome of ONE gate evaluation — ``PASS`` / ``FAIL`` / ``VOID`` (PI decision D7, 2026-07-28).

    ``VOID`` is the third verdict the repo lacked: the run produced a number, but the number cannot be
    interpreted as physics because the solver residual was too large a fraction of the signal being
    measured.  It is not a soft FAIL.  A FAIL says "the physics missed the target"; a VOID says
    "nothing was measured, so no claim of either kind is available from this run".

    Without it, solver artifacts get reported as physics, which has already happened here once: the
    **withdrawn 1.92 pN** SF tension passed its gate on an explicit relax that had not converged
    (residual/tension = 15.5%).  The live example is ``sf_implicit``, where the implicit solve dropped
    absolute residual ~1e5x but tension collapsed by a similar factor — so the RATIO is uninformative
    and neither PASS nor FAIL is honest.

    Distinct from :class:`QuantitativeClaim`, and the two are easy to conflate:

    * :class:`QuantitativeClaim` asks *may this NUMBER be quoted* — an artifact-level statement that
      coexists with a perfectly good run (``CONNECTED`` + ``BLOCKED`` is the success case).
    * :class:`GateVerdict` asks *is this RUN interpretable at all* — a gate-level statement.

    A ``VOID`` run should also carry ``QuantitativeClaim.BLOCKED``, but the converse does not hold.
    """

    PASS = "PASS"
    FAIL = "FAIL"
    VOID = "VOID"


@dataclass(frozen=True, slots=True)
class VoidCeiling:
    """A residual/signal ceiling, declared BEFORE the run that it judges.

    Declaring it after seeing the residual is gate-loosening — the run would choose the ceiling that
    lets it pass.  This type exists to make the ordering explicit in code: build it at configuration
    time, record it in the artifact, and hand it to :func:`classify_gate` afterwards.

    Attributes:
        ratio: The largest residual/signal fraction at which the gate's observable still means
            something.  Dimensionless, in ``(0, 1]``.
        rationale: Why this ceiling, recorded with it so the artifact carries the justification and
            not just the number.  Required.
    """

    ratio: float
    rationale: str

    def __post_init__(self) -> None:
        if not (0.0 < self.ratio <= 1.0):
            raise ValueError(
                f"void ceiling must be a fraction in (0, 1]; got {self.ratio!r}. A ceiling at or "
                "above 1 voids nothing: residual equal to signal would still be called measurable"
            )
        if not self.rationale.strip():
            raise ValueError("a void ceiling declared without a rationale is a bare magic number")


def classify_gate(
    *,
    passed: bool,
    residual: float,
    signal: float,
    ceiling: VoidCeiling,
) -> GateVerdict:
    """Return ``PASS`` / ``FAIL`` / ``VOID`` for one gate, with ``VOID`` overriding both others.

    The override is the entire point.  A run whose residual swamps its signal has not measured
    anything, so it must not be able to report a PASS — that is precisely how the withdrawn 1.92 pN
    was published.

    Args:
        passed: Whether the gate's own criterion was met, evaluated by the caller.
        residual: The solver residual, in the same units as ``signal``.
        signal: The magnitude of what the gate is measuring, same units as ``residual``.
        ceiling: The pre-declared :class:`VoidCeiling`.

    Returns:
        ``VOID`` if ``residual/signal`` exceeds the ceiling or is not computable; otherwise ``PASS``
        when ``passed`` and ``FAIL`` when not.

    Raises:
        ValueError: If ``residual`` or ``signal`` is negative or not finite.

    Note:
        A zero signal yields ``VOID``, never ``PASS``: a residual/signal criterion cannot be evaluated
        without a signal.  A gate whose observable is *legitimately* zero — "tension is exactly zero
        while no head is bound" is a real and valuable gate here — is asserting an absolute zero, not a
        ratio, and must not be routed through this classifier.
    """
    import math

    if not math.isfinite(residual) or not math.isfinite(signal):
        raise ValueError(f"residual and signal must be finite; got {residual!r}, {signal!r}")
    if residual < 0.0 or signal < 0.0:
        raise ValueError(f"residual and signal are magnitudes; got {residual!r}, {signal!r}")
    if signal == 0.0:
        return GateVerdict.VOID
    if residual / signal > ceiling.ratio:
        return GateVerdict.VOID
    return GateVerdict.PASS if passed else GateVerdict.FAIL


#: The rungs in order, below-ladder dump markers first (PI 2026-07-22 ladder, D1-B axis split).
EVIDENCE_ORDER: tuple[EvidenceRung, ...] = tuple(EvidenceRung)
_EVIDENCE_RANK: dict[str, int] = {rung.value: index for index, rung in enumerate(EVIDENCE_ORDER)}
#: At or above this rank the label is on the real 8-state ladder rather than a wiring marker.
LADDER_FLOOR: int = _EVIDENCE_RANK[EvidenceRung.CONTRACTED.value]


def rung_rank(name: str | EvidenceRung | None) -> int:
    """Rank an evidence rung; unknown or missing sorts below every declared rung.

    Args:
        name: A rung value, an :class:`EvidenceRung`, or ``None``.

    Returns:
        The rung's index in :data:`EVIDENCE_ORDER`, or ``-1`` when it is not a declared rung.
    """
    return _EVIDENCE_RANK.get(str(name), -1)


@dataclass(frozen=True, slots=True)
class EvidenceLabel:
    """The two-axis label an artifact carries: a structural rung plus a quantitative status.

    ``basis`` is required and carries no default.  A rung asserted without naming the measurement
    that earned it is the exact failure this class exists to stop — ``ac_gate_b_cortex_motor_native``
    stamped ``"evidence": "CONNECTED"`` as a string constant, so the artifact recorded an intention
    rather than an observation.

    Deliberately NOT encoded: any cross-axis rule beyond the one below.  Whether, say, ``CONFIRMED``
    should require a native full-population run is a gate-contract question that D1-B did not decide,
    so encoding it here would be a gate change nobody signed.  The single encoded rule is definitional
    rather than a new gate: a below-ladder marker proves topology or census, never physics, so it
    cannot carry a confirmed magnitude.
    """

    rung: EvidenceRung
    quantitative: QuantitativeClaim = QuantitativeClaim.BLOCKED
    basis: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.rung, EvidenceRung):
            raise TypeError("rung must be an EvidenceRung, not a free string")
        if not isinstance(self.quantitative, QuantitativeClaim):
            raise TypeError("quantitative must be a QuantitativeClaim, not a free string")
        if not self.basis.strip():
            raise ValueError(
                f"evidence rung {self.rung.value!r} needs a basis naming what was measured"
            )
        if (
            self.quantitative is QuantitativeClaim.CONFIRMED
            and rung_rank(self.rung) < LADDER_FLOOR
        ):
            raise ValueError(
                f"{self.rung.value!r} is a wiring marker below the ladder and proves no physics, "
                "so it cannot carry a CONFIRMED magnitude"
            )

    @property
    def on_ladder(self) -> bool:
        """Whether this rung is on the 8-state ladder rather than a below-ladder wiring marker."""
        return rung_rank(self.rung) >= LADDER_FLOOR

    def as_artifact_fields(self) -> dict[str, str]:
        """Return the flat fields a run artifact records, both axes always present."""
        return {
            "evidence": self.rung.value,
            "quantitative_claim_status": self.quantitative.value,
            "evidence_basis": self.basis,
        }


@dataclass(frozen=True, slots=True)
class ComponentContract:
    """Build-time declaration of one state-owning component."""

    name: str
    role: ComponentRole
    representation: str
    solver: str
    owns_geometry: bool = True
    dynamically_evolving: bool = True
    has_events: bool = False

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("component name must be non-empty")
        if not self.representation.strip():
            raise ValueError(f"component {self.name!r} needs a representation")
        if not self.solver.strip():
            raise ValueError(f"component {self.name!r} needs a solver")
        if self.has_events and not self.dynamically_evolving:
            raise ValueError(
                f"component {self.name!r} cannot own events without evolving under the transaction"
            )


@dataclass(frozen=True, slots=True)
class ConnectorContract:
    """Build-time declaration of a bidirectional inter-component force path."""

    name: str
    family: ConnectorFamily
    component_a: str
    component_b: str
    kinetics: bool
    commit_on_accept: bool
    scope: ConnectorScope = ConnectorScope.INTER_COMPONENT
    endpoint_role_a: str = "any"
    endpoint_role_b: str = "any"
    chemistry_card: str | None = None
    mechanical_group: str | None = None
    bidirectional: bool = True
    adjoint_transfer_required: bool = True
    generation_required: bool = False
    remap_on_accept: bool = False
    blocks_sleep_refine: bool = False

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("connector name must be non-empty")
        if not self.component_a.strip() or not self.component_b.strip():
            raise ValueError(f"connector {self.name!r} needs two component endpoints")
        same_component = self.component_a == self.component_b
        if self.scope is ConnectorScope.INTER_COMPONENT and same_component:
            raise ValueError(f"inter-component connector {self.name!r} must join distinct components")
        if self.scope is ConnectorScope.INTERNAL and not same_component:
            raise ValueError(f"internal connector {self.name!r} must stay inside one component")
        if self.kinetics and not self.commit_on_accept:
            raise ValueError(
                f"kinetic connector {self.name!r} must commit only on an accepted physical step"
            )
        if not self.bidirectional:
            raise ValueError(f"connector {self.name!r} cannot be a one-way load")
        if not self.endpoint_role_a.strip() or not self.endpoint_role_b.strip():
            raise ValueError(f"connector {self.name!r} needs non-empty endpoint roles")
        if self.mechanical_group is not None and not self.mechanical_group.strip():
            raise ValueError(f"connector {self.name!r} has an empty mechanical group")
        if (self.remap_on_accept or self.blocks_sleep_refine) and not self.generation_required:
            raise ValueError(
                f"connector {self.name!r} cannot remap/lock topology without generation checks"
            )


#: The component that HOLDS the resting osmotic turgor Pi_0 (PI decision **D3, option C**, 2026-07-28).
#:
#: Variant A — the cortex is the pressurised envelope — is ratified and stays, because every native
#: resting number in the repo was measured on it (``ac/cell/assemble.py`` builds the pressure field at
#: ``PI_0_PA`` and the force manifest records that Pi_0 sets the entire resting cortical tension through
#: ``gamma = dP*R/2``).  What option C adds is that the membrane must then receive that load through
#: DECLARED CONNECTORS and never by carrying Pi_0 itself: co-location is not a connection, so a membrane
#: that independently applied Pi_0 would be the framework's double-count trap, counting one osmotic load
#: on two surfaces.  :meth:`CellArchitecture.osmotic_envelope_load_path` is the structural check that the
#: transmission path exists at all; it cannot check magnitudes, and does not claim to.
#:
#: This also settles the Variant A/B question open since 2026-07-16 without invalidating a measurement:
#: A and B differ in WHERE the pressure boundary sits, but once transmission is an explicit connector the
#: difference collapses to that connector's stiffness, which is a parameter, not a contract.
OSMOTIC_ENVELOPE: str = "cortex"

#: The surface the envelope must reach: the outer boundary Pi_0 ultimately loads.
OSMOTIC_LOADED_SURFACE: str = "membrane"


@dataclass(frozen=True, slots=True)
class CellArchitecture:
    """Validated component registry plus its explicit connector graph."""

    components: tuple[ComponentContract, ...]
    connectors: tuple[ConnectorContract, ...]

    def __post_init__(self) -> None:
        component_names = tuple(component.name for component in self.components)
        if len(set(component_names)) != len(component_names):
            raise ValueError("component names must be unique")
        connector_names = tuple(connector.name for connector in self.connectors)
        if len(set(connector_names)) != len(connector_names):
            raise ValueError("connector names must be unique")
        known = set(component_names)
        for connector in self.connectors:
            missing = {connector.component_a, connector.component_b} - known
            if missing:
                raise ValueError(
                    f"connector {connector.name!r} references unknown components {sorted(missing)}"
                )

    def component(self, name: str) -> ComponentContract:
        """Return a named component or raise a useful error."""
        for component in self.components:
            if component.name == name:
                return component
        raise KeyError(name)

    def port_kind(self, name: str) -> PortKind:
        """Return WHAT the named connector is, derived from its own declaration.

        The four rules are checked in an order that is not arbitrary and is asserted rather than
        assumed: an edge matching two rules would be classified by whichever ran first, silently, so
        a co-occurrence that the rules do not intend raises instead.

        Args:
            name: a declared connector name.

        Returns:
            Its :class:`PortKind`.

        Raises:
            KeyError: if no connector is declared under ``name``.
            ValueError: if the declaration matches two rules that must not co-occur — an internal
                element cannot also be a series half or reach a reservoir, because both of those
                require two distinct endpoints.
        """
        connector = next((c for c in self.connectors if c.name == name), None)
        if connector is None:
            raise KeyError(name)
        internal = connector.scope is ConnectorScope.INTERNAL
        reservoir = any(
            not self.component(endpoint).dynamically_evolving
            for endpoint in (connector.component_a, connector.component_b)
        )
        series = connector.mechanical_group is not None
        if internal and (reservoir or series):
            raise ValueError(
                f"connector {name!r} is declared INTERNAL but also "
                f"{'reaches a non-evolving endpoint' if reservoir else 'carries a mechanical group'}; "
                "both need two distinct endpoints, so the declaration is inconsistent"
            )
        if internal:
            return PortKind.INTERNAL_ELEMENT
        if reservoir:
            return PortKind.RESERVOIR_PORT
        if series:
            return PortKind.SERIES_HALF
        return PortKind.POWER_PORT

    def power_ports(self) -> tuple[str, ...]:
        """Return the connectors for which the force-balance invariant is meaningful.

        This is the honest answer to "how many connectors are there", for the question that usually
        means: how many places do two state-owning bodies exchange force. It is smaller than the
        declaration count, and it is the set a coupled solve has to converge.
        """
        return tuple(
            c.name for c in self.connectors if self.port_kind(c.name) is PortKind.POWER_PORT
        )

    def balance_cuts(self) -> dict[tuple[str, str], tuple[str, ...]]:
        """Return every two-body cut the force-balance gate can score, and what each cut tests.

        The gate's unit is a CUT, not a connector: :meth:`GlobalCellLedger.add_body_force` reduces each
        sub-body's whole force array, so what it verifies is that everything crossing the cut is an
        adjoint pair.  A cut therefore resolves a single connector only when that connector is the two
        components' ONLY power port — otherwise it tests them jointly and cannot attribute a failure.

        Measured on the reference architecture: 27 of the 28 power ports are the sole coupling between
        their endpoints, so the cut IS the connector for all but one.  The exception is
        ``cortex``/``membrane``, joined by ``membrane_erm_cortex`` and ``membrane_cortex_contact``,
        which is deliberate — the ERM tether is force-free in compression, so the contact edge is the
        only one that can carry an outward load, and the two are complementary rather than redundant.
        A cut across that pair scores both together, and an artifact must say so rather than name one.

        Returns:
            ``{(component_a, component_b): (connector names…)}``, endpoints sorted, listing only
            :attr:`PortKind.POWER_PORT` edges.  A one-element value is an attributable cut; a longer
            one is not.
        """
        cuts: dict[tuple[str, str], list[str]] = {}
        for connector in self.connectors:
            if self.port_kind(connector.name) is not PortKind.POWER_PORT:
                continue
            key = tuple(sorted((connector.component_a, connector.component_b)))
            cuts.setdefault(key, []).append(connector.name)  # type: ignore[arg-type]
        return {key: tuple(names) for key, names in sorted(cuts.items())}

    def attributable_cuts(self) -> dict[tuple[str, str], str]:
        """Return only the cuts that isolate ONE connector, mapped to that connector's name.

        These are the edges whose balance failure names itself.  Everything else needs a per-connector
        accumulator before a failure can be attributed, and pretending otherwise would put a
        connector's name on evidence that does not single it out.
        """
        return {
            key: names[0] for key, names in self.balance_cuts().items() if len(names) == 1
        }

    def neighbors(self, name: str) -> frozenset[str]:
        """Return components joined to ``name`` by registered connector families."""
        self.component(name)
        adjacent: set[str] = set()
        for connector in self.connectors:
            if connector.component_a == name:
                adjacent.add(connector.component_b)
            elif connector.component_b == name:
                adjacent.add(connector.component_a)
        return frozenset(adjacent)

    def force_path(self, source: str, target: str) -> tuple[str, ...] | None:
        """Return one component-level connector path, or ``None`` when no path exists."""
        self.component(source)
        self.component(target)
        queue: list[tuple[str, tuple[str, ...]]] = [(source, (source,))]
        visited = {source}
        while queue:
            current, path = queue.pop(0)
            if current == target:
                return path
            for neighbor in sorted(self.neighbors(current)):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, (*path, neighbor)))
        return None

    def osmotic_envelope_load_path(
        self,
        envelope: str = OSMOTIC_ENVELOPE,
        loaded_surface: str = OSMOTIC_LOADED_SURFACE,
    ) -> tuple[ConnectorContract, ...]:
        """Return the declared connectors that carry turgor from the envelope to the loaded surface.

        Args:
            envelope: The component holding Pi_0.  Defaults to :data:`OSMOTIC_ENVELOPE`.
            loaded_surface: The surface that must receive the load.  Defaults to
                :data:`OSMOTIC_LOADED_SURFACE`.

        Returns:
            Every connector joining the two, in declaration order.

        Raises:
            KeyError: If either component is not registered.
            ValueError: If no connector joins them — meaning the pressurised envelope could not load
                the surface at all, which is the structural form of the Variant A/B defect.

        Note:
            The path is direct by construction: turgor is a surface load, so a multi-hop route through
            a third component would be transmitting pressure through something that does not hold it.
            Both edges matter and neither is redundant — an ERM tether is force-free in compression
            ("a molecular linker is a tether, not a strut"), so tension alone cannot balance an
            outward pressure; the CONTACT edge is what carries it.
        """
        self.component(envelope)
        self.component(loaded_surface)
        pair = {envelope, loaded_surface}
        path = tuple(
            connector
            for connector in self.connectors
            if {connector.component_a, connector.component_b} == pair
        )
        if not path:
            raise ValueError(
                f"osmotic envelope {envelope!r} has no declared connector to {loaded_surface!r}: "
                "the pressurised envelope cannot load the surface, so turgor would have to be "
                "applied twice or not at all"
            )
        return path

    def mechanical_group(self, name: str) -> tuple[ConnectorContract, ...]:
        """Return semantic edges resolved by one composite mechanical joint."""
        grouped = tuple(connector for connector in self.connectors if connector.mechanical_group == name)
        if not grouped:
            raise KeyError(name)
        return grouped


def reference_cell_architecture() -> CellArchitecture:
    """Return the fixed top-level composition for the first connected cell-engine build.

    The graph declares allowed component-level load paths.  Runtime bond populations remain
    geometry- and kinetics-dependent; this function does not pre-bind molecular joints.

    Every component is declared ``has_events=False``: which components own an intra-component event
    kernel (cortex nucleation/severing, NMII assembly, …) is a per-component biology-phase decision
    that must arrive with a sourced rate (Magic-Number Block), so the flag is flipped then, not here.
    Inter-component events stay connector-owned (``ConnectorContract.kinetics``); those are enforced
    to expose ``propose_events`` today.
    """
    components = (
        ComponentContract(
            "membrane", ComponentRole.SURFACE_BODY, "fluid Helfrich surface FEM", "surface solve"
        ),
        ComponentContract(
            "cortex", ComponentRole.SURFACE_BODY,
            "explicit crosslinked F-actin network; gated surface condensation",
            "network mechanics or validated surface solve",
        ),
        ComponentContract(
            "cytosol", ComponentRole.FLUID_VOLUME, "Biot/Darcy finite-volume field", "field solve"
        ),
        ComponentContract(
            "nucleus", ComponentRole.CORE_BODY,
            "native lamina/chromatin FEM; gated modal reduction",
            "native or validated reduced body solve",
        ),
        ComponentContract(
            "sf_arc", ComponentRole.ACTIVE_LOAD_PATH, "active rod/cable graph", "graph mechanics"
        ),
        ComponentContract(
            "focal_adhesion", ComponentRole.ACTIVE_LOAD_PATH, "stochastic clutch graph", "joint KMC",
            owns_geometry=False,
        ),
        ComponentContract(
            "microtubule", ComponentRole.STRUCTURAL_RIG, "dynamic rod graph", "graph mechanics"
        ),
        ComponentContract(
            "intermediate_filament", ComponentRole.STRUCTURAL_RIG, "nonlinear cable graph", "graph mechanics"
        ),
        ComponentContract(
            "lamellipodium", ComponentRole.ACTIVE_LOAD_PATH,
            "local adaptive explicit Arp2/3 branched F-actin",
            "adaptive branch-graph mechanics",
        ),
        ComponentContract(
            "filopodium", ComponentRole.ACTIVE_LOAD_PATH,
            "local adaptive explicit bundled F-actin",
            "adaptive bundle-graph mechanics",
        ),
        ComponentContract(
            "nmii", ComponentRole.ACTIVE_LOAD_PATH,
            "explicit Stam-Hocky bipolar backbone + individual heads",
            "head-resolved internal mechanics plus graph MOTOR KMC",
        ),
        ComponentContract(
            "ecm", ComponentRole.ENVIRONMENT, "fiber/crosslink world", "environment mechanics"
        ),
        ComponentContract(
            "world_boundary", ComponentRole.ENVIRONMENT,
            "physiological far-field reference frame", "boundary constraint",
            owns_geometry=False, dynamically_evolving=False,
        ),
        ComponentContract(
            # 14th component. PI decision D5, option A (2026-07-28): DECLARED NOW, IMPLEMENTED AT T10.
            #
            # Read the second half of that sentence as load-bearing. Declaring this component does NOT
            # fix what it exists to fix. The argument for it was never magnitude (~1e-5 pN at 30 nm/s);
            # it is that an exterior medium is the physically correct regulariser of the cell's six
            # whole-body rigid modes and the only path that can LOAD them. A declaration reserves the
            # contract slot at one-line cost and buys exactly zero of that — the rigid modes stay
            # regularised by numerics until T10 lands a Stokes/Brinkman exterior (or at minimum a
            # boundary-integral surface mobility) behind this seam. Severity on the plan's own scale:
            # BLOCKS_CRAWL.
            #
            # The failure mode this guards against is already in the tree: sf_implicit.py:66-72's
            # regulariser numerically masks unbound minifilaments as free rigid bodies. Same class.
            #
            # T10 implementation carries ONE hard precondition (PI framework trap #4): when the medium
            # starts dragging, the per-node drag term `gamma_node = 6*pi*eta*R/Nc` must be removed in
            # the SAME change. Leaving both is a 2x error, not a modelling choice.
            "extracellular_medium", ComponentRole.ENVIRONMENT,
            "exterior Stokes/Brinkman medium (free fluid: k -> inf, so Darcy is wrong here)",
            "exterior fluid solve or boundary-integral surface mobility",
            owns_geometry=False, dynamically_evolving=False,
        ),
    )
    connectors = (
        ConnectorContract(
            "membrane_erm_cortex", ConnectorFamily.ERM, "membrane", "cortex", True, True,
            endpoint_role_a="membrane material point", endpoint_role_b="cortex material point",
            chemistry_card="ezrin_membrane_f_actin",
        ),
        ConnectorContract(
            "surface_porous_transfer", ConnectorFamily.IMMERSED_TRANSFER, "cortex", "cytosol", False, False,
            endpoint_role_a="cortical shell quadrature",
            endpoint_role_b="porous-fluid transfer stencil",
        ),
        ConnectorContract(
            "membrane_cytosol_boundary", ConnectorFamily.FLUID_BOUNDARY, "membrane", "cytosol", False, False,
            endpoint_role_a="live membrane face quadrature",
            endpoint_role_b="fluid boundary stencil",
        ),
        ConnectorContract(
            "nucleus_cytosol_boundary", ConnectorFamily.FLUID_BOUNDARY, "nucleus", "cytosol", False, False,
            endpoint_role_a="moving impermeable nuclear envelope quadrature",
            endpoint_role_b="interior fluid boundary stencil",
        ),
        ConnectorContract(
            "fa_actin_anchor", ConnectorFamily.ACTIN_ANCHOR, "sf_arc", "focal_adhesion", True, True,
            endpoint_role_a="ventral end or dorsal basal end", endpoint_role_b="FA actin-side state",
            chemistry_card="actin_talin_integrin", mechanical_group="alpha2beta1_collagen_series",
        ),
        ConnectorContract(
            "integrin_collagen_clutch", ConnectorFamily.FA_CLUTCH, "focal_adhesion", "ecm", True, True,
            endpoint_role_a="FA ligand-side state", endpoint_role_b="collagen ligand material point",
            chemistry_card="alpha2beta1_collagen", mechanical_group="alpha2beta1_collagen_series",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            "sf_cortex_transient", ConnectorFamily.TRANSIENT_ACTIN, "sf_arc", "cortex", True, True,
            endpoint_role_a="stress-fiber or arc material point",
            endpoint_role_b="cortical actin material point",
            chemistry_card="transient_actin_crosslink",
        ),
        ConnectorContract(
            "actin_cap_linc", ConnectorFamily.LINC, "sf_arc", "nucleus", True, True,
            endpoint_role_a="perinuclear actin-cap material point",
            endpoint_role_b="nuclear-envelope LINC actin-facing site",
            chemistry_card="nesprin_actin_linc",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            "mt_nucleus_linc", ConnectorFamily.LINC, "microtubule", "nucleus", True, True,
            endpoint_role_a="microtubule plus-end or motor-bound lattice site",
            endpoint_role_b="nuclear-envelope LINC microtubule-facing site",
            chemistry_card="microtubule_motor_linc",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            "mt_cortex_capture", ConnectorFamily.MOTOR, "microtubule", "cortex", True, True,
            endpoint_role_a="dynamic microtubule plus-end or lattice motor site",
            endpoint_role_b="cortical dynein capture site",
            chemistry_card="microtubule_cortical_dynein",
        ),
        ConnectorContract(
            "if_nucleus_linc", ConnectorFamily.LINC, "intermediate_filament", "nucleus", True, True,
            endpoint_role_a="intermediate-filament junction or material point",
            endpoint_role_b="nuclear-envelope LINC IF-facing site",
            chemistry_card="nesprin3_plectin_if_linc",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            "if_sf_plectin", ConnectorFamily.PLECTIN, "intermediate_filament", "sf_arc", True, True,
            endpoint_role_a="intermediate-filament material point",
            endpoint_role_b="stress-fiber or arc material point",
            chemistry_card="plectin_if_actin",
        ),
        ConnectorContract(
            "mt_sf_spectraplakin", ConnectorFamily.SPECTRAPLAKIN, "microtubule", "sf_arc", True, True,
            endpoint_role_a="microtubule lattice site",
            endpoint_role_b="stress-fiber or arc material point",
            chemistry_card="spectraplakin_mt_actin",
        ),
        ConnectorContract(
            "sf_cytosol_transfer", ConnectorFamily.IMMERSED_TRANSFER, "sf_arc", "cytosol", False, False,
            endpoint_role_a="stress-fiber or arc quadrature",
            endpoint_role_b="porous-fluid transfer stencil",
        ),
        ConnectorContract(
            "mt_cytosol_transfer", ConnectorFamily.IMMERSED_TRANSFER, "microtubule", "cytosol", False, False,
            endpoint_role_a="microtubule quadrature",
            endpoint_role_b="porous-fluid transfer stencil",
        ),
        ConnectorContract(
            "if_cytosol_transfer", ConnectorFamily.IMMERSED_TRANSFER, "intermediate_filament", "cytosol", False, False,
            endpoint_role_a="intermediate-filament quadrature",
            endpoint_role_b="porous-fluid transfer stencil",
        ),
        ConnectorContract(
            "lamellipodium_membrane_contact", ConnectorFamily.CONTACT,
            "lamellipodium", "membrane", True, True,
            endpoint_role_a="growing barbed-end material point",
            endpoint_role_b="membrane contact quadrature",
            chemistry_card="actin_membrane_brownian_ratchet_contact",
        ),
        ConnectorContract(
            "lamellipodium_cortex_seam", ConnectorFamily.TRANSIENT_ACTIN,
            "lamellipodium", "cortex", True, True,
            endpoint_role_a="lamellipodial network rear seam",
            endpoint_role_b="cortical actin material point",
            chemistry_card="transient_actin_crosslink",
        ),
        ConnectorContract(
            "lamellipodium_cytosol_transfer", ConnectorFamily.IMMERSED_TRANSFER,
            "lamellipodium", "cytosol", False, False,
            endpoint_role_a="active filament quadrature and barbed-end sink",
            endpoint_role_b="porous drag and G-actin transport stencil",
        ),
        ConnectorContract(
            "lamellipodium_nascent_fa", ConnectorFamily.ACTIN_ANCHOR,
            "lamellipodium", "focal_adhesion", True, True,
            endpoint_role_a="lamellipodial actin material point",
            endpoint_role_b="nascent FA actin-side state",
            chemistry_card="actin_talin_integrin_nascent",
            mechanical_group="lamellipodium_nascent_series",
        ),
        # PI 2026-08-09 (option (a)): a nascent adhesion is a SERIES, like the mature one.  Without an
        # ECM-side partner the actin anchor has only ONE geometric endpoint — `focal_adhesion` declares
        # `owns_geometry=False` — so there is nothing to pull against and no mechanical runtime can be
        # built for it.  The mature FA works precisely because `fa_actin_anchor` and
        # `integrin_collagen_clutch` share one mechanical group, so the spring runs actin -> collagen and
        # the geometry-less FA is the middle both edges are registered under.  Each nascent adhesion gets
        # its OWN group: folding them into `alpha2beta1_collagen_series` would make `bind_connector`
        # demand ONE runtime for the SF, lamellipodial and filopodial anchors, i.e. one spring shared by
        # three different adhesions.
        ConnectorContract(
            "lamellipodium_nascent_clutch", ConnectorFamily.FA_CLUTCH,
            "focal_adhesion", "ecm", True, True,
            endpoint_role_a="nascent FA ligand-side state",
            endpoint_role_b="collagen ligand material point",
            chemistry_card="alpha2beta1_collagen_nascent",
            mechanical_group="lamellipodium_nascent_series",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            "filopodium_membrane_tip", ConnectorFamily.CONTACT,
            "filopodium", "membrane", True, True,
            endpoint_role_a="bundled barbed-end tip",
            endpoint_role_b="membrane tip contact quadrature",
            chemistry_card="actin_membrane_brownian_ratchet_contact",
        ),
        ConnectorContract(
            "filopodium_cortex_root", ConnectorFamily.TRANSIENT_ACTIN,
            "filopodium", "cortex", True, True,
            endpoint_role_a="filopodium bundle root",
            endpoint_role_b="cortical actin material point",
            chemistry_card="formin_fascin_root_coupling",
        ),
        ConnectorContract(
            "filopodium_cytosol_transfer", ConnectorFamily.IMMERSED_TRANSFER,
            "filopodium", "cytosol", False, False,
            endpoint_role_a="bundle drag quadrature and barbed-end sink",
            endpoint_role_b="porous drag and G-actin transport stencil",
        ),
        ConnectorContract(
            "filopodium_nascent_fa", ConnectorFamily.ACTIN_ANCHOR,
            "filopodium", "focal_adhesion", True, True,
            endpoint_role_a="filopodium base or shaft actin material point",
            endpoint_role_b="nascent FA actin-side state",
            chemistry_card="actin_talin_integrin_nascent",
            mechanical_group="filopodium_nascent_series",
        ),
        # PI 2026-08-09 (option (a)) — see `lamellipodium_nascent_clutch` for why the series is required
        # and why each nascent adhesion needs its own group rather than joining the mature one.
        ConnectorContract(
            "filopodium_nascent_clutch", ConnectorFamily.FA_CLUTCH,
            "focal_adhesion", "ecm", True, True,
            endpoint_role_a="nascent FA ligand-side state",
            endpoint_role_b="collagen ligand material point",
            chemistry_card="alpha2beta1_collagen_nascent",
            mechanical_group="filopodium_nascent_series",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            "nmii_sf_motor", ConnectorFamily.MOTOR, "nmii", "sf_arc", True, True,
            endpoint_role_a="individual NMII head crossbridge",
            endpoint_role_b="live polar actin material coordinate",
            chemistry_card="nmii_head_actin_hill_bell",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            "nmii_cortex_motor", ConnectorFamily.MOTOR, "nmii", "cortex", True, True,
            endpoint_role_a="individual NMII head crossbridge",
            endpoint_role_b="live polar actin material coordinate",
            chemistry_card="nmii_head_actin_hill_bell",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            "nmii_lamellipodium_motor", ConnectorFamily.MOTOR,
            "nmii", "lamellipodium", True, True,
            endpoint_role_a="individual NMII head crossbridge",
            endpoint_role_b="live polar actin material coordinate",
            chemistry_card="nmii_head_actin_hill_bell",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            "nmii_filopodium_motor", ConnectorFamily.MOTOR,
            "nmii", "filopodium", True, True,
            endpoint_role_a="individual NMII head crossbridge",
            endpoint_role_b="live polar actin material coordinate",
            chemistry_card="nmii_head_actin_hill_bell",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            "dorsal_arc_crosslink", ConnectorFamily.TRANSIENT_ACTIN, "sf_arc", "sf_arc", True, True,
            scope=ConnectorScope.INTERNAL, endpoint_role_a="dorsal free end",
            endpoint_role_b="transverse arc material point",
        ),
        ConnectorContract(
            "ecm_crosslink", ConnectorFamily.FIBER_CROSSLINK, "ecm", "ecm", True, True,
            scope=ConnectorScope.INTERNAL, endpoint_role_a="collagen segment",
            endpoint_role_b="collagen segment", chemistry_card="collagen_crosslink",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            # RESOLVED 2026-08-09 (PI): `the-far-field-anchor-is-a-clamp-and-the-registry-calls-it-a-
            # connector` is ratified as a CONNECTOR, because its premise conflates two different fields.
            # `dynamically_evolving=False` on `world_boundary` governs POSITION — the far-field frame does
            # not move.  `adjoint_transfer_required=True` governs FORCE — the reaction is accumulated into
            # the environment reservoir rather than discarded.  A clamp would discard it, and then the
            # balance ledger's two channels would agree by OMISSION, which is the failure mode a one-sided
            # scatter is supposed to be caught by.  `medium_exterior` already names the far side "the
            # reservoir's half of the adjoint pair".  Same reasoning applies to `membrane_medium_traction`.
            # Fixed by `test_contracts::test_environment_boundaries_are_inert_in_position_not_in_force`.
            "ecm_far_field_anchor", ConnectorFamily.ENVIRONMENT_BOUNDARY,
            "ecm", "world_boundary", False, False,
            endpoint_role_a="finite collagen boundary material point",
            endpoint_role_b="physiological far-field reference frame",
        ),
        ConnectorContract(
            # 36th connector, the load path of the 14th component. DECLARED NOW, IMPLEMENTED AT T10
            # (PI D5-A, 2026-07-28) — see the `extracellular_medium` component note for why declaring
            # it does not yet regularise the rigid modes.
            #
            # This is the free face of the ASYMMETRIC world boundary the plan's baseline requires:
            # basal face = 2D collagen ECM + far-field anchor, free face = media. A symmetric boundary
            # can neither spread nor crawl, so this edge is not decoration.
            "membrane_medium_traction", ConnectorFamily.ENVIRONMENT_BOUNDARY,
            "membrane", "extracellular_medium", False, False,
            endpoint_role_a="live membrane surface quadrature",
            endpoint_role_b="exterior medium traction stencil",
        ),
        ConnectorContract(
            "membrane_ecm_contact", ConnectorFamily.CONTACT, "membrane", "ecm", False, False,
            endpoint_role_a="live membrane contact quadrature",
            endpoint_role_b="live collagen segment material point",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        # --- PI-ratified 2026-07-25 (ARCHITECTURE_COMPLETENESS_2026-07-25.md §2) -------------------
        # Three connectors the completeness audit found missing. All non-kinetic, all force-free at the
        # physiological gap, so none introduces a sourced magnitude beyond a numerical exclusion stiffness
        # derived from the existing WCA / derived-mobility pattern.
        ConnectorContract(
            # The cortex cannot PUSH the nucleus today: all three nucleus solid edges are tension tethers
            # (actin_cap_linc / mt_nucleus_linc / if_nucleus_linc) and the only volumetric edge carries
            # isotropic pressure, so the interaction was sign-unavailable rather than merely un-routed.
            # Measured consequence: cortex moved 39.8 nm mean / 751 nm max under 10% strain while nuclear
            # displacement was 0.0 nm and flattening 0.0% (INTERNAL_DISPLACEMENT_DIAGNOSIS_2026-07-16.md).
            "nucleus_cortex_contact", ConnectorFamily.CONTACT, "nucleus", "cortex", False, False,
            endpoint_role_a="live nuclear envelope contact quadrature",
            endpoint_role_b="cortex material point",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            # The compressive channel the ERM tether structurally cannot supply: erm_cortex_connector.py:34
            # states a compressed / broken / resting tether is force-free — "a molecular linker is a tether,
            # not a strut" — so membrane turgor cannot be balanced by cortex form-finding through ERM alone
            # (RESTING_BASELINE_DIAGNOSIS_2026-07-23c.md:137 calls that path "proven inert"). Specified but
            # never promoted into the tuple by SURFACE_BODY_PLAN.md:198-199,220,257,394: "no binding kinetics
            # … a surface contact coupler … separate from ERM so collision cannot manufacture an adhesion
            # bond … force-free at the physiological gap", with no no-slip constraint.
            # RATIFIED AS THE TURGOR TRANSMISSION PATH (PI decision D3, option C, 2026-07-28). The
            # osmotic envelope stays the CORTEX (Variant A, :data:`OSMOTIC_ENVELOPE`) so no measured
            # resting number is invalidated; what changes is that the membrane receives that load only
            # through this edge plus `membrane_erm_cortex`, and never by applying Pi_0 itself — which
            # would count one osmotic load on two surfaces. The two edges are complementary, not
            # redundant: ERM is a tension tether that is force-free in compression (see the note above),
            # so this CONTACT edge is the only one that can carry an outward pressure. Still
            # declared-only — this records WHICH connector owns the path, not that it runs.
            "membrane_cortex_contact", ConnectorFamily.CONTACT, "membrane", "cortex", False, False,
            endpoint_role_a="live membrane contact quadrature",
            endpoint_role_b="cortex material point",
            generation_required=True, remap_on_accept=True, blocks_sleep_refine=True,
        ),
        ConnectorContract(
            # nmii was the only geometry-owning, dynamically-evolving component with no fluid coupling: its
            # four MOTOR edges were its only edges, so an all-heads-unbound minifilament had zero force, zero
            # drag and zero connectors — a free rigid body whose null space sf_implicit.py:66-72's regularizer
            # was closing numerically. This is the physical fix for that.
            "nmii_cytosol_transfer", ConnectorFamily.IMMERSED_TRANSFER, "nmii", "cytosol", False, False,
            endpoint_role_a="minifilament backbone and head quadrature",
            endpoint_role_b="porous-fluid transfer stencil",
        ),
    )
    return CellArchitecture(components=components, connectors=connectors)
