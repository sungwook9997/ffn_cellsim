"""Ring A/B/C archetypes and the endpoint/owner grammar validator.

The census these tests assert against is READ FROM THE ENGINE SOURCE, never transcribed here.  The
only literal counts in this file are in :func:`test_canonical_census_matches_state_document`, which
exists precisely to catch the day the status document and the code disagree -- and it is written so
that the code wins.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from aleph.virtual_cell.archetypes import (
    ARCHETYPE_IDS,
    RING_A_ARCHETYPES,
    RING_BC_ARCHETYPES,
    ApplicabilityBlock,
    ArchetypeBuild,
    ArchetypeRing,
    DistributionFamily,
    MissingLaw,
    MissingLawVerdict,
    PriorDistribution,
    PriorStatus,
    build_all_archetypes,
    build_archetype,
)
from aleph.virtual_cell.contracts import CellStateManifest
from aleph.virtual_cell.grammar import (
    NAMESPACE_SEPARATOR,
    CellJunction,
    CellWiring,
    CensusUnavailableError,
    ComponentDeclaration,
    ComponentExtent,
    ConnectorDeclaration,
    ConnectorInstance,
    EngineGrammar,
    PopulationWiring,
    Severity,
    SharedComponent,
    ViolationCode,
    check_colocation,
    load_engine_grammar,
    namespaced_ref,
    parse_ref,
    shared_ref,
    validate_manifest,
    validate_population,
)


@pytest.fixture(scope="module")
def grammar() -> EngineGrammar:
    """Load the canonical census once for the module."""
    return load_engine_grammar()


@pytest.fixture(scope="module")
def builds(grammar: EngineGrammar) -> tuple[ArchetypeBuild, ...]:
    """Build every registered archetype once."""
    return build_all_archetypes(grammar=grammar)


# --------------------------------------------------------------------------------------------
# The census itself
# --------------------------------------------------------------------------------------------


def test_canonical_census_matches_state_document(grammar: EngineGrammar) -> None:
    """The census read from code is 14 components / 38 connectors, as STATE.md claims."""
    assert len(grammar.components) == 14
    assert len(grammar.connectors) == 38
    assert "engine/contracts.py" in grammar.source
    assert "engine/dispatch.py" in grammar.source


def test_every_connector_has_exactly_one_declared_owner(grammar: EngineGrammar) -> None:
    """Ownership is the exact-once dispatch claim; no edge is unowned or doubly owned."""
    owners = [connector.owner for connector in grammar.connectors]
    assert len(owners) == len(grammar.connectors)
    assert set(owners) <= set(grammar.facade_components)
    for connector in grammar.connectors:
        assert grammar.owner_of(connector.name) == connector.owner


def test_owner_is_an_endpoint_except_for_the_one_composite_joint(grammar: EngineGrammar) -> None:
    """35 of 36 owners own an endpoint; the exception is the composite FA series joint.

    The engine's ``ConnectorContract`` carries NO owner field, so "the owner is one of its
    endpoints" is not a rule the engine states -- it is a relation this lane measured.  It holds
    everywhere except ``fa_actin_anchor``, whose ``focal_adhesion`` endpoint declares
    ``owns_geometry=False`` and therefore has no facade of its own.
    """
    exceptions = []
    for connector in grammar.connectors:
        claimed = set(grammar.facade_components[connector.owner])
        if not (connector.endpoints & claimed):
            exceptions.append(connector)
    assert [connector.name for connector in exceptions] == ["fa_actin_anchor"]
    assert exceptions[0].mechanical_group == "alpha2beta1_collagen_series"
    assert not grammar.component("focal_adhesion").owns_geometry


def test_environment_components_are_the_shared_ones(grammar: EngineGrammar) -> None:
    """The role that must be one shared object is ``environment``."""
    assert grammar.environment_components == {"ecm", "world_boundary", "extracellular_medium"}


def test_census_load_does_not_import_warp_or_the_engine_package() -> None:
    """The census must be readable on a CPU-only lane; ``ac.engine.__init__`` initialises Warp.

    Run in a subprocess so an unrelated test that legitimately imports Warp cannot make this pass.
    """
    script = textwrap.dedent(
        """
        import sys
        from aleph.virtual_cell.grammar import load_engine_grammar
        from aleph.virtual_cell.archetypes import build_all_archetypes
        grammar = load_engine_grammar()
        build_all_archetypes(grammar=grammar)
        assert "warp" not in sys.modules, "grammar load initialised Warp"
        assert "aleph.engine" not in sys.modules, "grammar load imported the engine package"
        print(len(grammar.components), len(grammar.connectors))
        """
    )
    completed = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        cwd=Path(__file__).resolve().parents[4],
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.split() == ["14", "38"]


def test_missing_census_source_fails_loudly_with_no_fallback(tmp_path: Path) -> None:
    """A census that cannot be located must raise, never fall back to a transcribed list."""
    with pytest.raises(CensusUnavailableError, match="not found"):
        load_engine_grammar(contracts_path=tmp_path / "absent.py")


def test_unreadable_ownership_manifest_fails_loudly(tmp_path: Path) -> None:
    """A dispatch module without the ownership manifest stops validation."""
    stub = tmp_path / "dispatch.py"
    stub.write_text("def something_else() -> None:\n    return None\n", encoding="utf-8")
    with pytest.raises(CensusUnavailableError, match="canonical_facade_claims"):
        load_engine_grammar(dispatch_path=stub)


def test_contracts_module_without_architecture_fails_loudly(tmp_path: Path) -> None:
    """A contracts module that no longer declares the architecture stops validation."""
    stub = tmp_path / "contracts.py"
    stub.write_text("VALUE = 1\n", encoding="utf-8")
    with pytest.raises(CensusUnavailableError, match="reference_cell_architecture"):
        load_engine_grammar(contracts_path=stub)


# --------------------------------------------------------------------------------------------
# Manifest-level grammar checks
# --------------------------------------------------------------------------------------------


def _manifest(
    manifest_id: str, components: tuple[str, ...], connectors: tuple[str, ...]
) -> CellStateManifest:
    return CellStateManifest(
        manifest_id=manifest_id,
        species="Homo sapiens",
        lineage="fixture",
        identity={"kind": "fixture"},
        biological_state={"state": "fixture"},
        environment={"context": "fixture"},
        components=components,
        connectors=connectors,
    )


def test_manifest_with_an_undeclared_component_is_refused(grammar: EngineGrammar) -> None:
    report = validate_manifest(
        _manifest("bogus", ("membrane", "cortex", "sarcomere"), ("membrane_erm_cortex",)),
        grammar,
    )
    assert not report.ok
    assert ViolationCode.UNKNOWN_COMPONENT in report.codes


def test_manifest_with_an_undeclared_connector_is_refused(grammar: EngineGrammar) -> None:
    report = validate_manifest(
        _manifest("bogus", ("membrane", "cortex"), ("membrane_desmosome",)), grammar
    )
    assert not report.ok
    assert ViolationCode.UNKNOWN_CONNECTOR in report.codes


def test_connector_whose_endpoint_is_absent_is_refused(grammar: EngineGrammar) -> None:
    """``membrane_erm_cortex`` needs both membrane and cortex present."""
    report = validate_manifest(_manifest("bogus", ("membrane",), ("membrane_erm_cortex",)), grammar)
    assert not report.ok
    assert ViolationCode.ENDPOINT_NOT_PRESENT in report.codes


def test_co_location_in_a_manifest_is_never_connection(grammar: EngineGrammar) -> None:
    """Two components listed together with no connector are UNCONNECTED, not coupled."""
    report = validate_manifest(_manifest("bogus", ("membrane", "cortex"), ()), grammar)
    assert not report.ok
    unconnected = {
        finding.subject
        for finding in report.refusals
        if finding.code is ViolationCode.UNCONNECTED_COMPONENT
    }
    assert unconnected == {"membrane", "cortex"}


def _synthetic_grammar(**overrides: object) -> EngineGrammar:
    """Build a two-component grammar so the connector-shape rules can be exercised directly."""
    defaults: dict[str, object] = {
        "name": "a_b_link",
        "family": "contact",
        "endpoint_a": "alpha",
        "endpoint_b": "beta",
        "scope": "inter_component",
        "kinetics": False,
        "commit_on_accept": False,
        "bidirectional": True,
        "adjoint_transfer_required": True,
        "mechanical_group": None,
        "owner": "AlphaFacade",
    }
    defaults.update(overrides)
    components = tuple(
        ComponentDeclaration(name, "surface_body", "fixture", "fixture", True, True)
        for name in ("alpha", "beta")
    )
    return EngineGrammar.from_declarations(
        components=components,
        connectors=(ConnectorDeclaration(**defaults),),  # type: ignore[arg-type]
        facade_components={"AlphaFacade": ("alpha",), "GammaFacade": ("beta",)},
        source="synthetic fixture",
    )


def test_one_sided_connector_declaration_is_refused() -> None:
    """A connector that is not bidirectional violates Newton's 3rd law and is refused."""
    grammar = _synthetic_grammar(bidirectional=False)
    report = validate_manifest(_manifest("fixture", ("alpha", "beta"), ("a_b_link",)), grammar)
    assert not report.ok
    assert ViolationCode.ONE_SIDED_CONNECTOR in report.codes


def test_connector_without_declared_adjoint_closure_is_refused() -> None:
    grammar = _synthetic_grammar(adjoint_transfer_required=False)
    report = validate_manifest(_manifest("fixture", ("alpha", "beta"), ("a_b_link",)), grammar)
    assert not report.ok
    assert ViolationCode.ADJOINT_NOT_DECLARED in report.codes


def test_owner_that_owns_no_endpoint_and_no_group_is_refused() -> None:
    """Without a mechanical group, an owner that owns neither endpoint is a wiring failure."""
    grammar = _synthetic_grammar(owner="OrphanFacade")
    report = validate_manifest(_manifest("fixture", ("alpha", "beta"), ("a_b_link",)), grammar)
    assert not report.ok
    assert ViolationCode.OWNER_MISSING in report.codes


def test_owner_outside_the_endpoints_is_permitted_only_for_a_composite_group() -> None:
    """The measured exemption: a mechanical group's composite joint is owned once, elsewhere."""
    refused = _synthetic_grammar(owner="GammaFacade", endpoint_a="alpha", endpoint_b="beta")
    # GammaFacade owns beta, which IS an endpoint, so build a genuinely orphaned case instead.
    components = (
        *refused.components,
        ComponentDeclaration("gamma", "surface_body", "fixture", "fixture", True, True),
    )
    orphaned = EngineGrammar.from_declarations(
        components=components,
        connectors=(
            ConnectorDeclaration(
                "a_b_link",
                "contact",
                "alpha",
                "beta",
                "inter_component",
                False,
                False,
                True,
                True,
                None,
                "GammaFacade",
            ),
        ),
        facade_components={"AlphaFacade": ("alpha",), "GammaFacade": ("gamma",)},
        source="synthetic fixture",
    )
    report = validate_manifest(_manifest("fixture", ("alpha", "beta"), ("a_b_link",)), orphaned)
    assert ViolationCode.OWNER_NOT_ENDPOINT in report.codes
    assert not report.ok

    composite = EngineGrammar.from_declarations(
        components=components,
        connectors=(
            ConnectorDeclaration(
                "a_b_link",
                "contact",
                "alpha",
                "beta",
                "inter_component",
                False,
                False,
                True,
                True,
                "series_joint",
                "GammaFacade",
            ),
        ),
        facade_components={"AlphaFacade": ("alpha",), "GammaFacade": ("gamma",)},
        source="synthetic fixture",
    )
    grouped = validate_manifest(_manifest("fixture", ("alpha", "beta"), ("a_b_link",)), composite)
    assert ViolationCode.OWNER_NOT_ENDPOINT_COMPOSITE in grouped.codes
    assert grouped.ok
    assert grouped.expansions[0].severity is Severity.EXPANSION


# --------------------------------------------------------------------------------------------
# Co-location
# --------------------------------------------------------------------------------------------


def test_shared_array_without_a_connector_is_reported_unconnected() -> None:
    """The project's oldest invariant: co-location in a shared array is never a connection."""
    extents = (
        ComponentExtent("cell-a::cortex", (0.0, 0.0, 0.0), (8.0, 8.0, 8.0), array_id="node_xyz"),
        ComponentExtent("cell-a::nucleus", (2.0, 2.0, 2.0), (5.0, 5.0, 5.0), array_id="node_xyz"),
    )
    findings = check_colocation(extents, ())
    assert [finding.code for finding in findings] == [ViolationCode.CO_LOCATED_WITHOUT_CONNECTOR]
    assert "co-location is never connection" in findings[0].detail


def test_overlapping_extents_joined_by_a_declared_connector_are_clean() -> None:
    extents = (
        ComponentExtent("cell-a::membrane", (0.0, 0.0, 0.0), (8.0, 8.0, 8.0)),
        ComponentExtent("cell-a::cortex", (0.1, 0.1, 0.1), (7.9, 7.9, 7.9)),
    )
    instance = ConnectorInstance(
        "erm-0",
        "membrane_erm_cortex",
        "cell-a::membrane",
        "cell-a::cortex",
        owner="SurfaceBody",
    )
    assert check_colocation(extents, (instance,)) == ()


def test_disjoint_extents_are_not_reported() -> None:
    extents = (
        ComponentExtent("cell-a::cortex", (0.0, 0.0, 0.0), (8.0, 8.0, 8.0)),
        ComponentExtent("cell-b::cortex", (20.0, 0.0, 0.0), (28.0, 8.0, 8.0)),
    )
    assert check_colocation(extents, ()) == ()


def test_duplicate_extent_reference_is_rejected() -> None:
    extent = ComponentExtent("cell-a::cortex", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    with pytest.raises(ValueError, match="only one extent"):
        check_colocation((extent, extent), ())


# --------------------------------------------------------------------------------------------
# Population namespacing, sharing, and cross-cell wiring
# --------------------------------------------------------------------------------------------


def _two_cell_population(**overrides: object) -> PopulationWiring:
    cells = (
        CellWiring("cell-a", "endothelial", ("membrane", "cortex"), ("ecm",)),
        CellWiring("cell-b", "endothelial", ("membrane", "cortex"), ("ecm",)),
    )
    shared = (SharedComponent("ecm", "collagen-world-0", ("cell-a", "cell-b")),)
    connectors = (
        ConnectorInstance(
            "erm-a", "membrane_erm_cortex", "cell-a::membrane", "cell-a::cortex", "SurfaceBody"
        ),
        ConnectorInstance(
            "erm-b", "membrane_erm_cortex", "cell-b::membrane", "cell-b::cortex", "SurfaceBody"
        ),
    )
    payload: dict[str, object] = {
        "population_id": "two-cell",
        "cells": cells,
        "shared": shared,
        "connectors": connectors,
    }
    payload.update(overrides)
    return PopulationWiring(**payload)  # type: ignore[arg-type]


def test_two_cell_instances_cannot_collide_on_a_component_name(grammar: EngineGrammar) -> None:
    """Namespacing by cell ID is what keeps two cells' membranes and cortices distinct."""
    assert namespaced_ref("cell-a", "cortex") != namespaced_ref("cell-b", "cortex")
    assert parse_ref(namespaced_ref("cell-a", "cortex")) == ("cell-a", "cortex")
    assert parse_ref(shared_ref("collagen-world-0")) == (None, "collagen-world-0")
    report = validate_population(_two_cell_population(), grammar)
    assert report.ok, [finding.to_dict() for finding in report.refusals]


def test_duplicate_cell_id_is_refused(grammar: EngineGrammar) -> None:
    population = _two_cell_population(
        cells=(
            CellWiring("cell-a", "endothelial", ("membrane", "cortex"), ("ecm",)),
            CellWiring("cell-a", "fibroblast-myofibroblast", ("membrane", "cortex"), ("ecm",)),
        ),
        shared=(SharedComponent("ecm", "collagen-world-0", ("cell-a",)),),
        connectors=(),
    )
    report = validate_population(population, grammar)
    assert not report.ok
    assert ViolationCode.DUPLICATE_CELL_ID in report.codes
    assert ViolationCode.COMPONENT_NAME_COLLISION in report.codes


def test_connector_wired_across_two_cells_cortices_is_detected_and_refused(
    grammar: EngineGrammar,
) -> None:
    """Cell-cell coupling must go through a junction, never through two cells' cortices."""
    population = _two_cell_population(
        connectors=(
            ConnectorInstance(
                "bad-weld",
                "sf_cortex_transient",
                "cell-a::cortex",
                "cell-b::cortex",
                "StressFiberActor",
            ),
        )
    )
    report = validate_population(population, grammar)
    assert not report.ok
    codes = report.codes
    assert ViolationCode.CROSS_CELL_INTRACELLULAR_CONNECTOR in codes
    detail = next(
        finding.detail
        for finding in report.refusals
        if finding.code is ViolationCode.CROSS_CELL_INTRACELLULAR_CONNECTOR
    )
    assert "explicit junction" in detail


def test_cell_cell_coupling_through_a_junction_is_accepted_and_flagged_as_a_missing_law(
    grammar: EngineGrammar,
) -> None:
    """A junction is the sanctioned route -- and the census declares no junction law at all."""
    population = _two_cell_population(
        junctions=(
            CellJunction(
                "ve-cadherin-0",
                "law.ve_cadherin_adherens_junction",
                "cell-a",
                "cell-b",
                "cell-a::membrane",
                "cell-b::membrane",
            ),
        )
    )
    report = validate_population(population, grammar)
    assert report.ok
    assert ViolationCode.JUNCTION_LAW_NOT_IN_GRAMMAR in report.codes
    expansion = report.expansions[0]
    assert expansion.severity is Severity.EXPANSION
    assert "missing component / missing law project" in expansion.expansion


def test_junction_between_a_cell_and_itself_is_refused(grammar: EngineGrammar) -> None:
    population = _two_cell_population(
        junctions=(
            CellJunction(
                "self-junction",
                "law.ve_cadherin_adherens_junction",
                "cell-a",
                "cell-a",
                "cell-a::membrane",
                "cell-a::cortex",
            ),
        )
    )
    report = validate_population(population, grammar)
    assert not report.ok
    assert ViolationCode.JUNCTION_ENDPOINTS_SAME_CELL in report.codes


def test_shared_ecm_must_be_one_object_not_duplicated_per_cell(grammar: EngineGrammar) -> None:
    """A shared environment is ONE object referenced by many cells, never a per-cell copy."""
    population = _two_cell_population(
        cells=(
            CellWiring("cell-a", "endothelial", ("membrane", "cortex", "ecm")),
            CellWiring("cell-b", "endothelial", ("membrane", "cortex", "ecm")),
        ),
        shared=(),
    )
    report = validate_population(population, grammar)
    assert not report.ok
    duplicated = {
        finding.subject
        for finding in report.refusals
        if finding.code is ViolationCode.SHARED_COMPONENT_DUPLICATED
    }
    assert duplicated == {"cell-a:ecm", "cell-b:ecm"}


def test_shared_object_must_be_declared_by_the_population(grammar: EngineGrammar) -> None:
    population = _two_cell_population(shared=())
    report = validate_population(population, grammar)
    assert not report.ok
    assert ViolationCode.SHARED_COMPONENT_UNDECLARED in report.codes


def test_per_cell_component_cannot_be_declared_shared(grammar: EngineGrammar) -> None:
    """A cortex owns per-cell state; it can never be one object several cells reference."""
    population = _two_cell_population(
        cells=(
            CellWiring("cell-a", "endothelial", ("membrane",), ("cortex", "ecm")),
            CellWiring("cell-b", "endothelial", ("membrane", "cortex"), ("ecm",)),
        ),
        connectors=(),
    )
    report = validate_population(population, grammar)
    assert not report.ok
    assert ViolationCode.PER_CELL_COMPONENT_SHARED in report.codes


def test_bare_component_endpoint_is_refused_as_unnamespaced(grammar: EngineGrammar) -> None:
    population = _two_cell_population(
        connectors=(
            ConnectorInstance(
                "erm-a", "membrane_erm_cortex", "membrane", "cell-a::cortex", "SurfaceBody"
            ),
        )
    )
    report = validate_population(population, grammar)
    assert not report.ok
    assert ViolationCode.COMPONENT_NOT_NAMESPACED in report.codes
    assert NAMESPACE_SEPARATOR in next(
        finding.detail
        for finding in report.refusals
        if finding.code is ViolationCode.COMPONENT_NOT_NAMESPACED
    )


def test_connector_instance_endpoints_must_match_the_declared_edge(grammar: EngineGrammar) -> None:
    population = _two_cell_population(
        cells=(
            CellWiring("cell-a", "endothelial", ("membrane", "cortex", "nucleus"), ("ecm",)),
            CellWiring("cell-b", "endothelial", ("membrane", "cortex"), ("ecm",)),
        ),
        connectors=(
            ConnectorInstance(
                "wrong-edge",
                "membrane_erm_cortex",
                "cell-a::membrane",
                "cell-a::nucleus",
                "SurfaceBody",
            ),
        ),
    )
    report = validate_population(population, grammar)
    assert not report.ok
    assert ViolationCode.ENDPOINT_MISMATCH in report.codes


def test_connector_instance_owner_must_match_the_canonical_manifest(
    grammar: EngineGrammar,
) -> None:
    population = _two_cell_population(
        connectors=(
            ConnectorInstance(
                "erm-a", "membrane_erm_cortex", "cell-a::membrane", "cell-a::cortex", "ECMWorld"
            ),
        )
    )
    report = validate_population(population, grammar)
    assert not report.ok
    assert ViolationCode.OWNER_MISSING in report.codes


def test_one_sided_connector_instance_is_refused(grammar: EngineGrammar) -> None:
    population = _two_cell_population(
        connectors=(
            ConnectorInstance(
                "erm-a",
                "membrane_erm_cortex",
                "cell-a::membrane",
                "cell-a::cortex",
                "SurfaceBody",
                bidirectional=False,
                adjoint_declared=False,
            ),
        )
    )
    report = validate_population(population, grammar)
    assert not report.ok
    assert ViolationCode.ONE_SIDED_CONNECTOR in report.codes
    assert ViolationCode.ADJOINT_NOT_DECLARED in report.codes


def test_report_raise_for_refusals_names_every_refusal(grammar: EngineGrammar) -> None:
    report = validate_manifest(_manifest("bogus", ("membrane", "cortex"), ()), grammar)
    with pytest.raises(ValueError, match="not expressible in the canonical grammar"):
        report.raise_for_refusals()
    payload = report.to_dict()
    assert payload["ok"] is False
    assert payload["census_source"] == grammar.source


# --------------------------------------------------------------------------------------------
# Ring A archetypes
# --------------------------------------------------------------------------------------------


def test_registry_covers_ring_a_and_ring_bc() -> None:
    assert RING_A_ARCHETYPES == (
        "amoeboid-immune",
        "endothelial",
        "fibroblast-myofibroblast",
        "junctional-epithelial",
        "stem-ipsc-like",
        "pulsatile-contractile",
    )
    assert RING_BC_ARCHETYPES == (
        "erythrocyte",
        "neuron-growth-cone",
        "ciliated-cell",
        "cardiomyocyte",
    )
    assert (*RING_A_ARCHETYPES, *RING_BC_ARCHETYPES) == ARCHETYPE_IDS


def test_every_ring_a_archetype_validates_against_the_canonical_grammar(
    builds: tuple[ArchetypeBuild, ...], grammar: EngineGrammar
) -> None:
    ring_a = [build for build in builds if build.ring is ArchetypeRing.RING_A]
    assert len(ring_a) == len(RING_A_ARCHETYPES)
    for build in ring_a:
        assert build.supported
        assert build.manifest is not None
        assert build.verdict is None
        report = validate_manifest(build.manifest, grammar)
        assert report.ok, (build.archetype_id, [f.to_dict() for f in report.refusals])
        assert set(build.manifest.components) <= grammar.component_names
        assert set(build.manifest.connectors) <= grammar.connector_names
        assert grammar.isolated_components(build.manifest.components) == ()


def test_ring_a_manifests_are_distinct_and_hash_stably(
    builds: tuple[ArchetypeBuild, ...],
) -> None:
    hashes = {build.manifest.manifest_hash for build in builds if build.manifest is not None}
    assert len(hashes) == len(RING_A_ARCHETYPES)


def test_the_programme_is_not_renarrowed_to_two_cell_lines(
    builds: tuple[ArchetypeBuild, ...],
) -> None:
    """No archetype is TARGETED at a specific line; a line is a benchmark, not the target.

    The scope is the archetype identity itself -- id, species, lineage, identity block and the
    applicability text.  A PI-gap reason may name a line as a cautionary example (the resting turgor
    is exactly that), and that is the opposite of re-narrowing.
    """
    for build in builds:
        identity = [build.archetype_id, build.applicability.to_dict()]
        if build.manifest is not None:
            identity.append(
                {
                    "species": build.manifest.species,
                    "lineage": build.manifest.lineage,
                    "identity": dict(build.manifest.identity),
                    "biological_state": dict(build.manifest.biological_state),
                    "environment": dict(build.manifest.environment),
                }
            )
        if build.verdict is not None:
            identity.append(build.verdict.to_dict())
        blob = repr(identity).lower()
        for banned in ("mcf7", "mcf-7", "mda-mb-231", "mdamb231", "hela"):
            assert banned not in blob, (build.archetype_id, banned)


def test_every_ring_a_archetype_declares_missing_laws_and_unsupported_claims(
    builds: tuple[ArchetypeBuild, ...],
) -> None:
    """Ring A is supported, not complete: the applicability block is where the archetype lives."""
    for build in builds:
        if build.ring is not ArchetypeRing.RING_A:
            continue
        assert build.applicability.missing_laws
        assert build.applicability.unsupported_claims
        assert build.manifest is not None
        assert build.manifest.missing_laws == build.applicability.law_ids
        for law in build.applicability.missing_laws:
            assert law.blocked_by


def test_component_census_barely_separates_the_ring_a_archetypes(
    builds: tuple[ArchetypeBuild, ...],
) -> None:
    """A stated finding, asserted so it cannot rot: the component set is archetype-coarse.

    Three of the six Ring A archetypes take the FULL 14-component census, so what separates an
    endothelial cell from a fibroblast from an epithelial cell is entirely in ``missing_laws``.
    """
    sets = {
        build.archetype_id: frozenset(build.manifest.components)
        for build in builds
        if build.manifest is not None
    }
    full = max(sets.values(), key=len)
    at_full = {name for name, value in sets.items() if value == full}
    assert {"endothelial", "fibroblast-myofibroblast"} <= at_full
    assert len({frozenset(value) for value in sets.values()}) < len(sets)


def test_every_ring_a_prior_is_a_distribution_slot_and_currently_an_open_gap(
    builds: tuple[ArchetypeBuild, ...],
) -> None:
    """Priors are distributions, never point values; an unsourceable prior is a PI gap."""
    for build in builds:
        assert build.priors
        for prior in build.priors:
            assert prior.unit
            if prior.status is PriorStatus.PI_GAP:
                assert prior.parameters == {}
                assert prior.pi_gap_reason
                assert prior.source is None
            else:  # pragma: no cover - no sourced prior is registered today
                assert prior.source
        assert build.pi_gap_prior_ids == tuple(prior.prior_id for prior in build.priors)
        assert build.applicability.pi_gap_prior_ids == build.pi_gap_prior_ids


# --------------------------------------------------------------------------------------------
# Ring B/C refusals
# --------------------------------------------------------------------------------------------


def test_every_ring_bc_archetype_refuses_with_a_named_missing_law(
    builds: tuple[ArchetypeBuild, ...],
) -> None:
    """A successful Ring B/C entry produces NO manifest and names exactly what is missing."""
    ring_bc = [build for build in builds if build.ring is not ArchetypeRing.RING_A]
    assert [build.archetype_id for build in ring_bc] == list(RING_BC_ARCHETYPES)
    for build in ring_bc:
        assert not build.supported
        assert build.manifest is None
        assert build.verdict is not None
        assert build.verdict.refused is True
        assert build.verdict.expansion == (
            "archetype failure -> missing component / missing law project"
        )
        assert build.verdict.missing_laws
        assert build.verdict.missing_components
        assert build.verdict.unsupported_claims
        for law in build.verdict.missing_laws:
            assert law.law_id.startswith("law.")
            assert len(law.statement) > 40
            assert len(law.blocked_by) > 40


@pytest.mark.parametrize(
    ("archetype_id", "required_law"),
    [
        ("erythrocyte", "law.membrane_in_plane_shear_elasticity"),
        ("neuron-growth-cone", "law.processive_motor_cargo_transport"),
        ("ciliated-cell", "law.axonemal_dynein_sliding"),
        ("cardiomyocyte", "law.sarcomeric_registration"),
    ],
)
def test_ring_bc_refusal_names_the_defining_missing_law(
    archetype_id: str, required_law: str
) -> None:
    """The refusal must name a SPECIFIC law, not a generic failure."""
    build = build_archetype(archetype_id)
    assert build.verdict is not None
    assert required_law in build.verdict.law_ids


def test_ring_bc_refusals_quote_the_census_representation_that_blocks_them() -> None:
    """Each ``blocked_by`` justification cites what the engine actually declares."""
    grammar = load_engine_grammar()
    representations = {component.representation for component in grammar.components}
    erythrocyte = build_archetype("erythrocyte", grammar=grammar)
    assert erythrocyte.verdict is not None
    shear = next(
        law
        for law in erythrocyte.verdict.missing_laws
        if law.law_id == "law.membrane_in_plane_shear_elasticity"
    )
    assert grammar.component("membrane").representation in shear.blocked_by
    assert "fluid Helfrich surface FEM" in representations

    cardiomyocyte = build_archetype("cardiomyocyte", grammar=grammar)
    assert cardiomyocyte.verdict is not None
    register = next(
        law
        for law in cardiomyocyte.verdict.missing_laws
        if law.law_id == "law.sarcomeric_registration"
    )
    assert grammar.component("sf_arc").representation in register.blocked_by


def test_ring_bc_does_not_need_the_engine_census_to_refuse() -> None:
    """A refusal is a statement about the grammar's shape and must not depend on loading it."""
    build = build_archetype("ciliated-cell", grammar=None)
    assert not build.supported
    assert build.verdict is not None


def test_unknown_archetype_raises() -> None:
    with pytest.raises(KeyError, match="unknown archetype"):
        build_archetype("tardigrade")


# --------------------------------------------------------------------------------------------
# Contract types
# --------------------------------------------------------------------------------------------


def test_pi_gap_prior_cannot_carry_a_number() -> None:
    with pytest.raises(ValueError, match="invented value"):
        PriorDistribution(
            prior_id="cell_radius_um",
            quantity="radius",
            unit="um",
            status=PriorStatus.PI_GAP,
            parameters={"mean": 7.5},
            pi_gap_reason="unsourced",
        )


def test_pi_gap_prior_needs_a_reason() -> None:
    with pytest.raises(ValueError, match="pi_gap_reason"):
        PriorDistribution(
            prior_id="cell_radius_um",
            quantity="radius",
            unit="um",
            status=PriorStatus.PI_GAP,
        )


def test_sourced_prior_needs_a_family_parameters_and_a_source() -> None:
    with pytest.raises(ValueError, match="needs distribution parameters"):
        PriorDistribution(
            prior_id="cell_radius_um",
            quantity="radius",
            unit="um",
            status=PriorStatus.SOURCED,
            family=DistributionFamily.LOGNORMAL,
            source="fixture",
        )
    with pytest.raises(ValueError, match="source"):
        PriorDistribution(
            prior_id="cell_radius_um",
            quantity="radius",
            unit="um",
            status=PriorStatus.SOURCED,
            family=DistributionFamily.LOGNORMAL,
            parameters={"mu": 2.0, "sigma": 0.2},
        )
    sourced = PriorDistribution(
        prior_id="cell_radius_um",
        quantity="radius",
        unit="um",
        status=PriorStatus.SOURCED,
        family=DistributionFamily.LOGNORMAL,
        parameters={"mu": 2.0, "sigma": 0.2},
        source="test fixture, not literature",
    )
    assert sourced.to_dict()["parameters"] == {"mu": 2.0, "sigma": 0.2}


def test_a_refusal_without_a_named_law_is_rejected() -> None:
    """A generic failure is not a verdict."""
    with pytest.raises(ValueError, match="without naming a missing law"):
        MissingLawVerdict(
            archetype_id="something",
            ring=ArchetypeRing.RING_B,
            missing_components=("thing",),
            missing_laws=(),
            unsupported_claims=("anything",),
        )


def test_ring_a_cannot_emit_a_refusal_verdict() -> None:
    with pytest.raises(ValueError, match="does not emit a refusal"):
        MissingLawVerdict(
            archetype_id="endothelial",
            ring=ArchetypeRing.RING_A,
            missing_components=("cell_cell_junction",),
            missing_laws=(MissingLaw("law.x", "statement", "blocked because the census says so"),),
            unsupported_claims=("nothing",),
        )


def test_a_ring_bc_build_cannot_smuggle_in_a_manifest() -> None:
    """The refusal is the deliverable; forcing a manifest out of Ring B/C is rejected."""
    verdict = MissingLawVerdict(
        archetype_id="erythrocyte",
        ring=ArchetypeRing.RING_B,
        missing_components=("spectrin_membrane_skeleton",),
        missing_laws=(MissingLaw("law.x", "statement", "blocked because the census says so"),),
        unsupported_claims=("shear modulus",),
    )
    applicability = ApplicabilityBlock(
        missing_components=verdict.missing_components,
        missing_laws=verdict.missing_laws,
        unsupported_claims=verdict.unsupported_claims,
    )
    with pytest.raises(ValueError, match="refused but still produced a manifest"):
        ArchetypeBuild(
            archetype_id="erythrocyte",
            ring=ArchetypeRing.RING_B,
            supported=False,
            priors=(),
            applicability=applicability,
            manifest=_manifest("erythrocyte", ("membrane",), ()),
            verdict=verdict,
        )


def test_a_build_that_refuses_must_carry_a_verdict() -> None:
    applicability = ApplicabilityBlock((), (), ("nothing",))
    with pytest.raises(ValueError, match="a bare failure is not an answer"):
        ArchetypeBuild(
            archetype_id="erythrocyte",
            ring=ArchetypeRing.RING_B,
            supported=False,
            priors=(),
            applicability=applicability,
        )
