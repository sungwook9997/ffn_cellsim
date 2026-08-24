"""Structural gates for the component-first cell-engine architecture."""

from __future__ import annotations

import pytest

from aleph.engine import (
    EVIDENCE_ORDER,
    LADDER_FLOOR,
    OSMOTIC_ENVELOPE,
    OSMOTIC_LOADED_SURFACE,
    CellArchitecture,
    ComponentContract,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
    ConnectorScope,
    EvidenceLabel,
    EvidenceRung,
    GateVerdict,
    QuantitativeClaim,
    VoidCeiling,
    classify_gate,
    reference_cell_architecture,
    rung_rank,
)


def test_membrane_and_cortex_are_distinct_but_erm_connected() -> None:
    architecture = reference_cell_architecture()
    assert architecture.component("membrane") is not architecture.component("cortex")
    assert "cortex" in architecture.neighbors("membrane")
    connector = next(item for item in architecture.connectors if item.name == "membrane_erm_cortex")
    assert connector.family is ConnectorFamily.ERM
    assert connector.kinetics and connector.commit_on_accept
    assert "cytosol" in architecture.neighbors("membrane")
    boundary = next(
        item for item in architecture.connectors if item.name == "membrane_cytosol_boundary"
    )
    assert boundary.family is ConnectorFamily.FLUID_BOUNDARY


def test_load_path_reaches_environment_and_nucleus() -> None:
    architecture = reference_cell_architecture()
    assert architecture.force_path("sf_arc", "ecm") == (
        "sf_arc", "focal_adhesion", "ecm"
    )
    assert architecture.force_path("microtubule", "nucleus") == ("microtubule", "nucleus")
    assert architecture.force_path("intermediate_filament", "nucleus") == (
        "intermediate_filament", "nucleus"
    )


def test_fluid_boundaries_and_skeleton_transfers_are_explicit() -> None:
    architecture = reference_cell_architecture()
    connectors = {connector.name: connector for connector in architecture.connectors}
    for boundary_name in ("membrane_cytosol_boundary", "nucleus_cytosol_boundary"):
        boundary = connectors[boundary_name]
        assert boundary.family is ConnectorFamily.FLUID_BOUNDARY
        assert boundary.bidirectional and boundary.adjoint_transfer_required
    for transfer_name in (
        "surface_porous_transfer",
        "sf_cytosol_transfer",
        "mt_cytosol_transfer",
        "if_cytosol_transfer",
    ):
        transfer = connectors[transfer_name]
        assert transfer.family is ConnectorFamily.IMMERSED_TRANSFER
        assert transfer.bidirectional and transfer.adjoint_transfer_required


def test_microtubule_cortical_capture_is_a_dynamic_motor_joint() -> None:
    architecture = reference_cell_architecture()
    capture = next(
        connector for connector in architecture.connectors
        if connector.name == "mt_cortex_capture"
    )
    assert capture.family is ConnectorFamily.MOTOR
    assert capture.kinetics and capture.commit_on_accept
    assert capture.chemistry_card == "microtubule_cortical_dynein"


def test_protrusions_are_explicit_distinct_actors_with_dynamic_paths() -> None:
    architecture = reference_cell_architecture()
    lamellipodium = architecture.component("lamellipodium")
    filopodium = architecture.component("filopodium")
    assert lamellipodium is not filopodium
    assert "explicit" in lamellipodium.representation
    assert "explicit" in filopodium.representation
    assert architecture.force_path("lamellipodium", "ecm") == (
        "lamellipodium", "focal_adhesion", "ecm"
    )
    assert architecture.force_path("filopodium", "ecm") == (
        "filopodium", "focal_adhesion", "ecm"
    )


def test_ecm_has_collagen_specific_crosslinks_contact_and_far_field_closure() -> None:
    architecture = reference_cell_architecture()
    connectors = {connector.name: connector for connector in architecture.connectors}
    crosslink = connectors["ecm_crosslink"]
    assert crosslink.family is ConnectorFamily.FIBER_CROSSLINK
    assert crosslink.generation_required
    assert crosslink.remap_on_accept
    assert crosslink.blocks_sleep_refine
    boundary = connectors["ecm_far_field_anchor"]
    assert boundary.family is ConnectorFamily.ENVIRONMENT_BOUNDARY
    assert {boundary.component_a, boundary.component_b} == {"ecm", "world_boundary"}
    assert connectors["membrane_ecm_contact"].family is ConnectorFamily.CONTACT


def test_nmii_is_head_resolved_and_connects_to_every_actin_actor() -> None:
    architecture = reference_cell_architecture()
    nmii = architecture.component("nmii")
    assert "individual heads" in nmii.representation
    connectors = {connector.name: connector for connector in architecture.connectors}
    for name in (
        "nmii_sf_motor",
        "nmii_cortex_motor",
        "nmii_lamellipodium_motor",
        "nmii_filopodium_motor",
    ):
        motor = connectors[name]
        assert motor.family is ConnectorFamily.MOTOR
        assert motor.kinetics and motor.commit_on_accept
        assert motor.chemistry_card == "nmii_head_actin_hill_bell"
        assert motor.generation_required and motor.blocks_sleep_refine


def test_all_linc_edges_lock_live_endpoint_generations() -> None:
    architecture = reference_cell_architecture()
    linc_edges = tuple(
        connector for connector in architecture.connectors
        if connector.family is ConnectorFamily.LINC
    )
    assert {connector.name for connector in linc_edges} == {
        "actin_cap_linc", "mt_nucleus_linc", "if_nucleus_linc"
    }
    assert all(connector.generation_required for connector in linc_edges)
    assert all(connector.remap_on_accept for connector in linc_edges)
    assert all(connector.blocks_sleep_refine for connector in linc_edges)


def test_reference_cell_has_one_connected_fourteen_actor_graph() -> None:
    """The membership set and the count are deliberate TRIPWIRES, not redundancy.

    They are spelled out so that changing the composition cannot happen silently — a contract change
    has to break this test and be confirmed by a human. 13/32 -> 13/35 (PI 2026-07-25) -> 14/36 -> 14/38 (PI 2026-08-09, nascent-FA option (a):
    each nascent adhesion needs an ECM-side clutch to be a SERIES, because `focal_adhesion` owns no
    geometry and an actin anchor alone has nothing to pull against)
    (PI D5-A, 2026-07-28: `extracellular_medium`, declared only, no exterior solve until T10).
    """
    architecture = reference_cell_architecture()
    expected = {
        "membrane",
        "cortex",
        "cytosol",
        "nucleus",
        "sf_arc",
        "focal_adhesion",
        "microtubule",
        "intermediate_filament",
        "lamellipodium",
        "filopodium",
        "nmii",
        "ecm",
        "world_boundary",
        "extracellular_medium",
    }
    assert {component.name for component in architecture.components} == expected
    assert len(architecture.connectors) == 38
    assert all(
        architecture.force_path(component, "world_boundary") is not None
        for component in expected
    )


def test_concatenation_without_connector_is_not_connectivity() -> None:
    architecture = CellArchitecture(
        components=(
            ComponentContract("a", ComponentRole.SURFACE_BODY, "mesh", "surface"),
            ComponentContract("b", ComponentRole.ACTIVE_LOAD_PATH, "cable", "graph"),
        ),
        connectors=(),
    )
    assert architecture.force_path("a", "b") is None


def test_unknown_connector_endpoint_is_rejected() -> None:
    component = ComponentContract("a", ComponentRole.SURFACE_BODY, "mesh", "surface")
    connector = ConnectorContract(
        "bad", ConnectorFamily.ERM, "a", "missing", kinetics=True, commit_on_accept=True
    )
    with pytest.raises(ValueError, match="unknown components"):
        CellArchitecture((component,), (connector,))


def test_kinetic_connector_cannot_mutate_outside_acceptance() -> None:
    with pytest.raises(ValueError, match="accepted physical step"):
        ConnectorContract(
            "bad", ConnectorFamily.FA_CLUTCH, "fa", "ecm", kinetics=True, commit_on_accept=False
        )


def test_one_way_load_is_rejected() -> None:
    with pytest.raises(ValueError, match="one-way"):
        ConnectorContract(
            "bad", ConnectorFamily.LINC, "nucleus", "mt", kinetics=False, commit_on_accept=False,
            bidirectional=False,
        )


def test_internal_connector_stays_inside_component() -> None:
    connector = ConnectorContract(
        "dorsal_arc", ConnectorFamily.TRANSIENT_ACTIN, "sf", "sf", True, True,
        scope=ConnectorScope.INTERNAL,
    )
    assert connector.component_a == connector.component_b == "sf"
    with pytest.raises(ValueError, match="stay inside"):
        ConnectorContract(
            "bad", ConnectorFamily.TRANSIENT_ACTIN, "sf", "arc", True, True,
            scope=ConnectorScope.INTERNAL,
        )


def test_topology_remap_or_sleep_lock_requires_generation_checks() -> None:
    with pytest.raises(ValueError, match="without generation checks"):
        ConnectorContract(
            "bad", ConnectorFamily.CONTACT, "membrane", "ecm", False, False,
            remap_on_accept=True,
        )


def test_fa_semantic_edges_resolve_to_one_series_mechanical_group() -> None:
    architecture = reference_cell_architecture()
    group = architecture.mechanical_group("alpha2beta1_collagen_series")
    assert {connector.name for connector in group} == {
        "fa_actin_anchor", "integrin_collagen_clutch"
    }
    assert {connector.chemistry_card for connector in group} == {
        "actin_talin_integrin", "alpha2beta1_collagen"
    }


# ── D1-B: evidence is TWO axes — a structural rung and an orthogonal magnitude status ─────────────


def test_evidence_ladder_is_owned_by_the_engine_not_a_viz_script() -> None:
    """The rung vocabulary must be importable from the contract module.

    It previously existed only inside ``scripts/ac_architecture_dashboard.py``, which is why drivers
    hand-typed rungs as free strings and one stamped a rung it had not earned.
    """
    assert EVIDENCE_ORDER[LADDER_FLOOR] is EvidenceRung.CONTRACTED
    assert [rung.value for rung in EVIDENCE_ORDER[LADDER_FLOOR:]] == [
        "CONTRACTED", "SEAMED", "KERNEL_BOUND", "CUDA_UNIT",
        "CONNECTED", "NATIVE", "OPTIMISED", "PRODUCTION",
    ]
    # the two below-ladder markers prove wiring, not physics, so they sort under the ladder
    assert rung_rank(EvidenceRung.SCHEMA_PROOF) < LADDER_FLOOR
    assert rung_rank(EvidenceRung.CENSUS_WIRED) < LADDER_FLOOR
    assert rung_rank("NOT_A_RUNG") == -1


def test_rung_and_magnitude_are_independent_axes() -> None:
    """A high rung with a blocked magnitude is the SUCCESS case, not a contradiction (PI D1-B)."""
    label = EvidenceLabel(
        rung=EvidenceRung.CONNECTED,
        quantitative=QuantitativeClaim.BLOCKED,
        basis="connector dispatched; every motor parameter is a PI-GAP",
    )
    assert label.on_ladder
    assert label.as_artifact_fields() == {
        "evidence": "CONNECTED",
        "quantitative_claim_status": "BLOCKED",
        "evidence_basis": "connector dispatched; every motor parameter is a PI-GAP",
    }
    # BLOCKED is the default: an artifact that says nothing about magnitude claims nothing
    assert EvidenceLabel(EvidenceRung.CUDA_UNIT, basis="6 CUDA-gated tests").quantitative is (
        QuantitativeClaim.BLOCKED
    )


def test_a_rung_cannot_be_asserted_without_naming_what_measured_it() -> None:
    with pytest.raises(ValueError, match="needs a basis"):
        EvidenceLabel(rung=EvidenceRung.NATIVE, basis="   ")
    with pytest.raises(TypeError, match="not a free string"):
        EvidenceLabel(rung="CONNECTED", basis="a string rung is what D1-B removes")  # type: ignore[arg-type]


def test_below_ladder_marker_cannot_confirm_a_magnitude() -> None:
    """Definitional, not a new gate: a topology-only dump proves no physics."""
    with pytest.raises(ValueError, match="cannot carry a CONFIRMED magnitude"):
        EvidenceLabel(
            rung=EvidenceRung.SCHEMA_PROOF,
            quantitative=QuantitativeClaim.CONFIRMED,
            basis="synthetic topology graph",
        )


# ── D3-C: turgor stays on the cortex, and reaches the membrane only through declared connectors ──


def test_osmotic_envelope_loads_the_membrane_through_declared_connectors() -> None:
    architecture = reference_cell_architecture()
    assert OSMOTIC_ENVELOPE == "cortex"          # Variant A ratified (PI D3-C, 2026-07-28)
    assert OSMOTIC_LOADED_SURFACE == "membrane"
    path = architecture.osmotic_envelope_load_path()
    names = {connector.name for connector in path}
    assert names == {"membrane_erm_cortex", "membrane_cortex_contact"}
    # Both edges are load-bearing and neither is redundant: an ERM tether is force-free in
    # compression, so only the CONTACT edge can carry an outward pressure.
    families = {connector.name: connector.family for connector in path}
    assert families["membrane_erm_cortex"] is ConnectorFamily.ERM
    assert families["membrane_cortex_contact"] is ConnectorFamily.CONTACT
    # every edge on the path is bidirectional + adjoint, so Pi_0 cannot be transmitted one-way
    assert all(connector.bidirectional for connector in path)
    assert all(connector.adjoint_transfer_required for connector in path)


def test_an_envelope_with_no_connector_to_the_loaded_surface_is_rejected() -> None:
    """The structural form of the Variant A/B defect: pressure that cannot reach the surface."""
    architecture = reference_cell_architecture()
    with pytest.raises(ValueError, match="no declared connector"):
        architecture.osmotic_envelope_load_path(envelope="nucleus", loaded_surface="ecm")
    with pytest.raises(KeyError):
        architecture.osmotic_envelope_load_path(envelope="not_a_component")


# ── D7: VOID is a third gate verdict, and it OVERRIDES pass ───────────────────────────────────────


def test_void_overrides_pass_which_is_the_whole_point() -> None:
    """The withdrawn 1.92 pN reproduced: a gate that "passed" on an unconverged solve.

    Its residual/tension was 15.5%. Under a declared 5% ceiling that run is VOID, not PASS — it
    measured nothing, so neither a pass nor a fail is available from it.
    """
    ceiling = VoidCeiling(0.05, "residual above 5% of signal leaves no separable physics")
    assert classify_gate(passed=True, residual=15.5, signal=100.0, ceiling=ceiling) is GateVerdict.VOID
    assert classify_gate(passed=False, residual=15.5, signal=100.0, ceiling=ceiling) is GateVerdict.VOID
    # below the ceiling the caller's own criterion decides, unchanged
    assert classify_gate(passed=True, residual=1.0, signal=100.0, ceiling=ceiling) is GateVerdict.PASS
    assert classify_gate(passed=False, residual=1.0, signal=100.0, ceiling=ceiling) is GateVerdict.FAIL
    # exactly at the ceiling is still interpretable; strictly above is not
    assert classify_gate(passed=True, residual=5.0, signal=100.0, ceiling=ceiling) is GateVerdict.PASS
    assert classify_gate(passed=True, residual=5.001, signal=100.0, ceiling=ceiling) is GateVerdict.VOID


def test_a_ceiling_cannot_be_a_bare_number_and_cannot_void_nothing() -> None:
    with pytest.raises(ValueError, match="bare magic number"):
        VoidCeiling(0.05, "   ")
    # a ceiling at or above 1 would call a residual equal to its own signal "measurable"
    with pytest.raises(ValueError, match="voids nothing"):
        VoidCeiling(1.5, "too permissive to mean anything")
    with pytest.raises(ValueError, match="fraction in"):
        VoidCeiling(0.0, "a zero ceiling voids every run")


def test_zero_signal_is_void_never_pass() -> None:
    """A ratio criterion cannot be evaluated without a signal.

    Gates whose observable is legitimately zero — `nmii_sf_motor`'s "tension is exactly 0 while no
    head is bound" is a real one — assert an absolute zero and must not route through this classifier.
    """
    ceiling = VoidCeiling(0.05, "declared before the run")
    assert classify_gate(passed=True, residual=0.0, signal=0.0, ceiling=ceiling) is GateVerdict.VOID

    for bad in ((-1.0, 1.0), (1.0, -1.0), (float("nan"), 1.0), (1.0, float("inf"))):
        with pytest.raises(ValueError):
            classify_gate(passed=True, residual=bad[0], signal=bad[1], ceiling=ceiling)


def test_gate_verdict_and_quantitative_claim_are_different_questions() -> None:
    """Easy to conflate, so pinned: one judges the RUN, the other judges the NUMBER."""
    assert set(GateVerdict) == {GateVerdict.PASS, GateVerdict.FAIL, GateVerdict.VOID}
    assert "VOID" not in {claim.value for claim in QuantitativeClaim}
    assert "BLOCKED" not in {verdict.value for verdict in GateVerdict}


def test_environment_boundaries_are_inert_in_position_not_in_force() -> None:
    """RESOLVED 2026-08-09 (PI): the far-field anchor is a CONNECTOR, not a clamp.

    The proposal `the-far-field-anchor-is-a-clamp-and-the-registry-calls-it-a-connector` reads an inert
    far side as proof the edge is a boundary condition.  That conflates two different fields:
    ``dynamically_evolving`` governs POSITION and ``adjoint_transfer_required`` governs FORCE.  The
    environment frame does not move AND its reaction is accumulated into a reservoir rather than
    discarded — ``medium_exterior`` already calls the far side "the reservoir's half of the adjoint
    pair".  A clamp would discard the reaction, and then the balance ledger's two channels would agree
    by OMISSION, which is the failure a one-sided scatter is supposed to be caught by.

    This test is what the resolution rests on: flip it and the ratification is withdrawn.
    """
    architecture = reference_cell_architecture()
    environment_edges = [
        connector for connector in architecture.connectors
        if connector.family is ConnectorFamily.ENVIRONMENT_BOUNDARY
    ]
    assert {c.name for c in environment_edges} == {"ecm_far_field_anchor", "membrane_medium_traction"}
    for connector in environment_edges:
        far_side = architecture.component(connector.component_b)
        assert far_side.role is ComponentRole.ENVIRONMENT
        assert far_side.dynamically_evolving is False, "the environment frame does not move"
        assert far_side.owns_geometry is False
        assert connector.bidirectional is True, "force flows both ways even where position does not"
        assert connector.adjoint_transfer_required is True, "the reaction is accumulated, not discarded"
        assert connector.kinetics is False, "an environment boundary forms and breaks no bonds"
