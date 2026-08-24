"""Structural preflight for the suspended-round AFM cortical-tension experiment."""

from __future__ import annotations

from dataclasses import replace

import pytest

from aleph.engine.afm_cortical_tension import (
    AFMContactCalibration,
    AFMSweepDeclaration,
    C2Acceptance,
    CellGeometry,
    Holder,
    OsmoticSetpoint,
    Pi0Evidence,
    ResultClass,
    SphericalIndenterProtocol,
    afm_experiment_architecture,
)
from aleph.engine.connector_joints import force_kernel_for
from aleph.engine.contracts import ConnectorFamily, reference_cell_architecture


def _protocol() -> SphericalIndenterProtocol:
    return SphericalIndenterProtocol(
        target_protocol="example direct suspended-cell AFM protocol",
        target_cell_line="target-line",
        geometry=CellGeometry.SUSPENDED_ROUND,
        holder=Holder.EXTERIOR_STOKES,
        indenter_radius_um=5.0,
        exterior_viscosity_pa_s=1.0e-3,
        apparatus_source="example protocol methods section",
    )


def _direct_pi0() -> OsmoticSetpoint:
    return OsmoticSetpoint(
        value_pa=50.0,
        evidence=Pi0Evidence.DIRECT_TARGET_PROTOCOL,
        source="example target-protocol volume/osmotic measurement",
    )


def _contact() -> AFMContactCalibration:
    return AFMContactCalibration(
        contact_distance_um=0.05,
        stiffness_pn_per_um=100.0,
        provenance="example grid-refinement declaration, not a production value",
    )


def _c2(passed: bool = True) -> C2Acceptance:
    return C2Acceptance(
        build_commit="deadbeef",
        inner_converged=passed,
        outer_accepted=passed,
        committed_time_s=0.01 if passed else 0.0,
        max_projected_force_pn=0.2 if passed else 0.4,
        projected_force_tolerance_pn=0.21066859316597425,
        membrane_subdivisions=8,
        full_native_population=True,
        full_compartments=True,
        accepted_seeds=(0, 1, 2),
    )


def test_reference_census_is_unchanged_and_apparatus_adds_only_one_component_and_edge() -> None:
    reference = reference_cell_architecture()
    experiment = afm_experiment_architecture(reference)
    assert len(experiment.components) == len(reference.components) + 1
    assert len(experiment.connectors) == len(reference.connectors) + 1
    assert len(reference.components) == 14
    assert len(reference.connectors) == 38


def test_indenter_membrane_is_the_sixth_contact_and_uses_the_existing_unilateral_law() -> None:
    architecture = afm_experiment_architecture()
    contacts = [c for c in architecture.connectors if c.family is ConnectorFamily.CONTACT]
    assert len(contacts) == 6
    indenter = next(c for c in contacts if c.name == "indenter_membrane_contact")
    assert {indenter.component_a, indenter.component_b} == {"afm_indenter", "membrane"}
    assert indenter.bidirectional and indenter.adjoint_transfer_required
    assert not indenter.remap_on_accept and not indenter.generation_required
    assert force_kernel_for(indenter).__name__ == "_unilateral_contact_force_kernel"


def test_adherent_cell_and_fa_holder_are_different_experiments() -> None:
    kwargs = dict(
        target_protocol="x",
        target_cell_line="x",
        indenter_radius_um=5.0,
        exterior_viscosity_pa_s=1.0e-3,
        apparatus_source="x",
    )
    with pytest.raises(ValueError, match="SF/FA-dominated"):
        SphericalIndenterProtocol(
            geometry=CellGeometry.ADHERENT_SPREAD, holder=Holder.EXTERIOR_STOKES, **kwargs
        )
    with pytest.raises(ValueError, match="not focal adhesions"):
        SphericalIndenterProtocol(
            geometry=CellGeometry.SUSPENDED_ROUND, holder=Holder.FOCAL_ADHESION, **kwargs
        )


def test_hela_40_pa_proxy_is_mechanism_only_and_blocks_quantitative_enumeration() -> None:
    proxy = OsmoticSetpoint(
        value_pa=40.0,
        evidence=Pi0Evidence.CROSS_PROTOCOL_PROXY,
        source="Fischer-Friedrich 2014 HeLa proxy used outside its exact target protocol",
    )
    assert proxy.result_class is ResultClass.MECHANISM_DEMO
    sweep = AFMSweepDeclaration(
        protocol=_protocol(), pi0=proxy, contact=_contact(), c2=_c2(),
        speeds_um_s=(0.1,), depths_um=(0.0, 0.2), seeds=(0, 1, 2),
    )
    assert sweep.reaction_channel == 0
    assert any(blocker.startswith("PI0_PROXY") for blocker in sweep.blockers())
    with pytest.raises(RuntimeError, match="PI0_PROXY"):
        sweep.points()
    assert len(sweep.points(require_quantitative=False)) == 6


def test_contact_and_c2_are_hard_preflight_blockers() -> None:
    sweep = AFMSweepDeclaration(
        protocol=_protocol(), pi0=_direct_pi0(), contact=None, c2=_c2(False),
        speeds_um_s=(0.1,), depths_um=(0.1,), seeds=(0, 1, 2),
    )
    assert sweep.blockers() == (
        "ACTIVE_CORTEX_STATE_MISSING: no post-induction, equilibrated, dt-converged active-cortex state "
        "handoff exists; STATE (c) 3/17 therefore still block cortical-tension magnitudes",
        "CONTACT_CALIBRATION_MISSING: rest_um and stiffness are undeclared",
        "C2_FAILED: require accepted full-native/full-compartment mechanics, maxPF<=the unchanged predicate, "
        "membrane subdivision>=8, committed time, and >=3 accepted seeds",
    )
    with pytest.raises(RuntimeError, match="CONTACT_CALIBRATION_MISSING.*C2_FAILED"):
        sweep.points()


@pytest.mark.parametrize(
    "evidence",
    (
        {"membrane_subdivisions": 6},
        {"max_projected_force_pn": 0.22},
        {"full_native_population": False},
        {"full_compartments": False},
        {"accepted_seeds": (0,)},
    ),
)
def test_c2_cannot_pass_on_reduced_or_under_replicated_evidence(evidence: dict[str, object]) -> None:
    assert not replace(_c2(), **evidence).passed


def test_three_distinct_seeds_are_required_and_points_are_seed_replicates() -> None:
    with pytest.raises(ValueError, match="three distinct"):
        AFMSweepDeclaration(
            protocol=_protocol(), pi0=_direct_pi0(), contact=_contact(), c2=_c2(),
            speeds_um_s=(0.1,), depths_um=(0.1,), seeds=(0, 0, 1),
        )
    sweep = AFMSweepDeclaration(
        protocol=_protocol(), pi0=_direct_pi0(), contact=_contact(), c2=_c2(),
        speeds_um_s=(0.1, 0.2), depths_um=(0.0, 0.1), seeds=(3, 5, 8),
    )
    points = sweep.points(require_quantitative=False)
    assert len(points) == 12
    assert {point.seed for point in points} == {3, 5, 8}
