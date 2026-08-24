"""The tether forward operator, wired — and the two ways the wiring can lie.

Pure host arithmetic over declared arrays. No law is evaluated, no step is taken and no device is
touched, so nothing here is a physics measurement and nothing here may be quoted as one. What it
checks is that the seam between the membrane and :mod:`aleph.observe.tether` composes the two
tension terms it claims to, and that a plateau with no scatter cannot leave wearing an error bar.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from aleph.laws.membrane_surface import GAMMA_MCA_PN_UM, erm_rupture_force
from aleph.observe.operator import ApparentQuantity
from aleph.observe.tether import equilibrium_tether_force_pn
from aleph.scripts.ac_observe_tether_native import (
    REQUIRED_PROTOCOL_FIELDS,
    REQUIRED_STATE_FIELDS,
    tether_from_declaration,
)

# The resolved physiological membrane (`laws/compartments.py::resolve_membrane` defaults, KB-3.B1)
# and the recorded preload area of the native cell.
GAMMA_MEM = 10.0
KAPPA = 0.0828
A0_UM2 = 703.4902162858084


def _declaration(areas: list[float]) -> dict:
    """A complete declaration. Every protocol field is a string because every one is required."""
    return {
        # protocol
        "temperature": "310.0 K",
        "medium": "DMEM + 10% FBS, 37 C",
        "probe_geometry": "streptavidin bead, 3 um diameter, on a tipless cantilever",
        "control_mode": "constant-speed retraction",
        "rate": "0.5 um/s",
        "sampling_cadence": "0.01 s",
        "fit_window": "2-6 um extension, declared before the pull",
        "resolution_limit": "5 pN force noise floor",
        "specimen_identity": "simulated cell, aleph native membrane",
        "batch": "unit check, not a run",
        "citation_audit": "kappa/gamma_mem KB-3.B1.2/3.B1.1; W is the GAMMA_MCA_PN_UM literal",
        "applicability_statement": "unit check of the composition seam only",
        "evidence_role": "structural check; not evidence for any physical magnitude",
        "operator_version": "test",
        # state
        "membrane_area_um2": areas,
        "rest_area_um2": A0_UM2,
        "reservoir_capacity": 0.10,
        "reservoir_capacity_provenance": "swept point; NOT sourced for this cell type",
        "gamma_mem_pn_per_um": GAMMA_MEM,
        "k_area_pn_per_um": 2.35e5,
        "tau_lysis_pn_per_um": 5.0e3,
        "bending_rigidity_pn_um": KAPPA,
        "sample_interval_s": 0.01,
        "temperature_k": 310.0,
        "retraction_speed_um_per_s": 0.5,
    }


def test_the_runtime_already_computes_this_force_as_a_rupture_threshold():
    """``erm_rupture_force`` and the operator's forward map are the same closed form.

    This is the finding that reframes the wiring: the number is not absent from the runtime, it is
    present there stripped of its protocol. If the two ever diverge, one of them changed physics
    silently and the other is still being quoted.
    """
    assert erm_rupture_force(KAPPA, GAMMA_MEM, GAMMA_MCA_PN_UM) == pytest.approx(
        equilibrium_tether_force_pn(KAPPA, GAMMA_MEM + GAMMA_MCA_PN_UM), rel=1e-15
    )


def test_constant_area_is_reported_as_an_input_echo_not_a_measurement():
    """Today's runtime: constant gamma_mem, constant kappa, constant W -> no scatter, no error bar."""
    body = tether_from_declaration(_declaration([A0_UM2] * 32))

    assert body["verdict"] == "INPUT_ECHO_NOT_A_MEASUREMENT"
    quantity = body["apparent_quantity"]
    assert quantity["standard_error"] is None, "a zero-scatter plateau must not carry an error bar"
    assert body["inputs"]["bilayer_tension_is_constant"] is True
    # And the value is exactly the closed form the two constants imply -- which is the whole point.
    assert quantity["value"] == pytest.approx(
        2.0 * math.pi * math.sqrt(2.0 * KAPPA * (GAMMA_MEM + GAMMA_MCA_PN_UM)), rel=1e-12
    )


def test_a_varying_area_makes_the_force_move_and_earns_an_error_bar():
    """Past the reservoir capacity the K_A branch engages, sigma varies, and the operator measures.

    No edit to the driver is needed for this: it is the same code path, given a state that moves.
    """
    ceiling = A0_UM2 * 1.10
    areas = list(ceiling * (1.0 + 0.002 * np.sin(np.arange(64) * 0.7)))
    body = tether_from_declaration(_declaration(areas))

    # The bare declaration leaves eleven fields NOT_RECORDED -- a confession, not a failure -- so
    # the verdict is downgraded rather than refused. That downgrade IS the manifest being wired.
    assert body["verdict"] == "MEASURED_UNDER_INCOMPLETE_PROTOCOL"
    assert body["manifest_complete"] is False
    assert "pharmacology" in body["manifest_not_recorded_fields"]
    quantity = body["apparent_quantity"]
    assert quantity["standard_error"] is not None and quantity["standard_error"] > 0.0
    assert body["inputs"]["bilayer_tension_is_constant"] is False
    # The tether radius is the operator's own falsifiability check: outside ~10-50 nm the
    # (kappa, sigma_app) pair is not a cell membrane, whatever the force came out as.
    assert 0.005 < quantity["detail"]["mean_tether_radius_um"] < 0.060


def test_every_required_field_is_actually_required():
    """Drop any one field and the driver refuses, naming it. No field may be quietly defaulted."""
    for field in REQUIRED_PROTOCOL_FIELDS + REQUIRED_STATE_FIELDS:
        declaration = _declaration([A0_UM2] * 16)
        del declaration[field]
        with pytest.raises(ValueError, match=field):
            tether_from_declaration(declaration)


def test_too_few_frames_is_a_refusal_carried_on_the_record():
    """A refusal is a result: it is written, with its reason, and does not raise."""
    body = tether_from_declaration(_declaration([A0_UM2] * 4))

    assert body["verdict"] == "REFUSED"
    assert body["refusal"]["code"] == "INSUFFICIENT_SAMPLES"
    assert body["manifest_hash"], "a refusal still records the protocol it refused under"


def test_zero_standard_error_is_refused_for_every_operator_not_just_this_one():
    """The guard lives on ``ApparentQuantity``, so all seven operators inherit it."""
    with pytest.raises(ValueError, match="infinite precision"):
        ApparentQuantity(
            name="anything",
            operator_id="test",
            value=1.0,
            units="pN",
            standard_error=0.0,
            n_effective=10.0,
            manifest_hash="deadbeef",
            analysis_chain=("mean",),
        )


def test_declaring_the_last_eleven_fields_lifts_the_downgrade_and_moves_the_hash():
    """A complete protocol earns ``MEASURED`` — and hashes differently, which is the join key."""
    ceiling = A0_UM2 * 1.10
    areas = list(ceiling * (1.0 + 0.002 * np.sin(np.arange(64) * 0.7)))
    partial = tether_from_declaration(_declaration(areas))

    declaration = _declaration(areas)
    declaration.update(
        {
            "probe_location": "cell apex",
            "probe_direction": "outward normal",
            "adhesion_state": "biotin-streptavidin",
            "dwell": "1.0 s contact dwell",
            "loading_history": "first pull on this cell",
            "preconditioning": "none",
            "calibration": "forward operator; no transducer, so no spring constant exists to state",
            "osmotic_state": "isotonic",
            "substrate": "suspended, no substrate",
            "confinement": "none",
            "pharmacology": "vehicle only",
        }
    )
    body = tether_from_declaration(declaration)

    assert body["verdict"] == "MEASURED"
    assert body["manifest_complete"] is True
    assert body["manifest_not_recorded_fields"] == []
    assert body["manifest_hash"] != partial["manifest_hash"], (
        "two protocols that differ must not share a digest, or a number joins to the wrong record"
    )
