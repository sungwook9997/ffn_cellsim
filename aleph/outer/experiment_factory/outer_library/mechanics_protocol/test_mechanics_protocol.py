from __future__ import annotations

import json
from pathlib import Path

import jsonschema
import numpy as np

from .infer import route_protocol
from .train import _method_superclass, _probe_scaling, _wire_roles

HERE = Path(__file__).parent
RESULT = HERE.parents[3] / "outputs/outer/mechanics_protocol/model/mechanics_protocol_report.json"


def test_probe_scaling_recovers_quadratic_exponent() -> None:
    length = np.tile(np.asarray([1.0, 2.0, 3.0, 5.0]), 3)
    labels = np.repeat(np.arange(3), 4)
    groups = np.tile(np.asarray(["a", "a", "b", "b"]), 3)
    report = _probe_scaling(length, 7.0 * length ** 2, labels, groups)
    for cell_line in ("MCF-10A", "MCF-7", "MDA-MB-231"):
        assert np.isclose(report[cell_line]["raw_wire_level_log_log_exponent"], 2.0)
        assert np.isclose(report[cell_line]["predicted_eta_fold_1_to_5_um_from_raw_fit"], 25.0)


def test_date_split_is_disjoint_per_cell_line() -> None:
    labels = np.repeat(np.arange(3), 6)
    groups = np.tile(np.asarray(["a", "a", "b", "b", "c", "c"]), 3)
    roles = _wire_roles(labels, groups)
    for label in range(3):
        split_sets = [set(groups[(labels == label) & (roles == role)]) for role in ("fit", "calibration", "test")]
        assert all(not split_sets[i] & split_sets[j] for i in range(3) for j in range(i + 1, 3))


def test_protocol_superclasses_do_not_collapse_all_afm_or_magnetic_methods() -> None:
    assert _method_superclass("magnetic_rotational_spectroscopy") == _method_superclass("magnetic_bead_microrheometry")
    assert _method_superclass("afm_force_relaxation") == "local_afm"
    assert _method_superclass("dynamic_afm_parallel_plate") == "whole_cell_afm_confinement"
    assert _method_superclass("continuous_micropipette_aspiration") == "whole_cell_suction"


def test_committed_report_is_fail_closed() -> None:
    report = json.loads(RESULT.read_text(encoding="utf-8"))
    assert report["aleph_authority"] == "none"
    assert report["observable_evidence_eligible"] is False
    assert report["may_select_aleph_parameter"] is False
    assert report["architecture"]["gpu_used"] is False
    assert report["protocol_router"]["external_holdout"] is False
    assert "insufficient_external_holdout" in report["protocol_router"]["status"]
    assert report["wire_head"]["experimental_date_disjoint"]["split_group_overlap_count"] == 0


def test_catalog_never_pools_units_or_expands_aggregate_rows() -> None:
    catalog = json.loads((HERE / "protocol_catalog.json").read_text(encoding="utf-8"))
    schema = json.loads((HERE / "protocol_catalog.schema.json").read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(catalog)
    assert catalog["aleph_authority"] == "none"
    assert len(catalog["records"]) == 19
    assert len({record["source_id"] for record in catalog["records"]}) == 9
    assert set(catalog["sources"]) == {record["source_id"] for record in catalog["records"]}
    assert sum(source["repo_source_audit"] == "OK" for source in catalog["sources"].values()) == 5
    mcf7_viscosity_scopes = {
        record["mechanical_scope"] for record in catalog["records"]
        if record["cell_line"] == "MCF-7" and record["unit"] == "Pa s"
    }
    assert mcf7_viscosity_scopes == {
        "effective_cytoplasm_viscoelasticity", "cytoplasmic_solvent_microviscosity"
    }
    rules = " ".join(catalog["non_pooling_contract"])
    assert "never normalized into one target" in rules
    assert "never expanded into synthetic cells" in rules


def test_router_refuses_unseen_scope_instead_of_nearest_neighbor_transfer() -> None:
    checkpoint = RESULT.parent / "mechanics_protocol_router.npz"
    request = {
        "cell_state": "adherent_interphase_unsynchronized",
        "observable": "invented_global_mechanics_score",
        "probe_geometry": "internalized_wire",
        "probe_size_um": 3.0,
        "time_scale_s": 1.0,
        "unit": "arbitrary",
    }
    result = route_protocol(checkpoint, HERE / "protocol_catalog.json", request)
    assert result["refused"] is True
    assert result["nearest_scope_transfer_attempted"] is False
