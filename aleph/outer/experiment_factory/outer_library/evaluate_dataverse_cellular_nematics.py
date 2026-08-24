#!/usr/bin/env python3
"""Evaluate cellular-nematic force observations as a pre-sweep constraint."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np


DATASET = "dataverse-data2772-nih3t3-cellular-nematics"


def _paired(rows: list[dict[str, object]], observable: str) -> dict[str, object]:
    grouped: dict[str, dict[str, float]] = defaultdict(dict)
    for row in rows:
        if row["observable"] == observable:
            grouped[str(row["technical_replicate"])][str(row["condition"])] = float(row["value"])
    if len(grouped) != 12 or any(len(pair) != 2 for pair in grouped.values()):
        raise ValueError(f"paired Figure S9 structure changed: {observable}")
    before = np.asarray([grouped[key]["before_blebbistatin"] for key in sorted(grouped)])
    after = np.asarray([
        grouped[key]["after_blebbistatin_10uM_30min"] for key in sorted(grouped)
    ])
    difference = after - before
    rng = np.random.default_rng(2772)
    boot = np.empty(20_000)
    for index in range(len(boot)):
        selected = rng.integers(0, len(difference), len(difference))
        boot[index] = float(np.mean(difference[selected]))
    return {
        "before_FOV_mean": float(np.mean(before)),
        "after_FOV_mean": float(np.mean(after)),
        "pooled_after_to_before_ratio": float(np.mean(after) / np.mean(before)),
        "mean_paired_change_after_minus_before": float(np.mean(difference)),
        "technical_FOV_bootstrap_95ci": [
            float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))
        ],
        "FOVs_with_decrease": int(np.sum(difference < 0)),
        "FOV_count": len(difference),
        "exact_technical_sign_probability_all_same_direction": float(2.0 ** -12),
        "inference_warning": (
            "fields of view are nested technical replicates from one dated culture/batch; "
            "the bootstrap and sign probability are descriptive, not biological-replicate inference"
        ),
    }


def evaluate(manifest: Path, manifest_report: Path) -> dict[str, object]:
    rows = [json.loads(line) for line in manifest.read_text(encoding="utf-8").splitlines()]
    audit = json.loads(manifest_report.read_text(encoding="utf-8"))
    if len(rows) != 372 or {row["dataset_id"] for row in rows} != {DATASET}:
        raise ValueError("cellular-nematic manifest identity/cardinality changed")
    if {row["aleph_authority"] for row in rows} != {"none"}:
        raise ValueError("external cellular-nematic data acquired Aleph authority")
    if audit["tfm_msm"]["biological_culture_or_batch_count"] != 1:
        raise ValueError("Figure S9 biological hierarchy was inflated")
    contrasts = {
        name: _paired(rows, name) for name in (
            "mean_traction_magnitude", "median_traction_magnitude",
            "rms_traction_magnitude", "mean_maximum_principal_tension",
            "median_maximum_principal_tension",
        )
    }
    direction_reproduced = all(
        item["FOVs_with_decrease"] == 12 and item["pooled_after_to_before_ratio"] < 0.25
        for item in contrasts.values()
    )
    return {
        "schema": "aleph.outer_library.dataverse_data2772_evidence.v1",
        "dataset_id": DATASET,
        "lab_group": "ibec-trepat-guillamat-cellular-nematics",
        "task": "myosin_II_inhibition_force_and_tension_direction",
        "manifest": audit,
        "figure_S9_paired_technical_FOV_contrasts": contrasts,
        "published_direction_reproduced_from_source_tables": direction_reproduced,
        "source_role_discrepancy": {
            "detected": True,
            "t0_t1_reproduce_published_before_after_histograms": True,
            "t2_admitted": False,
            "resolution": "withhold t2 until authors or source metadata disambiguate its role",
        },
        "piv_representation": {
            "single_sequence": True,
            "frames": audit["piv"]["frame_count"],
            "valid_and_documented_interpolated_vectors": audit["piv"][
                "vectors_admitted_total"
            ],
            "external_holdout_eligible": False,
        },
        "design_audit": {
            "fields_treated_as_biological_replicates": False,
            "frames_treated_as_independent_samples": False,
            "independent_lab_holdout_within_source": False,
            "uncertainty_holdout_passed": False,
        },
        "pre_sweep_constraint": {
            "purpose": (
                "reduce and rank a future Aleph contractility sweep; this does not replace "
                "the required Aleph forward simulation"
            ),
            "mechanism_family": "actomyosin_contractility_and_load_transfer",
            "direction": "myosin_II_inhibition_reduces_traction_and_monolayer_tension",
            "required_Aleph_operators": [
                "geometry_matched_mean_traction_magnitude_Pa",
                "geometry_matched_maximum_principal_monolayer_tension_mN_per_m",
            ],
            "candidate_axes": [
                "aleph.vertical.sf_arc.MaterialCard.active_tension_pn",
                "aleph.vertical.nmii.MotorKinetics.crossbridge_stiffness",
                "aleph.vertical.nmii.MotorKinetics.attachment_rate",
                "aleph.vertical.focal_adhesion.ClutchCard.bond_stiffness_pn_per_um",
            ],
            "identifiability": (
                "TFM/MSM jointly constrain an effective contractility/load-transfer family; "
                "they do not separately identify motor, stress-fibre, cortex, or adhesion axes"
            ),
            "future_sweep_range_reduction_intended": True,
            "may_emit_numeric_parameter_range": False,
            "may_select_Aleph_parameter": False,
        },
        "observable_evidence_eligible": False,
        "independent_lab_holdout": False,
        "uncertainty_holdout_passed": False,
        "aleph_parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "limitations": [
            "one dated culture/batch with 12 technical fields of view",
            "PIV contains one longitudinal sequence",
            "no geometry-matched Aleph TFM/MSM observation operator has passed parity",
            "the t2 table role is unresolved and therefore excluded",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate(args.manifest, args.manifest_report)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
