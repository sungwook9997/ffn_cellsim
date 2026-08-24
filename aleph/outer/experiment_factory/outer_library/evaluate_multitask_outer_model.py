#!/usr/bin/env python3
"""Build the task-scoped Outer Library model contract.

This module deliberately does not turn heterogeneous evidence into one global
score.  Each task head receives only the authority earned by its own dataset/
lab split, uncertainty audit, and measurement operator.  The resulting JSON is
the runtime contract consumed by the Aleph sweep gate; it is not a physics
parameter predictor and cannot mutate Aleph.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


SCHEMA = "aleph.outer_library.multitask_model.v1"


def _require_none_authority(name: str, report: dict[str, Any]) -> None:
    if report.get("aleph_authority", "none") != "none":
        raise ValueError(f"{name} claims forbidden Aleph authority")
    if report.get("aleph_parameter_proposal_eligible", False):
        raise ValueError(f"{name} claims direct Aleph parameter proposal authority")


def _classification_metrics_pass(metrics: dict[str, Any]) -> bool:
    return all(float(metrics[name]) >= 0.70 for name in (
        "balanced_accuracy", "macro_f1", "auroc",
    ))


def _head(
    *,
    head_id: str,
    task: str,
    kind: str,
    modalities: list[str],
    split_contract: dict[str, Any],
    metrics: dict[str, Any],
    status: str,
    evidence_eligible: bool,
    uncertainty_eligible: bool,
    external_holdout: bool,
    pristine_confirmation: bool,
    limitations: list[str],
) -> dict[str, Any]:
    return {
        "head_id": head_id,
        "task": task,
        "kind": kind,
        "modalities": modalities,
        "split_contract": split_contract,
        "metrics": metrics,
        "status": status,
        "external_holdout": external_holdout,
        "pristine_confirmation": pristine_confirmation,
        "uncertainty_eligible": uncertainty_eligible,
        "observable_evidence_eligible": evidence_eligible,
        "parameter_proposal_eligible": False,
        "aleph_authority": "none",
        "limitations": limitations,
    }


def evaluate(
    calcium_report: dict[str, Any],
    piezo1_report: dict[str, Any],
    calcium_uncertainty_report: dict[str, Any],
    mechanism_report: dict[str, Any],
    migration_report: dict[str, Any],
    dryad_migration_report: dict[str, Any],
    fibrotic_report: dict[str, Any],
    mechano_osmotic_report: dict[str, Any],
    piv_report: dict[str, Any],
    cell_monolayer_report: dict[str, Any],
    temporal_cnn_report: dict[str, Any] | None = None,
    if_translocation_report: dict[str, Any] | None = None,
    hpa_if_report: dict[str, Any] | None = None,
    microglia_report: dict[str, Any] | None = None,
    cell_cycle_report: dict[str, Any] | None = None,
    mitochondrial_stress_report: dict[str, Any] | None = None,
    sciplex_design_report: dict[str, Any] | None = None,
    apoptosis_report: dict[str, Any] | None = None,
    sciplex_state_report: dict[str, Any] | None = None,
    senscout_report: dict[str, Any] | None = None,
    figshare_cell_death_report: dict[str, Any] | None = None,
    external_senescence_report: dict[str, Any] | None = None,
    senescence_confirmation_report: dict[str, Any] | None = None,
    tfm_tether_report: dict[str, Any] | None = None,
    piv_contractility_report: dict[str, Any] | None = None,
    nematic_tfm_report: dict[str, Any] | None = None,
    opto_force_report: dict[str, Any] | None = None,
    tfm_fret_report: dict[str, Any] | None = None,
    exif_emt_report: dict[str, Any] | None = None,
    mtrack_emt_report: dict[str, Any] | None = None,
    lincs_mcf10a_report: dict[str, Any] | None = None,
    karacosta_emt_report: dict[str, Any] | None = None,
    geo_emt_rna_report: dict[str, Any] | None = None,
    gse325309_emt_report: dict[str, Any] | None = None,
    pbmc_cross_lab_report: dict[str, Any] | None = None,
    pbmc_ifn_direction_report: dict[str, Any] | None = None,
    pbmc_final_confirmation_report: dict[str, Any] | None = None,
    idr0072_if_report: dict[str, Any] | None = None,
    lightmycells_report: dict[str, Any] | None = None,
    idr0168_decode_refusal: dict[str, Any] | None = None,
    mechanics_protocol_report: dict[str, Any] | None = None,
) -> dict[str, Any]:
    reports = {
        "calcium": calcium_report,
        "piezo1": piezo1_report,
        "calcium_uncertainty": calcium_uncertainty_report,
        "mechanism": mechanism_report,
        "migration": migration_report,
        "dryad_migration": dryad_migration_report,
        "fibrotic": fibrotic_report,
        "mechano_osmotic": mechano_osmotic_report,
        "piv": piv_report,
        "cell_monolayer": cell_monolayer_report,
    }
    if temporal_cnn_report is not None:
        reports["temporal_cnn"] = temporal_cnn_report
    if if_translocation_report is not None:
        reports["if_translocation"] = if_translocation_report
    if hpa_if_report is not None:
        reports["hpa_if"] = hpa_if_report
    if microglia_report is not None:
        reports["microglia"] = microglia_report
    if cell_cycle_report is not None:
        reports["cell_cycle"] = cell_cycle_report
    if mitochondrial_stress_report is not None:
        reports["mitochondrial_stress"] = mitochondrial_stress_report
    if sciplex_design_report is not None:
        reports["sciplex_design"] = sciplex_design_report
    if apoptosis_report is not None:
        reports["apoptosis"] = apoptosis_report
    if sciplex_state_report is not None:
        reports["sciplex_state"] = sciplex_state_report
    if senscout_report is not None:
        reports["senscout"] = senscout_report
    if figshare_cell_death_report is not None:
        reports["figshare_cell_death"] = figshare_cell_death_report
    if external_senescence_report is not None:
        reports["external_senescence"] = external_senescence_report
    if senescence_confirmation_report is not None:
        reports["senescence_confirmation"] = senescence_confirmation_report
    if tfm_tether_report is not None:
        reports["tfm_tether"] = tfm_tether_report
    if piv_contractility_report is not None:
        reports["piv_contractility"] = piv_contractility_report
    if nematic_tfm_report is not None:
        reports["nematic_tfm"] = nematic_tfm_report
    if opto_force_report is not None:
        reports["opto_force"] = opto_force_report
    if tfm_fret_report is not None:
        reports["tfm_fret"] = tfm_fret_report
    if exif_emt_report is not None:
        reports["exif_emt"] = exif_emt_report
    if mtrack_emt_report is not None:
        reports["mtrack_emt"] = mtrack_emt_report
    if lincs_mcf10a_report is not None:
        reports["lincs_mcf10a"] = lincs_mcf10a_report
    if karacosta_emt_report is not None:
        reports["karacosta_emt"] = karacosta_emt_report
    if geo_emt_rna_report is not None:
        reports["geo_emt_rna"] = geo_emt_rna_report
    if gse325309_emt_report is not None:
        reports["gse325309_emt"] = gse325309_emt_report
    if pbmc_cross_lab_report is not None:
        reports["pbmc_cross_lab"] = pbmc_cross_lab_report
    if pbmc_final_confirmation_report is not None:
        reports["pbmc_final_confirmation"] = pbmc_final_confirmation_report
    if idr0072_if_report is not None:
        reports["idr0072_if"] = idr0072_if_report
    if lightmycells_report is not None:
        reports["lightmycells"] = lightmycells_report
    if idr0168_decode_refusal is not None:
        reports["idr0168_decode_refusal"] = idr0168_decode_refusal
    if mechanics_protocol_report is not None:
        reports["mechanics_protocol"] = mechanics_protocol_report
    for name, report in reports.items():
        _require_none_authority(name, report)

    if calcium_report.get("sample_overlap_count") != 0:
        raise ValueError("calcium internal split contains sample overlap")
    if piezo1_report.get("sample_overlap_count") != 0:
        raise ValueError("Piezo1 cross-dataset split contains sample overlap")
    calcium_selection = piezo1_report["aligned_active_nanoswitch_task"][
        "selection_protocol"
    ]
    if "never used for selection" not in calcium_selection:
        raise ValueError("Piezo1 external labels are not sealed from model selection")
    uncertainty_selection = calcium_uncertainty_report["selection_protocol"]
    if "no EA.hy926 label tunes anything" not in uncertainty_selection:
        raise ValueError("calcium external labels are not sealed from uncertainty fitting")
    if temporal_cnn_report is not None:
        if "EA.hy926 labels never select" not in temporal_cnn_report[
            "selection_protocol"
        ]:
            raise ValueError("temporal CNN used external labels for model selection")
        if temporal_cnn_report["split"]["sample_overlap_count"] != 0:
            raise ValueError("temporal CNN train/external split contains sample overlap")

    migration_protocol = migration_report["protocol"]
    migration_labs = {
        migration_protocol["train_lab"],
        migration_protocol["validation_and_conformal_calibration_lab"],
        migration_protocol["external_test_lab"],
    }
    if len(migration_labs) != 3:
        raise ValueError("migration train/calibration/test labs must be disjoint")
    if migration_protocol["test_lab_labels_used_for_numeric_fitting"]:
        raise ValueError("migration external labels were used for numeric fitting")

    fibrotic_labs = {
        fibrotic_report["discovery"]["lab_group"],
        fibrotic_report["external_confirmatory"]["lab_group"],
    }
    if len(fibrotic_labs) != 2:
        raise ValueError("fibrotic discovery and confirmation labs must be disjoint")

    calcium_external = piezo1_report["aligned_active_nanoswitch_task"][
        "external_metrics"
    ]
    calcium_classification_passed = _classification_metrics_pass(calcium_external)
    calcium_uncertainty_passed = bool(
        calcium_uncertainty_report["uncertainty_holdout_passed"]
    )
    migration_classification_passed = bool(
        migration_report["external_lab_holdout_passed"]
        and _classification_metrics_pass(migration_report["external_lab3_metrics"])
    )
    migration_uncertainty_passed = bool(
        migration_report["conformal_uncertainty"]["uncertainty_holdout_passed"]
    )
    fibrotic_mechanism_passed = bool(
        fibrotic_report["independent_lab_mechanism_holdout_passed"]
        and fibrotic_report["observable_constraint_eligible"]
    )
    calcium_mechanism_passed = bool(
        mechanism_report["mechanism_holdout_passed"]
        and mechanism_report["observable_constraint_eligible"]
    )

    heads = [
        _head(
            head_id="calcium_cell_state",
            task=piezo1_report["task"],
            kind="cell_level_binary_classifier",
            modalities=["single_cell_calcium_time_series"],
            split_contract={
                "train_dataset": piezo1_report["training_dataset"],
                "validation_scope": "disjoint HeLa cells",
                "calibration_scope": "disjoint HeLa cells",
                "external_test_dataset": piezo1_report["external_dataset"],
                "external_test_lab_holdout": piezo1_report["external_lab_holdout"],
                "external_labels_used_for_selection": False,
                "split_unit": piezo1_report["split_unit"],
            },
            metrics={
                "internal_test": calcium_report["test"],
                "external_test": calcium_external,
                "external_ood": calcium_uncertainty_report["dataset_shift"],
                "external_conformal": calcium_uncertainty_report["conformal"],
                "trainable_temporal_CNN_diagnostic": temporal_cnn_report,
            },
            status=(
                "validated_external" if calcium_classification_passed
                else "blocked_external_classification"
            ),
            evidence_eligible=False,
            uncertainty_eligible=calcium_uncertainty_passed,
            external_holdout=True,
            pristine_confirmation=True,
            limitations=list(piezo1_report["limitations"]) + [
                "external rank signal and uncertainty coverage do not rescue failed per-cell classification",
                "the trained temporal CNN also failed external macro-F1 and is not an accepted head",
            ],
        ),
        _head(
            head_id="piezo1_calcium_mechanism",
            task=mechanism_report["task"],
            kind="matched_condition_contrast",
            modalities=["single_cell_calcium_time_series_derived_AUC"],
            split_contract={
                "discovery_dataset": mechanism_report["discovery"]["dataset"],
                "external_test_dataset": mechanism_report[
                    "external_confirmatory"
                ]["dataset"],
                "external_test_lab_holdout": mechanism_report[
                    "external_confirmatory"
                ]["lab_holdout"],
                "inference_scope": "equal_weight_predeclared_active_strata",
            },
            metrics={
                "discovery": mechanism_report["discovery"],
                "external_test": mechanism_report["external_confirmatory"],
            },
            status=(
                "validated_external_evidence_only"
                if calcium_mechanism_passed else "blocked_external_evidence"
            ),
            evidence_eligible=calcium_mechanism_passed,
            uncertainty_eligible=False,
            external_holdout=True,
            pristine_confirmation=True,
            limitations=mechanism_report["limitations"],
        ),
        _head(
            head_id="paired_migration_state",
            task=migration_report["task"],
            kind="matched_control_multivariate_classifier",
            modalities=["live_cell_morphology_time_series_derived"],
            split_contract={
                "train_lab": migration_protocol["train_lab"],
                "validation_and_calibration_lab": migration_protocol[
                    "validation_and_conformal_calibration_lab"
                ],
                "external_test_lab": migration_protocol["external_test_lab"],
                "test_labels_used_for_fitting": False,
                "inference_unit": migration_protocol["inference_unit"],
                "matched_control_required": True,
            },
            metrics={
                "validation": migration_report["validation_metrics"],
                "external_test": migration_report["external_lab3_metrics"],
                "external_conformal": migration_report["conformal_uncertainty"],
                "external_experiment_effect_ci95": migration_report[
                    "external_lab3_effect_bootstrap_95ci"
                ],
                "dataset_shift": migration_report["dataset_shift"],
            },
            status=(
                "validated_retrospective_external_evidence_only"
                if migration_classification_passed and migration_uncertainty_passed
                else "blocked_external_classification"
            ),
            evidence_eligible=(
                migration_classification_passed and migration_uncertainty_passed
            ),
            uncertainty_eligible=migration_uncertainty_passed,
            external_holdout=True,
            pristine_confirmation=bool(
                migration_report["pristine_sealed_confirmation"]
            ),
            limitations=migration_report["limitations"],
        ),
        _head(
            head_id="migration_mode_representation",
            task=dryad_migration_report["task"],
            kind="same_lab_grouped_representation_classifier",
            modalities=["single_cell_tracking_organizational_features"],
            split_contract={
                "split_group": dryad_migration_report["protocol"]["split_group"],
                "train_dataset": dryad_migration_report["protocol"][
                    "training_dataset"
                ],
                "test_dataset": dryad_migration_report["protocol"][
                    "retrospective_test_dataset"
                ],
                "independent_lab_holdout": False,
            },
            metrics={
                "validation": dryad_migration_report[
                    "validation_observation_metrics"
                ],
                "retrospective_tracked_cell": dryad_migration_report[
                    "retrospective_talin_tracked_cell_metrics"
                ],
                "conformal": dryad_migration_report["conformal_uncertainty"],
            },
            status="training_only_same_lab",
            evidence_eligible=False,
            uncertainty_eligible=False,
            external_holdout=False,
            pristine_confirmation=False,
            limitations=[
                "strong grouped retrospective transfer is confined to one lab",
                "author migration-mode labels do not define a universal cell state",
            ],
        ),
        _head(
            head_id="fibrotic_state_mechanism",
            task=fibrotic_report["task"],
            kind="cross_modality_condition_contrast",
            modalities=["immunofluorescence_derived", "bulk_RNA_seq_derived"],
            split_contract={
                "discovery_lab": fibrotic_report["discovery"]["lab_group"],
                "external_test_lab": fibrotic_report[
                    "external_confirmatory"
                ]["lab_group"],
                "external_dataset_holdout": fibrotic_report[
                    "independent_dataset_holdout"
                ],
                "external_lab_holdout": fibrotic_report[
                    "independent_lab_holdout"
                ],
                "inference_scope": "condition_level_marker_panel_contrast",
            },
            metrics={
                "discovery": fibrotic_report["discovery"],
                "external_test": fibrotic_report["external_confirmatory"],
            },
            status=(
                "validated_external_evidence_only"
                if fibrotic_mechanism_passed else "blocked_external_evidence"
            ),
            evidence_eligible=fibrotic_mechanism_passed,
            uncertainty_eligible=False,
            external_holdout=True,
            pristine_confirmation=True,
            limitations=fibrotic_report["limitations"],
        ),
        _head(
            head_id="mechano_osmotic_observable",
            task=mechano_osmotic_report["task"],
            kind="experiment_grouped_observable_contrast",
            modalities=[
                "FXm_RICM_cell_spreading_derived",
                "AFM_membrane_tether_force_derived",
            ],
            split_contract={
                "inference_unit": "author_independent_experiment",
                "independent_lab_holdout": False,
                "measurement_operator_required": "AFM_tether_operator_pending",
            },
            metrics={
                "initial_spreading_area_rate": mechano_osmotic_report[
                    "initial_spreading_area_rate"
                ],
                "initial_volume_flux": mechano_osmotic_report[
                    "initial_volume_flux"
                ],
                "AFM_membrane_tether_force": mechano_osmotic_report[
                    "AFM_membrane_tether_force"
                ],
            },
            status="evidence_only_operator_blocked",
            evidence_eligible=bool(
                mechano_osmotic_report["observable_constraint_eligible"]
            ),
            uncertainty_eligible=False,
            external_holdout=False,
            pristine_confirmation=False,
            limitations=mechano_osmotic_report["operator_assumptions"],
        ),
        _head(
            head_id="velocity_field_representation",
            task="PIV_and_collective_velocity_representation",
            kind="grouped_representation_encoder",
            modalities=["PIV_derived", "single_cell_tracking_velocity_field_derived"],
            split_contract={
                "piv_biological_replicates": piv_report[
                    "biological_replicate_counts"
                ],
                "collective_velocity_biological_replicates": cell_monolayer_report[
                    "raw_velocity_tensor"
                ]["biological_replicate_count"],
                "mandatory_collective_velocity_split_group": cell_monolayer_report[
                    "raw_velocity_tensor"
                ]["mandatory_split_group"],
                "external_holdout": False,
            },
            metrics={
                "piv": piv_report["metric_summaries"],
                "collective_velocity": cell_monolayer_report[
                    "field_descriptions"
                ],
            },
            status="training_only_no_biological_replication",
            evidence_eligible=False,
            uncertainty_eligible=False,
            external_holdout=False,
            pristine_confirmation=False,
            limitations=(
                piv_report["limitations"] + cell_monolayer_report["limitations"]
            ),
        ),
    ]

    if if_translocation_report is not None:
        dataset_holdout = if_translocation_report["dataset_holdout"]
        if dataset_holdout["independent_lab_holdout"]:
            raise ValueError(
                "shared-lineage BBBC plates cannot claim independent-lab holdout"
            )
        if not dataset_holdout["shared_provider_lineage"]:
            raise ValueError("BBBC shared provider lineage was erased")
        if if_translocation_report["may_select_aleph_parameter"]:
            raise ValueError("IF representation claims forbidden parameter authority")
        representation_passed = bool(
            dataset_holdout["translocation_representation_passed"]
        )
        uncertainty_passed = bool(
            if_translocation_report["uncertainty"]["uncertainty_holdout_passed"]
        )
        heads.append(
            _head(
                head_id="if_translocation_representation",
                task=if_translocation_report["task"],
                kind="paired_channel_IF_field_representation",
                modalities=["IF_FKHR_EGFP_DRAQ_field", "IF_NFkB_FITC_DAPI_field"],
                split_contract={
                    "train_dataset": dataset_holdout["train_dataset"],
                    "sealed_external_dataset": dataset_holdout[
                        "sealed_external_dataset"
                    ],
                    "dataset_holdout": dataset_holdout["dataset_holdout_present"],
                    "independent_lab_holdout": False,
                    "shared_provider_lineage": True,
                    "split_unit": if_translocation_report["observation_unit"],
                    "external_labels_used_for_selection": False,
                },
                metrics={
                    "contexts": if_translocation_report["contexts"],
                    "uncertainty": if_translocation_report["uncertainty"],
                    "dataset_ood": if_translocation_report["ood"],
                },
                status=(
                    "training_only_dataset_holdout_threshold_blocked"
                    if not representation_passed
                    else (
                        "training_only_dataset_holdout_uncertainty_blocked"
                        if not uncertainty_passed
                        else "training_only_independent_lab_holdout_blocked"
                    )
                ),
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=True,
                pristine_confirmation=True,
                limitations=if_translocation_report["limitations"],
            )
        )

    if hpa_if_report is not None:
        if hpa_if_report["independent_dataset_holdout"]:
            raise ValueError("HPA grouped split cannot claim an independent dataset")
        if hpa_if_report["independent_lab_holdout"]:
            raise ValueError("HPA grouped split cannot claim an independent lab")
        if hpa_if_report["split"]["gene_or_antibody_cross_split_overlap_count"]:
            raise ValueError("HPA IF split contains gene/antibody leakage")
        if hpa_if_report["split"]["test_labels_used_for_fitting_or_threshold_selection"]:
            raise ValueError("HPA test labels were used for fitting or selection")
        heads.append(
            _head(
                head_id="hpa_if_localization_representation",
                task=hpa_if_report["task"],
                kind="provider_embedding_multilabel_classifier",
                modalities=["immunofluorescence_subcellular_localization_image_embedding"],
                split_contract={
                    "dataset": hpa_if_report["dataset_id"],
                    "split_unit": hpa_if_report["split"]["unit"],
                    "gene_or_antibody_overlap_count": 0,
                    "sealed_grouped_test": True,
                    "independent_dataset_holdout": False,
                    "independent_lab_holdout": False,
                    "test_labels_used_for_selection": False,
                },
                metrics={
                    "validation": hpa_if_report["validation_metrics"],
                    "sealed_grouped_test": hpa_if_report[
                        "sealed_grouped_test_metrics"
                    ],
                    "conformal_diagnostic": hpa_if_report[
                        "conformal_diagnostic"
                    ],
                },
                status="training_only_grouped_same_provider",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=hpa_if_report["limitations"],
            )
        )

    if exif_emt_report is not None:
        if not exif_emt_report["representation_training_eligible"]:
            raise ValueError("ExIF EMT raw-image path is not training eligible")
        if (
            exif_emt_report["observable_evidence_eligible"]
            or exif_emt_report["independent_dataset_holdout"]
            or exif_emt_report["independent_lab_holdout"]
            or exif_emt_report["uncertainty_holdout_passed"]
            or exif_emt_report["OOD_refusal_validated"]
        ):
            raise ValueError("single-date ExIF EMT pilot claims external authority")
        if exif_emt_report["pre_sweep_role"]["may_emit_sweep_axis"]:
            raise ValueError("ExIF EMT pilot claims a sweep axis")
        heads.append(
            _head(
                head_id="raw_IF_EMT_condition_representation",
                task=exif_emt_report["task"],
                kind="field_grouped_raw_multichannel_IF_condition_representation",
                modalities=["raw_five_channel_4plex_IF_field"],
                split_contract=exif_emt_report["design_audit"],
                metrics={
                    "paired_well_column_transfer": exif_emt_report[
                        "paired_well_column_transfer_general_channels_only"
                    ],
                    "mean_paired_transfer_accuracy": exif_emt_report[
                        "mean_paired_transfer_accuracy"
                    ],
                    "mean_paired_transfer_macro_f1": exif_emt_report[
                        "mean_paired_transfer_macro_f1"
                    ],
                    "TGF_beta1_marker_direction_check": exif_emt_report[
                        "TGF_beta1_marker_direction_check"
                    ],
                },
                status="training_only_single_plate_treatment_labels_no_external_holdout",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=exif_emt_report["limitations"],
            )
        )

    if mtrack_emt_report is not None:
        if not mtrack_emt_report["representation_training_eligible"]:
            raise ValueError("M-TRACK EMT temporal path is not training eligible")
        if (
            mtrack_emt_report["observable_evidence_eligible"]
            or mtrack_emt_report["independent_dataset_holdout"]
            or mtrack_emt_report["independent_lab_holdout"]
            or mtrack_emt_report["uncertainty_holdout_passed"]
            or mtrack_emt_report["OOD_refusal_validated"]
        ):
            raise ValueError("single-trajectory M-TRACK pilot claims external authority")
        if mtrack_emt_report["pre_sweep_role"]["may_emit_sweep_axis"]:
            raise ValueError("M-TRACK EMT pilot claims a sweep axis")
        design = mtrack_emt_report["design_audit"]
        if (
            design["trajectory_split_groups"] != 1
            or design["frames_permitted_as_train_test_units"]
            or design["biological_replicates_released"]
            or design["untreated_control_released"]
        ):
            raise ValueError("M-TRACK one-trajectory grouping contract changed")
        heads.append(
            _head(
                head_id="live_VIM_RFP_EMT_temporal_representation",
                task=mtrack_emt_report["task"],
                kind="single_trajectory_raw_fluorescence_temporal_diagnostic",
                modalities=["live_VIM_RFP_single_trajectory_with_segmentation"],
                split_contract=design,
                metrics={
                    "temporal_descriptive_trends": mtrack_emt_report[
                        "temporal_descriptive_trends"
                    ],
                    "cross_lab_vimentin_direction_context": mtrack_emt_report[
                        "cross_lab_vimentin_direction_context"
                    ],
                },
                status="training_only_single_trajectory_no_control_direction_conflict",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=mtrack_emt_report["limitations"],
            )
        )

    if lincs_mcf10a_report is not None:
        role = lincs_mcf10a_report["pre_sweep_role"]
        design = lincs_mcf10a_report["design_audit"]
        if not lincs_mcf10a_report["collection_holdout_passed"]:
            raise ValueError("LINCS MCF10A collection holdout failed")
        if (
            lincs_mcf10a_report["independent_lab_holdout"]
            or lincs_mcf10a_report["OOD_refusal_validated"]
            or role["may_emit_sweep_axis"]
            or lincs_mcf10a_report["may_select_aleph_parameter"]
        ):
            raise ValueError("LINCS MCF10A collection evidence overclaims authority")
        if not design["technical_rows_aggregated_before_inference"]:
            raise ValueError(
                "LINCS MCF10A technical wells entered inference independently"
            )
        heads.append(
            _head(
                head_id="MCF10A_IF_ligand_state_collection_holdout",
                task="TGFB_plus_EGF_vs_EGF_MCF10A_IF_phenotype_state",
                kind="replicate_aggregated_IF_phenotype_nearest_centroid_classifier",
                modalities=["fixed_IF_author_well_level_population_summary"],
                split_contract=design,
                metrics={
                    "collection_holdout_classifier": lincs_mcf10a_report[
                        "collection_holdout_classifier"
                    ],
                    "paired_direction_holdout": lincs_mcf10a_report[
                        "paired_direction_holdout"
                    ],
                    "uncertainty_status": lincs_mcf10a_report["uncertainty_status"],
                },
                status=(
                    "training_only_same_lab_collection_holdout_"
                    "OOD_and_common_lab_blocked"
                ),
                evidence_eligible=lincs_mcf10a_report["observable_evidence_eligible"],
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=lincs_mcf10a_report["limitations"] + [
                    "this state constraint narrows a future Aleph sweep and does not replace it",
                ],
            )
        )

    if karacosta_emt_report is not None:
        design = karacosta_emt_report["design_audit"]
        direction = karacosta_emt_report["cross_lab_common_marker_direction"]
        role = karacosta_emt_report["pre_sweep_role"]
        if not (
            karacosta_emt_report["representation_training_eligible"]
            and karacosta_emt_report["directional_bridge_eligible"]
            and direction["independent_lab_common_marker_direction_reproduced"]
            and direction["dataset_and_lab_disjoint"]
        ):
            raise ValueError("Karacosta independent-lab EMT direction bridge failed")
        if (
            karacosta_emt_report["observable_evidence_eligible"]
            or karacosta_emt_report["independent_lab_inferential_holdout"]
            or karacosta_emt_report["uncertainty_holdout_passed"]
            or karacosta_emt_report["OOD_refusal_validated"]
            or role["may_emit_sweep_axis"]
            or karacosta_emt_report["may_select_aleph_parameter"]
        ):
            raise ValueError("Karacosta one-experiment source overclaims authority")
        if (
            design["mandatory_split_groups"] != 1
            or not design["single_cells_are_not_independent_biological_replicates"]
            or not design["cell_state_classifier_on_these_features_is_label_circular_and_forbidden"]
        ):
            raise ValueError("Karacosta pseudoreplication/circular-label guard changed")
        heads.append(
            _head(
                head_id="cross_modal_EMT_common_protein_direction_bridge",
                task=karacosta_emt_report["task"],
                kind="independent_lab_cross_modality_direction_filter",
                modalities=[
                    "single_cell_mass_cytometry_CyTOF_nontransformed",
                    "raw_five_channel_4plex_IF_field",
                ],
                split_contract=design,
                metrics={
                    "progressive_EMT": karacosta_emt_report[
                        "progressive_EMT_0d_to_10d_TGFB"
                    ],
                    "withdrawal_reversal": karacosta_emt_report[
                        "withdrawal_10d_TGFB_to_10dW"
                    ],
                    "cross_lab_common_marker_direction": direction,
                },
                status=(
                    "training_only_independent_lab_common_marker_direction_"
                    "biological_replication_and_OOD_blocked"
                ),
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=karacosta_emt_report["limitations"] + [
                    "this direction filter compresses a future Aleph sweep and cannot replace it",
                ],
            )
        )

    if geo_emt_rna_report is not None:
        design = geo_emt_rna_report["design_audit"]
        direction = geo_emt_rna_report["replicated_direction"]
        cross_assay = geo_emt_rna_report["cross_assay_bridge"]
        ood = geo_emt_rna_report["semantic_OOD_challenge"]
        role = geo_emt_rna_report["pre_sweep_role"]
        if not (
            geo_emt_rna_report["representation_training_eligible"]
            and geo_emt_rna_report["directional_pre_sweep_filter_eligible"]
            and geo_emt_rna_report["independent_study_direction_holdout"]
            and direction["all_ten_study_gene_directions_expected"]
            and direction["all_ten_between_group_ranges_completely_separated"]
            and direction["retrospective_direction_replication_passed"]
            and not direction["study_selection_was_outcome_blind"]
            and not direction["confirmatory_direction_threshold_0_05_passed"]
            and not geo_emt_rna_report["study_level_direction_inference_passed"]
            and cross_assay[
                "RNA_direction_agrees_with_independent_IF_and_CyTOF_protein_direction"
            ]
        ):
            raise ValueError("replicated GEO EMT direction constraint failed")
        if (
            design["studies"] != 5
            or design["independent_lab_groups"] != 5
            or not design["biological_replicates_not_expression_features"]
            or geo_emt_rna_report["numeric_uncertainty_promotion_passed"]
            or ood["OOD_refusal_model_validated"]
            or geo_emt_rna_report["OOD_refusal_validated"]
            or role["numeric_range"] is not None
            or role["may_emit_sweep_axis"]
            or geo_emt_rna_report["may_select_aleph_parameter"]
        ):
            raise ValueError("replicated GEO EMT source overclaims sweep authority")
        heads.append(
            _head(
                head_id="replicated_multi_study_EMT_RNA_direction_constraint",
                task=geo_emt_rna_report["task"],
                kind="replicated_cross_study_RNA_direction_filter",
                modalities=["bulk_RNA_expression_derived"],
                split_contract=design,
                metrics={
                    "endpoint_contrasts": geo_emt_rna_report["endpoint_contrasts"],
                    "replicated_direction": direction,
                    "cross_assay_bridge": cross_assay,
                    "semantic_OOD_challenge": ood,
                },
                status=(
                    "training_only_retrospective_replicated_direction_"
                    "preregistered_holdout_numeric_uncertainty_OOD_and_operator_blocked"
                ),
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=geo_emt_rna_report["limitations"],
            )
        )

    if gse325309_emt_report is not None:
        primary = gse325309_emt_report["frozen_primary_test"]
        design = gse325309_emt_report["design_audit"]
        role = gse325309_emt_report["pre_sweep_role"]
        if not (
            gse325309_emt_report["pristine_external_direction_holdout_passed"]
            and primary["passed"] and primary["exact_one_sided_p"] == 0.05
            and primary["directions"] == {"CDH1": "decrease", "VIM": "increase"}
            and design["biological_samples"] == 6
            and design["control_replicates"] == 3
            and design["TGFB1_replicates"] == 3
            and design["sample_exclusions"] == 0
            and design["protocol_frozen_before_expression_access"]
        ):
            raise ValueError("GSE325309 preregistered external EMT holdout failed")
        if (
            gse325309_emt_report["numeric_uncertainty_promotion_passed"]
            or gse325309_emt_report["OOD_refusal_validated"]
            or role["numeric_range"] is not None or role["may_emit_sweep_axis"]
            or gse325309_emt_report["may_select_aleph_parameter"]
        ):
            raise ValueError("GSE325309 direction holdout overclaims Aleph authority")
        heads.append(_head(
            head_id="preregistered_GSE325309_EMT_direction_holdout",
            task=gse325309_emt_report["task"],
            kind="preregistered_external_RNA_direction_holdout",
            modalities=["bulk_RNA_seq_author_processed"],
            split_contract=design,
            metrics={"frozen_primary_test": primary,
                     "secondary_diagnostics": gse325309_emt_report["secondary_diagnostics"]},
            status="validated_external_direction_only_numeric_OOD_and_Aleph_operator_blocked",
            evidence_eligible=False, uncertainty_eligible=False,
            external_holdout=True, pristine_confirmation=True,
            limitations=[
                "direction-only RNA holdout does not calibrate a protein or Aleph numeric scale",
                "one frozen context cannot validate semantic OOD refusal",
            ],
        ))

    if microglia_report is not None:
        if microglia_report["independent_biological_replicate_holdout"]:
            raise ValueError("BBBC054 has only one annotated biological replicate")
        if microglia_report["independent_lab_holdout"]:
            raise ValueError("BBBC054 cannot claim an independent lab holdout")
        if microglia_report["split"]["field_overlap_count"]:
            raise ValueError("BBBC054 contains field leakage")
        heads.append(
            _head(
                head_id="microglia_immune_morphophenotype_representation",
                task=microglia_report["task"],
                kind="brightfield_patch_multiclass_classifier",
                modalities=["brightfield_single_cell_patch_derived"],
                split_contract={
                    "split_unit": microglia_report["split"]["unit"],
                    "field_overlap_count": 0,
                    "annotated_biological_replicates": 1,
                    "independent_lab_holdout": False,
                    "pristine_confirmation": microglia_report["split"][
                        "pristine_confirmation"
                    ],
                },
                metrics={
                    "validation": microglia_report["validation_metrics"],
                    "retrospective_field_test": microglia_report[
                        "sealed_field_test_metrics"
                    ],
                    "conformal_diagnostic": microglia_report[
                        "split_conformal_diagnostic"
                    ],
                    "descriptive_LPS_timecourse": microglia_report[
                        "descriptive_LPS_timecourse"
                    ],
                },
                status="training_only_single_replicate_classifier_rejected",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=microglia_report["limitations"],
            )
        )

    if cell_cycle_report is not None:
        if cell_cycle_report["independent_lab_holdout"]:
            raise ValueError("BBBC048 same-provider split cannot claim an independent lab")
        if cell_cycle_report["split_sample_id_overlap"]:
            raise ValueError("BBBC048 image head contains cell split leakage")
        heads.append(
            _head(
                head_id="cell_cycle_image_representation",
                task="seven_phase_Jurkat_cell_cycle_classification",
                kind="raw_three_channel_image_MLP_classifier",
                modalities=["imaging_flow_cytometry_multichannel"],
                split_contract={
                    "split_unit": "cell_with_all_three_channels",
                    "train_selection_calibration_test_disjoint": True,
                    "independent_dataset_holdout": False,
                    "independent_lab_holdout": False,
                },
                metrics={
                    "same_provider_test": cell_cycle_report["test"],
                    "same_provider_conformal": cell_cycle_report["conformal_90"],
                    "best_selection_epoch": cell_cycle_report["best_selection_epoch"],
                },
                status="training_only_same_provider_rare_phase_limited",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=cell_cycle_report["notes"] + [
                    "marginal conformal coverage does not establish rare-phase conditional coverage",
                ],
            )
        )

    if mitochondrial_stress_report is not None:
        if mitochondrial_stress_report["biological_replicates_reported"]:
            raise ValueError("BBBC053 unexpectedly claims biological replicate identifiers")
        if mitochondrial_stress_report["independent_lab_holdout"]:
            raise ValueError("BBBC053 cannot claim an independent lab holdout")
        heads.append(
            _head(
                head_id="mitochondrial_stress_IF_representation",
                task="FCCP_mitochondrial_morphology_representation",
                kind="raw_TOM20_IF_representation",
                modalities=["TOM20_super_resolution_IF_image"],
                split_contract={
                    "source_images": mitochondrial_stress_report["image_count"],
                    "biological_replicates_reported": False,
                    "independent_lab_holdout": False,
                },
                metrics={
                    "condition_counts": mitochondrial_stress_report["condition_counts"],
                    "compact_image_shape": mitochondrial_stress_report["compact_image_shape"],
                    "morphology_feature_count": mitochondrial_stress_report["morphology_feature_count"],
                },
                status="training_only_no_biological_replication",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=[
                    "FCCP versus DMSO directory labels do not identify a universal stress latent",
                    "no independent biological replicate or external provider test is available",
                ],
            )
        )

    if sciplex_design_report is not None:
        if sciplex_design_report["unassigned_context_cells_may_enter_supervised_training"]:
            raise ValueError("unassigned sci-Plex cells cannot enter supervised training")
        if sciplex_state_report is not None:
            if not sciplex_design_report["training_ready"]:
                raise ValueError("sci-Plex state model requires a verified matrix receipt")
            if sciplex_state_report["split"]["cross_split_treatment_overlap_count"]:
                raise ValueError("sci-Plex state model leaks treatments across splits")
            if sciplex_state_report["split"][
                "test_treatments_used_for_fitting_selection_or_calibration"
            ]:
                raise ValueError("sci-Plex sealed treatments entered model selection")
            if (
                sciplex_state_report["independent_dataset_holdout"]
                or sciplex_state_report["independent_lab_holdout"]
                or sciplex_state_report["production_eligible"]
            ):
                raise ValueError("same-provider sci-Plex model claims external authority")
            state_status = "training_only_unseen_compound_holdout_same_provider"
            state_metrics = {
                "model": sciplex_state_report["model"],
                "sealed_unseen_treatment_test": sciplex_state_report[
                    "sealed_unseen_treatment_test"
                ],
                "conformal": sciplex_state_report["conformal_90"],
                "vehicle_ood": sciplex_state_report["vehicle_ood_diagnostic"],
                "replicate_consistency": sciplex_state_report["replicate_consistency"],
            }
            state_limitations = sciplex_state_report["limitations"]
        else:
            state_status = (
                "training_only_model_pending"
                if sciplex_design_report["training_ready"]
                else "training_only_matrix_pending"
            )
            state_metrics = {"design_audit": sciplex_design_report}
            state_limitations = [
                (
                    "context tensor and treatment-disjoint model training are pending"
                    if sciplex_design_report["training_ready"] else
                    "expression matrix verification and context tensor construction are pending"
                ),
                "drug and pathway metadata are perturbation context, not direct cell-state truth",
                "both replicates share one provider lab",
            ]
        heads.append(
            _head(
                head_id="sciplex_perturbation_state_representation",
                task="drug_dose_time_cell_type_transcriptomic_context_representation",
                kind="single_cell_RNA_context_encoder",
                modalities=["single_cell_RNA_seq"],
                split_contract={
                    "mandatory_context_group_count": sciplex_design_report[
                        "mandatory_context_group_count"
                    ],
                    "assigned_context_cells": sciplex_design_report[
                        "assigned_context_cell_count"
                    ],
                    "excluded_unassigned_context_cells": sciplex_design_report[
                        "unassigned_context_cell_count"
                    ],
                    "replicate_holdout_is_same_lab": True,
                    "independent_lab_holdout": False,
                },
                metrics=state_metrics,
                status=state_status,
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=state_limitations,
            )
        )

    if apoptosis_report is not None:
        if apoptosis_report["independent_lab_holdout"]:
            raise ValueError("S-BIAD2515 same-provider wells cannot claim an independent lab")
        if apoptosis_report["split"]["well_overlap_count"]:
            raise ValueError("S-BIAD2515 apoptosis head contains well leakage")
        if apoptosis_report["split"]["test_used_for_fitting_selection_or_calibration"]:
            raise ValueError("S-BIAD2515 sealed test wells entered fitting")
        heads.append(
            _head(
                head_id="apoptosis_AnnexinV_endpoint_representation",
                task=apoptosis_report["task"],
                kind="author_processed_live_morphology_ridge_regressor",
                modalities=[
                    "live_multichannel_morphology_to_AnnexinV_endpoint_derived"
                ],
                split_contract={
                    "split_unit": apoptosis_report["split"]["unit"],
                    "dose_balanced_fit_calibration_test": True,
                    "sealed_test_wells": len(
                        apoptosis_report["split"]["sealed_test_wells"]
                    ),
                    "independent_dataset_holdout": False,
                    "independent_lab_holdout": False,
                },
                metrics={
                    "sealed_same_provider_test": apoptosis_report[
                        "sealed_same_provider_test"
                    ],
                    "baselines": apoptosis_report["baselines_on_same_sealed_test"],
                    "conformal": apoptosis_report["conformal_90"],
                },
                status="training_only_same_provider_well_holdout",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=apoptosis_report["limitations"],
            )
        )

    if tfm_tether_report is not None:
        promotion = tfm_tether_report["promotion"]
        if promotion["observable_evidence_eligible"]:
            raise ValueError("operator-blocked TFM/tether data claim evidence authority")
        if promotion["may_emit_numeric_parameter_range"]:
            raise ValueError("TFM/tether head claims a forbidden numeric sweep range")
        if not promotion["future_sweep_range_reduction_intended"]:
            raise ValueError("TFM/tether head lost its pre-sweep purpose")
        heads.append(
            _head(
                head_id="TFM_tether_pre_sweep_observation_constraint",
                task=tfm_tether_report["task"],
                kind="descriptive_measurement_and_operator_audit",
                modalities=[
                    "traction_force_microscopy_derived",
                    "optical_tweezers_membrane_tether_derived",
                ],
                split_contract=tfm_tether_report["design_audit"],
                metrics=tfm_tether_report["descriptive_results"],
                status="training_only_operator_and_replication_blocked",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=tfm_tether_report["limitations"] + [
                    "this head exists to reduce a future Aleph sweep, not to replace it",
                ],
            )
        )

    if piv_contractility_report is not None:
        constraint = piv_contractility_report["pre_sweep_constraint"]
        if piv_contractility_report["observable_evidence_eligible"]:
            raise ValueError("small-n PIV/IF head claims external evidence authority")
        if constraint["may_emit_numeric_parameter_range"]:
            raise ValueError("PIV/IF head claims a forbidden numeric sweep range")
        if not constraint["future_sweep_range_reduction_intended"]:
            raise ValueError("PIV/IF head lost its pre-sweep purpose")
        if not piv_contractility_report["replicated_contractility_direction"]:
            raise ValueError("PIV/IF author-replicated intervention direction failed")
        heads.append(
            _head(
                head_id="C2C12_PIV_IF_supracellular_contractility_pre_sweep_constraint",
                task=piv_contractility_report["task"],
                kind="replicated_directional_PIV_IF_observation_constraint",
                modalities=[
                    "particle_image_velocimetry_and_live_actin_derived",
                    "immunofluorescence_fibronectin_orientation_derived",
                ],
                split_contract=piv_contractility_report["design_audit"],
                metrics={
                    "late_post_blebbistatin": piv_contractility_report[
                        "late_post_blebbistatin_contrasts_t_ge_42h"
                    ],
                    "fibronectin_orientation": piv_contractility_report[
                        "fibronectin_orientation_nematic_order"
                    ],
                    "uncertainty": piv_contractility_report["uncertainty"],
                },
                status="training_only_directional_constraint_operator_and_external_holdout_blocked",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=piv_contractility_report["limitations"] + [
                    "this directional constraint narrows a future Aleph sweep and does not replace it",
                ],
            )
        )

    if nematic_tfm_report is not None or opto_force_report is not None:
        selected = [
            report for report in (nematic_tfm_report, opto_force_report)
            if report is not None
        ]
        for report in selected:
            constraint = report["pre_sweep_constraint"]
            if report["observable_evidence_eligible"]:
                raise ValueError("technical-only force source claims evidence authority")
            if report["independent_lab_holdout"] or report["uncertainty_holdout_passed"]:
                raise ValueError("unresolved force-source hierarchy claims a holdout")
            if constraint["may_emit_numeric_parameter_range"]:
                raise ValueError("force pre-sweep source claims a numeric range")
            if not constraint["future_sweep_range_reduction_intended"]:
                raise ValueError("force source lost its Aleph pre-sweep purpose")
        if nematic_tfm_report is not None and not nematic_tfm_report[
            "published_direction_reproduced_from_source_tables"
        ]:
            raise ValueError("NIH3T3 TFM/MSM published direction failed reproduction")
        if opto_force_report is not None and not opto_force_report[
            "overall_direction"
        ]["usually_positive_in_both_axes"]:
            raise ValueError("opto-MDCK cortical direction failed reproduction")
        labs = {report["lab_group"] for report in selected}
        if piv_contractility_report is not None:
            labs.add("anseth-white-skillin-lab-cu-boulder")
        heads.append(
            _head(
                head_id="cross_cell_type_contractility_force_pre_sweep_constraint",
                task="cross_cell_type_actomyosin_contractility_direction",
                kind="cross_lab_directional_force_constraint",
                modalities=[
                    "traction_force_and_monolayer_stress_microscopy_derived",
                    "contour_model_from_TFM_MSM_and_actin_tracking",
                    "particle_image_velocimetry_and_immunofluorescence_derived",
                ],
                split_contract={
                    "lab_groups": sorted(labs),
                    "lab_group_count": len(labs),
                    "aligned_same_observable_holdout": False,
                    "technical_units_promoted_to_biological_replicates": False,
                },
                metrics={
                    "NIH3T3_myosin_inhibition": (
                        nematic_tfm_report[
                            "figure_S9_paired_technical_FOV_contrasts"
                        ] if nematic_tfm_report is not None else None
                    ),
                    "opto_MDCK_RhoA_activation": (
                        opto_force_report["overall_direction"]
                        if opto_force_report is not None else None
                    ),
                    "C2C12_long_term_myosin_inhibition": (
                        piv_contractility_report["replicated_contractility_direction"]
                        if piv_contractility_report is not None else None
                    ),
                },
                status="training_only_cross_lab_direction_operator_and_replication_blocked",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=[
                    limitation for report in selected
                    for limitation in report["limitations"]
                ] + [
                    "different cell types and observables support a mechanism-family direction, not a numeric cross-lab holdout",
                    "this constraint ranks a future Aleph sweep and does not replace forward simulation",
                ],
            )
        )

    if tfm_fret_report is not None:
        constraint = tfm_fret_report["pre_sweep_constraint"]
        if tfm_fret_report["observable_evidence_eligible"]:
            raise ValueError("replicate-unmapped TFM--FRET source claims evidence authority")
        if tfm_fret_report["independent_lab_holdout"]:
            raise ValueError("direction-only TFM--FRET source claims inferential holdout")
        if tfm_fret_report["uncertainty_holdout_passed"]:
            raise ValueError("TFM--FRET source claims unmeasured uncertainty coverage")
        if constraint["may_emit_numeric_parameter_range"]:
            raise ValueError("TFM--FRET source claims a numeric sweep range")
        if not tfm_fret_report["cross_lab_same_observable_direction"][
            "independent_lab_direction_reproduced"
        ]:
            raise ValueError("independent-lab stiffness--TFM direction was not reproduced")
        if not tfm_fret_report["author_correlation_categories_reproduced"]:
            raise ValueError("paired TFM--FRET correlation categories were not reproduced")
        heads.append(
            _head(
                head_id="cross_lab_TFM_FRET_load_transfer_pre_sweep_constraint",
                task=tfm_fret_report["task"],
                kind="cross_lab_monotonic_TFM_and_multiscale_load_transfer_constraint",
                modalities=[
                    "traction_force_microscopy_derived",
                    "FLIM_FRET_vinculin_tension_sensor_derived",
                    "paired_FA_TFM_FLIM_FRET_derived",
                ],
                split_contract=tfm_fret_report["design_audit"],
                metrics={
                    "cross_lab_same_observable_direction": tfm_fret_report[
                        "cross_lab_same_observable_direction"
                    ],
                    "stiffness_directions": tfm_fret_report["stiffness_directions"],
                    "paired_FA_cell_correlations": tfm_fret_report[
                        "paired_FA_cell_correlations"
                    ],
                },
                status="training_only_same_observable_direction_operator_and_replication_blocked",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=tfm_fret_report["limitations"] + [
                    "this multiscale constraint reduces a future Aleph sweep and does not replace it",
                ],
            )
        )

    if senscout_report is not None:
        if senscout_report["split"]["cross_split_image_group_overlap_count"]:
            raise ValueError("SenSCOUT morphology head leaks image groups")
        if senscout_report["split"]["cells_treated_as_independent_biological_replicates"]:
            raise ValueError("SenSCOUT cannot treat cells as biological replicates")
        if (
            senscout_report["independent_dataset_holdout"]
            or senscout_report["independent_lab_holdout"]
            or senscout_report["production_eligible"]
        ):
            raise ValueError("same-provider SenSCOUT model claims external authority")
        heads.append(
            _head(
                head_id="senscout_morphology_senescence_representation",
                task=senscout_report["task"],
                kind="curated_morphology_multitask_encoder",
                modalities=["morphology_to_senescence_biomarker_state"],
                split_contract={
                    "split_unit": "image_group",
                    "fit_repeat": "Bio1",
                    "sealed_same_provider_repeat": "Bio2",
                    "frozen_shift_repeat": "Bio3",
                    "image_group_overlap_count": 0,
                    "independent_dataset_holdout": False,
                    "independent_lab_holdout": False,
                },
                metrics={
                    "direct_bgal": senscout_report["direct_bgal_head"],
                    "continuous_biomarkers": senscout_report[
                        "continuous_biomarker_heads"
                    ],
                    "repeat_shift": senscout_report["repeat_shift_diagnostic"],
                },
                status="training_only_same_provider_biorepeat_holdout",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=senscout_report["limitations"],
            )
        )

    if figshare_cell_death_report is not None:
        split = figshare_cell_death_report["split"]
        if split["cross_split_compound_overlap_count"]:
            raise ValueError("Figshare cell-death model leaks compound identities")
        if split["test_compounds_used_for_fit_selection_or_calibration"]:
            raise ValueError("Figshare sealed compounds entered model selection")
        if (
            figshare_cell_death_report["independent_dataset_holdout"]
            or figshare_cell_death_report["independent_lab_holdout"]
            or figshare_cell_death_report["production_eligible"]
            or figshare_cell_death_report["observable_evidence_eligible"]
        ):
            raise ValueError("same-provider Figshare model claims external authority")
        heads.append(
            _head(
                head_id="programmed_cell_death_mechanism_representation",
                task=figshare_cell_death_report["task"],
                kind="compound_grouped_Cell_Painting_classifier",
                modalities=["DeepProfiler_Cell_Painting_well_profile"],
                split_contract={
                    "split_unit": split["unit"],
                    "compound_counts": split["compound_counts"],
                    "cross_split_compound_overlap_count": 0,
                    "well_rows_are_not_independent_training_units": True,
                    "independent_dataset_holdout": False,
                    "independent_lab_holdout": False,
                },
                metrics={
                    "model": figshare_cell_death_report["model"],
                    "sealed_unseen_compound_test": figshare_cell_death_report[
                        "sealed_unseen_compound_test"
                    ],
                },
                status="training_only_unseen_compound_same_provider_rejected",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=figshare_cell_death_report["limitations"],
            )
        )

    if external_senescence_report is not None:
        split = external_senescence_report["split"]
        if split["dataset_overlap_count"] or split["lab_overlap_count"]:
            raise ValueError("external senescence split contains dataset/lab leakage")
        if split["external_test_labels_used_for_fit_selection_or_calibration"]:
            raise ValueError("external senescence labels entered model selection")
        if not (
            external_senescence_report["independent_dataset_holdout"]
            and external_senescence_report["independent_lab_holdout"]
        ):
            raise ValueError("external senescence head lost its dataset/lab holdout")
        if senescence_confirmation_report is None:
            raise ValueError("external senescence head requires post-freeze confirmation")
        if not senescence_confirmation_report[
            "checkpoint_frozen_before_dataset_acquisition"
        ]:
            raise ValueError("senescence confirmation is not genuinely post-freeze")
        if senescence_confirmation_report[
            "model_fit_selection_calibration_or_retraining_on_GSE301164"
        ]:
            raise ValueError("GSE301164 entered model fitting or selection")
        confirmation = senescence_confirmation_report["combined"]
        directional_passed = bool(
            confirmation["directional_effect_evidence_eligible"]
            and confirmation["positive_pair_count"] == confirmation["pair_count"]
            and float(confirmation["exact_one_sided_sign_test_p"]) < 0.05
        )
        sample_state_refused = bool(
            senescence_confirmation_report["ood"]["rejection_fraction"] > 0
            or not external_senescence_report["production_eligible"]
        )
        if not sample_state_refused:
            raise ValueError("senescence absolute-state head bypassed OOD refusal")
        heads.append(
            _head(
                head_id="senescence_transcriptomic_direction",
                task="fibroblast_senescence_direction_across_labs_and_triggers",
                kind="predeclared_paired_directional_contrast",
                modalities=["bulk_RNA_seq_frozen_35_marker_panel"],
                split_contract={
                    "fit_selection_calibration_dataset": "GSE63577",
                    "first_external_lab_test": "GSE297406",
                    "postfreeze_confirmation_dataset": senescence_confirmation_report[
                        "dataset"
                    ],
                    "postfreeze_confirmation_lab": senescence_confirmation_report[
                        "lab"
                    ],
                    "confirmation_pair_count": confirmation["pair_count"],
                    "external_labels_used_for_model_selection": False,
                    "absolute_sample_state_use_refused": True,
                },
                metrics={
                    "first_external_lab_test": external_senescence_report[
                        "sealed_external_dataset_lab_test"
                    ],
                    "postfreeze_directional_confirmation": senescence_confirmation_report,
                },
                status=(
                    "validated_external_direction_only_absolute_state_OOD_blocked"
                    if directional_passed else
                    "blocked_external_direction_confirmation"
                ),
                evidence_eligible=directional_passed,
                uncertainty_eligible=False,
                external_holdout=True,
                pristine_confirmation=True,
                limitations=(
                    external_senescence_report["limitations"]
                    + senescence_confirmation_report["limitations"]
                ),
            )
        )

    if pbmc_cross_lab_report is not None:
        cell_type = pbmc_cross_lab_report["cell_type"]
        passes = cell_type["frozen_endpoint_passes"]
        if pbmc_cross_lab_report["overall_external_context_gate_pass"] != all(passes.values()):
            raise ValueError("PBMC aggregate gate disagrees with frozen endpoint checks")
        direction_passed = bool(
            pbmc_ifn_direction_report
            and pbmc_ifn_direction_report["all_three_paired_donor_differences_positive"]
            and not pbmc_ifn_direction_report["numeric_scale_transfer_permitted"]
        )
        final_confirmation_passed = bool(
            pbmc_final_confirmation_report
            and pbmc_final_confirmation_report["single_use_confirmation_consumed"]
            and pbmc_final_confirmation_report["all_frozen_endpoints_passed"]
        )
        if pbmc_final_confirmation_report is not None:
            if pbmc_final_confirmation_report["may_be_reused_as_pristine"]:
                raise ValueError("consumed PBMC confirmation claims pristine reuse")
            if (
                (pbmc_final_confirmation_report["status"] == "passed")
                != final_confirmation_passed
            ):
                raise ValueError("PBMC final confirmation status disagrees with frozen endpoints")
        heads.append(
            _head(
                head_id="PBMC_cell_type_IFN_pre_sweep_context",
                task="PBMC_cell_type_and_IFN_context_conditioning",
                kind="single_cell_RNA_multinomial_classifier",
                modalities=["single_cell_RNA_seq", "bulk_RNA_seq_direction_confirmation"],
                split_contract={
                    "fit_selection_unit": "donor",
                    "calibration_unit": "sealed_donor",
                    "internal_test_unit": "sealed_donor",
                    "external_test_dataset_lab": "GSE132044_Broad",
                    "external_labels_used_for_selection": False,
                    "IFN_direction_external_dataset_lab": "GSE72502_Genentech",
                    "uncertainty_calibration_dataset_lab": "GSE164378_NYGC",
                    "single_use_final_confirmation_dataset_lab": (
                        "GSE150728_Stanford_Wilk"
                        if pbmc_final_confirmation_report is not None else None
                    ),
                },
                metrics={
                    "sealed_internal_donor": cell_type["internal"]["sealed_internal_donor"],
                    "pristine_external_lab": cell_type["pristine_external_lab"],
                    "semantic_ood": cell_type["semantic_ood"],
                    "frozen_endpoint_passes": passes,
                    "IFN_direction_confirmation": pbmc_ifn_direction_report,
                    "postcalibration_single_use_confirmation": pbmc_final_confirmation_report,
                },
                status=(
                    "validated_external" if final_confirmation_passed
                    else "blocked_external_uncertainty"
                ),
                evidence_eligible=False,
                uncertainty_eligible=final_confirmation_passed,
                external_holdout=True,
                pristine_confirmation=True,
                limitations=[
                    (
                        "postcalibration Stanford confirmation passes every frozen endpoint"
                        if final_confirmation_passed else
                        "postcalibration Stanford confirmation passes macro-F1, ECE, and semantic OOD but fails marginal and class-conditional conformal coverage"
                    ),
                    "cell rows are not independent biological replicates",
                    "IFN direction is confirmed in three paired donors" if direction_passed else "IFN cross-lab direction confirmation is absent or failed",
                    "the head may rank a future context-conditioned sweep but cannot delete axes or emit numeric ranges",
                ],
            )
        )

    if idr0072_if_report is not None:
        if idr0168_decode_refusal is not None:
            if (
                idr0168_decode_refusal["status"]
                != "failed_source_integrity_and_permanently_consumed"
                or not idr0168_decode_refusal["single_use_confirmation_consumed"]
                or idr0168_decode_refusal["may_be_reused_as_pristine"]
                or idr0168_decode_refusal["model_predictions_made"] != 0
                or idr0168_decode_refusal["replacement_field_used"]
                or idr0168_decode_refusal["external_provider_validation_complete"]
                or idr0168_decode_refusal["Aleph_parameter_authority"] != "none"
            ):
                raise ValueError("IDR0168 source-integrity refusal was rewritten")
        endpoints = idr0072_if_report["endpoint_results"]
        if idr0072_if_report["all_applicable_endpoints_passed"] != all(
            endpoints.values()
        ):
            raise ValueError("IDR0072 aggregate result disagrees with frozen endpoints")
        if any(endpoints.values()):
            raise ValueError("IDR0072 diagnostic unexpectedly passed a frozen endpoint")
        if idr0072_if_report["status"] != "domain_shift_diagnostic_failed":
            raise ValueError("IDR0072 failed transfer status was rewritten")
        if (
            idr0072_if_report["pristine_holdout_authority"]
            or idr0072_if_report["uncertainty_authority"]
            or idr0072_if_report["may_emit_numeric_sweep_range"]
        ):
            raise ValueError("consumed IDR0072 diagnostic claims forbidden authority")
        heads.append(
            _head(
                head_id="raw_IF_subcellular_localization_transfer",
                task="cross_provider_subcellular_localization",
                kind="raw_IF_multilabel_CNN_external_transfer_diagnostic",
                modalities=["two_channel_confocal_IF"],
                split_contract={
                    "training_provider": "Human_Protein_Atlas",
                    "external_diagnostic_provider": "IDR0072_Subcellular_Protein_Atlas",
                    "evaluated_image_count": idr0072_if_report["evaluated_image_count"],
                    "pre_prediction_refusal_count": idr0072_if_report[
                        "pre_prediction_refusal_count"
                    ],
                    "external_labels_used_for_checkpoint_or_threshold_fitting": False,
                    "diagnostic_consumed_and_may_not_be_reused_as_pristine": True,
                    "schema_optics_and_source_resolution_were_informed_before_evaluation": True,
                    "single_use_confirmation_provider": (
                        "IDR0168_PUPS"
                        if idr0168_decode_refusal is not None else None
                    ),
                    "single_use_confirmation_consumed_before_prediction": bool(
                        idr0168_decode_refusal
                    ),
                },
                metrics={
                    "frozen_endpoint_results": endpoints,
                    "image": idr0072_if_report["image_metrics"],
                    "landmark_plate_group": idr0072_if_report[
                        "landmark_plate_group_metrics"
                    ],
                    "species_screen": idr0072_if_report["species_screen_metrics"],
                    "IDR0168_single_use_confirmation": idr0168_decode_refusal,
                },
                status="blocked_external_domain_shift",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=True,
                pristine_confirmation=False,
                limitations=[
                    "all five frozen transfer and uncertainty endpoints failed",
                    "IDR0072 is a consumed schema- and optics-informed diagnostic, not a pristine confirmation set",
                    "landmark/plate aggregation does not create independent biological replicates",
                    "the failed diagnostic may guide only qualitative architecture hypotheses and may not tune a replacement model",
                    (
                        "the frozen IDR0168 confirmation was consumed before prediction because one official archive failed gzip CRC and lacked native Zarr metadata"
                        if idr0168_decode_refusal is not None else
                        "a new untouched provider confirmation is still required"
                    ),
                ],
            )
        )

    if lightmycells_report is not None:
        passes = lightmycells_report["frozen_endpoint_passes"]
        if lightmycells_report["all_applicable_confirmation_endpoints_passed"] != all(
            passes.values()
        ):
            raise ValueError(
                "Light My Cells aggregate result disagrees with frozen endpoints"
            )
        if lightmycells_report["status"] != "blocked_external_spatial_transfer":
            raise ValueError("Light My Cells failed transfer status was rewritten")
        if passes["macro_study_Pearson_correlation"] or passes[
            "embedding_OOD_AUROC"
        ]:
            raise ValueError("Light My Cells failed endpoints were rewritten")
        if (
            lightmycells_report["observable_evidence_eligible"]
            or lightmycells_report["uncertainty_authority_for_this_spatial_task"]
            or lightmycells_report["may_emit_numeric_sweep_range"]
        ):
            raise ValueError("failed Light My Cells head claims forbidden authority")
        heads.append(
            _head(
                head_id="label_free_organelle_spatial_prediction",
                task="transmitted_light_to_organelle_fluorescence_spatial_prediction",
                kind="conditional_residual_UNet_external_transfer",
                modalities=["BF", "PC", "DIC", "organelle_fluorescence"],
                split_contract={
                    "fit_study_count": 16,
                    "selection_study_count": 4,
                    "calibration_study_count": 4,
                    "single_use_confirmation_study_count": lightmycells_report[
                        "confirmation_study_count"
                    ],
                    "study_overlap_permitted": False,
                    "all_confirmation_labs_independent": lightmycells_report[
                        "all_confirmation_labs_independent"
                    ],
                    "novel_experimenter_group_count": lightmycells_report[
                        "novel_experimenter_group_partial_holdout"
                    ]["group_count"],
                    "confirmation_consumed_and_may_not_be_reused_as_pristine": True,
                    "confirmation_labels_used_for_fitting_selection_or_calibration": False,
                },
                metrics={
                    "frozen_endpoint_values": lightmycells_report[
                        "frozen_endpoint_values"
                    ],
                    "frozen_endpoint_passes": passes,
                    "per_target": lightmycells_report["spatial_metrics"]["per_target"],
                    "uncertainty": lightmycells_report["uncertainty"]["per_target"],
                    "novel_experimenter_group_partial_holdout": lightmycells_report[
                        "novel_experimenter_group_partial_holdout"
                    ],
                },
                status="blocked_external_spatial_transfer",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=True,
                pristine_confirmation=True,
                limitations=[
                    "single-use cross-study confirmation missed Pearson and embedding-OOD endpoints",
                    "three of five confirmation experimenter groups overlap training experimenter groups",
                    "Tubulin and Actin confirmation support is below the frozen 25-acquisition rule",
                    "the consumed confirmation partition may not tune a replacement model",
                    "this visual observation head has no Aleph parameter or numeric sweep authority",
                ],
            )
        )

    if mechanics_protocol_report is not None:
        wire = mechanics_protocol_report["wire_head"]
        router = mechanics_protocol_report["protocol_router"]
        run_disjoint = wire["experimental_date_disjoint"]
        if run_disjoint["split_group_overlap_count"]:
            raise ValueError("mechanics wire head leaks experiment dates")
        if wire["production_eligible"] or router["external_holdout"]:
            raise ValueError("mechanics protocol head claims unsupported production/external authority")
        if mechanics_protocol_report["observable_evidence_eligible"]:
            raise ValueError("mechanics protocol representation claims observable evidence authority")
        heads.append(
            _head(
                head_id="mechanics_protocol_and_probe_scale_representation",
                task=mechanics_protocol_report["task"],
                kind="protocol_router_plus_wire_multihead_network",
                modalities=[
                    "cell_mechanics_protocol_metadata",
                    "magnetic_wire_dynamic_response",
                ],
                split_contract={
                    "wire_split_unit": "experimental_date",
                    "wire_split_group_overlap_count": 0,
                    "individual_wire_source_count": 1,
                    "catalog_source_count": mechanics_protocol_report["dataset"][
                        "catalog_source_count"
                    ],
                    "protocol_router_external_holdout": False,
                    "cross_scope_interpolation_permitted": False,
                },
                metrics={
                    "wire_experimental_date_disjoint": run_disjoint["metrics"],
                    "wire_random_interpolation_diagnostic": wire[
                        "random_wire_holdout_interpolation_only"
                    ],
                    "protocol_router": router,
                    "probe_size_dependence": mechanics_protocol_report[
                        "probe_size_dependence"
                    ],
                },
                status="training_only_operator_representation_cell_line_transfer_failed",
                evidence_eligible=False,
                uncertainty_eligible=False,
                external_holdout=False,
                pristine_confirmation=False,
                limitations=mechanics_protocol_report["refusals"],
            )
        )

    passed_external_evidence = [
        item["head_id"] for item in heads
        if item["observable_evidence_eligible"] and item["external_holdout"]
    ]
    blocked_classifiers = [
        item["head_id"] for item in heads
        if item["status"].startswith("blocked_external")
    ]
    training_only = [
        item["head_id"] for item in heads
        if item["status"].startswith("training_only")
    ]
    pristine_external_classifier_heads = [
        item["head_id"] for item in heads
        if item["kind"].endswith("classifier")
        and item["status"] == "validated_external"
        and item["pristine_confirmation"]
    ]
    production_eligible = bool(pristine_external_classifier_heads)
    if production_eligible:
        # Kept as an explicit invariant: adding one successful head must not
        # silently grant global cell-type/state production authority.
        production_eligible = all(
            item["status"] == "validated_external"
            for item in heads if item["kind"].endswith("classifier")
        )

    return {
        "schema": SCHEMA,
        "model_family": "task_scoped_multimodal_multihead",
        "architecture": {
            "shared_latent": False,
            "reason": (
                "modalities lack enough aligned samples for a leakage-safe shared latent; "
                "each head keeps its own measurement and split semantics"
            ),
            "fusion": "task_scoped_late_fusion_only",
            "head_count": len(heads),
            "trainable_raw_encoder_status": (
                "trained"
                if lightmycells_report is not None else
                "three_encoders_trained_external_IF_failed"
                if (
                    temporal_cnn_report is not None
                    and cell_cycle_report is not None
                    and idr0072_if_report is not None
                ) else "two_encoders_trained_no_external_acceptance"
                if temporal_cnn_report is not None and cell_cycle_report is not None
                else "one_temporal_encoder_trained_rejected_external"
                if temporal_cnn_report is not None else "not_yet_trained"
            ),
            "trained_encoder_heads": [
                head_id for head_id, present in (
                    ("calcium_cell_state", temporal_cnn_report is not None),
                    ("cell_cycle_image_representation", cell_cycle_report is not None),
                    (
                        "hpa_if_localization_representation",
                        idr0072_if_report is not None,
                    ),
                    (
                        "label_free_organelle_spatial_prediction",
                        lightmycells_report is not None,
                    ),
                    (
                        "sciplex_perturbation_state_representation",
                        sciplex_state_report is not None,
                    ),
                    (
                        "senscout_morphology_senescence_representation",
                        senscout_report is not None,
                    ),
                    (
                        "programmed_cell_death_mechanism_representation",
                        figshare_cell_death_report is not None,
                    ),
                    (
                        "senescence_transcriptomic_direction",
                        external_senescence_report is not None,
                    ),
                    ("PBMC_cell_type_IFN_pre_sweep_context", pbmc_cross_lab_report is not None),
                    (
                        "mechanics_protocol_and_probe_scale_representation",
                        mechanics_protocol_report is not None,
                    ),
                ) if present
            ],
            "gpu_required_for_current_report": False,
        },
        "heads": heads,
        "readiness": {
            "task_scoped_evidence_runtime_ready": bool(passed_external_evidence),
            "general_cell_type_state_model_ready": False,
            "production_eligible": production_eligible,
            "passed_external_evidence_heads": passed_external_evidence,
            "blocked_external_classifier_heads": blocked_classifiers,
            "training_only_heads": training_only,
            "pristine_external_classifier_heads": pristine_external_classifier_heads,
            "next_model_requirement": (
                ("improve sci-Plex unseen-compound pathway transfer; "
                 if sciplex_state_report is not None else
                 "finish the verified sci-Plex expression tensor; ")
                + "acquire pristine "
                "aligned external labs and biological replicates for cell cycle, "
                "mitochondrial stress, calcium, TFM, and PIV before any new "
                "encoder can earn evidence authority; replace scalar PBMC calibration "
                "with a preregistered multi-source/domain-conditional uncertainty model "
                "and reserve new calibration and confirmation labs; for IF, train a "
                "domain-robust multi-provider encoder on providers other than the consumed "
                + (
                    "IDR0072 diagnostic; IDR0168 was consumed by official-source integrity failure before prediction, so reserve another untouched provider for confirmation; "
                    if idr0168_decode_refusal is not None else
                    "IDR0072 diagnostic and reserve a new untouched provider for confirmation; "
                )
                + "for label-free organelle prediction, replace the consumed Light My Cells "
                "confirmation with a newly preregistered external provider and improve "
                "cross-study correlation plus semantic OOD separation"
            ),
        },
        "sweep_interface": {
            "status": "evidence_only",
            "admissible_head_ids": passed_external_evidence,
            "may_emit_numeric_parameter_range": False,
            "may_select_aleph_parameter": False,
            "may_mutate_physics": False,
            "required_downstream_gate": "aleph.outer_library.sweep_gate",
        },
        "aleph_authority": "none",
    }


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--calcium-report", type=Path, required=True)
    parser.add_argument("--piezo1-report", type=Path, required=True)
    parser.add_argument("--calcium-uncertainty-report", type=Path, required=True)
    parser.add_argument("--mechanism-report", type=Path, required=True)
    parser.add_argument("--migration-report", type=Path, required=True)
    parser.add_argument("--dryad-migration-report", type=Path, required=True)
    parser.add_argument("--fibrotic-report", type=Path, required=True)
    parser.add_argument("--mechano-osmotic-report", type=Path, required=True)
    parser.add_argument("--piv-report", type=Path, required=True)
    parser.add_argument("--cell-monolayer-report", type=Path, required=True)
    parser.add_argument("--temporal-cnn-report", type=Path)
    parser.add_argument("--if-translocation-report", type=Path)
    parser.add_argument("--hpa-if-report", type=Path)
    parser.add_argument("--microglia-report", type=Path)
    parser.add_argument("--cell-cycle-report", type=Path)
    parser.add_argument("--mitochondrial-stress-report", type=Path)
    parser.add_argument("--sciplex-design-report", type=Path)
    parser.add_argument("--apoptosis-report", type=Path)
    parser.add_argument("--sciplex-state-report", type=Path)
    parser.add_argument("--senscout-report", type=Path)
    parser.add_argument("--figshare-cell-death-report", type=Path)
    parser.add_argument("--external-senescence-report", type=Path)
    parser.add_argument("--senescence-confirmation-report", type=Path)
    parser.add_argument("--tfm-tether-report", type=Path)
    parser.add_argument("--piv-contractility-report", type=Path)
    parser.add_argument("--nematic-tfm-report", type=Path)
    parser.add_argument("--opto-force-report", type=Path)
    parser.add_argument("--tfm-fret-report", type=Path)
    parser.add_argument("--exif-emt-report", type=Path)
    parser.add_argument("--mtrack-emt-report", type=Path)
    parser.add_argument("--lincs-mcf10a-report", type=Path)
    parser.add_argument("--karacosta-emt-report", type=Path)
    parser.add_argument("--geo-emt-rna-report", type=Path)
    parser.add_argument("--gse325309-emt-report", type=Path)
    parser.add_argument("--pbmc-cross-lab-report", type=Path)
    parser.add_argument("--pbmc-ifn-direction-report", type=Path)
    parser.add_argument("--pbmc-final-confirmation-report", type=Path)
    parser.add_argument("--idr0072-if-report", type=Path)
    parser.add_argument("--idr0168-decode-refusal", type=Path)
    parser.add_argument("--lightmycells-report", type=Path)
    parser.add_argument("--mechanics-protocol-report", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(
        *[
            _read(getattr(args, name))
            for name in (
                "calcium_report",
                "piezo1_report",
                "calcium_uncertainty_report",
                "mechanism_report",
                "migration_report",
                "dryad_migration_report",
                "fibrotic_report",
                "mechano_osmotic_report",
                "piv_report",
                "cell_monolayer_report",
            )
        ],
        (_read(args.temporal_cnn_report) if args.temporal_cnn_report else None),
        (
            _read(args.if_translocation_report)
            if args.if_translocation_report else None
        ),
        (_read(args.hpa_if_report) if args.hpa_if_report else None),
        (_read(args.microglia_report) if args.microglia_report else None),
        (_read(args.cell_cycle_report) if args.cell_cycle_report else None),
        (
            _read(args.mitochondrial_stress_report)
            if args.mitochondrial_stress_report else None
        ),
        (_read(args.sciplex_design_report) if args.sciplex_design_report else None),
        (_read(args.apoptosis_report) if args.apoptosis_report else None),
        (_read(args.sciplex_state_report) if args.sciplex_state_report else None),
        (_read(args.senscout_report) if args.senscout_report else None),
        (
            _read(args.figshare_cell_death_report)
            if args.figshare_cell_death_report else None
        ),
        (
            _read(args.external_senescence_report)
            if args.external_senescence_report else None
        ),
        (
            _read(args.senescence_confirmation_report)
            if args.senescence_confirmation_report else None
        ),
        (_read(args.tfm_tether_report) if args.tfm_tether_report else None),
        (
            _read(args.piv_contractility_report)
            if args.piv_contractility_report else None
        ),
        (_read(args.nematic_tfm_report) if args.nematic_tfm_report else None),
        (_read(args.opto_force_report) if args.opto_force_report else None),
        (_read(args.tfm_fret_report) if args.tfm_fret_report else None),
        (_read(args.exif_emt_report) if args.exif_emt_report else None),
        (_read(args.mtrack_emt_report) if args.mtrack_emt_report else None),
        (_read(args.lincs_mcf10a_report) if args.lincs_mcf10a_report else None),
        (_read(args.karacosta_emt_report) if args.karacosta_emt_report else None),
        (_read(args.geo_emt_rna_report) if args.geo_emt_rna_report else None),
        (_read(args.gse325309_emt_report) if args.gse325309_emt_report else None),
        (_read(args.pbmc_cross_lab_report) if args.pbmc_cross_lab_report else None),
        (_read(args.pbmc_ifn_direction_report) if args.pbmc_ifn_direction_report else None),
        (
            _read(args.pbmc_final_confirmation_report)
            if args.pbmc_final_confirmation_report else None
        ),
        (_read(args.idr0072_if_report) if args.idr0072_if_report else None),
        (_read(args.lightmycells_report) if args.lightmycells_report else None),
        (
            _read(args.idr0168_decode_refusal)
            if args.idr0168_decode_refusal else None
        ),
        (
            _read(args.mechanics_protocol_report)
            if args.mechanics_protocol_report else None
        ),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result["readiness"], sort_keys=True))


if __name__ == "__main__":
    main()
