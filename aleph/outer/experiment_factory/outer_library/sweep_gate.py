#!/usr/bin/env python3
"""Evidence gate between Outer Library predictions and Aleph sweep proposals."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def evaluate(
    report: dict,
    mechanism_report: dict | None = None,
    tendon_report: dict | None = None,
    fibrotic_report: dict | None = None,
    uncertainty_report: dict | None = None,
    migration_report: dict | None = None,
    rock_context_report: dict | None = None,
    mechano_osmotic_report: dict | None = None,
    cell_monolayer_report: dict | None = None,
    multitask_model_report: dict | None = None,
    piv_contractility_report: dict | None = None,
    nematic_tfm_report: dict | None = None,
    opto_force_report: dict | None = None,
    tfm_fret_report: dict | None = None,
    exif_emt_report: dict | None = None,
    mtrack_emt_report: dict | None = None,
    lincs_mcf10a_report: dict | None = None,
    karacosta_emt_report: dict | None = None,
    geo_emt_rna_report: dict | None = None,
    gse325309_emt_report: dict | None = None,
    observation_operator_report: dict | None = None,
) -> dict:
    checks = {
        "macro_f1_at_least_0_70": report["test"]["macro_f1"] >= 0.70,
        "ece_at_most_0_15": report["test"]["ece_10_bin"] <= 0.15,
        "semantic_ood_auroc_at_least_0_70": report["semantic_ood"]["mean_auroc"] >= 0.70,
        "independent_dataset_task_holdout_present": (
            report["independent_dataset_task_holdout_count"] >= 1
        ),
        "independent_lab_task_holdout_present": (
            report["independent_lab_task_holdout_count"] >= 1
        ),
        "independent_dataset_task_holdout_passed": (
            report["independent_dataset_task_holdout_passed_count"] >= 1
        ),
        "independent_lab_task_holdout_passed": (
            report["independent_lab_task_holdout_passed_count"] >= 1
        ),
        "model_has_no_aleph_authority": report["aleph_authority"] == "none",
    }
    operator_by_observable = {}
    if observation_operator_report is not None:
        operator_by_observable = {
            item["external_observable"]: item
            for item in observation_operator_report.get("entries", [])
        }
        checks["observation_operator_registry_passed"] = bool(
            observation_operator_report.get("passed")
            and not observation_operator_report.get("may_emit_numeric_parameter_range")
            and not observation_operator_report.get("may_select_aleph_parameter")
        )
        checks["EMT_marker_operators_are_explicitly_unavailable"] = all(
            operator_by_observable.get(name, {}).get("status") == "unavailable"
            and not operator_by_observable[name].get("promotable_to_numeric_sweep")
            for name in (
                "CDH1_RNA_abundance", "VIM_RNA_abundance",
                "E_cadherin_protein_abundance", "Vimentin_protein_abundance",
            )
        )
    if uncertainty_report is not None:
        checks["external_lab_uncertainty_holdout_passed"] = bool(
            uncertainty_report["uncertainty_holdout_passed"]
        )
    if migration_report is not None:
        checks["independent_multisite_migration_external_evaluation_passed"] = bool(
            migration_report["external_lab_holdout_passed"]
        )
        checks["independent_multisite_migration_retrospective_uncertainty_passed"] = bool(
            migration_report["conformal_uncertainty"]["uncertainty_holdout_passed"]
        )
        checks["independent_multisite_migration_pristine_confirmation_present"] = bool(
            migration_report["pristine_sealed_confirmation"]
        )
    multitask_model_audit = None
    if multitask_model_report is not None:
        sweep_interface = multitask_model_report["sweep_interface"]
        safe_contract = bool(
            multitask_model_report["schema"]
            == "aleph.outer_library.multitask_model.v1"
            and multitask_model_report["aleph_authority"] == "none"
            and sweep_interface["status"] == "evidence_only"
            and not sweep_interface["may_emit_numeric_parameter_range"]
            and not sweep_interface["may_select_aleph_parameter"]
            and not sweep_interface["may_mutate_physics"]
            and all(
                not head["parameter_proposal_eligible"]
                and head["aleph_authority"] == "none"
                for head in multitask_model_report["heads"]
            )
        )
        checks["multitask_model_contract_is_non_authoritative"] = safe_contract
        admitted = set(sweep_interface["admissible_head_ids"])
        declared = {
            head["head_id"] for head in multitask_model_report["heads"]
            if head["observable_evidence_eligible"]
            and head["external_holdout"]
            and head["status"].startswith("validated")
        }
        checks["multitask_model_admits_only_externally_validated_evidence_heads"] = (
            admitted == declared
        )
        multitask_model_audit = {
            "model_family": multitask_model_report["model_family"],
            "task_scoped_evidence_runtime_ready": multitask_model_report[
                "readiness"
            ]["task_scoped_evidence_runtime_ready"],
            "general_cell_type_state_model_ready": multitask_model_report[
                "readiness"
            ]["general_cell_type_state_model_ready"],
            "production_eligible": multitask_model_report["readiness"][
                "production_eligible"
            ],
            "admitted_evidence_heads": sorted(admitted),
            "blocked_external_classifier_heads": multitask_model_report[
                "readiness"
            ]["blocked_external_classifier_heads"],
            "may_select_aleph_parameter": False,
        }
    if piv_contractility_report is not None:
        constraint = piv_contractility_report["pre_sweep_constraint"]
        checks["piv_contractility_constraint_is_non_authoritative"] = bool(
            piv_contractility_report["aleph_authority"] == "none"
            and not piv_contractility_report["observable_evidence_eligible"]
            and not piv_contractility_report["uncertainty_holdout_passed"]
            and not piv_contractility_report["independent_lab_holdout"]
            and not piv_contractility_report["aleph_parameter_proposal_eligible"]
            and constraint["future_sweep_range_reduction_intended"]
            and not constraint["may_emit_numeric_parameter_range"]
            and not constraint["may_select_Aleph_parameter"]
        )
    for name, force_report in (
        ("nematic_tfm", nematic_tfm_report),
        ("opto_force", opto_force_report),
        ("tfm_fret", tfm_fret_report),
    ):
        if force_report is None:
            continue
        constraint = force_report["pre_sweep_constraint"]
        checks[f"{name}_constraint_is_non_authoritative"] = bool(
            force_report["aleph_authority"] == "none"
            and not force_report["observable_evidence_eligible"]
            and not force_report["uncertainty_holdout_passed"]
            and not force_report["independent_lab_holdout"]
            and not force_report["aleph_parameter_proposal_eligible"]
            and constraint["future_sweep_range_reduction_intended"]
            and not constraint["may_emit_numeric_parameter_range"]
            and not constraint["may_select_Aleph_parameter"]
        )
    observable_channel = {
        "status": "not_evaluated",
        "constraints": [],
    }
    evidence_exclusions = []
    if rock_context_report is not None and rock_context_report[
        "context_heterogeneity_detected"
    ]:
        evidence_exclusions.append({
            "observable": "universal_ROCK_inhibition_migration_speed_direction",
            "status": "rejected_context_dependent",
            "reason": rock_context_report["interpretation"],
            "dryad_2D_direction": rock_context_report["dryad_H1299_PL_2D"]["direction"],
            "elife_3D_CDM_direction": rock_context_report["elife_NIH3T3_3D_CDM"]["direction"],
            "may_select_aleph_parameter": False,
        })
    if mechano_osmotic_report is not None:
        evidence_exclusions.append({
            "observable": "AFM_membrane_tether_force_as_direct_cortical_tension",
            "status": "rejected_operator_mismatch",
            "reason": (
                "Tether force is a pN-scale proxy for apparent membrane tension only "
                "under a constant-bending-rigidity tether model; it is not a direct "
                "measurement of Aleph cortical or Helfrich surface tension."
            ),
            "may_select_aleph_parameter": False,
        })
    if piv_contractility_report is not None:
        evidence_exclusions.append({
            "observable": "C2C12_blebbistatin_PIV_IF_direction_as_one_Aleph_parameter",
            "status": "rejected_pending_operator_and_aligned_external_lab",
            "reason": piv_contractility_report["pre_sweep_constraint"]["identifiability"],
            "replicated_author_direction": piv_contractility_report[
                "replicated_contractility_direction"
            ],
            "may_select_aleph_parameter": False,
        })
    if nematic_tfm_report is not None:
        evidence_exclusions.append({
            "observable": "NIH3T3_blebbistatin_TFM_MSM_as_one_Aleph_parameter",
            "status": "rejected_pending_biological_replication_and_forward_operator",
            "reason": nematic_tfm_report["pre_sweep_constraint"]["identifiability"],
            "published_direction_reproduced": nematic_tfm_report[
                "published_direction_reproduced_from_source_tables"
            ],
            "may_select_aleph_parameter": False,
        })
    if opto_force_report is not None:
        evidence_exclusions.append({
            "observable": "opto_MDCK_RhoA_contour_response_as_one_Aleph_parameter",
            "status": "rejected_pending_sample_hierarchy_and_forward_operator",
            "reason": opto_force_report["pre_sweep_constraint"]["identifiability"],
            "usually_positive_in_both_axes": opto_force_report[
                "overall_direction"
            ]["usually_positive_in_both_axes"],
            "may_select_aleph_parameter": False,
        })
    if mechanism_report is not None:
        eligible = bool(mechanism_report["observable_constraint_eligible"])
        external = mechanism_report["external_confirmatory"]
        observable_channel = {
            "status": "evidence_constraint_eligible" if eligible else "refused",
            "constraints": ([{
                "observable": "intracellular_calcium_response_AUC",
                "direction": "decreases_with_GsMTx4_under_active_PIEZO1_nanoswitch_stimulation",
                "external_equal_stratum_mean_difference": external[
                    "equal_stratum_mean_difference"
                ],
                "external_bootstrap_ci95": external["equal_stratum_bootstrap_ci95"],
                "scope": "HeLa discovery; EA.hy926 LooPINS/CaPINS confirmation",
                "authority": "evidence_only",
            }] if eligible else []),
            "may_select_aleph_parameter": False,
        }
    if fibrotic_report is not None and fibrotic_report["observable_constraint_eligible"]:
        external = fibrotic_report["external_confirmatory"]
        observable_channel["status"] = "evidence_constraint_eligible"
        observable_channel.setdefault("constraints", []).append({
            "observable": "TGF_beta_1_fibrotic_state_signature",
            "direction": "increases_with_TGF_beta_1_and_reverses_with_celastrol",
            "external_TGF_beta_1_bootstrap_ci95": external[
                "TGF_beta_1_minus_control"
            ]["bootstrap_ci95"],
            "external_celastrol_reversal_bootstrap_ci95": external[
                "celastrol_reversal_TGF_minus_TGF_plus_celastrol"
            ]["bootstrap_ci95"],
            "scope": "rat tendon IF discovery; independent human dermal fibroblast bulk-RNA confirmation",
            "authority": "evidence_only",
        })
        observable_channel["may_select_aleph_parameter"] = False
    if migration_report is not None and migration_report["external_lab_holdout_passed"]:
        observable_channel["status"] = "evidence_constraint_eligible"
        observable_channel.setdefault("constraints", []).append({
            "observable": "multivariate_live_cell_ROCK_inhibition_signature",
            "direction": "positive_projection_on_lab1_paired_perturbation_direction",
            "external_lab3_metrics": migration_report["external_lab3_metrics"],
            "external_lab3_conformal_uncertainty": migration_report[
                "conformal_uncertainty"
            ],
            "external_lab3_experiment_effect_bootstrap_ci95": migration_report[
                "external_lab3_effect_bootstrap_95ci"
            ],
            "scope": "HT1080; three labs; matched-control assay only",
            "confirmation_status": migration_report["protocol"][
                "confirmation_status"
            ],
            "identifiability": "pathway_perturbation_not_one_Aleph_parameter",
            "authority": "evidence_only",
        })
        observable_channel["may_select_aleph_parameter"] = False
    if (mechano_osmotic_report is not None
            and mechano_osmotic_report["observable_constraint_eligible"]):
        early = mechano_osmotic_report["AFM_membrane_tether_force"][
            "spreading_phase_30_to_90_min"
        ]
        observable_channel["status"] = "evidence_constraint_eligible"
        observable_channel.setdefault("constraints", []).append({
            "observable": "HeLa_mechano_osmotic_spreading_and_AFM_tether_signature",
            "direction": (
                "Y27632_increases_initial_area_rate_decreases_volume_flux_and_"
                "transiently_increases_membrane_tether_force"
            ),
            "area_rate_experiment_bootstrap_ci95": mechano_osmotic_report[
                "initial_spreading_area_rate"
            ]["independent_experiment_bootstrap_95ci"],
            "volume_flux_experiment_bootstrap_ci95": mechano_osmotic_report[
                "initial_volume_flux"
            ]["independent_experiment_bootstrap_95ci"],
            "early_tether_force_paired_experiment_bootstrap_ci95_pN": early[
                "paired_experiment_bootstrap_95ci_pN"
            ],
            "steady_state_tether_force_experiment_ci_includes_zero": True,
            "scope": "HeLa Kyoto; fibronectin spreading; 100 micromolar Y-27632",
            "identifiability": "multimechanism_drug_contrast_not_one_Aleph_parameter",
            "authority": "evidence_only",
        })
        observable_channel["may_select_aleph_parameter"] = False
    parameter_axis_candidates = []
    calibration_requests = []
    representation_training_audits = []
    if cell_monolayer_report is not None:
        safely_grouped = bool(
            cell_monolayer_report["representation_training_eligible"]
            and cell_monolayer_report["raw_velocity_tensor"][
                "biological_replicate_count"
            ] == 1
            and not cell_monolayer_report["independent_lab_holdout"]
            and not cell_monolayer_report["aleph_parameter_proposal_eligible"]
            and cell_monolayer_report["aleph_authority"] == "none"
        )
        representation_training_audits.append({
            "dataset": "dryad-d96w58-hdf-collective-velocity",
            "status": "admitted_training_only" if safely_grouped else "refused",
            "technical_vector_count": cell_monolayer_report[
                "raw_velocity_tensor"
            ]["technical_vectors"],
            "mandatory_split_group": cell_monolayer_report[
                "raw_velocity_tensor"
            ]["mandatory_split_group"],
            "may_count_as_external_holdout": False,
            "may_select_aleph_parameter": False,
            "reason": (
                "spatial vectors and longitudinal points are useful for grouped "
                "representation learning but the source reports no independent "
                "culture identifiers"
            ),
        })
    if piv_contractility_report is not None:
        representation_training_audits.append({
            "dataset": "dryad-08kprr59c-c2c12-supracontractility",
            "status": "admitted_training_direction_only",
            "biological_substrates": piv_contractility_report["manifest"][
                "author_independent_LCN_substrates"
            ],
            "replicated_contractility_direction": piv_contractility_report[
                "replicated_contractility_direction"
            ],
            "may_count_as_external_holdout": False,
            "may_select_aleph_parameter": False,
            "reason": (
                "exact-unit replicated PIV/IF direction is useful for pre-sweep "
                "ranking, but n=2 per intervention arm and no aligned external "
                "lab or Aleph observation operator prevent numeric calibration"
            ),
        })
        calibration_requests.append({
            "parameter": "mechanism_family:supracellular_contractility",
            "candidate_parameters": [
                "aleph.vertical.sf_arc.MaterialCard.active_tension_pn",
                "aleph.vertical.nmii.MotorKinetics.crossbridge_stiffness",
                "aleph.vertical.nmii.MotorKinetics.attachment_rate",
                "aleph.vertical.focal_adhesion.ClutchCard.bond_stiffness_pn_per_um",
            ],
            "status": "blocked_directional_constraint_only",
            "numeric_range": None,
            "required_forward_observables": piv_contractility_report[
                "pre_sweep_constraint"
            ]["required_Aleph_operators"],
            "required_before_unblock": [
                "implement unit- and geometry-matched Aleph order/correlation observation operators",
                "pass grouped biological-replicate validation in an aligned independent lab",
                "separate nonmuscle-myosin, stress-fibre, focal-adhesion, and ECM-remodelling contributions",
                "show monotonic replicated forward sensitivity over a declared safe range",
            ],
            "may_emit_sweep_axis": False,
        })
    if nematic_tfm_report is not None or opto_force_report is not None:
        reports = [
            item for item in (nematic_tfm_report, opto_force_report)
            if item is not None
        ]
        labs = {item["lab_group"] for item in reports}
        if piv_contractility_report is not None:
            labs.add("anseth-white-skillin-lab-cu-boulder")
        representation_training_audits.append({
            "dataset": "cross-cell-type-actomyosin-contractility",
            "status": "admitted_directional_pre_sweep_ranking_only",
            "lab_groups": sorted(labs),
            "cell_types": [
                "C2C12_mouse_myoblast_myotube",
                "NIH3T3_mouse_fibroblast",
                "opto_MDCK_canine_epithelial",
            ],
            "aligned_same_observable_holdout": False,
            "may_count_as_external_holdout": False,
            "may_select_aleph_parameter": False,
            "reason": (
                "three labs give directionally compatible intervention evidence, but "
                "different observables and unresolved biological hierarchies do not form "
                "a numeric cross-lab calibration"
            ),
        })
        required_operators = sorted({
            operator for item in reports
            for operator in item["pre_sweep_constraint"]["required_Aleph_operators"]
        })
        if piv_contractility_report is not None:
            required_operators += [
                operator for operator in piv_contractility_report[
                    "pre_sweep_constraint"
                ]["required_Aleph_operators"] if operator not in required_operators
            ]
        calibration_requests.append({
            "parameter": "mechanism_family:cross_cell_type_actomyosin_contractility",
            "candidate_parameters": [
                "aleph.vertical.sf_arc.MaterialCard.active_tension_pn",
                "aleph.vertical.nmii.MotorKinetics.crossbridge_stiffness",
                "aleph.vertical.nmii.MotorKinetics.attachment_rate",
                "aleph.vertical.focal_adhesion.ClutchCard.bond_stiffness_pn_per_um",
            ],
            "status": "blocked_directionally_ranked_but_not_numerically_calibrated",
            "numeric_range": None,
            "required_forward_observables": required_operators,
            "required_before_unblock": [
                "implement unit- and geometry-matched Aleph PIV, TFM, MSM, and contour operators",
                "acquire at least two independent biological cultures per intervention in an aligned external task",
                "pass dataset/lab holdout with calibrated uncertainty and OOD refusal",
                "show monotonic forward sensitivity and separate motor/cortex/stress-fibre/adhesion confounding",
            ],
            "may_emit_sweep_axis": False,
        })
    if tfm_fret_report is not None:
        constraint = tfm_fret_report["pre_sweep_constraint"]
        representation_training_audits.append({
            "dataset": tfm_fret_report["source"]["dataset_id"],
            "status": "admitted_cross_lab_direction_and_multiscale_identifiability_only",
            "lab_group": tfm_fret_report["lab_group"],
            "same_observable_external_direction_reproduced": tfm_fret_report[
                "cross_lab_same_observable_direction"
            ]["independent_lab_direction_reproduced"],
            "paired_TFM_FRET_categories_reproduced": tfm_fret_report[
                "author_correlation_categories_reproduced"
            ],
            "may_count_as_external_holdout": False,
            "may_select_aleph_parameter": False,
            "reason": (
                "mean traction in Pa reproduces the stiffness direction in an independent "
                "lab, while simultaneous FA traction/FRET exposes load-transfer "
                "nonidentifiability; released cells lack replicate mapping"
            ),
        })
        calibration_requests.append({
            "parameter": "mechanism_family:focal_adhesion_load_transfer",
            "candidate_parameters": constraint["candidate_parameters"],
            "status": "blocked_directionally_ranked_and_multiscale_identifiability_constrained",
            "numeric_range": None,
            "required_forward_observables": constraint["required_Aleph_operators"],
            "required_before_unblock": [
                "implement the author-equivalent cell-mask and focal-adhesion traction operators",
                "implement or validate a relative vinculin-tension observation bridge",
                "acquire released cell-to-biological-replicate mapping or a replicate-complete aligned lab",
                "pass calibrated dataset/lab uncertainty and OOD refusal",
                "separate active tension from clutch bond and anchor stiffness by forward sensitivity",
            ],
            "may_emit_sweep_axis": False,
        })
    if exif_emt_report is not None:
        role = exif_emt_report["pre_sweep_role"]
        safely_grouped = bool(
            exif_emt_report["aleph_authority"] == "none"
            and exif_emt_report["representation_training_eligible"]
            and not exif_emt_report["observable_evidence_eligible"]
            and not exif_emt_report["independent_dataset_holdout"]
            and not exif_emt_report["independent_lab_holdout"]
            and not exif_emt_report["uncertainty_holdout_passed"]
            and not exif_emt_report["OOD_refusal_validated"]
            and role["role"] == "future_state_conditioning_and_sweep_space_reduction"
            and role["numeric_range"] is None
            and not role["may_emit_sweep_axis"]
            and not exif_emt_report["may_select_aleph_parameter"]
        )
        checks["exif_emt_constraint_is_non_authoritative"] = safely_grouped
        representation_training_audits.append({
            "dataset": exif_emt_report["source"]["dataset_id"],
            "status": "admitted_training_only" if safely_grouped else "refused",
            "state_axis": exif_emt_report["state_axis"],
            "condition_labels_are_per_cell_ground_truth": False,
            "common_task_with_HPA_localization_labels": False,
            "future_sweep_role": (
                "condition the eventual Aleph sweep on EMT context and reject "
                "context-incompatible candidates; never replace the sweep"
            ),
            "candidate_mechanism_families": role["candidate_mechanism_families"],
            "may_count_as_external_holdout": False,
            "may_emit_numeric_parameter_range": False,
            "may_select_aleph_parameter": False,
            "reason": (
                "raw multichannel IF adds independent-provider EMT treatment "
                "representation, but one plate date and treatment-level labels "
                "cannot identify a physical parameter or validate uncertainty"
            ),
        })
    if mtrack_emt_report is not None:
        role = mtrack_emt_report["pre_sweep_role"]
        design = mtrack_emt_report["design_audit"]
        safely_grouped = bool(
            mtrack_emt_report["aleph_authority"] == "none"
            and mtrack_emt_report["representation_training_eligible"]
            and not mtrack_emt_report["observable_evidence_eligible"]
            and not mtrack_emt_report["independent_dataset_holdout"]
            and not mtrack_emt_report["independent_lab_holdout"]
            and not mtrack_emt_report["uncertainty_holdout_passed"]
            and not mtrack_emt_report["OOD_refusal_validated"]
            and design["trajectory_split_groups"] == 1
            and not design["frames_permitted_as_train_test_units"]
            and not design["biological_replicates_released"]
            and not design["untreated_control_released"]
            and role["role"] == "future_temporal_state_conditioning_and_sweep_space_reduction"
            and role["numeric_range"] is None
            and not role["may_emit_sweep_axis"]
            and not mtrack_emt_report["may_select_aleph_parameter"]
        )
        checks["mtrack_emt_constraint_is_non_authoritative"] = safely_grouped
        direction = mtrack_emt_report["cross_lab_vimentin_direction_context"]
        representation_training_audits.append({
            "dataset": mtrack_emt_report["source"]["dataset_id"],
            "status": "admitted_training_only" if safely_grouped else "refused",
            "state_axis": mtrack_emt_report["state_axis"],
            "single_trajectory": True,
            "may_count_frames_as_independent_samples": False,
            "may_count_as_external_holdout": False,
            "fixed_IF_vs_live_VIM_RFP_direction_agrees": (
                None if direction is None else direction["directions_agree"]
            ),
            "future_sweep_role": (
                "block intensity-only EMT sweep conditioning when live and fixed "
                "assays conflict; retain temporal morphology as training context"
            ),
            "candidate_mechanism_families": role["candidate_mechanism_families"],
            "may_emit_numeric_parameter_range": False,
            "may_select_aleph_parameter": False,
            "reason": (
                "one TGF-beta-treated live trajectory has no control or biological "
                "replicates and its late VIM-RFP intensity direction conflicts with "
                "the fixed-IF treatment contrast"
            ),
        })
    if lincs_mcf10a_report is not None:
        role = lincs_mcf10a_report["pre_sweep_role"]
        design = lincs_mcf10a_report["design_audit"]
        safely_bounded = bool(
            lincs_mcf10a_report["aleph_authority"] == "none"
            and lincs_mcf10a_report["collection_holdout_passed"]
            and lincs_mcf10a_report["observable_evidence_eligible"]
            and not lincs_mcf10a_report["independent_lab_holdout"]
            and not lincs_mcf10a_report["OOD_refusal_validated"]
            and design["technical_rows_aggregated_before_inference"]
            and design["C1_and_C2_are_distinct_collection_periods_same_OHSU_site"]
            and role["numeric_range"] is None
            and not role["may_emit_sweep_axis"]
            and not lincs_mcf10a_report["may_select_aleph_parameter"]
        )
        checks["lincs_mcf10a_constraint_is_non_authoritative"] = safely_bounded
        representation_training_audits.append({
            "dataset": lincs_mcf10a_report["source"]["dataset_id"],
            "status": (
                "admitted_collection_validated_training_only"
                if safely_bounded else "refused"
            ),
            "state_axis": lincs_mcf10a_report["state_axis"],
            "technical_wells_count_as_independent_biology": False,
            "collection_holdout_passed": lincs_mcf10a_report["collection_holdout_passed"],
            "may_count_as_independent_lab_holdout": False,
            "OOD_refusal_validated": False,
            "candidate_mechanism_families": role["candidate_mechanism_families"],
            "may_emit_numeric_parameter_range": False,
            "may_select_aleph_parameter": False,
            "reason": (
                "C1-to-C2 transfer and paired directions reproduce after biological-group "
                "aggregation, but both collections are OHSU and the TGFB+EGF dose conflict, "
                "semantic OOD, Aleph operators, and forward sensitivity remain unresolved"
            ),
        })
        calibration_requests.append({
            "parameter": "mechanism_family:EMT_spatial_dispersion_and_cell_cell_adhesion",
            "status": (
                "blocked_collection_validated_but_independent_lab_and_operator_missing"
            ),
            "numeric_range": None,
            "required_forward_observables": role["required_Aleph_operators"],
            "required_before_unblock": role["required_before_numeric_promotion"],
            "may_emit_sweep_axis": False,
        })
    if karacosta_emt_report is not None:
        direction = karacosta_emt_report["cross_lab_common_marker_direction"]
        role = karacosta_emt_report["pre_sweep_role"]
        design = karacosta_emt_report["design_audit"]
        safely_bounded = bool(
            karacosta_emt_report["aleph_authority"] == "none"
            and karacosta_emt_report["representation_training_eligible"]
            and karacosta_emt_report["directional_bridge_eligible"]
            and direction["dataset_and_lab_disjoint"]
            and direction["independent_lab_common_marker_direction_reproduced"]
            and not direction["numeric_scale_transfer_permitted"]
            and design["mandatory_split_groups"] == 1
            and design["single_cells_are_not_independent_biological_replicates"]
            and design["cell_state_classifier_on_these_features_is_label_circular_and_forbidden"]
            and not karacosta_emt_report["observable_evidence_eligible"]
            and not karacosta_emt_report["independent_lab_inferential_holdout"]
            and not karacosta_emt_report["uncertainty_holdout_passed"]
            and not karacosta_emt_report["OOD_refusal_validated"]
            and role["numeric_range"] is None
            and not role["may_emit_sweep_axis"]
            and not karacosta_emt_report["may_select_aleph_parameter"]
        )
        checks["karacosta_emt_direction_filter_is_non_authoritative"] = safely_bounded
        representation_training_audits.append({
            "dataset": karacosta_emt_report["source"]["dataset_id"],
            "status": (
                "admitted_independent_lab_common_marker_direction_filter_only"
                if safely_bounded else "refused"
            ),
            "state_axis": karacosta_emt_report["state_axis"],
            "raw_single_cells": design["raw_single_cell_events"],
            "biological_replicate_groups": design["mandatory_split_groups"],
            "single_cells_count_as_independent_biology": False,
            "label_circular_classifier_forbidden": True,
            "independent_lab_common_marker_direction_reproduced": direction[
                "independent_lab_common_marker_direction_reproduced"
            ],
            "may_count_as_inferential_external_holdout": False,
            "may_emit_numeric_parameter_range": False,
            "may_select_aleph_parameter": False,
            "future_sweep_role": role["newly_supported_action"],
            "reason": (
                "Stanford HCC827 CyTOF reproduces the UNSW A549 IF E-cadherin-down/"
                "Vimentin-up direction across dataset, lab, cell type, and modality, "
                "but the released CyTOF source exposes one biological experiment"
            ),
        })
        calibration_requests.append({
            "parameter": "mechanism_family:EMT_marker_direction_conditioning",
            "candidate_mechanism_families": role["candidate_mechanism_families"],
            "status": (
                "direction_filter_admitted_numeric_sweep_blocked_by_"
                "replication_OOD_operator_and_forward_sensitivity"
            ),
            "numeric_range": None,
            "newly_supported_action": role["newly_supported_action"],
            "required_forward_observables": role["required_Aleph_operators"],
            "required_before_unblock": role["required_before_numeric_promotion"],
            "may_emit_sweep_axis": False,
        })
    if geo_emt_rna_report is not None:
        direction = geo_emt_rna_report["replicated_direction"]
        cross_assay = geo_emt_rna_report["cross_assay_bridge"]
        ood = geo_emt_rna_report["semantic_OOD_challenge"]
        role = geo_emt_rna_report["pre_sweep_role"]
        design = geo_emt_rna_report["design_audit"]
        safely_bounded = bool(
            geo_emt_rna_report["aleph_authority"] == "none"
            and geo_emt_rna_report["representation_training_eligible"]
            and geo_emt_rna_report["directional_pre_sweep_filter_eligible"]
            and geo_emt_rna_report["independent_study_direction_holdout"]
            and direction["all_ten_study_gene_directions_expected"]
            and direction["all_ten_between_group_ranges_completely_separated"]
            and direction["independent_studies"] == 5
            and direction["independent_labs"] == 5
            and direction[
                "nominal_exact_one_sided_sign_test_p_for_five_retrospective_study_signatures"
            ]
            == 0.03125
            and direction["retrospective_direction_replication_passed"]
            and not direction["study_selection_was_outcome_blind"]
            and not direction["confirmatory_direction_threshold_0_05_passed"]
            and not geo_emt_rna_report["study_level_direction_inference_passed"]
            and design["samples"] == 72
            and design["biological_replicates_not_expression_features"]
            and design["cross_platform_absolute_scale_pooling_forbidden"]
            and cross_assay[
                "RNA_direction_agrees_with_independent_IF_and_CyTOF_protein_direction"
            ]
            and not cross_assay["RNA_protein_absolute_scale_transfer_permitted"]
            and ood["challenge_defined"]
            and not ood["OOD_refusal_model_validated"]
            and not geo_emt_rna_report["numeric_uncertainty_promotion_passed"]
            and not geo_emt_rna_report["OOD_refusal_validated"]
            and role["numeric_range"] is None
            and not role["may_emit_sweep_axis"]
            and not geo_emt_rna_report["may_select_aleph_parameter"]
        )
        checks["geo_emt_rna_direction_constraint_is_non_authoritative"] = safely_bounded
        representation_training_audits.append({
            "dataset": "geo-emt-rna-five-study-panel",
            "status": (
                "admitted_retrospective_replicated_study_direction_constraint_only"
                if safely_bounded else "refused"
            ),
            "studies": design["studies"],
            "independent_lab_groups": design["independent_lab_groups"],
            "biological_samples": design["samples"],
            "expression_features_count_as_independent_biology": False,
            "cross_platform_absolute_scale_transfer": False,
            "RNA_protein_absolute_scale_transfer": False,
            "study_level_sign_test_p": direction[
                "nominal_exact_one_sided_sign_test_p_for_five_retrospective_study_signatures"
            ],
            "study_selection_was_outcome_blind": False,
            "may_count_as_confirmatory_direction_holdout": False,
            "may_count_as_numeric_uncertainty_holdout": False,
            "may_emit_numeric_parameter_range": False,
            "may_select_aleph_parameter": False,
            "future_sweep_role": role["supported_action"],
            "reason": (
                "five independent GEO laboratories retrospectively reproduce CDH1-down/"
                "VIM-up with nominal p=0.03125, but outcome-informed study selection "
                "blocks confirmatory authority in addition to numeric uncertainty, OOD, "
                "and operator validation"
            ),
        })
        calibration_requests.append({
            "parameter": "mechanism_family:replicated_EMT_RNA_direction_conditioning",
            "status": (
                "retrospective_direction_constraint_admitted_numeric_sweep_blocked_by_"
                "preregistered_holdout_numeric_uncertainty_OOD_operator_and_sensitivity"
            ),
            "numeric_range": None,
            "newly_supported_action": role["supported_action"],
            "required_forward_observables": [
                "CDH1_RNA_or_declared_protein_surrogate_observation_operator",
                "VIM_RNA_or_declared_protein_surrogate_observation_operator",
            ],
            "current_operator_resolution": {
                name: operator_by_observable.get(name, {"status": "not_audited"})
                for name in ("CDH1_RNA_abundance", "VIM_RNA_abundance")
            },
            "required_before_unblock": role["remaining_blockers"],
            "may_emit_sweep_axis": False,
        })
    if gse325309_emt_report is not None:
        primary = gse325309_emt_report["frozen_primary_test"]
        design = gse325309_emt_report["design_audit"]
        role = gse325309_emt_report["pre_sweep_role"]
        safely_bounded = bool(
            gse325309_emt_report["aleph_authority"] == "none"
            and gse325309_emt_report["pristine_external_direction_holdout_passed"]
            and primary["passed"] and primary["exact_one_sided_p"] == 0.05
            and primary["exact_assignment_count"] == 20
            and primary["directions"] == {"CDH1": "decrease", "VIM": "increase"}
            and design["biological_samples"] == 6
            and design["sample_exclusions"] == 0
            and design["biological_replicate_is_test_unit"]
            and design["protocol_frozen_before_expression_access"]
            and not gse325309_emt_report["numeric_uncertainty_promotion_passed"]
            and not gse325309_emt_report["OOD_refusal_validated"]
            and role["numeric_range"] is None and not role["may_emit_sweep_axis"]
            and not gse325309_emt_report["may_select_aleph_parameter"]
        )
        checks["gse325309_preregistered_direction_holdout_is_non_authoritative"] = safely_bounded
        representation_training_audits.append({
            "dataset": "geo-gse325309-preregistered-a549-tgfb1-rnaseq",
            "status": "admitted_pristine_external_direction_holdout_only" if safely_bounded else "refused",
            "biological_samples": design["biological_samples"],
            "sample_exclusions": design["sample_exclusions"],
            "exact_one_sided_p": primary["exact_one_sided_p"],
            "pristine_confirmation": True,
            "may_count_as_confirmatory_direction_holdout": safely_bounded,
            "may_count_as_numeric_uncertainty_holdout": False,
            "may_emit_numeric_parameter_range": False,
            "may_select_aleph_parameter": False,
            "future_sweep_role": role["supported_action"],
        })
        calibration_requests.append({
            "parameter": "mechanism_family:preregistered_EMT_direction_confirmation",
            "status": "pristine_direction_holdout_passed_numeric_sweep_blocked_by_operator_OOD_and_scale",
            "numeric_range": None, "newly_supported_action": role["supported_action"],
            "required_forward_observables": [
                "CDH1_RNA_abundance", "VIM_RNA_abundance"
            ],
            "current_operator_resolution": {
                name: operator_by_observable.get(name, {"status": "not_audited"})
                for name in ("CDH1_RNA_abundance", "VIM_RNA_abundance")
            },
            "required_before_unblock": role["remaining_blockers"],
            "may_emit_sweep_axis": False,
        })
    if tendon_report is not None and tendon_report["observable_constraint_eligible"]:
        matrix = tendon_report["matrix_tension_contrasts"]
        hypotheses = (
            ("aleph.vertical.sf_arc.MaterialCard.active_tension_pn", "pN",
             "effective stress-fibre/arc active tension"),
            ("aleph.vertical.nmii.MotorKinetics.crossbridge_stiffness", "pN_per_um",
             "single-head motor force scale"),
            ("aleph.vertical.nmii.MotorKinetics.attachment_rate", "per_s",
             "motor bound-fraction kinetics"),
            ("aleph.vertical.focal_adhesion.ClutchCard.bond_stiffness_pn_per_um", "pN_per_um",
             "ligand-side load transmission"),
            ("aleph.vertical.focal_adhesion.ClutchCard.anchor_stiffness_pn_per_um", "pN_per_um",
             "actin-side load transmission"),
        )
        parameter_axis_candidates.extend({
            "parameter": parameter,
            "parameter_unit": unit,
            "mechanistic_role": role,
            "status": "blocked_pending_forward_sensitivity_calibration",
            "numeric_range": None,
            "evidence_relation": "higher matrix tension co-occurs with higher traction, alpha-SMA, and Acta2",
            "identifiability": "confounded_candidate_not_identified_by_observational_contrast",
            "calibration_contract": "aleph.outer_library.forward_sensitivity.v1",
            "may_emit_sweep_axis": False,
        } for parameter, unit, role in hypotheses)
        calibration_requests.append({
            "parameter": "aleph.vertical.sf_arc.MaterialCard.active_tension_pn",
            "candidate_parameters": [parameter for parameter, _, _ in hypotheses],
            "evaluator": "evaluate_sweep_calibration.py",
            "required_forward_observables": [
                "tissue_traction_force_per_cell",
                "nuclear_circularity",
            ],
            "external_target_constraints": {
                name: {
                    "mean_difference": matrix[name]["mean_difference"],
                    "bootstrap_ci95": matrix[name]["bootstrap_ci95"],
                }
                for name in (
                    "tissue_traction_force_per_cell",
                    "alpha_SMA_relative_fluorescence",
                    "Acta2_relative_expression",
                )
            },
            "required_before_unblock": [
                "one-at-a-time Aleph forward sensitivity is monotonic over a declared safe range",
                "traction measurement operator has matched units and geometry",
            ] + ([] if (
                fibrotic_report is not None
                and fibrotic_report["independent_lab_mechanism_holdout_passed"]
            ) else ["at least one independent-lab traction/state holdout passes"]),
        })
    if (mechano_osmotic_report is not None
            and mechano_osmotic_report["observable_constraint_eligible"]):
        parameter_axis_candidates.append({
            "parameter": "requested_nonexistent_axis:membrane_cortex_attachment_energy_pn_per_um",
            "parameter_unit": "pN_per_um",
            "mechanistic_role": "membrane-cortex attachment energy W in sigma_app = sigma_bilayer + W",
            "status": "blocked_missing_physics_term_and_exact_nonidentifiability",
            "numeric_range": None,
            "evidence_relation": (
                "Y-27632 transiently increases AFM tether force during spreading, "
                "but has no experiment-level steady-state tether-force effect"
            ),
            "identifiability": (
                "one tether force identifies only sigma_bilayer + W; no current "
                "Aleph axis represents W and neither term is separately identified"
            ),
            "calibration_contract": "aleph.outer_library.forward_sensitivity.v1",
            "may_emit_sweep_axis": False,
        })
        calibration_requests.append({
            "parameter": "requested_nonexistent_axis:membrane_cortex_attachment_energy_pn_per_um",
            "candidate_parameters": [
                "requested_nonexistent_axis:membrane_cortex_attachment_energy_pn_per_um",
                "aleph.vertical.membrane.MembraneMaterial.tension_pn_per_um",
                "aleph.vertical.membrane.MembraneMaterial.bending_rigidity_pn_um",
            ],
            "evaluator": "evaluate_sweep_calibration.py",
            "status": "blocked_missing_W_axis_and_second_identifying_observable",
            "required_forward_observables": [
                "membrane_tether_plateau_force:pN",
                "initial_spreading_area_rate_square_micrometre_per_minute",
                "initial_volume_flux_cubic_micrometre_per_minute",
            ],
            "external_target_constraints": {
                "AFM_membrane_tether_force_spreading_phase_Y27632_minus_control_pN": (
                    mechano_osmotic_report["AFM_membrane_tether_force"][
                        "spreading_phase_30_to_90_min"
                    ]["paired_experiment_bootstrap_95ci_pN"]
                ),
                "initial_spreading_area_rate_Y27632_minus_control": (
                    mechano_osmotic_report["initial_spreading_area_rate"][
                        "independent_experiment_bootstrap_95ci"
                    ]
                ),
                "initial_volume_flux_Y27632_minus_control": (
                    mechano_osmotic_report["initial_volume_flux"][
                        "independent_experiment_bootstrap_95ci"
                    ]
                ),
            },
            "required_before_unblock": [
                "PI decides whether membrane-cortex attachment energy W enters Aleph or remains out of scope",
                "add an independent observable that separates sigma_bilayer from W",
                "reproduce the transient-versus-steady-state contrast in matched geometry",
                "show that actomyosin, ion-flux, and membrane-reservoir effects are not collapsed into one tension parameter",
                "pass an independent-lab confirmatory experiment",
            ],
        })
    # All optional evidence contracts must be included before the final gate
    # decision. Computing this near the initial core checks would let a later
    # failed non-authority audit appear in failed_checks without changing the
    # overall status.
    passed = all(checks.values())
    return {
        "status": "eligible_for_non_authoritative_proposal" if passed else "refused",
        "checks": checks,
        "failed_checks": [name for name, value in checks.items() if not value],
        "may_mutate_physics": False,
        "may_write_decision_authority": False,
        "observable_evidence_channel": observable_channel,
        "evidence_exclusions": evidence_exclusions,
        "context_heterogeneity_audit": rock_context_report,
        "mechano_osmotic_audit": mechano_osmotic_report,
        "piv_contractility_audit": piv_contractility_report,
        "nematic_tfm_audit": nematic_tfm_report,
        "opto_force_audit": opto_force_report,
        "tfm_fret_audit": tfm_fret_report,
        "exif_emt_audit": exif_emt_report,
        "mtrack_emt_audit": mtrack_emt_report,
        "geo_emt_rna_audit": geo_emt_rna_report,
        "gse325309_emt_audit": gse325309_emt_report,
        "observation_operator_audit": observation_operator_report,
        "parameter_axis_candidates": parameter_axis_candidates,
        "calibration_requests": calibration_requests,
        "representation_training_audits": representation_training_audits,
        "multitask_model_audit": multitask_model_audit,
        "uncertainty_audit": uncertainty_report,
        "next_action": (
            "emit bounded candidate sweep for human/physics review" if passed
            else (
                "run unit-matched Aleph forward sensitivity calibration; separately improve the per-cell classifier"
                if fibrotic_report is not None
                and fibrotic_report["independent_lab_mechanism_holdout_passed"]
                else "acquire aligned independent-lab data and rerun grouped holdout"
            )
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mechanism-report", type=Path)
    parser.add_argument("--tendon-report", type=Path)
    parser.add_argument("--fibrotic-report", type=Path)
    parser.add_argument("--uncertainty-report", type=Path)
    parser.add_argument("--migration-report", type=Path)
    parser.add_argument("--rock-context-report", type=Path)
    parser.add_argument("--mechano-osmotic-report", type=Path)
    parser.add_argument("--cell-monolayer-report", type=Path)
    parser.add_argument("--multitask-model-report", type=Path)
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
    parser.add_argument("--observation-operator-report", type=Path)
    args = parser.parse_args()
    mechanism = (
        json.loads(args.mechanism_report.read_text(encoding="utf-8"))
        if args.mechanism_report else None
    )
    result = evaluate(
        json.loads(args.model_report.read_text(encoding="utf-8")),
        mechanism,
        (json.loads(args.tendon_report.read_text(encoding="utf-8"))
         if args.tendon_report else None),
        (json.loads(args.fibrotic_report.read_text(encoding="utf-8"))
         if args.fibrotic_report else None),
        (json.loads(args.uncertainty_report.read_text(encoding="utf-8"))
         if args.uncertainty_report else None),
        (json.loads(args.migration_report.read_text(encoding="utf-8"))
         if args.migration_report else None),
        (json.loads(args.rock_context_report.read_text(encoding="utf-8"))
         if args.rock_context_report else None),
        (json.loads(args.mechano_osmotic_report.read_text(encoding="utf-8"))
         if args.mechano_osmotic_report else None),
        (json.loads(args.cell_monolayer_report.read_text(encoding="utf-8"))
         if args.cell_monolayer_report else None),
        (json.loads(args.multitask_model_report.read_text(encoding="utf-8"))
         if args.multitask_model_report else None),
        (json.loads(args.piv_contractility_report.read_text(encoding="utf-8"))
         if args.piv_contractility_report else None),
        (json.loads(args.nematic_tfm_report.read_text(encoding="utf-8"))
         if args.nematic_tfm_report else None),
        (json.loads(args.opto_force_report.read_text(encoding="utf-8"))
         if args.opto_force_report else None),
        (json.loads(args.tfm_fret_report.read_text(encoding="utf-8"))
         if args.tfm_fret_report else None),
        (json.loads(args.exif_emt_report.read_text(encoding="utf-8"))
         if args.exif_emt_report else None),
        (json.loads(args.mtrack_emt_report.read_text(encoding="utf-8"))
         if args.mtrack_emt_report else None),
        (json.loads(args.lincs_mcf10a_report.read_text(encoding="utf-8"))
         if args.lincs_mcf10a_report else None),
        (json.loads(args.karacosta_emt_report.read_text(encoding="utf-8"))
         if args.karacosta_emt_report else None),
        (json.loads(args.geo_emt_rna_report.read_text(encoding="utf-8"))
         if args.geo_emt_rna_report else None),
        (json.loads(args.gse325309_emt_report.read_text(encoding="utf-8"))
         if args.gse325309_emt_report else None),
        (json.loads(args.observation_operator_report.read_text(encoding="utf-8"))
         if args.observation_operator_report else None),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
