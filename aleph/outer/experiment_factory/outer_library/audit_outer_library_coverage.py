#!/usr/bin/env python3
"""Adversarially audit Outer Library coverage without inflating row counts.

This report treats time points and image-derived observables as observations,
not independent biological samples.  It also distinguishes a perturbation
condition from a biological cell-state ontology label.
"""

from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path
from typing import Any


SCHEMA = "aleph.outer_library.adversarial_coverage_audit.v1"
UNKNOWN_UNIT_TOKENS = (
    "unknown", "unspecified", "author_processed", "author_scaled", "author_normalized"
)
UNKNOWN_REPLICATE_TOKENS = (
    "unknown", "not_reported", "unreported", "not_released", "not_available",
)

# Conservative manual ontology: a dataset counts only when its author labels or
# assay explicitly instantiate the axis. Row volume never upgrades lab/provider
# independence.
STATE_AXES: dict[str, dict[str, Any]] = {
    "mechanosensitive_calcium": {
        "datasets": [
            "zenodo-18495219-hela-mechanotransduction",
            "zenodo-19890985-eahy926-piezo1",
            "zenodo-8215150-mdck-mechanocalcium",
        ],
        "status": "covered_but_external_classifier_failed",
    },
    "migration_and_ROCK_cytoskeletal_state": {
        "datasets": [
            "scilifelab-21407402-multisite-live-cell-2d",
            "dryad-9jh6m-ROCK-modulation",
            "elife-71032-nih3t3-3d-cdm-migration",
        ],
        "status": "covered_context_dependent",
    },
    "fibrotic_myofibroblast_state": {
        "datasets": [
            "nature-2026-tendon-mechanoculture",
            "geo-gse226374-tgfb-celastrol-rnaseq",
        ],
        "status": "covered_condition_level_two_labs",
    },
    "mechano_osmotic_membrane_state": {
        "datasets": ["elife-72381-hela-mechano-osmotic"],
        "status": "single_lab_operator_and_identifiability_blocked",
    },
    "supracellular_contractility_and_ECM_alignment": {
        "datasets": [
            "dryad-08kprr59c-c2c12-supracontractility",
            "dataverse-data2772-nih3t3-cellular-nematics",
            "dryad-sj3tx9683-opto-mdck-force-propagation",
        ],
        "status": "three_lab_cross_cell_type_direction_operator_and_biological_replication_blocked",
    },
    "focal_adhesion_load_transfer": {
        "datasets": [
            "zenodo-7432971-tfm-tether",
            "zenodo-14692589-tfm-fret",
        ],
        "status": "two_lab_same_unit_stiffness_direction_and_TFM_FRET_identifiability_replication_blocked",
    },
    "subcellular_localization_state": {
        "datasets": [
            "hpa-v25.1-subcellular-image-embeddings",
            "bbbc013-v1-fkhr-translocation",
            "bbbc014-v1-nfkb-translocation",
        ],
        "status": "training_signal_present_no_independent_provider_head",
    },
    "cell_cycle_and_proliferation": {
        "datasets": ["bbbc048-v1-jurkat-cell-cycle"],
        "status": "author_labelled_single_provider_training_only",
    },
    "apoptosis_and_cell_death": {
        "datasets": ["biad2515-hela-staurosporine-apoptosis"],
        "status": "author_AnnexinV_endpoint_single_provider_well_holdout_training_only",
    },
    "senescence": {
        "datasets": [
            "geo-gse63577-replicative-senescence",
            "geo-gse297406-replicative-senescence",
            "geo-gse282425-doxorubicin-senescence",
            "geo-gse301164-postfreeze-senescence-confirmation",
        ],
        "status": "directionally_confirmed_across_labs_absolute_state_OOD_blocked",
    },
    "EMT_epithelial_mesenchymal_state": {
        "datasets": [
            "figshare-26500210-exif-a549-emt-4plex",
            "github-mtrack-a549-tgfb4ng-xy01",
            "lincs-mcf10a-ligand-if-phenotypes",
            "nature-2019-karacosta-hcc827-emt-cytof",
        ],
        "status": "ten_lab_sources_with_preregistered_GSE325309_direction_holdout_passed_but_no_numeric_OOD_or_Aleph_operator",
    },
    "immune_activation_and_inflammation": {
        "datasets": ["bbbc054-v1-LPS-microglia-timecourse"],
        "status": "single_annotated_replicate_training_only",
    },
    "mitochondrial_metabolic_stress": {
        "datasets": ["bbbc053-v1-CAD-FCCP-mitochondria"],
        "status": "FCCP_single_provider_training_only_no_biological_replicate_ids",
    },
    "hypoxia_state": {"datasets": [], "status": "missing"},
    "differentiation_and_stemness": {"datasets": [], "status": "missing"},
    "DNA_damage_and_genotoxic_stress": {"datasets": [], "status": "missing"},
    "infection_response": {"datasets": [], "status": "missing"},
}

MODALITY_GAPS = {
    "western_blot": "missing",
    "single_cell_RNA_seq": "missing",
    "spatial_transcriptomics": "missing",
    "proteomics_phosphoproteomics": "missing",
    "metabolomics": "missing",
    "flow_cytometry_or_CyTOF": "imaging_flow_cytometry_plus_90066_cell_HCC827_CyTOF_one_released_biological_experiment",
    "electron_microscopy": "missing",
    "raw_multichannel_IF_for_trainable_encoder": "HPA_pilot_plus_independent_provider_ExIF_raw_five_channel_pilot_and_MCF10A_author_well_phenotypes_but_no_raw_common_task_holdout",
    "live_fluorescence_state_trajectory": "paper_linked_MTRACK_A549_VIM_RFP_single_treated_trajectory_no_control_or_biological_replicates",
    "traction_force_microscopy": "exact_Pa_two_lab_stiffness_direction_plus_mN_per_m_NIH3T3_source_added_but_replication_mapping_and_Aleph_operator_blocked",
    "PIV_collective_velocity": "C2C12_replicates_plus_NIH3T3_sequence_present_aligned_external_same_task_holdout_still_missing",
    "RT_qPCR": "present_condition_level",
    "bulk_RNA_seq": "present_small_marker_panel",
    "live_multichannel_apoptosis_imaging": "author_processed_profiles_present_raw_archive_indexed_not_local",
}


def audit(results: Path) -> dict[str, Any]:
    manifests = sorted(results.glob("*_observations.jsonl"))
    datasets: set[str] = set()
    labs: set[str] = set()
    samples: set[tuple[str, str]] = set()
    observations = 0
    unit_unknown = 0
    biological_replicate_unknown = 0
    per_modality: dict[str, dict[str, set[Any] | int]] = {}
    state_labels: collections.Counter[str] = collections.Counter()
    for manifest in manifests:
        with manifest.open(encoding="utf-8") as handle:
            for line in handle:
                row = json.loads(line)
                observations += 1
                dataset = row["dataset_id"]
                lab = row["lab_group"]
                sample = (dataset, row["sample_id"])
                modality = row["modality"]
                datasets.add(dataset)
                labs.add(lab)
                samples.add(sample)
                state_labels[row["cell_state"]] += 1
                unit_unknown += any(token in row["unit"].lower() for token in UNKNOWN_UNIT_TOKENS)
                replicate = row["biological_replicate"].strip().lower()
                biological_replicate_unknown += (
                    not replicate
                    or any(token in replicate for token in UNKNOWN_REPLICATE_TOKENS)
                )
                info = per_modality.setdefault(
                    modality,
                    {"observations": 0, "datasets": set(), "labs": set(), "samples": set()},
                )
                info["observations"] = int(info["observations"]) + 1
                info["datasets"].add(dataset)  # type: ignore[union-attr]
                info["labs"].add(lab)  # type: ignore[union-attr]
                info["samples"].add(sample)  # type: ignore[union-attr]

    hpa_receipt_path = results / "hpa_if_embedding_receipt.json"
    hpa_receipt = json.loads(hpa_receipt_path.read_text()) if hpa_receipt_path.exists() else None
    if hpa_receipt:
        datasets.add(hpa_receipt["dataset_id"])
        labs.add(hpa_receipt["lab_group"])

    modality_rows: dict[str, Any] = {}
    for name, info in sorted(per_modality.items()):
        modality_rows[name] = {
            "observations": info["observations"],
            "sample_count": len(info["samples"]),  # type: ignore[arg-type]
            "dataset_count": len(info["datasets"]),  # type: ignore[arg-type]
            "lab_group_count": len(info["labs"]),  # type: ignore[arg-type]
        }

    missing_axes = [name for name, item in STATE_AXES.items() if item["status"] == "missing"]
    return {
        "schema": SCHEMA,
        "scope": "canonical_manifests_plus_HPA_embedding_receipt",
        "counts": {
            "canonical_observations": observations,
            "canonical_samples": len(samples),
            "datasets_including_HPA_embedding_source": len(datasets),
            "conservative_lab_groups_including_HPA": len(labs),
            "observations_per_sample": observations / len(samples),
            "HPA_embedding_images_separate_from_physical_observations": (
                hpa_receipt["image_count"] if hpa_receipt else 0
            ),
        },
        "provenance_defects": {
            "unknown_or_author_processed_unit_observations": unit_unknown,
            "unknown_unit_fraction": unit_unknown / observations,
            "biological_replicate_unreported_observations": biological_replicate_unknown,
            "biological_replicate_unreported_fraction": biological_replicate_unknown / observations,
        },
        "per_modality": modality_rows,
        "state_axis_coverage": STATE_AXES,
        "missing_state_axis_count": len(missing_axes),
        "missing_state_axes": missing_axes,
        "modality_gap_status": MODALITY_GAPS,
        "adversarial_findings": [
            f"{observations:,} observation rows are not {observations:,} independent samples; time points and derived observables inflate row count",
            "cell_state currently mixes perturbation, phenotype, motion mode, and assay context; it is not a controlled ontology",
            "HPA contributes 81k images but only one provider/annotation lineage and cannot validate universal state transfer",
            "strong HPA macro AUROC coexists with modest macro F1, exposing rare-label and threshold failure",
            "cross-species force datasets cannot be pooled as one traction-state target without a declared measurement bridge",
            "missing units and biological replicate IDs block several otherwise large imaging datasets from inferential authority",
            "no paired multi-omic samples justify a shared latent space; task-scoped late fusion remains mandatory",
            "ExIF EMT labels are author treatment conditions, not per-cell EMT ground truth, and one plate date cannot support uncertainty or OOD authority",
            "M-TRACK contributes ten correlated frames from one treated trajectory; its late VIM-RFP intensity conflicts with the ExIF fixed-IF treatment direction and therefore blocks intensity-only sweep conditioning",
            "LINCS MCF10A C1-to-C2 transfer is a same-site collection holdout, not an independent-lab holdout; its TGFB condition also carries an unresolved companion-EGF dose conflict",
            "Karacosta HCC827 CyTOF reproduces E-cadherin-down/Vimentin-up across an independent lab, cell line, and modality, but 90,066 cells belong to one released biological experiment and CCAST labels are circular with the six clustering markers",
            "Five biologically replicated GEO studies retrospectively reproduce CDH1-down/VIM-up at nominal study-level one-sided p=0.03125, but outcome-informed literature selection blocks confirmatory authority; RNA magnitude is not a protein or Aleph parameter scale and one OOD context cannot establish general refusal",
            "Protocol-frozen GSE325309 passes its exact 20-assignment CDH1/VIM composite at p=0.05 after zero sample exclusions, establishing a pristine direction holdout but not numeric cross-assay calibration or OOD authority",
        ],
        "priority_queue": [
            {
                "priority": 1,
                "target": "common-task independent-provider raw IF localization holdout",
                "reason": "ExIF adds independent-provider raw IF state training, but its EMT treatment task cannot validate the HPA localization head",
            },
            {
                "priority": 2,
                "target": "independent cell-cycle/apoptosis panels plus replicated protein EMT, culture-context OOD, and controlled trajectories",
                "reason": "GSE325309 now supplies a preregistered outcome-blind CDH1/VIM direction holdout, but protein-scale uncertainty, multi-context OOD, and Aleph observation operators remain unresolved; M-TRACK is still one uncontrolled trajectory",
            },
            {
                "priority": 3,
                "target": "single-cell transcriptomic state source with donor/batch hierarchy",
                "reason": "adds molecular state coverage and independent donor holdouts absent from bulk marker panels",
            },
            {
                "priority": 4,
                "target": "aligned independent-lab TFM/PIV with biological replicate IDs and exact units",
                "reason": "a second lab now reproduces the same-unit stiffness-to-mean-traction direction and adds paired TFM-FRET, but released cell-to-replicate mapping, aligned intervention validation, and direct Aleph force/flow operators remain missing",
            },
            {
                "priority": 5,
                "target": "WB/proteomics mechanism confirmation",
                "reason": "orthogonal protein-level validation is currently absent",
            },
        ],
        "authority": "none",
        "production_eligible": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.results)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"missing_state_axes": report["missing_state_axes"], "provenance_defects": report["provenance_defects"]}, sort_keys=True))


if __name__ == "__main__":
    main()
