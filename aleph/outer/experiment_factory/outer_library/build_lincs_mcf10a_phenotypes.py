#!/usr/bin/env python3
"""Build replicate-safe LINCS MCF10A IF phenotype observations."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mcf10a_xlsx import read_first_sheet
from schema import Observation


DATASET_ID = "lincs-mcf10a-ligand-if-phenotypes"
LAB_GROUP = "OHSU_MEP_LINCS_MCF10A"
ACCESSION = "10.1038/s42003-022-03975-9"
PHENOTYPE_FILE = "42003_2022_3975_MOESM3_ESM.xlsx"
METADATA_FILE = "42003_2022_3975_MOESM25_ESM.xlsx"
OBSERVABLES = {
    "Well_Cell_Count": ("well_cell_count", "cells_per_well"),
    "Intensity_IntegratedIntensity_EdU": (
        "EdU_integrated_intensity", "author_fluorescence_integrated_intensity_arbitrary_unit"
    ),
    "Intensity_MeanIntensity_EdU": (
        "EdU_mean_intensity", "author_fluorescence_mean_intensity_arbitrary_unit"
    ),
    "Edu_Positive_Proportion": ("EdU_positive_proportion", "fraction"),
    "G2m_Proportion": ("G2M_proportion", "fraction"),
    "Intensity_IntegratedIntensity_KRT5": (
        "KRT5_integrated_intensity", "author_fluorescence_integrated_intensity_arbitrary_unit"
    ),
    "Intensity_MeanIntensity_KRT5": (
        "KRT5_mean_cytoplasmic_intensity", "author_fluorescence_mean_intensity_arbitrary_unit"
    ),
    "Normalized_Second_Neighbor_Dist": (
        "normalized_second_neighbor_distance", "dimensionless_ratio_to_random_expectation"
    ),
    "Mean_Cells_per_Cluster": ("mean_cells_per_cluster", "cells_per_cluster"),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(source_dir: Path, receipt_path: Path, manifest_path: Path,
          tensor_path: Path) -> dict[str, object]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("schema") != "aleph.outer_library.lincs_mcf10a_acquisition.v1":
        raise ValueError("LINCS MCF10A acquisition schema changed")
    records = {item["local_name"]: item for item in receipt["files"]}
    phenotype_path, metadata_path = source_dir / PHENOTYPE_FILE, source_dir / METADATA_FILE
    for path in (phenotype_path, metadata_path):
        record = records[path.name]
        if not path.is_file() or path.stat().st_size != record["size"] or _sha256(path) != record["sha256"]:
            raise ValueError(f"LINCS MCF10A source changed: {path.name}")
    headers, phenotype_rows = read_first_sheet(phenotype_path)
    metadata_headers, metadata_rows = read_first_sheet(metadata_path)
    if headers != ["specimenName", "ligand", "experimentalTimePoint", "collection", "replicate", "WellIndex", *OBSERVABLES]:
        raise ValueError("LINCS phenotype headers changed")
    if metadata_headers[:19] != [
        "specimenID", "specimenName", "experimentalCondition", "ligand",
        "experimentalTimePoint", "replicate", "collection", "ligandDose",
        "timePointUnit", "species", "sex", "cellLine", "isCellLine",
        "cellType", "ligandDoseUnit", "secondLigand", "secondLigandDose",
        "secondLigandDoseUnit", "studySite",
    ]:
        raise ValueError("LINCS metadata headers changed")
    metadata_by_name = {str(row["specimenName"]): row for row in metadata_rows}
    if len(metadata_by_name) != len(metadata_rows):
        raise ValueError("LINCS specimen metadata names are not unique")
    if any(str(row["specimenName"]) not in metadata_by_name for row in phenotype_rows):
        raise ValueError("LINCS phenotype specimen misses metadata")
    observations, tensor_rows, row_ids, split_groups = [], [], [], []
    phenotype_sha, metadata_sha = _sha256(phenotype_path), _sha256(metadata_path)
    for row in phenotype_rows:
        name = str(row["specimenName"])
        meta = metadata_by_name[name]
        if any(str(row[key]) != str(meta[key]) for key in ("ligand", "experimentalTimePoint", "collection", "replicate")):
            raise ValueError(f"LINCS phenotype/metadata mismatch for {name}")
        source_row = int(row["_source_row"])
        technical_id = f"{name}:xlsx_row_{source_row}:WellIndex_{row['WellIndex']}"
        row_ids.append(technical_id)
        split_group = f"{row['collection']}:{row['replicate']}"
        split_groups.append(split_group)
        values = []
        for source_name, (observable, unit) in OBSERVABLES.items():
            raw = row[source_name]
            value = np.nan if raw in (None, "NA") else float(raw)
            values.append(value)
            if np.isnan(value):
                continue
            observation = Observation(
                observation_id=f"lincs-mcf10a-row{source_row}-{observable}",
                dataset_id=DATASET_ID, lab_group=LAB_GROUP,
                provider="MEP-LINCS/OHSU via Communications Biology supplementary data",
                accession=ACCESSION, license_id="cc-by-4.0",
                source_sha256=phenotype_sha,
                source_locator=f"{PHENOTYPE_FILE}:Sup Data 1 Phenotypes!row{source_row}:{source_name}",
                sample_id=technical_id,
                biological_replicate=split_group,
                technical_replicate=f"author_well_row_{source_row}_WellIndex_{row['WellIndex']}",
                condition=(
                    f"{meta['ligand']}_{meta['ligandDose']}{meta['ligandDoseUnit']}"
                    f"_plus_{meta['secondLigand']}_{meta['secondLigandDose']}{meta['secondLigandDoseUnit']}"
                    f"_{meta['experimentalTimePoint']}{meta['timePointUnit']}"
                ),
                cell_type="MCF10A_human_female_normal_breast_epithelium",
                cell_state="author_ligand_treatment_condition_not_per_cell_state_ground_truth",
                modality="fixed_IF_author_well_level_population_summary",
                observable=observable, value=value, unit=unit,
                unit_source_sha256=metadata_sha,
                unit_source_locator=f"{METADATA_FILE}:MDD_sample_annotations:{name}",
                coordinate_name="treatment_time",
                coordinate_value=float(meta["experimentalTimePoint"]),
                coordinate_unit=str(meta["timePointUnit"]),
            )
            observation.validate()
            observations.append(observation)
        tensor_rows.append(values)
    X = np.asarray(tensor_rows, dtype=np.float64)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text("".join(json.dumps(row.as_dict(), sort_keys=True) + "\n" for row in observations))
    np.savez_compressed(
        tensor_path, X=X, feature_names=np.asarray(list(OBSERVABLES), dtype="U64"),
        specimen_name=np.asarray([str(row["specimenName"]) for row in phenotype_rows], dtype="U32"),
        ligand=np.asarray([str(row["ligand"]) for row in phenotype_rows], dtype="U8"),
        time_h=np.asarray([int(row["experimentalTimePoint"]) for row in phenotype_rows], dtype=np.int16),
        collection=np.asarray([str(row["collection"]) for row in phenotype_rows], dtype="U2"),
        replicate=np.asarray([str(row["replicate"]) for row in phenotype_rows], dtype="U1"),
        well_index=np.asarray([int(row["WellIndex"]) for row in phenotype_rows], dtype=np.int16),
        source_row=np.asarray([int(row["_source_row"]) for row in phenotype_rows], dtype=np.int16),
        sample_id=np.asarray(row_ids, dtype="U64"), split_group=np.asarray(split_groups, dtype="U8"),
    )
    return {
        "schema": "aleph.outer_library.lincs_mcf10a_if_phenotypes.v1",
        "source": {"dataset_id": DATASET_ID, "lab_group": LAB_GROUP, "paper_doi": ACCESSION},
        "counts": {
            "author_well_rows": len(phenotype_rows), "canonical_observations": len(observations),
            "biological_replicate_groups": len(set(split_groups)), "collections": len(set(row["collection"] for row in phenotype_rows)),
        },
        "tensor_shape": list(X.shape), "feature_names": list(OBSERVABLES),
        "replicate_contract": {
            "split_unit": "collection_plus_replicate_biological_group",
            "technical_unit": "source_row_because_WellIndex_repeats_across_unreleased_plates",
            "rows_may_be_split_independently": False,
            "paper_reports_eight_biological_replicates_but_released_IF_table_contains_seven_groups": True,
        },
        "metadata_conflict": {
            "paper_methods_TGFB_companion_EGF_ng_per_mL": 10.0,
            "supplementary_metadata_TGFB_companion_EGF_ng_per_mL": 20.0,
            "resolved": False,
        },
        "manifest_sha256": _sha256(manifest_path), "tensor_sha256": _sha256(tensor_path),
        "aleph_authority": "none", "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build(args.source_dir, args.receipt, args.manifest, args.tensor)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
