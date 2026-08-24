#!/usr/bin/env python3
"""Build compound-disjoint well profiles from the Figshare Cell Painting corpus."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np


SOURCE_FILE = "aggregate_profiles_DP_aggregated.parquet"
SOURCE_SHA256 = "01d10a24f2fa448cb560425ce41712ff264ce701fab82a5ae5aaba273f677844"


def build(
    aggregate_dir: Path,
    extraction_receipt: Path,
    design_audit_path: Path,
    output: Path,
    report_path: Path,
) -> dict[str, object]:
    try:
        import pyarrow.parquet as parquet
    except ImportError as error:
        raise RuntimeError("PyArrow is required only for the one-time Parquet-to-NPZ conversion") from error
    receipt = json.loads(extraction_receipt.read_text(encoding="utf-8"))
    source_record = receipt.get("aggregate_files", {}).get(SOURCE_FILE, {})
    if source_record.get("sha256") != SOURCE_SHA256 or receipt.get("archive_paths_safe") is not True:
        raise ValueError("Figshare aggregate extraction receipt is not verified")
    audit = json.loads(design_audit_path.read_text(encoding="utf-8"))
    if audit.get("cross_split_compound_overlap_count") != 0:
        raise ValueError("Figshare design audit contains compound leakage")
    assignment = {
        compound: (label, role)
        for label, compounds in audit["compound_assignments"].items()
        for compound, role in compounds.items()
    }
    path = aggregate_dir / SOURCE_FILE
    table = parquet.read_table(path)
    metadata = [
        "moa", "Metadata_Plate", "Metadata_Well", "Metadata_Site",
        "Metadata_cmpdConc", "label", "Metadata_cmpdName",
    ]
    features = [name for name in table.column_names if name.startswith("Feature_")]
    if len(features) != 672 or table.num_rows != 6_626:
        raise ValueError("Figshare DeepProfiler aggregate schema changed")
    meta = {name: table[name].to_pylist() for name in metadata}
    values = np.column_stack([
        table[name].to_numpy(zero_copy_only=False) for name in features
    ]).astype(np.float32)
    if not np.isfinite(values).all():
        raise ValueError("Figshare aggregate profiles contain non-finite values")
    groups: dict[tuple[str, str], list[int]] = defaultdict(list)
    for index, key in enumerate(zip(meta["Metadata_Plate"], meta["Metadata_Well"], strict=True)):
        groups[key].append(index)

    X, plates, wells, compounds, labels, roles, concentrations, site_counts = [], [], [], [], [], [], [], []
    for (plate, well), indices in sorted(groups.items()):
        compound = str(meta["Metadata_cmpdName"][indices[0]])
        label = str(meta["moa"][indices[0]])
        concentration = float(meta["Metadata_cmpdConc"][indices[0]])
        if any(
            meta[field][index] != meta[field][indices[0]]
            for field in ("moa", "Metadata_cmpdName", "Metadata_cmpdConc")
            for index in indices
        ):
            raise ValueError(f"metadata changes within well {plate}/{well}")
        if compound not in assignment or assignment[compound][0] != label:
            raise ValueError(f"compound/label absent from frozen audit: {compound}/{label}")
        X.append(values[indices].mean(axis=0))
        plates.append(plate); wells.append(well); compounds.append(compound); labels.append(label)
        roles.append(assignment[compound][1]); concentrations.append(concentration); site_counts.append(len(indices))
    X_array = np.asarray(X, dtype=np.float32)
    if X_array.shape != (954, 672) or set(compounds) != set(assignment):
        raise ValueError(f"Figshare well tensor contract changed: {X_array.shape}")
    role_compounds = {role: set(np.asarray(compounds)[np.asarray(roles) == role]) for role in set(roles)}
    if sum(bool(role_compounds[left] & role_compounds[right]) for left in role_compounds for right in role_compounds if left < right):
        raise ValueError("compound leaked across generated roles")
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output, X=X_array, feature_name=np.asarray(features), plate=np.asarray(plates),
        well=np.asarray(wells), compound=np.asarray(compounds), moa=np.asarray(labels),
        compound_split=np.asarray(roles), concentration=np.asarray(concentrations, dtype=np.float32),
        site_count=np.asarray(site_counts, dtype=np.int8),
    )
    report = {
        "schema": "aleph.outer_library.figshare_cell_death_tensor.v1",
        "source_doi": "10.17044/scilifelab.28202864.v2",
        "source_aggregate_sha256": SOURCE_SHA256,
        "source_profile_family": "DeepProfiler_aggregate_profiles",
        "site_rows": int(table.num_rows), "well_rows": int(X_array.shape[0]),
        "feature_count": int(X_array.shape[1]), "compound_count": len(set(compounds)),
        "plate_count": len(set(plates)),
        "split_compound_counts": {role: len(values) for role, values in sorted(role_compounds.items())},
        "split_well_counts": dict(sorted(Counter(roles).items())),
        "label_well_counts": dict(sorted(Counter(labels).items())),
        "site_count_range": [min(site_counts), max(site_counts)],
        "split_unit": "compound_identity",
        "profile_unit": "plate_well_mean_of_available_sites",
        "cross_split_compound_overlap_count": 0,
        "single_cell_profiles_extracted": False,
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "eligible_role": "same_provider_unseen_compound_programmed_cell_death_training",
        "same_endpoint_as_biad2515_annexin_v": False,
        "aleph_authority": "none", "may_select_aleph_parameter": False,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aggregate-dir", type=Path, required=True)
    parser.add_argument("--extraction-receipt", type=Path, required=True)
    parser.add_argument("--design-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = build(args.aggregate_dir, args.extraction_receipt, args.design_audit, args.output, args.report)
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
