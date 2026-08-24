#!/usr/bin/env python3
"""Create a chemical-identity-disjoint split for programmed cell death.

The authors' published splits are well-disjoint, which is appropriate for
within-compound prediction, but most compounds occur in multiple roles.  A
mechanism generalization claim therefore needs the stricter split emitted here.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
import tarfile

from download_figshare_cell_death import PROGRAMMED_DEATH_LABELS, _digest


QC_SHA256 = "fc44402f50fe6d79be6ae28effe6c99a9e98c50719c05ef9fad1ace928a752b1"
SPLITS_SHA256 = "91adadc4cdebfb9c50ee917ccf4d26e972d5f43b7f495d54c0fe84e7a78f7f76"
ROLE_ORDER = ("train", "selection", "calibration", "sealed_test")


def _rank(label: str, compound: str) -> bytes:
    return hashlib.sha256(
        f"aleph-figshare-cell-death-compound-split-v1:{label}:{compound}".encode()
    ).digest()


def _assign(label: str, compounds: set[str]) -> dict[str, str]:
    ordered = sorted(compounds, key=lambda compound: _rank(label, compound))
    if len(ordered) < 2:
        raise ValueError(f"{label} has fewer than two distinct compounds")
    assignment = {compound: "train" for compound in ordered}
    assignment[ordered[-1]] = "sealed_test"
    if len(ordered) >= 3:
        assignment[ordered[-2]] = "calibration"
    if len(ordered) >= 4:
        assignment[ordered[-3]] = "selection"
    return assignment


def audit(source_dir: Path) -> dict[str, object]:
    qc_path = source_dir / "qc_df.csv"
    split_path = source_dir / "splits.tar.gz"
    if _digest(qc_path, "sha256") != QC_SHA256:
        raise ValueError("Figshare cell-death QC table checksum changed")
    if _digest(split_path, "sha256") != SPLITS_SHA256:
        raise ValueError("Figshare cell-death author split checksum changed")
    with qc_path.open(newline="", encoding="utf-8-sig") as handle:
        rows = [
            row for row in csv.DictReader(handle)
            if row["moa"] in PROGRAMMED_DEATH_LABELS
        ]

    label_compounds: dict[str, set[str]] = defaultdict(set)
    compound_labels: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        label = row["moa"]
        compound = row["Metadata_cmpdName"]
        label_compounds[label].add(compound)
        compound_labels[compound].add(label)
    conflicts = {
        compound: sorted(labels)
        for compound, labels in compound_labels.items() if len(labels) != 1
    }
    if conflicts:
        raise ValueError(f"compounds have conflicting author MoA labels: {conflicts}")

    assignments: dict[str, dict[str, str]] = {}
    for label in sorted(PROGRAMMED_DEATH_LABELS):
        assignments[label] = _assign(label, label_compounds[label])
    compound_role = {
        compound: role
        for label in assignments.values() for compound, role in label.items()
    }
    role_compounds = {
        role: sorted(compound for compound, value in compound_role.items() if value == role)
        for role in ROLE_ORDER
    }
    role_rows = Counter(compound_role[row["Metadata_cmpdName"]] for row in rows)
    role_wells = {
        role: len({
            (row["Metadata_Plate"], row["Metadata_Well"])
            for row in rows if compound_role[row["Metadata_cmpdName"]] == role
        })
        for role in ROLE_ORDER
    }
    role_labels = {
        role: sorted({
            row["moa"] for row in rows
            if compound_role[row["Metadata_cmpdName"]] == role
        })
        for role in ROLE_ORDER
    }

    author_split_audit: dict[str, object] = {}
    with tarfile.open(split_path, "r:gz") as archive:
        for member in sorted(name for name in archive.getnames() if name.endswith(".csv")):
            extracted = archive.extractfile(member)
            if extracted is None:
                raise ValueError(f"cannot read {member}")
            author_rows = list(csv.DictReader(io.TextIOWrapper(extracted, encoding="utf-8")))
            roles_by_compound: dict[str, set[str]] = defaultdict(set)
            for row in author_rows:
                roles_by_compound[row["Metadata_cmpdName"]].add(row["split"])
            leaking = {
                compound: sorted(roles)
                for compound, roles in roles_by_compound.items() if len(roles) > 1
            }
            author_split_audit[Path(member).stem] = {
                "compound_count": len(roles_by_compound),
                "compounds_crossing_author_roles": len(leaking),
                "compounds_in_all_three_author_roles": sum(
                    len(roles) == 3 for roles in leaking.values()
                ),
                "well_disjoint_but_not_compound_disjoint": bool(leaking),
            }

    role_sets = [set(value) for value in role_compounds.values()]
    overlap = sum(
        len(role_sets[left] & role_sets[right])
        for left in range(len(role_sets)) for right in range(left + 1, len(role_sets))
    )
    return {
        "schema": "aleph.outer_library.figshare_cell_death_design_audit.v1",
        "dataset_id": "figshare-28202864-mcf7-programmed-cell-death",
        "eligible_qc_site_rows": len(rows),
        "author_label_count": len(label_compounds),
        "compound_count": len(compound_role),
        "label_compound_counts": {
            label: len(compounds) for label, compounds in sorted(label_compounds.items())
        },
        "compound_assignments": assignments,
        "split_compound_counts": {
            role: len(compounds) for role, compounds in role_compounds.items()
        },
        "split_qc_site_counts": {role: role_rows[role] for role in ROLE_ORDER},
        "split_well_counts": role_wells,
        "split_label_coverage": role_labels,
        "cross_split_compound_overlap_count": overlap,
        "sealed_test_has_all_six_labels": len(role_labels["sealed_test"]) == 6,
        "selection_missing_rare_labels": sorted(
            PROGRAMMED_DEATH_LABELS - set(role_labels["selection"])
        ),
        "calibration_missing_rare_labels": sorted(
            PROGRAMMED_DEATH_LABELS - set(role_labels["calibration"])
        ),
        "rare_label_conditional_calibration_authority": False,
        "author_well_split_leakage_audit": author_split_audit,
        "split_policy": "author_MoA_stratified_compound_identity_disjoint_v1",
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.source_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        key: report[key] for key in (
            "compound_count", "split_compound_counts",
            "cross_split_compound_overlap_count"
        )
    }, sort_keys=True))


if __name__ == "__main__":
    main()
