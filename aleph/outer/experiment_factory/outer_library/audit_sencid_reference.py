#!/usr/bin/env python3
"""Audit the public SenCID model snapshot without executing pickle artifacts."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
import subprocess


EXPECTED_COMMIT = "f93cba19b9d44e84f89504ecc6aaeb108ce8e508"
MODEL_NAMES = tuple(
    [name for sid in range(1, 7) for name in (f"SID{sid}.pkl", f"SID{sid}_L.pkl")]
    + ["recommend_model.pkl"]
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _feature_sets(prediction_source: Path) -> dict[int, list[str]]:
    tree = ast.parse(prediction_source.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "feature_dict"
            for target in node.targets
        ):
            value = ast.literal_eval(node.value)
            if sorted(value) != list(range(1, 7)):
                raise ValueError("SenCID SID feature dictionary changed")
            return {int(key): list(features) for key, features in value.items()}
    raise ValueError("SenCID SID feature dictionary not found")


def audit(source_dir: Path) -> dict[str, object]:
    commit = subprocess.run(
        ["git", "-C", str(source_dir), "rev-parse", "HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    if commit != EXPECTED_COMMIT:
        raise ValueError(f"SenCID source commit changed: {commit}")
    license_text = (source_dir / "LICENSE").read_text(encoding="utf-8")
    if not license_text.startswith("MIT License"):
        raise ValueError("SenCID license text changed")

    model_dir = source_dir / "SenCID" / "model"
    model_files = {
        name: {
            "bytes": (model_dir / name).stat().st_size,
            "sha256": _sha256(model_dir / name),
            "executed_during_audit": False,
        }
        for name in MODEL_NAMES
    }
    features = _feature_sets(source_dir / "SenCID" / "Pred.py")
    seneset = [
        line.strip() for line in
        (source_dir / "SenCID" / "resource" / "seneset.txt").read_text().splitlines()
        if line.strip()
    ]
    tracked = subprocess.run(
        ["git", "-C", str(source_dir), "ls-files"],
        check=True, capture_output=True, text=True,
    ).stdout.splitlines()
    data_like = [
        name for name in tracked
        if Path(name).suffix.lower() in {".csv", ".h5", ".h5ad", ".mtx", ".tsv"}
    ]
    training_code_present = any(
        "train" in Path(name).name.lower() and Path(name).suffix == ".py"
        for name in tracked
    )
    return {
        "schema": "aleph.outer_library.sencid_reference_audit.v1",
        "source": {
            "repository": "https://github.com/JackieHanLab/SenCID",
            "commit": commit,
            "license": "MIT",
            "paper_doi": "10.1016/j.cmet.2024.03.009",
            "paper_reported_training_scope": {
                "samples": 602, "studies": 52, "cell_types": 30, "senescence_identities": 6,
            },
        },
        "model_artifacts": model_files,
        "pickle_execution_policy": "never_execute_during_corpus_audit",
        "sid_feature_counts": {str(key): len(value) for key, value in features.items()},
        "unique_sid_features": len(set().union(*(set(value) for value in features.values()))),
        "seneset_gene_count": len(seneset),
        "tracked_file_count": len(tracked),
        "tracked_training_matrix_files": data_like,
        "training_code_present": training_code_present,
        "study_level_training_manifest_present": False,
        "independent_retraining_possible_from_repository_alone": False,
        "preprocessing_risks": [
            "default path depends on DCA and TensorFlow versions older than the Aleph environment",
            "Zcol fits a scaler independently within every sample across genes",
            "joblib pickle artifacts are executable serialization and are not loaded by this audit",
        ],
        "eligible_role": "frozen_reference_expert_pending_untouched_dataset_evaluation",
        "external_validation_eligible": False,
        "reason_not_external_validation": (
            "the public repository ships models but not the 52-study training matrix, labels, "
            "or study-level split manifest needed to exclude overlap and reproduce fitting"
        ),
        "may_select_aleph_parameter": False,
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.source_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "sid_feature_counts", "unique_sid_features", "seneset_gene_count",
        "independent_retraining_possible_from_repository_alone"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()
