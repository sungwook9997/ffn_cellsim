import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import train


def _row(domain, class_index, split, ordinal, *, gold=False, group=None):
    source = f"doi:10.0000/{class_index}-{split}-{ordinal}{'-gold' if gold else ''}"
    return {
        "authority_status": "proposed",
        "source_family_id": source,
        "leakage_group_id": group or f"lg:{class_index:02x}{split}{ordinal}",
        "split": split,
        "evaluation_cohort": "frozen_300_gold_candidate_not_ground_truth" if gold else "corpus",
        "supervision_policy": {"gold_candidate_is_ground_truth": False},
        "views": {
            "text_metadata": {
                "present": True, "domain": domain,
                "acquisition_pass": "first" if ordinal % 2 else "second",
                "assays": [f"assay_{class_index}"], "cell_types": [f"cell_{class_index % 4}"],
            },
            "numeric": {
                "present": ordinal % 3 != 0, "record_count": ordinal + 1,
                "observation_count": 3 + class_index, "candidate_target_count": class_index + 1,
                "ambiguous_target_mask_count": ordinal, "quantities": {f"q{class_index % 3}": 2},
                "roles": {"reported_measurement_candidate": 1},
            },
            "visual_metadata": {
                "present": ordinal % 2 == 0, "figure_count": ordinal + 1,
                "linked_figure_count": ordinal, "modality_tags": {f"visual_{class_index % 3}": 1},
            },
            "visual_asset_receipt": {"present": ordinal % 2 == 0, "acquired_asset_count": ordinal + 1},
            "visual_descriptor": {
                "present": ordinal % 2 == 0, "descriptor_count": ordinal + 1,
                "numeric_aggregates": {"entropy": {"mean": class_index + 0.5, "std": 0.1}},
                "descriptor_refs": [{"descriptor_sha256": "a" * 64}],
            },
            "dataset_manifest": {
                "present": ordinal % 4 == 0, "accessions": [f"GEO:GSE{class_index + 1000}"],
                "manifest_ids": [f"dataset:{class_index}"], "leakage_components": [f"ds:{class_index}"],
                "modality_hints": {f"hint_{class_index % 2}": 1},
            },
        },
    }


def _fixture(path):
    rows = []
    for index, domain in enumerate(train.DOMAINS):
        rows.append(_row(domain, index, "train", 1))
        rows.append(_row(domain, index, "train", 2))
        rows.append(_row(domain, index, "validation", 3))
        rows.append(_row(domain, index, "test", 4))
        rows.append(_row(domain, index, "train", 5, gold=True))
    path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows), encoding="utf-8")
    return rows


def test_feature_map_excludes_target_and_content_identifiers():
    row = _row(train.DOMAINS[0], 0, "train", 2)
    features = train.feature_map(row, train.VIEWS["metadata_numeric_visual_dataset"])
    assert not any(train.DOMAINS[0] in key for key in features)
    assert not any("sha256" in key or "descriptor_refs" in key for key in features)
    assert "visual:descriptor:entropy:mean" in features
    assert "missing:dataset_manifest" in features


def test_split_audit_kills_one_group_across_two_splits(tmp_path):
    path = tmp_path / "projection.jsonl"
    rows = [
        _row(train.DOMAINS[0], 0, "train", 1, group="lg:same"),
        _row(train.DOMAINS[0], 0, "test", 2, group="lg:same"),
    ]
    path.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")
    audit = train.audit_splits(rows, path)
    assert audit["leakage_split_violations"] == 1
    assert audit["errors"] >= 1


def test_three_seed_ablation_is_reproducible_and_gold_is_audit_only(tmp_path):
    projection = tmp_path / "projection.jsonl"
    _fixture(projection)
    first = tmp_path / "first"
    second = tmp_path / "second"
    one = train.run(projection, first, seeds=(11, 29, 47), max_epochs=25, patience=5)
    two = train.run(projection, second, seeds=(11, 29, 47), max_epochs=25, patience=5)
    assert one["audit"]["leakage_split_violations"] == 0
    assert one["metrics"]["seeds"] == [11, 29, 47]
    assert set(one["metrics"]["views"]) == set(train.VIEWS)
    assert one["gold_audit"]["sources"] == len(train.DOMAINS)
    assert one["gold_audit"]["performance_metrics_computed"] is False
    assert one["gold_audit"]["used_for_validation_model_selection"] is False
    assert one["hashes"] == two["hashes"]
    assert (first / "weights.npz").read_bytes() == (second / "weights.npz").read_bytes()


def test_metrics_include_support_confusion_balanced_accuracy_and_ece():
    y = np.asarray([0, 1, 1], dtype=np.int64)
    probability = np.full((3, len(train.DOMAINS)), 1e-6)
    probability[0, 0] = 0.9
    probability[1, 1] = 0.8
    probability[2, 0] = 0.7
    probability /= probability.sum(axis=1, keepdims=True)
    result = train.classification_metrics(y, probability)
    assert result["support"] == 3
    assert len(result["confusion_matrix"]) == len(train.DOMAINS)
    assert 0 <= result["balanced_accuracy"] <= 1
    assert 0 <= result["ece_10_bin"] <= 1


def test_projection_summary_digest_mismatch_refuses(tmp_path):
    projection = tmp_path / "projection.jsonl"
    _fixture(projection)
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({"projection_sha256": "0" * 64, "leakage_split_violations": 0}))
    with pytest.raises(RuntimeError, match="does not bind"):
        train.run(projection, tmp_path / "output", projection_summary=summary, seeds=(11,), max_epochs=2, patience=1)


def test_cli_forces_deterministic_blas_before_numpy_import(tmp_path):
    projection = tmp_path / "projection.jsonl"
    _fixture(projection)
    outputs = [tmp_path / "subprocess_one", tmp_path / "subprocess_two"]
    environment = dict(os.environ)
    # The trainer must override a caller's multi-thread settings before importing NumPy.
    environment.update({
        "OPENBLAS_NUM_THREADS": "8", "OMP_NUM_THREADS": "8", "MKL_NUM_THREADS": "8",
        "VECLIB_MAXIMUM_THREADS": "8", "NUMEXPR_NUM_THREADS": "8",
    })
    for output in outputs:
        subprocess.run(
            [
                sys.executable, str(Path(train.__file__)), "--projection", str(projection),
                "--output-dir", str(output), "--max-epochs", "8", "--patience", "3",
            ],
            check=True, capture_output=True, text=True, env=environment,
        )
    assert (outputs[0] / "SHA256SUMS.json").read_bytes() == (outputs[1] / "SHA256SUMS.json").read_bytes()
    for artifact in json.loads((outputs[0] / "SHA256SUMS.json").read_text()):
        assert (outputs[0] / artifact).read_bytes() == (outputs[1] / artifact).read_bytes()
