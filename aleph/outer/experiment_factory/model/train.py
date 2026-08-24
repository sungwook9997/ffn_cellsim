#!/usr/bin/env python
"""Deterministic NumPy weak broad-route baselines over the frozen learning projection."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import tempfile
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

# Dense reductions differed by up to 3.33e-16 under multi-threaded BLAS.  The CLI is
# a reproducibility artifact, so fix its reduction order before NumPy initializes.
for _thread_variable in (
    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import numpy as np


DOMAINS = [
    "adhesion_ecm_mechanics", "cell_state_mechanics", "cortical_membrane_pressure",
    "cytoskeleton_motor_rheology", "division_morphogenesis_tissue", "flow_shear",
    "measurement_inference", "mechanotransduction", "migration_invasion",
    "nuclear_mechanics", "organelle_general_mechanics", "osmotic_volume_poroelasticity",
]
VIEWS = {
    "metadata": ("metadata",),
    "metadata_numeric": ("metadata", "numeric"),
    "metadata_numeric_visual": ("metadata", "numeric", "visual"),
    "metadata_numeric_visual_dataset": ("metadata", "numeric", "visual", "dataset"),
}
SEEDS = (11, 29, 47)
SCHEMA = "aleph.external_training.weak_route_baseline.v1"


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def atomic(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        handle.write(data)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def load_rows(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def audit_splits(rows: list[dict[str, Any]], projection_path: Path) -> dict[str, Any]:
    errors: Counter[str] = Counter()
    sources = set()
    group_splits: dict[str, set[str]] = defaultdict(set)
    split_counts: Counter[str] = Counter()
    assignment_rows = []
    for row in rows:
        source = row["source_family_id"]
        if source in sources:
            errors["duplicate_source"] += 1
        sources.add(source)
        split = row["split"]
        if split not in {"train", "validation", "test"}:
            errors["invalid_split"] += 1
        group_splits[row["leakage_group_id"]].add(split)
        split_counts[split] += 1
        assignment_rows.append((source, row["leakage_group_id"], split))
        if row.get("authority_status") != "proposed":
            errors["authority_promotion"] += 1
        if row.get("supervision_policy", {}).get("gold_candidate_is_ground_truth") is not False:
            errors["gold_truth_policy"] += 1
    errors["leakage_split"] += sum(len(splits) != 1 for splits in group_splits.values())
    return {
        "schema": "aleph.external_training.model_split_audit.v1",
        "sources": len(rows),
        "unique_sources": len(sources),
        "leakage_groups": len(group_splits),
        "leakage_split_violations": errors["leakage_split"],
        "split_source_counts": dict(sorted(split_counts.items())),
        "split_assignment_sha256": digest(canonical(sorted(assignment_rows))),
        "projection_file_sha256": digest(projection_path.read_bytes()),
        "errors": sum(errors.values()),
        "error_classes": dict(sorted(errors.items())),
    }


def _flatten_numeric(prefix: str, value: Any, output: dict[str, float]) -> None:
    if isinstance(value, bool):
        output[prefix] = float(value)
    elif isinstance(value, (int, float)) and math.isfinite(float(value)):
        output[prefix] = math.log1p(max(float(value), 0.0))
    elif isinstance(value, dict):
        for key, child in sorted(value.items()):
            _flatten_numeric(f"{prefix}:{key}", child, output)


def feature_map(row: dict[str, Any], active: tuple[str, ...]) -> dict[str, float]:
    views = row["views"]
    result: dict[str, float] = {}
    metadata = views["text_metadata"]
    result[f"metadata:pass:{metadata['acquisition_pass']}"] = 1.0
    assays = metadata.get("assays") or []
    cells = metadata.get("cell_types") or []
    for value in assays:
        result[f"metadata:assay:{value}"] = 1.0
    for value in cells:
        result[f"metadata:cell:{value}"] = 1.0
    result["missing:metadata_assays"] = float(not assays)
    result["missing:metadata_cell_types"] = float(not cells)

    if "numeric" in active:
        numeric = views["numeric"]
        result["missing:numeric"] = float(not numeric.get("present"))
        if numeric.get("present"):
            for key in ("record_count", "observation_count", "candidate_target_count", "ambiguous_target_mask_count"):
                _flatten_numeric(f"numeric:{key}", numeric.get(key, 0), result)
            for family in ("quantities", "roles"):
                for key, value in sorted((numeric.get(family) or {}).items()):
                    _flatten_numeric(f"numeric:{family}:{key}", value, result)

    if "visual" in active:
        visual = views["visual_metadata"]
        descriptor = views["visual_descriptor"]
        asset = views.get("visual_asset_receipt") or {"present": False}
        result["missing:visual_metadata"] = float(not visual.get("present"))
        result["missing:visual_descriptor"] = float(not descriptor.get("present"))
        result["missing:visual_asset_receipt"] = float(not asset.get("present"))
        if asset.get("present"):
            _flatten_numeric("visual:acquired_asset_count", asset.get("acquired_asset_count", 0), result)
        if visual.get("present"):
            for key in ("figure_count", "linked_figure_count", "available_asset_count"):
                _flatten_numeric(f"visual:{key}", visual.get(key, 0), result)
            for key, value in sorted((visual.get("modality_tags") or {}).items()):
                _flatten_numeric(f"visual:modality:{key}", value, result)
        if descriptor.get("present"):
            _flatten_numeric("visual:descriptor_count", descriptor.get("descriptor_count", 0), result)
            # The final projection may add source-level numeric aggregates. Content identifiers and
            # digests are deliberately never feature-hashed.
            aggregates = descriptor.get("aggregates") or descriptor.get("numeric_aggregates") or {}
            _flatten_numeric("visual:descriptor", aggregates, result)
            result["missing:visual_descriptor_values"] = float(not bool(aggregates))

    if "dataset" in active:
        dataset = views["dataset_manifest"]
        result["missing:dataset_manifest"] = float(not dataset.get("present"))
        if dataset.get("present"):
            _flatten_numeric("dataset:accession_count", len(dataset.get("accessions") or []), result)
            _flatten_numeric("dataset:manifest_count", len(dataset.get("manifest_ids") or []), result)
            _flatten_numeric("dataset:component_count", len(dataset.get("leakage_components") or []), result)
            for key, value in sorted((dataset.get("modality_hints") or {}).items()):
                _flatten_numeric(f"dataset:modality:{key}", value, result)
    return result


class Vectorizer:
    def __init__(self, active: tuple[str, ...]):
        self.active = active
        self.columns: list[str] = []
        self.mean: np.ndarray | None = None
        self.scale: np.ndarray | None = None

    def fit(self, rows: list[dict[str, Any]]) -> "Vectorizer":
        self.columns = sorted({key for row in rows for key in feature_map(row, self.active)})
        raw = self._raw(rows)
        self.mean = raw.mean(axis=0)
        scale = raw.std(axis=0)
        self.scale = np.where(scale > 1e-12, scale, 1.0)
        return self

    def _raw(self, rows: list[dict[str, Any]]) -> np.ndarray:
        index = {name: i for i, name in enumerate(self.columns)}
        matrix = np.zeros((len(rows), len(self.columns)), dtype=np.float64)
        for row_index, row in enumerate(rows):
            for key, value in feature_map(row, self.active).items():
                column = index.get(key)
                if column is not None:
                    matrix[row_index, column] = value
        return matrix

    def transform(self, rows: list[dict[str, Any]]) -> np.ndarray:
        if self.mean is None or self.scale is None:
            raise RuntimeError("vectorizer is not fitted")
        raw = self._raw(rows)
        standardized = (raw - self.mean) / self.scale
        return np.concatenate([np.ones((len(rows), 1), dtype=np.float64), standardized], axis=1)


def labels(rows: list[dict[str, Any]]) -> np.ndarray:
    index = {label: i for i, label in enumerate(DOMAINS)}
    try:
        return np.asarray([index[row["views"]["text_metadata"]["domain"]] for row in rows], dtype=np.int64)
    except KeyError as exc:
        raise ValueError(f"unknown weak routing label: {exc}") from exc


def confusion_matrix(y_true: np.ndarray, y_pred: np.ndarray, classes: int) -> np.ndarray:
    matrix = np.zeros((classes, classes), dtype=np.int64)
    np.add.at(matrix, (y_true, y_pred), 1)
    return matrix


def classification_metrics(y_true: np.ndarray, probabilities: np.ndarray) -> dict[str, Any]:
    predictions = probabilities.argmax(axis=1)
    matrix = confusion_matrix(y_true, predictions, len(DOMAINS))
    supports = matrix.sum(axis=1)
    tp = np.diag(matrix).astype(np.float64)
    precision = np.divide(tp, matrix.sum(axis=0), out=np.zeros_like(tp), where=matrix.sum(axis=0) != 0)
    recall = np.divide(tp, supports, out=np.zeros_like(tp), where=supports != 0)
    f1 = np.divide(2 * precision * recall, precision + recall, out=np.zeros_like(tp), where=(precision + recall) != 0)
    macro_f1 = float(f1.mean())
    micro_f1 = float(tp.sum() / max(matrix.sum(), 1))
    balanced = float(recall.mean())
    confidence = probabilities.max(axis=1)
    correct = predictions == y_true
    ece = 0.0
    bins = np.linspace(0.0, 1.0, 11)
    for left, right in zip(bins[:-1], bins[1:]):
        mask = (confidence >= left) & (confidence < right if right < 1 else confidence <= right)
        if mask.any():
            ece += float(mask.mean()) * abs(float(correct[mask].mean()) - float(confidence[mask].mean()))
    per_class = {
        label: {"support": int(supports[i]), "precision": float(precision[i]), "recall": float(recall[i]), "f1": float(f1[i])}
        for i, label in enumerate(DOMAINS)
    }
    return {
        "macro_f1": macro_f1,
        "micro_f1": micro_f1,
        "balanced_accuracy": balanced,
        "ece_10_bin": ece,
        "support": int(len(y_true)),
        "per_class": per_class,
        "confusion_matrix": matrix.tolist(),
    }


def softmax(logits: np.ndarray) -> np.ndarray:
    shifted = logits - logits.max(axis=1, keepdims=True)
    values = np.exp(shifted)
    return values / values.sum(axis=1, keepdims=True)


def train_linear(
    x_train: np.ndarray, y_train: np.ndarray, x_validation: np.ndarray, y_validation: np.ndarray,
    *, seed: int, max_epochs: int, patience: int, learning_rate: float, l2: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    rng = np.random.default_rng(seed)
    weights = rng.normal(0.0, 0.01, size=(x_train.shape[1], len(DOMAINS)))
    counts = np.bincount(y_train, minlength=len(DOMAINS)).astype(np.float64)
    if np.any(counts == 0):
        raise ValueError("every weak route must have train support")
    class_weights = len(y_train) / (len(DOMAINS) * counts)
    sample_weights = class_weights[y_train]
    m = np.zeros_like(weights)
    v = np.zeros_like(weights)
    best = weights.copy()
    best_score = -1.0
    best_epoch = 0
    stale = 0
    for epoch in range(1, max_epochs + 1):
        probabilities = softmax(x_train @ weights)
        gradient_logits = probabilities
        gradient_logits[np.arange(len(y_train)), y_train] -= 1.0
        gradient_logits *= sample_weights[:, None] / sample_weights.sum()
        gradient = x_train.T @ gradient_logits + l2 * np.vstack([np.zeros((1, weights.shape[1])), weights[1:]])
        m = 0.9 * m + 0.1 * gradient
        v = 0.999 * v + 0.001 * gradient * gradient
        m_hat = m / (1 - 0.9 ** epoch)
        v_hat = v / (1 - 0.999 ** epoch)
        weights -= learning_rate * m_hat / (np.sqrt(v_hat) + 1e-8)
        score = classification_metrics(y_validation, softmax(x_validation @ weights))["macro_f1"]
        if score > best_score + 1e-12:
            best_score = score
            best_epoch = epoch
            best = weights.copy()
            stale = 0
        else:
            stale += 1
            if stale >= patience:
                break
    return best, {
        "best_epoch": best_epoch,
        "epochs_run": epoch,
        "validation_macro_f1": best_score,
        "class_weights": {DOMAINS[i]: float(class_weights[i]) for i in range(len(DOMAINS))},
    }


def aggregate_seed_metrics(values: list[dict[str, Any]]) -> dict[str, Any]:
    keys = ("macro_f1", "micro_f1", "balanced_accuracy", "ece_10_bin")
    return {
        key: {
            "mean": float(np.mean([value[key] for value in values])),
            "std": float(np.std([value[key] for value in values], ddof=0)),
            "min": float(np.min([value[key] for value in values])),
            "max": float(np.max([value[key] for value in values])),
        }
        for key in keys
    }


def gold_representation_audit(
    gold_rows: list[dict[str, Any]], vectorizers: dict[str, Vectorizer], primary_weights: dict[str, np.ndarray],
) -> dict[str, Any]:
    coverage = {
        name: sum(row["views"][name]["present"] for row in gold_rows)
        for name in ("numeric", "visual_metadata", "visual_asset_receipt", "visual_descriptor", "dataset_manifest")
    }
    predictions = {}
    for view, vectorizer in vectorizers.items():
        probability = softmax(vectorizer.transform(gold_rows) @ primary_weights[view])
        distribution = Counter(DOMAINS[index] for index in probability.argmax(axis=1))
        predictions[view] = {
            "predicted_route_distribution": dict(sorted(distribution.items())),
            "mean_max_probability": float(probability.max(axis=1).mean()),
        }
    return {
        "schema": "aleph.external_training.gold_candidate_representation_audit.v1",
        "authority_status": "proposed",
        "sources": len(gold_rows),
        "view_coverage": coverage,
        "prediction_representation_only": predictions,
        "performance_metrics_computed": False,
        "used_for_training": False,
        "used_for_validation_model_selection": False,
        "used_for_test_metrics": False,
        "limitation": "gold_candidate is a sampling cohort, not label truth",
    }


def run(
    projection: Path, output_dir: Path, *, projection_summary: Path | None = None,
    seeds: tuple[int, ...] = SEEDS, max_epochs: int = 300, patience: int = 35,
    learning_rate: float = 0.03, l2: float = 1e-4,
) -> dict[str, Any]:
    rows = load_rows(projection)
    audit = audit_splits(rows, projection)
    if audit["errors"]:
        raise RuntimeError(f"split audit failed: {audit}")
    if projection_summary is not None:
        summary = json.loads(projection_summary.read_text(encoding="utf-8"))
        if summary["projection_sha256"] != audit["projection_file_sha256"]:
            raise RuntimeError("projection summary SHA does not bind the supplied projection")
        if summary["leakage_split_violations"] != 0:
            raise RuntimeError("projection summary reports leakage violations")
    corpus_all = [row for row in rows if row["evaluation_cohort"] == "corpus"]
    unmapped_excluded = Counter(
        row["split"] for row in corpus_all
        if row["views"]["text_metadata"]["domain"] not in DOMAINS
    )
    corpus = [
        row for row in corpus_all
        if row["views"]["text_metadata"]["domain"] in DOMAINS
    ]
    gold_rows = [row for row in rows if row["evaluation_cohort"].startswith("frozen_300")]
    partitions = {split: [row for row in corpus if row["split"] == split] for split in ("train", "validation", "test")}
    if any(not partition for partition in partitions.values()):
        raise ValueError("all corpus-only partitions must be non-empty")
    y = {split: labels(partition) for split, partition in partitions.items()}
    metrics: dict[str, Any] = {
        "schema": SCHEMA,
        "authority_status": "proposed",
        "target": "weak_12_class_discovery_route_only",
        "seeds": list(seeds),
        "corpus_only_split_counts": {key: len(value) for key, value in partitions.items()},
        "unmapped_weak_route_excluded_counts": {
            split: unmapped_excluded[split] for split in ("train", "validation", "test")
        },
        "views": {},
    }
    vectorizers: dict[str, Vectorizer] = {}
    primary_weights: dict[str, np.ndarray] = {}
    weights_artifact: dict[str, np.ndarray] = {}
    feature_schemas = {}
    for view, active in VIEWS.items():
        vectorizer = Vectorizer(active).fit(partitions["train"])
        vectorizers[view] = vectorizer
        matrices = {split: vectorizer.transform(partition) for split, partition in partitions.items()}
        seed_rows = []
        for seed in seeds:
            weights, training = train_linear(
                matrices["train"], y["train"], matrices["validation"], y["validation"],
                seed=seed, max_epochs=max_epochs, patience=patience,
                learning_rate=learning_rate, l2=l2,
            )
            test_metrics = classification_metrics(y["test"], softmax(matrices["test"] @ weights))
            seed_rows.append({"seed": seed, "training": training, "test": test_metrics})
            weights_artifact[f"{view}__seed_{seed}"] = weights
            if seed == seeds[0]:
                primary_weights[view] = weights
        metrics["views"][view] = {
            "active_feature_groups": list(active),
            "feature_count_including_intercept": matrices["train"].shape[1],
            "per_seed": seed_rows,
            "test_aggregate": aggregate_seed_metrics([row["test"] for row in seed_rows]),
        }
        feature_schemas[view] = {
            "columns": vectorizer.columns,
            "train_mean": vectorizer.mean.tolist(),
            "train_scale": vectorizer.scale.tolist(),
        }
    gold_audit = gold_representation_audit(gold_rows, vectorizers, primary_weights)
    output_dir.mkdir(parents=True, exist_ok=True)
    atomic(output_dir / "split_audit.json", json.dumps(audit, indent=2, sort_keys=True).encode() + b"\n")
    atomic(output_dir / "metrics.json", json.dumps(metrics, indent=2, sort_keys=True).encode() + b"\n")
    atomic(output_dir / "feature_schemas.json", json.dumps(feature_schemas, indent=2, sort_keys=True).encode() + b"\n")
    atomic(output_dir / "gold_candidate_representation_audit.json", json.dumps(gold_audit, indent=2, sort_keys=True).encode() + b"\n")
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as handle:
        for name in sorted(weights_artifact):
            payload = io.BytesIO()
            np.save(payload, weights_artifact[name], allow_pickle=False)
            entry = zipfile.ZipInfo(f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_DEFLATED
            entry.external_attr = 0o600 << 16
            handle.writestr(entry, payload.getvalue())
    atomic(output_dir / "weights.npz", archive.getvalue())
    config = {
        "schema": SCHEMA, "authority_status": "proposed", "seeds": list(seeds),
        "max_epochs": max_epochs, "patience": patience, "learning_rate": learning_rate,
        "l2": l2, "class_weighting": "balanced inverse train support",
        "early_stopping": "validation macro-F1 only", "test_used_for_model_selection": False,
        "gold_candidate_used_for_model_selection": False,
        "projection_sha256": audit["projection_file_sha256"],
        "split_assignment_sha256": audit["split_assignment_sha256"],
    }
    atomic(output_dir / "config.json", json.dumps(config, indent=2, sort_keys=True).encode() + b"\n")
    hashes = {
        name: digest((output_dir / name).read_bytes())
        for name in ("config.json", "feature_schemas.json", "gold_candidate_representation_audit.json", "metrics.json", "split_audit.json", "weights.npz")
    }
    atomic(output_dir / "SHA256SUMS.json", json.dumps(hashes, indent=2, sort_keys=True).encode() + b"\n")
    return {"metrics": metrics, "audit": audit, "gold_audit": gold_audit, "hashes": hashes}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--projection-summary", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-epochs", type=int, default=300)
    parser.add_argument("--patience", type=int, default=35)
    args = parser.parse_args(argv)
    result = run(
        args.projection, args.output_dir, projection_summary=args.projection_summary,
        max_epochs=args.max_epochs, patience=args.patience,
    )
    print(json.dumps({"split_audit": result["audit"], "artifact_hashes": result["hashes"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
