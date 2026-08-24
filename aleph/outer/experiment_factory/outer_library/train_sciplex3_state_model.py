#!/usr/bin/env python3
"""Train a treatment-disjoint sci-Plex3 pathway-state representation head."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np


SPLIT_NAMES = ("train", "selection", "calibration", "test", "vehicle_ood")


def _softmax(logits: np.ndarray) -> np.ndarray:
    shifted = logits - logits.max(axis=1, keepdims=True)
    exponential = np.exp(shifted)
    return exponential / exponential.sum(axis=1, keepdims=True)


def _forward(
    X: np.ndarray, params: dict[str, np.ndarray]
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    hidden = np.maximum(0.0, X @ params["W_shared"] + params["b_shared"])
    logits = {
        name: hidden @ params[f"W_{name}"] + params[f"b_{name}"]
        for name in ("pathway", "cell_type", "time")
    }
    return hidden, logits


def _metrics(
    probability: np.ndarray,
    labels: np.ndarray,
    class_names: np.ndarray,
    eligible: np.ndarray | None = None,
) -> dict[str, object]:
    prediction = probability.argmax(axis=1)
    class_count = len(class_names)
    confusion = np.zeros((class_count, class_count), dtype=np.int64)
    np.add.at(confusion, (labels, prediction), 1)
    if eligible is None:
        eligible = confusion.sum(axis=1) > 0
    f1: list[float | None] = []
    recall: list[float | None] = []
    for label in range(class_count):
        tp = int(confusion[label, label])
        fp = int(confusion[:, label].sum() - tp)
        fn = int(confusion[label].sum() - tp)
        f1.append(float(2 * tp / max(1, 2 * tp + fp + fn)) if eligible[label] else None)
        recall.append(float(tp / max(1, tp + fn)) if eligible[label] else None)
    confidence = probability.max(axis=1)
    correct = prediction == labels
    ece = 0.0
    for lower in np.linspace(0.0, 0.9, 10):
        mask = (confidence >= lower) & (confidence < lower + 0.1 + 1e-12)
        if mask.any():
            ece += float(mask.mean() * abs(correct[mask].mean() - confidence[mask].mean()))
    eligible_f1 = [value for value in f1 if value is not None]
    eligible_recall = [value for value in recall if value is not None]
    return {
        "sample_count": int(len(labels)),
        "accuracy": float(correct.mean()),
        "balanced_accuracy": float(np.mean(eligible_recall)),
        "macro_f1": float(np.mean(eligible_f1)),
        "ece_10_bin": ece,
        "eligible_class_count": int(eligible.sum()),
        "per_class_f1": {
            str(name): value for name, value in zip(class_names, f1, strict=True)
        },
        "confusion_matrix_true_by_predicted": confusion.tolist(),
    }


def _aggregate_treatments(
    logits: np.ndarray,
    labels: np.ndarray,
    treatments: np.ndarray,
    mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    aggregate_logits: list[np.ndarray] = []
    aggregate_labels: list[int] = []
    names: list[str] = []
    for treatment in sorted(set(treatments[mask])):
        rows = mask & (treatments == treatment)
        values = set(labels[rows].tolist())
        if len(values) != 1:
            raise ValueError(f"pathway labels change within treatment {treatment}")
        aggregate_logits.append(logits[rows].mean(axis=0))
        aggregate_labels.append(values.pop())
        names.append(str(treatment))
    return (
        np.asarray(aggregate_logits, dtype=np.float64),
        np.asarray(aggregate_labels, dtype=np.int64),
        np.asarray(names, dtype="U16"),
    )


def _temperature(
    logits: np.ndarray, labels: np.ndarray
) -> tuple[float, np.ndarray]:
    grid = np.linspace(0.35, 4.0, 147)
    losses = []
    for value in grid:
        probability = _softmax(logits / value)
        losses.append(
            float(-np.log(probability[np.arange(len(labels)), labels] + 1e-12).mean())
        )
    selected = float(grid[int(np.argmin(losses))])
    return selected, np.asarray(losses)


def _replicate_consistency(
    hidden: np.ndarray, context: np.ndarray
) -> dict[str, object]:
    # Match on every declared experimental factor except the replicate.
    groups: dict[tuple[str, str, str, str], dict[str, list[int]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for index, row in enumerate(context):
        key = (str(row[0]), str(row[2]), str(row[3]), str(row[4]))
        groups[key][str(row[1])].append(index)
    similarities: list[float] = []
    for replicates in groups.values():
        if set(replicates) != {"rep1", "rep2"}:
            continue
        left = hidden[replicates["rep1"]].mean(axis=0)
        right = hidden[replicates["rep2"]].mean(axis=0)
        denominator = np.linalg.norm(left) * np.linalg.norm(right)
        if denominator > 0:
            similarities.append(float(left @ right / denominator))
    return {
        "matched_context_pair_count": len(similarities),
        "mean_latent_cosine_similarity": float(np.mean(similarities)),
        "median_latent_cosine_similarity": float(np.median(similarities)),
        "same_lab_technical_or_biological_status": "same_provider_replicates_not_external_labs",
    }


def train(
    tensor_path: Path,
    checkpoint_path: Path,
    *,
    epochs: int = 80,
    hidden_width: int = 96,
    selected_gene_count: int = 768,
) -> dict[str, object]:
    with np.load(tensor_path, allow_pickle=False) as data:
        required = {
            "log1p_cpm", "context", "context_field", "context_metadata",
            "context_metadata_field", "gene_symbol", "treatment_split",
        }
        if not required.issubset(data.files):
            raise ValueError(f"sci-Plex3 tensor misses {sorted(required - set(data.files))}")
        expression = np.asarray(data["log1p_cpm"], dtype=np.float32)
        context = np.asarray(data["context"]).astype(str)
        context_fields = np.asarray(data["context_field"]).astype(str)
        metadata = np.asarray(data["context_metadata"]).astype(str)
        metadata_fields = np.asarray(data["context_metadata_field"]).astype(str)
        gene_symbols = np.asarray(data["gene_symbol"]).astype(str)
        split = np.asarray(data["treatment_split"], dtype=np.int64)
    if context_fields.tolist() != ["cell_type", "replicate", "time_point", "treatment", "dose"]:
        raise ValueError("sci-Plex3 context field order changed")
    metadata_index = {name: i for i, name in enumerate(metadata_fields)}
    needed_metadata = {"vehicle", "pathway_level_1"}
    if not needed_metadata.issubset(metadata_index):
        raise ValueError("sci-Plex3 tensor lacks pathway/vehicle metadata")
    if expression.shape != (len(context), len(gene_symbols)) or not np.isfinite(expression).all():
        raise ValueError(f"invalid sci-Plex3 expression tensor {expression.shape}")
    if set(np.unique(split)) != {0, 1, 2, 3, 4}:
        raise ValueError("sci-Plex3 split must contain four supervised partitions and vehicle OOD")

    pathway_text = metadata[:, metadata_index["pathway_level_1"]]
    vehicle = metadata[:, metadata_index["vehicle"]] == "TRUE"
    if not np.array_equal(vehicle, split == 4):
        raise ValueError("vehicle metadata and OOD split disagree")
    pathway_names = np.asarray(sorted(set(pathway_text[~vehicle])), dtype="U64")
    pathway_index = {name: i for i, name in enumerate(pathway_names)}
    pathway = np.asarray(
        [pathway_index.get(name, -1) for name in pathway_text], dtype=np.int64
    )
    cell_type_names = np.asarray(sorted(set(context[:, 0])), dtype="U16")
    cell_type_index = {name: i for i, name in enumerate(cell_type_names)}
    cell_type = np.asarray([cell_type_index[name] for name in context[:, 0]], dtype=np.int64)
    time_names = np.asarray(sorted(set(context[:, 2])), dtype="U16")
    time_index = {name: i for i, name in enumerate(time_names)}
    time_label = np.asarray([time_index[name] for name in context[:, 2]], dtype=np.int64)
    treatments = context[:, 3]
    treatment_split_sets = {
        name: set(treatments[split == index]) for index, name in enumerate(SPLIT_NAMES)
    }
    overlap = sum(
        bool(treatment_split_sets[left] & treatment_split_sets[right])
        for i, left in enumerate(SPLIT_NAMES)
        for right in SPLIT_NAMES[i + 1:]
    )
    if overlap:
        raise ValueError("sci-Plex3 treatment leaked across splits")

    train_mask = split == 0
    # Remove cell-line/replicate/time baselines using fitting compounds only.
    # This asks the pathway head to learn perturbational response rather than
    # the much larger constitutive cell-line signal.  Sealed compounds never
    # contribute to a nuisance mean.
    residual_expression = expression.copy()
    nuisance_group_count = 0
    for cell_type_value in cell_type_names:
        for replicate_value in sorted(set(context[:, 1])):
            for time_value in time_names:
                group = (
                    (context[:, 0] == cell_type_value)
                    & (context[:, 1] == replicate_value)
                    & (context[:, 2] == time_value)
                )
                fitting_group = group & train_mask
                if fitting_group.any():
                    residual_expression[group] -= expression[fitting_group].mean(axis=0)
                    nuisance_group_count += 1
    if nuisance_group_count != len(set(map(tuple, context[:, :3]))):
        raise ValueError("a sci-Plex nuisance group lacks fitting compounds")
    variance = residual_expression[train_mask].var(axis=0)
    selected_gene_count = min(selected_gene_count, expression.shape[1])
    selected_genes = np.argsort(variance, kind="stable")[-selected_gene_count:]
    X = residual_expression[:, selected_genes]
    mean = X[train_mask].mean(axis=0)
    scale = X[train_mask].std(axis=0)
    scale[scale < 1e-6] = 1.0
    X = ((X - mean) / scale).astype(np.float32)

    rng = np.random.default_rng(139944)
    output_sizes = {
        "pathway": len(pathway_names),
        "cell_type": len(cell_type_names),
        "time": len(time_names),
    }
    params: dict[str, np.ndarray] = {
        "W_shared": (
            rng.standard_normal((selected_gene_count, hidden_width))
            * np.sqrt(2.0 / selected_gene_count)
        ).astype(np.float32),
        "b_shared": np.zeros(hidden_width, dtype=np.float32),
    }
    for name, size in output_sizes.items():
        params[f"W_{name}"] = (
            rng.standard_normal((hidden_width, size)) * np.sqrt(2.0 / hidden_width)
        ).astype(np.float32)
        params[f"b_{name}"] = np.zeros(size, dtype=np.float32)
    moments = {key: np.zeros_like(value) for key, value in params.items()}
    velocities = {key: np.zeros_like(value) for key, value in params.items()}

    # Each compound receives equal total pathway weight before class balancing.
    treatment_context_count = Counter(treatments[train_mask].tolist())
    pathway_treatment_count = Counter(
        {
            treatment: int(pathway[train_mask & (treatments == treatment)][0])
            for treatment in treatment_context_count
        }.values()
    )
    pathway_weight = np.asarray(
        [
            np.sqrt(len(treatment_context_count) / (len(pathway_names) * pathway_treatment_count[i]))
            for i in range(len(pathway_names))
        ],
        dtype=np.float32,
    )
    pathway_weight = np.minimum(pathway_weight, 8.0)
    cell_counts = np.bincount(cell_type[train_mask], minlength=len(cell_type_names))
    time_counts = np.bincount(time_label[train_mask], minlength=len(time_names))
    cell_weight = np.sqrt(train_mask.sum() / (len(cell_type_names) * cell_counts))
    time_weight = np.sqrt(train_mask.sum() / (len(time_names) * time_counts))
    task_weight = {"pathway": 1.0, "cell_type": 0.2, "time": 0.2}

    eligible_pathway = np.asarray(
        [all(np.any((split == part) & (pathway == label)) for part in range(4))
         for label in range(len(pathway_names))],
        dtype=bool,
    )
    best_params = {key: value.copy() for key, value in params.items()}
    best_selection_macro_f1 = -1.0
    best_epoch = 0
    train_indices = np.flatnonzero(train_mask)
    step = 0
    for epoch in range(epochs):
        for batch in np.array_split(rng.permutation(train_indices), max(1, len(train_indices) // 128)):
            xb = X[batch]
            hidden, logits = _forward(xb, params)
            labels_by_task = {
                "pathway": pathway[batch],
                "cell_type": cell_type[batch],
                "time": time_label[batch],
            }
            row_weights = {
                "pathway": np.asarray(
                    [
                        pathway_weight[label] / treatment_context_count[str(treatment)]
                        for label, treatment in zip(pathway[batch], treatments[batch], strict=True)
                    ],
                    dtype=np.float32,
                ),
                "cell_type": cell_weight[cell_type[batch]],
                "time": time_weight[time_label[batch]],
            }
            gradients: dict[str, np.ndarray] = {}
            dhidden = np.zeros_like(hidden)
            for name in ("pathway", "cell_type", "time"):
                labels = labels_by_task[name]
                probability = _softmax(logits[name])
                dlogits = probability
                dlogits[np.arange(len(batch)), labels] -= 1.0
                weights = row_weights[name]
                dlogits *= (task_weight[name] * weights / weights.sum())[:, None]
                gradients[f"W_{name}"] = hidden.T @ dlogits + 1e-4 * params[f"W_{name}"]
                gradients[f"b_{name}"] = dlogits.sum(axis=0)
                dhidden += dlogits @ params[f"W_{name}"].T
            dhidden[hidden <= 0] = 0
            gradients["W_shared"] = xb.T @ dhidden + 1e-4 * params["W_shared"]
            gradients["b_shared"] = dhidden.sum(axis=0)
            step += 1
            for key in params:
                moments[key] = 0.9 * moments[key] + 0.1 * gradients[key]
                velocities[key] = 0.999 * velocities[key] + 0.001 * gradients[key] ** 2
                mhat = moments[key] / (1.0 - 0.9 ** step)
                vhat = velocities[key] / (1.0 - 0.999 ** step)
                params[key] -= 7e-4 * mhat / (np.sqrt(vhat) + 1e-8)
        _, all_logits = _forward(X, params)
        selection_logits, selection_labels, _ = _aggregate_treatments(
            all_logits["pathway"], pathway, treatments, split == 1
        )
        score = _metrics(
            _softmax(selection_logits), selection_labels, pathway_names, eligible_pathway
        )["macro_f1"]
        if float(score) > best_selection_macro_f1:
            best_selection_macro_f1 = float(score)
            best_epoch = epoch + 1
            best_params = {key: value.copy() for key, value in params.items()}

    hidden, logits = _forward(X, best_params)
    calibration_logits, calibration_labels, calibration_treatments = _aggregate_treatments(
        logits["pathway"], pathway, treatments, split == 2
    )
    temperature, _ = _temperature(calibration_logits, calibration_labels)
    calibration_probability = _softmax(calibration_logits / temperature)
    conformal_scores = 1.0 - calibration_probability[
        np.arange(len(calibration_labels)), calibration_labels
    ]
    rank = min(
        len(conformal_scores) - 1,
        int(np.ceil((len(conformal_scores) + 1) * 0.90)) - 1,
    )
    conformal_threshold = float(np.partition(conformal_scores, rank)[rank])

    test_logits, test_labels, test_treatments = _aggregate_treatments(
        logits["pathway"], pathway, treatments, split == 3
    )
    test_probability = _softmax(test_logits / temperature)
    treatment_metrics = _metrics(
        test_probability, test_labels, pathway_names, eligible_pathway
    )
    test_sets = test_probability >= (1.0 - conformal_threshold)
    treatment_coverage = float(
        test_sets[np.arange(len(test_labels)), test_labels].mean()
    )
    test_context_probability = _softmax(logits["pathway"][split == 3] / temperature)
    context_metrics = _metrics(
        test_context_probability, pathway[split == 3], pathway_names, eligible_pathway
    )

    calibration_hidden = hidden[split == 2]
    train_hidden = hidden[train_mask]
    hidden_mean = train_hidden.mean(axis=0)
    hidden_scale = train_hidden.std(axis=0)
    hidden_scale[hidden_scale < 1e-6] = 1.0
    calibration_distance = np.sqrt(
        np.mean(((calibration_hidden - hidden_mean) / hidden_scale) ** 2, axis=1)
    )
    ood_rank = min(
        len(calibration_distance) - 1,
        int(np.ceil((len(calibration_distance) + 1) * 0.95)) - 1,
    )
    ood_threshold = float(np.partition(calibration_distance, ood_rank)[ood_rank])
    vehicle_distance = np.sqrt(
        np.mean(((hidden[vehicle] - hidden_mean) / hidden_scale) ** 2, axis=1)
    )

    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        checkpoint_path,
        **best_params,
        selected_gene_index=selected_genes,
        selected_gene_symbol=gene_symbols[selected_genes],
        input_mean=mean,
        input_scale=scale,
        pathway_name=pathway_names,
        cell_type_name=cell_type_names,
        time_name=time_names,
        temperature=np.asarray(temperature),
        conformal_threshold=np.asarray(conformal_threshold),
        ood_hidden_mean=hidden_mean,
        ood_hidden_scale=hidden_scale,
        ood_threshold=np.asarray(ood_threshold),
    )
    return {
        "schema": "aleph.outer_library.sciplex3_state_model.v1",
        "task": "author_annotated_drug_pathway_state_from_context_aggregated_scRNAseq",
        "model": {
            "family": f"multitask_MLP_{selected_gene_count}_{hidden_width}_pathway_celltype_time",
            "input_representation": "fit_only_celltype_replicate_time_residual_log1p_CPM",
            "nuisance_group_count": nuisance_group_count,
            "selected_gene_count": selected_gene_count,
            "pathway_class_count": len(pathway_names),
            "cell_type_class_count": len(cell_type_names),
            "time_class_count": len(time_names),
            "selection_epoch": best_epoch,
            "selection_treatment_macro_f1": best_selection_macro_f1,
            "temperature_selected_on_calibration_treatments": temperature,
        },
        "split": {
            "unit": "compound_treatment_id",
            "treatment_counts": {
                name: len(treatment_split_sets[name]) for name in SPLIT_NAMES
            },
            "cross_split_treatment_overlap_count": overlap,
            "test_treatments_used_for_fitting_selection_or_calibration": False,
            "rare_pathway_labels_are_train_only": True,
            "evaluable_pathway_labels": pathway_names[eligible_pathway].tolist(),
        },
        "sealed_unseen_treatment_test": {
            "treatment_level": treatment_metrics,
            "context_level_nonindependent_diagnostic": context_metrics,
            "treatment_ids": test_treatments.tolist(),
        },
        "conformal_90": {
            "unit": "held_out_compound_treatment",
            "calibration_treatment_count": len(calibration_treatments),
            "test_treatment_count": len(test_treatments),
            "threshold": conformal_threshold,
            "test_marginal_coverage": treatment_coverage,
            "test_mean_set_size": float(test_sets.sum(axis=1).mean()),
            "class_conditional_authority": False,
        },
        "vehicle_ood_diagnostic": {
            "vehicle_context_count": int(vehicle.sum()),
            "distance_threshold_calibrated_on_heldout_drug_contexts": ood_threshold,
            "vehicle_context_rejection_fraction": float((vehicle_distance > ood_threshold).mean()),
            "authority": "diagnostic_only_single_vehicle_treatment_same_lab",
        },
        "replicate_consistency": _replicate_consistency(hidden, context),
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "status": "training_only_unseen_compound_holdout_same_provider",
        "production_eligible": False,
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
        "limitations": [
            "pathway labels are vendor/author compound annotations, not direct single-cell state truth",
            "contexts average many cells and are the training units; cell rows are not treated as independent samples",
            "compound holdout tests chemical generalization inside one study, not a new dataset or laboratory",
            "vehicle OOD has one treatment identity and cannot establish general OOD authority",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--epochs", type=int, default=80)
    args = parser.parse_args()
    report = train(args.tensor, args.checkpoint, epochs=args.epochs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["sealed_unseen_treatment_test"]["treatment_level"], sort_keys=True))


if __name__ == "__main__":
    main()
