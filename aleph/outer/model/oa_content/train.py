#!/usr/bin/env python
"""Train same-split title-only and bounded OA JATS-content neural baselines on CPU."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import time
from collections import Counter
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
V1_SPEC = importlib.util.spec_from_file_location("external_model_v1", HERE.parent / "train.py")
v1 = importlib.util.module_from_spec(V1_SPEC)
V1_SPEC.loader.exec_module(v1)


def digest_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def aggregate_digest(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths):
        digest.update(path.name.encode())
        digest.update(bytes.fromhex(digest_file(path)))
    return digest.hexdigest()


def split_name(family: str, config: dict) -> str:
    value = int.from_bytes(hashlib.sha256(f"{config['seed']}:{family}".encode()).digest()[:8], "big") / 2**64
    if value < config["split"]["train"]:
        return "train"
    if value < config["split"]["train"] + config["split"]["validation"]:
        return "validation"
    return "test"


def load_oa_records(repo_root: Path, config: dict, compute_aggregate: bool = True) -> tuple[list[dict], dict]:
    corpus_root = repo_root / "corpus" / "external_training"
    catalog, _ = v1.load_records(corpus_root)
    metadata = {record["source_family_id"]: record for record in catalog}
    derived_root = repo_root / config["derived_root"]
    files = sorted(derived_root.glob("*.chunks.jsonl"))
    empty_files = []
    invalid_files = []
    by_family = {}
    selected_classes = tuple(config["features"]["section_token_limits"])
    coverage = Counter()
    for path in files:
        rows = []
        with path.open(encoding="utf-8") as handle:
            rows = [json.loads(line) for line in handle if line.strip()]
        if not rows:
            empty_files.append(path.name)
            continue
        families = {v1.canonical_source_family(row.get("source_family_id", "")) for row in rows}
        if len(families) != 1:
            invalid_files.append(path.name)
            continue
        family = next(iter(families))
        if family in by_family:
            invalid_files.append(path.name)
            continue
        section_tokens = {section: [] for section in selected_classes}
        for row in rows:
            section = row.get("section_class")
            if section not in section_tokens:
                continue
            if row.get("text"):
                section_tokens[section].extend(v1.TOKEN_RE.findall(row["text"].lower()))
        bounded = {}
        for section, tokens in section_tokens.items():
            limit = config["features"]["section_token_limits"][section]
            bounded[section] = tokens[:limit]
            if bounded[section]:
                coverage[section] += 1
        if family not in metadata:
            invalid_files.append(path.name)
            continue
        base = dict(metadata[family])
        base["section_tokens"] = bounded
        base["derived_file"] = path.name
        by_family[family] = base
    usable = [record for record in by_family.values() if any(record["section_tokens"].values())]
    usable.sort(key=lambda record: record["source_family_id"])
    return usable, {
        "chunk_files": len(files),
        "nonempty_unique_source_families": len(by_family),
        "usable_source_families": len(usable),
        "empty_files": empty_files,
        "invalid_or_duplicate_files": invalid_files,
        "section_family_coverage": dict(sorted(coverage.items())),
        "derived_aggregate_sha256": aggregate_digest(files) if compute_aggregate else None,
    }


def stable_bucket(token: str, width: int) -> int:
    return int.from_bytes(hashlib.sha256(token.encode()).digest()[:8], "big") % width


def feature_ids(record: dict, config: dict, include_content: bool) -> np.ndarray:
    width = config["features"]["hash_buckets"]
    minimum = config["features"]["minimum_token_length"]
    tokens = []
    for prefix, text in (("title", record["title"]), ("journal", record["journal"])):
        tokens.extend(f"{prefix}:{token}" for token in v1.TOKEN_RE.findall(text.lower()) if len(token) >= minimum)
    if record["year"]:
        bucket = record["year"] // config["features"]["year_bucket_width"] * config["features"]["year_bucket_width"]
        tokens.append(f"year:{bucket}")
    if include_content:
        for section, values in record["section_tokens"].items():
            tokens.extend(f"content:{section}:{token}" for token in values if len(token) >= minimum)
    return np.asarray(sorted({stable_bucket(token, width) for token in tokens}), dtype=np.int32)


def dense_batch(features: list[np.ndarray], indices: np.ndarray, width: int) -> np.ndarray:
    x = np.zeros((len(indices), width), dtype=np.float32)
    for row, index in enumerate(indices):
        ids = features[int(index)]
        x[row, ids] = np.float32(1.0 / math.sqrt(max(1, len(ids))))
    return x


def predict(features, indices, weights, batch_size):
    w1, b1, w2, b2 = weights
    predictions, embeddings = [], []
    for start in range(0, len(indices), batch_size):
        batch = indices[start:start + batch_size]
        x = dense_batch(features, batch, w1.shape[0])
        hidden = np.maximum(x @ w1 + b1, 0)
        predictions.append(np.argmax(hidden @ w2 + b2, axis=1))
        embeddings.append(hidden)
    return np.concatenate(predictions), np.concatenate(embeddings)


def macro_f1(y_true, y_pred, n_classes):
    scores = []
    for cls in range(n_classes):
        tp = int(np.sum((y_true == cls) & (y_pred == cls)))
        fp = int(np.sum((y_true != cls) & (y_pred == cls)))
        fn = int(np.sum((y_true == cls) & (y_pred != cls)))
        denominator = 2 * tp + fp + fn
        scores.append(0.0 if denominator == 0 else 2 * tp / denominator)
    return float(np.mean(scores)), scores


def fit(features, y, splits, labels, config, seed_offset):
    started = time.perf_counter()
    rng = np.random.default_rng(config["seed"] + seed_offset)
    width = config["features"]["hash_buckets"]
    hidden = config["model"]["hidden_width"]
    classes = len(labels)
    w1 = (rng.standard_normal((width, hidden)) * math.sqrt(2 / width)).astype(np.float32)
    b1 = np.zeros(hidden, dtype=np.float32)
    w2 = (rng.standard_normal((hidden, classes)) * math.sqrt(2 / hidden)).astype(np.float32)
    b2 = np.zeros(classes, dtype=np.float32)
    params = [w1, b1, w2, b2]
    moments = [np.zeros_like(param) for param in params]
    variances = [np.zeros_like(param) for param in params]
    best = None
    best_f1 = -1.0
    stale = step = 0
    history = []
    for epoch in range(config["model"]["epochs"]):
        order = rng.permutation(splits["train"])
        losses = []
        for start in range(0, len(order), config["model"]["batch_size"]):
            idx = order[start:start + config["model"]["batch_size"]]
            x = dense_batch(features, idx, width)
            pre = x @ w1 + b1
            h = np.maximum(pre, 0)
            logits = h @ w2 + b2
            logits -= logits.max(axis=1, keepdims=True)
            probs = np.exp(logits)
            probs /= probs.sum(axis=1, keepdims=True)
            losses.append(float(-np.log(np.maximum(probs[np.arange(len(idx)), y[idx]], 1e-12)).mean()))
            probs[np.arange(len(idx)), y[idx]] -= 1
            probs /= len(idx)
            dh = (probs @ w2.T) * (pre > 0)
            gradients = [x.T @ dh + config["model"]["l2"] * w1, dh.sum(axis=0), h.T @ probs + config["model"]["l2"] * w2, probs.sum(axis=0)]
            step += 1
            for i, (param, gradient) in enumerate(zip(params, gradients)):
                moments[i] = 0.9 * moments[i] + 0.1 * gradient
                variances[i] = 0.999 * variances[i] + 0.001 * gradient * gradient
                param -= config["model"]["learning_rate"] * (moments[i] / (1 - 0.9**step)) / (np.sqrt(variances[i] / (1 - 0.999**step)) + 1e-8)
        val_pred, _ = predict(features, splits["validation"], params, config["model"]["batch_size"])
        val_f1, _ = macro_f1(y[splits["validation"]], val_pred, classes)
        history.append({"epoch": epoch + 1, "train_cross_entropy": float(np.mean(losses)), "validation_macro_f1": val_f1})
        if val_f1 > best_f1 + 1e-6:
            best_f1 = val_f1
            best = [param.copy() for param in params]
            stale = 0
        else:
            stale += 1
            if stale >= config["model"]["patience"]:
                break
    test_pred, test_emb = predict(features, splits["test"], best, config["model"]["batch_size"])
    train_pred, train_emb = predict(features, splits["train"], best, config["model"]["batch_size"])
    test_y = y[splits["test"]]
    f1, per_class = macro_f1(test_y, test_pred, classes)
    retrieval = v1.retrieval_metrics(test_emb, test_y, train_emb, y[splits["train"]], config["retrieval"]["top_k"])
    return best, {
        "macro_f1": f1,
        "accuracy": float(np.mean(test_pred == test_y)),
        "retrieval": retrieval,
        "best_validation_macro_f1": best_f1,
        "epochs_completed": len(history),
        "history": history,
        "per_domain": {label: {"f1": per_class[i], "test_support": int(np.sum(test_y == i))} for i, label in enumerate(labels)},
        "runtime_seconds": time.perf_counter() - started,
    }


def run(config_path: Path) -> dict:
    config = json.loads(config_path.read_text())
    repo_root = config_path.resolve().parents[4]
    records, audit = load_oa_records(repo_root, config)
    labels = sorted({record["label"] for record in records})
    label_to_id = {label: index for index, label in enumerate(labels)}
    for record in records:
        record["split"] = split_name(record["source_family_id"], config)
    y = np.asarray([label_to_id[record["label"]] for record in records], dtype=np.int64)
    splits = {name: np.asarray([i for i, record in enumerate(records) if record["split"] == name], dtype=np.int64) for name in ("train", "validation", "test")}
    missing_support = {name: sorted(set(labels) - {records[i]["label"] for i in indices}) for name, indices in splits.items()}
    if any(missing_support.values()):
        raise RuntimeError(f"split_missing_label_support:{missing_support}")
    title_features = [feature_ids(record, config, False) for record in records]
    content_features = [feature_ids(record, config, True) for record in records]
    title_weights, title_metrics = fit(title_features, y, splits, labels, config, 0)
    content_weights, content_metrics = fit(content_features, y, splits, labels, config, 1)
    output = config_path.parent
    np.savez_compressed(output / "weights_title_only.npz", w1=title_weights[0], b1=title_weights[1], w2=title_weights[2], b2=title_weights[3], labels=np.asarray(labels))
    np.savez_compressed(output / "weights_oa_content.npz", w1=content_weights[0], b1=content_weights[1], w2=content_weights[2], b2=content_weights[3], labels=np.asarray(labels))
    with (output / "splits.jsonl").open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps({"source_family_id": record["source_family_id"], "derived_file": record["derived_file"], "label": record["label"], "split": record["split"]}, sort_keys=True) + "\n")
    metrics = {
        "schema_version": "1.0.0",
        "authority": config["authority"],
        "limitations": ["Discovery domains are weak labels, not reviewed biology.", "OA text is proposed evidence and cannot train physical parameters or promote authority.", "The 2,701 chunk-file count contains one empty file; only 2,695 families have bounded target-section text."],
        "audit": audit,
        "splits": {name: len(indices) for name, indices in splits.items()},
        "label_count": len(labels),
        "label_counts_by_split": {name: dict(sorted(Counter(records[i]["label"] for i in indices).items())) for name, indices in splits.items()},
        "title_only": title_metrics,
        "oa_content_plus_metadata": content_metrics,
        "delta": {"macro_f1": content_metrics["macro_f1"] - title_metrics["macro_f1"], "accuracy": content_metrics["accuracy"] - title_metrics["accuracy"], "retrieval_hit_rate_at_5": content_metrics["retrieval"]["hit_rate_at_k"] - title_metrics["retrieval"]["hit_rate_at_k"], "retrieval_precision_at_5": content_metrics["retrieval"]["precision_at_k"] - title_metrics["retrieval"]["precision_at_k"]},
    }
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    receipt_names = ["config.json", "train.py", "splits.jsonl", "metrics.json", "weights_title_only.npz", "weights_oa_content.npz"]
    (output / "SHA256SUMS.json").write_text(json.dumps({name: digest_file(output / name) for name in receipt_names}, indent=2, sort_keys=True) + "\n")
    return metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    metrics = run(args.config.resolve())
    print(json.dumps({"usable": metrics["audit"]["usable_source_families"], "splits": metrics["splits"], "title_macro_f1": metrics["title_only"]["macro_f1"], "content_macro_f1": metrics["oa_content_plus_metadata"]["macro_f1"], "delta": metrics["delta"]}, indent=2))


if __name__ == "__main__":
    main()
