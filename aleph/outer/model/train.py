#!/usr/bin/env python
"""CPU-only weak-supervision metadata tagger/retriever.

This is an external catalogue experiment, not Aleph learning or physics training.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


TOKEN_RE = re.compile(r"[a-z0-9]+")


def canonical_source_family(source_family_id: str) -> str:
    """Canonicalise external identifier namespaces before deduplication and splitting."""
    namespace, separator, identifier = source_family_id.strip().partition(":")
    if not separator:
        return source_family_id.strip().lower()
    return f"{namespace.lower()}:{identifier.lower()}"


def canonical_domain(domain: str) -> str:
    if domain.startswith("2010s_"):
        domain = domain[6:]
    if domain == "stem_cell_differentiation_states":
        return "stem_cell_states"
    return domain


HISTORICAL_RULES = (
    ("nucleus_chromatin", ("nucle", "chromatin", "lamin")),
    ("adhesion_traction", ("adhesion", "traction", "integrin", "focal adhesion")),
    ("cortex_membrane_pressure", ("cortex", "cortical", "membrane tension", "bleb", "surface tension")),
    ("osmotic_poroelastic_volume", ("osmotic", "cell volume", "poroelastic", "hydraulic")),
    ("ecm_invasion", ("extracellular matrix", "matrix stiffness", "invasion", "collagen")),
    ("protrusion_migration", ("migration", "motility", "lamell", "filopod")),
    ("cell_division_cytokinesis", ("mitosis", "cytokinesis", "cell division")),
    ("cytoskeleton_networks", ("actin", "microtub", "cytoskeleton", "myosin")),
    ("mechanotransduction_signalling", ("mechanotrans", "yap", "taz", "force signaling")),
    ("mechanical_measurement_modalities", ("atomic force", "micropipette", "optical tweezer", "rheolog")),
    ("fluid_shear_endothelium", ("fluid shear", "shear stress", "endothel")),
    ("developmental_morphogenesis", ("morphogenesis", "development", "embry")),
    ("cancer_invasion_mechanics", ("cancer", "tumor", "tumour", "metasta")),
)


def historical_domain(record: dict) -> str:
    text = (record.get("title") or "").lower()
    for domain, needles in HISTORICAL_RULES:
        if any(needle in text for needle in needles):
            return domain
    return "general_cell_biology"


def seed_domain(record: dict) -> str:
    tags = " ".join(record.get("topic_tags", ())).lower()
    title = (record.get("title") or "").lower()
    merged = f"{tags} {title}"
    for domain, needles in HISTORICAL_RULES:
        if any(needle in merged for needle in needles):
            return domain
    if any(x in merged for x in ("simulation_based_inference", "network", "gnn", "neural", "cell_painting", "multimodal")):
        return "observation_inference_methods"
    return "general_cell_biology"


def load_records(root: Path) -> tuple[list[dict], dict]:
    snapshots = sorted((root / "snapshots").glob("*.jsonl"))
    rows: dict[str, dict] = {}
    duplicate_rows = 0
    for path in snapshots:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                record = json.loads(line)
                family = canonical_source_family(record["source_family_id"])
                if family in rows:
                    duplicate_rows += 1
                    continue
                domain = record["discovery_domain"]
                label = historical_domain(record) if domain.startswith("historical_") else canonical_domain(domain)
                rows[family] = {
                    "source_family_id": family,
                    "title": record.get("title") or "",
                    "journal": record.get("journal") or "",
                    "year": int(record.get("year") or 0),
                    "label": label,
                    "label_source": "historical_title_rule" if domain.startswith("historical_") else "discovery_query",
                }
    seed_path = root / "sources" / "seed_primary_sources.jsonl"
    seed_overlaps = 0
    with seed_path.open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            family = canonical_source_family(record["source_family_id"])
            seed_overlaps += int(family in rows)
            rows[family] = {
                "source_family_id": family,
                "title": record.get("title") or "",
                "journal": "",
                "year": int(record.get("year") or 0),
                "label": seed_domain(record),
                "label_source": "seed_topic_rule",
            }
    ordered = [rows[key] for key in sorted(rows)]
    def line_count(path: Path) -> int:
        with path.open(encoding="utf-8") as handle:
            return sum(1 for _ in handle)

    return ordered, {
        "snapshot_files": [str(p.relative_to(root.parent.parent.parent)) for p in snapshots],
        "snapshot_rows": sum(line_count(p) for p in snapshots),
        "duplicate_snapshot_rows": duplicate_rows,
        "seed_rows": line_count(seed_path),
        "seed_overlaps": seed_overlaps,
        "unique_source_families": len(ordered),
    }


def split_name(family: str, seed: int, train: float, validation: float) -> str:
    digest = hashlib.sha256(f"{seed}:{family}".encode()).digest()
    value = int.from_bytes(digest[:8], "big") / 2**64
    if value < train:
        return "train"
    if value < train + validation:
        return "validation"
    return "test"


def token_strings(record: dict, cfg: dict) -> list[str]:
    minimum = cfg["minimum_token_length"]
    result = []
    for prefix, text in ((cfg["title_prefix"], record["title"]), (cfg["journal_prefix"], record["journal"])):
        result.extend(f"{prefix}:{token}" for token in TOKEN_RE.findall(text.lower()) if len(token) >= minimum)
    width = cfg["year_bucket_width"]
    if record["year"]:
        result.append(f"year:{record['year'] // width * width}")
    return sorted(set(result))


def feature_ids(record: dict, cfg: dict) -> np.ndarray:
    buckets = cfg["hash_buckets"]
    return np.asarray(sorted({int.from_bytes(hashlib.sha256(t.encode()).digest()[:8], "big") % buckets for t in token_strings(record, cfg)}), dtype=np.int32)


def dense_batch(features: list[np.ndarray], indices: np.ndarray, width: int) -> np.ndarray:
    x = np.zeros((len(indices), width), dtype=np.float32)
    for row, index in enumerate(indices):
        ids = features[int(index)]
        if ids.size:
            x[row, ids] = np.float32(1.0 / math.sqrt(ids.size))
    return x


def macro_f1(y_true: np.ndarray, y_pred: np.ndarray, n_classes: int) -> tuple[float, list[float]]:
    values = []
    for cls in range(n_classes):
        tp = int(np.sum((y_true == cls) & (y_pred == cls)))
        fp = int(np.sum((y_true != cls) & (y_pred == cls)))
        fn = int(np.sum((y_true == cls) & (y_pred != cls)))
        denom = 2 * tp + fp + fn
        values.append(0.0 if denom == 0 else 2 * tp / denom)
    return float(np.mean(values)), values


def predict(features: list[np.ndarray], indices: np.ndarray, w1: np.ndarray, b1: np.ndarray, w2: np.ndarray, b2: np.ndarray, batch_size: int) -> tuple[np.ndarray, np.ndarray]:
    outputs, embeddings = [], []
    for start in range(0, len(indices), batch_size):
        batch = indices[start : start + batch_size]
        x = dense_batch(features, batch, w1.shape[0])
        hidden = np.maximum(x @ w1 + b1, 0)
        outputs.append(np.argmax(hidden @ w2 + b2, axis=1))
        embeddings.append(hidden)
    return np.concatenate(outputs), np.concatenate(embeddings)


def retrieval_metrics(query_emb: np.ndarray, query_y: np.ndarray, corpus_emb: np.ndarray, corpus_y: np.ndarray, k: int) -> dict:
    q = query_emb / np.maximum(np.linalg.norm(query_emb, axis=1, keepdims=True), 1e-12)
    c = corpus_emb / np.maximum(np.linalg.norm(corpus_emb, axis=1, keepdims=True), 1e-12)
    hits = precisions = 0.0
    for start in range(0, len(q), 128):
        scores = q[start : start + 128] @ c.T
        top = np.argpartition(scores, -k, axis=1)[:, -k:]
        relevant = corpus_y[top] == query_y[start : start + len(top), None]
        hits += float(np.sum(np.any(relevant, axis=1)))
        precisions += float(np.sum(np.mean(relevant, axis=1)))
    return {"hit_rate_at_k": hits / len(q), "precision_at_k": precisions / len(q), "queries": len(q), "k": k}


def train(config_path: Path) -> dict:
    started = time.perf_counter()
    config = json.loads(config_path.read_text())
    model_dir = config_path.parent
    corpus_root = model_dir.parent
    records, input_summary = load_records(corpus_root)
    labels = sorted({r["label"] for r in records})
    label_to_id = {label: i for i, label in enumerate(labels)}
    for r in records:
        r["split"] = split_name(r["source_family_id"], config["seed"], config["split"]["train"], config["split"]["validation"])
    features = [feature_ids(r, config["features"]) for r in records]
    y = np.asarray([label_to_id[r["label"]] for r in records], dtype=np.int64)
    splits = {name: np.asarray([i for i, r in enumerate(records) if r["split"] == name], dtype=np.int64) for name in ("train", "validation", "test")}

    rng = np.random.default_rng(config["seed"])
    width = config["features"]["hash_buckets"]
    hidden_width = config["model"]["hidden_width"]
    n_classes = len(labels)
    w1 = (rng.standard_normal((width, hidden_width)) * math.sqrt(2 / width)).astype(np.float32)
    b1 = np.zeros(hidden_width, dtype=np.float32)
    w2 = (rng.standard_normal((hidden_width, n_classes)) * math.sqrt(2 / hidden_width)).astype(np.float32)
    b2 = np.zeros(n_classes, dtype=np.float32)
    params = [w1, b1, w2, b2]
    m = [np.zeros_like(p) for p in params]
    v = [np.zeros_like(p) for p in params]
    best = None
    best_val = -1.0
    stale = 0
    step = 0
    history = []
    batch_size = config["model"]["batch_size"]
    lr = config["model"]["learning_rate"]
    l2 = config["model"]["l2"]
    for epoch in range(config["model"]["epochs"]):
        order = rng.permutation(splits["train"])
        losses = []
        for start in range(0, len(order), batch_size):
            idx = order[start : start + batch_size]
            x = dense_batch(features, idx, width)
            hpre = x @ w1 + b1
            h = np.maximum(hpre, 0)
            logits = h @ w2 + b2
            logits -= logits.max(axis=1, keepdims=True)
            probs = np.exp(logits)
            probs /= probs.sum(axis=1, keepdims=True)
            losses.append(float(-np.log(np.maximum(probs[np.arange(len(idx)), y[idx]], 1e-12)).mean()))
            probs[np.arange(len(idx)), y[idx]] -= 1
            probs /= len(idx)
            grads = [x.T @ ((probs @ w2.T) * (hpre > 0)) + l2 * w1,
                     np.sum((probs @ w2.T) * (hpre > 0), axis=0),
                     h.T @ probs + l2 * w2,
                     np.sum(probs, axis=0)]
            step += 1
            for j, (param, grad) in enumerate(zip(params, grads)):
                m[j] = 0.9 * m[j] + 0.1 * grad
                v[j] = 0.999 * v[j] + 0.001 * grad * grad
                param -= lr * (m[j] / (1 - 0.9**step)) / (np.sqrt(v[j] / (1 - 0.999**step)) + 1e-8)
        val_pred, _ = predict(features, splits["validation"], w1, b1, w2, b2, batch_size)
        val_f1, _ = macro_f1(y[splits["validation"]], val_pred, n_classes)
        history.append({"epoch": epoch + 1, "train_cross_entropy": float(np.mean(losses)), "validation_macro_f1": val_f1})
        if val_f1 > best_val + 1e-6:
            best_val = val_f1
            best = [p.copy() for p in params]
            stale = 0
        else:
            stale += 1
            if stale >= config["model"]["patience"]:
                break
    w1, b1, w2, b2 = best

    train_pred, train_emb = predict(features, splits["train"], w1, b1, w2, b2, batch_size)
    test_pred, test_emb = predict(features, splits["test"], w1, b1, w2, b2, batch_size)
    test_y = y[splits["test"]]
    test_f1, per_class = macro_f1(test_y, test_pred, n_classes)
    majority = Counter(y[splits["train"]]).most_common(1)[0][0]
    majority_pred = np.full_like(test_y, majority)
    majority_f1, _ = macro_f1(test_y, majority_pred, n_classes)

    # Simple non-neural nearest-centroid baseline in the same fixed hashed metadata space.
    centroids = np.zeros((n_classes, width), dtype=np.float32)
    for cls in range(n_classes):
        members = splits["train"][y[splits["train"]] == cls]
        for start in range(0, len(members), batch_size):
            centroids[cls] += dense_batch(features, members[start : start + batch_size], width).sum(axis=0)
        centroids[cls] /= max(1, len(members))
    centroids /= np.maximum(np.linalg.norm(centroids, axis=1, keepdims=True), 1e-12)
    centroid_pred = []
    test_fixed = []
    for start in range(0, len(splits["test"]), batch_size):
        xb = dense_batch(features, splits["test"][start : start + batch_size], width)
        test_fixed.append(xb)
        centroid_pred.append(np.argmax(xb @ centroids.T, axis=1))
    centroid_pred = np.concatenate(centroid_pred)
    centroid_f1, _ = macro_f1(test_y, centroid_pred, n_classes)

    query_limit = config["retrieval"]["max_test_queries"]
    query_count = min(query_limit, len(test_emb))
    query_sel = np.arange(query_count)
    neural_retrieval = retrieval_metrics(test_emb[query_sel], test_y[query_sel], train_emb, y[splits["train"]], config["retrieval"]["top_k"])
    fixed_train = np.concatenate([dense_batch(features, splits["train"][s:s+batch_size], width) for s in range(0, len(splits["train"]), batch_size)])
    fixed_test = np.concatenate(test_fixed)
    fixed_retrieval = retrieval_metrics(fixed_test[query_sel], test_y[query_sel], fixed_train, y[splits["train"]], config["retrieval"]["top_k"])

    weights_path = model_dir / "weights.npz"
    np.savez_compressed(weights_path, w1=w1, b1=b1, w2=w2, b2=b2, labels=np.asarray(labels))
    split_path = model_dir / "splits.jsonl"
    with split_path.open("w", encoding="utf-8") as handle:
        for r in records:
            handle.write(json.dumps({k: r[k] for k in ("source_family_id", "label", "label_source", "split")}, sort_keys=True) + "\n")
    label_counts = {split: dict(sorted(Counter(records[i]["label"] for i in indices).items())) for split, indices in splits.items()}
    metrics = {
        "schema_version": "1.0.0",
        "authority": config["authority"],
        "limitations": [
            "Labels are weak discovery/query provenance and title rules, not reviewed biological ground truth.",
            "Metadata-only metrics do not validate article quality, mechanism, trajectories, or physical parameters.",
            "This external model does not unseal Aleph MTG-PN and cannot satisfy R5."
        ],
        "input": input_summary,
        "splits": {name: int(len(indices)) for name, indices in splits.items()},
        "label_count": n_classes,
        "label_counts_by_split": label_counts,
        "training": {"epochs_completed": len(history), "best_validation_macro_f1": best_val, "history": history},
        "tagging": {
            "majority_macro_f1": majority_f1,
            "nearest_centroid_macro_f1": centroid_f1,
            "neural_macro_f1": test_f1,
            "neural_accuracy": float(np.mean(test_pred == test_y)),
            "per_domain": {label: {"f1": per_class[i], "test_support": int(np.sum(test_y == i))} for i, label in enumerate(labels)}
        },
        "retrieval": {"fixed_hashed_baseline": fixed_retrieval, "neural_hidden": neural_retrieval},
        "runtime": {"device": "cpu", "numpy_version": np.__version__, "seconds": time.perf_counter() - started},
    }
    metrics_path = model_dir / "metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    receipt_targets = [config_path, Path(__file__), weights_path, split_path, metrics_path]
    receipts = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in receipt_targets}
    (model_dir / "SHA256SUMS.json").write_text(json.dumps(receipts, indent=2, sort_keys=True) + "\n")
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    metrics = train(args.config.resolve())
    print(json.dumps({"families": metrics["input"]["unique_source_families"], "splits": metrics["splits"], "neural_macro_f1": metrics["tagging"]["neural_macro_f1"], "neural_retrieval": metrics["retrieval"]["neural_hidden"], "seconds": metrics["runtime"]["seconds"]}, indent=2))


if __name__ == "__main__":
    main()
