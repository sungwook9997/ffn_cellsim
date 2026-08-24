#!/usr/bin/env python3
"""Train protocol-routing and magnetic-wire mechanics heads from external data.

The two heads have deliberately different authority.  The protocol router maps
measurement metadata to a non-pooling observation family.  The wire head learns
only within Dessard 2024 and is evaluated both on random held-out wires and on
entire experimental-date holdouts.  Neither is an Aleph parameter predictor.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

CELL_LINES = np.asarray(["MCF-10A", "MCF-7", "MDA-MB-231"], dtype="U16")
WIRE_FEATURES = ("wire_length_um", "wire_diameter_um", "omega_c_rad_s", "wire_aspect_effective")


def _softmax(logits: np.ndarray) -> np.ndarray:
    shifted = logits - logits.max(axis=1, keepdims=True)
    values = np.exp(shifted)
    return values / values.sum(axis=1, keepdims=True)


def _metrics(probability: np.ndarray, label: np.ndarray, names: np.ndarray) -> dict[str, Any]:
    prediction = probability.argmax(axis=1)
    confusion = np.zeros((len(names), len(names)), dtype=np.int64)
    np.add.at(confusion, (label, prediction), 1)
    recalls = []
    f1s = []
    per_class = {}
    for index, name in enumerate(names):
        tp = int(confusion[index, index])
        fp = int(confusion[:, index].sum() - tp)
        fn = int(confusion[index].sum() - tp)
        recall = tp / max(1, tp + fn)
        f1 = 2 * tp / max(1, 2 * tp + fp + fn)
        recalls.append(recall)
        f1s.append(f1)
        per_class[str(name)] = {"recall": recall, "f1": f1, "support": int(confusion[index].sum())}
    return {
        "sample_count": int(len(label)),
        "accuracy": float(np.mean(prediction == label)),
        "balanced_accuracy": float(np.mean(recalls)),
        "macro_f1": float(np.mean(f1s)),
        "per_class": per_class,
        "confusion_matrix_true_by_predicted": confusion.tolist(),
    }


def _wire_roles(labels: np.ndarray, groups: np.ndarray) -> np.ndarray:
    roles = np.empty(len(labels), dtype="U12")
    for label in range(len(CELL_LINES)):
        indices = np.flatnonzero(labels == label)
        ordered = sorted(set(groups[indices].tolist()))
        if len(ordered) < 3:
            raise ValueError(f"{CELL_LINES[label]} lacks three experimental dates")
        roles[indices] = "fit"
        roles[indices[groups[indices] == ordered[-2]]] = "calibration"
        roles[indices[groups[indices] == ordered[-1]]] = "test"
    return roles


def _random_roles(labels: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    roles = np.full(len(labels), "fit", dtype="U12")
    for label in range(len(CELL_LINES)):
        indices = np.flatnonzero(labels == label)
        rng.shuffle(indices)
        n_test = max(1, int(round(0.2 * len(indices))))
        n_calibration = max(1, int(round(0.2 * len(indices))))
        roles[indices[:n_test]] = "test"
        roles[indices[n_test:n_test + n_calibration]] = "calibration"
    return roles


def _fit_wire_model(
    raw_x: np.ndarray,
    log_eta: np.ndarray,
    label: np.ndarray,
    roles: np.ndarray,
    seed: int,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    fit = roles == "fit"
    calibration = roles == "calibration"
    test = roles == "test"
    x = np.log(raw_x)
    mean = x[fit].mean(axis=0)
    scale = x[fit].std(axis=0)
    scale[scale < 1e-8] = 1.0
    x = (x - mean) / scale
    eta_mean = float(log_eta[fit].mean())
    eta_scale = float(log_eta[fit].std())
    eta_target = (log_eta - eta_mean) / eta_scale
    one_hot = np.eye(len(CELL_LINES))[label]
    rng = np.random.default_rng(seed)
    params = {
        "W_shared": rng.normal(0.0, 0.18, (x.shape[1], 16)),
        "b_shared": np.zeros(16),
        "W_class": rng.normal(0.0, 0.18, (16, len(CELL_LINES))),
        "b_class": np.zeros(len(CELL_LINES)),
        "W_eta": rng.normal(0.0, 0.18, (16, 1)),
        "b_eta": np.zeros(1),
    }
    for epoch in range(6_000):
        hidden = np.tanh(x @ params["W_shared"] + params["b_shared"])
        probability = _softmax(hidden @ params["W_class"] + params["b_class"])
        eta_prediction = (hidden @ params["W_eta"] + params["b_eta"]).ravel()
        dz = (probability - one_hot) / fit.sum()
        dz[~fit] = 0.0
        deta = 0.35 * (eta_prediction - eta_target) / fit.sum()
        deta[~fit] = 0.0
        d_w_class = hidden.T @ dz + 0.003 * params["W_class"]
        d_b_class = dz.sum(axis=0)
        d_w_eta = hidden.T @ deta[:, None] + 0.003 * params["W_eta"]
        d_b_eta = np.asarray([deta.sum()])
        d_hidden = dz @ params["W_class"].T + deta[:, None] @ params["W_eta"].T
        d_pre = d_hidden * (1.0 - hidden * hidden)
        d_w_shared = x.T @ d_pre + 0.003 * params["W_shared"]
        d_b_shared = d_pre.sum(axis=0)
        learning_rate = 0.025 if epoch < 3_000 else 0.01
        for name, gradient in (
            ("W_shared", d_w_shared), ("b_shared", d_b_shared),
            ("W_class", d_w_class), ("b_class", d_b_class),
            ("W_eta", d_w_eta), ("b_eta", d_b_eta),
        ):
            params[name] -= learning_rate * gradient

    hidden = np.tanh(x @ params["W_shared"] + params["b_shared"])
    probability = _softmax(hidden @ params["W_class"] + params["b_class"])
    predicted_log_eta = (hidden @ params["W_eta"] + params["b_eta"]).ravel() * eta_scale + eta_mean
    split_metrics = {}
    for name, mask in (("fit", fit), ("calibration", calibration), ("test", test)):
        residual = predicted_log_eta[mask] - log_eta[mask]
        denominator = np.sum((log_eta[mask] - log_eta[mask].mean()) ** 2)
        split_metrics[name] = {
            "cell_line": _metrics(probability[mask], label[mask], CELL_LINES),
            "log_viscosity": {
                "mae_log_pa_s": float(np.mean(np.abs(residual))),
                "median_absolute_percent_error": float(np.median(np.abs(np.exp(residual) - 1.0))),
                "r2": float(1.0 - np.sum(residual ** 2) / denominator) if denominator else None,
            },
        }
    checkpoint = {
        **params,
        "input_mean": mean,
        "input_scale": scale,
        "eta_log_mean": np.asarray(eta_mean),
        "eta_log_scale": np.asarray(eta_scale),
        "feature_name": np.asarray(WIRE_FEATURES),
        "class_name": CELL_LINES,
    }
    return checkpoint, split_metrics


def _ols_slope(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    design = np.column_stack([np.ones(len(x)), x])
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    prediction = design @ coefficients
    r2 = 1.0 - float(np.sum((y - prediction) ** 2) / np.sum((y - y.mean()) ** 2))
    return float(coefficients[1]), r2


def _probe_scaling(
    length: np.ndarray, eta: np.ndarray, labels: np.ndarray, groups: np.ndarray
) -> dict[str, Any]:
    rng = np.random.default_rng(10392024)
    results = {}
    for label, name in enumerate(CELL_LINES):
        mask = labels == label
        x = np.log(length[mask])
        y = np.log(eta[mask])
        group = groups[mask]
        slope, r2 = _ols_slope(x, y)
        unique_groups = sorted(set(group.tolist()))
        boot = []
        for _ in range(4_000):
            selected = rng.choice(unique_groups, size=len(unique_groups), replace=True)
            indices = np.concatenate([np.flatnonzero(group == value) for value in selected])
            if len(set(x[indices].round(6))) > 1:
                boot.append(_ols_slope(x[indices], y[indices])[0])
        interval = np.quantile(boot, [0.025, 0.975])
        results[str(name)] = {
            "n_wires": int(mask.sum()),
            "n_experiment_dates": len(unique_groups),
            "raw_wire_level_log_log_exponent": slope,
            "experiment_date_cluster_bootstrap_95_interval": interval.tolist(),
            "raw_wire_level_r2": r2,
            "predicted_eta_fold_1_to_5_um_from_raw_fit": float(5.0 ** slope),
            "published_binned_exponent": {"MCF-10A": 2.2, "MCF-7": 2.2, "MDA-MB-231": 2.5}[str(name)],
            "published_standard_error": 0.2,
        }
    return results


def _tokenize(record: dict[str, Any]) -> list[str]:
    fields = ("probe_geometry", "observable", "unit")
    tokens = []
    for field in fields:
        tokens.extend(f"{field}:{token}" for token in re.findall(r"[a-z0-9]+", str(record[field]).lower()))
    state = str(record["cell_state"])
    tokens.append("surface:suspended" if state.startswith("suspended") else "surface:adherent")
    tokens.append("time:unknown" if record["time_scale_s"] is None else (
        "time:subsecond" if record["time_scale_s"] < 1 else "time:seconds_or_longer"
    ))
    return sorted(set(tokens))


def _method_superclass(method: str) -> str:
    if method in {"magnetic_rotational_spectroscopy", "magnetic_bead_microrheometry"}:
        return "internal_or_local_magnetic_rheology"
    if method == "optical_tweezers_brownian_probe":
        return "optical_probe_microrheology"
    if method == "continuous_micropipette_aspiration":
        return "whole_cell_suction"
    if method in {"dynamic_afm_parallel_plate", "afm_parallel_plate_equilibrium"}:
        return "whole_cell_afm_confinement"
    if method == "electrodeformation_relaxation":
        return "noncontact_whole_cell_deformation"
    return "local_afm"


def _protocol_router(catalog: dict[str, Any]) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    records = catalog["records"]
    vocabulary = np.asarray(sorted({token for record in records for token in _tokenize(record)}), dtype="U80")
    token_index = {token: index for index, token in enumerate(vocabulary)}
    x = np.zeros((len(records), len(vocabulary) + 2), dtype=np.float64)
    for row, record in enumerate(records):
        for token in _tokenize(record):
            x[row, token_index[token]] = 1.0
        x[row, -2] = math.log10(float(record["probe_size_um"]) + 1e-6)
        x[row, -1] = -6.0 if record["time_scale_s"] is None else math.log10(float(record["time_scale_s"]) + 1e-6)
    method_text = np.asarray([_method_superclass(record["method_family"]) for record in records])
    scope_text = np.asarray([record["mechanical_scope"] for record in records])
    state_text = np.asarray(["suspended" if record["cell_state"].startswith("suspended") else "adherent" for record in records])
    heads = {"method": method_text, "scope": scope_text, "surface_state": state_text}
    encoded = {}
    names = {}
    for head, values in heads.items():
        names[head] = np.asarray(sorted(set(values.tolist())), dtype="U64")
        lookup = {name: index for index, name in enumerate(names[head])}
        encoded[head] = np.asarray([lookup[value] for value in values], dtype=np.int64)

    rng = np.random.default_rng(20260810)
    width = 24
    params: dict[str, np.ndarray] = {
        "W_shared": rng.normal(0.0, 0.15, (x.shape[1], width)), "b_shared": np.zeros(width)
    }
    for head in heads:
        params[f"W_{head}"] = rng.normal(0.0, 0.15, (width, len(names[head])))
        params[f"b_{head}"] = np.zeros(len(names[head]))
    for _epoch in range(5_000):
        hidden = np.tanh(x @ params["W_shared"] + params["b_shared"])
        d_hidden = np.zeros_like(hidden)
        for head in heads:
            probability = _softmax(hidden @ params[f"W_{head}"] + params[f"b_{head}"])
            dz = (probability - np.eye(len(names[head]))[encoded[head]]) / len(records)
            d_w = hidden.T @ dz + 0.005 * params[f"W_{head}"]
            d_b = dz.sum(axis=0)
            d_hidden += dz @ params[f"W_{head}"].T
            params[f"W_{head}"] -= 0.025 * d_w
            params[f"b_{head}"] -= 0.025 * d_b
        d_pre = d_hidden * (1.0 - hidden * hidden)
        params["W_shared"] -= 0.025 * (x.T @ d_pre + 0.005 * params["W_shared"])
        params["b_shared"] -= 0.025 * d_pre.sum(axis=0)
    hidden = np.tanh(x @ params["W_shared"] + params["b_shared"])
    training_metrics = {
        head: _metrics(_softmax(hidden @ params[f"W_{head}"] + params[f"b_{head}"]), encoded[head], names[head])
        for head in heads
    }
    source_counts = Counter(record["source_id"] for record in records)
    class_source_counts = {
        head: {
            str(name): len({records[i]["source_id"] for i in np.flatnonzero(values == name)})
            for name in names[head]
        }
        for head, values in heads.items()
    }
    checkpoint = {
        **params, "vocabulary": vocabulary,
        **{f"class_name_{head}": value for head, value in names.items()},
    }
    report = {
        "record_count": len(records),
        "source_count": len(source_counts),
        "source_record_counts": dict(sorted(source_counts.items())),
        "training_resubstitution_diagnostic": training_metrics,
        "independent_source_counts_by_class": class_source_counts,
        "external_holdout": False,
        "status": "routing_representation_only_insufficient_external_holdout",
        "reason": "several method/scope classes occur in only one independent paper; training fit is not a generalization metric",
        "runtime_policy": "exact ontology route or abstain; no interpolation across mechanical_scope or unit",
    }
    return checkpoint, report


def _svg_confusion(path: Path, confusion: list[list[int]], title: str) -> None:
    matrix = np.asarray(confusion)
    maximum = max(1, int(matrix.max()))
    cells = []
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = int(matrix[row, column])
            opacity = 0.12 + 0.78 * value / maximum
            cells.append(
                f'<rect x="{160 + 90 * column}" y="{70 + 70 * row}" width="86" height="66" fill="#2563eb" opacity="{opacity:.3f}"/>'
                f'<text x="{203 + 90 * column}" y="{110 + 70 * row}" text-anchor="middle" font-size="20">{value}</text>'
            )
    labels = "".join(
        f'<text x="{203 + 90 * i}" y="{300}" text-anchor="middle" font-size="13">{name}</text>'
        f'<text x="150" y="{110 + 70 * i}" text-anchor="end" font-size="13">{name}</text>'
        for i, name in enumerate(CELL_LINES)
    )
    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="520" height="340" viewBox="0 0 520 340">'
        '<rect width="100%" height="100%" fill="white"/>'
        f'<text x="260" y="30" text-anchor="middle" font-size="18" font-weight="bold">{title}</text>'
        + "".join(cells) + labels
        + '<text x="300" y="330" text-anchor="middle" font-size="13">predicted</text>'
        + '<text x="28" y="175" text-anchor="middle" font-size="13" transform="rotate(-90 28 175)">true</text></svg>'
    )
    path.write_text(svg, encoding="utf-8")


def train(tensor_path: Path, catalog_path: Path, output_dir: Path) -> dict[str, Any]:
    with np.load(tensor_path, allow_pickle=False) as data:
        feature_names = np.asarray(data["feature_name"]).astype(str)
        measurement = np.asarray(data["measurement"], dtype=np.float64)
        labels = np.asarray(data["cell_label"], dtype=np.int64)
        groups = np.asarray(data["experiment_group"]).astype(str)
    if measurement.shape != (196, 10) or feature_names[:5].tolist() != [
        "wire_length_um", "wire_diameter_um", "omega_c_rad_s", "wire_aspect_effective", "viscosity_pa_s"
    ]:
        raise ValueError("Dessard normalized tensor contract changed")
    raw_x = measurement[:, :4]
    log_eta = np.log(measurement[:, 4])
    run_roles = _wire_roles(labels, groups)
    random_roles = _random_roles(labels, 20260810)
    run_checkpoint, run_metrics = _fit_wire_model(raw_x, log_eta, labels, run_roles, 20260810)
    _, random_metrics = _fit_wire_model(raw_x, log_eta, labels, random_roles, 20260811)
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    router_checkpoint, router_report = _protocol_router(catalog)
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_dir / "dessard_wire_multihead.npz", **run_checkpoint)
    np.savez_compressed(output_dir / "mechanics_protocol_router.npz", **router_checkpoint)
    probe_scaling = _probe_scaling(measurement[:, 0], measurement[:, 4], labels, groups)
    test_confusion = run_metrics["test"]["cell_line"]["confusion_matrix_true_by_predicted"]
    _svg_confusion(output_dir / "run_disjoint_cell_line_confusion.svg", test_confusion, "Experimental-date holdout")
    report = {
        "schema": "aleph.outer_library.mechanics_protocol_model.v1",
        "task": "protocol_conditioned_cell_mechanics_representation",
        "dataset": {
            "individual_wire_source": "Dessard2024_NanoscaleAdv",
            "wire_count": 196,
            "experiment_date_count": len(set(groups.tolist())),
            "catalog_record_count": len(catalog["records"]),
            "catalog_source_count": len(set(record["source_id"] for record in catalog["records"])),
            "repo_source_audit_counts": dict(sorted(Counter(
                source["repo_source_audit"] for source in catalog["sources"].values()
            ).items())),
        },
        "architecture": {
            "wire_head": "4_input_16_tanh_shared_3_class_plus_log_viscosity_regression",
            "protocol_router": "hashed_protocol_tokens_24_tanh_three_heads",
            "device": "cpu",
            "gpu_used": False,
        },
        "wire_head": {
            "input_features": list(WIRE_FEATURES),
            "random_wire_holdout_interpolation_only": random_metrics,
            "experimental_date_disjoint": {
                "split_group_overlap_count": 0,
                "roles_by_cell_line": {
                    str(name): {
                        role: sorted(set(groups[(labels == index) & (run_roles == role)].tolist()))
                        for role in ("fit", "calibration", "test")
                    }
                    for index, name in enumerate(CELL_LINES)
                },
                "metrics": run_metrics,
            },
            "status": "within_paper_representation_only_failed_robust_cell_line_transfer",
            "independent_lab_holdout": False,
            "production_eligible": False,
        },
        "protocol_router": router_report,
        "probe_size_dependence": {
            "external_mrs_wire_length_scaling": probe_scaling,
            "decisive_afm_probe_radius_sweep": {
                "poroelastic_prediction": "tau proportional to a^2",
                "fold_change_a_1_to_5_um": 25.0,
                "probe_independent_clock_prediction": "tau proportional to a^0",
                "fold_change_a_1_to_5_um_probe_independent": 1.0,
                "status": "preregistered_discriminator_not_yet_an_external_afm_radius_dataset",
                "warning": "Dessard wire length scaling is evidence that apparent viscosity is probe-scale dependent; wire length is not AFM contact radius and does not validate tau(a).",
            },
        },
        "capabilities": [
            "route known protocols to non-poolable measurement and mechanical-scope heads",
            "learn the deterministic dynamic-response-to-effective-viscosity mapping from individual wires",
            "quantify how a cell-line classifier degrades when an entire acquisition date is withheld",
            "attach a probe-scale hypothesis to a future relaxation curve without choosing an Aleph parameter",
        ],
        "refusals": [
            "no cross-unit global mechanics score",
            "no claim of general cell-state classification from one-lab magnetic-wire data",
            "no substitution of MRS wire length for AFM contact radius",
            "no engine validation, gate verdict, parameter selection, or physics mutation",
        ],
        "observable_evidence_eligible": False,
        "may_select_aleph_parameter": False,
        "aleph_authority": "none",
    }
    (output_dir / "mechanics_protocol_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    report = train(args.tensor, args.catalog, args.output_dir)
    print(json.dumps({
        "wire_run_disjoint_test": report["wire_head"]["experimental_date_disjoint"]["metrics"]["test"],
        "protocol_router_status": report["protocol_router"]["status"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
