#!/usr/bin/env python3
"""Inference interface for the external mechanics protocol network."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .train import _softmax, _tokenize


def route_protocol(
    checkpoint_path: Path, catalog_path: Path, request: dict[str, Any]
) -> dict[str, Any]:
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    surface = "suspended" if str(request["cell_state"]).startswith("suspended") else "adherent"
    exact_records = [record for record in catalog["records"] if (
        record["probe_geometry"] == request["probe_geometry"]
        and record["observable"] == request["observable"]
        and record["unit"] == request["unit"]
        and ("suspended" if record["cell_state"].startswith("suspended") else "adherent") == surface
    )]
    if not exact_records:
        return {
            "schema": "aleph.outer_library.mechanics_protocol_inference.v1",
            "refused": True,
            "reason": "no exact probe-geometry/observable/unit/surface-state route in the external catalog",
            "nearest_scope_transfer_attempted": False,
            "aleph_authority": "none",
        }
    with np.load(checkpoint_path, allow_pickle=False) as data:
        vocabulary = np.asarray(data["vocabulary"]).astype(str)
        params = {name: np.asarray(data[name]) for name in data.files}
    token_index = {token: index for index, token in enumerate(vocabulary)}
    x = np.zeros((1, len(vocabulary) + 2), dtype=np.float64)
    for token in _tokenize(request):
        if token in token_index:
            x[0, token_index[token]] = 1.0
    x[0, -2] = np.log10(float(request["probe_size_um"]) + 1e-6)
    x[0, -1] = -6.0 if request["time_scale_s"] is None else np.log10(float(request["time_scale_s"]) + 1e-6)
    hidden = np.tanh(x @ params["W_shared"] + params["b_shared"])
    heads = {}
    for head in ("method", "scope", "surface_state"):
        names = np.asarray(params[f"class_name_{head}"]).astype(str)
        probability = _softmax(hidden @ params[f"W_{head}"] + params[f"b_{head}"])[0]
        heads[head] = {
            "prediction": str(names[int(probability.argmax())]),
            "probabilities": {str(name): float(value) for name, value in zip(names, probability, strict=True)},
        }
    return {
        "schema": "aleph.outer_library.mechanics_protocol_inference.v1",
        "refused": False,
        "matched_record_ids": [record["record_id"] for record in exact_records],
        "heads": heads,
        "authority": "ontology_routing_diagnostic_only",
        "external_holdout": False,
        "may_compare_values_across_scopes": False,
        "aleph_authority": "none",
    }


def infer_wire(checkpoint_path: Path, values: list[float]) -> dict[str, Any]:
    with np.load(checkpoint_path, allow_pickle=False) as data:
        params = {name: np.asarray(data[name]) for name in data.files}
    x = (np.log(np.asarray(values, dtype=np.float64))[None, :] - params["input_mean"]) / params["input_scale"]
    hidden = np.tanh(x @ params["W_shared"] + params["b_shared"])
    probability = _softmax(hidden @ params["W_class"] + params["b_class"])[0]
    log_eta = float((hidden @ params["W_eta"] + params["b_eta"])[0, 0] * params["eta_log_scale"] + params["eta_log_mean"])
    names = np.asarray(params["class_name"]).astype(str)
    return {
        "schema": "aleph.outer_library.mechanics_wire_inference.v1",
        "effective_viscosity_pa_s": float(np.exp(log_eta)),
        "effective_viscosity_authority": "Dessard_MRS_operator_interpolation_only",
        "cell_line_probabilities_diagnostic": {str(name): float(value) for name, value in zip(names, probability, strict=True)},
        "cell_line_prediction_refused": True,
        "cell_line_refusal_reason": "experimental-date-disjoint balanced accuracy was 0.385",
        "independent_lab_holdout": False,
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    route = subparsers.add_parser("route")
    route.add_argument("--checkpoint", type=Path, required=True)
    route.add_argument("--catalog", type=Path, required=True)
    route.add_argument("--request", type=Path, required=True)
    wire = subparsers.add_parser("wire")
    wire.add_argument("--checkpoint", type=Path, required=True)
    wire.add_argument("--length-um", type=float, required=True)
    wire.add_argument("--diameter-um", type=float, required=True)
    wire.add_argument("--omega-c-rad-s", type=float, required=True)
    wire.add_argument("--effective-aspect", type=float, required=True)
    args = parser.parse_args()
    if args.command == "route":
        request = json.loads(args.request.read_text(encoding="utf-8"))
        result = route_protocol(args.checkpoint, args.catalog, request)
    else:
        result = infer_wire(args.checkpoint, [args.length_um, args.diameter_um, args.omega_c_rad_s, args.effective_aspect])
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
