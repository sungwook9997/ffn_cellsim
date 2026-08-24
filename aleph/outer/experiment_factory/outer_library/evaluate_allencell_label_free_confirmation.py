#!/usr/bin/env python3
"""Evaluate frozen v2/v3 models on Allen raw hiPSC images once."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import zoom
import tifffile
import torch

from build_lightmycells_tensor import normalize_resize
from evaluate_lightmycells_confirmation import coverage_report
from train_lightmycells_model import group_metrics, predict, sha256
from train_lightmycells_model_v2 import (
    ConditionalResidualUNet,
    baseline_comparison,
)
from train_lightmycells_model_v3 import DomainInvariantUNet, leave_study_out_ood


def read_channel(path: Path, channel: int) -> tuple[np.ndarray, str]:
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        if series.axes != "ZCYX":
            raise ValueError(f"unsupported Allen axes {series.axes}: {path}")
        z_count, channel_count = series.shape[:2]
        if not 0 <= channel < channel_count:
            raise ValueError(f"channel {channel} outside {series.shape}: {path}")
        array = np.stack([
            tif.pages[z * channel_count + channel].asarray()
            for z in range(z_count)
        ])
        return array, tif.ome_metadata or ""


def dna_channel(ome: str) -> int:
    tags = re.findall(r"<Channel\s+[^>]*>", ome)
    for index, tag in enumerate(tags):
        if "H3342" in tag or "Hoechst" in tag:
            return index
    raise ValueError("FOV OME metadata has no Hoechst DNA channel")


def physical_size_z(ome: str) -> float:
    match = re.search(r'PhysicalSizeZ="([^"]+)"', ome)
    if match is None:
        raise ValueError("FOV OME metadata has no PhysicalSizeZ")
    return float(match.group(1))


def crop_yx_with_zero_padding(
    volume: np.ndarray, y0: int, y1: int, x0: int, x1: int,
) -> np.ndarray:
    result = np.zeros(
        (volume.shape[0], y1 - y0, x1 - x0), dtype=volume.dtype
    )
    source_y0, source_y1 = max(0, y0), min(volume.shape[1], y1)
    source_x0, source_x1 = max(0, x0), min(volume.shape[2], x1)
    target_y0, target_x0 = source_y0 - y0, source_x0 - x0
    result[
        :, target_y0:target_y0 + source_y1 - source_y0,
        target_x0:target_x0 + source_x1 - source_x0,
    ] = volume[:, source_y0:source_y1, source_x0:source_x1]
    return result


def endpoint_summary(
    prediction: np.ndarray, truth: np.ndarray, plates: np.ndarray,
    mean_map: np.ndarray, denominators: dict[str, float],
) -> dict[str, object]:
    targets = np.asarray(["Mitochondria"] * len(truth))
    metrics = group_metrics(prediction, truth, targets, plates)
    baseline = baseline_comparison(
        prediction, truth, targets, plates, {"Mitochondria": mean_map}
    )
    row = metrics["per_target"]["Mitochondria"]
    values = {
        "Pearson": row["macro_study_Pearson"],
        "global_SSIM": row["macro_study_global_SSIM"],
        "relative_MAE_improvement": baseline["Mitochondria"][
            "relative_MAE_improvement"
        ],
    }
    ratios = {key: values[key] / denominators[key] for key in values}
    return {
        "values": values, "endpoint_ratios": ratios,
        "minimum_endpoint_ratio": min(ratios.values()),
        "all_spatial_endpoints_passed": all(value >= 1 for value in ratios.values()),
        "metrics": metrics, "baseline": baseline,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--training-tensor", type=Path, required=True)
    parser.add_argument("--v2-checkpoint", type=Path, required=True)
    parser.add_argument("--v3-checkpoint", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    acquisition = json.loads(args.acquisition.read_text())
    if acquisition["protocol_sha256"] != sha256(args.protocol):
        raise ValueError("acquisition does not pin evaluation protocol")
    if acquisition["manifest_sha256"] != sha256(args.manifest):
        raise ValueError("acquisition does not pin candidate manifest")
    receipt = {row["path"]: row for row in acquisition["files"]}
    for relative, row in receipt.items():
        path = args.data_root / relative
        if path.stat().st_size != row["size_bytes"] or sha256(path) != row["sha256"]:
            raise ValueError(f"Allen payload changed after acquisition: {relative}")
    rows = [json.loads(line) for line in args.manifest.read_text().splitlines()]
    inputs, truths, geometry = [], [], []
    for number, row in enumerate(rows, 1):
        crop_path = args.data_root / row["crop_raw"]
        fov_path = args.data_root / row["fov_path"]
        crop_dna, _ = read_channel(crop_path, 0)
        target, _ = read_channel(crop_path, 2)
        brightfield_index = int(row["ChannelNumberBrightfield"])
        brightfield, ome = read_channel(fov_path, brightfield_index)
        fov_dna, _ = read_channel(fov_path, dna_channel(ome))
        z0, z1, y0, y1, x0, x1 = ast.literal_eval(row["roi"])
        if crop_dna.shape != (z1 - z0, y1 - y0, x1 - x0):
            raise ValueError(f"ROI/crop shape mismatch for CellId {row['CellId']}")
        factor = physical_size_z(ome) / float(
            ast.literal_eval(row["scale_micron"])[0]
        )
        bf_iso = zoom(
            crop_yx_with_zero_padding(brightfield, y0, y1, x0, x1),
            (factor, 1, 1), order=1,
            prefilter=False,
        )[z0:z1]
        dna_iso = zoom(
            crop_yx_with_zero_padding(fov_dna, y0, y1, x0, x1),
            (factor, 1, 1), order=1,
            prefilter=False,
        )[z0:z1]
        if bf_iso.shape != target.shape or dna_iso.shape != crop_dna.shape:
            raise ValueError(f"isotropic ROI shape mismatch for CellId {row['CellId']}")
        difference = np.abs(dna_iso.astype(np.int64) - crop_dna.astype(np.int64))
        middle = (len(target) - 1) // 2
        input64, _ = normalize_resize(bf_iso[middle], 64)
        truth64, _ = normalize_resize(target[middle], 64)
        inputs.append(input64); truths.append(truth64)
        geometry.append({
            "CellId": row["CellId"], "FOVId": row["FOVId"],
            "PlateId": row["PlateId"],
            "DNA_exact_array_equal": bool(np.array_equal(dna_iso, crop_dna)),
            "DNA_uint16_MAE": float(difference.mean()),
            "DNA_exact_pixel_fraction": float((difference == 0).mean()),
            "input_plane_raw_standard_deviation": float(bf_iso[middle].std()),
            "target_plane_raw_standard_deviation": float(target[middle].std()),
        })
        print(json.dumps({"decoded": number, "total": len(rows)}, sort_keys=True), flush=True)
    inputs = np.asarray(inputs, dtype=np.float16)
    truths = np.asarray(truths, dtype=np.float16)
    plates = np.asarray([row["PlateId"] for row in rows])
    target_index = np.ones(len(rows), dtype=np.int64)
    modality_index = np.zeros(len(rows), dtype=np.int64)
    indices = np.arange(len(rows))

    v2 = torch.load(args.v2_checkpoint, map_location="cpu", weights_only=True)
    model2 = ConditionalResidualUNet(int(v2["base_channels"]))
    model2.load_state_dict(v2["state_dict"])
    prediction2, embedding2 = predict(
        model2, inputs, truths, target_index, modality_index, indices,
        32, 4, 3,
    )
    v3 = torch.load(args.v3_checkpoint, map_location="cpu", weights_only=True)
    spec = v3["model"]
    model3 = DomainInvariantUNet(
        int(spec["base_channels"]), int(spec["condition_channels"]),
        int(spec["study_count"]), float(spec["gradient_reversal_strength"]),
    )
    model3.load_state_dict(v3["state_dict"])
    prediction3, embedding3 = predict(
        model3, inputs, truths, target_index, modality_index, indices,
        32, 4, 3,
    )
    denominators = protocol["evaluation"]["endpoint_denominators"]
    v2_summary = endpoint_summary(
        prediction2, truths.astype(np.float32), plates,
        v2["fit_mean_maps"]["Mitochondria"].numpy(), denominators,
    )
    v3_summary = endpoint_summary(
        prediction3, truths.astype(np.float32), plates,
        v3["fit_mean_maps"]["Mitochondria"].numpy(), denominators,
    )
    coverage = coverage_report(
        prediction3, truths.astype(np.float32),
        np.asarray(["Mitochondria"] * len(rows)), plates,
        v3["calibration_thresholds"], ["Mitochondria"],
    )["required_target_macro_study_pixel_interval_coverage"]
    with np.load(args.training_tensor, allow_pickle=False) as data:
        train_inputs, train_truths = data["inputs"], data["targets"]
        train_roles, train_studies = data["role"], data["study_id"]
        train_targets, train_modalities = data["target_channel"], data["input_modality"]
    target_order, modality_order = v3["target_order"], v3["modality_order"]
    target_lookup = {name: index for index, name in enumerate(target_order)}
    modality_lookup = {name: index for index, name in enumerate(modality_order)}
    train_target_index = np.asarray([target_lookup[x] for x in train_targets])
    train_modality_index = np.asarray([modality_lookup[x] for x in train_modalities])
    fit = np.flatnonzero(train_roles == "fit")
    _, fit_embedding = predict(
        model3, train_inputs, train_truths, train_target_index,
        train_modality_index, fit, 48, 4, 3,
    )
    ood = leave_study_out_ood(
        fit_embedding, train_studies[fit], embedding3, plates,
    )
    exact_geometry = all(row["DNA_exact_array_equal"] for row in geometry)
    evaluation = protocol["evaluation"]
    comparison_improvement = (
        v3_summary["minimum_endpoint_ratio"] - v2_summary["minimum_endpoint_ratio"]
    )
    diagnostic_passes = {
        "v3_spatial_endpoints": v3_summary["all_spatial_endpoints_passed"],
        "v3_pixel_interval_coverage": evaluation[
            "v3_pixel_interval_coverage_range"
        ][0] <= coverage <= evaluation["v3_pixel_interval_coverage_range"][1],
        "v3_embedding_OOD_AUROC": ood["study_macro_AUROC"]
        >= evaluation["v3_embedding_OOD_AUROC_minimum"],
        "v3_minimum_ratio_improvement_over_v2": comparison_improvement
        >= evaluation["v3_minimum_endpoint_ratio_improvement_over_v2"],
    }
    args.tensor.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.tensor, inputs=inputs, targets=truths,
        v2_predictions=prediction2.astype(np.float16),
        v3_predictions=prediction3.astype(np.float16), plates=plates,
    )
    figure_indices = np.linspace(0, len(rows) - 1, 8, dtype=int)
    fig, axes = plt.subplots(8, 5, figsize=(11, 17), constrained_layout=True)
    titles = ("Brightfield", "TOMM20 truth", "v2", "v3", "|v3-truth|")
    for column, title in enumerate(titles):
        axes[0, column].set_title(title)
    for axis_row, index in zip(axes, figure_indices, strict=True):
        images = (
            inputs[index], truths[index], prediction2[index], prediction3[index],
            np.abs(prediction3[index] - truths[index]),
        )
        for axis, image in zip(axis_row, images, strict=True):
            axis.imshow(image, cmap="gray", vmin=0, vmax=1)
            axis.axis("off")
        axis_row[0].set_ylabel(f"Plate {plates[index]}")
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure, dpi=180)
    plt.close(fig)
    report = {
        "schema": "aleph.outer_library.allencell_label_free_confirmation_result.v1",
        "protocol_sha256": sha256(args.protocol),
        "acquisition_sha256": sha256(args.acquisition),
        "candidate_count": len(rows), "plate_count": len(set(plates.tolist())),
        "geometry": {
            "all_DNA_exact_array_equal": exact_geometry,
            "rows": geometry,
            "mean_DNA_uint16_MAE": float(np.mean([row["DNA_uint16_MAE"] for row in geometry])),
        },
        "preregistered_confirmation": {
            "status": "evaluated" if exact_geometry else "refused_geometry_exact_match",
            "all_endpoints_passed": bool(exact_geometry and all(diagnostic_passes.values())),
        },
        "post_refusal_non_authoritative_diagnostic": {
            "v2": v2_summary, "v3": v3_summary,
            "v3_pixel_interval_coverage": coverage,
            "v3_embedding_OOD": ood,
            "v3_minimum_endpoint_ratio_improvement_over_v2": comparison_improvement,
            "passes": diagnostic_passes,
            "all_diagnostic_endpoints_passed": all(diagnostic_passes.values()),
        },
        "tensor_sha256": sha256(args.tensor),
        "figure_sha256": sha256(args.figure),
        "single_use_confirmation_consumed": True,
        "production_eligible": False,
        "observable_evidence_eligible": False,
        "uncertainty_authority": False,
        "aleph_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "confirmation_status": report["preregistered_confirmation"]["status"],
        "v2": v2_summary["values"], "v3": v3_summary["values"],
        "v3_minus_v2_minimum_ratio": comparison_improvement,
        "v3_OOD_AUROC": ood["study_macro_AUROC"],
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
