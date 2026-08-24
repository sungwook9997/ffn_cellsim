#!/usr/bin/env python3
"""Evaluate the frozen HPA raw-IF CNN on the frozen IDR0072 FLEX manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np
from PIL import Image
import tifffile
import torch

from train_hpa_raw_if_cnn import RawIFCNN, metrics, multilabel_ece, predict, sigmoid

tifffile.logger().setLevel(logging.CRITICAL)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def normalize(value: str) -> str:
    return re.sub(r"[^0-9a-z]+", "", value.casefold())


def filter_tuple(value: str) -> tuple[int, int] | None:
    match = re.fullmatch(r"\s*(\d+)\s*/\s*(\d+)\s*", value)
    return (int(match.group(1)), int(match.group(2))) if match else None


def channel_planes(
    path: Path, image_name: str, contract: dict, screen_id: int,
    render_contract: dict | None = None,
) -> tuple[np.ndarray, np.ndarray, dict]:
    field_match = re.search(r"Field\s*#\s*(\d+)", image_name, re.IGNORECASE)
    if not field_match:
        raise ValueError(f"cannot resolve field from {image_name}")
    field = int(field_match.group(1))
    with tifffile.TiffFile(path) as tif:
        if 65200 not in tif.pages[0].tags:
            raise ValueError(f"missing FLEX XML tag: {path}")
        root = ET.fromstring(tif.pages[0].tags[65200].value)
        wavelength = {}
        for source in root.findall(".//LightSource"):
            value = source.findtext("Wavelength")
            if value is not None:
                wavelength[source.attrib["ID"]] = int(float(value))
        combinations = {}
        for combination in root.findall(".//LightSourceCombination"):
            combinations[combination.attrib["ID"]] = {
                wavelength[node.attrib["ID"]]
                for node in combination.findall("LightSourceRef")
                if node.attrib.get("ID") in wavelength
            }
        filters = {}
        for combination in root.findall(".//FilterCombination"):
            filters[combination.attrib["ID"]] = {
                normalize(node.attrib.get("ID", "")): node.attrib.get("Filter", "")
                for node in combination.findall("SliderRef")
            }
        candidates = {"nucleus": [], "target_protein": []}
        selected_wavelengths: set[int] = set()
        selected_metadata = []
        optical = contract["field_and_channel_contract"]
        selected_images = [
            image for image in root.findall(".//Image")
            if int(image.findtext("Sublayout", "-1")) == field
        ]
        if render_contract is not None:
            selected_images.sort(key=lambda image: int(image.attrib["BufferNo"]))
            if len(selected_images) != 2:
                raise ValueError(
                    f"screen-render contract requires two images for {path.name} field {field}, found {len(selected_images)}"
                )
            roles = render_contract[str(screen_id)]["ordered_roles"]
            for image, role in zip(selected_images, roles, strict=True):
                buffer_no = int(image.attrib["BufferNo"])
                candidates[role].append(buffer_no)
                selected_metadata.append({
                    "buffer": buffer_no,
                    "camera": normalize(image.findtext("CameraRef", "")),
                    "channel_name": image.findtext("ChannelName", ""),
                    "role": role,
                    "authority": "official_screen_render_channel_order",
                })
        else:
          for image in selected_images:
            buffer_no = int(image.attrib["BufferNo"])
            if buffer_no < 0:
                raise ValueError(f"negative buffer {buffer_no} in {path.name}")
            camera = normalize(image.findtext("CameraRef", ""))
            light_id = image.findtext("LightSourceCombinationRef", "")
            filter_id = image.findtext("FilterCombinationRef", "")
            light_values = combinations.get(light_id, set())
            selected_wavelengths.update(light_values)
            camera_number = "camera1" if camera.endswith(("cam1", "camera1")) else "camera3" if camera.endswith(("cam3", "camera3")) else ""
            emission = filter_tuple(filters.get(filter_id, {}).get(camera_number, ""))
            role = None
            target = optical["target_protein"]
            nucleus = optical["nucleus"]
            if camera_number == "camera1" and emission == (target["emission_center_nm"], target["emission_bandwidth_nm"]):
                role = "target_protein"
            if camera_number == "camera3" and emission is not None and emission[0] == nucleus["emission_center_nm"] and emission[1] in nucleus["allowed_emission_bandwidth_nm"]:
                if role is not None:
                    raise ValueError(f"ambiguous optical role in {path.name} buffer {buffer_no}")
                role = "nucleus"
            if role is not None:
                candidates[role].append(buffer_no)
            selected_metadata.append({"buffer": buffer_no, "camera": camera, "emission_filter": emission, "light_sources_nm": sorted(light_values), "role": role})
        required_lights = render_contract is not None or (
            optical["target_protein"]["required_light_source_nm"] in selected_wavelengths
            and any(value in selected_wavelengths for value in optical["nucleus"]["allowed_light_source_nm"])
        )
        if not required_lights or not candidates["nucleus"] or not candidates["target_protein"]:
            raise ValueError(f"field/optics contract failed for {path.name}: lights={selected_wavelengths}, candidates={candidates}")
        arrays = {}
        for role, pages in candidates.items():
            stack = np.stack([np.asarray(tif.pages[index].asarray(), dtype=np.float32) for index in pages])
            arrays[role] = np.max(stack, axis=0)
    return arrays["nucleus"], arrays["target_protein"], {"field": field, "pages": candidates, "metadata": selected_metadata}


def scale_resize(value: np.ndarray) -> np.ndarray:
    if value.ndim != 2 or not np.isfinite(value).all():
        raise ValueError("FLEX plane is not a finite 2D array")
    low, high = np.percentile(value, [1.0, 99.0])
    if not high > low:
        scaled = np.zeros_like(value, dtype=np.float32)
    else:
        scaled = np.clip((value - low) / (high - low), 0.0, 1.0).astype(np.float32)
    side = min(scaled.shape)
    top = (scaled.shape[0] - side) // 2
    left = (scaled.shape[1] - side) // 2
    crop = np.ascontiguousarray(scaled[top:top + side, left:left + side], dtype=np.float32)
    result = np.asarray(Image.fromarray(crop).resize((128, 128), resample=Image.Resampling.BOX), dtype=np.float32)
    if not np.isfinite(result).all() or result.min() < 0 or result.max() > 1:
        raise ValueError("resized FLEX channel violates finite [0,1]")
    return result


def grouped_metrics(rows: list[dict], logits: np.ndarray, y: np.ndarray, temperature: float, cutoffs: np.ndarray) -> dict:
    groups: dict[str, list[int]] = {}
    for index, row in enumerate(rows):
        groups.setdefault(row["split_group"], []).append(index)
    group_logits, group_y = [], []
    for indices in groups.values():
        labels = y[indices]
        if not np.all(labels == labels[0]):
            raise ValueError("landmark/plate group crosses canonical labels")
        group_logits.append(logits[indices].mean(axis=0))
        group_y.append(labels[0])
    probability = sigmoid(np.stack(group_logits) / temperature)
    result = metrics(np.stack(group_y), probability, cutoffs)
    result["group_count"] = len(groups)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--acquisition-report", type=Path, required=True)
    parser.add_argument("--image-dir", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    protocol_wrapper = json.loads(args.protocol.read_text(encoding="utf-8"))
    render_contract = None
    minimums = None
    if "base_corrected_source_protocol" in protocol_wrapper:
        if not protocol_wrapper["frozen_after_full_xml_metadata_audit_but_before_any_pixel_array_access_or_prediction"]:
            raise ValueError("screen-render protocol was not frozen")
        corrected_path = args.protocol.parent / protocol_wrapper["base_corrected_source_protocol"]["path"]
        corrected = json.loads(corrected_path.read_text(encoding="utf-8"))
        base_path = args.protocol.parent / corrected["base_optical_protocol"]["path"]
        if sha256(base_path) != corrected["base_optical_protocol"]["sha256"]:
            raise ValueError("base optical protocol hash mismatch")
        protocol = json.loads(base_path.read_text(encoding="utf-8"))
        expected_manifest = corrected["corrected_source_contract"]["flex_manifest_sha256"]
        render_contract = protocol_wrapper["screen_render_channel_contract"]
        minimums = protocol_wrapper["post_refusal_minimums"]
    elif "base_optical_protocol" in protocol_wrapper:
        if not protocol_wrapper["frozen_after_source_path_collision_detection_but_before_corrected_acquisition_or_prediction"]:
            raise ValueError("corrected source protocol was not frozen")
        base_path = args.protocol.parent / protocol_wrapper["base_optical_protocol"]["path"]
        if sha256(base_path) != protocol_wrapper["base_optical_protocol"]["sha256"]:
            raise ValueError("base optical protocol hash mismatch")
        protocol = json.loads(base_path.read_text(encoding="utf-8"))
        expected_manifest = protocol_wrapper["corrected_source_contract"]["flex_manifest_sha256"]
    else:
        protocol = protocol_wrapper
        if not protocol["frozen_after_human_and_mouse_metadata_probes_but_before_pixel_array_access_or_prediction"]:
            raise ValueError("final optical protocol was not frozen")
        expected_manifest = protocol["source_contract"]["flex_manifest_sha256"]
    if sha256(args.manifest) != expected_manifest:
        raise ValueError("FLEX manifest differs from frozen protocol")
    if sha256(args.checkpoint) != protocol["frozen_model"]["checkpoint_sha256"]:
        raise ValueError("HPA checkpoint differs from frozen protocol")
    rows = [json.loads(line) for line in args.manifest.read_text(encoding="utf-8").splitlines() if line]
    acquisition = json.loads(args.acquisition_report.read_text(encoding="utf-8"))
    receipts = {int(item["image_id"]): item for item in acquisition["files"]}
    if set(receipts) != {int(row["image_id"]) for row in rows}:
        raise ValueError("acquisition report does not exactly cover the FLEX manifest")
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    classes = checkpoint["classes"]
    class_index = {name: index for index, name in enumerate(classes)}
    images, labels, optical_audit, used_rows, refusals = [], [], [], [], []
    for position, row in enumerate(rows, start=1):
        receipt = receipts[int(row["image_id"])]
        path = args.image_dir / receipt["path"]
        expected_size = receipt.get("retained_bytes", receipt.get("bytes"))
        if path.stat().st_size != expected_size or sha256(path) != receipt["sha256"]:
            raise ValueError(f"FLEX acquisition mismatch: {path}")
        try:
            nucleus, target, audit = channel_planes(
                path, row["image_name"], protocol, int(row["screen_id"]), render_contract
            )
        except ValueError as error:
            refusals.append({"image_id": row["image_id"], "path": receipt["path"], "reason": str(error)})
            continue
        images.append(np.stack((scale_resize(nucleus), scale_resize(target))))
        label = np.zeros(len(classes), dtype=np.float32)
        label[class_index[row["canonical_class"]]] = 1.0
        labels.append(label)
        used_rows.append(row)
        optical_audit.append({"image_id": row["image_id"], **audit})
        if position % 100 == 0 or position == len(rows):
            print(json.dumps({"preprocessed": position, "images": len(rows)}), flush=True)
    if minimums is not None:
        retained_classes = {row["canonical_class"] for row in used_rows}
        retained_proteins = {row["gene_symbol"] for row in used_rows}
        if len(used_rows) < minimums["images"] or len(retained_classes) < minimums["canonical_classes"] or len(retained_proteins) < minimums["landmark_proteins"]:
            raise ValueError("post-refusal support fell below the frozen minimums")
    x = np.stack(images).astype(np.float32)
    y = np.stack(labels).astype(np.float32)
    model = RawIFCNN(len(classes))
    model.load_state_dict(checkpoint["state_dict"])
    torch.set_num_threads(min(8, max(1, torch.get_num_threads())))
    logits = predict(model, torch.from_numpy(x), batch_size=32)
    temperature = float(checkpoint["temperature"])
    probability = sigmoid(logits / temperature)
    cutoffs = np.asarray(checkpoint["decision_thresholds"], dtype=np.float64)
    image_metrics = metrics(y, probability, cutoffs)
    image_metrics["multilabel_ece_15_bin"] = multilabel_ece(y, probability)
    conformal = np.asarray(checkpoint["conformal_thresholds"], dtype=np.float64)
    included = (1.0 - probability) <= conformal
    true_positive = y == 1
    image_metrics["positive_label_conformal_coverage"] = float(included[true_positive].mean())
    image_metrics["mean_prediction_set_size"] = float(included.sum(axis=1).mean())
    group_result = grouped_metrics(used_rows, logits, y, temperature, cutoffs)
    species = {}
    for screen_id in (2952, 2953):
        mask = np.asarray([row["screen_id"] == screen_id for row in used_rows])
        result = metrics(y[mask], probability[mask], cutoffs)
        result["multilabel_ece_15_bin"] = multilabel_ece(y[mask], probability[mask])
        species[str(screen_id)] = result
    endpoints = protocol["frozen_endpoints"]
    endpoint_results = {
        "macro_average_precision": image_metrics["macro_average_precision"] >= endpoints["supported_class_image_level_macro_average_precision_minimum"],
        "macro_F1": image_metrics["macro_f1"] >= endpoints["supported_class_image_level_macro_F1_minimum"],
        "ECE": image_metrics["multilabel_ece_15_bin"] <= endpoints["multilabel_ECE_maximum"],
        "positive_conformal_coverage": image_metrics["positive_label_conformal_coverage"] >= endpoints["positive_label_conformal_coverage_minimum"],
        "mean_prediction_set_size": image_metrics["mean_prediction_set_size"] <= endpoints["mean_prediction_set_size_maximum"],
    }
    args.tensor.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.tensor, logits=logits.astype(np.float32), labels=y,
        image_ids=np.asarray([row["image_id"] for row in used_rows], dtype=np.int64),
        screen_ids=np.asarray([row["screen_id"] for row in used_rows], dtype=np.int32),
        classes=np.asarray(classes),
    )
    report = {
        "schema": "aleph.outer_library.idr0072_if_external_evaluation.v1",
        "status": "domain_shift_diagnostic_passed" if all(endpoint_results.values()) else "domain_shift_diagnostic_failed",
        "protocol_sha256": sha256(args.protocol),
        "manifest_sha256": sha256(args.manifest),
        "checkpoint_sha256": sha256(args.checkpoint),
        "tensor_sha256": sha256(args.tensor),
        "classes": classes,
        "supported_classes": protocol["evaluation_contract"]["supported_primary_classes"],
        "input_image_count": len(rows),
        "evaluated_image_count": len(used_rows),
        "pre_prediction_refusal_count": len(refusals),
        "pre_prediction_refusals": refusals,
        "image_metrics": image_metrics,
        "landmark_plate_group_metrics": group_result,
        "species_screen_metrics": species,
        "endpoint_results": endpoint_results,
        "all_applicable_endpoints_passed": all(endpoint_results.values()),
        "optical_audit": optical_audit,
        "pristine_holdout_authority": False,
        "uncertainty_authority": False,
        "aleph_parameter_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "image_macro_AP": image_metrics["macro_average_precision"], "image_macro_F1": image_metrics["macro_f1"], "ECE": image_metrics["multilabel_ece_15_bin"], "coverage": image_metrics["positive_label_conformal_coverage"], "mean_set_size": image_metrics["mean_prediction_set_size"]}, sort_keys=True))


if __name__ == "__main__":
    main()
