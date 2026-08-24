#!/usr/bin/env python3
"""Build raw-image and morphology tensors for BBBC053 FCCP stress IF."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile

import numpy as np


def _features(image: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    low, high = np.quantile(image, (0.01, 0.995))
    scale = max(1.0, float(high - low))
    normalized = np.clip((image.astype(np.float32) - low) / scale, 0.0, 1.0)
    gx = np.diff(normalized, axis=1)
    gy = np.diff(normalized, axis=0)
    q10, q25, q50, q75, q90 = np.quantile(normalized, (0.1, 0.25, 0.5, 0.75, 0.9))
    feature = np.asarray([
        normalized.mean(), normalized.std(), q10, q25, q50, q75, q90,
        np.abs(gx).mean(), np.abs(gy).mean(),
        np.sqrt(gx[:-1, :] ** 2 + gy[:, :-1] ** 2).mean(),
        (normalized > 0.25).mean(), (normalized > 0.50).mean(),
        (normalized > 0.75).mean(),
    ], dtype=np.float32)
    return normalized, feature


def build(archive_path: Path, manifest_path: Path, tensor_path: Path) -> dict[str, object]:
    try:
        from PIL import Image
    except ImportError as error:
        raise RuntimeError("Pillow is required only for BBBC053 TIFF decoding") from error
    rows, compact, features, labels, sample_ids = [], [], [], [], []
    with zipfile.ZipFile(archive_path) as archive:
        if archive.testzip() is not None:
            raise ValueError("BBBC053 ZIP CRC failure")
        for member in sorted(archive.namelist()):
            match = re.fullmatch(
                r"FCCP/(DMSO|FCCP)/MAX_(?:DMSO|FCCP)_ \((\d+)\) - Deconvolved\.tif",
                member,
            )
            if not match:
                raise ValueError(f"unexpected BBBC053 archive member {member!r}")
            condition, image_number = match.group(1), int(match.group(2))
            with Image.open(io.BytesIO(archive.read(member))) as image:
                if image.mode != "I;16B" or image.size != (2048, 2048):
                    raise ValueError(f"BBBC053 TIFF contract changed: {image.mode} {image.size}")
                pixels = np.asarray(image, dtype=np.uint16)
                normalized, feature = _features(pixels)
                reduced = Image.fromarray((normalized * 255).astype(np.uint8)).resize(
                    (64, 64), resample=Image.Resampling.BOX
                )
                compact.append(np.asarray(reduced, dtype=np.uint8))
                features.append(feature)
            label = int(condition == "FCCP")
            sample_id = f"BBBC053:{condition}:image{image_number:02d}"
            labels.append(label)
            sample_ids.append(sample_id)
            rows.append({
                "schema": "aleph.outer_library.observation.v1",
                "observation_id": sample_id,
                "dataset_id": "bbbc053-v1-CAD-FCCP-mitochondria",
                "lab_group": "Vitriol_Lab_Augusta_University",
                "sample_id": sample_id,
                "split_group": sample_id,
                "modality": "TOM20_super_resolution_IF_image",
                "cell_type": "murine_Cath_a_differentiated_CAD",
                "cell_state": "FCCP_mitochondrial_stress" if label else "DMSO_vehicle",
                "condition": condition,
                "unit": "categorical_author_directory_label",
                "value": label,
                "biological_replicate": "not_reported",
                "source_locator": member,
                "license_id": "cc-by-nc-sa-3.0",
                "training_role": "same_provider_mitochondrial_stress_supervision",
                "external_lab_holdout": False,
                "aleph_authority": "none",
                "may_select_aleph_parameter": False,
            })
    if len(rows) != 58 or sum(labels) != 29:
        raise ValueError(f"BBBC053 condition cardinality changed: {len(rows)}, FCCP={sum(labels)}")
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows), encoding="utf-8"
    )
    np.savez_compressed(
        tensor_path,
        compact_image=np.stack(compact),
        X=np.stack(features),
        label=np.asarray(labels, dtype=np.uint8),
        sample_id=np.asarray(sample_ids, dtype="U30"),
    )
    return {
        "schema": "aleph.outer_library.bbbc053_mito_stress.v1",
        "image_count": len(rows),
        "condition_counts": {"DMSO": labels.count(0), "FCCP": labels.count(1)},
        "source_image_shape": [2048, 2048],
        "compact_image_shape": [58, 64, 64],
        "morphology_feature_count": int(np.stack(features).shape[1]),
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "tensor_sha256": hashlib.sha256(tensor_path.read_bytes()).hexdigest(),
        "biological_replicates_reported": False,
        "independent_lab_holdout": False,
        "status": "same_provider_training_only",
        "aleph_authority": "none",
        "may_select_aleph_parameter": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build(args.archive, args.manifest, args.tensor)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
