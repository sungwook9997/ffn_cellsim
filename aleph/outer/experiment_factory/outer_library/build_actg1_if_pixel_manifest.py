#!/usr/bin/env python3
"""Decode a public author-labelled IF composite into image-level observables.

The source object is a rendered RGB TIFF despite its ``.czi`` filename.  Red is
author-labelled ACTG1 and green is author-labelled occludin.  The image contains
a mixed WT/ACTG1-KO culture, so this adapter deliberately does *not* assign a
genotype to pixels or pretend that tiles are independent cells/replicates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from schema import Observation


SOURCE_SHA256 = "51d283dcbbcc426e7a544a6241fb6a532570e00b547b0736af4368deff2351ee"
ACCESSION = "10.6084/m9.figshare.28309337"
DATASET_ID = "figshare-28309337-actg1-occludin-if"
LAB_GROUP = "mauperin-citi-2025"


def _otsu_threshold(channel: np.ndarray) -> float:
    values = np.asarray(channel, dtype=np.uint8).reshape(-1)
    hist = np.bincount(values, minlength=256).astype(np.float64)
    probability = hist / hist.sum()
    omega = np.cumsum(probability)
    means = np.cumsum(probability * np.arange(256, dtype=np.float64))
    total_mean = means[-1]
    denominator = omega * (1.0 - omega)
    between = np.zeros(256, dtype=np.float64)
    valid = denominator > 0
    between[valid] = (total_mean * omega[valid] - means[valid]) ** 2 / denominator[valid]
    return float(np.argmax(between))


def pixel_metrics(rgb: np.ndarray) -> dict[str, tuple[float, str]]:
    image = np.asarray(rgb)
    if image.ndim != 3 or image.shape[2] != 3 or image.dtype != np.uint8:
        raise ValueError(f"expected uint8 RGB image, found {image.shape} {image.dtype}")
    result: dict[str, tuple[float, str]] = {}
    channels = {"ACTG1_red": image[..., 0], "occludin_green": image[..., 1]}
    thresholds: dict[str, float] = {}
    for name, channel in channels.items():
        threshold = _otsu_threshold(channel)
        thresholds[name] = threshold
        result[f"{name}_intensity_mean"] = (float(channel.mean()), "uint8_intensity")
        result[f"{name}_intensity_std"] = (float(channel.std()), "uint8_intensity")
        result[f"{name}_intensity_p90"] = (float(np.quantile(channel, 0.90)), "uint8_intensity")
        result[f"{name}_intensity_p99"] = (float(np.quantile(channel, 0.99)), "uint8_intensity")
        result[f"{name}_otsu_threshold"] = (threshold, "uint8_intensity")
        result[f"{name}_foreground_fraction"] = (
            float(np.mean(channel > threshold)),
            "fraction",
        )

    red = channels["ACTG1_red"].astype(np.float64)
    green = channels["occludin_green"].astype(np.float64)
    red_fg = red > thresholds["ACTG1_red"]
    green_fg = green > thresholds["occludin_green"]
    union = red_fg | green_fg
    if union.any() and red[union].std() > 0 and green[union].std() > 0:
        pearson = float(np.corrcoef(red[union], green[union])[0, 1])
    else:
        pearson = 0.0
    result["ACTG1_occludin_foreground_pearson"] = (pearson, "dimensionless")
    result["ACTG1_occludin_double_positive_fraction"] = (
        float(np.mean(red_fg & green_fg)),
        "fraction",
    )
    result["ACTG1_positive_with_occludin_fraction"] = (
        float(np.sum(red_fg & green_fg) / max(np.sum(red_fg), 1)),
        "fraction",
    )
    result["occludin_positive_with_ACTG1_fraction"] = (
        float(np.sum(red_fg & green_fg) / max(np.sum(green_fg), 1)),
        "fraction",
    )
    result["rendered_blue_channel_nonzero_fraction"] = (
        float(np.mean(image[..., 2] > 0)),
        "fraction",
    )
    return result


def _load_rgb(path: Path) -> np.ndarray:
    try:
        from PIL import Image
    except ImportError as exc:  # pragma: no cover - environment-specific refusal
        raise RuntimeError(
            "Pillow is required only for the IF pixel adapter; use the bundled "
            "workspace dependency Python recorded in README.md"
        ) from exc
    with Image.open(path) as image:
        if image.format != "TIFF" or image.n_frames != 1:
            raise ValueError(f"expected one-frame rendered TIFF, found {image.format}")
        return np.asarray(image.convert("RGB"), dtype=np.uint8)


def build(source_image: Path) -> list[Observation]:
    payload = source_image.read_bytes()
    if hashlib.sha256(payload).hexdigest() != SOURCE_SHA256:
        raise ValueError("IF source content hash mismatch")
    rgb = _load_rgb(source_image)
    metrics = pixel_metrics(rgb)
    rows: list[Observation] = []
    for observable, (value, unit) in metrics.items():
        locator = f"tiff:frame-0:all-pixels::{observable}"
        identifier = hashlib.sha256(f"{SOURCE_SHA256}:{locator}".encode()).hexdigest()[:20]
        row = Observation(
            observation_id=f"actg1-if:{identifier}",
            dataset_id=DATASET_ID,
            lab_group=LAB_GROUP,
            provider="Figshare",
            accession=ACCESSION,
            license_id="cc-by-4.0",
            source_sha256=SOURCE_SHA256,
            source_locator=locator,
            sample_id=f"{DATASET_ID}:figure-4a-bottom:representative-image",
            biological_replicate="representative_image_single",
            technical_replicate="rendered_RGB_composite",
            condition="mixed_WT_and_ACTG1_KO_clone_2G3",
            cell_type="MDCK epithelial cells",
            cell_state="confluent_mixed_culture",
            modality="immunofluorescence_pixel_derived",
            observable=observable,
            value=value,
            unit=unit,
        )
        row.validate()
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-image", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build(args.source_image)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")
    print(json.dumps({"observations": len(rows), "images": 1}, sort_keys=True))


if __name__ == "__main__":
    main()
