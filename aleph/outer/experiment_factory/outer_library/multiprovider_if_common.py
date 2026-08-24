"""Shared frozen preprocessing and label mapping for multiprovider raw IF."""

from __future__ import annotations

import re
from typing import Any

import numpy as np
from PIL import Image


def normalized_label(value: str) -> str:
    """Apply the exact normalization frozen in the multiprovider protocol."""
    return re.sub(r"[-\s]+", "_", value.casefold().strip())


def map_author_label(value: str, ontology: list[dict[str, Any]]) -> list[str]:
    """Map one author label with exact-match precedence, then substring fallback."""
    exact = normalized_label(value)
    exact_matches = [
        item["canonical"]
        for item in ontology
        if any(normalized_label(token) == exact for token in item["author_tokens"])
    ]
    if exact_matches:
        return exact_matches
    lowered = value.casefold()
    return [
        item["canonical"]
        for item in ontology
        if any(token.casefold() in lowered for token in item["author_tokens"])
    ]


def map_author_labels(
    values: list[str], ontology: list[dict[str, Any]]
) -> list[str]:
    return sorted({label for value in values for label in map_author_label(value, ontology)})


def normalize_crop_resize(
    array: np.ndarray, side: int
) -> tuple[np.ndarray, bool]:
    """p1-p99 normalize, center-crop square, BOX resize, and uint8 encode."""
    values = np.asarray(array, dtype=np.float32)
    if values.ndim != 2 or not np.isfinite(values).all():
        raise ValueError("raw IF channel must be a finite two-dimensional array")
    low, high = np.percentile(values, [1.0, 99.0])
    constant = bool(not np.isfinite(low) or not np.isfinite(high) or high <= low)
    if constant:
        normalized = np.zeros(values.shape, dtype=np.float32)
    else:
        normalized = np.clip((values - low) / (high - low), 0.0, 1.0)
    crop_side = min(normalized.shape)
    top = (normalized.shape[0] - crop_side) // 2
    left = (normalized.shape[1] - crop_side) // 2
    crop = np.ascontiguousarray(
        normalized[top : top + crop_side, left : left + crop_side]
    )
    resized = np.asarray(
        Image.fromarray(crop).resize((side, side), Image.Resampling.BOX),
        dtype=np.float32,
    )
    if not np.isfinite(resized).all() or resized.min() < 0 or resized.max() > 1:
        raise ValueError("normalized IF channel left finite [0,1] range")
    return np.rint(resized * 255.0).astype(np.uint8), constant
