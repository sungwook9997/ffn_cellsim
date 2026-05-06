"""V2 imaging measurement-protocol gate — canonical area extractor reproducing PI's CSV ``Area_px``.

Paired lock: ``docs/v2/v2_imaging_measurement_protocol_gate_locked.md``.
Paired sanity gate: ``docs/v2/v2_imaging_measurement_protocol_gate_sanity_gate.md``.

This module is the **source-pipeline reproduction** of PI's analysis at
``/Users/sw1/Desktop/spread_analysis/spread_infer.py``. It is the single
authoritative implementation of the algorithm that produced the
``Area_px`` column in ``data/experimental/260313_*.csv``. The 260313
dataset is a multi-cell / spheroid-stage validation reservoir
(V2-3 / V2-5 territory), NOT a Layer 1 single-cell anchor — this gate
validates the **measurement protocol** (top-down spheroid projection
area), not Layer 1 single-cell physics.

The earlier framing as "V2 Layer-2 Unit 1 / imaging→SingleCellState
loader" was retracted in commit ``f4e09aa`` after PI's phase-scope
correction (260313 = spheroid-level, not single-cell).

PI source = ``spread_infer.py:661`` (``build_row_from_saved_mask``):

1. Read saved mask PNG (``cv2.IMREAD_GRAYSCALE``); operator-final state
   (auto-detected contour OR manual polygon OR re-reviewed contour) is
   already rasterised into the saved mask, so the mask is the
   authoritative artefact.
2. ``cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)``.
3. ``select_central_contour(contours, w, h)`` (``spread_infer.py:312``)
   selects the most-centred non-border contour with cost
   ``dist + 0.0005 / (area + 1e-6)`` under the constants
   ``CENTER_MAX_DIST_FRAC = 0.45``, ``BORDER_MARGIN_PX = 12``,
   ``REJECT_TOUCH_BORDER = True``.
4. ``Area_px = cv2.contourArea(cnt)`` — Green's theorem polygon area on
   the integer-vertex contour, **NOT** a connected-component pixel count.
5. ``Area_um2 = Area_px * PIXEL_AREA_UM2`` with
   ``PIXEL_AREA_UM2 = 4.194304 = 2.048**2`` (Hard Rule 11 inline
   unit-chain: pixel_count × pixel_size_um² = μm²; pixel side 2.048 μm).

This is reproduction, NOT calibration. Hard Rule 1 (``forbidden_use:
physics_parameter_fitting``) still applies; the algorithm is read-only
toward PI's data.

Tolerance contract
==================

Under the source-pipeline contract on identical OpenCV major version,
the algorithm is byte-exact. The measurement-protocol gate lock pins the contract
threshold at **1.5% relative residual** as a cross-version drift
budget; observed residual on cv2 4.13.0 (validated 2026-05-06) is
**0.0** on every non-excluded row. The 1.5% slack is reserved for
``cv2`` cross-version contour discretisation, **NOT** for fitting.

PI-directed provenance exclusion: Lam4 Position(31)
====================================================

By PI directive (resolution of ``id=2278``, design-discussion
``id=2304`` Codex guardrail), Lam4 Position(31) is excluded from the
protocol-identification set; reported separately. The original
24.79% mismatch flagged in design-discussion ``id=2278`` was
specifically a ``largest_cc area`` vs ``Area_px`` measurement-protocol
mismatch — under canonical ``select_central_contour`` reproduction,
Pos31 also reproduces at 0.0 residual, but the provenance exclusion
is preserved per directive (``provenance, NOT tolerance tuning``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import cv2
import numpy as np


# Source constants — must match spread_infer.py:55-57. Changing any of
# these breaks byte-exact reproduction of PI's CSV Area_px.
BORDER_MARGIN_PX = 12
REJECT_TOUCH_BORDER = True
CENTER_MAX_DIST_FRAC = 0.45
SELECT_CENTRAL_CONTOUR_AREA_COST_WEIGHT = 0.0005
PIXEL_AREA_UM2 = 4.194304  # = 2.048**2 μm² per pixel; Hard Rule 11 inline unit-chain.

# PI-directed provenance exclusion set (id=2278 resolution, Codex
# guardrail id=2304). Keyed by (group, series). Excluded rows are
# surfaced with provenance_status=excluded_by_pi_measurement_provenance
# and MUST NOT contribute to thresholds, calibration, or rates.
PI_DIRECTED_PROVENANCE_EXCLUSIONS: frozenset[tuple[str, str]] = frozenset({
    ("Lam4", "Position(31)"),
})

PROVENANCE_STATUS_NORMAL = "normal"
PROVENANCE_STATUS_EXCLUDED = "excluded_by_pi_measurement_provenance"
PROVENANCE_STATUS_MASK_MISSING = "mask_missing"
PROVENANCE_STATUS_MASK_READ_FAIL = "mask_read_fail"
PROVENANCE_STATUS_NO_CONTOUR = "no_central_contour"


@dataclass(frozen=True, slots=True)
class AreaExtractionResult:
    """One row of source-pipeline reproduction.

    ``provenance_status`` is the discipline channel for downstream
    validators. Pos31 is surfaced (not silently filtered) and downstream
    code MUST refuse to use ``excluded_by_pi_measurement_provenance``
    rows for tolerance derivation, calibration, or acceptance windows.
    """

    group: str
    series: str
    frame: int
    filename: str
    csv_area_px: float
    repro_area_px: float | None
    rel_residual: float | None
    provenance_status: str

    @property
    def is_excluded_by_pi(self) -> bool:
        return self.provenance_status == PROVENANCE_STATUS_EXCLUDED


def contour_centroid(cnt) -> tuple[float, float] | None:
    M = cv2.moments(cnt)
    if M["m00"] == 0:
        return None
    return (M["m10"] / M["m00"], M["m01"] / M["m00"])


def touches_border(cnt, w: int, h: int, margin_px: int = BORDER_MARGIN_PX) -> bool:
    x, y, bw, bh = cv2.boundingRect(cnt)
    if x <= margin_px or y <= margin_px:
        return True
    if (x + bw) >= (w - margin_px) or (y + bh) >= (h - margin_px):
        return True
    return False


def select_central_contour(contours, img_w: int, img_h: int):
    """Mirror of ``spread_infer.py:312``.

    Picks the most-centred non-border contour by minimising
    ``cost = dist + 0.0005 / (area + 1e-6)``. If no contour passes the
    border / max-distance filter, falls back to the closest-to-center
    contour (without border check), matching PI source line 345-357.
    """
    if not contours:
        return None
    cx0 = img_w / 2.0
    cy0 = img_h / 2.0
    diag = math.sqrt(img_w * img_w + img_h * img_h)
    max_dist = diag * CENTER_MAX_DIST_FRAC

    best = None
    best_cost = 1e18

    for cnt in contours:
        area = cv2.contourArea(cnt)
        if area <= 0:
            continue
        c = contour_centroid(cnt)
        if c is None:
            continue
        cx, cy = c
        dist = math.sqrt((cx - cx0) ** 2 + (cy - cy0) ** 2)
        if dist > max_dist:
            continue
        if REJECT_TOUCH_BORDER and touches_border(cnt, img_w, img_h):
            continue
        cost = dist + SELECT_CENTRAL_CONTOUR_AREA_COST_WEIGHT * (1.0 / (area + 1e-6))
        if cost < best_cost:
            best_cost = cost
            best = cnt

    if best is None:
        for cnt in contours:
            area = cv2.contourArea(cnt)
            if area <= 0:
                continue
            c = contour_centroid(cnt)
            if c is None:
                continue
            cx, cy = c
            dist = math.sqrt((cx - cx0) ** 2 + (cy - cy0) ** 2)
            if dist < best_cost:
                best_cost = dist
                best = cnt
    return best


def extract_area_px_from_mask(mask_path: Path) -> tuple[float | None, str]:
    """Run the canonical algorithm on one saved mask PNG.

    Returns ``(area_px, status)`` where ``area_px`` is None on failure
    paths (mask missing, mask read fail, no central contour). The
    ``status`` is one of ``PROVENANCE_STATUS_*`` (excluding
    ``EXCLUDED`` — that one is keyed by ``(group, series)`` and
    applied at row level by ``classify_provenance``).
    """
    if not mask_path.exists():
        return None, PROVENANCE_STATUS_MASK_MISSING
    mask = cv2.imread(str(mask_path), cv2.IMREAD_GRAYSCALE)
    if mask is None:
        return None, PROVENANCE_STATUS_MASK_READ_FAIL
    cnts, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    cnt = select_central_contour(cnts, mask.shape[1], mask.shape[0])
    if cnt is None:
        return None, PROVENANCE_STATUS_NO_CONTOUR
    return float(cv2.contourArea(cnt)), PROVENANCE_STATUS_NORMAL


def classify_provenance(group: str, series: str) -> str:
    """Apply PI-directed provenance exclusion at the (group, series) level.

    Exclusion is surfaced as a status flag, NOT a silent filter.
    Downstream validators MUST refuse to use rows with status
    ``excluded_by_pi_measurement_provenance`` for tolerance
    derivation, calibration, rates, or acceptance windows.
    """
    if (group, series) in PI_DIRECTED_PROVENANCE_EXCLUSIONS:
        return PROVENANCE_STATUS_EXCLUDED
    return PROVENANCE_STATUS_NORMAL


def reproduce_one_row(
    group: str,
    series: str,
    frame: int,
    filename: str,
    csv_area_px: float,
    mask_path: Path,
) -> AreaExtractionResult:
    """Compute one source-pipeline reproduction with provenance discipline."""
    repro_area, mask_status = extract_area_px_from_mask(mask_path)

    pi_status = classify_provenance(group, series)
    # If a row would otherwise be EXCLUDED but the mask itself is missing
    # or unreadable, the missing/fail status is reported (more specific)
    # — but PI exclusion is still recorded via a side flag accessible
    # through the (group, series) lookup. For simplicity here, the
    # missing/fail path takes priority; downstream still has the PI
    # exclusion set as a frozen module-level constant for cross-check.
    if mask_status != PROVENANCE_STATUS_NORMAL:
        return AreaExtractionResult(
            group=group,
            series=series,
            frame=frame,
            filename=filename,
            csv_area_px=csv_area_px,
            repro_area_px=None,
            rel_residual=None,
            provenance_status=mask_status,
        )

    if csv_area_px == 0:
        rel = 0.0 if repro_area == 0 else float("nan")
    else:
        rel = abs(repro_area - csv_area_px) / abs(csv_area_px)

    return AreaExtractionResult(
        group=group,
        series=series,
        frame=frame,
        filename=filename,
        csv_area_px=csv_area_px,
        repro_area_px=repro_area,
        rel_residual=rel,
        provenance_status=pi_status,
    )


def aggregate_normal_residuals(results: Iterable[AreaExtractionResult]) -> dict:
    """Summary statistics over the non-excluded set.

    By contract, the measurement-protocol-gate acceptance threshold (1.5% relative
    residual) is computed here over normal rows ONLY — Pos31 and
    other PI-excluded rows are reported in a separate aggregate so
    they cannot contribute to the acceptance decision.
    """
    normal: list[float] = []
    excluded: list[float] = []
    other_status_counts: dict[str, int] = {}
    for r in results:
        if r.rel_residual is None:
            other_status_counts[r.provenance_status] = (
                other_status_counts.get(r.provenance_status, 0) + 1
            )
            continue
        if r.provenance_status == PROVENANCE_STATUS_NORMAL:
            normal.append(r.rel_residual)
        elif r.provenance_status == PROVENANCE_STATUS_EXCLUDED:
            excluded.append(r.rel_residual)
        else:
            other_status_counts[r.provenance_status] = (
                other_status_counts.get(r.provenance_status, 0) + 1
            )
    arr = np.array(normal, dtype=np.float64) if normal else np.empty(0)
    arr_ex = np.array(excluded, dtype=np.float64) if excluded else np.empty(0)
    return {
        "rows_normal": int(arr.size),
        "rows_excluded_by_pi": int(arr_ex.size),
        "other_status_counts": other_status_counts,
        "normal_max_rel_residual": float(arr.max()) if arr.size else None,
        "normal_mean_rel_residual": float(arr.mean()) if arr.size else None,
        "normal_p99_rel_residual": (
            float(np.quantile(arr, 0.99)) if arr.size else None
        ),
        "excluded_max_rel_residual": float(arr_ex.max()) if arr_ex.size else None,
        "excluded_mean_rel_residual": float(arr_ex.mean()) if arr_ex.size else None,
    }


# Measurement-protocol-gate acceptance contract — fixed at lock seal time, NOT tuned
# to observed residuals. The tolerance is justified as a cv2 cross-
# version drift budget (PI source uses cv2 4.x; this impl validated
# under cv2 4.13.0 at 0.0 observed). Lowering this threshold based on
# observation is permitted; raising it is a contract change requiring
# a new PI directive.
ACCEPTANCE_REL_RESIDUAL_THRESHOLD = 0.015
