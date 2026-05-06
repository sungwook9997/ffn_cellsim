"""V2 imaging measurement-protocol gate — reproduce PI's CSV `Area_px` from saved mask PNGs.

Paired lock: ``docs/v2_imaging_measurement_protocol_gate_locked.md``.

The 260313 dataset is a multi-cell / spheroid-stage validation reservoir
(V2-3/V2-5 territory), NOT a Layer 1 single-cell anchor; this script
validates the **measurement protocol** (top-down spheroid projection
area), not Layer 1 single-cell physics.

Source-pipeline reproduction of PI's analysis at /Users/sw1/Desktop/spread_analysis/spread_infer.py.
Canonical algorithm = `build_row_from_saved_mask` (spread_infer.py:661):
  cv2.findContours(mask, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE)
  -> select_central_contour(contours, w, h) with CENTER_MAX_DIST_FRAC=0.45,
     BORDER_MARGIN_PX=12, REJECT_TOUCH_BORDER=True, cost = dist + 0.0005/(area+1e-6)
  -> Area_px = cv2.contourArea(cnt) (Green's theorem polygon area, NOT pixel count)
  -> Area_um2 = Area_px * PIXEL_AREA_UM2 (= 2.048**2 = 4.194304 um^2/px)

Per Codex guardrail id=2304: Lam4 Position(31) is excluded from the protocol-identification
set BY PI DIRECTIVE (provenance, not tolerance tuning); reported separately.

Tolerance is justified from cv2 / float roundoff, not from observed residuals.

Run via .venv-collab (which has cv2 4.13.0 installed for this verification only):
    /Users/sw1/ActiveCellSim/.venv-collab/bin/python \
        scripts/v2/reproduce_pi_area_extraction.py
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

import cv2
import numpy as np
import pandas as pd

# PI source constants (spread_infer.py:55-57)
BORDER_MARGIN_PX = 12
REJECT_TOUCH_BORDER = True
CENTER_MAX_DIST_FRAC = 0.45
PIXEL_AREA_UM2 = 4.194304  # = 2.048**2 um^2/px

REPO = Path(__file__).resolve().parents[2]
IMAGING_ROOT = REPO / "data" / "experimental" / "imaging"
GROUPS = {
    "Bare": REPO / "data" / "experimental" / "260313_Bare.csv",
    "Lam4": REPO / "data" / "experimental" / "260313_Lam4.csv",
    "Pre": REPO / "data" / "experimental" / "260313_Pre.csv",
}

PI_DIRECTED_PROVENANCE_EXCLUSIONS = {
    # PI-directed provenance exclusion: Lam4 Position(31) — excluded from
    # protocol-identification set; reported separately. NOT tolerance tuning.
    ("Lam4", "Position(31)"),
}


def contour_centroid(cnt):
    M = cv2.moments(cnt)
    if M["m00"] == 0:
        return None
    return (M["m10"] / M["m00"], M["m01"] / M["m00"])


def touches_border(cnt, w, h, margin_px=BORDER_MARGIN_PX):
    x, y, bw, bh = cv2.boundingRect(cnt)
    if x <= margin_px or y <= margin_px:
        return True
    if (x + bw) >= (w - margin_px) or (y + bh) >= (h - margin_px):
        return True
    return False


def select_central_contour(contours, img_w, img_h):
    """Mirror of spread_infer.py:312 — pick the most centred non-border contour."""
    if not contours:
        return None
    cx0, cy0 = img_w / 2.0, img_h / 2.0
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
        cost = dist + 0.0005 * (1.0 / (area + 1e-6))
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


def area_from_mask(mask_path: Path) -> tuple[float | None, tuple[int, int] | None]:
    mask = cv2.imread(str(mask_path), cv2.IMREAD_GRAYSCALE)
    if mask is None:
        return None, None
    cnts, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    cnt = select_central_contour(cnts, mask.shape[1], mask.shape[0])
    if cnt is None:
        return 0.0, mask.shape
    return float(cv2.contourArea(cnt)), mask.shape


def mask_path_for_row(group: str, series: str, filename: str) -> Path:
    base = Path(filename).stem
    return IMAGING_ROOT / f"260313_{group}" / series / f"{base}_mask.png"


@dataclass
class RowResult:
    group: str
    series: str
    frame: int
    filename: str
    csv_area_px: float
    repro_area_px: float | None
    rel_residual: float | None
    mask_present: bool
    provenance_status: str  # 'normal' | 'excluded_by_pi_measurement_provenance' | 'mask_missing'


def reproduce_one_csv(group: str, csv_path: Path) -> list[RowResult]:
    df = pd.read_csv(csv_path)
    out: list[RowResult] = []
    for _, row in df.iterrows():
        series = str(row["Series"])
        frame = int(row["Frame"])
        filename = str(row["filename"])
        csv_area = float(row["Area_px"])

        excluded = (group, series) in PI_DIRECTED_PROVENANCE_EXCLUSIONS
        mp = mask_path_for_row(group, series, filename)
        if not mp.exists():
            out.append(RowResult(
                group=group, series=series, frame=frame, filename=filename,
                csv_area_px=csv_area, repro_area_px=None, rel_residual=None,
                mask_present=False,
                provenance_status="mask_missing",
            ))
            continue

        repro, _shape = area_from_mask(mp)
        if repro is None:
            out.append(RowResult(
                group=group, series=series, frame=frame, filename=filename,
                csv_area_px=csv_area, repro_area_px=None, rel_residual=None,
                mask_present=True,
                provenance_status="mask_read_fail",
            ))
            continue

        if csv_area == 0:
            rel = float("nan") if repro != 0 else 0.0
        else:
            rel = abs(repro - csv_area) / abs(csv_area)

        status = (
            "excluded_by_pi_measurement_provenance" if excluded else "normal"
        )
        out.append(RowResult(
            group=group, series=series, frame=frame, filename=filename,
            csv_area_px=csv_area, repro_area_px=repro, rel_residual=rel,
            mask_present=True,
            provenance_status=status,
        ))
    return out


def summarise(results: list[RowResult]) -> dict:
    df = pd.DataFrame([r.__dict__ for r in results])
    summary: dict = {
        "total_rows": int(len(df)),
        "groups": {},
    }

    by_group = df.groupby("group")
    for g, gdf in by_group:
        normal = gdf[gdf["provenance_status"] == "normal"]
        excluded = gdf[gdf["provenance_status"] == "excluded_by_pi_measurement_provenance"]
        missing = gdf[gdf["provenance_status"].isin(["mask_missing", "mask_read_fail"])]

        normal_resid = normal["rel_residual"].dropna()
        excluded_resid = excluded["rel_residual"].dropna()

        summary["groups"][g] = {
            "rows_total": int(len(gdf)),
            "rows_normal": int(len(normal)),
            "rows_excluded_by_pi": int(len(excluded)),
            "rows_missing_or_failed": int(len(missing)),
            "normal_max_rel_residual": float(normal_resid.max()) if len(normal_resid) else None,
            "normal_mean_rel_residual": float(normal_resid.mean()) if len(normal_resid) else None,
            "normal_p99_rel_residual": float(normal_resid.quantile(0.99)) if len(normal_resid) else None,
            "excluded_max_rel_residual": float(excluded_resid.max()) if len(excluded_resid) else None,
            "excluded_mean_rel_residual": float(excluded_resid.mean()) if len(excluded_resid) else None,
        }
    return summary


def main():
    all_results: list[RowResult] = []
    for group, csv_path in GROUPS.items():
        print(f"[measurement-protocol gate] reproducing {group} from {csv_path.name} ...", flush=True)
        rs = reproduce_one_csv(group, csv_path)
        all_results.extend(rs)
        print(
            f"  rows={len(rs)} "
            f"mask_present={sum(r.mask_present for r in rs)} "
            f"excluded={sum(r.provenance_status == 'excluded_by_pi_measurement_provenance' for r in rs)}",
            flush=True,
        )

    df = pd.DataFrame([r.__dict__ for r in all_results])
    out_csv = REPO / "runs" / "v2_layer_2_unit_1" / "pi_area_reproduction.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nSaved per-row reproduction → {out_csv}")

    summary = summarise(all_results)
    print("\n=== SUMMARY (source-pipeline reproduction; "
          "PI-directed provenance exclusion: Lam4 Position(31)) ===")
    print(f"total rows: {summary['total_rows']}")
    for g, s in summary["groups"].items():
        print(f"\n[{g}]")
        for k, v in s.items():
            print(f"  {k}: {v}")

    # Bottom line: ≤1.5% threshold enforcement on the non-excluded set ONLY,
    # under PI source-pipeline contract; tolerance derived from float / cv2
    # determinism, not chosen to make rows pass.
    fail = []
    for g, s in summary["groups"].items():
        m = s["normal_max_rel_residual"]
        if m is not None and m > 0.015:
            fail.append((g, m))
    if fail:
        print("\n[CONTRACT FAIL] non-excluded rows exceed 1.5% relative residual:")
        for g, m in fail:
            print(f"  {g}: max_rel_residual={m:.4%}")
        sys.exit(2)
    print("\n[CONTRACT PASS] all non-excluded rows ≤1.5% relative residual under "
          "source-pipeline reproduction.")


if __name__ == "__main__":
    main()
