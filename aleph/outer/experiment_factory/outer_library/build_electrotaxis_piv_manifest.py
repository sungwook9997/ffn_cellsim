#!/usr/bin/env python3
"""Normalize real epithelial-tissue PIV fields into frame-level observables."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from schema import Observation


ACCESSION = "10.5281/zenodo.5880000"
DATASET_ID = "zenodo-5880000-mdck-electrotaxis-piv"
OFFICIAL_MD5 = {
    "control": "7350858b9bbf1f40aed956d3128e1aa0",
    "stimulated": "a72dad4ee010500eeb70e9b916adc196",
}
OFFICIAL_SIZE = {
    "control": 153_729_717,
    "stimulated": 159_083_657,
}


def _digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _field_metrics(u: np.ndarray, v: np.ndarray) -> dict[str, tuple[float, str]]:
    u = np.asarray(u, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    if u.shape != v.shape or u.ndim != 2:
        raise ValueError(f"invalid PIV component shapes: {u.shape}, {v.shape}")
    valid = np.isfinite(u) & np.isfinite(v)
    if not valid.any():
        raise ValueError("PIV frame has no finite vectors")
    uv = u[valid]
    vv = v[valid]
    speed = np.hypot(uv, vv)
    moving = speed > 0
    if moving.any():
        unit_x = uv[moving] / speed[moving]
        unit_y = vv[moving] / speed[moving]
        polar_order = np.hypot(unit_x.mean(), unit_y.mean())
        nematic_real = np.mean(unit_x**2 - unit_y**2)
        nematic_imag = np.mean(2 * unit_x * unit_y)
        nematic_order = np.hypot(nematic_real, nematic_imag)
    else:
        polar_order = 0.0
        nematic_order = 0.0
    return {
        "mean_velocity_x": (float(uv.mean()), "source_velocity_unit_unknown"),
        "mean_velocity_y": (float(vv.mean()), "source_velocity_unit_unknown"),
        "mean_speed": (float(speed.mean()), "source_velocity_unit_unknown"),
        "rms_speed": (float(np.sqrt(np.mean(speed**2))), "source_velocity_unit_unknown"),
        "speed_p90": (float(np.quantile(speed, 0.9)), "source_velocity_unit_unknown"),
        "velocity_polar_order": (float(polar_order), "dimensionless"),
        "velocity_nematic_order": (float(nematic_order), "dimensionless"),
        "valid_vector_fraction": (float(valid.mean()), "fraction"),
    }


def build(source_mat: Path, condition: str) -> list[Observation]:
    if condition not in OFFICIAL_MD5:
        raise ValueError(f"unsupported author condition: {condition}")
    if source_mat.stat().st_size != OFFICIAL_SIZE[condition]:
        raise ValueError(f"official size mismatch for {source_mat}")
    if _digest(source_mat, "md5") != OFFICIAL_MD5[condition]:
        raise ValueError(f"official MD5 mismatch for {source_mat}")
    sha256 = _digest(source_mat, "sha256")
    payload = loadmat(source_mat, squeeze_me=True, struct_as_record=False)
    frames = np.asarray(payload["S"], dtype=object).reshape(-1)
    if len(frames) != 60:
        raise ValueError(f"expected 60 PIV frames, found {len(frames)}")
    rows: list[Observation] = []
    for frame_index, frame in enumerate(frames):
        metrics = _field_metrics(frame.u_filtered, frame.v_filtered)
        for metric, (value, unit) in metrics.items():
            locator = (
                f"mat:S[{frame_index}].u_filtered+S[{frame_index}].v_filtered::{metric}"
            )
            obs_digest = hashlib.sha256(f"{sha256}:{locator}".encode()).hexdigest()[:20]
            row = Observation(
                observation_id=f"electrotaxis-piv:{obs_digest}",
                dataset_id=DATASET_ID,
                lab_group="zenodo-5880000-authors",
                provider="Zenodo",
                accession=ACCESSION,
                license_id="cc-by-4.0",
                source_sha256=sha256,
                source_locator=locator,
                sample_id=(
                    f"{DATASET_ID}:{condition}:representative-tissue:frame-{frame_index}"
                ),
                biological_replicate="representative_tissue_single",
                technical_replicate="not_reported_single_representative_tissue",
                condition=condition,
                cell_type="MDCK-II epithelial tissue",
                cell_state=(
                    "three_hour_electrical_stimulation" if condition == "stimulated"
                    else "unstimulated_control"
                ),
                modality="particle_image_velocimetry_derived",
                observable=metric,
                value=value,
                unit=unit,
                coordinate_name="frame_index",
                coordinate_value=float(frame_index),
                coordinate_unit="index",
            )
            row.validate()
            rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--control", type=Path, required=True)
    parser.add_argument("--stimulated", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build(args.control, "control") + build(args.stimulated, "stimulated")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")
    print(json.dumps({"observations": len(rows), "frames": len({r.sample_id for r in rows})}))


if __name__ == "__main__":
    main()
