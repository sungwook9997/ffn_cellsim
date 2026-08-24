#!/usr/bin/env python3
"""Normalize verified cellular-nematic PIV and paired TFM/MSM observations."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from schema import Observation


DATASET = "dataverse-data2772-nih3t3-cellular-nematics"
LAB = "ibec-trepat-guillamat-cellular-nematics"
ACCESSION = "10.34810/DATA2772"
PAPER = "10.1126/science.adz9174"
PIV_SHA256 = "d036631d56037a7197047967a1e94d3a755cf806561715e275b3b7bff1b2d71e"
TFM_MD5 = "4dbe7f6bd1187c1c6d82f4cf608d2473"


def _observation(
    *, source_sha256: str, locator: str, sample_id: str,
    biological_replicate: str, technical_replicate: str, condition: str,
    state: str, modality: str, observable: str, value: float, unit: str,
    unit_locator: str, coordinate_name: str, coordinate_value: float,
    coordinate_unit: str,
) -> Observation:
    identity = hashlib.sha256(
        f"{source_sha256}:{locator}:{observable}".encode()
    ).hexdigest()[:20]
    row = Observation(
        observation_id=f"{DATASET}:{identity}", dataset_id=DATASET,
        lab_group=LAB, provider="CORA Research Data Repository Dataverse",
        accession=ACCESSION, license_id="cc0-1.0",
        source_sha256=source_sha256, source_locator=locator,
        sample_id=sample_id, biological_replicate=biological_replicate,
        technical_replicate=technical_replicate, condition=condition,
        cell_type="NIH3T3_mouse_fibroblast", cell_state=state,
        modality=modality, observable=observable, value=float(value), unit=unit,
        unit_source_sha256=source_sha256, unit_source_locator=unit_locator,
        coordinate_name=coordinate_name, coordinate_value=coordinate_value,
        coordinate_unit=coordinate_unit,
    )
    row.validate()
    return row


def _load_receipt(receipt_path: Path) -> dict[str, object]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("schema") != "aleph.outer_library.dataverse_data2772_acquisition.v1":
        raise ValueError("cellular-nematic acquisition schema changed")
    if receipt.get("accession") != ACCESSION:
        raise ValueError("cellular-nematic accession changed")
    if receipt["metadata"].get("license") != "cc0-1.0":
        raise ValueError("cellular-nematic source is no longer verified CC0")
    sources = receipt["sources"]
    if sources["Fig1G_PIVdata.rar"]["sha256"] != PIV_SHA256:
        raise ValueError("PIV archive hash changed")
    if sources["FigS09.rar"]["md5"] != TFM_MD5:
        raise ValueError("TFM/MSM archive hash changed")
    return receipt


def _piv_rows(
    source_root: Path, source_sha256: str, tensor_output: Path | None,
) -> tuple[list[Observation], dict[str, object]]:
    piv_root = source_root / "Fig1G_PIVdata" / "PIV_CollA"
    files = sorted(
        piv_root.glob("PIV_*.txt"), key=lambda path: int(path.stem.split("_")[1])
    )
    if [int(path.stem.split("_")[1]) for path in files] != list(range(1, 64)):
        raise ValueError("expected exactly PIV_1 through PIV_63")
    tensors: list[np.ndarray] = []
    rows: list[Observation] = []
    flags: dict[int, int] = {}
    for frame_index, path in enumerate(files):
        table = np.loadtxt(path, dtype=np.float64)
        if table.shape != (3_844, 16) or not np.isfinite(table).all():
            raise ValueError(f"unexpected PIV table: {path.name} {table.shape}")
        for flag, count in zip(*np.unique(table[:, 15].astype(int), return_counts=True)):
            flags[int(flag)] = flags.get(int(flag), 0) + int(count)
        # The plugin's own format contract says its vector plot uses columns
        # x,y,ux1,uy1,mag1. ux0/uy0/mag0 are the prior-iteration interpolation,
        # not the plotted result. Flag 999 is a documented interpolation and is
        # retained with provenance; undocumented composite flags remain withheld.
        valid = np.isin(table[:, 15], (0, 999))
        final = table[valid][:, [0, 1, 2, 3, 4, 15]]
        final[:, 2:5] *= 0.66  # 0.66 um/px and one-minute frame interval.
        frame_column = np.full((len(final), 1), frame_index, dtype=np.float64)
        tensors.append(np.concatenate((frame_column, final), axis=1))
        ux, uy, speed = final[:, 2], final[:, 3], final[:, 4]
        metrics = {
            "mean_cell_sheet_speed": (np.mean(speed), "micrometre_per_minute"),
            "median_cell_sheet_speed": (np.median(speed), "micrometre_per_minute"),
            "rms_cell_sheet_speed": (np.sqrt(np.mean(speed ** 2)), "micrometre_per_minute"),
            "flow_polar_coherence": (
                np.hypot(np.mean(ux), np.mean(uy)) / np.mean(speed),
                "dimensionless",
            ),
        }
        for observable, (value, unit) in metrics.items():
            rows.append(_observation(
                source_sha256=source_sha256,
                locator=(f"rar:Fig1G_PIVdata.rar;file=PIV_CollA/{path.name};"
                         "admitted_flags=0,999"),
                sample_id=f"{DATASET}:figure1G-single-sequence",
                biological_replicate="figure1G_biological_identity_not_reported",
                technical_replicate=f"PIV_frame_{frame_index:02d}",
                condition="collagenase_induced_in_plane_contraction_sequence",
                state="detaching_contractile_cellular_nematic",
                modality="particle_image_velocimetry_derived",
                observable=observable, value=float(value), unit=unit,
                unit_locator=(
                    "00_Readme.txt;Fig1G resolution=0.66 um/px;time interval=1 min;"
                    "ImageJ PIV plotted columns ux1,uy1,mag1"
                ),
                coordinate_name="time_from_sequence_start",
                coordinate_value=float(frame_index), coordinate_unit="minute",
            ))
    tensor = np.concatenate(tensors, axis=0)
    admitted = flags.get(0, 0) + flags.get(999, 0)
    if tensor.shape != (admitted, 7):
        raise ValueError("PIV tensor/vector count mismatch")
    if tensor_output is not None:
        tensor_output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            tensor_output, vectors=tensor,
            columns=np.asarray(["frame", "x_px", "y_px", "ux1_um_per_min",
                                "uy1_um_per_min", "speed1_um_per_min", "piv_flag"]),
        )
    return rows, {
        "sequence_count": 1,
        "frame_count": len(files),
        "vectors_total": int(sum(flags.values())),
        "vectors_admitted_flag_0": flags.get(0, 0),
        "vectors_admitted_interpolated_flag_999": flags.get(999, 0),
        "vectors_admitted_total": admitted,
        "vectors_excluded_by_flag": {
            str(flag): count for flag, count in sorted(flags.items())
            if flag not in (0, 999)
        },
        "split_contract": "all 63 frames remain one sequence and one split group",
    }


def _read_stress(path: Path) -> dict[str, float]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "Mask", "T_x (Pa)", "T_y (Pa)", "sigma_1 h (Pa m)",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"TFM/MSM column contract changed: {path}")
    mask = np.asarray([float(row["Mask"]) for row in rows]) == 1.0
    if int(mask.sum()) < 1_000:
        raise ValueError(f"unexpectedly small tissue mask: {path}")
    tx = np.asarray([float(row["T_x (Pa)"]) for row in rows])
    ty = np.asarray([float(row["T_y (Pa)"]) for row in rows])
    tension = np.asarray([float(row["sigma_1 h (Pa m)"]) for row in rows]) * 1_000.0
    traction = np.hypot(tx, ty)
    return {
        "mean_traction_magnitude": float(np.mean(traction[mask])),
        "median_traction_magnitude": float(np.median(traction[mask])),
        "rms_traction_magnitude": float(np.sqrt(np.mean(traction[mask] ** 2))),
        "mean_maximum_principal_tension": float(np.mean(tension[mask])),
        "median_maximum_principal_tension": float(np.median(tension[mask])),
    }


def _tfm_rows(source_root: Path, source_sha256: str) -> tuple[list[Observation], dict[str, object]]:
    root = source_root / "FigS09_verified" / "tfm_msm"
    rows: list[Observation] = []
    phases = ((0, "before_blebbistatin"), (1, "after_blebbistatin_10uM_30min"))
    units = {
        "mean_traction_magnitude": "Pa",
        "median_traction_magnitude": "Pa",
        "rms_traction_magnitude": "Pa",
        "mean_maximum_principal_tension": "mN_per_m",
        "median_maximum_principal_tension": "mN_per_m",
    }
    for fov in range(1, 13):
        folder = root / f"f{fov}"
        params = (folder / "Params_Files.csv").read_text(encoding="utf-8")
        if "20230609_Results" not in params or "n_times,4" not in params:
            raise ValueError(f"Figure S9 paired-FOV contract changed: f{fov}")
        for phase_index, condition in phases:
            matches = list((folder / "Stresses").glob(f"Stresses_*_t_{phase_index}.csv"))
            if len(matches) != 1:
                raise ValueError(f"missing unique paired stress table: f{fov} t{phase_index}")
            metrics = _read_stress(matches[0])
            for observable, value in metrics.items():
                rows.append(_observation(
                    source_sha256=source_sha256,
                    locator=(f"rar:FigS09.rar;nested=2 DEFECTS/TFM & MSM data.rar;"
                             f"folder=f{fov};file={matches[0].name}"),
                    sample_id=f"{DATASET}:20230609-figureS9-monolayer-culture",
                    biological_replicate="20230609_single_culture_or_batch",
                    technical_replicate=f"field_of_view_{fov:02d}",
                    condition=condition,
                    state=("contractile_two_defect_cellular_nematic" if phase_index == 0
                           else "myosin_II_inhibited_two_defect_cellular_nematic"),
                    modality="traction_force_and_monolayer_stress_microscopy_derived",
                    observable=observable, value=value, unit=units[observable],
                    unit_locator=(
                        "FigS09 Stresses CSV header: T_x,T_y in Pa; sigma_1 h in Pa m;"
                        "converted exactly by 1000 Pa m per mN/m"
                    ),
                    coordinate_name="author_figure_S9_phase",
                    coordinate_value=float(phase_index), coordinate_unit="index",
                ))
    # f13-f30 lack the complete masks/stress triplet used for n=12.  t2 is also
    # withheld: its numerical direction contradicts the published before/after
    # histogram and the archive provides no unambiguous role for that table.
    return rows, {
        "technical_fields_of_view": 12,
        "biological_culture_or_batch_count": 1,
        "before_table_index": 0,
        "after_table_index": 1,
        "withheld_table_index": 2,
        "withheld_reason": (
            "t0 and t1 reproduce Figure S9 before/after distributions; t2 does not, "
            "and its role is not unambiguously documented"
        ),
        "f13_to_f30_withheld_reason": (
            "not members of the complete n=12 paired mask/stress set used in Figure S9"
        ),
        "split_contract": (
            "12 fields are technical replicates nested in one dated culture/batch and "
            "must never cross train/test splits"
        ),
    }


def build(
    source_root: Path, receipt_path: Path, tensor_output: Path | None = None,
) -> tuple[list[Observation], dict[str, object]]:
    receipt = _load_receipt(receipt_path)
    piv_sha = receipt["sources"]["Fig1G_PIVdata.rar"]["sha256"]
    tfm_sha = receipt["sources"]["FigS09.rar"]["sha256"]
    piv, piv_audit = _piv_rows(source_root, piv_sha, tensor_output)
    tfm, tfm_audit = _tfm_rows(source_root, tfm_sha)
    rows = piv + tfm
    if len(rows) != 372 or len({row.observation_id for row in rows}) != len(rows):
        raise ValueError("unexpected cellular-nematic observation cardinality")
    report = {
        "schema": "aleph.outer_library.dataverse_data2772_manifest.v1",
        "dataset_id": DATASET, "accession": ACCESSION, "paper_doi": PAPER,
        "observations": len(rows),
        "samples": len({row.sample_id for row in rows}),
        "biological_replicates_claimed": 0,
        "piv": piv_audit, "tfm_msm": tfm_audit,
        "modalities": sorted({row.modality for row in rows}),
        "aleph_authority": "none",
    }
    return rows, report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--tensor-output", type=Path)
    args = parser.parse_args()
    rows, report = build(args.source_root, args.receipt, args.tensor_output)
    args.output.write_text(
        "".join(json.dumps(row.as_dict(), sort_keys=True) + "\n" for row in rows),
        encoding="utf-8",
    )
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
