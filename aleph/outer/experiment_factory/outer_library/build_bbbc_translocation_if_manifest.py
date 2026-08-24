#!/usr/bin/env python3
"""Build field-grouped IF features from the public BBBC013/BBBC014 plates.

Each well remains one sample.  The two channels are paired before any feature is
emitted, so neither pixels nor channels can leak across a split.  The adapter
uses only author-provided plate labels; image-derived measurements are
observables, never cell-type or cell-state labels.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import struct
from typing import Iterable
from zipfile import ZipFile

import numpy as np
from scipy import ndimage

from schema import Observation


BBBC013_ARCHIVE_SHA256 = (
    "c059b569d96f70ad5626fad144867e6ece4353622119c46a8af8f9794f1e7985"
)
BBBC013_PLATEMAP_SHA256 = (
    "e8db6666271d47962fa7d2abfa3ea965352b8e87bee461f2983d0f667bc7ff08"
)
BBBC014_ARCHIVE_SHA256 = (
    "223c3230c164c8295f4626f82e4bd9fe394d8e9339a4f82faf697fdd9b835341"
)
BBBC014_PLATEMAP_SHA256 = (
    "c3c4e645f9d16d2d8f0304038824cb7daccfd3a385cf10fce2a1b8ee395be266"
)

FEATURE_NAMES = (
    "signal_intensity_mean",
    "signal_intensity_std",
    "signal_intensity_p90",
    "signal_foreground_fraction",
    "nuclear_object_count",
    "nuclear_mask_fraction",
    "nuclear_object_median_area",
    "signal_nuclear_mean",
    "signal_perinuclear_mean",
    "signal_nuclear_perinuclear_log_ratio",
    "signal_nuclear_perinuclear_standardized_difference",
    "dna_signal_foreground_pearson",
)

FEATURE_UNITS = {
    "signal_intensity_mean": "uint8_intensity",
    "signal_intensity_std": "uint8_intensity",
    "signal_intensity_p90": "uint8_intensity",
    "signal_foreground_fraction": "fraction",
    "nuclear_object_count": "count",
    "nuclear_mask_fraction": "fraction",
    "nuclear_object_median_area": "pixel2",
    "signal_nuclear_mean": "uint8_intensity",
    "signal_perinuclear_mean": "uint8_intensity",
    "signal_nuclear_perinuclear_log_ratio": "dimensionless",
    "signal_nuclear_perinuclear_standardized_difference": "dimensionless",
    "dna_signal_foreground_pearson": "dimensionless",
}


@dataclass(frozen=True, slots=True)
class PlateSpec:
    dataset_id: str
    accession: str
    archive_sha256: str
    platemap_sha256: str
    provider: str
    lab_group: str
    image_shape: tuple[int, int]
    signal_channel: int
    dna_channel: int
    modality: str


SPECS = {
    "BBBC013": PlateSpec(
        dataset_id="bbbc013-v1-fkhr-translocation",
        accession="BBBC013v1",
        archive_sha256=BBBC013_ARCHIVE_SHA256,
        platemap_sha256=BBBC013_PLATEMAP_SHA256,
        provider="Ilya Ravkin via Broad Bioimage Benchmark Collection",
        # Shared provider lineage is deliberately one lab group.  Acquisition
        # site/instrument differences do not establish independent provenance.
        lab_group="ravkin-bbbc-translocation",
        image_shape=(640, 640),
        signal_channel=1,
        dna_channel=2,
        modality="IF_FKHR_EGFP_DRAQ_field",
    ),
    "BBBC014": PlateSpec(
        dataset_id="bbbc014-v1-nfkb-translocation",
        accession="BBBC014v1",
        archive_sha256=BBBC014_ARCHIVE_SHA256,
        platemap_sha256=BBBC014_PLATEMAP_SHA256,
        provider="Ilya Ravkin via Broad Bioimage Benchmark Collection",
        lab_group="ravkin-bbbc-translocation",
        image_shape=(1024, 1360),
        signal_channel=2,
        dna_channel=1,
        modality="IF_NFkB_FITC_DAPI_field",
    ),
}


@dataclass(frozen=True, slots=True)
class WellLabel:
    row: str
    column: int
    image_index: int
    treatment: str
    dose_text: str
    dose_value: float
    dose_unit: str
    cell_type: str
    cell_state: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def decode_bmp(payload: bytes) -> np.ndarray:
    """Decode an uncompressed 8/24/32-bit BMP into one uint8 intensity plane."""
    if len(payload) < 54 or payload[:2] != b"BM":
        raise ValueError("not a BMP image")
    pixel_offset = struct.unpack_from("<I", payload, 10)[0]
    dib_size = struct.unpack_from("<I", payload, 14)[0]
    if dib_size < 40:
        raise ValueError(f"unsupported BMP DIB header size {dib_size}")
    width, signed_height = struct.unpack_from("<ii", payload, 18)
    planes, bits_per_pixel = struct.unpack_from("<HH", payload, 26)
    compression = struct.unpack_from("<I", payload, 30)[0]
    if width <= 0 or signed_height == 0 or planes != 1 or compression != 0:
        raise ValueError("only uncompressed, positive-width BMP planes are supported")
    height = abs(signed_height)
    channels = {8: 1, 24: 3, 32: 4}.get(bits_per_pixel)
    if channels is None:
        raise ValueError(f"unsupported BMP bit depth {bits_per_pixel}")
    row_bytes = ((width * channels + 3) // 4) * 4
    required = pixel_offset + row_bytes * height
    if required > len(payload):
        raise ValueError("truncated BMP pixel array")
    packed = np.frombuffer(payload, dtype=np.uint8, count=row_bytes * height, offset=pixel_offset)
    packed = packed.reshape(height, row_bytes)[:, : width * channels]
    if signed_height > 0:
        packed = packed[::-1]
    if bits_per_pixel == 8:
        indexes = packed[:, :width]
        palette_start = 14 + dib_size
        palette_bytes = pixel_offset - palette_start
        if palette_bytes >= 4:
            palette = np.frombuffer(
                payload,
                dtype=np.uint8,
                count=(palette_bytes // 4) * 4,
                offset=palette_start,
            ).reshape(-1, 4)
            if int(indexes.max(initial=0)) >= len(palette):
                raise ValueError("BMP palette index out of range")
            # BMP palettes are B,G,R,0.  Preserve a grayscale palette exactly;
            # otherwise return displayed luminance for a single-channel assay.
            bgr = palette[:, :3].astype(np.float64)
            if np.array_equal(bgr[:, 0], bgr[:, 1]) and np.array_equal(bgr[:, 1], bgr[:, 2]):
                lookup = palette[:, 0]
            else:
                lookup = np.rint(0.114 * bgr[:, 0] + 0.587 * bgr[:, 1] + 0.299 * bgr[:, 2]).astype(np.uint8)
            return np.asarray(lookup[indexes], dtype=np.uint8)
        return np.asarray(indexes, dtype=np.uint8)
    bgr = packed.reshape(height, width, channels)[..., :3].astype(np.float64)
    return np.rint(0.114 * bgr[..., 0] + 0.587 * bgr[..., 1] + 0.299 * bgr[..., 2]).astype(np.uint8)


def _otsu(channel: np.ndarray) -> float:
    values = np.asarray(channel, dtype=np.uint8).reshape(-1)
    hist = np.bincount(values, minlength=256).astype(np.float64)
    probability = hist / max(hist.sum(), 1.0)
    omega = np.cumsum(probability)
    means = np.cumsum(probability * np.arange(256, dtype=np.float64))
    denominator = omega * (1.0 - omega)
    between = np.zeros(256, dtype=np.float64)
    valid = denominator > 0
    between[valid] = (means[-1] * omega[valid] - means[valid]) ** 2 / denominator[valid]
    return float(np.argmax(between))


def field_features(dna: np.ndarray, signal: np.ndarray) -> dict[str, float]:
    """Return scale-aware field features without manufacturing cell labels."""
    dna = np.asarray(dna)
    signal = np.asarray(signal)
    if dna.shape != signal.shape or dna.ndim != 2:
        raise ValueError("paired DNA and signal images must be equally shaped 2D planes")
    if dna.dtype != np.uint8 or signal.dtype != np.uint8:
        raise ValueError("paired IF planes must be uint8")

    dna_threshold = _otsu(dna)
    raw_nuclei = dna > dna_threshold
    raw_nuclei = ndimage.binary_opening(raw_nuclei, structure=np.ones((3, 3), dtype=bool))
    labels, count = ndimage.label(raw_nuclei)
    areas = np.bincount(labels.reshape(-1), minlength=count + 1)
    min_area = max(9, int(round(dna.size * 1.0e-5)))
    # Five percent is permissive enough for low-cell-count fields and synthetic
    # QA images while still removing whole-field illumination islands.
    max_area = max(min_area + 1, int(round(dna.size * 5.0e-2)))
    keep = (areas >= min_area) & (areas <= max_area)
    keep[0] = False
    border = np.unique(
        np.concatenate((labels[0], labels[-1], labels[:, 0], labels[:, -1]))
    )
    keep[border] = False
    nuclei = keep[labels]
    kept_areas = areas[keep]
    if not nuclei.any():
        raise ValueError("DNA plane produced no valid nuclear foreground objects")

    equivalent_radius = float(np.median(np.sqrt(kept_areas / np.pi)))
    ring_width = float(np.clip(0.75 * equivalent_radius, 2.0, 20.0))
    distance = ndimage.distance_transform_edt(~nuclei)
    perinuclear = (distance > 0.0) & (distance <= ring_width)
    signal_f = signal.astype(np.float64)
    dna_f = dna.astype(np.float64)
    nuclear_mean = float(signal_f[nuclei].mean())
    perinuclear_mean = float(signal_f[perinuclear].mean())
    scale = float(signal_f.std())
    signal_threshold = _otsu(signal)
    joint = nuclei | perinuclear
    if joint.any() and dna_f[joint].std() > 0 and signal_f[joint].std() > 0:
        pearson = float(np.corrcoef(dna_f[joint], signal_f[joint])[0, 1])
    else:
        pearson = 0.0
    eps = 0.5  # half of one uint8 intensity bin
    values = {
        "signal_intensity_mean": float(signal_f.mean()),
        "signal_intensity_std": scale,
        "signal_intensity_p90": float(np.quantile(signal_f, 0.90)),
        "signal_foreground_fraction": float(np.mean(signal > signal_threshold)),
        "nuclear_object_count": float(len(kept_areas)),
        "nuclear_mask_fraction": float(nuclei.mean()),
        "nuclear_object_median_area": float(np.median(kept_areas)),
        "signal_nuclear_mean": nuclear_mean,
        "signal_perinuclear_mean": perinuclear_mean,
        "signal_nuclear_perinuclear_log_ratio": float(
            np.log((nuclear_mean + eps) / (perinuclear_mean + eps))
        ),
        "signal_nuclear_perinuclear_standardized_difference": float(
            (nuclear_mean - perinuclear_mean) / max(scale, eps)
        ),
        "dna_signal_foreground_pearson": pearson,
    }
    if tuple(values) != FEATURE_NAMES or not np.isfinite(list(values.values())).all():
        raise ValueError("invalid IF feature vector")
    return values


def _plate_values(path: Path, expected_sha256: str) -> list[str]:
    if _sha256(path) != expected_sha256:
        raise ValueError(f"plate map content hash mismatch: {path}")
    values = [item.strip() for item in path.read_text(encoding="utf-8").splitlines() if item.strip()]
    if not values or not values[0].startswith("DESCRIPTION") or len(values[1:]) != 96:
        raise ValueError(f"expected one 96-well plate map: {path}")
    return values[1:]


def plate_labels(accession: str, platemap: Path) -> list[WellLabel]:
    spec = SPECS[accession]
    values = _plate_values(platemap, spec.platemap_sha256)
    labels: list[WellLabel] = []
    for index, dose_text in enumerate(values, start=1):
        row = chr(ord("A") + (index - 1) // 12)
        column = (index - 1) % 12 + 1
        if accession == "BBBC013":
            treatment = "Wortmannin" if row <= "D" else "LY294002"
            cell_type = "U2OS human osteosarcoma cells"
            dose_value = float(dose_text)
            dose_unit = "nM"
            cell_state = "PI3K_PKB_inhibition_dose_response"
        else:
            treatment = "TNF-alpha"
            cell_type = (
                "MCF7 human breast adenocarcinoma cells"
                if row <= "D"
                else "A549 human alveolar basal epithelial cells"
            )
            if not dose_text.startswith("10^"):
                raise ValueError(f"unexpected BBBC014 dose label {dose_text!r}")
            dose_value = float(10.0 ** float(dose_text[3:]))
            dose_unit = "g/mL"
            cell_state = "TNF_alpha_dose_response"
        labels.append(
            WellLabel(
                row=row,
                column=column,
                image_index=index,
                treatment=treatment,
                dose_text=dose_text,
                dose_value=dose_value,
                dose_unit=dose_unit,
                cell_type=cell_type,
                cell_state=cell_state,
            )
        )
    return labels


_MEMBER_PATTERNS = {
    "BBBC013": re.compile(
        r"(?:^|/)Channel(?P<channel>[12])-(?P<index>\d{2})-"
        r"(?P<row>[A-H])-(?P<column>\d{2})\.BMP$",
        re.IGNORECASE,
    ),
    "BBBC014": re.compile(
        r"(?:^|/)Channel (?P<channel>[12])-(?P<index>\d{2})-"
        r"(?P<row>[A-H])-(?P<column>\d{2})-00\.Bmp$",
        re.IGNORECASE,
    ),
}


def _members(archive: ZipFile, accession: str) -> dict[tuple[int, int], str]:
    result: dict[tuple[int, int], str] = {}
    pattern = _MEMBER_PATTERNS[accession]
    for name in archive.namelist():
        match = pattern.search(name)
        if not match:
            continue
        channel = int(match.group("channel"))
        index = int(match.group("index"))
        row = match.group("row").upper()
        column = int(match.group("column"))
        expected = (ord(row) - ord("A")) * 12 + column
        if expected != index:
            raise ValueError(f"archive member index disagrees with well: {name}")
        key = (index, channel)
        if key in result:
            raise ValueError(f"duplicate archive member for well/channel: {key}")
        result[key] = name
    if len(result) != 192:
        raise ValueError(f"expected 192 paired-channel images, found {len(result)}")
    return result


def build_plate(
    accession: str,
    source_zip: Path,
    platemap: Path,
) -> tuple[list[Observation], dict[str, np.ndarray]]:
    spec = SPECS[accession]
    if _sha256(source_zip) != spec.archive_sha256:
        raise ValueError(f"{accession} image archive content hash mismatch")
    labels = plate_labels(accession, platemap)
    rows: list[Observation] = []
    matrix: list[list[float]] = []
    sample_ids: list[str] = []
    plate_rows: list[str] = []
    columns: list[int] = []
    treatments: list[str] = []
    cell_types: list[str] = []
    doses: list[float] = []
    dose_texts: list[str] = []
    with ZipFile(source_zip) as archive:
        members = _members(archive, accession)
        for label in labels:
            signal_name = members[(label.image_index, spec.signal_channel)]
            dna_name = members[(label.image_index, spec.dna_channel)]
            signal = decode_bmp(archive.read(signal_name))
            dna = decode_bmp(archive.read(dna_name))
            if signal.shape != spec.image_shape or dna.shape != spec.image_shape:
                raise ValueError(
                    f"{accession} image shape mismatch at {label.row}{label.column:02d}: "
                    f"{signal.shape}, {dna.shape}"
                )
            features = field_features(dna, signal)
            well = f"{label.row}{label.column:02d}"
            sample_id = f"{spec.dataset_id}:plate-1:{well}"
            condition = (
                f"{label.treatment};author_dose={label.dose_text} {label.dose_unit}"
            )
            for feature_name, value in features.items():
                locator = (
                    f"zip:{signal_name}+{dna_name}|platemap:{platemap.name}:"
                    f"sha256={spec.platemap_sha256}:well={well}:"
                    f"value={label.dose_text}|feature:{feature_name}"
                )
                identifier = hashlib.sha256(
                    f"{spec.archive_sha256}:{locator}".encode("utf-8")
                ).hexdigest()[:20]
                observation = Observation(
                    observation_id=f"{accession.lower()}-if:{identifier}",
                    dataset_id=spec.dataset_id,
                    lab_group=spec.lab_group,
                    provider=spec.provider,
                    accession=spec.accession,
                    license_id="cc-by-3.0",
                    source_sha256=spec.archive_sha256,
                    source_locator=locator,
                    sample_id=sample_id,
                    biological_replicate="not_reported_single_plate",
                    technical_replicate=f"author_replica_row_{label.row}",
                    condition=condition,
                    cell_type=label.cell_type,
                    cell_state=label.cell_state,
                    modality=spec.modality,
                    observable=feature_name,
                    value=float(value),
                    unit=FEATURE_UNITS[feature_name],
                    coordinate_name="treatment_dose",
                    coordinate_value=label.dose_value,
                    coordinate_unit=label.dose_unit,
                )
                observation.validate()
                rows.append(observation)
            matrix.append([features[name] for name in FEATURE_NAMES])
            sample_ids.append(sample_id)
            plate_rows.append(label.row)
            columns.append(label.column)
            treatments.append(label.treatment)
            cell_types.append(label.cell_type)
            doses.append(label.dose_value)
            dose_texts.append(label.dose_text)
    tensor = {
        "X": np.asarray(matrix, dtype=np.float64),
        "feature_names": np.asarray(FEATURE_NAMES, dtype="U64"),
        "dataset_id": np.asarray([spec.dataset_id] * 96, dtype="U64"),
        "lab_group": np.asarray([spec.lab_group] * 96, dtype="U64"),
        "sample_id": np.asarray(sample_ids, dtype="U96"),
        "plate_row": np.asarray(plate_rows, dtype="U1"),
        "column": np.asarray(columns, dtype=np.int16),
        "treatment": np.asarray(treatments, dtype="U32"),
        "cell_type": np.asarray(cell_types, dtype="U64"),
        "dose": np.asarray(doses, dtype=np.float64),
        "dose_text": np.asarray(dose_texts, dtype="U24"),
        "source_sha256": np.asarray([spec.archive_sha256] * 96, dtype="U64"),
        "platemap_sha256": np.asarray([spec.platemap_sha256] * 96, dtype="U64"),
    }
    return rows, tensor


def build(
    bbbc013_zip: Path,
    bbbc013_platemap: Path,
    bbbc014_zip: Path,
    bbbc014_platemap: Path,
) -> tuple[list[Observation], dict[str, np.ndarray]]:
    rows_13, tensor_13 = build_plate("BBBC013", bbbc013_zip, bbbc013_platemap)
    rows_14, tensor_14 = build_plate("BBBC014", bbbc014_zip, bbbc014_platemap)
    if tuple(tensor_13) != tuple(tensor_14):
        raise ValueError("IF plate tensor schemas disagree")
    tensor = {
        key: np.concatenate((tensor_13[key], tensor_14[key]), axis=0)
        if tensor_13[key].ndim > 0 and key != "feature_names"
        else tensor_13[key]
        for key in tensor_13
    }
    return rows_13 + rows_14, tensor


def _write_manifest(path: Path, rows: Iterable[Observation]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row.as_dict(), sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bbbc013-zip", type=Path, required=True)
    parser.add_argument("--bbbc013-platemap", type=Path, required=True)
    parser.add_argument("--bbbc014-zip", type=Path, required=True)
    parser.add_argument("--bbbc014-platemap", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tensor-output", type=Path, required=True)
    args = parser.parse_args()
    rows, tensor = build(
        args.bbbc013_zip,
        args.bbbc013_platemap,
        args.bbbc014_zip,
        args.bbbc014_platemap,
    )
    _write_manifest(args.output, rows)
    args.tensor_output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.tensor_output, **tensor)
    print(
        json.dumps(
            {
                "datasets": 2,
                "fields": int(tensor["X"].shape[0]),
                "features": int(tensor["X"].shape[1]),
                "observations": len(rows),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
