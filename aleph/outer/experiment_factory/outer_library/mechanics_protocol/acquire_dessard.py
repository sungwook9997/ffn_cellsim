#!/usr/bin/env python3
"""Acquire and normalize the Dessard 2024 magnetic-wire source workbook.

The workbook is author-released supplementary data.  This builder downloads the
Europe PMC supplementary archive, checks the publisher-declared MD5 embedded in
the article XML, and parses OOXML with the standard library.  It never rewrites
the source measurements.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

import numpy as np

PMCID = "PMC10929591"
DOI = "10.1039/D4NA00003J"
ARCHIVE_URL = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles"
MEMBER = "NA-006-D4NA00003J-s001.xlsx"
EXPECTED_MD5 = "cedad5da75da899061f8176986c3ba05"
LICENSE = "CC BY-NC 3.0"
CELL_LINES = ("MCF-10A", "MCF-7", "MDA-MB-231")
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
REL_NS = {"r": "http://schemas.openxmlformats.org/package/2006/relationships"}


def _download(url: str) -> bytes:
    with urlopen(url, timeout=120) as response:
        return response.read()


def _column_index(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference)
    if letters is None:
        raise ValueError(f"invalid cell reference {reference}")
    value = 0
    for letter in letters.group(0):
        value = value * 26 + ord(letter) - ord("A") + 1
    return value - 1


def _strings(archive: ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return ["".join(node.text or "" for node in item.findall(".//m:t", NS)) for item in root]


def _sheet_paths(archive: ZipFile) -> dict[str, str]:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {item.attrib["Id"]: item.attrib["Target"] for item in relationships.findall("r:Relationship", REL_NS)}
    result = {}
    for sheet in workbook.findall(".//m:sheet", NS):
        relation = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
        target = targets[relation].lstrip("/")
        result[sheet.attrib["name"]] = target if target.startswith("xl/") else f"xl/{target}"
    return result


def _rows(archive: ZipFile, path: str, shared: list[str]) -> list[list[object | None]]:
    root = ET.fromstring(archive.read(path))
    result: list[list[object | None]] = []
    for row in root.findall(".//m:sheetData/m:row", NS):
        values: list[object | None] = [None] * 15
        for cell in row.findall("m:c", NS):
            index = _column_index(cell.attrib["r"])
            if index >= len(values):
                continue
            value_node = cell.find("m:v", NS)
            if value_node is None or value_node.text is None:
                continue
            raw = value_node.text
            kind = cell.attrib.get("t")
            if kind == "s":
                values[index] = shared[int(raw)]
            elif kind == "str":
                values[index] = raw
            else:
                try:
                    values[index] = float(raw)
                except ValueError:
                    values[index] = raw
        result.append(values)
    return result


def _date_label(value: object | None) -> str | None:
    if isinstance(value, str):
        match = re.search(r"(?:Run\s+)?(\d{1,2}-[A-Za-z]{3}-\d{2})", value)
        if match:
            return datetime.strptime(match.group(1), "%d-%b-%y").date().isoformat()
    if isinstance(value, float) and 40_000 < value < 60_000:
        return (datetime(1899, 12, 30) + timedelta(days=value)).date().isoformat()
    return None


def _normalize(workbook: bytes) -> dict[str, np.ndarray]:
    records: list[dict[str, object]] = []
    # ZipFile accepts a seekable in-memory buffer; keeping it local also makes
    # it impossible for normalization to mutate the downloaded workbook.
    from io import BytesIO
    with ZipFile(BytesIO(workbook)) as archive:
        shared = _strings(archive)
        paths = _sheet_paths(archive)
        if set(paths) != set(CELL_LINES):
            raise ValueError(f"unexpected workbook sheets {sorted(paths)}")
        for label, sheet_name in enumerate(CELL_LINES):
            rows = _rows(archive, paths[sheet_name], shared)
            current_date: str | None = None
            pending: list[int] = []
            sheet_record_indices: list[int] = []
            for row_number, row in enumerate(rows, start=1):
                new_date = _date_label(row[0])
                if new_date is not None:
                    current_date = new_date
                    for index in pending:
                        records[index]["experiment_group"] = current_date
                    pending.clear()
                numeric = []
                for column in (2, 3, 4, 6, 8, 9, 11, 12, 13, 14):
                    value = row[column]
                    numeric.append(float(value) if isinstance(value, (int, float)) else np.nan)
                if not np.isfinite(numeric[:5]).all():
                    continue
                record = {
                    "cell_line": sheet_name,
                    "cell_label": label,
                    "experiment_group": current_date or "pending_date",
                    "source_row": row_number,
                    "measurement_id": str(row[1] or f"row-{row_number}"),
                    "values": numeric,
                }
                records.append(record)
                sheet_record_indices.append(len(records) - 1)
                if current_date is None:
                    pending.append(len(records) - 1)
            if pending:
                raise ValueError(f"{sheet_name} has measurements without an experiment date")

    values = np.asarray([record["values"] for record in records], dtype=np.float64)
    if values.shape != (196, 10):
        raise ValueError(f"Dessard row contract changed: {values.shape}, expected (196, 10)")
    return {
        "cell_line": np.asarray([record["cell_line"] for record in records], dtype="U16"),
        "cell_label": np.asarray([record["cell_label"] for record in records], dtype=np.int8),
        "experiment_group": np.asarray([record["experiment_group"] for record in records], dtype="U10"),
        "source_row": np.asarray([record["source_row"] for record in records], dtype=np.int16),
        "measurement_id": np.asarray([record["measurement_id"] for record in records], dtype="U64"),
        "feature_name": np.asarray([
            "wire_length_um", "wire_diameter_um", "omega_c_rad_s", "wire_aspect_effective",
            "viscosity_pa_s", "viscosity_uncertainty_pa_s", "theta_b_low_rad",
            "theta_b_high_rad", "elastic_modulus_pa", "elastic_modulus_uncertainty_pa",
        ], dtype="U40"),
        "measurement": values,
    }


def acquire(output_dir: Path, archive_bytes: bytes | None = None) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = archive_bytes if archive_bytes is not None else _download(ARCHIVE_URL)
    archive_path = output_dir / f"{PMCID}-supplementary.zip"
    archive_path.write_bytes(payload)
    from io import BytesIO
    with ZipFile(BytesIO(payload)) as archive:
        workbook = archive.read(MEMBER)
    md5 = hashlib.md5(workbook).hexdigest()  # nosec B303: publisher integrity identifier
    if md5 != EXPECTED_MD5:
        raise ValueError(f"source workbook MD5 changed: {md5}")
    workbook_path = output_dir / MEMBER
    workbook_path.write_bytes(workbook)
    tensor = _normalize(workbook)
    tensor_path = output_dir / "dessard2024_mrs_wire_data.npz"
    np.savez_compressed(tensor_path, **tensor)
    receipt = {
        "schema": "aleph.outer_library.mechanics_source_receipt.v1",
        "pmcid": PMCID,
        "doi": DOI,
        "license": LICENSE,
        "source_url": ARCHIVE_URL,
        "archive_member": MEMBER,
        "publisher_declared_md5": EXPECTED_MD5,
        "workbook_sha256": hashlib.sha256(workbook).hexdigest(),
        "normalized_row_count": int(len(tensor["cell_label"])),
        "cell_line_counts": {name: int(np.sum(tensor["cell_line"] == name)) for name in CELL_LINES},
        "measurement_columns": tensor["feature_name"].tolist(),
        "raw_measurements_modified": False,
        "aleph_authority": "none",
    }
    (output_dir / "source_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(acquire(args.output_dir), sort_keys=True))


if __name__ == "__main__":
    main()
