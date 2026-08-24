#!/usr/bin/env python3
"""Verify and extract author spreadsheet annotations for Dryad 08kprr59c.

The release uses cell fill as data: green identifies the maximum-density row
used for temporal alignment, while red identifies measurements the authors
excluded.  This audit reads OOXML directly so those labels survive a tabular
export without making a spreadsheet library part of the runtime.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from xml.etree import ElementTree as ET
import zipfile


NS = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
ACCESSION = "10.5061/dryad.08kprr59c"
ARTICLE = "10.1126/sciadv.adn0235"
FILES = {
    "Data_Fig2.xlsx": {
        "file_id": 3_082_702,
        "size": 36_747,
        "sha256": "d8de16a23e53801796a6802895230f1682167cd5e770c4e64e58e2a9732e2d44",
    },
    "Data_Fig4.xlsx": {
        "file_id": 3_082_704,
        "size": 25_014,
        "sha256": "07c9606199fb01cad84e9cfd16aa898fb99aaf054de9fde7e02d16b60e89b471",
    },
    "README.md": {
        "file_id": 3_082_804,
        "size": 4_586,
        "sha256": "770549fcd3daaf6a102db8b76fd05bd091063cf41cbab30850da830630b99766",
    },
}
EXPECTED_GREEN = {
    "N21", "S21", "X21", "AC21", "AH21", "AN21", "AR21", "AV21", "AZ21",
}
EXPECTED_RED_COUNTS = {"Data_Fig2.xlsx": 81, "Data_Fig4.xlsx": 31}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _shared_strings(archive: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return [
        "".join(node.text or "" for node in item.findall(".//x:t", NS))
        for item in root.findall("x:si", NS)
    ]


def _value(cell: ET.Element, shared: list[str]) -> str | None:
    if cell.attrib.get("t") == "inlineStr":
        return "".join(node.text or "" for node in cell.findall(".//x:t", NS))
    node = cell.find("x:v", NS)
    if node is None:
        return None
    if cell.attrib.get("t") == "s":
        return shared[int(node.text)]
    return node.text


def _annotations(path: Path) -> dict[str, list[dict[str, object]]]:
    with zipfile.ZipFile(path) as archive:
        shared = _shared_strings(archive)
        styles = ET.fromstring(archive.read("xl/styles.xml"))
        fills_node = styles.find("x:fills", NS)
        xfs_node = styles.find("x:cellXfs", NS)
        if fills_node is None or xfs_node is None:
            raise ValueError(f"missing OOXML style table: {path.name}")
        fills = list(fills_node)
        xfs = list(xfs_node)
        style_fill = [int(xf.attrib.get("fillId", "0")) for xf in xfs]
        fill_attributes: list[dict[str, str]] = []
        for fill in fills:
            color = fill.find("x:patternFill/x:fgColor", NS)
            fill_attributes.append({} if color is None else dict(color.attrib))
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        result = {"author_excluded": [], "maximum_density_alignment": []}
        for cell in sheet.findall(".//x:c", NS):
            style = int(cell.attrib.get("s", "0"))
            fill = style_fill[style]
            color = fill_attributes[fill]
            record = {
                "cell": cell.attrib["r"], "style_id": style,
                "fill_id": fill, "fill": color, "raw_value": _value(cell, shared),
            }
            if color.get("rgb") == "FFFF0000":
                result["author_excluded"].append(record)
            elif color.get("theme") == "9" and color.get("tint") == "0.39997558519241921":
                result["maximum_density_alignment"].append(record)
        return result


def audit(source_dir: Path) -> dict[str, object]:
    sources: dict[str, object] = {}
    annotations: dict[str, object] = {}
    for name, expected in FILES.items():
        path = source_dir / name
        if not path.is_file() or path.stat().st_size != expected["size"]:
            raise ValueError(f"source size mismatch: {name}")
        digest = _sha256(path)
        if digest != expected["sha256"]:
            raise ValueError(f"source SHA-256 mismatch: {name}: {digest}")
        sources[name] = {
            **expected,
            "api_metadata_url": f"https://datadryad.org/api/v2/files/{expected['file_id']}",
        }
        if name.endswith(".xlsx"):
            annotation = _annotations(path)
            if len(annotation["author_excluded"]) != EXPECTED_RED_COUNTS[name]:
                raise ValueError(f"red author-exclusion annotation count changed: {name}")
            annotations[name] = annotation
    green = {
        item["cell"]
        for item in annotations["Data_Fig2.xlsx"]["maximum_density_alignment"]
    }
    if green != EXPECTED_GREEN:
        raise ValueError(f"maximum-density annotations changed: {sorted(green)}")
    if annotations["Data_Fig4.xlsx"]["maximum_density_alignment"]:
        raise ValueError("unexpected maximum-density annotation in Fig. 4")
    return {
        "schema": "aleph.outer_library.dryad_08kprr59c_source_audit.v1",
        "dataset_id": "dryad-08kprr59c-c2c12-supracontractility",
        "accession": ACCESSION,
        "paper_doi": ARTICLE,
        "provider": "Dryad Digital Repository",
        "license": "cc0-1.0",
        "dataset_version_id": 289_530,
        "dataset_api_url": "https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.08kprr59c",
        "sources": sources,
        "annotations": annotations,
        "author_annotation_contract": {
            "red_fill": "datapoint excluded due to stage motion artifact or inability to link PIV analysis to subsequent frame",
            "green_fill": "maximum cell density used for time-series alignment",
            "style_is_data": True,
        },
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.source_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "green_alignment_cells": len(report["annotations"]["Data_Fig2.xlsx"]["maximum_density_alignment"]),
        "red_excluded_cells": sum(len(v["author_excluded"]) for v in report["annotations"].values()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
