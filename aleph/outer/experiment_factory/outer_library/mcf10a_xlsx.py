"""Dependency-free reader for the two pinned LINCS MCF10A XLSX tables."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET
import zipfile


NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def _column_index(reference: str) -> int:
    match = re.match(r"[A-Z]+", reference)
    if match is None:
        raise ValueError(f"invalid XLSX cell reference {reference!r}")
    result = 0
    for character in match.group():
        result = result * 26 + ord(character) - 64
    return result - 1


def read_first_sheet(path: Path) -> tuple[list[str], list[dict[str, Any]]]:
    """Return first-sheet rows while preserving source row numbers."""
    with zipfile.ZipFile(path) as archive:
        strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            strings = [
                "".join(node.text or "" for node in item.findall(".//m:t", NS))
                for item in root.findall("m:si", NS)
            ]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    raw_rows: list[tuple[int, dict[int, Any]]] = []
    for row in sheet.findall(".//m:sheetData/m:row", NS):
        values: dict[int, Any] = {}
        for cell in row.findall("m:c", NS):
            index = _column_index(cell.attrib["r"])
            value_node = cell.find("m:v", NS)
            value: Any = None if value_node is None else value_node.text
            if cell.attrib.get("t") == "s" and value is not None:
                value = strings[int(value)]
            elif cell.attrib.get("t") == "inlineStr":
                value = "".join(
                    node.text or "" for node in cell.findall(".//m:t", NS)
                )
            elif value is not None:
                try:
                    number = float(value)
                    value = int(number) if number.is_integer() else number
                except ValueError:
                    pass
            values[index] = value
        raw_rows.append((int(row.attrib["r"]), values))
    if not raw_rows:
        raise ValueError("empty XLSX sheet")
    width = max(raw_rows[0][1]) + 1
    headers = [str(raw_rows[0][1].get(index, "")) for index in range(width)]
    rows = []
    for source_row, values in raw_rows[1:]:
        if values.get(0) in (None, ""):
            continue
        record = {headers[index]: values.get(index) for index in range(width)}
        record["_source_row"] = source_row
        rows.append(record)
    return headers, rows
