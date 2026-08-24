"""Dependency-light streaming reader for the pinned Karacosta Source Data XLSX.

The workbook is 79 MB and contains 90,066 raw CyTOF events.  Loading it as a
generic office document creates unnecessary object overhead, so this module
reads only the OOXML parts required by the corpus builder.  It is deliberately
read-only and accepts no formulas, macros, or external relationships.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
DOC_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
NS = {"m": MAIN_NS, "r": DOC_REL_NS}
REL_NS = {"p": PKG_REL_NS}
CELL = f"{{{MAIN_NS}}}c"
ROW = f"{{{MAIN_NS}}}row"
MAX_SHARED_STRINGS = 10_000
MAX_SHEETS = 32
MAX_COLUMNS = 128
MAX_ROWS = 200_000


def _column_index(reference: str) -> int:
    match = re.match(r"[A-Z]+", reference)
    if match is None:
        raise ValueError(f"invalid OOXML cell reference: {reference}")
    result = 0
    for char in match.group(0):
        result = result * 26 + ord(char) - 64
    return result - 1


def _safe_member(name: str) -> str:
    normalized = name.lstrip("/")
    if normalized.startswith("../") or "/../" in normalized:
        raise ValueError("unsafe OOXML relationship target")
    return normalized


@dataclass(frozen=True)
class SheetInfo:
    name: str
    member: str
    max_row: int
    max_column: int


class Workbook:
    def __init__(self, path: Path):
        self.path = path
        with zipfile.ZipFile(path) as archive:
            members = set(archive.namelist())
            required = {
                "xl/workbook.xml", "xl/_rels/workbook.xml.rels",
                "xl/sharedStrings.xml",
            }
            if not required.issubset(members):
                raise ValueError("Karacosta workbook OOXML structure changed")
            if any(name.lower().endswith("vbaproject.bin") for name in members):
                raise ValueError("macro-enabled workbook is not accepted")
            self.shared_strings = self._read_shared_strings(archive)
            self.sheets = self._read_sheets(archive, members)

    @staticmethod
    def _read_shared_strings(archive: zipfile.ZipFile) -> tuple[str, ...]:
        root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        values = tuple(
            "".join(node.text or "" for node in item.findall(".//m:t", NS))
            for item in root.findall("m:si", NS)
        )
        if len(values) > MAX_SHARED_STRINGS:
            raise ValueError("unexpectedly large shared-string table")
        return values

    @staticmethod
    def _read_sheets(
        archive: zipfile.ZipFile, members: set[str]
    ) -> dict[str, SheetInfo]:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(
            archive.read("xl/_rels/workbook.xml.rels")
        )
        targets = {
            relation.attrib["Id"]: _safe_member(relation.attrib["Target"])
            for relation in relationships.findall("p:Relationship", REL_NS)
        }
        sheets: dict[str, SheetInfo] = {}
        for sheet in workbook.findall("m:sheets/m:sheet", NS):
            name = sheet.attrib["name"]
            relation_id = sheet.attrib[f"{{{DOC_REL_NS}}}id"]
            target = targets[relation_id]
            member = target if target.startswith("xl/") else f"xl/{target}"
            if member not in members:
                raise ValueError(f"missing worksheet part for {name}")
            dimension_ref = None
            with archive.open(member) as part:
                for _, element in ET.iterparse(part, events=("start",)):
                    if element.tag == f"{{{MAIN_NS}}}dimension":
                        dimension_ref = element.attrib.get("ref")
                        break
                    if element.tag == f"{{{MAIN_NS}}}sheetData":
                        break
            if dimension_ref is None:
                max_row = max_column = 0
            else:
                final = dimension_ref.split(":")[-1]
                row_match = re.search(r"[0-9]+$", final)
                if row_match is None:
                    raise ValueError(f"invalid worksheet dimension for {name}")
                max_row = int(row_match.group(0))
                max_column = _column_index(final) + 1
            if max_row > MAX_ROWS or max_column > MAX_COLUMNS:
                raise ValueError(f"worksheet limits exceeded for {name}")
            sheets[name] = SheetInfo(name, member, max_row, max_column)
        if not sheets or len(sheets) > MAX_SHEETS:
            raise ValueError("unexpected worksheet count")
        return sheets

    def _value(self, cell: ET.Element) -> str | None:
        if cell.find("m:f", NS) is not None:
            raise ValueError("formula cells are not accepted as source observations")
        kind = cell.attrib.get("t")
        if kind == "inlineStr":
            return "".join(
                node.text or "" for node in cell.findall(".//m:t", NS)
            )
        value = cell.find("m:v", NS)
        if value is None:
            return None
        if kind == "s":
            return self.shared_strings[int(value.text or "0")]
        return value.text

    def iter_rows(self, sheet_name: str):
        info = self.sheets[sheet_name]
        with zipfile.ZipFile(self.path) as archive, archive.open(info.member) as part:
            for _, element in ET.iterparse(part, events=("end",)):
                if element.tag != ROW:
                    continue
                row_index = int(element.attrib["r"])
                cells: dict[int, str | None] = {}
                for cell in element:
                    if cell.tag != CELL:
                        continue
                    column = _column_index(cell.attrib["r"])
                    if column >= MAX_COLUMNS:
                        raise ValueError("worksheet column limit exceeded")
                    cells[column] = self._value(cell)
                yield row_index, cells
                element.clear()

    def records_after_header(
        self, sheet_name: str, first_header: str,
        unnamed_first_column: str | None = None,
    ):
        header: list[str | None] | None = None
        header_row = 0
        for row_index, cells in self.iter_rows(sheet_name):
            width = max(cells, default=-1) + 1
            values = [cells.get(index) for index in range(width)]
            if header is None and first_header in values:
                header = values
                if header and header[0] is None and unnamed_first_column:
                    header[0] = unnamed_first_column
                header_row = row_index
                continue
            if header is None or row_index <= header_row or not cells:
                continue
            yield row_index, {
                name: cells.get(index)
                for index, name in enumerate(header)
                if name is not None
            }
