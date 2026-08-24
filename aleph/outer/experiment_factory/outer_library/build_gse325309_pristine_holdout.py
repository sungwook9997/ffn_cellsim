#!/usr/bin/env python3
"""Build the frozen GSE325309 CDH1/VIM holdout without generic office objects."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import numpy as np

from schema import Observation


NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
EXPECTED_HEADER = [
    "gene_id", "chr", "start", "end", "strand", "gene_name", "EntrezID",
    "transcript_id", "Description", "trans_type", "FPKM.con_1", "FPKM.con_2",
    "FPKM.con_3", "FPKM.T_1", "FPKM.T_2", "FPKM.T_3",
]
SAMPLES = [
    ("GSM9600356", "control", "1", 0.0),
    ("GSM9600357", "control", "2", 0.0),
    ("GSM9600358", "control", "3", 0.0),
    ("GSM9600359", "TGFB1_10ng_mL_48h", "1", 48.0),
    ("GSM9600360", "TGFB1_10ng_mL_48h", "2", 48.0),
    ("GSM9600361", "TGFB1_10ng_mL_48h", "3", 48.0),
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _shared_strings(archive: ZipFile) -> list[str]:
    result = []
    with archive.open("xl/sharedStrings.xml") as handle:
        for _, element in ET.iterparse(handle, events=("end",)):
            if element.tag == NS + "si":
                result.append("".join(node.text or "" for node in element.iter(NS + "t")))
                element.clear()
    if len(result) != 124_798:
        raise ValueError("GSE325309 shared-string inventory changed")
    return result


def _column(reference: str) -> str:
    return re.match(r"[A-Z]+", reference).group(0)  # type: ignore[union-attr]


def _rows(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with ZipFile(path) as archive:
        shared = _shared_strings(archive)
        header: list[str] = []
        markers: dict[str, dict[str, str]] = {}
        with archive.open("xl/worksheets/sheet1.xml") as handle:
            for _, row in ET.iterparse(handle, events=("end",)):
                if row.tag != NS + "row":
                    continue
                values: dict[str, str] = {}
                for cell in row.findall(NS + "c"):
                    node = cell.find(NS + "v")
                    if node is None or node.text is None:
                        continue
                    value = shared[int(node.text)] if cell.get("t") == "s" else node.text
                    values[_column(cell.get("r", ""))] = value
                if row.get("r") == "1":
                    header = [values.get(chr(ord("A") + i), "") for i in range(16)]
                gene = values.get("F")
                if gene in {"CDH1", "VIM"}:
                    if gene in markers:
                        raise ValueError(f"duplicate exact gene symbol in GSE325309: {gene}")
                    markers[gene] = values
                row.clear()
    if header != EXPECTED_HEADER or set(markers) != {"CDH1", "VIM"}:
        raise ValueError("GSE325309 frozen workbook layout/gene contract changed")
    expected_ids = {"CDH1": ("ENSG00000039068", "999,51343"), "VIM": ("ENSG00000026025", "7431")}
    for gene, (ensembl, entrez) in expected_ids.items():
        if (markers[gene].get("A"), markers[gene].get("G")) != (ensembl, entrez):
            raise ValueError(f"GSE325309 gene identity changed: {gene}")
    return header, markers


def build(source: Path, receipt_path: Path, protocol_path: Path, manifest_path: Path, tensor_path: Path) -> dict[str, object]:
    receipt = json.loads(receipt_path.read_text())
    protocol = json.loads(protocol_path.read_text())
    if receipt.get("schema") != "aleph.outer_library.gse325309_pristine_acquisition.v1":
        raise ValueError("GSE325309 receipt schema changed")
    if _sha256(source) != receipt["file"]["sha256"] or _sha256(protocol_path) != receipt["freeze_proof"]["protocol_sha256"]:
        raise ValueError("GSE325309 source or frozen protocol changed")
    _, markers = _rows(source)
    columns = ("K", "L", "M", "N", "O", "P")
    genes = ("CDH1", "VIM")
    X = np.asarray([[float(markers[gene][column]) for gene in genes] for column in columns], dtype=np.float64)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8") as handle:
        for sample_index, (sample_id, condition, replicate, hours) in enumerate(SAMPLES):
            for gene_index, gene in enumerate(genes):
                observation = Observation(
                    observation_id=f"GSE325309:{sample_id}:{gene}",
                    dataset_id="geo-gse325309-preregistered-a549-tgfb1-rnaseq",
                    lab_group="Lv_Huang_Jiang_Labs_Sun_Yat_sen_University",
                    provider="NCBI Gene Expression Omnibus", accession="GSE325309",
                    license_id="geo_public_no_explicit_license",
                    source_sha256=receipt["file"]["sha256"],
                    source_locator=f"GSE325309_processed_data.xlsx:sheet1:{gene}:{columns[sample_index]}",
                    sample_id=f"GSE325309:{sample_id}",
                    biological_replicate=f"GSE325309:{condition}:biological_replicate_{replicate}",
                    technical_replicate="none_reported", condition=condition,
                    cell_type="A549_human_lung_adenocarcinoma_cell_line",
                    cell_state="author_TGFB1_conditioned_or_control",
                    modality="bulk_RNA_seq_author_processed", observable=f"{gene}_RNA_expression",
                    value=float(X[sample_index, gene_index]), unit="FPKM",
                    coordinate_name="TGFB1_exposure_time", coordinate_value=hours,
                    coordinate_unit="hour",
                )
                observation.validate()
                handle.write(json.dumps(observation.as_dict(), sort_keys=True) + "\n")
    np.savez_compressed(
        tensor_path, X_FPKM=X, gene=np.asarray(genes),
        sample_id=np.asarray([row[0] for row in SAMPLES]),
        condition=np.asarray([row[1] for row in SAMPLES]),
        biological_replicate=np.asarray([row[2] for row in SAMPLES]),
    )
    return {
        "schema": "aleph.outer_library.gse325309_pristine_manifest.v1",
        "source": receipt, "protocol": protocol,
        "counts": {"biological_samples": 6, "observations": 12, "genes": 2,
                   "control_replicates": 3, "TGFB1_replicates": 3},
        "tensor_shape": [6, 2], "sample_exclusions": 0,
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build(args.source, args.receipt, args.protocol, args.manifest, args.tensor)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
