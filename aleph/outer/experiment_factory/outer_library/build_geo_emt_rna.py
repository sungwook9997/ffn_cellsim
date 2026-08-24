#!/usr/bin/env python3
"""Build a study-separated CDH1/VIM EMT RNA tensor and canonical manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile

import numpy as np

from schema import Observation


PROBES = {"CDH1": "201131_s_at", "VIM": "201426_s_at"}
GSE49644_PROBES = {"CDH1": "201131_s_at", "VIM": "1555938_x_at"}
STUDIES = {
    "GSE17708": {"dataset": "geo-gse17708-a549-tgfb-timecourse", "lab": "Kuick_Lab_University_of_Michigan"},
    "GSE42373": {"dataset": "geo-gse42373-a549-2d-3d-tgfb-tnfa", "lab": "Cieslik_Lab_University_of_Virginia"},
    "GSE125369": {"dataset": "geo-gse125369-a549-long-tgfb3-rnaseq", "lab": "Haley_Lab_Stony_Brook_University"},
    "GSE69667": {"dataset": "geo-gse69667-a549-tgfb-timecourse-rnaseq", "lab": "Xue_Lab_Shanghai_Advanced_Research_Institute"},
    "GSE49644": {"dataset": "geo-gse49644-nsclc-long-tgfb-emt", "lab": "Sun_Setttleman_Lab_Genentech"},
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix(
    path: Path, probes: dict[str, str] = PROBES
) -> tuple[list[str], list[str], dict[str, list[float]]]:
    titles: list[str] = []
    accessions: list[str] = []
    values: dict[str, list[float]] = {}
    wanted = set(probes.values())
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            elif line.startswith("!Sample_geo_accession"):
                accessions = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            elif line.startswith('"'):
                row = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                if row[0] in wanted:
                    values[row[0]] = [float(item) for item in row[1:]]
    if not titles or len(titles) != len(accessions) or set(values) != wanted:
        raise ValueError(f"incomplete GEO matrix contract: {path.name}")
    if any(len(row) != len(titles) for row in values.values()):
        raise ValueError(f"GEO matrix sample axis changed: {path.name}")
    return titles, accessions, values


def _gse17708(title: str) -> tuple[str, str, float]:
    replicate = re.search(r"experiment #(\d+)$", title)
    if replicate is None:
        raise ValueError(f"GSE17708 title changed: {title}")
    time = re.search(r"for ([0-9.]+) h", title)
    hours = float(time.group(1)) if time else 0.0
    return ("untreated" if hours == 0 else f"TGFB1_{hours:g}h", replicate.group(1), hours)


def _gse42373(title: str) -> tuple[str, str, float]:
    match = re.fullmatch(r"a549-(2d|3d)-(control|treated)-(\d+)", title)
    if match is None:
        raise ValueError(f"GSE42373 title changed: {title}")
    culture, treatment, replicate = match.groups()
    return f"{culture}_{treatment}", replicate, 96.0 if treatment == "treated" else 0.0


def _gse69667_samples(source_dir: Path) -> list[dict[str, object]]:
    times = (0, 6, 12, 24, 36, 48, 72, 96)
    result: list[dict[str, object]] = []
    for replicate in (1, 2):
        path = source_dir / f"GSE69667_RNASEQ_TPM_rep{replicate}.txt.gz"
        source_hash = _sha256(path)
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["GENE", "ENSG", *(str(hour) for hour in times)]:
                raise ValueError(f"GSE69667 TPM header changed: rep{replicate}")
            genes = {
                row["GENE"]: row for row in reader if row["GENE"] in PROBES
            }
        if set(genes) != set(PROBES):
            raise ValueError(f"GSE69667 marker inventory changed: rep{replicate}")
        for time_index, hours in enumerate(times):
            sample_number = 1_706_458 + 2 * time_index + (replicate - 1)
            sample_id = f"GSM{sample_number}"
            condition = "untreated" if hours == 0 else f"TGFB1_{hours}h"
            for gene in PROBES:
                result.append({
                    "study": "GSE69667", "sample_id": sample_id,
                    "condition": condition, "replicate": str(replicate),
                    "hours": float(hours), "gene": gene,
                    "value": float(genes[gene][str(hours)]), "unit": "TPM",
                    "source_hash": source_hash,
                    "locator": f"{path.name}:{gene}:{hours}h",
                    "cell_type": "A549_human_lung_adenocarcinoma_cell_line",
                })
    return result


def _validate_gse49644_annotation(path: Path) -> None:
    found: dict[str, tuple[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["symbol"] in GSE49644_PROBES:
                gene = row["symbol"]
                if gene in found:
                    raise ValueError(f"GSE49644 exact-symbol marker is not unique: {gene}")
                found[gene] = (row["ProbeID"], row["entrez.gene"])
    expected = {
        "CDH1": (GSE49644_PROBES["CDH1"], "999"),
        "VIM": (GSE49644_PROBES["VIM"], "7431"),
    }
    if found != expected:
        raise ValueError(f"GSE49644 exact-symbol annotation changed: {found}")


def _gse49644_samples(source_dir: Path) -> list[dict[str, object]]:
    path = source_dir / "GSE49644_series_matrix.txt.gz"
    _validate_gse49644_annotation(
        source_dir / "GSE49644_annotated_results.txt.gz"
    )
    titles, sample_ids, matrix = _matrix(path, GSE49644_PROBES)
    source_hash = _sha256(path)
    group_counts: dict[tuple[str, str], int] = {}
    result: list[dict[str, object]] = []
    for index, (title, sample_id) in enumerate(zip(titles, sample_ids)):
        match = re.fullmatch(r"\d+\.(A549|HCC827|NCI-H358)\.(Parental|EMT)", title)
        if match is None:
            raise ValueError(f"GSE49644 title changed: {title}")
        cell_line, state = match.groups()
        key = (cell_line, state)
        group_counts[key] = group_counts.get(key, 0) + 1
        condition = f"{cell_line}_{'parental' if state == 'Parental' else 'TGFB_EMT_3week'}"
        for gene, probe in GSE49644_PROBES.items():
            result.append({
                "study": "GSE49644", "sample_id": sample_id,
                "condition": condition, "replicate": str(group_counts[key]),
                "hours": 0.0 if state == "Parental" else 504.0,
                "gene": gene, "value": matrix[probe][index],
                "unit": "GEO_series_matrix_log2_expression",
                "source_hash": source_hash,
                "locator": f"{path.name}:{probe}:{sample_id}",
                "cell_type": f"{cell_line}_human_NSCLC_cell_line",
            })
    if set(group_counts.values()) != {3} or len(group_counts) != 6:
        raise ValueError(f"GSE49644 replicate hierarchy changed: {group_counts}")
    return result


def _rsem_samples(path: Path) -> list[tuple[str, str, str, dict[str, float]]]:
    members = {
        "GSM3572736_A549_EMT_L4.LB9.rsem.genes.results.txt.gz": ("GSM3572736", "TGFB3_21d", "LB9"),
        "GSM3572737_A549_EMT_L4.LB10.rsem.genes.results.txt.gz": ("GSM3572737", "TGFB3_21d", "LB10"),
        "GSM3572738_CNTR_L4.LB13.rsem.genes.results.txt.gz": ("GSM3572738", "control", "LB13"),
        "GSM3572739_CNTR_L4.LB16.rsem.genes.results.txt.gz": ("GSM3572739", "control", "LB16"),
    }
    result = []
    with tarfile.open(path) as archive:
        if set(archive.getnames()) != set(members):
            raise ValueError("GSE125369 RAW member inventory changed")
        for name, metadata in members.items():
            extracted = archive.extractfile(name)
            if extracted is None:
                raise ValueError(f"missing archive member: {name}")
            genes = {}
            with gzip.open(io.BytesIO(extracted.read()), "rt", encoding="utf-8") as handle:
                for row in csv.reader(handle, delimiter="\t"):
                    if row[0] in PROBES:
                        genes[row[0]] = float(row[2]) * 1_000_000.0
            if set(genes) != set(PROBES):
                raise ValueError(f"GSE125369 marker rows changed: {name}")
            result.append((*metadata, genes))
    return result


def _validate_probe_annotation(path: Path) -> None:
    observed: dict[str, set[str]] = {gene: set() for gene in PROBES}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) >= 3 and row[2] in observed:
                observed[row[2]].add(row[0])
    for gene, probe in PROBES.items():
        if probe not in observed[gene]:
            raise ValueError(f"GPL570 no longer maps {probe} to {gene}")


def build(source_dir: Path, receipt_path: Path, manifest_path: Path, tensor_path: Path) -> dict[str, object]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("schema") != "aleph.outer_library.geo_emt_rna_acquisition.v1":
        raise ValueError("GEO EMT receipt schema changed")
    for name, contract in receipt["files"].items():
        path = source_dir / name
        if not path.is_file() or path.stat().st_size != contract["size"] or _sha256(path) != contract["sha256"]:
            raise ValueError(f"GEO EMT source changed: {name}")
    _validate_probe_annotation(source_dir / "GPL570.annot.gz")

    rows: list[dict[str, object]] = []
    for accession, parser in (("GSE17708", _gse17708), ("GSE42373", _gse42373)):
        path = source_dir / f"{accession}_series_matrix.txt.gz"
        titles, sample_ids, matrix = _matrix(path)
        source_hash = _sha256(path)
        for index, (title, sample_id) in enumerate(zip(titles, sample_ids)):
            condition, replicate, hours = parser(title)
            for gene, probe in PROBES.items():
                rows.append({"study": accession, "sample_id": sample_id, "condition": condition,
                             "replicate": replicate, "hours": hours, "gene": gene,
                             "value": matrix[probe][index], "unit": "GEO_series_matrix_log2_expression",
                             "source_hash": source_hash, "locator": f"{path.name}:{probe}:{sample_id}",
                             "cell_type": "A549_human_lung_adenocarcinoma_cell_line"})
    raw = source_dir / "GSE125369_RAW.tar"
    source_hash = _sha256(raw)
    for sample_id, condition, replicate, genes in _rsem_samples(raw):
        for gene, value in genes.items():
            rows.append({"study": "GSE125369", "sample_id": sample_id, "condition": condition,
                         "replicate": replicate, "hours": 504.0 if condition != "control" else 0.0,
                         "gene": gene, "value": value, "unit": "RSEM_normalized_fraction_times_1e6",
                         "source_hash": source_hash, "locator": f"GSE125369_RAW.tar:{sample_id}:{gene}",
                         "cell_type": "A549_human_lung_adenocarcinoma_cell_line"})
    rows.extend(_gse69667_samples(source_dir))
    rows.extend(_gse49644_samples(source_dir))
    if len(rows) != 144 or len({(r["study"], r["sample_id"]) for r in rows}) != 72:
        raise ValueError("GEO EMT sample/observation count changed")

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8") as handle:
        for row in rows:
            study = str(row["study"])
            meta = STUDIES[study]
            observation = Observation(
                observation_id=f"{study}:{row['sample_id']}:{row['gene']}",
                dataset_id=meta["dataset"], lab_group=meta["lab"],
                provider="NCBI Gene Expression Omnibus", accession=study,
                license_id="geo_public_no_explicit_license", source_sha256=str(row["source_hash"]),
                source_locator=str(row["locator"]), sample_id=f"{study}:{row['sample_id']}",
                biological_replicate=f"{study}:{row['condition']}:biological_replicate_{row['replicate']}",
                technical_replicate="none_reported", condition=str(row["condition"]),
                cell_type=str(row["cell_type"]),
                cell_state="author_conditioned_EMT_context", modality="bulk_RNA_expression_derived",
                observable=f"{row['gene']}_RNA_expression", value=float(row["value"]),
                unit=str(row["unit"]), coordinate_name="TGFB_exposure_time",
                coordinate_value=float(row["hours"]), coordinate_unit="hour",
            )
            observation.validate()
            handle.write(json.dumps(observation.as_dict(), sort_keys=True) + "\n")

    np.savez_compressed(
        tensor_path,
        study=np.asarray([r["study"] for r in rows], dtype="U12"),
        sample_id=np.asarray([r["sample_id"] for r in rows], dtype="U16"),
        condition=np.asarray([r["condition"] for r in rows], dtype="U32"),
        biological_replicate=np.asarray([r["replicate"] for r in rows], dtype="U8"),
        exposure_hours=np.asarray([r["hours"] for r in rows], dtype=np.float64),
        gene=np.asarray([r["gene"] for r in rows], dtype="U4"),
        value=np.asarray([r["value"] for r in rows], dtype=np.float64),
        unit=np.asarray([r["unit"] for r in rows], dtype="U40"),
        cell_type=np.asarray([r["cell_type"] for r in rows], dtype="U48"),
    )
    counts = {study: len({r["sample_id"] for r in rows if r["study"] == study}) for study in STUDIES}
    return {
        "schema": "aleph.outer_library.geo_emt_rna_manifest.v1",
        "counts": {"studies": 5, "independent_lab_groups": 5, "samples": 72,
                   "observations": 144, "samples_by_study": counts},
        "features": ["CDH1_RNA_expression", "VIM_RNA_expression"],
        "probe_contract": {
            "GPL570_default": PROBES,
            "GSE49644_author_annotated_exact_symbol_unique_rows": GSE49644_PROBES,
            "selection": (
                "deterministic canonical _s_at probes except GSE49644, whose author-"
                "annotated release omits 201426_s_at and contains one exact-symbol/"
                "Entrez VIM row; no replicate or alternate probe was selected by outcome"
            ),
        },
        "split_contract": {"unit": "biological_sample_nested_within_study",
                           "study_is_outer_holdout_unit": True,
                           "RNA_features_are_not_independent_replicates": True},
        "scale_contract": {"within_study_effects_only": True,
                           "cross_platform_absolute_values_may_not_be_combined": True,
                           "RNA_and_protein_absolute_values_may_not_be_equated": True},
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tensor", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build(args.source_dir, args.receipt, args.manifest, args.tensor)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
