#!/usr/bin/env python3
"""Acquire the frozen IDR0168 single-use confirmation archive subset."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import concurrent.futures
import csv
import io
import json
from pathlib import Path
import time
import urllib.request

from train_hpa_raw_if_cnn import sha256


USER_AGENT = "Project-Aleph-research/1.0"


def request(url: str, *, method: str = "GET", headers: dict[str, str] | None = None):
    values = {"User-Agent": USER_AGENT}
    if headers:
        values.update(headers)
    return urllib.request.Request(url, method=method, headers=values)


def download_small(url: str, path: Path, expected_size: int, expected_sha: str) -> None:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(request(url), timeout=180) as response:
            data = response.read()
        path.write_bytes(data)
    if path.stat().st_size != expected_size or sha256(path) != expected_sha:
        raise ValueError(f"official metadata hash or size changed: {path.name}")


def head_size(url: str) -> int:
    with urllib.request.urlopen(request(url, method="HEAD"), timeout=180) as response:
        if response.status != 200:
            raise ValueError(f"unexpected HEAD status {response.status}: {url}")
        value = response.headers.get("Content-Length")
    if value is None or int(value) <= 0:
        raise ValueError(f"missing positive Content-Length: {url}")
    return int(value)


def download_archive(url: str, path: Path, expected_size: int) -> dict[str, object]:
    if path.exists() and path.stat().st_size == expected_size:
        return {
            "local_name": path.name,
            "size_bytes": expected_size,
            "sha256": sha256(path),
            "download_reused": True,
        }
    if path.exists():
        raise ValueError(f"completed archive has wrong size: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(6):
        start = partial.stat().st_size if partial.exists() else 0
        if start > expected_size:
            partial.unlink()
            start = 0
        headers = {"Range": f"bytes={start}-"} if start else {}
        try:
            with urllib.request.urlopen(
                request(url, headers=headers), timeout=300
            ) as response:
                append = bool(start and response.status == 206)
                mode = "ab" if append else "wb"
                if start and not append:
                    start = 0
                with partial.open(mode) as handle:
                    while chunk := response.read(8 * 1024 * 1024):
                        handle.write(chunk)
            if partial.stat().st_size != expected_size:
                raise IOError(
                    f"archive size {partial.stat().st_size} != {expected_size}"
                )
            partial.replace(path)
            return {
                "local_name": path.name,
                "size_bytes": expected_size,
                "sha256": sha256(path),
                "download_reused": False,
            }
        except Exception as exc:
            last_error = exc
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"failed archive download after retries: {url}: {last_error}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if not protocol["frozen_before_any_idr0168_pixel_access_or_model_prediction"]:
        raise ValueError("IDR0168 protocol is not frozen before pixel access")
    source = protocol["official_source"]
    metadata_specs = (
        (
            "annotation.csv",
            source["annotation_csv_url"],
            source["annotation_csv_size_bytes"],
            source["annotation_csv_sha256"],
        ),
        (
            "filePaths.tsv",
            source["file_paths_tsv_url"],
            source["file_paths_tsv_size_bytes"],
            source["file_paths_tsv_sha256"],
        ),
        (
            "study.txt",
            source["study_metadata_url"],
            source["study_metadata_size_bytes"],
            source["study_metadata_sha256"],
        ),
    )
    for name, url, size, digest in metadata_specs:
        download_small(url, args.cache_dir / name, int(size), digest)

    annotation_path = args.cache_dir / "annotation.csv"
    with annotation_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        headers = reader.fieldnames or []
        rows = list(reader)
    contract = protocol["metadata_contract"]
    if len(rows) != contract["exact_row_count"]:
        raise ValueError("IDR0168 annotation row count changed")
    missing = sorted(set(contract["required_columns"]) - set(headers))
    if missing:
        raise ValueError(f"IDR0168 required columns missing: {missing}")
    background = [row for row in rows if not row["Comment [Gene Symbol]"].strip()]
    admitted = [row for row in rows if row["Comment [Gene Symbol]"].strip()]
    if len(background) != contract["background_row_count"]:
        raise ValueError("IDR0168 background-row count changed")
    if len(admitted) != contract["nonbackground_row_count"]:
        raise ValueError("IDR0168 nonbackground-row count changed")
    cell_lines = sorted({row["Characteristics [Cell Line]"] for row in admitted})
    genes = sorted({row["Comment [Gene Symbol]"] for row in admitted})
    if cell_lines != sorted(contract["exact_cell_lines"]):
        raise ValueError("IDR0168 cell-line set changed")
    if genes != sorted(contract["exact_gene_symbols"]):
        raise ValueError("IDR0168 gene set changed")
    for row in admitted:
        gene = row["Comment [Gene Symbol]"]
        expected_channels = contract["required_channel_string_template"].format(
            gene_symbol=gene
        )
        if row["Channels"] != expected_channels:
            raise ValueError(f"IDR0168 channel contract changed: {row['Image Name']}")
        if row["Comment [Image File Type]"] != contract["required_image_file_type"]:
            raise ValueError(f"IDR0168 non-raw source row: {row['Image Name']}")

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in admitted:
        grouped[(row["Comment [Gene Symbol]"], row["Characteristics [Cell Line]"])].append(row)
    selected = [min(values, key=lambda row: row["Image Name"]) for values in grouped.values()]
    selected.sort(key=lambda row: (
        row["Comment [Gene Symbol]"], row["Characteristics [Cell Line]"]
    ))
    selection = protocol["single_use_selection"]
    if len(selected) != selection["expected_selected_stack_count"]:
        raise ValueError("IDR0168 selected stack count changed")
    if len(grouped) != selection["expected_group_count"]:
        raise ValueError("IDR0168 group count changed")

    file_paths = set()
    with (args.cache_dir / "filePaths.tsv").open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) != 3:
                raise ValueError("IDR0168 filePaths TSV row width changed")
            file_paths.add(row[2])
    if any(row["Image Name"] not in file_paths for row in selected):
        raise ValueError("selected IDR0168 image missing from pinned filePaths TSV")

    base_url = source["archive_base_url"]
    archive_rows = []
    for row in selected:
        submitted = row["Image File"]
        if not submitted.endswith(".zarr.zip"):
            raise ValueError(f"unexpected submitted archive suffix: {submitted}")
        archive_name = submitted.removesuffix(".zip") + ".tgz"
        archive_rows.append({
            "cell_line": row["Characteristics [Cell Line]"],
            "gene_symbol": row["Comment [Gene Symbol]"],
            "gene_identifier": row["Comment [Gene Identifier]"],
            "uniprot_id": row["Comment [UniProt ID]"],
            "image_name": row["Image Name"],
            "source_name": row["Source Name"],
            "channels": row["Channels"],
            "archive_name": archive_name,
            "url": base_url + archive_name,
        })
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        sizes = list(pool.map(lambda row: head_size(row["url"]), archive_rows))
    for row, size in zip(archive_rows, sizes):
        row["size_bytes"] = size
    total_bytes = sum(sizes)
    if total_bytes != selection["expected_selected_archive_total_bytes"]:
        raise ValueError(
            f"IDR0168 selected archive total changed: {total_bytes}"
        )

    archive_dir = args.cache_dir / "archives"
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(
                download_archive,
                row["url"],
                archive_dir / row["archive_name"],
                int(row["size_bytes"]),
            ): row
            for row in archive_rows
        }
        completed = 0
        for future in concurrent.futures.as_completed(futures):
            row = futures[future]
            row.update(future.result())
            completed += 1
            print(
                json.dumps({
                    "downloaded_or_verified": completed,
                    "total": len(archive_rows),
                    "archive": row["archive_name"],
                }, sort_keys=True),
                flush=True,
            )

    report = {
        "schema": "aleph.outer_library.idr0168_confirmation_acquisition.v1",
        "protocol_sha256": sha256(args.protocol),
        "metadata_hashes_verified": True,
        "metadata_schema_verified": True,
        "channel_strings_verified": True,
        "row_count": len(rows),
        "selected_stack_count": len(archive_rows),
        "selected_group_count": len(grouped),
        "selected_gene_counts": dict(sorted(Counter(
            row["gene_symbol"] for row in archive_rows
        ).items())),
        "selected_cell_line_counts": dict(sorted(Counter(
            row["cell_line"] for row in archive_rows
        ).items())),
        "selected_archive_total_bytes": total_bytes,
        "archives": sorted(archive_rows, key=lambda row: row["archive_name"]),
        "archive_body_objects_accessed": len(archive_rows),
        "pixel_archives_opened_for_decode": 0,
        "model_predictions_made": 0,
        "external_provider_validation_complete": False,
        "Aleph_parameter_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "selected_stacks": len(archive_rows),
        "bytes": total_bytes,
        "report": str(args.report),
        "pixels_decoded": 0,
        "predictions": 0,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
