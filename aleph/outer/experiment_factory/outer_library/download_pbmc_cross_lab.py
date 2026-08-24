#!/usr/bin/env python3
"""Acquire pinned public PBMC sources and only the required benchmark labels."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import struct
import tarfile
import urllib.request
import zlib
from pathlib import Path
from typing import Any


GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
FILES = {
    "GSE96583_RAW.tar": f"{GEO_BASE}/GSE96nnn/GSE96583/suppl/GSE96583_RAW.tar",
    "GSE96583_batch1.total.tsne.df.tsv.gz": f"{GEO_BASE}/GSE96nnn/GSE96583/suppl/GSE96583_batch1.total.tsne.df.tsv.gz",
    "GSE96583_batch2.total.tsne.df.tsv.gz": f"{GEO_BASE}/GSE96nnn/GSE96583/suppl/GSE96583_batch2.total.tsne.df.tsv.gz",
    "GSE96583_batch2.genes.tsv.gz": f"{GEO_BASE}/GSE96nnn/GSE96583/suppl/GSE96583_batch2.genes.tsv.gz",
    "GSE132044_pbmc_hg38_cell.tsv.gz": f"{GEO_BASE}/GSE132nnn/GSE132044/suppl/GSE132044_pbmc_hg38_cell.tsv.gz",
    "GSE132044_pbmc_hg38_gene.tsv.gz": f"{GEO_BASE}/GSE132nnn/GSE132044/suppl/GSE132044_pbmc_hg38_gene.tsv.gz",
    "GSE132044_pbmc_hg38_count_matrix.mtx.gz": f"{GEO_BASE}/GSE132nnn/GSE132044/suppl/GSE132044_pbmc_hg38_count_matrix.mtx.gz",
    "GSE72502_PBMC_RPKMs.xls.gz": f"{GEO_BASE}/GSE72nnn/GSE72502/suppl/GSE72502_PBMC_RPKMs.xls.gz",
    "GSE72502_series_matrix.txt.gz": f"{GEO_BASE}/GSE72nnn/GSE72502/matrix/GSE72502_series_matrix.txt.gz",
    # Original author cell annotations exposed by the public Broad Single Cell
    # Portal study API.  The secondary Zenodo benchmark labels below are kept
    # only for a compatibility audit; their row counts do not match the GEO
    # expression matrix and they must not be used as classifier targets.
    "SCP424_Harmony_TSNE_CellType.json": (
        "https://singlecell.broadinstitute.org/single_cell/api/v1/studies/"
        "SCP424/clusters/Harmony%20TSNE?annotation_name=CellType&"
        "annotation_type=group&annotation_scope=study&subsample=all&consensus="
    ),
}
ZENODO_RECORD = "https://zenodo.org/api/records/3357167"
ZENODO_ARCHIVE = {
    "record_id": 3357167,
    "title": "A comparison of automatic cell identification methods for single-cell RNA-sequencing data",
    "size": 3671466589,
    "checksum": "md5:b799a660b8bcaf5f3580a9b6f9372e5b",
}
ZENODO_LABEL_MEMBERS = (
    "Inter-dataset/PbmcBench/10Xv2/10Xv2_pbmc1Labels.csv",
    "Inter-dataset/PbmcBench/10Xv3/10Xv3_pbmc1Labels.csv",
    "Inter-dataset/PbmcBench/CEL-Seq/CL_pbmc1Labels.csv",
    "Inter-dataset/PbmcBench/Drop-Seq/DR_pbmc1Labels.csv",
    "Inter-dataset/PbmcBench/inDrop/iD_pbmc1Labels.csv",
    "Inter-dataset/PbmcBench/Smart-Seq2/SM2_pbmc1Labels.csv",
    "Inter-dataset/PbmcBench/Seq-Well/SW_pbmc1Labels.csv",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size:
        return
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Mozilla/5.0 (compatible; Project-Aleph/1)",
            "Accept": "application/json,*/*;q=0.8",
        },
    )
    with urllib.request.urlopen(request, timeout=120) as response, path.open("wb") as out:
        while block := response.read(1024 * 1024):
            out.write(block)


def _range(url: str, start: int, end: int) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"Range": f"bytes={start}-{end}", "User-Agent": "Project-Aleph/1"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        data = response.read()
    expected = end - start + 1
    if len(data) != expected:
        raise ValueError(f"range response length {len(data)} != {expected}")
    return data


def _zip_index(url: str, size: int) -> dict[str, dict[str, int]]:
    tail_start = max(0, size - 65557)
    tail = _range(url, tail_start, size - 1)
    eocd_position = tail.rfind(b"PK\x05\x06")
    if eocd_position < 0:
        raise ValueError("ZIP end-of-central-directory not found")
    eocd = struct.unpack_from("<4s4H2LH", tail, eocd_position)
    central_size, central_offset = int(eocd[5]), int(eocd[6])
    central = _range(url, central_offset, central_offset + central_size - 1)
    entries: dict[str, dict[str, int]] = {}
    position = 0
    while position < len(central):
        values = struct.unpack_from("<4s6H3L5H2L", central, position)
        if values[0] != b"PK\x01\x02":
            raise ValueError(f"invalid central-directory signature at {position}")
        method = int(values[4])
        compressed_size, uncompressed_size = int(values[8]), int(values[9])
        name_length, extra_length, comment_length = map(int, values[10:13])
        local_offset = int(values[-1])
        start = position + 46
        name = central[start:start + name_length].decode("utf-8")
        entries[name] = {
            "method": method,
            "compressed_size": compressed_size,
            "uncompressed_size": uncompressed_size,
            "local_offset": local_offset,
        }
        position = start + name_length + extra_length + comment_length
    return entries


def _extract_remote_member(url: str, entry: dict[str, int], output: Path) -> None:
    offset = entry["local_offset"]
    header = _range(url, offset, offset + 29)
    values = struct.unpack("<4s5H3L2H", header)
    if values[0] != b"PK\x03\x04":
        raise ValueError("invalid local ZIP header")
    name_length, extra_length = int(values[-2]), int(values[-1])
    data_start = offset + 30 + name_length + extra_length
    compressed = _range(
        url, data_start, data_start + entry["compressed_size"] - 1
    )
    if entry["method"] == 0:
        payload = compressed
    elif entry["method"] == 8:
        payload = zlib.decompress(compressed, -15)
    else:
        raise ValueError(f"unsupported ZIP compression method {entry['method']}")
    if len(payload) != entry["uncompressed_size"]:
        raise ValueError("remote ZIP member size mismatch")
    output.write_bytes(payload)


def acquire(output_dir: Path, protocol: Path) -> dict[str, Any]:
    protocol_payload = json.loads(protocol.read_text(encoding="utf-8"))
    if not protocol_payload.get("frozen_before_expression_or_cell_annotation_access"):
        raise ValueError("protocol was not frozen before access")
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for name, url in FILES.items():
        path = output_dir / name
        _download(url, path)
        if name.endswith(".gz"):
            with gzip.open(path, "rb") as handle:
                handle.read(1)
        if name.endswith(".tar"):
            with tarfile.open(path) as archive:
                tar_members = [member.name for member in archive.getmembers()]
        else:
            tar_members = None
        records.append({
            "name": name, "url": url, "size_bytes": path.stat().st_size,
            "sha256": _sha256(path), "container_members": tar_members,
        })

    label_records = []
    label_dir = output_dir / "zenodo_3357167_labels"
    label_dir.mkdir(exist_ok=True)
    missing_members = [
        member for member in ZENODO_LABEL_MEMBERS
        if not (label_dir / Path(member).name).exists()
        or not (label_dir / Path(member).name).stat().st_size
    ]
    archive_url = ""
    index: dict[str, dict[str, int]] = {}
    if missing_members:
        with urllib.request.urlopen(ZENODO_RECORD, timeout=120) as response:
            zenodo = json.load(response)
        archive = next(
            item for item in zenodo["files"]
            if item["key"] == "scRNAseq_Benchmark_datasets.zip"
        )
        archive_url = archive["links"]["self"]
        if int(archive["size"]) != ZENODO_ARCHIVE["size"]:
            raise ValueError("pinned Zenodo archive size changed")
        if archive["checksum"] != ZENODO_ARCHIVE["checksum"]:
            raise ValueError("pinned Zenodo archive checksum changed")
        index = _zip_index(archive_url, int(archive["size"]))
    for member in ZENODO_LABEL_MEMBERS:
        path = label_dir / Path(member).name
        if member in missing_members:
            if member not in index:
                raise ValueError(f"required Zenodo label member missing: {member}")
            _extract_remote_member(archive_url, index[member], path)
        label_records.append({
            "member": member, "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        })
    return {
        "schema": "aleph.outer_library.pbmc_cross_lab_acquisition.v1",
        "protocol_sha256": _sha256(protocol),
        "protocol_frozen_commit": "b49f9ed",
        "official_files": records,
        "secondary_benchmark_label_source": {
            "record": ZENODO_RECORD,
            "record_id": ZENODO_ARCHIVE["record_id"],
            "title": ZENODO_ARCHIVE["title"],
            "archive_size_bytes": ZENODO_ARCHIVE["size"],
            "archive_checksum": ZENODO_ARCHIVE["checksum"],
            "selected_members": label_records,
            "full_archive_downloaded": False,
            "classifier_target_authority": "none_pending_compatibility_audit",
        },
        "aleph_authority": "none",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    result = acquire(args.output_dir, args.protocol)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "official_files": len(result["official_files"]),
        "selected_label_members": len(result["secondary_benchmark_label_source"]["selected_members"]),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
