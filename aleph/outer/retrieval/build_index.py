#!/usr/bin/env python
"""Build a compact, text-free receipt index for locally derived OA chunks."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        handle.write(payload)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def build(output_dir: Path = HERE) -> dict:
    config = json.loads((HERE / "config.json").read_text())
    derived_root = REPO_ROOT / config["derived_root"]
    model_dir = HERE.parent / "model" / "oa_content"
    model_rows = [json.loads(line) for line in (model_dir / "splits.jsonl").read_text().splitlines()]
    rag_result_path = HERE.parent / "ragtagcag" / "results" / "full_ingest_2026-08-05.json"
    rag_result = json.loads(rag_result_path.read_text())
    rows = []
    classes: Counter[str] = Counter()
    total_chunks = 0
    for model_row in sorted(model_rows, key=lambda row: row["source_family_id"]):
        chunk_path = derived_root / model_row["derived_file"]
        if not chunk_path.is_file():
            continue
        manifest_path = chunk_path.with_name(chunk_path.name.replace(".chunks.jsonl", ".manifest.json"))
        manifest = json.loads(manifest_path.read_text())
        chunk_digest = sha256(chunk_path)
        if manifest["chunks_file_sha256"] != chunk_digest:
            raise RuntimeError(f"derived_manifest_digest_mismatch:{chunk_path.name}")
        if manifest["authority_status"] != "proposed":
            raise RuntimeError(f"non_proposed_derived_article:{chunk_path.name}")
        section_counts: Counter[str] = Counter()
        with chunk_path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                chunk = json.loads(line)
                if chunk["source_family_id"].lower() != model_row["source_family_id"].lower():
                    raise RuntimeError(f"chunk_source_family_mismatch:{chunk_path.name}")
                if chunk["authority_status"] != "proposed":
                    raise RuntimeError(f"non_proposed_chunk:{chunk_path.name}")
                section_counts[chunk["section_class"]] += 1
        count = sum(section_counts.values())
        if count != manifest["chunks"]:
            raise RuntimeError(f"chunk_count_mismatch:{chunk_path.name}")
        total_chunks += count
        classes.update(section_counts)
        rows.append({
            "authority_status": "proposed",
            "chunks": count,
            "chunks_file": chunk_path.name,
            "chunks_file_sha256": chunk_digest,
            "pmcid": manifest["pmcid"],
            "section_class_counts": dict(sorted(section_counts.items())),
            "source_family_id": model_row["source_family_id"],
            "split": model_row["split"],
            "weak_domain": model_row["label"],
        })
    index_bytes = b"".join(
        (json.dumps(row, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode()
        for row in rows
    )
    index_path = output_dir / "index.jsonl"
    atomic_write(index_path, index_bytes)
    receipt = {
        "authority_status": "proposed",
        "chunk_count": total_chunks,
        "index_sha256": hashlib.sha256(index_bytes).hexdigest(),
        "local_source_families": len(rows),
        "model_content_weights_sha256": json.loads((model_dir / "SHA256SUMS.json").read_text())["weights_oa_content.npz"],
        "model_source_families": len(model_rows),
        "ragtagcag": {
            "ledger_head_hash": rag_result["first_ingest"]["ledger_head_hash"],
            "projection_hash": rag_result["first_ingest"]["projection_hash"],
            "source_families": rag_result["input"]["source_families"],
            "validation_receipt_sha256": sha256(rag_result_path),
        },
        "schema": "aleph.external_training.retrieval_index_receipt.v1",
        "section_class_counts": dict(sorted(classes.items())),
        "text_committed": False,
    }
    atomic_write(output_dir / "receipt.json", (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode())
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    print(json.dumps(build(args.output_dir.resolve()), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
