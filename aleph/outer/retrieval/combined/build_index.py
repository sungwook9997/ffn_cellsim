#!/usr/bin/env python
"""Build the compact text-free two-pass OA retrieval index."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[3]


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
    model_dir = HERE.parents[1] / "model" / "combined_content"
    split_rows = [json.loads(line) for line in (model_dir / "splits.jsonl").read_text().splitlines()]
    model_by_family = {row["source_family_id"]: row for row in split_rows}
    if len(model_by_family) != len(split_rows):
        raise RuntimeError("combined_model_duplicate_source_family")
    rag_path = HERE.parents[1] / "ragtagcag" / "second_pass" / "combined_ingest_2026-08-05.json"
    rag = json.loads(rag_path.read_text())
    rows = []
    section_totals: Counter[str] = Counter()
    pass_counts: Counter[str] = Counter()
    chunks_total = 0
    for pass_name, relative_root in sorted(config["derived_roots"].items()):
        derived_root = REPO_ROOT / relative_root
        for chunk_path in sorted(derived_root.glob("*.chunks.jsonl")):
            manifest_path = chunk_path.with_name(chunk_path.name.replace(".chunks.jsonl", ".manifest.json"))
            manifest = json.loads(manifest_path.read_text())
            chunk_digest = sha256(chunk_path)
            if manifest["chunks_file_sha256"] != chunk_digest:
                raise RuntimeError(f"derived_manifest_digest_mismatch:{pass_name}:{chunk_path.name}")
            if manifest["authority_status"] != "proposed":
                raise RuntimeError(f"non_proposed_derived_article:{pass_name}:{chunk_path.name}")
            section_counts: Counter[str] = Counter()
            with chunk_path.open(encoding="utf-8") as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    chunk = json.loads(line)
                    if chunk["source_family_id"].lower() != manifest["source_family_id"].lower():
                        raise RuntimeError(f"chunk_source_family_mismatch:{pass_name}:{chunk_path.name}")
                    if chunk["authority_status"] != "proposed":
                        raise RuntimeError(f"non_proposed_chunk:{pass_name}:{chunk_path.name}")
                    section_counts[chunk["section_class"]] += 1
            count = sum(section_counts.values())
            if count != manifest["chunks"]:
                raise RuntimeError(f"chunk_count_mismatch:{pass_name}:{chunk_path.name}")
            model_row = model_by_family.get(manifest["source_family_id"].lower())
            if model_row is not None and model_row["pass"] != pass_name:
                raise RuntimeError(f"model_pass_mismatch:{manifest['source_family_id']}")
            chunks_total += count
            pass_counts[pass_name] += 1
            section_totals.update(section_counts)
            rows.append({
                "authority_status": "proposed",
                "chunks": count,
                "chunks_file": chunk_path.name,
                "chunks_file_sha256": chunk_digest,
                "neural_routable": model_row is not None,
                "pass": pass_name,
                "pmcid": manifest["pmcid"],
                "proposed_coarse_domain": model_row["label"] if model_row else None,
                "section_class_counts": dict(sorted(section_counts.items())),
                "source_family_id": manifest["source_family_id"].lower(),
                "split": model_row["split"] if model_row else None,
                "weak_label": model_row["weak_label"] if model_row else None,
            })
    rows.sort(key=lambda row: (row["source_family_id"], row["pass"]))
    families = [row["source_family_id"] for row in rows]
    if len(families) != len(set(families)):
        raise RuntimeError("cross_pass_local_source_family_overlap")
    index_bytes = b"".join(
        (json.dumps(row, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode()
        for row in rows
    )
    atomic_write(output_dir / "index.jsonl", index_bytes)
    model_sums = json.loads((model_dir / "SHA256SUMS.json").read_text())
    receipt = {
        "authority_status": "proposed",
        "chunk_count": chunks_total,
        "index_sha256": hashlib.sha256(index_bytes).hexdigest(),
        "local_source_families": len(rows),
        "model": {
            "coarse_taxonomy_sha256": model_sums["taxonomy.json"],
            "neural_routable_source_families": len(model_by_family),
            "splits_sha256": model_sums["splits.jsonl"],
            "weights_sha256": model_sums["weights_combined_content.npz"],
        },
        "pass_source_families": dict(sorted(pass_counts.items())),
        "ragtagcag": {
            "combined_source_families": rag["families"]["combined_unique"],
            "ledger_head_hash": rag["first"]["ledger_head_hash"],
            "projection_hash": rag["first"]["projection_hash"],
            "snapshot_manifest_hash": rag["first"]["snapshot_manifest_hash"],
            "validation_receipt_sha256": sha256(rag_path),
        },
        "schema": "aleph.external_training.combined_retrieval_index_receipt.v1",
        "section_class_counts": dict(sorted(section_totals.items())),
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
