"""Deterministic combined-corpus validation over Aleph's public harness."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from aleph.harness import DataRoot, ObjectStore
from aleph.harness.snapshot import resolve
from aleph.outer.ragtagcag.ingest import verify_idempotent_ingest


class CombinedValidationError(RuntimeError):
    """The combined projection did not preserve its evidence contract."""


def _records(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, dict):
                raise CombinedValidationError(f"{path}:{line_number}: expected an object")
            rows.append(value)
    return rows


def _family(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        raise CombinedValidationError("record lacks source_family_id")
    prefix, separator, identifier = text.partition(":")
    if separator and prefix.lower() in {"doi", "pmid", "pmcid", "arxiv"}:
        return f"{prefix.lower()}:{identifier.strip().lower()}"
    return text


def _families(paths: Sequence[Path]) -> set[str]:
    return {_family(row.get("source_family_id")) for path in paths for row in _records(path)}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _receipt_summary(paths: Sequence[Path]) -> dict[str, Any]:
    rows = [row for path in paths for row in _records(path)]
    statuses = Counter(str(row.get("status", "unknown")) for row in rows)
    return {
        "files": len(paths),
        "records": len(rows),
        "status_counts": dict(sorted(statuses.items())),
        "successful_bytes": sum(
            int(row.get("bytes", 0)) for row in rows if row.get("status") == "success"
        ),
        "file_sha256": {_portable(path): _sha256(path) for path in sorted(paths)},
    }


def _portable(path: Path) -> str:
    marker = "aleph/outer/"
    text = str(path)
    index = text.find(marker)
    return text[index:] if index >= 0 else path.name


def _audit_receipt_objects(data_root: DataRoot, snapshot_id: str) -> dict[str, int]:
    snapshot = resolve(data_root, snapshot_id)
    store = ObjectStore.open(data_root)
    status_counts: Counter[str] = Counter()
    receipt_objects = 0
    payload_present = 0
    for entity in snapshot.entities.values():
        if entity.attrs.get("authority_status") != "proposed" or entity.attrs.get("status") != "proposed":
            raise CombinedValidationError(f"authority promotion in {entity.subject_id}")
        for summary in entity.attrs.get("acquisition_receipts", []):
            digest = summary["receipt_object"]
            receipt = store.get_json(digest, verify=True)
            receipt_objects += 1
            status_counts[str(receipt.get("status", "unknown"))] += 1
            payload_digest = summary.get("payload_sha256")
            if summary.get("payload_present"):
                if not payload_digest or not store.exists(payload_digest) or not store.verify(payload_digest):
                    raise CombinedValidationError(f"missing or corrupt payload for {entity.subject_id}")
                payload_present += 1
    return {
        "receipt_objects": receipt_objects,
        "receipt_success": status_counts["success"],
        "receipt_failed": status_counts["failed"],
        "receipt_unknown": sum(status_counts.values())
        - status_counts["success"]
        - status_counts["failed"],
        "payloads_present_and_verified": payload_present,
    }


def run_combined_validation(
    *,
    initial_manifests: Iterable[str | Path],
    second_pass_manifest: str | Path,
    receipt_paths: Iterable[str | Path],
    data_root: DataRoot,
) -> dict[str, Any]:
    """Run twice and return a compact, text-free validation receipt."""

    initial = tuple(Path(path).resolve() for path in initial_manifests)
    second = Path(second_pass_manifest).resolve()
    receipts = tuple(Path(path).resolve() for path in receipt_paths)
    initial_families = _families(initial)
    second_families = _families((second,))
    combined_families = initial_families | second_families

    result = verify_idempotent_ingest(
        (*initial, second), data_root, receipt_paths=receipts
    )
    first = result["first"]
    second_run = result["second"]
    if first["source_families"] != len(combined_families):
        raise CombinedValidationError("adapter family count disagrees with independent preflight")
    if second_run["events_appended"] != 0:
        raise CombinedValidationError("second run appended events")

    receipt_audit = _audit_receipt_objects(data_root, first["snapshot_id"])
    expected_receipts = _receipt_summary(receipts)
    if receipt_audit["receipt_objects"] != expected_receipts["records"]:
        raise CombinedValidationError("TAG receipt-object count disagrees with receipt inputs")
    if first["missing_objects"] or first["cag_lint_violations"]:
        raise CombinedValidationError("missing object or CAG lint violation")

    return {
        "schema": "aleph.external_training.ragtagcag_combined.v1",
        "authority_status": "proposed",
        "families": {
            "initial_unique": len(initial_families),
            "second_pass_candidate_unique": len(second_families),
            "overlap": len(initial_families & second_families),
            "new_unique": len(second_families - initial_families),
            "combined_unique": len(combined_families),
        },
        "receipts": expected_receipts,
        "object_audit": receipt_audit,
        "first": first,
        "second": second_run,
        "idempotent": True,
        "committed_payloads": False,
    }
