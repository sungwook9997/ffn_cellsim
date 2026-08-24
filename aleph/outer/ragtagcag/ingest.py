"""Ingest external-source manifests through Aleph's public harness APIs.

The adapter deliberately contains no second database or authority model.  A
source import is a RAG event, the current source projection is a TAG snapshot,
and the prose attached to the import is checked by the shipped CAG lint.  Every
record remains ``proposed``; importing literature cannot promote an Aleph claim.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from aleph.harness import CAG, DataRoot, Event, Ledger, ObjectStore, create_snapshot
from aleph.harness.ledger import canonical_json
from aleph.harness.paths import repo_root
from aleph.harness.snapshot import Snapshot, latest

ACTOR = "external-corpus-ragtagcag-adapter"
AUTHORITY_EPOCH = "R0"
EVENT_IMPORTED = "source.metadata.imported"
EVENT_REVISED = "source.metadata.revised"
EVENT_RETRACTED = "source.metadata.retracted"


class CorpusIngestError(RuntimeError):
    """The input cannot be represented without losing provenance or authority."""


@dataclass(frozen=True)
class IngestResult:
    input_files: int
    input_records: int
    source_families: int
    duplicate_records: int
    events_appended: int
    ledger_events: int
    objects: int
    acquisition_receipts: int
    acquisition_payloads_present: int
    acquisition_payloads_unavailable: int
    retracted_entities: int
    snapshot_id: str
    snapshot_manifest_hash: str
    projection_hash: str
    ledger_head_hash: str
    chain_verified: bool
    missing_objects: int
    cag_lint_violations: int

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def _read_json_records(path: Path) -> list[dict[str, Any]]:
    if path.suffix.lower() == ".jsonl":
        rows: list[dict[str, Any]] = []
        with path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                value = json.loads(line)
                if not isinstance(value, dict):
                    raise CorpusIngestError(f"{path}:{line_number}: record is not an object")
                rows.append(value)
        return rows
    value = json.loads(path.read_text("utf-8"))
    if isinstance(value, dict):
        return [value]
    if isinstance(value, list) and all(isinstance(item, dict) for item in value):
        return list(value)
    raise CorpusIngestError(f"{path}: expected a JSON object, object list, or JSONL")


def _normalise_family(value: Any) -> str:
    family = str(value or "").strip()
    if not family:
        raise CorpusIngestError("source record has no source_family_id")
    prefix, separator, identifier = family.partition(":")
    if separator and prefix.lower() in {"doi", "pmid", "pmcid", "arxiv"}:
        return f"{prefix.lower()}:{identifier.strip().lower()}"
    return family


def _subject_id(family: str) -> str:
    token = hashlib.sha256(family.encode("utf-8")).hexdigest()[:24]
    return f"source/external-{token}"


def _record_rank(path: Path, record: Mapping[str, Any]) -> tuple[int, int, str]:
    # Hand-curated seed records carry explicit access/licence evidence and may
    # fill or replace provider-only values.  No rank changes authority.
    curated = 1 if "sources" in path.parts else 0
    populated = sum(value not in (None, "", [], {}) for value in record.values())
    return curated, populated, canonical_json(record)


def _portable_path(path: Path) -> str:
    try:
        return str(path.relative_to(repo_root()))
    except ValueError:
        return str(path)


def _canonical_source(
    family: str, variants: Sequence[tuple[Path, Mapping[str, Any]]]
) -> dict[str, Any]:
    for path, record in variants:
        authority = record.get("authority_status", "proposed")
        if authority != "proposed":
            raise CorpusIngestError(
                f"{path}: {family} has authority_status={authority!r}; only proposed imports are allowed"
            )

    ordered = sorted(variants, key=lambda item: _record_rank(item[0], item[1]))
    merged: dict[str, Any] = {}
    for _path, record in ordered:
        for key, value in record.items():
            if value not in (None, "", [], {}):
                merged[key] = value

    discovery_domains = sorted(
        {
            str(record["discovery_domain"])
            for _path, record in variants
            if record.get("discovery_domain")
        }
    )
    source_files = sorted({_portable_path(path) for path, _record in variants})
    provider_oa = any(bool(record.get("is_open_access_provider_flag")) for _, record in variants)

    access_route = merged.get("access_route")
    if access_route is None:
        access_route = "open_access" if provider_oa else "metadata_only"
    license_status = merged.get("license_status") or merged.get("provider_license") or "unverified"
    redistribution = merged.get("redistribution", "unknown")

    merged.update(
        {
            "source_family_id": family,
            "authority_status": "proposed",
            "access_route": access_route,
            "license_status": license_status,
            "redistribution": redistribution,
            "discovery_domains": discovery_domains,
            "input_files": source_files,
            "input_record_count": len(variants),
        }
    )
    return {key: merged[key] for key in sorted(merged)}


def _explicitly_retracted(record: Mapping[str, Any]) -> bool:
    if record.get("retracted") is True:
        return True
    statuses = {
        str(record.get(key, "")).strip().lower()
        for key in ("status", "source_status", "retraction_status", "quality_status")
    }
    if statuses & {"retracted", "withdrawn", "voided"}:
        return True
    publication_types = {str(item).strip().lower() for item in record.get("publication_types", [])}
    return "retracted publication" in publication_types


def _payload_path(receipt: Mapping[str, Any], receipt_path: Path) -> Path | None:
    for key in ("payload_path", "local_path", "object_path"):
        raw = receipt.get(key)
        if not raw:
            continue
        candidate = Path(str(raw)).expanduser()
        if not candidate.is_absolute():
            beside_receipt = receipt_path.parent / candidate
            candidate = beside_receipt if beside_receipt.exists() else repo_root() / candidate
        return candidate
    return None


def _receipt_digest(receipt: Mapping[str, Any]) -> str | None:
    for key in ("payload_sha256", "sha256", "object_digest"):
        raw = receipt.get(key)
        if raw:
            digest = str(raw).lower().removeprefix("sha256:")
            if len(digest) != 64 or any(ch not in "0123456789abcdef" for ch in digest):
                raise CorpusIngestError(f"invalid receipt payload digest: {raw!r}")
            return digest
    return None


def _load_inputs(
    manifest_paths: Sequence[Path], receipt_paths: Sequence[Path]
) -> tuple[dict[str, list[tuple[Path, Mapping[str, Any]]]], int, dict[str, list[tuple[Path, dict[str, Any]]]]]:
    variants: dict[str, list[tuple[Path, Mapping[str, Any]]]] = defaultdict(list)
    total = 0
    for path in sorted(manifest_paths):
        for record in _read_json_records(path):
            family = _normalise_family(record.get("source_family_id"))
            variants[family].append((path, record))
            total += 1

    receipts: dict[str, list[tuple[Path, dict[str, Any]]]] = defaultdict(list)
    for path in sorted(receipt_paths):
        for receipt in _read_json_records(path):
            family = _normalise_family(receipt.get("source_family_id"))
            receipts[family].append((path, receipt))
    unknown = sorted(set(receipts) - set(variants))
    if unknown:
        raise CorpusIngestError(
            f"{len(unknown)} acquisition receipt family/families are absent from manifests: {unknown[:3]}"
        )
    return variants, total, receipts


def _projection_hash(snapshot: Snapshot) -> str:
    body = {
        subject_id: {
            "attrs": entity.attrs,
            "objects": entity.objects,
            "retracted": entity.retracted,
            "event_types": entity.event_types,
        }
        for subject_id, entity in sorted(snapshot.entities.items())
    }
    return hashlib.sha256(canonical_json(body).encode("utf-8")).hexdigest()


def _snapshot_at_current_head(data_root: DataRoot, ledger: Ledger, store: ObjectStore) -> Snapshot:
    head = ledger.head()
    if head is None:
        raise CorpusIngestError("cannot create a corpus snapshot from an empty ledger")
    current = latest(data_root)
    if current is not None and current.ledger_head_id == head.id:
        return current
    return create_snapshot(
        data_root,
        f"external-training-corpus-proposed-head-{head.id}",
        ledger=ledger,
        store=store,
        authority_epoch=AUTHORITY_EPOCH,
        actor=ACTOR,
        record_event=False,
        notes={
            "authority_status": "proposed",
            "use": "literature evidence only; no Aleph claim promotion",
        },
    )


def ingest_corpus(
    manifest_paths: Iterable[str | Path],
    data_root: DataRoot,
    *,
    receipt_paths: Iterable[str | Path] = (),
) -> IngestResult:
    """Ingest manifests and optional acquisition receipts into ``data_root``.

    Identical source metadata is idempotent.  A changed canonical record appends
    a revision event; an explicit retraction appends a retraction event and can
    never erase the earlier history.
    """

    manifests = [Path(path).resolve() for path in manifest_paths]
    receipts_in = [Path(path).resolve() for path in receipt_paths]
    if not manifests:
        raise CorpusIngestError("at least one manifest is required")
    for path in manifests + receipts_in:
        if not path.is_file():
            raise FileNotFoundError(path)

    variants, input_records, receipts = _load_inputs(manifests, receipts_in)
    data_root.ensure()
    receipt_count = sum(len(items) for items in receipts.values())
    payloads_present = 0
    payloads_unavailable = 0
    appended = 0

    cag = CAG.load()
    authority_statement = (
        "External literature records retain proposed evidence status and recorded access and "
        "licence classes."
    )
    lint = cag.check_statement_against_invariants(authority_statement)
    if lint:
        raise CorpusIngestError("CAG lint rejected adapter authority statement: " + "; ".join(v.describe() for v in lint))

    with Ledger.open(data_root) as ledger:
        store = ObjectStore.open(data_root)
        for family in sorted(variants):
            canonical = _canonical_source(family, variants[family])
            object_digests: list[str] = [
                store.put_json(
                    canonical,
                    kind="external-source-metadata",
                    origin=";".join(canonical["input_files"]),
                    extra={"source_family_id": family, "authority_status": "proposed"},
                )
            ]

            receipt_summaries: list[dict[str, Any]] = []
            for receipt_path, receipt in sorted(
                receipts.get(family, []), key=lambda item: (str(item[0]), canonical_json(item[1]))
            ):
                receipt_object = store.put_json(
                    receipt,
                    kind="external-acquisition-receipt",
                    origin=str(receipt_path),
                    extra={"source_family_id": family},
                )
                object_digests.append(receipt_object)
                payload_digest = _receipt_digest(receipt)
                payload_path = _payload_path(receipt, receipt_path)
                payload_available = False
                if payload_digest is not None and store.exists(payload_digest):
                    if not store.verify(payload_digest):
                        raise CorpusIngestError(f"stored acquisition payload is corrupt: {payload_digest}")
                    payload_available = True
                    object_digests.append(payload_digest)
                elif payload_digest is not None and payload_path is not None and payload_path.is_file():
                    stored = store.put(
                        payload_path,
                        kind="external-source-payload",
                        origin=str(payload_path),
                        extra={"source_family_id": family, "redistribution": receipt.get("redistribution", "unknown")},
                    )
                    if stored != payload_digest:
                        raise CorpusIngestError(
                            f"{receipt_path}: payload digest says {payload_digest}, bytes hash to {stored}"
                        )
                    payload_available = True
                    object_digests.append(stored)
                if payload_digest is not None:
                    if payload_available:
                        payloads_present += 1
                    else:
                        payloads_unavailable += 1
                receipt_summaries.append(
                    {
                        "receipt_object": receipt_object,
                        "payload_sha256": payload_digest,
                        "payload_present": payload_available,
                        "access_route": receipt.get("access_route", canonical["access_route"]),
                        "license_status": receipt.get("license_status", canonical["license_status"]),
                        "redistribution": receipt.get("redistribution", canonical["redistribution"]),
                    }
                )

            subject_id = _subject_id(family)
            previous = ledger.events(subject_id=subject_id, descending=True, limit=1)
            previous_payload = previous[0].payload if previous else {}
            was_retracted = any(
                record.event_type.endswith((".retracted", ".withdrawn", ".voided"))
                for record in ledger.events(subject_id=subject_id)
            )
            payload = {
                "name": canonical.get("title", family),
                "source_family_id": family,
                "status": "proposed",
                "authority_status": "proposed",
                "access_route": canonical["access_route"],
                "license_status": canonical["license_status"],
                "redistribution": canonical["redistribution"],
                "objects": sorted(set(object_digests)),
                "metadata_object": object_digests[0],
                "acquisition_receipts": receipt_summaries,
                "retracted": was_retracted or _explicitly_retracted(canonical),
            }
            comparison = dict(payload)
            if previous_payload == comparison:
                continue
            event_type = EVENT_RETRACTED if payload["retracted"] else (EVENT_REVISED if previous else EVENT_IMPORTED)
            snapshot_date = str(canonical.get("snapshot_date") or canonical.get("access_verified_at") or "2026-08-05")
            timestamp = f"{snapshot_date[:10]}T00:00:00.000000+00:00"
            ledger.append(Event(event_type, subject_id, ACTOR, payload, ts=timestamp))
            appended += 1

        verification = ledger.verify_chain()
        if not verification.ok:
            raise CorpusIngestError(verification.describe())
        snapshot = _snapshot_at_current_head(data_root, ledger, store)
        corrupt = store.verify_all()
        if corrupt:
            raise CorpusIngestError(f"object verification failed for {len(corrupt)} digest(s)")
        if snapshot.missing_objects:
            raise CorpusIngestError(f"snapshot references {len(snapshot.missing_objects)} missing object(s)")
        retracted_count = sum(entity.retracted for entity in snapshot.entities.values())
        head = ledger.head()
        assert head is not None
        return IngestResult(
            input_files=len(manifests),
            input_records=input_records,
            source_families=len(variants),
            duplicate_records=input_records - len(variants),
            events_appended=appended,
            ledger_events=ledger.count(),
            objects=sum(1 for _ in store.iter_digests()),
            acquisition_receipts=receipt_count,
            acquisition_payloads_present=payloads_present,
            acquisition_payloads_unavailable=payloads_unavailable,
            retracted_entities=retracted_count,
            snapshot_id=snapshot.id,
            snapshot_manifest_hash=snapshot.manifest_hash,
            projection_hash=_projection_hash(snapshot),
            ledger_head_hash=head.row_hash,
            chain_verified=verification.ok,
            missing_objects=len(snapshot.missing_objects),
            cag_lint_violations=len(lint),
        )


def verify_idempotent_ingest(
    manifest_paths: Iterable[str | Path],
    data_root: DataRoot,
    *,
    receipt_paths: Iterable[str | Path] = (),
) -> dict[str, Any]:
    manifests = tuple(manifest_paths)
    receipts = tuple(receipt_paths)
    first = ingest_corpus(manifests, data_root, receipt_paths=receipts)
    second = ingest_corpus(manifests, data_root, receipt_paths=receipts)
    stable = (
        second.events_appended == 0
        and first.ledger_events == second.ledger_events
        and first.snapshot_manifest_hash == second.snapshot_manifest_hash
        and first.projection_hash == second.projection_hash
        and first.ledger_head_hash == second.ledger_head_hash
    )
    if not stable:
        raise CorpusIngestError("second ingest changed events, ledger head, or TAG projection")
    return {"first": first.as_dict(), "second": second.as_dict(), "idempotent": True}
