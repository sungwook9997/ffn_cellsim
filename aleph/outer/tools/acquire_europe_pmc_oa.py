#!/usr/bin/env python3
"""Acquire provider-flagged Europe PMC open-access full-text XML.

This tool intentionally refuses candidates unless the discovery record has both
``is_open_access_provider_flag == true`` and a syntactically valid PMCID.  It
uses Europe PMC's official REST fullTextXML endpoint, validates the response as
XML, stores payloads by SHA-256 outside Git, and atomically checkpoints a
versioned receipt JSONL plus a machine-readable summary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import re
import signal
import sys
import tempfile
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


ENDPOINT = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
PMCID_RE = re.compile(r"PMC[1-9][0-9]*\Z")
USER_AGENT = "Project-Aleph-OA-Acquisition/1.0 (polite research client; Europe PMC REST)"
RECEIPT_SCHEMA = "aleph.external_training.oa_acquisition_receipt.v1"
SUMMARY_SCHEMA = "aleph.external_training.oa_acquisition_summary.v1"
TERMINATE = False


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp_name, path)
    except BaseException:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise


def canonical_json(record: dict[str, Any]) -> bytes:
    return (json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n").encode("utf-8")


def read_jsonl(paths: Iterable[Path]) -> list[dict[str, Any]]:
    candidates: dict[str, dict[str, Any]] = {}
    for path in paths:
        with path.open("r", encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, 1):
                if not line.strip():
                    continue
                try:
                    record = json.loads(line)
                except json.JSONDecodeError as exc:
                    raise ValueError(f"{path}:{line_number}: invalid JSON: {exc}") from exc
                if record.get("is_open_access_provider_flag") is not True:
                    continue
                pmcid = str(record.get("pmcid") or "").upper()
                if not PMCID_RE.fullmatch(pmcid):
                    continue
                record = dict(record)
                record["pmcid"] = pmcid
                record["candidate_input"] = path.as_posix()
                existing = candidates.get(pmcid)
                if existing is None:
                    candidates[pmcid] = record
                    continue
                # Deterministically retain the richer duplicate discovery row.
                old_score = sum(value not in (None, "", [], {}) for value in existing.values())
                new_score = sum(value not in (None, "", [], {}) for value in record.values())
                if (new_score, record.get("source_family_id", "")) > (
                    old_score,
                    existing.get("source_family_id", ""),
                ):
                    candidates[pmcid] = record
    return [candidates[key] for key in sorted(candidates, key=lambda value: int(value[3:]))]


def load_receipts(path: Path) -> dict[str, dict[str, Any]]:
    receipts: dict[str, dict[str, Any]] = {}
    if not path.exists():
        return receipts
    with path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            try:
                receipt = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{line_number}: invalid receipt JSON: {exc}") from exc
            pmcid = receipt.get("pmcid")
            if not isinstance(pmcid, str) or not PMCID_RE.fullmatch(pmcid):
                raise ValueError(f"{path}:{line_number}: invalid receipt PMCID")
            receipts[pmcid] = receipt
    return receipts


def strip_namespace(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def validate_and_describe_xml(payload: bytes) -> dict[str, Any]:
    if not payload.lstrip().startswith(b"<"):
        raise ValueError("response is not XML")
    try:
        root = ET.fromstring(payload)
    except ET.ParseError as exc:
        raise ValueError(f"malformed XML: {exc}") from exc
    root_name = strip_namespace(root.tag)
    article_nodes = [node for node in root.iter() if strip_namespace(node.tag) == "article"]
    if root_name not in {"article", "pmc-articleset"} or not article_nodes:
        raise ValueError(f"unexpected XML root/content: {root_name}")

    license_types: set[str] = set()
    license_hrefs: set[str] = set()
    license_texts: list[str] = []
    for node in root.iter():
        if strip_namespace(node.tag) != "license":
            continue
        license_type = node.attrib.get("license-type")
        if license_type:
            license_types.add(license_type.strip())
        for key, value in node.attrib.items():
            if strip_namespace(key) == "href" and value:
                license_hrefs.add(value.strip())
        text = " ".join(" ".join(node.itertext()).split())
        if text:
            license_texts.append(text[:1000])
    return {
        "xml_root": root_name,
        "xml_article_count": len(article_nodes),
        "license_types": sorted(license_types),
        "license_hrefs": sorted(license_hrefs),
        "license_text": " | ".join(license_texts[:3]) or None,
    }


def retrieve(url: str, timeout: float, retries: int, backoff: float) -> tuple[bytes, int, str, int]:
    last_error: BaseException | None = None
    for attempt in range(retries + 1):
        request = urllib.request.Request(
            url,
            headers={"Accept": "application/xml,text/xml", "User-Agent": USER_AGENT},
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                status = int(response.status)
                media_type = response.headers.get_content_type()
                payload = response.read()
                return payload, status, media_type, attempt + 1
        except urllib.error.HTTPError as exc:
            last_error = exc
            # Missing/restricted records will not improve inside a single run.
            if exc.code in {400, 401, 403, 404, 410}:
                break
            retry_after = exc.headers.get("Retry-After") if exc.headers else None
            try:
                wait = float(retry_after) if retry_after else backoff * (2**attempt)
            except ValueError:
                wait = backoff * (2**attempt)
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last_error = exc
            wait = backoff * (2**attempt)
        if attempt < retries:
            time.sleep(wait + random.uniform(0.0, min(0.25, backoff)))
    assert last_error is not None
    raise last_error


def store_object(object_root: Path, payload: bytes) -> tuple[str, Path, bool]:
    digest = hashlib.sha256(payload).hexdigest()
    target = object_root / digest
    if target.exists():
        if target.stat().st_size != len(payload):
            raise RuntimeError(f"content-address collision/size mismatch for {digest}")
        return digest, target, False
    atomic_write(target, payload)
    return digest, target, True


def summary_for(
    candidates: list[dict[str, Any]], receipts: dict[str, dict[str, Any]], started_at: str
) -> dict[str, Any]:
    statuses = Counter(receipt.get("status", "unknown") for receipt in receipts.values())
    successes = [receipt for receipt in receipts.values() if receipt.get("status") == "success"]
    digests = {receipt.get("sha256") for receipt in successes if receipt.get("sha256")}
    total_bytes = sum(int(receipt.get("bytes") or 0) for receipt in successes)
    elapsed = sum(float(receipt.get("elapsed_seconds") or 0.0) for receipt in receipts.values())
    retrieved_times = sorted(
        datetime.fromisoformat(receipt["retrieved_at"].replace("Z", "+00:00"))
        for receipt in receipts.values()
        if receipt.get("retrieved_at")
    )
    acquisition_span = (
        (retrieved_times[-1] - retrieved_times[0]).total_seconds() if len(retrieved_times) > 1 else 0.0
    )
    provider_licenses = Counter(str(receipt.get("provider_license") or "UNKNOWN") for receipt in successes)
    http_failures = Counter(
        str(receipt.get("http_status") or "NO_HTTP_STATUS")
        for receipt in receipts.values()
        if receipt.get("status") != "success"
    )
    return {
        "schema": SUMMARY_SCHEMA,
        "generated_at": utc_now(),
        "run_started_at": started_at,
        "eligible_unique_pmcids": len(candidates),
        "receipt_count": len(receipts),
        "status_counts": dict(sorted(statuses.items())),
        "success_count": len(successes),
        "failure_count": sum(count for status, count in statuses.items() if status != "success"),
        "unique_payload_count": len(digests),
        "duplicate_payload_count": len(successes) - len(digests),
        "validated_xml_bytes": total_bytes,
        "request_elapsed_seconds_sum": round(elapsed, 6),
        "mean_request_seconds": round(elapsed / len(receipts), 6) if receipts else 0.0,
        "successful_megabytes": round(total_bytes / 1_000_000, 6),
        "acquisition_first_retrieved_at": (
            retrieved_times[0].isoformat().replace("+00:00", "Z") if retrieved_times else None
        ),
        "acquisition_last_retrieved_at": (
            retrieved_times[-1].isoformat().replace("+00:00", "Z") if retrieved_times else None
        ),
        "acquisition_span_seconds": acquisition_span,
        "successes_per_wall_second": (
            round(len(successes) / acquisition_span, 6) if acquisition_span > 0 else 0.0
        ),
        "provider_license_counts": dict(sorted(provider_licenses.items())),
        "xml_embedded_license_count": sum(bool(receipt.get("license_text")) for receipt in successes),
        "failure_http_status_counts": dict(sorted(http_failures.items())),
    }


def checkpoint(
    receipt_path: Path,
    summary_path: Path,
    candidates: list[dict[str, Any]],
    receipts: dict[str, dict[str, Any]],
    started_at: str,
) -> None:
    ordered = b"".join(canonical_json(receipts[key]) for key in sorted(receipts))
    atomic_write(receipt_path, ordered)
    atomic_write(summary_path, canonical_json(summary_for(candidates, receipts, started_at)))


def failure_receipt(candidate: dict[str, Any], url: str, started: float, exc: BaseException) -> dict[str, Any]:
    status_code = exc.code if isinstance(exc, urllib.error.HTTPError) else None
    return {
        "schema": RECEIPT_SCHEMA,
        "status": "failed",
        "pmcid": candidate["pmcid"],
        "source_family_id": candidate.get("source_family_id"),
        "doi": candidate.get("doi"),
        "title": candidate.get("title"),
        "candidate_input": candidate.get("candidate_input"),
        "provider_oa_flag": True,
        "provider_license": candidate.get("provider_license"),
        "retrieval_url": url,
        "retrieved_at": utc_now(),
        "http_status": status_code,
        "bytes": 0,
        "sha256": None,
        "media_type": None,
        "license_types": [],
        "license_hrefs": [],
        "license_text": None,
        "access_route": "open_access_provider_flag_but_payload_unavailable",
        "redistribution": "none_no_payload",
        "error_type": type(exc).__name__,
        "error": str(exc)[:1000],
        "elapsed_seconds": round(time.monotonic() - started, 6),
    }


def acquire_candidate(candidate: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    """Acquire one candidate; delay in the worker to bound each request lane."""
    pmcid = candidate["pmcid"]
    url = ENDPOINT.format(pmcid=pmcid)
    started = time.monotonic()
    try:
        payload, status, media_type, attempts = retrieve(
            url, args.timeout_seconds, args.retries, args.backoff_seconds
        )
        xml_details = validate_and_describe_xml(payload)
        digest, _object_path, object_created = store_object(args.object_root, payload)
        return {
            "schema": RECEIPT_SCHEMA,
            "status": "success",
            "pmcid": pmcid,
            "source_family_id": candidate.get("source_family_id"),
            "doi": candidate.get("doi"),
            "title": candidate.get("title"),
            "candidate_input": candidate.get("candidate_input"),
            "provider_oa_flag": True,
            "provider_license": candidate.get("provider_license"),
            "retrieval_url": url,
            "retrieved_at": utc_now(),
            "http_status": status,
            "request_attempts": attempts,
            "media_type": media_type,
            "sha256": digest,
            "bytes": len(payload),
            "object_path": f"data/external_training/objects/{digest}",
            "object_created": object_created,
            "access_route": "open_access",
            "redistribution": "subject_to_recorded_source_license",
            "elapsed_seconds": round(time.monotonic() - started, 6),
            **xml_details,
        }
    except (OSError, ValueError, RuntimeError, urllib.error.URLError) as exc:
        return failure_receipt(candidate, url, started, exc)
    finally:
        if args.delay_seconds:
            time.sleep(args.delay_seconds)


def stop_on_signal(signum: int, _frame: object) -> None:
    global TERMINATE
    TERMINATE = True
    print(f"received signal {signum}; checkpointing after current request", file=sys.stderr)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", type=Path, required=True, help="candidate JSONL (repeatable)")
    parser.add_argument("--object-root", type=Path, required=True, help="ignored content-addressed object directory")
    parser.add_argument("--receipts", type=Path, required=True, help="versioned acquisition receipt JSONL")
    parser.add_argument("--summary", type=Path, required=True, help="versioned acquisition summary JSON")
    parser.add_argument("--limit", type=int, default=0, help="maximum new attempts; 0 means all")
    parser.add_argument("--delay-seconds", type=float, default=0.15, help="minimum delay after each request")
    parser.add_argument("--timeout-seconds", type=float, default=45.0)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--backoff-seconds", type=float, default=1.0)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    parser.add_argument("--workers", type=int, default=1, help="bounded request lanes; use conservatively")
    parser.add_argument("--retry-failures", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.limit < 0 or args.delay_seconds < 0 or args.timeout_seconds <= 0:
        parser.error("limit/delay/timeout values are invalid")
    if args.retries < 0 or args.backoff_seconds < 0 or args.checkpoint_every <= 0 or args.workers <= 0:
        parser.error("retry/backoff/checkpoint values are invalid")
    return args


def main() -> int:
    args = parse_args()
    candidates = read_jsonl(args.input)
    receipts = load_receipts(args.receipts)
    pending = [
        candidate
        for candidate in candidates
        if candidate["pmcid"] not in receipts
        or (args.retry_failures and receipts[candidate["pmcid"]].get("status") != "success")
    ]
    if args.limit:
        pending = pending[: args.limit]
    print(
        json.dumps(
            {
                "eligible_unique_pmcids": len(candidates),
                "existing_receipts": len(receipts),
                "pending_attempts": len(pending),
                "dry_run": args.dry_run,
            },
            sort_keys=True,
        ),
        flush=True,
    )
    if args.dry_run:
        return 0

    signal.signal(signal.SIGINT, stop_on_signal)
    signal.signal(signal.SIGTERM, stop_on_signal)
    started_at = utc_now()
    completed = 0
    try:
        pending_iter = iter(pending)
        with ThreadPoolExecutor(max_workers=args.workers, thread_name_prefix="europe-pmc") as executor:
            active: dict[Future[dict[str, Any]], str] = {}
            for _ in range(args.workers):
                candidate = next(pending_iter, None)
                if candidate is None:
                    break
                active[executor.submit(acquire_candidate, candidate, args)] = candidate["pmcid"]
            while active:
                done, _ = wait(active, return_when=FIRST_COMPLETED)
                for future in done:
                    pmcid = active.pop(future)
                    receipts[pmcid] = future.result()
                    completed += 1
                    if completed % args.checkpoint_every == 0:
                        checkpoint(args.receipts, args.summary, candidates, receipts, started_at)
                        counts = Counter(item.get("status") for item in receipts.values())
                        print(
                            f"checkpoint attempted={completed} total={len(receipts)} status={dict(counts)}",
                            flush=True,
                        )
                    if not TERMINATE:
                        candidate = next(pending_iter, None)
                        if candidate is not None:
                            active[executor.submit(acquire_candidate, candidate, args)] = candidate["pmcid"]
    finally:
        checkpoint(args.receipts, args.summary, candidates, receipts, started_at)
    print(json.dumps(summary_for(candidates, receipts, started_at), sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
