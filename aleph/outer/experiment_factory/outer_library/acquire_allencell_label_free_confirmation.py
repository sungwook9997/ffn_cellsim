#!/usr/bin/env python3
"""Acquire the frozen Allen raw-image confirmation payloads with receipts."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import urllib.request


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def acquire(url: str, relative: str, root: Path) -> dict[str, object]:
    output = root / relative
    output.parent.mkdir(parents=True, exist_ok=True)
    error: Exception | None = None
    for attempt in range(5):
        try:
            head = urllib.request.Request(
                url, method="HEAD", headers={"User-Agent": "Project-Aleph-research/1.0"}
            )
            with urllib.request.urlopen(head, timeout=120) as response:
                size = int(response.headers["Content-Length"])
                etag = response.headers.get("ETag")
                last_modified = response.headers.get("Last-Modified")
            if not output.exists() or output.stat().st_size != size:
                part = output.with_suffix(output.suffix + ".part")
                offset = part.stat().st_size if part.exists() else 0
                if offset > size:
                    part.unlink()
                    offset = 0
                headers = {"User-Agent": "Project-Aleph-research/1.0"}
                if offset:
                    headers["Range"] = f"bytes={offset}-"
                request = urllib.request.Request(
                    url, headers=headers
                )
                with urllib.request.urlopen(request, timeout=600) as response:
                    if offset and response.status != 206:
                        part.unlink()
                        raise ValueError(f"source refused range resume for {relative}")
                    mode = "ab" if offset else "wb"
                    with part.open(mode) as handle:
                        while True:
                            chunk = response.read(4 * 1024 * 1024)
                            if not chunk:
                                break
                            handle.write(chunk)
                if part.stat().st_size != size:
                    raise ValueError(f"short payload for {relative}")
                part.replace(output)
            return {
                "path": relative, "url": url, "size_bytes": size,
                "etag": etag, "last_modified": last_modified,
                "sha256": sha256(output),
            }
        except Exception as exc:
            error = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"failed Allen acquisition {relative}: {error}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    root = args.protocol.parent
    pinned = protocol["candidate_artifacts"]
    checks = {
        "metadata_audit_protocol_sha256": root / "allencell_label_free_metadata_audit_protocol.json",
        "metadata_audit_report_sha256": root / "results/allencell_label_free_metadata_audit.json",
        "candidate_manifest_sha256": args.manifest,
    }
    for key, path in checks.items():
        if sha256(path) != pinned[key]:
            raise ValueError(f"pinned Allen artifact changed: {path}")
    model_paths = {
        "v2": {
            "protocol_sha256": root / "lightmycells_model_training_protocol_v2.json",
            "report_sha256": root / "results/lightmycells_model_training_report_v2.json",
            "checkpoint_sha256": root / "results/lightmycells_model_checkpoint_v2.pt",
        },
        "v3": {
            "protocol_sha256": root / "lightmycells_model_training_protocol_v3.json",
            "report_sha256": root / "results/lightmycells_model_training_report_v3.json",
            "checkpoint_sha256": root / "results/lightmycells_model_checkpoint_v3.pt",
        },
    }
    for version, paths in model_paths.items():
        for key, path in paths.items():
            if sha256(path) != protocol["frozen_models"][version][key]:
                raise ValueError(f"pinned {version} artifact changed: {path}")
    rows = [json.loads(line) for line in args.manifest.read_text().splitlines()]
    if len(rows) != pinned["candidate_count"]:
        raise ValueError("candidate count changed after protocol freeze")
    base = protocol["source_contract"]["base_url"]
    payloads = sorted({row[key] for row in rows for key in ("crop_raw", "fov_path")})
    receipts = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(acquire, base + relative, relative, args.data_root): relative
            for relative in payloads
        }
        for future in as_completed(futures):
            receipts.append(future.result())
            print(json.dumps({
                "completed": len(receipts), "total": len(payloads),
                "path": futures[future],
            }, sort_keys=True), flush=True)
    receipts.sort(key=lambda row: row["path"])
    report = {
        "schema": "aleph.outer_library.allencell_label_free_acquisition.v1",
        "protocol_sha256": sha256(args.protocol),
        "manifest_sha256": sha256(args.manifest),
        "candidate_count": len(rows),
        "payload_count": len(receipts),
        "payload_bytes": sum(int(row["size_bytes"]) for row in receipts),
        "files": receipts,
        "image_payloads_accessed": len(receipts),
        "image_pixels_decoded": 0,
        "single_use_confirmation_consumed": True,
        "model_predictions_made": 0,
        "external_confirmation_complete": False,
        "aleph_parameter_authority": "none",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
