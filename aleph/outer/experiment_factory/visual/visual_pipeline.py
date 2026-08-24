from __future__ import annotations

import argparse
import gzip
import hashlib
import html
import json
import math
import re
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from aleph.outer.experiment_factory.schema.id_helpers import observation_id


SCHEMA = "aleph.external_training.visual_observation_candidate.v1"
RECEIPT_SCHEMA = "aleph.external_training.visual_asset_receipt.v1"
XLINK = "{http://www.w3.org/1999/xlink}href"
USER_AGENT = "Project-Aleph-research-corpus/1.0 (polite OA asset audit)"

MODALITY_RULES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("tfm_map", (r"traction force microscopy", r"traction map", r"traction stress")),
    ("piv_vector_field", (r"particle image velocimetry", r"\bpiv\b", r"velocity vector field", r"flow field")),
    ("afm", (r"atomic force microscop", r"\bafm\b", r"force.indentation")),
    ("western_blot_or_gel", (r"western blot", r"immunoblot", r"\bblot\b", r"\bgel electrophoresis\b")),
    ("pcr_or_qpcr_curve", (r"\bq?pcr\b", r"polymerase chain reaction", r"amplification curve", r"melting curve")),
    ("immunofluorescence", (r"immunofluorescen", r"immunostain", r"fluorescently labeled", r"fluorescence image")),
    ("fluorescence_microscopy", (r"confocal", r"fluorescence microscop", r"wide.?field", r"two.?photon")),
    ("morphology", (r"morpholog", r"cell shape", r"phase.?contrast", r"bright.?field", r"micrograph")),
    ("time_series", (r"time.?lapse", r"over time", r"time series", r"live.?cell")),
)


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def stable_id(prefix: str, *parts: object) -> str:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return f"aleph:{prefix}:{hashlib.sha256(payload).hexdigest()[:16]}"


def clean_text(element: ET.Element | None) -> str:
    if element is None:
        return ""
    return " ".join("".join(element.itertext()).split())


def modality_tags(caption: str) -> list[dict[str, Any]]:
    folded = caption.casefold()
    found: list[dict[str, Any]] = []
    for tag, rules in MODALITY_RULES:
        hits = [rule for rule in rules if re.search(rule, folded)]
        if hits:
            found.append({"tag": tag, "basis": "caption_rule", "rule_hits": len(hits), "authority_status": "proposed"})
    return found


def iter_receipts(paths: Iterable[Path]) -> Iterable[tuple[dict[str, Any], str]]:
    for path in paths:
        acquisition_pass = "second_pass" if "second_pass" in path.parts else "initial"
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                row = json.loads(line)
                if row.get("status") == "success":
                    yield row, acquisition_pass


def figure_records(receipt: dict[str, Any], acquisition_pass: str) -> Iterable[dict[str, Any]]:
    root = ET.parse(receipt["object_path"]).getroot()
    source_id = str(receipt["source_family_id"]).casefold()
    for ordinal, fig in enumerate(root.findall(".//fig")):
        figure_id = fig.get("id") or f"ordinal-{ordinal}"
        caption = clean_text(fig.find("caption"))
        label = clean_text(fig.find("label"))
        graphics = fig.findall(".//graphic")
        experiment_candidate_id = stable_id("experiment-candidate", source_id, f"fig[{figure_id}]", ordinal)
        locator = f"fig[id={figure_id}]"
        observation_id = stable_id("observation-candidate", experiment_candidate_id, "visual_figure", locator, 0)
        graphic_rows = []
        for graphic_ordinal, graphic in enumerate(graphics):
            href = graphic.get(XLINK)
            graphic_rows.append({
                "asset_candidate_id": stable_id("visual-asset-candidate", observation_id, graphic_ordinal),
                "ordinal": graphic_ordinal,
                "exact_locator": f"{locator}/graphic[{graphic_ordinal}]",
                "href": href,
                "href_sha256": sha256_bytes((href or "").encode("utf-8")),
            })
        yield {
            "schema": SCHEMA,
            "authority_status": "proposed",
            "source_family_id": source_id,
            "pmcid": receipt["pmcid"],
            "acquisition_pass": acquisition_pass,
            "source_xml_sha256": receipt["sha256"],
            "experiment_candidate_id": experiment_candidate_id,
            "experiment_link_state": "figure_scoped_candidate_pending_schema_reconciliation",
            "observation_candidate_id": observation_id,
            "observation_type": "visual_figure",
            "figure_id": figure_id,
            "figure_ordinal": ordinal,
            "exact_locator": locator,
            "figure_label_sha256": sha256_bytes(label.encode("utf-8")),
            "caption_sha256": sha256_bytes(caption.encode("utf-8")),
            "caption_bytes": len(caption.encode("utf-8")),
            "graphics": graphic_rows,
            "proposed_modality_tags": modality_tags(caption),
            "asset_state": "not_in_annotation_selection",
            "scientific_value_state": "refused_without_calibration",
        }


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with path.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as handle:
            for row in rows:
                handle.write(canonical_json(row))
    else:
        with path.open("wb") as handle:
            for row in rows:
                handle.write(canonical_json(row))


def selected_source_ids(path: Path) -> set[str]:
    values: set[str] = set()
    for row in read_jsonl(path):
        source = row.get("source_family_id") or row.get("source_id")
        if source:
            values.add(str(source).casefold())
    return values


def build_index(args: argparse.Namespace) -> None:
    records: list[dict[str, Any]] = []
    for receipt, acquisition_pass in iter_receipts(args.receipts):
        records.extend(figure_records(receipt, acquisition_pass))
    records.sort(key=lambda row: (row["source_family_id"], row["figure_ordinal"]))
    selected = selected_source_ids(args.selection) if args.selection else set()
    for row in records:
        if row["source_family_id"] in selected:
            row["asset_state"] = "selected_pending_acquisition"
    write_jsonl(args.output, records)
    summary = {
        "schema": "aleph.external_training.visual_index_summary.v1",
        "authority_status": "proposed",
        "records": len(records),
        "sources": len({r["source_family_id"] for r in records}),
        "figures": len({(r["source_family_id"], r["figure_id"]) for r in records}),
        "graphic_href_present": sum(bool(g["href"]) for r in records for g in r["graphics"]),
        "selected_sources_requested": len(selected),
        "selection_sha256": sha256_bytes(args.selection.read_bytes()) if args.selection else None,
        "selected_sources_with_figures": len(selected & {r["source_family_id"] for r in records}),
        "selected_records": sum(r["asset_state"] == "selected_pending_acquisition" for r in records),
        "modality_record_counts": dict(sorted(Counter(tag["tag"] for r in records for tag in r["proposed_modality_tags"]).items())),
        "index_sha256": sha256_bytes(args.output.read_bytes()),
        "caption_text_committed": False,
        "scientific_values_from_pixels": False,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_bytes(canonical_json(summary))


def request_bytes(url: str, timeout: float, attempts: int = 3) -> tuple[bytes | None, int | None, str | None, str | None]:
    last_status: int | None = None
    last_error: str | None = None
    for attempt in range(attempts):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "text/html,image/*"})
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read(), response.status, response.headers.get_content_type(), None
        except urllib.error.HTTPError as exc:
            last_status, last_error = exc.code, f"http_{exc.code}"
            if exc.code not in {429, 500, 502, 503, 504}:
                break
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = type(exc).__name__
        if attempt + 1 < attempts:
            time.sleep(0.5 * (2**attempt))
    return None, last_status, None, last_error


def html_asset_map(page: bytes) -> dict[str, str]:
    text = page.decode("utf-8", errors="replace")
    urls = re.findall(r'(?:src|href)=["\'](https://cdn\.ncbi\.nlm\.nih\.gov/pmc/blobs/[^"\']+)["\']', text)
    result: dict[str, str] = {}
    for raw in urls:
        url = html.unescape(raw)
        result.setdefault(urllib.parse.unquote(url.rsplit("/", 1)[-1]), url)
    return result


def resolve_asset_url(mapping: dict[str, str], href: str | None) -> str | None:
    if not href:
        return None
    base = urllib.parse.unquote(href.rsplit("/", 1)[-1])
    if base in mapping:
        return mapping[base]
    stem = base.rsplit(".", 1)[0].casefold()
    alternatives = sorted(
        ((name, url) for name, url in mapping.items() if name.rsplit(".", 1)[0].casefold() == stem),
        key=lambda item: ((item[0].rsplit(".", 1)[-1].casefold() not in {"jpg", "jpeg", "png", "tif", "tiff"}), item[0]),
    )
    return alternatives[0][1] if alternatives else None


def acquire_asset(
    row: dict[str, Any],
    graphic: dict[str, Any],
    mapping: dict[str, str],
    page: bytes | None,
    page_url: str,
    page_status: int | None,
    page_media: str | None,
    page_error: str | None,
    timeout: float,
    object_dir: Path,
) -> dict[str, Any]:
    href = graphic.get("href")
    asset_url = resolve_asset_url(mapping, href)
    payload = None
    status = None
    media_type = None
    error = None
    if not href:
        state, error = "missing_graphic_href", "xml_graphic_href_absent"
    elif not page:
        state, error = "article_page_unavailable", page_error
    elif not asset_url:
        state, error = "asset_url_unresolved", "href_not_found_on_official_article_page"
    else:
        payload, status, media_type, error = request_bytes(asset_url, timeout)
        if payload and media_type and media_type.startswith("image/"):
            state = "acquired"
        elif payload:
            state, error, payload = "unexpected_media_type", f"media_type:{media_type}", None
        else:
            state = "asset_unavailable"
    digest = sha256_bytes(payload) if payload else None
    if payload and digest:
        target = object_dir / digest
        if not target.exists():
            target.write_bytes(payload)
    return {
        "schema": RECEIPT_SCHEMA,
        "authority_status": "proposed",
        "observation_candidate_id": row["observation_candidate_id"],
        "asset_candidate_id": graphic["asset_candidate_id"],
        "exact_locator": graphic["exact_locator"],
        "source_family_id": row["source_family_id"],
        "pmcid": row["pmcid"],
        "graphic_href_sha256": graphic["href_sha256"],
        "asset_state": state,
        "official_article_page": page_url,
        "article_page_http_status": page_status,
        "article_page_media_type": page_media,
        "retrieval_url": asset_url,
        "http_status": status,
        "media_type": media_type,
        "bytes": len(payload) if payload else 0,
        "sha256": digest,
        "local_object": f"data/external_training/experiment_factory/visual/objects/{digest}" if digest else None,
        "license_inheritance": "source_article_OA_receipt_subject_to_recorded_license",
        "error": error,
    }


def acquire(args: argparse.Namespace) -> None:
    records = read_jsonl(args.index)
    selected = selected_source_ids(args.selection)
    by_pmcid: dict[str, list[dict[str, Any]]] = {}
    for row in records:
        if row["source_family_id"] in selected:
            by_pmcid.setdefault(row["pmcid"], []).append(row)
    args.object_dir.mkdir(parents=True, exist_ok=True)
    retryable = {"article_page_unavailable", "asset_url_unresolved", "asset_unavailable"}
    prior = read_jsonl(args.output) if args.output.exists() else []
    for prior_row in prior:
        prior_row.pop("elapsed_seconds", None)
    receipts: list[dict[str, Any]] = [row for row in prior if row["asset_state"] not in retryable]
    completed = {row["asset_candidate_id"] for row in receipts}
    for article_number, pmcid in enumerate(sorted(by_pmcid)):
        page_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
        ranked_rows = sorted(
            by_pmcid[pmcid],
            key=lambda row: (not bool(row["proposed_modality_tags"]), row["figure_ordinal"], row["figure_id"]),
        )
        active_ids = {row["observation_candidate_id"] for row in ranked_rows[: args.max_figures_per_source]}
        pending: list[tuple[dict[str, Any], dict[str, Any]]] = []
        deferred: list[dict[str, Any]] = []
        for row in ranked_rows:
            fallback = [{
                "asset_candidate_id": stable_id("visual-asset-candidate", row["observation_candidate_id"], 0),
                "href": None,
                "href_sha256": sha256_bytes(b""),
                "exact_locator": f'{row["exact_locator"]}/graphic[0]',
            }]
            graphics = row["graphics"] or fallback
            for graphic_number, graphic in enumerate(graphics):
                if graphic["asset_candidate_id"] in completed:
                    continue
                if row["observation_candidate_id"] in active_ids and graphic_number == 0:
                    pending.append((row, graphic))
                else:
                    deferred.append({
                        "schema": RECEIPT_SCHEMA,
                        "authority_status": "proposed",
                        "observation_candidate_id": row["observation_candidate_id"],
                        "asset_candidate_id": graphic["asset_candidate_id"],
                        "exact_locator": graphic["exact_locator"],
                        "source_family_id": row["source_family_id"],
                        "pmcid": pmcid,
                        "graphic_href_sha256": graphic["href_sha256"],
                        "asset_state": "deferred_per_source_figure_budget",
                        "official_article_page": page_url,
                        "article_page_http_status": None,
                        "article_page_media_type": None,
                        "retrieval_url": None,
                        "http_status": None,
                        "media_type": None,
                        "bytes": 0,
                        "sha256": None,
                        "local_object": None,
                        "license_inheritance": "source_article_OA_receipt_subject_to_recorded_license",
                        "error": f"deterministic_max_figures_per_source:{args.max_figures_per_source}",
                    })
        receipts.extend(deferred)
        completed.update(row["asset_candidate_id"] for row in deferred)
        page = None
        page_status = None
        page_media = None
        page_error = None
        mapping: dict[str, str] = {}
        if pending:
            for page_attempt in range(1):
                page, page_status, page_media, page_error = request_bytes(page_url, args.timeout, attempts=1)
                mapping = html_asset_map(page or b"")
                if mapping:
                    break
                page_error = page_error or "official_article_page_has_no_asset_links"
        def fetch(pair: tuple[dict[str, Any], dict[str, Any]]) -> dict[str, Any]:
            return acquire_asset(
                pair[0], pair[1], mapping, page, page_url, page_status, page_media,
                page_error, args.timeout, args.object_dir,
            )
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            new_receipts = list(pool.map(fetch, pending))
        receipts.extend(new_receipts)
        completed.update(row["asset_candidate_id"] for row in new_receipts)
        receipts.sort(key=lambda record: record["asset_candidate_id"])
        write_jsonl(args.output, receipts)
        if article_number + 1 < len(by_pmcid):
            time.sleep(args.delay)
    receipts.sort(key=lambda row: row["asset_candidate_id"])
    write_jsonl(args.output, receipts)


def decode_rgb(path: Path, size: int = 128) -> np.ndarray:
    cmd = [
        "/opt/homebrew/bin/ffmpeg", "-v", "error", "-i", str(path),
        "-vf", f"scale={size}:{size}:force_original_aspect_ratio=decrease,pad={size}:{size}:(ow-iw)/2:(oh-ih)/2:color=white",
        "-f", "rawvideo", "-pix_fmt", "rgb24", "-",
    ]
    proc = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    expected = size * size * 3
    if len(proc.stdout) != expected:
        raise ValueError(f"decoded_byte_count:{len(proc.stdout)}")
    return np.frombuffer(proc.stdout, dtype=np.uint8).reshape(size, size, 3)


def probe_image(path: Path) -> dict[str, Any]:
    cmd = [
        "/opt/homebrew/bin/ffprobe", "-v", "error", "-select_streams", "v:0",
        "-show_entries", "stream=width,height,pix_fmt,nb_frames", "-of", "json", str(path),
    ]
    payload = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
    streams = json.loads(payload).get("streams", [])
    if not streams:
        raise ValueError("no_visual_stream")
    stream = streams[0]
    frames = stream.get("nb_frames")
    return {
        "width": int(stream["width"]),
        "height": int(stream["height"]),
        "source_pixel_format": stream.get("pix_fmt"),
        "source_frame_count": int(frames) if isinstance(frames, str) and frames.isdigit() else None,
    }


def runs(mask: np.ndarray) -> list[tuple[int, int]]:
    padded = np.pad(mask.astype(np.int8), (1, 1))
    changes = np.flatnonzero(np.diff(padded))
    return [(int(changes[i]), int(changes[i + 1])) for i in range(0, len(changes), 2)]


def panels(gray: np.ndarray) -> list[dict[str, float]]:
    size = gray.shape[0]
    cuts: dict[str, list[int]] = {"x": [], "y": []}
    for axis, means, stds in (
        ("x", gray.mean(axis=0), gray.std(axis=0)),
        ("y", gray.mean(axis=1), gray.std(axis=1)),
    ):
        mask = (means >= 248.0) & (stds <= 3.0)
        for start, stop in runs(mask):
            midpoint = (start + stop - 1) // 2
            if start > 0 and stop < size and stop - start >= 2 and 8 <= midpoint <= size - 9:
                cuts[axis].append(midpoint)
    xs = [0] + cuts["x"][:3] + [size]
    ys = [0] + cuts["y"][:3] + [size]
    boxes = []
    for y0, y1 in zip(ys, ys[1:]):
        for x0, x1 in zip(xs, xs[1:]):
            if x1 - x0 >= 12 and y1 - y0 >= 12:
                boxes.append({"x0": x0 / size, "y0": y0 / size, "x1": x1 / size, "y1": y1 / size})
    return boxes[:16] if len(boxes) > 1 else []


def descriptor(path: Path) -> dict[str, Any]:
    probe = probe_image(path)
    width, height = probe["width"], probe["height"]
    rgb = decode_rgb(path)
    gray = np.rint(0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2]).astype(np.uint8)
    gx = np.abs(np.diff(gray.astype(np.int16), axis=1))
    gy = np.abs(np.diff(gray.astype(np.int16), axis=0))
    histogram = np.bincount(gray.ravel(), minlength=256).astype(np.float64)
    probs = histogram[histogram > 0] / gray.size
    channel_std = rgb.reshape(-1, 3).std(axis=0)
    q = np.quantile(gray, [0.05, 0.25, 0.5, 0.75, 0.95], method="linear")
    return {
        "descriptor_version": "cpu_ffmpeg_128_v1",
        "width": width,
        "height": height,
        "aspect_ratio": round(width / height, 8),
        "decoded_channels": 3,
        "source_pixel_format": probe["source_pixel_format"],
        "source_frame_count": probe["source_frame_count"],
        "colorfulness_channel_std_mean": round(float(channel_std.mean()), 6),
        "gray_mean": round(float(gray.mean()), 6),
        "gray_std": round(float(gray.std()), 6),
        "gray_quantiles_05_25_50_75_95": [round(float(x), 6) for x in q],
        "entropy_bits": round(float(-(probs * np.log2(probs)).sum()), 6),
        "edge_density_absdiff_gt_24": round(float(((gx > 24).mean() + (gy > 24).mean()) / 2), 8),
        "panel_candidates": panels(gray),
        "ocr_state": "unavailable_no_local_ocr_backend",
        "scientific_value_state": "refused_without_calibration",
    }


def describe(args: argparse.Namespace) -> None:
    receipts = read_jsonl(args.receipts)
    digests = sorted({receipt["sha256"] for receipt in receipts if receipt.get("asset_state") == "acquired" and receipt.get("sha256")})
    def compute(digest: str) -> tuple[str, dict[str, Any] | None, str | None]:
        try:
            return digest, descriptor(args.object_dir / digest), None
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            return digest, None, f"{type(exc).__name__}:{exc}"
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        computed = {digest: (value, error) for digest, value, error in pool.map(compute, digests)}
    output: list[dict[str, Any]] = []
    for receipt in receipts:
        row = {
            "schema": "aleph.external_training.visual_descriptor.v1",
            "authority_status": "proposed",
            "source_family_id": receipt["source_family_id"],
            "pmcid": receipt["pmcid"],
            "observation_candidate_id": receipt["observation_candidate_id"],
            "asset_candidate_id": receipt["asset_candidate_id"],
            "descriptor_id": stable_id("visual-descriptor", receipt["asset_candidate_id"], receipt.get("sha256")),
            "asset_sha256": receipt.get("sha256"),
            "descriptor_state": "not_acquired",
            "descriptor": None,
            "error": None,
        }
        if receipt.get("asset_state") == "acquired" and receipt.get("sha256"):
            value, error = computed[receipt["sha256"]]
            if value is not None:
                row["descriptor"] = value
                row["descriptor_state"] = "described"
            else:
                row["descriptor_state"] = "decode_failed"
                row["error"] = error
        output.append(row)
    output.sort(key=lambda row: row["asset_candidate_id"])
    write_jsonl(args.output, output)


def audit(args: argparse.Namespace) -> None:
    index = read_jsonl(args.index)
    receipts = read_jsonl(args.receipts)
    descriptors = read_jsonl(args.descriptors)
    ids = [r["observation_candidate_id"] for r in index]
    receipt_ids = [r["asset_candidate_id"] for r in receipts]
    descriptor_ids = [r["asset_candidate_id"] for r in descriptors]
    acquired = [r for r in receipts if r["asset_state"] == "acquired"]
    asset_tags: dict[str, list[str]] = {}
    selected_observations: set[str] = set()
    selected_asset_ids: set[str] = set()
    for visual in index:
        tags = [tag["tag"] for tag in visual["proposed_modality_tags"]]
        for graphic in visual["graphics"]:
            asset_tags[graphic["asset_candidate_id"]] = tags
        if visual["asset_state"] == "selected_pending_acquisition":
            selected_observations.add(visual["observation_candidate_id"])
            if visual["graphics"]:
                selected_asset_ids.update(graphic["asset_candidate_id"] for graphic in visual["graphics"])
            else:
                selected_asset_ids.add(stable_id("visual-asset-candidate", visual["observation_candidate_id"], 0))
    digest_failures = 0
    for row in acquired:
        path = args.object_dir / row["sha256"]
        if not path.exists() or sha256_bytes(path.read_bytes()) != row["sha256"]:
            digest_failures += 1
    summary = {
        "schema": "aleph.external_training.visual_audit.v1",
        "authority_status": "proposed",
        "index_records": len(index),
        "index_unique_ids": len(set(ids)),
        "receipt_records": len(receipts),
        "receipt_unique_ids": len(set(receipt_ids)),
        "selected_asset_identities": len(selected_asset_ids),
        "receipt_id_set_equal_selected_assets": set(receipt_ids) == selected_asset_ids,
        "missing_selected_asset_receipts": len(selected_asset_ids - set(receipt_ids)),
        "extra_asset_receipts": len(set(receipt_ids) - selected_asset_ids),
        "descriptor_records": len(descriptors),
        "descriptor_unique_ids": len(set(descriptor_ids)),
        "acquisition_states": dict(sorted(Counter(r["asset_state"] for r in receipts).items())),
        "descriptor_states": dict(sorted(Counter(r["descriptor_state"] for r in descriptors).items())),
        "acquired_sources": len({r["source_family_id"] for r in acquired}),
        "receipt_sources": len({r["source_family_id"] for r in receipts}),
        "selected_figures": len(selected_observations),
        "figures_with_acquired_asset": len({r["observation_candidate_id"] for r in acquired}),
        "acquired_bytes": sum(r["bytes"] for r in acquired),
        "acquired_unique_objects": len({r["sha256"] for r in acquired}),
        "acquired_modality_receipt_counts": dict(sorted(Counter(tag for r in acquired for tag in asset_tags.get(r["asset_candidate_id"], [])).items())),
        "described_images_with_panel_candidates": sum(bool(row.get("descriptor", {}).get("panel_candidates")) for row in descriptors if row.get("descriptor")),
        "panel_candidate_boxes": sum(len(row.get("descriptor", {}).get("panel_candidates", [])) for row in descriptors if row.get("descriptor")),
        "official_asset_url_violations": sum(
            bool(r.get("retrieval_url")) and not str(r["retrieval_url"]).startswith("https://cdn.ncbi.nlm.nih.gov/pmc/blobs/")
            for r in receipts
        ),
        "non_proposed_records": sum(r.get("authority_status") != "proposed" for r in index + receipts + descriptors),
        "object_digest_failures": digest_failures,
        "receipt_descriptor_id_set_equal": set(receipt_ids) == set(descriptor_ids),
        "duplicate_index_ids": len(ids) - len(set(ids)),
        "duplicate_receipt_ids": len(receipt_ids) - len(set(receipt_ids)),
        "index_sha256": sha256_bytes(args.index.read_bytes()),
        "receipts_sha256": sha256_bytes(args.receipts.read_bytes()),
        "descriptors_sha256": sha256_bytes(args.descriptors.read_bytes()),
        "scientific_values_extracted_from_pixels": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical_json(summary))


def normalized_figure_id(locator_id: str) -> str | None:
    match = re.fullmatch(r"figure:(.+?)/caption", locator_id, flags=re.IGNORECASE)
    return match.group(1).casefold() if match else None


def link_experiments(args: argparse.Namespace) -> None:
    links: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for record in read_jsonl(args.figure_join_index):
        source_id = str(record.get("source_family_id", "")).casefold()
        figure_id = normalized_figure_id(str(record.get("figure_locator_id", "")))
        if source_id and figure_id:
            links.setdefault((source_id, figure_id), []).append(record)
    output: list[dict[str, Any]] = []
    matched_experiments: set[str] = set()
    for visual in read_jsonl(args.index):
        key = (visual["source_family_id"].casefold(), visual["figure_id"].casefold())
        candidates = links.get(key, [])
        hash_matched = [row for row in candidates if row.get("caption_text_sha256") == visual["caption_sha256"]]
        experiment_ids = sorted({row["experiment_record_id"] for row in hash_matched})
        matched_experiments.update(experiment_ids)
        output.append({
            "schema": "aleph.external_training.visual_experiment_link.v1",
            "authority_status": "proposed",
            "source_family_id": visual["source_family_id"],
            "figure_id": visual["figure_id"],
            "exact_locator": visual["exact_locator"],
            "visual_observation_candidate_id": visual["observation_candidate_id"],
            "figure_scoped_experiment_candidate_id": visual["experiment_candidate_id"],
            "linked_experiment_record_ids": experiment_ids,
            "canonical_visual_observation_ids": [
                observation_id(exp_id, "image", visual["exact_locator"], ordinal)
                for ordinal, exp_id in enumerate(experiment_ids)
            ],
            "link_state": (
                "linked_exact_source_figure_and_caption_hash" if experiment_ids
                else "caption_hash_mismatch_retained_figure_scoped_candidate" if candidates
                else "unmatched_retained_figure_scoped_candidate"
            ),
            "asset_candidate_ids": [graphic["asset_candidate_id"] for graphic in visual["graphics"]],
            "caption_sha256": visual["caption_sha256"],
        })
    write_jsonl(args.output, output)
    summary = {
        "schema": "aleph.external_training.visual_experiment_link_summary.v1",
        "authority_status": "proposed",
        "visual_figures": len(output),
        "linked_visual_figures": sum(bool(row["linked_experiment_record_ids"]) for row in output),
        "unmatched_visual_figures": sum(not row["linked_experiment_record_ids"] for row in output),
        "linked_experiment_records": len(matched_experiments),
        "source_figure_keys_in_extraction": len(links),
        "caption_hash_mismatch_visual_figures": sum(row["link_state"].startswith("caption_hash_mismatch") for row in output),
        "figure_join_input_sha256": sha256_bytes(args.figure_join_index.read_bytes()),
        "link_basis": "exact_casefolded_source_family_id_and_JATS_figure_id_and_caption_sha256",
        "output_sha256": sha256_bytes(args.output.read_bytes()),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_bytes(canonical_json(summary))


def merge_receipts(args: argparse.Namespace) -> None:
    priority = {
        "acquired": 5,
        "unexpected_media_type": 4,
        "asset_url_unresolved": 3,
        "asset_unavailable": 3,
        "article_page_unavailable": 2,
        "missing_graphic_href": 2,
        "deferred_per_source_figure_budget": 1,
    }
    chosen: dict[str, dict[str, Any]] = {}
    for path in args.input:
        for row in read_jsonl(path):
            key = row["asset_candidate_id"]
            previous = chosen.get(key)
            if previous is None or priority.get(row["asset_state"], 0) > priority.get(previous["asset_state"], 0):
                chosen[key] = row
    write_jsonl(args.output, (chosen[key] for key in sorted(chosen)))


def finalize_receipts(args: argparse.Namespace) -> None:
    receipts = {row["asset_candidate_id"]: row for row in read_jsonl(args.receipts)}
    for visual in read_jsonl(args.index):
        if visual["asset_state"] != "selected_pending_acquisition":
            continue
        fallback = [{
            "asset_candidate_id": stable_id("visual-asset-candidate", visual["observation_candidate_id"], 0),
            "href": None,
            "href_sha256": sha256_bytes(b""),
            "exact_locator": f'{visual["exact_locator"]}/graphic[0]',
        }]
        for graphic in visual["graphics"] or fallback:
            if graphic["asset_candidate_id"] in receipts:
                continue
            receipts[graphic["asset_candidate_id"]] = {
                "schema": RECEIPT_SCHEMA,
                "authority_status": "proposed",
                "observation_candidate_id": visual["observation_candidate_id"],
                "asset_candidate_id": graphic["asset_candidate_id"],
                "exact_locator": graphic["exact_locator"],
                "source_family_id": visual["source_family_id"],
                "pmcid": visual["pmcid"],
                "graphic_href_sha256": graphic["href_sha256"],
                "asset_state": "official_asset_unresolved_after_bounded_attempt",
                "official_article_page": f'https://pmc.ncbi.nlm.nih.gov/articles/{visual["pmcid"]}/',
                "article_page_http_status": None,
                "article_page_media_type": None,
                "retrieval_url": None,
                "http_status": None,
                "media_type": None,
                "bytes": 0,
                "sha256": None,
                "local_object": None,
                "license_inheritance": "source_article_OA_receipt_subject_to_recorded_license",
                "error": "bounded_official_PMC_attempt_completed_without_resolved_asset",
            }
    write_jsonl(args.output, (receipts[key] for key in sorted(receipts)))


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(dest="command", required=True)
    build = commands.add_parser("build-index")
    build.add_argument("--receipts", type=Path, action="append", required=True)
    build.add_argument("--selection", type=Path)
    build.add_argument("--output", type=Path, required=True)
    build.add_argument("--summary", type=Path, required=True)
    build.set_defaults(func=build_index)
    get = commands.add_parser("acquire")
    get.add_argument("--index", type=Path, required=True)
    get.add_argument("--selection", type=Path, required=True)
    get.add_argument("--object-dir", type=Path, required=True)
    get.add_argument("--output", type=Path, required=True)
    get.add_argument("--delay", type=float, default=0.2)
    get.add_argument("--timeout", type=float, default=30.0)
    get.add_argument("--workers", type=int, choices=range(1, 5), default=3)
    get.add_argument("--max-figures-per-source", type=int, choices=range(1, 9), default=2)
    get.set_defaults(func=acquire)
    desc = commands.add_parser("describe")
    desc.add_argument("--receipts", type=Path, required=True)
    desc.add_argument("--object-dir", type=Path, required=True)
    desc.add_argument("--output", type=Path, required=True)
    desc.add_argument("--workers", type=int, choices=range(1, 5), default=3)
    desc.set_defaults(func=describe)
    check = commands.add_parser("audit")
    check.add_argument("--index", type=Path, required=True)
    check.add_argument("--receipts", type=Path, required=True)
    check.add_argument("--descriptors", type=Path, required=True)
    check.add_argument("--object-dir", type=Path, required=True)
    check.add_argument("--output", type=Path, required=True)
    check.set_defaults(func=audit)
    link = commands.add_parser("link-experiments")
    link.add_argument("--index", type=Path, required=True)
    link.add_argument("--figure-join-index", type=Path, required=True)
    link.add_argument("--output", type=Path, required=True)
    link.add_argument("--summary", type=Path, required=True)
    link.set_defaults(func=link_experiments)
    merge = commands.add_parser("merge-receipts")
    merge.add_argument("--input", type=Path, action="append", required=True)
    merge.add_argument("--output", type=Path, required=True)
    merge.set_defaults(func=merge_receipts)
    finish = commands.add_parser("finalize-receipts")
    finish.add_argument("--index", type=Path, required=True)
    finish.add_argument("--receipts", type=Path, required=True)
    finish.add_argument("--output", type=Path, required=True)
    finish.set_defaults(func=finalize_receipts)
    return root


def main() -> None:
    args = parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
