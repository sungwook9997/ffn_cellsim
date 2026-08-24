#!/usr/bin/env python3
"""Deterministic, resumable stdlib JATS extraction for local OA objects."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

SCHEMA = "aleph.external_training.jats_chunk.v1"
EXTRACTOR_VERSION = "1"
MAX_CHARS = 1800
TAG_DIMENSIONS = (
    "cell_type_state",
    "measurement_modality",
    "mechanics_observable",
    "perturbation",
)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def digest_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def normalized_text(element: ET.Element | None) -> str:
    if element is None:
        return ""
    return " ".join("".join(element.itertext()).split())


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        handle.write(data)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def load_vocabulary(path: Path) -> dict[str, dict[str, tuple[re.Pattern[str], ...]]]:
    raw = json.loads(path.read_text("utf-8"))["dimensions"]
    compiled: dict[str, dict[str, tuple[re.Pattern[str], ...]]] = {}
    for dimension in TAG_DIMENSIONS:
        compiled[dimension] = {}
        for tag, phrases in raw[dimension].items():
            compiled[dimension][tag] = tuple(
                re.compile(r"(?<![a-z0-9])" + re.escape(phrase.lower()) + r"(?![a-z0-9])")
                for phrase in phrases
            )
    return compiled


def infer_tags(text: str, vocabulary: dict[str, dict[str, tuple[re.Pattern[str], ...]]]) -> dict[str, list[str]]:
    lowered = text.lower()
    return {
        dimension: sorted(
            tag for tag, patterns in vocabulary[dimension].items()
            if any(pattern.search(lowered) for pattern in patterns)
        )
        for dimension in TAG_DIMENSIONS
    }


def classify_section(title: str, sec_type: str, kind: str = "section") -> str:
    if kind == "figure":
        return "figure_caption"
    if kind == "table":
        return "table_caption"
    value = f"{sec_type} {title}".lower()
    rules = (
        ("abstract", ("abstract", "summary")),
        ("introduction", ("intro", "background")),
        ("methods", ("method", "material", "experimental procedure")),
        ("results", ("result", "finding")),
        ("discussion", ("discussion",)),
        ("conclusion", ("conclusion", "concluding")),
        ("supplementary", ("supplement", "appendix")),
    )
    for classification, needles in rules:
        if any(needle in value for needle in needles):
            return classification
    return "other"


def text_without_nested_sections(section: ET.Element) -> str:
    pieces: list[str] = []

    def visit(node: ET.Element, root: bool = False) -> None:
        if not root and local_name(node.tag) == "sec":
            return
        if local_name(node.tag) in {"fig", "table-wrap", "ref-list"}:
            return
        if node.text:
            pieces.append(node.text)
        for child in node:
            visit(child)
            if child.tail:
                pieces.append(child.tail)

    visit(section, root=True)
    return " ".join(" ".join(pieces).split())


def split_chunks(text: str, max_chars: int = MAX_CHARS) -> list[str]:
    text = " ".join(text.split())
    if not text:
        return []
    sentences = re.split(r"(?<=[.!?])\s+", text)
    chunks: list[str] = []
    current = ""
    for sentence in sentences:
        if len(sentence) > max_chars:
            if current:
                chunks.append(current)
                current = ""
            chunks.extend(sentence[i:i + max_chars] for i in range(0, len(sentence), max_chars))
        elif not current:
            current = sentence
        elif len(current) + 1 + len(sentence) <= max_chars:
            current += " " + sentence
        else:
            chunks.append(current)
            current = sentence
    if current:
        chunks.append(current)
    return chunks


def descendants(root: ET.Element, name: str) -> Iterable[ET.Element]:
    return (element for element in root.iter() if local_name(element.tag) == name)


def first_descendant(root: ET.Element, name: str) -> ET.Element | None:
    return next(descendants(root, name), None)


def extract_payload(payload: bytes, receipt: dict[str, Any], vocabulary: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, int]]:
    root = ET.fromstring(payload)
    if local_name(root.tag) not in {"article", "article-set"}:
        raise ValueError(f"unsupported XML root: {local_name(root.tag)}")
    article = root if local_name(root.tag) == "article" else first_descendant(root, "article")
    if article is None:
        raise ValueError("no article element")
    pmcid = receipt["pmcid"]
    family = receipt["source_family_id"]
    payload_sha = receipt["sha256"]
    article_locator = f"pmcid:{pmcid}"
    records: list[dict[str, Any]] = []

    def emit(text: str, locator: str, title: str, section_class: str) -> None:
        for ordinal, chunk_text in enumerate(split_chunks(text)):
            base = {
                "schema": SCHEMA,
                "authority_status": "proposed",
                "source_family_id": family,
                "pmcid": pmcid,
                "payload_sha256": payload_sha,
                "article_locator": article_locator,
                "section_locator": locator,
                "section_title": title,
                "section_class": section_class,
                "chunk_ordinal": ordinal,
                "text_sha256": digest_bytes(chunk_text.encode("utf-8")),
                "text": chunk_text,
                "tags": infer_tags(chunk_text, vocabulary),
            }
            hash_basis = dict(base)
            hash_basis.pop("text")
            base["chunk_sha256"] = digest_bytes(canonical_bytes(hash_basis))
            records.append(base)

    abstract = first_descendant(article, "abstract")
    if abstract is not None:
        emit(normalized_text(abstract), "abstract:0", "Abstract", "abstract")

    sections = list(descendants(article, "sec"))
    for index, section in enumerate(sections):
        title_element = next((child for child in section if local_name(child.tag) == "title"), None)
        title = normalized_text(title_element)
        identifier = section.attrib.get("id") or str(index)
        emit(
            text_without_nested_sections(section),
            f"section:{identifier}",
            title,
            classify_section(title, section.attrib.get("sec-type", "")),
        )

    figures = list(descendants(article, "fig"))
    for index, figure in enumerate(figures):
        caption = first_descendant(figure, "caption")
        emit(normalized_text(caption), f"figure:{figure.attrib.get('id') or index}", "", "figure_caption")

    tables = list(descendants(article, "table-wrap"))
    for index, table in enumerate(tables):
        caption = first_descendant(table, "caption")
        emit(normalized_text(caption), f"table:{table.attrib.get('id') or index}", "", "table_caption")

    counts = {"sections": len(sections) + (1 if abstract is not None else 0), "figures": len(figures), "tables": len(tables)}
    return records, counts


def successful_receipts(path: Path) -> list[dict[str, Any]]:
    by_family: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            receipt = json.loads(line)
            if receipt.get("status") != "success":
                continue
            if receipt.get("schema") != "aleph.external_training.oa_acquisition_receipt.v1":
                raise ValueError(f"line {line_number}: unsupported receipt schema")
            by_family[receipt["source_family_id"]] = receipt
    return [by_family[key] for key in sorted(by_family)]


def process(receipts_path: Path, object_root: Path, derived_root: Path, vocabulary_path: Path, rebuild: bool = False) -> dict[str, Any]:
    vocabulary = load_vocabulary(vocabulary_path)
    receipts = successful_receipts(receipts_path)
    articles_dir = derived_root / "articles"
    articles_dir.mkdir(parents=True, exist_ok=True)
    manifests: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    totals = Counter()
    tag_totals = {dimension: Counter() for dimension in TAG_DIMENSIONS}
    skipped = 0

    for receipt in receipts:
        payload_hash = receipt["sha256"]
        local_id = digest_bytes(receipt["source_family_id"].encode("utf-8"))[:24]
        chunk_path = articles_dir / f"{local_id}.chunks.jsonl"
        manifest_path = articles_dir / f"{local_id}.manifest.json"
        existing = None
        if manifest_path.is_file() and not rebuild:
            try:
                existing = json.loads(manifest_path.read_text("utf-8"))
            except (OSError, json.JSONDecodeError):
                existing = None
        if existing and existing.get("payload_sha256") == payload_hash and existing.get("extractor_version") == EXTRACTOR_VERSION and chunk_path.is_file():
            if digest_bytes(chunk_path.read_bytes()) == existing.get("chunks_file_sha256"):
                manifests.append(existing)
                skipped += 1
                continue
        try:
            object_path = Path(receipt["object_path"])
            if not object_path.is_absolute():
                object_path = object_root / object_path
            payload = object_path.read_bytes()
            if digest_bytes(payload) != payload_hash:
                raise ValueError("payload SHA-256 mismatch")
            records, counts = extract_payload(payload, receipt, vocabulary)
            chunk_data = b"".join(canonical_bytes(record) + b"\n" for record in records)
            atomic_write(chunk_path, chunk_data)
            article_tag_counts = {dimension: sorted({tag for record in records for tag in record["tags"][dimension]}) for dimension in TAG_DIMENSIONS}
            manifest = {
                "schema": "aleph.external_training.jats_article_digest.v1",
                "authority_status": "proposed",
                "extractor_version": EXTRACTOR_VERSION,
                "source_family_id": receipt["source_family_id"],
                "pmcid": receipt["pmcid"],
                "payload_sha256": payload_hash,
                "payload_bytes": len(payload),
                "chunks": len(records),
                "chunks_file_sha256": digest_bytes(chunk_data),
                **counts,
                "tags_present": article_tag_counts,
            }
            atomic_write(manifest_path, canonical_bytes(manifest) + b"\n")
            manifests.append(manifest)
        except (OSError, ET.ParseError, ValueError, KeyError, TypeError) as exc:
            failures.append({"source_family_id": str(receipt.get("source_family_id", "UNKNOWN")), "failure_class": type(exc).__name__})

    for manifest in manifests:
        totals["payload_bytes"] += manifest["payload_bytes"]
        totals["chunks"] += manifest["chunks"]
        totals["sections"] += manifest["sections"]
        totals["figures"] += manifest["figures"]
        totals["tables"] += manifest["tables"]
        for dimension in TAG_DIMENSIONS:
            for tag in manifest["tags_present"][dimension]:
                tag_totals[dimension][tag] += 1
    digest_manifests = sorted(manifests, key=lambda item: item["source_family_id"])
    repeat_hash = digest_bytes(b"".join(canonical_bytes(item) + b"\n" for item in digest_manifests))
    return {
        "schema": "aleph.external_training.jats_extraction_summary.v1",
        "authority_status": "proposed",
        "extractor_version": EXTRACTOR_VERSION,
        "receipt_successes": len(receipts),
        "xml_parsed": len(manifests),
        "failures": len(failures),
        "failure_classes": dict(sorted(Counter(item["failure_class"] for item in failures).items())),
        "resumed_articles": skipped,
        **dict(totals),
        "tagged_articles": {dimension: sum(1 for item in manifests if item["tags_present"][dimension]) for dimension in TAG_DIMENSIONS},
        "tag_counts": {dimension: dict(sorted(counter.items())) for dimension, counter in tag_totals.items()},
        "article_manifest_set_sha256": repeat_hash,
        "failed_source_family_sha256": digest_bytes(canonical_bytes(sorted(item["source_family_id"] for item in failures))),
    }


def main() -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipts", type=Path, required=True)
    parser.add_argument("--object-root", type=Path, default=Path("."))
    parser.add_argument("--derived-root", type=Path, required=True)
    parser.add_argument("--vocabulary", type=Path, default=here / "tag_vocabulary.json")
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()
    summary = process(args.receipts, args.object_root, args.derived_root, args.vocabulary, args.rebuild)
    data = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.summary:
        atomic_write(args.summary, data.encode("utf-8"))
    print(data, end="")
    return 0 if summary["failures"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())

