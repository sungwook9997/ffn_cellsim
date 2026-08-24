#!/usr/bin/env python3
"""Acquire and validate IDR0006 confirmation annotations without pixels."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re
import time
import urllib.request

from train_hpa_raw_if_cnn import sha256


def normalized(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.casefold())


def download(url: str, path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(path.suffix + ".part")
    error: Exception | None = None
    for attempt in range(5):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "Project-Aleph-research/1.0"}
            )
            with urllib.request.urlopen(request, timeout=180) as response, part.open(
                "wb"
            ) as handle:
                while chunk := response.read(1024 * 1024):
                    handle.write(chunk)
            part.replace(path)
            return
        except Exception as exc:
            error = exc
            part.unlink(missing_ok=True)
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"failed IDR0006 metadata download: {error}")


def resolve_unique(headers: list[str], aliases: dict[str, list[str]]) -> tuple[dict[str, str], dict[str, list[str]]]:
    canonical = {header: normalized(header) for header in headers}
    resolved, ambiguous = {}, {}
    for semantic, values in aliases.items():
        accepted = {normalized(value) for value in values}
        matches = [header for header, value in canonical.items() if value in accepted]
        if len(matches) == 1:
            resolved[semantic] = matches[0]
        elif len(matches) > 1:
            ambiguous[semantic] = matches
    return resolved, ambiguous


def phenotype_class(value: str, mapping: dict[str, list[str]]) -> str | None:
    candidate = normalized(value)
    for canonical, values in mapping.items():
        if candidate in {normalized(item) for item in values}:
            return canonical
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    if not protocol["frozen_before_idr0006_annotation_csv_parquet_api_or_pixel_access"]:
        raise ValueError("IDR0006 protocol was not frozen before metadata access")
    source = protocol["official_source"]
    csv_path = args.cache_dir / "idr0006-screenA-annotation.csv"
    download(source["annotation_csv_url"], csv_path)
    lower, upper = source["annotation_csv_expected_repository_size_range_bytes"]
    if not lower <= csv_path.stat().st_size <= upper:
        raise ValueError("IDR0006 annotation CSV size outside frozen repository range")
    with csv_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        headers = reader.fieldnames or []
        rows = list(reader)
    gate = protocol["metadata_gate"]
    resolved, ambiguous = resolve_unique(
        headers, gate["required_unique_column_aliases"]
    )
    missing = sorted(
        set(gate["required_unique_column_aliases"]) - set(resolved) - set(ambiguous)
    )
    long_matches = [
        header for header in headers
        if normalized(header) in {
            normalized(value) for value in gate["accepted_phenotype_schema_A_long_column_aliases"]
        }
    ]
    wide_lookup = {
        normalized(value): value
        for value in gate["accepted_phenotype_schema_B_wide_exact_author_columns"]
    }
    wide_matches = [header for header in headers if normalized(header) in wide_lookup]
    if len(long_matches) == 1 and not wide_matches:
        phenotype_schema = "A_long"
    elif not long_matches and wide_matches:
        phenotype_schema = "B_wide"
    else:
        phenotype_schema = None

    ORF_classes: dict[str, set[str]] = defaultdict(set)
    ORF_genes: dict[str, set[str]] = defaultdict(set)
    ORF_rows: dict[str, int] = Counter()
    unmapped_phenotypes: Counter[str] = Counter()
    positive_tokens = {
        normalized(value) for value in gate["positive_value_tokens_for_wide_schema"]
    }
    if not missing and not ambiguous and phenotype_schema is not None:
        for row in rows:
            ORF = row[resolved["ORF_identifier"]].strip()
            gene = row[resolved["gene_symbol"]].strip()
            if not ORF:
                continue
            ORF_rows[ORF] += 1
            if gene:
                ORF_genes[ORF].add(gene)
            if phenotype_schema == "A_long":
                values = [
                    item.strip() for item in re.split(r"[;,|]", row[long_matches[0]])
                    if item.strip()
                ]
                for value in values:
                    mapped = phenotype_class(value, gate["phenotype_mapping"])
                    if mapped is None:
                        unmapped_phenotypes[value] += 1
                    else:
                        ORF_classes[ORF].add(mapped)
            else:
                for header in wide_matches:
                    if normalized(row[header].strip()) in positive_tokens:
                        mapped = phenotype_class(wide_lookup[normalized(header)], gate["phenotype_mapping"])
                        if mapped is None:
                            unmapped_phenotypes[wide_lookup[normalized(header)]] += 1
                        else:
                            ORF_classes[ORF].add(mapped)
    mixed = sorted(ORF for ORF, values in ORF_classes.items() if len(values) > 1)
    class_ORFs = {
        canonical: sorted(
            ORF for ORF, values in ORF_classes.items()
            if values == {canonical}
        )
        for canonical in gate["phenotype_mapping"]
    }
    schema_passed = not missing and not ambiguous and phenotype_schema is not None
    support_passed = (
        len(class_ORFs["nuclear_interior"]) >= gate["minimum_unique_nuclear_interior_ORFs"]
        and len(class_ORFs["nuclear_envelope"]) >= gate["minimum_unique_nuclear_envelope_ORFs"]
    )
    passed = schema_passed and support_passed
    report = {
        "schema": "aleph.outer_library.idr0006_confirmation_metadata_gate.v1",
        "protocol_sha256": sha256(args.protocol),
        "annotation_csv_url": source["annotation_csv_url"],
        "annotation_csv_local_name": csv_path.name,
        "annotation_csv_size_bytes": csv_path.stat().st_size,
        "annotation_csv_sha256": sha256(csv_path),
        "row_count": len(rows),
        "headers": headers,
        "resolved_required_columns": resolved,
        "ambiguous_required_columns": ambiguous,
        "missing_required_columns": missing,
        "long_phenotype_column_matches": long_matches,
        "wide_phenotype_column_matches": wide_matches,
        "selected_phenotype_schema": phenotype_schema,
        "unique_ORF_count": len(ORF_rows),
        "unique_gene_symbol_count": len({gene for values in ORF_genes.values() for gene in values}),
        "class_unique_ORF_counts": {
            key: len(value) for key, value in class_ORFs.items()
        },
        "mixed_interior_and_envelope_ORF_count": len(mixed),
        "mixed_interior_and_envelope_ORFs": mixed,
        "unmapped_phenotype_counts": dict(sorted(unmapped_phenotypes.items())),
        "schema_gate_passed": schema_passed,
        "class_support_gate_passed": support_passed,
        "metadata_gate_passed": passed,
        "status": (
            "metadata_passed_channel_gate_pending"
            if passed else "refused_before_IDR0006_pixel_access"
        ),
        "IDR0006_parquet_requests": 0,
        "IDR0006_channel_metadata_requests": 0,
        "IDR0006_image_objects_downloaded": 0,
        "IDR0006_pixels_accessed": 0,
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
        "rows": len(rows), "ORFs": len(ORF_rows),
        "phenotype_schema": phenotype_schema,
        "class_ORFs": report["class_unique_ORF_counts"],
        "metadata_gate_passed": passed, "IDR0006_pixels": 0,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
