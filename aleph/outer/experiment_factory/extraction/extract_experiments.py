#!/usr/bin/env python3
"""Deterministic, resumable, evidence-local ExperimentRecord candidate extractor.

Only explicit reported numbers with units are emitted.  Lexical context can label a
candidate but never supplies a missing value or relation.  Full records remain in an
ignored local data root; committed outputs are text-free digests and aggregates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Iterable

SCHEMA = "aleph.external_training.experiment_record.v1"
EXTRACTOR_VERSION = "1"
AUTHORITY_STATUS = "proposed"
NUMBER = r"[+-]?(?:\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?|\.\d+)(?:[eE][+-]?\d+)?"


@dataclass(frozen=True)
class UnitDef:
    pattern: str
    quantity: str
    canonical_unit: str
    factor: str
    offset: str = "0"


# Domain-canonical units preserve useful magnitudes.  Conversion is exact Decimal
# arithmetic over the reported decimal token; no physical value is inferred.
UNIT_DEFS = (
    UnitDef(r"pN[·*]?s/(?:µ|μ|u)m", "drag_coefficient", "pN*s/um", "1"),
    UnitDef(r"(?:µ|μ|u)m/(?:min)", "speed", "um/s", "0.01666666666666666666666666667"),
    UnitDef(r"(?:µ|μ|u)m/s", "speed", "um/s", "1"),
    UnitDef(r"nm/s", "speed", "um/s", "0.001"),
    UnitDef(r"kPa", "pressure_or_modulus", "Pa", "1000"),
    UnitDef(r"MPa", "pressure_or_modulus", "Pa", "1000000"),
    UnitDef(r"Pa", "pressure_or_modulus", "Pa", "1"),
    UnitDef(r"nN", "force", "pN", "1000"),
    UnitDef(r"pN", "force", "pN", "1"),
    UnitDef(r"mN", "force", "pN", "1000000000"),
    UnitDef(r"(?<![a-zA-Z])N", "force", "pN", "1000000000000"),
    UnitDef(r"nm", "length", "um", "0.001"),
    UnitDef(r"(?:µ|μ|u)m", "length", "um", "1"),
    UnitDef(r"mm", "length", "um", "1000"),
    UnitDef(r"cm", "length", "um", "10000"),
    UnitDef(r"ms", "time", "s", "0.001"),
    UnitDef(r"min(?:ute)?s?", "time", "s", "60"),
    UnitDef(r"h(?:ours?|r)?", "time", "s", "3600"),
    UnitDef(r"days?", "time", "s", "86400"),
    UnitDef(r"s(?:ec(?:ond)?s?)?", "time", "s", "1"),
    UnitDef(r"°C", "temperature", "K", "1", "273.15"),
    UnitDef(r"mM", "concentration", "mol/L", "0.001"),
    UnitDef(r"(?:µ|μ|u)M", "concentration", "mol/L", "0.000001"),
    UnitDef(r"nM", "concentration", "mol/L", "0.000000001"),
    UnitDef(r"ng/(?:mL|ml)", "mass_concentration", "g/L", "0.000001"),
    UnitDef(r"(?:µ|μ|u)g/(?:mL|ml)", "mass_concentration", "g/L", "0.001"),
    UnitDef(r"mg/(?:mL|ml)", "mass_concentration", "g/L", "1"),
    UnitDef(r"%", "percentage", "%", "1"),
    UnitDef(r"(?i:-?fold)", "fold_change", "fold", "1"),
)

UNIT_ALT = "|".join(f"(?:{item.pattern})" for item in UNIT_DEFS)
MEASURE_RE = re.compile(
    rf"(?P<value>{NUMBER})\s*(?:±\s*(?P<uncertainty>{NUMBER})\s*)?(?P<unit>{UNIT_ALT})(?![A-Za-z])",
)
N_RE = re.compile(r"(?<![A-Za-z])(?:n|N)\s*=\s*(\d{1,5})(?!\d)")
REPLICATE_RE = re.compile(r"(?P<n>\d{1,4})\s+(?P<kind>biological|technical|independent)\s+replicates?", re.I)

ASSAYS = {
    "traction_force_microscopy": ("traction force microscopy", "tfm"),
    "particle_image_velocimetry": ("particle image velocimetry", "piv"),
    "atomic_force_microscopy": ("atomic force microscopy", "afm"),
    "micropipette_aspiration": ("micropipette aspiration",),
    "optical_tweezers": ("optical tweezer", "optical trap"),
    "magnetic_tweezers": ("magnetic tweezer",),
    "immunofluorescence": ("immunofluorescence", "immunofluorescent", "confocal microscopy"),
    "western_blot": ("western blot", "immunoblot"),
    "pcr": ("polymerase chain reaction", "qpcr", "rt-pcr", "pcr"),
    "flow_cytometry": ("flow cytometry", "facs"),
    "live_cell_imaging": ("live-cell imaging", "live cell imaging", "time-lapse microscopy"),
    "rheology": ("rheology", "rheometer",),
    "indentation": ("indentation", "nanoindentation"),
}

CELL_TYPES = {
    "fibroblast": ("fibroblast",), "epithelial": ("epithelial", "keratinocyte"),
    "endothelial": ("endothelial",), "immune": ("macrophage", "lymphocyte", "neutrophil", "t cell", "b cell"),
    "stem_progenitor": ("stem cell", "ipsc", "progenitor"), "muscle": ("myocyte", "cardiomyocyte", "myoblast"),
    "neural_glial": ("neuron", "astrocyte", "microglia", "glial"), "cancer": ("carcinoma", "tumor cell", "cancer cell"),
    "organoid": ("organoid",), "red_blood_cell": ("erythrocyte", "red blood cell", "rbc"),
}

STATES = {
    "migrating": ("migrating", "migration"), "mitotic": ("mitotic", "mitosis"),
    "apoptotic": ("apoptotic", "apoptosis"), "differentiating": ("differentiating", "differentiation"),
    "hypoxic": ("hypoxic", "hypoxia"), "quiescent": ("quiescent", "quiescence"),
    "activated": ("activated", "activation"), "diseased": ("disease", "diabetic", "fibrotic"),
}

CONDITIONS = {
    "control": ("control", "vehicle", "untreated", "wild-type", "wild type"),
    "treatment": ("treated", "treatment", "inhibitor", "agonist", "drug"),
    "genetic": ("knockout", "knockdown", "silencing", "overexpression", "mutant", "crispr"),
    "substrate": ("substrate", "matrix", "hydrogel", "collagen", "fibronectin"),
    "environment": ("hypoxia", "temperature", "shear", "confinement", "osmotic", "stretch", "compression"),
}


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        handle.write(data)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def lname(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def text(element: ET.Element | None) -> str:
    return "" if element is None else " ".join("".join(element.itertext()).split())


def descendants(root: ET.Element, name: str) -> Iterable[ET.Element]:
    return (node for node in root.iter() if lname(node.tag) == name)


def first(root: ET.Element, name: str) -> ET.Element | None:
    return next(descendants(root, name), None)


def lexical_labels(value: str, vocabulary: dict[str, tuple[str, ...]]) -> list[str]:
    lowered = value.lower()
    return sorted(label for label, terms in vocabulary.items() if any(term in lowered for term in terms))


def unit_definition(unit: str) -> UnitDef | None:
    for definition in UNIT_DEFS:
        if re.fullmatch(definition.pattern, unit):
            return definition
    return None


def decimal_string(value: Decimal) -> str:
    result = format(value.normalize(), "f")
    return "0" if result in {"-0", ""} else result


def normalize_measure(value: str, unit: str) -> dict[str, Any]:
    definition = unit_definition(unit)
    if definition is None:
        return {"reported": unit, "normalized": None, "status": "unsupported"}
    try:
        number = Decimal(value.replace(",", ""))
        normalized = number * Decimal(definition.factor) + Decimal(definition.offset)
    except InvalidOperation:
        return {"reported": unit, "normalized": None, "status": "invalid_decimal"}
    return {
        "reported": unit,
        "normalized": definition.canonical_unit,
        "status": "exact_lexical_conversion",
        "normalized_value": decimal_string(normalized),
        "quantity": definition.quantity,
    }


def extract_observations(evidence_text: str, evidence: dict[str, Any], observation_role: str) -> list[dict[str, Any]]:
    observations: list[dict[str, Any]] = []
    sample_match = N_RE.search(evidence_text)
    replicate_match = REPLICATE_RE.search(evidence_text)
    replicate_kind = replicate_match.group("kind").lower() if replicate_match else None
    replicate_n = int(replicate_match.group("n")) if replicate_match else None
    replicate = {
        "biological_n": replicate_n if replicate_kind == "biological" else None,
        "technical_n": replicate_n if replicate_kind == "technical" else None,
        "unit_of_replication": f"reported_n={sample_match.group(1)};type_unspecified" if sample_match and not replicate_match else None,
        "aggregation": "unknown",
        "independence_status": "reported" if replicate_kind in {"biological", "technical"} else ("ambiguous" if sample_match or replicate_match else "unverified"),
    }
    for ordinal, match in enumerate(MEASURE_RE.finditer(evidence_text)):
        raw_value = match.group("value").replace(",", "")
        unit = match.group("unit")
        # SI seconds are lowercase.  Uppercase S occurs frequently in gene names
        # (notably 18 S rRNA) and is not accepted as a time unit.
        if unit == "S":
            continue
        if unit == "N":
            context = evidence_text[max(0, match.start() - 60):match.end() + 60]
            after_unit = evidence_text[match.end():match.end() + 5]
            if (not re.search(r"\b(?:force|load|newton|tension|thrust|traction)\b", context, re.I)
                    or re.match(r"\s*(?:/|m\s*[−-]?1)", after_unit)):
                continue
        normalized = normalize_measure(raw_value, unit)
        following = evidence_text[match.end():match.end() + 40]
        ambiguous_micro_molar = (
            unit in {"µm", "μm", "um"}
            and bool(re.match(r"\s+(?:[A-Z]{2,}[A-Z0-9-]*|[A-Za-z]+-\d+)", following))
        )
        if ambiguous_micro_molar:
            normalized = {
                "reported": unit, "normalized": None, "status": "ambiguous_case_or_typesetting",
                "normalized_value": None, "quantity": "length_or_concentration",
            }
        uncertainty: list[dict[str, Any]] = []
        if match.group("uncertainty"):
            uncertainty.append({
                "type": "unknown",
                "lower": None,
                "upper": None,
                "value": float(match.group("uncertainty").replace(",", "")),
                "unit": unit,
                "level": None,
            })
        converted = normalized.get("normalized_value")
        is_converted = normalized.get("normalized") != unit or normalized.get("normalized_value") != raw_value
        unit_record = {
            "reported": unit,
            "normalized": normalized.get("normalized"),
            "normalization_status": "unknown" if ambiguous_micro_molar else ("converted" if is_converted else "exact"),
            "conversion_formula": (
                f"normalized = reported * {unit_definition(unit).factor} + {unit_definition(unit).offset}"
                if is_converted and not ambiguous_micro_molar and unit_definition(unit) is not None else None
            ),
            "conversion_evidence": [evidence] if is_converted and not ambiguous_micro_molar else [],
        }
        obs_basis = f"{evidence['document_sha256']}|{evidence['locator_id']}|{evidence['text_or_asset_sha256']}|{match.start()}|{match.end()}|{ordinal}"
        observations.append({
            "observation_id": f"aleph:observation:{sha(obs_basis.encode())[:16]}",
            "observation_type": "scalar",
            "protocol_id": "protocol:pending",
            "condition_ids": ["condition:pending"],
            "measured_entity": f"{observation_role}:{normalized.get('quantity')}",
            "value_state": "extracted",
            "value": {"reported_value": raw_value, "normalized_value": converted},
            "unit": unit_record,
            "uncertainty": uncertainty,
            "replicates": replicate,
            "visual_payload": None,
            "evidence": [{**evidence, "char_start": match.start(), "char_end": match.end()}],
            "extraction": {
                "status": "ambiguous" if uncertainty or ambiguous_micro_molar else "machine_candidate",
                "confidence": 0.5 if ambiguous_micro_molar else (0.9 if uncertainty else 1.0),
                "method": "explicit_number_and_unit_v1",
                "ambiguity_reason": ("micrometre token followed by compound-like name; lowercase typesetting may mean micromolar"
                                     if ambiguous_micro_molar else ("plus/minus semantics not reported in evidence unit" if uncertainty else None)),
                "candidates": (["length", "concentration"] if ambiguous_micro_molar else (["sd", "se", "unknown"] if uncertainty else [])),
            },
        })
    return observations


def acquisition_receipts(paths: list[tuple[str, Path]]) -> list[dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    for acquisition_pass, path in paths:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                row = json.loads(line)
                if row.get("status") == "success":
                    row = dict(row)
                    row["receipt_sha256"] = sha(canonical_bytes(row))
                    row["acquisition_pass"] = acquisition_pass
                    rows[row["source_family_id"].lower()] = row
    return [rows[key] for key in sorted(rows)]


def publication_year(article: ET.Element) -> int | None:
    years = []
    for date in descendants(article, "pub-date"):
        year = first(date, "year")
        if year is not None and text(year).isdigit():
            years.append(int(text(year)))
    return years[0] if years else None


def section_units(article: ET.Element) -> Iterable[tuple[list[str], str, str, str]]:
    """Yield (path, locator_type, locator_id, normalized evidence text)."""
    abstract = first(article, "abstract")
    if abstract is not None and text(abstract):
        yield ["Abstract"], "section", "abstract:0", text(abstract)

    def walk(section: ET.Element, parent: list[str], ordinal: list[int]) -> Iterable[tuple[list[str], str, str, str]]:
        ordinal[0] += 1
        title_node = next((child for child in section if lname(child.tag) == "title"), None)
        title = text(title_node) or "Untitled section"
        path = parent + [title]
        locator = section.attrib.get("id") or str(ordinal[0])
        paragraphs = [child for child in section if lname(child.tag) in {"p", "list"}]
        for index, paragraph in enumerate(paragraphs):
            value = text(paragraph)
            if value:
                yield path, "section", f"section:{locator}/paragraph:{index}", value
        for child in section:
            if lname(child.tag) == "sec":
                yield from walk(child, path, ordinal)

    body = first(article, "body")
    if body is not None:
        ordinal = [0]
        parent = {child: node for node in body.iter() for child in node}
        for section in (child for child in body.iter() if lname(child.tag) == "sec"):
            ancestor = parent.get(section)
            has_section_ancestor = False
            while ancestor is not None and ancestor is not body:
                if lname(ancestor.tag) == "sec":
                    has_section_ancestor = True
                    break
                ancestor = parent.get(ancestor)
            if not has_section_ancestor:
                yield from walk(section, [], ordinal)

    for table_index, table in enumerate(descendants(article, "table-wrap")):
        table_id = table.attrib.get("id") or str(table_index)
        caption = text(first(table, "caption"))
        rows = list(descendants(table, "tr"))
        headers: list[str] = []
        for row_index, row in enumerate(rows):
            cells = [text(cell) for cell in row if lname(cell.tag) in {"td", "th"}]
            if row_index == 0:
                headers = cells
            row_context = " | ".join(part for part in [caption, *headers, *cells] if part)
            if row_context:
                yield ["Table", caption[:120]], "table", f"table:{table_id}/row:{row_index}", row_context

    for figure_index, figure in enumerate(descendants(article, "fig")):
        caption = text(first(figure, "caption"))
        if caption:
            figure_id = figure.attrib.get("id") or str(figure_index)
            yield ["Figure"], "figure", f"figure:{figure_id}/caption", caption


def explicit_comparison(value: str, evidence: dict[str, Any]) -> list[dict[str, Any]]:
    relation = None
    if re.search(r"\b(?:versus|vs\.?|compared with|compared to)\b", value, re.I):
        relation = "explicit_comparison_mentioned"
    elif re.search(r"\b(?:increased|decreased|higher|lower)\s+than\b", value, re.I):
        relation = "explicit_directional_comparison"
    # The schema requires resolved arm references.  A lexical comparison without
    # resolved arms is therefore not emitted as a comparison object.
    return []


def ontology(label: str | None, status: str = "candidate") -> dict[str, Any]:
    return {"label": label, "identifier": None, "normalization_status": status, "candidates": []}


def condition_factor_type(label: str) -> str:
    return {
        "control": "control", "treatment": "pharmacological", "genetic": "genetic",
        "substrate": "substrate_ecm", "environment": "other",
    }.get(label, "unknown")


def modality(label: str) -> str:
    return {
        "atomic_force_microscopy": "afm", "particle_image_velocimetry": "piv",
        "rheology": "other", "indentation": "other",
    }.get(label, label)


def record_for_unit(receipt: dict[str, Any], year: int | None, section_path: list[str], locator_type: str,
                    locator_id: str, evidence_text: str) -> dict[str, Any] | None:
    assays = lexical_labels(evidence_text, ASSAYS)
    cells = lexical_labels(evidence_text, CELL_TYPES)
    states = lexical_labels(evidence_text, STATES)
    conditions = lexical_labels(evidence_text, CONDITIONS)
    # Restrict ordinary prose to primary experiment-bearing evidence classes.
    # Background/introduction/discussion numbers often restate other studies and
    # are not automatic ExperimentRecords.  Abstracts require explicit study verbs.
    path_text = " ".join(section_path).lower()
    if locator_type == "section":
        if locator_id == "abstract:0":
            if not re.search(r"\b(?:we|this study|our (?:results|experiments?))\b.*\b(?:measured|tested|examined|investigated|quantified|assessed|found|showed|demonstrated)\b", evidence_text, re.I):
                return None
        elif not any(term in path_text for term in ("method", "material", "protocol", "experimental", "result", "finding")):
            return None
    # An explicit biological or experimental context is required. This rejects e.g.
    # non-biological PIV papers unless another biological signal is present.
    if not cells and not assays and not states:
        return None
    evidence_hash = sha(evidence_text.encode("utf-8"))
    schema_locator_type = "paragraph" if locator_type == "section" else locator_type
    evidence = {
        "document_sha256": receipt["sha256"],
        "section_path": section_path,
        "locator_type": schema_locator_type,
        "locator_id": locator_id,
        "text_or_asset_sha256": evidence_hash,
        "char_start": None,
        "char_end": None,
    }
    is_methods = any(term in path_text for term in ("method", "material", "protocol", "experimental"))
    observation_role = "protocol_parameter" if is_methods else "reported_measurement_candidate"
    observations = extract_observations(evidence_text, evidence, observation_role)
    if not observations:
        return None
    record_basis = f"{receipt['source_family_id'].lower()}|{locator_type}|{locator_id}|{evidence_hash}"
    ambiguity = []
    if not cells:
        ambiguity.append("biological_system_not_explicit_in_evidence_unit")
    if not assays:
        ambiguity.append("protocol_not_explicit_in_evidence_unit")
    if any(item["uncertainty"] for item in observations):
        ambiguity.append("plus_minus_semantics_unspecified")
    if any(item["extraction"]["ambiguity_reason"] and "micrometre token" in item["extraction"]["ambiguity_reason"] for item in observations):
        ambiguity.append("micrometre_or_micromolar_typesetting_ambiguous")
    record_id = f"aleph:experiment:{sha(record_basis.encode())[:16]}"
    system_id = f"system:{record_id.rsplit(':', 1)[-1]}"
    condition_id = f"condition:{record_id.rsplit(':', 1)[-1]}"
    protocol_labels = assays or ["unknown"]
    protocol_ids = [f"protocol:{record_id.rsplit(':', 1)[-1]}:{index}" for index in range(len(protocol_labels))]
    for observation in observations:
        observation["protocol_id"] = protocol_ids[0]
        observation["condition_ids"] = [condition_id]
    factors = [{
        "factor_type": condition_factor_type(condition), "name": condition, "value": None,
        "unit": None, "duration": None, "status": "extracted",
    } for condition in conditions]
    if not factors:
        factors = [{"factor_type": "unknown", "name": None, "value": None, "unit": None, "duration": None, "status": "missing"}]
    system_cell = cells[0] if cells else None
    return {
        "schema": SCHEMA,
        "authority_status": AUTHORITY_STATUS,
        "review_state": "gold_candidate",
        "experiment_record_id": record_id,
        "source": {
            "source_family_id": receipt["source_family_id"].lower(),
            "pmcid": receipt.get("pmcid"),
            "acquisition_pass": receipt["acquisition_pass"],
            "publication_year": year,
            "primary_domain": "unmapped",
        },
        "biological_systems": [{
            "biological_system_id": system_id,
            "organism": ontology(None, "missing"), "tissue": ontology(None, "missing"),
            "cell_type": ontology(system_cell, "candidate" if system_cell else "missing"),
            "cell_line": ontology(None, "missing"),
            "cell_states": [ontology(state) for state in states],
            "genotype": None,
            "donor_attributes": {"sex": None, "age": None, "disease": None, "passage": None},
        }],
        "conditions": [{
            "condition_id": condition_id, "biological_system_id": system_id,
            "role": "control" if "control" in conditions else ("treatment" if any(c in conditions for c in ("treatment", "genetic")) else "unknown"),
            "factors": factors,
            "environment": {"media": None, "temperature": None, "substrate": None, "geometry": None},
        }],
        "protocols": [{
            "protocol_id": pid, "modality": modality(label), "instrument": None,
            "calibration": None, "analysis_pipeline": None, "evidence": [evidence],
        } for pid, label in zip(protocol_ids, protocol_labels)],
        "observations": observations,
        "comparisons": explicit_comparison(evidence_text, evidence),
        "provenance": {
            "source_receipt_sha256": receipt["receipt_sha256"],
            "candidate_record_sha256": evidence_hash,
            "annotation_template_version": "1.0.0",
            "annotator_ids": [], "adjudicator_id": None,
            "created_at": "2026-08-05T00:00:00Z", "updated_at": "2026-08-05T00:00:00Z",
            "record_sha256": None,
        },
        # Local-only bounded evidence makes human review possible. It is excluded
        # from compact committed indices and never promoted to authority.
        "_local_evidence_text": evidence_text,
        "_local_ambiguity_flags": ambiguity,
    }


def extract_article(payload: bytes, receipt: dict[str, Any]) -> list[dict[str, Any]]:
    root = ET.fromstring(payload)
    article = root if lname(root.tag) == "article" else first(root, "article")
    if article is None:
        raise ValueError("no article element")
    year = publication_year(article)
    records = []
    for args in section_units(article):
        candidate = record_for_unit(receipt, year, *args)
        if candidate is not None:
            records.append(candidate)
    records.sort(key=lambda item: item["experiment_record_id"])
    return records


def text_free_record(record: dict[str, Any]) -> dict[str, Any]:
    copy = dict(record)
    copy.pop("_local_evidence_text", None)
    copy.pop("_local_ambiguity_flags", None)
    return copy


def load_domain_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    result = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            result[row["source_family_id"].lower()] = row["label"]
    return result


def process(receipt_inputs: list[tuple[str, Path]], object_root: Path, output_root: Path,
            rebuild: bool = False, limit: int | None = None, domain_map_path: Path | None = None) -> dict[str, Any]:
    receipts = acquisition_receipts(receipt_inputs)
    domain_map = load_domain_map(domain_map_path)
    if limit is not None:
        receipts = receipts[:limit]
    records_dir = output_root / "records"
    evidence_dir = output_root / "evidence"
    manifests_dir = output_root / "manifests"
    records_dir.mkdir(parents=True, exist_ok=True)
    evidence_dir.mkdir(parents=True, exist_ok=True)
    manifests_dir.mkdir(parents=True, exist_ok=True)
    manifests = []
    failures = []
    resumed = 0
    for receipt in receipts:
        local_id = sha(receipt["source_family_id"].lower().encode())[:24]
        record_path = records_dir / f"{local_id}.jsonl"
        manifest_path = manifests_dir / f"{local_id}.json"
        if manifest_path.is_file() and record_path.is_file() and not rebuild:
            existing = json.loads(manifest_path.read_text("utf-8"))
            if (existing.get("payload_sha256") == receipt["sha256"] and
                    existing.get("extractor_version") == EXTRACTOR_VERSION and
                    sha(record_path.read_bytes()) == existing.get("records_file_sha256")):
                manifests.append(existing)
                resumed += 1
                continue
        try:
            object_path = Path(receipt["object_path"])
            if not object_path.is_absolute():
                object_path = object_root / object_path
            payload = object_path.read_bytes()
            if sha(payload) != receipt["sha256"]:
                raise ValueError("payload SHA-256 mismatch")
            records = extract_article(payload, receipt)
            domain = domain_map.get(receipt["source_family_id"].lower(), "unmapped")
            for record in records:
                record["source"]["primary_domain"] = domain
            schema_records = [text_free_record(record) for record in records]
            record_data = b"".join(canonical_bytes(record) + b"\n" for record in schema_records)
            evidence_data = b"".join(canonical_bytes({
                "experiment_record_id": record["experiment_record_id"],
                "text_sha256": sha(record["_local_evidence_text"].encode()),
                "text": record["_local_evidence_text"],
                "ambiguity_flags": record["_local_ambiguity_flags"],
            }) + b"\n" for record in records)
            atomic_write(record_path, record_data)
            evidence_path = evidence_dir / f"{local_id}.jsonl"
            atomic_write(evidence_path, evidence_data)
            assays = sorted({p["modality"] for r in records for p in r["protocols"] if p["modality"] != "unknown"})
            cell_types = sorted({b["cell_type"]["label"] for r in records for b in r["biological_systems"] if b["cell_type"]["label"]})
            observations = [o for r in records for o in r["observations"]]
            role_counts = Counter(o["measured_entity"].split(":", 1)[0] for o in observations)
            quantity_counts = Counter(o["measured_entity"].split(":", 1)[-1] for o in observations)
            manifest = {
                "schema": "aleph.external_training.experiment_article_manifest.v1",
                "authority_status": AUTHORITY_STATUS,
                "extractor_version": EXTRACTOR_VERSION,
                "source_family_id": receipt["source_family_id"].lower(),
                "pmcid": receipt.get("pmcid"),
                "acquisition_pass": receipt["acquisition_pass"],
                "domain": domain,
                "payload_sha256": receipt["sha256"],
                "records": len(records),
                "observations": len(observations),
                "normalized_observations": sum(o["unit"]["normalization_status"] in {"exact", "converted"} for o in observations),
                "ambiguous_records": sum(bool(r["_local_ambiguity_flags"]) for r in records),
                "table_records": sum(r["protocols"][0]["evidence"][0]["locator_type"] == "table" for r in records),
                "figure_records": sum(r["protocols"][0]["evidence"][0]["locator_type"] == "figure" for r in records),
                "records_with_comparison": sum(bool(r["comparisons"]) for r in records),
                "observations_with_uncertainty": sum(bool(o["uncertainty"]) for o in observations),
                "observations_with_sample_size_or_replicates": sum(o["replicates"]["independence_status"] in {"reported", "ambiguous"} for o in observations),
                "observation_role_counts": dict(sorted(role_counts.items())),
                "quantity_counts": dict(sorted(quantity_counts.items())),
                "assays": assays,
                "cell_types": cell_types,
                "records_file_sha256": sha(record_data),
                "evidence_file_sha256": sha(evidence_data),
                "record_set_sha256": sha(record_data),
            }
            atomic_write(manifest_path, canonical_bytes(manifest) + b"\n")
            manifests.append(manifest)
        except (OSError, ValueError, KeyError, ET.ParseError, json.JSONDecodeError) as exc:
            failures.append({"source_family_id": receipt.get("source_family_id"), "failure_class": type(exc).__name__})

    assay_counts: Counter[str] = Counter()
    cell_counts: Counter[str] = Counter()
    pass_counts: dict[str, Counter[str]] = defaultdict(Counter)
    domain_counts: dict[str, Counter[str]] = defaultdict(Counter)
    role_counts: Counter[str] = Counter()
    quantity_counts: Counter[str] = Counter()
    totals: Counter[str] = Counter()
    for item in manifests:
        for key in ("records", "observations", "normalized_observations", "ambiguous_records", "table_records", "figure_records",
                    "records_with_comparison", "observations_with_uncertainty", "observations_with_sample_size_or_replicates"):
            totals[key] += item[key]
            pass_counts[item["acquisition_pass"]][key] += item[key]
            domain_counts[item["domain"]][key] += item[key]
        if item["records"]:
            totals["articles_with_records"] += 1
            pass_counts[item["acquisition_pass"]]["articles_with_records"] += 1
            domain_counts[item["domain"]]["articles_with_records"] += 1
        domain_counts[item["domain"]]["input_articles"] += 1
        role_counts.update(item["observation_role_counts"])
        quantity_counts.update(item["quantity_counts"])
        assay_counts.update(item["assays"])
        cell_counts.update(item["cell_types"])
    ordered = sorted(manifests, key=lambda item: item["source_family_id"])
    summary = {
        "schema": "aleph.external_training.experiment_extraction_summary.v1",
        "authority_status": AUTHORITY_STATUS,
        "extractor_version": EXTRACTOR_VERSION,
        "input_articles": len(receipts),
        "processed_articles": len(manifests),
        "resumed_articles": resumed,
        "failures": len(failures),
        "failure_classes": dict(sorted(Counter(f["failure_class"] for f in failures).items())),
        **dict(totals),
        "unit_normalization_success_rate": (totals["normalized_observations"] / totals["observations"] if totals["observations"] else 0.0),
        "record_ambiguity_rate": (totals["ambiguous_records"] / totals["records"] if totals["records"] else 0.0),
        "article_coverage_rate": (totals["articles_with_records"] / len(receipts) if receipts else 0.0),
        "by_acquisition_pass": {key: dict(counter) for key, counter in sorted(pass_counts.items())},
        "by_domain": {key: dict(counter) for key, counter in sorted(domain_counts.items())},
        "observation_role_counts": dict(sorted(role_counts.items())),
        "quantity_counts": dict(sorted(quantity_counts.items())),
        "assay_article_counts": dict(sorted(assay_counts.items())),
        "cell_type_article_counts": dict(sorted(cell_counts.items())),
        "manifest_set_sha256": sha(b"".join(canonical_bytes(item) + b"\n" for item in ordered)),
        "failed_family_set_sha256": sha(canonical_bytes(sorted(f["source_family_id"] for f in failures))),
    }
    atomic_write(output_root / "summary.json", json.dumps(summary, indent=2, sort_keys=True).encode() + b"\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", action="append", nargs=2, metavar=("PASS", "PATH"), required=True)
    parser.add_argument("--object-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--domain-map", type=Path)
    args = parser.parse_args()
    summary = process([(name, Path(path)) for name, path in args.receipt], args.object_root, args.output_root,
                      args.rebuild, args.limit, args.domain_map)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["failures"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
