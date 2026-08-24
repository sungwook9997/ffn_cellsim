from __future__ import annotations

import gzip
import json
from pathlib import Path

from .build_projection import build, digest, split_for_group


def write_jsonl(path: Path, values: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(json.dumps(value) + "\n" for value in values))


def fixture(tmp_path: Path) -> dict[str, Path]:
    article = tmp_path / "articles.jsonl"
    write_jsonl(article, [{
        "source_family_id": f"doi:s{i}", "payload_sha256": ("a" if i < 2 else "b") * 64,
        "record_set_sha256": str(i) * 64, "domain": "mechanics", "assays": ["afm"],
        "cell_types": ["fibroblast"], "acquisition_pass": "first",
    } for i in range(3)])
    gold = tmp_path / "gold.jsonl"
    write_jsonl(gold, [{"source_family_id": "doi:s0"}])
    registry = tmp_path / "registry.jsonl"
    write_jsonl(registry, [
        {"source_family_id": "doi:s1", "provider": "GEO", "accession": "GSE1", "leakage_component_id": "ds:x", "modality_hints": ["omics"]},
        {"source_family_id": "doi:s2", "provider": "GEO", "accession": "GSE1", "leakage_component_id": "ds:x", "modality_hints": ["omics"]},
    ])
    manifests = tmp_path / "manifests.jsonl"
    write_jsonl(manifests, [{"leakage_component_id": "ds:x", "manifest_id": "dataset:x"}])
    records = tmp_path / "records"
    write_jsonl(records / "one.jsonl", [{
        "source": {"source_family_id": "doi:s0"}, "experiment_record_id": "aleph:experiment:0000000000000000",
        "biological_systems": [{"cell_type": {"normalization_status": "candidate"}}],
        "protocols": [{"modality": "afm"}],
        "observations": [
            {"observation_id": "aleph:observation:0000000000000000", "measured_entity": "reported_measurement_candidate:force",
             "value": {"reported_value": "1", "normalized_value": "1"}, "unit": {"reported": "pN"},
             "extraction": {"status": "machine_candidate"}, "evidence": [{"text_or_asset_sha256": "c" * 64}]},
            {"observation_id": "aleph:observation:1111111111111111", "measured_entity": "reported_measurement_candidate:length_or_concentration",
             "value": {"reported_value": "2", "normalized_value": None}, "unit": {"reported": "um"},
             "extraction": {"status": "ambiguous"}, "evidence": [{"text_or_asset_sha256": "d" * 64}]},
        ],
    }])
    visual = tmp_path / "visual.jsonl.gz"
    with gzip.open(visual, "wt") as handle:
        handle.write(json.dumps({"source_family_id": "doi:s0", "observation_candidate_id": "v:1", "proposed_modality_tags": [{"tag": "afm"}]}) + "\n")
    links = tmp_path / "links.jsonl.gz"
    with gzip.open(links, "wt") as handle:
        handle.write(json.dumps({"source_family_id": "doi:s0", "linked_experiment_record_ids": ["e:1"]}) + "\n")
    assets = tmp_path / "assets.jsonl"
    write_jsonl(assets, [{"source_family_id": "doi:s0", "sha256": "e" * 64}])
    descriptors = tmp_path / "descriptors.jsonl"
    write_jsonl(descriptors, [
        {"source_family_id": "doi:s0", "descriptor_state": "described", "descriptor_id": "vd:1",
         "asset_sha256": "e" * 64, "observation_candidate_id": "v:1",
         "descriptor": {"width": 100, "height": 50, "entropy_bits": 3.5, "panel_candidates": [{}, {}], "source_pixel_format": "rgb24"}},
        {"source_family_id": "doi:s1", "descriptor_state": "not_acquired", "descriptor_id": "vd:2", "descriptor": None},
    ])
    return {"article_index": article, "gold_candidates": gold, "records_root": records,
            "visual_index": visual, "visual_links": links, "asset_receipts": assets,
            "visual_descriptors": descriptors,
            "dataset_registry": registry, "dataset_manifests": manifests}


def test_split_is_stable() -> None:
    assert split_for_group("lg:x") == split_for_group("lg:x")


def test_projection_joins_and_masks_without_leakage(tmp_path: Path) -> None:
    args = fixture(tmp_path)
    output, summary, local = tmp_path / "out.jsonl", tmp_path / "summary.json", tmp_path / "numeric.jsonl"
    result = build(output=output, summary_path=summary, local_numeric_output=local, **args)
    values = [json.loads(line) for line in output.read_text().splitlines()]
    assert result["sources"] == 3
    # payload and accession edges connect all three sources.
    assert result["leakage_groups"] == 1
    assert result["leakage_split_violations"] == 0
    assert len({row["split"] for row in values}) == 1
    assert values[0]["evaluation_cohort"] == "frozen_300_gold_candidate_not_ground_truth"
    numeric = [json.loads(line) for line in local.read_text().splitlines()]
    assert [row["target_mask"] for row in numeric] == [True, False]
    assert result["ambiguous_observations_promoted_to_labels"] == 0
    assert values[0]["views"]["visual_descriptor"]["numeric_aggregates"]["width"]["mean"] == 100.0
    assert values[0]["views"]["visual_metadata"]["present"] is True
    assert values[0]["views"]["visual_asset_receipt"]["present"] is True
    assert values[1]["missing_modality_mask"]["visual_descriptor"] is True
    assert '"text":' not in output.read_text().lower()


def test_optional_inputs_can_be_missing_and_replay_is_identical(tmp_path: Path) -> None:
    args = fixture(tmp_path)
    output, summary = tmp_path / "out.jsonl", tmp_path / "summary.json"
    first = build(article_index=args["article_index"], output=output, summary_path=summary)
    data = output.read_bytes()
    second = build(article_index=args["article_index"], output=output, summary_path=summary)
    assert data == output.read_bytes()
    assert first["projection_sha256"] == second["projection_sha256"]
    assert "visual_index" in second["missing_optional_inputs"]
