import json
from pathlib import Path

import numpy as np

from aleph.outer.experiment_factory.visual.visual_pipeline import (
    descriptor,
    figure_records,
    html_asset_map,
    modality_tags,
    normalized_figure_id,
    panels,
    resolve_asset_url,
    sha256_bytes,
    stable_id,
)


def test_ids_are_stable_and_domain_separated():
    assert stable_id("experiment-candidate", "doi:x", "fig[F1]", 0) == stable_id("experiment-candidate", "doi:x", "fig[F1]", 0)
    assert stable_id("experiment-candidate", "doi:x", "fig[F1]", 0) != stable_id("observation-candidate", "doi:x", "fig[F1]", 0)


def test_modality_rules_are_proposed_and_do_not_infer_values():
    tags = modality_tags("Confocal immunofluorescence and a TFM traction map over time")
    names = {tag["tag"] for tag in tags}
    assert {"immunofluorescence", "fluorescence_microscopy", "tfm_map"} <= names
    assert all(tag["authority_status"] == "proposed" for tag in tags)


def test_html_mapping_accepts_only_official_ncbi_cdn_assets():
    page = b'''<img src="https://cdn.ncbi.nlm.nih.gov/pmc/blobs/aa/1/hash/fig1.jpg">
    <img src="https://publisher.example/fig2.jpg">'''
    assert html_asset_map(page) == {"fig1.jpg": "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/aa/1/hash/fig1.jpg"}


def test_asset_resolution_uses_exact_then_same_stem_image_fallback():
    mapping = {"fig1.jpg": "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/a/fig1.jpg"}
    assert resolve_asset_url(mapping, "fig1.jpg") == mapping["fig1.jpg"]
    assert resolve_asset_url(mapping, "fig1.gif") == mapping["fig1.jpg"]
    assert resolve_asset_url(mapping, "publisher-only.gif") is None


def test_panel_candidates_split_on_bright_internal_gutter():
    image = np.full((128, 128), 30, dtype=np.uint8)
    image[:, 62:66] = 255
    boxes = panels(image)
    assert len(boxes) == 2
    assert boxes[0]["x1"] <= boxes[1]["x0"]


def test_no_panel_claim_for_uniform_or_border_only_image():
    assert panels(np.full((128, 128), 255, dtype=np.uint8)) == []
    image = np.full((128, 128), 20, dtype=np.uint8)
    image[:, :6] = 255
    assert panels(image) == []


def test_sha256_is_exact_bytes_not_filename(tmp_path: Path):
    first = b"same filename, first payload"
    second = b"same filename, second payload"
    assert sha256_bytes(first) != sha256_bytes(second)


def test_missing_graphic_is_typed_and_caption_text_is_not_stored(tmp_path: Path):
    xml = tmp_path / "article.xml"
    xml.write_text("<article><body><fig id='F1'><label>Figure 1</label><caption><p>Western blot</p></caption></fig></body></article>")
    receipt = {"object_path": str(xml), "source_family_id": "doi:X", "pmcid": "PMC1", "sha256": "a" * 64}
    record = next(iter(figure_records(receipt, "fixture")))
    assert record["graphics"] == []
    assert "Western blot" not in json.dumps(record)
    assert record["scientific_value_state"] == "refused_without_calibration"


def test_descriptor_replay_is_exact_and_refuses_values(tmp_path: Path):
    ppm = tmp_path / "fixture.ppm"
    header = b"P6\n2 2\n255\n"
    ppm.write_bytes(header + bytes([0, 0, 0, 255, 0, 0, 0, 255, 0, 255, 255, 255]))
    first = descriptor(ppm)
    second = descriptor(ppm)
    assert first == second
    assert first["width"] == 2 and first["height"] == 2
    assert first["scientific_value_state"] == "refused_without_calibration"


def test_figure_locator_normalization_is_exact_not_substring():
    assert normalized_figure_id("figure:Fig-2A/caption") == "fig-2a"
    assert normalized_figure_id("section:Fig-2A/caption") is None
    assert normalized_figure_id("figure:Fig-2A/body") is None
