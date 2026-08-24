from pathlib import Path
import gzip
import hashlib
import io
import json
import sys
import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from multiprovider_if_common import map_author_label, normalize_crop_resize

from build_tfm_manifest import build
from build_nuclear_manifest import build as build_nuclear
from build_mechanotransduction_manifest import build as build_mechanotransduction
from sweep_gate import evaluate
from external_inference import infer as infer_external
from evaluate_external_inference import evaluate as evaluate_external_inference
from resolve_idr0072_if_flex_sources_v2 import expected_filename as idr0072_expected_filename
from latent_promotion import evaluate as evaluate_latent_promotion
from audit_observation_operators import audit as audit_observation_operators
from audit_sciplex3_design import audit as audit_sciplex3_design
from build_sciplex3_state_tensor import _triplet_chunks, _treatment_splits
from train_sciplex3_state_model import train as train_sciplex3_state_model
from download_sciplex3_matrix import _verify_gzip_stream
from download_sciplex3_metadata import _range as download_resumable_range
from audit_gse250041_senescence_design import audit as audit_gse250041_senescence_design
from download_gse250041_matrices import _matrix_market_header
from build_gse250041_senescence_tensor import (
    RNA_MARKERS as GSE250041_RNA_MARKERS,
    build as build_gse250041_senescence_tensor,
)
from verify_senscout_archive import _safe_member as senscout_safe_member
from train_senscout_morphology_model import (
    _partition_groups as senscout_partition_groups,
    _transform as senscout_transform,
)
from extract_figshare_cell_death_aggregates import _safe as figshare_archive_path_safe
from evaluate_gse301164_confirmation import _load_counts as load_gse301164_counts
from download_zenodo_tfm_figure_tables import _safe as zenodo_tfm_archive_path_safe
from evaluate_zenodo_tfm_tension import evaluate as evaluate_zenodo_tfm_tension
from download_zenodo_tfm_fret_source import _safe as zenodo_tfm_fret_path_safe
from build_zenodo_tfm_fret_manifest import build as build_zenodo_tfm_fret
from evaluate_zenodo_tfm_fret import evaluate as evaluate_zenodo_tfm_fret
from build_figshare_exif_emt_manifest import build as build_figshare_exif_emt
from evaluate_figshare_exif_emt import evaluate as evaluate_figshare_exif_emt
from download_figshare_exif_emt_subset import _safe as figshare_exif_path_safe
from build_mtrack_emt_manifest import (
    build as build_mtrack_emt,
    decode_gray16_png as decode_mtrack_mask,
    decode_uncompressed_tiff as decode_mtrack_tiff,
)
from evaluate_mtrack_emt import evaluate as evaluate_mtrack_emt
from build_lincs_mcf10a_phenotypes import build as build_lincs_mcf10a
from evaluate_lincs_mcf10a_phenotypes import evaluate as evaluate_lincs_mcf10a
from build_karacosta_emt_cytof import build as build_karacosta_emt_cytof
from evaluate_karacosta_emt_cytof import evaluate as evaluate_karacosta_emt_cytof
from build_geo_emt_rna import build as build_geo_emt_rna
from evaluate_geo_emt_rna import evaluate as evaluate_geo_emt_rna
from build_gse325309_pristine_holdout import build as build_gse325309_pristine
from evaluate_gse325309_pristine_holdout import evaluate as evaluate_gse325309_pristine
from audit_dryad_supracontractility_sources import (
    audit as audit_dryad_supracontractility_sources,
)
from build_dryad_supracontractility_manifest import (
    build as build_dryad_supracontractility,
)
from evaluate_dryad_supracontractility import (
    evaluate as evaluate_dryad_supracontractility,
)
from build_dataverse_cellular_nematics_manifest import (
    build as build_dataverse_cellular_nematics,
)
from evaluate_dataverse_cellular_nematics import (
    evaluate as evaluate_dataverse_cellular_nematics,
)
from audit_dryad_force_propagation_source import (
    audit as audit_dryad_force_propagation_source,
)
from build_dryad_force_propagation_manifest import (
    build as build_dryad_force_propagation,
)
from evaluate_dryad_force_propagation import (
    evaluate as evaluate_dryad_force_propagation,
)
from audit_sencid_reference import audit as audit_sencid_reference
from theory_runtime import (
    TheoryRefusal,
    persistent_random_walk_msd_2d,
    standard_linear_solid_relaxation,
    tether_apparent_tension,
    tether_force,
)
from train_calcium_baseline import train
from build_mdck_calcium_manifest import build as build_mdck
from build_eahy_calcium_manifest import build as build_eahy
from evaluate_piezo1_transfer import evaluate as evaluate_transfer
from evaluate_mechanism_contrasts import evaluate as evaluate_contrasts
from build_electrotaxis_piv_manifest import build as build_piv, OFFICIAL_SIZE as PIV_SIZE
from evaluate_electrotaxis_piv import evaluate as evaluate_piv
from build_actg1_if_pixel_manifest import pixel_metrics, build as build_if_pixels
from build_tendon_multimodal_manifest import build as build_tendon_multimodal
from evaluate_tendon_multimodal import evaluate as evaluate_tendon_multimodal
from build_gse226374_manifest import build as build_gse226374
from evaluate_fibrotic_state_transfer import evaluate as evaluate_fibrotic_transfer
from evaluate_sweep_calibration import (
    ALLOWED_AXES,
    SCHEMA as FORWARD_SCHEMA,
    evaluate as evaluate_sweep_calibration,
)
from evaluate_calcium_uncertainty import evaluate as evaluate_calcium_uncertainty
from evaluate_multisite_migration import evaluate as evaluate_multisite_migration
from build_multisite_migration_manifest import (
    PAPER_SELECTED_OBSERVABLES as MULTISITE_OBSERVABLES,
    _condition as multisite_condition,
    _unit as multisite_unit,
)
from build_dryad_migration_manifest import (
    BEHAVIORAL_FEATURES as DRYAD_BEHAVIORAL_FEATURES,
    SOURCE_SHA256 as DRYAD_SOURCE_SHA256,
    build as build_dryad_migration,
    build_tensor as build_dryad_tensor,
)
from evaluate_dryad_migration import evaluate as evaluate_dryad_migration
from build_elife_71032_manifest import build as build_elife_71032
from evaluate_rock_context_transfer import evaluate as evaluate_rock_context_transfer
from build_elife_72381_mechano_osmotic_manifest import (
    FIG2_ARCHIVE_SHA256 as ELIFE_72381_FIG2_SHA256,
    FIG3_ARCHIVE_SHA256 as ELIFE_72381_FIG3_SHA256,
    build as build_elife_72381,
)
from evaluate_mechano_osmotic_evidence import evaluate as evaluate_mechano_osmotic
from build_cell_monolayer_velocity_manifest import (
    BIOLOGICAL_GROUP as MONOLAYER_BIOLOGICAL_GROUP,
    build as build_cell_monolayer_velocity,
    build_tensor as build_cell_monolayer_tensor,
)
from evaluate_cell_monolayer_velocity import evaluate as evaluate_cell_monolayer_velocity
from evaluate_multitask_outer_model import (
    SCHEMA as MULTITASK_MODEL_SCHEMA,
    evaluate as evaluate_multitask_outer_model,
)
from train_calcium_temporal_cnn import train as train_calcium_temporal_cnn
from predict_calcium_temporal_cnn import (
    load_checkpoint as load_calcium_temporal_checkpoint,
    predict as predict_calcium_temporal_cnn,
)
from evaluate_piezo1_transfer import _metrics as calcium_binary_metrics
from build_bbbc_translocation_if_manifest import (
    BBBC013_ARCHIVE_SHA256,
    BBBC013_PLATEMAP_SHA256,
    BBBC014_ARCHIVE_SHA256,
    BBBC014_PLATEMAP_SHA256,
    FEATURE_NAMES as BBBC_IF_FEATURE_NAMES,
    build as build_bbbc_translocation_if,
    decode_bmp,
    field_features as bbbc_if_field_features,
    plate_labels as bbbc_plate_labels,
    _sha256 as bbbc_sha256,
)
from evaluate_bbbc_translocation_if import evaluate as evaluate_bbbc_translocation_if
from build_hpa_if_embedding_tensor import (
    OFFICIAL_ARCHIVE_SHA256 as HPA_ARCHIVE_SHA256,
    OFFICIAL_ARCHIVE_SIZE as HPA_ARCHIVE_SIZE,
    build as build_hpa_if_embedding,
)
from evaluate_hpa_if_localization import evaluate as evaluate_hpa_if_localization
from audit_outer_library_coverage import audit as audit_outer_library_coverage
from build_hpa_raw_image_subset import select as select_hpa_raw_image_subset
from build_bbbc054_microglia_manifest import (
    ANNOTATION_SHA256 as BBBC054_ANNOTATION_SHA256,
    ARCHIVE_SHA256 as BBBC054_ARCHIVE_SHA256,
    FEATURE_NAMES as BBBC054_FEATURE_NAMES,
    build as build_bbbc054_microglia,
    decode_uncompressed_tiff as decode_bbbc054_tiff,
)
from evaluate_bbbc054_microglia import evaluate as evaluate_bbbc054_microglia
from build_biad2515_apoptosis import build as build_biad2515_apoptosis
from download_figshare_cell_death import acquire as acquire_figshare_cell_death
from audit_figshare_cell_death_design import audit as audit_figshare_cell_death_design
from evaluate_biad2515_apoptosis import evaluate as evaluate_biad2515_apoptosis


PILOT = Path("data/external_training/experiment_factory/datasets/pilot")
WORKBOOK_EXPORT = Path(
    "data/external_training/experiment_factory/outer_library/nuclear_workbook.json"
)
MECH_WORKBOOK_EXPORT = Path(
    "data/external_training/experiment_factory/outer_library/mechanotransduction_workbooks.json"
)
TENDON_WORKBOOK_EXPORT = Path(
    "data/external_training/experiment_factory/outer_library/tendon_source_data_1.json"
)
GSE226374_COUNTS = Path(
    "data/external_training/experiment_factory/outer_library/downloads/gse226374/"
    "GSE226374_normalized_counts.txt.gz"
)
DRYAD_MIGRATION_ARCHIVE = Path(
    "data/external_training/experiment_factory/outer_library/downloads/dryad_9jh6m/"
    "Shafqat-Abbasi_Source_Datasets.zip"
)
ELIFE_71032_EXPORT = Path(
    "data/external_training/experiment_factory/outer_library/elife_71032_figure4.json"
)
ELIFE_72381_CONTROL_EXPORT = Path(
    "data/external_training/experiment_factory/outer_library/"
    "elife_72381_control_dVdt_dAdt.json"
)
ELIFE_72381_FIG2_ZIP = Path(
    "data/external_training/experiment_factory/outer_library/downloads/"
    "elife_72381/elife-72381-fig2-data1-v2.zip"
)
ELIFE_72381_FIG3_ZIP = Path(
    "data/external_training/experiment_factory/outer_library/downloads/"
    "elife_72381/elife-72381-fig3-data1-v2.zip"
)
ELIFE_72381_Y27_EXPORT = Path(
    "data/external_training/experiment_factory/outer_library/"
    "elife_72381_Y27_dVdt_dAdt.json"
)
ELIFE_72381_TETHER_EXPORT = Path(
    "data/external_training/experiment_factory/outer_library/"
    "elife_72381_tether_force_all.json"
)
ELIFE_72381_DATA_SUM_EXPORT = Path(
    "data/external_training/experiment_factory/outer_library/"
    "elife_72381_data_sum.json"
)
CELL_MONOLAYER_SOURCE = Path(
    "data/external_training/experiment_factory/outer_library/downloads/"
    "cell_monolayer_velocity"
)
BBBC013_PLATEMAP = Path(
    "data/external_training/experiment_factory/outer_library/downloads/"
    "bbbc013/BBBC013_v1_platemap_all.txt"
)
BBBC013_IMAGES = Path(
    "data/external_training/experiment_factory/outer_library/downloads/"
    "bbbc013/BBBC013_v1_images_bmp.zip"
)
BBBC014_PLATEMAP = Path(
    "data/external_training/experiment_factory/outer_library/downloads/"
    "bbbc014/BBBC014_v1_platemap_all.txt"
)
BBBC014_IMAGES = Path(
    "data/external_training/experiment_factory/outer_library/downloads/"
    "bbbc014/BBBC014_v1_images.zip"
)
HPA_IF_ARCHIVE = Path(
    "data/external_training/experiment_factory/outer_library/downloads/"
    "hpa_v25_1/subcell_image_umap_features.tsv.zip"
)
BBBC054_REPLICATE1 = Path(
    "data/external_training/experiment_factory/outer_library/downloads/"
    "bbbc054/Replicate_1.zip"
)
BBBC054_ANNOTATION = Path(
    "data/external_training/experiment_factory/outer_library/downloads/"
    "bbbc054/Replicate1annotation.csv"
)
SCIPLEX3_PDATA = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "sciplex3/GSM4150378_sciPlex3_pData.txt.gz"
)
SCIPLEX3_CELL_ANNOTATIONS = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "sciplex3/GSM4150378_sciPlex3_A549_MCF7_K562_screen_cell.annotations.txt.gz"
)
SCIPLEX3_MATRIX_RECEIPT = HERE / "results" / "sciplex3_matrix_receipt.json"
BIAD2515_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "biad2515"
)
GSE250041_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "gse250041"
)
FIGSHARE_CELL_DEATH_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "figshare_28202864"
)
SENCID_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "sencid/SenCID-code"
)
DRYAD_SUPRACONTRACTILITY_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "dryad_08kprr59c"
)
DRYAD_SUPRACONTRACTILITY_EXPORT = Path(
    "data/aleph/outer/experiment_factory/outer_library/"
    "dryad_08kprr59c_workbooks.json"
)
DATAVERSE_CELLULAR_NEMATICS_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "dataverse_data2772"
)
DRYAD_FORCE_PROPAGATION_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "dryad_sj3tx9683"
)
ZENODO_TFM_FRET_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "zenodo_14692589/Source-Data.zip"
)
FIGSHARE_EXIF_EMT_RAW = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "figshare_exif_emt/raw_exif_emt"
)
MTRACK_EMT_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "mtrack_emt_xy01"
)
LINCS_MCF10A_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/downloads/"
    "lincs_mcf10a"
)
KARACOSTA_EMT_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/"
    "karacosta_emt_cytof"
)
GEO_EMT_RNA_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/geo_emt_rna"
)
GSE325309_SOURCE = Path(
    "data/aleph/outer/experiment_factory/outer_library/gse325309/"
    "GSE325309_processed_data.xlsx"
)


def _result(name: str) -> dict:
    return json.loads((HERE / "results" / name).read_text(encoding="utf-8"))


def _multitask_inputs() -> list[dict]:
    return [
        _result(name) for name in (
            "calcium_baseline_report.json",
            "piezo1_transfer_report.json",
            "calcium_uncertainty_report.json",
            "mechanism_contrast_report.json",
            "multisite_migration_report.json",
            "dryad_migration_report.json",
            "fibrotic_state_transfer_report.json",
            "mechano_osmotic_evidence_report.json",
            "electrotaxis_piv_report.json",
            "cell_monolayer_velocity_report.json",
            "calcium_temporal_cnn_report.json",
            "bbbc_translocation_if_report.json",
            "hpa_if_localization_report.json",
            "bbbc054_microglia_report.json",
            "bbbc048_mlp_report.json",
            "bbbc053_mito_stress_report.json",
            "sciplex3_design_audit.json",
            "biad2515_apoptosis_report.json",
        )
    ]


def require(*paths: Path) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        pytest.skip(f"local licensed source data not present: {missing}")


def test_tfm_manifest_is_paired_and_non_authoritative():
    require(PILOT)
    rows = build(PILOT)
    assert len(rows) == 1_610
    assert len({row.sample_id for row in rows}) == 327
    assert {row.aleph_authority for row in rows} == {"none"}
    assert {row.label_authority for row in rows} == {"author_reported"}
    assert {row.modality for row in rows} == {"traction_force_microscopy_derived"}
    assert {row.observable: row.unit for row in rows} == {
        "projected_area": "square_micrometre",
        "mean_total_traction_force": "nN",
        "mean_abs_traction_force_x": "nN",
        "mean_abs_traction_force_y": "nN",
        "major_axis_length": "micrometre",
        "minor_axis_length": "micrometre",
    }
    assert all(row.unit_source_sha256 for row in rows)
    assert all("#Fig2;panel:E" in row.unit_source_locator for row in rows)


def test_tfm_manifest_never_uses_source_hash_as_sample_split():
    require(PILOT)
    rows = build(PILOT)
    sample_to_source = {}
    for row in rows:
        sample_to_source.setdefault(row.sample_id, row.source_sha256)
        assert sample_to_source[row.sample_id] == row.source_sha256
    assert len(sample_to_source) == 327


def test_nuclear_manifest_has_exact_cell_provenance():
    require(WORKBOOK_EXPORT)
    rows = build_nuclear(WORKBOOK_EXPORT)
    assert len(rows) > 1_000
    assert len({row.observation_id for row in rows}) == len(rows)
    assert {row.condition for row in rows} == {"WT", "BIX", "CHT", "DZNep"}
    assert all(row.source_locator.startswith("xlsx:") for row in rows)
    assert {row.aleph_authority for row in rows} == {"none"}


def test_nuclear_manifest_is_not_falsely_paired_across_figure_5_columns():
    require(WORKBOOK_EXPORT)
    rows = build_nuclear(WORKBOOK_EXPORT)
    figure5 = [row for row in rows if ":figure5:" in row.sample_id]
    assert figure5
    assert all(len({row.observable for row in figure5 if row.sample_id == sample}) == 1
               for sample in {row.sample_id for row in figure5})


def test_mechanotransduction_manifest_preserves_time_coordinates():
    require(MECH_WORKBOOK_EXPORT)
    rows = build_mechanotransduction(MECH_WORKBOOK_EXPORT)
    calcium = [row for row in rows if row.observable == "intracellular_calcium_F_over_F0"]
    assert len(calcium) == 33_900
    assert len({row.sample_id for row in calcium}) == 150
    assert {row.condition for row in calcium} == {
        "untreated_control", "Piezo1_inhibited_GsMTx4", "TRPV4_inhibited_HC-067047"
    }
    assert {row.coordinate_unit for row in calcium} == {"s"}
    assert all(row.aleph_authority == "none" for row in rows)


def test_mechanotransduction_fret_cell_type_is_author_header_label():
    require(MECH_WORKBOOK_EXPORT)
    rows = build_mechanotransduction(MECH_WORKBOOK_EXPORT)
    fret = [row for row in rows if row.modality == "live_cell_FRET_imaging_derived"]
    assert {row.cell_type for row in fret} == {"HeLa", "microvilli-deficient HeLa"}


def test_calcium_split_is_cell_grouped_and_gate_refuses_missing_lab_holdout():
    manifest = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "mechanotransduction_observations.jsonl"
    )
    require(manifest)
    report = train(manifest)
    assert report["sample_count"] == 150
    assert report["split"] == {"train": 90, "validation": 30, "test": 30}
    assert report["sample_overlap_count"] == 0
    gate = evaluate(report)
    assert gate["status"] == "refused"
    assert not gate["checks"]["independent_lab_task_holdout_present"]
    assert gate["may_mutate_physics"] is False


def test_mdck_external_calcium_is_author_aggregate_not_fake_single_cells():
    source = Path(
        "data/external_training/experiment_factory/outer_library/downloads/"
        "mdck_model_data_fit.zip"
    )
    require(source)
    rows = build_mdck(source)
    assert len(rows) == 826
    assert len({row.sample_id for row in rows}) == 7
    assert {row.biological_replicate for row in rows} == {"author_aggregate"}
    assert {row.cell_type for row in rows} == {"MDCK-II"}


def test_eahy_manifest_is_independent_lab_and_cell_grouped():
    source = Path(
        "data/external_training/experiment_factory/outer_library/downloads/"
        "piezo_nanoswitch/Figure_4_A_ii.csv"
    )
    require(source)
    rows = build_eahy(source)
    assert len(rows) == 62_996
    assert len({row.sample_id for row in rows}) == 625
    assert {row.lab_group for row in rows} == {"zenodo-19890985-authors"}
    assert all(row.aleph_authority == "none" for row in rows)


def test_transfer_evaluation_has_no_cross_lab_sample_overlap():
    train_manifest = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "mechanotransduction_observations.jsonl"
    )
    external_manifest = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "eahy_calcium_observations.jsonl"
    )
    require(train_manifest, external_manifest)
    report = evaluate_transfer(
        train_manifest,
        external_manifest,
    )
    assert report["external_lab_holdout"] is True
    assert report["sample_overlap_count"] == 0
    assert report["metrics"]["sample_count"] == 584
    aligned = report["aligned_active_nanoswitch_task"]
    assert aligned["external_metrics"]["sample_count"] == 311
    assert aligned["external_metrics"]["auroc"] >= 0.70
    assert aligned["external_metrics"]["balanced_accuracy"] >= 0.70
    assert report["task_holdout_passed"] is False


def test_mechanism_contrast_can_pass_without_granting_parameter_authority():
    discovery = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "mechanotransduction_observations.jsonl"
    )
    external = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "eahy_calcium_observations.jsonl"
    )
    require(discovery, external)
    report = evaluate_contrasts(discovery, external, draws=2_000)
    assert report["mechanism_holdout_passed"] is True
    assert report["external_confirmatory"]["direction_consistency"] == "2/2"
    assert report["observable_constraint_eligible"] is True
    assert report["aleph_parameter_proposal_eligible"] is False

    model_report = train(discovery)
    gate = evaluate(model_report, report)
    assert gate["observable_evidence_channel"]["status"] == "evidence_constraint_eligible"
    assert gate["parameter_axis_candidates"] == []
    assert gate["status"] == "refused"


def test_real_piv_manifest_uses_frame_groups_and_no_invented_velocity_unit():
    root = Path(
        "data/external_training/experiment_factory/outer_library/downloads/"
        "electrotaxis_piv"
    )
    control = root / "control_tissue_piv.mat"
    stimulated = root / "stimulated_tissue_piv.mat"
    require(control, stimulated)
    if (control.stat().st_size != PIV_SIZE["control"]
            or stimulated.stat().st_size != PIV_SIZE["stimulated"]):
        pytest.skip("official PIV downloads are present but incomplete")
    rows = build_piv(control, "control") + build_piv(stimulated, "stimulated")
    assert len(rows) == 960
    assert len({row.sample_id for row in rows}) == 120
    assert {row.condition for row in rows} == {"control", "stimulated"}
    assert {row.technical_replicate for row in rows} == {
        "not_reported_single_representative_tissue"
    }
    velocity_rows = [row for row in rows if "velocity" in row.observable or "speed" in row.observable]
    assert velocity_rows
    assert all(
        row.unit in {"source_velocity_unit_unknown", "dimensionless"}
        for row in velocity_rows
    )
    assert all(row.aleph_authority == "none" for row in rows)


def test_piv_evaluator_refuses_frame_pseudoreplication(tmp_path):
    manifest = tmp_path / "piv.jsonl"
    rows = []
    for condition, offset in (("control", 0.0), ("stimulated", 1.0)):
        for frame in range(60):
            for metric in (
                "mean_velocity_x", "mean_velocity_y", "mean_speed", "rms_speed",
                "speed_p90", "velocity_polar_order", "velocity_nematic_order",
                "valid_vector_fraction",
            ):
                rows.append({
                    "condition": condition,
                    "sample_id": f"{condition}:frame:{frame}",
                    "biological_replicate": "representative_tissue_single",
                    "observable": metric,
                    "value": offset + frame / 100.0,
                })
    manifest.write_text(
        "".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8"
    )
    report = evaluate_piv(manifest)
    assert report["frame_counts"] == {"control": 60, "stimulated": 60}
    assert report["biological_replicate_counts"] == {"control": 1, "stimulated": 1}
    assert report["metric_summaries"]["mean_speed"][
        "stimulated_minus_control_frame_mean"
    ] == pytest.approx(1.0)
    assert report["metric_summaries"]["mean_speed"]["inferential_ci95"] is None
    assert report["external_holdout_passed"] is False
    assert report["observable_constraint_eligible"] is False
    assert report["aleph_parameter_proposal_eligible"] is False


def test_if_pixel_metrics_are_deterministic_and_do_not_require_cell_labels():
    image = np.zeros((8, 8, 3), dtype=np.uint8)
    image[:4, :, 0] = 200
    image[:, :4, 1] = 100
    metrics = pixel_metrics(image)
    assert metrics["ACTG1_red_foreground_fraction"] == (0.5, "fraction")
    assert metrics["occludin_green_foreground_fraction"] == (0.5, "fraction")
    assert metrics["ACTG1_occludin_double_positive_fraction"] == (0.25, "fraction")


def test_if_pixel_manifest_keeps_mixed_culture_as_one_image_replicate():
    source = PILOT / "51d283dcbbcc426e7a544a6241fb6a532570e00b547b0736af4368deff2351ee"
    require(source)
    pytest.importorskip("PIL")
    rows = build_if_pixels(source)
    assert len({row.sample_id for row in rows}) == 1
    assert {row.biological_replicate for row in rows} == {"representative_image_single"}
    assert {row.condition for row in rows} == {"mixed_WT_and_ACTG1_KO_clone_2G3"}
    assert {row.modality for row in rows} == {"immunofluorescence_pixel_derived"}
    assert all(row.aleph_authority == "none" for row in rows)


def _grayscale_bmp(image: np.ndarray, *, top_down: bool = False) -> bytes:
    import struct

    image = np.asarray(image, dtype=np.uint8)
    height, width = image.shape
    row_bytes = ((width + 3) // 4) * 4
    palette = b"".join(bytes((value, value, value, 0)) for value in range(256))
    pixel_offset = 14 + 40 + len(palette)
    stored = image if top_down else image[::-1]
    pixels = b"".join(
        row.tobytes() + bytes(row_bytes - width) for row in stored
    )
    file_size = pixel_offset + len(pixels)
    header = b"BM" + struct.pack("<IHHI", file_size, 0, 0, pixel_offset)
    dib = struct.pack(
        "<IiiHHIIiiII",
        40,
        width,
        -height if top_down else height,
        1,
        8,
        0,
        len(pixels),
        0,
        0,
        256,
        0,
    )
    return header + dib + palette + pixels


def test_bbbc_bmp_decoder_preserves_orientation_palette_and_padding():
    image = np.asarray([[0, 1, 2], [250, 251, 252]], dtype=np.uint8)
    assert np.array_equal(decode_bmp(_grayscale_bmp(image)), image)
    assert np.array_equal(decode_bmp(_grayscale_bmp(image, top_down=True)), image)


def test_bbbc_if_features_pair_channels_before_emitting_field_observables():
    yy, xx = np.mgrid[:96, :96]
    nuclei = ((xx - 30) ** 2 + (yy - 48) ** 2 < 9 ** 2) | (
        (xx - 66) ** 2 + (yy - 48) ** 2 < 8 ** 2
    )
    dna = np.full((96, 96), 5, dtype=np.uint8)
    dna[nuclei] = 220
    signal = np.full((96, 96), 20, dtype=np.uint8)
    signal[nuclei] = 180
    features = bbbc_if_field_features(dna, signal)
    assert tuple(features) == BBBC_IF_FEATURE_NAMES
    assert features["nuclear_object_count"] == 2.0
    assert features["signal_nuclear_perinuclear_log_ratio"] > 1.0
    assert features["signal_nuclear_perinuclear_standardized_difference"] > 0.0


def test_bbbc_author_plate_maps_preserve_wells_doses_and_cell_types():
    require(BBBC013_PLATEMAP, BBBC014_PLATEMAP)
    assert hashlib.sha256(BBBC013_PLATEMAP.read_bytes()).hexdigest() == BBBC013_PLATEMAP_SHA256
    assert hashlib.sha256(BBBC014_PLATEMAP.read_bytes()).hexdigest() == BBBC014_PLATEMAP_SHA256
    labels_13 = bbbc_plate_labels("BBBC013", BBBC013_PLATEMAP)
    labels_14 = bbbc_plate_labels("BBBC014", BBBC014_PLATEMAP)
    assert len(labels_13) == len(labels_14) == 96
    assert labels_13[0].treatment == "Wortmannin"
    assert labels_13[48].treatment == "LY294002"
    assert labels_14[0].cell_type.startswith("MCF7")
    assert labels_14[48].cell_type.startswith("A549")
    assert labels_14[5].dose_text == "10^-9.523"
    assert labels_14[11].dose_value == pytest.approx(1e-13)


def test_bbbc_dataset_holdout_never_masquerades_as_independent_lab(tmp_path):
    feature_names = np.asarray(BBBC_IF_FEATURE_NAMES, dtype="U64")
    score_index = list(BBBC_IF_FEATURE_NAMES).index(
        "signal_nuclear_perinuclear_log_ratio"
    )
    matrices, datasets, rows, treatments, cell_types, doses = [], [], [], [], [], []
    for dataset_index, dataset_id in enumerate(
        ("bbbc013-v1-fkhr-translocation", "bbbc014-v1-nfkb-translocation")
    ):
        matrix = np.zeros((96, 12), dtype=np.float64)
        for index in range(96):
            row = chr(ord("A") + index // 12)
            column = index % 12
            if dataset_index == 0:
                treatment = "Wortmannin" if row <= "D" else "LY294002"
                cell_type = "U2OS"
                dose = float(column)
                target = column / 11.0
            else:
                treatment = "TNF-alpha"
                cell_type = "MCF7" if row <= "D" else "A549"
                dose = 10.0 ** (-7.0 - 6.0 * column / 11.0)
                target = 1.0 - column / 11.0
            row_offset = (ord(row) - ord("A")) * 0.002
            matrix[index] = np.linspace(0.1, 1.2, 12) * target
            matrix[index, score_index] = target + row_offset
            datasets.append(dataset_id)
            rows.append(row)
            treatments.append(treatment)
            cell_types.append(cell_type)
            doses.append(dose)
        matrices.append(matrix)
    X = np.concatenate(matrices)
    tensor = tmp_path / "bbbc_if_rows.npz"
    np.savez_compressed(
        tensor,
        X=X,
        feature_names=feature_names,
        dataset_id=np.asarray(datasets, dtype="U64"),
        lab_group=np.asarray(["ravkin-bbbc-translocation"] * 192, dtype="U64"),
        sample_id=np.asarray([f"sample-{i}" for i in range(192)], dtype="U32"),
        plate_row=np.asarray(rows, dtype="U1"),
        column=np.tile(np.arange(1, 13), 16),
        treatment=np.asarray(treatments, dtype="U32"),
        cell_type=np.asarray(cell_types, dtype="U64"),
        dose=np.asarray(doses, dtype=np.float64),
        dose_text=np.asarray([str(value) for value in doses], dtype="U24"),
        source_sha256=np.asarray(["0" * 64] * 192, dtype="U64"),
        platemap_sha256=np.asarray(["1" * 64] * 192, dtype="U64"),
    )
    report = evaluate_bbbc_translocation_if(tensor, bootstrap_draws=50)
    assert report["dataset_holdout"]["dataset_holdout_present"] is True
    assert report["dataset_holdout"]["translocation_representation_passed"] is True
    assert report["dataset_holdout"]["independent_lab_holdout"] is False
    assert report["aleph_authority"] == "none"
    assert report["may_select_aleph_parameter"] is False


def test_bbbc_real_archives_build_intact_well_grouped_fields():
    require(BBBC013_IMAGES, BBBC013_PLATEMAP, BBBC014_IMAGES, BBBC014_PLATEMAP)
    assert bbbc_sha256(BBBC013_IMAGES) == BBBC013_ARCHIVE_SHA256
    assert bbbc_sha256(BBBC014_IMAGES) == BBBC014_ARCHIVE_SHA256
    observations, tensor = build_bbbc_translocation_if(
        BBBC013_IMAGES,
        BBBC013_PLATEMAP,
        BBBC014_IMAGES,
        BBBC014_PLATEMAP,
    )
    assert tensor["X"].shape == (192, 12)
    assert len(observations) == 192 * 12
    assert len({row.sample_id for row in observations}) == 192
    assert {row.lab_group for row in observations} == {"ravkin-bbbc-translocation"}
    assert {row.biological_replicate for row in observations} == {
        "not_reported_single_plate"
    }
    assert {row.modality for row in observations} == {
        "IF_FKHR_EGFP_DRAQ_field",
        "IF_NFkB_FITC_DAPI_field",
    }
    assert all(row.label_authority == "author_reported" for row in observations)
    assert all(row.aleph_authority == "none" for row in observations)


def test_tendon_multimodal_manifest_preserves_force_units_and_unpaired_nuclei():
    require(TENDON_WORKBOOK_EXPORT)
    rows = build_tendon_multimodal(TENDON_WORKBOOK_EXPORT)
    force = [row for row in rows if row.modality == "cantilever_post_force_derived"]
    nuclei = [row for row in rows if row.modality == "nuclear_morphology_derived"]
    qpcr = [row for row in rows if row.modality == "RT_qPCR_derived"]
    assert len(force) == 89
    assert {row.unit for row in force} == {"nN_per_cell"}
    assert {row.cell_type for row in force} == {
        "rat tail tendon stromal cells",
        "human tendon fibroblasts",
        "human dermal fibroblasts",
    }
    assert nuclei
    assert {row.biological_replicate for row in nuclei} == {"unknown"}
    assert len(qpcr) == 210
    human_qpcr = [row for row in qpcr if row.cell_type == "human tendon fibroblasts"]
    assert len(human_qpcr) == 132
    assert all(row.biological_replicate.startswith("author_donor_column:")
               for row in human_qpcr)
    assert {row.unit for row in qpcr} == {
        "source_relative_expression_unit_unspecified"
    }
    assert len({row.sample_id for row in rows}) == len(rows)
    assert all(row.aleph_authority == "none" for row in rows)


def test_tendon_multimodal_evidence_passes_cell_type_not_lab_holdout():
    manifest = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "tendon_multimodal_observations.jsonl"
    )
    require(manifest)
    report = evaluate_tendon_multimodal(manifest, draws=2_000)
    assert report["multimodal_primary_coherence_passed"] is True
    assert report["human_fibrosis_force"]["cell_type_holdout_passed"] is True
    assert report["human_fibrosis_force"]["independent_lab_holdout"] is False
    assert report["independent_lab_task_holdout_passed"] is False
    assert report["aleph_parameter_proposal_eligible"] is False
    assert report["matrix_tension_contrasts"]["nuclear_circularity"][
        "evidence_eligible"
    ] is False

    model_report_path = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "calcium_baseline_report.json"
    )
    require(model_report_path)
    gate = evaluate(
        json.loads(model_report_path.read_text(encoding="utf-8")),
        tendon_report=report,
    )
    candidate = gate["parameter_axis_candidates"][0]
    assert len(gate["parameter_axis_candidates"]) == 5
    assert candidate["parameter"].endswith("MaterialCard.active_tension_pn")
    assert candidate["numeric_range"] is None
    assert candidate["may_emit_sweep_axis"] is False
    assert gate["status"] == "refused"


def test_gse226374_is_independent_lab_marker_panel_not_single_cells():
    require(GSE226374_COUNTS)
    rows = build_gse226374(GSE226374_COUNTS)
    assert len(rows) == 40
    assert len({row.sample_id for row in rows}) == 8
    assert {row.modality for row in rows} == {"bulk_RNA_seq_derived"}
    assert {row.lab_group for row in rows} == {"leask-lab-gse226374"}
    assert {row.condition for row in rows} == {
        "vehicle_control", "TGF_beta_1", "celastrol_only",
        "TGF_beta_1_plus_celastrol",
    }
    assert all(row.aleph_authority == "none" for row in rows)


def test_fibrotic_state_transfers_across_lab_species_and_modality():
    tendon_report = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "tendon_multimodal_report.json"
    )
    external_manifest = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "gse226374_observations.jsonl"
    )
    require(tendon_report, external_manifest)
    report = evaluate_fibrotic_transfer(tendon_report, external_manifest, draws=2_000)
    assert report["independent_lab_holdout"] is True
    assert report["independent_lab_mechanism_holdout_passed"] is True
    assert report["cell_state_classifier_holdout_passed"] is False
    assert report["aleph_parameter_proposal_eligible"] is False
    model_report_path = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "calcium_baseline_report.json"
    )
    require(model_report_path)
    gate = evaluate(
        json.loads(model_report_path.read_text(encoding="utf-8")),
        tendon_report=json.loads(tendon_report.read_text(encoding="utf-8")),
        fibrotic_report=report,
    )
    requirements = gate["calibration_requests"][0]["required_before_unblock"]
    assert not any("independent-lab" in requirement for requirement in requirements)
    assert gate["parameter_axis_candidates"][0]["may_emit_sweep_axis"] is False


def _synthetic_forward_calibration(observable_unit="nN_per_cell"):
    levels = ((0.0, 5.0), (1.0, 50.0), (2.0, 110.0), (3.0, 150.0))
    return {
        "schema": "aleph.outer_library.forward_sensitivity.v1",
        "engine": "Project_Aleph",
        "physics_commit": "synthetic-test-only",
        "scenario": "synthetic-test-only",
        "observable_operator": "synthetic-test-only",
        "observable_operator_status": "native_validated",
        "observable_operator_validation_sha256": "a" * 64,
        "axis": "aleph.vertical.sf_arc.MaterialCard.active_tension_pn",
        "axis_unit": "pN",
        "observable": "tissue_traction_force_per_cell",
        "observable_unit": observable_unit,
        "reference_axis_value": 0.0,
        "expected_direction": "increases",
        "runs": [
            {
                "run_id": f"level-{level}-seed-{seed}",
                "axis_value": level,
                "observable_value": response + offset,
                "converged": True,
            }
            for level, response in levels
            for seed, offset in enumerate((-1.0, 0.0, 1.0))
        ],
    }


def test_sweep_calibration_emits_only_inside_replicated_simulated_hull():
    tendon_report = json.loads(Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "tendon_multimodal_report.json"
    ).read_text(encoding="utf-8"))
    result = evaluate_sweep_calibration(
        _synthetic_forward_calibration(), tendon_report
    )
    assert result["status"] == "eligible_for_non_authoritative_bounded_axis"
    assert result["may_emit_sweep_axis"] is True
    assert result["within_simulated_hull"] is True
    assert 1.0 < result["numeric_axis_range"][0] < 2.0
    assert 2.0 < result["numeric_axis_range"][1] < 3.0
    assert result["may_mutate_physics"] is False


def test_sweep_calibration_refuses_traction_total_as_force_per_cell():
    tendon_report = json.loads(Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "tendon_multimodal_report.json"
    ).read_text(encoding="utf-8"))
    result = evaluate_sweep_calibration(
        _synthetic_forward_calibration(observable_unit="pN"), tendon_report
    )
    assert result["status"] == "refused"
    assert result["numeric_axis_range"] is None
    assert "observable_unit_matches_experiment" in result["failed_checks"]


def test_sweep_calibration_refuses_unvalidated_similarly_named_operator():
    tendon_report = json.loads(Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "tendon_multimodal_report.json"
    ).read_text(encoding="utf-8"))
    calibration = _synthetic_forward_calibration()
    calibration["observable_operator_status"] = "candidate_native"
    result = evaluate_sweep_calibration(calibration, tendon_report)
    assert result["status"] == "refused"
    assert "observable_operator_is_validated" in result["failed_checks"]
    assert result["may_emit_sweep_axis"] is False


def test_sweep_calibration_preserves_typed_upstream_refusal():
    calibration = _synthetic_forward_calibration()
    calibration["runs"] = []
    calibration["refusal"] = {
        "refused": True,
        "code": "UPSTREAM_DEFECT_OPEN",
        "reason": "load-path precondition failed",
        "operator_id": "aleph.observe.traction.TractionPerCellOperator/v1",
        "detail": {"n_offending_rows": 2},
    }
    result = evaluate_sweep_calibration(calibration, {})
    assert result["status"] == "refused"
    assert result["refusal_is_well_formed"] is True
    assert result["upstream_refusal"] == calibration["refusal"]
    assert result["may_emit_sweep_axis"] is False


def test_tether_observable_cannot_be_misfit_to_existing_tension_axis():
    report = json.loads(Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "mechano_osmotic_evidence_report.json"
    ).read_text(encoding="utf-8"))
    calibration = _synthetic_forward_calibration(observable_unit="pN")
    calibration["observable"] = "membrane_tether_plateau_force"
    result = evaluate_sweep_calibration(calibration, report)
    assert result["status"] == "refused"
    assert "observable_axis_pair_is_identifiable" in result["failed_checks"]
    assert result["may_emit_sweep_axis"] is False


def test_forward_sensitivity_schema_and_evaluator_allow_the_same_axes():
    schema = json.loads((HERE / "forward_sensitivity.schema.json").read_text(
        encoding="utf-8"
    ))
    assert schema["properties"]["schema"]["const"] == FORWARD_SCHEMA
    assert set(schema["properties"]["axis"]["enum"]) == set(ALLOWED_AXES)
    assert set(schema["properties"]["observable"]["enum"]) == {
        "tissue_traction_force_per_cell",
        "membrane_tether_plateau_force",
    }


def test_external_calcium_mondrian_conformal_restores_training_only_coverage():
    train_manifest = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "mechanotransduction_observations.jsonl"
    )
    external_manifest = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "eahy_calcium_observations.jsonl"
    )
    require(train_manifest, external_manifest)
    report = evaluate_calcium_uncertainty(train_manifest, external_manifest)
    assert report["sample_overlap_count"] == 0
    assert report["split"] == {"train": 50, "validation": 20, "calibration": 30}
    assert report["conformal"]["external_count"] == 311
    assert report["conformal"]["external_coverage"] >= 0.90
    assert report["conformal"]["singleton_fraction"] >= 0.30
    assert report["conformal"]["singleton_accuracy"] >= 0.70
    assert report["conformal"]["pooled_diagnostic"]["external_coverage"] < 0.90
    assert report["dataset_shift"]["external_vs_calibration_auroc"] >= 0.70
    assert report["uncertainty_holdout_passed"] is True
    assert report["aleph_parameter_proposal_eligible"] is False
    model_report = json.loads(Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "calcium_baseline_report.json"
    ).read_text(encoding="utf-8"))
    gate = evaluate(model_report, uncertainty_report=report)
    assert gate["checks"]["external_lab_uncertainty_holdout_passed"] is True
    assert gate["status"] == "refused"


def _write_synthetic_multisite_manifest(path: Path, *, reverse_lab3: bool = False) -> None:
    rows = []
    observables = ("cell_speed", "persistence", "cell_area")
    treatment_effect = np.asarray((4.0, 2.0, -1.0))
    for lab_index, lab in enumerate(("lab_1", "lab_2", "lab_3"), start=1):
        lab_offset = np.asarray((64.0, -32.0, 16.0)) * lab_index
        for experiment_index in range(9):
            experiment = (
                f"{lab}:person_{experiment_index // 3}:experiment_{experiment_index}"
            )
            experiment_offset = np.asarray((0.5, -0.25, 0.125)) * experiment_index
            for condition, label in (
                ("control", 0), ("ROCK_inhibited", 1)
            ):
                effect = treatment_effect * label
                if reverse_lab3 and lab == "lab_3":
                    effect *= -1
                for technical_index in range(3):
                    technical = f"technical_{technical_index}"
                    technical_noise = np.asarray((0.0625, -0.03125, 0.015625)) * (
                        technical_index - 1
                    )
                    vector = lab_offset + experiment_offset + effect + technical_noise
                    for time_index in range(4):
                        time_noise = np.asarray((0.015625, 0.0078125, -0.00390625)) * (
                            time_index - 1.5
                        )
                        for observable, value in zip(observables, vector + time_noise):
                            rows.append({
                                "lab_group": lab,
                                "biological_replicate": experiment,
                                "technical_replicate": technical,
                                "condition": condition,
                                "observable": observable,
                                "value": float(value),
                                "time_index": time_index,
                            })
    path.write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
        encoding="utf-8",
    )


def test_multisite_migration_uses_a_sealed_third_lab_and_experiment_inference(tmp_path):
    manifest = tmp_path / "multisite.jsonl"
    _write_synthetic_multisite_manifest(manifest)
    report = evaluate_multisite_migration(manifest, draws=2_000)
    assert report["protocol"]["train_lab"] == "lab_1"
    assert report["protocol"]["validation_and_conformal_calibration_lab"] == "lab_2"
    assert report["protocol"]["decision_threshold"] == 0.0
    assert report["protocol"]["external_test_lab"] == "lab_3"
    assert report["protocol"]["test_lab_labels_used_for_numeric_fitting"] is False
    assert report["pristine_sealed_confirmation"] is False
    assert report["experiment_counts"] == {"lab_1": 9, "lab_2": 9, "lab_3": 9}
    assert report["person_counts"] == {"lab_1": 3, "lab_2": 3, "lab_3": 3}
    assert report["external_lab3_metrics"]["balanced_accuracy"] >= 0.70
    assert report["conformal_uncertainty"]["coverage"] >= 0.90
    assert report["conformal_uncertainty"]["uncertainty_holdout_passed"] is True
    assert report["external_lab3_effect_bootstrap_95ci"][0] > 0
    assert report["external_lab_holdout_passed"] is True
    assert report["aleph_parameter_proposal_eligible"] is False


def test_multisite_migration_rejects_reversed_sealed_lab_effect(tmp_path):
    manifest = tmp_path / "multisite_reversed.jsonl"
    _write_synthetic_multisite_manifest(manifest, reverse_lab3=True)
    report = evaluate_multisite_migration(manifest, draws=2_000)
    assert report["external_lab3_effect_bootstrap_95ci"][1] < 0
    assert report["external_lab_holdout_passed"] is False
    assert report["aleph_parameter_proposal_eligible"] is False


def test_multisite_manifest_declares_only_author_supported_labels_and_units():
    assert len(MULTISITE_OBSERVABLES) == 18
    assert multisite_condition("C") == ("control", "untreated_migrating_HT1080")
    assert multisite_condition("T") == (
        "ROCK_inhibited", "ROCK_inhibited_migrating_HT1080"
    )
    assert multisite_unit("Cells_Area") == "square_micrometre"
    assert multisite_unit("Cells_Perimeter") == "micrometre"
    assert multisite_unit("Cells_ICS") == "author_scaled_source_unit"


def test_real_multisite_manifest_preserves_independent_experiments_and_sealed_lab():
    manifest = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "multisite_migration_observations.jsonl"
    )
    require(manifest)
    report = evaluate_multisite_migration(manifest, draws=2_000)
    assert report["observable_count"] == 18
    assert report["experiment_counts"] == {"lab_1": 9, "lab_2": 9, "lab_3": 9}
    assert report["person_counts"] == {"lab_1": 3, "lab_2": 3, "lab_3": 3}
    assert report["external_lab3_metrics"]["auroc"] >= 0.98
    assert report["external_lab3_metrics"]["balanced_accuracy"] >= 0.74
    assert report["external_lab3_metrics"]["macro_f1"] >= 0.72
    assert report["conformal_uncertainty"]["coverage"] >= 0.90
    assert report["conformal_uncertainty"]["singleton_fraction"] >= 0.60
    assert report["conformal_uncertainty"]["uncertainty_holdout_passed"] is True
    assert report["external_lab3_effect_bootstrap_95ci"][0] > 0
    assert report["external_lab_holdout_passed"] is True
    assert report["aleph_parameter_proposal_eligible"] is False
    model_report_path = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "calcium_baseline_report.json"
    )
    require(model_report_path)
    gate = evaluate(
        json.loads(model_report_path.read_text(encoding="utf-8")),
        migration_report=report,
    )
    assert gate["checks"][
        "independent_multisite_migration_external_evaluation_passed"
    ] is True
    assert gate["checks"][
        "independent_multisite_migration_retrospective_uncertainty_passed"
    ] is True
    assert gate["checks"][
        "independent_multisite_migration_pristine_confirmation_present"
    ] is False
    constraints = gate["observable_evidence_channel"]["constraints"]
    assert constraints[-1]["observable"] == "multivariate_live_cell_ROCK_inhibition_signature"
    assert constraints[-1]["authority"] == "evidence_only"
    assert gate["status"] == "refused"
    assert gate["may_mutate_physics"] is False


def test_dryad_migration_manifest_preserves_date_cell_hierarchy_and_author_scope():
    require(DRYAD_MIGRATION_ARCHIVE)
    rows = build_dryad_migration(DRYAD_MIGRATION_ARCHIVE)
    assert len(rows) == 22_336
    assert len({row.sample_id for row in rows}) == 367
    assert len({row.biological_replicate for row in rows}) == 33
    assert len({row.dataset_id for row in rows}) == 3
    assert {row.source_sha256 for row in rows} == {DRYAD_SOURCE_SHA256}
    assert {row.lab_group for row in rows} == {"stromblad-lab-karolinska"}
    assert {row.aleph_authority for row in rows} == {"none"}
    assert all(row.technical_replicate.startswith("tracked_cell:") for row in rows)
    assert all(row.biological_replicate.startswith("author_experimental_date:") for row in rows)


def test_dryad_row_tensor_is_groupable_and_excludes_behavior_from_mode_model(tmp_path):
    require(DRYAD_MIGRATION_ARCHIVE)
    tensor = tmp_path / "dryad_rows.npz"
    summary = build_dryad_tensor(DRYAD_MIGRATION_ARCHIVE, tensor)
    assert summary["rows"] == 22_647
    assert summary["tracked_cells"] == 367
    assert summary["features"] == 60
    report = evaluate_dryad_migration(tensor, draws=2_000)
    protocol = report["protocol"]
    assert protocol["label_defining_behavioral_features_excluded"] == sorted(
        DRYAD_BEHAVIORAL_FEATURES
    )
    assert protocol["selected_organizational_feature_count"] == 55
    assert report["retrospective_talin_tracked_cell_metrics"]["auroc"] >= 0.98
    assert report["retrospective_talin_tracked_cell_metrics"]["balanced_accuracy"] >= 0.90
    assert report["conformal_uncertainty"]["coverage"] >= 0.90
    assert report["conformal_uncertainty"]["singleton_fraction"] >= 0.40
    assert report["same_lab_cross_perturbation_transfer_passed"] is True
    assert report["independent_lab_holdout"] is False
    assert report["pristine_lab4_confirmation"] is False
    assert report["aleph_parameter_proposal_eligible"] is False


def test_elife_71032_manifest_keeps_cells_distinct_from_unknown_experiments():
    require(ELIFE_71032_EXPORT)
    rows = build_elife_71032(ELIFE_71032_EXPORT)
    assert len(rows) == 106
    assert len({row.sample_id for row in rows}) == 53
    assert {row.observable for row in rows} == {
        "persistent_speed", "speed_oscillation_period"
    }
    assert {row.biological_replicate for row in rows} == {"not_reported"}
    assert {row.lab_group for row in rows} == {
        "riveline-godeau-strasbourg-curie"
    }
    assert {row.aleph_authority for row in rows} == {"none"}


def test_rock_context_transfer_rejects_a_universal_speed_constraint(tmp_path):
    require(ELIFE_71032_EXPORT)
    elife_manifest = tmp_path / "elife.jsonl"
    elife_manifest.write_text(
        "".join(
            json.dumps(row.as_dict(), sort_keys=True) + "\n"
            for row in build_elife_71032(ELIFE_71032_EXPORT)
        ),
        encoding="utf-8",
    )
    dryad_manifest = Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "dryad_migration_observations.jsonl"
    )
    require(dryad_manifest)
    report = evaluate_rock_context_transfer(dryad_manifest, elife_manifest, draws=2_000)
    assert report["dryad_H1299_PL_2D"]["direction"] == "increases"
    assert report["elife_NIH3T3_3D_CDM"]["direction"] == "decreases"
    assert report["context_heterogeneity_detected"] is True
    assert report["independent_lab_dataset_present"] is True
    assert report["independent_lab_mechanism_holdout_passed"] is False
    assert report["observable_constraint_eligible"] is False
    assert report["aleph_parameter_proposal_eligible"] is False
    model_report = json.loads(Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "calcium_baseline_report.json"
    ).read_text(encoding="utf-8"))
    gate = evaluate(model_report, rock_context_report=report)
    assert gate["evidence_exclusions"][0]["status"] == "rejected_context_dependent"
    assert gate["evidence_exclusions"][0]["may_select_aleph_parameter"] is False
    assert gate["status"] == "refused"


def test_elife_72381_manifest_preserves_cells_inside_author_experiments():
    require(
        ELIFE_72381_CONTROL_EXPORT,
        ELIFE_72381_Y27_EXPORT,
        ELIFE_72381_TETHER_EXPORT,
        ELIFE_72381_FIG2_ZIP,
        ELIFE_72381_FIG3_ZIP,
    )
    rows = build_elife_72381(
        ELIFE_72381_FIG2_ZIP,
        ELIFE_72381_FIG3_ZIP,
        ELIFE_72381_CONTROL_EXPORT,
        ELIFE_72381_Y27_EXPORT,
        ELIFE_72381_TETHER_EXPORT,
    )
    spreading = [
        row for row in rows
        if row.observable in {"initial_spreading_area_rate", "initial_volume_flux"}
    ]
    assert spreading
    assert all(row.unit_source_sha256 == ELIFE_72381_FIG2_SHA256 for row in spreading)
    assert {row.unit_source_locator for row in spreading} == {
        "zip:Figure2/Figure2D/data_sum.xlsx!Sheet1!B1",
        "zip:Figure2/Figure2D/data_sum.xlsx!Sheet1!E1",
    }
    assert len(rows) == 944
    assert len({row.sample_id for row in rows}) == 629
    assert {row.source_sha256 for row in rows} == {
        ELIFE_72381_FIG2_SHA256,
        ELIFE_72381_FIG3_SHA256,
    }
    spreading = [
        row for row in rows
        if row.modality == "FXm_RICM_cell_spreading_derived"
    ]
    control = [row for row in spreading if row.condition == "vehicle_control"]
    treated = [
        row for row in spreading
        if row.condition == "Y27632_100_micromolar"
    ]
    assert len({row.sample_id for row in control}) == 194
    assert len({row.sample_id for row in treated}) == 121
    assert len({row.biological_replicate for row in control}) == 3
    assert len({row.biological_replicate for row in treated}) == 4
    tether = [
        row for row in rows
        if row.modality == "AFM_membrane_tether_force_derived"
    ]
    assert len(tether) == 314
    assert {row.unit for row in tether} == {"pN"}
    assert all(row.aleph_authority == "none" for row in rows)


def test_mechano_osmotic_gate_uses_experiment_effects_but_refuses_direct_tension(tmp_path):
    require(
        ELIFE_72381_CONTROL_EXPORT,
        ELIFE_72381_Y27_EXPORT,
        ELIFE_72381_TETHER_EXPORT,
        ELIFE_72381_FIG2_ZIP,
        ELIFE_72381_FIG3_ZIP,
        ELIFE_72381_DATA_SUM_EXPORT,
    )
    manifest = tmp_path / "elife_72381.jsonl"
    manifest.write_text(
        "".join(
            json.dumps(row.as_dict(), sort_keys=True) + "\n"
            for row in build_elife_72381(
                ELIFE_72381_FIG2_ZIP,
                ELIFE_72381_FIG3_ZIP,
                ELIFE_72381_CONTROL_EXPORT,
                ELIFE_72381_Y27_EXPORT,
                ELIFE_72381_TETHER_EXPORT,
            )
        ),
        encoding="utf-8",
    )
    report = evaluate_mechano_osmotic(
        manifest, ELIFE_72381_DATA_SUM_EXPORT, draws=2_000
    )
    assert report["initial_spreading_area_rate"][
        "independent_experiment_bootstrap_95ci"
    ][0] > 0
    assert report["initial_volume_flux"][
        "independent_experiment_bootstrap_95ci"
    ][1] < 0
    early = report["AFM_membrane_tether_force"]["spreading_phase_30_to_90_min"]
    steady = report["AFM_membrane_tether_force"]["steady_state_spread_4_to_5_hour"]
    assert early["paired_experiment_bootstrap_95ci_pN"][0] > 0
    assert steady["paired_experiment_bootstrap_95ci_pN"][0] < 0
    assert steady["paired_experiment_bootstrap_95ci_pN"][1] > 0
    assert report["observable_constraint_eligible"] is True
    reconciliation = report["source_reconciliation"]
    assert reconciliation[
        "control_reconciles_after_one_explicit_author_summary_exclusion"
    ] is True
    assert reconciliation["Y27632_volume_flux_reconciles"] is True
    assert reconciliation["Y27632_area_rate_raw_table_minus_summary_mean"] == pytest.approx(
        0.676071981487446
    )
    assert report["tether_force_is_direct_cortical_tension"] is False
    assert report["aleph_parameter_proposal_eligible"] is False

    model_report = json.loads(Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "calcium_baseline_report.json"
    ).read_text(encoding="utf-8"))
    gate = evaluate(model_report, mechano_osmotic_report=report)
    assert gate["status"] == "refused"
    assert gate["observable_evidence_channel"]["constraints"][-1][
        "observable"
    ] == "HeLa_mechano_osmotic_spreading_and_AFM_tether_signature"
    assert gate["evidence_exclusions"][-1]["status"] == "rejected_operator_mismatch"
    candidate = gate["parameter_axis_candidates"][-1]
    assert candidate["parameter"] == (
        "requested_nonexistent_axis:membrane_cortex_attachment_energy_pn_per_um"
    )
    assert candidate["status"] == "blocked_missing_physics_term_and_exact_nonidentifiability"
    assert candidate["numeric_range"] is None
    assert candidate["may_emit_sweep_axis"] is False
    assert gate["calibration_requests"][-1]["evaluator"] == "evaluate_sweep_calibration.py"
    assert gate["may_mutate_physics"] is False


def test_cell_monolayer_velocity_manifest_keeps_technical_points_grouped(tmp_path):
    require(CELL_MONOLAYER_SOURCE)
    rows = build_cell_monolayer_velocity(CELL_MONOLAYER_SOURCE)
    assert len(rows) == 3_154
    assert len({row.sample_id for row in rows}) == 832
    assert {row.cell_type for row in rows} == {"human dermal fibroblasts (hdF)"}
    assert {row.biological_replicate for row in rows} == {
        MONOLAYER_BIOLOGICAL_GROUP
    }
    assert {row.lab_group for row in rows} == {"valentine-lab-ucsb"}
    assert {row.aleph_authority for row in rows} == {"none"}
    tensor = tmp_path / "cell_monolayer_velocity.npz"
    result = build_cell_monolayer_tensor(CELL_MONOLAYER_SOURCE, tensor)
    assert result == {
        "vectors": 7_645,
        "features": 5,
        "split_group": "valentine-lab-ucsb",
        "biological_group_count": 1,
    }
    payload = np.load(tensor, allow_pickle=False)
    assert payload["vectors"].shape == (7_645, 5)
    assert np.isnan(payload["vectors"][:, 4]).sum() == 1
    assert set(payload["split_group"].tolist()) == {"valentine-lab-ucsb"}


def test_cell_monolayer_velocity_report_refuses_false_holdout_and_fak_causality(tmp_path):
    require(CELL_MONOLAYER_SOURCE)
    manifest = tmp_path / "cell_monolayer_velocity.jsonl"
    manifest.write_text(
        "".join(
            json.dumps(row.as_dict(), sort_keys=True) + "\n"
            for row in build_cell_monolayer_velocity(CELL_MONOLAYER_SOURCE)
        ),
        encoding="utf-8",
    )
    tensor = tmp_path / "cell_monolayer_velocity.npz"
    build_cell_monolayer_tensor(CELL_MONOLAYER_SOURCE, tensor)
    report = evaluate_cell_monolayer_velocity(manifest, tensor)
    assert report["raw_velocity_tensor"]["technical_vectors"] == 7_645
    assert report["raw_velocity_tensor"]["biological_replicate_count"] == 1
    assert report["representation_training_eligible"] is True
    assert report["independent_lab_holdout"] is False
    assert report["uncertainty_holdout_passed"] is False
    assert report["observable_constraint_eligible"] is False
    assert report["aleph_parameter_proposal_eligible"] is False
    assert report["FAKi_descriptive_contrasts"]["cell_number_density"][
        "FAKi_minus_control_timepoint_mean"
    ] < -100
    assert report["FAKi_descriptive_contrasts"]["ensemble_mean_speed"][
        "inferential_ci95"
    ] is None
    model_report = json.loads(Path(
        "aleph/outer/experiment_factory/outer_library/results/"
        "calcium_baseline_report.json"
    ).read_text(encoding="utf-8"))
    gate = evaluate(model_report, cell_monolayer_report=report)
    audit = gate["representation_training_audits"][0]
    assert audit["status"] == "admitted_training_only"
    assert audit["technical_vector_count"] == 7_645
    assert audit["may_count_as_external_holdout"] is False
    assert audit["may_select_aleph_parameter"] is False
    assert gate["status"] == "refused"


def test_multitask_outer_model_keeps_task_specific_validation_authority():
    report = evaluate_multitask_outer_model(*_multitask_inputs())
    assert report["schema"] == MULTITASK_MODEL_SCHEMA
    assert report["architecture"]["shared_latent"] is False
    assert report["architecture"]["trainable_raw_encoder_status"] == (
        "two_encoders_trained_no_external_acceptance"
    )
    assert report["architecture"]["trained_encoder_heads"] == [
        "calcium_cell_state", "cell_cycle_image_representation"
    ]
    assert report["readiness"]["task_scoped_evidence_runtime_ready"] is True
    assert report["readiness"]["general_cell_type_state_model_ready"] is False
    assert report["readiness"]["production_eligible"] is False
    heads = {head["head_id"]: head for head in report["heads"]}
    assert heads["calcium_cell_state"]["status"] == "blocked_external_classification"
    assert heads["calcium_cell_state"]["uncertainty_eligible"] is True
    assert heads["calcium_cell_state"]["observable_evidence_eligible"] is False
    cnn = heads["calcium_cell_state"]["metrics"][
        "trainable_temporal_CNN_diagnostic"
    ]
    assert cnn["architecture"]["filters_are_trainable"] is True
    assert cnn["architecture_accepted"] is False
    assert heads["paired_migration_state"]["status"] == (
        "validated_retrospective_external_evidence_only"
    )
    assert heads["paired_migration_state"]["pristine_confirmation"] is False
    assert heads["fibrotic_state_mechanism"]["observable_evidence_eligible"] is True
    assert heads["velocity_field_representation"]["status"] == (
        "training_only_no_biological_replication"
    )
    assert heads["hpa_if_localization_representation"]["status"] == (
        "training_only_grouped_same_provider"
    )
    assert heads["hpa_if_localization_representation"][
        "observable_evidence_eligible"
    ] is False
    assert heads["microglia_immune_morphophenotype_representation"]["status"] == (
        "training_only_single_replicate_classifier_rejected"
    )
    assert heads["cell_cycle_image_representation"]["status"] == (
        "training_only_same_provider_rare_phase_limited"
    )
    assert heads["mitochondrial_stress_IF_representation"]["status"] == (
        "training_only_no_biological_replication"
    )
    assert heads["sciplex_perturbation_state_representation"]["status"] == (
        "training_only_model_pending"
    )
    assert heads["apoptosis_AnnexinV_endpoint_representation"]["status"] == (
        "training_only_same_provider_well_holdout"
    )
    assert heads["apoptosis_AnnexinV_endpoint_representation"]["metrics"][
        "sealed_same_provider_test"
    ]["r2"] > 0.87
    assert all(not head["parameter_proposal_eligible"] for head in heads.values())
    assert set(report["sweep_interface"]["admissible_head_ids"]) == {
        "piezo1_calcium_mechanism",
        "paired_migration_state",
        "fibrotic_state_mechanism",
    }
    assert report["sweep_interface"]["may_select_aleph_parameter"] is False


def test_hpa_if_official_embedding_archive_builds_without_gene_antibody_leakage(tmp_path):
    require(HPA_IF_ARCHIVE)
    assert HPA_IF_ARCHIVE.stat().st_size == HPA_ARCHIVE_SIZE
    tensor = tmp_path / "hpa_if.npz"
    receipt = build_hpa_if_embedding(HPA_IF_ARCHIVE, tensor)
    assert receipt["archive_sha256"] == HPA_ARCHIVE_SHA256
    assert receipt["image_count"] == 81_007
    assert receipt["gene_count"] == 13_366
    assert receipt["antibody_count"] == 15_700
    assert receipt["gene_or_antibody_cross_split_overlap_count"] == 0
    with np.load(tensor, allow_pickle=False) as data:
        assert data["X"].shape == (81_007, 128)
        assert set(data["split"].astype(str)) == {
            "train", "validation", "calibration", "test"
        }


def test_hpa_if_grouped_holdout_is_training_only_not_external_authority():
    tensor = HERE / "results" / "hpa_if_embedding_rows.npz"
    require(tensor)
    report = evaluate_hpa_if_localization(tensor)
    assert report == _result("hpa_if_localization_report.json")
    assert report["split"]["gene_or_antibody_cross_split_overlap_count"] == 0
    assert report["sealed_grouped_test_metrics"]["macro_auroc"] >= 0.90
    assert report["sealed_grouped_test_metrics"]["micro_f1"] >= 0.75
    assert report["independent_dataset_holdout"] is False
    assert report["independent_lab_holdout"] is False
    assert report["conformal_diagnostic"]["eligible_for_uncertainty_authority"] is False
    assert report["may_select_aleph_parameter"] is False


def test_adversarial_coverage_audit_does_not_confuse_rows_with_samples():
    report = audit_outer_library_coverage(HERE / "results")
    assert report == _result("adversarial_coverage_audit.json")
    assert report["counts"]["canonical_observations"] == 329_303
    assert report["counts"]["canonical_samples"] == 193_832
    assert report["counts"]["datasets_including_HPA_embedding_source"] == 38
    assert report["counts"]["conservative_lab_groups_including_HPA"] == 37
    assert report["counts"]["HPA_embedding_images_separate_from_physical_observations"] == 81_007
    assert report["provenance_defects"]["unknown_or_author_processed_unit_observations"] == 29_055
    assert report["provenance_defects"]["unknown_unit_fraction"] > 0.08
    assert report["provenance_defects"]["biological_replicate_unreported_observations"] == 238_071
    assert report["provenance_defects"]["biological_replicate_unreported_fraction"] > 0.72
    assert "cell_cycle_and_proliferation" not in report["missing_state_axes"]
    assert report["state_axis_coverage"]["cell_cycle_and_proliferation"]["status"] == (
        "author_labelled_single_provider_training_only"
    )
    assert "EMT_epithelial_mesenchymal_state" not in report["missing_state_axes"]
    assert report["state_axis_coverage"]["EMT_epithelial_mesenchymal_state"][
        "status"
    ] == "ten_lab_sources_with_preregistered_GSE325309_direction_holdout_passed_but_no_numeric_OOD_or_Aleph_operator"
    assert "immune_activation_and_inflammation" not in report["missing_state_axes"]
    assert "mitochondrial_metabolic_stress" not in report["missing_state_axes"]
    assert "hypoxia_state" in report["missing_state_axes"]
    assert "apoptosis_and_cell_death" not in report["missing_state_axes"]
    assert "senescence" not in report["missing_state_axes"]
    assert report["state_axis_coverage"]["senescence"]["status"] == (
        "directionally_confirmed_across_labs_absolute_state_OOD_blocked"
    )
    assert report["state_axis_coverage"][
        "supracellular_contractility_and_ECM_alignment"
    ]["status"] == (
        "three_lab_cross_cell_type_direction_operator_and_biological_replication_blocked"
    )
    assert report["state_axis_coverage"]["apoptosis_and_cell_death"]["status"] == (
        "author_AnnexinV_endpoint_single_provider_well_holdout_training_only"
    )
    assert report["production_eligible"] is False


def test_sciplex3_design_preserves_context_groups_and_unassigned_cells():
    require(SCIPLEX3_PDATA, SCIPLEX3_CELL_ANNOTATIONS, SCIPLEX3_MATRIX_RECEIPT)
    report = audit_sciplex3_design(
        SCIPLEX3_PDATA,
        SCIPLEX3_CELL_ANNOTATIONS,
        json.loads(SCIPLEX3_MATRIX_RECEIPT.read_text(encoding="utf-8")),
    )
    assert report == _result("sciplex3_design_audit.json")
    assert report["cell_count"] == 799_317
    assert report["assigned_context_cell_count"] == 762_795
    assert report["unassigned_context_cell_count"] == 36_522
    assert report["mandatory_context_group_count"] == 4_896
    assert report["replicate_holdout_is_same_lab_not_external_lab"] is True
    assert report["expression_matrix_acquired"] is True
    assert report["training_ready"] is True
    assert report["matrix_evidence"]["sha256"] == (
        "7d632716aa6ed0fc1780996003d0abc440ff78340609bceaaa4c8ade9345d00a"
    )
    verified = audit_sciplex3_design(
        SCIPLEX3_PDATA,
        SCIPLEX3_CELL_ANNOTATIONS,
        {
            "schema": "aleph.outer_library.sciplex3_matrix_receipt.v1",
            "accession": "GSM4150378",
            "size_bytes": 3_096_224_606,
            "sha256": "a" * 64,
            "gzip_stream_verified": True,
            "decompressed_line_count": 2_000_000,
            "aleph_authority": "none",
        },
    )
    assert verified["expression_matrix_acquired"] is True
    assert verified["training_ready"] is True
    assert verified["reason_not_training_ready"] is None


def test_reactome_state_panel_is_source_derived_and_sparse_stream_is_chunk_safe():
    panel = _result("../reactome_state_panel.json")
    assert panel["pathway_count"] == 9
    assert panel["panel_gene_symbol_count"] == 2_228
    assert panel["manual_marker_additions"] == []
    assert {item["reactome_stable_id"] for item in panel["pathways"]} == {
        "R-HSA-1640170", "R-HSA-109581", "R-HSA-2559583",
        "R-HSA-2173791", "R-HSA-1234174", "R-HSA-73894",
        "R-HSA-381119", "R-HSA-168249", "R-HSA-913531",
    }
    assert all(len(item["participants_response_sha256"]) == 64 for item in panel["pathways"])
    payload = b"3 1 3\n17 1 2\n46 2 4\n49 8 1\n"
    parsed = np.vstack(list(_triplet_chunks(io.BytesIO(payload), block_bytes=7)))
    assert np.array_equal(parsed, np.asarray([[3, 1, 3], [17, 1, 2], [46, 2, 4], [49, 8, 1]]))
    contexts, metadata = [], []
    for pathway, count in (("common", 7), ("rare", 3)):
        for index in range(count):
            treatment = f"{pathway}-{index}"
            contexts.append(("A549", "rep1", "24", treatment, "10"))
            metadata.append(("FALSE", treatment, pathway, "level2", treatment, "target", "path"))
    contexts.append(("A549", "rep1", "24", "vehicle", "0"))
    metadata.append(("TRUE", "vehicle", "Vehicle", "Vehicle", "DMSO", "Vehicle", "Vehicle"))
    treatment_split, evaluable = _treatment_splits(contexts, metadata)
    assert treatment_split["vehicle"] == 4
    assert set(treatment_split[treatment] for treatment in treatment_split if treatment.startswith("common")) == {0, 1, 2, 3}
    assert set(treatment_split[treatment] for treatment in treatment_split if treatment.startswith("rare")) == {0}
    assert evaluable == ["common"]


def test_sciplex3_state_model_preserves_treatment_holdout_and_non_authority(tmp_path):
    contexts, metadata, splits, expression = [], [], [], []
    rng = np.random.default_rng(139944)
    for pathway_index, pathway in enumerate(("pathway_A", "pathway_B")):
        for split, suffix in enumerate(("train", "selection", "calibration", "test")):
            treatment = f"{pathway}_{suffix}"
            for cell_index, cell_type in enumerate(("A549", "MCF7")):
                for replicate in ("rep1", "rep2"):
                    contexts.append((cell_type, replicate, "24", treatment, "10"))
                    metadata.append(("FALSE", pathway))
                    splits.append(split)
                    signal = np.zeros(12, dtype=np.float32)
                    signal[pathway_index * 3:(pathway_index + 1) * 3] = 3.0
                    signal[6 + cell_index * 2:8 + cell_index * 2] = 1.5
                    expression.append(signal + rng.normal(0.0, 0.05, 12))
    for cell_type in ("A549", "MCF7"):
        for replicate in ("rep1", "rep2"):
            contexts.append((cell_type, replicate, "24", "vehicle", "0"))
            metadata.append(("TRUE", "Vehicle"))
            splits.append(4)
            expression.append(rng.normal(0.0, 0.05, 12))
    tensor = tmp_path / "sciplex_synthetic.npz"
    np.savez_compressed(
        tensor,
        log1p_cpm=np.asarray(expression, dtype=np.float32),
        context=np.asarray(contexts, dtype="U24"),
        context_field=np.asarray(
            ["cell_type", "replicate", "time_point", "treatment", "dose"]
        ),
        context_metadata=np.asarray(metadata, dtype="U24"),
        context_metadata_field=np.asarray(["vehicle", "pathway_level_1"]),
        gene_symbol=np.asarray([f"G{i}" for i in range(12)]),
        treatment_split=np.asarray(splits, dtype=np.uint8),
    )
    checkpoint = tmp_path / "sciplex_model.npz"
    report = train_sciplex3_state_model(
        tensor, checkpoint, epochs=3, hidden_width=8, selected_gene_count=8
    )
    assert report["split"]["cross_split_treatment_overlap_count"] == 0
    assert report["split"]["test_treatments_used_for_fitting_selection_or_calibration"] is False
    assert report["conformal_90"]["calibration_treatment_count"] == 2
    assert report["conformal_90"]["test_treatment_count"] == 2
    assert report["independent_lab_holdout"] is False
    assert report["aleph_authority"] == "none"
    assert checkpoint.is_file()


def test_sciplex3_gzip_verifier_counts_bytes_and_final_unterminated_row(tmp_path):
    import gzip
    terminated = tmp_path / "terminated.gz"
    unterminated = tmp_path / "unterminated.gz"
    with gzip.open(terminated, "wb") as handle:
        handle.write(b"1 1 2\n3 1 4\n")
    with gzip.open(unterminated, "wb") as handle:
        handle.write(b"1 1 2\n3 1 4")
    assert _verify_gzip_stream(terminated) == (12, 2)
    assert _verify_gzip_stream(unterminated) == (11, 2)


def test_range_download_resumes_inside_a_part_after_connection_reset(tmp_path, monkeypatch):
    calls = []

    class Response:
        def __init__(self, start, payload, fail):
            self.headers = {"Content-Range": f"bytes {start}-9/10"}
            self.payload = payload
            self.fail = fail
            self.read_count = 0

        def __enter__(self):
            return self

        def __exit__(self, *_):
            return False

        def read(self, _size):
            self.read_count += 1
            if self.read_count == 1:
                return self.payload
            if self.fail:
                raise ConnectionResetError("synthetic reset")
            return b""

    def fake_urlopen(request, timeout):
        assert timeout == 90
        header = request.get_header("Range")
        calls.append(header)
        start = int(header.removeprefix("bytes=").split("-")[0])
        if len(calls) == 1:
            return Response(start, b"abcd", True)
        return Response(start, b"efghij", False)

    monkeypatch.setattr("download_sciplex3_metadata.urllib.request.urlopen", fake_urlopen)
    monkeypatch.setattr("download_sciplex3_metadata.time.sleep", lambda _: None)
    destination = tmp_path / "part000"
    download_resumable_range("https://example.invalid/data", 0, 9, destination, 10)
    assert destination.read_bytes() == b"abcdefghij"
    assert calls == ["bytes=0-9", "bytes=4-9"]
    assert not (tmp_path / "part000.partial").exists()


def test_gse250041_senescence_is_direct_label_but_not_replicated_holdout():
    require(GSE250041_SOURCE)
    report = audit_gse250041_senescence_design(GSE250041_SOURCE)
    assert report["barcode_counts"] == {
        "untreated_proliferating": 8_664,
        "IR_10_Gy_day_10_senescent": 4_949,
    }
    assert report["feature_type_counts"] == {
        "Gene Expression": 36_601, "Antibody Capture": 8
    }
    assert report["feature_axis_identical_across_conditions"] is True
    assert report["independent_dataset_from_sciplex3"] is True
    assert report["independent_lab_from_sciplex3"] is True
    assert report["single_cells_are_not_independent_biological_replicates"] is True
    assert report["external_validation_eligible"] is False
    assert report["matrix_training_ready"] is False
    assert report["aleph_authority"] == "none"


def test_gse250041_matrix_header_verifier_reads_declared_dimensions(tmp_path):
    import gzip
    matrix = tmp_path / "matrix.mtx.gz"
    with gzip.open(matrix, "wt", encoding="ascii") as handle:
        handle.write("%%MatrixMarket matrix coordinate integer general\n")
        handle.write("% synthetic fixture\n36609 8664 2\n1 1 3\n36609 8664 1\n")
    assert _matrix_market_header(matrix) == (36_609, 8_664, 2)


def test_gse250041_audit_rejects_untyped_matrix_receipt(tmp_path):
    require(GSE250041_SOURCE)
    receipt = tmp_path / "receipt.json"
    receipt.write_text('{"schema":"wrong","accession":"GSE250041"}\n')
    with pytest.raises(ValueError, match="receipt schema"):
        audit_gse250041_senescence_design(GSE250041_SOURCE, receipt)


def test_gse250041_tensor_keeps_capture_confounding_explicit(tmp_path):
    import gzip
    source = tmp_path / "source"
    source.mkdir()
    adt = ["CD112", "HLA_A_B_C", "CD44", "CD54", "CD26", "CD49a", "CD73", "CD109"]
    feature_rows = [
        (f"ENSG{index:011d}", marker, "Gene Expression")
        for index, marker in enumerate(GSE250041_RNA_MARKERS, start=1)
    ] + [(name, f"{name}_TotalSeqB", "Antibody Capture") for name in adt]
    for prefix, cells in (("Proliferating", 2), ("Senescent", 1)):
        with gzip.open(source / f"GSE250041_{prefix}_features.tsv.gz", "wt") as handle:
            for row in feature_rows:
                handle.write("\t".join(row) + "\n")
        with gzip.open(source / f"GSE250041_{prefix}_barcodes.tsv.gz", "wt") as handle:
            for cell in range(cells):
                handle.write(f"cell-{cell}\n")
        entries = [
            (feature, cell, 1 + (feature + cell) % 3)
            for feature in range(1, len(feature_rows) + 1)
            for cell in range(1, cells + 1)
        ]
        with gzip.open(source / f"GSE250041_{prefix}_matrix.mtx.gz", "wt") as handle:
            handle.write("%%MatrixMarket matrix coordinate integer general\n")
            handle.write(f"{len(feature_rows)} {cells} {len(entries)}\n")
            for entry in entries:
                handle.write(" ".join(str(value) for value in entry) + "\n")
    output = tmp_path / "tensor.npz"
    report = build_gse250041_senescence_tensor(source, output, tmp_path / "report.json")
    tensor = np.load(output)
    assert tensor["X"].shape == (3, len(GSE250041_RNA_MARKERS) + 8)
    assert tensor["author_condition_label"].tolist() == [0, 0, 1]
    assert report["split_group"] == "capture_id"
    assert report["valid_dataset_split_count"] == 0
    assert report["classification_metric_authority"].startswith("none_")
    assert report["external_validation_eligible"] is False


def test_senscout_archive_member_guard_rejects_traversal_and_absolute_paths():
    assert senscout_safe_member("SenSCOUT-main/Data/features.csv") is True
    assert senscout_safe_member("../escape.csv") is False
    assert senscout_safe_member("/absolute/model.h5") is False
    assert senscout_safe_member("folder\\windows_escape.csv") is False


def test_senscout_group_split_and_transform_do_not_leak_fit_statistics():
    groups = np.asarray([f"image-{index // 3}" for index in range(90)])
    split = senscout_partition_groups(groups)
    assert set(split.tolist()) == {0, 1, 2}
    for name in set(groups.tolist()):
        assert len(set(split[groups == name].tolist())) == 1
    assert not (
        set(groups[split == 0]) & set(groups[split == 1])
        or set(groups[split == 0]) & set(groups[split == 2])
        or set(groups[split == 1]) & set(groups[split == 2])
    )

    raw = np.asarray([[-3.0, 0.0], [1.0, 2.0], [1e9, 1e9]])
    transformed, mean, scale = senscout_transform(raw, np.asarray([True, True, False]))
    expected = np.sign(raw[:2]) * np.log1p(np.abs(raw[:2]))
    assert np.array_equal(mean, expected.mean(axis=0))
    assert np.array_equal(scale, expected.std(axis=0))
    assert np.isfinite(transformed).all()


def test_senscout_model_is_same_provider_candidate_not_sweep_authority():
    report = _result("senscout_morphology_model_report.json")
    assert report["split"]["unit"] == "image_group"
    assert report["split"]["cross_split_image_group_overlap_count"] == 0
    assert report["split"]["cells_treated_as_independent_biological_replicates"] is False
    assert report["independent_lab_holdout"] is False
    assert report["production_eligible"] is False
    assert report["may_select_aleph_parameter"] is False
    assert report["direct_bgal_head"]["Bio2_image_group_prevalence"][
        "conformal_90_coverage"
    ] < 0.9
    assert "pending_external_dataset_lab_validation" in report["eventual_sweep_role"]

    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(), senscout_report=report
    )
    head = next(
        item for item in multitask["heads"]
        if item["head_id"] == "senscout_morphology_senescence_representation"
    )
    assert head["status"] == "training_only_same_provider_biorepeat_holdout"
    assert head["observable_evidence_eligible"] is False
    assert head["parameter_proposal_eligible"] is False
    assert head["external_holdout"] is False
    assert "senscout_morphology_senescence_representation" in multitask[
        "architecture"
    ]["trained_encoder_heads"]


def test_sencid_is_reference_expert_not_reproducible_external_authority():
    require(SENCID_SOURCE)
    report = audit_sencid_reference(SENCID_SOURCE)
    assert report["source"]["paper_reported_training_scope"] == {
        "samples": 602, "studies": 52, "cell_types": 30, "senescence_identities": 6,
    }
    assert report["sid_feature_counts"] == {
        "1": 19, "2": 13, "3": 186, "4": 87, "5": 83, "6": 224,
    }
    assert report["unique_sid_features"] == 445
    assert report["seneset_gene_count"] == 1290
    assert report["tracked_training_matrix_files"] == []
    assert report["study_level_training_manifest_present"] is False
    assert report["independent_retraining_possible_from_repository_alone"] is False
    assert all(
        item["executed_during_audit"] is False
        for item in report["model_artifacts"].values()
    )
    assert report["external_validation_eligible"] is False
    assert report["aleph_authority"] == "none"


def test_biad2515_apoptosis_endpoint_is_well_disjoint_and_terminal_leakage_free(tmp_path):
    require(BIAD2515_SOURCE)
    manifest = tmp_path / "observations.jsonl"
    tensor = tmp_path / "rows.npz"
    tensor_report = build_biad2515_apoptosis(
        BIAD2515_SOURCE, manifest, tensor
    )
    assert tensor_report["live_feature_count"] == 2_336
    assert tensor_report["terminal_feature_count_excluded_from_inputs"] == 530
    assert tensor_report["terminal_AnnexinV_features_present_in_X"] is False
    assert tensor_report["split_well_counts"] == {
        "fit": 10, "calibration": 10, "sealed_test": 10
    }
    checkpoint = tmp_path / "ridge.npz"
    report = evaluate_biad2515_apoptosis(tensor, checkpoint)
    assert report["split"]["well_overlap_count"] == 0
    assert report["split"]["test_used_for_fitting_selection_or_calibration"] is False
    assert report["sealed_same_provider_test"]["r2"] > 0.87
    assert report["sealed_same_provider_test"]["mae"] < 0.10
    assert report["baselines_on_same_sealed_test"]["image_model_beats_dose_only_mae"] is True
    assert report["conformal_90"]["sealed_test_coverage"] >= 0.90
    assert report["independent_lab_holdout"] is False
    assert report["production_eligible"] is False


def test_figshare_cell_death_design_is_verified_and_not_endpoint_conflated(tmp_path):
    require(FIGSHARE_CELL_DEATH_SOURCE)
    receipt = acquire_figshare_cell_death(FIGSHARE_CELL_DEATH_SOURCE)
    assert receipt["license_id"] == "CC-BY-4.0"
    assert receipt["plate_count"] == 8
    assert receipt["metadata_site_rows"] == 21_271
    assert receipt["qc_site_rows"] == 15_430
    assert len(receipt["programmed_cell_death_labels"]) == 6
    assert {
        item["unique_plate_wells"]
        for item in receipt["author_split_summary"].values()
    } == {887, 894, 913}
    assert receipt["independent_from_biad2515_provider"] is True
    assert receipt["same_endpoint_as_biad2515_annexin_v"] is False
    assert receipt["feature_archive_ready_for_training"] is False
    assert receipt["aleph_authority"] == "none"
    assert receipt["may_select_aleph_parameter"] is False


def test_figshare_cell_death_model_split_is_compound_disjoint():
    require(FIGSHARE_CELL_DEATH_SOURCE)
    report = audit_figshare_cell_death_design(FIGSHARE_CELL_DEATH_SOURCE)
    assert report["compound_count"] == 55
    assert report["split_compound_counts"] == {
        "train": 40, "selection": 4, "calibration": 5, "sealed_test": 6
    }
    assert report["cross_split_compound_overlap_count"] == 0
    assert report["sealed_test_has_all_six_labels"] is True
    assert report["selection_missing_rare_labels"] == [
        "necroptosis inducer", "pyroptosis inducer"
    ]
    assert report["calibration_missing_rare_labels"] == ["pyroptosis inducer"]
    assert report["rare_label_conditional_calibration_authority"] is False
    author = report["author_well_split_leakage_audit"]
    assert all(item["well_disjoint_but_not_compound_disjoint"] for item in author.values())
    assert min(item["compounds_crossing_author_roles"] for item in author.values()) >= 45
    assert report["aleph_authority"] == "none"


def test_bbbc048_raw_image_head_is_disjoint_calibrated_and_training_only():
    tensor_receipt = _result("bbbc048_image_tensor_receipt.json")
    report = _result("bbbc048_mlp_report.json")
    assert tensor_receipt["tensor_shape"] == [32_266, 3, 16, 16]
    assert tensor_receipt["split_group_overlap"] == 0
    assert set(tensor_receipt["split_class_counts"]) == {
        "train", "selection", "calibration", "test"
    }
    assert report["test"]["macro_f1"] > 0.60
    assert report["test"]["balanced_accuracy"] > 0.65
    assert report["conformal_90"]["test_marginal_coverage"] >= 0.89
    assert report["conformal_90"]["rare_class_conditional_coverage_authority"] is False
    assert report["independent_lab_holdout"] is False
    assert report["production_eligible"] is False


def test_bbbc053_mitochondrial_stress_is_raw_but_not_external_authority():
    report = _result("bbbc053_mito_stress_report.json")
    assert report["condition_counts"] == {"DMSO": 29, "FCCP": 29}
    assert report["compact_image_shape"] == [58, 64, 64]
    assert report["morphology_feature_count"] == 13
    assert report["biological_replicates_reported"] is False
    assert report["independent_lab_holdout"] is False
    assert report["status"] == "same_provider_training_only"


def test_hpa_raw_subset_is_split_preserving_and_diverse():
    tensor = HERE / "results" / "hpa_if_embedding_rows.npz"
    require(tensor)
    manifest = select_hpa_raw_image_subset(tensor)
    committed = _result("hpa_raw_if_subset_manifest.json")
    assert manifest == committed
    assert manifest["image_count"] == 512
    assert manifest["download_file_count"] == 2_048
    assert manifest["cell_line_count"] >= 25
    assert manifest["author_location_label_count"] == 35
    assert manifest["gene_antibody_component_count"] == 512
    assert manifest["gene_antibody_component_overlap_across_splits"] == 0
    assert manifest["independent_lab_holdout"] is False


def test_bbbc054_author_labels_and_brightfield_features_are_exact_cell_grouped():
    require(BBBC054_REPLICATE1, BBBC054_ANNOTATION)
    observations, tensor = build_bbbc054_microglia(
        BBBC054_REPLICATE1, BBBC054_ANNOTATION
    )
    assert len(observations) == 58_186
    assert tensor["X"].shape == (58_186, len(BBBC054_FEATURE_NAMES))
    assert np.isfinite(tensor["X"]).all()
    assert set(tensor["label"].astype(str)) == {"amoeboid", "ramified", "round"}
    assert set(tensor["field"]) == set(range(60))
    assert len(set(tensor["sample_id"].astype(str))) == 58_186
    assert {row.biological_replicate for row in observations} == {"replicate_1"}
    assert {row.aleph_authority for row in observations} == {"none"}


def test_bbbc054_tiff_decoder_preserves_official_16bit_field():
    require(BBBC054_REPLICATE1)
    import zipfile
    with zipfile.ZipFile(BBBC054_REPLICATE1) as archive:
        image = decode_bbbc054_tiff(archive.read("Replicate 1/IMG_20x_0.tif"))
    assert image.shape == (1080, 1280)
    assert image.dtype == np.float32
    assert float(image.std()) > 100.0


def test_bbbc054_microglia_head_reports_failure_without_false_external_authority():
    tensor = HERE / "results" / "bbbc054_microglia_rows.npz"
    require(tensor)
    report = evaluate_bbbc054_microglia(tensor)
    assert report == _result("bbbc054_microglia_report.json")
    assert report["split"]["field_overlap_count"] == 0
    assert report["split"]["pristine_confirmation"] is False
    assert report["sealed_field_test_metrics"]["macro_f1"] < 0.50
    assert report["split_conformal_diagnostic"]["mean_prediction_set_size"] > 2.0
    assert report["independent_biological_replicate_holdout"] is False
    assert report["may_select_aleph_parameter"] is False


def test_temporal_cnn_is_reproducible_but_rejected_by_external_holdout(tmp_path):
    train_manifest = HERE / "results" / "mechanotransduction_observations.jsonl"
    external_manifest = HERE / "results" / "eahy_calcium_observations.jsonl"
    require(train_manifest, external_manifest)
    checkpoint = tmp_path / "calcium_temporal_cnn_checkpoint.npz"
    report = train_calcium_temporal_cnn(
        train_manifest, external_manifest, checkpoint
    )
    committed = _result("calcium_temporal_cnn_report.json")
    assert report == committed
    assert report["architecture"]["parameter_count"] == 114
    assert report["architecture"]["filters_are_trainable"] is True
    assert report["selected_validation_metrics"]["macro_f1"] >= 0.94
    assert report["external_retrospective_metrics"]["auroc"] >= 0.77
    assert report["external_retrospective_metrics"]["macro_f1"] < 0.70
    assert report["external_embedding_ood"]["external_vs_validation_auroc"] < 0.70
    assert report["external_embedding_ood"]["ood_detected"] is False
    label_shift = report["unlabelled_label_shift_diagnostic"]
    assert label_shift["actual_prior_used_for_fitting_or_selection"] is False
    assert label_shift["estimated_external_prior_without_labels"][1] < 0.50
    assert label_shift["retrospective_actual_external_prior"][1] > 0.80
    assert label_shift["adaptation_accepted"] is False
    assert report["split"]["pristine_external_confirmation"] is False
    assert report["task_holdout_passed"] is False
    assert report["architecture_accepted"] is False
    assert report["aleph_authority"] == "none"
    assert hashlib.sha256(checkpoint.read_bytes()).hexdigest() == committed[
        "checkpoint"
    ]["sha256"]


def test_temporal_cnn_checkpoint_has_read_only_probability_and_ood_inference():
    checkpoint = HERE / "results" / "calcium_temporal_cnn_checkpoint.npz"
    manifest = HERE / "results" / "eahy_calcium_observations.jsonl"
    require(checkpoint, manifest)
    model, info = load_calcium_temporal_checkpoint(checkpoint)
    assert info["schema"] == "aleph.outer_library.calcium_temporal_cnn.v1"
    assert info["filter_count"] == 8
    assert model["class_centroids"].shape == (2, 16)
    allowed = (
        "LooPINS (+GsMTx4)", "LooPINS (+Mag)",
        "CaPINS (+GsMTx4)", "CaPINS (+Mag)",
    )
    rows = predict_calcium_temporal_cnn(
        checkpoint, manifest, allowed_conditions=allowed
    )
    assert len(rows) == 311
    probability = np.asarray([
        [row["probability_non_inhibited_reference"],
         row["probability_Piezo1_inhibited"]]
        for row in rows
    ])
    labels = np.asarray([
        int("+GsMTx4" in row["condition"]) for row in rows
    ])
    assert np.allclose(probability.sum(axis=1), 1.0, rtol=0.0, atol=1e-15)
    assert calcium_binary_metrics(probability, labels) == _result(
        "calcium_temporal_cnn_report.json"
    )["external_retrospective_metrics"]
    assert all(np.isfinite(row["embedding_ood_distance"]) for row in rows)
    assert {row["aleph_authority"] for row in rows} == {"none"}
    assert {row["may_select_aleph_parameter"] for row in rows} == {False}
    assert not any("parameter" in row and row["may_select_aleph_parameter"]
                   for row in rows)


def test_multitask_outer_model_rejects_split_leakage_and_claimed_authority():
    inputs = _multitask_inputs()
    inputs[4]["protocol"]["external_test_lab"] = inputs[4]["protocol"]["train_lab"]
    with pytest.raises(ValueError, match="must be disjoint"):
        evaluate_multitask_outer_model(*inputs)

    inputs = _multitask_inputs()
    inputs[6]["aleph_parameter_proposal_eligible"] = True
    with pytest.raises(ValueError, match="direct Aleph parameter"):
        evaluate_multitask_outer_model(*inputs)


def test_multitask_model_schema_and_sweep_gate_preserve_non_authority():
    schema = _result("../multitask_outer_model.schema.json")
    assert schema["properties"]["schema"]["const"] == MULTITASK_MODEL_SCHEMA
    assert schema["properties"]["aleph_authority"]["const"] == "none"
    assert schema["properties"]["sweep_interface"]["properties"][
        "may_select_aleph_parameter"
    ]["const"] is False
    multitask = evaluate_multitask_outer_model(*_multitask_inputs())
    gate = evaluate(
        _result("calcium_baseline_report.json"),
        multitask_model_report=multitask,
    )
    assert gate["checks"]["multitask_model_contract_is_non_authoritative"] is True
    assert gate["checks"][
        "multitask_model_admits_only_externally_validated_evidence_heads"
    ] is True
    assert gate["multitask_model_audit"]["production_eligible"] is False
    assert gate["multitask_model_audit"]["may_select_aleph_parameter"] is False
    assert gate["status"] == "refused"

    tampered = evaluate_multitask_outer_model(*_multitask_inputs())
    calcium = next(
        head for head in tampered["heads"]
        if head["head_id"] == "calcium_cell_state"
    )
    calcium["observable_evidence_eligible"] = True
    tampered["sweep_interface"]["admissible_head_ids"].append(
        "calcium_cell_state"
    )
    tampered_gate = evaluate(
        _result("calcium_baseline_report.json"),
        multitask_model_report=tampered,
    )
    assert tampered_gate["checks"][
        "multitask_model_admits_only_externally_validated_evidence_heads"
    ] is False
    assert tampered_gate["status"] == "refused"


def test_pbmc_preregistered_external_lab_gate_fails_only_transfer_uncertainty():
    tensor = _result("pbmc_cross_lab_tensor_report.json")
    model = _result("pbmc_cross_lab_model_report.json")
    direction = _result("pbmc_ifn_direction_report.json")
    assert tensor["training"]["cells"] == 24_413
    assert tensor["training"]["common_frozen_features"] == 508
    assert tensor["external"]["matched_and_retained_cells"] == 22_594
    assert tensor["secondary_benchmark_label_audit"]["decision"] == (
        "rejected_as_classifier_targets"
    )
    passes = model["cell_type"]["frozen_endpoint_passes"]
    assert passes == {
        "external_conformal_coverage": False,
        "external_ece": True,
        "external_macro_f1": True,
        "sealed_internal_macro_f1": True,
        "semantic_ood_auroc": True,
    }
    assert model["overall_external_context_gate_pass"] is False
    assert model["authority"]["may_emit_numeric_sweep_range"] is False
    assert direction["all_three_paired_donor_differences_positive"] is True
    assert direction["numeric_scale_transfer_permitted"] is False


def test_pbmc_head_remains_pre_sweep_ranking_only_after_external_failure():
    report = evaluate_multitask_outer_model(
        *_multitask_inputs(),
        pbmc_cross_lab_report=_result("pbmc_cross_lab_model_report.json"),
        pbmc_ifn_direction_report=_result("pbmc_ifn_direction_report.json"),
    )
    head = next(
        item for item in report["heads"]
        if item["head_id"] == "PBMC_cell_type_IFN_pre_sweep_context"
    )
    assert head["status"] == "blocked_external_uncertainty"
    assert head["external_holdout"] is True
    assert head["pristine_confirmation"] is True
    assert head["observable_evidence_eligible"] is False
    assert head["parameter_proposal_eligible"] is False
    assert report["readiness"]["general_cell_type_state_model_ready"] is False
    assert "PBMC_cell_type_IFN_pre_sweep_context" in report["readiness"][
        "blocked_external_classifier_heads"
    ]


def test_pbmc_disjoint_calibration_and_final_lab_keep_uncertainty_blocked():
    calibration = _result("pbmc_gse164378_calibration_report.json")
    confirmation = _result("wilk_pbmc_confirmation_report.json")
    assert calibration["classifier_weights_bit_identical"] is True
    assert calibration["source_parameter_sha256"] == calibration[
        "output_parameter_sha256"
    ]
    assert len(calibration["calibration"]["distinct_donors"]) == 8
    assert confirmation["single_use_confirmation_consumed"] is True
    assert confirmation["may_be_reused_as_pristine"] is False
    assert confirmation["population"]["healthy_donors"] == [
        "H1", "H2", "H3", "H4", "H5", "H6"
    ]
    assert confirmation["primary_metrics"]["macro_f1"] >= 0.70
    assert confirmation["primary_metrics"]["ece_15_bin"] <= 0.15
    assert confirmation["semantic_ood"]["auroc"] >= 0.70
    assert confirmation["primary_metrics"]["conformal"][
        "marginal_coverage"
    ] < 0.90
    assert confirmation["all_frozen_endpoints_passed"] is False

    report = evaluate_multitask_outer_model(
        *_multitask_inputs(),
        pbmc_cross_lab_report=_result("pbmc_cross_lab_model_report.json"),
        pbmc_ifn_direction_report=_result("pbmc_ifn_direction_report.json"),
        pbmc_final_confirmation_report=confirmation,
    )
    head = next(
        item for item in report["heads"]
        if item["head_id"] == "PBMC_cell_type_IFN_pre_sweep_context"
    )
    assert head["status"] == "blocked_external_uncertainty"
    assert head["uncertainty_eligible"] is False
    assert head["metrics"]["postcalibration_single_use_confirmation"][
        "status"
    ] == "failed_and_permanently_retired"


def test_multitask_contract_accepts_trained_sciplex_only_as_same_provider_head():
    inputs = _multitask_inputs()
    design = json.loads(json.dumps(inputs[16]))
    design["expression_matrix_acquired"] = True
    design["training_ready"] = True
    design["reason_not_training_ready"] = None
    design["matrix_evidence"] = {
        "size_bytes": 3_096_224_606,
        "sha256": "a" * 64,
        "decompressed_line_count": 2_000_000,
        "gzip_stream_verified": True,
    }
    inputs[16] = design
    state = {
        "schema": "aleph.outer_library.sciplex3_state_model.v1",
        "model": {"family": "synthetic_contract_test"},
        "split": {
            "cross_split_treatment_overlap_count": 0,
            "test_treatments_used_for_fitting_selection_or_calibration": False,
        },
        "sealed_unseen_treatment_test": {"treatment_level": {"macro_f1": 0.5}},
        "conformal_90": {"test_marginal_coverage": 0.9},
        "vehicle_ood_diagnostic": {"authority": "diagnostic_only"},
        "replicate_consistency": {"same_lab_technical_or_biological_status": "same_provider"},
        "independent_dataset_holdout": False,
        "independent_lab_holdout": False,
        "production_eligible": False,
        "aleph_authority": "none",
        "limitations": ["same-provider contract test"],
    }
    report = evaluate_multitask_outer_model(*inputs, sciplex_state_report=state)
    head = next(
        item for item in report["heads"]
        if item["head_id"] == "sciplex_perturbation_state_representation"
    )
    assert head["status"] == "training_only_unseen_compound_holdout_same_provider"
    assert head["external_holdout"] is False
    assert head["observable_evidence_eligible"] is False
    assert head["parameter_proposal_eligible"] is False
    assert "sciplex_perturbation_state_representation" in report["architecture"][
        "trained_encoder_heads"
    ]
    assert report["readiness"]["general_cell_type_state_model_ready"] is False


def test_theory_runtime_roundtrips_and_rejects_invalid_domains():
    force = tether_force(0.082, 54.0)["tether_force_pn"]
    recovered = tether_apparent_tension(force, 0.082)[
        "apparent_tension_pn_per_um"
    ]
    assert recovered == pytest.approx(54.0, rel=1e-15)
    assert standard_linear_solid_relaxation(
        [0.0, 2.0], 5.0, 1.0, 2.0
    )["response"][0] == 5.0
    assert persistent_random_walk_msd_2d(
        [0.0], 0.5, 3.0
    )["mean_squared_displacement_um2"][0] == 0.0
    with pytest.raises(TheoryRefusal):
        tether_force(0.0, 1.0)


def test_external_tether_inference_preserves_exact_ambiguity():
    result = infer_external({
        "modality": "AFM_membrane_tether_force",
        "value": 18.697398,
        "unit": "pN",
        "bending_rigidity_pn_um": 0.082,
        "provenance": {"dataset": "worked-example"},
    })
    assert result["refused"] is False
    assert result["may_emit_Aleph_parameter"] is False
    assert result["ambiguity_sets"][0]["may_choose_one_point"] is False
    assert set(result["ambiguity_sets"][0]["bounds_under_nonnegative_components"]) == {
        "sigma_bilayer_pn_per_um", "attachment_energy_density_pn_per_um"
    }


def test_latent_promotion_is_temporary_pre_sweep_gate_not_terminal_redesign():
    phase = evaluate_external_inference()
    assert phase["objective_change"] is False
    assert phase["execution_order_change"] is True
    assert "Aleph_parameter_sweeps" in phase["final_program_objective"]

    inference = infer_external({
        "modality": "AFM_membrane_tether_force", "value": 18.7, "unit": "pN",
        "bending_rigidity_pn_um": 0.082, "provenance": {"dataset": "example"},
    })
    blocked = evaluate_latent_promotion(
        inference,
        external_dataset_holdout_passed=False,
        external_lab_holdout_passed=False,
        uncertainty_calibrated=False,
    )
    assert blocked["stage"] == "external_training_or_validation_only"
    assert blocked["bounded_sweep_proposal_eligible"] is False
    assert blocked["numeric_parameter_range"] is None

    resolved = json.loads(json.dumps(inference))
    for ambiguity in resolved["ambiguity_sets"]:
        ambiguity["resolved_by_independent_observable"] = True
    mapping = {
        "operator_status": "native_validated",
        "operator_id": "aleph.observe.test_only",
        "aleph_observable": "AFM_membrane_tether_force_pN",
        "external_unit": "pN", "aleph_unit": "pN",
        "unit_conversion_verified": True, "geometry_matched": True,
        "validation_report_sha256": "a" * 64,
    }
    promoted = evaluate_latent_promotion(
        resolved,
        external_dataset_holdout_passed=True,
        external_lab_holdout_passed=True,
        uncertainty_calibrated=True,
        operator_mapping=mapping,
        forward_sensitivity_passed=True,
    )
    assert promoted["stage"] == "eligible_for_bounded_sweep_proposal"
    assert promoted["bounded_sweep_proposal_eligible"] is True
    assert promoted["may_select_Aleph_parameter"] is False


def test_observation_operator_registry_blocks_semantic_shortcuts():
    registry = json.loads((HERE / "observation_operator_registry.json").read_text(
        encoding="utf-8"
    ))
    report = audit_observation_operators(registry, HERE.parents[3])
    assert report["passed"] is True
    assert report["promotable_external_observables"] == []
    assert set(report["unavailable_external_observables"]) >= {
        "CDH1_RNA_abundance", "VIM_RNA_abundance",
        "E_cadherin_protein_abundance", "Vimentin_protein_abundance",
    }
    assert report["candidate_native_external_observables"] == [
        "effective_laplace_surface_tension",
        "tissue_traction_force_per_cell",
    ]
    by_name = {item["external_observable"]: item for item in report["entries"]}
    assert by_name["VIM_RNA_abundance"]["operator_id"] is None
    assert by_name["tissue_traction_force_per_cell"]["source_symbol_present"] is True

    tampered = json.loads(json.dumps(registry))
    tampered["entries"][0]["status"] = "native_validated"
    failed = audit_observation_operators(tampered, HERE.parents[3])
    assert failed["passed"] is False
    assert "validated_statuses_have_complete_validation" in failed["failed_checks"]


def test_latent_promotion_rejects_unvalidated_operator_name():
    inference = infer_external({
        "modality": "AFM_membrane_tether_force", "value": 18.7, "unit": "pN",
        "bending_rigidity_pn_um": 0.082, "provenance": {"dataset": "example"},
    })
    for ambiguity in inference["ambiguity_sets"]:
        ambiguity["resolved_by_independent_observable"] = True
    mapping = {
        "operator_status": "candidate_native",
        "operator_id": "aleph.vertical.some_similarly_named_callable",
        "aleph_observable": "tether_force", "external_unit": "pN",
        "aleph_unit": "pN", "unit_conversion_verified": True,
        "geometry_matched": True, "validation_report_sha256": "a" * 64,
    }
    result = evaluate_latent_promotion(
        inference, external_dataset_holdout_passed=True,
        external_lab_holdout_passed=True, uncertainty_calibrated=True,
        operator_mapping=mapping, forward_sensitivity_passed=True,
    )
    assert result["operator_constraint_ready"] is False
    assert result["bounded_sweep_proposal_eligible"] is False
    assert "Aleph_observation_operator_mapping_complete" in result["failed_checks"]


def test_sweep_gate_resolves_preregistered_emt_markers_as_unavailable():
    registry = json.loads((HERE / "observation_operator_registry.json").read_text(
        encoding="utf-8"
    ))
    operator_report = audit_observation_operators(registry, HERE.parents[3])
    gate = evaluate(
        _result("calcium_baseline_report.json"),
        gse325309_emt_report=_result("gse325309_pristine_evaluation.json"),
        observation_operator_report=operator_report,
    )
    assert gate["checks"]["observation_operator_registry_passed"] is True
    assert gate["checks"]["EMT_marker_operators_are_explicitly_unavailable"] is True
    request = next(
        item for item in gate["calibration_requests"]
        if item["parameter"] == "mechanism_family:preregistered_EMT_direction_confirmation"
    )
    assert request["current_operator_resolution"]["CDH1_RNA_abundance"][
        "status"
    ] == "unavailable"
    assert request["numeric_range"] is None
    assert request["may_emit_sweep_axis"] is False


def test_figshare_programmed_cell_death_is_compound_grouped_training_only():
    report = _result("figshare_cell_death_model_report.json")
    assert report["split"]["unit"] == "compound_identity"
    assert report["split"]["cross_split_compound_overlap_count"] == 0
    assert report["split"]["well_rows_are_not_independent_training_units"] is True
    assert report["sealed_unseen_compound_test"]["metrics"]["compound_count"] == 6
    assert report["sealed_unseen_compound_test"]["conformal_90"]["coverage"] < 0.9
    assert report["observable_evidence_eligible"] is False
    assert report["may_select_aleph_parameter"] is False
    assert figshare_archive_path_safe("features/aggregate.parquet") is True
    assert figshare_archive_path_safe("../escape.parquet") is False
    assert figshare_archive_path_safe("/absolute.parquet") is False
    assert figshare_archive_path_safe("features\\escape.parquet") is False

    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(), figshare_cell_death_report=report
    )
    head = next(
        item for item in multitask["heads"]
        if item["head_id"] == "programmed_cell_death_mechanism_representation"
    )
    assert head["status"] == "training_only_unseen_compound_same_provider_rejected"
    assert head["external_holdout"] is False
    assert head["parameter_proposal_eligible"] is False


def test_gse301164_library_size_uses_all_genes_not_only_frozen_markers(tmp_path):
    path = tmp_path / "counts.tsv.gz"
    header = [
        "Geneid", "Chr", "Start", "End", "Strand", "Length", "gene_name",
        "gene_type", "sample_1", "sample_2",
    ]
    rows = [
        ["ENSG1", "chr1", "1", "2", "+", "2", "MARK", "protein_coding", "10", "10"],
        ["ENSG2", "chr1", "3", "4", "+", "2", "HOUSE", "protein_coding", "90", "990"],
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    counts, library, samples = load_gse301164_counts(path, np.asarray(["MARK"]))
    assert np.array_equal(samples, np.asarray(["sample_1", "sample_2"]))
    assert np.array_equal(counts[:, 0], np.asarray([10.0, 10.0]))
    assert np.array_equal(library, np.asarray([100.0, 1000.0]))


def test_external_senescence_grants_direction_only_and_refuses_absolute_state():
    model = _result("external_senescence_model_report.json")
    confirmation = _result("gse301164_postfreeze_confirmation.json")
    assert model["split"]["dataset_overlap_count"] == 0
    assert model["split"]["lab_overlap_count"] == 0
    assert model["sealed_external_dataset_lab_test"]["conformal_90"]["coverage"] == 0.125
    assert model["sealed_external_dataset_lab_test"]["ood_rejection_fraction"] == 1.0
    assert confirmation["checkpoint_frozen_before_dataset_acquisition"] is True
    assert confirmation["model_fit_selection_calibration_or_retraining_on_GSE301164"] is False
    assert confirmation["combined"] == {
        "directional_effect_evidence_eligible": True,
        "exact_one_sided_sign_test_p": 0.015625,
        "pair_count": 6,
        "positive_pair_count": 6,
    }
    assert confirmation["ood"]["rejection_fraction"] == 1.0
    assert confirmation["absolute_state_classifier_promoted"] is False
    assert confirmation["may_select_aleph_parameter"] is False

    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(),
        external_senescence_report=model,
        senescence_confirmation_report=confirmation,
    )
    head = next(
        item for item in multitask["heads"]
        if item["head_id"] == "senescence_transcriptomic_direction"
    )
    assert head["status"] == "validated_external_direction_only_absolute_state_OOD_blocked"
    assert head["observable_evidence_eligible"] is True
    assert head["uncertainty_eligible"] is False
    assert head["split_contract"]["absolute_sample_state_use_refused"] is True
    assert head["parameter_proposal_eligible"] is False
    assert multitask["sweep_interface"]["may_emit_numeric_parameter_range"] is False
    assert multitask["sweep_interface"]["may_select_aleph_parameter"] is False


def test_zenodo_tfm_tether_is_a_future_sweep_constraint_not_a_replacement():
    manifest = HERE / "results" / "zenodo_7432971_tfm_tether_observations.jsonl"
    report = evaluate_zenodo_tfm_tension(manifest)
    assert report["purpose"].startswith(
        "accumulate externally grounded observables that can reduce a future Aleph"
    )
    assert report["design_audit"]["cells_treated_as_biological_replicates"] is False
    assert report["design_audit"]["neuron_stiffness_partly_confound_with_acquisition_date"] is True
    fibroblast = report["descriptive_results"][
        "fibroblast_mean_traction_by_stiffness"
    ]
    assert fibroblast["PAA_shear_modulus_1_kPa"]["cell_measurement_count"] == 30
    assert fibroblast["PAA_shear_modulus_10_kPa"]["median"] > fibroblast[
        "PAA_shear_modulus_1_kPa"
    ]["median"]
    assert set(
        item["direction"] for item in report["descriptive_results"][
            "fibroblast_stiffness_direction_by_batch"
        ].values()
    ) == {"increases"}
    assert report["observation_operator_audit"][
        "time_averaged_mean_traction_stress"
    ]["unit_identity"] == "1 pN/um2 == 1 Pa"
    assert report["promotion"]["current_executable_Aleph_sweep_available"] is False
    assert report["promotion"]["future_sweep_range_reduction_intended"] is True
    assert report["promotion"]["observable_evidence_eligible"] is False
    assert report["promotion"]["numeric_range"] is None
    assert report["may_select_aleph_parameter"] is False

    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(), tfm_tether_report=report
    )
    head = next(
        item for item in multitask["heads"]
        if item["head_id"] == "TFM_tether_pre_sweep_observation_constraint"
    )
    assert head["status"] == "training_only_operator_and_replication_blocked"
    assert head["parameter_proposal_eligible"] is False
    assert head["head_id"] not in multitask["sweep_interface"]["admissible_head_ids"]
    assert multitask["sweep_interface"]["may_emit_numeric_parameter_range"] is False
    gate = evaluate(
        _result("calcium_baseline_report.json"),
        multitask_model_report=multitask,
    )
    assert gate["status"] == "refused"

    assert zenodo_tfm_archive_path_safe("data_zenodo/Figure Data/Figure2_Data.xlsx")
    assert not zenodo_tfm_archive_path_safe("../escape.xlsx")
    assert not zenodo_tfm_archive_path_safe("/absolute.xlsx")
    assert not zenodo_tfm_archive_path_safe("folder\\escape.xlsx")


def test_zenodo_tfm_fret_preserves_cells_and_nested_focal_adhesions():
    require(ZENODO_TFM_FRET_SOURCE)
    rows, manifest_report = build_zenodo_tfm_fret(
        ZENODO_TFM_FRET_SOURCE,
        HERE / "results" / "zenodo_14692589_acquisition.json",
    )
    assert len(rows) == 847
    assert len({row.sample_id for row in rows}) == 358
    assert len({row.observation_id for row in rows}) == 847
    assert manifest_report["cell_level_TFM_cells"] == 89
    assert manifest_report["cell_level_FLIM_FRET_cells"] == 266
    assert manifest_report["paired_focal_adhesions"] == 246
    assert manifest_report["cell_to_replicate_mapping_released"] is False
    assert all("mapping_not_released" in row.biological_replicate for row in rows)
    paired = [row for row in rows if row.modality == "paired_FA_TFM_FLIM_FRET_derived"]
    assert len(paired) == 492
    assert {row.unit for row in paired} == {"author_normalized_fraction"}
    assert {row.aleph_authority for row in rows} == {"none"}
    assert zenodo_tfm_fret_path_safe("Source-Data/SD_Fig2c.prism")
    assert not zenodo_tfm_fret_path_safe("../escape.prism")
    assert not zenodo_tfm_fret_path_safe("folder\\escape.prism")


def test_zenodo_tfm_fret_adds_cross_lab_direction_without_numeric_promotion():
    report = evaluate_zenodo_tfm_fret(
        HERE / "results" / "zenodo_14692589_tfm_fret_observations.jsonl",
        HERE / "results" / "zenodo_14692589_manifest_report.json",
        HERE / "results" / "zenodo_7432971_tfm_tether_observations.jsonl",
    )
    traction = report["stiffness_directions"]["TFM"]
    assert traction["median_13_minus_4.5_Pa"] == pytest.approx(57.0870205)
    assert traction["increases_with_stiffness"] is True
    assert report["stiffness_directions"]["VinTS_FLIM_FRET"][
        "FRET_decreases_and_vinculin_tension_increases"
    ] is True
    assert report["author_correlation_categories_reproduced"] is True
    correlations = report["paired_FA_cell_correlations"]
    assert correlations["directly_correlated_cell"]["FRET_on_traction_slope"] < 0
    assert correlations["inversely_correlated_cell"]["FRET_on_traction_slope"] > 0
    assert abs(correlations["uncorrelated_cell"]["pearson_correlation"]) < 0.1
    assert report["cross_lab_same_observable_direction"][
        "independent_lab_direction_reproduced"
    ] is True
    assert report["design_audit"]["independent_lab_inferential_holdout"] is False
    assert report["pre_sweep_constraint"]["numeric_range"] is None
    assert report["pre_sweep_constraint"]["may_emit_numeric_parameter_range"] is False

    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(), tfm_fret_report=report
    )
    head = next(
        item for item in multitask["heads"]
        if item["head_id"] == "cross_lab_TFM_FRET_load_transfer_pre_sweep_constraint"
    )
    assert head["observable_evidence_eligible"] is False
    assert head["parameter_proposal_eligible"] is False
    assert head["head_id"] not in multitask["sweep_interface"]["admissible_head_ids"]
    gate = evaluate(
        _result("calcium_baseline_report.json"),
        multitask_model_report=multitask,
        tfm_fret_report=report,
    )
    assert gate["checks"]["tfm_fret_constraint_is_non_authoritative"] is True
    request = next(
        item for item in gate["calibration_requests"]
        if item["parameter"] == "mechanism_family:focal_adhesion_load_transfer"
    )
    assert request["numeric_range"] is None
    assert request["may_emit_sweep_axis"] is False


def test_figshare_exif_emt_raw_if_is_grouped_pre_sweep_training_only(tmp_path):
    receipt = HERE / "results" / "figshare_exif_emt_acquisition.json"
    require(FIGSHARE_EXIF_EMT_RAW, receipt)
    manifest = tmp_path / "exif_emt_observations.jsonl"
    tensor = tmp_path / "exif_emt_rows.npz"
    manifest_report = build_figshare_exif_emt(
        FIGSHARE_EXIF_EMT_RAW, receipt, manifest, tensor
    )
    assert manifest_report["counts"] == {
        "field_samples": 48,
        "well_split_groups": 48,
        "raw_channel_files": 240,
        "marker_panels": 8,
        "conditions": 3,
    }
    assert manifest_report["biological_replicates_reported"] is False
    assert manifest_report["aleph_authority"] == "none"
    with np.load(tensor, allow_pickle=False) as data:
        assert data["X"].shape == (48, 60)
        assert data["compact_image"].shape == (48, 5, 64, 64)
        assert len(set(data["split_group"].astype(str))) == 48

    report_path = tmp_path / "manifest_report.json"
    report_path.write_text(json.dumps(manifest_report), encoding="utf-8")
    report = evaluate_figshare_exif_emt(tensor, report_path)
    assert report["state_axis_status"] == (
        "author_treatment_labelled_raw_IF_training_only"
    )
    assert report["representation_training_eligible"] is True
    assert report["observable_evidence_eligible"] is False
    assert report["independent_lab_holdout"] is False
    assert report["pre_sweep_role"]["numeric_range"] is None
    assert report["pre_sweep_role"]["may_emit_sweep_axis"] is False
    assert report["may_select_aleph_parameter"] is False

    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(), exif_emt_report=report
    )
    head = next(
        item for item in multitask["heads"]
        if item["head_id"] == "raw_IF_EMT_condition_representation"
    )
    assert head["observable_evidence_eligible"] is False
    assert head["parameter_proposal_eligible"] is False
    assert head["head_id"] not in multitask["sweep_interface"]["admissible_head_ids"]
    gate = evaluate(
        _result("calcium_baseline_report.json"),
        multitask_model_report=multitask,
        exif_emt_report=report,
    )
    assert gate["checks"]["exif_emt_constraint_is_non_authoritative"] is True
    audit = next(
        item for item in gate["representation_training_audits"]
        if item["dataset"] == "figshare-26500210-exif-a549-emt-4plex"
    )
    assert audit["status"] == "admitted_training_only"
    assert audit["may_emit_numeric_parameter_range"] is False
    assert audit["may_select_aleph_parameter"] is False
    assert figshare_exif_path_safe("Images/raw/ExIF_EMT/A01_0000_DAPI.tif")
    assert not figshare_exif_path_safe("../escape.tif")
    assert not figshare_exif_path_safe("folder\\escape.tif")


def test_mtrack_emt_is_one_conflicting_temporal_diagnostic_not_a_holdout(tmp_path):
    receipt = HERE / "results" / "mtrack_emt_acquisition.json"
    require(MTRACK_EMT_SOURCE, receipt)
    manifest = tmp_path / "mtrack_observations.jsonl"
    tensor = tmp_path / "mtrack_rows.npz"
    manifest_report = build_mtrack_emt(
        MTRACK_EMT_SOURCE, receipt, manifest, tensor
    )
    assert manifest_report["counts"] == {
        "timepoints": 10,
        "trajectory_split_groups": 1,
        "biological_replicates": 0,
    }
    assert manifest_report["all_frames_form_one_split_group"] is True
    assert manifest_report["untreated_control_released_in_pilot"] is False
    first_mask = decode_mtrack_mask(
        (MTRACK_EMT_SOURCE / "seg_a549_tgf4ngxy01t003c2.png").read_bytes()
    )
    first_image = decode_mtrack_tiff(
        (MTRACK_EMT_SOURCE / "a549_tgf4ngxy01t003c1.tif").read_bytes()
    )
    assert first_mask.shape == first_image.shape == (1952, 1952)
    assert np.unique(first_mask).tolist() == [0, 1, 2, 3, 4, 5]
    assert (float(first_image.min()), float(first_image.max())) == (347.0, 748.0)

    manifest_report_path = tmp_path / "mtrack_manifest_report.json"
    manifest_report_path.write_text(json.dumps(manifest_report), encoding="utf-8")
    report = evaluate_mtrack_emt(
        tensor,
        manifest_report_path,
        HERE / "results" / "figshare_exif_emt_evaluation.json",
    )
    assert report["design_audit"]["frames_permitted_as_train_test_units"] is False
    assert report["independent_dataset_holdout"] is False
    assert report["independent_lab_holdout"] is False
    assert report["cross_lab_vimentin_direction_context"]["directions_agree"] is False
    assert report["pre_sweep_role"]["numeric_range"] is None
    assert report["pre_sweep_role"]["may_emit_sweep_axis"] is False

    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(), mtrack_emt_report=report
    )
    head = next(
        item for item in multitask["heads"]
        if item["head_id"] == "live_VIM_RFP_EMT_temporal_representation"
    )
    assert head["observable_evidence_eligible"] is False
    assert head["external_holdout"] is False
    assert head["parameter_proposal_eligible"] is False
    gate = evaluate(
        _result("calcium_baseline_report.json"),
        multitask_model_report=multitask,
        mtrack_emt_report=report,
    )
    assert gate["checks"]["mtrack_emt_constraint_is_non_authoritative"] is True
    audit = next(
        item for item in gate["representation_training_audits"]
        if item["dataset"] == "github-mtrack-a549-tgfb4ng-xy01"
    )
    assert audit["fixed_IF_vs_live_VIM_RFP_direction_agrees"] is False
    assert audit["may_count_frames_as_independent_samples"] is False
    assert audit["may_emit_numeric_parameter_range"] is False


def test_lincs_mcf10a_collection_holdout_preserves_biological_hierarchy(tmp_path):
    receipt = HERE / "results" / "lincs_mcf10a_acquisition.json"
    require(LINCS_MCF10A_SOURCE, receipt)
    manifest = tmp_path / "lincs_mcf10a_observations.jsonl"
    tensor = tmp_path / "lincs_mcf10a_rows.npz"
    manifest_report = build_lincs_mcf10a(
        LINCS_MCF10A_SOURCE, receipt, manifest, tensor
    )
    assert manifest_report["counts"] == {
        "author_well_rows": 159,
        "canonical_observations": 1263,
        "biological_replicate_groups": 7,
        "collections": 2,
    }
    assert manifest_report["replicate_contract"][
        "rows_may_be_split_independently"
    ] is False
    assert manifest_report["metadata_conflict"]["resolved"] is False
    with np.load(tensor, allow_pickle=False) as data:
        assert data["X"].shape == (159, 9)
        assert len(set(data["split_group"].astype(str))) == 7
        duplicate = np.flatnonzero(
            data["specimen_name"].astype(str) == "TGFB_48_C2_C"
        )
        assert len(duplicate) == 3
        assert len(set(data["sample_id"][duplicate].astype(str))) == 3

    manifest_report_path = tmp_path / "lincs_mcf10a_manifest_report.json"
    manifest_report_path.write_text(json.dumps(manifest_report), encoding="utf-8")
    report = evaluate_lincs_mcf10a(tensor, manifest_report_path)
    assert report["collection_holdout_passed"] is True
    assert [
        (row["accuracy"], row["AUROC"])
        for row in report["collection_holdout_classifier"]
    ] == [(1.0, 1.0), (pytest.approx(5 / 6), 1.0)]
    for feature in report["paired_direction_holdout"].values():
        assert feature["C1_discovery_C2_holdout_reproduced"] is True
        assert all(
            collection["bootstrap_CI_excludes_zero_in_expected_direction"]
            for collection in feature["collections"].values()
        )
    assert report["OOD_refusal_validated"] is False
    assert report["semantic_OOD_diagnostic"][0][
        "distance_OOD_AUROC"
    ] == pytest.approx(0.78125)
    assert report["semantic_OOD_diagnostic"][1][
        "C2_in_distribution_acceptance"
    ] == pytest.approx(1 / 3)
    assert report["independent_lab_holdout"] is False
    assert report["pre_sweep_role"]["numeric_range"] is None
    assert report["pre_sweep_role"]["may_emit_sweep_axis"] is False

    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(), lincs_mcf10a_report=report
    )
    head = next(
        item for item in multitask["heads"]
        if item["head_id"] == "MCF10A_IF_ligand_state_collection_holdout"
    )
    assert head["observable_evidence_eligible"] is True
    assert head["external_holdout"] is False
    assert head["head_id"] not in multitask["sweep_interface"]["admissible_head_ids"]
    gate = evaluate(
        _result("calcium_baseline_report.json"),
        multitask_model_report=multitask,
        lincs_mcf10a_report=report,
    )
    assert gate["checks"]["lincs_mcf10a_constraint_is_non_authoritative"] is True
    request = next(
        item for item in gate["calibration_requests"]
        if item["parameter"] == "mechanism_family:EMT_spatial_dispersion_and_cell_cell_adhesion"
    )
    assert request["numeric_range"] is None
    assert request["may_emit_sweep_axis"] is False


def test_karacosta_cytof_build_preserves_one_experiment_and_author_states(tmp_path):
    receipt = HERE / "results" / "karacosta_emt_cytof_acquisition.json"
    require(KARACOSTA_EMT_SOURCE, receipt)
    manifest = tmp_path / "karacosta_observations.jsonl"
    tensor = tmp_path / "karacosta_rows.npz"
    report = build_karacosta_emt_cytof(
        KARACOSTA_EMT_SOURCE, receipt, manifest, tensor
    )
    assert report["counts"] == {
        "raw_single_cell_events": 90_066,
        "raw_author_retained_state_events": 87_354,
        "raw_author_unretained_cluster_events": 2_712,
        "transformed_downsampled_state_cells": 29_066,
        "canonical_observations": 90_066,
        "canonical_samples": 90_066,
        "biological_replicate_groups": 1,
        "timepoints": 7,
    }
    assert report["replicate_contract"][
        "single_cells_may_be_split_as_biological_replicates"
    ] is False
    assert report["replicate_contract"][
        "released_Source_Data_maps_that_independent_replicate"
    ] is False
    assert report["label_contract"][
        "state_classification_on_these_features_would_be_label_circular"
    ] is True
    assert report["transformed_state_counts"] == {
        "E1": 2394, "E2": 3284, "E3": 1316, "M": 9609, "MET": 3498,
        "pEMT1": 3359, "pEMT2": 2943, "pEMT3": 2663,
    }
    with np.load(tensor, allow_pickle=False) as data:
        assert data["raw_X"].shape == (90_066, 6)
        assert data["state_X"].shape == (29_066, 6)
        assert len(set(data["split_group"].astype(str))) == 1
        assert set(data["state_label"].astype(str)) == {
            "E1", "E2", "E3", "pEMT1", "pEMT2", "pEMT3", "M", "MET",
        }


def test_karacosta_common_marker_direction_filters_but_cannot_emit_sweep():
    report = _result("karacosta_emt_cytof_evaluation.json")
    progressive = report["progressive_EMT_0d_to_10d_TGFB"]
    withdrawal = report["withdrawal_10d_TGFB_to_10dW"]
    direction = report["cross_lab_common_marker_direction"]
    assert progressive["E_cadherin_down_Vimentin_up_reproduced"] is True
    assert withdrawal["partial_marker_reversal_reproduced"] is True
    assert direction["dataset_and_lab_disjoint"] is True
    assert direction["independent_lab_common_marker_direction_reproduced"] is True
    assert direction["numeric_scale_transfer_permitted"] is False
    assert report["independent_lab_direction_holdout"] is True
    assert report["independent_lab_inferential_holdout"] is False
    assert report["uncertainty_holdout_passed"] is False
    assert report["OOD_refusal_validated"] is False
    assert report["pre_sweep_role"]["numeric_range"] is None
    assert report["pre_sweep_role"]["may_emit_sweep_axis"] is False

    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(), karacosta_emt_report=report
    )
    head = next(
        item for item in multitask["heads"]
        if item["head_id"] == "cross_modal_EMT_common_protein_direction_bridge"
    )
    assert head["external_holdout"] is False
    assert head["observable_evidence_eligible"] is False
    assert head["parameter_proposal_eligible"] is False
    assert head["head_id"] not in multitask["sweep_interface"]["admissible_head_ids"]
    gate = evaluate(
        _result("calcium_baseline_report.json"),
        multitask_model_report=multitask,
        karacosta_emt_report=report,
    )
    assert gate["checks"][
        "karacosta_emt_direction_filter_is_non_authoritative"
    ] is True
    audit = next(
        item for item in gate["representation_training_audits"]
        if item["dataset"] == "nature-2019-karacosta-hcc827-emt-cytof"
    )
    assert audit["status"] == (
        "admitted_independent_lab_common_marker_direction_filter_only"
    )
    assert audit["biological_replicate_groups"] == 1
    request = next(
        item for item in gate["calibration_requests"]
        if item["parameter"] == "mechanism_family:EMT_marker_direction_conditioning"
    )
    assert request["numeric_range"] is None
    assert request["may_emit_sweep_axis"] is False


def test_geo_emt_rna_preserves_biological_samples_and_study_boundaries(tmp_path):
    receipt = HERE / "results" / "geo_emt_rna_acquisition.json"
    require(GEO_EMT_RNA_SOURCE, receipt)
    manifest = tmp_path / "geo_emt_observations.jsonl"
    tensor = tmp_path / "geo_emt_rows.npz"
    report = build_geo_emt_rna(
        GEO_EMT_RNA_SOURCE, receipt, manifest, tensor
    )
    assert report["counts"] == {
        "studies": 5,
        "independent_lab_groups": 5,
        "samples": 72,
        "observations": 144,
        "samples_by_study": {
            "GSE17708": 26, "GSE42373": 8, "GSE125369": 4,
            "GSE69667": 16, "GSE49644": 18,
        },
    }
    assert report["split_contract"]["study_is_outer_holdout_unit"] is True
    assert report["split_contract"]["RNA_features_are_not_independent_replicates"] is True
    rows = [json.loads(line) for line in manifest.read_text().splitlines()]
    assert len({row["sample_id"] for row in rows}) == 72
    assert len({row["lab_group"] for row in rows}) == 5
    assert {row["aleph_authority"] for row in rows} == {"none"}


def test_geo_emt_rna_replicates_direction_but_does_not_replace_sweep():
    report = evaluate_geo_emt_rna(
        HERE / "results" / "geo_emt_rna_rows.npz",
        HERE / "results" / "geo_emt_rna_manifest_report.json",
        HERE / "results" / "karacosta_emt_cytof_evaluation.json",
    )
    replicated = report["replicated_direction"]
    assert replicated["all_ten_study_gene_directions_expected"] is True
    assert replicated["all_ten_between_group_ranges_completely_separated"] is True
    assert replicated[
        "nominal_exact_one_sided_sign_test_p_for_five_retrospective_study_signatures"
    ] == 0.03125
    assert replicated["retrospective_direction_replication_passed"] is True
    assert replicated["study_selection_was_outcome_blind"] is False
    assert replicated["confirmatory_direction_threshold_0_05_passed"] is False
    assert report["study_level_direction_inference_passed"] is False
    assert report["numeric_uncertainty_promotion_passed"] is False
    assert report["within_GSE49644_cross_cell_line_transfer"][
        "all_four_marker_directions_agree"
    ] is True
    assert report["cross_assay_bridge"][
        "RNA_direction_agrees_with_independent_IF_and_CyTOF_protein_direction"
    ] is True
    assert report["semantic_OOD_challenge"]["challenge_defined"] is True
    assert report["semantic_OOD_challenge"]["diagnostic_passed"] is True
    assert report["semantic_OOD_challenge"]["OOD_refusal_model_validated"] is False
    assert report["pre_sweep_role"]["numeric_range"] is None
    assert report["pre_sweep_role"]["may_emit_sweep_axis"] is False
    assert report["may_select_aleph_parameter"] is False


def test_geo_emt_rna_enters_multitask_and_gate_only_as_pre_sweep_constraint():
    geo = _result("geo_emt_rna_evaluation.json")
    karacosta = _result("karacosta_emt_cytof_evaluation.json")
    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(),
        karacosta_emt_report=karacosta,
        geo_emt_rna_report=geo,
    )
    head = next(
        item for item in multitask["heads"]
        if item["head_id"] == "replicated_multi_study_EMT_RNA_direction_constraint"
    )
    without_geo = evaluate_multitask_outer_model(
        *_multitask_inputs(), karacosta_emt_report=karacosta
    )
    assert multitask["architecture"]["head_count"] == (
        without_geo["architecture"]["head_count"] + 1
    )
    assert head["split_contract"]["studies"] == 5
    assert head["split_contract"]["independent_lab_groups"] == 5
    assert head["split_contract"]["biological_replicates_not_expression_features"] is True
    assert head["observable_evidence_eligible"] is False
    assert head["uncertainty_eligible"] is False
    assert head["external_holdout"] is False
    assert head["parameter_proposal_eligible"] is False
    assert head["head_id"] not in multitask["sweep_interface"]["admissible_head_ids"]

    gate = evaluate(
        _result("calcium_baseline_report.json"),
        multitask_model_report=multitask,
        karacosta_emt_report=karacosta,
        geo_emt_rna_report=geo,
    )
    assert gate["checks"][
        "geo_emt_rna_direction_constraint_is_non_authoritative"
    ] is True
    audit = next(
        item for item in gate["representation_training_audits"]
        if item["dataset"] == "geo-emt-rna-five-study-panel"
    )
    assert audit["status"] == "admitted_retrospective_replicated_study_direction_constraint_only"
    assert audit["study_level_sign_test_p"] == 0.03125
    assert audit["study_selection_was_outcome_blind"] is False
    assert audit["may_count_as_confirmatory_direction_holdout"] is False
    assert audit["may_count_as_numeric_uncertainty_holdout"] is False
    request = next(
        item for item in gate["calibration_requests"]
        if item["parameter"] == "mechanism_family:replicated_EMT_RNA_direction_conditioning"
    )
    assert request["numeric_range"] is None
    assert request["may_emit_sweep_axis"] is False
    assert gate["status"] == "refused"


def test_geo_emt_rna_overclaim_is_rejected_by_multitask_and_gate():
    geo = json.loads(json.dumps(_result("geo_emt_rna_evaluation.json")))
    geo["numeric_uncertainty_promotion_passed"] = True
    with pytest.raises(ValueError, match="overclaims sweep authority"):
        evaluate_multitask_outer_model(
            *_multitask_inputs(), geo_emt_rna_report=geo
        )

    core = json.loads(json.dumps(_result("calcium_baseline_report.json")))
    core["test"]["macro_f1"] = 0.8
    core["test"]["ece_10_bin"] = 0.1
    core["semantic_ood"]["mean_auroc"] = 0.8
    core["independent_dataset_task_holdout_count"] = 1
    core["independent_lab_task_holdout_count"] = 1
    core["independent_dataset_task_holdout_passed_count"] = 1
    core["independent_lab_task_holdout_passed_count"] = 1
    core["aleph_authority"] = "none"
    assert evaluate(core)["status"] == "eligible_for_non_authoritative_proposal"
    gate = evaluate(core, geo_emt_rna_report=geo)
    assert gate["checks"][
        "geo_emt_rna_direction_constraint_is_non_authoritative"
    ] is False
    assert gate["status"] == "refused"


def test_gse325309_holdout_protocol_is_frozen_before_expression_access():
    protocol = json.loads(
        (HERE / "gse325309_preregistered_holdout_protocol.json").read_text()
    )
    assert protocol["frozen_before_expression_file_download_or_inspection"] is True
    source = protocol["source_metadata"]
    controls = set(source["control_biological_replicates"])
    treated = set(source["TGFB1_biological_replicates"])
    assert len(controls) == len(treated) == 3
    assert controls.isdisjoint(treated)
    assert protocol["frozen_feature_contract"]["gene_symbols"] == ["CDH1", "VIM"]
    assert protocol["frozen_feature_contract"]["sample_exclusion"] == "none"
    primary = protocol["frozen_primary_test"]
    assert primary["null_distribution"] == (
        "exact_all_20_assignments_of_3_of_6_samples_to_TGFB1"
    )
    assert primary["alpha"] == 0.05
    authority = protocol["authority_contract"]
    assert authority["numeric_Aleph_range_permitted"] is False
    assert authority["may_emit_sweep_axis"] is False
    assert authority["may_select_Aleph_parameter"] is False


def test_gse325309_pristine_builder_preserves_all_six_biological_samples(tmp_path):
    receipt = HERE / "results" / "gse325309_pristine_acquisition.json"
    protocol = HERE / "gse325309_preregistered_holdout_protocol.json"
    require(GSE325309_SOURCE, receipt)
    manifest = tmp_path / "gse325309.jsonl"
    tensor = tmp_path / "gse325309.npz"
    report = build_gse325309_pristine(
        GSE325309_SOURCE, receipt, protocol, manifest, tensor
    )
    assert report["counts"] == {
        "biological_samples": 6, "observations": 12, "genes": 2,
        "control_replicates": 3, "TGFB1_replicates": 3,
    }
    assert report["sample_exclusions"] == 0
    rows = [json.loads(line) for line in manifest.read_text().splitlines()]
    assert len({row["sample_id"] for row in rows}) == 6
    assert len({row["biological_replicate"] for row in rows}) == 6
    assert {row["aleph_authority"] for row in rows} == {"none"}


def test_gse325309_pristine_exact_test_passes_direction_only():
    report = evaluate_gse325309_pristine(
        HERE / "results" / "gse325309_pristine_rows.npz",
        HERE / "results" / "gse325309_pristine_manifest_report.json",
    )
    primary = report["frozen_primary_test"]
    assert primary["directions"] == {"CDH1": "decrease", "VIM": "increase"}
    assert primary["exact_assignment_count"] == 20
    assert primary["exact_one_sided_p"] == 0.05
    assert primary["passed"] is True
    assert report["pristine_external_direction_holdout_passed"] is True
    assert report["numeric_uncertainty_promotion_passed"] is False
    assert report["OOD_refusal_validated"] is False
    assert report["pre_sweep_role"]["numeric_range"] is None
    assert report["pre_sweep_role"]["may_emit_sweep_axis"] is False


def test_gse325309_pristine_holdout_enters_model_and_gate_without_numeric_authority():
    report = _result("gse325309_pristine_evaluation.json")
    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(), gse325309_emt_report=report
    )
    head = next(item for item in multitask["heads"] if item["head_id"] ==
                "preregistered_GSE325309_EMT_direction_holdout")
    assert head["external_holdout"] is True
    assert head["pristine_confirmation"] is True
    assert head["observable_evidence_eligible"] is False
    assert head["uncertainty_eligible"] is False
    assert head["parameter_proposal_eligible"] is False
    assert head["head_id"] not in multitask["sweep_interface"]["admissible_head_ids"]
    gate = evaluate(
        _result("calcium_baseline_report.json"),
        multitask_model_report=multitask,
        gse325309_emt_report=report,
    )
    assert gate["checks"][
        "gse325309_preregistered_direction_holdout_is_non_authoritative"
    ] is True
    audit = next(item for item in gate["representation_training_audits"]
                 if item["dataset"] == "geo-gse325309-preregistered-a549-tgfb1-rnaseq")
    assert audit["status"] == "admitted_pristine_external_direction_holdout_only"
    assert audit["may_count_as_confirmatory_direction_holdout"] is True
    assert audit["may_count_as_numeric_uncertainty_holdout"] is False
    request = next(item for item in gate["calibration_requests"] if item["parameter"] ==
                   "mechanism_family:preregistered_EMT_direction_confirmation")
    assert request["numeric_range"] is None
    assert request["may_emit_sweep_axis"] is False


def test_dryad_supracontractility_preserves_author_styles_and_replicates():
    require(
        DRYAD_SUPRACONTRACTILITY_SOURCE,
        DRYAD_SUPRACONTRACTILITY_EXPORT,
    )
    audit = audit_dryad_supracontractility_sources(
        DRYAD_SUPRACONTRACTILITY_SOURCE
    )
    assert sum(
        len(item["author_excluded"])
        for item in audit["annotations"].values()
    ) == 112
    assert len(audit["annotations"]["Data_Fig2.xlsx"][
        "maximum_density_alignment"
    ]) == 9
    rows, report = build_dryad_supracontractility(
        DRYAD_SUPRACONTRACTILITY_EXPORT,
        HERE / "results" / "dryad_08kprr59c_source_audit.json",
    )
    assert len(rows) == 2_155
    assert len({row.sample_id for row in rows}) == 16
    assert len({row.biological_replicate for row in rows}) == 16
    assert report["author_excluded_cells_not_emitted"] == 112
    assert {row.aleph_authority for row in rows} == {"none"}
    assert {row.modality for row in rows} == {
        "particle_image_velocimetry_and_live_actin_derived",
        "immunofluorescence_fibronectin_orientation_derived",
    }


def test_dryad_supracontractility_is_directional_pre_sweep_only():
    manifest = HERE / "results" / "dryad_08kprr59c_supracontractility_observations.jsonl"
    report = evaluate_dryad_supracontractility(manifest)
    assert report["replicated_contractility_direction"] is True
    late = report["late_post_blebbistatin_contrasts_t_ge_42h"]
    assert late["orientation_order_parameter"]["group_A_minus_group_B"] > 0.28
    assert late["nematic_correlation_length"]["group_A_minus_group_B"] > 390
    assert late["velocity_correlation_length"]["group_A_minus_group_B"] > 54
    fibronectin = report["fibronectin_orientation_nematic_order"][
        "control_minus_blebbistatin"
    ]
    assert fibronectin["group_A_minus_group_B"] > 0.53
    assert fibronectin["exact_randomization_two_sided_p"] == pytest.approx(1 / 3)
    assert report["pre_sweep_constraint"][
        "future_sweep_range_reduction_intended"
    ] is True
    assert report["pre_sweep_constraint"]["may_emit_numeric_parameter_range"] is False
    assert report["observable_evidence_eligible"] is False
    assert report["independent_lab_holdout"] is False
    assert report["aleph_parameter_proposal_eligible"] is False
    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(), piv_contractility_report=report
    )
    head = next(
        item for item in multitask["heads"]
        if item["head_id"] == (
            "C2C12_PIV_IF_supracellular_contractility_pre_sweep_constraint"
        )
    )
    assert head["status"] == "training_only_directional_constraint_operator_and_external_holdout_blocked"
    assert head["observable_evidence_eligible"] is False
    assert head["head_id"] not in multitask["sweep_interface"]["admissible_head_ids"]
    assert multitask["sweep_interface"]["may_emit_numeric_parameter_range"] is False
    gate = evaluate(
        _result("calcium_baseline_report.json"),
        multitask_model_report=multitask,
        piv_contractility_report=report,
    )
    assert gate["checks"]["piv_contractility_constraint_is_non_authoritative"] is True
    audit = next(
        item for item in gate["representation_training_audits"]
        if item["dataset"] == "dryad-08kprr59c-c2c12-supracontractility"
    )
    assert audit["status"] == "admitted_training_direction_only"
    request = next(
        item for item in gate["calibration_requests"]
        if item["parameter"] == "mechanism_family:supracellular_contractility"
    )
    assert request["numeric_range"] is None
    assert request["may_emit_sweep_axis"] is False


def test_dataverse_cellular_nematics_preserves_sequence_and_FOV_hierarchy(tmp_path):
    require(DATAVERSE_CELLULAR_NEMATICS_SOURCE)
    rows, report = build_dataverse_cellular_nematics(
        DATAVERSE_CELLULAR_NEMATICS_SOURCE,
        HERE / "results" / "dataverse_data2772_acquisition.json",
        tmp_path / "piv_vectors.npz",
    )
    assert len(rows) == 372
    assert report["biological_replicates_claimed"] == 0
    assert report["piv"]["sequence_count"] == 1
    assert report["piv"]["frame_count"] == 63
    assert report["piv"]["vectors_admitted_total"] == 239_150
    assert report["piv"]["vectors_excluded_by_flag"] == {"-22": 78, "9999": 2_944}
    assert report["tfm_msm"]["technical_fields_of_view"] == 12
    assert report["tfm_msm"]["biological_culture_or_batch_count"] == 1
    assert report["tfm_msm"]["withheld_table_index"] == 2
    assert {row.aleph_authority for row in rows} == {"none"}


def test_cross_cell_type_force_data_ranks_future_sweep_without_replacing_it():
    nematic = evaluate_dataverse_cellular_nematics(
        HERE / "results" / "dataverse_data2772_cellular_nematic_observations.jsonl",
        HERE / "results" / "dataverse_data2772_manifest_report.json",
    )
    traction = nematic["figure_S9_paired_technical_FOV_contrasts"][
        "mean_traction_magnitude"
    ]
    tension = nematic["figure_S9_paired_technical_FOV_contrasts"][
        "mean_maximum_principal_tension"
    ]
    assert traction["before_FOV_mean"] == pytest.approx(136.88581948980973)
    assert traction["after_FOV_mean"] == pytest.approx(30.178018661492075)
    assert tension["before_FOV_mean"] == pytest.approx(25.402471159298596)
    assert tension["after_FOV_mean"] == pytest.approx(5.1825546655008505)
    assert traction["FOVs_with_decrease"] == 12
    assert tension["FOVs_with_decrease"] == 12
    assert nematic["source_role_discrepancy"]["t2_admitted"] is False
    assert nematic["pre_sweep_constraint"]["may_emit_numeric_parameter_range"] is False

    opto = evaluate_dryad_force_propagation(
        HERE / "results" / "dryad_sj3tx9683_force_propagation_observations.jsonl",
        HERE / "results" / "dryad_sj3tx9683_manifest_report.json",
        HERE / "results" / "dryad_sj3tx9683_source_audit.json",
    )
    assert opto["overall_direction"] == {
        "usually_positive_in_both_axes": True,
        "x_positive_objects": 47, "x_total": 68,
        "y_positive_objects": 46, "y_total": 68,
    }
    assert opto["design_audit"]["objects_treated_as_biological_replicates"] is False
    assert opto["design_audit"]["stats_fullstim_admitted"] is False

    c2c12 = _result("dryad_08kprr59c_supracontractility_evidence.json")
    multitask = evaluate_multitask_outer_model(
        *_multitask_inputs(), piv_contractility_report=c2c12,
        nematic_tfm_report=nematic, opto_force_report=opto,
    )
    head = next(
        item for item in multitask["heads"]
        if item["head_id"] == "cross_cell_type_contractility_force_pre_sweep_constraint"
    )
    assert head["split_contract"]["lab_group_count"] == 3
    assert head["split_contract"]["aligned_same_observable_holdout"] is False
    assert head["observable_evidence_eligible"] is False
    assert head["parameter_proposal_eligible"] is False
    gate = evaluate(
        _result("calcium_baseline_report.json"),
        multitask_model_report=multitask,
        piv_contractility_report=c2c12,
        nematic_tfm_report=nematic,
        opto_force_report=opto,
    )
    request = next(
        item for item in gate["calibration_requests"]
        if item["parameter"] == "mechanism_family:cross_cell_type_actomyosin_contractility"
    )
    assert request["numeric_range"] is None
    assert request["may_emit_sweep_axis"] is False
    assert len(request["required_forward_observables"]) >= 4


def test_dryad_force_propagation_rejects_mismatched_stats_and_pickle_execution(tmp_path):
    require(DRYAD_FORCE_PROPAGATION_SOURCE)
    audit = audit_dryad_force_propagation_source(DRYAD_FORCE_PROPAGATION_SOURCE)
    assert audit["csv_contract"]["raw_fullstim_group_counts"] == {
        "AR1to1d": 22, "AR1to2d": 30, "AR2to1d": 16,
    }
    assert audit["csv_contract"]["stats_fullstim_group_counts"] == {
        "AR1to1d": 35, "AR1to1s": 15,
    }
    assert audit["csv_contract"]["stats_fullstim_matches_raw_fullstim"] is False
    assert len(audit["untrusted_pickle_members_not_executed"]) == 9
    audit_path = tmp_path / "audit.json"
    audit_path.write_text(json.dumps(audit), encoding="utf-8")
    rows, report = build_dryad_force_propagation(
        DRYAD_FORCE_PROPAGATION_SOURCE, audit_path
    )
    assert len(rows) == 476
    assert report["samples"] == 1
    assert report["stats_fullstim_csv_admitted"] is False
    assert report["untrusted_dat_pickles_executed"] is False


def test_committed_summary_reconciles_every_generated_manifest_and_tensor():
    results = HERE / "results"
    summary = json.loads((results / "summary.json").read_text(encoding="utf-8"))
    observations = 0
    datasets, labs, samples = set(), set(), set()
    modalities: dict[str, int] = {}
    for filename, expected_digest in summary["generated_manifest_sha256"].items():
        path = results / filename
        require(path)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_digest
        for line in path.read_text(encoding="utf-8").splitlines():
            row = json.loads(line)
            observations += 1
            datasets.add(row["dataset_id"])
            labs.add(row["lab_group"])
            samples.add((row["dataset_id"], row["sample_id"]))
            modalities[row["modality"]] = modalities.get(row["modality"], 0) + 1
    assert summary["dataset_count"] == len(datasets)
    assert summary["lab_group_count"] == len(labs)
    assert summary["observation_count"] == observations
    assert summary["sample_count"] == len(samples)
    assert summary["modalities"] == modalities
    for filename, expected_digest in summary["generated_tensor_sha256"].items():
        path = results / filename
        require(path)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_digest
    for filename, expected_digest in summary["generated_report_sha256"].items():
        path = results / filename
        require(path)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_digest
    assert summary["authority"] == "none"
    assert summary["production_eligible"] is False


def test_idr0072_corrected_flex_manifest_is_one_to_one_and_frozen():
    assert idr0072_expected_filename({"screen_id": 2952, "well_name": "m20"}) == "013020001.flex"
    assert idr0072_expected_filename({"screen_id": 2953, "well_name": "A13"}) == "001013000.flex"
    with pytest.raises(ValueError):
        idr0072_expected_filename({"screen_id": 2952, "well_name": "Q25"})

    collision = _result("idr0072_if_flex_resolution_collision.json")
    assert collision["resolved_manifest_status"] == "invalid_retired_before_prediction"
    assert collision["input_image_ids"] == 2594
    assert collision["unique_source_urls"] == 81
    assert collision["model_predictions_made"] == 0

    report = _result("idr0072_if_flex_resolution_v2.json")
    manifest_path = HERE / "results" / "idr0072_if_flex_manifest_v2.jsonl"
    assert hashlib.sha256(manifest_path.read_bytes()).hexdigest() == report["resolved_manifest_sha256"]
    rows = [json.loads(line) for line in manifest_path.read_text(encoding="utf-8").splitlines()]
    assert len(rows) == report["resolved_rows"] == 2594
    assert len({row["image_id"] for row in rows}) == 2594
    assert len({row["source_url"] for row in rows}) == 2594
    assert len({row["local_relative_path"] for row in rows}) == 2594
    assert all(row["flex_filename"] == idr0072_expected_filename(row) for row in rows)

    wrapper = json.loads((HERE / "idr0072_if_corrected_source_protocol.json").read_text(encoding="utf-8"))
    base = HERE / wrapper["base_optical_protocol"]["path"]
    assert hashlib.sha256(base.read_bytes()).hexdigest() == wrapper["base_optical_protocol"]["sha256"]
    assert wrapper["corrected_source_contract"]["flex_manifest_sha256"] == report["resolved_manifest_sha256"]
    assert wrapper["authority_contract"]["pristine_holdout_authority"] is False
    assert wrapper["authority_contract"]["aleph_parameter_authority"] == "none"


def test_idr0072_external_if_failure_cannot_gain_authority():
    report = _result("idr0072_if_evaluation.json")
    assert report["input_image_count"] == 2594
    assert report["evaluated_image_count"] == 2593
    assert report["pre_prediction_refusal_count"] == 1
    assert report["status"] == "domain_shift_diagnostic_failed"
    assert report["all_applicable_endpoints_passed"] is False
    assert report["endpoint_results"] == {
        "ECE": False,
        "macro_F1": False,
        "macro_average_precision": False,
        "mean_prediction_set_size": False,
        "positive_conformal_coverage": False,
    }
    assert report["pristine_holdout_authority"] is False
    assert report["uncertainty_authority"] is False
    assert report["aleph_parameter_authority"] == "none"
    assert report["may_emit_numeric_sweep_range"] is False
    tensor = HERE / "results" / "idr0072_if_predictions.npz"
    assert hashlib.sha256(tensor.read_bytes()).hexdigest() == report["tensor_sha256"]


def test_idr0072_failed_transfer_enters_multitask_only_as_consumed_diagnostic():
    baseline = evaluate_multitask_outer_model(*_multitask_inputs())
    report = evaluate_multitask_outer_model(
        *_multitask_inputs(), idr0072_if_report=_result("idr0072_if_evaluation.json")
    )
    assert report["architecture"]["head_count"] == (
        baseline["architecture"]["head_count"] + 1
    )
    assert report["architecture"]["trainable_raw_encoder_status"] == (
        "three_encoders_trained_external_IF_failed"
    )
    head = {item["head_id"]: item for item in report["heads"]}[
        "raw_IF_subcellular_localization_transfer"
    ]
    assert head["status"] == "blocked_external_domain_shift"
    assert head["external_holdout"] is True
    assert head["pristine_confirmation"] is False
    assert head["uncertainty_eligible"] is False
    assert head["observable_evidence_eligible"] is False
    assert head["parameter_proposal_eligible"] is False
    assert head["split_contract"][
        "diagnostic_consumed_and_may_not_be_reused_as_pristine"
    ] is True
    assert "raw_IF_subcellular_localization_transfer" in report["readiness"][
        "blocked_external_classifier_heads"
    ]
    assert report["readiness"]["general_cell_type_state_model_ready"] is False


def test_idr0168_source_integrity_refusal_is_consumed_without_predictions():
    protocol_path = HERE / "idr0168_multiprovider_if_confirmation_protocol.json"
    acquisition_path = HERE / "results" / "idr0168_confirmation_acquisition.json"
    refusal = _result("idr0168_confirmation_decode_refusal.json")
    assert refusal["protocol_sha256"] == hashlib.sha256(
        protocol_path.read_bytes()
    ).hexdigest()
    assert refusal["acquisition_sha256"] == hashlib.sha256(
        acquisition_path.read_bytes()
    ).hexdigest()
    assert refusal["status"] == (
        "failed_source_integrity_and_permanently_consumed"
    )
    assert refusal["source_integrity"]["gzip_integrity_passed"] is False
    assert refusal["source_integrity"]["native_Zarr_metadata_present"] is False
    assert refusal["source_integrity"][
        "all_sampled_remote_ranges_match_acquired_object"
    ] is True
    assert len(refusal["source_integrity"]["remote_range_checks"]) == 5
    assert refusal["tensor_written"] is False
    assert refusal["model_predictions_made"] == 0
    assert refusal["replacement_field_used"] is False
    assert refusal["single_use_confirmation_consumed"] is True
    assert refusal["may_be_reused_as_pristine"] is False
    assert refusal["Aleph_parameter_authority"] == "none"


def test_idr0168_refusal_keeps_external_if_head_blocked_and_rejects_rewrite():
    refusal = _result("idr0168_confirmation_decode_refusal.json")
    report = evaluate_multitask_outer_model(
        *_multitask_inputs(),
        idr0072_if_report=_result("idr0072_if_evaluation.json"),
        idr0168_decode_refusal=refusal,
    )
    head = {item["head_id"]: item for item in report["heads"]}[
        "raw_IF_subcellular_localization_transfer"
    ]
    assert head["status"] == "blocked_external_domain_shift"
    assert head["metrics"]["IDR0168_single_use_confirmation"][
        "model_predictions_made"
    ] == 0
    assert head["split_contract"][
        "single_use_confirmation_consumed_before_prediction"
    ] is True
    assert head["observable_evidence_eligible"] is False
    assert head["parameter_proposal_eligible"] is False
    assert head["head_id"] not in report["sweep_interface"][
        "admissible_head_ids"
    ]

    tampered = json.loads(json.dumps(refusal))
    tampered["model_predictions_made"] = 1
    with pytest.raises(ValueError, match="source-integrity refusal was rewritten"):
        evaluate_multitask_outer_model(
            *_multitask_inputs(),
            idr0072_if_report=_result("idr0072_if_evaluation.json"),
            idr0168_decode_refusal=tampered,
        )


def test_lightmycells_metadata_freezes_cross_study_roles_without_confirmation_pixels():
    protocol_path = HERE / "lightmycells_preregistered_protocol.json"
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    report = _result("lightmycells_metadata_audit.json")
    assert report["protocol_sha256"] == hashlib.sha256(
        protocol_path.read_bytes()
    ).hexdigest()
    assert protocol["access_history_before_freeze"]["image_pixels_opened"] == 0
    assert report["study_count"] == 30
    assert report["listed_file_count"] == 56_984
    assert report["listed_total_bytes"] == 109_287_061_290
    assert report["matched_acquisition_set_count"] == 2_574
    assert report["unmatched_file_count"] == 0
    assert report["matched_channel_file_counts"] == {
        "Actin": 27,
        "BF": 41_213,
        "DIC": 3_499,
        "Mitochondria": 1_819,
        "Nucleus": 2_533,
        "PC": 7_670,
        "Tubulin": 223,
    }
    assert report["confirmation_pixels_accessed"] == 0
    assert report["model_predictions_made"] == 0
    assert report["aleph_authority"] == "none"
    confirmation = report["split_summary"]["single_use_confirmation"]
    assert confirmation["study_count"] == 6
    assert confirmation["acquisition_set_count"] == 102
    assert confirmation["target_acquisition_set_counts"]["Nucleus"] == 102
    assert confirmation["target_acquisition_set_counts"]["Mitochondria"] == 91
    assert confirmation["target_acquisition_set_counts"]["Tubulin"] == 23
    assert confirmation["target_acquisition_set_counts"]["Actin"] == 10

    manifest = HERE / "results" / "lightmycells_metadata_manifest.jsonl"
    rows = [json.loads(line) for line in manifest.read_text().splitlines()]
    assert len(rows) == 2_574
    assert hashlib.sha256(manifest.read_bytes()).hexdigest() == report["manifest_sha256"]
    assert len({(row["study_id"], row["acquisition_set_id"]) for row in rows}) == len(rows)
    role_names = ("fit", "selection", "calibration", "single_use_confirmation")
    by_role = {
        role: {row["study_id"] for row in rows if row["role"] == role}
        for role in role_names
    }
    assert sum(len(studies) for studies in by_role.values()) == 30
    assert not any(
        by_role[left] & by_role[right]
        for index, left in enumerate(role_names)
        for right in role_names[index + 1:]
    )
    assert all(row["pixel_content_accessed"] is False for row in rows)
    assert all(row["aleph_authority"] == "none" for row in rows)


def test_lightmycells_v1_pair_resolution_is_compact_and_records_position_refusals():
    report = _result("lightmycells_pair_resolution.json")
    assert report["pair_count"] == 3_253
    assert report["unique_selected_payload_count"] == 5_971
    assert report["selected_payload_bytes"] == 6_929_551_641
    assert report["header_request_count"] == 6_518
    assert report["header_error_count"] == 2_783
    assert report["refusal_count"] == 841
    assert report["confirmation_study_count_skipped_before_header_access"] == 6
    assert report["confirmation_header_requests"] == 0
    assert report["confirmation_pixels_accessed"] == 0
    assert report["model_predictions_made"] == 0
    assert report["aleph_authority"] == "none"
    assert sum(
        "OME Plane lacks PositionZ" in item["reason"]
        for item in report["header_errors"]
    ) == 2_782
    assert sum(
        "not a TIFF byte order marker" in item["reason"]
        for item in report["header_errors"]
    ) == 1
    pair_path = HERE / "results" / "lightmycells_spatial_pairs.jsonl"
    assert hashlib.sha256(pair_path.read_bytes()).hexdigest() == report["output_sha256"]
    rows = [json.loads(line) for line in pair_path.read_text().splitlines()]
    assert len(rows) == 3_253
    assert {row["role"] for row in rows} == {"fit", "selection", "calibration"}
    assert all(row["confirmation_pixels_accessed"] is False for row in rows)


def test_lightmycells_v2_fallback_recovers_pairs_without_touching_confirmation():
    wrapper_path = HERE / "lightmycells_pair_resolution_protocol.json"
    wrapper = json.loads(wrapper_path.read_text())
    base = HERE / wrapper["base_protocol"]["path"]
    diagnostic = HERE / wrapper["source_informed_diagnostic"]["path"]
    assert hashlib.sha256(base.read_bytes()).hexdigest() == wrapper["base_protocol"]["sha256"]
    assert hashlib.sha256(diagnostic.read_bytes()).hexdigest() == wrapper[
        "source_informed_diagnostic"
    ]["sha256"]
    report = _result("lightmycells_pair_resolution_v2.json")
    assert report["resolution_protocol_sha256"] == hashlib.sha256(
        wrapper_path.read_bytes()
    ).hexdigest()
    assert report["resolution_fallback_enabled"] is True
    assert report["pair_count"] == 4_375
    assert report["unique_selected_payload_count"] == 7_933
    assert report["selected_payload_bytes"] == 12_678_835_256
    assert report["alignment_method_pair_counts"] == {
        "lower_median_z_fallback": 1_065,
        "nearest_PositionZ": 3_253,
        "singleton_input": 57,
    }
    assert report["role_pair_counts"] == {
        "calibration": 569,
        "fit": 3_744,
        "selection": 62,
    }
    assert report["refusal_count"] == 1
    assert report["refusals"][0]["study_id"] == "Study_25"
    assert report["refusals"][0]["target"] == "Nucleus"
    assert report["confirmation_header_requests"] == 0
    assert report["confirmation_pixels_accessed"] == 0
    pair_path = HERE / "results" / "lightmycells_spatial_pairs_v2.jsonl"
    assert hashlib.sha256(pair_path.read_bytes()).hexdigest() == report["output_sha256"]
    rows = [json.loads(line) for line in pair_path.read_text().splitlines()]
    assert len(rows) == 4_375
    assert all(row["role"] != "single_use_confirmation" for row in rows)
    assert all(row["aleph_authority"] == "none" for row in rows)
    assert sum(row["exact_position_alignment"] for row in rows) == 3_253


def test_lightmycells_integrity_filter_acquisition_and_tensor_are_frozen():
    protocol_path = HERE / "lightmycells_source_integrity_protocol.json"
    protocol = json.loads(protocol_path.read_text())
    report = _result("lightmycells_source_integrity_filter.json")
    assert report["integrity_protocol_sha256"] == hashlib.sha256(
        protocol_path.read_bytes()
    ).hexdigest()
    assert report["input_pair_count"] == 4_375
    assert report["admitted_pair_count"] == 4_374
    assert report["refused_pair_count"] == 1
    assert report["unique_admitted_payload_count"] == 7_932
    assert report["admitted_payload_bytes"] == 12_672_726_679
    assert report["model_predictions_made"] == 0
    assert report["confirmation_pixels_accessed"] == 0
    receipt = _result("lightmycells_pair_acquisition.json")
    assert receipt["integrity_protocol_sha256"] == hashlib.sha256(
        protocol_path.read_bytes()
    ).hexdigest()
    assert receipt["integrity_report_sha256"] == hashlib.sha256(
        (HERE / "results" / "lightmycells_source_integrity_filter.json").read_bytes()
    ).hexdigest()
    assert receipt["file_count"] == 7_932
    assert receipt["total_bytes"] == 12_672_726_679
    assert len({item["path"] for item in receipt["files"]}) == 7_932
    assert all(len(item["sha256"]) == 64 for item in receipt["files"])
    assert receipt["confirmation_file_count"] == 0
    assert receipt["confirmation_pixels_accessed"] == 0

    tensor_path = HERE / "results" / "lightmycells_spatial_tensor.npz"
    tensor_report = _result("lightmycells_spatial_tensor_report.json")
    assert hashlib.sha256(tensor_path.read_bytes()).hexdigest() == tensor_report[
        "tensor_sha256"
    ]
    assert tensor_report["pair_count"] == 4_374
    assert tensor_report["shape"] == [4_374, 64, 64]
    assert tensor_report["role_counts"] == {
        "calibration": 569,
        "fit": 3_744,
        "selection": 61,
    }
    assert tensor_report["study_count"] == 24
    assert tensor_report["acquisition_set_count"] == 2_472
    assert tensor_report["confirmation_pair_count"] == 0
    with np.load(tensor_path, allow_pickle=False) as data:
        assert data["inputs"].shape == data["targets"].shape == (4_374, 64, 64)
        assert data["inputs"].dtype == data["targets"].dtype == np.float16
        assert np.isfinite(data["inputs"]).all()
        assert np.isfinite(data["targets"]).all()
        assert 0 <= data["inputs"].min() <= data["inputs"].max() <= 1
        assert 0 <= data["targets"].min() <= data["targets"].max() <= 1
        assert "single_use_confirmation" not in set(data["role"].tolist())


def test_lightmycells_v1_model_is_frozen_but_not_external_confirmation():
    protocol_path = HERE / "lightmycells_model_training_protocol.json"
    report = _result("lightmycells_model_training_report.json")
    checkpoint = HERE / "results" / "lightmycells_model_checkpoint.pt"
    assert report["protocol_sha256"] == hashlib.sha256(
        protocol_path.read_bytes()
    ).hexdigest()
    assert report["checkpoint_sha256"] == hashlib.sha256(
        checkpoint.read_bytes()
    ).hexdigest()
    assert report["architecture"] == {
        "embedding_dimension": 32,
        "family": "conditional_small_UNet",
        "parameter_count": 26_641,
        "selected_base_channels": 16,
        "selected_epoch": 9,
    }
    history = [
        (
            row["selection_required_target_macro_study_MAE"],
            candidate["base_channels"],
            row["epoch"],
        )
        for candidate in report["candidate_selection"]
        for row in candidate["history"]
    ]
    assert min(history) == (0.12536347564309835, 16, 9)
    selection = report["selection_required_target_macro"]
    assert selection["Pearson_correlation"] < 0.35
    assert selection["global_SSIM"] < 0.25
    assert selection["relative_MAE_improvement_over_fit_mean_target"] < 0.10
    assert report["uncertainty_calibration"]["Nucleus"]["pixel_coverage"] > 0.95
    assert report["uncertainty_calibration"]["Mitochondria"]["pixel_coverage"] > 0.95
    assert report["independent_dataset_or_lab_confirmation_performed"] is False
    assert report["confirmation_pairs_seen"] == 0
    assert report["production_eligible"] is False
    assert report["observable_evidence_eligible"] is False
    assert report["uncertainty_authority"] is False
    assert report["aleph_authority"] == "none"


def test_lightmycells_single_use_confirmation_fails_without_authority():
    acquisition = _result("lightmycells_confirmation_acquisition.json")
    tensor = _result("lightmycells_confirmation_tensor_report.json")
    report = _result("lightmycells_confirmation_evaluation.json")
    assert acquisition["study_count"] == 6
    assert acquisition["pair_count"] == 226
    assert acquisition["target_pair_counts"] == {
        "Actin": 10,
        "Mitochondria": 91,
        "Nucleus": 102,
        "Tubulin": 23,
    }
    assert acquisition["integrity_refused_pair_count"] == 0
    assert acquisition["all_confirmation_labs_independent"] is False
    assert len(acquisition["novel_confirmation_experimenter_groups"]) == 2
    assert tensor["confirmation_consumed"] is True
    assert tensor["confirmation_pair_count"] == 226
    assert tensor["input_constant_image_count"] == 0
    assert tensor["target_constant_image_count"] == 0
    assert report["status"] == "blocked_external_spatial_transfer"
    assert report["all_applicable_confirmation_endpoints_passed"] is False
    assert report["frozen_endpoint_passes"] == {
        "embedding_OOD_AUROC": False,
        "macro_study_Pearson_correlation": False,
        "macro_study_SSIM": True,
        "pixel_interval_coverage": True,
        "relative_MAE_improvement_over_fit_mean_target": True,
    }
    values = report["frozen_endpoint_values"]
    assert values["macro_study_Pearson_correlation"] < 0.35
    assert values["embedding_OOD_AUROC"] < 0.70
    assert report["observable_evidence_eligible"] is False
    assert report["uncertainty_authority_for_this_spatial_task"] is False
    assert report["general_cell_state_production_eligible"] is False
    assert report["aleph_authority"] == "none"
    assert report["may_emit_numeric_sweep_range"] is False


def test_lightmycells_failed_confirmation_enters_multitask_as_blocked_head():
    baseline = evaluate_multitask_outer_model(*_multitask_inputs())
    report = evaluate_multitask_outer_model(
        *_multitask_inputs(),
        lightmycells_report=_result("lightmycells_confirmation_evaluation.json"),
    )
    assert report["architecture"]["head_count"] == (
        baseline["architecture"]["head_count"] + 1
    )
    head = {item["head_id"]: item for item in report["heads"]}[
        "label_free_organelle_spatial_prediction"
    ]
    assert head["status"] == "blocked_external_spatial_transfer"
    assert head["external_holdout"] is True
    assert head["pristine_confirmation"] is True
    assert head["observable_evidence_eligible"] is False
    assert head["uncertainty_eligible"] is False
    assert head["parameter_proposal_eligible"] is False
    assert head["head_id"] not in report["sweep_interface"]["admissible_head_ids"]
    assert head["head_id"] in report["readiness"][
        "blocked_external_classifier_heads"
    ]


def test_lightmycells_v2_model_passes_selection_before_confirmation_is_opened():
    protocol_path = HERE / "lightmycells_model_training_protocol_v2.json"
    report = _result("lightmycells_model_training_report_v2.json")
    checkpoint = HERE / "results" / "lightmycells_model_checkpoint_v2.pt"
    assert report["protocol_sha256"] == hashlib.sha256(
        protocol_path.read_bytes()
    ).hexdigest()
    assert report["checkpoint_sha256"] == hashlib.sha256(
        checkpoint.read_bytes()
    ).hexdigest()
    assert report["architecture"] == {
        "base_channels": 16,
        "embedding_dimension": 64,
        "family": "conditional_two_level_residual_UNet",
        "parameter_count": 244_001,
        "selected_epoch": 16,
    }
    scores = [
        (row["selection"]["minimum_endpoint_ratio"], row["epoch"])
        for row in report["history"]
    ]
    assert max(scores, key=lambda item: (item[0], -item[1])) == (
        report["selection_endpoints"]["minimum_endpoint_ratio"], 16
    )
    assert report["selection_endpoints"]["all_selection_endpoints_passed"] is True
    assert report["selection_endpoints"]["values"] == {
        "Pearson": 0.558005190481734,
        "global_SSIM": 0.47546305594295896,
        "relative_MAE_improvement": 0.30090364736140096,
    }
    for target in ("Nucleus", "Mitochondria"):
        calibration = report["uncertainty_calibration"][target]
        assert calibration["status"] == "calibrated_same_source_studies"
        assert 0.85 <= calibration["pixel_coverage"] <= 0.95
    assert report["confirmation_open_condition_passed"] is True
    assert report["confirmation_pairs_seen"] == 0
    assert report["production_eligible"] is False
    assert report["observable_evidence_eligible"] is False
    assert report["uncertainty_authority"] is False
    assert report["aleph_authority"] == "none"


def test_lightmycells_v3_protocol_seals_consumed_confirmation():
    protocol = json.loads(
        (HERE / "lightmycells_model_training_protocol_v3.json").read_text()
    )
    consumed = protocol["consumed_confirmation"]
    assert protocol["scope"].startswith("non_disease_label_free")
    assert consumed["may_be_used_for_gradient_updates"] is False
    assert consumed["may_be_used_for_epoch_or_hyperparameter_selection"] is False
    assert consumed["may_be_re_evaluated_by_v3"] is False
    assert protocol["training_tensor"]["confirmation_pairs"] == 0
    assert protocol["optimization"][
        "selection_or_calibration_rows_used_for_gradient_updates"
    ] is False
    assert protocol["authority_contract"] == {
        "v3_success_before_new_confirmation": "training_only_checkpoint_candidate",
        "external_confirmation_authority": False,
        "general_cell_state_production_authority": False,
        "aleph_parameter_authority": "none",
        "may_emit_numeric_sweep_range": False,
    }


def test_allencell_label_free_metadata_audit_selects_raw_normal_hipsc_only():
    protocol_path = HERE / "allencell_label_free_metadata_audit_protocol.json"
    report = _result("allencell_label_free_metadata_audit.json")
    manifest_path = HERE / "results" / "allencell_label_free_candidate_manifest.jsonl"
    rows = [json.loads(line) for line in manifest_path.read_text().splitlines()]
    assert report["protocol_sha256"] == hashlib.sha256(
        protocol_path.read_bytes()
    ).hexdigest()
    assert report["selected_manifest_sha256"] == hashlib.sha256(
        manifest_path.read_bytes()
    ).hexdigest()
    assert report["range_request_count"] == 64
    assert report["range_bytes_read"] == 128 * 1024 * 1024
    assert report["distinct_eligible_FOVs"] >= 32
    assert report["selected_candidate_count"] == len(rows) == 32
    assert len({row["FOVId"] for row in rows}) == 32
    assert {row["structure_name"] for row in rows} == {"TOMM20"}
    assert {row["cell_stage"] for row in rows} == {"M0"}
    assert {row["outlier"] for row in rows} == {"No"}
    assert report["image_headers_accessed"] == 0
    assert report["image_pixels_accessed"] == 0
    assert report["model_predictions_made"] == 0
    assert report["external_confirmation_complete"] is False
    assert report["aleph_parameter_authority"] == "none"


def test_allencell_label_free_v3_improves_but_external_gate_remains_failed():
    protocol_path = HERE / "allencell_label_free_confirmation_protocol.json"
    acquisition = _result("allencell_label_free_acquisition.json")
    report = _result("allencell_label_free_confirmation_evaluation.json")
    assert acquisition["protocol_sha256"] == hashlib.sha256(
        protocol_path.read_bytes()
    ).hexdigest()
    assert acquisition["payload_count"] == 64
    assert acquisition["single_use_confirmation_consumed"] is True
    assert report["preregistered_confirmation"] == {
        "status": "refused_geometry_exact_match",
        "all_endpoints_passed": False,
    }
    diagnostic = report["post_refusal_non_authoritative_diagnostic"]
    assert diagnostic["v2"]["values"] == {
        "Pearson": 0.1306182146147445,
        "global_SSIM": 0.10221794914475538,
        "relative_MAE_improvement": -0.6029442803664952,
    }
    assert diagnostic["v3"]["values"] == {
        "Pearson": 0.3678856266210219,
        "global_SSIM": 0.35902407890358756,
        "relative_MAE_improvement": 0.022969473599104666,
    }
    assert diagnostic["v3_minimum_endpoint_ratio_improvement_over_v2"] > 6.25
    assert diagnostic["v3_pixel_interval_coverage"] > 0.95
    assert diagnostic["v3_embedding_OOD"]["study_macro_AUROC"] < 0.45
    assert diagnostic["all_diagnostic_endpoints_passed"] is False
    assert report["production_eligible"] is False
    assert report["observable_evidence_eligible"] is False
    assert report["uncertainty_authority"] is False
    assert report["aleph_authority"] == "none"
    assert report["may_emit_numeric_sweep_range"] is False


def test_multiprovider_if_mapping_uses_exact_match_precedence():
    protocol = json.loads((HERE / "multiprovider_if_preregistered_protocol.json").read_text())
    ontology = protocol["coarse_ontology_in_order"]
    assert map_author_label("ER", ontology) == ["endoplasmic_reticulum"]
    assert map_author_label("Membrane", ontology) == ["plasma_membrane_or_junction"]
    assert map_author_label("Nuclear membrane", ontology) == ["nuclear_envelope"]
    assert map_author_label("cell_contact", ontology) == [
        "plasma_membrane_or_junction"
    ]


def test_multiprovider_if_preprocessing_center_crops_before_resize():
    rectangular = np.asarray([[10, 0, 0, 10], [10, 0, 0, 10]], dtype=np.float32)
    encoded, constant = normalize_crop_resize(rectangular, 8)
    assert encoded.shape == (8, 8)
    assert encoded.dtype == np.uint8
    assert constant is False
    assert np.all(encoded == 0)


def test_multiprovider_if_development_tensors_are_disjoint_and_non_authoritative():
    protocol_path = HERE / "multiprovider_if_preregistered_protocol.json"
    protocol_sha = hashlib.sha256(protocol_path.read_bytes()).hexdigest()
    assert protocol_sha == "41dc73848e0208dd1b4c7295e2740dd18accb07e26d0f9e2cc1155e59684f77a"
    hpa_report = _result("hpa_multiprovider_tensor_report.json")
    hpa_tensor = HERE / "results" / "hpa_multiprovider_tensor.npz"
    assert hpa_report["protocol_sha256"] == protocol_sha
    assert hpa_report["tensor_sha256"] == hashlib.sha256(hpa_tensor.read_bytes()).hexdigest()
    assert hpa_report["shape"] == [2_858, 2, 128, 128]
    assert hpa_report["split_component_counts"] == {
        "calibration": 384,
        "selection": 381,
        "test": 56,
        "train": 2_037,
    }
    assert hpa_report["component_overlap_across_splits"] == 0
    assert hpa_report["excluded_unmapped_image_count"] == 14
    assert hpa_report["HPA_is_external_confirmation"] is False
    assert hpa_report["may_emit_numeric_sweep_range"] is False
    with np.load(hpa_tensor, allow_pickle=False) as data:
        assert data["images"].shape == (2_858, 2, 128, 128)
        assert data["images"].dtype == np.uint8
        assert data["labels"].shape == (2_858, 9)
        assert len(set(data["component_id"].tolist())) == 2_858

    opencell_report = _result("opencell_multiprovider_tensor_report.json")
    opencell_tensor = HERE / "results" / "opencell_multiprovider_tensor.npz"
    assert opencell_report["protocol_sha256"] == protocol_sha
    assert opencell_report["tensor_sha256"] == hashlib.sha256(
        opencell_tensor.read_bytes()
    ).hexdigest()
    assert opencell_report["shape"] == [2_260, 2, 128, 128]
    assert opencell_report["target_count"] == 1_130
    assert opencell_report["OpenCell_is_pristine_confirmation"] is False
    assert opencell_report["may_emit_numeric_sweep_range"] is False


def test_allencell_confirmation_is_refused_before_pixel_access():
    protocol_path = HERE / "multiprovider_if_preregistered_protocol.json"
    amendment_path = HERE / "multiprovider_if_allencell_metadata_amendment.json"
    amendment = json.loads(amendment_path.read_text())
    report = _result("allencell_multiprovider_metadata_gate.json")
    assert amendment["base_protocol_sha256"] == hashlib.sha256(
        protocol_path.read_bytes()
    ).hexdigest()
    assert report["metadata_amendment_sha256"] == hashlib.sha256(
        amendment_path.read_bytes()
    ).hexdigest()
    assert report["metadata_gate_passed"] is False
    assert report["status"] == "refused_before_Allen_image_access"
    candidate = report["candidate_reports_until_selection"][0]
    assert candidate["resolved_columns"] == {
        "cell_id": "CellId",
        "image_object_path": "crop_raw",
        "tagged_structure": "structure_name",
    }
    assert candidate["missing_semantic_fields"] == ["tagged_gene_or_protein"]
    assert report["Allen_image_objects_downloaded"] == 0
    assert report["Allen_pixels_accessed"] == 0
    assert report["model_predictions_made"] == 0


def test_multiprovider_if_model_is_frozen_but_blocked_by_development_ood():
    protocol_path = HERE / "multiprovider_if_model_training_protocol.json"
    checkpoint = HERE / "results" / "multiprovider_if_model_checkpoint.pt"
    report = _result("multiprovider_if_model_report.json")
    assert report["protocol_sha256"] == hashlib.sha256(
        protocol_path.read_bytes()
    ).hexdigest()
    assert report["checkpoint_sha256"] == hashlib.sha256(
        checkpoint.read_bytes()
    ).hexdigest()
    assert report["architecture"] == {
        "class_count": 9,
        "embedding_dimension": 64,
        "family": "domain_balanced_residual_CNN",
        "parameter_count": 349_321,
    }
    history = [
        (row["selection"]["score"], -row["epoch"]) for row in report["training"]["history"]
    ]
    assert max(history) == (report["training"]["selected_score"], -80)
    assert report["training"]["epochs_executed"] == 80
    assert report["training"]["selected_epoch"] == 80
    assert report["selection"]["HPA"]["macro_average_precision"] > 0.64
    assert report["selection"]["OpenCell"]["macro_average_precision"] > 0.89
    assert report["selection"]["HPA"]["every_supported_class_AP_above_prevalence"] is True
    assert report["selection"]["OpenCell"]["every_supported_class_AP_above_prevalence"] is True
    assert report["development_open_gates"] == {
        "HPA_selection_macro_AP_minimum": True,
        "OpenCell_selection_macro_AP_minimum": True,
        "OpenCell_semantic_OOD_AUROC_minimum": False,
        "every_supported_provider_class_AP_above_prevalence": True,
        "minimum_provider_prevalence_normalized_AP_minimum": True,
    }
    assert report["OOD"]["development_AUROC"] < 0.35
    assert report["all_development_open_gates_passed"] is False
    assert report["next_external_confirmation_pixels_may_be_opened"] is False
    assert report["Allen_pixels_seen"] == 0
    assert report["external_provider_validation_complete"] is False
    assert report["production_eligible"] is False
    assert report["uncertainty_authority"] is False
    assert report["observable_evidence_eligible"] is False
    assert report["Aleph_parameter_authority"] == "none"
    assert report["may_emit_numeric_sweep_range"] is False


def test_multiprovider_if_cross_provider_calibration_remains_limited():
    report = _result("multiprovider_if_model_report.json")
    eligible = {
        name
        for name, value in report["calibration"]["positive_conformal"].items()
        if value["cross_provider_calibration_eligible"]
    }
    assert eligible == {
        "cytoskeleton_or_adhesion",
        "cytosol",
        "endoplasmic_reticulum",
        "golgi",
        "nuclear_interior",
        "plasma_membrane_or_junction",
        "vesicular_organelle",
    }
    assert report["calibration"]["internal_provider_checks"]["HPA"][
        "positive_label_coverage"
    ] >= 0.90
    assert report["calibration"]["internal_provider_checks"]["OpenCell"][
        "positive_label_coverage"
    ] >= 0.90
    assert report["calibration"]["external_uncertainty_authority"] is False


def test_multiprovider_if_ood_failure_is_label_semantic_not_visual_novelty():
    report = _result("multiprovider_if_ood_failure_diagnostic.json")
    assert report["analysis_role"] == "post-freeze_failure_diagnosis_only"
    assert report["may_change_frozen_gate_or_checkpoint"] is False
    assert report["semantic_OOD_high_grade_author_annotation_counts"] == {
        "big_aggregates": 3,
        "nuclear_punctae": 10,
    }
    assert report["semantic_OOD_nearest_frozen_coarse_class_counts"] == {
        "nuclear_interior": 10,
        "vesicular_organelle": 3,
    }
    assert report["semantic_OOD_score_summary"]["median"] < report[
        "in_domain_score_summary"
    ]["median"]
    assert report["reproduced_development_OOD_AUROC"] < 0.35
    assert report["frozen_development_gate_remains_failed"] is True
    assert report["next_external_confirmation_pixels_may_be_opened"] is False
    assert report["production_eligible"] is False
    assert report["Aleph_parameter_authority"] == "none"


def test_multiprovider_if_imaging_ood_v2_passes_without_external_authority():
    protocol_path = HERE / "multiprovider_if_ood_v2_protocol.json"
    head = HERE / "results" / "multiprovider_if_ood_v2_head.npz"
    report = _result("multiprovider_if_ood_v2_report.json")
    assert report["protocol_sha256"] == hashlib.sha256(
        protocol_path.read_bytes()
    ).hexdigest()
    assert report["quality_head_sha256"] == hashlib.sha256(head.read_bytes()).hexdigest()
    assert report["v1_semantic_OOD_gate_remains_failed"] is True
    assert report["semantic_novelty_identifiable_from_image_alone"] is False
    assert report["role_unit_counts"] == {
        "OOD_calibration": {"HPA": 160, "OpenCell": 116},
        "OOD_fit": {"HPA": 512, "OpenCell": 512},
        "OOD_selection": {"HPA": 160, "OpenCell": 132},
    }
    assert report["cross_role_unit_overlap_count"] == 0
    assert report["final_development_test_unit_counts"] == {
        "HPA": 381,
        "OpenCell": 157,
    }
    final = report["final_v2_development_test"]
    assert final["minimum_provider_by_family_AUROC"] > 0.86
    assert final["macro_provider_by_family_AUROC"] > 0.97
    assert final["clean_false_refusal_rate"]["HPA"] < 0.04
    assert final["clean_false_refusal_rate"]["OpenCell"] < 0.02
    assert all(report["development_gates"].values())
    assert report["all_v2_development_gates_passed"] is True
    assert report["replacement_external_provider_may_be_preregistered_and_opened_once"] is True
    assert report["external_provider_pixels_accessed"] == 0
    assert report["external_provider_validation_complete"] is False
    assert report["external_OOD_authority"] is False
    assert report["production_eligible"] is False
    assert report["uncertainty_authority"] is False
    assert report["observable_evidence_eligible"] is False
    assert report["Aleph_parameter_authority"] == "none"
    assert report["may_emit_numeric_sweep_range"] is False


def test_multitask_model_admits_fail_closed_mechanics_protocol_head():
    mechanics = json.loads((
        HERE.parents[2]
        / "outputs/outer/mechanics_protocol/model/mechanics_protocol_report.json"
    ).read_text(encoding="utf-8"))
    baseline = evaluate_multitask_outer_model(*_multitask_inputs())
    report = evaluate_multitask_outer_model(
        *_multitask_inputs(), mechanics_protocol_report=mechanics
    )
    assert report["architecture"]["head_count"] == baseline["architecture"]["head_count"] + 1
    head = next(
        item for item in report["heads"]
        if item["head_id"] == "mechanics_protocol_and_probe_scale_representation"
    )
    assert head["status"] == "training_only_operator_representation_cell_line_transfer_failed"
    assert head["split_contract"]["wire_split_group_overlap_count"] == 0
    assert head["observable_evidence_eligible"] is False
    assert head["parameter_proposal_eligible"] is False
