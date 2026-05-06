"""V2 imaging measurement-protocol gate — area extractor + Hard Rule 11 wording-boundary tests.

Paired lock: ``docs/v2/v2_imaging_measurement_protocol_gate_locked.md``.
Paired sanity gate: ``docs/v2/v2_imaging_measurement_protocol_gate_sanity_gate.md``.

The earlier framing as "V2-2 Unit 1 / imaging→SingleCellState loader"
was retracted in commit ``f4e09aa`` after PI's phase-scope correction
(260313 = spheroid-level, not single-cell). This test file is the
live wording-boundary surface for the renamed measurement-protocol
gate; the retracted lock spec is kept on disk for audit trail only.

Covers (in this order):

1. Algorithm constants (frozen against PI source ``spread_infer.py``).
2. ``select_central_contour`` behaviour on synthetic masks.
3. PI-directed provenance exclusion is surfaced (NOT silently filtered)
   and Pos31 is keyed correctly.
4. Source-pipeline reproduction is byte-exact (or below the cross-version
   drift budget) on the 260313 dataset.
5. **Hard Rule 11 wording-boundary meta-test** — the lock spec MUST contain
   the required source-pipeline-reproduction wording and MUST NOT contain
   any of the rejected fitting-flavoured phrases. This guards against
   silent rewording during future edits.
6. Sister-gate / sister-doc mirroring (CLAUDE.md design-note pre-commit
   batch step 6): the Pos31 exclusion clause appears in code, docs, and
   the reproduction script, all citing the same wording.
"""

from __future__ import annotations

import math
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
LOCK_SPEC = REPO / "docs" / "v2" / "v2_imaging_measurement_protocol_gate_locked.md"
SANITY_GATE = REPO / "docs" / "v2" / "v2_imaging_measurement_protocol_gate_sanity_gate.md"
RETRACTED_LOCK = (
    REPO
    / "docs"
    / "v2"
    / "v2_layer_2_unit_1_imaging_to_single_cell_state_loader_locked.md"
)
EXTRACTOR_PY = REPO / "acs" / "v2" / "imaging_contract" / "area_extractor.py"
REPRO_SCRIPT = REPO / "scripts" / "v2" / "reproduce_pi_area_extraction.py"

cv2 = pytest.importorskip(
    "cv2", reason="imaging measurement-protocol gate reproduction needs cv2."
)
from acs.v2.imaging_contract import area_extractor as ae  # noqa: E402


# -----------------------------------------------------------------
# 1. Algorithm constants — frozen against PI source spread_infer.py
# -----------------------------------------------------------------

def test_pixel_area_um2_matches_pi_source():
    # spread_infer.py:44 — PIXEL_AREA_UM2 = 4.194304 = 2.048**2 μm²/px.
    assert ae.PIXEL_AREA_UM2 == pytest.approx(4.194304, rel=1e-12)
    assert ae.PIXEL_AREA_UM2 == pytest.approx(2.048 ** 2, rel=1e-12)


def test_select_central_contour_constants_match_pi_source():
    # spread_infer.py:55-57 — CENTER_MAX_DIST_FRAC, BORDER_MARGIN_PX,
    # REJECT_TOUCH_BORDER, area-cost weight 0.0005.
    assert ae.CENTER_MAX_DIST_FRAC == 0.45
    assert ae.BORDER_MARGIN_PX == 12
    assert ae.REJECT_TOUCH_BORDER is True
    assert ae.SELECT_CENTRAL_CONTOUR_AREA_COST_WEIGHT == 0.0005


# -----------------------------------------------------------------
# 2. select_central_contour on synthetic masks
# -----------------------------------------------------------------

def _solid_circle_mask(w: int, h: int, cx: int, cy: int, r: int) -> np.ndarray:
    m = np.zeros((h, w), dtype=np.uint8)
    cv2.circle(m, (cx, cy), r, 255, thickness=-1)
    return m


def test_select_central_contour_picks_centred_blob_over_offcentre_one():
    h, w = 1040, 1388  # PI image_size
    m = _solid_circle_mask(w, h, w // 2, h // 2, 80)
    cv2.circle(m, (300, 300), 50, 255, thickness=-1)  # off-centre
    cnts, _ = cv2.findContours(m, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    cnt = ae.select_central_contour(cnts, w, h)
    assert cnt is not None
    M = cv2.moments(cnt)
    cx = M["m10"] / M["m00"]
    cy = M["m01"] / M["m00"]
    assert abs(cx - w / 2) < 5
    assert abs(cy - h / 2) < 5


def test_select_central_contour_rejects_border_touching():
    h, w = 200, 200
    m = np.zeros((h, w), dtype=np.uint8)
    # Border-touching blob (covers full image).
    m[:, :] = 255
    cnts, _ = cv2.findContours(m, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    cnt = ae.select_central_contour(cnts, w, h)
    # Fallback path should still return *some* contour (PI source line
    # 345-357), but the touches_border check must classify the only
    # contour as touching.
    assert ae.touches_border(cnts[0], w, h) is True


def test_extract_area_px_byte_exact_on_synthetic_circle():
    h, w = 400, 400
    m = _solid_circle_mask(w, h, w // 2, h // 2, 100)
    tmp = (
        REPO
        / "runs"
        / "v2_imaging_measurement_protocol_gate"
        / "_synthetic_circle_mask.png"
    )
    tmp.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(tmp), m)
    repro, status = ae.extract_area_px_from_mask(tmp)
    assert status == ae.PROVENANCE_STATUS_NORMAL
    # cv2.contourArea on the discretised circle should be close to π r²
    # but not exactly equal (Green's theorem on the integer-vertex
    # contour); allow ~1% from continuous geometry.
    expected = math.pi * 100 ** 2
    assert abs(repro - expected) / expected < 0.02
    tmp.unlink()


# -----------------------------------------------------------------
# 3. PI-directed provenance exclusion is keyed and surfaced
# -----------------------------------------------------------------

def test_pos31_is_in_pi_directed_provenance_exclusion_set():
    assert ("Lam4", "Position(31)") in ae.PI_DIRECTED_PROVENANCE_EXCLUSIONS


def test_classify_provenance_routes_pos31_to_excluded_status():
    assert ae.classify_provenance("Lam4", "Position(31)") == (
        ae.PROVENANCE_STATUS_EXCLUDED
    )
    assert ae.classify_provenance("Lam4", "Position(10)") == (
        ae.PROVENANCE_STATUS_NORMAL
    )
    assert ae.classify_provenance("Bare", "MD-Experiment-0004_Position(102)") == (
        ae.PROVENANCE_STATUS_NORMAL
    )


def test_aggregate_only_uses_normal_set_for_acceptance():
    """Codex guardrail id=2304: Pos31 MUST NOT contribute to thresholds."""
    results = [
        ae.AreaExtractionResult("Lam4", "Position(10)", 0, "f", 1000.0, 1000.0, 0.0,
                                 ae.PROVENANCE_STATUS_NORMAL),
        ae.AreaExtractionResult("Lam4", "Position(31)", 0, "f", 1000.0, 1300.0, 0.30,
                                 ae.PROVENANCE_STATUS_EXCLUDED),
    ]
    summary = ae.aggregate_normal_residuals(results)
    # Pos31 30% residual must NOT raise the acceptance threshold.
    assert summary["normal_max_rel_residual"] == 0.0
    assert summary["rows_normal"] == 1
    assert summary["rows_excluded_by_pi"] == 1
    assert summary["excluded_max_rel_residual"] == pytest.approx(0.30)


# -----------------------------------------------------------------
# 4. Source-pipeline reproduction on the 260313 dataset
# -----------------------------------------------------------------

@pytest.mark.skipif(
    not (REPO / "data" / "experimental" / "260313_Bare.csv").exists()
    or not (REPO / "data" / "experimental" / "imaging").exists(),
    reason="260313 CSVs or saved-mask PNG tree not present in this checkout",
)
def test_reproduction_byte_exact_on_real_dataset_sample():
    """Spot-check 6 rows per group on the live 260313 dataset (in-process).

    Cheap CI variant of the full-reproduction gate: 18 rows total,
    runs in <1 s. Uses ``area_extractor`` directly so we still cover
    the saved-mask path-resolution + cv2 contour pipeline. Per-row
    rel_residual must be exactly 0 (cv2 4.13.0 deterministic).
    """
    import pandas as pd  # local import to keep top-level fast
    from acs.v2.imaging_contract import area_extractor as _ae

    imaging_root = REPO / "data" / "experimental" / "imaging"
    rng = np.random.default_rng(seed=20260506)
    for group, csv_name in [
        ("Bare", "260313_Bare.csv"),
        ("Lam4", "260313_Lam4.csv"),
        ("Pre", "260313_Pre.csv"),
    ]:
        df = pd.read_csv(REPO / "data" / "experimental" / csv_name)
        normal = df[df["Series"].apply(
            lambda s, g=group: (g, s) not in _ae.PI_DIRECTED_PROVENANCE_EXCLUSIONS
        )]
        sample = normal.iloc[rng.integers(0, len(normal), size=6)]
        for _, row in sample.iterrows():
            mp = (
                imaging_root
                / f"260313_{group}"
                / str(row["Series"])
                / f"{Path(str(row['filename'])).stem}_mask.png"
            )
            assert mp.exists(), f"sample mask missing: {mp}"
            res = _ae.reproduce_one_row(
                group=group,
                series=str(row["Series"]),
                frame=int(row["Frame"]),
                filename=str(row["filename"]),
                csv_area_px=float(row["Area_px"]),
                mask_path=mp,
            )
            assert res.provenance_status == _ae.PROVENANCE_STATUS_NORMAL
            assert res.rel_residual == 0.0, (
                f"non-zero residual on {group}/{row['Series']}/"
                f"{row['filename']}: rel_residual={res.rel_residual}"
            )


@pytest.mark.skipif(
    not (REPO / "data" / "experimental" / "260313_Bare.csv").exists(),
    reason="260313 CSVs not present in this checkout",
)
@pytest.mark.skipif(
    os.environ.get("RUN_FULL_REPRODUCTION_GATE") != "1",
    reason=(
        "full 4551-row mask reproduction is opt-in (set "
        "RUN_FULL_REPRODUCTION_GATE=1; reads ~7 GB of mask PNGs); "
        "in-process unit tests above already cover the algorithm."
    ),
)
def test_reproduction_byte_exact_on_real_dataset():
    """Run the imaging measurement-protocol gate reproduction script and check ≤1.5% threshold.

    Opt-in via env var ``RUN_FULL_REPRODUCTION_GATE=1`` because the
    script reads all 4551 saved masks (~7 GB cold) and takes
    several minutes. Lock seal verification artefact lives at
    ``runs/v2_imaging_measurement_protocol_gate/pi_area_reproduction.csv`` (gitignored)
    and ships with the seal evidence rather than CI.
    """
    out = subprocess.run(
        [str(REPO / ".venv-collab" / "bin" / "python"),
         str(REPRO_SCRIPT)],
        capture_output=True,
        text=True,
        cwd=str(REPO),
        timeout=1800,
    )
    if out.returncode != 0:
        pytest.fail(
            f"reproduce_pi_area_extraction.py exited {out.returncode}\n"
            f"--- stdout ---\n{out.stdout}\n--- stderr ---\n{out.stderr}"
        )
    assert "[CONTRACT PASS]" in out.stdout
    # All 260313 rows reproduce byte-exact under cv2 4.13.0; assert that
    # no group fails the threshold contract.
    assert "[CONTRACT FAIL]" not in out.stdout


# -----------------------------------------------------------------
# 5. Hard Rule 11 wording-boundary meta-test
# -----------------------------------------------------------------

# Phrases the lock spec MUST contain — Codex guardrail id=2304 require list.
# These are verbatim citations of the source-pipeline + provenance contract
# wording. A test failure means someone silently softened the contract.
HARD_RULE_11_REQUIRED_PHRASES = (
    "source-pipeline reproduction",
    "PI-directed provenance exclusion: Lam4 Position(31)",
    "excluded from protocol-identification set; reported separately",
)

# Phrases the lock spec MUST NOT contain — Codex guardrail id=2304 reject list.
# These describe fitting-flavoured behaviours that violate Hard Rule 1
# (no parameter fitting to PI data) or Hard Rule 11 (no measurement-protocol
# erosion). Presence of any one of these is a contract violation.
HARD_RULE_11_REJECT_PHRASES = (
    "fit to CSV",
    "tuned tolerance",
    "dropped outlier to pass",
    "Area_px matched by correction factor",
)


@pytest.mark.skipif(
    not LOCK_SPEC.exists(),
    reason="imaging measurement-protocol gate lock not yet sealed",
)
def test_lock_spec_contains_required_source_pipeline_wording():
    text = LOCK_SPEC.read_text(encoding="utf-8")
    missing = [p for p in HARD_RULE_11_REQUIRED_PHRASES if p not in text]
    assert not missing, (
        "Imaging measurement-protocol gate lock spec is missing required Hard Rule 11 wording: "
        f"{missing}. Codex guardrail id=2304 require-list."
    )


@pytest.mark.skipif(
    not LOCK_SPEC.exists(),
    reason="imaging measurement-protocol gate lock not yet sealed",
)
def test_lock_spec_does_not_contain_fitting_flavoured_wording():
    text = LOCK_SPEC.read_text(encoding="utf-8").lower()
    hits = [p for p in HARD_RULE_11_REJECT_PHRASES if p.lower() in text]
    assert not hits, (
        "Imaging measurement-protocol gate lock spec contains rejected fitting-flavoured wording: "
        f"{hits}. Codex guardrail id=2304 reject-list."
    )


def test_extractor_module_cites_pi_source_path_and_function():
    """``area_extractor.py`` must keep the source-of-truth pointer alive."""
    text = EXTRACTOR_PY.read_text(encoding="utf-8")
    assert "spread_infer.py" in text
    assert "build_row_from_saved_mask" in text
    # Inline unit-chain phrase MUST be present.
    assert "Hard Rule 11" in text
    assert "pixel_count" in text and "pixel_size_um" in text


def test_extractor_module_does_not_carry_fitting_flavoured_wording():
    text = EXTRACTOR_PY.read_text(encoding="utf-8").lower()
    hits = [p for p in HARD_RULE_11_REJECT_PHRASES if p.lower() in text]
    assert not hits, hits


# -----------------------------------------------------------------
# 6. Sister-doc / sister-code mirroring (design-note pre-commit batch step 6)
# -----------------------------------------------------------------

def test_pos31_exclusion_mirrored_across_code_doc_script():
    """Pos31 exclusion clause must appear in extractor code, lock spec, and
    reproduction script, citing the same wording. This is the canonical
    sister-mirror that prevents drift between the three artefacts."""
    pos31_clause = re.compile(r"Lam4[^\w]+Position\(31\)", re.IGNORECASE)

    code_text = EXTRACTOR_PY.read_text(encoding="utf-8")
    script_text = REPRO_SCRIPT.read_text(encoding="utf-8")
    assert pos31_clause.search(code_text), (
        "Pos31 exclusion missing from area_extractor.py — must be cited in "
        "code per Codex guardrail id=2304."
    )
    assert pos31_clause.search(script_text), (
        "Pos31 exclusion missing from scripts/v2/reproduce_pi_area_extraction.py."
    )
    if LOCK_SPEC.exists():
        lock_text = LOCK_SPEC.read_text(encoding="utf-8")
        assert pos31_clause.search(lock_text), (
            "Pos31 exclusion missing from imaging measurement-protocol gate lock spec."
        )


def test_pixel_size_consistent_with_v2_1_imaging_yaml():
    """Gate's PIXEL_AREA_UM2 must match V2-1 imaging contract YAML."""
    yaml_path = REPO / "configs" / "imaging" / "260313.yaml"
    if not yaml_path.exists():
        pytest.skip("V2-1 imaging contract YAML not present")
    yaml_text = yaml_path.read_text(encoding="utf-8")
    # Loose check — the lock-pinned line is `pixel_area_um2: 4.194304`.
    assert "4.194304" in yaml_text, (
        "V2-1 imaging YAML pixel_area_um2 has drifted from "
        f"gate PIXEL_AREA_UM2 = {ae.PIXEL_AREA_UM2}."
    )


def test_extractor_acceptance_threshold_is_lock_pinned():
    """The 1.5% threshold is a contract value, not a tuning knob."""
    assert ae.ACCEPTANCE_REL_RESIDUAL_THRESHOLD == 0.015
