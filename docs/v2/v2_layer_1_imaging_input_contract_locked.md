# V2 Layer 1 Imaging Input Contract — Locked Design

**Date**: 2026-05-06 KST
**Authors**: Claude + Codex design-discussion (3-round adversarial lock,
PI id=809 + id=1985 aggressive debate posture)
**Source unit**: design-discussion `topic=v2-layer-1-imaging-input-contract`,
MCP id 2237–2252 (rounds 1–3 + seal ack)
**PI directives**:
- `id=2187` Option 1 pivot — vertical 5-layer framework completion
- `id=2188` PI A안 imaging input choice
- `id=2190` "260313만 사용해" (scope to single date)
- `id=2232` PI Q1=(a) raw imaging stack 별도 보관, Q3=(a) Mac-only V2-1/V2-2,
  Q4=(b) 추가 imaging 가능
- `id=2235` "동의" PI Q5 = 70/30 condition-stratified scene-level + seed=20260506
- `id=1985` "건강하고 활발한 토론"
- `id=2114` external audit: framework-completion-vs-plumbing discipline
**Tier**: **A-tier full debate** (3 rounds + 17 substantive Y trace).
**First Layer-framework-completion deliverable**: this is V2-1 of
`docs/v2/10_dev_roadmap_v2.md`. Subsequent V2-2 Layer 1/2 framework
completion design cycles depend on this contract.
**PI ratify status**: full delegation per PI id=939/1008/2235. impl-work
uses this for the V2-1 contract implementation.

---

## 0. Scope + framework-completion framing

### What this unit IS

- **V2-1 imaging input contract** for PI 260313 single date dataset
- Canonical file layout, required metadata, metric registry, contract
  validation tests, calibration/validation split policy
- Pure read-only loader contract — does NOT mutate PI data
- **First grounding** of v2's "image-constrained" identity in real PI
  data (audit `id=2114` plumbing-vs-mechanism gap closure)

### What this unit is NOT

- ❌ NOT a v2 dynamics parameter fitting pipeline (Hard Rule 1; PI
  data is calibration/validation reference only, NEVER physics
  parameter input)
- ❌ NOT a Cellpose/SAM segmentation pipeline (PI segmentation is
  authoritative input — `model_choice=saved, used_thresh=saved`)
- ❌ NOT raw imaging stack consumer (PI Q1=(a): raw stack 별도 보관,
  V2-2 후반 또는 V2-3에서 별 contract)
- ❌ NOT a multi-channel pipeline (PI Q4=(b): future channels
  possible, not in V2-1)
- ❌ NOT a Win box / GPU pipeline (PI Q3=(a): Mac-only V2-1/V2-2)
- ❌ NOT cross-condition / cross-date validation (Lam1 from 260126
  dropped per PI id=2190; only 260313 Lam4)
- ❌ NOT a frame-level data split (Y2 P0-2: same-cell leakage risk;
  scene-level only)

### Wording boundary (audit `id=2114` framework-completion discipline)

This is a **framework-completion deliverable**, NOT plumbing scale-up.
First V2-1 deliverable per `docs/v2/10_dev_roadmap_v2.md`. Lock seal is
"V2-1 imaging input contract design lock sealed" (unit-level), with
Layer 1 framework completion still pending V2-2 deliverables.

PI-facing milestone wording: **"V2-1 imaging input contract design lock
sealed (V2-1 of V2-1/V2-2 Layer 1 framework completion path)"**. NOT
"Layer 1 milestone".

---

## 1. Final Lock Summary

### YAML contract spec (Y1-Y17 incorporated)

```yaml
# configs/imaging/260313.yaml
contract_id: pi_260313_v1
contract_version: 1
roadmap_layer: v2_1_imaging_input
canonical_root: data/experimental/imaging/

groups:
  260313_Bare:
    condition: Bare                     # Y9: explicit, NOT derived from group_id
    imaging_date: 2026-03-13
    path: 260313_Bare/
    csv: 260313_Bare.csv
    naming_convention: "MD-Experiment-0004_Position(NN)"
    n_scenes: 8
  260313_Lam4:
    condition: Lam4
    imaging_date: 2026-03-13
    path: 260313_Lam4/
    csv: 260313_Lam4.csv
    naming_convention: "Position(NN)"   # Lam4 distinct (no MD- prefix)
    n_scenes: 26
  260313_Pre:
    condition: Pre
    imaging_date: 2026-03-13
    path: 260313_Pre/
    csv: 260313_Pre.csv
    naming_convention: "MD-Experiment-0004_Position(NN)"
    n_scenes: 25

# Y1 P0-1 critical: pixel size 2.048 μm/px (NOT 1.024); audit estimate
# was 4× wrong. Validated against CSV via sqrt(Area_um2 / Area_px).
pixel_size_um_xy: [2.048, 2.048]
pixel_area_um2: 4.194304
pixel_size_validation:                  # Y12
  method: sqrt_area_um2_over_area_px
  rtol: 1.0e-6
  atol_um: 1.0e-9

# Y16: frame_interval_min is first-class top-level field; the test
# "missing frame interval raises" is unambiguous because of this.
frame_interval_min: 60

# Y11 C3 critical: per-scene frame count variation; Bare 8×83=664
# expected vs actual 602 reveals not all scenes have full 83 frames.
# Use min/max + monotonic time validation (NOT exact equality).
# (Frame-count semantics = Python `len(frames_in_scene)`, NOT `max(Frame)`;
#  see lock amendment 2026-05-06 — original lock had `max(Frame)` = 82/83/82
#  off-by-one, impl-work id=2268 caught this BEFORE first execution.)
expected_timepoints:
  min_frames_per_scene: 1
  max_frames_per_scene_per_group:
    260313_Bare: 83
    260313_Lam4: 84
    260313_Pre: 83
  expected_frame_interval_min: 60       # consistency check vs top-level
  require_monotonic_time_min: true
  require_monotonic_frame: true

# Y17: image_size_px axis order explicit (NumPy/PIL convention can swap)
image_size_px_width_height: [1388, 1040]
image_dtype: uint8

channels:                               # Y6: typed (not just count)
  - name: segmentation_mask
    artifact: segmentation_mask
    dtype: uint8
    source: pi_inhouse_segmentation_output

metric_registry:
  - name: Area_um2
    csv_column: Area_um2
    units: um2
    csv_columns_required: [Area_um2]
  - name: A_over_A0
    csv_column: A_over_A0
    units: dimensionless
    csv_columns_required: [A_over_A0]
  - name: Solidity
    csv_column: Solidity
    units: dimensionless
    csv_columns_required: [Solidity]
  - name: Circularity
    csv_column: Circularity
    units: dimensionless
    csv_columns_required: [Circularity]
  - name: EffectiveRadius_um
    csv_column: EffectiveRadius_um
    units: um
    csv_columns_required: [EffectiveRadius_um]
  - name: Perimeter_um                   # Y6 P0-1 sister: Perimeter_px * pixel_size_um
    csv_column: Perimeter_px
    units: um
    csv_columns_required: [Perimeter_px]
    transform: multiply_by_pixel_size_um  # NOT area-derived

# Y2 P0-2 critical: CSV `filename` is authoritative; NEVER derive
# t-index from Frame index (Bare uses t002/t005/t008 stride 3,
# Lam4 mixed). Loader uses CSV column to resolve mask paths.
filename_resolution:
  csv_authoritative: true
  mask_suffix: "_mask.png"
  overlay_suffix: "_overlay.png"
  overlay_required: false                # diagnostic only

calibration_validation_split:           # Y14: pipeline-only role clarity
  policy: condition_stratified_scene_level
  ratio: [0.7, 0.3]
  seed: 20260506
  manifest: configs/imaging/generated/260313_split_seed20260506.json
  calibration_use: pipeline_loader_metric_dashboard_development_only
  forbidden_use: physics_parameter_fitting

sidecars:                               # Y7: machine-readable
  model2_framewise:
    Bare: _MODEL2_FRAMEWISE_260313_Bare/
    Lam4: _MODEL2_FRAMEWISE_260313_Lam4/
    Pre: _MODEL2_FRAMEWISE_260313_Pre/
    allowed_use: reduced_surrogate_reference_only
    forbidden_use: dynamics_parameter_fitting

provenance:
  pi_decisions:
    raw_imaging_available: true          # PI Q1=(a)
    raw_path: null                       # PI provides later
    sync_policy: mac_v2_1_v2_2__win_sync_v2_3_plus  # PI Q3=(a)
    future_channels_possible: true       # PI Q4=(b)
    split_decision: condition_stratified_70_30  # PI Q5 동의
  imaging_date: 2026-03-13
  segmentation_provenance: "PI in-house, model_choice=saved, used_thresh=saved"
  csv_authoritative: true                # Y2 P0-2 reminder
```

### Module structure

```
acs/v2/
├── data_contract.py                    # EXISTING — reused (Y3)
│                                       #   ImagingDatasetSpec, MetricSpec, V2DataContract,
│                                       #   ArtifactKind validators
├── imaging_contract/                   # NEW (Y3)
│   ├── __init__.py                     # public exports
│   ├── loader.py                       # load_imaging_contract(path) -> ImagingInputContract
│   │                                   # validate_imaging_contract(contract) -> ImagingInputContractValidationResult
│   ├── split_builder.py                # generate 70/30 manifest (Y8/Y13)
│   │                                   # CLI: python -m acs.v2.imaging_contract.split_builder configs/imaging/260313.yaml
│   └── path_resolver.py                # resolve_mask_path / resolve_overlay_path (Y10)
│                                       # naming convention per-group dispatch
configs/imaging/
├── 260313.yaml                         # NEW (this contract)
└── generated/
    └── 260313_split_seed20260506.json  # NEW, generated by split_builder, committed,
                                        # manually edit FORBIDDEN (Y13)
data/experimental/imaging/              # symlinks (gitignored, already created)
├── 260313_Bare/  -> /Users/sw1/Desktop/Spread_data/.../spread_output_260313_Bare/
├── 260313_Lam4/  -> ...
└── 260313_Pre/   -> ...
tests/
└── test_v2_imaging_input_contract.py   # NEW, 13 tests per §4
```

### Loader API (Y5 two-level)

```python
from pathlib import Path
from typing import Optional

from acs.v2.data_contract import V2DataContract  # reuse existing


@dataclass(frozen=True, slots=True)
class ImagingInputContract:
    """Loaded V2-1 imaging input contract from YAML."""
    contract_id: str
    contract_version: int
    canonical_root: Path
    groups: dict[str, GroupSpec]
    pixel_size_um_xy: tuple[float, float]
    frame_interval_min: float                       # Y16 first-class
    image_size_px_width_height: tuple[int, int]     # Y17 explicit
    channels: tuple[ChannelSpec, ...]
    metric_registry: tuple[MetricSpec, ...]
    filename_resolution: FilenameResolutionSpec
    calibration_validation_split: SplitSpec
    sidecars: SidecarSpec
    provenance: ProvenanceSpec


@dataclass(frozen=True, slots=True)
class ImagingInputContractValidationResult:
    """Typed validation result (Y5: PI-facing 'what was checked')."""
    contract_id: str
    n_groups_validated: int
    n_scenes_validated: int
    n_csv_rows_validated: int
    n_mask_paths_resolved: int
    n_overlay_paths_resolved: int                   # diagnostic
    pixel_size_consistency_max_deviation_um: float
    split_manifest_validated: bool
    sidecar_paths_validated: bool
    all_passed: bool
    failure_messages: tuple[str, ...]               # empty if all_passed


def load_imaging_contract(path: Path) -> ImagingInputContract:
    """Load YAML + light validation. Raises on missing required fields."""
    ...


def validate_imaging_contract(
    contract: ImagingInputContract,
) -> ImagingInputContractValidationResult:
    """Heavy validation: 13 fail-closed tests per §4. Returns typed result."""
    ...
```

### Filename resolver (Y10 explicit)

```python
def resolve_mask_path(group_root: Path, csv_row: dict) -> Path:
    """Y10 explicit resolver — CSV authoritative; NEVER derive t-index from Frame."""
    series_dir = group_root / csv_row["Series"]
    stem = Path(csv_row["filename"]).stem
    return series_dir / f"{stem}_mask.png"


def resolve_overlay_path(group_root: Path, csv_row: dict) -> Optional[Path]:
    """Optional/diagnostic. Returns None if missing (NOT error)."""
    series_dir = group_root / csv_row["Series"]
    stem = Path(csv_row["filename"]).stem
    candidate = series_dir / f"{stem}_overlay.png"
    return candidate if candidate.exists() else None
```

### Split builder (Y13 generated, manually edit forbidden)

```python
"""acs/v2/imaging_contract/split_builder.py

CLI: python -m acs.v2.imaging_contract.split_builder configs/imaging/260313.yaml

Generates configs/imaging/generated/<contract_id>_split_seed<N>.json
deterministically from YAML. Manifest is committed; manual edit FORBIDDEN.
Regenerate after any YAML change.
"""

def build_split_manifest(yaml_path: Path) -> dict:
    """Generate split manifest from YAML.

    Output structure:
    {
      "contract_id": "pi_260313_v1",
      "seed": 20260506,
      "policy": "condition_stratified_scene_level",
      "manifest_version": 1,
      "generator": "acs.v2.imaging_contract.split_builder",
      "generated_at_iso": null,                # NOT timestamped (Y13: HB#5/Item5 sister, reproducibility)
      "input_yaml": "configs/imaging/260313.yaml",
      "input_yaml_sha256": "...",              # tamper detection
      "groups": {
        "260313_Bare": {
          "calibration": ["MD-Experiment-0004_Position(94)", ...],
          "validation": ["MD-Experiment-0004_Position(118)", ...]
        },
        ...
      }
    }

    Deterministic: same YAML + same seed → exact same manifest bytes.
    """
    ...
```

### Forbidden in V2-1 contract

- ❌ **Frame-level split** (Y2/PI directive `id=2233/2235`: same-cell leakage risk)
- ❌ **Pixel size 1.024 μm/px** (Y1 P0-1: audit estimate wrong; actual 2.048)
- ❌ **Derive t-index from Frame index** (Y2 P0-2: CSV `filename` authoritative;
  Bare uses stride 3, Lam4 mixed)
- ❌ **Group-id string-split for condition** (Y9: explicit `condition` field)
- ❌ **Exact `frames_per_scene` enforcement** (Y11 C3: Bare 8×83=664 vs actual
  602 → impossible constraint; use min/max + monotonic)
- ❌ **Manual edit of generated split manifest** (Y13: regenerate via
  `split_builder` only; `input_yaml_sha256` enforces consistency)
- ❌ **PI data physics parameter fitting** (Y14 + Hard Rule 1; calibration role
  is `pipeline_loader_metric_dashboard_development_only`)
- ❌ **`_MODEL2_FRAMEWISE_<group>/` sidecar consumed for v2 dynamics parameter
  fitting** (Y7: `forbidden_use: dynamics_parameter_fitting`; allowed only as
  `reduced_surrogate_reference_only` for V2-6 future)
- ❌ **Silent t-index swap or naming convention assumption** (Y10: explicit
  resolver per group; Bare/Pre = `MD-Experiment-0004_Position(NN)`,
  Lam4 = `Position(NN)`)
- ❌ **Image axis order ambiguity** (Y17: `image_size_px_width_height` explicit)
- ❌ **Frame interval inconsistency** (Y16: top-level `frame_interval_min`
  + `expected_timepoints.expected_frame_interval_min` consistency)
- ❌ Cross-condition (Lam1 vs Lam4) cross-validation (PI id=2190: 260126
  dropped; only Lam4 from 260313)
- ❌ Cross-date validation (PI id=2190: 260313 single date only)

---

## 2. Reasoned-acceptance trace (Y1–Y17)

The lock converged after 3 rounds (id 2237–2252) under PI directive
id=1985 aggressive debate posture, with audit `id=2114`/`id=2210`/`id=2231`
direct CSV/disk inspection driving P0 corrections.

**Y1 (P0-1 CRITICAL pixel size correction)**: Audit estimate
`voxel_size_um: 1.024` was wrong. Codex direct CSV inspection
(`Area_um2/Area_px = 162214.7072 / 38675.0 = 4.194304 μm²/px →
sqrt = 2.048 μm/px`) caught the 4× error before lock. All derived
metrics (Perimeter conversion, dashboard scales) would have been
4× corrupted.

**Y2 (P0-2 CRITICAL CSV authoritative, NOT derive t-index from Frame)**:
Audit estimate `t001/t002/...` stride 1 was wrong. Bare uses
`t002/t005/t008` stride 3; Lam4 mixed. Loader uses CSV `Series` +
`filename` columns directly; never invents t-index from Frame index.

**Y3 (Q1 reuse `acs/v2/data_contract.py` + new `imaging_contract/`
package)**: Existing types (`ImagingDatasetSpec`, `MetricSpec`,
`V2DataContract`, `ArtifactKind`) reused; PI-260313-specific loader
in new package. Avoids duplicate contract types.

**Y4 (Q2 YAML at `configs/imaging/260313.yaml`)**: Not in
`data/experimental/imaging/` (gitignored); under `configs/`
sister-pattern. User-facing YAML + typed Python loader.

**Y5 (Q3 two-level loader API)**: `load_imaging_contract` +
`validate_imaging_contract` (typed result). `validate_*` returns typed
result for PI-facing reports ("what was checked"); `load_*` calls
validate by default.

**Y6 (Q4 12 validation tests)**: 5 roadmap fail-closed (missing
voxel/frame interval/channels/empty metrics/cal-val overlap) + 7
additions (group scene counts match CSV, row count consistency, Series
dirs exist, CSV filename → mask resolution, pixel size CSV-derived
cross-check, Perimeter conversion via pixel size, sidecar marker).

**Y7 (Q5 machine-readable sidecar handling)**: `_MODEL2_FRAMEWISE_<group>/`
sidecars exposed as metadata with `allowed_use: reduced_surrogate_reference_only`,
`forbidden_use: dynamics_parameter_fitting`. Loader does NOT load
model2 parameters into main metric registry. No AST guard in first lock
(brittle).

**Y8 (Q6 NEW materialized split manifest)**: 70/30 split generated by
`split_builder.py` from YAML, materialized to
`configs/imaging/generated/260313_split_seed20260506.json`. Committed
(small file). Prevents accidental changes from directory ordering / RNG
implementation drift.

**Y9 (C1 separate `condition` from `group_id`)**: YAML `condition` field
explicit; not derived by string-splitting `group_id`. Future multi-date
or non-PI dataset support natural.

**Y10 (C2 explicit mask/overlay filename resolver)**: `resolve_mask_path`
+ `resolve_overlay_path` defined as `series_dir / f"{stem}_<suffix>"`
with stem from CSV `filename` column. Mask required (test 9
fail-closed); overlay diagnostic-only.

**Y11 (C3 CRITICAL per-scene frame count variation)**: Codex direct
CSV inspection: Bare 8×83=664 expected vs actual 602 → not all scenes
have full 83 frames. YAML uses `min_frames_per_scene: 1` +
`max_frames_per_scene_per_group` + monotonic Frame/Time_min validation.
My initial YAML would have locked impossible constraint.

**Y11 lock amendment 2026-05-06 03:25 KST (impl-work id=2268)**: original
lock max values 82/83/82 used `max(Frame)` (last index) instead of Python
`len(frames_in_scene)` (count). Actual ground-truth `len()` is 83/84/83.
impl-work caught this before first execution per CLAUDE.md gate-edit
prohibition (refused to silently bump YAML to clear test). Amendment
shifts integers only; field name + Python loader semantics unchanged.
Off-by-one is uniform across all 3 conditions, confirming audit-time
counting error not per-condition data anomaly.

**Y12 (C4 pixel-size tolerance unit-specific)**: `pixel_size_validation`
declares `rtol: 1.0e-6` + `atol_um: 1.0e-9`. Loader compares CSV-derived
`sqrt(Area_um2/Area_px)` to declared `pixel_size_um_xy` row-wise via
`np.isclose`. Rounding variation handled via lock amendment, not inline
relaxation.

**Y13 (C5 split manifest = generated output, manually edit FORBIDDEN)**:
Manifest is OUTPUT of `split_builder.py` from YAML, NOT hand-authored.
`input_yaml_sha256` field for tamper detection. Lock §1 mandates
"manually edit FORBIDDEN — regenerate via CLI". Reproducibility
guaranteed across runs / hosts.

**Y14 (C6 calibration role clarity)**: `calibration_use:
pipeline_loader_metric_dashboard_development_only` + `forbidden_use:
physics_parameter_fitting`. Hard Rule 1 strengthened in YAML metadata.
"Calibration" name is standard ML/statistics; this contract restricts
scope explicitly.

**Y15 (C7 test 7 corrected + test 13 NEW)**: Test 7 no longer requires
exact `frames_per_scene` equality (Y11 C3 reason). New test 13: split
manifest covers every CSV scene exactly once across cal+val per group;
no frame-level entries; `input_yaml_sha256` matches current YAML.

**Y16 (Codex SEAL guard 1 — frame_interval_min first-class)**: Top-level
`frame_interval_min: 60` field; `expected_timepoints.expected_frame_interval_min`
references for consistency. "Missing frame interval raises" test
unambiguous because of top-level field.

**Y17 (Codex SEAL guard 2 — image_size_px axis order explicit)**:
`image_size_px_width_height: [1388, 1040]` with explicit name. NumPy/PIL
report shape as (height, width); validation must compare correctly,
not accidentally swap axes.

---

## 3. Sanity Gate (impl writes in module docstring or sibling sanity doc)

1. **Units**:
   - `pixel_size_um_xy`: μm/pixel (P0-1 derivation:
     `sqrt(Area_um2/Area_px) = sqrt(4.194304) = 2.048`)
   - `pixel_area_um2`: μm² (= `pixel_size_um_xy[0] · pixel_size_um_xy[1]`)
   - `frame_interval_min`: minutes (60)
   - `image_size_px_width_height`: pixels (1388 × 1040 width × height)
   - `Area_um2`, `EffectiveRadius_um`, `Perimeter_um`: μm² / μm
   - `A_over_A0`, `Solidity`, `Circularity`: dimensionless

2. **Boundary**:
   - Empty CSV (no rows): validation 7 raises "no scenes found"
   - Missing mask file: validation 9 fail-closed `MaskNotFoundError`
   - Missing overlay (optional): NOT error; diagnostic count
   - Invalid YAML structure: validation 1-4 raise specific
     `MissingFieldError`
   - Pixel size mismatch beyond tolerance: validation 10 raises with
     row index + observed vs declared values
   - Per-scene frame count > declared max: validation 7 raises
   - Non-monotonic Frame or Time_min within scene: validation 7 raises
   - Cal/val overlap: validation 5 raises with overlapping scene IDs
   - Sidecar dir missing: validation 12 raises `SidecarPathError`
   - Manifest YAML SHA mismatch: validation 13 raises
     `ManifestStaleError`

3. **Conservation**:
   - Loader is pure read-only — no PI data mutation
   - Split manifest deterministic: same YAML + seed → exact same bytes
   - CSV row count = sum over groups (validation 7)
   - Cal ∪ val = all scenes per group, cal ∩ val = ∅ (validations 5+13)

4. **Numerical**:
   - `np.float64` for pixel size, area, ratios
   - `np.isclose(rtol=1.0e-6, atol=1.0e-9)` for pixel size CSV
     consistency (Y12)
   - `int` for scene/frame/row counts
   - Hash: SHA-256 hex for `input_yaml_sha256` (Y13 tamper detection)

5. **Sign**:
   - All counts ≥ 0 (validation rejects negative)
   - Pixel size > 0 (validation 1)
   - Frame interval > 0 (validation 2)
   - Probabilities `ratio: [0.7, 0.3]` sum to 1.0 within tolerance

6. **Measurement-protocol** (Hard Rule 11, central anchor):
   - **P0-1 pixel size derivation**: declared `pixel_size_um_xy` validated
     against CSV-derived `sqrt(Area_um2/Area_px)` row-wise. Audit
     estimate (1.024 μm/px) was 4× wrong; lock requires CSV consistency
     check (validation 10) to catch any future audit drift.
   - **P0-2 CSV authoritative**: NEVER derive t-index from Frame index;
     CSV `Series` + `filename` columns are source of truth. Bare uses
     stride 3 (`t002/t005/t008/...`); Lam4 mixed; uniform stride
     assumption silently skips frames.
   - **Y2 mask path resolution**: explicit `resolve_mask_path` /
     `resolve_overlay_path`; NEVER assume t-index pattern; mask required
     (validation 9), overlay diagnostic-only.
   - **Y9 condition explicit**: `condition` field separate from
     `group_id`; future multi-date contracts cannot ambiguate condition
     by string-splitting.
   - **Y11 frame count variation**: per-scene frame counts vary (Bare
     8×83=664 expected vs actual 602); validation uses `min/max +
     monotonic` instead of exact equality.
   - **Y13 split manifest discipline**: generated output, manual edit
     FORBIDDEN; `input_yaml_sha256` enforces YAML-manifest consistency
     (HB#5/Item5 sister-pattern: no wall-clock timestamp for
     reproducibility).
   - **Y14 calibration role clarity**: `pipeline_loader_metric_dashboard_development_only`
     + `forbidden_use: physics_parameter_fitting`. Hard Rule 1
     strengthened in YAML metadata; calibration scope restricted.
   - **Y17 image axis order explicit**: NumPy `arr.shape == (height,
     width)` vs PIL/typical convention `(width, height)` — naming
     `image_size_px_width_height` prevents silent axis swap.
   - **Audit `id=2114` framework-completion discipline**: this is V2-1
     deliverable, NOT plumbing scale-up. Lock seal is "V2-1 imaging
     input contract design lock sealed (V2-1 of V2-1/V2-2 Layer 1
     framework completion path)"; NOT "Layer 1 milestone" or "Layer 1
     framework complete".

---

## 4. Test Catalog (13 tests)

### Required field presence (4)

1. `test_imaging_contract_missing_pixel_size_raises` — Y1+Y16 sister:
   YAML without `pixel_size_um_xy` raises `MissingFieldError`
2. `test_imaging_contract_missing_frame_interval_raises` — Y16: YAML
   without top-level `frame_interval_min` raises `MissingFieldError`
3. `test_imaging_contract_missing_channels_raises` — YAML without
   `channels` raises `MissingFieldError`
4. `test_imaging_contract_empty_metric_registry_or_name_raises` — empty
   `metric_registry` or any metric with empty `name` raises

### Calibration/validation split (1)

5. `test_imaging_contract_calibration_validation_scene_overlap_raises`
   — same scene appearing in both cal and val for any group raises
   `SplitOverlapError` (Y8 sister)

### Group/scene count consistency (3)

6. `test_imaging_contract_declared_group_scene_counts_match_csv_unique_series`
   — Declared `n_scenes` per group matches `len(set(csv["Series"]))`
   exactly
7. `test_imaging_contract_per_scene_monotonic_frame_time_max_count`
   (Y11 corrected): each scene has monotonic `Frame` AND `Time_min`;
   total CSV row count matches sum across groups; max frame count per
   scene ≤ `max_frames_per_scene_per_group`. NOT exact equality.
8. `test_imaging_contract_every_csv_series_directory_exists` — for each
   group, every unique `csv["Series"]` value resolves to existing
   directory under `group_root`

### Filename resolution (1)

9. `test_imaging_contract_every_csv_filename_resolves_to_mask_png`
   (Y10): for each row, `resolve_mask_path(group_root, row)` returns
   path to existing `.png` file. Overlay optional.

### Pixel size + transform consistency (2)

10. `test_imaging_contract_pixel_size_matches_sqrt_area_um2_over_area_px_within_tolerance`
    (Y1+Y12): for rows with `Area_px > 0`,
    `np.isclose(sqrt(Area_um2/Area_px), pixel_size_um_xy[0],
    rtol=1e-6, atol_um=1e-9)`. Both x and y axes (Y17 axis order
    explicit).
11. `test_imaging_contract_perimeter_um_uses_pixel_size_not_area_derived`
    (Y6 P0-1 sister): `Perimeter_um = Perimeter_px * pixel_size_um`,
    NOT derived from `EffectiveRadius_um` or `Area_um2`. Loader
    output has correct values for sample rows.

### Sidecar marker (1)

12. `test_imaging_contract_sidecar_model2_framewise_marked_reference_only_not_for_fitting`
    (Y7): YAML loads sidecars with `allowed_use:
    reduced_surrogate_reference_only` and `forbidden_use:
    dynamics_parameter_fitting`; loader exposes paths as metadata only
    (NOT loaded into metric registry).

### Split manifest discipline (1)

13. `test_imaging_contract_split_manifest_covers_every_scene_exactly_once_no_frame_entries_yaml_sha_matches`
    (Y13+Y15 NEW): generated manifest covers `cal ∪ val = all scenes`
    per group, `cal ∩ val = ∅` per group, no frame-level entries (only
    scene-level), `input_yaml_sha256` matches current YAML SHA-256.

---

## 5. Files

- `docs/v2/v2_layer_1_imaging_input_contract_locked.md` (this file)
- `acs/v2/data_contract.py` (EXISTING — reused, no changes per Y3)
- `acs/v2/imaging_contract/__init__.py` (NEW)
- `acs/v2/imaging_contract/loader.py` (NEW): `load_imaging_contract`,
  `validate_imaging_contract`, `ImagingInputContract`,
  `ImagingInputContractValidationResult`, plus dataclasses
  (`GroupSpec`, `ChannelSpec`, `FilenameResolutionSpec`, `SplitSpec`,
  `SidecarSpec`, `ProvenanceSpec`)
- `acs/v2/imaging_contract/split_builder.py` (NEW): `build_split_manifest`
  + CLI entry (`python -m acs.v2.imaging_contract.split_builder`)
- `acs/v2/imaging_contract/path_resolver.py` (NEW): `resolve_mask_path`,
  `resolve_overlay_path`, naming convention dispatchers per group
- `configs/imaging/260313.yaml` (NEW, this contract spec verbatim)
- `configs/imaging/generated/260313_split_seed20260506.json` (NEW,
  generated by split_builder, committed, manual edit FORBIDDEN)
- `tests/v2/test_v2_imaging_input_contract.py` (NEW, 13 tests per §4)
- `data/experimental/imaging/{260313_Bare,260313_Lam4,260313_Pre}/`
  (already created, gitignored, symlinks to PI Spread_data/)
- `.gitignore` (already updated, commit `2111c2a`)

---

## 6. References

- PI directives:
  - `id=2187` Option 1 pivot — vertical 5-layer framework completion
  - `id=2188` PI A안 imaging input choice
  - `id=2190` "260313만 사용해" (single date scope)
  - `id=2232` Q1=(a)/Q3=(a)/Q4=(b)
  - `id=2235` "동의" Q5 70/30 condition-stratified scene-level + seed=20260506
  - `id=1985` aggressive design debate posture
  - `id=2114` external audit (plumbing-vs-mechanism gap)
- v2 roadmap source: `docs/v2/10_dev_roadmap_v2.md` Phase V2-1
- v2 vision source: `docs/v2/00_project_vision_v2.md` ("image-constrained" identity)
- 5-layer architecture source: `docs/01_model_architecture.md`
- Sister locks (verified at source for Y3 reuse):
  - `acs/v2/data_contract.py` (`ArtifactKind`, `ImagingDatasetSpec`,
    `MetricSpec`, `V2DataContract`, validators)
- Audit-driven framing source:
  - `id=2114` plumbing-vs-mechanism gap
  - `id=2210` 9-group inspection (deprecated by `id=2231`)
  - `id=2231` 260313 single-date inspection (corrected scope)
- File ops verified at source:
  - `data/experimental/imaging/{260313_Bare,260313_Lam4,260313_Pre}/`
    symlinks created
  - `.gitignore` updated commit `2111c2a` push
- Direct CSV inspection (P0-1 + P0-2 + Y11):
  - `260313_Bare.csv` row 0: `Area_um2=162214.7072, Area_px=38675.0,
    sqrt = 2.048 μm/px`
  - Bare scene t-index: `t002, t005, t008, ...` (stride 3)
  - Lam4 scene t-index: `t001, t002, t005, t008, ...` (mixed)
  - Bare expected `8 × 83 = 664` rows vs actual `602` rows
- Hard Rule 1 (no PI data fitting): CLAUDE.md
- Hard Rule 11 (sim-experiment measurement matching): CLAUDE.md
- Magic-Number Block: CLAUDE.md (Y14 calibration role clarity
  enforcement)
- HB#5/Item 5 sweep harness sister-pattern (Y13 manifest no-timestamp
  for reproducibility):
  `docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md`,
  `docs/v2/v2_item_5_sweep_harness_locked.md`

---

## 7. Cross-room dispatch

This file is the design-team input to implementation-work for:

1. impl Claude/Codex implements per §1 module structure + §4 13-test
   catalog. **Skip Sanity Gate doc** (A-tier full debate done in
   design; impl direct + Codex single review = milestone).
2. impl-side single review pass on commit covering:
   - **Y1 P0-1**: pixel size 2.048 μm/px (NOT 1.024); validation 10
     enforces CSV cross-check
   - **Y2 P0-2**: CSV `Series` + `filename` authoritative; NEVER
     derive t-index from Frame; explicit `resolve_mask_path` /
     `resolve_overlay_path`
   - **Y3 module reuse**: `acs/v2/data_contract.py` reused (no parallel
     contract types); new `acs/v2/imaging_contract/` package only
   - **Y11 per-scene frame count variation**: validation 7 uses min/max
     + monotonic, NOT exact equality
   - **Y13 split manifest discipline**: generated by split_builder CLI,
     committed, manual edit FORBIDDEN, `input_yaml_sha256` enforced
   - **Y14 calibration role clarity**: YAML `forbidden_use:
     physics_parameter_fitting` enforced in loader docstring;
     calibration set NOT used for v2 dynamics parameter input
   - **Y16 frame_interval_min top-level**: top-level field + consistency
     with `expected_timepoints`
   - **Y17 image axis order explicit**: `image_size_px_width_height`
     name unambiguous; validation 10 + 11 use correct axis
3. On PASS: V2-1 imaging input contract design lock sealed milestone
   (NOT Layer 1 framework completion milestone — V2-2 deliverables
   pending).
4. After milestone: design-discussion can open V2-2 first deliverable
   design unit. Candidates per `docs/v2/10_dev_roadmap_v2.md`:
   - V2-2 Layer 1 single-cell core completion (boundary representation,
     projected metric trajectories, validation dashboard)
   - V2-2 Layer 2 boundary biology (lamellipodia/filopodia event
     extraction from imaging, FA state machine dynamics)

A-tier debate cycle MCP id 2237–2252 (rounds 1-3 + Codex SEAL ack
with Y16/Y17 wording guards).
