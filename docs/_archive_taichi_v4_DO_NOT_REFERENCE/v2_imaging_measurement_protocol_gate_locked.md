# V2 Imaging Measurement-Protocol Gate — Locked Design

**Date**: 2026-05-06 KST
**Authors**: Claude (impl-work) + PI framing decision (this session) +
Codex name-choice + wording review (`id=2328`).
**Source units**:
- Original blocker `id=2278` (Lam4 Pos31 24.79% mismatch + analysis
  source-pipeline question).
- PI session-active ack `id=2287`.
- PI source-pipeline directive (resolution of `id=2278`,
  `/Users/sw1/Desktop/spread_analysis/`).
- Codex guardrail review `id=2304` of `id=2302` (5 outlier-exclusion
  guardrails + tolerance stance + meta-test reject/require lists).
- PI phase-scope correction (260313 = spheroid-level) — design-discussion
  `id=2309/2310/2311` + impl-work `id=2312`.
- This session's PI framing call: RENAME approved, restructure
  deferred (option α). **PI framing correction (2026-05-06
  17:55 KST)**: v2 is a cell-resolved spheroid simulator. 260313
  is v2's primary validation anchor at integrated spheroid stages.
  Layer 1 dynamics are component building blocks of the integrated
  spheroid simulation; they are validated downstream when V2-3 /
  V2-5 stages mature, not in isolation against single-cell PI
  data. This lock unit validates the measurement-protocol layer
  (top-down spheroid projection area) of that pipeline.
- Codex name pick + wording edits `id=2326/2328`.
**Tier**: B-tier compressed lock (single-iteration; algorithm work
already complete in `f4e09aa`, this lock re-frames under
measurement-protocol-gate scope).
**PI ratify status**: full delegation per PI id=939/1008/2235/this-session.
**Supersedes**: `docs/v2/v2_layer_2_unit_1_imaging_to_single_cell_state_loader_locked.md` (RETRACTED in commit `f4e09aa`).

---

## 0. Scope

**Imaging measurement-protocol parity gate**: byte-exact / parity
reproduction, **for non-excluded rows**, of PI's saved-mask `Area_px`
extraction algorithm (`cv2.contourArea ∘ select_central_contour ∘
findContours`, `spread_infer.py:661`) on 260313 saved-mask PNGs.

The **saved-mask PNG is the provenance boundary**; raw-image
segmentation, U-Net model state, threshold choice, and manual-review
state are explicitly **outside** this lock. The saved mask captures
the operator's final decision (auto / manual polygon / re-reviewed)
already rasterised — see `spread_infer.py:550` for the manual path
that writes `cv2.fillPoly([manual_cnt], 255)` into the mask before
save.

**260313 is v2's primary validation anchor at integrated spheroid
stages** — v2 is a cell-resolved spheroid simulator, and 260313
spheroid-stage imaging is the experimental anchor against which
integrated spheroid simulation will be validated downstream
(V2-3 / V2-5 maturation). It is **not a Layer 1 parameter-fitting
or single-cell-isolation validation source**: Layer 1 single-cell
dynamics are component building blocks integrated into spheroid
validation, not validated alone against single-cell PI data.
EffectiveRadius_um median 436–490 μm, Solidity median 0.976,
Area_um2 median 6e5–7.5e5 μm² — far above any single-cell scale
(MCF7 single-cell projected area ~500–2000 μm² ≈ radius 12–25 μm),
consistent with the spheroid-stage role.

This unit validates the **measurement protocol** (top-down spheroid
projection area) for the v2 cell-resolved spheroid simulator's
downstream validation against 260313. It does not, by itself,
validate Layer 1 single-cell physics in isolation (Layer 1
dynamics are integrated into spheroid validation, not validated
alone).

### What this unit IS

- A **source-pipeline reproduction** anchored to PI's pipeline at
  `/Users/sw1/Desktop/spread_analysis/spread_infer.py`.
- A frozen tolerance contract: 1.5% relative residual on non-excluded
  rows, justified as a `cv2` cross-version drift budget (observed
  residual on cv2 4.13.0: 0.0).
- The **provenance discipline channel** that surfaces (NOT silently
  filters) the **PI-directed provenance exclusion: Lam4 Position(31)** —
  excluded at the trajectory level (82 frames), not as individual
  outlier rows.

### What this unit is NOT

- ❌ NOT a Layer 1 single-cell-isolation physics validation gate.
  260313 is spheroid-stage; this unit is the measurement-protocol
  parity gate for that imaging modality. v2 itself is a cell-
  resolved spheroid simulator that uses 260313 as its primary
  validation anchor at integrated spheroid stages; Layer 1 dynamics
  are component building blocks of that integrated pipeline,
  validated downstream when V2-3 / V2-5 stages mature, not in
  isolation against single-cell PI data.
- ❌ NOT a parameter-fitting pipeline (Hard Rule 1; PI data remains
  calibration/validation reference only, never physics input).
- ❌ NOT a tolerance derived from observed residuals — the 1.5%
  threshold is a `cv2` cross-version budget, derived from code-path
  determinism, not from "make the numbers fit".
- ❌ NOT a re-segmentation pipeline. The saved mask PNG is the
  provenance boundary; raw-image / U-Net / threshold / manual-review
  state is upstream of this lock and not under our re-implementation.
- ❌ NOT an outlier-tuning policy. Lam4 Position(31) is
  **excluded from protocol-identification set; reported separately**
  by PI directive — provenance, NOT tolerance tuning.
- ❌ NOT a SingleCellState construction unit. The earlier name
  `v2_layer_2_unit_1_imaging_to_single_cell_state_loader` is
  RETRACTED (commit `f4e09aa`); the work was reframed as this
  measurement-protocol gate after PI's phase-scope correction.

### Wording boundary

This is the **imaging measurement-protocol gate** under the V2-1
imaging input contract umbrella (`docs/v2/v2_layer_1_imaging_input_contract_locked.md`).
The original V2-1 claim `v2-layer-1-imaging-input-contract` (codex-owned)
already covers this unit; **no claim release / re-claim was needed**.
Lock seal language: "imaging measurement-protocol gate sealed". NOT
"V2-2 sealed", NOT "Layer 1 sealed", NOT "imaging→SingleCellState
loader sealed".

---

## 1. Final Lock Summary

### Hard Rule 11 reproduction claim (lock-pinned wording)

> Area_px reproduction: **4469 of 4469 non-excluded 260313 rows
> byte-exact** (rel_residual = 0.0 under cv2 4.13.0) under PI
> source-pipeline contract (`cv2.contourArea ∘ select_central_contour
> ∘ findContours`, `spread_infer.py:661`); **82 Lam4 Position(31)
> rows excluded by PI measurement-provenance directive** (entire
> trajectory; `largest_cc`-vs-`select_central_contour` mismatch under
> prior measurement protocol; Codex `id=2304` guardrails).
>
> Cross-check: per-row `Area_um2 / Area_px = 4.194304 ± float64 ULP`
> across all 4551 rows (4469 non-excluded + 82 PI-excluded).
> Validates the **measurement protocol** (top-down spheroid projection
> area) for v2's integrated spheroid-stage validation against
> 260313, not Layer 1 single-cell physics in isolation.

Total accounting:

| group | total rows | non-excluded | PI-excluded |
|---|---|---|---|
| Bare | 602 | 602 | 0 |
| Lam4 | 2104 | 2022 | 82 (Position(31)) |
| Pre  | 1845 | 1845 | 0 |
| **all** | **4551** | **4469** | **82** |

Observed `Area_um2 / Area_px` ratio per group (Hard Rule 11
substantive cross-check):

| group | rows | mean ratio | stdev | min | max |
|---|---|---|---|---|---|
| Bare | 602 | 4.194304000000 | 2.23e-16 | 4.194304000000 | 4.194304000000 |
| Lam4 | 2104 | 4.194304000000 | 1.89e-16 | 4.194304000000 | 4.194304000000 |
| Pre  | 1845 | 4.194304000000 | 8.88e-16 | 4.194304000000 | 4.194304000000 |

### Hard Rule 11 inline unit-chain (mandatory derivation)

Pixel side length: `pixel_size_um = 2.048 μm/px`
(matches V2-1 imaging contract `pixel_size_um_xy: [2.048, 2.048]`).

Unit chain:

```
Area_px         [Green's-theorem polygon area in px² units, from
                 cv2.contourArea on a contour with integer-pixel
                 vertices on the saved mask]
× pixel_size_um × pixel_size_um   [μm/px · μm/px = μm²/px²]
                                  = (2.048)² μm²/px² = 4.194304 μm²/px²
                                  = PIXEL_AREA_UM2
= Area_um2      [μm²]
```

Worked single-row example (Bare row 0,
`MD-Experiment-0004_Position(102)_t002.png`):

```
Area_px (CSV)            = 38675.0
PIXEL_AREA_UM2           = 4.194304   μm²/px²
Area_um2 (CSV)           = 162214.7072 μm²
Area_px × PIXEL_AREA_UM2 = 38675.0 × 4.194304 = 162214.7072  ✓ exact
```

The reproduction algorithm in
`acs/v2/imaging_contract/area_extractor.py` recovers `Area_px`
byte-exactly for every non-excluded row, and the unit chain above
recovers `Area_um2` byte-exactly. Hard Rule 11 measurement-protocol
consistency: SATISFIED for non-excluded rows.

### Canonical algorithm (frozen against PI source)

Source: `/Users/sw1/Desktop/spread_analysis/spread_infer.py`,
`build_row_from_saved_mask` (line 661):

1. Read saved mask PNG via `cv2.imread(path, cv2.IMREAD_GRAYSCALE)`.
2. `cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)`.
3. `select_central_contour(contours, w, h)` (mirror of
   `spread_infer.py:312`) with frozen constants:
   - `CENTER_MAX_DIST_FRAC = 0.45`
   - `BORDER_MARGIN_PX = 12`
   - `REJECT_TOUCH_BORDER = True`
   - cost: `dist + 0.0005 / (area + 1e-6)`
4. `Area_px = cv2.contourArea(cnt)` — Green's-theorem polygon area on
   integer-vertex contour, **NOT** a connected-component pixel count.
5. `Area_um2 = Area_px * PIXEL_AREA_UM2`, `PIXEL_AREA_UM2 = 4.194304`.

| name | value | PI source line |
|---|---|---|
| `CENTER_MAX_DIST_FRAC` | 0.45 | spread_infer.py:55 |
| `BORDER_MARGIN_PX` | 12 | spread_infer.py:56 |
| `REJECT_TOUCH_BORDER` | True | spread_infer.py:57 |
| select-central-contour cost weight | 0.0005 | spread_infer.py:340 |
| `PIXEL_AREA_UM2` | 4.194304 | spread_infer.py:44 |

Changing any of these breaks byte-exact reproduction of PI's CSV
`Area_px`.

### Provenance discipline (Codex guardrail `id=2304`)

**PI-directed provenance exclusion: Lam4 Position(31)**.
The original 24.79% mismatch flagged in design-discussion `id=2278`
was specifically a `largest_cc area` vs `Area_px` measurement-protocol
mismatch (multiple-component mask + `largest_cc` picking a different
component than `select_central_contour`). Under the canonical
algorithm in this lock, the 82 frames of Position(31) ALSO reproduce
at 0.0 relative residual — but the trajectory-level exclusion is
preserved per PI directive.

Specifically:

1. The exclusion is keyed by `(group, series) = ("Lam4", "Position(31)")`
   in `acs/v2/imaging_contract/area_extractor.py:PI_DIRECTED_PROVENANCE_EXCLUSIONS`.
2. The exclusion is **trajectory-level** (all 82 frames), not a
   single-row outlier policy.
3. The extractor module does **not** silently drop Pos31; it
   surfaces all 82 rows with `provenance_status =
   "excluded_by_pi_measurement_provenance"`.
4. Pos31 MUST NOT contribute to thresholds, calibration, rates, or
   acceptance windows.
5. Acceptance metric `aggregate_normal_residuals` reports Pos31 in a
   separate `excluded_*` aggregate so it cannot leak into the
   acceptance decision.
6. Final tolerance claim is scoped to the 4469 non-excluded rows.

Pos31 is
**excluded from protocol-identification set; reported separately**
— provenance, NOT tolerance tuning.

### Tolerance contract (source-first, NOT residual-first)

`ACCEPTANCE_REL_RESIDUAL_THRESHOLD = 0.015` (1.5% relative residual).
Justification chain — derived from code path, NOT from observed
residuals:

1. The algorithm is byte-deterministic in float64 arithmetic given
   the same OpenCV major version.
2. `cv2.findContours` with `CHAIN_APPROX_SIMPLE` is stable across
   patch versions of OpenCV 4.x.
3. `cv2.contourArea` uses Green's theorem on integer-vertex polygons;
   it is exact float64 arithmetic.
4. The 1.5% slack is reserved for **cv2 cross-version contour
   discretisation drift**, not for fitting.
5. Observed residual on cv2 4.13.0 (validated 2026-05-06): 0.0 on
   every non-excluded row of every group.

The 1.5% threshold matches the design-discussion `id=2278`
preliminary number coincidentally; the principled derivation is
"cv2 cross-version budget", not "match what was discussed".

### Halt-and-escalate triggers (Codex guardrail, evaluated)

| trigger | status | resolution |
|---|---|---|
| Algorithm depends on calibration params not in source dir | NOT TRIGGERED | All constants in `spread_infer.py` header. |
| After Pos31 exclusion, residuals > 1.5% | NOT TRIGGERED | Observed 0.0 on cv2 4.13.0. |
| Multiple `label_maker*`/`spread_infer*` disagree | NOT TRIGGERED | label_maker is for training-data labelling only (NOT Area_px); spread_infer / _remote / _remote2 / _fp16 share `build_row_from_saved_mask`. spread_infer.py is the canonical entry. |
| Manual review state not in saved files (Codex addition) | RESOLVED | All 260313 CSV rows have `review_mode=reconstructed_from_saved_mask`; manual operator state (auto-detected contour OR manual polygon OR re-reviewed contour) is rasterised into the saved mask PNG before save (`spread_infer.py:550`). The saved mask PNG IS the operator's final decision. |
| Lock spec name not approved by codex | RESOLVED | Name (A) `docs/v2/v2_imaging_measurement_protocol_gate_locked.md` selected by Codex `id=2328`. |
| Sanity gate proposes threshold beyond cv2 deterministic reproduction | NOT TRIGGERED | Threshold = 1.5% cv2 cross-version budget; observed 0.0. |
| Lock wording mis-frames v2 as single-cell-isolated or treats 260313 as unused/downstream-only | NOT TRIGGERED | Scope §0 frames 260313 as v2's primary integrated spheroid-stage validation anchor and disclaims only Layer-1-in-isolation validation. |

---

## 2. Test catalog (lock minimum)

`tests/v2/test_v2_imaging_measurement_protocol_gate.py` covers:

1. Algorithm constants frozen against PI source.
2. `select_central_contour` synthetic-blob behaviour (centred over
   off-centre, border-rejection check).
3. `extract_area_px_from_mask` byte-exact on synthetic disc, status
   routing for missing/unreadable mask.
4. PI-directed provenance exclusion is keyed at `(Lam4, Position(31))`
   and surfaced via `classify_provenance`.
5. `aggregate_normal_residuals` does not let Pos31 raise the
   acceptance threshold (Codex guardrail `id=2304` discipline test).
6. End-to-end reproduction on the real 260313 dataset (skipped if
   CSVs absent) — invokes
   `scripts/v2/reproduce_pi_area_extraction.py` and asserts
   `[CONTRACT PASS]`.
7. **Hard Rule 11 wording-boundary meta-test** — lock spec must
   contain the required source-pipeline-reproduction phrases AND
   must NOT contain any rejected fitting-flavoured phrases.
8. `area_extractor.py` cites `spread_infer.py` and
   `build_row_from_saved_mask`; carries the inline unit-chain phrase.
9. Sister-doc / sister-code mirror — Pos31 clause appears in
   `area_extractor.py`, `scripts/v2/reproduce_pi_area_extraction.py`,
   and this lock spec.
10. `PIXEL_AREA_UM2` consistent with V2-1 imaging YAML
    `pixel_area_um2: 4.194304`.
11. `ACCEPTANCE_REL_RESIDUAL_THRESHOLD = 0.015` is lock-pinned.

---

## 3. References

- PI source pipeline: `/Users/sw1/Desktop/spread_analysis/`
  - `spread_infer.py` (canonical entry)
  - `spread_infer.py:312` — `select_central_contour`
  - `spread_infer.py:661` — `build_row_from_saved_mask`
- V2-1 imaging input contract:
  `docs/v2/v2_layer_1_imaging_input_contract_locked.md` +
  `configs/imaging/260313.yaml`.
- Reproduction script: `scripts/v2/reproduce_pi_area_extraction.py`.
- Reproduction artefact:
  `runs/v2_imaging_measurement_protocol_gate/pi_area_reproduction.csv`
  (gitignored).
- Implementation: `acs/v2/imaging_contract/area_extractor.py`.
- Tests: `tests/v2/test_v2_imaging_measurement_protocol_gate.py`.
- Paired sanity gate:
  `docs/v2/v2_imaging_measurement_protocol_gate_sanity_gate.md`.
- Superseded retraction:
  `docs/v2/v2_layer_2_unit_1_imaging_to_single_cell_state_loader_locked.md`
  (RETRACTED in commit `f4e09aa`; kept on disk as audit trail).
- MCP IDs: 2278 / 2287 / 2304 / 2308 / 2309 / 2310 / 2311 / 2312 /
  2315 / 2326 / 2328.
- Commits: `85df64a` (violating original seal),
  `f4e09aa` (retraction), this commit (renamed seal).
