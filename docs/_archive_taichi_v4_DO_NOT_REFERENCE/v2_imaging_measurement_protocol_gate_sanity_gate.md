# V2 Imaging Measurement-Protocol Gate — Sanity Gate

**Date**: 2026-05-06 KST
**Paired lock**: `docs/v2/v2_imaging_measurement_protocol_gate_locked.md`
**Protocol**: CLAUDE.md "Sanity Gate Protocol" (mandatory before
first execution; the V2-1 imaging input contract was the prior
sealed sister gate).

This is a **measurement-protocol parity gate**, not a physics-numerics
gate. The Sanity Gate items below are interpreted in that context:
"first execution" means the first time the reproduction algorithm is
run against the 4551 260313 rows; the algorithm itself contains no
new physics or new discretisation, only a `cv2`-deterministic
reproduction of PI's source pipeline.

---

## What this gate measures

- **Primary**: byte-exact / parity reproduction, **for non-excluded
  rows**, of PI's saved-mask `Area_px` extraction algorithm
  (`cv2.contourArea ∘ select_central_contour ∘ findContours`,
  `spread_infer.py:661`) on 260313 saved-mask PNGs.
- **Cross-check**: per-row `Area_um2 / Area_px = 4.194304 ± float64
  ULP` across all 4551 rows (Hard Rule 11 substantive cross-check).
- **Provenance discipline**: 82 Lam4 Position(31) rows surfaced with
  `provenance_status = "excluded_by_pi_measurement_provenance"`
  (entire trajectory; PI directive); they do not contribute to
  thresholds or acceptance.

Total accounting: 4551 rows = 4469 non-excluded + 82 PI-excluded.

---

## 1. Dimensional analysis

- Inputs: saved-mask PNG (`uint8`, shape `[H, W]`).
- Algorithm: `cv2.findContours` returns integer-pixel-vertex contours.
  `cv2.contourArea` applies Green's theorem to those contours,
  yielding an area in **px² units** (i.e., the area of the integer-
  vertex polygon in pixel coordinates).
- Conversion: `Area_um2 = Area_px × PIXEL_AREA_UM2`,
  `PIXEL_AREA_UM2 = 4.194304 μm²/px²` (= `pixel_size_um² = (2.048 μm/px)²`).
- Non-dimensional numbers: not applicable (this is a measurement-
  reproduction protocol, no transport/wave/elasticity equations).
- CFL / time-step: not applicable (no time integration).
- **Dimensional check: PASS** — Area_px [px²] × [μm²/px²] = Area_um2
  [μm²].

## 2. Boundary cases

- **Empty mask** (all zeros): `findContours` returns no contours;
  `select_central_contour` returns `None`;
  `extract_area_px_from_mask` returns
  `(None, PROVENANCE_STATUS_NO_CONTOUR)`. No NaN / divide-by-zero.
- **Mask file missing**: returns
  `(None, PROVENANCE_STATUS_MASK_MISSING)` without raising.
- **Mask file unreadable** (corrupt PNG): returns
  `(None, PROVENANCE_STATUS_MASK_READ_FAIL)`.
- **Single full-image foreground** (touching all four borders): the
  primary path rejects the contour via `REJECT_TOUCH_BORDER`;
  fallback path (mirror of `spread_infer.py:345-357`) picks the
  closest-to-center contour without border check, matching PI source
  behaviour.
- **Multiple components**: `select_central_contour` minimises
  `dist + 0.0005/(area + 1e-6)` — the same cost PI uses (line 340),
  preventing border-touching or off-centre artefacts from being
  picked when a centred valid contour exists. The Pos31 24.79%
  prior mismatch was caused by `largest_cc` (NOT this algorithm)
  picking the wrong component on multi-component frames.
- **Zero CSV `Area_px`**: handled in
  `reproduce_one_row` — relative residual reported as 0.0 if both
  reproduce to 0, NaN if reproduction differs.

## 3. Conservation invariants

This gate has no conservation laws of its own — it is a reproduction
contract, not a dynamic system. The relevant invariant is
**measurement-protocol consistency**: per-row
`Area_um2 / Area_px = 4.194304` is constant to within float64 ULP
across all 4551 rows. **Verified empirically** (see lock spec §1
ratio table).

## 4. Numerical sanity

- Float precision: `cv2.contourArea` returns `float` (float64);
  CSV stores float values; reading via pandas keeps float64; the
  reproduction stays in float64. No f32 hot path.
- `cv2` version drift: lock spec pins observed residual = 0.0 under
  cv2 4.13.0; the 1.5% acceptance threshold is reserved as a
  **cross-version contour-discretisation budget**, NOT a fitting
  knob. If a future cv2 minor release perturbs `findContours` /
  `contourArea` discretisation by > 1.5% on any non-excluded row,
  the gate fails — the response is to investigate the cv2 change
  and, if the change is a real bugfix, update the lock pin (NOT to
  raise the threshold).
- `cv2` deterministic across runs: yes — `findContours` and
  `contourArea` are deterministic on a given input mask.
- Mask resolution: saved masks are at the same resolution PI's
  pipeline writes them; the V2-1 imaging contract YAML records
  `image_size_px_width_height: [1388, 1040]` and
  `pixel_area_um2: 4.194304`. The reproduction makes no assumption
  about resolution — it reads `mask.shape` and uses that for the
  `select_central_contour` width/height arguments.

## 5. Sign / sense check

- `cv2.contourArea`: returns non-negative area for a non-self-
  intersecting closed contour. Sign is meaningful only for
  open-contour / orientation tests, which do not apply here.
- `select_central_contour` cost: `dist + 0.0005/(area + 1e-6)` —
  smaller is better (more central, larger). Sign agrees with the
  intent ("most centred non-border contour").
- `provenance_status` is a discrete-valued discipline channel; no
  sign semantics.

## 6. Measurement-protocol consistency (Sanity Gate item)

The reproduction algorithm exactly matches the PI-source measurement
protocol on the saved mask. The substantive cross-check is per-row
`Area_um2 / Area_px = 4.194304` constancy across all 4551 rows
(verified to float64 ULP). The cross-check is performed against the
CSV's own `Area_um2` column, not just `Area_px` — Hard Rule 11
satisfied.

The "off-proof" measurement protocol items (CLAUDE.md Sanity Gate
6) have been audited:

- **Sample domain**: every row of every CSV in `data/experimental/
  260313_*.csv`, with mask path resolved via the V2-1 imaging
  contract `filename_resolution: csv_authoritative`. No band
  averaging, no integration, no off-proof sampling — the gate
  evaluates the algorithm at the same integer-pixel domain that the
  CSV was generated on.
- **Off-proof contributions**: none — the algorithm is the same
  function on the same input; no auxiliary measurement has been
  introduced.

---

## Failure conditions

The gate fails if any of the following hold:

1. Any non-excluded row has `rel_residual > 0.015` after running
   `scripts/v2/reproduce_pi_area_extraction.py`.
2. Any non-excluded row has `rel_residual > 0.0` AND no `cv2`
   version-drift root cause has been documented (i.e., a new mode
   of failure that is not "cv2 patch perturbation").
3. The Hard Rule 11 cross-check ratio `Area_um2 / Area_px`
   deviates from `4.194304` beyond float64 ULP (~1e-15) on any row.
4. The reproduction script exits non-zero (mask read failure,
   missing file, etc.).
5. Pos31 rows are silently dropped by the extractor module
   (provenance-discipline violation; Codex `id=2304` guardrail).
6. Pos31 rows raise the acceptance threshold (Pos31 leaks into
   `aggregate_normal_residuals.normal_max_rel_residual`).
7. Lock spec wording fails the wording-boundary meta-test (reject
   phrases present, or required phrases missing).

## Pre-execution status

- Algorithm: complete (`acs/v2/imaging_contract/area_extractor.py`,
  layer-agnostic; surviving from the retracted commit `85df64a`
  through `f4e09aa` retraction; algorithm itself is valid).
- Reproduction artefact:
  `runs/v2_imaging_measurement_protocol_gate/pi_area_reproduction.csv`
  (gitignored;
  4551 rows produced, all `rel_residual = 0.0`, all
  `provenance_status` either `normal` or
  `excluded_by_pi_measurement_provenance`).
- Test catalog: `tests/v2/test_v2_imaging_measurement_protocol_gate.py`
  (16/17 measurement-protocol tests PASS, 1 opt-in full-reproduction
  subprocess test SKIPPED unless `RUN_FULL_REPRODUCTION_GATE=1`).
- V2-1 regression: 21/21 PASS (no break from the algorithm work).
- Combined pre-seal command:
  `.venv-collab/bin/python -m pytest tests/v2/test_v2_imaging_input_contract.py tests/v2/test_v2_imaging_measurement_protocol_gate.py`
  reported 37 passed, 1 skipped before this path cleanup.

---

## References

- Paired lock: `docs/v2/v2_imaging_measurement_protocol_gate_locked.md`.
- PI source: `/Users/sw1/Desktop/spread_analysis/spread_infer.py`
  (`build_row_from_saved_mask` line 661, `select_central_contour`
  line 312, `compute_shape_metrics` line 634).
- V2-1 imaging input contract:
  `docs/v2/v2_layer_1_imaging_input_contract_locked.md`.
- Implementation: `acs/v2/imaging_contract/area_extractor.py`.
- Reproduction script: `scripts/v2/reproduce_pi_area_extraction.py`.
- MCP IDs: 2278 / 2287 / 2304 / 2308 / 2309 / 2310 / 2311 / 2312 /
  2315 / 2326 / 2328.
