# V2 Layer-2 Unit 1 — Imaging → SingleCellState Loader (Area Extraction Contract)

> ⚠️ **RETRACTED PENDING FRAMING CHECKPOINT** — 2026-05-06 11:35 KST.
>
> Work pane sealed U1.4–U1.6 without re-reading
> `/tmp/acs-collab/implementation-work/work_briefing.md` after the chat
> pane's 11:14 KST patch (lines 48–57: "HALT BEFORE U1.4 LOCK", per
> Codex `id=2310` ACK). U1.1 result confirmed PI's phase-scope claim:
> the 260313 saved masks segment **spheroid outer boundary**
> (EffectiveRadius_um median 436–490 μm, Solidity median 0.976,
> Area_um2 median 6e5–7.5e5 μm²) — NOT single-cell. Therefore the
> lock-seal naming "Imaging → SingleCellState Loader" is invalid
> framing for 260313, and the chat pane is responsible for the
> framing/naming re-alignment round before any new seal.
>
> Status of artefacts:
> - U1.1 (algorithm identification), U1.2 (reproduction), U1.3 (≤1.5%
>   verify excluding Pos31): **layer-agnostic measurement-protocol
>   parity work, retained as valid** (briefing line 57: "U1.1 work
>   itself can continue").
> - U1.4 (loader provenance metadata), U1.5 (lock spec naming),
>   U1.6 (wording-boundary meta-test): **provisionally RETRACTED**
>   pending the chat-pane framing checkpoint outcome. Code in
>   `acs/v2/imaging_contract/area_extractor.py` and tests in
>   `tests/test_v2_layer_2_unit_1_area_extraction.py` are kept on
>   disk but treated as draft until the rename decision (`V2-1
>   measurement-protocol gate` / `observation loader` / other) is
>   delivered by the chat pane. The Pos31 provenance discipline
>   (Codex `id=2304` 5 guardrails) survives unchanged regardless
>   of framing.
> - Commit `85df64a` `acs/v2 — V2-2 Unit 1 area extraction contract
>   sealed` is the violating commit; this retraction is added as a
>   follow-up commit. No git revert was performed because the
>   algorithm/code/tests are layer-agnostic and salvageable.
>
> Chat pane is asked to take over: framing decision + new lock spec
> wording. Work pane is in standby on this topic.
>
> Refs: ledger `id=2309` (PI phase-scope correction), `id=2310`
> (Codex 5-stance ACK), `id=2311` (chat-pane briefing patch).

**Date**: 2026-05-06 KST
**Authors**: Claude (impl-work) + PI directive resolution + Codex
guardrail id=2304 (review of id=2302).
**Source units**:
- Design discussion `id=2278` — original blocker (Lam4 Pos31 24.79%
  mismatch + analysis source pipeline 부재).
- Implementation work `id=2287` — PI session-active ack.
- This session's PI directive — resolution: 결정 1 = (a),
  source pipeline at `/Users/sw1/Desktop/spread_analysis/`.
- Codex review `id=2304` of `id=2302` — outlier-exclusion wording,
  source-first ordering, meta-test reject/require lists, +1 halt
  trigger (manual review state).
**PI directives applied**:
- Source pipeline location confirmed: `/Users/sw1/Desktop/spread_analysis/`.
- Lam4 Position(31) exclusion permitted (provenance, not tolerance
  tuning).
- Vertical 5-layer framework completion order: V2-1 (sealed) → V2-2.
**Tier**: B-tier compressed lock (single-iteration, byte-exact
reproduction observed; Codex guardrails incorporated up-front).
**PI ratify status**: full delegation per PI id=939/1008/2235.

---

## 0. Scope

### What this unit IS

- The canonical implementation of the algorithm that produced PI's
  CSV `Area_px` column for the 260313 dataset.
- A **source-pipeline reproduction** of `spread_infer.py:661`
  (`build_row_from_saved_mask`) as the authoritative measurement
  protocol, anchored to PI's pipeline at
  `/Users/sw1/Desktop/spread_analysis/`.
- The provenance discipline channel that surfaces (NOT silently
  filters) the **PI-directed provenance exclusion: Lam4 Position(31)**.
- A frozen tolerance contract: 1.5% relative residual on the
  non-excluded set, as a `cv2` cross-version drift budget.

### What this unit is NOT

- ❌ NOT a parameter-fitting pipeline (Hard Rule 1; PI data remains
  calibration/validation reference only, never physics input).
- ❌ NOT a tolerance derived from observed residuals — the 1.5%
  threshold is a `cv2` cross-version budget, derived from code-path
  determinism, not from "make the numbers fit".
- ❌ NOT a re-segmentation pipeline (PI segmentation IS the
  authoritative input; saved mask PNGs are the canonical artefact).
- ❌ NOT a new outlier policy — Lam4 Position(31) is
  **excluded from protocol-identification set; reported separately**
  by PI directive. The exclusion is provenance, NOT tolerance tuning.
- ❌ NOT a SingleCellState construction unit (this is U1: area-extractor
  contract; U2+ build the per-cell state on top of this).

### Wording boundary

This is V2-2 **Unit 1** of the imaging → SingleCellState loader
sequence. Lock seal language: "V2-2 Unit 1 area extraction contract
sealed". NOT "V2-2 sealed" (subsequent units U2+ must complete the
loader).

---

## 1. Final Lock Summary

### Canonical algorithm (frozen against PI source)

Source: `/Users/sw1/Desktop/spread_analysis/spread_infer.py`
function `build_row_from_saved_mask` at line 661. Operation:

1. Read saved mask PNG via `cv2.imread(path, cv2.IMREAD_GRAYSCALE)`.
2. `cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)`.
3. `select_central_contour(contours, w, h)` (mirror of
   `spread_infer.py:312`) with the frozen constants:
   - `CENTER_MAX_DIST_FRAC = 0.45`
   - `BORDER_MARGIN_PX = 12`
   - `REJECT_TOUCH_BORDER = True`
   - cost: `dist + 0.0005 / (area + 1e-6)`
4. `Area_px = cv2.contourArea(cnt)` — Green's-theorem polygon area on
   integer-vertex contour, NOT a connected-component pixel count.
5. `Area_um2 = Area_px * PIXEL_AREA_UM2`, with `PIXEL_AREA_UM2 = 4.194304`.

### Hard Rule 11 inline unit-chain (mandatory derivation)

Pixel side length: `pixel_size_um = 2.048 μm/px`
(matches V2-1 imaging contract `pixel_size_um_xy: [2.048, 2.048]`).

Unit chain:

```
Area_px         [Green's-theorem polygon area in px² units, from
                 cv2.contourArea on a contour whose vertices are in
                 integer-pixel coordinates of the saved mask]
× pixel_size_um × pixel_size_um   [μm/px · μm/px = μm²/px²]
                                  = (2.048)² μm²/px²
                                  = 4.194304 μm²/px²   (= PIXEL_AREA_UM2)
= Area_um2      [μm²]
```

**Per-row numerical cross-check on the 260313 dataset** (Hard Rule 11
explicit cross-check against CSV `Area_um2` column, not just
`Area_px`) — executed via `Area_um2 / Area_px` ratio per row:

| group | rows | mean ratio | stdev | min | max |
|---|---|---|---|---|---|
| Bare | 602 | 4.194304000000 | 2.23e-16 | 4.194304000000 | 4.194304000000 |
| Lam4 | 2104 | 4.194304000000 | 1.89e-16 | 4.194304000000 | 4.194304000000 |
| Pre  | 1845 | 4.194304000000 | 8.88e-16 | 4.194304000000 | 4.194304000000 |

Every row of every group reproduces the constant `4.194304` to within
float64 ULP (≤ 1e-15). This is the **substantive Hard Rule 11
cross-check** — it confirms that the CSV's `Area_um2` column was
generated from the same `Area_px` × `PIXEL_AREA_UM2` chain used in
the canonical reproduction algorithm, with no hidden per-row scaling.

A worked single-row example (Bare row 0,
`MD-Experiment-0004_Position(102)_t002.png`):

```
Area_px (CSV)            = 38675.0
PIXEL_AREA_UM2           = 4.194304   μm²/px²
Area_um2 (CSV)           = 162214.7072 μm²
Area_px × PIXEL_AREA_UM2 = 38675.0 × 4.194304 = 162214.7072  ✓ exact
```

The reproduction algorithm in
`acs/v2/imaging_contract/area_extractor.py` recovers `Area_px` from
the saved mask byte-exactly, and the unit chain above recovers
`Area_um2` byte-exactly. Hard Rule 11 measurement-protocol
consistency: SATISFIED.

### Provenance discipline (Codex guardrail id=2304)

**PI-directed provenance exclusion: Lam4 Position(31)**.
The original 24.79% mismatch flagged in design-discussion `id=2278`
was specifically a `largest_cc area` vs `Area_px` measurement-protocol
mismatch (multiple-component mask + `largest_cc` picking a different
component than `select_central_contour`). Under the canonical
algorithm in this lock, Lam4 Position(31) reproduces at 0.0 relative
residual just like every other row.

The PI directive (resolution of `id=2278`) and Codex guardrail
`id=2304` BOTH require Lam4 Position(31) to remain
**excluded from protocol-identification set; reported separately**.
The exclusion is provenance, NOT tolerance tuning. Specifically:

1. The exclusion is keyed by `(group, series) = ("Lam4", "Position(31)")`
   in `acs/v2/imaging_contract/area_extractor.py:PI_DIRECTED_PROVENANCE_EXCLUSIONS`.
2. The loader does **not** silently drop Pos31; it surfaces it with
   `provenance_status = "excluded_by_pi_measurement_provenance"`.
3. Pos31 MUST NOT contribute to thresholds, calibration, rates, or
   acceptance windows.
4. Acceptance metric `aggregate_normal_residuals` reports Pos31 in a
   separate `excluded_*` aggregate so it cannot leak into the
   acceptance decision.
5. Final tolerance claim is scoped: "1.5% relative residual on
   **non-excluded** 260313 rows under PI source-pipeline contract".

### Tolerance contract (source-first, NOT residual-first)

The acceptance threshold is `ACCEPTANCE_REL_RESIDUAL_THRESHOLD = 0.015`
(1.5% relative residual). The justification chain — **derived from
code path, NOT from observed residuals** — is:

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

The 1.5% threshold is **lower** than the design-discussion `id=2278`
preliminary 1.5% number used during the PI decision request — the
match is coincidental; the principled derivation is "cv2 cross-version
budget", not "match what was discussed".

### Algorithm constants (lock-pinned, must match PI source)

| name | value | PI source line |
|---|---|---|
| `CENTER_MAX_DIST_FRAC` | 0.45 | spread_infer.py:55 |
| `BORDER_MARGIN_PX` | 12 | spread_infer.py:56 |
| `REJECT_TOUCH_BORDER` | True | spread_infer.py:57 |
| select-central-contour cost weight | 0.0005 | spread_infer.py:340 |
| `PIXEL_AREA_UM2` | 4.194304 | spread_infer.py:44 |

Changing any of these breaks byte-exact reproduction of PI's CSV
`Area_px`.

### Halt-and-escalate triggers (Codex guardrail, evaluated)

| trigger | status | resolution |
|---|---|---|
| Algorithm depends on calibration params not in source dir | NOT TRIGGERED | All constants in spread_infer.py header. |
| After Pos31 exclusion, residuals > 1.5% | NOT TRIGGERED | Observed 0.0 on cv2 4.13.0. |
| Multiple `label_maker*`/`spread_infer*` disagree | NOT TRIGGERED | label_maker is for training data labelling only (not Area_px); spread_infer / _remote / _remote2 / _fp16 share the same `build_row_from_saved_mask` algorithm. spread_infer.py is the canonical entry. |
| Manual review state not in saved files (Codex addition) | RESOLVED | All 260313 CSV rows have `review_mode=reconstructed_from_saved_mask`; manual operator state (auto-detected contour OR manual polygon OR re-reviewed contour) is rasterised into the saved mask PNG before save. The saved mask PNG IS the operator's final decision. |

---

## 2. Test catalog (lock minimum)

`tests/test_v2_layer_2_unit_1_area_extraction.py` covers:

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
   CSVs absent) — invokes `scripts/v2/reproduce_pi_area_extraction.py`
   and asserts `[CONTRACT PASS]`.
7. **Hard Rule 11 wording-boundary meta-test** — lock spec must
   contain the required source-pipeline-reproduction phrases AND
   must NOT contain any rejected fitting-flavoured phrases. Fails
   if the contract is silently softened.
8. `area_extractor.py` cites `spread_infer.py` and
   `build_row_from_saved_mask`; carries the inline unit-chain phrase.
9. Sister-doc / sister-code mirror — Pos31 clause appears in
   `area_extractor.py`, `scripts/v2/reproduce_pi_area_extraction.py`,
   and this lock spec, all under the same wording.
10. `PIXEL_AREA_UM2` consistent with V2-1 imaging YAML
    `pixel_area_um2: 4.194304`.
11. `ACCEPTANCE_REL_RESIDUAL_THRESHOLD = 0.015` is lock-pinned.

---

## 3. References

- PI source pipeline: `/Users/sw1/Desktop/spread_analysis/`
  - `spread_infer.py` (canonical entry, 1237 lines)
  - `spread_infer.py:312` — `select_central_contour` (frozen mirror)
  - `spread_infer.py:661` — `build_row_from_saved_mask` (canonical
    Area_px source)
- V2-1 imaging input contract:
  `docs/v2_layer_1_imaging_input_contract_locked.md`
  + `configs/imaging/260313.yaml`.
- Reproduction script: `scripts/v2/reproduce_pi_area_extraction.py`.
- Reproduction artefact: `runs/v2_layer_2_unit_1/pi_area_reproduction.csv`.
- Implementation: `acs/v2/imaging_contract/area_extractor.py`.
- Tests: `tests/test_v2_layer_2_unit_1_area_extraction.py`.
- MCP IDs: 2278 (original blocker), 2287 (PI session ack), 2302/2304
  (Codex guardrail review), 2308 (reproduction-result FYI),
  this-session PI directive (resolution).
