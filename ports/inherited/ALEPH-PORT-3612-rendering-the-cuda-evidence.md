# ALEPH-PORT-3612 — the per-channel parity matrix: rendering what a single number per case erases

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3612` |
| Lane | `46143f30` Track V2 (visualisation) |
| Status | `PROPOSED` |
| Written | `2026-08-01` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `CONCEPT_ONLY` |
| Exists because | Lane G7 drove 49 law cases on an RTX A5000 and wrote 228 KB of JSON. The evidence is real and **nobody can read it without opening JSON.** That is the convenience half of the PI's standing instruction failing. Worse, the natural summary — one pass/fail per case — *erases the port's central finding*, because `WRONG_ENERGY_HALF` and `WRONG_FROZEN_NORMAL` are caught in one channel while remaining **bit-identical to the true kernel in another**. |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
from aleph.viz.evidence import (
    DECADE_HIGH,
    DECADE_LOW,
    EXACT_LANE,
    CHANNEL_STYLE,
    THIN_HEADROOM,
    ChannelState,
    ChannelStyle,
    CaseRow,
    DescentWindow,
    EvidenceBundle,
    EvidenceError,
    FieldRow,
    ParityEvidence,
    ScatterEvidence,
    discover_artefacts,
    load_bundle,
    log_position,
    read_descent_window,
    read_parity,
    read_scatter,
)

from aleph.viz.evidence_figure import (
    EXIT_OK,
    EXIT_REFUSED,
    EXIT_FAULT,
    evidence_figure,
    main,
    render_evidence_svg,
)
```

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none.** `/Users/sw1/ffn_cellsim` was not read for this entry. |
| Source commit | n/a |
| Source path | n/a |
| Source symbol(s) | n/a |
| Read from | Aleph's own run artefacts, written by `scripts/unattended_runner.py` (`ALEPH-PORT-3607`) on host `GBook`, fetched read-only over ssh. |
| Working tree == commit? | n/a |

The inputs are `~/.aleph_runner/artefacts/{j1_law_parity,j2_scatter_determinism,g6_descent_window}.json`,
schema `aleph.unattended_runner.summary/1` and the per-job shapes `ALEPH-PORT-3607` §6 defines. This
entry reads them and **writes none**; it is a viewer, and `aleph/viz/**` may not import
`validation/**` or run a device.

## 3. Why source-derived porting beats clean-room

It does not, and nothing is ported. There is no `ffn_cellsim` symbol behind this: the reference
visualisation has no law-parity harness, no per-field ULP, and no wrong-kernel enumeration to render,
so there is nothing there to derive from. This is written clean-room against artefacts Aleph produced,
and the entry exists because `ports/TEMPLATE.md` §12 and the acceptance discipline are worth applying
to a renderer whose failure mode is *silent* — see §9.

## 4. Physical or mathematical law represented

**None.** What crosses the boundary is a **reading convention**, and saying so is the honest answer.
Three parts of it are load-bearing enough to derive rather than assume:

**(a) A channel is the unit of judgement, not a case.** A law case yields `k` fields (an energy and
one or two force arrays). `ALEPH-PORT-3605` M9 and `-3606` §13a M8 both record gates that were silent;
this entry records the dual — a gate that is *loud on one field and silent on another within the same
case*. Define, for a wrong variant `w` of family `F` with true variant `T`, and field `f`:

    blind(w, f)  <=>  ulp(w, f) == ulp(T, f)  exactly

Equality of the reported ULP is not a proxy for bit-identity — it *is* it at the resolution the
artefact carries, because ULP is a deterministic function of the two arrays and the reference scale,
and the reference scale is a property of `F` and the state, not of the variant. A variant that moved
any element of `f` would move `max_abs` and hence `ulp`. §14.2 states what this does not establish.

**(b) ULP is logarithmic.** The measured range is `[0.0159, 2.5166e7]` — **nine decades**. A linear
axis maps every real result onto the first pixel and the whole figure becomes a picture of
`WRONG_ENERGY_HALF`. The axis is therefore `log10`, with one exception that is not a detail:

    ulp == 0.0  is an EXACT result, not a small one.

`log10(0)` is `-inf`, and clamping zero into the bottom decade would draw an *exact* agreement — the
thing `ALEPH-PORT-3604` shipped as its headline — as merely a very good one. Exact zero gets its own
lane, drawn left of the axis with a visible break, and `log_position` refuses to place it on the
continuum.

**(c) Margin is a ratio, and there are two of them.**

    pressure(f) = ulp(f) / effective_ulp_budget(f)          how close this result is to its gate
    headroom(f) = roundoff_ulp_bound(f) / effective_ulp_budget(f)   how far the gate is under round-off

A budget with `headroom < 2` is *thin*: it sits less than a binary order under the round-off bound its
own ledger derived, so a modest change of mesh or state moves the gate through the data. `THIN_HEADROOM
= 2.0`, and the figure marks those gates rather than leaving the reader to divide two columns.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `ulp` | float32 ULP of the field's reference scale | dimensionless | `[0, inf)`; `0.0` is exact |
| `effective_ulp_budget` | float32 ULP | dimensionless | `(0, inf)` |
| `roundoff_ulp_bound` | float32 ULP | dimensionless | `(0, inf)` |
| `pressure`, `headroom` | ratio | dimensionless | `[0, inf)` |
| `slack_multiple_of_u32` | multiples of `u32 = 2^-24` | dimensionless | `{0, 1, 2, 4, 8}` as measured |

Singular and boundary cases, each with the behaviour required:

- **`ulp == 0.0`** — placed in `EXACT_LANE`, never on the log continuum. `log_position` **raises**
  `EvidenceError` if asked to place it, so a caller cannot get a plausible number by accident.
- **`ulp` below the bottom decade** (`10**DECADE_LOW`) — clamped to the axis start **and flagged**, so
  a clamped mark is never read as a measured one.
- **a field present for one variant of a family and absent for another** — emitted as a `FieldRow`
  with `state == ABSENT` and a reason. It is **not** skipped: `render.py` established that an owner
  that cannot be drawn is listed, and a field that is not in the artefact gets the same treatment.
- **a family with no `TRUE` variant** — every `blind` verdict for that family is `None` (unknown), not
  `False`. There is no baseline to be bit-identical to, and reporting "not blind" would be a claim.
- **`effective_ulp_budget == 0.0`** — `pressure` is undefined, not `inf`; carried as `None`.
  `ALEPH-PORT-3604` declares budgets of exactly `0.0` for its all-separated contact, so this is a real
  state and not a hypothetical.
- **an artefact absent from disk** — the bundle records it in `missing` with a reason and the figure
  draws the panel as absent. It never silently omits a panel.

Invariants that must hold, each with the test that asserts it:

- **I1.** Every case in the artefact is either a row of the matrix or named in `absent` with a reason,
  and the two sets are disjoint and cover — the partition `render.py` enforces for owners.
  `tests/viz/test_evidence.py::test_every_case_is_either_a_row_or_named_absent`
- **I2.** Every wrong kernel that is uncaught on both devices is either in the artefact's
  `blind_on_both_devices` list (**declared**) or rendered as a fault (**undeclared**). The two never
  collapse. `tests/viz/test_evidence.py::test_an_undeclared_blind_variant_is_a_fault_and_a_declared_one_is_not`
- **I3.** The CPU and CUDA columns are read from different device sub-dicts.
  `tests/viz/test_evidence.py::test_the_two_device_columns_are_not_the_same_column`
- **I4.** `ulp == 0.0` never receives a log position.
  `tests/viz/test_evidence.py::test_exact_zero_is_not_placed_on_the_log_axis`
- **I5.** The rendered SVG contains one mark per field row, and the count is non-zero.
  `tests/viz/test_evidence_figure.py::test_the_figure_draws_one_mark_per_field_row_and_the_set_is_not_empty`

## 6. Source evidence class and known retractions

No source, so no source claim to inherit. The **inputs** carry evidence classes and this entry
propagates them rather than restating them:

- `j1_law_parity.json` declares `evidence_rung: GPU_UNIT` and three `not_established` sentences,
  including *"Nothing here is evidence about a cell. A ULP figure is STRUCTURAL evidence about an
  implementation."* The figure prints that rung and those sentences **on the figure**, not in a
  docstring, because a picture travels further than the artefact it came from.
- `j1`'s `gate_comparison.note` records that `erm_tether_step[WRONG_NEWTON_MIRROR]` is a **declared**
  blind spot (`ALEPH-PORT-3604` §9a: bit-identical to the true kernel at every input by an IEEE-754
  identity). That is a retraction-shaped fact and it is rendered as declared, per I2.
- Looked for retractions in: `ALEPH-PORT-3601`…`-3610` §14 sections, `docs/ACTIVE_SESSIONS.md` lane
  closing blocks, and the artefacts' own `findings` / `not_established` arrays. `j2`, `j3`, `j4` carry
  `findings: []` and `verdict: PASS`; `g6_descent_window.json` carries the **disagreement**
  `ALEPH-PORT-3609` §14.6 raised and did not settle, and the figure shows it as a disagreement.

## 7. Independent oracle or derivation

A renderer's oracle is not a physical law — it is that **the picture and the artefact say the same
thing**. Four checks, none of which is "the SVG parses":

1. **Recount the gate comparison from the per-case rows and require it to equal the artefact's own
   `gate_comparison` block.** The artefact reports `caught_on_cuda = 38`, `caught_on_cpu = 38`,
   `wrong_kernels_total = 39`, and two empty set differences. The reader recomputes all five from the
   49 per-case records it actually drew, and disagreement is a refusal. This is the strongest control
   in the lane: it fails if the reader drops a case, double-counts one, or reads one device twice.
2. **The blind set, recomputed.** Recompute `{w : ulp(w, f) == ulp(TRUE, f) for every f}` and require
   the uncaught-on-both set to be a subset of it — a variant can be uncaught without being bit-identical
   (it moved but stayed inside its budget), but a variant bit-identical on every field **cannot** be
   caught, and if one is, the reader is misaligned.
3. **The log axis is an isotonic map.** `log_position` is strictly increasing on `(0, inf)`, checked
   on the artefact's own 130 measured values sorted — the property that makes "further right is worse"
   readable at a glance, and the one a bug in decade arithmetic destroys silently.
4. **The exact figures survive the round trip.** The worst per-field ULP, the case count and the
   converged fraction printed on the figure are re-derived from the artefact and compared to the
   values the runner recorded (`j4`: `12/12`, `1.5748e-05`).

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/viz/test_evidence.py::test_the_recounted_gate_comparison_equals_the_artefacts_own` | Recounting from the drawn rows reproduces `caught_on_cuda == 38`, `caught_on_cpu == 38`, `wrong_kernels_total == 39` and both empty set differences, on the real `j1` artefact. |
| Positive | `tests/viz/test_evidence.py::test_wrong_energy_half_is_over_budget_on_energy_and_bit_identical_on_both_forces` | The named case's `energy_pn_um` is `OVER_BUDGET` at `8.3886e6` ULP while `forces_a_pn` and `forces_b_pn` are `BLIND_TO_TRUE` — the finding, asserted per channel. |
| Positive | `tests/viz/test_evidence.py::test_wrong_frozen_normal_moves_one_force_field_and_nothing_else` | `cortex_membrane_contact[WRONG_FROZEN_NORMAL]`: `forces_cortex_pn` `OVER_BUDGET`, `energy_pn_um` **and** `forces_membrane_pn` `BLIND_TO_TRUE`. |
| Positive | `tests/viz/test_evidence.py::test_the_log_axis_is_strictly_increasing_over_the_measured_values` | `log_position` is isotonic on the 130 measured ULP values. |
| Positive | `tests/viz/test_evidence.py::test_a_thin_headroom_budget_is_marked_and_a_roomy_one_is_not` | A gate with `headroom < 2` is flagged; one with headroom above it is not — asserted on a non-empty set of each. |
| Positive | `tests/viz/test_evidence_figure.py::test_the_figure_draws_one_mark_per_field_row_and_the_set_is_not_empty` | Mark count equals row count, and both exceed zero. |
| Positive | `tests/viz/test_evidence_figure.py::test_the_descent_window_panel_shows_the_intersection_and_the_prior_cpu_window` | The `{1.0}` device window, the `{2.0}` prior CPU window and their **disagreement** all appear. |

## 9. Deliberately failing negative control

A renderer's negative controls have to attack the two ways it fails *quietly*: drawing nothing and
reporting success, and asserting over an empty set.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/viz/test_evidence.py::test_a_reader_that_drops_a_case_is_refused` | An artefact with one case removed is **refused**, because the recount no longer matches `gate_comparison`. This is the control that kills "drew nothing, exited 0". |
| Negative (must fail) | `tests/viz/test_evidence.py::test_an_undeclared_blind_variant_is_a_fault_and_a_declared_one_is_not` | Moving `WRONG_NEWTON_MIRROR` out of `blind_on_both_devices` turns it into a rendered fault; leaving it in does not. |
| Negative (must fail) | `tests/viz/test_evidence.py::test_the_two_device_columns_are_not_the_same_column` | Feeding the *same* device dict as both CPU and CUDA is refused. Without this, a reader that read `cuda:0` twice would report "no cross-device difference" forever and look perfect. |
| Negative (must fail) | `tests/viz/test_evidence.py::test_exact_zero_is_not_placed_on_the_log_axis` | `log_position(0.0)` raises rather than returning the axis origin. |
| Negative (must fail) | `tests/viz/test_evidence.py::test_a_family_with_no_true_variant_reports_unknown_not_not_blind` | Removing a family's `TRUE` case makes every blind verdict `None`, never `False`. |
| Negative (must fail) | `tests/viz/test_evidence_figure.py::test_an_empty_matrix_is_refused_and_never_written_as_a_figure` | A parity artefact with zero cases exits `EXIT_REFUSED` and writes **no** file. |
| Negative (must fail) | `tests/viz/test_evidence_figure.py::test_a_missing_artefact_is_drawn_as_absent_and_never_omitted` | With `g6_descent_window.json` absent the panel is rendered as absent with a reason; the figure does not quietly lose a panel. |
| Vacuity | `tests/viz/test_evidence.py::test_every_control_set_in_this_module_is_non_empty` | Each set the controls above assert over (`wrong kernels`, `blind field rows`, `thin budgets`, `exact-zero rows`) has a **stated, asserted cardinality**, so `np.all([])` cannot pass one of them. `ALEPH-PORT-3606` §13a note 3 is why this test exists. |

## 10. Numerical and precision envelope

Working precision is float64 throughout, on numbers that arrived as JSON doubles. The reader performs
**no** arithmetic on the physics: `ulp`, `roundoff_ulp_bound` and `effective_ulp_budget` are read, not
recomputed, so this entry cannot disagree with the runner about a measurement — only about how it is
drawn. The two derived ratios (`pressure`, `headroom`) are single divisions, and `pressure` is
**cross-checked** against the artefact's own `pressure` field to `1e-9` relative where present; a
disagreement is a refusal, because it means the reader has paired a ULP with the wrong budget.

Bit-identity (§4a) is tested with `==` on float64 values that came from the same producer in the same
process — the comparison `ALEPH-PORT-3604` used for the exact zero, for the same reason: exactness is
the evidence, and `isclose` would convert an exact statement into an approximate one.

`log10` is applied only to values `> 0` (I4). Below `10**DECADE_LOW = 10**-2` the mark is clamped and
carries `clamped=True`; the smallest measured value is `0.0159`, so clamping is live in this data and
not a dead branch.

## 11. Production-backend residency and transfer

**This never runs on the production backend.** It is a viewer: host-only, NumPy-free on the hot path
(the matrix is built from Python floats), and it holds no device array at any point. `aleph/viz/**`
is forbidden a module-scope accelerator import by `tests/viz/test_import_boundaries.py`, which globs
the package and therefore covers both new modules without being edited.

Transfer is one direction and happens **outside** the process where possible: `discover_artefacts`
prefers a local `~/.aleph_runner/artefacts/`, and falls back to `scp` from `gbook` in a subprocess. No
GPU is touched, no `warp` is imported, and `scripts/gpu_preflight.py` is irrelevant here because
nothing is scheduled — reading a JSON file a run already wrote is not a run.

## 12. Comments and docstrings to discard

No source prose exists to discard, since nothing was ported. What must **not** appear in this module's
output, and what replaces it:

- **A single verdict per case.** Replaced by the per-channel matrix — the whole point of §4a.
- **"PASS" as a headline.** The artefact's `verdict: PASS` is shown *attributed to the artefact*, next
  to the `GPU_UNIT` rung and the three `not_established` sentences, never as the figure's own claim.
- **Any ULP figure without its budget.** A bare `326.66` is unreadable; every mark is drawn against the
  gate it was judged by.
- **Any implied statement about a cell.** `j1`'s own `not_established[0]` is rendered verbatim.
- **"38/39 caught" without the 39th.** The declared blind spot is named on the figure with its ledger
  citation, because a reader who cannot see why one is uncaught will assume a defect.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | `2026-08-01`. `tests/viz`: **329 passed** (was 277 before this entry), zero skipped. Every control named in §8–§9 resolves to a real test, checked by name. Mutation study: **13 planted defects, 13 killed**, second pass from a bytecode tree verified empty under `aleph/viz/` and `tests/viz/`. The command runs end to end with **no arguments**, fetching the artefacts off `gbook` over ssh; the fetched files are **byte-identical** to the copies committed under `docs/results/2026-08-01-cuda-evidence/artefacts/`, so the committed evidence is the run. Refusals exit `2` with the reason. **Zero GPU jobs run.** |
| Reviewer | **Agent-proposed. Unratified.** No `decided_by` field appears in this entry by construction. |
| Rollback | Delete `aleph/viz/evidence.py`, `aleph/viz/evidence_figure.py` and their two test files, and revert the append-only names in `aleph/viz/__init__.py`. Nothing else imports them — `figure.py` is untouched. The CUDA evidence returns to being readable only as JSON. |

### 13a. The mutation study

| # | Planted defect | Killed by |
|---|---|---|
| M1 | the renderer draws no marks and reports success | `test_evidence_figure.py::test_the_figure_draws_one_mark_per_field_row_and_the_set_is_not_empty` |
| M2 | both device columns read from `cuda:0` | `test_evidence.py::test_the_two_device_columns_are_not_the_same_column` |
| M3 | one case silently dropped from the matrix | `test_evidence.py::test_the_recounted_gate_comparison_equals_the_artefacts_own` |
| M4 | `BLIND_TO_TRUE` collapsed into `WITHIN_BUDGET` | `test_evidence.py::test_wrong_energy_half_is_over_budget_on_energy_and_bit_identical_on_both_forces` |
| M5 | `log_position` returns the axis origin for exact zero | `test_evidence.py::test_exact_zero_is_not_placed_on_the_log_axis` |
| M6 | a **linear** axis instead of a logarithmic one | `test_evidence.py::test_the_axis_is_logarithmic_and_not_merely_increasing` |
| M7 | every blind variant treated as declared | `test_evidence.py::test_an_undeclared_blind_variant_is_a_fault_and_a_declared_one_is_not` |
| M8 | an absent field skipped instead of listed | `test_evidence.py::test_a_field_present_for_one_variant_and_missing_for_another_is_rendered_absent` |
| M9 | the CPU column's blind verdict judged against the **CUDA** baseline | `test_evidence.py::test_the_blind_set_is_the_same_on_both_devices` |
| M10 | no gate ever counted as thin | `test_evidence.py::test_a_thin_headroom_budget_is_marked_and_a_roomy_one_is_not` |
| M11 | a refusal still writes the figure and exits 0 | `test_evidence_figure.py::test_an_empty_matrix_is_refused_and_never_written_as_a_figure` |
| M12 | bit-identity by `isclose` instead of `==` | `test_evidence.py::test_a_nearly_identical_field_is_not_reported_as_bit_identical` |
| M13 | the recount check disabled | `test_evidence.py::test_a_gate_comparison_that_disagrees_with_the_rows_is_refused` |

**Two of these were found by the study rather than confirmed by it, and both are worth the space.**

**M9 was a real defect in this module, not a planted one.** The first working version judged the CPU
column's bit-identity against the **CUDA** true-kernel baseline. It produced a plausible figure: 95 of
130 channels differ between devices even for a true kernel, so almost every CPU field read "not
blind" and the finding quietly vanished from half the picture. It was caught by asking what the CPU
column was being compared *to* — a cross-device comparison wearing a same-device label. The control
that now holds it is that the blind set is **identical on both devices** (34 channels), which is a
strong statement precisely because the underlying ULP figures are not.

**M12 survived the first pass.** Swapping `==` for `math.isclose(rel_tol=1e-3)` changed **no verdict
on this artefact**, because the blind fields are exactly equal and every other field differs by orders
of magnitude. §10 claims that exactness is the evidence; the control set did not back that claim,
because the looseness was never exercised. The repair drives the boundary directly — a field nudged
one part in a million away from the true kernel's must read *not* bit-identical — and it kills M12.
This is the sixth entry in the project's running series on controls that are silent, and its shape is
new: **the control was not vacuous and not wrong; it was never exercised at the boundary it claimed
to defend.**

## 14. Honest limits

1. **This entry establishes nothing about a cell, and nothing about physics.** It is a reading of a
   structural measurement. `j1`'s own `not_established` is reproduced on the figure precisely so the
   picture cannot outrun it.
2. **Bit-identity is inferred from equal ULP, not from the arrays.** §4a argues why equal `ulp` implies
   equal arrays *for two variants of one family on one state*, and the argument is sound only because
   the reference scale is a property of the family and state. The arrays themselves are not in the
   artefact, so this is `INHERITED_UNVERIFIED` at the array level. A future runner that emitted a hash
   per field would upgrade it, and that is a one-line change to `ALEPH-PORT-3607`'s job — **not made
   here**, because `scripts/**` belongs to another lane.
3. **`j5_precision_modes.json` and `j3_tension_identity.json` are read but not rendered in this pass.**
   Both are per-mode tables that deserve their own panel; drawing them badly would be worse than the
   bundle naming them as present-but-not-drawn, which is what it does. Named here rather than left to
   be discovered as a gap.
4. **The figure is a projection selected for legibility.** With 130 field rows the labels are elided at
   the right edge; the SVG carries `<title>` on every mark so the full text is recoverable on hover, and
   `--json` prints the complete table. The picture is not the record — `ALEPH-PORT-3602` set that limit
   and it holds here.
5. **The descent-window panel renders a disagreement this lane cannot resolve.** `ALEPH-PORT-3609`
   §14.6 measured `{1.0}` on CUDA against `-3608`'s derived-and-ratified `2·u32`, named three candidate
   repairs and took none. The figure shows both windows and marks them as disagreeing. **It does not
   pick one**, and a renderer that did would be settling a PI question in a picture.
6. **No GPU was used and none may be.** Every number here was produced by Lane G7's runner on
   2026-07-31; this lane re-reads them. If the artefacts are stale relative to the kernels, the figure
   is stale, and the only thing it can do about that is print `finished_utc` and the device identity on
   the figure — which it does.
7. **The ssh fallback trusts the `gbook` host alias.** If the alias resolves elsewhere, the figure
   describes another machine's run. The host recorded *in the artefact* is printed on the figure, so a
   mismatch is visible, but it is not enforced.
8. **This run contains no exact-zero channel, so `EXACT` is a state with no data behind it.**
   `ALEPH-PORT-3604`'s exact zeros live in its all-separated contact state, which `j1` does not drive.
   The style, the lane and the refusal in `log_position` are all exercised by unit controls, but
   **no control asserts over `ParityEvidence.exact_rows`**, because such a control would pass
   vacuously today. `test_every_control_set_in_this_module_is_non_empty` asserts that set is
   *empty* and says why, so the day a run does produce one, the assertion fails and somebody writes
   the control rather than the existing ones silently starting to cover it. This is a live instance
   of `ALEPH-PORT-3606` §13a note 3 handled in advance instead of after.
9. **The figure renders the CUDA column's channel state and carries the CPU column as a second mark
   only.** `FieldRow.state_cpu` is computed, controlled and available in `--json`, but the row's
   colour is the CUDA verdict. On this artefact the two never disagree — the blind set is identical
   and no variant is caught on one device and missed on the other — so nothing is currently hidden.
   **If a future run produces a disagreement, the row colour would show only one side of it.** The
   summary band counts such cases and the faults list would name them, so it cannot pass unnoticed,
   but the per-row encoding is not yet two-sided. Named here rather than left to be discovered.
10. **`j1` reports one position-precision mode (`global_f32`).** The per-mode comparison that
    `ALEPH-PORT-3605` §13.4 and `-3606` §13.4 disagree about lives in `j5_precision_modes.json`,
    which this pass reads and does not render (§14.3). The figure therefore says nothing about
    option D of the position-precision proposal, and should not be read as evidence for or against
    it.
