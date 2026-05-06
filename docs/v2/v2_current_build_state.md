# V2 Current Build State

Last updated: 2026-05-06 KST (post-Option-A+B repo restructure merge `0533422`).

This is the canonical handoff page for restarting ActiveCellSim simulation
construction work in Claude/Codex workrooms. It is intentionally compact:
read this first, then progressively load the referenced docs and code.

## Executive State

ActiveCellSim new development is on **v2**, an image-constrained,
cell-resolved mechanobiology framework. The v1 continuum/MPM spheroid-first
prototype is preserved as historical/reduced-surrogate work and should not
receive new biology unless explicitly requested.

Top-level layout (post-Option-A+B):

- `acs/v2/` — v2 code home (no v1/v2 directory split inside `acs/`; the
  v1 vs v2 split lives at `acs/v2/` vs the rest of `acs/`).
- `docs/v1/` — frozen v1 docs (roadmap, outcome reports, sanity gates).
- `docs/v2/` — active v2 docs (vision, roadmap, locked specs, sanity gates,
  briefs, this file).
- `tests/v1/` and `tests/v2/` — corresponding test split.
- `STRUCTURE.md` — root file→bucket map.
- `docs/claude_codex_log.md` — shared Claude/Codex ledger.

Pointers:

- v2 roadmap: `docs/v2/10_dev_roadmap_v2.md`
- v2 vision/framing: `docs/v2/00_project_vision_v2.md`
- v1 frozen reference: `docs/v1/v1_continuum_backup.md`

## Recent seals and milestones (since 2026-05-03)

- 2026-05-04 — `acs/v2/metrics.py` pure metric registry added
  (`75824af`); spheroid CSV metrics extension (`6e48f18`).
- 2026-05-05 — V2-1 imaging input contract sealed
  (`docs/v2/v2_layer_1_imaging_input_contract_locked.md`,
  reconciliation in `b463616`, `40ed049`, `79376e7`).
- 2026-05-06 — V2 imaging measurement-protocol gate sealed
  (`docs/v2/v2_imaging_measurement_protocol_gate_locked.md`,
  HEAD-side cleanup `49e246c`). 260313 reframed as multi-cell /
  spheroid-stage validation reservoir for V2-3/V2-5; **NOT** a
  Layer 1 single-cell calibration source. Hard Rule 11 contract.
- 2026-05-06 — Option A+B repo structure cleanup merged (`0533422`):
  docs split to `docs/v1/` + `docs/v2/`, tests split to
  `tests/v1/` + `tests/v2/`, root `STRUCTURE.md` added, path-reference
  rewrites throughout. No `acs/v2/` code moves.
- 2026-05-06 (active) — Layer 1 Sanity-Gate reference anchors draft
  (`docs/v2_layer_1_sanity_gate_anchors_draft.md`, untracked at root,
  to be relocated to `docs/v2/...` when promoted) cross-posted to
  `design-discussion` as a debate seed (`mcp_msg:2500`). Promotion
  blocked on Q1 (boundary representation), Q2 (single-cell dataset),
  Q3 (FA-channel inclusion in first runnable Layer 1 build).

The full set of locked specs and sanity gates lives in `docs/v2/`. As of
2026-05-06 the surface includes (non-exhaustive): V2-1 imaging input
contract; V2 imaging measurement-protocol gate; adhesion network
layer; protrusion-coupled FA (63b); closed-loop ECM phased plan;
focal adhesion dynamics result-typed; hard blockers 1/2/3/4/5 (each
with brief + locked + sanity gate); item 5 sweep harness; P1
active-contour derivation; phase D no-op scaffolding; phase E
composition + step 2; phase F minimal cell motility pilot; phase F
proper-1 FA rate response; sister-gate-mirror audit. New code work
should always cross-check the relevant locked spec + sanity gate
before editing.

## Code state pointers (as of 0533422)

`acs/v2/` is no longer at "schema seed" stage. Active modules:

- Contract / typing: `data_contract.py`, `single_cell.py`,
  `measurement_boundary.py`, `metrics.py` (pure registry).
- Imaging contract / I/O: `imaging_contract/` (area_extractor,
  loader, split_builder).
- Active contour / boundary: `active_contour.py`,
  `active_contour_harness.py`, `dynamics/active_contour.py`.
- Protrusion / FA / adhesion: `protrusion.py`, `focal_adhesion.py`,
  `adhesion_network_state.py`, `dynamics/focal_adhesion.py`,
  `dynamics/adhesion_network_dynamics.py`,
  `dynamics/protrusion_coupled_focal_adhesion.py`,
  `dynamics/fa_rate_response.py` (**WIP HALTED — implementation
  coherent but untested; lock §4/§5 test catalog is not implemented;
  Phase F proper resume requires explicit PI directive after Option 1
  pivot. See `docs/v2/v2_phase_f_proper_1_fa_rate_response_locked.md`
  STATUS banner and commit `a1d97e8`. Halt approved post-merge in
  claude-work mcp_msg:2550, codex mcp_msg:2551.**),
  `dynamics/fa_to_ecm_scattering.py`,
  `dynamics/ecm_to_fa_bias.py`.
- ECM: `ecm_substrate.py`, `ecm_open_loop_harness.py`,
  `dynamics/ecm_constitutive_response.py`,
  `dynamics/ecm_lyapunov_metric.py`,
  `dynamics/ecm_open_loop.py`.
- Closed-loop / phased composition:
  `dynamics/closed_loop_phase_d.py`,
  `dynamics/closed_loop_phase_e.py`,
  `dynamics/closed_loop_phase_e_sweep.py`,
  `phase_e_v2_pilot_runner.py`,
  `phase_f_minimal_motility_pilot_runner.py`,
  `dynamics/phase_f_minimal_motility.py`.
- Other: `cell_cluster.py`, `junction.py`, `output/frame_dump.py`,
  `viz/`.

Test surface (`tests/v2/`): 718 passed + 1 skipped on `0533422`
(`.venv-collab/bin/python -m pytest tests/v2`, ~21 min on Mac
without GPU). The single skip is the opt-in 4551-row imaging
measurement-protocol gate (set `RUN_FULL_REPRODUCTION_GATE=1` to
run the full subprocess, requires ~7 GB of mask PNGs).

## Active Roadmap

Follow `docs/v2/10_dev_roadmap_v2.md`, not the original v1 roadmap.

Phase status (revised against actual code/lock state, not the original
roadmap text):

- **Phase V2-0**: complete. v1 frozen, v2 package open, AGENTS/CLAUDE
  point to v2.
- **Phase V2-1** (imaging data contract): **sealed** at
  `v2_layer_1_imaging_input_contract_locked.md`. Canonical layout,
  metadata fields, segmentation provenance, calibration/validation
  split all addressed in lock; metric registry is implemented in
  `acs/v2/metrics.py`. Example YAML/JSON contract present.
- **Phase V2-2** (single-cell minimal mechanobiology): partially
  implemented. Single-cell state object exists; protrusion + focal
  adhesion + adhesion-network + ECM bias + closed-loop phase D/E/F
  dynamics modules exist with locked specs and sanity gates. Active
  contour boundary is in place. Layer 1 Sanity-Gate reference anchor
  contract is **draft only** (working reference at
  `docs/v2_layer_1_sanity_gate_anchors_draft.md`); promotion gated on
  design-discussion Q1-Q3.
- **Phase V2-3** (cell-resolved multi-cell assembly): not implemented.
  260313 reservoir is reserved as one possible multi-cell calibration
  set per the imaging measurement-protocol gate.
- **Phase V2-4** (ECM fiber network) and **Phase V2-5** (brute-force
  spheroid reference): not implemented.

When in doubt, do not jump directly to spheroid dynamics. The Layer 1
single-cell anchor work is still incomplete (Q1-Q3 unresolved).

## Modeling Rules That Matter Most

- Measurement modality must match the experiment or image source.
  Top-down projected area means top-down projected polygon/convex-hull
  measurement, not substrate-contact area (Hard Rule 11).
- Separate parameter roles explicitly: literature-fixed,
  fitted/calibrated, exploratory, diagnostic. v2 allows
  image-constrained calibration only when the parameter is
  interpretable and validation is held out.
- No scalar fudge factors. Any empirical constant must pass the
  Magic-Number Block in `AGENTS.md` / `CLAUDE.md`.
- Any new physics or numerics module must pass and record the Sanity
  Gate Protocol before first execution.
- Literature anchors for Layer 1/2 are in
  `docs/v2/10_dev_roadmap_v2.md`; do not use biology papers to
  justify numerical discretization choices.
- 260313 PI imaging is **multi-cell / spheroid-stage** (V2-3/V2-5
  reservoir); it is **not** a Layer 1 single-cell calibration source.
- GPU portability remains mandatory. Do not hard-code device IDs or
  paths.

## Workroom Routing

Use the two-room layout:

- `design-discussion`: model design, literature anchors, sanity-gate
  review, measurement protocol, and unresolved modeling choices.
- `implementation-work`: code edits, tests, data pipeline work,
  benchmarks, and execution dispatch.

Each room has exactly one Claude session and one Codex session. There
is no separate chat/work LLM split anymore.

## Bootstrap protocol for fresh sessions

Per `CLAUDE.md`'s "Session bootstrap (acs-collab MCP)" section, every
fresh session — including `claude --continue` resumes, brand-new
desktop-app threads, and any thread where the assistant has no
in-memory record of recent collab traffic — calls
`mcp__acs-collab__bootstrap(limit=50)` before processing the first
user prompt.

After bootstrap, load (in order, only what's relevant):

1. `AGENTS.md` (Codex) or `CLAUDE.md` (Claude).
2. `docs/v2/00_project_vision_v2.md`.
3. `docs/v2/10_dev_roadmap_v2.md`.
4. The locked spec(s) + sanity gate(s) for the unit at hand.
5. Relevant files under `acs/v2/`.
6. Recent entries from `docs/claude_codex_log.md`.

Do not rely on MCP history alone for the full simulation plan.

## Next Concrete Work

Good next tasks for `implementation-work`:

- Promote Layer 1 Sanity-Gate reference anchors (`docs/v2_layer_1_sanity_gate_anchors_draft.md`)
  by closing promotion criteria (a)-(f) — this requires
  design-discussion Q1-Q3 resolution + JCR cross-check + path-stability
  sweep + PI sign-off.
- Watch for and address Codex correctness-review BLOCKERs on existing
  V2-2 dynamics modules; the dynamics surface is large and many of
  the locked specs are recent.
- Cross-check that every `acs/v2/dynamics/*.py` module is satisfying
  its corresponding sanity gate doc before any new run.

Good next tasks for `design-discussion`:

- Resolve Layer 1 Q1: boundary representation choice (vertex /
  phase-field / hybrid) — this gates how anchors 1.4 (per-event
  protrusion) and 1.5 (FA per-state) are operationally applicable.
- Resolve Layer 1 Q2: select a single-cell imaging dataset (260313 is
  ineligible per the Hard Rule 11 contract) — without this the Layer 1
  Sanity Gate cannot be calibrated even literature-fixed-only.
- Resolve Layer 1 Q3: include FA state machine in first runnable
  Layer 1 build, or defer to Layer 1+? The TFM/paxillin channel
  requirement falls out of this decision.

When (γ) draft is promoted, this page should reference the locked
spec by path under `docs/v2/...` rather than the working draft path.
