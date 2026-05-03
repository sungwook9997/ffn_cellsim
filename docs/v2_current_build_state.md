# V2 Current Build State

Last updated: 2026-05-03 KST.

This is the canonical handoff page for restarting ActiveCellSim simulation
construction work in Claude/Codex workrooms. It is intentionally compact:
read this first, then progressively load the referenced docs and code.

## Executive State

ActiveCellSim new development is on **v2**, an image-constrained,
cell-resolved mechanobiology framework. The v1 continuum/MPM spheroid-first
prototype is preserved as historical/reduced-surrogate work and should not
receive new biology unless explicitly requested.

Current implementation state:

- v2 code home: `acs/v2/`
- roadmap: `docs/10_dev_roadmap_v2.md`
- vision/framing: `docs/00_project_vision_v2.md`
- v1 frozen reference: `docs/v1_continuum_backup.md`
- shared Claude/Codex ledger: `docs/claude_codex_log.md`

The current v2 code is still at the schema/contract seed stage:

- `acs/v2/data_contract.py` defines imaging dataset specs, metric specs, and
  a minimal single-cell contract.
- `acs/v2/single_cell.py` defines single-cell state, protrusion events, focal
  adhesion observations, validation, and projected polygon area.
- `tests/test_v2_contract.py` covers basic contract validation and single-cell
  projected-area/nested-state checks.

No dynamics solver should be assumed to exist yet.

## Active Roadmap

Follow `docs/10_dev_roadmap_v2.md`, not the original v1 roadmap.

Current phases:

- Phase V2-0: mostly complete. V1 is frozen/labeled; v2 package is open.
- Phase V2-1: in progress. Imaging data contract exists, but canonical data
  layout, metadata examples, metric registry, and example YAML/JSON contract
  are still missing.
- Phase V2-2: started only as state schemas. Single-cell dynamics, boundary
  representation, event extraction, focal adhesion state machine, polarity /
  traction representation, and dashboard are not implemented.
- Phase V2-3 and later are future work and must build on V2-1/V2-2 decisions.

When in doubt, do not jump directly to spheroid dynamics. Establish the
single-cell and imaging-contract ground truth first.

## Modeling Rules That Matter Most

- Measurement modality must match the experiment or image source. Top-down
  projected area means top-down projected polygon/convex-hull measurement, not
  substrate-contact area.
- Separate parameter roles explicitly: literature-fixed, fitted/calibrated,
  exploratory, diagnostic. v2 allows image-constrained calibration only when
  the parameter is interpretable and validation is held out.
- No scalar fudge factors. Any empirical constant must pass the Magic-Number
  Block in `AGENTS.md` / `CLAUDE.md`.
- Any new physics or numerics module must pass and record the Sanity Gate
  Protocol before first execution.
- Literature anchors for Layer 1/2 are in `docs/10_dev_roadmap_v2.md`; do not
  use biology papers to justify numerical discretization choices.
- GPU portability remains mandatory. Do not hard-code device IDs or paths.

## Workroom Routing

Use the two-room layout:

- `design-discussion`: model design, literature anchors, sanity-gate review,
  measurement protocol, and unresolved modeling choices.
- `implementation-work`: code edits, tests, data pipeline work, benchmarks,
  and execution dispatch.

Each room has exactly one Claude session and one Codex session. There is no
separate chat/work LLM split anymore.

## Known Context Gap Fixed By This Page

On 2026-05-03, Claude audited its memory/MCP state and reported that its
persistent memory mostly contained workroom/collaboration infrastructure, not
the actual v2 simulation-building context. The v2 docs and `acs/v2/` files
were present on disk, but not loaded into the active Claude context.

Therefore, any fresh simulation-construction turn should begin with this page,
then load:

1. `AGENTS.md`
2. `docs/00_project_vision_v2.md`
3. `docs/10_dev_roadmap_v2.md`
4. relevant files under `acs/v2/`
5. relevant recent entries from `docs/claude_codex_log.md`

Do not rely on MCP history alone for the full simulation plan.

## Next Concrete Implementation Work

Good next tasks for `implementation-work`:

- add an example v2 data-contract YAML/JSON file and loader tests
- define a metric registry for single-cell morphology metrics
- add a dedicated `acs/v2/metrics.py` with projected area, perimeter,
  roughness, and validation-safe naming
- draft but do not execute a V2-2 single-cell active-contour dynamics module
  until its Sanity Gate is written

Good next tasks for `design-discussion`:

- resolve the exact Layer 1 boundary representation: active contour,
  vertex model, phase-field, or hybrid
- define the minimum imaging channels required before fitting any protrusion
  or focal-adhesion rates
- pre-register which observables validate protrusion, polarity, traction, and
  focal adhesion dynamics
