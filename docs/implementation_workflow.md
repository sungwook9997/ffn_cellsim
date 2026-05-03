# Implementation Workflow SOP — Claude (main) + Codex (review)

**Status:** active operating model for `implementation-work` workroom
(2026-05-03, MCP ids: PI 567/572/575/578, Claude 570/580/582/585,
Codex 568/574/576/579/583).
**Scope:** ActiveCellSim v2 (`acs/v2/`, `docs/10_dev_roadmap_v2.md`).
**Mirrors:** `docs/workroom_layout.md` (room split SOP), `CLAUDE.md`
(Sanity Gate, Magic-Number Block, Hard Rules).

## 1. Roles

| Agent | Responsibility | Forbidden |
|---|---|---|
| **Claude** (coding main) | v2 module design + implementation; Sanity Gate authoring; smoke tests (≤2 GB, dev.yaml); commit hygiene; production handoff drafting | Approving its own diff for production; modifying gate tolerances to make a run pass; introducing magic numbers without all three Magic-Number Block tests recorded |
| **Codex** (review/debug) | Independent diff review; Sanity Gate 6-check re-verification; conservation/edge-case probing; additional test authoring; reproducibility audit (commit/config/env/host/GPU); magic-number block enforcement | Re-implementing Claude's modules from scratch; silently fixing Claude's diff (always reply through MCP); rubber-stamping without running the smoke test it is reviewing |
| **PI** | Design decisions; Sanity Gate failure resolution; production-run authorization; experimental-data overlay direction | None — but PI is in `design-discussion` by default; do not page PI for routine review work |

PI is the only authority for: design changes, gate-contract changes,
production-run go/no-go, deviation from any `CLAUDE.md` Hard Rule.

## 2. Agent-to-agent MCP discipline

The `tmux_relay` runs with `--include-llm-messages` in both rooms
(see `tools/collab_mcp/tmux_relay.py` `route_targets`, line 218).
A message with `to=claude` or `to=codex` lands in the partner's tmux
pane wrapped exactly like a PI message (`: mcp_msg id=...`). Use this.

**Direct channel** (`to=claude` / `to=codex`):
- Implementation/review questions, diff feedback, smoke results,
  Sanity-Gate cross-check, test additions, "is this clean enough to
  hand off?" questions.
- These do **not** wait for PI. The relay wakes the partner immediately.

**PI escalation** (`to=pi` with `status=open-question` or
`decision-needed`):
- Sanity Gate FAIL on a physics/numerics module (per `CLAUDE.md`
  failure-handling rule — surface with at least three concrete options).
- Magic-number test 3 = yes, or tests 1/2 = no.
- Gate-contract change request (never edit gate tolerances inline).
- Production-run authorization.
- Disagreement between Claude and Codex that does not converge after
  one review cycle.
- Anything that touches `data/experimental/*` semantics.

**Visibility-only to PI** (`to=pi` with `status=FYI`):
- Cycle complete + handoff ready (one line + commit hash).
- Smoke test result one-liner if PI asked for status.

Default: keep PI off the timeline for routine review.

## 3. The 4-step cycle

A "unit of work" is one logically-coherent module change (e.g. one
physics term, one analysis function, one config field). Bigger
changes are split into multiple units.

### Step 1 — Code (Claude)
1. Pick the unit from `docs/10_dev_roadmap_v2.md` or PI directive.
2. Write code under `acs/v2/...`.
3. Pre-execution checks (mandatory before first run):
   - Sanity Gate 6 checks recorded as docstring section or
     `*_sanity.md` sibling. FAIL → halt, escalate.
   - Magic-Number Block answered for any new constant. Test 3 = yes
     or 1/2 = no → halt, escalate.
4. Commit. Atomic. Message names the module + cycle id.

### Step 2 — Self-check (Claude)
1. Run unit tests: `pytest tests/<relevant>` (must pass).
2. Run dev smoke: `ACS_GPU_BACKEND=auto python -m acs.v2.<entrypoint>
   --config configs/dev.yaml --frames 10` (≤ 2 GB, ≤ a few minutes).
3. Capture: peak VRAM (`nvidia-smi` snapshot), wall-clock,
   conservation diagnostic (mass/momentum drift), Sanity Gate
   confirmation lines.
4. Post a single MCP message `to=codex` with: commit hash, smoke
   metrics, Sanity Gate PASS lines, any concerning observation. No
   long log dumps — link to file path if needed.

### Step 3 — Review (Codex)
1. Read the diff (`git show <hash>` or `git diff <prev>..<hash>`).
2. Re-verify Sanity Gate 6 checks against the diff. Walk the
   measurement-protocol consistency check (Sanity Gate item 6) — do
   not accept "exact at the peak" proofs that are evaluated off-peak.
3. Independently reason the dimensional balance for any new term
   (Hard Rule 10: same unit basis for any order-of-magnitude
   comparison).
4. Verify measurement-protocol matches experimental modality where
   relevant (Hard Rule 11: PI A/A₀ = top-down xy convex hull, not
   substrate-contact area).
5. Spot-run the smoke or write an additional probe test if any check
   is uncertain. Record exact command in the review reply.
6. Reply MCP `to=claude`. Status: `ack` (accept), `open-question`
   (need clarification), or `decision-needed` (Sanity / magic-number
   / gate-contract issue → also CC `to=pi`).

### Step 4 — Converge
- Codex `ack` → Claude marks the unit complete. Both ids recorded in
  the unit's commit message trailer (`Reviewed-By: codex (mcp=<id>)`).
- Codex `open-question` → Claude addresses, posts new commit hash, go
  to step 3.
- Codex `decision-needed` → both wait for PI.

A cycle that bounces step 3 ↔ step 1 more than twice without
converging is itself a `decision-needed` — escalate to PI rather
than spiraling.

## 4. Sanity Gate and Magic-Number Block

Both are defined in `CLAUDE.md`. This SOP enforces *when* they are
applied:

- **Before first execution** of any new physics/numerics module —
  Sanity Gate 6 checks recorded.
- **Before any new tunable constant lands in code** — all three
  Magic-Number Block tests answered. Test 3 = yes or tests 1/2 = no
  → blocker.
- **At Codex review** — both are re-verified independently. Codex
  treats Claude's gate record as a claim, not a guarantee.
- **At gate failure** — never edit the gate (Hard Rule). Surface
  with three options, wait for PI direction.

## 5. Production handoff format

PI's desktop runs production. Implementation-work hands off via a
single brief file at `runs/<UTC-yyyymmdd-HHMM>_<topic>.md`:

```
# Production handoff — <topic>

## Reproducibility
- git commit: <40-char hash>
- git status --short: <output, must be empty for production>
- branch: <branch name>
- repo path: <absolute path on origin host>

## Configuration
- config: configs/production_16gb.yaml  # or _24gb.yaml
- entry command: ACS_GPU_BACKEND=cuda python -m acs.v2.<entrypoint> --config <path>
- expected wall-clock: <hours>
- expected VRAM peak: <GB>
- frame interval: <minutes>

## Target environment
- target host: <hostname>            # e.g. Gbook2 via `ssh win`
- nvidia-smi: <name, memory.total, memory.free, driver>
- ACS_GPU_BACKEND: <cuda|vulkan|opengl|metal|cpu|auto>
- key env vars: <other relevant TI_*, OMP_*, etc.>

## Outputs
- output path: <absolute path on target>
- key metrics: <list>
- Sanity Gate (re-confirmed at handoff): 6/6 PASS — see <path>

## Sign-off
- Claude: mcp id=<id>
- Codex: mcp id=<id>
```

### Clean-tree rule (production / sweep)
`git status --short` MUST be empty. Dirty tree → handoff refused;
commit or stash & branch first. Reason: production runs are recorded
against a commit hash; uncommitted changes invalidate that record.

### Dev-smoke exception
Dev smoke / debug runs MAY proceed dirty, but the handoff (if any)
records `git diff --stat` and the touched files explicitly, and is
flagged "not promotable to production". Promotion requires re-running
on a clean commit.

### Always-blocker
Any dirty file under `data/experimental/` is a hard blocker
regardless of run type. Per `CLAUDE.md` Hard Rule, those CSVs are
read-only.

## 6. Sanity Gate / Magic-Number escalation

When Claude (step 1 or 2) or Codex (step 3) hits a FAIL:
1. Stop further code work on the module.
2. Post `to=pi` with `status=decision-needed` containing:
   - what failed and on which check / test
   - the dimensional or conservation reasoning behind the FAIL
   - **at least three concrete options** (e.g. reduce Δt, switch to
     implicit, switch to overdamped; or: derive principled term,
     drop the feature, surface to literature lookup)
   - which option the implementer recommends and why
3. Wait for PI direction. Never silently work around.

Example anti-patterns to never repeat (from `CLAUDE.md`):
- `csf_kappa_scale=0.05` (April 2026): magic number chosen to make
  gate pass.
- `g_star/κ` per-volume vs per-area mismatch (Path C, 2026-04-29):
  unit-mismatched dimensional comparison.
- v13 `_build_curvature` band-averaged vs analytical-peak proof
  (April 2026): off-protocol measurement.

## 7. Heavy-load mode (parallel sessions)

When PI opens additional sessions for sweep throughput:
- Claude responsibility: split the sweep config (`configs/sweep_*.yaml`
  or condition × repeat matrix) into N batches by parameter, not by
  random hash — keep each session's batch reproducible standalone.
- Codex responsibility: review remains single-Codex by default.
  Reviews are sequential per unit, not per session.
- Heavy-load handoff records a sweep manifest
  `runs/<UTC>_sweep_<topic>.md` listing every per-session
  reproducibility block from §5.
- Heavy-load mode terminates on PI signal; revert to single
  Claude + single Codex per CLAUDE.md "1Claude+1Codex per room"
  invariant (commit `1e77485`).

## 8. References

- `CLAUDE.md` — vision, Hard Rules, Sanity Gate, Magic-Number Block.
- `docs/00_project_vision_v2.md` — v2 framing (image-constrained,
  cell-resolved).
- `docs/10_dev_roadmap_v2.md` — current implementation roadmap.
- `docs/workroom_layout.md` — 2-room workroom SOP.
- `tools/collab_mcp/tmux_relay.py` — bilingual relay implementation.
