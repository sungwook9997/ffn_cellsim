# V2 Phase 1 Project alpha — P0 Morning Summary

**Generated**: 2026-05-04 ~01:38 KST.
**Audience**: PI ↔ design-discussion handoff.
**Status**: P0 (Plan §13.1) **complete and committed**. P1 not started;
derivation document required first per Plan §13.2 / Codex MCP id=896 / 899.

---

## 1. P0 deliverables (committed)

All commits live on `main` between `5568a21` (pre-α backup tag) and current
HEAD `f04e194`. Backup tag `backup/pre-project-alpha-2026-05-04T00-30Z`
captures the pre-overnight state. Tag `v2-alpha-checkpoint-1` (set by Codex
gamma backup) → `75824af` (end of Cycle B).

| # | Commit | Cycle | Headline | Files | LOC | Tests |
|---|---|---|---|---|---|---|
| 1 | `64332b4` | P-1 | `requirements-dev.txt` for `.venv-collab` | 1 | +42 | smoke-OK |
| 2 | `e61f3a3` | A | data_contract.py — 9-member `ArtifactKind` (incl. `ECM_FIELD`); artifact-aware `MetricSpec`; duplicate key `(name, modality)` | 3 | +437/-54 | 16 |
| 3 | `75824af` | B | `metrics.py` registry + `SingleCellState` slow-biology hooks | 5 | +519/-2 | +14+15 |
| 4 | `bab037b` | C | 5 schema-only validators: `ECMSubstrateState`, canonical `ProtrusionEvent` (Plan §5.4 fields), canonical `FocalAdhesionState` (xy 2D, no-traction-without-attachment), `JunctionState` (canonical pair), `CellClusterState` (symmetric neighbor graph) | 13 | +1432/-71 | 66 |
| 5 | `8b00f57` | D | HDF5 `frame_dump` writer/reader with **full schema fidelity**: every optional field round-trips, including slow-biology hooks, neighbor graph, role labels, lineage refs, full protrusion / FA / junction metadata. Atomic write via `mkstemp` + `os.replace`. | 4 | +1218 | 13 |
| 6 | `f04e194` | E | Stub 3D viz: matplotlib Agg PNG + self-contained HTML/SVG. Headless-SSH friendly. **Confocal-quality Blender pipeline is a separate design unit**, not in alpha. | 4 | +395 | 5 |

Codex side: `4808848` plan-docs commit and tag `v2-alpha-checkpoint-1` were
created during the gamma-backup window.

**Total v2 test count**: **217 passed, 2 expected Taichi skips**
(`tests/test_physics_conservation.py`, `tests/test_taichi_gpu.py` — both
require Taichi which is intentionally not in `.venv-collab`).

## 2. Plan §15 morning-acceptance checklist

| Criterion | Status | Note |
|---|---|---|
| Consolidated plan present | ✅ | `docs/v2_phase1_plan_consolidated.md` (committed in `4808848`) |
| P0 schemas / contracts / frame-dump / viz drafted or implemented | ✅ | Implemented + tested |
| No PI experimental data used | ✅ | Hard Rule. No reads from `data/experimental/*` introduced. |
| No hidden numeric defaults | ✅ | All optional/sentinel constants are named module-level (`_LIFETIME_CONSISTENCY_TOL_S`, `_ORIENTATION_SYMMETRY_TOL`, `_ORIENTATION_BOUND`, `_DEFAULT_FIGSIZE_INCHES`, `_HEIGHT_TINT_SCALE_UM`). No physics tunables in P0; the listed constants are validation/render numerical tolerances, named and documented, not fitted. |
| No gate tolerance edits | ✅ | None |
| Blocked P1 reported, not forced through | ⏳ | P1 not yet attempted — see §3 below |

## 3. P1 status — design-team input requested

P1 = `ActiveContourState` + cortex / area-restoring numerics
(Plan §13.2 / §6.1). Plan §11 and §13.2 require a Sanity Gate that
derives:

1. unit chain for vertices (μm), force (nN), mobility (μm·s/nN), time (s)
2. timestep bound from implemented stiffness and mobility
   (`dt < c · 1/max(μ·K_eff)`)
3. boundary cases (N<3, NaN, self-intersection, tiny edges)
4. area drift and symmetry preservation
5. force sign/sense check (cortex contractile, area restoring)
6. measurement protocol for the exported `MeasurementBoundary`

before any executable physics is written. Plan §11 already labels
`K_A`, `substrate drag coefficient`, `dt_cell` as
`requires Sanity Gate derivation before execution`. Until those are
literature-anchored or first-principles-derived, no executable physics
should land in the repo per CLAUDE.md / AGENTS.md Hard Rules.

**Open questions for design-discussion** (the questions PI is taking
back to design):

A. **Cortex / line tension convention**:
   - **Need from design**: choose a cortex convention with explicit
     units, anchored to IF≥15 literature (Chugh et al. 2017 Nat Cell
     Biol, Salbreux–Charras–Paluch 2012 review, or another anchor
     they prefer). Candidate conventions to choose from:
       (i) **line-tension energy**: `E = γ · L` with γ in nN
       (force per unit boundary length integrated to give a line
       tension), so the per-vertex force is the variational derivative
       `F_v = γ · κ · n̂ · ℓ_v` where `ℓ_v` is the vertex control
       length (sum of half adjacent edges). γ ranges Chugh-anchored.
       (ii) **2.5D cortex sheet**: `E = γ_s · area_cortex` with
       γ_s in nN/μm; per-vertex force has an additional integration
       step against the cell's height profile.
       (iii) other convention they propose.
   - The Sanity Gate derivation **must show** the per-vertex force
     comes out in nN, then `v_v = μ · F_v` is in μm/s, then `dt`
     follows. Plan §6.1 only sets the qualitative direction
     ("shape-restoring, contractile"), not the convention.

B. **Area-restoring convention**:
   - **Need from design**: choose how the area term enters the
     vertex equation, anchored to literature. The three usual
     conventions are:
       (i) **2.5D pressure** (cell with fixed height H): area energy
       `E = ½ K_A · (A - A_target)² / A_target` with K_A in nN·μm
       (so `dE/dA` is in nN/μm and `dE/dx_v = K_A · ΔA · ∂A/∂x_v`
       gives nN per vertex once `∂A/∂x_v` carries μm).
       (ii) **line-pressure**: area term is integrated as a uniform
       outward pressure per boundary segment, with K_A in nN/μm.
       (iii) **bulk pressure** (3D-like): K_A in nN/μm² × volume,
       requires height H to be a state variable.
   - The chosen convention **must** include an explicit derivation
     showing `force_per_vertex [nN]` and the segment- or
     control-length integration that gets there. Plan §11 labels
     K_A as "requires Sanity Gate derivation before execution"
     precisely because the convention choice changes the unit chain.

C. **Mobility / substrate drag**:
   - Plan §11 labels substrate drag coefficient
     `requires Sanity Gate derivation before execution`. Question for
     design: is Phase 1 a uniform overdamped boundary-vertex
     `μ_uniform · v = F_total`, or a viscous-like drag against the
     substrate that depends on FA distribution?

D. **Timestep**:
   - Plan §6.1 forbids fixed `dt = 10 ms` without derivation.
   - Once γ_c, K_A, μ are picked above, `dt` follows from the explicit
     Euler stability bound. Design should confirm whether implicit
     stepping is on the table for alpha (would relax the bound but
     adds a linear-solve dependency; CLAUDE.md "halt and surface to PI"
     would apply).

E. **Vertex-density adaptive scheme**:
   - Plan §4.1 requires baseline boundary spacing + refinement around
     active filopodia. P1 alpha does not include filopodia dynamics
     (that is staged after FA/ECM). For the alpha cortex/area solver,
     is a **fixed vertex spacing** acceptable (refinement disabled)
     until the filopodia module lands? This dramatically simplifies
     the gate.

F. **Boundary / substrate contact**:
   - Free boundary (substrate flat, no walls) for the alpha symmetric
     polygon relaxation test? Plan §13.2 asks for "symmetric polygon
     relaxation". Confirm.

If design provides anchors for A–F, the Sanity Gate derivation can be
written without ambiguity. If any of A–F stays unresolved, Codex's
recommended path (id=896) is to surface as `status=blocker` to PI
with at least three concrete options rather than write executable
physics under uncertainty.

## 4. Recommended next steps

1. PI takes §3 questions A–F to design-discussion.
2. design-discussion returns either (a) a closed answer to each, or
   (b) deferral notes with explicit "not for alpha" labels.
3. Implementation-work then writes
   `docs/v2_p1_active_contour_sanity_gate.md` covering the six items.
   Codex review re-verifies units (`μ·K_eff` must reduce to 1/s),
   per-vertex force from `K_A`, and grid-invariance of `dt`
   (id=899).
4. If gate PASSES → `acs/v2/active_contour.py` (state) +
   `acs/v2/dynamics/active_contour.py` (cortex + area) lands in two
   small commits with the Sanity Gate docstring already in place.
5. If gate FAILS → blocker message to PI with three options
   (Plan-compliant). No partial physics committed.

## 5. Operational notes

- Overnight auto-commit policy in effect (PI mcp id=790). High-risk
  items (Sanity Gate FAIL, magic-number exception, gate tolerance
  change, production/sweep, destructive op, experimental data risk)
  still go to PI approval queue.
- ssh win RTX A5000 ≤ 3 hours: PI approval-exempt (PI mcp id=845).
  Currently unused; readiness only.
- Live tree dirty files outside α scope: none in `acs/v2/` or
  `tests/`. The collab tooling files
  (`tools/collab_mcp/{room.py, setup_tmux_workroom.py}`) and
  `tests/test_collab_room_upload.py` were dirty at various points and
  are managed separately.
