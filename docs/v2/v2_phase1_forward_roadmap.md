# ActiveCellSim v2 Phase 1 Forward Roadmap

> ⚠️ **CURRENT STATUS NOTE — 2026-05-06**
> (post-Option-A+B restructure audit cycle, claude-work
> `mcp_msg:2652`, codex approval `mcp_msg:2655`).
>
> This roadmap is a 2026-05-04 forward-looking design/ratification
> artifact, **not current implementation routing**.
>
> - PI directives after this roadmap superseded the "pending PI
>   ratification/update" status, including Option 1 pivot
>   (`mcp_msg:2187`) and V2-1 Q5 split decision (`mcp_msg:2235`).
> - The Closed-Loop ECM Gate 6-item blocker list below has since
>   been addressed through HB#1-5, Phase D no-op, Phase E v1,
>   Phase E v2 step 2, and Item 5 sweep work. Use the audit chain
>   and the current build-state doc for current status.
> - Closed-loop Tier 1 progressed faster than the original
>   projected timeline; older week-scale projections below are
>   **historical planning context**.
> - For current state, see
>   `docs/v2/v2_current_build_state.md` (refreshed 2026-05-06,
>   commit `da581cd`).
>
> Known caveats that remain current:
>
> - F-proper-1 is WIP HALTED (commit `8a24d89`).
> - Phase E v2 Item 5 sweep is surfaced as a B-tier sister
>   extension (commit `0f5e96f`).

Date: 2026-05-04 KST
Status: Claude/Codex design-discussion lock, pending PI ratification/update
Source discussion: `design-discussion`, topic `v2-layer-2-forward-roadmap`

## One-Line Recommendation

Phase 1 should remain a mechanics/ECM/protrusion/FA/junction/contact-inhibition
small-cluster reference engine. Planning headline: 7-8 weeks median, 6 weeks
stretch, 9-10 weeks conservative.

Executable cytokinesis and slow-biology dynamics should stay outside Phase 1.
Phase 1 keeps those fields schema-only / design-note only unless PI explicitly
overrides.

## Timeline

| Milestone | Stretch | Median | Conservative |
|---|---:|---:|---:|
| P0 alpha: data contract, metrics, schemas, frame dump, stub viz, tests | 24h | 24h | 36h |
| P1 active-contour gate attempt: cortex+area, Sanity Gate first | 24h | 72h | 1 week |
| Separated dynamics modes: active contour stable, ECM open-loop preflight, static FA, protrusion API force-disabled | 72h | 1 week | 2 weeks |
| Closed-loop Tier 1 single-cell: protrusion force, FA-protrusion coupling, FA to ECM, ECM bias | 2 weeks | 3 weeks | 4 weeks |
| First 2-cell contact + junction | 3 weeks | 4 weeks | 5 weeks |
| N=5-20 Tier 2 cluster + contact inhibition | 5 weeks | 6-7 weeks | 9 weeks |
| Phase 1 completion candidate | 6 weeks | 7-8 weeks | 9-10 weeks |
| Buffer / publication polish | - | - | 10+ weeks |

## Phase 1 In Scope

- Active-contour mechanics: cortex, area, protrusion force, FA traction.
- `ProtrusionEventScheduler` as one production API:
  - first gate: force output disabled for event/geometry/adaptive-vertex checks
  - second gate: force output enabled after active-contour timestep implications
    are derived.
- FA molecular clutch split into:
  - 6.3a static-boundary FA traction preflight
  - 6.3b protrusion-coupled FA maturation.
- ECM substrate dynamics:
  - open-loop preflight first
  - closed-loop only after the six-item gate below.
- Cell-cell junction and contact inhibition.
- HDF5 closure dump production format.
- Diagnostic visualization:
  - P0 stub
  - richer diagnostic 3D viewer around week 2
  - Blender diagnostic prototype around week 4
  - publication-polish rendering after Phase 1.
- Small clusters N=5/10/20 as stability permits.
- Cytokinesis design-note only, no executable implementation requirement.

## Phase 1 Out / Schema-Only

- Cell-cycle dynamics: schema-only.
- Cytokinesis execution: design-note only in Phase 1; executable work deferred
  to Phase 1.5 / Phase 2.
- YAP/TAZ feedback dynamics: schema-only optional field.
- O2/nutrient/necrosis/apoptosis dynamics.
- Full signaling pathways such as Notch/Wnt/TGF-beta.
- Publication-quality confocal Blender pass.

## Closed-Loop ECM Gate

Closed-loop ECM activation requires all six items before first execution:

1. Response monotonicity under prescribed traction.
2. Saturation behavior under repeated traction.
3. No response when traction is zero.
4. Bounded feedback in a single-cell loop.
5. Sensitivity sweep over grid spacing and `dt_ecm`.
6. Rule 10 unit-chain proof for traction-density storage and comparisons:
   `nN`, `um2`, `nN*s/um2`, and kPa / force-per-area comparisons before any
   ECM force or remodeling update.

Failure of any item blocks closed-loop ECM; open-loop ECM remains allowed.

## Top Risks

1. P1 active-contour timestep bound. If the Sanity Gate cannot derive a clean
   timestep from force and mobility, stop at P0 and report blocker/options.
2. Closed-loop ECM feedback. FA traction, ECM memory/alignment, and protrusion
   bias form a circular dependency that can saturate, oscillate, or hide
   arbitrary clipping.
3. N=20 cluster stability. If N=10 contact graph and junction maturation are
   unstable, do not force N=20.

## Process Note

This roadmap came from an adversarial design discussion under PI guidance
(`id=809` / `id=814`): Claude opened with an aggressive linear activation
roadmap, Codex challenged open-loop vs closed-loop ECM, dynamics-vs-schema
cycle-time extrapolation, and scope creep. Claude challenged missing
cytokinesis and protrusion staging. The final lock removed executable
slow-biology/cytokinesis from Phase 1 and changed the planning headline to a
7-8 week median.

The process worked better than the previous fast-convergence pattern because
it surfaced and corrected two late scope leaks before PI reporting:

- executable cytokinesis inside Phase 1
- conditional slow-biology activation inside Phase 1

Future design units should keep this pattern: adversarial opening positions,
explicit concessions, unresolved-disagreement list before lock, and a final
scope-leak pass.
