# ffn_cellsim — agent charter pointer

**This file is NOT a charter. `/CLAUDE.md` is the single project charter — read it fully at minute 0
of every session (Codex, Claude Code, or any other agent) and follow it verbatim.** This file exists
only so an agent that boots from `AGENTS.md` is redirected there instead of acting on a stale copy.

Do not re-expand this file into a parallel charter — two charters drift, and the drift is silent. The
247-line 2026-07-16 body (archived verbatim at `aleph/docs/v2_audit/_historical/AGENTS_2026-07-16.md`)
still described `ac/` as a pre-canonical scaffold and still carried the **retired** unified-network double-count
framing nine days after PI reframed it (audit R3, `aleph/docs/v2_audit/AUDIT_WHOLE_REPO_2026-07-25.md:79`).
If a rule changes, change it in `CLAUDE.md`.

## The four things a stale boot gets wrong (all are HARD rules; full text in `CLAUDE.md`)

1. **`ac/engine/` is the CANONICAL composition layer** (13 components / 32 connectors; SoT
   `aleph/docs/v2_audit/cell_engine/`) — not a scaffold. `ac/cell/` is the running incumbent being
   bound under it, and is **feature-frozen** (PI 2026-07-22): instrumentation / parameter guards /
   P0-blocker fixes only, no new biology.
2. **Disjoint populations, PI-reframed 2026-07-22 (SF-separate).** Cortex, SF/arc, lamellipodium and
   filopodium are **separate state-owning components, each owning a DISJOINT filament population**;
   every physical filament belongs to exactly **ONE** component. There is no unified label-blind
   network, so "don't double-count one filament under multiple labels" is the wrong mental model.
3. **Co-location is NEVER a connection.** Sharing an array is not a weld. Connectors are the only
   mechanical connections, and each is bidirectional + adjoint (Newton's 3rd law); kinetic connectors
   commit only on an accepted physical step, under one device-resident accepted-step transaction.
4. **Warp CUDA GPU is the only simulation runtime.** HOOMD is never imported or executed — not for
   production, dev, parity, fallback, benchmarks, or validation. The retired HOOMD tree was
   DELETED 2026-07-29 and survives only in git history (`1a9ded66^`). Validate and debug at the **full native population + compartments**
   at physiological setpoints; coarse/stripped configurations cannot support a conclusion.

## Codex-session operational note (still current)

Codex output is **advisory** — debugging ideas and independent review, **not verification authority**
(`aleph/docs/v2_audit/SUCCESSOR_HANDOFF_2026-07-25b.md:135`). A finding is adjudicated by measurement
on the native run, not accepted because Codex said it. PI-gated items are never changed autonomously:
gate contracts and thresholds, magic-number triggers, parameter sourcing, mechanistic-model decisions,
and the frozen `ac/cell` contract. Write the precise proposal and surface it to the PI instead.

Repo map: `/STRUCTURE.md`. Current `ac/engine` status: `aleph/docs/v2_audit/cell_engine/ROLLING_ROADMAP.md`.
