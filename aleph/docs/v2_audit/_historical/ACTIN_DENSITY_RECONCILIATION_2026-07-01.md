---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# Cortical actin filament count/density — conflict reconciliation

> **SUPERSEDING PI DECISION 2026-07-16 (new engine I0-A):** production cortex = **70,686 active F-actin
> filaments** at the MCF7 reference (`100 µm⁻² × 4π(7.5 µm)²`). The 38,000 and ×40 values below remain
> historical provenance only and are not production options. Non-cortical filaments are additional,
> compartment-specific density×geometry populations tracked by global unique IDs and the GPU-memory ledger;
> no global coarse-graining factor is applied to fit hardware.

**Date:** 2026-07-01  **Trigger:** the FF Stage-6Q lit-faithful-density work surfaced a conflict between
CLAUDE.md's "~38,000 native filaments" and KB-3.18's "~100/µm²" cortical actin density. This doc
reconciles the three project values, records their provenance, sets the FF working convention, and
isolates the one genuinely open (source-dependent) decision for the PI.

## The three conflicting values

| value | implied density @ MCF7 R=7.5µm (707µm²) | provenance | source-grounded? |
|---|---|---|---|
| **~38,000 native** (CLAUDE.md, Plan v2 §3 H.3) | ~54/µm² (or ~30/µm² @ R=10µm) | ratified ×40 coarse-graining ANCHOR (1000 effective = 38000/38) | modeling convention (plan), not a lit density measurement |
| **~100/µm²** (KB-3.18) | → N ≈ 70,686 | KB-3.18 "Cortex molecular composition (5-component)" | ✅ SOURCED — 18 SE edges, 9 citations (Salbreux2012_TCB, ChughPaluch2018_JCS, Charras2008_BJ, Ferrer2008_PNAS, Murrell2015, Gittes1993, …); see internal inconsistency below |
| **N~1.3e4** (KB-3.18, same row) | ~18/µm² | KB-3.18 companion figure | ⚠️ internally inconsistent with the SAME row's "100/µm²" (100×707 = 70,686 ≠ 13,000) |

Key facts established (git + KB audit, 2026-07-01):
- CLAUDE.md's "~38,000" is a **ratified coarse-graining anchor**, not an actin-density claim — it is the
  native count the ×40 mesoscale divides down to (~1000). It is used project-wide (GPU native pilots
  `--n-fil 38000`, `CORTICAL_TENSION_RECORD`, Plan v2). **Changing "38,000" retroactively edits a ratified
  Plan decision AND breaks the self-consistent ×40 ratio (38000/38 ≈ 1000).**
- KB-3.18 **IS sourced** — **18 `edges` to `source_evidence`** + 9 citations (Salbreux2012_TCB,
  ChughPaluch2018_JCS, Murrell2015_NRMCB, Charras2008_BJ, Ferrer2008_PNAS, Gittes1993_JCB, Isambert1995_JBC,
  Goldmann2002, PollardBorisy2003_Cell). ⚠️ Methodology note: the `edges` table keys on **Notion UUIDs**, so
  `edges WHERE src_id LIKE '%3.18%'` returns [] (a false negative on the label); querying by the row UUID
  `372120daec5d813795efc28ee677345f` returns the 18 edges. The real defect is **internal, not sourcing**:
  the same row's "100/µm²" and "N~1.3e4" are **mutually inconsistent** (100/µm² × 707µm² = 70,686 ≠ 13,000;
  13,000 ⇒ 18/µm²).

## Resolution

1. **FF working convention (in code, DONE, Stage 6Q):** the FF native cortex uses the KB-3.18 areal density
   **~100/µm² → N=70,686 at MCF7 R=7.5µm** (`cortex_assembly.CortexParams`, `gamma_floor.PROD_N_FIL`,
   `architecture_spec.CORTEX`). Rationale: an areal density is the physically transferable quantity (the
   count follows from density × the cell's actual area), and 100/µm² is the highest of the project's values
   so it is the least likely to under-count. This is the FF line's explicit choice, documented here.

2. **Historical note, superseded 2026-07-16:** CLAUDE.md formerly preserved "~38,000 native" as the ×40
   coarse-graining anchor. The PI has now performed the required plan revision: 38,000/×40 is historical-only,
   while 70,686 is the new-engine production cortical count.

3. **The γ-floor is UNAFFECTED either way.** The active cortical tension is myosin-bound (γ_myo channel);
   the actin filament COUNT does not drive it (FF_STAGE6Q: γ_myo is set by engaged NMII density × per-motor
   force, not the actin scaffold count). So this conflict is a citation/consistency matter, not a physics
   lever — it cannot rescue or deepen the floor.

## The ONE open item (PI + Notion SoT edit)

**KB-3.18's internal inconsistency ("100/µm²" vs "N~1.3e4") is a real KB bug that needs a Notion fix.** The
sources needed to settle it are ALREADY linked to KB-3.18 (Salbreux2012_TCB, ChughPaluch2018_JCS,
Charras2008_BJ, Ferrer2008_PNAS + 5 more) — this is NOT a "missing source" problem. Per the hard rule
("never hand-edit the vault/duckdb — fix Notion and refresh"), it is not hand-fixed here. Recommended action:

- Check which cortical actin areal density the already-linked primaries actually support (Salbreux2012 /
  Chugh-Paluch2018 / Charras2008), then correct the KB-3.18 row in Notion so the density and N are mutually
  consistent at the cell's real area (100/µm² × 707µm² ⇒ ~70,686, so the "N~1.3e4" companion is the figure
  in error if 100/µm² is retained). **PI-gated** (SoT edit).
- Until the Notion row is corrected, treat "N~1.3e4" as the erroneous companion and 100/µm² as the intended
  density (which is what the FF line uses). The floor is unaffected either way.

## Status

- ✅ Conflict documented + provenance established; FF working convention fixed (100/µm², in code).
- ✅ PI superseded the historical 38,000/×40 anchor; AGENTS/CLAUDE/master now bind 70,686 production cortex.
- ⏳ KB-3.18 Notion correction (resolve "100/µm²" ↔ "N~1.3e4" against its already-linked sources) = PI + SoT edit (flagged, not hand-fixed).

> Correction (2026-07-01, post adversarial review): an earlier draft wrongly stated KB-3.18 is "unsourced /
> no SE edge" — that was a `LIKE '%3.18%'` false negative (edges key on Notion UUIDs). KB-3.18 has 18 SE
> edges + 9 citations; the genuine defect is the row's INTERNAL 100/µm²↔N~1.3e4 inconsistency, corrected above.

Related: [[project-gamma-floor-likely-deficit]] (floor is myosin-bound, actin count irrelevant),
FF_STAGE6Q_LITFAITHFUL_DENSITY_2026-07-01.md, KB-3.18 (Notion SoT).
