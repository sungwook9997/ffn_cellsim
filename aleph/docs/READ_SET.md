# READ SET — what a fresh session should read, and why

Moved out of `STATE.md` on 2026-07-28: the list of paths has to stay in `STATE.md` (its validator
asserts every one exists on disk), but the paragraph explaining each did not, and the paragraphs were
the expensive part. `STATE.md` keeps the ordered paths and one clause each; the reasoning is here.

**⚠ This file is status.** Which documents are worth reading changes as documents are archived,
superseded or found to mislead. If an entry here names a file that no longer exists, this file is the
defect — `STATE.md`'s list is the one under test.

---

## (d) Minimal correct read set for a fresh session

<!-- STATE-READSET:BEGIN -->

1. **`ROADMAP.md` — READ THIS FIRST, it is where the stage gates now live.** Moved out of `CLAUDE.md` on 2026-07-28 because it was 46% of the charter's boot cost and was current state wearing the clothes of a rule. It states what each stage must PROVE and never whether it has been proved — that is (b) below, which wins if the two disagree. Its own header lists which of its contents change and what to verify each against; note in particular that **the gate wording is under a PROPOSED, unratified amendment**.
2. **`aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md` — the per-stage working detail.** The PI-directed reframe of 2026-07-25: purpose becomes inferring per-cell-type / per-cell-state molecular parameters (parameters are OUTPUTS), the gate becomes mechanical CONNECTEDNESS rather than magnitude, and the baseline becomes a cell in media adherent on 2D collagen. It carries the eight open PI decisions, the three measurements that open six tracks, the data-representation rule (fields/graphs/modal coefficients, never summary scalars), and five verifier-confirmed code defects. Track-level detail is its companion `aleph/docs/v2_audit/COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md` (145 KB — read the track you are working, not the whole file); raw audit output for all of it is `aleph/docs/v2_audit/_raw_2026-07-25/` (4.0 MB, cite rather than re-run).
3. `CLAUDE.md` — the charter, and it is now one reframe behind this file. Read past three known defects: the retired-then-remandated handoff/4-role protocol, the stale KB scale figures, and the H.1→H.7 program it still names as active (h7/h8/h9 have no closeout and the chain has had no work since 2026-06-11).
4. `aleph/docs/v2_audit/cell_engine/ROLLING_ROADMAP.md` — the accurate per-component status, and the only document that says the GATE A/B cortex is coarse-grained.
5. `aleph/docs/v2_audit/SUCCESSOR_HANDOFF_2026-07-25b.md` (later, 16:35) and `aleph/docs/v2_audit/SUCCESSOR_HANDOFF_TO_OPUS5_2026-07-25.md` (12:01) — read the later, scan the earlier for anything not carried forward; two same-day handoffs is itself the finding.
6. `aleph/docs/v2_audit/cell_engine/CELL_ENGINE_ARCHITECTURE.md` + `aleph/docs/v2_audit/cell_engine/COMPONENT_BUILD_MATRIX.md` — the component/connector contract and the 8-state ladder (the matrix is staler than the roadmap; roadmap wins).
7. `aleph/docs/v2_audit/PI_GAP_EVIDENCE_CARDS_2026-07-25.md` + `aleph/components/motor/params_i0b3.yaml` — what is blocked on PI; the YAML is the reference implementation of a recorded parameter gap.
8. `aleph/docs/v2_audit/CORTEX_STRUCTURE_AUDIT_2026-07-23.md` + `aleph/docs/v2_audit/CORTEX_MESH_FIDELITY_2026-07-24.md` — the two documents that decide which older results are still citable. Without these a fresh reader re-quotes the γ-floor.
9. `aleph/docs/v2_audit/PARAM_PROVENANCE_AUDIT_2026-07-24.md` — the SOURCED / DERIVED / CONVENIENCE / PI-GAP classification, including turgor Π₀ = 40 Pa named as the worst convenience default.
10. `aleph/outputs/tag_kb/run_audit_report.md` + `aleph/outputs/tag_kb/param_audit_report.md` + `aleph/outputs/tag_kb/source_audit_report.md` — the three disk-grounded verdicts, read with (b)'s caveats about what the manifests still omit.

<!-- STATE-READSET:END -->

~288 KB measured 2026-07-25, of which `source_audit_report.md` is 91 KB (skim it). Most likely to mislead, in order:
`PHASE_0_CLOSEOUT.md` + `PHASE_0_3_DECISIONS.md` (CLAUDE.md's §Active roadmap now names both historical, but they
still read as current), the 12 HOOMD-era H briefs, `NEW_ENGINE_BUILD_PLAN_2026-07-16.md` (64 KB, invites booting from
it, superseded in six days), `FF_ENGINE_REVALIDATION_LADDER_2026-07-16.md`, `FF_STAGE6O_NATIVE_GAMMA_2026-07-01.md`,
`README.md` + gh-pages, `results_manifest.yaml` — a green gate reads as coverage. **Archived 2026-07-28: `v2_audit` root 261 → 52** — the 209 that pre-date the reframe AND no live file cites now sit in `_historical/` with `superseded_by` front-matter saying why.

