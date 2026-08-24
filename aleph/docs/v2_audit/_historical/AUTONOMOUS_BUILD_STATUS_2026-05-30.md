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

# Autonomous build status — 2026-05-30 (overnight /loop, PI asleep)

Durable handoff for PI wake. Branch `phase1/h3-cortex`, boot `a5afe2b`.

## Committed + independently verified (additive, non-governed only)

| Commit | What | Verification |
|---|---|---|
| `f36327b` | Mechanism audit + FA-integration design + PI-decision queue + production hardening (checkpoint/restart + integrity guard + KU-3.5 driver wiring) | regression 309/14/0; checkpoint fresh→RESUMED→recompute verified |
| `58b9f98` | Viz unlock: GSD trajectory writer + npz→GSD + fresnel renderer | render verified on-machine (journal-quality curvature-colored cortex, AO+colorbar+scale bar) |
| `0f655f6` | FA wiring S0–S2 into Cell.build (substrate anchor + integrin catch + fa_actin_clutch), additive default-off | 326/15/0, off-path bit-for-bit |
| `ec8a4e6` | PI-queue: A7 (S5 finding) + corrected per-occurrence citation plan | doc-only |
| `8ba3a00` | **S5 tag-space unification** — FA coexists with myosin/xlink/lamellipodium; unblocks KU-3.5 v4 build path | 328/15/0; off-path noop 5✓, H.4 standalone 18✓, KU-2.x catch-slip 26✓, integrator/ diff empty |
| (this commit) | Per-occurrence "Funk 2022" phantom-citation fix (27 occurrences, 6 files): capping→Li/Bieling 2022 eLife, abortive branching→Funk 2021 Nat Commun; δ_cap units typo pN→nm + PI-pending note. Doc/comment-only. | regression unchanged 328/15/0; yaml values (δ_cap 3e-10, abortive 500 Pa) byte-intact; AFINES file uncorrupted (897 lines, 0 embedded prefixes) |

**Known follow-up (left intentionally):** `AFINES_ALGORITHM_NOTES.md` has one residual DOI line `[Funk et al. 2022 eLife](https://doi.org/10.7554/eLife.73145)` — the DOI is the REAL Li/Bieling 2022 eLife paper, only the author-name string is wrong ("Funk"→"Li"). Trivial author-name fix, deferred (not attempted on a degraded tool channel late-session to avoid a bad edit).

## Notion (closeout protocol)
- Milestone child page created under Dev Logs (2026-05-30 H.4α strengthening pass).
- Dev Logs board: inline pointer prepended.
- RAG Index: reconciliation note appended (v1-era Unit-5 = collective/spheroid vs ffn_cellsim KU-5.x = lamellipodium); **KB restructure staged for PI, not auto-applied** (shared scientific contract + large-page corruption risk).

## What is DONE vs what needs PI

**Done autonomously** = everything fully verifiable without a gate-contract / frozen-integrator / magic-number decision: the audit, the viz pipeline, the FA wiring chain (S0–S2 + S5), production hardening, and the docs/Notion. KU-3.5 v4 *build path* is unblocked (FA+myosin coexist, short warm-up clean).

**Needs PI (the v4 production run + remaining co-requisites are gated on these):**
- **B1** global `dt = min(τ)` reconciliation (frozen-integrator).
- **B2** `equilibrate_no_shear` prelude in Cell.build (frozen-integrator) — also the step that brings the cortex onto the substrate at real geometry so clutch bonds actually form (today they don't at the physical 1.5 µm capture radius; the test uses a wide override).
- **A3** KU-5.1 density band ratification; **A4** KU-3.5 interphase-vs-metaphase regime.
- Co-requisites (membrane object S7, actin turnover, enclosed-volume pressure) are designed (FA design §4/§7) but deliberately NOT built tonight: they are larger new-physics modules whose target gates (KU-5.1, KU-3.5 v4) are PI-gated, so building them now yields unverifiable scaffold. Recommend building each as additive default-off once A3/A4/B1/B2 are decided.

See `PI_DECISION_QUEUE_2026-05-30.md` for the full A/B/C/D item list + recommended order.

## Honest notes
- Two false-alarm detours, both caught by verifying before acting: (1) briefly mis-suspected a Syncthing/.git race — it's a non-issue (`.stignore` excludes `.git`); (2) almost did a blanket "Funk 2022"→Li/Bieling replace — caught that it's a phantom mapping to TWO papers by context, fixed the plan to per-occurrence. No work lost.
- The persistent Bash shell's stdin got messy from early heredoc commit messages → switched to `git commit -F <file>`. Late-session the direct Read/Bash channel degraded (stale/abbreviated output), so the citation fixes were delegated to a fresh-channel subagent rather than hand-edited. All commits confirmed via clean fresh-shell `git log`.
