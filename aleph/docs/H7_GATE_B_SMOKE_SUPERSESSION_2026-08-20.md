---
kb_record:
  topic: H7-gate-b-smoke-artifacts
  gate: VG-H7-gate-b
  status: authoritative
  authoritative_as_of: 2026-08-20
  aliases: [h7 gateB smoke, gate-b plumbing, gate-b emacheck, gate-b probe]
  supersedes:
    - RUN-h7-h7-gateb-relaxed-myooff.cond
    - RUN-h7-h7-gateb-rigid-myooff.cond
    - RUN-h7-h7-gateb-plumbing
    - RUN-h7-h7-gateb-plumbing2
    - RUN-h7-h7-gateb-emacheck
    - RUN-h7-h7-gbook-probe
  conclusion: >-
    Six H.7 Gate-B smoke RunResult rows are SUPERSEDED by the condition-matrix run set
    RUN-h7-h7-gateb-cm-*, and their artifacts survive in no branch. Retrieval must prefer the
    cm rows; the six are retained as record that the runs happened, NOT as citable evidence,
    because nothing can reproduce them. A seventh row, RUN-h7-gate-b-probe-seed1, is NOT
    superseded here — it has no successor and no git trace at all, and is flagged for PI review
    rather than swept in with the six.
---

# H.7 Gate-B smoke artifacts — superseded, and absent from every surviving branch

**Scope.** This record is about SIX RUN ARTIFACTS and nothing else. It makes no claim about
H.7 Gate-B itself, about `VG-H7-gate-b`'s verdict, or about any tension magnitude. The gate is
named in the front-matter only so the traversal can reach this note from the gate those rows hang
off — not because this note scores it.

## What was found (2026-08-20)

Reconciling the Notion `RunResult` table against disk turned up 19 rows with no corresponding
file. Seventeen of them were explained by path: the tree was renamed `ffn_sim/` -> `aleph/` on
2026-08-09 and reorganised again on 2026-08-20 (43 top-level entries -> 8), while the rows still
carry the old prefix. Those files are all present under `aleph/outputs/_archive/`.

**These six are different: the files are genuinely gone.**

| Row | was |
|---|---|
| `RUN-h7-h7-gateb-relaxed-myooff.cond` | `h7/production/h7_gateB_relaxed_myoOFF.cond.json` |
| `RUN-h7-h7-gateb-rigid-myooff.cond` | `h7/production/h7_gateB_rigid_myoOFF.cond.json` |
| `RUN-h7-h7-gateb-plumbing` | `h7/production/h7_gateB_plumbing.json` |
| `RUN-h7-h7-gateb-plumbing2` | `h7/production/h7_gateB_plumbing2.json` |
| `RUN-h7-h7-gateb-emacheck` | `h7/production/h7_gateB_emacheck.json` |
| `RUN-h7-h7-gbook-probe` | `h7/production/h7_gbook_probe.json` |

**Why they are gone, established rather than assumed.** All six were ADDED by `07a75b81`
("H.7 WIP cleanup: Gate-B smoke artifacts + SimuCell3D foundry-engine exploration", 2026-06-09).
That commit is **not an ancestor of `HEAD`**, and `git branch -a --contains 07a75b81` returns
nothing: it went with a branch deleted in the 45-branches-to-3 consolidation. So there is no
deletion commit to point at — the files were never on the line that survived.

This is not a wholesale directory loss. `aleph/outputs/_archive/h7/production/` holds 34 json
files, including the five `h7_gateB_cm_*` rows. Only these six are missing.

## Why the cm set is the successor

The naming carries it, and the survivors confirm it. `plumbing`, `plumbing2`, `emacheck` and
`gbook_probe` are scaffolding — plumbing checks, an EMA check, and a machine probe — run while the
condition matrix was being built. The two `.cond` rows have exact 1:1 successors that differ only
by the `_cm_` infix:

    h7_gateB_relaxed_myoOFF.cond.json  ->  h7_gateB_cm_relaxed_myoOFF.cond.json   (on disk)
    h7_gateB_rigid_myoOFF.cond.json    ->  h7_gateB_cm_rigid_myoOFF.cond.json     (on disk)

The matrix that replaced them is complete on disk and in Notion: `cm_relaxed_myoOFF`,
`cm_relaxed_myoON`, `cm_rigid_myoOFF`, `cm_rigid_myoON`, plus `cm_SUMMARY`.

## What was decided, and what was refused

**PI 2026-08-20 — option c + a: supersede AND label, archive nothing.**

- **(c)** the six are superseded by the cm set — this front-matter block, which is the only place
  `supersession.py` will read a supersession from (`kb_record:` front-matter, or an explicit
  `@supersession` marker in a DecisionLedger row; never inferred from dates or text).
- **(a)** each row's `Artifact Path` in Notion now reads
  `MISSING (branch lost, verified 2026-08-20) — was <original path>`, with the evidence above in
  its `Notes`.

**Archiving them was the option NOT taken, deliberately.** A row whose artifact cannot be produced
is a fact worth keeping: deleting it means the next person re-runs this same investigation from
scratch and finds the same nothing. What must not happen is the row being CITED as evidence, and
that is what the MISSING label and the supersession edge exist to prevent.

## ⚠ The seventh row is not part of this

`RUN-h7-gate-b-probe-seed1` is deliberately absent from `supersedes:` above. Unlike the six, it
appears **nowhere in git history** — not as an added file, not in any branch, and a content search
(`git log -S`) finds it only inside `tag_kb` Notion-snapshot commits, i.e. only as a string in the
KB's own snapshot of itself. Its siblings `gate_b_probe.json` and `gate_b_probe_seed2.json` both
exist. The likeliest explanation is a row harvested from a transient file that was never committed.

It has been labelled MISSING like the others but is **not** declared superseded, because there is
nothing to supersede it with, and asserting a successor that does not exist is the failure this
whole layer is built to prevent. **PI item: keep as a flagged unreproducible row, or archive it.**
