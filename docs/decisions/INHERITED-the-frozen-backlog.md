# INHERITED — the decision backlog both trees brought, frozen

> **Status: AGENT-PROPOSED as an index.** Decision **D** of
> [`PROPOSAL-the-merge-and-the-rename.md`](PROPOSAL-the-merge-and-the-rename.md). This file names
> what is waiting and where it lives. It does not decide any of it, and **nothing here is a queue to
> pick work from.**

| Field | Value |
|---|---|
| Written | 2026-08-09, session `5ccfd28d` |
| Items | **47** — 35 Aleph proposals, 7 numbered Aleph decisions, 5 `STATE.md` (e) rows |
| Promoted | **4**, and they are one question |
| Bodies | The 35 + 7 stay in `Project_Aleph/docs/decisions/`; that tree is read-only and is not being demoted while `b4fad06b` is open |

---

## 0. Why one index and not forty-seven files

Copying 35 proposal documents into this tree would be inheriting a decision backlog wholesale, which
is the thing the merge spent the day *not* doing with code. The bodies are two directories away and
they are not going anywhere. What this file provides is the one thing an index gives that a
directory does not: **which of them are blocking, stated once, so a session does not read all
forty-seven to find out.**

**Most of them are questions belonging to a stage nobody has reached.** *Is the nuclear envelope
impenetrable*, *Darcy or Stokes for the cytosol*, *does the chromatin series need interior degrees of
freedom* — those are stage-3 physics. Answering them now would be answering before the thing that
makes the answer checkable exists.

## 1. Promoted — the four that block C, and they are one question wearing four hats

`STATE.md` (e) 1, `ALEPH-DQ-104`, and two Aleph proposals are all asking **when a step counts as
accepted**:

| id | where | what it asks |
|---|---|---|
| `STATE.md` (e) 1 | this tree | the held gate amendment — stationarity instead of convergence. PI HELD it 2026-07-28, and a defect surfaced after: stage 1's inner gate accepts on `balance_ok_d`, which `sf_motor_slice.py:65` says tests **adjoint wiring, not convergence** |
| `ALEPH-DQ-104` | `Project_Aleph` | the acceptance predicate — *"the one I deliberately left open"* |
| `PROPOSAL-what-does-converged-mean-for-a-cell-with-motors-in-it` | `Project_Aleph` | the same question, from the physics side |
| `PROPOSAL-the-engine-must-work-before-a-sweep-means-anything` | `Project_Aleph` | the same question, from the programme side |

**Decide them together or not at all.** Answering one without the others produces three records that
disagree about the same predicate, which is how `STATE.md` (c) got seventeen rows.

The prerequisite is not a decision, it is a measurement `ROADMAP.md` already specifies and costs one
GPU-day: **nobody knows what a converged native run costs, because no converged native run exists.**
172 s/step is a converging *attempt*.

## 2. Frozen — Aleph's 35 proposals

Bodies: `Project_Aleph/docs/decisions/PROPOSAL-<slug>.md`.

**Census and ontology** — `census-cannot-finalize` · `both-census-blockers-are-two-artefact-contradictions` ·
`the-census-counts-wiring-not-delivery` · `the-census-does-not-distinguish-cell-types` ·
`a-composite-group-is-one-connector-object` · `declared-parameter-multiplicity` ·
`focal-adhesion-endpoint-shape` · `the-far-field-anchor-is-a-clamp-and-the-registry-calls-it-a-connector` ·
`substrate-adhesion-belongs-to-an-owner` · `cortex-ecm-state-key-collision` ·
`cortex-is-a-filament-graph-not-the-envelope`

**Stage-3 physics** — `darcy-or-stokes-for-the-cytosol` · `the-nuclear-envelope-is-not-impenetrable` ·
`chromatin-series-needs-interior-degrees-of-freedom` · `solid-to-fluid-is-two-books-not-one` ·
`a-bond-has-no-formation-criterion-and-no-rupture-criterion` ·
`a-contact-is-not-a-bond-and-should-not-be-paired-once` · `the-cortex-has-no-sourced-excluded-volume` ·
`the-physiological-band-has-no-source` · `cortex-fluid-immersion-cost` · `only-the-cortex-has-a-resolution` ·
`steric-is-inert-and-the-real-index-cost-is-elsewhere` · `multirate-needs-a-registration-contract-not-a-cheaper-cortex` ·
`the-fixture-is-44-percent-placement-prestrain`

**Solver and numerics** — `the-descent-slack-is-derived-from-energy-not-from-a-length` ·
`position-precision-for-the-cuda-port`

**Governance** — `a-certificate-is-information-not-a-brake` · `caught-means-two-different-things` ·
`declaration-open-computation-sealed` · `resolve-the-seal-contradiction` ·
`p2b-lenses-are-not-aleph-lenses` · `plan-to-completion` · `two-session-scope-division`

⚠ **`two-session-scope-division` may already be answered differently here.** Sessions on this machine
can now message each other directly, which removes the constraint that proposal was written under.
`ownership.yaml`'s re-cut lanes are this tree's answer; reconcile rather than inherit.

## 3. Frozen — the numbered Aleph decisions

`ALEPH-DQ-101` initial scope · **`-102` first observation modality** · `-103` smallest latent target ·
`-104` acceptance predicate *(promoted, §1)* · `-105` repository and data root · `-107` production
runtime backend.

**`-105` and `-107` are answered by the merge and should be closed rather than carried.** The
repository is this tree, renamed; the production backend is Warp-CUDA on the shared workstation,
which decision B2 and the charter now state.

**`-102` is C-1 in disguise** — *"which observation does the sweep target"* is the same question as
*"which experiment is the engine finished against"*. Its standing candidate is the passive membrane
fluctuation spectrum. It is not promoted only because C-1 is the form this tree will decide it in.

## 4. Frozen — this tree's `STATE.md` (e)

1. the gate amendment — **promoted**, §1
2. the blebbistatin analogue needs a different knob; do not re-run the `k_xb` sweep
3. boot budget 46,000 B — **closed by this session.** The budget was not raised; `CLAUDE.md`'s dead
   GPU section (gbook, the lease, the queue) was removed and the boot context now sits at 45,688 B
4. `STRUCTURE.md` names directories that no longer exist — **closed.** It was rewritten from scratch
5. the master plan's IDENTITY claim — **closed by the rename.** That is what the rename *was*

## 5. What may be picked up from here

Nothing, except §1. Everything else is thawed by the stage that needs it, and a session that finds
itself blocked on a frozen item should say so rather than answer it — the freeze is load-bearing.
Three of these were written by a lane that had measured the problem, and re-deciding them from a
one-line summary in this index would lose exactly the evidence that made them worth writing.
