"""The inner network: learned from this engine's own accepted runs. **Stage 2. Not open yet.**

What goes here
--------------
The surrogate and its sequential design — the thing that decides which native run to spend a GPU-day
on next, and the emulator that stands in for the ones not worth spending it on. Downstream of that,
the **Inner Library**: accepted native states and the observables derived from them, admitted only
when the run converged *and* carries an accepted digest.

Why it is not "after the sweep"
-------------------------------
It reads as a post-processing stage and it is not. `ROADMAP.md`'s arithmetic: an 8 x 5 grid is 14.9
years before the x4 for cell states, so **there is no sweep for a surrogate to be downstream of**.
The surrogate is one of the four multipliers that make the sweep exist at all, beside raising `dt`
under a physical predicate, cutting dimension by Fisher, and replacing the grid with the analytic
gradient this tree already validated. Together: ~390,000 runs -> ~100-300.

What blocks it
--------------
**C-1.** Until the one experiment the sweep targets is declared — name, protocol, observable,
uncertainty — a surrogate has no output to emulate and Fisher has no observation to be informative
about. Opening this package before that produces a model of an unnamed quantity.

⚠ And a warning about the arithmetic above: `ROADMAP.md`'s budget figures that rest on an
extrapolated exponent are **re-derivable, not citable** — a `steps^2.08` extrapolation in a planning
document measured **0.97** on the GPU, and `ALEPH-PORT-3665`'s IMEX driver already returns 28.3x
whole-cell at 10x `dt`. Re-derive before sizing anything on them.

`aleph/virtual_cell/` holds candidate material — `surrogate.py`, `sweep.py`, `reduction.py`,
`tensor_train.py`, all tested, none wired to the engine. Triage them here when this opens, per
decision A4 as revised.

The inherited side of the specification is `docs/design/INHERITED-sealed-declaration-layer.md`
(`represent/` — blocks, cuts, families, index, plan, qualification) and
`docs/design/INHERITED-ALEPH-TN-LANGUAGE-SPEC.md` for the tensor-network language it is written in.
Both are `AGENT-PROPOSED`: a check, never a veto.
"""
