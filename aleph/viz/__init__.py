"""Rendering and figures. Above `laws/`, so it may look down at everything and is looked at by none.

Two sets of modules with different origins, deliberately not merged into one.

**`viz_*.py` — eight modules from this tree's `ff/`.** They were living *in* the law library, which is
the same defect `Project_Aleph`'s `runtime/law_cases.py` had: a thing that must see everything, filed
among the things it sees. Four of them import `archive.gamma_floor`, which is legal one layer up and
was a violation one layer down.

**Everything else — twenty modules from `Project_Aleph/aleph/viz/`**, carried 2026-08-09. A scene
graph, a renderer, frames, replay, a publisher, transport, channels, capture, cell and evidence
figures. The other session's porting plan flagged this as the package *"with a visible product — the
PI has been reading interactive cell figures out of it"*, and this session had it backwards until
that note landed: the first pass proposed keeping this tree's eight loose scripts and dropping the
system.

Three of the twenty do not import here
--------------------------------------
`boundary_figure.py`, `cpu_source.py` and `scene.py` import `aleph.runtime` —
`Transaction`, `Pipeline`, `participant`, `ForceWorkLedger`, `WarpBackend`. That is the layer this
tree **replaced** with `engine/` (decision A2), so there is nothing to point them at without either
carrying a second runtime or rewriting them against `engine/`'s contracts.

They are carried unmodified and unimported, like `outer/ragtagcag/`. Rewriting them is a port with a
ledger entry, not an import fix, and the adaptation is small but not mechanical: `engine/` has
`CellActor` and `CellWorldTransaction` where Aleph had `Pipeline` and `Transaction`, and the two
disagree about who owns the schedule. `tests/architecture/test_outer_imports.py` pins the quarantine at
exactly three so it cannot grow, and so nobody resolves it by importing a runtime that was removed on
purpose.

The seventeen that do import need no adaptation: they take arrays and write pictures.
"""
