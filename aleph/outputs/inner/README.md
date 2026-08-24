# `outputs/inner/` — empty on purpose

The inner network's artifacts land here when it opens. It has not opened, so this directory holds
nothing but this file, and that is the correct state rather than an oversight.

`aleph/inner/__init__.py` carries the reason in full. In one line: the inner surrogate is **Stage 2
and blocked on C-1** — until the one experiment the sweep targets is declared (name, protocol,
observable, uncertainty) a surrogate has no output to emulate, and opening the package early
"produces a model of an unnamed quantity."

**Why the directory exists at all.** `outputs/` is now split by layer — `ac/` for the engine,
`inner/` here, `outer/` for the external-research layer, `_archive/` for superseded lines. An absent
`inner/` reads as *forgotten*; an empty one that says why reads as *not yet*. The distinction is the
whole point of the directory.

Candidate material waiting for triage is in `aleph/virtual_cell/` — `surrogate.py`, `sweep.py`,
`reduction.py`, `tensor_train.py`, all tested, none wired to the engine.
