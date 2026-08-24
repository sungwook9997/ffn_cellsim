# The three layers, and what actually separates them

Written 2026-08-20, session `d380041d`, at PI instruction *"엔진/내부신경망 데이터/외부신경망 데이터도
구분이 안되어있고"*.

`CLAUDE.md` names three layers: the mechanistic Warp/CUDA instrument, **a network learned from its own
accepted runs (`inner/`)**, and **one learned from external research (`outer/`)**. The names existed. What
did not exist was any place a session could look and see which of the three a given directory belonged to
— `aleph/outputs/` in particular held 43 top-level entries across four eras with nothing marking the split.

> **This file is not `STRUCTURE.md`.** That one is PI-authored and read-only during unit work
> (`ownership.yaml:166`), and it still describes the pre-2026-08-20 `outputs/` tree. Where the two
> disagree about `outputs/`, this file is the newer statement and `STRUCTURE.md` needs a PI-approved
> refresh — listed at the bottom.

## The split, as it now stands on disk

| layer | source of truth it learns from | code | artifacts |
|---|---|---|---|
| **engine** | nothing — it *is* the mechanism | `laws/` `engine/` `components/` `world/` `viz/` `observe/` `units/` `validation/` | `outputs/ac/` |
| **inner** | this engine's own ACCEPTED runs | `inner/` (empty — Stage 2, blocked on C-1) | `outputs/inner/` (empty, with a README saying why) |
| **outer** | external literature and experiment | `outer/` | `outputs/outer/` |
| — | superseded lines, kept not deleted | `archive/` `dcm/` | `outputs/_archive/` |

Live KB infrastructure — `outputs/tag_kb/`, `outputs/obsidian_rag_full/`, `references/` — is outer-layer
material by content and **was deliberately left where it is**: 80+ files reference those paths and
`make kb-check` reads them. Moving working machinery to make a listing tidier is the wrong trade.

## What enforces it — a test, not this document

`aleph/tests/architecture/test_layer_directions.py` holds the import graph's direction in its `LAYERS`
table. Two rows carry the layer separation:

    outer  MUST NOT import  engine, components, laws, dcm, inner, ...
    inner  MUST NOT import  components, outer, ...

**`inner` was added on 2026-08-20 while the package still holds nothing but `__init__.py`**, which is the
cheap moment to do it. That test's own docstring records why: the cycle it was written for came from
"four good decisions", each defensible in its own diff, and a counter saw fifty-three imports where a
reviewer saw one. A direction declared before a package opens costs one line.

`laws` and `engine` are deliberately **not** forbidden to `inner`: a surrogate must be able to name the
observables it emulates and the units they carry, and forbidding that forces a copy. What it may not do
is reach into a specific component's state, or into the other network.

## Why `inner/` being empty is a state and not a gap

`aleph/inner/__init__.py` is the authority. Short version: the surrogate is Stage 2 and blocked on **C-1**
— until the experiment the sweep targets is declared (name, protocol, observable, uncertainty) there is
nothing for a surrogate to emulate, and opening it early "produces a model of an unnamed quantity."
Candidate material sits in `aleph/virtual_cell/` (`surrogate.py`, `sweep.py`, `reduction.py`,
`tensor_train.py` — all tested, none wired), to be triaged into `inner/` when it opens.

## ⚠ Open to the PI

1. **`STRUCTURE.md` is stale on `outputs/`** and may not be edited without PI sign-off. Its `outputs/`
   rows describe the flat 43-entry tree that no longer exists.
2. **`.gitignore:160** (`aleph/outputs/**/*.npz`) makes every raw run dump unversioned by policy, while
   the charter says a run's data is committed beside its record. It cost the GATE-B row its only artifact
   until 2026-08-20. The rule was NOT changed here — that is a PI call.
