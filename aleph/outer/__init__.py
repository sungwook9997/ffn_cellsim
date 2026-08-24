"""The outer network: learned from external research. **Stage 4. Open, and no longer empty.**

53,602 lines of pipeline arrived 2026-08-09 from `Project_Aleph`'s `codex/external-training-corpus`
branch at `803ae3ec` — acquisition, extraction, the experiment factory, retrieval, verification,
review, orchestration and a model directory. **This package was a docstring pointing at another
repository for six hours before anyone noticed**, which is recorded in the merge record because the
thing that hid it was a session registry row nobody opened.

Layout
------
`acquisition/` and `ACQUISITION_POLICY.md` · `extraction/` · `experiment_factory/` (196 modules, the
bulk) · `retrieval/` · `review/`, `review_oa/` · `verification/` · `orchestration/` · `model/` ·
`queues/`, `snapshots/`, `sources/`, `tools/` · `QUALITY_POLICY.md` · `manifest.schema.json`.

**`ragtagcag/` does not import here and that is deliberate** — six of its modules want
`aleph.harness`, `Project_Aleph`'s own KB substrate, and decision B4 keeps this tree's. See its
README; the quarantine is pinned by `tests/architecture/test_outer_imports.py`.

**The 280 MB of built tensors are not in git.** They are outputs, gitignored under
`aleph/outputs/outer/`, with `PROVENANCE.md` recording what they are and how to rebuild them. About
12.7 GB of raw acquisition data was never in git anywhere and still lives in the source worktree —
⚠ which means `Project_Aleph` cannot be deleted until that data has somewhere else to be.

Why this one is not blocked
---------------------------
Its inputs are literature, omics and experiment, not simulation, and the knowledge base those come
from already exists in this tree (Notion contract-graph + DuckDB TAG + the Obsidian projection). So
stage 4 runs in parallel with stage 1 and never touches the engine — which is also why it is its own
ownership lane.

Three layers, per `docs/design/INHERITED-MULTIMODAL_OUTER_LIBRARY_MTGPN_REPORT.md`
-----------------------------------------------------------------------------
1. **Experimental Observation Compiler** — IF, WB, PCR, PIV, TFM, AFM, FRET, FRAP, microscopy and
   direct force/tension data into one observation language carrying units, coordinates, uncertainty
   and protocol. `SOURCE_UNIT_UNKNOWN` is a value it must be able to record: a published series that
   arrived without a unit is a fact about the source, not a gap to fill in.
2. **MTG-PN** — seven stages, not three. H1-H3 pool, G relates, L brings measurements in, P answers,
   V says how much to believe the answer and when to refuse. Conflating the last four with the
   hierarchy is how "the encoder" ends up meaning both "pools segments into filaments" and "decides
   whether to refuse".
3. The **native simulator**, which the learner does not replace. Out of distribution, it runs.

What travels here as a document, not as code
--------------------------------------------
`learn/` and `represent/` (12,400 lines, 15 files, zero numpy/torch imports — verified) are
**computation-sealed by ratification**: no array import, no float literal, every compute entry point
raises. There is no code to port, only a specification, and decision A4 carries it across as one —
`docs/design/INHERITED-sealed-declaration-layer.md`, every module docstring and public symbol, with
the arithmetic that does not exist left out because it does not exist.

The specification is `AGENT-PROPOSED`. That means it cannot be used to **reject** an implementation;
it can be used to **check** one. Where this tree already has something on the same axis, this tree's
is what runs — §1a of the merge record, where `virtual_cell/observation_operator.py` turned out to be
the stronger of the two.

`aleph/virtual_cell/` holds this tree's own attempt — `observation_operator.py` (1,285 lines, 884 of
tests) is **stronger** than Aleph's on the axis they share, and the one idea it lacks is the
apparent-vs-true split with `analysis_chain`. That narrow port is `ALEPH-PORT-4002`, unwritten.
"""
