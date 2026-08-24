# `ragtagcag/` — does not import in this tree, and should not be made to yet

**880 lines, 8 modules. Six of them import `aleph.harness`, which does not exist here.** That is not
an oversight in the migration; it is the migration surfacing a decision.

## What it depends on

`Project_Aleph/aleph/harness/` — 3,831 lines of contract-graph, vault, object store, ledger and
snapshot machinery (`CAG`, `DataRoot`, `Ledger`, `ObjectStore`, `Snapshot`, `create_snapshot`,
`canonical_json`). It is that tree's own knowledge-base substrate.

## Why it was not brought with it

**Decision B4**: the knowledge base stays this tree's. Notion contract-graph as the single source of
truth, with a DuckDB + TAG query engine and an Obsidian projection materialised *from* it. Aleph's
`harness/` is a second answer to the same question, and the name says it — `ragtagcag` is RAG + TAG +
CAG, which is what `aleph/outputs/tag_kb/` already is here.

Bringing `harness/` to make these six files import would install two knowledge bases in one tree, and
then every later question about a source has two places to be answered from. That is the shape the
whole merge was against.

## What actually has to happen, and it is not a port

**Adjudicate `ragtagcag/` against `outputs/tag_kb/`, per module.** This tree's KB has a Notion
contract-graph with `source_audit.verdict` per source, a BM25 content layer, and code→contract edges
that come only from an explicit `Implements: KU-x.y` docstring line. Aleph's has a CAG with snapshots
and a canonical-JSON ledger. They are not the same design and neither is obviously the better one:

- **`tag_kb` is stronger on provenance** — every source carries an audit verdict, and a claim may not
  be cited in a deliverable until its verdict is OK.
- **`ragtagcag` is stronger on immutability** — content-addressed object store, snapshots, an
  append-only ledger with canonical JSON. `tag_kb` has no equivalent; it is refreshed, not versioned.

The honest outcome is probably neither wholesale. **It is `ALEPH-PORT-4002`+ and it is unwritten.**

## Until then

The modules stay, unimported and unmodified. They are not stubbed, not deleted, and not made to work
by pulling a second KB in behind them. `pyproject.toml`'s `testpaths` is `["aleph/tests"]`, so
`ragtagcag/tests/` is not collected and the broken imports cost nothing at rest —
`tests/architecture/test_outer_imports.py` pins that arrangement so it cannot rot into a surprise.

**Do not "fix" these six imports by adding `aleph/harness/`.** If that is the answer, it is the
answer to a port-ledger entry that argued it, not to an import error.
