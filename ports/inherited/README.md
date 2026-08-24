# inherited — 135 port records from `Project_Aleph`, immutable

These are the ledger entries `Project_Aleph` wrote while porting **out of** `ffn_cellsim`, between
2026-07-30 and 2026-08-09. Decision **B5** carried them over with the mechanism, and

> ⚠ **135, not 136.** Every count of this set — B5's, and `PORTING_PLAN_FOR_THE_MERGE.md`'s own
> "136 entries" — came from `ls ports/ledger/*.md`, which counts `INDEX.md` as an entry. The index
> is not a port record. Corrected here and in the decision record.

`PORTING_PLAN_FOR_THE_MERGE.md` §7.2 asked whether they should travel at all.

They travel. They do **not** live in `ports/ledger/`, and the distinction is the point.

## Why they are not in the active ledger

`ports/ledger/` holds ports **into this tree**, and
`aleph/tests/ports/test_port_discipline.py` enforces that each one names controls that exist here,
cites a 40-character source sha, and declares exactly one status. Copied into that directory, these
135 fail on six counts at once — and every one of those failures is a **category error**, not a
defect in the entry:

| the check | why it fires, and why that is wrong here |
|---|---|
| named controls resolve to real tests | they name tests under `Project_Aleph/tests/`, which is a different tree. The tests exist; they are not here |
| section 2 records a source commit | these cite `ffn_cellsim` at revisions, in the format that tree's template asked for, not this one's |
| exactly one recognised status | Aleph's status vocabulary is close to this tree's and not identical |
| no placeholder tokens | some carry `TBD` in a section this tree's template does not have |

Holding a record of a port into another repository to this repository's port discipline would say
something false about the record. What these entries are is **evidence about why Aleph's contracts
exist**, which is exactly what `PORTING_PLAN_FOR_THE_MERGE.md` §4 calls them: the argument, kept so
the rules do not look arbitrary once the tree that produced them is an archive.

## The rules that do apply

1. **Immutable.** Do not edit an entry. If one is wrong, the correction goes in a new entry in
   `ports/ledger/` that cites it. Editing a record of a decision that was taken is not a correction,
   it is a different history.
2. **Not a precedent.** A rule stated here has whatever authority its own PI record gives it, and no
   more from having been carried. Several are `AGENT-PROPOSED` and unratified — including the
   representation and MTG-PN specifications that `aleph/outer/` refers to.
3. **Read them before re-deriving.** `-3657` through `-3665` in particular carry measurements about
   cells rather than about Aleph — sourced filament-length distributions, the microtubule count as a
   length density, the NF2007 semi-implicit step at 708x, the IMEX driver at 28.3x whole-cell. Those
   are facts this tree can use, and `docs/decisions/PROPOSAL-the-merge-and-the-rename.md` §1b acts on
   several of them.

## Numbering

Inherited IDs run up to `ALEPH-PORT-3665` and are not renumbered. **Ports in the new direction —
`Project_Aleph` into this tree — start at `ALEPH-PORT-4001`**, so the ID alone says which way the
code moved. The gap between the two bands is deliberate.

`INDEX.md` is Aleph's generated index of these entries, carried as-is. It indexes this directory, not
`ports/ledger/`.
