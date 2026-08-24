# The chromatin/lamina fix had no evidence until the listing found it

**2026-08-21, session `d380041d`.** A small verification, recorded because it did not exist.

Earlier tonight two builders were changed to state their own topology in a form a consumer can read:

* `aleph/world/build/chromatin.py` — added `"n_strands": 1, "nodes_per_strand": int(n_sub)`, a
  restatement of the builder's own claim.
* `aleph/world/build/lamina.py` — added `"n_segments": 0`, `"topology_is_built": False` and a
  `topology_note` to `plan_lamina`'s return.

Both were committed with an argument and **no before/after measurement**. The exporter was believed
to have been dropping them; nothing showed that it had been, or that it had stopped.

## What `--list` found on its first run

`cell_app --list` exists to answer *"what can I open"*. Its first real output was:

```
ac/world_phase1/cell_current.alephcell   95 MB   12 populations   5.1 h
ac/world_phase1/cell_full.alephcell      94 MB   10 populations  11.7 h
```

Two files of the same cell, two populations apart. `--info` on each:

```
cell_full     (01:15, PRE-fix)   10 pops   chromatin absent
                                           lamina    absent
cell_current  (07:46, POST-fix)  12 pops   chromatin 552 nodes / 551 segments
                                           lamina    2043 nodes /   0 segments
```

**The two populations the older export lacks are exactly the two that were fixed.** They are back,
and the count of recovered chromatin segments — **551** — is the number the fix predicted.

⚠ **And `lamina` comes back with 0 segments, which is the point rather than a failure.** Its
topology genuinely is not built; the fix was to make the builder *say so* instead of returning a
census that a consumer read as "nothing to export". The app names it: *"lamina (2,043 nodes, NO
TOPOLOGY — not drawn)"*. **A population that is present and empty is different from one that is
absent, and before tonight the export could not tell them apart.**

## Why this is written down

The fix was argued and committed with no measurement of the thing it fixed. That is the shape this
session has spent the night cataloguing, one layer up: **a change believed to work, with nothing that
would have shown it did not.** The before/after existed on disk the whole time — two exports six
hours apart — and nothing compared them until a tool built for a different reason listed them side
by side.

⚠ It is also a caution about `cell_full.alephcell`: it is the **pre-fix** state and a reader opening
it sees a cell with no chromatin and no lamina. It is not deleted, because it is the only copy of the
"before", and it is named here so nobody reads it as a current cell.
