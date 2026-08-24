# The run host's commit stamp was 101 commits stale, and nothing in any artifact could have said so

**2026-08-21, ~04:50 KST.** Found while trying to give two new `STATE.md` rows the build commit
`check_state_md.py` demands. Both rows said `2026-08-21 world/` and named no hash. The checker was
right to refuse them, and the reason they had no hash is the finding.

## 1. What was measured

The GPU host `sungwook@100.110.26.26` holds `~/ffn/world` as a **deployed tree with no `.git`**. Its
`COMMIT_STAMP` file read `b9de49cd` — *"feat(gpu,world): a warm executor inside one hold"*, committed
**2026-08-20 18:00:02**, the moment Slurm job 72 started.

By 04:50 the local branch was at `bb4cc80f`. Between them: **101 commits.**

Every native number produced tonight after the stamp was written — the balance-gate series, γ on the
resting path, the thermostat runs, the τ windows, the seven-population step cost — ran on code that
the stamp does not describe.

## 2. The stamp was not merely old — the tree was mixed

Hashing the host's 51 `aleph/world/*.py` files against candidate commits:

| commit | files identical |
|---|---:|
| `HEAD` (`bb4cc80f`) | 41 / 51 |
| `3f29f638` | 41 / 51 |
| `a1bdcd4d` | 40 / 51 |
| `b9de49cd` **(what the stamp claimed)** | **5 / 51** |

**No commit describes the tree.** Ten files matched neither `HEAD` nor the working tree, because
deployment here is **per-file and by hand**: files were pushed as each run needed them, and three
that were fixed tonight and never pushed —

| file | fixed locally | host copy |
|---|---|---|
| `bond.py` | 01:21 (`7f847af7`) | 01:14, 492 lines vs 540 |
| `step.py` | 03:58 (`93e84b78`) | 00:56, 478 lines vs 486 |
| `geometry.py` | 03:59 (`a1bdcd4d`) | 01:14, 292 lines vs 322 |

— were still the pre-fix versions when the 04:40 seven-population run executed.

### 2.1 ⚠ Correction — the three extra files were not stale copies, they were a FLATTENED directory

I first recorded them as *"superseded top-level copies"*. Session A checked its own name against the
repository and found the better answer: **`aleph/world/nmii_cortex_crossbridge.py` and
`aleph/world/membrane_erm_cortex.py` have never existed in this repository at all** — `git log --all
--follow` returns **zero commits** for either path, while `world/families/` versions of both have
three and two. The hand deployment did not leave old copies behind; **it wrote files to the wrong
place**, flattening a subdirectory into its parent.

That is a worse defect than the one it replaces, and in a specific way: a stale copy is at least a
file someone once wrote at that path. A flattened one is a path that **no commit can ever describe**,
sitting in the parent package where `import aleph.world.nmii_cortex_crossbridge` resolves. Nothing
imported them tonight. Nothing was stopping it.

## 3. The guard existed, was tested, and was wired to nothing

`aleph/scripts/run_provenance.py` was written for exactly this. Its module docstring says the
surviving half of the archived `ffn_gpu.py` is the part that checks *"that the tree on the run host
actually holds that code"*, and warns that **"a record whose commit does not describe what ran cannot
be corrected afterwards: the run is over and the number is in a file."** `remote_mismatches`' own
docstring names the mechanism outright: **"that tree is updated by hand."**

Callers of `import_closure` / `remote_mismatches`, in the whole repository:

* `aleph/archive/scripts/ffn_gpu.py` — **archived**
* `aleph/tests/coordination/test_ffn_gpu_build_stamp_guard.py` — **a test**

**Zero production drivers.** The check was documented, implemented, unit-tested, and never called. It
is not that the hazard was unknown; it is that knowing it was stored somewhere that does not run.

## 4. And the reused code had the same defect class in it

Wiring `stamp()` in, the new Sanity Gate immediately failed on a case that had been in the tree for
weeks: `_closure_digest([])` returned `sha256:e3b0c442…` — **the digest of nothing**, a perfectly
well-formed hex string. Every driver whose import closure failed to resolve would have stamped that
same value, and any two such records would compare **equal**: *"this is the same code"*, asserted from
having read none of it.

Fixed at the shared call site (`return ""` on an empty closure), which is where all callers route
through, rather than in the one function that tripped over it.

This is the fourth member of tonight's family — the zero-step acceptance pass, the empty-dict
placement pass, the constant-γ stationarity pass, and now the empty-closure digest. **The shape is
always the same: a check whose subject can be empty, and whose verdict does not depend on whether it
was.**

## 5. What changed

* `run_provenance.stamp(driver)` — new entry point that runs **on the run host, from inside the
  driver**, and returns `build_commit`, `closure_digest`, `closure_files`, `session`.
* Wired into `world_phase4_native.py` and `world_phase1_native.py`, at the top of the record.
* `closure_digest` is the load-bearing field, **not** `build_commit`. The host has no `.git`, so the
  commit falls back to `COMMIT_STAMP` and is *labelled* `(COMMIT_STAMP, hand-deployed)` — exactly the
  field that went stale. The digest is over the bytes that were imported, so it cannot.
* The host tree was re-synced: **48 / 48 `world/` files now identical to `bb4cc80f`**, the three
  superseded shadows removed, `COMMIT_STAMP` rewritten.

## 6. What this does NOT fix, and what it costs

⚠ **Tonight's numbers are not retroactively attributable.** The runs are over. The best available
statement for any pre-05:00 record is *"host tree as hashed at 04:50, which matches no commit"*, and
that is what they should carry — not a commit hash chosen afterwards because it is close.

⚠ **The two `STATE.md` rows still have no build commit**, and `check_state_md.py` still refuses them.
That refusal is correct and is being left in place rather than satisfied with a plausible hash. **The
fix is to re-run the measurement on the synced tree**, so the record stamps a digest that a local
commit can be matched against — not to write a hash that would make the checker quiet.

⚠ The commit is still hand-written on the host. `stamp()` labels it as such; it does not make it
true. A deploy that writes `COMMIT_STAMP` as part of the same action that copies the files is the
actual fix, and is **PI queue item 10**.


## 7. What a stale file does to a PREDICATE, which is worse

Session B's contribution, and it generalises past this incident. Everything above is about a
**number** being attributed to the wrong code. But `world/contact.py`'s `check_complementarity()` is
a **predicate**, and the failure mode is not the same shape:

> If the host holds a stale `contact.py`, the gate can **PASS against a different tolerance than the
> committed one**, and the artifact records only `PASS`.

A wrong number is at least visibly a number, and a reader can compare it to something. A `PASS`
carries no trace of what it was measured against, so a stale predicate is unfalsifiable after the
fact — there is nothing in the record to disagree with. **Any driver that first writes
`check_complementarity()`'s dict into a record must stamp beside it**, or which tolerance produced
the PASS is unrecoverable.

## 8. Whose numbers this actually touched

Asked directly, four sessions checked rather than assumed, and the answers were not uniform:

* **Lane A, Lane B, session 3b** — never submitted to the host at all. Nothing to re-run. Lane A's
  only quoted figure (26,807 characters of emitted CUDA C++) is a dev-Mac codegen result, re-checked
  at `c41d4e95` and unchanged; it is compile evidence, not a physics number.
* **Lane D** — clean on its own files, and then flagged the one thing I had not: the 12-population
  `n_outside = 0` run at 03:39 **is inside the window**, and its `sf_arc` portion depends on files
  that were two hours old at the time. Their words: *"the content is probably identical to the commit
  — and 'probably' is exactly what `closure_digest` replaces."* PHASE 1 was re-run stamped rather
  than quoted.

**Every session that could have been affected checked its own window before answering.** That is the
part of this worth keeping: the guard being absent was recoverable; a session assuming its numbers
were fine would not have been.


## 9. The re-run, and what it demonstrated on the way

PHASE 1 was re-run stamped on the synced tree. **Result: the flagged run is CONFIRMED, not retracted.**

    census live: node 4,591,262  segment 4,313,157  angle3 4,241,350
                 angle4 614,400  face 409,600  strand 2,443  grid_cell 250,047
    placement:   12 populations, n_outside = 0
    provenance:  closure_digest sha256:0742a2fa…, argv recorded

That census is **bit-identical** to `phase1_twelve.json` and `phase1_refactored.json`, the records the
03:39 run wrote. The code that ran during the stale window built the same cell the synced tree builds.

### 9.1 And the first attempt at the re-run demonstrated the argv gap live

The first re-run came back with **4,558,554** nodes and **11** populations, not 4,591,262 and 12. The
difference decomposed exactly:

    4,588,662 − 32,708 (nmii) + 2,600 (sf_arc) = 4,558,554        exact, to the node

For about a minute this looked like the stale-deployment finding landing a second time. It is not.
`--nmii-heads-per-side` is **UNRATIFIED** (PI queue item 2) and the driver builds NMII *only when the
flag is passed*; I had not passed it. **The difference was the invocation, and the old record could
not say so, because it did not store its invocation** — which is the exact gap closed in `763632d1`
an hour earlier, arriving to demonstrate itself.

Re-run with `--nmii-heads-per-side 30` — the value recovered from the old record's own `nmii` block,
which *was* stored — the census reproduces exactly.

**Two different failures wear the same face**: *the number changed, so the code must have been stale.*
One was real (§2), one was an unrecorded flag. Without `argv` in the record there is no way to tell
them apart, and the wrong one is the one you would act on.
