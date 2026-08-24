# FINDING — eight checks in one day that returned green for a reason unrelated to what they claim

> **`AGENT-WRITTEN`.** No `decided_by`. This is a measurement about the tree's verification surface,
> not a decision about its physics. It is filed here rather than in `STATE.md` (c) because (c) is a
> list of **claims that may not be quoted**, and this is a **class of defect** — the difference
> matters: a (c) entry retires a number, and nothing here retires a number. What it does is tell you
> which numbers were never protected in the first place.

| Field | Value |
|---|---|
| Written | 2026-08-09 23:1x KST, session `5ccfd28d`, lane `contracts` |
| Instances | **eight**, found by three sessions, none of them looking for this |
| Asked for by | `97fb3916`: *"If your (c) framing has room for a pattern rather than an instance, this is it."* |

---

## The eight

| | the check | what it actually tested |
|---|---|---|
| 1 | `require_complete=True` on the composed world | Satisfied only ever by `_CensusRuntime` stubs and `np.random.default_rng(23)` geometry. **Every vertical slice runs `require_complete=False`.** The contract has never been met by real runtimes, and physics has only ever run under the relaxed one |
| 2 | `ECMWorld.crosslink_connector` — name, internal endpoints, adjoint requirement, all validated | **Nothing satisfied the Protocol.** The tier-(a) row's *"zero crosslink bonds"* was a missing implementation, not a scoping choice. → `STATE.md` (c) 18 |
| 3 | `mechanical_joint_count` | Defaults to 1 **with no producer**. A census reporting a number nobody computes |
| 4 | `resolve_fa_series_group` | Keys rivalry on the component pair, so distinct adhesions collapse into one group |
| 5 | `assert_component_state_disjoint` | Compares device pointers for **exact equality**, so two *overlapping slice views* of one buffer pass the disjointness check |
| 6 | `ownership.yaml`'s overlap rule | The loader refuses two lanes declaring the **same glob string**. The property wanted is *no two lanes cover the same path*. `aleph/engine/co[^r]*.py` loaded cleanly and excluded nothing |
| 8 | `--declared-commit` on `ac_resting_residual_curve.py` | The flag parses, and `observation_artifact` takes `declared_commit` and forwards it — **the call does not pass it.** The record then reports `"no commit was declared"`, which reads as caller error. A parameter that exists at both ends and is unconnected in the middle. **Introduced by making the record schema-compliant**: the re-emission moved the stamp to a path with its own parameter, and the caller did not know to fill it |
| 7 | `residual_curve.json` declaring `record: "run-record@2"` | **It does not satisfy that schema.** `t0` and `timing` are REQUIRED — the docstring says so twice and gives the reason — and the record carries neither, nor `build`, `config_sha256` or `evidence`. A compliant one next to it carries all of them. Self-declaration stands in for conformance, with nothing checking |

Instances 1–2 and 5 are `97fb3916`'s and the Card-5 session's; 2 landed in (c) this evening. 6 and 7
are this session's, and **7 was found while trying to write the (b) row that 3 of these instances
argue for** — which is the pattern demonstrating itself: the check that would have caught it is
(b)'s own rule, *no row without a build commit*, and the artifact's self-declared schema is what a
reader would otherwise have trusted instead of checking.

## The smallest one, first, because it needs no code to understand

**A test that is not collected is a check returning green for a reason unrelated to what it claims.**

`pyproject.toml`'s `testpaths` was an allowlist of subdirectories, and it had gone stale: sixteen
files sat directly under `aleph/tests/` and **none of them was collected** — 120 test functions, 45
in `test_sanity_gate.py` alone. Four of the sixteen were CLI gate scripts wearing a `test_` name,
with a module-level `sys.exit` that would have aborted the whole pytest session the moment anything
collected them. The suite was green. It had been green for as long as the list had been stale.

There is no smaller version of the pattern than that, and it is the one to hold in mind for the
other six: **the check ran, the check passed, and what it examined was not what its name says.**

## The shape, stated once

**Each check names a property, and tests a proxy that is true more often than the property is.**

- exact pointer equality is a proxy for disjointness — and it is *implied by* disjointness, not
  equivalent to it, so every overlapping-but-not-identical case passes;
- a Protocol's field being declared and validated is a proxy for the Protocol being implemented;
- a default is a proxy for a value;
- a component pair is a proxy for an adhesion's identity.

None of these is a careless test. Each is what you write when you know what you mean and reach for
the nearest thing you can compare. **The failure is not laxity, it is that the cheap comparison
agrees with the expensive one on every case anybody happened to try.**

## Why eight in one day, and why not before

They surfaced because three things happened at once that had not before: the composed world was
driven with real runtimes rather than stubs (`97fb3916`), a session went looking for the shared-array
aliasing that `STATE.md` (c) 4 describes (Card-5), and **collection was widened from an allowlist of
directories to the whole test tree**, which added 410 tests that had never run.

The third is the `testpaths` widening described above — which is why the surfacing is evidence about
**how these are found**, not about how many remain. Each of the seven was turned up by somebody doing
something else: driving the world with real runtimes, hunting the shared-array aliasing, splitting a
lane, rendering a row.

## What this does not license

**It does not retire any result.** No number moves because of this document. What it licenses is a
*prior*: when a check in this tree reports satisfaction, the question *"satisfied by what?"* has
returned a surprising answer eight times in one day, and it is cheap to ask.

It also does not license a sweep. Auditing every predicate in the tree for proxy-vs-property would
be a large piece of work with no bound and no gate, and the entries above were each found by
somebody doing something else. **The right response is the prior, not a programme.**

## What was actually done about it

Nothing structural, deliberately. Four of the eight are `engine`'s to fix, one is `card5`'s, and
two are this session's — 6 is fixed (`tests/architecture/test_ownership_lanes.py` is the check the
loader cannot do) and 7 is reported below rather than patched, because a record's schema is not a
thing a reader should repair on the record's behalf. This
session's contribution is the record, `STATE.md` (c) 18 for the one that reached a quotable row, and
two guards written today that are the same shape as the defect and were written to be resistant to
it:

- `tests/architecture/test_layer_directions.py` — its **vacuity control** asserts the permitted
  direction is non-empty, because a scanner that matched nothing would pass every forbidden-direction
  test while checking nothing;
- `tests/ports/test_port_discipline.py` — `_defined_test_functions` parses rather than greps, after
  the grep version counted a **filename in a docstring** as a test definition and let a ledger entry
  name a control nobody had written.

Both of those were caught by their own vacuity controls before they shipped. That is the only
defence in this document that generalises: **every check gets one test that fails if the check stops
checking.**
