# ALEPH-PORT-4001 — a band without a citation is somebody's memory of a band

| Field | Value |
|---|---|
| Port ID | `ALEPH-PORT-4001` |
| Lane | `LN merge` (session `5ccfd28d`) |
| Status | `PROPOSED` |
| Written | `2026-08-09` — before the code, and the code is written clean-room per §3 |
| Port class | `CONCEPT_ONLY` |

---

## 1. API this entry authorises

```python
from aleph.units.band import (          # `aleph.units.band` after the rename
    LiteratureBand,
    OutOfBandError,
    band_guard,
)
```

Nothing else. In particular this entry does **not** authorise a dimensional-quantity type, a unit
algebra, or a conversion table — this tree already has a unit system (`ff/units.py`, the pN·µm·s
nondimensionalisation) and a second one would be a second answer to "what are these numbers".

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/Project_Aleph` (READ ONLY to this tree) |
| Source commit | `9ef077ba72200a725e01930332646d09900b63f1` (`9ef077ba`) |
| Source path | `aleph/units/guard.py` |
| Source symbol(s) | `LiteratureRange`, `RangeGuardFailure`, `guard_range`, `guard_all_in_range`, `guard_positive` |
| Read from | working tree — **stated explicitly** |
| Working tree == commit? | `yes` — `git status --short aleph/units/` reported nothing at the time of reading |

`Project_Aleph` has an open session (`b4fad06b`) but it owns `aleph/vertical/**` and
`aleph/scenarios/**`, not `aleph/units/**`, and the path was clean when read.

## 3. Why source-derived porting beats clean-room — **it does not, and that is the finding**

The template asks for the specific thing that would be expensive to re-derive. For this symbol the
honest answer is: **nothing.** A dataclass holding `lo`, `hi`, a unit string and a citation, with a
method that raises when a value falls outside, is what any competent author writes when asked for a
range check. Porting it would be the illegitimate reason the template names — *"it already works"*.

**What is not obvious, and is the entire asset, is a pair of design commitments:**

1. **The citation is a required constructor argument, not an optional one.** A bound with no source
   is somebody's memory of a bound, and it cannot be argued with — when the guard fires, the reader
   has no way to decide whether the band is wrong or the run is.
2. **The guard raises and returns the value untouched; it never clamps.** A quantity leaving its
   supported band is a finding. Clamping destroys the finding and leaves a run that looks healthy.

Those two are contracts, so they travel as contracts. **The code is written here from this entry,
with the source closed.** Port class is therefore `CONCEPT_ONLY` and the source is prior art.

This tree has been getting the second commitment wrong in the other direction: `STATE.md` (c) 2 and
`ff/ENGINE.md` describe headline numbers that were produced, quoted, and only later found to be
outside what supported them. A guard that fires at the point of use is the mechanism whose absence
that record documents.

## 4. Physical or mathematical law represented

**None. This is not a law and it must not be mistaken for one.** What crosses the boundary is an
argument about provenance: a numeric acceptance band is a triple `(interval, unit, source)`, and the
first two without the third are not checkable by a reader.

The predicate itself is set membership on a closed or open interval of the extended reals, with
non-finite values outside the band by definition. That is not physics; it is the definition of the
band.

## 5. Units, domains, singular cases, invariants

| Quantity | Unit | SI unit | Domain |
|---|---|---|---|
| `lo`, `hi` | declared per band, as a string | whatever the band declares | `lo <= hi`, both finite |
| the guarded value | must match the band's declared unit — **by the caller's discipline, not by a type** | — | any float |

The unit is a **string carried for the error message**, not a checked dimension. Stating that
plainly is part of the port: a reader who believes the unit is enforced will pass a value in the
wrong unit and get a green guard. Enforcing it needs a dimensional quantity type, which §1 excludes.

Singular and boundary cases:

- `NaN` — **outside every band**, including one spanning `(-inf, inf)`. It is not a small error.
- `+inf` / `-inf` — outside every band with a finite bound on that side.
- `lo == hi` — a legal degenerate band; only `inclusive=True` can then contain anything.
- `lo > hi` — **malformed**, refused at construction, not at first use.
- an empty or whitespace-only citation — **refused at construction.**

Invariants:

- **I1.** A guard that passes returns its input unchanged, `is`-identical for arrays. Asserted by
  `test_guard_returns_the_value_untouched`.
- **I2.** A guard never clamps. There is no code path that returns a boundary value in place of an
  out-of-band one. Asserted by `test_out_of_band_raises_and_does_not_clamp`.
- **I3.** A band cannot exist without a citation. Asserted by `test_a_band_without_a_citation_is_refused`.
- **I4.** `NaN` is outside every band. Asserted by `test_nan_is_outside_even_an_infinite_band`.

## 6. Source evidence class and known retractions

The source symbol carries no physical claim, so it has no evidence class to inherit; it is
infrastructure. Where this was looked for:

- `Project_Aleph/ports/ledger/` — grep for `units/guard`: no entry, i.e. the symbol was authored in
  Aleph rather than ported into it, so there is no upstream ledger to inherit.
- `Project_Aleph/docs/decisions/` — no `PROPOSAL-*` or `ALEPH-PD-*` names it.
- `Project_Aleph/README.md` and `PLAN.md` — not named.

It is reachable from Aleph's package and has live tests there. **Those tests are not evidence for
this port** and were not read; §7 is.

## 7. Independent oracle or derivation

The predicate is set membership, so the oracle is the definition rather than a reference
implementation: for a closed band, `contains(x) ⟺ lo ≤ x ≤ hi ∧ isfinite(x)`. The controls assert
that identity directly on a grid that straddles both bounds, plus the four singular cases in §5,
**without importing anything from `Project_Aleph`.**

There is no appeal to agreement with the source anywhere in the controls. If this port matched the
source and disagreed with §5, §5 wins.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `aleph/tests/units/test_band.py::test_membership_matches_the_definition` | For a band `[1e-2, 1.0]` and 2,001 points spanning `[-1, 2]`, `contains(x)` agrees with `1e-2 <= x <= 1.0` at **every** point, and `guard` raises on exactly the complement |
| Positive | `…::test_guard_returns_the_value_untouched` | I1 — the scalar guard returns an equal float and the array guard returns the same object |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `…::test_out_of_band_raises_and_does_not_clamp` | I2 — a value one ULP above `hi` raises `OutOfBandError`, and the exception message carries the citation. Fails if the guard is ever changed to clamp |
| Negative (must fail) | `…::test_a_band_without_a_citation_is_refused` | I3 — `LiteratureBand(..., citation="  ")` raises at construction |
| Negative (must fail) | `…::test_nan_is_outside_even_an_infinite_band` | I4 — the case a naive `lo <= x <= hi` silently passes as `False` but a naive `not (x < lo or x > hi)` silently passes as `True` |
| Vacuity | `…::test_the_positive_control_would_notice_a_broken_band` | A deliberately widened band makes `test_membership_matches_the_definition`'s comparison disagree, proving the comparison is load-bearing |

## 10. Numerical and precision envelope

Working precision is Python `float` (IEEE-754 binary64) throughout; there is no accumulation, so no
accumulation precision. The comparisons are exact — `<=` on binary64 — so the controls assert
**exact** agreement rather than a tolerance, and a tolerance here would be a way of not noticing a
boundary bug. The one-ULP negative control uses `math.nextafter(hi, inf)`, which is the smallest
input that must be refused; anything looser would pass a guard that had been widened by rounding.

Arrays are `np.float64`; a non-float array is converted by `np.asarray(..., dtype=float)`, and an
input that cannot convert raises rather than being coerced to `NaN`.

Outside the stated domain the behaviour is refusal, never silent degradation — that is I2.

## 11. Production-backend residency and transfer

**This code never runs on the production backend.** It is a construction-time and analysis-time
guard: bands are checked when a value is declared, when a run's inputs are assembled, and when a
result is read back between accepted steps. It allocates no device memory, launches no kernel, and
introduces **no authoritative per-step host state**, which is the condition this tree's charter sets.

A guard inside the inner mechanical loop would be exactly that forbidden thing, so `band_guard` is
not to be called from a kernel or from any function reachable from one. Nothing enforces that yet;
it is named in §14.

## 12. Comments and docstrings to discard

Nothing is copied, so nothing has to be stripped — that is the practical benefit of §3's answer.
The prose here is written against this tree's own record: `STATE.md` (c) 2 and `ff/ENGINE.md`'s
withdrawn headlines are the failures the contract exists to prevent, and they are this tree's
failures, not the source's.

Specifically **not** carried across: the source's module name (`guard.py` — too general for a module
that guards exactly one thing), its `RangeGuardFailure` / `guard_range` spellings, its `AlephUnitsError`
base class, and its dependency on a `dimension` module this entry excludes.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** `PROPOSED` until the controls in §8/§9 run green in a full-suite run |
| Reviewer | Agent-proposed (session `5ccfd28d`). Unratified |
| Rollback | Delete `aleph/units/` and `aleph/tests/units/`. Nothing imports it yet, so nothing breaks — which is also the reason it is cheap to accept and cheap to be wrong about |

## 14. Honest limits

1. **The unit string is not enforced.** §5 says so; a caller passing nm into a µm band gets a green
   guard. Fixing that needs a dimensional quantity type, which is a separate entry that does not
   exist.
2. **Nothing calls this yet.** Its value is entirely prospective until the first band is declared
   with a real citation, and a guard nobody calls has the same effect as no guard.
3. **No mechanism stops it being called from the inner loop**, where it would violate §11. Named,
   not solved.
4. **The citation is a free-text string, not a KB reference.** This tree has a Notion contract-graph
   with `source_audit.verdict` per source, and a band ought to cite a row in it rather than prose.
   That coupling is deliberately not built here: it would make `units` depend on `kb`, and the layer
   order in `STRUCTURE.md` §2 puts `units` at the bottom. The right shape is probably a separate
   `CitedBand` in `kb/`, and it is not designed.
5. **`inclusive=False` is implemented but has no caller and no literature motivation.** It is
   carried because a half-open band is the natural spelling for a strict physical bound, and it is
   tested, but it is unexercised in practice.
