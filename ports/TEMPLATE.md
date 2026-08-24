# ALEPH-PORT-NNN — <one-line name of the thing being ported>

<!--
  PORT LEDGER ENTRY TEMPLATE.

  Copy this file to ports/ledger/<ID>-<slug>.md and fill it in BEFORE the code lands.
  Every heading below is mandatory. `tests/ports/test_port_discipline.py` checks that each
  one is present and that none of them still contains a placeholder (TODO / TBD / FIXME / ???).

  Two rules that decide whether this entry is worth writing at all:

    * The unit of approval is a function, kernel, dataclass, or small cohesive module.
      A directory is never the unit of approval. If you cannot name the symbol, stop.
    * Passing the source's tests is not evidence. Matching the source numerically is not
      evidence. This entry must stand on a control this tree owns.

  Delete this comment block when you fill the template in.
-->

| Field | Value |
|---|---|
| Port ID | `ALEPH-PORT-NNN` |
| Lane | `LN <name>` |
| Status | `PROPOSED` \| `AUDITED` \| `ACCEPTED` \| `REJECTED` |
| Written | `YYYY-MM-DD` (before the code) |
| Port class | `RE-DERIVED` \| `SOURCE-DERIVED` \| `CONCEPT_ONLY` \| `DEFECT_STUDY` |

---

## 1. API this entry authorises

The exact public surface. Module path, symbol names, signatures. Nothing outside this list is
covered by this entry.

```python
from aleph.<subpackage>.<module> import (
    ...,
)
```

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/Project_Aleph` (READ ONLY to this tree) |
| Source commit | `<full 40-char sha>` (`<short>`) |
| Source path | `aleph/.../<file>.py` |
| Source symbol(s) | `<function / class / kernel names actually read>` |
| Read from | `git show <sha>:<path>` \| working tree — **state which** |
| Working tree == commit? | `yes` / `no` — verified with `git diff <sha> -- <path>` |

> The last two rows are not bureaucracy, and this tree is the reason the rule exists. At
> `be0e5876` its own working tree carried thirty-one uncommitted changes. A ledger entry that
> cites a commit but was written from the working tree is citing a revision it never read.
> `Project_Aleph` has a live session in it, so this check is not a formality there either.

## 3. Why source-derived porting beats clean-room

State the specific thing that would be expensive or unlikely to re-derive from nothing. A
convention any competent author would arrive at independently is **not** a reason to port — it is a
reason to write it clean-room and cite the source as prior art, if at all.

Legitimate: an empirical failure mode discovered by a real run; a branch structure whose case
analysis is easy to get subtly wrong; a numerically conditioned form.

Illegitimate: "it already works", "it is well tested there", "it saves time".

If clean-room wins, say so, set `Status: REJECTED`, and keep the entry.

## 4. Physical or mathematical law represented

The law, restated in this tree's vocabulary, derived here rather than quoted. If there is no law —
if what crosses is a convention, a schema, or an argument — say so and name what it actually is.

## 5. Units, domains, singular cases, invariants

| Quantity | Unit | SI unit | Domain |
|---|---|---|---|

Singular and boundary cases, each with the required behaviour:

- ...

Invariants, each with the test that asserts it:

- **I1.** ...

## 6. Source evidence class and known retractions

What the source's artifacts claim for this code and at what evidence class. Any retraction,
supersession, defect record, or self-labelled `BLOCKED` / `UNVERIFIED` status. Whether the symbol is
reachable from the source's production path or is dead code. Whether it has live tests there and —
separately — whether those tests exercise the law or only the plumbing.

"No retraction found" is acceptable only if you looked; say where.

## 7. Independent oracle or derivation

The thing this tree checks the port against, which is **not** the source. A closed form, an
analytically integrable case, a conservation law, a symmetry, an exact identity, a convergence rate
with a known exponent, or an independent implementation written from §4 with the source closed.

If the only available check is agreement with `Project_Aleph`, this port cannot be accepted.

## 8. Positive control

A named test that must pass, and the quantitative statement it makes. Not "the function runs".

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/.../test_x.py::test_...` | ... |

## 9. Deliberately failing negative control

A named test asserting that a **wrong** input, parameter, or model is rejected — and that fails if
the rejection stops working. A negative control that cannot fail is decoration.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/.../test_x.py::test_...` | ... |

## 10. Numerical and precision envelope

Working precision; accumulation precision; the tolerance each control is asserted at and why that
number and not a looser one; conditioning as written; the input range over which the tolerance
holds; what happens outside it (refusal, never silent degradation).

## 11. Production-backend residency and transfer

Where the data lives in production (host / device), what transfers and when, what precision each
copy carries, whether any host round-trip is required per step. **This tree forbids authoritative
per-step host state**, so a port that introduces one is refused here rather than argued. If the code
never runs on the production backend (an oracle, a schema, a diagnostic), say so — that is a
residency answer too, and it is what keeps `validation/**` out of the runtime.

## 12. Comments and docstrings to discard

Source prose that must not survive: foreign vocabulary, internal document references, gate names,
branch names, status labels, working-tree assumptions, dead-module citations, and any claim the
source asserts without evidence. List them so a reviewer can check they are gone, and state what
replaces them.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | `<what passed, at what tolerance, on what date>` |
| Reviewer | `<who reviewed; agent-proposed entries say so and remain unratified>` |
| Rollback | `<the exact change that removes this port and what breaks if it does>` |

## 14. Honest limits

What this entry does **not** establish. Anything carried at `UNVERIFIED`, `ASSUMED`, or
`INHERITED_UNVERIFIED`. Anything believed but not tested. Anything needing a GPU, a PI decision, or
a literature reading that has not happened.

An entry with an empty §14 is an entry whose author did not look hard enough.
