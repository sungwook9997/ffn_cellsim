# The distinction the balance-gate finding established was already written down — in a frozen layer

**Status:** FINDING. Read from source, not measured. Nothing here changes a result; it changes where
(e) 1's wording should start.

---

## 1. What tonight established, at length

`BALANCE_GATE_MEASURED_2026-08-21.md`: `ledger.py`'s acceptance predicate tests **adjoint closure** —
that every force appears with both signs, to within the Higham summation bound — and this is a
DIFFERENT property from convergence. It cannot fail for physics, and it accepted a run whose residual
grew four orders of magnitude while its own tolerance grew 86×.

## 2. ⚠ `aleph/virtual_cell/admissibility.py` has said so since before this engine existed

```python
ADMISSIBILITY_CLAUSE_ORDER = (
    AdmissibilityClause.STATE_FINITE,
    AdmissibilityClause.ADJOINT_CLOSURE,     # <- separate
    AdmissibilityClause.INNER_CONVERGED,     # <- separate
    AdmissibilityClause.ITERATION_BUDGET,
)
```

**Four clauses, and adjoint closure is not convergence.** More than that — the same module routes a
failure of each to a different KIND of work:

| clause | where a failure routes |
|---|---|
| `ADJOINT_CLOSURE` | *"**archetype** failure → missing component / missing law project — a one-sided or sign-flipped scatter means a connector law is missing or miswired"* |
| `INNER_CONVERGED` | *"**native cost** failure → solver / preconditioner / adaptive design"* |

**An architecture problem and a numerics problem.** Collapsing them into one gate does not only lose a
distinction — **it routes every failure to the wrong desk.** A run that fails to converge gets read as
a wiring bug; a genuinely miswired connector gets read as a solver that needs a bigger budget. The
2026-08-15 result — *"more iterations make it worse"* — is exactly a case that was routed as a budget
problem and was not one.

## 3. Nothing outside `virtual_cell` computes either clause

Searched: the only references outside that package are its own tests. `aleph/engine/`, `aleph/world/`,
`aleph/laws/` and `aleph/components/` contain neither name.

So the vocabulary exists, is correct, has a per-clause routing table — and **is imported by nothing.**
`aleph/virtual_cell/` is `archive-readonly` in `ownership.yaml` and 30,024 lines, and `STATE.md` (a)
records it as *"동결(archive-readonly), 배선 0"*.

## 4. ⚠ And `world/` may not import it

`tests/architecture/test_layer_directions.py` forbids exactly this edge:

```python
("world", ("engine", "components", "dcm", "virtual_cell", "scripts", "validation", "archive"))
```

So `world/step.py` **cannot** use `AdmissibilityClause`, and tonight's `AcceptancePredicate` had to
restate a vocabulary that already existed one layer over, without being able to reference it.

**This is the shape `units/provenance.py` already solved once.** That file was relocated out of
`engine/forces_manifest.py` so `world/` could read `PROVENANCE` without importing the port source —
and the reason recorded there was that *a vocabulary at the bottom of the layer stack that everything
reads is exactly what one lane must not be able to change alone.*

**The clause vocabulary is the same kind of thing and is in the wrong place.**

## 5. What this means for (e) 1

The amendment does not need to invent a clause structure. It needs to **use one that already exists
and is better than what the engine has**, and the concrete blocker is a layer rule rather than a
design question.

**For the PI, one decision:** relocate the admissibility clause vocabulary — the enum and its routing
table, not the oracle — to a layer `world/` may read, as `units/provenance.py` was. That is a move,
not a rewrite, and it is the same move for the same reason.

⚠ **Not done here.** `aleph/virtual_cell/**` is `archive-readonly` — *"Touch only to report, never to
fix"* — and relocating a vocabulary out of it is a governance change, not a session's call. This is
the report.

## 6. What is NOT claimed

That the four clauses are the right four, or that `INNER_CONVERGED` as `virtual_cell` defines it is
what (e) 1 should adopt. **Neither has been read against tonight's measurements**, and one of them —
`INNER_CONVERGED` — is precisely the thing (e) 1 says is undecidable, so it cannot simply be inherited.
What is claimed is narrower and sturdier: **the separation is already recorded, with a routing
consequence, in a place the engine cannot reach.**
