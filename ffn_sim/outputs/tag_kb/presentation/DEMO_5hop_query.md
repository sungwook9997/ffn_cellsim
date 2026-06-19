# Live demo — the 5-hop provenance query

*Reproducible: `python ffn_sim/outputs/tag_kb/presentation/demo_5hop_query.py`
(in-memory DuckDB, touches no real file).*

**The question a wet-lab reviewer actually asks:** *"This run failed its
validation gate — what physical constant is it testing, and which paper does
that number ultimately come from? And is that citation even verified?"*

A keyword/RAG search cannot answer this: the answer is not in any single chunk of
text, it lives in the **relationships** between a run, a gate, a contract, a
parameter, and a source. The contract-graph answers it with **one** query that
walks five hops:

```
run_result → validation_gate → model_contract → parameter → source_evidence
```

### The query (the `exec` core of the TAG loop)

```sql
SELECT r.run_id, r.outcome, g.vg_id, g.status AS gate, mc.mc_id,
       p.title AS parameter, p.value, p.unit,
       se.citation_key, se.verdict AS citation_verdict
FROM run_result r
JOIN edges e1 ON e1.src_id = r.id  AND e1.dst_type = 'validation_gate'
JOIN validation_gate g ON g.id = e1.dst_id
JOIN edges e2 ON e2.src_id = g.id  AND e2.dst_type = 'model_contract'
JOIN model_contract mc ON mc.id = e2.dst_id
JOIN edges e3 ON e3.src_id = mc.id AND e3.dst_type = 'parameter'
JOIN parameter p ON p.id = e3.dst_id
JOIN edges e4 ON e4.src_id = p.id  AND e4.dst_type = 'source_evidence'
JOIN source_evidence se ON se.id = e4.dst_id
WHERE lower(g.status) = 'failing';
```

### The result

| run_id | outcome | vg_id | gate | mc_id | parameter | value | unit | citation_key | citation_verdict |
|---|---|---|---|---|---|---|---|---|---|
| RUN-H7-001 | FAIL: gamma ~1000x under band | VG-H7-gate-a | failing | MC-U5-cortical-tension | gamma_cortex | 0.5e-3 | N/m | Chugh2017_NatCellBiol | CHECK |

### Reading it in plain language

Run **RUN-H7-001** (FAIL: gamma ~1000x under band) failed validation gate **VG-H7-gate-a**.
That gate tests contract **MC-U5-cortical-tension**, which constrains the constant
**gamma_cortex = 0.5e-3 N/m** — and that number is anchored to
**Chugh2017_NatCellBiol**, whose citation audit verdict is **CHECK** (i.e. *not*
confirmed-OK). So the talk's claim is concrete: from a red CI run we trace, in one
hop-chain, all the way back to the paper — and we can see that this particular
constant is riding a citation that still needs verification.

> In the full system the same SQL is what `tag_query.py` generates from the plain
> English question (the `syn` step) and then narrates with citations (the `gen`
> step). RAG would return "relevant-sounding" paragraphs; TAG returns *this row*.
