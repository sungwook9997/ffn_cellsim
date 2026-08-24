"""Live demo of the signature 5-hop provenance query — the capability RAG and
vanilla TAG cannot provide.

Builds a small but FAITHFUL contract-graph in a throwaway DuckDB (the real
schema: 8 node tables + an edges triple table, exactly as notion_to_duckdb.py
materializes it) and runs the signature query that walks
    run_result -> validation_gate -> model_contract -> parameter -> source_evidence
in a single SQL statement — "from this result, which paper does it ultimately
rest on, and is that citation verified?".

This is a reproducible, deterministic presentation artifact: it touches NO real
file (in-memory DuckDB) and writes a transcript to DEMO_5hop_query.md. It is the
"exec" core of the TAG loop (Fig 2); tag_query.py wraps the same SQL with an
LLM syn (NL->SQL) and gen (rows->answer+citations) step.

Run::

    conda activate ffn_sim
    python aleph/outputs/tag_kb/presentation/demo_5hop_query.py
"""
from __future__ import annotations

from pathlib import Path

import duckdb

OUT = Path(__file__).parent / "DEMO_5hop_query.md"


def build(con) -> None:
    """A faithful miniature of the live contract-graph (real column names)."""
    con.execute("CREATE TABLE source_evidence (id TEXT, citation_key TEXT, uid TEXT, "
                "doi TEXT, source_type TEXT, verdict TEXT)")
    con.executemany("INSERT INTO source_evidence VALUES (?,?,?,?,?,?)", [
        ("se_chugh", "Chugh2017_NatCellBiol", "SE201", "10.1038/ncb3525",
         "Direct measurement", "CHECK"),
        ("se_bell", "Bell1978_Science", "SE119", "10.1126/science.347575",
         "Model-derived", "OK"),
    ])
    con.execute("CREATE TABLE knowledge_claim (id TEXT, title TEXT, kb_id TEXT, status TEXT)")
    con.executemany("INSERT INTO knowledge_claim VALUES (?,?,?,?)", [
        ("kc35", "Cortex surface tension band", "KU-3.5", "accepted"),
    ])
    con.execute("CREATE TABLE model_contract (id TEXT, title TEXT, mc_id TEXT, status TEXT)")
    con.executemany("INSERT INTO model_contract VALUES (?,?,?,?)", [
        ("mc_cortex", "Cortical tension contract", "MC-U5-cortical-tension", "active"),
    ])
    con.execute("CREATE TABLE parameter (id TEXT, title TEXT, value TEXT, unit TEXT)")
    con.executemany("INSERT INTO parameter VALUES (?,?,?,?)", [
        ("p_gamma", "gamma_cortex", "0.5e-3", "N/m"),
    ])
    con.execute("CREATE TABLE validation_gate (id TEXT, title TEXT, vg_id TEXT, status TEXT)")
    con.executemany("INSERT INTO validation_gate VALUES (?,?,?,?)", [
        ("vg_a", "H7 Gate-A cortical tension band", "VG-H7-gate-a", "failing"),
    ])
    con.execute("CREATE TABLE run_result (id TEXT, title TEXT, run_id TEXT, commit TEXT, outcome TEXT)")
    con.executemany("INSERT INTO run_result VALUES (?,?,?,?,?)", [
        ("rr1", "H7 gamma production", "RUN-H7-001", "ac8add0", "FAIL: gamma ~1000x under band"),
    ])
    con.execute("CREATE TABLE code_mapping (id TEXT, title TEXT, path TEXT)")
    con.execute("CREATE TABLE decision_ledger (id TEXT, title TEXT)")
    con.execute("CREATE TABLE edges (src_id TEXT, src_type TEXT, rel TEXT, dst_id TEXT, dst_type TEXT)")
    con.executemany("INSERT INTO edges VALUES (?,?,?,?,?)", [
        ("rr1", "run_result", "tests", "vg_a", "validation_gate"),
        ("vg_a", "validation_gate", "validates", "mc_cortex", "model_contract"),
        ("mc_cortex", "model_contract", "constrains", "p_gamma", "parameter"),
        ("mc_cortex", "model_contract", "derived_from", "kc35", "knowledge_claim"),
        ("p_gamma", "parameter", "anchored_by", "se_chugh", "source_evidence"),
    ])


FIVE_HOP = """\
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
WHERE lower(g.status) = 'failing';"""


def md_table(cols, rows) -> str:
    head = "| " + " | ".join(cols) + " |"
    sep = "|" + "|".join(["---"] * len(cols)) + "|"
    body = "\n".join("| " + " | ".join(str(v) for v in r) + " |" for r in rows)
    return f"{head}\n{sep}\n{body}"


def main() -> None:
    con = duckdb.connect(":memory:")
    build(con)
    cur = con.execute(FIVE_HOP)
    cols = [d[0] for d in cur.description]
    rows = cur.fetchall()

    print(FIVE_HOP)
    print()
    for r in rows:
        print(dict(zip(cols, r)))

    doc = f"""# Live demo — the 5-hop provenance query

*Reproducible: `python aleph/outputs/tag_kb/presentation/demo_5hop_query.py`
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
{FIVE_HOP}
```

### The result

{md_table(cols, rows)}

### Reading it in plain language

Run **{rows[0][0]}** ({rows[0][1]}) failed validation gate **{rows[0][2]}**.
That gate tests contract **{rows[0][4]}**, which constrains the constant
**{rows[0][5]} = {rows[0][6]} {rows[0][7]}** — and that number is anchored to
**{rows[0][8]}**, whose citation audit verdict is **{rows[0][9]}** (i.e. *not*
confirmed-OK). So the talk's claim is concrete: from a red CI run we trace, in one
hop-chain, all the way back to the paper — and we can see that this particular
constant is riding a citation that still needs verification.

> In the full system the same SQL is what `tag_query.py` generates from the plain
> English question (the `syn` step) and then narrates with citations (the `gen`
> step). RAG would return "relevant-sounding" paragraphs; TAG returns *this row*.
"""
    OUT.write_text(doc)
    print(f"\nwrote {OUT}")
    con.close()


if __name__ == "__main__":
    main()
